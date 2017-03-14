#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <atomic>

#include "md.h"

#include "coords.h"
#include "coords_io.h"
#include "scon_chrono.h"
#include "scon_vect.h"
#include "scon_utility.h"


// Logging Functions

std::ostream& md::operator<<(std::ostream &strm, trace_data const &d)
{
  strm << std::right << std::setw(10) << d.i;
  strm << std::right << std::setw(20) << std::fixed << std::setprecision(3) << d.T;
  strm << std::right << std::setw(20) << std::fixed << std::setprecision(5) << d.P;
  strm << std::right << std::setw(20) << std::fixed << std::setprecision(5) << d.Ek;
  strm << std::right << std::setw(20) << std::fixed << std::setprecision(5) << d.Ep;
  strm << std::right << std::setw(20) << std::fixed << std::setprecision(5) << d.Ek + d.Ep;
  if (d.snapshot > 0u) strm << std::right << std::setw(20) << d.snapshot;
  else strm << std::right << std::setw(20) << '-';
  for (auto iae : d.Eia)
  {
    strm << std::right << std::setw(20) << std::fixed << std::setprecision(5) << iae;
  }
  strm << '\n';
  return strm;
}

void md::trace_writer::operator() (md::trace_data const & d)
{
	static std::atomic<bool> x(true);
	if (x && Config::get().general.verbosity > 1)
  {
    *strm << std::right << std::setw(10u) << "It";
    *strm << std::right << std::setw(20u) << "T";
    *strm << std::right << std::setw(20u) << "P";
    *strm << std::right << std::setw(20u) << "E_kin";
    *strm << std::right << std::setw(20u) << "E_pot";
    *strm << std::right << std::setw(20u) << "E_tot";
    *strm << std::right << std::setw(20u) << "Snapsh.";
    if (d.Eia.size() > 1u)
    {
      for (auto i : scon::index_range(d.Eia))
      {
        *strm << std::right << std::setw(20u) << "E_ia(# " << i << ')';
      }
    }
    *strm << '\n';
    x = false;
  }
  *strm << d;
}

namespace
{
  inline std::size_t gap(std::size_t const st, std::size_t const sn)
  {
    if (sn == 0u) return 0u;
    return st / std::min(st, sn); 
  }
}


md::Logger::Logger(coords::Coordinates &coords, std::size_t snap_offset) :
  snap_buffer(coords::make_buffered_cartesian_log(coords, "_MD_SNAP",
    Config::get().md.max_snap_buffer, snap_offset, Config::get().md.optimize_snapshots)),
  data_buffer(scon::offset_call_buffer<trace_data>(50u, Config::get().md.trackoffset, 
    trace_writer{ coords::output::filename("_MD_TRACE", ".txt").c_str() })),
  snapnum()
{
}

bool md::Logger::operator()(std::size_t const iter, coords::float_type const T,
  coords::float_type const P, coords::float_type const Ek, coords::float_type const Ep,
  std::vector<coords::float_type> const Eia, coords::Representation_3D const & x)
{
  if (iter % 5000u == 0)
  {
    if (std::isnan(T) || std::isnan(P) || std::isnan(Ek) || std::isnan(Ep))
    {
      std::cout << "NaN in simulation, please check your input options" << std::endl;
      throw std::logic_error("NaN in simulation, please check your input options");
    }
  }
  
  return data_buffer(trace_data(Eia, T, Ek, Ep, P, iter, snap_buffer(x) ? ++snapnum : 0u));
}

// Serialization Helper Function

static inline void append_to_buffer(std::vector<char> & buffer, void const * data, std::size_t const bytes)
{
  if (bytes > 0 && data)
  {
    buffer.resize(buffer.size() + bytes);
    std::vector<char>::iterator ins_pos_it(buffer.end() - bytes);
    memcpy(&(*ins_pos_it), data, bytes);
  }
}

static inline std::istream::pos_type istream_size_left(std::istream & S)
{
  std::istream::pos_type const T(S.tellg());
  S.seekg(0, std::istream::end);
  std::istream::pos_type const R(S.tellg());
  S.seekg(T, std::istream::beg);
  return R;
}


// Traces 


md::simulation::simulation(coords::Coordinates& coord_object) :
  coordobj(coord_object), 
  logging(coord_object, gap(Config::get().md.num_steps, Config::get().md.num_snapShots)),
  P(coord_object.xyz()), P_old(coord_object.xyz()),
  F(coord_object.g_xyz()), F_old(coord_object.g_xyz()),
  V(coord_object.xyz().size()), M(coord_object.xyz().size()),
  M_total(0.0), E_kin(0.0), T(Config::get().md.T_init), temp(0.0), dt(Config::get().md.timeStep),
  freedom(0), snapGap(0), C_geo(), C_mass(), press(0.0),
  nht(), rattle_bonds(), window(), restarted(true)
{
  std::sort(Config::set().md.heat_steps.begin(), Config::set().md.heat_steps.end());
}





void md::simulation::run(bool const restart)
{
  // Print initial info if required due to verbosity
  if (Config::get().general.verbosity > 0U) print_init_info();
  // Get first energy & derivatives
  coordobj.g();
  // get restarted flag
  restarted = restart;
  std::size_t iteration(0U);
  // (Re)initialize all values if simulation is to be fresh
  if (restarted)
  {
    nht = nose_hoover();
    T = Config::get().md.T_init;
    init();
	// remove rotation and translation of the molecule (only if no biased potential is applied)
	if (Config::get().md.set_active_center == 0) 
	{
		tune_momentum();
	}
  }

  // Read simulation data from resume file if required
  if (Config::get().md.resume)
  {
    std::ifstream restart_stream(
      (Config::get().general.outputFilename + "_MD_restart.cbf").c_str(),
      std::ifstream::in | std::ifstream::binary
      );
    scon::binary_stream<std::istream> bs{ restart_stream };
    if (bs >> iteration) bs >> *this;
  }
  // start Simulation
  integrate(false, iteration);
}

// Perform Umbrella Sampling run if requested
void md::simulation::umbrella_run(bool const restart) {

  // Generate umbrella sampling output file
  std::ofstream ofs;
  ofs.open("umbrella.txt", std::ofstream::out);
  std::size_t steps;
  steps = Config::get().md.num_steps;

  //General md initialization and config
  if (Config::get().general.verbosity > 0U)  print_init_info();
  coordobj.g();
  restarted = restart;
  if (restarted)
  {
    nht = nose_hoover();
    T = Config::get().md.T_init;
    init();
    // use nvect3d function to eliminate translation and rotation
    //V.eliminate_trans_rot(coordobj.xyz(), M);
    // use built in md function (identical)
    tune_momentum();
  }
  // Set kinetic Energy
  updateEkin();
  //run equilibration
  Config::set().md.num_steps = Config::get().md.usequil;
  integrate(false);
  // clear output vector and start production simulation
  udatacontainer.clear();
  // run production
  Config::set().md.num_steps = steps;
  integrate(false);
  //write output
  for (std::size_t i = 0; i < udatacontainer.size(); i++) {
    if (i% Config::get().md.usoffset == 0) {
      ofs << i << "   " << udatacontainer[i] << std::endl;
    }
  }
  ofs.close();
}


// Perform MD with the requested integrator
void md::simulation::integrate(bool fep, std::size_t const k_init)
{
  switch (Config::get().md.integrator)
  {
    case config::md_conf::integrators::BEEMAN:
      { // Beeman integrator
        beemanintegrator(fep, k_init);
        break;
      }
    default:
      { // Velocity verlet integrator
        velocity_verlet(fep, k_init);
      }
  }
}


// print initial MD conditions before simulation starts
void md::simulation::print_init_info(void)
{
  double const psec(Config::get().md.timeStep*static_cast<double>(Config::get().md.num_steps));
  std::cout << "Molecular dynamics simulation, limited to " << Config::get().md.num_steps;
  std::cout << " steps (representing a time step of " << Config::get().md.timeStep;
  std::cout << " ps each) started, simulating ";
  if (psec >= 500000.0) std::cout << (psec / 1000000.0) << " microseconds.\n";
  else if (psec >= 500) std::cout << (psec / 1000.0) << " nanoseconds.\n";
  else std::cout << psec << " picoseconds.\n";
  if (!Config::get().md.heat_steps.empty())
  {
    std::cout << "Temperature changes: ";
    std::cout << "[Step 0]:[" << Config::get().md.T_init << "K]";
    for (auto heatstep : Config::get().md.heat_steps)
    {
      std::cout << " -> [Step " << heatstep.offset << "]:[" << heatstep.raise << "K]";
    }
    std::cout << '\n';
  }
  std::cout << "Nose-Hoover thermostat is " << (Config::get().md.hooverHeatBath ? "active." : "inactive.") << '\n';
  if (Config::get().md.spherical.use)
  {
    std::cout << "Spherical boundaries will be applied.\n";
    std::cout << "Inner Sphere: (R: " << Config::get().md.spherical.r_inner << ", F: ";
    std::cout << Config::get().md.spherical.f1 << ", E: " << Config::get().md.spherical.e1 << ");";
    std::cout << "Outer Sphere: (R: " << Config::get().md.spherical.r_outer << ", F: ";
    std::cout << Config::get().md.spherical.f2 << ", E: " << Config::get().md.spherical.e2 << ")" << std::endl;
  }
  if (Config::get().md.rattle.use)
  {
    const std::size_t nr = Config::get().md.rattle.specified_rattle.size();
    if (Config::get().md.rattle.all) std::cout << "All covalent hydrogen bonds will be fixed\n";
    else if (nr > 0)
    {
      std::cout << "The following covalent hydrogen bonds will be fixed: \n";
      for (auto const & bond : Config::get().md.rattle.specified_rattle)
      {
        std::cout << "[" << bond.a << ":" << bond.b << "] ";
      }
      std::cout << std::endl;
    }
  }
  if (Config::get().md.resume)
  {
    std::cout << "The simulation state will be reverted to the state of the binary resume file '";
    std::cout << std::string(Config::get().general.outputFilename).append("_MD_restart.cbf") << "'." << std::endl;
  }
}

// generate random number for initial velocity distribution
static inline double ldrand(void)
{
  auto dist01 = std::uniform_real_distribution<double>{};
  return std::log(std::max(scon::random::threaded_rand(dist01), 1.0e-20));
}

// set up constraints for H-X bonds if requested
// ideal bond lengths are taken from the foce field parameter file
// specified by RATTpar in the INPUTFILE
void md::simulation::rattlesetup(void)
{
  config::md_conf::config_rattle::rattle_constraint_bond rctemp;
  //open rattle par file
  std::ifstream ifs;
  ifs.open(Config::get().md.rattle.ratpar.c_str());
  if (ifs.good()) std::cout << "Opened file for RATTLE parameters successfully." << std::endl;
  if (!ifs.good()) {
    std::cout << "Couldn't open file for RATTLE parameters. Check your input" << std::endl;
    throw;
  }
  // temp vars and vectors;
  struct trat {
    std::size_t ia, ib;
    double ideal;
  };
  struct ratatoms {
    std::size_t ia, ga;
  };
  ratatoms RTA;
  trat TRAT;
  std::vector<trat> temprat;
  std::vector<ratatoms> ratoms;
  std::size_t tia = std::size_t(), tib = std::size_t();
  char buffer[150];
  std::string bufferc;
  std::size_t found;
  std::vector < std::string > tokens;
  // read parameter file and extract ideal bond lengths
  while (!ifs.eof())
  {
    ifs.getline(buffer, 150);
    bufferc = buffer;
    found = bufferc.find("bond");
    if (found != std::string::npos)
    {
      tokens.clear();
      std::istringstream iss(bufferc);
      std::copy(std::istream_iterator <std::string>(iss), std::istream_iterator <std::string>(), std::back_inserter <std::vector <std::string >>(tokens));
      TRAT.ia = atoi(tokens[1].c_str());
      TRAT.ib = atoi(tokens[2].c_str());
      TRAT.ideal = atof(tokens[4].c_str());
      temprat.push_back(TRAT);
    }
  }// end of file check
  ifs.close();
  ifs.open(Config::get().md.rattle.ratpar.c_str());
  // read parameter file and get atoms
  while (!ifs.eof())
  {
    ifs.getline(buffer, 150);
    bufferc = buffer;
    found = bufferc.find("atom");
    if (found != std::string::npos)
    {
      tokens.clear();
      std::istringstream iss(bufferc);
      std::copy(std::istream_iterator <std::string>(iss), std::istream_iterator <std::string>(), std::back_inserter <std::vector <std::string >>(tokens));
      RTA.ia = atoi(tokens[1].c_str());
      RTA.ga = atoi(tokens[2].c_str());
      ratoms.push_back(RTA);
    }
  }
  ifs.close();
  // Generate vector with bonds which are to be constraint
  const std::size_t N = coordobj.size();
  //loop over all atoms
  for (std::size_t i = 0; i < N; i++) {
	  if (Config::get().md.rattle.all == true)  // all H-bonds are constraint
	  {
		  if (coordobj.atoms(i).number() == 1) //check if atom is hydrogen
		  {
			  rctemp.a = i;
			  rctemp.b = coordobj.atoms(i).bonds(0);
			  //loop over param vector
			  for (unsigned j = 0; j < ratoms.size(); j++)
			  {
				  if (ratoms[j].ia == coordobj.atoms(rctemp.a).energy_type()) tia = ratoms[j].ga;
				  if (ratoms[j].ia == coordobj.atoms(rctemp.b).energy_type()) tib = ratoms[j].ga;
			  }
			  for (unsigned k = 0; k < temprat.size(); k++)
			  {
				  // found matching parameters -> get ideal bond distance
				  if ((tia == temprat[k].ia || tia == temprat[k].ib) && (tib == temprat[k].ia || tib == temprat[k].ib))
				  {
					  rctemp.len = temprat[k].ideal;
				  }
			  }
			  rattle_bonds.push_back(rctemp);
		  }
	  }
	  else   // if MDrattle = 2 i.e. only spcified H-atoms are constrained
	  {
		  if (coordobj.atoms(i).number() == 1) //check if atom is hydrogen
		  {
			  rctemp.a = i;
			  rctemp.b = coordobj.atoms(i).bonds(0);
			  for (auto s : Config::get().md.rattle.specified_rattle)
			  {
				  if (s.a == rctemp.a)   // if H-atom is in the rattlebond list 
				  {
					  //loop over param vector
					  for (unsigned j = 0; j < ratoms.size(); j++)
					  {
						  if (ratoms[j].ia == coordobj.atoms(rctemp.a).energy_type()) tia = ratoms[j].ga;
						  if (ratoms[j].ia == coordobj.atoms(rctemp.b).energy_type()) tib = ratoms[j].ga;
					  }
					  for (unsigned k = 0; k < temprat.size(); k++)
					  {
						  // found matching parameters -> get ideal bond distance
						  if ((tia == temprat[k].ia || tia == temprat[k].ib) && (tib == temprat[k].ia || tib == temprat[k].ib))
						  {
							  rctemp.len = temprat[k].ideal;
						  }
					  }
					  rattle_bonds.push_back(rctemp);
				  }
			  }
		  }
	  }
    
  }
}

// Initialization of MD parameters and functions
void md::simulation::init(void)
{
	if (coordobj.validate_bonds() == false)  // test if there are broken bonds in the structure and save them
	{
		if (Config::get().general.verbosity > 1U)
		{
			std::cout << "Warning! Broken bonds in your structure even before the simulation starts! Atom numbers: \n";
			for (auto b : coordobj.broken_bonds)
			{
				std::cout << b[0] << " and " << b[1] << "\n";
			}
			throw std::logic_error("Please check your input structure.");
		}
	}
  using std::abs;
  std::size_t const N = coordobj.size();
  static double const twopi = 2.0*md::PI;
  double const twokbT = 2.0*md::kB*Config::get().md.T_init;
  M_total = 0;
  V.resize(N);
  F.resize(N);
  M.resize(N);
  P_start = coordobj.xyz();
  auto dist01 = std::uniform_real_distribution<double>{0,1};
  for (std::size_t i = 0; i<N; ++i)
  {
    // Get Atom Mass
    M[i] = coordobj.atoms(i).mass();
    // Sum total Mass
    M_total += M[i];
    // Set initial velocities to zero if fixed
    if (coordobj.atoms(i).fixed() || !(abs(Config::get().md.T_init) > 0.0)) V[i] = coords::Cartesian_Point(0);
    // initialize random velocities otherwise
    else
    { //ratioRoot = sqrt ( 2*kB*T / m )
      double const ratio(twopi*std::sqrt(twokbT / M[i]));
      V[i].x() = (std::sqrt(std::fabs(-2.0*ldrand())) *
        std::cos(scon::random::threaded_rand(dist01)*ratio));
      V[i].y() = (std::sqrt(std::fabs(-2.0*ldrand())) *
        std::cos(scon::random::threaded_rand(dist01)*ratio));
      V[i].z() = (std::sqrt(std::fabs(-2.0*ldrand())) *
        std::cos(scon::random::threaded_rand(dist01)*ratio));

      if (Config::get().general.verbosity > 4U) std::cout << "Initial Velocity of " << i << " is " << V[i] << " with rr: " << ratio << std::endl;
    }
    // sum position vectors for geometrical center
    C_geo += coordobj.xyz(i);
  }
  // get degrees of freedom
  freedom = 3U * N;
  if (Config::get().md.T_init != 0.0)
  {
	  //calculate temperature
	  double tempfactor(2.0 / (freedom*md::R));
	  updateEkin();
	  temp = E_kin * tempfactor;
	  // scale temperature
	  for (std::size_t i = 0; i < N; ++i)
	  {
		  if (!coordobj.atoms(i).fixed())
		  {
			  V[i].x() = V[i].x()*pow(Config::get().md.T_init / temp, 0.5);
			  V[i].y() = V[i].y()*pow(Config::get().md.T_init / temp, 0.5);
			  V[i].z() = V[i].z()*pow(Config::get().md.T_init / temp, 0.5);
		  }  
	  }
  }
  
  // Set up rattle vector for constraints
  if (Config::get().md.rattle.use == true)
  {
	  rattlesetup();
  }
  // constraint degrees of freedom
  if (Config::get().md.rattle.use == true) freedom -= rattle_bonds.size();
  // periodics and isothermal cases
  if (Config::get().md.hooverHeatBath == true)
  {
    if (Config::get().energy.periodic == true) freedom -= 3;
    else freedom -= 6;
  }
  if (Config::get().general.verbosity > 2U) std::cout << "Degrees of freedom: " << freedom << std::endl;
  // calc number of steps between snapshots (gap between snapshots)
  snapGap = gap(Config::get().md.num_steps, Config::get().md.num_snapShots);
  // scale geometrical center vector by number of atoms
  C_geo /= static_cast<double>(N);
  // call center of Mass method from coordinates object
  C_mass = coordobj.center_of_mass();

  // things for biased potentials
  inner_atoms.clear();    // in case of more than one MD, e.g. for an FEP calculation
  atoms_movable.clear();    
  if (Config::get().md.set_active_center == 1)
  {
	  distances = init_active_center(0);   //calculate initial active center and distances to active center
	  
	  for (int i(0U); i < N; ++i)  // determine which atoms are moved
	  {
		  if (distances[i] <= Config::get().md.outer_cutoff)
		  {
			  atoms_movable.push_back(i);
		  }
		  else   // set velocities of atoms that are not moved to zero
		  {
			  V[i] = coords::Cartesian_Point(0, 0, 0);
		  }
		  if (distances[i] <= Config::get().md.inner_cutoff)  //determine atoms inside inner cutoff
		  {                                                   // for temperature calculation
			  inner_atoms.push_back(i);
		  }
	  }
  }
  else   // if no active site is specified: all atoms are moved
  {
	  for (int i(0U); i < N; ++i)
	  {
		  atoms_movable.push_back(i);
	  }
  }
}

// If FEP calculation is requested: calculate lambda values for each window
// and print the scaling factors for van-der-Waals and electrostatics for each window
void md::simulation::fepinit(void)
{
  if (Config::get().general.verbosity > 1)   // print warning if vdw-softcore potential is shifted so strongly that is has no minimum
  {
    if (Config::get().general.energy_interface == config::interface_types::CHARMM22 || Config::get().general.energy_interface == config::interface_types::AMBER)
    {
      if (Config::get().fep.ljshift > 1.0 / ((1 - Config::get().fep.dlambda)*(1 - Config::get().fep.dlambda)))
      {
        std::cout << "WARNING! You should choose a smaller value for vshift!\n\n";
      }
    }
    else if (Config::get().general.energy_interface == config::interface_types::OPLSAA)
    {
      if (Config::get().fep.ljshift > 2.0 / ((1 - Config::get().fep.dlambda)*(1 - Config::get().fep.dlambda)))
      {
        std::cout << "WARNING! You should choose a smaller value for vshift!\n\n";
      }
    }
  }
  init();
  // center and temp var
  coordobj.move_all_by(-coordobj.center_of_geometry());
  double linear, dlin, linel;
  double  increment, tempo, tempo2, diff;
  Config::set().md.fep = true;
  FEPsum = 0.0;
  // increment electrostatics
  linear = 1.0 - Config::get().fep.eleccouple;
  dlin = linear / Config::get().fep.dlambda;
  linel = 1.0 / dlin;
  // fill vector with scaling increments
  increment = Config::get().fep.lambda / Config::get().fep.dlambda;
  //std::cout << Config::get().fep.lambda << "   " << Config::get().fep.dlambda << std::endl;
  std::cout << "Number of FEP windows:  " << increment << std::endl;
  coordobj.fep.window.resize(std::size_t(increment));
  coordobj.fep.window[0].step = 0;
  // Calculate lambda and dlambda for Electrostatics
  for (std::size_t i = 0; i< coordobj.fep.window.size(); i++) {

    tempo = i * Config::get().fep.dlambda;
    tempo2 = i *  Config::get().fep.dlambda + Config::get().fep.dlambda;
    if (tempo <= Config::get().fep.eleccouple)
    {
      coordobj.fep.window[i].ein = 0.0;
    }
    else
    {
      diff = std::abs(tempo - Config::get().fep.eleccouple);
      coordobj.fep.window[i].ein = (diff / Config::get().fep.dlambda) * linel;
      if (coordobj.fep.window[i].ein > 1.0) coordobj.fep.window[i].ein = 1.0;
    }
    if (tempo >= (1 - Config::get().fep.eleccouple))
    {
      coordobj.fep.window[i].eout = 0.0;
    }
    else {
      diff = tempo / Config::get().fep.dlambda;
      coordobj.fep.window[i].eout = 1.0 - (diff * linel);
    }
    if (tempo2 <= Config::get().fep.eleccouple) {
      coordobj.fep.window[i].dein = 0.0;
    }
    else {
      diff = std::abs(tempo2 - Config::get().fep.eleccouple);
      coordobj.fep.window[i].dein = (diff / Config::get().fep.dlambda) * linel;
      if (coordobj.fep.window[i].dein > 1.0) coordobj.fep.window[i].dein = 1.0;
    }
    if (tempo2 >= (1 - Config::get().fep.eleccouple)) {
      coordobj.fep.window[i].deout = 0.0;
    }
    else {
      diff = tempo2 / Config::get().fep.dlambda;
      coordobj.fep.window[i].deout = 1.0 - (diff * linel);
    }
    // calculate lambda and dlambda for van-der-Waals
    if (Config::get().fep.vdwcouple == 0) {
      coordobj.fep.window[i].vin = 1.0;
      coordobj.fep.window[i].dvin = 1.0;
      coordobj.fep.window[i].vout = 1.0;
      coordobj.fep.window[i].dvout = 1.0;
    }
    else if (Config::get().fep.vdwcouple == 1.0) {
      coordobj.fep.window[i].vout = 1.0 - i * Config::get().fep.dlambda;
      coordobj.fep.window[i].vin = i * Config::get().fep.dlambda;

      coordobj.fep.window[i].dvout = 1.0 - (i+1) * Config::get().fep.dlambda;
      coordobj.fep.window[i].dvin = (i+1) * Config::get().fep.dlambda;
    }
    else 
	{
		coordobj.fep.window[i].vout = (1.0 - i * Config::get().fep.dlambda)/ Config::get().fep.vdwcouple;
		coordobj.fep.window[i].vin = (i * Config::get().fep.dlambda)/ Config::get().fep.vdwcouple;

		coordobj.fep.window[i].dvout = (1.0 - (i + 1) * Config::get().fep.dlambda)/ Config::get().fep.vdwcouple;
		coordobj.fep.window[i].dvin = ((i + 1) * Config::get().fep.dlambda)/ Config::get().fep.vdwcouple;
		if (coordobj.fep.window[i].vout > 1)
		{
			coordobj.fep.window[i].vout = 1;
		}
		if (coordobj.fep.window[i].vin > 1)
		{
			coordobj.fep.window[i].vin = 1;
		}
		if (coordobj.fep.window[i].dvout > 1)
		{
			coordobj.fep.window[i].dvout = 1;
		}
		if (coordobj.fep.window[i].dvin > 1)
		{
			coordobj.fep.window[i].dvin = 1;
		}
      
    }
  }// end of loop
  // clear FEP output vector and print lambvda values
  std::ofstream fepclear;
  fepclear.open("alchemical.txt");
  fepclear.close();
  coordobj.fep.window[0].step = 0;
  std::cout << "FEP Coupling Parameters:" << std::endl;
  std::cout << std::endl;
  std::cout << "Van-der-Waals Coupling: " << std::endl;
  for (std::size_t i = 0; i < coordobj.fep.window.size(); i++)
  {
    std::cout << std::setw(8) << coordobj.fep.window[i].vout << std::setw(8) << coordobj.fep.window[i].dvout << std::setw(8) << coordobj.fep.window[i].vin << std::setw(8) << coordobj.fep.window[i].dvin << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Electrostatic Coupling:" << std::endl;
  for (std::size_t i = 0; i < coordobj.fep.window.size(); i++)
  {
    std::cout << std::setw(8) << coordobj.fep.window[i].eout << std::setw(8) << coordobj.fep.window[i].deout << std::setw(8) << coordobj.fep.window[i].ein << std::setw(8) << coordobj.fep.window[i].dein << std::endl;
  }


  // Check for PME flag, if yes allocate needed memeory for the grids
  /*if (Config::get().energy.pme == true)
  {
    std::cout << "PME electrostatics for FEP calculation. Setting up memory" << std::endl;
    coordobj.getfepinfo();
    // generate grids for incoming atoms, outgoing atoms and rest of the system (3 = x,y,z; size = atoms in the respective system)  
    coordobj.pme.pmetemp.initingrid.Allocate(3, coordobj.pme.pmetemp.fepi.size());
    coordobj.pme.pmetemp.initoutgrid.Allocate(3, coordobj.pme.pmetemp.fepo.size());
    coordobj.pme.pmetemp.initallgrid.Allocate(3, coordobj.pme.pmetemp.fepa.size());
    // Allocate memory for parallel execution of the PME calculation in FEP simulations
    coordobj.pme.pmetemp.parallelinpme.Allocate(coordobj.pme.pmetemp.fepi.size(), coordobj.pme.pmetemp.rgridtotal);
    coordobj.pme.pmetemp.paralleloutpme.Allocate(coordobj.pme.pmetemp.fepo.size(), coordobj.pme.pmetemp.rgridtotal);
    coordobj.pme.pmetemp.parallelallpme.Allocate(coordobj.pme.pmetemp.fepa.size(), coordobj.pme.pmetemp.rgridtotal);
  }*/
}

//Calculation of ensemble average and free energy change for each step if FEP calculation is performed
// calculation can be improved if at every step the current averages are stored
// currently calculation is performed at the end of each window
void md::simulation::freecalc()
{
  std::size_t iterator(0U), k(0U);
  // set conversion factors (conv) and constants (boltzmann, avogadro)
  double de_ensemble, temp_avg, boltz = 1.3806488E-23, avogad = 6.022E23, conv = 4184.0;
  // calculate ensemble average for current window
  for (std::size_t i = 0; i < coordobj.fep.fepdata.size(); i++) 
  {                // for every conformation in window
    iterator += 1;
    k = 0;
    de_ensemble = temp_avg = 0.0;
    coordobj.fep.fepdata[i].de_ens = exp(-1 / (boltz*coordobj.fep.fepdata[i].T)*conv*coordobj.fep.fepdata[i].dE / avogad);
    for (k = 0; k <= i; k++) {
      temp_avg += coordobj.fep.fepdata[k].T;
      de_ensemble += coordobj.fep.fepdata[k].de_ens;
    }
    de_ensemble = de_ensemble / iterator;
    temp_avg = temp_avg / iterator;
    coordobj.fep.fepdata[i].dG = -1 * std::log(de_ensemble)*temp_avg*boltz*avogad / conv;
  }// end of main loop
  // calculate final free energy change for the current window
  this->FEPsum += coordobj.fep.fepdata[coordobj.fep.fepdata.size() - 1].dG;
}

// write the output FEP calculations
// backward transformations are not supported anymore!!!!
void md::simulation::freewrite(std::size_t i)
{

  std::ofstream fep("alchemical.txt", std::ios_base::app);
  std::ofstream res("FEP_Results.txt", std::ios_base::app);
  // standard forward output
  if (i*Config::get().fep.dlambda == 0 && this->prod == false)
  {
    res << std::fixed << std::right << std::setprecision(4) << std::setw(10) << "0" << std::setw(10) << "0" << std::endl;
  }
  // equilibration is performed
  if (this->prod == false) {
    fep << "Equilibration for Lambda =  " << i * Config::get().fep.dlambda << 
      "  and dLambda =  " << (i * Config::get().fep.dlambda) + Config::get().fep.dlambda << std::endl;
  }
  // production run is performed
  else if (this->prod == true) {
    fep << "Starting new data collection with values:  " << 
      i * Config::get().fep.dlambda << "   " << (i * Config::get().fep.dlambda) + Config::get().fep.dlambda << std::endl;
  }
  // write output to alchemical.txt
  for (std::size_t k = 0; k < coordobj.fep.fepdata.size(); k++) {
    if (k%Config::get().fep.freq == 0) {
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].e_c_l1;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].e_c_l2;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].e_vdw_l1;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].e_vdw_l2;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].T;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].dE;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].dG;
      fep << std::endl;
    }
	if (Config::get().general.verbosity > 3u)
	{
		std::cout << "Coulomb: " << coordobj.fep.fepdata[k].e_c_l2 - coordobj.fep.fepdata[k].e_c_l1 << " ,vdW: " << coordobj.fep.fepdata[k].e_vdw_l2 - coordobj.fep.fepdata[k].e_vdw_l1 << "\n";
	}
  }
  // at the end of production data in alchemical.txt sum up the results and print the before the new window starts
  if (this->prod == true) {
    fep << "Free energy change for the current window:  ";
    fep << coordobj.fep.fepdata[coordobj.fep.fepdata.size() - 1].dG << std::endl;
    fep << "Total free energy change until current window:  " << FEPsum << std::endl;
    //fep <<  FEPsum << std::endl;
    fep << "End of collection. Increasing lambda value" << std::endl;
    res << std::fixed << std::right << std::setprecision(4) << std::setw(10) << (i * Config::get().fep.dlambda) + Config::get().fep.dlambda << std::setw(10) << FEPsum << std::endl;
  }
}

// perform FEP calculation if requested
void md::simulation::feprun()
{
    for (std::size_t i(0U); i < coordobj.fep.window.size(); ++i)  //for every window
    {
      std::cout << "Lambda:  " << i * Config::get().fep.dlambda << std::endl;
      coordobj.fep.window[0U].step = static_cast<int>(i);
      coordobj.fep.fepdata.clear();
      // equilibration run for window i
      Config::set().md.num_steps = Config::get().fep.equil;
      integrate(true);
      // write output for equlibration and clear fep vector
      this->prod = false;
      freewrite(i);
      coordobj.fep.fepdata.clear();
      // production run for window i
      Config::set().md.num_steps = Config::get().fep.steps;
      integrate(true);
      this->prod = true;
      // calculate free energy change for window and write output
      freecalc();
      freewrite(i);
    }// end of main window loop
}

// eliminate translation and rotation of the system at the beginning of a MD simulation
// can also be performed at the end of every MD step
void md::simulation::tune_momentum(void)
{
  std::size_t const N = coordobj.size();
  coords::Cartesian_Point momentum_linear, momentum_angular, mass_vector, velocity_angular;
  // 3x3 Matrix (moment of inertia tensor)
  scon::c3<scon::c3<coords::float_type>> InertiaTensor;
  // calculate system movement
  for (std::size_t i = 0; i < N; ++i)
  {
    // center of mass and linear and angular momentum
    mass_vector += coordobj.xyz(i) * M[i];
    momentum_linear += V[i] * M[i];
    momentum_angular += cross(coordobj.xyz(i), V[i])*M[i];
  }
  // scale by total mass
  momentum_linear /= M_total;
  mass_vector /= M_total;
  // MVxLM * M 
  momentum_angular -= cross(mass_vector, momentum_linear)*M_total;
  // momentum of inertia from each component
  double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
  for (std::size_t i = 0; i < N; ++i)
  {
    coords::Cartesian_Point r(coordobj.xyz(i) - mass_vector);
    xx += r.x()*r.x()*M[i];
    xy += r.x()*r.y()*M[i];
    xz += r.x()*r.z()*M[i];
    yy += r.y()*r.y()*M[i];
    yz += r.y()*r.z()*M[i];
    zz += r.z()*r.z()*M[i];
  }
  // Set up Inertia Tensor
  InertiaTensor.x() = scon::c3<coords::float_type>(yy + zz, -xy, -xz);
  InertiaTensor.y() = scon::c3<coords::float_type>(-xy, xx + zz, -yz);
  InertiaTensor.z() = scon::c3<coords::float_type>(-xz, -yz, xx + yy);
  // Invert inertia tensor
  scon::invert(InertiaTensor);
  // get angular velocity
  velocity_angular.x() = (dot(InertiaTensor.x(), momentum_angular));
  velocity_angular.y() = (dot(InertiaTensor.y(), momentum_angular));
  velocity_angular.z() = (dot(InertiaTensor.z(), momentum_angular));
  // remove angular and linear momentum
  for (std::size_t i = 0; i < N; ++i)
  {
    V[i] -= momentum_linear;
    coords::Cartesian_Point r(coordobj.xyz(i) - mass_vector);
    V[i] -= cross(velocity_angular, r);
  }
  if (Config::get().general.verbosity > 3u) std::cout << "Tuned momentum \n";
}

// call function for spherical boundary conditions
void md::simulation::boundary_adjustments()
{
  if (Config::get().md.spherical.use)
  {
    if (Config::get().general.verbosity > 3u) std::cout << "Adjusting boundary conditions.\n";
    spherical_adjust();
  }
}

// adjust gradients on atoms if spherical boundary conditions are used
void md::simulation::spherical_adjust()
{
  
  std::size_t const N = coordobj.size();
  for (std::size_t i(0U); i < N; ++i)
  {
    // distance from geometrical center
    coords::Cartesian_Point r = coordobj.xyz(i) - C_geo;
    double const L = geometric_length(r);
    // if distance is greater than inner radius of spherical boundary -> adjust atom
    //std::cout << coordobj.xyz(i) << "   " << C_geo << std::endl;
    //std::cout << L << "   " << Config::get().md.spherical.r_inner << std::endl;
    if (L > Config::get().md.spherical.r_inner)
    {
      // E = f1*d^e1 -> deriv: e1*f1*d^e1-1 [/L: normalization]
      coords::Cartesian_Point sp_g;
      sp_g = r*(Config::get().md.spherical.e1 * Config::get().md.spherical.f1 * std::pow(L - Config::get().md.spherical.r_inner, Config::get().md.spherical.e1 - 1) / L);
      coordobj.add_sp_gradients(i, sp_g);
    }
  }
}

// calculate kinetic energy of the system
// we use the tensor formulation of the kinetic energy
// tensor formulation is needed for anisotopic box shapes if constant pressure is applied
void md::simulation::updateEkin(void)
{
  // initialize tensor to zero
  using TenVal = coords::Tensor::value_type;
  TenVal z = {0.0,0.0,0.0};
  E_kin_tensor.fill(z);
  auto const N = V.size();
  // calculate contribution to kinetic energy for each atom
  for (std::size_t i = 0; i < N; ++i)
  {
	  auto const fact = 0.5 * M[i] / convert;
	  E_kin_tensor[0][0] += fact * V[i].x() * V[i].x();
	  E_kin_tensor[1][0] += fact * V[i].x() * V[i].y();
	  E_kin_tensor[2][0] += fact * V[i].x() * V[i].z();
	  E_kin_tensor[0][1] += fact * V[i].y() * V[i].x();
	  E_kin_tensor[1][1] += fact * V[i].y() * V[i].y();
	  E_kin_tensor[2][1] += fact * V[i].y() * V[i].z();
	  E_kin_tensor[0][2] += fact * V[i].z() * V[i].x();
	  E_kin_tensor[1][2] += fact * V[i].z() * V[i].y();
	  E_kin_tensor[2][2] += fact * V[i].z() * V[i].z();
  }
  // calculate total kinetic energy by the trace of the tensor
  E_kin = E_kin_tensor[0][0] + E_kin_tensor[1][1] + E_kin_tensor[2][2];
  if (Config::get().general.verbosity > 4u)
  {
    std::cout << "Updating kinetic Energy from " << E_kin_tensor[0][0] << ", "
      << E_kin_tensor[1][1] << ", " << E_kin_tensor[2][2] << " to " << E_kin << '\n';
  }
}

void md::simulation::updateEkin_some_atoms(std::vector<int> atom_list)
{
	// initialize tensor to zero
	using TenVal = coords::Tensor::value_type;
	TenVal z = { 0.0,0.0,0.0 };
	E_kin_tensor.fill(z);
	// calculate contribution to kinetic energy for each atom
	for (auto i : atom_list)
	{
		auto const fact = 0.5 * M[i] / convert;
		E_kin_tensor[0][0] += fact * V[i].x() * V[i].x();
		E_kin_tensor[1][0] += fact * V[i].x() * V[i].y();
		E_kin_tensor[2][0] += fact * V[i].x() * V[i].z();
		E_kin_tensor[0][1] += fact * V[i].y() * V[i].x();
		E_kin_tensor[1][1] += fact * V[i].y() * V[i].y();
		E_kin_tensor[2][1] += fact * V[i].y() * V[i].z();
		E_kin_tensor[0][2] += fact * V[i].z() * V[i].x();
		E_kin_tensor[1][2] += fact * V[i].z() * V[i].y();
		E_kin_tensor[2][2] += fact * V[i].z() * V[i].z();
	}
	// calculate total kinetic energy by the trace of the tensor
	E_kin = E_kin_tensor[0][0] + E_kin_tensor[1][1] + E_kin_tensor[2][2];
	if (Config::get().general.verbosity > 4u)
	{
		std::cout << "Updating kinetic Energy of some atoms from " << E_kin_tensor[0][0] << ", "
			<< E_kin_tensor[1][1] << ", " << E_kin_tensor[2][2] << " to " << E_kin << '\n';
	}
}

// apply pressure corrections if constant pressure simulation is performed
void md::simulation::berendsen(double const & time)
{
  // temp variables
  const std::size_t N = coordobj.size();
  double fac, scale, volume;
  double ptensor[3][3], aniso[3][3], anisobox[3][3], dimtemp[3][3];
  // get volume and pressure scaling factor
  volume = Config::get().energy.pb_box.x() *  Config::get().energy.pb_box.y() *  Config::get().energy.pb_box.z();
  fac = presc / volume;
  // pressure for ISOTROPIC boxes
  if (Config::get().energy.isotropic == true) {
    ptensor[0][0] = fac * 2.0 * (E_kin_tensor[0][0] - coordobj.virial()[0][0]);
    ptensor[1][1] = fac * 2.0 * (E_kin_tensor[1][1] - coordobj.virial()[1][1]);
    ptensor[2][2] = fac * 2.0 * (E_kin_tensor[2][2] - coordobj.virial()[2][2]);
    press = (ptensor[0][0] + ptensor[1][1] + ptensor[2][2]) / 3.0;  // pressure in bar or atm ???
    // Berendsen scaling for isotpropic boxes
    scale = std::pow((1.0 + (time*Config::get().md.pcompress / Config::get().md.pdelay)*(press - Config::get().md.ptarget)), 0.3333333333333);
    // Adjust box dimensions
    Config::set().energy.pb_box.x() = Config::get().energy.pb_box.x() * scale;
    Config::set().energy.pb_box.y() = Config::get().energy.pb_box.y() * scale;
    Config::set().energy.pb_box.z() = Config::get().energy.pb_box.z() * scale;
    // scale atomic coordinates
    for (size_t i = 0; i < N; ++i)
    {
      coordobj.scale_atom_by(i, scale);
    }
  }
  //  pressure for ANISOTROPIC boxes
  else {
    // get pressure tensor for anisotropic systems
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        ptensor[j][i] = fac * (2.0*E_kin_tensor[j][i] - coordobj.virial()[j][i]);
      }
    }
    // get isotropic pressure value
    press = (ptensor[0][0] + ptensor[1][1] + ptensor[2][2]) / 3.0;
    // Anisotropic scaling factors
    scale = time * Config::get().md.pcompress / (3.0*Config::get().md.pdelay);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (j == i) aniso[j][i] = 1.0 + scale * (ptensor[i][i] - Config::get().md.ptarget);
        else aniso[j][i] = scale * ptensor[j][i];
      }
    }
    // modify anisotropic box dimensions (only for cube or octahedron)
    // the 1.0 and  0.0 are placeholder for the scaling factors is monoclinic or triclinic boxes are introduced
    dimtemp[0][0] = Config::get().energy.pb_box.x();
    dimtemp[1][0] = 0.0;
    dimtemp[2][0] = 0.0;
    dimtemp[0][1] = Config::get().energy.pb_box.y() * 1.0;
    dimtemp[1][1] = Config::get().energy.pb_box.y() * 0.0;
    dimtemp[2][1] = 0.0;
    dimtemp[0][2] = Config::get().energy.pb_box.z() * 0.0;
    dimtemp[1][2] = Config::get().energy.pb_box.z() * (1.0 - 1.0 * 1.0) / 1.0;
    dimtemp[2][2] = Config::get().energy.pb_box.z() * sqrt(1.0 * 1.0 - 0.0 - 0.0);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        anisobox[j][i] = 0.0;
        for (int k = 0; k < 3; k++) {
          anisobox[j][i] = anisobox[j][i] + aniso[j][k] * dimtemp[k][i];
        }
      }
    }
    // scale box dimensions for anisotropic case
    Config::set().energy.pb_box.x() = std::sqrt((anisobox[0][0] * anisobox[0][0] + anisobox[1][0] * anisobox[1][0] + anisobox[2][0] * anisobox[2][0]));
    Config::set().energy.pb_box.y() = std::sqrt((anisobox[0][1] * anisobox[0][1] + anisobox[1][1] * anisobox[1][1] + anisobox[2][1] * anisobox[2][1]));
    Config::set().energy.pb_box.z() = std::sqrt((anisobox[0][2] * anisobox[0][2] + anisobox[1][2] * anisobox[1][2] + anisobox[2][2] * anisobox[2][2]));
    // scale atomic coordinates for anisotropic case
    for (size_t i = 0; i < N; ++i)
    {
      coordobj.set_atom_aniso(i, aniso, true);
    }
  }//end of isotropic else
}

// determine target temperature for heating
bool md::simulation::heat(std::size_t const step, bool fep)
{
	if (Config::get().md.heat_steps.size() == 0)  // no temperature control
	{
		return false;
	}
	else  // temperature control
	{
		config::md_conf::config_heat last;
		last.raise = Config::get().md.T_init;
		if (fep == true)   // keep constant temperature in FEP calculation
		{
			T = Config::get().md.T_final;
			return true;
		}
		for (auto const & heatstep : Config::get().md.heat_steps)  
		{
			if (heatstep.offset >= step)    // find first heatstep after current step
			{
				double const delta((heatstep.raise - last.raise) / static_cast<double>(heatstep.offset - last.offset));
				T += delta;  // adjust target temperature
				return true; // exit function
			}
			last = heatstep; // last heatstep before current step
		}
		if (step > Config::get().md.heat_steps[Config::get().md.heat_steps.size()-1].offset)
		{             // after last heatstep: keep final temperature
			T = Config::get().md.T_final;
			return true;
		}
		else
		{
			std::cout << "This should not happen!\n";
			throw std::exception();
		}
	}
}

// Nose-Hover thermostat. Variable names and implementation are identical to the book of
// Frenkel and Smit, Understanding Molecular Simulation, Appendix E
void md::simulation::nose_hoover_thermostat(void)
{
  double tempscale(0.0);
  std::size_t const N(coordobj.size());
  double const delt(Config::get().md.timeStep),
    d2(delt / 2.0), d4(d2 / 2.0), d8(d4 / 2.0),
    TR(T*md::R), fTR(freedom*TR);
  nht.G2 = (nht.Q1*nht.v1*nht.v1 - TR) / nht.Q2;
  nht.v2 += nht.G2*d4;
  nht.v1 *= exp(-nht.v2*d8);
  nht.G1 = (2.0*E_kin - fTR) / nht.Q1;
  nht.v1 += nht.G1*d4;
  nht.v1 *= exp(-nht.v2*d8);
  nht.x1 += nht.v1*d2;
  nht.x2 += nht.v2*d2;
  tempscale = exp(-nht.v1*d2);
  for (std::size_t i(0U); i < N; ++i) V[i] *= tempscale;
  E_kin *= tempscale*tempscale;
  if (Config::get().general.verbosity > 4u)
  {
    std::cout << "Nose-Hoover-Adjustment; Scaling factor: " << tempscale << '\n';
  }
  nht.v1 *= exp(-nht.v2*d8);
  nht.G1 = (2.0*E_kin - fTR) / nht.Q1;
  nht.v1 += nht.G1*d4;
  nht.v1 *= exp(-nht.v2*d8);
  nht.G2 = (nht.Q1*nht.v1*nht.v1 - TR) / nht.Q2;
  nht.v2 += nht.G2*d4;
}

double md::simulation::tempcontrol(bool thermostat, bool half)
{
	std::size_t const N = this->coordobj.size();  // total number of atoms
	double tempfactor(2.0 / (freedom*md::R));     // factor for calculation of temperature from kinetic energy  
	double temp, temp2=0, factor;     // current temperature before and after the temperature scaling, scaling factor

	if (thermostat)   // apply nose-hoover thermostat
	{
		if (Config::get().general.verbosity > 3 && half)
		{
			std::cout << "hoover halfstep\n";
		}
		else if (Config::get().general.verbosity > 3)
		{
			std::cout << "hoover fullstep\n";
		}
		updateEkin();
		nose_hoover_thermostat();
		temp2 = E_kin * tempfactor;
	}
	else // no thermostat -> direct scaling
	{
		if (Config::get().md.set_active_center == 1)
		{     // calculate temperature only for atoms inside inner cutoff
			updateEkin_some_atoms(inner_atoms); // kinetic energy of inner atoms
			size_t dof = 3u * inner_atoms.size();
			double T_factor = (2.0 / (dof*md::R));
			temp = E_kin*T_factor;           // temperature of inner atoms
			factor = std::sqrt(T / temp);    // temperature scaling factor
			for (auto i: atoms_movable) V[i] *= factor;   // new velocities (for all atoms that have a velocity)
			if (half == false)
			{
				updateEkin_some_atoms(inner_atoms);
				temp2 = E_kin * T_factor;                   // new temperature of inner atoms
				updateEkin();            // kinetic energy
			}	
		}
		else
		{
			updateEkin();
			temp = E_kin * tempfactor;      // temperature before
			factor = std::sqrt(T / temp);
			for (size_t i(0U); i < N; ++i) V[i] *= factor;  // new velocities
			if (half == false)
			{
				updateEkin();
				temp2 = E_kin * tempfactor;     // temperatures after
			}
		}
		
		
		if (Config::get().general.verbosity > 3 && half)
		{
			std::cout << "half step: desired temp: " << T << " current temp: " << temp << " factor: " << factor << "\n";
		}
		else if (Config::get().general.verbosity > 3)
		{
			std::cout << "full step: desired temp: " << T << " current temp: " << temp << " factor: " << factor << "\n";
		}		
	}
	return temp2;
}

std::vector<double> md::simulation::init_active_center(int counter)
{
	std::size_t const N = this->coordobj.size();                // total number of atoms
	config::molecular_dynamics const & CONFIG(Config::get().md);

	auto split = std::max(std::min(std::size_t(CONFIG.num_steps / 100u), size_t(10000u)), std::size_t{ 100u });

	std::vector<double> distances;
	std::vector<coords::Cartesian_Point> coords_act_center;
	for (auto & atom_number : Config::get().md.active_center)
	{
		if (atom_number > 0 && atom_number <= N)
		{
			coords_act_center.push_back(coordobj.xyz(atom_number - 1));  //(-1) because atom count in tinker starts with 1, not with 0
		}
		else
		{
			std::cout << "ERROR: Atom number " << atom_number << " for active site not valid!!!\n";
			throw std::exception();
		}
	}

	coords::Cartesian_Point summe_coords_act_center; //calculate geometrical center of active site
	for (auto & atom_coords : coords_act_center)
	{
		summe_coords_act_center += atom_coords;
	}
	coords::Cartesian_Point C_geo_act_center = summe_coords_act_center / double(coords_act_center.size());

	if (Config::get().general.verbosity > 2 && counter % split == 0)
	{
		std::cout << "Coordinates of active site: " << C_geo_act_center << "\n";
	}

	for (std::size_t i(0U); i < N; ++i)  // calculate distance to active center for every atom
	{
		coords::Cartesian_Point coords_atom = coordobj.xyz(i);
		double dist_x = C_geo_act_center.x() - coords_atom.x();
		double dist_y = C_geo_act_center.y() - coords_atom.y();
		double dist_z = C_geo_act_center.z() - coords_atom.z();
		double distance = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
		distances.push_back(distance);
		if (Config::get().general.verbosity > 3)
		{
			std::cout << "Atom " << i + 1 << ": Distance to active center: " << distance << "\n";
		}
	}
	return distances;
}

coords::Cartesian_Point md::simulation::adjust_velocities(int atom_number, double inner_cutoff, double outer_cutoff)
{
	double distance = distances[atom_number];
	coords::Cartesian_Point velocity = V[atom_number];

	if (distance > inner_cutoff)  // adjust velocities between inner and outer cutoff
	{
		velocity = velocity - velocity * ((distance - inner_cutoff) / (outer_cutoff - inner_cutoff));
		return velocity;
	}
	else if (distance <= inner_cutoff)  // normal movement inside inner cutoff
	{
		inner_atoms.push_back(atom_number); // save atoms for temperature calculation 
		return velocity;
	}
	else    // should not happen
	{ 
		if (distance > outer_cutoff)
		{
			velocity = coords::Cartesian_Point(0, 0, 0);
			std::cout << "This should not happen.\n"; //because velocities for these atoms are not calculated
			return velocity;
		}
		else
		{
			std::cout << "ERROR: really strange distance for atom " << atom_number << ": " << distance << "\n";
			exit(EXIT_FAILURE);
		}
	}
}

void md::simulation::restart_broken()
{
	if (Config::get().general.verbosity > 2U)
	{
		std::cout<<"Start again...\n";
	}
	coordobj.set_xyz(P_start);   // set all positions to start
	static double const twopi = 2.0*md::PI;
	double const twokbT = 2.0*md::kB*T;
	auto dist01 = std::uniform_real_distribution<double>{ 0,1 };
	std::size_t const N = coordobj.size();
	for (int i(0U); i < N; ++i)     // set random velocities around temperature T
	{
		double const ratio(twopi*std::sqrt(twokbT / M[i]));
		V[i].x() = (std::sqrt(std::fabs(-2.0*ldrand())) *
			std::cos(scon::random::threaded_rand(dist01)*ratio));
		V[i].y() = (std::sqrt(std::fabs(-2.0*ldrand())) *
			std::cos(scon::random::threaded_rand(dist01)*ratio));
		V[i].z() = (std::sqrt(std::fabs(-2.0*ldrand())) *
			std::cos(scon::random::threaded_rand(dist01)*ratio));
	}
	if (Config::get().md.set_active_center == 0)  
	{
		tune_momentum();     // remove translation and rotation
	}
}


// Velcoity verlet integrator
void md::simulation::velocity_verlet(bool fep, std::size_t k_init)
{
  scon::chrono::high_resolution_timer integration_timer;

  config::molecular_dynamics const & CONFIG(Config::get().md);

  // prepare tracking

  std::size_t const N = this->coordobj.size();
  // set average pressure to zero
  double p_average(0.0);
  // constant values
  double const
	  //velofactor(-0.5*dt*md::convert),
	  dt_2(0.5*dt);
    
  //inner and outer cutoff for biased potential
  double inner_cutoff = Config::get().md.inner_cutoff;
  double outer_cutoff = Config::get().md.outer_cutoff;

  if (Config::get().general.verbosity > 0U)
  {
    std::cout << "Saving " << std::size_t(snapGap > 0 ? (CONFIG.num_steps - k_init) / snapGap : 0);
    std::cout << " snapshots (" << Config::get().md.num_snapShots << " in config)\n";
  }
  // Main MD Loop
  auto split = std::max(std::min(std::size_t(CONFIG.num_steps / 100u ), size_t(10000u)), std::size_t{ 100u });
  for (std::size_t k(k_init); k < CONFIG.num_steps; ++k)
  {
    bool const HEATED(heat(k,fep));
    if (Config::get().general.verbosity > 1u && k % split == 0 && k > 1)
    {
      std::cout << k << " of " << CONFIG.num_steps << " steps completed\n";
    }
	
    // apply half step temperature corrections
	if (CONFIG.hooverHeatBath || HEATED)
	{
		temp = tempcontrol(CONFIG.hooverHeatBath, true);
	}

    // save old coordinates
    P_old = coordobj.xyz();
    // Calculate new positions and half step velocities
    for (auto i : atoms_movable)
    {
      // calc acceleration
      coords::Cartesian_Point const acceleration(coordobj.g_xyz(i)*md::negconvert / M[i]);

	  // update veloctiy
	  V[i] += acceleration*dt_2;

	  inner_atoms.clear();
	  if (Config::get().md.set_active_center == 1)  //adjustment of velocities by distance to active center
	  { 
		  V[i] = adjust_velocities(static_cast<int>(i), inner_cutoff, outer_cutoff);
	  }

      if (Config::get().general.verbosity > 4)
      {
        std::cout << "Move " << i << " by " << (V[i] * dt) 
          << " with g " << coordobj.g_xyz(i) << ", V: " << V[i] << std::endl;
      }

      // update coordinates
      coordobj.move_atom_by(i, V[i] * dt);
    }

	if (coordobj.validate_bonds() == false)  // look if all bonds are okay and save those which aren't 
	{
		if (Config::get().general.verbosity > 1U)
		{                                          // give warning if there are broken bonds
			std::cout << "Warning! Broken bonds between atoms...\n";
			for (auto b : coordobj.broken_bonds)
			{
				std::cout << b[0] << " and " << b[1] << ", distance: "<<b[2]<< "\n";
			}
		}
		if (Config::get().md.broken_restart == 1)
		{
			restart_broken();   // if desired: set simulation to original positions and random velocities
		}
	}
	if (Config::get().md.set_active_center == 1 && Config::get().md.adjustment_by_step == 1) 
	{
		distances = init_active_center(static_cast<int>(k));  //calculate active center and new distances to active center for every step
		atoms_movable.clear();            // determine again which atoms are moved
		for (int i(0U); i < N; ++i)
		{
			if (distances[i] <= outer_cutoff)
			{
				atoms_movable.push_back(i);
			}
			else
			{
				V[i] = coords::Cartesian_Point(0, 0, 0);
			}
		}
	}
    // Apply first part of RATTLE constraints if requested
    if (CONFIG.rattle.use) rattle_pre();
	
    // calculate new energy & gradients
    coordobj.g();
    // Apply umbrella potential if umbrella sampling is used
    if (CONFIG.umbrella == true)
    {
      coordobj.ubias(udatacontainer);
    }
    // refine nonbondeds if refinement is required due to configuration
    if (CONFIG.refine_offset != 0 && (k + 1U) % CONFIG.refine_offset == 0)
    {
      if (Config::get().general.verbosity > 3U) std::cout << "Refining structure/nonbondeds.\n";
      coordobj.energy_update(true);
    }
    // If spherical boundaries are used apply boundary potential
    boundary_adjustments();
    // add new acceleration and calculate full step velocities
	inner_atoms.clear();
    for (auto i: atoms_movable)
    {
      coords::Cartesian_Point const acceleration(coordobj.g_xyz(i)*md::negconvert / M[i]);
      V[i] += acceleration*dt_2;
	  if (Config::get().md.set_active_center == 1)   //adjustment of velocities by distance to active center
	  {
		  V[i] = adjust_velocities(static_cast<int>(i), inner_cutoff, outer_cutoff);
	  }
    }
	if (Config::get().general.verbosity > 3 && Config::get().md.set_active_center == 1)
	{
		std::cout << "number of atoms around active site: " << inner_atoms.size() << "\n";
	}
    // Apply full step RATTLE constraints
    if (CONFIG.rattle.use) rattle_post();
    // Apply full step temperature adjustments
	if (CONFIG.hooverHeatBath || HEATED)
	{
		temp = tempcontrol(CONFIG.hooverHeatBath, false);
	}
	else  // calculate E_kin and T if no temperature control is active
	{
		double tempfactor(2.0 / (freedom*md::R));
		updateEkin();
		temp = E_kin * tempfactor;
	}
    // Apply pressure adjustments
    if (CONFIG.pressure)
    {
      berendsen(dt);
    }
    // save temperature for FEP
    if (Config::get().md.fep)
    {
      coordobj.fep.fepdata.back().T = temp;
    }
    // if requested remove translation and rotation of the system
    if (Config::get().md.veloScale) tune_momentum();

    // Logging / Traces

    if (CONFIG.track)
    {
      std::vector<coords::float_type> iae;
      if (coordobj.interactions().size() > 1)
      {
        iae.reserve(coordobj.interactions().size());
        for (auto const & ia : coordobj.interactions()) iae.push_back(ia.energy);
      }
      logging(k, temp, press, E_kin, coordobj.pes().energy, iae, coordobj.xyz());
    }

    // Serialize to binary file if required.
    if (k > 0 && Config::get().md.restart_offset > 0 && k % Config::get().md.restart_offset == 0)
    {
      write_restartfile(k);
    }
    // add up pressure value
    p_average += press;
  }
  // calculate average pressure over whle simulation time
    p_average /= CONFIG.num_steps;
  if (Config::get().general.verbosity > 2U)
  {
    std::cout << "Average pressure: " << p_average << std::endl;
    //auto integration_time = integration_timer();
    std::cout << "Velocity-Verlet integration took " << integration_timer << '\n';
  }
}

void md::simulation::beemanintegrator(bool fep, std::size_t k_init)
{
	scon::chrono::high_resolution_timer integration_timer;

	config::molecular_dynamics const & CONFIG(Config::get().md);

	std::vector<coords::Cartesian_Point> F_old;

	std::size_t const N = this->coordobj.size();
	// set average pressure to zero
	double p_average(0.0);

	//inner and outer cutoff for biased potential
    double inner_cutoff = Config::get().md.inner_cutoff;
    double outer_cutoff = Config::get().md.outer_cutoff;

	if (Config::get().general.verbosity > 0U)
	{
		std::cout << "Saving " << std::size_t(snapGap > 0 ? (CONFIG.num_steps - k_init) / snapGap : 0);
		std::cout << " snapshots (" << Config::get().md.num_snapShots << " in config)\n";
	}
	// Main MD Loop
	auto split = std::max(std::min(std::size_t(CONFIG.num_steps / 100u), size_t(10000u)), std::size_t{ 100u });
	for (std::size_t k(k_init); k < CONFIG.num_steps; ++k)
	{
		if (k == 0)    // set F(t-dt) for first step to F(t)
		{
			for (size_t i = 0u; i < N; ++i)
			{
				F_old.push_back(coordobj.g_xyz(i));
			}
		}

		bool const HEATED(heat(k,fep));
		if (Config::get().general.verbosity > 1u && k % split == 0 && k > 1)
		{
			std::cout << k << " of " << CONFIG.num_steps << " steps completed\n";
		}

		// apply half step temperature corrections
		if (CONFIG.hooverHeatBath || HEATED)
		{
			temp = tempcontrol(CONFIG.hooverHeatBath, true);
		}

		// save old coordinates
		P_old = coordobj.xyz();

		// Calculate new positions and half step velocities
		for (auto i: atoms_movable)
		{
			// calc acceleration
			coords::Cartesian_Point const acceleration(coordobj.g_xyz(i)*md::negconvert / M[i]);
			coords::Cartesian_Point const acceleration_old(F_old[i] * md::negconvert / M[i]);

			// update veloctiy
			V[i] += acceleration*(2.0 / 3.0)*dt - acceleration_old*(1.0 / 6.0)*dt;
			
			inner_atoms.clear();
			if (Config::get().md.set_active_center == 1)  //adjustment of velocities by distance to active center
			{
				V[i] = adjust_velocities(static_cast<int>(i), inner_cutoff, outer_cutoff);
			}

			if (Config::get().general.verbosity > 4)
			{
				std::cout << "Move " << i << " by " << (V[i] * dt)
					<< " with g " << coordobj.g_xyz(i) << ", V: " << V[i] << std::endl;
			}
			// update coordinates
			coordobj.move_atom_by(i, V[i] * dt);
		}
		if (coordobj.validate_bonds() == false)  // look if all bonds are okay and save those which aren't 
		{
			if (Config::get().general.verbosity > 1U)
			{                                            // give warning if there are broken bonds
				std::cout << "Warning! Broken bonds between atoms...\n";
				for (auto b : coordobj.broken_bonds)
				{
					std::cout << b[0] << " and " << b[1] << ", distance: " << b[2] << "\n";
				}
				if (Config::get().md.broken_restart == 1)
				{         // if desired: set simulation to original positions and random velocities
					restart_broken();     
				}
			}
		}
		if (Config::get().md.set_active_center == 1 && Config::get().md.adjustment_by_step == 1)
		{
			distances = init_active_center(static_cast<int>(k)); 
			atoms_movable.clear();            // determine again which atoms are moved
			for (int i(0U); i < N; ++i)
			{
				if (distances[i] <= outer_cutoff)
				{
					atoms_movable.push_back(i);
				}
				else
				{
					V[i] = coords::Cartesian_Point(0, 0, 0);
				}
			}
		}
		// Apply first part of RATTLE constraints if requested
		if (CONFIG.rattle.use) rattle_pre();

		for (size_t i = 0u; i < N; ++i)   // save F(t) as F_old
		{
			F_old[i] = coordobj.g_xyz(i);
		}
		// calculate new energy & gradients -> F(t+dt)
		coordobj.g();

		// Apply umbrella potential if umbrella sampling is used
		if (CONFIG.umbrella == true)
		{
			coordobj.ubias(udatacontainer);
		}
		// refine nonbondeds if refinement is required due to configuration
		if (CONFIG.refine_offset != 0 && (k + 1U) % CONFIG.refine_offset == 0)
		{
			if (Config::get().general.verbosity > 3U) std::cout << "Refining structure/nonbondeds.\n";
			coordobj.energy_update(true);
		}
		// If spherical boundaries are used apply boundary potential
		boundary_adjustments();

		// add new acceleration and calculate full step velocities
		inner_atoms.clear();
		for (auto i: atoms_movable)
		{
			coords::Cartesian_Point const acceleration_new(coordobj.g_xyz(i)*md::negconvert / M[i]);
			coords::Cartesian_Point const acceleration(F_old[i] * md::negconvert / M[i]);
			V[i] += acceleration_new*(1.0 / 3.0)*dt + acceleration*(1.0 / 6.0)*dt;
			if (Config::get().md.set_active_center == 1)   //adjustment of velocities by distance to active center
			{
				V[i] = adjust_velocities(static_cast<int>(i), inner_cutoff, outer_cutoff);
			}
		}
		if (Config::get().general.verbosity > 3 && Config::get().md.set_active_center == 1)
		{
			std::cout << "number of atoms around active site: " << inner_atoms.size() << "\n";
		}

		// Apply full step RATTLE constraints
		if (CONFIG.rattle.use) rattle_post();
		// Apply full step temperature adjustments
		if (CONFIG.hooverHeatBath || HEATED)
		{
			temp = tempcontrol(CONFIG.hooverHeatBath, false);
		}
		else  // calculate E_kin and T if no temperature control is active
		{
			double tempfactor(2.0 / (freedom*md::R));
			updateEkin();
			temp = E_kin * tempfactor;
		}
		// Apply pressure adjustments
		if (CONFIG.pressure)
		{
			berendsen(dt);
		}
		// save temperature for FEP
		if (Config::get().md.fep)
		{
			coordobj.fep.fepdata.back().T = temp;
		}
		// if requested remove translation and rotation of the system
		if (Config::get().md.veloScale) tune_momentum();

		// Logging / Traces

		if (CONFIG.track)
		{
			std::vector<coords::float_type> iae;
			if (coordobj.interactions().size() > 1)
			{
				iae.reserve(coordobj.interactions().size());
				for (auto const & ia : coordobj.interactions()) iae.push_back(ia.energy);
			}
			logging(k, temp, press, E_kin, coordobj.pes().energy, iae, coordobj.xyz());
		}

			// Serialize to binary file if required.
			if (k > 0 && Config::get().md.restart_offset > 0 && k % Config::get().md.restart_offset == 0)
			{
				write_restartfile(k);
			}
			// add up pressure value
			p_average += press;
		}
		// calculate average pressure over whle simulation time
		p_average /= CONFIG.num_steps;
		if (Config::get().general.verbosity > 2U)
		{
			std::cout << "Average pressure: " << p_average << std::endl;
			//auto integration_time = integration_timer();
			std::cout << "Beeman integration took " << integration_timer << '\n';
		}
	}


 //First part of the RATTLE algorithm to constrain H-X bonds ( half step)
void md::simulation::rattle_pre(void)
{
  std::size_t niter(0U);
  bool done = false;
  // main loop till convergence is reached
  do
  {
    niter += 1;
    done = true;
    // loop over all constraint bonds
    for (auto const & rattlebond : rattle_bonds)
    {
      // get bond distance
      coords::Cartesian_Point const d(coordobj.xyz(rattlebond.b) - coordobj.xyz(rattlebond.a));
      // difference of distance square to optimal distance square
      double const delta = rattlebond.len*rattlebond.len - dot(d, d);
      // if delta is greater than tolerance -> rattle
      if (std::abs(delta) > Config::get().md.rattle.tolerance)
      {
        done = false;
        coords::Cartesian_Point const d_old(P_old[rattlebond.b] - P_old[rattlebond.a]);
        double inv_ma = 1.0 / M[rattlebond.a];
        double inv_mb = 1.0 / M[rattlebond.b];
        //calculate lagrange multiplier
        coords::Cartesian_Point rattle(1.25*delta / (2.0 * (inv_ma + inv_mb) * dot(d, d_old)));
        rattle *= d_old;
        // update half step positions
        coordobj.move_atom_by(rattlebond.a, -(rattle*inv_ma));
        coordobj.move_atom_by(rattlebond.b, rattle*inv_mb);
        inv_ma /= Config::get().md.timeStep;
        inv_mb /= Config::get().md.timeStep;
        // update half step velocities
        V[rattlebond.a] -= rattle*inv_ma;
        V[rattlebond.b] += rattle*inv_mb;
      }
    }
  } while (niter < Config::get().md.rattle.num_iter && done == false);
}
 //second part of RATTLE to constrain H.X bonds (full step)
void md::simulation::rattle_post(void)
{
  // initialize some local variables neede for the internal virial tensor
  std::size_t niter(0U);
  bool done = false;
  double tfact = 2.0 / (dt*1.0e-3*convert);
  double xvir, yvir, zvir;
  double vxx, vyx, vzx, vyy, vzy, vzz;
  std::array<std::array<double, 3>, 3> part_v;
  part_v[0][0] = part_v[0][1] = part_v[0][2] = part_v[1][0] = part_v[1][1] = part_v[1][2] = part_v[2][0] = part_v[2][1] = part_v[2][2] = 0.0;
  // main loop till convergence is reached
  do
  {
    niter += 1;
    done = true;
    // loop over all constraint bonds
    for (auto const & rattlebond : rattle_bonds)
    {
      // get bond distance
      coords::Cartesian_Point const d(coordobj.xyz(rattlebond.b) - coordobj.xyz(rattlebond.a));
      double inv_ma = 1.0 / M[rattlebond.a];
      double inv_mb = 1.0 / M[rattlebond.b];
      // calculate lagrange multiplier
      double const lagrange = -dot(V[rattlebond.b] - V[rattlebond.a], d) / (rattlebond.len*rattlebond.len*(inv_ma + inv_mb));
      if (std::fabs(lagrange) > Config::get().md.rattle.tolerance)
      {
        done = false;
        coords::Cartesian_Point rattle(d*lagrange*1.25);
        // update full step velocities
        V[rattlebond.a] -= rattle*inv_ma;
        V[rattlebond.b] += rattle*inv_mb;
        // compute RATTLE contributions to the internal virial tensor
        xvir = rattle.x() * tfact;
        yvir = rattle.y() * tfact;
        zvir = rattle.z() * tfact;
        vxx = d.x() * xvir;
        vyx = d.y() * xvir;
        vzx = d.z() * xvir;
        vyy = d.y() * yvir;
        vzy = d.z() * yvir;
        vzz = d.z() * zvir;
        part_v[0][0] -= vxx;
        part_v[1][0] -= vyx;
        part_v[2][0] -= vzx;
        part_v[0][1] -= vyx;
        part_v[1][1] -= vyy;
        part_v[2][1] -= vzy;
        part_v[0][2] -= vzx;
        part_v[1][2] -= vzy;
        part_v[2][2] -= vzz;
      }
    }
  } while (niter < Config::get().md.rattle.num_iter && done == false);
  // Add RATTLE contributions to the internal virial tensor
  coordobj.add_to_virial(part_v);
}


void md::simulation::write_restartfile(std::size_t const k)
{
  // Create binary buffer and copy data into vector char
  scon::binary_stream<std::vector<char>> buffer;
  buffer.v.reserve(scon::binary_size(k, *this));
  buffer << k << *this;
  // Create ofstream and insert buffer into stream
  std::ofstream restart_stream(
    (Config::get().general.outputFilename + "_MD_restart.cbf").c_str(),
    std::ofstream::out | std::ofstream::binary
    );
  restart_stream.write(buffer.v.data(), buffer.v.size());
}

//double md::barostat::berendsen::operator()(double const time, 
//  coords::Representation_3D & p,
//  coords::Tensor const & Ek_T, 
//  coords::Tensor const & Vir_T, 
//  coords::Cartesian_Point & box)
//{
//
//  double const volume = std::abs(box.x()) *  std::abs(box.y()) *  std::abs(box.z());
//  
//  double const fac = presc / volume;
//
//  double press = 0.0;
//
//  // ISOTROPIC BOX
//  if (isotropic == true) 
//  {
//    press = fac * ( (2.0 * Ek_T[0][0] - Vir_T[0][0]) + 
//                    (2.0 * Ek_T[1][1] - Vir_T[1][1]) + 
//                    (2.0 * Ek_T[2][2] - Vir_T[2][2])   ) / 3.0;
//    // Berendsen scaling for isotpropic boxes
//    double const scale = std::cbrt(1.0 + (time * compress / delay) * (press - target));
//    // Adjust box dimensions
//    box *= scale;
//    // scale atomic coordinates
//    p *= scale;
//  }
//  // ANISOTROPIC BOX
//  else 
//  {
//    coords::Tensor ptensor, aniso, anisobox, dimtemp;
//    // get pressure tensor for anisotropic systems
//    for (int i = 0; i < 3; i++) 
//    {
//      for (int j = 0; j < 3; j++) 
//      {
//        ptensor[j][i] = fac * (2.0*Ek_T[j][i] - Vir_T[j][i]);
//        //std::cout << "Ptensor:" << ptensor[j][i] << std::endl;
//      }
//    }
//    // get isotropic pressure value
//    press = (ptensor[0][0] + ptensor[1][1] + ptensor[2][2]) / 3.0;
//    // Anisotropic scaling factors
//    double const scale = time * compress / (3.0 * delay);
//    for (int i = 0; i < 3; i++) 
//    {
//      for (int j = 0; j < 3; j++) 
//      {
//        if (j == i)
//        {
//          aniso[i][i] = 1.0 + scale * (ptensor[i][i] - target);
//        }
//        else
//        {
//          aniso[j][i] = scale * ptensor[j][i];
//        }
//      }
//    }
//    // modify anisotropic box dimensions (only for cube or octahedron)
//    dimtemp[0][0] = box.x();
//    dimtemp[1][0] = 0.0;
//    dimtemp[2][0] = 0.0;
//    dimtemp[0][1] = box.y();
//    dimtemp[1][1] = 0.0;
//    dimtemp[2][1] = 0.0;
//    dimtemp[0][2] = 0.0;
//    dimtemp[1][2] = 0.0;
//    dimtemp[2][2] = box.z();
//
//    for (int i = 0; i < 3; i++) 
//    {
//      for (int j = 0; j < 3; j++) 
//      {
//        anisobox[j][i] = 0.0;
//        for (int k = 0; k < 3; k++) 
//        {
//          anisobox[j][i] = anisobox[j][i] + aniso[j][k] * dimtemp[k][i];
//        }
//      }
//    }
//
//    // scale box dimensions for anisotropic case
//
//    box.x() = std::sqrt((anisobox[0][0] * anisobox[0][0] + 
//      anisobox[1][0] * anisobox[1][0] + 
//      anisobox[2][0] * anisobox[2][0]));
//
//    box.y() = std::sqrt((anisobox[0][1] * anisobox[0][1] + 
//      anisobox[1][1] * anisobox[1][1] + 
//      anisobox[2][1] * anisobox[2][1]));
//
//    box.z() = std::sqrt((anisobox[0][2] * anisobox[0][2] + 
//      anisobox[1][2] * anisobox[1][2] + 
//      anisobox[2][2] * anisobox[2][2]));
//
//    // scale atomic coordinates for anisotropic case
//    for (auto & pos : p)
//    {
//      auto px = aniso[0][0] * pos.x() + aniso[0][1] * pos.y() + aniso[0][2] * pos.z();
//      //pos.x() = px;
//      auto py = aniso[1][0] * pos.x() + aniso[1][1] * pos.y() + aniso[1][2] * pos.z();
//      //pos.y() = py;
//      auto pz = aniso[2][0] * pos.x() + aniso[2][1] * pos.y() + aniso[2][2] * pos.z();
//      pos.x() = px;
//      pos.y() = py;
//      pos.z() = pz;
//    }
//
//  }//end of isotropic else
//  return press;
//}
