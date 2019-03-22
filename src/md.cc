#include "md.h"

// Logging Functions

std::ostream& md::operator<<(std::ostream &strm, trace_data const &d)
{
  strm  << d.i<<",";
  strm  << std::fixed << std::setprecision(3) << d.T << ",";
  strm << std::fixed << std::setprecision(5) << d.P << ",";
  strm   << std::fixed << std::setprecision(5) << d.Ek << ",";
  strm   << std::fixed << std::setprecision(5) << d.Ep << ",";
  strm   << std::fixed << std::setprecision(5) << d.Ek + d.Ep << ",";
  if (d.snapshot > 0u) strm <<  d.snapshot;
  else strm << std::right << std::setw(20) << '-';
  for (auto iae : d.Eia)
  {
    strm << "," << std::fixed << std::setprecision(5) << iae;
  }
  strm << '\n';
  return strm;
}

void md::trace_writer::operator() (md::trace_data const & d)
{
  static std::atomic<bool> x(true);
  if (x && Config::get().general.verbosity > 1)
  {
    *strm  << "It,";
    *strm << "T,";
    *strm  << "P,";
    *strm  << "kin.En.,";
    *strm  << "pot.En.,";
    *strm  << "tot.En.,";
    *strm  << "Snapsh.";
    if (d.Eia.size() > 1u)
    {
      for (auto i : scon::index_range(d.Eia))
      {
        *strm << "," << "E_ia(# " << i << ')';
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
    trace_writer{ coords::output::filename("_MD_TRACE", ".csv").c_str() })),
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
  M_total(0.0), E_kin(0.0), T(Config::get().md.T_init), temp(0.0), press(0.0), dt(Config::get().md.timeStep),
  freedom(0), snapGap(0), C_geo(), C_mass(), 
  nht(), rattle_bonds(), window(), restarted(true)
{
  std::sort(Config::set().md.heat_steps.begin(), Config::set().md.heat_steps.end());

  if (Config::get().md.ana_pairs.size() > 0)
  {  
    for (auto p : Config::get().md.ana_pairs)
    {                           // create atom pairs to analyze, fetch information and save
      ana_pair ap(p[0], p[1]);
      ap.symbol_a = coordobj.atoms(ap.a).symbol();
      ap.symbol_b = coordobj.atoms(ap.b).symbol();
      ap.name_a = ap.symbol_a + std::to_string(ap.a + 1);
      ap.name_b = ap.symbol_b + std::to_string(ap.b + 1);
      ap.legend = ap.name_a + "-" + ap.name_b;
      ana_pairs.push_back(ap);
    }
  }
}


void md::simulation::run(bool const restart)
{
  // Print initial info if required due to verbosity
  if (Config::get().general.verbosity > 1U) 
    print_init_info();
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
    if (bs >> iteration)
      bs >> *this;
  }
  // start Simulation
  integrate(false, iteration);
}

// Perform Umbrella Sampling run if requested
void md::simulation::umbrella_run(bool const restart) {
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
    tune_momentum(); // eliminate translation and rotation
  }
  // Set kinetic Energy
  updateEkin();
  //run equilibration
  Config::set().md.num_steps = Config::get().md.usequil;
  integrate(false);
  udatacontainer.clear();
  // run production
  Config::set().md.num_steps = steps;
  integrate(false);

  // write data into file
  std::ofstream ofs;
  ofs.open("umbrella.txt");
  auto &&number_of_restraints = Config::get().coords.bias.udist.size() + Config::get().coords.bias.utors.size() + Config::get().coords.bias.ucombs.size();
  for (auto s{ 0u }; s < udatacontainer.size()/number_of_restraints; ++s)  // for every step 
  {
    ofs << s << "   ";                                              // stepnumber
    for (auto b{ 0u }; b < number_of_restraints; ++b) {
      ofs << udatacontainer[b + number_of_restraints * s] << "  ";  // value(s)
    }
    ofs << "\n";
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
    integrator(fep, k_init, true);
    break;
  }
  default:
  { // Velocity verlet integrator
    integrator(fep, k_init, false);
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
      std::cout << "The following distances will be fixed: \n";
      for (auto const & bond : Config::get().md.rattle.specified_rattle)
      {
        std::cout << "[" << bond.a+1 << ":" << bond.b+1 << "] ";
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

// check if the two atom of a rattlepair are bonded with each other
void md::simulation::check_rattlepair_for_bond(config::md_conf::config_rattle::rattle_constraint_bond &rctemp)
{
	bool is_a_bond = false;
	for (auto bonding_partner : coordobj.atoms(rctemp.a).bonds())
	{
		if (bonding_partner == rctemp.b)
		{
			is_a_bond = true;
			break;
		}
	}
	if (is_a_bond == false)
	{
		throw std::runtime_error("No bond between " + std::to_string(rctemp.a + 1) + " and " + std::to_string(rctemp.b + 1) + ". It doesn't make sense to use a bonding parameter.");
	}
}

// set up constraints for H-X bonds if requested
// ideal bond lengths are taken from the foce field parameter file
// specified by RATTpar in the INPUTFILE
void md::simulation::rattlesetup(void)
{
  config::md_conf::config_rattle::rattle_constraint_bond rctemp;  // temporary rattlebond

	// temp vars and vectors;
	struct trat {
		std::size_t ia, ib;
		double ideal;
	};
	struct ratatoms {
		std::size_t ia, ga;
	};
	std::vector<trat> temprat;
	std::vector<ratatoms> ratoms;
	std::size_t tia = std::size_t(), tib = std::size_t();
  
	if (Config::get().md.rattle.use_paramfile)    //open rattle par file and look for atomtypes and distances
	{
		std::ifstream ifs;
		ifs.open(Config::get().general.paramFilename.c_str());
		if (ifs.good()) std::cout << "Opened file for RATTLE parameters successfully." << std::endl;
		if (!ifs.good()) {
			std::cout << "Couldn't open file for RATTLE parameters. Check your input" << std::endl;
			throw;
		}
		
		// more temporary variables
		ratatoms RTA;
		trat TRAT;
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
		ifs.open(Config::get().general.paramFilename.c_str());
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
	}

  // Generate vector with bonds which are to be constraint
  const std::size_t N = coordobj.size();

	if (Config::get().md.rattle.all == true)  // all H-bonds are constraint
	{
		for (std::size_t i = 0; i < N; i++) //loop over all atoms
		{
			if (coordobj.atoms(i).number() == 1) //check if atom is hydrogen
			{
				rctemp.a = i;                            // hydrogen atom
				rctemp.b = coordobj.atoms(i).bonds(0);   // bonding partner of hydrogen atom
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
	else   // if MDrattle = 2 i.e. only specified distances are constrained
	{
		if (Config::get().md.rattle.use_paramfile == false && Config::get().md.rattle.dists.size() != Config::get().md.rattle.specified_rattle.size())
		{
			throw std::runtime_error("Wrong number of atom distances given for rattlepairs!");
		}

		for (auto i=0u; i<Config::get().md.rattle.specified_rattle.size(); ++i)
		{
			rctemp = Config::get().md.rattle.specified_rattle[i];

			if (Config::get().md.rattle.use_paramfile)   // if distances are to be taken from parameterfile
			{
				// check if atoms share a bond
				check_rattlepair_for_bond(rctemp);

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
						break;
					}
				}
			}

			else   // get rattledist from inputfile
			{
				rctemp.len = Config::get().md.rattle.dists[i];
			}
			rattle_bonds.push_back(rctemp);
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
  // things for biased potentials
  inner_atoms.clear();    // in case of more than one MD, e.g. for an FEP calculation
  movable_atoms.clear();
  if (Config::get().md.set_active_center == 1)
  {
    distances = init_active_center(0);   //calculate initial active center and distances to active center

    for (auto i(0U); i < N; ++i)  // determine which atoms are moved
    {
      if (distances[i] <= Config::get().md.outer_cutoff)
      {
        movable_atoms.push_back(i);
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
  else if (Config::get().coords.fixed.size() != 0)  // if atoms are fixed 
  {
    for (auto i(0U); i < N; ++i)  // determine which atoms are moved
    {
      if (is_in(i, Config::get().coords.fixed) == false)
      {
        movable_atoms.push_back(i);
      }
    }
  }
  else   // if no active site is specified: all atoms are moved
  {
    for (auto i(0U); i < N; ++i)
    {
      movable_atoms.push_back(i);
    }
  }

  // temperature stuff
  M_total = 0;
  V.resize(N);
  F.resize(N);
  M.resize(N);
  P_start = coordobj.xyz();

  std::default_random_engine generator(static_cast<unsigned> (time(0)));  // generates random numbers
  auto dist01 = std::normal_distribution<double>{ 0,1 }; // normal distribution with mean=0 and standard deviation=1
  for (std::size_t i = 0; i < N; ++i)
  {
    // Get Atom Mass
    M[i] = coordobj.atoms(i).mass();
    // Sum total Mass
    M_total += M[i];
    // Set initial velocities to zero if fixed or not movable
    if (coordobj.atoms(i).fixed() || !(abs(Config::get().md.T_init) > 0.0) || std::find(movable_atoms.begin(), movable_atoms.end(), i) == movable_atoms.end())
    {
      V[i] = coords::Cartesian_Point(0);
    }
    // initialize random velocities otherwise (Maxwell Boltzmann distribution, see http://research.chem.psu.edu/shsgroup/chem647/newNotes/node6.html)
    else
    {
      V[i].x() = dist01(generator) * std::sqrt(kB*T / M[i]);
      V[i].y() = dist01(generator) * std::sqrt(kB*T / M[i]);
      V[i].z() = dist01(generator) * std::sqrt(kB*T / M[i]);

      if (Config::get().general.verbosity > 4U) std::cout << "Initial Velocity of " << i << " is " << V[i] << "\n";
    }
    // sum position vectors for geometrical center
    C_geo += coordobj.xyz(i);
  }
  // get degrees of freedom
  freedom = 3U * N;

  if (Config::get().md.rattle.use == true)
  {
    rattlesetup();                   // Set up rattle vector for constraints
		freedom -= rattle_bonds.size();  // constraint degrees of freedom
  }
  
  // periodics and isothermal cases
  if (Config::get().md.hooverHeatBath == true)
  {
    if (Config::get().periodics.periodic == true)
      freedom -= 3;
    else
      freedom -= 6;
  }
  if (Config::get().general.verbosity > 2U)
    std::cout << "Degrees of freedom: " << freedom << std::endl;

  // calc number of steps between snapshots (gap between snapshots)
  snapGap = gap(Config::get().md.num_steps, Config::get().md.num_snapShots);
  // scale geometrical center vector by number of atoms
  C_geo /= static_cast<double>(N);
  // call center of Mass method from coordinates object
  C_mass = coordobj.center_of_mass();

  if (Config::get().md.analyze_zones == true)
  {
    zones = find_zones();  // find atoms for every zone
  }
}

/**function that fills zones with atoms*/
std::vector<md::zone> md::simulation::find_zones()
{
  std::vector<zone> zones;
  double zone_width = Config::get().md.zone_width;

  // calculate initial active site and distances to active center
  distances = init_active_center(0);   

  // find out number of zones (maximum distance to active site)
  std::vector<double>::iterator it = std::max_element(distances.begin(), distances.end());
  double max_dist = *it;
  int number_of_zones = std::ceil(max_dist / zone_width);
  zones.resize(number_of_zones);

  // create zones, fill them with atoms and create legend
  for (auto i = 0u; i < distances.size(); i++)
  {
    int zone = std::floor(distances[i] / zone_width);
    zones[zone].atoms.push_back(i);
  }
  for (auto i = 0u; i < zones.size(); i++)
  {
    zones[i].legend = std::to_string(int(i * zone_width)) + " to " + std::to_string(int((i + 1)*zone_width)); 
  }

  // output
  if (Config::get().general.verbosity > 2)
  {
    std::cout << "Zones:\n";
    for (auto i = 0u; i < zones.size(); i++)
    {
      std::cout << zones[i].atoms.size() << " atoms in zone from " << i * zone_width << " to " << (i + 1)*zone_width << "\n";
    }
  }

  return zones;
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
  double  increment;
  Config::set().md.fep = true;
  FEPsum = 0.0;
  FEPsum_back = 0.0;
  FEPsum_SOS = 0.0;
  FEPsum_BAR = 0.0;
  // fill vector with scaling increments
  increment = Config::get().fep.lambda / Config::get().fep.dlambda;
  //std::cout << Config::get().fep.lambda << "   " << Config::get().fep.dlambda << std::endl;
  std::cout << "Number of FEP windows:  " << increment << std::endl;
  coordobj.fep.window.resize(std::size_t(increment)+1); //because of backward transformation one window more necessary
  coordobj.fep.window[0].step = 0;

  // calculate all lambda values for every window
  for (auto i = 0u; i < coordobj.fep.window.size(); i++) {

    double lambda = i * Config::get().fep.dlambda;  // lambda
    if (lambda < Config::get().fep.eleccouple) coordobj.fep.window[i].ein = 0;
    else coordobj.fep.window[i].ein = 1 - (1-lambda)/(1- Config::get().fep.eleccouple);
    if (lambda > Config::get().fep.vdwcouple) coordobj.fep.window[i].vin = 1;
    else coordobj.fep.window[i].vin = lambda / Config::get().fep.vdwcouple;
    if (lambda > 1-Config::get().fep.eleccouple) coordobj.fep.window[i].eout = 0;
    else coordobj.fep.window[i].eout = 1 - (lambda) / (1 - Config::get().fep.eleccouple);
    if (lambda < 1-Config::get().fep.vdwcouple) coordobj.fep.window[i].vout = 1;
    else coordobj.fep.window[i].vout = (1-lambda) / Config::get().fep.vdwcouple;

    double dlambda = (i+1) * Config::get().fep.dlambda;  // lambda + dlambda
    if (dlambda < Config::get().fep.eleccouple) coordobj.fep.window[i].dein = 0;
    else if (dlambda > 1) coordobj.fep.window[i].dein = 1;
    else coordobj.fep.window[i].dein = 1 - (1 - dlambda) / (1 - Config::get().fep.eleccouple);
    if (dlambda > Config::get().fep.vdwcouple) coordobj.fep.window[i].dvin = 1;
    else coordobj.fep.window[i].dvin = dlambda / Config::get().fep.vdwcouple;
    if (dlambda > 1 - Config::get().fep.eleccouple) coordobj.fep.window[i].deout = 0;
    else coordobj.fep.window[i].deout = 1 - (dlambda) / (1 - Config::get().fep.eleccouple);
    if (dlambda < 1-Config::get().fep.vdwcouple) coordobj.fep.window[i].dvout = 1;
    else if (dlambda > 1) coordobj.fep.window[i].dvout = 0;
    else coordobj.fep.window[i].dvout = (1 - dlambda) / Config::get().fep.vdwcouple;

    double mlambda = (int)(i-1) * Config::get().fep.dlambda;  // lambda - dlambda
    if (mlambda < Config::get().fep.eleccouple) coordobj.fep.window[i].mein = 0;
    else coordobj.fep.window[i].mein = 1 - (1 - mlambda) / (1 - Config::get().fep.eleccouple);
    if (mlambda > Config::get().fep.vdwcouple) coordobj.fep.window[i].mvin = 1;
    else if (mlambda < 0) coordobj.fep.window[i].mvin = 0;
    else coordobj.fep.window[i].mvin = mlambda / Config::get().fep.vdwcouple;
    if (mlambda > 1 - Config::get().fep.eleccouple) coordobj.fep.window[i].meout = 0;
    else if (mlambda < 0) coordobj.fep.window[i].meout = 1;
    else coordobj.fep.window[i].meout = 1 - (mlambda) / (1 - Config::get().fep.eleccouple);
    if (mlambda < 1-Config::get().fep.vdwcouple) coordobj.fep.window[i].mvout = 1;
    else coordobj.fep.window[i].mvout = (1 - mlambda) / Config::get().fep.vdwcouple;
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
    std::cout << std::setprecision(4) << std::setw(8) << coordobj.fep.window[i].mvout << std::setw(8) << coordobj.fep.window[i].vout << std::setw(8) << coordobj.fep.window[i].dvout << std::setw(8) << coordobj.fep.window[i].mvin << std::setw(8) << coordobj.fep.window[i].vin << std::setw(8) << coordobj.fep.window[i].dvin << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Electrostatic Coupling:" << std::endl;
  for (std::size_t i = 0; i < coordobj.fep.window.size(); i++)
  {
    std::cout << std::setprecision(4) << std::setw(8) << coordobj.fep.window[i].meout << std::setw(8) << coordobj.fep.window[i].eout << std::setw(8) << coordobj.fep.window[i].deout << std::setw(8) << coordobj.fep.window[i].mein << std::setw(8) << coordobj.fep.window[i].ein << std::setw(8) << coordobj.fep.window[i].dein << std::endl;
  }

}

// Calculation of ensemble average and free energy change for each step if FEP calculation is performed
// calculation can be improved if at every step the current averages are stored
// currently calculation is performed at the end of each window
void md::simulation::freecalc()
{
  std::size_t iterator(0U), k(0U);
  // set conversion factors (conv) and constants (boltzmann, avogadro)
  double de_ensemble, de_ensemble_back, de_ensemble_half=0, de_ensemble_back_half=0, temp_avg, boltz = 1.3806488E-23, avogad = 6.022E23, conv = 4184.0;
  // calculate ensemble average for current window
  for (std::size_t i = 0; i < coordobj.fep.fepdata.size(); i++)
  {                // for every conformation in window
    iterator += 1;
    k = 0;
    de_ensemble = de_ensemble_back = temp_avg = 0.0;
    double exponent = -1 / (boltz*coordobj.fep.fepdata[i].T)*conv*coordobj.fep.fepdata[i].dE / avogad;
    double exponent_back = 1 / (boltz*coordobj.fep.fepdata[i].T)*conv*coordobj.fep.fepdata[i].dE_back / avogad;
    coordobj.fep.fepdata[i].de_ens = exp(exponent);
    coordobj.fep.fepdata[i].de_ens_back = exp(exponent_back);
    double de_ens_half = exp(exponent/2);
    double de_ens_back_half = exp(exponent_back/2);
    for (k = 0; k <= i; k++) {
      temp_avg += coordobj.fep.fepdata[k].T;
      de_ensemble += coordobj.fep.fepdata[k].de_ens;
      de_ensemble_back += coordobj.fep.fepdata[k].de_ens_back;
    }
    de_ensemble = de_ensemble / iterator;
    de_ensemble_back = de_ensemble_back / iterator;
    de_ensemble_half += de_ens_half;
    de_ensemble_back_half += de_ens_back_half;
    temp_avg = temp_avg / iterator;
    coordobj.fep.fepdata[i].dG = -1 * std::log(de_ensemble)*temp_avg*boltz*avogad / conv;
    coordobj.fep.fepdata[i].dG_back = std::log(de_ensemble_back)*temp_avg*boltz*avogad / conv;
  }// end of main loop

  de_ensemble_half = de_ensemble_half / coordobj.fep.fepdata.size();
  de_ensemble_back_half = de_ensemble_back_half / coordobj.fep.fepdata.size();
  if (de_ensemble_back_half == 1)  dG_SOS = 0;
  else dG_SOS = -1 * std::log(de_ensemble_v_SOS / de_ensemble_back_half)*temp_avg*boltz*avogad / conv;
  de_ensemble_v_SOS = de_ensemble_half; //de_ensemble_half is needed for SOS-calculation in next step

  // calculate final free energy change for the current window
  this->FEPsum += coordobj.fep.fepdata[coordobj.fep.fepdata.size() - 1].dG;
  this->FEPsum_back += coordobj.fep.fepdata[coordobj.fep.fepdata.size() - 1].dG_back;
  this->FEPsum_SOS += dG_SOS;
}

void md::simulation::bar(int current_window)
{
   double boltz = 1.3806488E-23, avogad = 6.022E23, conv = 4184.0;
   double w, w_back;  // weighting function
   double dG_BAR = dG_SOS;  // start value for iteration
   double c; // constant C

   if (Config::get().general.verbosity > 3)
   {
     std::cout << "Start solution of BAR equation from dG_SOS: " << dG_SOS << "\n";
   }

	 auto count_iterations{ 0u };
   do    // iterative solution for BAR equation
   {
     c = dG_BAR;
     double ensemble = 0;
     double ensemble_back = 0;
     double temp_avg = 0;  // average temperature
     for (std::size_t i = 0; i < coordobj.fep.fepdata.size(); i++) // for every conformation in window
     {
       w = 2 / (exp(1 / (boltz*coordobj.fep.fepdata[i].T)*conv*((coordobj.fep.fepdata[i].dE - c) / 2) / avogad) + exp(-1 / (boltz*coordobj.fep.fepdata[i].T)*conv*((coordobj.fep.fepdata[i].dE - c) / 2) / avogad));
       double ens = w * exp(-1 / (boltz*coordobj.fep.fepdata[i].T)*conv*(coordobj.fep.fepdata[i].dE / 2) / avogad);
       w_back = 2 / (exp(1 / (boltz*coordobj.fep.fepdata[i].T)*conv*((coordobj.fep.fepdata[i].dE_back - c) / 2) / avogad) + exp(-1 / (boltz*coordobj.fep.fepdata[i].T)*conv*((coordobj.fep.fepdata[i].dE_back - c) / 2) / avogad));
       double ens_back = w_back * exp(1 / (boltz*coordobj.fep.fepdata[i].T)*conv*(coordobj.fep.fepdata[i].dE_back / 2) / avogad);
       ensemble += ens;
       ensemble_back += ens_back;
       temp_avg += coordobj.fep.fepdata[i].T;
     }
     ensemble = ensemble / coordobj.fep.fepdata.size();       // calculate averages
     ensemble_back = ensemble_back / coordobj.fep.fepdata.size();
     temp_avg = temp_avg / coordobj.fep.fepdata.size();

     if (current_window == 0)  dG_BAR = 0;                            // calculate dG for current window
     else dG_BAR = -1 * std::log(de_ensemble_v_BAR / ensemble_back)*temp_avg*boltz*avogad / conv;
     if (Config::get().general.verbosity > 3)
     {
       std::cout << "dG_BAR: " << dG_BAR << "\n";
     }
     de_ensemble_v_BAR = ensemble;  // this is needed in next step
     
		 count_iterations += 1;
   } while (fabs(c - dG_BAR) > 0.001 && count_iterations < 5000);  // 0.001 = convergence threshold and 5000 = maximum number of iterations (maybe later define by user?)
   this->FEPsum_BAR += dG_BAR;
}

// write the output FEP calculations
void md::simulation::freewrite(int i)
{
  std::ofstream fep("alchemical.txt", std::ios_base::app);
  std::ofstream res("FEP_Results.txt", std::ios_base::app);

  // standard forward output
  if (i*Config::get().fep.dlambda == 0 && this->prod == false)
  {
    res << std::fixed << std::right << std::setprecision(4) << std::setw(10) << "0" << std::setw(10) << "0";
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
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].e_c_l0;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].e_c_l1;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].e_c_l2;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].e_vdw_l0;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].e_vdw_l1;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].e_vdw_l2;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].T;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].dE;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].dG;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].dE_back;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.fep.fepdata[k].dG_back;
      fep << std::endl;
    }
    if (Config::get().general.verbosity > 3u)
    {
      std::cout << "Coulomb: " << coordobj.fep.fepdata[k].e_c_l2 - coordobj.fep.fepdata[k].e_c_l1 << ", vdW: " << coordobj.fep.fepdata[k].e_vdw_l2 - coordobj.fep.fepdata[k].e_vdw_l1 << "\n";
    }
  }

  // at the end of production data in alchemical.txt sum up the results and print the before the new window starts
  if (this->prod == true) {
    fep << "Free energy change for the current window:  ";
    fep << coordobj.fep.fepdata[coordobj.fep.fepdata.size() - 1].dG << std::endl;
    fep << "Total free energy change until current window:  " << FEPsum << std::endl;

    res << std::fixed << std::right << std::setprecision(4) << std::setw(10) << FEPsum_back << std::right << std::setprecision(4) << std::setw(10) << FEPsum_SOS<< std::right << std::setprecision(4) << std::setw(10) << FEPsum_BAR << std::endl;
    double rounded = std::stod(std::to_string(i * Config::get().fep.dlambda)); // round necessary when having a number of windows that is can't be expressed exactly in decimal numbers
    if (rounded < 1) {
      res << std::fixed << std::right << std::setprecision(4) << std::setw(10) << (i * Config::get().fep.dlambda) + Config::get().fep.dlambda << std::setw(10) << FEPsum;
    }
  }
}
#ifdef USE_PYTHON
std::string md::simulation::get_pythonpath()
{
  std::string pythonpaths_str = Py_GetPath();
  std::string path;
#ifdef __unix__
  std::vector<std::string> pythonpaths = split(pythonpaths_str, ':');
#elif defined(_WIN32) || defined(WIN32)
  std::vector<std::string> pythonpaths = split(pythonpaths_str, ';');
#endif
  path = "import sys\n";
  for (auto p : pythonpaths)  //keep pythonpath of system
  {
    path += "sys.path.append('" + p + "')\n";
  }
  path += "sys.path.append('" + get_python_modulepath("matplotlib") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("fractions") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("csv") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("atexit") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("calendar") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("Tkinter") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("FileDialog") + "')\n";
  return path;
}

std::vector<double> md::simulation::fepanalyze(std::vector<double> dE_pots, int window)
{
  std::string add_path = get_pythonpath();

  PyObject *modul, *funk, *prm, *ret, *pValue;

  // create python list with dE_pot_bac
  PyObject *E_pot_backs = PyList_New(coordobj.fep.fepdata.size());
  for (std::size_t k = 0; k < coordobj.fep.fepdata.size(); k++) {
    pValue = PyFloat_FromDouble(coordobj.fep.fepdata[k].dE_back);
    PyList_SetItem(E_pot_backs, k, pValue);
  }

  if (window > 0)   // no output for 0th window
  {
    // create python list with dE_pot from last run
    PyObject *E_pots = PyList_New(coordobj.fep.fepdata.size());   
    for (std::size_t k = 0; k < coordobj.fep.fepdata.size(); k++) {
      pValue = PyFloat_FromDouble(dE_pots[k]);
      PyList_SetItem(E_pots, k, pValue);
    }

    PySys_SetPath((char*)"./python_modules"); //set path
    const char *c = add_path.c_str();  //add paths pythonpath
    PyRun_SimpleString(c);

    modul = PyImport_ImportModule("FEP_analysis"); //import module 
    if (modul)
    {
      funk = PyObject_GetAttrString(modul, "plot_histograms_and_calculate_overlap"); //create function
      prm = Py_BuildValue("OOi", E_pots, E_pot_backs, window); //give parameters
      ret = PyObject_CallObject(funk, prm);  //call function with parameters
      std::string result_str = PyString_AsString(ret); //convert result to a C++ string
      if (result_str == "error")
      {
        std::cout << "An error occured during running python module 'FEP_analysis'\n";
      }
      else  // python function was successfull
      {
        float result = std::stof(result_str);  // convert result to float
        std::ofstream overlap("overlap.txt", std::ios_base::app);
        overlap << "Window " << window << ": " << result * 100 << " %\n";
      }
    }
    else
    {
      std::cout << "Error: module 'FEP_analysis' not found!\n";
      std::exit(0);
    }
    //delete PyObjects
    Py_DECREF(prm);
    Py_DECREF(ret);
    Py_DECREF(funk);
    Py_DECREF(modul);
    Py_DECREF(pValue);
    Py_DECREF(E_pots);
    Py_DECREF(E_pot_backs);
  }

  dE_pots.clear();  // save dE_pot for next run and return them
  for (std::size_t k = 0; k < coordobj.fep.fepdata.size(); k++) {  
    dE_pots.push_back(coordobj.fep.fepdata[k].dE);
  }
  return dE_pots;
}

void md::simulation::plot_distances(std::vector<ana_pair> &pairs)
{
  write_dists_into_file(pairs);

  std::string add_path = get_pythonpath();

  PyObject *modul, *funk, *prm, *ret, *pValue;

  // create python list with legends
  PyObject *legends = PyList_New(pairs.size());
  for (std::size_t k = 0; k < pairs.size(); k++) {
    pValue = PyString_FromString(pairs[k].legend.c_str());
    PyList_SetItem(legends, k, pValue);
  }

  // create a python list that contains a list with distances for every atom pair that is to be analyzed
  PyObject *distance_lists = PyList_New(pairs.size());
  int counter = 0;
  for (auto a : pairs)
  {
    PyObject *dists = PyList_New(a.dists.size());
    for (std::size_t k = 0; k < a.dists.size(); k++) {
      pValue = PyFloat_FromDouble(a.dists[k]);
      PyList_SetItem(dists, k, pValue);
    }
    PyList_SetItem(distance_lists, counter, dists);
    counter += 1;
  }

  PySys_SetPath((char*)"./python_modules"); //set path
  const char *c = add_path.c_str();  //add paths pythonpath
  PyRun_SimpleString(c);

  modul = PyImport_ImportModule("MD_analysis"); //import module 
  if (modul)
  {
    funk = PyObject_GetAttrString(modul, "plot_dists"); //create function
    prm = Py_BuildValue("(OO)", legends, distance_lists); //give parameters
    ret = PyObject_CallObject(funk, prm);  //call function with parameters
    std::string result_str = PyString_AsString(ret); //convert result to a C++ string
    if (result_str == "error")
    {
      std::cout << "An error occured during running python module 'MD_analysis'\n";
    }
  }
  else
  {
    std::cout << "Error: module 'MD_analysis' not found!\n";
    std::exit(0);
  }
  //delete PyObjects
  Py_DECREF(prm);
  Py_DECREF(ret);
  Py_DECREF(funk);
  Py_DECREF(modul);
  Py_DECREF(pValue);
  Py_DECREF(legends);
  Py_DECREF(distance_lists);
}

/**function to plot temperatures for all zones*/
void md::simulation::plot_zones()
{
  write_zones_into_file();

  std::string add_path = get_pythonpath();

  PyObject *modul, *funk, *prm, *ret, *pValue;

  // create python list with legends
  PyObject *legends = PyList_New(zones.size());
  for (std::size_t k = 0; k < zones.size(); k++) {
    pValue = PyString_FromString(zones[k].legend.c_str());
    PyList_SetItem(legends, k, pValue);
  }

  // create a python list that contains a list with temperatures for every zone
  PyObject *temp_lists = PyList_New(zones.size());
  int counter = 0;
  for (auto z : zones)
  {
    PyObject *temps = PyList_New(z.temperatures.size());
    for (std::size_t k = 0; k < z.temperatures.size(); k++) {
      pValue = PyFloat_FromDouble(z.temperatures[k]);
      PyList_SetItem(temps, k, pValue);
    }
    PyList_SetItem(temp_lists, counter, temps);
    counter += 1;
  }

  PySys_SetPath((char*)"./python_modules"); //set path
  const char *c = add_path.c_str();  //add paths pythonpath
  PyRun_SimpleString(c);

  modul = PyImport_ImportModule("MD_analysis"); //import module 
  if (modul)
  {
    funk = PyObject_GetAttrString(modul, "plot_zones"); //create function
    prm = Py_BuildValue("(OO)", legends, temp_lists); //give parameters
    ret = PyObject_CallObject(funk, prm);  //call function with parameters
    std::string result_str = PyString_AsString(ret); //convert result to a C++ string
    if (result_str == "error")
    {
      std::cout << "An error occured during running python module 'MD_analysis'\n";
    }
  }
  else
  {
    std::cout << "Error: module 'MD_analysis' not found!\n";
    std::exit(0);
  }
  //delete PyObjects
  Py_DECREF(prm);
  Py_DECREF(ret);
  Py_DECREF(funk);
  Py_DECREF(modul);
  Py_DECREF(pValue);
  Py_DECREF(legends);
  Py_DECREF(temp_lists);
}

#endif

void md::simulation::write_zones_into_file()
{
  std::ofstream zonefile;
  zonefile.open("zones.csv");

  zonefile << "Steps";                               // write headline
  for (auto &z : zones) zonefile << "," << z.legend;
  zonefile << "\n";

  for (auto i = 0u; i < Config::get().md.num_steps; ++i)   // for every MD step
  {
    zonefile << i + 1;
    for (auto &z : zones) zonefile << "," << z.temperatures[i];  // write a line with temperatures
    zonefile << "\n";
  }
  zonefile.close();
}

void md::simulation::write_dists_into_file(std::vector<ana_pair>& pairs)
{
  std::ofstream distfile;
  distfile.open("distances.csv");

  distfile << "Steps";                               // write headline
  for (auto &p : pairs) distfile << "," << p.legend;
  distfile << "\n";

  for (auto i = 0u; i < Config::get().md.num_steps; ++i)   // for every MD step
  {
    distfile << i + 1;
    for (auto &p : pairs) distfile << "," << p.dists[i];  // write a line with distances
    distfile << "\n";
  }
  distfile.close();

  for (auto &p : pairs) p.dists.clear();   // after writing: delete vector with distances
}

// perform FEP calculation if requested
void md::simulation::feprun()
{
  if (Config::get().fep.analyze)
  {
    std::remove("overlap.txt");
  }
  std::vector<double> dE_pots;

  for (auto i(0U); i < coordobj.fep.window.size(); ++i)  //for every window
  {
    std::cout << "Lambda:  " << i * Config::get().fep.dlambda << "\n";
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
    if (Config::get().fep.bar == true) bar(i);
    freewrite(i);

    if (Config::get().fep.analyze)
    {
#ifdef USE_PYTHON
      dE_pots = fepanalyze(dE_pots, i);
#else
      std::cout << "Analyzing is not possible without python!\n";
#endif
    }
  }
}// end of main window loop


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
  TenVal z = { 0.0,0.0,0.0 };
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
    std::cout << "New kinetic Energy is " <<E_kin<<" with x, y, z = " << E_kin_tensor[0][0] << ", "
      << E_kin_tensor[1][1] << ", " << E_kin_tensor[2][2] << '\n';
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
    std::cout << "New kinetic Energy is " << E_kin << " with x, y, z = " << E_kin_tensor[0][0] << ", "
      << E_kin_tensor[1][1] << ", " << E_kin_tensor[2][2] << '\n';
  }
}

// apply pressure corrections if constant pressure simulation is performed
void md::simulation::berendsen(double const time)
{
  // temp variables
  const std::size_t N = coordobj.size();
  double fac, scale, volume;
  double ptensor[3][3], aniso[3][3], anisobox[3][3], dimtemp[3][3];
  // get volume and pressure scaling factor
  volume = Config::get().periodics.pb_box.x() *  Config::get().periodics.pb_box.y() *  Config::get().periodics.pb_box.z();
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
    Config::set().periodics.pb_box.x() = Config::get().periodics.pb_box.x() * scale;
    Config::set().periodics.pb_box.y() = Config::get().periodics.pb_box.y() * scale;
    Config::set().periodics.pb_box.z() = Config::get().periodics.pb_box.z() * scale;
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
    dimtemp[0][0] = Config::get().periodics.pb_box.x();
    dimtemp[1][0] = 0.0;
    dimtemp[2][0] = 0.0;
    dimtemp[0][1] = Config::get().periodics.pb_box.y() * 1.0;
    dimtemp[1][1] = Config::get().periodics.pb_box.y() * 0.0;
    dimtemp[2][1] = 0.0;
    dimtemp[0][2] = Config::get().periodics.pb_box.z() * 0.0;
    dimtemp[1][2] = Config::get().periodics.pb_box.z() * (1.0 - 1.0 * 1.0) / 1.0;
    dimtemp[2][2] = Config::get().periodics.pb_box.z() * sqrt(1.0 * 1.0 - 0.0 - 0.0);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        anisobox[j][i] = 0.0;
        for (int k = 0; k < 3; k++) {
          anisobox[j][i] = anisobox[j][i] + aniso[j][k] * dimtemp[k][i];
        }
      }
    }
    // scale box dimensions for anisotropic case
    Config::set().periodics.pb_box.x() = std::sqrt((anisobox[0][0] * anisobox[0][0] + anisobox[1][0] * anisobox[1][0] + anisobox[2][0] * anisobox[2][0]));
    Config::set().periodics.pb_box.y() = std::sqrt((anisobox[0][1] * anisobox[0][1] + anisobox[1][1] * anisobox[1][1] + anisobox[2][1] * anisobox[2][1]));
    Config::set().periodics.pb_box.z() = std::sqrt((anisobox[0][2] * anisobox[0][2] + anisobox[1][2] * anisobox[1][2] + anisobox[2][2] * anisobox[2][2]));
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
    if (step > Config::get().md.heat_steps[Config::get().md.heat_steps.size() - 1].offset)
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

// Nose-Hover thermostat for inner atoms. Variable names and implementation are identical to the book of
// Frenkel and Smit, Understanding Molecular Simulation, Appendix E
double md::simulation::nose_hoover_thermostat_some_atoms(std::vector<int> active_atoms)
{
  double tempscale(0.0);
  int freedom_some = 3U * active_atoms.size();
  if (Config::get().periodics.periodic == true)
    freedom_some -= 3;
  else
    freedom_some -= 6;
  double const delt(Config::get().md.timeStep),
    d2(delt / 2.0), d4(d2 / 2.0), d8(d4 / 2.0),
    TR(T*md::R), fTR(freedom_some*TR);
  nht.G2 = (nht.Q1*nht.v1*nht.v1 - TR) / nht.Q2;
  nht.v2 += nht.G2*d4;
  nht.v1 *= exp(-nht.v2*d8);
  nht.G1 = (2.0*E_kin - fTR) / nht.Q1;
  nht.v1 += nht.G1*d4;
  nht.v1 *= exp(-nht.v2*d8);
  nht.x1 += nht.v1*d2;
  nht.x2 += nht.v2*d2;
  tempscale = exp(-nht.v1*d2);
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
  return tempscale;
}

double md::simulation::tempcontrol(bool thermostat, bool half)
{
  std::size_t const N = this->coordobj.size();  // total number of atoms
  double tempfactor(2.0 / (freedom*md::R));     // factor for calculation of temperature from kinetic energy  
  double temp1, temp2 = 0, factor;     // current temperature before and after the temperature scaling, scaling factor

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
    if (Config::get().md.set_active_center == 1)  // if biased potential
    {
      size_t dof = 3u * inner_atoms.size();
      double T_factor = (2.0 / (dof*md::R));
      updateEkin_some_atoms(inner_atoms);           // calculate kinetic energy of inner atoms
      factor = nose_hoover_thermostat_some_atoms(inner_atoms);     // calculate temperature scaling factor
      for (auto i : movable_atoms) V[i] *= factor;   // new velocities (for all atoms that have a velocity)
      temp2 = E_kin * T_factor;     // new temperature (only inner atoms)
      if (half == false)
      {
        updateEkin();  // new kinetic energy (whole molecule)
      }
    }
    else if (Config::get().coords.fixed.size() != 0)  // if fixed atoms
    {
      size_t dof = 3u * movable_atoms.size();
      double T_factor = (2.0 / (dof*md::R));
      updateEkin_some_atoms(movable_atoms);           // calculate kinetic energy of movable atoms
      factor = nose_hoover_thermostat_some_atoms(movable_atoms);     // calculate temperature scaling factor
      for (auto i : movable_atoms) V[i] *= factor;   // new velocities (for all atoms that have a velocity)
      temp2 = E_kin * T_factor;     // new temperature (only inner atoms)
      if (half == false)
      {
        updateEkin();  // new kinetic energy (whole molecule)
      }
    }
    else  // "normal" nose-hoover thermostat
    {
      updateEkin();
      nose_hoover_thermostat();
      temp2 = E_kin * tempfactor;
    }
  }
  else // no thermostat -> direct scaling
  {
    if (Config::get().md.set_active_center == 1)
    {     // calculate temperature only for atoms inside inner cutoff
      updateEkin_some_atoms(inner_atoms); // kinetic energy of inner atoms
      size_t dof = 3u * inner_atoms.size();
      double T_factor = (2.0 / (dof*md::R));
      temp1 = E_kin*T_factor;           // temperature of inner atoms
      factor = std::sqrt(T / temp1);    // temperature scaling factor
      for (auto i : movable_atoms) V[i] *= factor;   // new velocities (for all atoms that have a velocity)
      if (half == false)
      {
        updateEkin_some_atoms(inner_atoms);
        temp2 = E_kin * T_factor;                   // new temperature of inner atoms
        updateEkin();            // kinetic energy
      }
    }
    else if (Config::get().coords.fixed.size() != 0)
    {     // calculate temperature only for atoms inside inner cutoff
      updateEkin_some_atoms(movable_atoms); // kinetic energy of inner atoms
      size_t dof = 3u * movable_atoms.size();
      double T_factor = (2.0 / (dof*md::R));
      temp1 = E_kin * T_factor;           // temperature of inner atoms
      factor = std::sqrt(T / temp1);    // temperature scaling factor
      for (auto i : movable_atoms) V[i] *= factor;   // new velocities (for all atoms that have a velocity)
      if (half == false)
      {
        updateEkin_some_atoms(movable_atoms);
        temp2 = E_kin * T_factor;                   // new temperature of inner atoms
        updateEkin();            // kinetic energy
      }
    }
    else
    {
      updateEkin();
      temp1 = E_kin * tempfactor;      // temperature before
      factor = std::sqrt(T / temp1);
      for (size_t i(0U); i < N; ++i) V[i] *= factor;  // new velocities
      if (half == false)
      {
        updateEkin();
        temp2 = E_kin * tempfactor;     // temperatures after
      }
    }


    if (Config::get().general.verbosity > 3 && half)
    {
      std::cout << "half step: desired temp: " << T << " current temp: " << temp1 << " factor: " << factor << "\n";
    }
    else if (Config::get().general.verbosity > 3)
    {
      std::cout << "full step: desired temp: " << T << " current temp: " << temp1 << " factor: " << factor << "\n";
    }
  }
  return temp2;
}

std::vector<double> md::simulation::init_active_center(int counter)
{
  std::size_t const N = this->coordobj.size();                // total number of atoms
  config::molecular_dynamics const & CONFIG(Config::get().md);

  auto split = std::max(std::min(std::size_t(CONFIG.num_steps / 100u), size_t(10000u)), std::size_t{ 100u });

  std::vector<double> dists;
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
    dists.push_back(distance);
    if (Config::get().general.verbosity > 3)
    {
      std::cout << "Atom " << i + 1 << ": Distance to active center: " << distance << "\n";
    }
  }
  return dists;
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
    std::cout << "Start again...\n";
  }
  coordobj.set_xyz(P_start);   // set all positions to start

  std::default_random_engine generator(static_cast<unsigned> (time(0)));  // generates random numbers
  auto dist01 = std::normal_distribution<double>{ 0,1 }; // normal distribution with mean=0 and standard deviation=1
  std::size_t const N = coordobj.size();
  for (auto i(0U); i < N; ++i)     // set random velocities around temperature T
  {
    V[i].x() = dist01(generator) * std::sqrt(kB*T / M[i]);
    V[i].y() = dist01(generator) * std::sqrt(kB*T / M[i]);
    V[i].z() = dist01(generator) * std::sqrt(kB*T / M[i]);
  }
  if (Config::get().md.set_active_center == 0)
  {
    tune_momentum();     // remove translation and rotation
  }
}


// Velcoity verlet integrator
void md::simulation::integrator(bool fep, std::size_t k_init, bool beeman)
{
  scon::chrono::high_resolution_timer integration_timer;

  config::molecular_dynamics const & CONFIG(Config::get().md);

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
  auto split = std::max(std::min(std::size_t(CONFIG.num_steps / 100u), size_t(10000u)), std::size_t{ 100u });
  for (std::size_t k(k_init); k < CONFIG.num_steps; ++k)
  {
    if (k == 0 && beeman == true)    // set F(t-dt) for first step to F(t)
    {
      for (size_t i = 0u; i < N; ++i)
      {
        F_old.push_back(coordobj.g_xyz(i));
      }
    }

    if (Config::get().general.verbosity > 3u)
    {
      std::cout << k << " of " << CONFIG.num_steps << " steps completed\n";
    }
    else if (Config::get().general.verbosity > 1u && k % split == 0 && k > 1)
    {
      std::cout << k << " of " << CONFIG.num_steps << " steps completed\n";
    }

    bool const HEATED(heat(k, fep));
    if (Config::get().md.temp_control == true)
    {
      // apply half step temperature corrections
      if (CONFIG.hooverHeatBath || HEATED)
      {
        temp = tempcontrol(CONFIG.hooverHeatBath, true);
      }
    }
    
    // save old coordinates
    P_old = coordobj.xyz();
    // Calculate new positions and half step velocities
    for (auto i : movable_atoms)
    {
      coords::Cartesian_Point acceleration;
      if (beeman == false)  //velocity-verlet
      {
        acceleration = coordobj.g_xyz(i)*md::negconvert / M[i];
        V[i] += acceleration*dt_2;
      }
      else  //beeman
      {
        acceleration = coordobj.g_xyz(i)*md::negconvert / M[i];
        coords::Cartesian_Point const acceleration_old(F_old[i] * md::negconvert / M[i]);
        V[i] += acceleration*(2.0 / 3.0)*dt - acceleration_old*(1.0 / 6.0)*dt;
      }

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
          std::cout << b[0]+1 << " and " << b[1]+1 << ", distance: " << b[2] << "\n";
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
      movable_atoms.clear();            // determine again which atoms are moved
      inner_atoms.clear();
      for (auto i(0U); i < N; ++i)
      {
        if (distances[i] < inner_cutoff)
        {
          inner_atoms.push_back(i);
        }
        if (distances[i] <= outer_cutoff)
        {
          movable_atoms.push_back(i);
        }
        else
        {
          V[i] = coords::Cartesian_Point(0, 0, 0);
        }
      }
    }
    // Apply first part of RATTLE constraints if requested
    if (CONFIG.rattle.use) rattle_pre();

    if (beeman == true)
    {
      for (size_t i = 0u; i < N; ++i)   // save F(t) as F_old
      {
        F_old[i] = coordobj.g_xyz(i);
      }
    }

    // calculate new energy & gradients
    coordobj.g();

    // Apply umbrella potential if umbrella sampling is used
    if (CONFIG.umbrella == true)
    {
      // apply biases and fill udatacontainer with values for restrained coordinates
      coordobj.ubias(udatacontainer);
    }
    // refine nonbondeds if refinement is required due to configuration
    if (CONFIG.refine_offset != 0 && (k + 1U) % CONFIG.refine_offset == 0)
    {
      if (Config::get().general.verbosity > 3U)
        std::cout << "Refining structure/nonbondeds.\n";
      coordobj.energy_update(true);
    }
    // If spherical boundaries are used apply boundary potential
    boundary_adjustments();
    // add new acceleration and calculate full step velocities
    for (auto i : movable_atoms)
    {
      if (beeman == false)  // velocity verlet
      {
        coords::Cartesian_Point const acceleration(coordobj.g_xyz(i)*md::negconvert / M[i]);
        V[i] += acceleration*dt_2;
      }
      else  // beeman
      {
        coords::Cartesian_Point const acceleration_new(coordobj.g_xyz(i)*md::negconvert / M[i]);
        coords::Cartesian_Point const acceleration(F_old[i] * md::negconvert / M[i]);
        V[i] += acceleration_new*(1.0 / 3.0)*dt + acceleration*(1.0 / 6.0)*dt;
      }

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
    if (Config::get().md.temp_control == true)
    {
      if (CONFIG.hooverHeatBath || HEATED)
      {
        temp = tempcontrol(CONFIG.hooverHeatBath, false);
      }
      else  // calculate E_kin and T if no temperature control is active (just nothing switched on)
      {
        double tempfactor(2.0 / (freedom*md::R));
        updateEkin();
        temp = E_kin * tempfactor;
      }
    }
    else  // calculate E_kin and T if no temperature control is active (switched off by MDtemp_control)
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

    // calculate distances that should be analyzed
    if (ana_pairs.size() > 0)
    {
      for (auto &p : ana_pairs)
      {
        p.dists.push_back(dist(coordobj.xyz(p.a), coordobj.xyz(p.b)));
      }
    }

    // calculate average temperature for every zone
    if (Config::get().md.analyze_zones == true)
    {
      for (auto &z : zones)
      {
        updateEkin_some_atoms(z.atoms);           
        int dof = 3u * z.atoms.size();
        z.temperatures.push_back(E_kin * (2.0 / (dof*md::R)));
      }
    }
  }  // end loop for every MD step

#ifdef USE_PYTHON
  // plot distances from MD analyzing
  if (Config::get().md.ana_pairs.size() > 0) plot_distances(ana_pairs);

  // plot average temperatures of every zone
  if (Config::get().md.analyze_zones == true) plot_zones();
#else
  if (Config::get().md.ana_pairs.size() > 0)
  {
    std::cout << "Plotting is not possible without python!\n";
    write_dists_into_file(ana_pairs);
  }
  if (Config::get().md.analyze_zones == true)
  {
    std::cout << "Plotting is not possible without python!\n";
    write_zones_into_file();
  }
#endif

  // calculate average pressure over whole simulation time
  p_average /= CONFIG.num_steps;
  if (Config::get().general.verbosity > 2U)
  {
    std::cout << "Average pressure: " << p_average << "\n";
    if (beeman == false)
    {
      std::cout << "Velocity-Verlet integration took " << integration_timer << '\n';
    }
    else
    {
      std::cout << "Beeman integration took " << integration_timer << '\n';
    }
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
        coords::Cartesian_Point rattle(delta / (2.0 * (inv_ma + inv_mb) * dot(d, d_old)));
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
        coords::Cartesian_Point rattle(d*lagrange);
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
