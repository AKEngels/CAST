#include "md.h"
#ifndef CAST_PSEUDO_RNG_DEBUG
#define CAST_PSEUDO_RNG_DEBUG
#endif

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
    md::nose_hoover_2chained& nht = this->thermostat.nht_2chained;
    nht = nose_hoover_2chained();
    nht.setQ1(Config::set().md.nosehoover_Q);
    nht.setQ2(Config::set().md.nosehoover_Q);
    this->thermostat.nht_v2 = md::nose_hoover_arbitrary_length(std::vector<double>(Config::get().md.nosehoover_chainlength, Config::get().md.nosehoover_Q));

    desired_temp = Config::get().md.T_init;
    init();
    // remove rotation and translation of the molecule if desired (only if no biased potential is applied)
    if (Config::get().md.set_active_center == 0 && Config::get().md.veloScale)
    {
      removeTranslationalAndRotationalMomentumOfWholeSystem();
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
  double const psec(Config::get().md.timeStep * static_cast<double>(Config::get().md.num_steps));
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
  const bool thermostat_active = Config::get().md.thermostat_algorithm != 0u;
  std::cout << "Thermostat is " << (thermostat_active ? "active." : "inactive.") << '\n';
  if (thermostat_active)
  {
    if (Config::get().md.thermostat_algorithm == config::molecular_dynamics::thermostat_algorithms::TWO_NOSE_HOOVER_CHAINS)
    {
      std::cout << "Thermostat is " << "\"2 Nose-Hoover chains\". Implementation according to Frenkel & Smith, \"Understanding Molecular Simulation\", 2002, Appendix E." << '\n';
      std::cout << "Thermostat Mass Q is Q1=Q2=" << std::to_string(Config::get().md.nosehoover_Q) << ".\n";
    }
    else if (Config::get().md.thermostat_algorithm == config::molecular_dynamics::thermostat_algorithms::ARBITRARY_CHAIN_LENGTH_NOSE_HOOVER)
    {
      std::cout << "Thermostat is " << "\"Arbitrary Nose-Hoover chains\". Implementation according to:\n";
      std::cout << "Glenn J. Martyna, Mark E. Tuckerman, Douglas J. Tobias & Michael L. Klein (1996) Explicit reversible integrators for extended systems dynamics, Molecular Physics, 87:5,1117 - 1157, DOI : 10.1080 / 00268979600100761\n";
      std::cout << "Chainlength is: " << std::to_string(this->thermostat.nht_v2.chainlength) << '\n';
      std::cout << "Thermostat Masses Q are: { ";
      for (auto i : this->thermostat.nht_v2.masses_param_Q) 
        std::cout << std::to_string(i) << " , ";
      std::cout << " }.\n";
    }
    else if (Config::get().md.thermostat_algorithm == config::molecular_dynamics::thermostat_algorithms::BERENDSEN)
    {
      std::cout << "Thermostat is " << "\"BERENDSEN\". Implementation according to:\n";
      std::cout << "Adv. Polym. Sci. (2005) 173:105-149 DOI:10.1007 / b99427 - Thermostat Algorithms for Molecular Dynamics Simulations\n";
      std::cout << "Parameter t_B is: " << std::to_string(this->thermostat.berendsen_tB) << '\n';
    }
    else if (Config::get().md.thermostat_algorithm == config::molecular_dynamics::thermostat_algorithms::HOOVER_EVANS)
    {
      std::cout << "Thermostat is " << "\"HOOVER-EVANS\". Implementation according to:\n";
      std::cout << "Adv. Polym. Sci. (2005) 173:105-149 DOI:10.1007 / b99427 - Thermostat Algorithms for Molecular Dynamics Simulations\n";
    }
  }
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
      for (auto const& bond : Config::get().md.rattle.specified_rattle)
      {
        std::cout << "[" << bond.a + 1 << ":" << bond.b + 1 << "] ";
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


// Initialization of MD parameters and functions
void md::simulation::init(void)
{
  if (coordobj.validate_bonds() == false)  // test if there are broken bonds in the structure and save them
  {
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << "Warning! Broken bonds in your structure even before the simulation starts! Atom numbers: \n";
      for (auto b : coordobj.getBrokenBonds())
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

#ifdef CAST_PSEUDO_RNG_DEBUG
  // via https://stackoverflow.com/questions/15500621/c-c-algorithm-to-produce-same-pseudo-random-number-sequences-from-same-seed-on
  std::mt19937 generator(0); // Fixed seed of 0
#else
  std::default_random_engine generator(static_cast<unsigned> (time(0)));  // generates random numbers
#endif
  auto dist01 = std::normal_distribution<double>{ 0,1 }; // normal distribution with mean=0 and standard deviation=1

  for (std::size_t i = 0; i < N; ++i)
  {
    // Get Atom Mass
    M[i] = coordobj.atoms(i).mass();
    // Sum total Mass
    M_total += M[i];
    // Set initial velocities to zero if fixed or not movable
    if (coordobj.atoms(i).fixed() || abs(Config::get().md.T_init) == 0.0 || std::find(movable_atoms.begin(), movable_atoms.end(), i) == movable_atoms.end())
    {
      V[i] = coords::Cartesian_Point(0.);
    }
    // initialize random velocities otherwise (Maxwell Boltzmann distribution, see http://research.chem.psu.edu/shsgroup/chem647/newNotes/node6.html)
    else
    {
      V[i].x() = dist01(generator) * std::sqrt(gasconstant_R_1 * desired_temp / M[i]); // yields velocities in Angst/ps
      V[i].y() = dist01(generator) * std::sqrt(gasconstant_R_1 * desired_temp / M[i]); // yields velocities in Angst/ps
      V[i].z() = dist01(generator) * std::sqrt(gasconstant_R_1 * desired_temp / M[i]); // yields velocities in Angst/ps

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
    if ((Config::get().md.set_active_center == 1 || Config::get().coords.fixed.size() != 0)
      && Config::get().md.thermostat_algorithm != config::molecular_dynamics::thermostat_algorithms::VELOCITY_RESCALING)
    {
      std::cout << "ERROR! Currently NVT ensembles using thermostats are only available for either: \n";
      std::cout << " - RATTLE constraints\n - Fixed atoms\nor\n - Active Center Dynamics\n";
      std::cout << "No combination of these may be specified.";
      throw(std::runtime_error("Aborting..."));
    }
  }
  else if (Config::get().md.set_active_center == 1 && Config::get().coords.fixed.size() != 0)
  {
    std::cout << "ERROR! Currently NVT ensembles using thermostats are only available for either: \n";
    std::cout << " - RATTLE constraints\n - Fixed atoms\nor\n - Active Center Dynamics\n";
    std::cout << "No combination of these may be specified.";
    throw(std::runtime_error("Aborting..."));
  }
  if (Config::get().md.thermostat_algorithm == config::molecular_dynamics::thermostat_algorithms::HOOVER_EVANS)
  {
    if (Config::get().md.set_active_center == 1 || Config::get().coords.fixed.size() != 0 || Config::get().md.rattle.use || Config::get().md.heat_steps.size() != 1u
      || Config::get().md.T_final != Config::get().md.T_init)
    {
      std::cout << "ERROR!\nDon't use the Hoover-Evans thermostat with constraints or varying temperature.\n";
      std::cout << "Options causing this error might be fixed atoms, active center dynamics, RATTLE constraints or heat ramping.\nPlease change your inputfile. Aborting." << std::endl;
      throw std::runtime_error("Aborting...");
    }
  }

  // periodics and isothermal cases
  if (Config::get().md.thermostat_algorithm != config::molecular_dynamics::thermostat_algorithms::VELOCITY_RESCALING)
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

  if (Config::get().md.analyze_zones == true) zones = md_analysis::find_zones(this);  // find atoms for every zone

  if (abs(Config::get().md.T_init) == 0.0 && Config::get().md.temp_control 
  && !Config::get().md.thermostat_algorithm == config::molecular_dynamics::thermostat_algorithms::VELOCITY_RESCALING && false)
  {
    throw std::runtime_error("Velocity rescaling only works with non-zero starting temperature. Use a low value, like 5 K, instead. Exiting!");;
  }
  this->thermostat = md::thermostat_data(md::nose_hoover_arbitrary_length(std::vector<double>(Config::get().md.nosehoover_chainlength, 
    Config::get().md.nosehoover_Q)), md::nose_hoover_2chained(Config::get().md.nosehoover_Q));
}


// eliminate translation and rotation of the system at the beginning of a MD simulation
// can also be performed at the end of every MD step
void md::simulation::removeTranslationalAndRotationalMomentumOfWholeSystem(void)
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
    momentum_angular += cross(coordobj.xyz(i), V[i]) * M[i];
  }
  // scale by total mass
  momentum_linear /= M_total;
  mass_vector /= M_total;
  // MVxLM * M 
  momentum_angular -= cross(mass_vector, momentum_linear) * M_total;
  // momentum of inertia from each component
  double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
  for (std::size_t i = 0; i < N; ++i)
  {
    coords::Cartesian_Point r(coordobj.xyz(i) - mass_vector);
    xx += r.x() * r.x() * M[i];
    xy += r.x() * r.y() * M[i];
    xz += r.x() * r.z() * M[i];
    yy += r.y() * r.y() * M[i];
    yz += r.y() * r.z() * M[i];
    zz += r.z() * r.z() * M[i];
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
  if (Config::get().general.verbosity >= 5u) std::cout << "Eliminated Translation&Rotation of whole system.\n";
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
      sp_g = r * (Config::get().md.spherical.e1 * Config::get().md.spherical.f1 * std::pow(L - Config::get().md.spherical.r_inner, Config::get().md.spherical.e1 - 1) / L);
      coordobj.add_sp_gradients(i, sp_g);
    }
  }
}

// Calculates current kinetic energy from velocities
coords::float_type md::simulation::getEkin(std::vector<std::size_t> atom_list) const
{
  // initialize tensor to zero
  using TenVal = coords::Tensor::value_type;
  const TenVal z = { 0.0,0.0,0.0 };
  coords::Tensor cE_kin_tensor = coords::Tensor(this->E_kin_tensor);
  cE_kin_tensor.fill(z);
  // calculate contribution to kinetic energy for each atom
  for (auto i : atom_list)
  {
    auto const fact = 0.5 * M[i] / convert;
    cE_kin_tensor[0][0] += fact * V[i].x() * V[i].x();
    cE_kin_tensor[1][0] += fact * V[i].x() * V[i].y();
    cE_kin_tensor[2][0] += fact * V[i].x() * V[i].z();
    cE_kin_tensor[0][1] += fact * V[i].y() * V[i].x();
    cE_kin_tensor[1][1] += fact * V[i].y() * V[i].y();
    cE_kin_tensor[2][1] += fact * V[i].y() * V[i].z();
    cE_kin_tensor[0][2] += fact * V[i].z() * V[i].x();
    cE_kin_tensor[1][2] += fact * V[i].z() * V[i].y();
    cE_kin_tensor[2][2] += fact * V[i].z() * V[i].z();
  }
  // calculate total kinetic energy by the trace of the tensor
  coords::float_type cE_kin = cE_kin_tensor[0][0] + cE_kin_tensor[1][1] + cE_kin_tensor[2][2];
  if (Config::get().general.verbosity > 4u)
  {
    std::cout << "New kinetic Energy is " << cE_kin << "kcal/mol with E_kin(x), (y), (z) = " << cE_kin_tensor[0][0] << ", "
      << cE_kin_tensor[1][1] << ", " << cE_kin_tensor[2][2] << "." << '\n';
  }
  return cE_kin;
}

// Calculates current kinetic energy from velocities
void md::simulation::updateEkin(std::vector<std::size_t> atom_list)
{
  this->E_kin = getEkin(atom_list);
}

// apply pressure corrections if constant pressure simulation is performed
void md::simulation::berendsen(double const time)
{
  // temp variables
  const std::size_t N = coordobj.size();
  double fac, scale, volume;
  double ptensor[3][3], aniso[3][3], anisobox[3][3], dimtemp[3][3];
  // get volume and pressure scaling factor
  volume = Config::get().periodics.pb_box.x() * Config::get().periodics.pb_box.y() * Config::get().periodics.pb_box.z();
  fac = presc / volume;
  // pressure for ISOTROPIC boxes
  if (Config::get().energy.isotropic == true) {
    ptensor[0][0] = fac * 2.0 * (E_kin_tensor[0][0] - coordobj.virial()[0][0]);
    ptensor[1][1] = fac * 2.0 * (E_kin_tensor[1][1] - coordobj.virial()[1][1]);
    ptensor[2][2] = fac * 2.0 * (E_kin_tensor[2][2] - coordobj.virial()[2][2]);
    press = (ptensor[0][0] + ptensor[1][1] + ptensor[2][2]) / 3.0;  // pressure in bar or atm ???
    // Berendsen scaling for isotpropic boxes
    scale = std::pow((1.0 + (time * Config::get().md.pcompress / Config::get().md.pdelay) * (press - Config::get().md.ptarget)), 0.3333333333333);
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
        ptensor[j][i] = fac * (2.0 * E_kin_tensor[j][i] - coordobj.virial()[j][i]);
      }
    }
    // get isotropic pressure value
    press = (ptensor[0][0] + ptensor[1][1] + ptensor[2][2]) / 3.0;
    // Anisotropic scaling factors
    scale = time * Config::get().md.pcompress / (3.0 * Config::get().md.pdelay);
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

std::vector<double> md::simulation::init_active_center(int counter)
{
  std::size_t const N = this->coordobj.size();                // total number of atoms
  config::molecular_dynamics const& CONFIG(Config::get().md);

  auto split = std::max(std::min(std::size_t(CONFIG.num_steps / 100u), size_t(10000u)), std::size_t{ 100u });

  std::vector<double> dists;
  std::vector<coords::Cartesian_Point> coords_act_center;
  for (auto& atom_number : Config::get().md.active_center)
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
  for (auto& atom_coords : coords_act_center)
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
    double distance = sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
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

#ifdef CAST_PSEUDO_RNG_DEBUG
  // via https://stackoverflow.com/questions/15500621/c-c-algorithm-to-produce-same-pseudo-random-number-sequences-from-same-seed-on
  std::mt19937 generator(0); // Fixed seed of 0
#else
  std::default_random_engine generator(static_cast<unsigned> (time(0)));  // generates random numbers
#endif
  auto dist01 = std::normal_distribution<double>{ 0,1 }; // normal distribution with mean=0 and standard deviation=1
  std::size_t const N = coordobj.size();
  for (auto i(0U); i < N; ++i)     // set random velocities around temperature T
  {
    V[i].x() = dist01(generator) * std::sqrt(gasconstant_R_1 * desired_temp / M[i]); // yields velocities in Angst/ps
    V[i].y() = dist01(generator) * std::sqrt(gasconstant_R_1 * desired_temp / M[i]); // yields velocities in Angst/ps
    V[i].z() = dist01(generator) * std::sqrt(gasconstant_R_1 * desired_temp / M[i]); // yields velocities in Angst/ps
  }
  if (Config::get().md.set_active_center == 0 && Config::get().md.veloScale)
  {
    removeTranslationalAndRotationalMomentumOfWholeSystem();     // remove translation and rotation
  }
}


// Velcoity verlet integrator
void md::simulation::integrator(bool fep, std::size_t k_init, bool beeman)
{
  scon::chrono::high_resolution_timer integration_timer;

  config::molecular_dynamics const& CONFIG(Config::get().md);

  std::size_t const N = this->coordobj.size();
  // set average pressure to zero
  double p_average(0.0);
  // constant values
  double const
    //velofactor(-0.5*dt*md::convert),
    dt_2(0.5 * dt);

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



    // save old coordinates
    P_old = coordobj.xyz();
    // Calculate new positions and half step velocities
    for (auto i : movable_atoms)
    {
      coords::Cartesian_Point acceleration;
      if (beeman == false)  //velocity-verlet
      {
        acceleration = coordobj.g_xyz(i) * md::negconvert / M[i];
        V[i] += acceleration * dt_2;
      }
      else  //beeman
      {
        acceleration = coordobj.g_xyz(i) * md::negconvert / M[i];
        coords::Cartesian_Point const acceleration_old(F_old[i] * md::negconvert / M[i]);
        V[i] += acceleration * (2.0 / 3.0) * dt - acceleration_old * (1.0 / 6.0) * dt;
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
    }
    //Fetching target temperature
    bool const is_not_microcanonical = determine_current_desired_temperature(k, fep);
    if (CONFIG.temp_control == true && is_not_microcanonical)
    {
      // apply half step temperature corrections
      if (CONFIG.thermostat_algorithm == config::molecular_dynamics::thermostat_algorithms::ARBITRARY_CHAIN_LENGTH_NOSE_HOOVER
        || CONFIG.thermostat_algorithm == config::molecular_dynamics::thermostat_algorithms::TWO_NOSE_HOOVER_CHAINS)
      {
        this->instantaneous_temp = tempcontrol(CONFIG.thermostat_algorithm, true);
      }
    }
    for (auto i : movable_atoms)
    {  // update coordinates
      coordobj.move_atom_by(i, V[i] * dt);
    }

    if (coordobj.validate_bonds() == false)  // look if all bonds are okay and save those which aren't 
    {
      if (Config::get().general.verbosity > 1U)
      {                                          // give warning if there are broken bonds
        std::cout << "Warning! Broken bonds between atoms...\n";
        for (auto b : coordobj.getBrokenBonds())
        {
          std::cout << b[0] + 1 << " and " << b[1] + 1 << ", distance: " << b[2] << "\n";
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
        coords::Cartesian_Point const acceleration(coordobj.g_xyz(i) * md::negconvert / M[i]);
        V[i] += acceleration * dt_2;
      }
      else  // beeman
      {
        coords::Cartesian_Point const acceleration_new(coordobj.g_xyz(i) * md::negconvert / M[i]);
        coords::Cartesian_Point const acceleration(F_old[i] * md::negconvert / M[i]);
        V[i] += acceleration_new * (1.0 / 3.0) * dt + acceleration * (1.0 / 6.0) * dt;
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
    if (CONFIG.temp_control == true && is_not_microcanonical)
    {
      this->instantaneous_temp = tempcontrol(CONFIG.thermostat_algorithm, false);
    }
    else  // calculate E_kin and T if no temperature control is active (switched off by MDtemp_control)
    {
      const double tempfactor(2.0 / (freedom * md::R));
      updateEkin(range(N));            // kinetic energy
      this->instantaneous_temp = E_kin * tempfactor;
    }

    // Apply pressure adjustments
    if (CONFIG.pressure)
    {
      berendsen(dt);
    }
    // save temperature for FEP
    if (Config::get().md.fep)
    {
      coordobj.getFep().fepdata.back().T = this->instantaneous_temp;
    }
    // if requested remove translation and rotation of the system
    if (Config::get().md.veloScale) removeTranslationalAndRotationalMomentumOfWholeSystem();

    // Logging / Traces

    if (CONFIG.track)
    {
      std::vector<coords::float_type> iae;
      if (coordobj.interactions().size() > 1)
      {
        iae.reserve(coordobj.interactions().size());
        for (auto const& ia : coordobj.interactions()) iae.push_back(ia.energy);
      }
      logging(k, this->instantaneous_temp, press, E_kin, coordobj.pes().energy, iae, coordobj.xyz());
    }

    // Serialize to binary file if required.
    if (k > 0 && Config::get().md.restart_offset > 0 && k % Config::get().md.restart_offset == 0)
    {
      write_restartfile(k);
    }
    // add up pressure value
    p_average += press;

    // get info for analysis
    md_analysis::add_analysis_info(this);
  }

  // log analysis info
  md_analysis::write_and_plot_analysis_info(this);

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
