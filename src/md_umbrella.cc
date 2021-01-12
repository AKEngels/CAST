#include "md.h"
#include "md_umbrella.h"

// Perform Umbrella Sampling run if requested
void md::simulation::umbrella_run(bool const restart) {
  std::size_t steps;
  steps = Config::get().md.num_steps;

  //General md initialization and config
  if (Config::get().general.verbosity > 0U) print_init_info();
  coordobj.g();
  restarted = restart;
  if (restarted)
  {
    md::nose_hoover_2chained& nht = this->thermostat.nht_2chained;
    nht = nose_hoover_2chained();
    desired_temp = Config::get().md.T_init;
    init();
    removeTranslationalAndRotationalMomentumOfWholeSystem(); // eliminate translation and rotation
  }
  // if PMF-IC: create spline function
  if (Config::get().coords.umbrella.pmf_ic.use) create_uspline();
  // Set kinetic Energy
  updateEkin(range(coordobj.size()));            // kinetic energy
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
  auto&& number_of_restraints = Config::get().coords.bias.udist.size() + Config::get().coords.bias.utors.size()
    + Config::get().coords.bias.ucombs.size() + Config::get().coords.bias.uangles.size();
  for (auto s{ 0u }; s < udatacontainer.size() / number_of_restraints; ++s)  // for every step 
  {
    ofs << s << "   ";                                              // stepnumber
    for (auto b{ 0u }; b < number_of_restraints; ++b) {
      ofs << udatacontainer[b + number_of_restraints * s] << "  ";  // value(s)
    }
    ofs << "\n";
  }
  ofs.close();
}



// The Coords Object's functions
bool md::CoordinatesUBIAS::validate_bonds()
{
  bool status = true;
  broken_bonds.clear();
  auto const N = size();
  for (auto i = 0u; i < N; ++i)  // for every atom i
  {
    for (auto const& bound : atoms(i).bonds())  // for every atom b that is bound to i
    {
      auto const pos_1 = xyz(i);
      auto const pos_2 = xyz(bound);
      double const L(scon::geometric_length(pos_1 - pos_2));
      if (L < 0.3 || L > 5.0)  // test if bondlength between i and b is reasonable
      {
        status = false;
        if (i < bound)   // save all bonds with strange bondlengths in broken_bonds
        {
          std::vector<float> bond;
          bond.push_back(i);
          bond.push_back(bound);
          bond.push_back(L);
          broken_bonds.push_back(bond);
        }
      }
    }
  }
  return status;
}

void md::simulation::create_uspline()
{
  if (!file_exists(Config::get().coords.umbrella.pmf_ic.prepfile_name))
    throw std::runtime_error("File for getting spline function not found.");

  std::ifstream input(Config::get().coords.umbrella.pmf_ic.prepfile_name);
  std::string line;
  std::vector<std::string> linestr;
  std::getline(input, line);        // discard first line

  if (Config::get().coords.umbrella.pmf_ic.indices_xi.size() == 1)   // one-dimensional
  {
    std::vector<double> zs;
    std::vector<double> deltaEs;

    while (!input.eof())
    {
      std::getline(input, line);
      if (line == "") break;
      linestr = split(line, ',');
      zs.emplace_back(mapping::xi_to_z(std::stod(linestr[0]), Config::get().coords.umbrella.pmf_ic.xi0[0], Config::get().coords.umbrella.pmf_ic.L[0]));
      deltaEs.emplace_back(std::stod(linestr[4]));
    }
    Spline1D s;
    s.fill(zs, deltaEs);
    umbrella_spline = std::make_unique<Spline1D>(s);
  }

  else           // two-dimensional
  {
    std::vector<std::pair<double, double>> zs;
    std::vector<double> deltaEs;

    while (!input.eof())
    {
      std::getline(input, line);
      if (line == "") break;
      linestr = split(line, ',');
      double z1 = mapping::xi_to_z(std::stod(linestr[0]), Config::get().coords.umbrella.pmf_ic.xi0[0], Config::get().coords.umbrella.pmf_ic.L[0]);
      double z2 = mapping::xi_to_z(std::stod(linestr[2]), Config::get().coords.umbrella.pmf_ic.xi0[1], Config::get().coords.umbrella.pmf_ic.L[1]);
      zs.emplace_back(std::make_pair(z1, z2));
      deltaEs.emplace_back(std::stod(linestr[6]));
    }
    Spline2D s;
    s.fill(zs, deltaEs);
    umbrella_spline = std::make_unique<Spline2D>(s);
  }
}
