#include "md.h"
#include "md_umbrella.h"
#pragma once

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
    md::nose_hoover_2chained& nht = this->thermostat.nht_2chained;
    nht = nose_hoover_2chained();
    desired_temp = Config::get().md.T_init;
    init();
    removeTranslationalAndRotationalMomentumOfWholeSystem(); // eliminate translation and rotation
  }
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
      double const L(scon::geometric_length(xyz(i) - xyz(bound)));
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
