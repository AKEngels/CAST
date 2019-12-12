/**
CAST 3
Purpose: Tests for forcefield (at the moment only OPLSAA)

@author Susanne Sauer
@version 1.0
*/

// TODO: add tests for an example with improper dihedral!!! (compare to tinker)

#ifdef GOOGLE_MOCK

#include "../../energy_int_aco.h"
#include "../../tinker_parameters.h"
#include "../../coords_io.h"
#include <gtest/gtest.h>

TEST(forcefield, test_total_energy)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  double energy = coords.e();

  ASSERT_NEAR(energy, 6.5344, 0.0001);     // from tinker
}

TEST(forcefield, test_bonds_number_and_energy)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  energy::interfaces::aco::aco_ff y(&coords);

  // number of bonds
  auto number_of_bonds = y.refined.bonds().size();
  ASSERT_EQ(number_of_bonds, 14);          // from tinker
   
  // energy
  y.e();
  double energy = y.part_energy[energy::interfaces::aco::aco_ff::types::BOND];
  ASSERT_NEAR(energy, 0.6187, 0.0001);     // from tinker
}

TEST(forcefield, test_angles_number_and_energy)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  energy::interfaces::aco::aco_ff y(&coords);

  // number of angles
  auto number_of_angles = y.refined.angles().size();
  ASSERT_EQ(number_of_angles, 25);          // from tinker

  // energy
  y.e();
  double energy = y.part_energy[energy::interfaces::aco::aco_ff::types::ANGLE];
  ASSERT_NEAR(energy, 0.5885, 0.0001);     // from tinker
}

TEST(forcefield, test_torsions_number_and_energy)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  energy::interfaces::aco::aco_ff y(&coords);

  // number of torsion
  auto number_of_torsions = y.refined.torsions().size();
  ASSERT_EQ(number_of_torsions, 30);          // from tinker

  // energy
  y.e();
  double energy = y.part_energy[energy::interfaces::aco::aco_ff::types::TORSION];
  ASSERT_NEAR(energy, 1.1390, 0.0001);     // from tinker
}

TEST(forcefield, test_vdw_energy)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  energy::interfaces::aco::aco_ff y(&coords);

  // energy
  y.e();
  double energy = y.part_energy[energy::interfaces::aco::aco_ff::types::VDW];
  ASSERT_NEAR(energy, 0.4742, 0.0001);     // from tinker
}

TEST(forcefield, test_coulomb_energy)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  energy::interfaces::aco::aco_ff y(&coords);

  // energy
  y.e();
  double energy = y.part_energy[energy::interfaces::aco::aco_ff::types::CHARGE];
  ASSERT_NEAR(energy, 3.7139, 0.0001);     // from tinker
}

#endif
