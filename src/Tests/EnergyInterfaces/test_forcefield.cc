/**
CAST 3
Purpose: Tests for forcefield (at the moment only OPLSAA)
         Notice: charge and vdW interactions are done in file "test_forcefield_non-bonded.cc"

@author Susanne Sauer
@version 1.0
*/

// TODO: add tests for an example with improper dihedral!!! (compare to tinker)

#ifdef GOOGLE_MOCK

#include "../../energy_int_aco.h"
#include "../../tinker_parameters.h"
#include "../../coords_io.h"
#include <gtest/gtest.h>

/**the functions defined in this class are protected in base class so they can't be used directly*/
class ClassToChangeExternalCharges : public energy::interface_base
{
public:
  /**function to set external charges*/
  static void set_external_charges(std::vector<energy::PointCharge> const& new_ext_charges) {
    interface_base::set_external_charges(new_ext_charges);
  }
  /**function that deletes all external charges*/
  static void clear_external_charges() {
    interface_base::clear_external_charges();
  }
};

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
  y.update();  // initialization of interface

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
  y.update();  // initialization of interface

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
  y.update();  // initialization of interface

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
  y.update();  // initialization of interface

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
  y.update();  // initialization of interface

  // energy
  y.e();
  double energy = y.part_energy[energy::interfaces::aco::aco_ff::types::CHARGE];
  ASSERT_NEAR(energy, 3.7139, 0.0001);     // from tinker
}

TEST(forcefield, test_total_gradients)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  double energy = coords.g();
  ASSERT_NEAR(energy, 6.5344, 0.0001);     // from tinker

  // this gradient is just taken from CAST, hoping it is correct
  coords::Representation_3D expected_grad = {
    coords::r3(5.98878, 6.38636, -2.48425),
    coords::r3(1.90015, -8.80451, 7.23128),
    coords::r3(-2.39712, 0.892749, -2.16194),
    coords::r3(0.236449, 0.439926, 2.3236),
    coords::r3(2.26011, -1.62977, -1.7127),
    coords::r3(-5.18402, -10.0051, -7.16721),
    coords::r3(-1.68023, 2.80079, 2.88759),
    coords::r3(0.0844947, -0.00660548, -3.70373),
    coords::r3(-14.0849, 8.07328, -0.916618),
    coords::r3(1.05412, 1.57882, 3.83898),
    coords::r3(3.11387, 3.40486, -1.48488),
    coords::r3(0.972518, -22.0797, 25.3103),
    coords::r3(1.17297, -0.632528, -2.32377),
    coords::r3(1.05517, -0.952779, 3.21229),
    coords::r3(5.50761, 20.5343, -22.8489)
  };

  ASSERT_TRUE(is_nearly_equal(expected_grad, coords.g_xyz(), 0.0001));
}

TEST(forcefield, test_bonded_gradients)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  energy::interfaces::aco::aco_ff y(&coords);
  y.update();  // initialization of interface
  y.g();

  // this gradient is just taken from CAST, hoping it is correct
  coords::Representation_3D expected_grad = {
    coords::r3(4.13377, 6.31019, -1.8893),
    coords::r3(0.103714, -7.73059, 3.94103),
    coords::r3(-3.14763, 0.108121, -1.02539),
    coords::r3(0.00285454, -0.216821, 2.70685),
    coords::r3(1.21699, -2.4538, -1.1475),
    coords::r3(-0.890542, -7.40685, -2.98018),
    coords::r3(-2.23064, 3.73573, 1.63896),
    coords::r3(-0.187317, 0.241709, -4.25739),
    coords::r3(-7.38174, 4.29063, 0.303439),
    coords::r3(0.245427, 1.01147, 4.85198),
    coords::r3(2.24705, 4.0521, -2.36037),
    coords::r3(-1.9497, -19.5603, 22.9003),
    coords::r3(-0.381172, -0.950617, -2.77704),
    coords::r3(-0.827766, -2.28676, 1.948),
    coords::r3(9.0467, 20.8558, -21.8534),
  };
  auto calculated_grad = y.part_grad[energy::interfaces::aco::aco_ff::types::BOND];

  ASSERT_TRUE(is_nearly_equal(expected_grad, calculated_grad, 0.0001));
}

TEST(forcefield, test_angle_gradients)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  energy::interfaces::aco::aco_ff y(&coords);
  y.update();  // initialization of interface
  y.g();

  // this gradient is just taken from CAST, hoping it is correct
  coords::Representation_3D expected_grad = {
    coords::r3(-1.32223, 0.649877, -0.0594146),
    coords::r3(1.43172, -0.548423, 2.82758),
    coords::r3(0.355399, 0.340698, -1.05504),
    coords::r3(-0.263786, 0.223375, 0.0181707),
    coords::r3(-1.25065, -0.355016, -0.567236),
    coords::r3(-2.9449, -2.3828, -3.55351),
    coords::r3(0.4832, -0.0675362, 0.811576),
    coords::r3(0.00331078, -0.691372, -0.0393976),
    coords::r3(-2.06179, 4.45587, -1.67498),
    coords::r3(1.50154, 0.514236, -0.183152),
    coords::r3(-0.714459, 0.823454, 0.733484),
    coords::r3(2.51735, -3.01884, 2.4552),
    coords::r3(2.13593, -0.29333, -0.192764),
    coords::r3(2.33816, 0.302754, 1.34896),
    coords::r3(-2.20879, 0.0470505, -0.869473),
  };
  auto calculated_grad = y.part_grad[energy::interfaces::aco::aco_ff::types::ANGLE];

  ASSERT_TRUE(is_nearly_equal(expected_grad, calculated_grad, 0.0001));
}

TEST(forcefield, test_torsion_gradients)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  energy::interfaces::aco::aco_ff y(&coords);
  y.update();  // initialization of interface
  y.g();

  // this gradient is just taken from CAST, hoping it is correct
  coords::Representation_3D expected_grad = {
    coords::r3(0.263949, -0.392796, -0.656521),
    coords::r3(0.147911, -0.0849742, 0.400714),
    coords::r3(-0.0166528, 0.0294108, 0.0542201),
    coords::r3(0.044002, -0.0266906, -0.00218435),
    coords::r3(-0.0877141, 0.0101436, -0.114717),
    coords::r3(0.00216677, -0.173661, -1.06443),
    coords::r3(-0.256713, -0.287581, 0.306103),
    coords::r3(-0.112701, 0.859533, 0.0537577),
    coords::r3(-0.18497, -0.228971, 0.405684),
    coords::r3(-0.531684, 0.47203, -0.0715075),
    coords::r3(0.710963, -0.15112, 0.4174),
    coords::r3(-0.107265, 0.579247, 0.0103564),
    coords::r3(-0.580086, -0.331866, 0.193224),
    coords::r3(0.586404, 0.128942, 0.400546),
    coords::r3(0.122388, -0.401647, -0.332645),
  };
  auto calculated_grad = y.part_grad[energy::interfaces::aco::aco_ff::types::TORSION];

  ASSERT_TRUE(is_nearly_equal(expected_grad, calculated_grad, 0.0001));
}

TEST(forcefield, test_total_energy_with_external_charges)
{
  // setting some random external charges, only scaled_charge is used during calculation
  ClassToChangeExternalCharges::set_external_charges({ { -4, 5, -2, 0.75, double() }, { 5, -4, 2, 0.5, double() } });

  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  double energy = coords.e();

  ASSERT_NEAR(energy, 7.3448092340770454, 0.00001);     // just taken from CAST hoping it is correct
  ClassToChangeExternalCharges::clear_external_charges();
}

TEST(forcefield, test_total_gradients_with_external_charges)
{
  // setting some random external charges, only scaled_charge is used during calculation
  ClassToChangeExternalCharges::set_external_charges( {{-4, 5, -2, 0.75, double()}, {5, -4, 2, 0.5, double()}});

  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  double energy = coords.g();
  ASSERT_NEAR(energy, 7.3448092340770454, 0.00001);     // just taken from CAST hoping it is correct

  // this gradient is just taken from CAST, hoping it is correct
  coords::Representation_3D expected_grad = {
    coords::r3(6.06483, 4.87113, -1.70728),
    coords::r3(2.62259, -10.299, 8.08931),
    coords::r3(-2.29696, 1.50429, -2.42054),
    coords::r3(0.253274, 0.777142, 2.04098),
    coords::r3(2.25547, -1.24591, -1.82545),
    coords::r3(-4.50893, -10.6314, -6.6368),
    coords::r3(-2.13526, 3.74566, 1.92882),
    coords::r3(-0.449449, 1.06077, -3.8962),
    coords::r3(-14.4192, 8.47596, -1.0907),
    coords::r3(0.818244, 1.76258, 3.51958),
    coords::r3(2.57517, 3.6693, -1.73185),
    coords::r3(2.21904, -22.1024, 25.5906),
    coords::r3(1.02722, -0.412484, -2.30879),
    coords::r3(1.03842, -0.822801, 3.12971),
    coords::r3(4.15491, 20.6454, -22.8468),
  };

  ASSERT_TRUE(is_nearly_equal(expected_grad, coords.g_xyz(), 0.0001));
  ClassToChangeExternalCharges::clear_external_charges();
}

TEST(forcefield, test_total_energy_with_external_charges_is_sum)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  double energy_without_extCharges = coords.e();     
  ASSERT_NEAR(energy_without_extCharges, 6.5344, 0.0001);     // from tinker, see above

  // setting some random external charges, only scaled_charge is used during calculation
  ClassToChangeExternalCharges::set_external_charges({ {-4, 5, -2, 0.75, double()}, {5, -4, 2, 0.5, double()} });

  energy::interfaces::aco::aco_ff y(&coords);
  y.update();  // initialization of interface
  double energy_with_extCharges = y.e();   
  ASSERT_NEAR(energy_with_extCharges, 7.3448092340770454, 0.00001);     // see above

  double energy_extCharges = y.part_energy[energy::interfaces::aco::aco_ff::types::EXTERNAL_CHARGES];

  ASSERT_NEAR(energy_without_extCharges + energy_extCharges, energy_with_extCharges, 0.0000000001);
  ClassToChangeExternalCharges::clear_external_charges();
}

TEST(forcefield, test_total_gradients_with_external_charges_is_sum)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  double energy_without_extCharges = coords.g();
  ASSERT_NEAR(energy_without_extCharges, 6.5344, 0.0001);     // from tinker, see above
  auto gradients_without_extCharges = coords.g_xyz();

  // setting some random external charges, only scaled_charge is used during calculation
  ClassToChangeExternalCharges::set_external_charges({ {-4, 5, -2, 0.75, double()}, {5, -4, 2, 0.5, double()} });

  energy::interfaces::aco::aco_ff y(&coords);
  y.update();  // initialization of interface
  double energy_with_extCharges = y.g();
  ASSERT_NEAR(energy_with_extCharges, 7.3448092340770454, 0.00001);     // see above
  auto gradients_with_extCharges = coords.g_xyz();

  auto gradients_extCharges = y.part_grad[energy::interfaces::aco::aco_ff::types::EXTERNAL_CHARGES];

  ASSERT_TRUE(is_nearly_equal(gradients_without_extCharges + gradients_extCharges, gradients_with_extCharges, 0.0000000001));
  ClassToChangeExternalCharges::clear_external_charges();
}

#endif
