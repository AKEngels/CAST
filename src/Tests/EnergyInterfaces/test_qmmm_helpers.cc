/**
CAST 3
Purpose: Tests stuff for QM/MM interfaces

@author Susanne Sauer
@version 1.0
*/

#ifdef GOOGLE_MOCK
#include "../../coords_io.h"
#include "../../tinker_parameters.h"
#include "../../qmmm_helperfunctions.h"
#include "../../energy_int_aco.h"
#include <gtest/gtest.h>

// tests use the test system butanol.arc

TEST(qmmm, test_number_of_link_atoms)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };

  auto linkatoms = qmmm_helpers::create_link_atoms(qm_indizes, &coords, tp, { 85 });
  ASSERT_EQ(linkatoms.size(), 1);
}

TEST(qmmm, test_position_of_link_atom)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };

  auto linkatoms = qmmm_helpers::create_link_atoms(qm_indizes, &coords, tp, { 85 });
  auto x = linkatoms[0].position.x();
  auto y = linkatoms[0].position.y();
  auto z = linkatoms[0].position.z();

  ASSERT_FLOAT_EQ(x, -1.94975);
  ASSERT_FLOAT_EQ(y, 2.014036);
  ASSERT_FLOAT_EQ(z, -0.19086525);
}

TEST(qmmm, test_gradient_of_link_atom)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };

  auto linkatom = qmmm_helpers::create_link_atoms(qm_indizes, &coords, tp, { 85 })[0];

  coords::r3 g_la, g_qm, g_mm;
  g_la.x() = 4.81805;
  g_la.y() = 4.25433;
  g_la.z() = 0.482165;

  qmmm_helpers::calc_link_atom_grad(linkatom, g_la, &coords, g_qm, g_mm);

  auto x_qm = g_qm.x();
  auto y_qm = g_qm.y();
  auto z_qm = g_qm.z();

  auto x_mm = g_mm.x();
  auto y_mm = g_mm.y();
  auto z_mm = g_mm.z();

  ASSERT_FLOAT_EQ(x_qm, 4.75778);
  ASSERT_FLOAT_EQ(y_qm, 1.3922169);
  ASSERT_FLOAT_EQ(z_qm, 1.1590456);
  ASSERT_FLOAT_EQ(x_mm, 0.060268376);
  ASSERT_FLOAT_EQ(y_mm, 2.862113);
  ASSERT_FLOAT_EQ(z_mm, -0.6768806);
}

TEST(qmmm, test_get_mm_atoms)
{
  Config::set().energy.qmmm.qm_systems = { { 5,8,9,10,11,12,13,14 } };
  auto mm_atoms = qmmm_helpers::get_mm_atoms(15);

  std::vector<size_t> mm_atoms_ideal = { 0, 1, 2, 3, 4, 6, 7 };
  ASSERT_EQ(mm_atoms, mm_atoms_ideal);
}

TEST(qmmm, test_make_new_indices_for_qm_atoms)
{
  std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };
  auto new_indizes = qmmm_helpers::make_new_indices(qm_indizes, 15);

  std::vector<size_t> new_ideal;
  new_ideal.resize(15);

  new_ideal[5] = 0;
  new_ideal[8] = 1;
  new_ideal[9] = 2;
  new_ideal[10] = 3;
  new_ideal[11] = 4;
  new_ideal[12] = 5;
  new_ideal[13] = 6;
  new_ideal[14] = 7;

  ASSERT_EQ(new_indizes, new_ideal);
}

TEST(qmmm, test_make_new_indices_for_mm_atoms)
{
  std::vector<size_t> mm_indizes = { 0, 1, 2, 3, 4, 6, 7 };
  auto new_indizes = qmmm_helpers::make_new_indices(mm_indizes, 15);

  std::vector<size_t> new_ideal;
  new_ideal.resize(15);

  new_ideal[0] = 0;
  new_ideal[1] = 1;
  new_ideal[2] = 2;
  new_ideal[3] = 3;
  new_ideal[4] = 4;
  new_ideal[6] = 5;
  new_ideal[7] = 6;

  ASSERT_EQ(new_indizes, new_ideal);
}

TEST(qmmm, test_small_coords_for_mm_system)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  std::vector<size_t> mm_indizes = { 0, 1, 2, 3, 4, 6, 7 };
  auto new_indizes = qmmm_helpers::make_new_indices(mm_indizes, 15);

  auto small_coords = qmmm_helpers::make_small_coords(&coords, mm_indizes, new_indizes, config::interface_types::OPLSAA, "MM System: ");

  ASSERT_EQ(small_coords.size(), 7);
  ASSERT_EQ(small_coords.atoms(1).bonds().size(), 3);  // atom 2 which has the QM/MM bond in original system now has only 3 bonding partners
}

TEST(qmmm, test_small_coords_for_qm_system)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };
  auto new_indizes = qmmm_helpers::make_new_indices(qm_indizes, 15);

  auto linkatoms = qmmm_helpers::create_link_atoms(qm_indizes, &coords, tp, { 85 });
  auto small_coords = qmmm_helpers::make_small_coords(&coords, qm_indizes, new_indizes, config::interface_types::DFTB, "QM system: ", false, linkatoms);

  ASSERT_EQ(small_coords.size(), 9);  // +1 link atom
  ASSERT_EQ(small_coords.atoms(8).bonds()[0], 0);  // link atom has one bonding partner, this is atom 6 in original system, first atom in small_coords object
}

TEST(qmmm, test_get_index_of_qmcenter)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };

  auto index_of_QMcenter = qmmm_helpers::get_index_of_QM_center(0, qm_indizes, &coords);
  ASSERT_EQ(index_of_QMcenter, 8);
}

TEST(qmmm, test_add_external_charges_M1)
{
  Config::set().energy.qmmm.zerocharge_bonds = 1;    // default

  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };
  auto linkatoms = qmmm_helpers::create_link_atoms(qm_indizes, &coords, tp, { 85 });

  auto charges = coords.energyinterface()->charges();
  std::vector<size_t> all_indizes = range(coords.size());

  std::vector<int> result;
  qmmm_helpers::add_external_charges(qm_indizes, charges, all_indizes, linkatoms, result, &coords, 8);
  ASSERT_EQ(result.size(), 6);
}

TEST(qmmm, test_add_external_charges_M2)
{
  Config::set().energy.qmmm.zerocharge_bonds = 2;

  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };
  auto linkatoms = qmmm_helpers::create_link_atoms(qm_indizes, &coords, tp, { 85 });

  auto charges = coords.energyinterface()->charges();
  std::vector<size_t> all_indizes = range(coords.size());

  std::vector<int> result;
  qmmm_helpers::add_external_charges(qm_indizes, charges, all_indizes, linkatoms, result, &coords, 8);
  ASSERT_EQ(result.size(), 3);
}

TEST(qmmm, test_add_external_charges_M3)
{
  Config::set().energy.qmmm.zerocharge_bonds = 3;

  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  std::vector<size_t> qm_indizes = { 0, 1, 2, 3, 4, 6, 7 };   // switch QM and MM atoms
  auto linkatoms = qmmm_helpers::create_link_atoms(qm_indizes, &coords, tp, { 85 });

  auto charges = coords.energyinterface()->charges();
  std::vector<size_t> all_indizes = range(coords.size());

  std::vector<int> result;
  qmmm_helpers::add_external_charges(qm_indizes, charges, all_indizes, linkatoms, result, &coords, 8);
  ASSERT_EQ(result.size(), 1);
}

TEST(qmmm, test_add_external_charges_onlySE)
{
  Config::set().energy.qmmm.zerocharge_bonds = 1;    // default

  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read("test_files/butanol.arc"));

  tinker::parameter::parameters tp;
  tp.from_file("test_files/oplsaa.prm");

  std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };
  std::vector<size_t> qmse_indizes = { 1,5,6,7,8,9,10,11,12,13,14 };      // SE atoms: CH2 group next to QM region
  std::vector<size_t> all_indizes = range(coords.size());
  auto linkatoms = qmmm_helpers::create_link_atoms(qm_indizes, &coords, tp, { 85 });

  auto all_charges = coords.energyinterface()->charges();
  std::vector<double> charges;
  for (auto i : qmse_indizes)
  {
    charges.push_back(all_charges[i]);
  }

  std::vector<int> result;
  qmmm_helpers::add_external_charges(qm_indizes, charges, qmse_indizes, linkatoms, result, &coords, 8);
  ASSERT_EQ(result.size(), 2);  // only charges for SE atoms (7 and 8)
}

#endif
