/**
CAST 3
Purpose: Tests stuff for QM/MM interfaces

@author Susanne Sauer
@version 1.0
*/

#ifdef GOOGLE_MOCK
#include "../coords_io.h"
#include "../tinker_parameters.h"
#include "../qmmm_helperfunctions.h"
#include "../energy_int_aco.h"
#include "gtest/gtest.h"

TEST(qmmm, test_number_of_link_atoms)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
	coords::Coordinates coords(ci->read("test_files/butanol.arc"));

	tinker::parameter::parameters tp;
	tp.from_file("test_files/oplsaa.prm");

	std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };

	auto linkatoms = qmmm_helpers::create_link_atoms(&coords, qm_indizes, tp);
	ASSERT_EQ(linkatoms.size(), 1);
}

TEST(qmmm, test_position_of_link_atom)
{
	std::unique_ptr<coords::input::format> ci(coords::input::new_format());
	coords::Coordinates coords(ci->read("test_files/butanol.arc"));

	tinker::parameter::parameters tp;
	tp.from_file("test_files/oplsaa.prm");

	std::vector<size_t> qm_indizes = { 5,8,9,10,11,12,13,14 };

	auto linkatoms = qmmm_helpers::create_link_atoms(&coords, qm_indizes, tp);
	auto x = linkatoms[0].position.x();
	auto y = linkatoms[0].position.y();
	auto z = linkatoms[0].position.z();

	ASSERT_FLOAT_EQ(x, -1.94975);
	ASSERT_FLOAT_EQ(y, 2.014036);
	ASSERT_FLOAT_EQ(z, -0.19086525);
}

#endif
