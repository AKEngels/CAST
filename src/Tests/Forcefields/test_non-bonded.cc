/**
CAST 3
Purpose: Tests energy and gradient functions for non-bonding atom pairs

@author Susanne Sauer
@version 1.0
*/


#ifdef GOOGLE_MOCK

#include "../../energy_int_aco.h"
#include "../../tinker_parameters.h"
#include "../../coords_io.h"
#include <gtest/gtest.h>

TEST(energy_calculations, test_eQ)
{
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type e_cb = y.eQ(10.0, 0.666666);
	ASSERT_EQ(e_cb, 6.66666);
}

TEST(energy_calculations, test_gQ)
{
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dQ;
	y.gQ(10.0, 0.6666666666666, dQ);
	ASSERT_FLOAT_EQ(dQ, -4.4444447);
}

TEST(energy_calculations, test_gQ_fep_energy_without_shift)
{
	Config::set().fep.cshift = 0;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dQ;
	coords::float_type e_cb = y.gQ_fep(10.0, 1.5, 1, dQ);
	ASSERT_FLOAT_EQ(e_cb, 6.6666665);
}

TEST(energy_calculations, test_gQ_fep_gradient_without_shift)
{
	Config::set().fep.cshift = 0;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dQ;
	y.gQ_fep(10.0, 1.5, 1, dQ);
	ASSERT_FLOAT_EQ(dQ, -4.4444447);
}

TEST(energy_calculations, test_gQ_fep_energy_with_shift)
{
	Config::set().fep.cshift = 4;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dQ;
	coords::float_type e_cb = y.gQ_fep(10.0, 1.5, 0.2, dQ);
	ASSERT_FLOAT_EQ(e_cb, 1.2794334);
}

TEST(energy_calculations, test_gQ_fep_gradient_with_shift)
{
	Config::set().fep.cshift = 4;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dQ;
	y.gQ_fep(10.0, 1.5, 0.2, dQ);
	ASSERT_FLOAT_EQ(dQ, -0.66588628);
}

TEST(energy_calculations, test_eV_Rmin)
{
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type e_vdw = y.eV< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 0.666666666);
	ASSERT_FLOAT_EQ(e_vdw, 203.3198041);
}

TEST(energy_calculations, test_eV_sigma)
{
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type e_vdw = y.eV< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 0.666666666);
	ASSERT_FLOAT_EQ(e_vdw, 259.506361);
}

TEST(energy_calculations, test_gV_Rmin)
{
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dV;
	y.gV< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 0.666666666, dV);
	ASSERT_FLOAT_EQ(dV, -2076.05);
}

TEST(energy_calculations, test_gV_sigma)
{
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dV;
	y.gV< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 0.666666666, dV);
	ASSERT_FLOAT_EQ(dV, -2300.7971);
}

TEST(energy_calculations, test_gV_fep_Rmin_energy_without_shift)
{
	Config::set().fep.ljshift = 0;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dV;
	coords::float_type e_vdw = y.gV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 1.5, 1, dV);
	ASSERT_FLOAT_EQ(e_vdw, 203.3198041);
}

TEST(energy_calculations, test_gV_fep_Rmin_gradient_without_shift)
{
	Config::set().fep.ljshift = 0;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dV;
	y.gV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 1.5, 1, dV);
	ASSERT_FLOAT_EQ(dV, -2076.05);
}

TEST(energy_calculations, test_gV_fep_Rmin_energy_with_shift)
{
	Config::set().fep.ljshift = 4;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dV;
	coords::float_type e_vdw = y.gV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 1.5, 0.2, dV);
	ASSERT_FLOAT_EQ(e_vdw, -1.1941416);
}

TEST(energy_calculations, test_gV_fep_Rmin_gradient_with_shift)
{
	Config::set().fep.ljshift = 4;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dV;
	y.gV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 1.5, 0.2, dV);
	ASSERT_FLOAT_EQ(dV, 0.24112479);
}

TEST(energy_calculations, test_gV_fep_sigma_energy_without_shift)
{
	Config::set().fep.ljshift = 0;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dV;
	coords::float_type e_vdw = y.gV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 1.5, 1, dV);
	ASSERT_FLOAT_EQ(e_vdw, 259.506361);
}

TEST(energy_calculations, test_gV_fep_sigma_gradient_without_shift)
{
	Config::set().fep.ljshift = 0;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dV;
	y.gV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 1.5, 1, dV);
	ASSERT_FLOAT_EQ(dV, -2300.7971);
}

TEST(energy_calculations, test_gV_fep_sigma_energy_with_shift)
{
	Config::set().fep.ljshift = 4;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dV;
	coords::float_type e_vdw = y.gV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 1.5, 0.2, dV);
	ASSERT_FLOAT_EQ(e_vdw, -0.46367568);
}

TEST(energy_calculations, test_gV_fep_sigma_gradient_with_shift)
{
	Config::set().fep.ljshift = 4;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type dV;
	y.gV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 1.5, 0.2, dV);
	ASSERT_FLOAT_EQ(dV, 0.05119307);
}

TEST(energy_calculations, test_gQV_Rmin)
{
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	y.g_QV< ::tinker::parameter::radius_types::R_MIN>(10, 10, 2, 0.66666666, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 6.6666665);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, 203.3198041);
	ASSERT_FLOAT_EQ(dE, -1386.9967);    //gradient 
}

TEST(energy_calculations, test_gQV_SIGMA)
{
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	y.g_QV< ::tinker::parameter::radius_types::SIGMA>(10., 10., 2., 0.66666666, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 6.6666665);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, 259.506361);
	ASSERT_FLOAT_EQ(dE, -1536.8277);    //gradient 
}

TEST(energy_calculations, test_gQV_Rmin_fep_without_shift)
{
	Config::set().fep.ljshift = 0;
	Config::set().fep.cshift = 0;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	y.g_QV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 10, 2, 1.5, 1, 1, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 6.6666665);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, 203.3198041);
	ASSERT_FLOAT_EQ(dE, -1386.9967);    //gradient 
}

TEST(energy_calculations, test_gQV_sigma_fep_without_shift)
{
	Config::set().fep.ljshift = 0;
	Config::set().fep.cshift = 0;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	y.g_QV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 10, 2, 1.5, 1, 1, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 6.6666665);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, 259.506361);
	ASSERT_FLOAT_EQ(dE, -1536.8277);    //gradient 
}

TEST(energy_calculations, test_gQV_Rmin_fep_with_shift)
{
	Config::set().fep.ljshift = 4;
	Config::set().fep.cshift = 4;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	y.g_QV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 10, 2, 1.5, 0.2, 0.2, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 1.2794334);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, -1.1941416);
	ASSERT_FLOAT_EQ(dE, -0.2831743);    //gradient 
}

TEST(energy_calculations, test_gQV_SIGMA_fep_with_shift)
{
	Config::set().fep.ljshift = 4;
	Config::set().fep.cshift = 4;
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	y.g_QV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 10, 2, 1.5, 0.2, 0.2, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 1.2794334);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, -0.46367568);
	ASSERT_FLOAT_EQ(dE, -0.4097954);    //gradient 
}

TEST(energy_calculations, test_g_nb)
{
	std::unique_ptr<coords::input::format> ci(coords::input::new_format());
	coords::Coordinates coords(ci->read("test_files/butanol.arc"));

	tinker::parameter::parameters tp;
	tp.from_file("test_files/oplsaa.prm");

	energy::interfaces::aco::aco_ff y(&coords);
	y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC].assign(coords.size(), coords::Cartesian_Point(0.0, 0.0, 0.0));
	y.g_nb<::tinker::parameter::radius_types::SIGMA>();

	// energy values are compared with those from tinker
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::CHARGE], 3.7139, 0.0001);
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::VDW], 0.4742, 0.0001);

	// this gradient is just taken from CAST, hoping it is correct
	coords::Representation_3D expected_grad = {
		coords::r3(2.91329 , -0.180913 , 0.120984),
	  coords::r3(0.216802 , -0.440519 , 0.0619539),
	  coords::r3(0.411764 , 0.41452 , -0.135731),
	  coords::r3(0.453379 , 0.460062 , -0.399239),
	  coords::r3(2.38148 , 1.1689 , 0.116752),
	  coords::r3(-1.35075 , -0.041838 , 0.43091),
	  coords::r3(0.323918 , -0.579819 , 0.130951),
	  coords::r3(0.381201 , -0.416475 , 0.539308),
	  coords::r3(-4.45636 , -0.44424 , 0.0492423),
	  coords::r3(-0.16116 , -0.418916 , -0.758331),
	  coords::r3(0.870319 , -1.31957 , -0.275401),
	  coords::r3(0.512139 , -0.0798732 , -0.0556127),
	  coords::r3(-0.00170722 , 0.943285 , 0.452813),
	  coords::r3(-1.04163 , 0.902286 , -0.485219),
	  coords::r3(-1.45269 , 0.0331112 , 0.20662)
	};

	ASSERT_TRUE(is_nearly_equal(expected_grad, y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC], 0.00001));
}

#endif
