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
	coords::float_type e_cb = energy::interfaces::aco::aco_ff::eQ(10.0, 0.666666);
	ASSERT_EQ(e_cb, 6.66666);
}

TEST(energy_calculations, test_gQ)
{
	coords::float_type dQ;
	energy::interfaces::aco::aco_ff::gQ(10.0, 0.6666666666666, dQ);
	ASSERT_FLOAT_EQ(dQ, -4.4444447);
}

TEST(energy_calculations, test_gQ_fep_energy_without_shift)
{
	Config::set().fep.cshift = 0;
	coords::float_type dQ;
	coords::float_type e_cb = energy::interfaces::aco::aco_ff::gQ_fep(10.0, 1.5, 1, dQ);
	ASSERT_FLOAT_EQ(e_cb, 6.6666665);
}

TEST(energy_calculations, test_gQ_fep_gradient_without_shift)
{
	Config::set().fep.cshift = 0;
	coords::float_type dQ;
	energy::interfaces::aco::aco_ff::gQ_fep(10.0, 1.5, 1, dQ);
	ASSERT_FLOAT_EQ(dQ, -4.4444447);
}

TEST(energy_calculations, test_gQ_fep_energy_with_shift)
{
	Config::set().fep.cshift = 4;
	coords::float_type dQ;
	coords::float_type e_cb = energy::interfaces::aco::aco_ff::gQ_fep(10.0, 1.5, 0.2, dQ);
	ASSERT_FLOAT_EQ(e_cb, 1.2794334);
}

TEST(energy_calculations, test_gQ_fep_gradient_with_shift)
{
	Config::set().fep.cshift = 4;
	coords::float_type dQ;
	energy::interfaces::aco::aco_ff::gQ_fep(10.0, 1.5, 0.2, dQ);
	ASSERT_FLOAT_EQ(dQ, -0.66588628);
}

TEST(energy_calculations, test_eV_Rmin)
{
	coords::float_type e_vdw = energy::interfaces::aco::aco_ff::eV< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 0.666666666);
	ASSERT_FLOAT_EQ(e_vdw, 203.3198041);
}

TEST(energy_calculations, test_eV_sigma)
{
	coords::float_type e_vdw = energy::interfaces::aco::aco_ff::eV< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 0.666666666);
	ASSERT_FLOAT_EQ(e_vdw, 259.506361);
}

TEST(energy_calculations, test_gV_Rmin)
{
	coords::float_type dV;
	energy::interfaces::aco::aco_ff::gV< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 0.666666666, dV);
	ASSERT_FLOAT_EQ(dV, -2076.05);
}

TEST(energy_calculations, test_gV_sigma)
{
	coords::float_type dV;
	energy::interfaces::aco::aco_ff::gV< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 0.666666666, dV);
	ASSERT_FLOAT_EQ(dV, -2300.7971);
}

TEST(energy_calculations, test_gV_fep_Rmin_energy_without_shift)
{
	Config::set().fep.ljshift = 0;
	coords::float_type dV;
	coords::float_type e_vdw = energy::interfaces::aco::aco_ff::gV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 1.5, 1, dV);
	ASSERT_FLOAT_EQ(e_vdw, 203.3198041);
}

TEST(energy_calculations, test_gV_fep_Rmin_gradient_without_shift)
{
	Config::set().fep.ljshift = 0;
	coords::float_type dV;
	energy::interfaces::aco::aco_ff::gV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 1.5, 1, dV);
	ASSERT_FLOAT_EQ(dV, -2076.05);
}

TEST(energy_calculations, test_gV_fep_Rmin_energy_with_shift)
{
	Config::set().fep.ljshift = 4;
	coords::float_type dV;
	coords::float_type e_vdw = energy::interfaces::aco::aco_ff::gV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 1.5, 0.2, dV);
	ASSERT_FLOAT_EQ(e_vdw, -1.1941416);
}

TEST(energy_calculations, test_gV_fep_Rmin_gradient_with_shift)
{
	Config::set().fep.ljshift = 4;
	coords::float_type dV;
	energy::interfaces::aco::aco_ff::gV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 2.0, 1.5, 0.2, dV);
	ASSERT_FLOAT_EQ(dV, 0.24112479);
}

TEST(energy_calculations, test_gV_fep_sigma_energy_without_shift)
{
	Config::set().fep.ljshift = 0;
	coords::float_type dV;
	coords::float_type e_vdw = energy::interfaces::aco::aco_ff::gV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 1.5, 1, dV);
	ASSERT_FLOAT_EQ(e_vdw, 259.506361);
}

TEST(energy_calculations, test_gV_fep_sigma_gradient_without_shift)
{
	Config::set().fep.ljshift = 0;
	coords::float_type dV;
	energy::interfaces::aco::aco_ff::gV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 1.5, 1, dV);
	ASSERT_FLOAT_EQ(dV, -2300.7971);
}

TEST(energy_calculations, test_gV_fep_sigma_energy_with_shift)
{
	Config::set().fep.ljshift = 4;
	coords::float_type dV;
	coords::float_type e_vdw = energy::interfaces::aco::aco_ff::gV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 1.5, 0.2, dV);
	ASSERT_FLOAT_EQ(e_vdw, -0.46367568);
}

TEST(energy_calculations, test_gV_fep_sigma_gradient_with_shift)
{
	Config::set().fep.ljshift = 4;
	coords::float_type dV;
	energy::interfaces::aco::aco_ff::gV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 2.0, 1.5, 0.2, dV);
	ASSERT_FLOAT_EQ(dV, 0.05119307);
}

TEST(energy_calculations, test_gQV_Rmin)
{
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	energy::interfaces::aco::aco_ff::g_QV< ::tinker::parameter::radius_types::R_MIN>(10, 10, 2, 0.66666666, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 6.6666665);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, 203.3198041);
	ASSERT_FLOAT_EQ(dE, -1386.9967);    //gradient 
}

TEST(energy_calculations, test_gQV_SIGMA)
{
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	energy::interfaces::aco::aco_ff::g_QV< ::tinker::parameter::radius_types::SIGMA>(10., 10., 2., 0.66666666, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 6.6666665);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, 259.506361);
	ASSERT_FLOAT_EQ(dE, -1536.8277);    //gradient 
}

TEST(energy_calculations, test_gQV_Rmin_fep_without_shift)
{
	Config::set().fep.ljshift = 0;
	Config::set().fep.cshift = 0;
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	energy::interfaces::aco::aco_ff::g_QV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 10, 2, 1.5, 1, 1, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 6.6666665);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, 203.3198041);
	ASSERT_FLOAT_EQ(dE, -1386.9967);    //gradient 
}

TEST(energy_calculations, test_gQV_sigma_fep_without_shift)
{
	Config::set().fep.ljshift = 0;
	Config::set().fep.cshift = 0;
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	energy::interfaces::aco::aco_ff::g_QV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 10, 2, 1.5, 1, 1, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 6.6666665);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, 259.506361);
	ASSERT_FLOAT_EQ(dE, -1536.8277);    //gradient 
}

TEST(energy_calculations, test_gQV_Rmin_fep_with_shift)
{
	Config::set().fep.ljshift = 4;
	Config::set().fep.cshift = 4;
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	energy::interfaces::aco::aco_ff::g_QV_fep< ::tinker::parameter::radius_types::R_MIN>(10, 10, 2, 1.5, 0.2, 0.2, e_c, e_v, dE);
	ASSERT_FLOAT_EQ(e_c, 1.2794334);    //same as in e_QV
	ASSERT_FLOAT_EQ(e_v, -1.1941416);
	ASSERT_FLOAT_EQ(dE, -0.2831743);    //gradient 
}

TEST(energy_calculations, test_gQV_SIGMA_fep_with_shift)
{
	Config::set().fep.ljshift = 4;
	Config::set().fep.cshift = 4;
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	coords::float_type dE;
	energy::interfaces::aco::aco_ff::g_QV_fep< ::tinker::parameter::radius_types::SIGMA>(10, 10, 2, 1.5, 0.2, 0.2, e_c, e_v, dE);
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

TEST(energy_calculations, test_g_nb_cutoff)
{
	Config::set().energy.cutoff = 5.0;
	Config::set().energy.switchdist = 3.0;

	std::unique_ptr<coords::input::format> ci(coords::input::new_format());
	coords::Coordinates coords(ci->read("test_files/butanol.arc"));

	tinker::parameter::parameters tp;
	tp.from_file("test_files/oplsaa.prm");

	energy::interfaces::aco::aco_ff y(&coords);
	y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC].assign(coords.size(), coords::Cartesian_Point(0.0, 0.0, 0.0));
	y.g_nb<::tinker::parameter::radius_types::SIGMA>();

	// energy values are taken from CAST
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::CHARGE], 2.37532, 0.0001);
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::VDW], 0.630211, 0.0001);

	// this gradient is just taken from CAST, hoping it is correct
	coords::Representation_3D expected_grad = {
	  coords::r3(2.66575 , -0.079883 , 0.0324658),
    coords::r3(0.181848 , -0.404577 , 0.08333),
    coords::r3(0.643512 , 0.425222 , -0.106661),
    coords::r3(0.560345 , 0.456812 , -0.41929),
    coords::r3(2.45419 , 1.19877 , 0.134055),
    coords::r3(-1.62682 , 0.0883042 , 0.535632),
    coords::r3(0.452109 , -0.744505 , 0.173214),
    coords::r3(0.498707 , -0.5261 , 0.660312),
    coords::r3(-4.55026 , -0.450489 , 0.0339252),
    coords::r3(-0.209483 , -0.386559 , -0.844672),
    coords::r3(1.00701 , -1.3501 , -0.391613),
    coords::r3(0.366643 , -0.333948 , -0.105486),
    coords::r3(0.213947 , 1.09561 , 0.49338),
    coords::r3(-0.885455 , 1.01218 , -0.543312),
    coords::r3(-1.77205 , -0.000741292 , 0.264718)
	};

	ASSERT_TRUE(is_nearly_equal(expected_grad, y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC], 0.00001));
}

TEST(energy_calculations, test_g_nb_periodics)
{
	Config::set().energy.cutoff = 5.0;
	Config::set().energy.switchdist = 3.0;
	Config::set().periodics.periodic = true;
	Config::set().periodics.pb_box = { 10.0, 10.0, 10.0 };

	std::unique_ptr<coords::input::format> ci(coords::input::new_format());
	coords::Coordinates coords(ci->read("test_files/butanol.arc"));

	tinker::parameter::parameters tp;
	tp.from_file("test_files/oplsaa.prm");

	energy::interfaces::aco::aco_ff y(&coords);
	y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC].assign(coords.size(), coords::Cartesian_Point(0.0, 0.0, 0.0));
	y.g_nb<::tinker::parameter::radius_types::SIGMA>();

	// energy values are taken from CAST
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::CHARGE], 2.3940066036485312, 0.0000000001);
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::VDW], 0.62811548650874593, 0.0000000001);

	// this gradient is just taken from CAST, hoping it is correct
	coords::Representation_3D expected_grad = {
	  coords::r3(2.66575 , -0.079883 , 0.0324658),
    coords::r3(0.181848 , -0.404577 , 0.08333),
    coords::r3(0.705656 , 0.468727 , -0.157258),
    coords::r3(0.560345 , 0.456812 , -0.41929),
    coords::r3(2.45419 , 1.19877 , 0.134055),
    coords::r3(-1.62682 , 0.0883042 , 0.535632),
    coords::r3(0.452109 , -0.744505 , 0.173214),
    coords::r3(0.498707 , -0.5261 , 0.660312),
    coords::r3(-4.55026 , -0.450489 , 0.0339252),
    coords::r3(-0.209483 , -0.386559 , -0.844672),
    coords::r3(1.00701 , -1.3501 , -0.391613),
    coords::r3(-0.00275321 , -0.311504 , -0.0851016),
    coords::r3(0.213947 , 1.09561 , 0.49338),
    coords::r3(-0.885455 , 1.01218 , -0.543312),
    coords::r3(-1.4648 , -0.0666901 , 0.294931)
	};

	ASSERT_TRUE(is_nearly_equal(expected_grad, y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC], 0.00001));

	// reset config
	Config::set().energy.cutoff = std::numeric_limits<double>::max();
	Config::set().energy.switchdist = std::numeric_limits<double>::max() - 4;
	Config::set().periodics.periodic = false;
}

TEST(energy_calculations, test_g_nb_fep)
{
	std::unique_ptr<coords::input::format> ci(coords::input::new_format());
	coords::Coordinates coords(ci->read("test_files/ethan_FEP.arc"));

	tinker::parameter::parameters tp;
	tp.from_file("test_files/oplsaa.prm");

	// configuration stuff for FEP
	Config::set().md.fep = true;
	coords.getFep().window.resize(1);
	coords.getFep().window[coords.getFep().window[0].step].vout = 0.3;
	coords.getFep().window[coords.getFep().window[0].step].vin = 0.7;
	coords.getFep().window[coords.getFep().window[0].step].eout = 0.3;
	coords.getFep().window[coords.getFep().window[0].step].ein = 0.7;
	Config::set().fep.ljshift = 4;
	Config::set().fep.cshift = 4;

	energy::interfaces::aco::aco_ff y(&coords);
	y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC].assign(coords.size(), coords::Cartesian_Point(0.0, 0.0, 0.0));
	y.g_nb<::tinker::parameter::radius_types::SIGMA>();

	// energy values are taken from CAST
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::CHARGE], 1.9507770300556686, 0.0000000001);
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::VDW], -0.10399513896139564, 0.0000000001);

	// this gradient is just taken from CAST, hoping it is correct
	coords::Representation_3D expected_grad = {
		coords::r3(0 , 0 , 0),
    coords::r3(0 , 0 , 0),
    coords::r3(-0.108651 , 0.1153 , -0.0179165),
    coords::r3(-0.198244 , 0.0575492 , 0.00268284),
    coords::r3(-0.11509 , 0.103354 , 0.0540107),
    coords::r3(0 , 0 , 0),
    coords::r3(0.104161 , -0.190225 , 0.0409208),
    coords::r3(0.111436 , -0.158585 , -0.10326),
    coords::r3(0.0376066 , 0.0361289 , 0.0249788),
    coords::r3(0.0592395 , 0.0176594 , 0.00636159),
    coords::r3(0.0393275 , 0.0421409 , -0.00662291),
    coords::r3(-0.0503409 , -0.00816881 , -0.00375402),
    coords::r3(0.120556 , -0.015154 , 0.00259893)
	};

	ASSERT_TRUE(is_nearly_equal(expected_grad, y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC], 0.00001));
}

TEST(energy_calculations, test_g_nb_fep_diff_window)
{
	std::unique_ptr<coords::input::format> ci(coords::input::new_format());
	coords::Coordinates coords(ci->read("test_files/ethan_FEP.arc"));

	tinker::parameter::parameters tp;
	tp.from_file("test_files/oplsaa.prm");

	// configuration stuff for FEP
	Config::set().md.fep = true;
	coords.getFep().window.resize(1);
	coords.getFep().window[coords.getFep().window[0].step].vout = 0.8;
	coords.getFep().window[coords.getFep().window[0].step].vin = 0.2;
	coords.getFep().window[coords.getFep().window[0].step].eout = 0.8;
	coords.getFep().window[coords.getFep().window[0].step].ein = 0.2;
	Config::set().fep.ljshift = 4;
	Config::set().fep.cshift = 4;

	energy::interfaces::aco::aco_ff y(&coords);
	y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC].assign(coords.size(), coords::Cartesian_Point(0.0, 0.0, 0.0));
	y.g_nb<::tinker::parameter::radius_types::SIGMA>();

	// energy values are taken from CAST
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::CHARGE], 1.9511434909014964, 0.0000000001);
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::VDW], -0.087852330065593523, 0.0000000001);

	// this gradient is just taken from CAST, hoping it is correct
	coords::Representation_3D expected_grad = {
    coords::r3(0 , 0 , 0),
    coords::r3(0 , 0 , 0),
    coords::r3(-0.0269377 , 0.0274031 , -0.00656142),
    coords::r3(-0.0403453 , 0.0100977 , -7.99348e-05),
    coords::r3(-0.0282304 , 0.0227154 , 0.014776),
    coords::r3(0 , 0 , 0),
    coords::r3(-0.184335 , -0.240822 , 0.0359661),
    coords::r3(-0.175923 , -0.210819 , -0.138969),
    coords::r3(0.134056 , 0.147149 , 0.0777395),
    coords::r3(0.287224 , 0.103537 , 0.0323126),
    coords::r3(0.142217 , 0.169632 , -0.00443942),
    coords::r3(-0.142112 , -0.0248893 , -0.0114263),
    coords::r3(0.0343861 , -0.00400397 , 0.00068158)
	};

	ASSERT_TRUE(is_nearly_equal(expected_grad, y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC], 0.00001));
}

TEST(energy_calculations, test_g_nb_fep_cutoff)
{
	Config::set().energy.cutoff = 5.0;
	Config::set().energy.switchdist = 3.0;

	std::unique_ptr<coords::input::format> ci(coords::input::new_format());
	coords::Coordinates coords(ci->read("test_files/ethan_FEP.arc"));

	tinker::parameter::parameters tp;
	tp.from_file("test_files/oplsaa.prm");

	// configuration stuff for FEP
	Config::set().md.fep = true;
	coords.getFep().window.resize(1);
	coords.getFep().window[coords.getFep().window[0].step].vout = 0.3;
	coords.getFep().window[coords.getFep().window[0].step].vin = 0.7;
	coords.getFep().window[coords.getFep().window[0].step].eout = 0.3;
	coords.getFep().window[coords.getFep().window[0].step].ein = 0.7;
	Config::set().fep.ljshift = 4;
	Config::set().fep.cshift = 4;

	energy::interfaces::aco::aco_ff y(&coords);
	y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC].assign(coords.size(), coords::Cartesian_Point(0.0, 0.0, 0.0));
	y.g_nb<::tinker::parameter::radius_types::SIGMA>();

	// energy values are taken from CAST
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::CHARGE], 0.94546443146117964, 0.0000000001);
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::VDW], -0.10334100850688555, 0.0000000001);

	// this gradient is just taken from CAST, hoping it is correct
	coords::Representation_3D expected_grad = {
		coords::r3(0 , 0 , 0),
    coords::r3(0 , 0 , 0),
    coords::r3(-0.142275 , 0.148968 , -0.0265193),
    coords::r3(-0.243393 , 0.0697826 , 0.002701),
    coords::r3(-0.150283 , 0.130839 , 0.07262),
    coords::r3(0 , 0 , 0),
    coords::r3(0.124825 , -0.243775 , 0.0540208),
    coords::r3(0.133946 , -0.203267 , -0.134253),
    coords::r3(0.0512833 , 0.048941 , 0.0341772),
    coords::r3(0.0784282 , 0.023734 , 0.00851389),
    coords::r3(0.0536236 , 0.0571385 , -0.00928494),
    coords::r3(-0.067975 , -0.0114413 , -0.00513211),
    coords::r3(0.161821 , -0.0209211 , 0.00315653)
	};

	ASSERT_TRUE(is_nearly_equal(expected_grad, y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC], 0.00001));
}

TEST(energy_calculations, test_g_nb_fep_periodic)
{
	Config::set().energy.cutoff = 2.5;
	Config::set().energy.switchdist = 2.0;
	Config::set().periodics.periodic = true;
	Config::set().periodics.pb_box = { 5.0, 5.0, 5.0 };

	std::unique_ptr<coords::input::format> ci(coords::input::new_format());
	coords::Coordinates coords(ci->read("test_files/ethan_FEP.arc"));

	tinker::parameter::parameters tp;
	tp.from_file("test_files/oplsaa.prm");

	// configuration stuff for FEP
	Config::set().md.fep = true;
	coords.getFep().window.resize(1);
	coords.getFep().window[coords.getFep().window[0].step].vout = 0.3;
	coords.getFep().window[coords.getFep().window[0].step].vin = 0.7;
	coords.getFep().window[coords.getFep().window[0].step].eout = 0.3;
	coords.getFep().window[coords.getFep().window[0].step].ein = 0.7;
	Config::set().fep.ljshift = 4;
	Config::set().fep.cshift = 4;

	energy::interfaces::aco::aco_ff y(&coords);
	y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC].assign(coords.size(), coords::Cartesian_Point(0.0, 0.0, 0.0));
	y.g_nb<::tinker::parameter::radius_types::SIGMA>();

	// energy values are taken from CAST
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::CHARGE], 0.087038394018264173, 0.0000000001);
	ASSERT_NEAR(y.part_energy[energy::interfaces::aco::aco_ff::types::VDW], 0.081725523622103696, 0.0000000001);

	// this gradient is just taken from CAST, hoping it is correct
	coords::Representation_3D expected_grad = {
    coords::r3(0 , 0 , 0),
    coords::r3(0 , 0 , 0),
    coords::r3(0 , 0 , 0),
    coords::r3(0.356422 , -0.152445 , -0.0690268),
    coords::r3(5.14282e-05 , -7.60584e-05 , -1.39899e-05),
    coords::r3(0 , 0 , 0),
    coords::r3(0.000363331 , 0.000366175 , -0.00015517),
    coords::r3(0.00149822 , -9.69581e-05 , -0.000270525),
    coords::r3(0 , 0 , 0),
    coords::r3(-0.0751791 , -0.0262904 , 0.00247108),
    coords::r3(-7.45337e-05 , -0.000127216 , -2.74298e-05),
    coords::r3(0.0743824 , 0.0259395 , -0.00249897),
    coords::r3(-0.357464 , 0.15273 , 0.0695218)
	};

	ASSERT_TRUE(is_nearly_equal(expected_grad, y.part_grad[energy::interfaces::aco::aco_ff::types::VDWC], 0.00001));

	// reset config
	Config::set().energy.cutoff = std::numeric_limits<double>::max();
	Config::set().energy.switchdist = std::numeric_limits<double>::max() - 4;
	Config::set().periodics.periodic = false;
}

#endif
