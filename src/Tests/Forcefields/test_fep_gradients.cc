/**
CAST 3
Purpose: Tests energy and gradient functions for non-bonding atom pairs

@author Susanne Sauer
@version 1.0
*/


#ifdef GOOGLE_MOCK

#include "../../energy_int_aco.h"
#include "../../tinker_parameters.h"
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

TEST(energy_calculations, test_eQV_Rmin)  
{
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
  coords::float_type e_c(0.0);
  coords::float_type e_v(0.0);
  y.e_QV<::tinker::parameter::radius_types::R_MIN>(10, 10, 2, 0.66666666, e_c, e_v);
	ASSERT_FLOAT_EQ(e_c, 6.6666665);   
	ASSERT_FLOAT_EQ(e_v, 203.3198041);
}

TEST(energy_calculations, test_eQV_SIGMA)
{
	energy::interfaces::aco::aco_ff y(&coords::Coordinates());
	coords::float_type e_c(0.0);
	coords::float_type e_v(0.0);
	y.e_QV< ::tinker::parameter::radius_types::SIGMA>(10., 10., 2., 0.66666666, e_c, e_v);
	ASSERT_FLOAT_EQ(e_c, 6.6666665);   
	ASSERT_FLOAT_EQ(e_v, 259.506361);
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

#endif
