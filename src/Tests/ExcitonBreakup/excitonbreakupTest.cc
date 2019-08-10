/**
CAST 3
configurationTest
Purpose: Tests functions in exciton_breakup.h

@author Dustin Kaiser
@version 1.0
*/


#ifdef GOOGLE_MOCK
#include <gtest/gtest.h>

#include <iostream>
#include "../../exciton_breakup_old.h"
#include "../../exciton_breakup.h"

// These tests currently depend on inputfiles

TEST(XB_throws, when_number_of_molecules_dont_match)
{
  std::ofstream schwerpunkt;
  schwerpunkt.open("_tmp_xbtest_massCenterTest.txt", std::ios::out);
  std::size_t gesamtanzahl = 15u;
  schwerpunkt << gesamtanzahl << "\n\n\n";
  schwerpunkt.close();

  Config::set().exbreak.pscnumber = 3;
  Config::set().exbreak.nscnumber = 1;
  EXPECT_ANY_THROW(XB::exciton_breakup(5, 12, Config::get().exbreak.interfaceorientation, "_tmp_xbtest_massCenterTest.txt",
    Config::get().exbreak.nscpairrates, Config::get().exbreak.pscpairexrates, Config::get().exbreak.pscpairchrates, Config::get().exbreak.pnscpairrates));
  
  EXPECT_ANY_THROW(XB::ExcitonBreakup("_tmp_xbtest_massCenterTest.txt",Config::get().exbreak.nscpairrates, Config::get().exbreak.pscpairexrates, \
    Config::get().exbreak.pscpairchrates, Config::get().exbreak.pnscpairrates));


  std::remove("_tmp_xbtest_massCenterTest.txt");
}

TEST(XB_correctly, reads_files_and_stores_raw_data)
{
  std::ofstream file;
  file.open("_tmp_xbtest_massCenterTest.txt", std::ios::out);
  file << "6\n\n1    2.8443712    2.7002553    1.1455320\n2    -2.2574961    6.8744065   24.1365451\n3    1.2881799    9.4028842    3.0530986\n4   -1.2150799   15.1905763   22.3438163\n5   -1.4737771   19.7372300    3.4026387\n6   -0.9276370   21.1267348   20.8234782" << "\n\n\n";
  file.close();
  file.open("_tmp_xbtest_nSC_homodimer.txt", std::ios::out);
  file << "4   5   0.001\n5   6   0.005\n4   6   0.01\n\n";
  file.close();
  file.open("_tmp_xbtest_pscpair_exrates.txt", std::ios::out);
  file << "1   2   0.01\n2   3   0.05\n1   3   0.1\n\n";
  file.close();
  file.open("_tmp_xbtest_pscpair_chargerates.txt", std::ios::out);
  file << "1   2   0.005\n2   3   0.004\n1   3   0.001\n\n";
  file.close();
  file.open("_tmp_xbtest_heterodimer.txt", std::ios::out);
  file << "1   4   0.0008   0.01\n2   5   0.0007   0.015\n\n\n";
  file.close();

  XB::exciton_breakup(3, 3, Config::get().exbreak.interfaceorientation, "_tmp_xbtest_massCenterTest.txt",
    "_tmp_xbtest_nSC_homodimer.txt", "_tmp_xbtest_homodimer_exciton.txt", "_tmp_xbtest_homodimer_ladung.txt", "_tmp_xbtest_heterodimer.txt");

  Config::set().exbreak.pscnumber = 3;
  Config::set().exbreak.nscnumber = 3;
  XB::ExcitonBreakup xb("_tmp_xbtest_massCenterTest.txt","_tmp_xbtest_nSC_homodimer.txt", "_tmp_xbtest_pscpair_exrates.txt", "_tmp_xbtest_pscpair_chargerates.txt", "_tmp_xbtest_heterodimer.txt");

  ASSERT_EQ(xb.x, std::vector<double>({ 0.0, 2.8443712 , -2.2574961 , 1.2881799 ,-1.2150799,-1.4737771,-0.9276370 }));
  ASSERT_EQ(xb.y, std::vector<double>({ 0.0, 2.7002553 , 6.8744065 , 9.4028842 ,15.1905763,19.7372300, 21.1267348 }));
  ASSERT_EQ(xb.z, std::vector<double>({ 0.0, 1.1455320 , 24.1365451 , 3.0530986 ,22.3438163,3.4026387,20.8234782 }));

  std::size_t gesamtanzahl = 6u;

  // Diese Matrix ist gesamtzahl_x_gesamtzahl und beinhaltet allerdings nur Kopplungden der P-SCs, könnte also kleiner gemacht werden...
  ASSERT_EQ(xb.coupling_exciton.at(1).at(2), 0.01);
  ASSERT_EQ(xb.coupling_exciton.at(3).at(2), 0.05);
  ASSERT_EQ(xb.coupling_exciton.at(4).at(4), 0.00);

  ASSERT_EQ(xb.coupling_ladung.at(1).at(2), 0.005);
  ASSERT_EQ(xb.coupling_ladung.at(3).at(2), 0.004);
  ASSERT_EQ(xb.coupling_ladung.at(4).at(4), 0.00);

  ASSERT_EQ(xb.coupling_ct.at(1).at(4), 0.0008);
  ASSERT_EQ(xb.coupling_ct.at(2).at(5), 0.0007);
  ASSERT_EQ(xb.coupling_ct.at(4).at(4), 0.00);

  ASSERT_EQ(xb.coupling_rek.at(1).at(4), 0.01);
  ASSERT_EQ(xb.coupling_rek.at(2).at(5), 0.015);
  ASSERT_EQ(xb.coupling_rek.at(5).at(5), 0.00);

  ASSERT_EQ(xb.coupling_fulleren.at(4).at(5), 0.001);
  ASSERT_EQ(xb.coupling_fulleren.at(5).at(6), 0.005);
  ASSERT_EQ(xb.coupling_fulleren.at(5).at(5), 0.00);
  ASSERT_EQ(xb.coupling_fulleren.at(2).at(2), 0.00);

  std::remove("_tmp_xbtest_massCenterTest.txt");
  std::remove("_tmp_xbtest_pscpair_chargerates.txt");
  std::remove("_tmp_xbtest_pscpair_exrates.txt");
  std::remove("_tmp_xbtest_heterodimer.txt");
}


#endif
