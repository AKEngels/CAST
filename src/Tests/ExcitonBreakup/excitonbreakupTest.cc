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
#include "../../exciton_breakup.h"

void setupTestFiles()
{
  std::ofstream file;
  file.open("_tmp_xbtest_massCenterTest.txt", std::ios::out);
  file << "6\n\n1    0.0000    0.0000    0.0000\n2    5.0000    0.0000   0.0000\n3    10.0000    0.0000    0.0000\n4   15.0000   0.0000   0.0000\n5   20.0000   0.0000    0.0000\n6   25.0000   0.0000   0.0000" << "\n\n\n";
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
}

void cleanupTestFiles()
{
  std::remove("_tmp_xbtest_massCenterTest.txt");
  std::remove("_tmp_xbtest_pscpair_chargerates.txt");
  std::remove("_tmp_xbtest_pscpair_exrates.txt");
  std::remove("_tmp_xbtest_heterodimer.txt");
}

TEST(XB_throws, when_number_of_molecules_dont_match)
{
  std::ofstream schwerpunkt;
  schwerpunkt.open("_tmp_xbtest_massCenterTest.txt", std::ios::out);
  std::size_t gesamtanzahl = 15u;
  schwerpunkt << gesamtanzahl << "\n\n\n";
  schwerpunkt.close();

  Config::set().exbreak.pscnumber = 3;
  Config::set().exbreak.nscnumber = 1;
  
  EXPECT_ANY_THROW(XB::ExcitonBreakup("_tmp_xbtest_massCenterTest.txt",Config::get().exbreak.nscpairrates, Config::get().exbreak.pscpairexrates, \
    Config::get().exbreak.pscpairchrates, Config::get().exbreak.pnscpairrates));


  std::remove("_tmp_xbtest_massCenterTest.txt");
}

TEST(XB_correctly, reads_files_and_stores_raw_data)
{
  cleanupTestFiles();
  setupTestFiles();

  Config::set().exbreak.pscnumber = 3;
  Config::set().exbreak.nscnumber = 3;
  XB::ExcitonBreakup xb("_tmp_xbtest_massCenterTest.txt","_tmp_xbtest_nSC_homodimer.txt", "_tmp_xbtest_pscpair_exrates.txt", "_tmp_xbtest_pscpair_chargerates.txt", "_tmp_xbtest_heterodimer.txt");

  ASSERT_EQ(xb.x, std::vector<double>({ 0.0, 0.0, 5.0 , 10.0, 15.0, 20.0, 25.0 }));
  ASSERT_EQ(xb.y, std::vector<double>({ 0.0, 0.0 ,0.0 , 0.0 , 0.0 , 0.0 , 0.0  }));
  ASSERT_EQ(xb.z, std::vector<double>({ 0.0,0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0  }));

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

  cleanupTestFiles();
}

TEST(XB_correctly, identifies_startingpoints_independent_Of_orientation)
{
  setupTestFiles();


  Config::set().exbreak.pscnumber = 3;
  Config::set().exbreak.nscnumber = 3;
  XB::ExcitonBreakup xb("_tmp_xbtest_massCenterTest.txt", "_tmp_xbtest_nSC_homodimer.txt", "_tmp_xbtest_pscpair_exrates.txt", "_tmp_xbtest_pscpair_chargerates.txt", "_tmp_xbtest_heterodimer.txt");

  std::size_t numPoints = 0u;
  std::vector <std::size_t> vecOfStartingPoints;
  vecOfStartingPoints = xb.calculateStartingpoints('x', numPoints, 0.5);

  // Diese Matrix ist gesamtzahl_x_gesamtzahl und beinhaltet allerdings nur Kopplungden der P-SCs, könnte also kleiner gemacht werden...
  ASSERT_EQ(vecOfStartingPoints.at(0) , 0u);
  ASSERT_EQ(vecOfStartingPoints.at(1) , 1u);
  ASSERT_EQ(vecOfStartingPoints.at(2) , 2u);
  ASSERT_EQ(vecOfStartingPoints.at(3) , 0u);
  ASSERT_EQ(vecOfStartingPoints.size(), 4u);
  ASSERT_EQ(numPoints, 2u);

  std::remove("_tmp_xbtest_massCenterTest.txt");
  std::ofstream file;
  file.open("_tmp_xbtest_massCenterTest.txt", std::ios::out);
  file << "6\n\n1    0.0000    0.0000    0.0000\n2    5.0000    0.0000   0.0000\n3    10.0000    0.0000    0.0000\n4   -5.0000   0.0000   0.0000\n5   -10.0000   0.0000    0.0000\n6   -15.0000   0.0000   0.0000" << "\n\n\n";
  file.close();

  XB::ExcitonBreakup xb2("_tmp_xbtest_massCenterTest.txt", "_tmp_xbtest_nSC_homodimer.txt", "_tmp_xbtest_pscpair_exrates.txt", "_tmp_xbtest_pscpair_chargerates.txt", "_tmp_xbtest_heterodimer.txt");
  vecOfStartingPoints = xb2.calculateStartingpoints('x', numPoints, 0.5);
  ASSERT_EQ(vecOfStartingPoints.at(0), 0u);
  ASSERT_EQ(vecOfStartingPoints.at(1), 2u);
  ASSERT_EQ(vecOfStartingPoints.at(2), 3u);
  ASSERT_EQ(vecOfStartingPoints.at(3), 0u);
  ASSERT_EQ(vecOfStartingPoints.size(), 4u);
  ASSERT_EQ(numPoints, 2u);
  cleanupTestFiles();
}


#endif
