/**
CAST 3
pmf_ic_prep.h
Purpose: preparation for PMF-IC calculation (enhanced umbrella)

@author Susanne Sauer
@version 1.0
*/
#pragma once
#include"coords.h"
#include"coords_io.h"

/**class to perform preparation for PMF-IC*/
class pmf_ic_prep
{
public:

  /**constructor
  @param coords: coordinates object
  @param ci: pointer to input format that contains all given structures
  @param outfilename: name of the outputfile
  @param splinefilename: name of the file where spline is written to*/
  pmf_ic_prep(coords::Coordinates& coords, coords::input::format& ci, std::string const& outfilename, std::string const& splinefilename);

  /**function that performs PMF-IC preparation*/
  void run();

private:

  // input variables
  coords::Coordinates coordobj;
  coords::input::format* coord_input;
  std::string outfilename;
  std::string splinefilename;

  // variables that are calculated during PMF-IC
  std::vector<double> xis;
  std::vector<double> zs;
  std::vector<double> E_HLs;
  std::vector<double> E_LLs;
  std::vector<double> deltaEs;

  /**calculates values for xi (reaction coordinate), z (mapped reaction coordinate) and E_HL (high level energy) for every structure
  stores them into member variables xis, zs and E_HLs respectively*/
  void calc_xis_zs_and_E_HLs();
  /**calculates low level energies for every structure and stores them into E_LLs*/
  void calc_E_LLs();
  /**calculates energy differences HL - LL for every structure and stores them into deltaEs*/
  void calc_deltaEs();
  /**writes outputfile*/
  void write_to_file();
  /**writes splinefile*/
  void write_spline();
};