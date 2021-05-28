/**
CAST 3
pmf_ic_prep.h
Purpose: preparation for PMF-IC calculation (enhanced umbrella)

@author Susanne Sauer
@version 1.0
*/
#pragma once
#include <optional>

#include"coords.h"
#include"coords_io.h"

#include "spline.h"

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
  /**dimension of the PMF (can be 1 or 2)*/
  std::size_t dimension;

  // variables that are calculated during PMF-IC

  /**xi value for every structure (only used in 1D)*/
  std::vector<double> xis;  
  /**z value for every structure (only used in 1D)*/
  std::vector<double> zs;
  /**high level energy for every structure*/
  std::vector<double> E_HLs;
  /**low level energy for every structure*/
  std::vector<double> E_LLs;
  /**energy difference for every structure*/
  std::vector<double> deltaEs;
  // only for 2D
  std::vector < std::pair<double, double>> xi_2d; // used instead of 'xis'
  std::vector < std::pair<double, double>> z_2d;  // used instead of 'zs'

  /** Mapper for first reaction coordinate */
  XiToZMapper mapper1;
  /** Mapper for second reaction coordinate, not always needed */
  std::optional<XiToZMapper> mapper2;

  /**calculates values for xi (reaction coordinate), z (mapped reaction coordinate) and E_HL (high level energy) for every structure
  stores them into member variables xis, zs and E_HLs respectively (for 1D) or in xi_2d, z_2d and E_HLs (for 2D)*/
  void calc_xis_zs_and_E_HLs();
  /**calculates low level energies for every structure and stores them into E_LLs*/
  void calc_E_LLs();
  /**calculates energy differences HL - LL for every structure and stores them into deltaEs*/
  void calc_deltaEs();
  /**writes outputfile*/
  void write_to_file();
  /**writes splinefile for 1d spline*/
  void write_spline_1d();
  /**writes splinefile for 2d spline*/
  void write_spline_2d();
};