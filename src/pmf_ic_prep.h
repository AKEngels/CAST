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

#include "InternalCoordinates/PrimitiveInternalCoordinates.h"

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
  /**high level energy for every structure*/
  std::vector<double> E_HLs;
  /**low level energy for every structure*/
  std::vector<double> E_LLs;
  /**energy difference for every structure*/
  std::vector<double> deltaEs;
  /** energy gradient wrt reaction coordinate (difference between HL and LL method) */
  std::vector<double> grad_Es;
  // only for 2D
  std::vector < std::pair<double, double>> xi_2d; // used instead of 'xis'

  std::unique_ptr<internals::PrimitiveInternalCoordinates> ic_system_;
  InternalCoordinates::InternalCoordinate* rc_;
  std::size_t rc_index_;

  /**calculates values for xi (reaction coordinate), and E_HL (high level energy) for every structure
  stores them into member variables xis or xi_2d (depending on number of dimenstions) and E_HLs, E_LLs, deltaEs and grad_Es*/
  void calc_energies();
  /**writes outputfile*/
  void write_to_file();
  /**writes splinefile for 1d spline*/
  void write_spline_1d();
  /**writes splinefile for 2d spline*/
  void write_spline_2d();

  double calc_gradient_difference(coords::Representation_3D const& xyz,
                                  coords::Gradients_3D const& grad_hl,
                                  coords::Gradients_3D const& grad_ll);
};