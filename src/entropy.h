/**
CAST 3
entropy.h
Purpose:
Specific algorithms for calculations of entropy.
Currently stable and tested:
Karplus, Knapp marginal
Currently stable although dubious results:
Hnizdo, Hnizdo marginal, Knapp

@author Dustin Kaiser
@version 1.0
*/

#pragma once
#include "matop.h"
#include "alignment.h"

namespace entropy
{
  /**
  * Outputs the !SQUARED! next-neighbor distance in eucledean space
  * of the k-nearest neighbor to the query-Point. Rows are dimensions
  * of the data points, columns are the actual data points.
  *
  * @param dimension_in Dimensionality of the search (needs to
  * be identical to row_querypts.size().
  * @param k_in The k-th nearest neighbor will be searched.
  * @param row_querypts std::vector of rows (ie dimensions)
  * used in the search (example {1, 4, 5}).
  * @param col_querypt Columns index of the query Points.
  */
  float_type knn_distance(
    Matrix_Class const& input, size_t
    const& dimension_in, size_t const& k_in,
    std::vector<size_t>& row_querypts,
    size_t const& col_querypt,
    coords::float_type* buffer = nullptr);


  float_type maximum_norm_knn_distance(
    Matrix_Class const& input, size_t
    const& dimension_in, size_t const& k_in,
    std::vector<size_t>& row_querypts,
    size_t const& col_querypt,
    coords::float_type* buffer = nullptr);

  /**
  * Outputs the !SQUARED! next-neighbor distance in eucledean space
  * of the k-nearest neighbor to the query-Point. Rows are dimensions
  * of the data points, columns are the actual data points.
  *
  * @param dimension_in Dimensionality of the search (will simply use
  * subsequent rows starting from "row_query_Pt").
  * @param k_in The k-th nearest neighbor will be searched.
  * @param row_querypt Row index of the query Point.
  * @param col_querypt Columns index of the query Point.
  */
  float_type knn_distance(Matrix_Class const& input, size_t const& dimension_in, size_t const& k_in, size_t const& row_querypt, size_t const& col_querypt, coords::float_type* buffer = nullptr);


  /**
  * Class used for Entropy calculations based
  * on MD trajectories
  * For sample use see task "ENTROPY"
  */
  class TrajectoryMatrixRepresentation
  {
  public:
    /**
    * Constructor, subsequently calls (in this order):
    * -> generateCoordinateMatrix
    */
    TrajectoryMatrixRepresentation(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords);

    /**
    * Generates matrix representation of
    * a MD trajectory
    */
    void generateCoordinateMatrix(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords);

    /**
    * Performs entropy calculation according to Hnizdo et al.
    * see: (DOI: 10.1002/jcc.20589)
    * If rows > 10 this might be really computationally expensive
    * Basically, many "next-neighbor"-searches are performed using a
    * primitive brute-force approach.
    *
    * @param k: value of k for k nearest neighboor search
    */
    float_type hnizdo(size_t const k = 5 );

    /**
    * Performs entropy calculation according to Hnizdo et al.
    * Sums marginal (== 1-dimensional) entropy, i.e. calculating
    * values for each "row" independently and then building sum.
    * see: (DOI: 10.1002/jcc.20589)
    * Marginal Entropy is upper limit to actual entropy which is
    * computationally expensive ( "hnizdo()" )
    *
    * @param k: value of k for k nearest neighboor search
    */
    float_type hnizdo_marginal(size_t const k = 5);

    /**
     * Karplus entropy,
     * basically never works (by design, its a really old approach)
     * see: DOI 10.1021/ma50003a019
     */
    float_type karplus();

    /**
    * Performs entropy calculation according to Knapp et al. without corrections
    * for anharmonicity or Mutual Information.
    * Quasi-Harmonic-Approximation used.
    * This method is almost identical to Schlitter's approach.
    * Gives a strict upper limit to the actual entropy.
    * see: (Genome Inform. 2007;18:192-205.)
    *
    */
    float_type knapp_marginal(float_type const temperatureInKelvin = 300.0, bool removeDOF = false);

    /**
    * Performs entropy calculation according to Knapp et al. with corrections
    * for anharmonicity or Mutual Information.
    * Quasi-Harmonic-Approximation used.
    * Hnizdo's entropy estimator is used to calculate the corrections.
    * Gives a strict upper limit to the actual entropy.
    * see: (Genome Inform. 2007;18:192-205.)
    *
    */
    float_type knapp(float_type const temperatureInKelvin = 300.0, size_t const k = 5, bool removeDOF = false);

    /**
    * Performs entropy calculation according to Schlitter
    * Quasi-Harmonic-Approximation used.
    * Gives a strict upper limit to the actual entropy.
    * see: (doi:10.1016/0009-2614(93)89366-P)
    *
    */
    float_type schlitter(float_type const temperatureInKelvin = 300.0);

  private:
    // This matrix is massweightend when cartesians are used
    Matrix_Class coordsMatrix;
  };
}