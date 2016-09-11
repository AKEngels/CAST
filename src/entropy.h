#pragma once
#include "matop.h"
namespace entrop
{
  class TrajectoryMatrixRepresentation
  {
  public:
    TrajectoryMatrixRepresentation(std::unique_ptr<coords::input::format>& ci, coords::Coordinates& coords);

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
    float_type TrajectoryMatrixRepresentation::knapp_marginal(float_type const temperatureInKelvin = 300.0, bool removeDOF = false);

    /**
    * Performs entropy calculation according to Knapp et al. with corrections
    * for anharmonicity or Mutual Information.
    * Quasi-Harmonic-Approximation used.
    * Hnizdo's entropy estimator is used to calculate the corrections.
    * Gives a strict upper limit to the actual entropy.
    * see: (Genome Inform. 2007;18:192-205.)
    *
    */
    float_type TrajectoryMatrixRepresentation::knapp(float_type const temperatureInKelvin = 300.0, size_t const k = 5, bool removeDOF = false);

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