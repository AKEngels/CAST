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

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <omp.h>
#include <limits>
#include <string>
#include <cmath>
#include <tuple>

#include "constants.h"
#include "histogram.h"
#include "Scon/scon_chrono.h"
#include "matop.h"
#include "alignment.h"
#include "kahan_summation.h"

/////////////////
// Some constants
/////////////////
using namespace constants;
using namespace mathFunctions;

// VIA http://www.cplusplus.com/forum/general/209784/
// Computes the distance between two std::vectors
template <typename T>
T	eucledeanDistance(const std::vector<T>& a, std::vector<T> b = std::vector<T>{})
{

  if (b == std::vector<T>{})
  {
    for (unsigned int i = 0u; i < a.size(); i++)
      b.push_back(T(0.));
  }
  std::vector<T>	auxiliary;

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(auxiliary),//
    [](T element1, T element2) {return pow((element1 - element2), 2); });

  return std::sqrt(std::accumulate(auxiliary.begin(), auxiliary.end(), 0.0));

} // end template vectors_distance

// Returns the ardakani Corrected NN distance, see PHYSICAL REVIEW E 83, 051121 (2011)
float_type ardakaniCorrection1D(float_type const& globMin, float_type const& globMax, float_type const& currentPoint, float_type const& NNdistance);

template<typename T>
float_type ardakaniCorrectionGeneralizedEucledeanNorm(std::vector<T> const& globMin, std::vector<T> const& globMax, std::vector<T> const& currentPoint, T const& NNdistance)
{
#ifdef _DEBUG
  if (!(globMin.size() == globMax.size() && globMax.size() == currentPoint.size()))
  {
    throw std::runtime_error("Size Mismatch in N Dimensional Eucledean Ardakani Correction. Aborting.");
  }
#endif
  std::vector<T> radiiOfHyperEllipsoid(globMin.size());
  for (unsigned int i = 0u; i < radiiOfHyperEllipsoid.size(); i++)
  {
    double min = std::min(currentPoint.at(i) + NNdistance * 0.5, globMax.at(i));
    double max = std::max(currentPoint.at(i) - NNdistance * 0.5, globMin.at(i));
    radiiOfHyperEllipsoid.at(i) = min - max;
  }

  // Getting determinant of Mahalanobis distance for calculation of
  // Hyperelipsoid volume. Determinant is product of all eigenvalues. 
  // Eigenvalues of Hyperelipsoid quartic matrix are radius^-2
  //T determinant = 1.;
  //for (unsigned int i = 0u; i < radiiOfHyperEllipsoid.size(); i++)
  //{
  //  determinant *= std::pow(radiiOfHyperEllipsoid.at(i), -2);
  //}
  //
  //// For volume of hyperelipsoid https://math.stackexchange.com/questions/332391/volume-of-hyperellipsoid
  //// Is radius of hypersphere volume 1? I think, because we have all scaling the eigenvalues
  //// But I am not sure...
  //const T V_d = 2. / radiiOfHyperEllipsoid.size() * (std::pow(pi, radiiOfHyperEllipsoid.size() / 2.) / tgamma(radiiOfHyperEllipsoid.size() / 2.)) * std::pow(1, radiiOfHyperEllipsoid.size());
  //const T volumeOfHyperellipsoid = V_d * std::sqrt(determinant);
  //const T equivalentRadiusOfHypersphere = std::pow(volumeOfHyperellipsoid * tgamma(radiiOfHyperEllipsoid.size() / 2. + 1) / std::pow(pi, radiiOfHyperEllipsoid.size() / 2.), 1. / radiiOfHyperEllipsoid.size());
  //return equivalentRadiusOfHypersphere;

  T equivalentRadiusOfHypersphere = T(1.);
  for (unsigned int i = 0u; i < radiiOfHyperEllipsoid.size(); i++)
  {
    equivalentRadiusOfHypersphere *= radiiOfHyperEllipsoid.at(i);
  }
  return std::pow(equivalentRadiusOfHypersphere, 1. / double(radiiOfHyperEllipsoid.size()));
}

template<typename T>
float_type ardakaniCorrectionGeneralizedMaximumNorm(std::vector<T> const& globMin, std::vector<T> const& globMax, std::vector<T> const& currentPoint, T const& NNdistance)
{
#ifdef _DEBUG
  if (!(globMin.size() == globMax.size() && globMax.size() == currentPoint.size()))
  {
    throw std::runtime_error("Size Mismatch in N Dimensional Maximum Norm Ardakani Correction. Aborting.");
  }
#endif
  std::vector<T> radiiOfBox(globMin.size());
  for (unsigned int i = 0u; i < radiiOfBox.size(); i++)
  {
    double min = std::min(currentPoint.at(i) + NNdistance * 0.5, globMax.at(i));
    double max = std::max(currentPoint.at(i) - NNdistance * 0.5, globMin.at(i));
    radiiOfBox.at(i) = min - max;
  }

  T volumeOfBox = 1.;
  for (unsigned int i = 0u; i < radiiOfBox.size(); i++)
  {
    volumeOfBox *= radiiOfBox.at(i);
  }

  const T equivalentRadiusOfHyperSquare = std::pow(volumeOfBox, 1. / radiiOfBox.size());
  return equivalentRadiusOfHyperSquare;

}

void scalePCACoordinatesForQuasiHarmonicTreatment(Matrix_Class& modes, float_type const& temperatureInK);

enum kNN_NORM
{
  EUCLEDEAN = 0,
  MAXIMUM
};

enum kNN_FUNCTION
{
  GORIA = 0,
  LOMBARDI,
  HNIZDO
};

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
  float_type knn_distance_eucl_squared(
    Matrix_Class const& input, 
    size_t const& dimension_in, 
    size_t const& k_in,
    std::vector<size_t> const& row_querypts,
    size_t const& col_querypt,
    coords::float_type* buffer = nullptr);


  float_type maximum_norm_knn_distance(
    Matrix_Class const& input, 
    size_t const& dimension_in, 
    size_t const& k_in,
    std::vector<size_t> const& row_querypts,
    size_t const& col_querypt,
    coords::float_type* buffer = nullptr);


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
    float_type knapp_marginal(bool removeDOF = false);

    /**
    * Performs entropy calculation according to Schlitter
    * Quasi-Harmonic-Approximation used.
    * Gives a strict upper limit to the actual entropy.
    * see: (doi:10.1016/0009-2614(93)89366-P)
    *
    */
    float_type schlitter(float_type const temperatureInKelvin = 300.0);

    Matrix_Class const& getCoordsMatrix(void) const
    {
      return this->coordsMatrix;
    }

    void setCoordsMatrix(Matrix_Class const& in)
    {
      this->coordsMatrix = in;
    }

  private:
    // This matrix is massweightend when cartesians are used
    Matrix_Class coordsMatrix;
  };
}



// entropyobj
// Contains a matrix with draws from a distribution and 
// the number of draws(iter), dimensionality(dimension)
// and standard deviation (sigma)
// 
class entropyobj
{
public:
  size_t numberOfDraws, dimension;
  Matrix_Class drawMatrix;
  std::vector<size_t> subDims;

  entropyobj(Matrix_Class const& drawMatrix_, size_t dimensions_, size_t numberOfDraws_)
    : numberOfDraws(numberOfDraws_), dimension(dimensions_)
  {
    this->drawMatrix = drawMatrix_;
    this->dimension = dimensions_;
    this->numberOfDraws = numberOfDraws_;
  }

  entropyobj(entropy::TrajectoryMatrixRepresentation const& traj)
  {
    this->drawMatrix = traj.getCoordsMatrix();
    this->dimension = traj.getCoordsMatrix().rows();
    this->numberOfDraws = traj.getCoordsMatrix().cols();
    transpose(this->drawMatrix);
  }
};

// Calculated entropy object calculates estiamted entropy
class calculatedentropyobj : public entropyobj
{
public:
  size_t kNN;
  double mean, standardDeviation;
  double empiricalNormalDistributionEntropy;
  Matrix_Class pcaModes;

  calculatedentropyobj(size_t k_, entropyobj const& obj) :
    entropyobj(obj),
    kNN(k_),
    mean(std::numeric_limits<double>::quiet_NaN()),
    standardDeviation(std::numeric_limits<double>::quiet_NaN()),
    empiricalNormalDistributionEntropy(std::numeric_limits<double>::quiet_NaN()),
    pcaModes(Matrix_Class(0u, 0u))
    {

  }

  void empiricalGaussianEntropy()
  {
    std::cout << "Commencing empirical gaussian entropy calculation." << std::endl;
    const unsigned int dimensionality = this->subDims != std::vector<size_t>() ? this->subDims.size() : this->dimension;

    std::cout << "Dimensionality: " << dimensionality << std::endl;
    standardDeviation = 0.0;
    mean = 0.0;
    if (this->dimension == 1)
    {
      // Calculate Mean
      for (unsigned int n = 0; n < numberOfDraws; ++n)
        mean += drawMatrix(n, 0);
      mean /= numberOfDraws;

      // Calculate standard Deviation
      for (unsigned int n = 0; n < numberOfDraws; ++n)
        standardDeviation += (drawMatrix(n, 0) - mean) * (drawMatrix(n, 0) - mean);

      standardDeviation /= numberOfDraws;
      standardDeviation = sqrt(standardDeviation);
      std::cout << "Mean of Empirical Gaussian: " << mean << "\n";
      std::cout << "Standard Deviation of Empricial Gaussian: " << standardDeviation << std::endl;

      empiricalNormalDistributionEntropy = log(standardDeviation * sqrt(2. * pi * e));
    }
    else
    {
      Matrix_Class cov_matr, meanPerDim;

      if (this->subDims != std::vector<size_t>())
      {
        cov_matr = Matrix_Class(this->subDims.size(), this->subDims.size(), 0.);
        meanPerDim = Matrix_Class(this->subDims.size(), 1u, std::numeric_limits<double>::quiet_NaN());
      }
      else
      {
        cov_matr = Matrix_Class(this->dimension, this->dimension, 0.);
        meanPerDim = Matrix_Class(this->dimension, 1u, std::numeric_limits<double>::quiet_NaN());
      }

      //const unsigned int dimensionality = std::max(this->dimension, this->subDims.size());
      for (unsigned int dim = 0u; dim < dimensionality; dim++)
      {
        const unsigned int currentDimension = this->subDims.size() > 0 ? this->subDims.at(dim) : dim;
        double meanThisDim = 0u;
        for (unsigned int i = 0u; i < this->numberOfDraws; i++)
        {
          meanThisDim += this->drawMatrix(i, currentDimension);
        }
        meanThisDim /= double(this->numberOfDraws);
        meanPerDim(dim, 0u) = meanThisDim;

        double sum = 0.;
        for (unsigned int i = 0u; i < this->numberOfDraws; i++)
        {
          sum += std::pow(this->drawMatrix(i, currentDimension) - meanThisDim, 2);
        }
        sum /= double(this->numberOfDraws);
        cov_matr(dim, dim) = sum;

        for (unsigned int dim2 = dim + 1; dim2 < dimensionality; dim2++)
        {
          const unsigned int currentDimension2 = this->subDims.size() > 0 ? this->subDims.at(dim2) : dim2;
          if (meanPerDim(dim2, 0u) != meanPerDim(dim2, 0u))
          {
            double mean = 0u;
            for (unsigned int i = 0u; i < this->numberOfDraws; i++)
            {
              mean += this->drawMatrix(i, currentDimension2);
            }
            mean /= double(this->numberOfDraws);
            meanPerDim(dim2, 0u) = mean;
          }

          double sum2 = 0.;
          for (unsigned int i = 0u; i < this->numberOfDraws; i++)
          {
            sum2 += (this->drawMatrix(i, currentDimension2) - meanPerDim(dim2, 0u)) * (this->drawMatrix(i, currentDimension) - meanPerDim(dim, 0u));
          }
          sum2 /= double(this->numberOfDraws);
          cov_matr(dim, dim2) = sum2;
          cov_matr(dim2, dim) = sum2;
        }
      }

      Matrix_Class eigenval, eigenvec;
      //const double determinant = cov_matr.determ();

      const double gaussentropy = 0.5 * log(cov_matr.determ()) + double(dimensionality) / 2. * std::log(2. * ::constants::pi * ::constants::e);
      empiricalNormalDistributionEntropy = gaussentropy;
      std::cout << "Empirical gaussian entropy (statistical): " << gaussentropy << std::endl;
      //Covariance Matrix


      //Calculate PCA Frequencies in quasi-harmonic approximation and Entropy in SHO approximation; provides upper limit of entropy
      Matrix_Class eigenvalues;
      Matrix_Class eigenvectors;
      std::tie(eigenvalues, eigenvectors) = cov_matr.eigensym(true);
      Matrix_Class pca_frequencies(eigenvalues.rows(), 1u);
      Matrix_Class alpha_i(pca_frequencies.rows(), 1u);
      Matrix_Class quantum_entropy(pca_frequencies.rows(), 1u);
      float_type entropy_sho = 0;
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "Gaussian Entropy in qQH-approximation is estimated from eigenvalues of draw matrix.\n";
        std::cout << "NOTICE: The draw-matrix will be transformed to obtain thermodynamically valid units." << std::endl;
      }
      for (std::size_t i = 0; i < eigenvalues.rows(); i++)
      {
        if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
        {
          // Thermodynamic transformation is applied to eigenvalues here
          //pca_frequencies(i, 0u) = sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp / eigenvalues(i, 0u));
          alpha_i(i, 0u) = 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp) * sqrt(eigenvalues(i, 0u)));
          quantum_entropy(i, 0u) = ((alpha_i(i, 0u) / (exp(alpha_i(i, 0u)) - 1)) - log(1 - exp(-1 * alpha_i(i, 0u)))) * 1.380648813 * 6.02214129 * 0.239005736;
          entropy_sho += quantum_entropy(i, 0u);
          if (Config::get().general.verbosity >= 4)
          {
            std::cout << "MODE " << i << " - Entropy: " << quantum_entropy(i, 0u) << " cal / (mol * K)" << std::endl;
          }
        }
      }
      std::cout << "Entropy in qQH-approximation from PCA-Modes: " << entropy_sho << " cal / (mol * K)" << std::endl;
    }
  }

  void histogramProbabilityDensity(size_t numberOFBins, std::string filename, std::vector<size_t> dimensionsToBeUsed = std::vector<size_t>())
  {
    using namespace histo;
    std::cout << "Writing histogrammed probability density to file " << filename << "." << std::endl;

    std::vector<size_t> dimensionsToBeUsedInHistogramming = dimensionsToBeUsed;
    if (dimensionsToBeUsedInHistogramming == std::vector<size_t>())
    {
      for (unsigned int i = 0u; i < this->dimension; i++)
        dimensionsToBeUsedInHistogramming.push_back(i);
    }
    DimensionalHistogram<float_type> *histograms_p;
    const size_t histogramBins = numberOFBins;

    histograms_p = new DimensionalHistogram<float_type>((size_t)dimensionsToBeUsedInHistogramming.size(), histogramBins);

    //Filling the histogram
    for (size_t j = 0u; j < this->drawMatrix.rows(); j++)
    {
      std::vector<float_type> temp(dimensionsToBeUsedInHistogramming.size());
      for (size_t i = 0u; i < dimensionsToBeUsedInHistogramming.size(); i++)
      {
        temp[i] = drawMatrix(j, dimensionsToBeUsedInHistogramming[i]);
      }
      histograms_p->add_value(temp);
    }

    histograms_p->distribute();
    histograms_p->writeProbabilityDensity("entropytrails_" + filename);
    histograms_p->writeAuxilaryData("auxfile_entropytrails_" + filename);
    delete histograms_p;
  }

  // Full dimensional nearest neighbor entropy computation (without MI-Expansion)
  double calculateNN(const kNN_NORM norm, bool const& ardakaniCorrection, const kNN_FUNCTION func = kNN_FUNCTION::HNIZDO)
  {
    std::cout << "Commencing NNEntropy calculation." << std::endl;

    const unsigned int dimensionality = this->subDims != std::vector<size_t>() ? this->subDims.size() : this->dimension;

    std::cout << "Dimensionality: " << dimensionality << std::endl;
    Matrix_Class dimPurgedDrawMatrix = Matrix_Class(static_cast<uint_type>(this->drawMatrix.rows()), static_cast<uint_type>(dimensionality));
    if (this->subDims != std::vector<size_t>())
    {
      for (size_t i = 0u; i < this->subDims.size(); i++)
      {
        for (size_t j = 0u; j < numberOfDraws; j++)
        {
          dimPurgedDrawMatrix(j, i) = drawMatrix(j, this->subDims.at(i));
        }
      }
    }
    else
    {
      dimPurgedDrawMatrix = drawMatrix;
    }

    //Neccessarry
    transpose(dimPurgedDrawMatrix);
    //transpose(drawMatrix);

    Matrix_Class copytemp = dimPurgedDrawMatrix;
    Matrix_Class eucl_kNN_distances(1u, numberOfDraws, 0.);
    Matrix_Class maxnorm_kNN_distances(1u, numberOfDraws, 0.);
    Matrix_Class eucl_kNN_distances_ardakani_corrected(1u, numberOfDraws, 0.);
    Matrix_Class maxnorm_kNN_distances_ardakani_corrected(1u, numberOfDraws, 0.);

    scon::chrono::high_resolution_timer timer;
    //std::function<std::vector<double>(std::vector<double> const& x)> PDFtemporary = this->probdens.function();
    std::vector<float_type> ardakaniCorrection_minimumValueInDataset(dimensionality, std::numeric_limits<float_type>::max());
    std::vector<float_type> ardakaniCorrection_maximumValueInDataset(dimensionality, -std::numeric_limits<float_type>::max());
    if (ardakaniCorrection)
    {
      for (size_t j = 0; j < drawMatrix.cols(); j++)
      {
        for (unsigned int i = 0u; i < dimensionality; i++)
        {
          if (ardakaniCorrection_minimumValueInDataset.at(i) > drawMatrix(i, j))
            ardakaniCorrection_minimumValueInDataset.at(i) = drawMatrix(i, j);

          if (ardakaniCorrection_maximumValueInDataset.at(i) < drawMatrix(i, j))
            ardakaniCorrection_maximumValueInDataset.at(i) = drawMatrix(i, j);
        }
      }
    }

#ifdef _OPENMP
#pragma omp parallel firstprivate(copytemp, ardakaniCorrection_minimumValueInDataset, ardakaniCorrection_maximumValueInDataset ) \
    shared(eucl_kNN_distances,maxnorm_kNN_distances, eucl_kNN_distances_ardakani_corrected, maxnorm_kNN_distances_ardakani_corrected)
    {
#endif
      float_type* buffer = new float_type[kNN];
#ifdef _OPENMP
      auto const n_omp = static_cast<std::ptrdiff_t>(numberOfDraws);

#pragma omp for
      for (std::ptrdiff_t i = 0; i < n_omp; ++i)
#else
      for (size_t i = 0u; i < numberOfDraws; i++)
#endif
      {
        std::vector<size_t> rowQueryPts;
        for (unsigned int currentDim = 0u; currentDim < dimensionality; currentDim++)
        {
          rowQueryPts.push_back(currentDim);
        }

        std::vector<double> current;
        if (ardakaniCorrection)
        {
          for (unsigned int j = 0u; j < dimensionality; j++)
            current.push_back(copytemp(j, i));
        }

        if (norm == kNN_NORM::EUCLEDEAN)
        {
          const float_type holdNNdistanceEucl = sqrt(entropy::knn_distance_eucl_squared(copytemp, dimensionality, kNN, rowQueryPts, i, buffer));
          eucl_kNN_distances(0, i) = holdNNdistanceEucl;

          if (ardakaniCorrection)
          {
            eucl_kNN_distances_ardakani_corrected(0, i) = ardakaniCorrectionGeneralizedEucledeanNorm(ardakaniCorrection_minimumValueInDataset,
              ardakaniCorrection_maximumValueInDataset, current, holdNNdistanceEucl);
          }
        }
        else // norm == kNN_NORM::MAX
        {
          const float_type holdNNdistanceMax = entropy::maximum_norm_knn_distance(copytemp, dimensionality, kNN, rowQueryPts, i, buffer);
          maxnorm_kNN_distances(0, i) = holdNNdistanceMax;

          if (ardakaniCorrection)
          {
            maxnorm_kNN_distances_ardakani_corrected(0, i) = ardakaniCorrectionGeneralizedMaximumNorm(ardakaniCorrection_minimumValueInDataset,
              ardakaniCorrection_maximumValueInDataset, current, holdNNdistanceMax);
          }
        }

      }
      delete[] buffer;
#ifdef _OPENMP
    }
#endif
    double returnValue = std::numeric_limits<double>::quiet_NaN();

    // Eucledean ArdakaniSum
    if (norm == kNN_NORM::EUCLEDEAN && ardakaniCorrection)
    {
      KahanAccumulation<double> kahan_acc_eucl_ardakani_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_eucl_ardakani_sum = KahanSum(kahan_acc_eucl_ardakani_sum, log(eucl_kNN_distances_ardakani_corrected(0, i)));

      double ardakaniSum = kahan_acc_eucl_ardakani_sum.sum / double(numberOfDraws);
      ardakaniSum *= double(dimensionality);
      ardakaniSum += log(pow(pi, double(dimensionality) / 2.) / (tgamma(0.5 * dimensionality + 1)));

      ardakaniSum -= digammal(double(kNN));

      double ardakaniSum_lombardi = ardakaniSum + digammal(double(numberOfDraws));
      double ardakaniSum_goria = ardakaniSum + log(double(numberOfDraws - 1.));
      double ardakaniSum_hnizdo = ardakaniSum + log(double(numberOfDraws));



      if (func == kNN_FUNCTION::LOMBARDI)
        returnValue = ardakaniSum_lombardi;
      else if (func == kNN_FUNCTION::GORIA)
        returnValue = ardakaniSum_goria;
      else if (func == kNN_FUNCTION::HNIZDO)
        returnValue = ardakaniSum_hnizdo;
      else
        throw std::runtime_error("Critical Error in NN Entropy.");
    }
    else if (norm == kNN_NORM::MAXIMUM && ardakaniCorrection)
    {
      // Maximum Norm ArdakaniSum
      KahanAccumulation<double> kahan_acc_max_ardakani_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_max_ardakani_sum = KahanSum(kahan_acc_max_ardakani_sum, log(maxnorm_kNN_distances_ardakani_corrected(0, i)));

      double maxArdakaniEntropy = kahan_acc_max_ardakani_sum.sum / double(numberOfDraws);
      maxArdakaniEntropy *= double(dimensionality);
      maxArdakaniEntropy += log(pow(2., dimensionality));




      maxArdakaniEntropy -= digammal(double(kNN));

      double ardakaniSum_lombardi = maxArdakaniEntropy + digammal(double(numberOfDraws));
      double ardakaniSum_goria = maxArdakaniEntropy + log(double(numberOfDraws - 1.));
      double ardakaniSum_hnizdo = maxArdakaniEntropy + log(double(numberOfDraws));


      if (func == kNN_FUNCTION::LOMBARDI)
        returnValue = ardakaniSum_lombardi;
      else if (func == kNN_FUNCTION::GORIA)
        returnValue = ardakaniSum_goria;
      else if (func == kNN_FUNCTION::HNIZDO)
        returnValue = ardakaniSum_hnizdo;
      else
        throw std::runtime_error("Critical Error in NN Entropy.");
    }
    else if (norm == kNN_NORM::EUCLEDEAN && !ardakaniCorrection)
    {
      // ENTROPY according to Hnzido
      KahanAccumulation<double> kahan_acc_eucl_sum;

      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_eucl_sum = KahanSum(kahan_acc_eucl_sum, log(eucl_kNN_distances(0, i)));

      double sum = kahan_acc_eucl_sum.sum;
      sum /= double(numberOfDraws);
      sum *= double(dimensionality);

      sum = sum + log(pow(pi, double(dimensionality) / 2.) / (tgamma(0.5 * dimensionality + 1)));

      sum -= digammal(double(kNN));

      double sum_lombardi = sum + digammal(double(numberOfDraws));
      double sum_goria = sum + log(double(numberOfDraws - 1.));
      double sum_hnizdo = sum + log(double(numberOfDraws));


      if (func == kNN_FUNCTION::LOMBARDI)
        returnValue = sum_lombardi;
      else if (func == kNN_FUNCTION::GORIA)
        returnValue = sum_goria;
      else if (func == kNN_FUNCTION::HNIZDO)
        returnValue = sum_hnizdo;
      else
        throw std::runtime_error("Critical Error in NN Entropy.");

    }
    else if (norm == kNN_NORM::MAXIMUM && !ardakaniCorrection)
    {
      KahanAccumulation<double> kahan_acc_max_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_max_sum = KahanSum(kahan_acc_max_sum, log(maxnorm_kNN_distances(0, i)));
      double maxNormSum = kahan_acc_max_sum.sum;

      //
      maxNormSum = maxNormSum / double(numberOfDraws);
      maxNormSum *= double(dimensionality);
      maxNormSum += log(pow(2., dimensionality));

      maxNormSum -= digammal(double(kNN));

      double sum_lombardi = maxNormSum + digammal(double(numberOfDraws));
      double sum_goria = maxNormSum + log(double(numberOfDraws - 1.));
      double sum_hnizdo = maxNormSum + log(double(numberOfDraws));


      if (func == kNN_FUNCTION::LOMBARDI)
        returnValue = sum_lombardi;
      else if (func == kNN_FUNCTION::GORIA)
        returnValue = sum_goria;
      else if (func == kNN_FUNCTION::HNIZDO)
        returnValue = sum_hnizdo;
      else
        throw std::runtime_error("Critical Error in NN Entropy.");
    }

    if (Config::get().general.verbosity >= 3)
    {
      std::cout << "NN Entropy";
      if (ardakaniCorrection)
        std::cout << " with Ardakani-Correction";
      if (norm == kNN_NORM::EUCLEDEAN)
        std::cout << " with L2 norm";
      else
        std::cout << " with Lmax norm";
      if (func == kNN_FUNCTION::LOMBARDI)
        std::cout << " using Lombardi's function";
      else if (func == kNN_FUNCTION::GORIA)
        std::cout << " using Goria's function";
      else if (func == kNN_FUNCTION::HNIZDO)
        std::cout << " using Hnizdo's function";
      std::cout << ": " << returnValue << "." << std::endl;
      std::cout << "NN Calculation took " << timer << " ." << std::endl;
    }
    //Neccessarry
    //transpose(drawMatrix);
    return returnValue;
  }

  // pca Transformation is applied to the Draw-Matrix and the resulting PCA eigenvalues und eigenvectors are stored.
  double pcaTransformDraws(Matrix_Class & eigenvaluesPCA, Matrix_Class & eigenvectorsPCA, bool removeDOF)
  {
    Matrix_Class input(this->drawMatrix);
    transpose(input);

    Matrix_Class cov_matr = Matrix_Class{ transpose(input) };
    cov_matr = cov_matr - Matrix_Class(input.cols(), input.cols(), 1.) * cov_matr / static_cast<float_type>(input.cols());
    cov_matr = transpose(cov_matr) * cov_matr;
    cov_matr *= (1.f / static_cast<float_type>(input.cols()));
    Matrix_Class eigenvalues;
    Matrix_Class eigenvectors;
    float_type cov_determ = 0.;
    int cov_rank = cov_matr.rank();
    std::tie(eigenvalues, eigenvectors) = cov_matr.eigensym(true);


    //Remove Eigenvalues that should be zero if cov_matr is singular
    if ((cov_rank < (int)eigenvalues.rows()) || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90))
    {
      std::cout << "Notice: covariance matrix is singular, attempting to fix by truncation of Eigenvalues.\n";
      std::cout << "Details: rank of covariance matrix is " << cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
      if (removeDOF)
      {
        size_t temp = std::max(6, int((cov_matr.rows() - cov_rank)));
        eigenvalues.shed_rows(eigenvalues.rows() - temp, eigenvalues.rows() - 1u);
        eigenvectors.shed_cols(eigenvectors.cols() - temp, eigenvectors.cols() - 1u);
      }
      else
      {
        eigenvalues.shed_rows((cov_rank), eigenvalues.rows() - 1u);
        eigenvectors.shed_cols((cov_rank), eigenvectors.cols() - 1u);
      }
    }
    else if (removeDOF)
    {
      eigenvectors.shed_cols(0, 5);
      eigenvalues.shed_rows(0, 5);
    }

    eigenvaluesPCA = eigenvalues;
    eigenvectorsPCA = eigenvectors;

    //Calculate PCA Frequencies in quasi-harmonic approximation and Entropy in SHO approximation; provides upper limit of entropy
    Matrix_Class pca_frequencies(eigenvalues.rows(), 1u);
    Matrix_Class alpha_i(pca_frequencies.rows(), 1u);
    Matrix_Class quantum_entropy(pca_frequencies.rows(), 1u);
    float_type entropy_sho = 0;
    for (std::size_t i = 0; i < eigenvalues.rows(); i++)
    {
      if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
      {
        pca_frequencies(i, 0u) = sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp / eigenvalues(i, 0u));
        alpha_i(i, 0u) = 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * Config::get().entropy.entropy_temp) * sqrt(eigenvalues(i, 0u)));
        quantum_entropy(i, 0u) = ((alpha_i(i, 0u) / (exp(alpha_i(i, 0u)) - 1)) - log(1 - exp(-1 * alpha_i(i, 0u)))) * 1.380648813 * 6.02214129 * 0.239005736;
        entropy_sho += quantum_entropy(i, 0u);
      }
    }
    std::cout << "Entropy in qQH-approximation from PCA-Modes: " << entropy_sho << " cal / (mol * K)" << std::endl;
    const unsigned int dimensionality = this->subDims != std::vector<size_t>() ? this->subDims.size() : this->dimension;
    std::cout << "Dimensionality: " << dimensionality << std::endl;


    //Corrections for anharmonicity and M.I.
    // I. Create PCA-Modes matrix
    Matrix_Class eigenvectors_t(transpose(eigenvectorsPCA));

    Matrix_Class input2(this->drawMatrix);
    transpose(input2);

    this->pcaModes = Matrix_Class(eigenvectors_t * input2);

    return entropy_sho;
  }

  /**
  * Performs entropy calculation according to Knapp et al. with corrections
  * for anharmonicity or Mutual Information.
  * Quasi-Harmonic-Approximation used.
  * Hnizdo's entropy estimator is used to calculate the corrections.
  * Gives a strict upper limit to the actual entropy.
  * see: (Genome Inform. 2007;18:192-205.)
  *
  */
  double numataCorrectionsFromMI(size_t orderOfCorrection, Matrix_Class & eigenvaluesPCA,
    const double temperatureInK, const kNN_NORM norm, const kNN_FUNCTION func, const bool removeNegativeMI = true, const float_type anharmonicityCutoff = 0.007)
  {
    const unsigned int dimensionality = this->subDims != std::vector<size_t>() ? this->subDims.size() : this->dimension;

    std::cout << "Dimensionality: " << dimensionality << std::endl;
    scon::chrono::high_resolution_timer timer;
    Matrix_Class pca_modes = this->pcaModes;

    Matrix_Class entropy_anharmonic(pca_modes.rows(), 1u, 0.);

    Matrix_Class statistical_entropy(pca_modes.rows(), 1u, 0.);
    Matrix_Class classical_entropy(pca_modes.rows(), 1u, 0.);


    // Modify PCA modes as the PCA eigenvalues have been modified. This is not detailed in the original paper
    // but sensible and reasonable to obtain valid values.
    scalePCACoordinatesForQuasiHarmonicTreatment(pca_modes, temperatureInK);

    const Matrix_Class storeDrawMatrix = this->drawMatrix;

    const size_t storeDim = this->dimension;
    this->dimension = pca_modes.rows();

    this->drawMatrix = pca_modes;
    this->calculateNN_MIExpansion(orderOfCorrection, norm, func, false);
    this->drawMatrix = storeDrawMatrix;
    this->dimension = storeDim;

    Matrix_Class pca_frequencies(eigenvaluesPCA.rows(), 1u);
    Matrix_Class alpha_i(pca_frequencies.rows(), 1u);
    Matrix_Class quantum_entropy(pca_frequencies.rows(), 1u);
    float_type entropy_sho = 0;
    for (std::size_t i = 0; i < eigenvaluesPCA.rows(); i++)
    {
      pca_frequencies(i, 0u) = sqrt(1.380648813 * 10e-23 * temperatureInK / eigenvaluesPCA(i, 0u));
      alpha_i(i, 0u) = 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * temperatureInK) * sqrt(eigenvaluesPCA(i, 0u)));
      quantum_entropy(i, 0u) = ((alpha_i(i, 0u) / (exp(alpha_i(i, 0u)) - 1)) - log(1 - exp(-1 * alpha_i(i, 0u)))); // This was wrong (remove this comment only for github commit)
      if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
      {
        entropy_sho += quantum_entropy(i, 0u) * 1.380648813 * 6.02214129 * 0.239005736;
      }
    }

    std::vector<calcBuffer> tempMIs = this->calculatedMIs;

    Matrix_Class entropy_kNN(pca_modes.rows(), 1u, 0.);
    for (auto&& element : tempMIs)
    {
      for (size_t i = 1u; i <= orderOfCorrection; i++)
      {
        if (element.dim == i)
        {
          bool isOneModeNotInClassicalLimit = false;
          for (auto&& rowNr : element.rowIdent)
          {
            if (pca_frequencies(rowNr, 0u) > (temperatureInK * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34)))
            {
              isOneModeNotInClassicalLimit = true;
            }
          }
          if (isOneModeNotInClassicalLimit)
          {
            element.entropyValue = 0.0;
          }

        }

        if (element.dim == 1u)
        {
          entropy_kNN(element.rowIdent.at(0), 0u) = element.entropyValue;
        }
      }
    }


    for (size_t i = 0; i < entropy_kNN.rows(); i++)
    {
      if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
      {
        //These are in units S/k_B (therefore: not multiplied by k_B)
        statistical_entropy(i, 0u) = -1.0 * (log(alpha_i(i, 0u)) -/*this might be plus or minus?!*/ log(sqrt(2. * 3.14159265358979323846 * 2.71828182845904523536)));
        classical_entropy(i, 0u) = -1.0 * (log(alpha_i(i, 0u)) - 1.); // should this be +1??? // The formula written HERE NOW is correct, there is a sign error in the original pape rof Knapp/numata
        entropy_anharmonic(i, 0u) = statistical_entropy(i, 0u) - (-1.0) * entropy_kNN(i, 0u);

        // Debug output for developers
        if (Config::get().general.verbosity >= 4)
        {
          std::cout << "---------------------" << std::endl;
          std::cout << "--- Units: S/k_B ---" << std::endl;
          std::cout << "Mode " << i << ": entropy kNN: " << -1.0*entropy_kNN(i, 0u) << "\n";
          std::cout << "Mode " << i << ": entropy anharmonic correction: " << entropy_anharmonic(i, 0u) << "\n";
          std::cout << "Mode " << i << ": classical entropy: " << classical_entropy(i, 0u) << "\n";
          std::cout << "Mode " << i << ": statistical entropy: " << statistical_entropy(i, 0u) << "\n";
          std::cout << "Mode " << i << ": quantum entropy: " << quantum_entropy(i, 0u) << "\n";
          std::cout << "Mode " << i << ": pca freq: " << pca_frequencies(i, 0u) << "\n";
          std::cout << "Mode " << i << ": alpha (dimensionless, standard deviation): " << alpha_i(i, 0u) << "\n";
          std::cout << "Mode " << i << ": standard deviation in mw-pca-units: " << sqrt(eigenvaluesPCA(i, 0u)) << std::endl;
          std::cout << "---------------------" << std::endl;
        }

        if (pca_frequencies(i, 0u) < (temperatureInK * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34)))
        {
          if (std::abs(entropy_anharmonic(i, 0u) / quantum_entropy(i, 0u)) < anharmonicityCutoff)
          {
            entropy_anharmonic(i, 0u) = 0.0;
            std::cout << "Notice: PCA-Mode " << i << " not corrected for anharmonicity (value too small: " << std::abs(entropy_anharmonic(i, 0u) / quantum_entropy(i, 0u)) << " < " << anharmonicityCutoff << ").\n";
          }
          else
          {
            if (entropy_anharmonic(i, 0u) < 0.0)
            {
              std::cout << "Notice: PCA-Mode " << i << " not corrected for anharmonicity (estimated correction is < 0.0 : " << entropy_anharmonic(i, 0u) << ").\n";
              entropy_anharmonic(i, 0u) = 0.0;
            }
          }
        }
        else
        {
          std::cout << "Notice: PCA-Mode " << i << " not corrected since it is not within the classical limit (PCA-Freq needs to be smaller than ";
          std::cout << (temperatureInK * 1.380648813 * 10e-23 / (1.05457172647 * 10e-34)) << ", but is " << pca_frequencies(i, 0u) << "; equipartition is not a valid assumption).\n";
          entropy_anharmonic(i, 0u) = 0.0;
        }

        // Change dimensionless entropy to cal / K * mol
        entropy_anharmonic(i, 0u) *= 1.380648813 * 6.02214129 * 0.239005736;
      }
    }

    // III. Calculate Difference of Entropies
    double delta_entropy = 0;
    for (size_t i = 0; i < entropy_anharmonic.rows(); i++)
    {
      if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
      {
        delta_entropy += entropy_anharmonic(i, 0u);
      }
      if (i == entropy_anharmonic.rows() - 1u)
        std::cout << "Correction for entropy (anharmonicity, order 1): " << delta_entropy << " cal / (mol * K)\n";
    }

    float_type higher_order_entropy = 0.;
    size_t countNegativeSecondOrderMIs = 0u;
    float_type sumOfNegativeSecondOrderMIs = 0.;
    for (size_t i = 2u; i <= orderOfCorrection; i++)
    {
      for (auto&& element : tempMIs)
      {
        if (element.dim == i && i != 1u)
        {
          if (element.dim == 2 && element.entropyValue < 0.)
          {
            if (!removeNegativeMI)
            {
              higher_order_entropy += element.entropyValue * 1.380648813 * 6.02214129 * 0.239005736;
            }
            countNegativeSecondOrderMIs++;
            if (Config::get().general.verbosity >= 4)
            {
              std::cout << "Notice: Negative 2nd order MI for modes " << element.rowIdent.at(0u) << " and " << element.rowIdent.at(1u);
              std::cout << ": " << element.entropyValue * 1.380648813 * 6.02214129 * 0.239005736 << " cal / (mol * K)\n";
            }
            sumOfNegativeSecondOrderMIs += element.entropyValue * 1.380648813 * 6.02214129 * 0.239005736;
          }
          else
          {
            higher_order_entropy += element.entropyValue * 1.380648813 * 6.02214129 * 0.239005736;
          }
        }
      }
      std::cout << "Correction for higher order entropy (mutual information, up to order " << i << "): " << higher_order_entropy << " cal / (mol * K)" << std::endl;
      if (i == 2u)
      {
        std::cout << "Counted " << countNegativeSecondOrderMIs << " negative second order MI terms with a summed value of " << sumOfNegativeSecondOrderMIs << " cal / (mol * K)" << std::endl;
        std::cout << "--- This is not good and should not happen ---\n";
        if (removeNegativeMI)
        {
          std::cout << "--- Negative 2nd order MI terms are set to zero and thus ignored ---\n";
        }
      }
    }

    std::cout << "Correction for entropy: " << delta_entropy + higher_order_entropy << " cal / (mol * K)\n";
    std::cout << "Entropy after correction: " << entropy_sho - delta_entropy - higher_order_entropy << " cal / (mol * K)" << std::endl;

    this->drawMatrix = storeDrawMatrix;
    this->dimension = storeDim;


    std::cout << "NN Calculation took " << timer << " ." << std::endl;

    return entropy_sho - delta_entropy;
  }

  // Calculates the full dimensional Mutual Information Expansion of the entropy up to order "N"
  double calculateNN_MIExpansion(const size_t order_N, const kNN_NORM norm,
    const kNN_FUNCTION func, bool const& ardakaniCorrection)
  {
    scon::chrono::high_resolution_timer timer;

    double miEntropy = 0.;
    std::vector<calcBuffer> buffer;

    // Permutate all var combinatins of size "dim".
    // Calculate jointEntropies of each combination
    // Calculate HMI from joint Entropies up to order dim
    // Sum HMIs accordingly

    std::vector<size_t> rowPts;
    if (this->subDims != std::vector<size_t>())
      rowPts = this->subDims;
    else
    {
      for (size_t i = 0u; i < this->dimension; i++)
        rowPts.push_back(i);
    }
    //transpose(drawMatrix);
    helperFKT(norm, func, ardakaniCorrection, order_N, buffer, rowPts);
    transpose(drawMatrix);
    //Insert unique
    for (auto&& element : buffer)
    {
      if (std::find(this->calculatedMIs.begin(), this->calculatedMIs.end(), element) == this->calculatedMIs.end())
      {
        this->calculatedMIs.push_back(element);
      }
    }


    if (Config::get().general.verbosity >= 4)
    {
      std::cout << "--- MI-EXPANSION INFO ---" << "\n";
      std::cout << "--- Order: << " << order_N << "---" << "\n";
    }
    for (size_t i = 0u; i < buffer.size(); i++)
    {
      if (Config::get().general.verbosity >= 4)
      {
        std::cout << "Contribution " << i << ":  Dimensionality " << buffer.at(i).dim << "\n";
        for (size_t j = 0u; j < buffer.at(i).dim; j++)
        {
          std::cout << "Mode " << buffer.at(i).rowIdent.at(j) << "\n";
        }
        std::cout << "Value: " << std::pow(-1., buffer.at(i).dim + 1.) * buffer.at(i).entropyValue << std::endl;
      }
      miEntropy += std::pow(-1., buffer.at(i).dim + 1.) * buffer.at(i).entropyValue;
    }
    if (Config::get().general.verbosity >= 4)
    {
      std::cout << "--- MI-EXPANSION INFO DONE ---" << "\n";
    }

    storeMI mi;
    mi.order = order_N;
    mi.value = miEntropy;
    this->calculatedMIEs.push_back(mi);

    if (Config::get().general.verbosity >= 3)
    {
      std::cout << "NN Calculation took " << timer << " ." << std::endl;
      std::cout << "NN Value: " << miEntropy << " ." << std::endl;

    }

    return miEntropy;

  }

private:

  struct storeMI
  {
    size_t order = 0u;
    double value = 0u;
  };

  std::vector<storeMI> calculatedMIEs;

  struct calcBuffer
  {
    size_t dim = 0u;
    std::vector<size_t> rowIdent;
    double entropyValue = 0u;

    inline bool operator == (const calcBuffer &b) const
    {
      return this->dim == b.dim && this->rowIdent == b.rowIdent && this->entropyValue == b.entropyValue;
    }

    inline bool operator != (const calcBuffer &b) const
    {
      return !((*this) == b);
    }
  };

  std::vector<calcBuffer> calculatedMIs;

  // Calculates the NN entropy for the specified row indicices
  double calculateNNsubentropy(const kNN_NORM norm, const kNN_FUNCTION func, bool const& ardakaniCorrection, std::vector<size_t> const& rowIndices)
  {
    //transpose(drawMatrix);

    Matrix_Class copytemp = drawMatrix;
    Matrix_Class eucl_kNN_distances(1u, numberOfDraws, 0.);
    Matrix_Class maxnorm_kNN_distances(1u, numberOfDraws, 0.);
    Matrix_Class eucl_kNN_distances_ardakani_corrected(1u, numberOfDraws, 0.);
    Matrix_Class maxnorm_kNN_distances_ardakani_corrected(1u, numberOfDraws, 0.);

    std::vector<float_type> ardakaniCorrection_minimumValueInDataset(rowIndices.size(), std::numeric_limits<float_type>::max());
    std::vector<float_type> ardakaniCorrection_maximumValueInDataset(rowIndices.size(), -std::numeric_limits<float_type>::max());

    if (ardakaniCorrection)
    {
      for (size_t j = 0; j < drawMatrix.cols(); j++)
      {
        for (unsigned int i = 0u; i < rowIndices.size(); i++)
        {
          if (ardakaniCorrection_minimumValueInDataset.at(i) > drawMatrix(rowIndices.at(i), j))
            ardakaniCorrection_minimumValueInDataset.at(i) = drawMatrix(rowIndices.at(i), j);

          if (ardakaniCorrection_maximumValueInDataset.at(i) < drawMatrix(rowIndices.at(i), j))
            ardakaniCorrection_maximumValueInDataset.at(i) = drawMatrix(rowIndices.at(i), j);
        }
      }
    }

#ifdef _OPENMP
#pragma omp parallel firstprivate(copytemp, ardakaniCorrection_minimumValueInDataset, ardakaniCorrection_maximumValueInDataset ) \
    shared(eucl_kNN_distances,maxnorm_kNN_distances, eucl_kNN_distances_ardakani_corrected, maxnorm_kNN_distances_ardakani_corrected)
    {
#endif
      float_type* buffer = new float_type[kNN];
#ifdef _OPENMP
      auto const n_omp = static_cast<std::ptrdiff_t>(numberOfDraws);

#pragma omp for
      for (std::ptrdiff_t i = 0; i < n_omp; ++i)
#else
      for (size_t i = 0u; i < numberOfDraws; i++)
#endif
      {
        std::vector<size_t> rowQueryPts;
        for (unsigned int currentDim = 0u; currentDim < rowIndices.size(); currentDim++)
        {
          rowQueryPts.push_back(rowIndices.at(currentDim));
        }

        std::vector<double> current;
        if (ardakaniCorrection)
        {
          for (unsigned int currentDim = 0u; currentDim < rowIndices.size(); currentDim++)
            current.push_back(copytemp(rowIndices.at(currentDim), i));
        }

        if (norm == kNN_NORM::EUCLEDEAN)
        {
          const float_type holdNNdistanceEucl = sqrt(entropy::knn_distance_eucl_squared(copytemp, rowIndices.size(), kNN, rowQueryPts, i, buffer));
          eucl_kNN_distances(0, i) = holdNNdistanceEucl;

          if (ardakaniCorrection)
          {
            eucl_kNN_distances_ardakani_corrected(0, i) = ardakaniCorrectionGeneralizedEucledeanNorm(ardakaniCorrection_minimumValueInDataset,
              ardakaniCorrection_maximumValueInDataset, current, holdNNdistanceEucl);
          }
        }
        else if (norm == kNN_NORM::MAXIMUM)
        {
          const float_type holdNNdistanceMax = entropy::maximum_norm_knn_distance(copytemp, rowIndices.size(), kNN, rowQueryPts, i, buffer);
          maxnorm_kNN_distances(0, i) = holdNNdistanceMax;

          if (ardakaniCorrection)
          {
            maxnorm_kNN_distances_ardakani_corrected(0, i) = ardakaniCorrectionGeneralizedMaximumNorm(ardakaniCorrection_minimumValueInDataset,
              ardakaniCorrection_maximumValueInDataset, current, holdNNdistanceMax);
          }
        }
        else
        {
          throw std::runtime_error("Unknown norm in kNN calculation. Exiting.\n");
        }


      }
#ifdef _OPENMP
    }
#endif

    //transpose(drawMatrix);

    // Eucledean ArdakaniSum
    if (norm == kNN_NORM::EUCLEDEAN && ardakaniCorrection)
    {
      KahanAccumulation<double> kahan_acc_eucl_ardakani_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_eucl_ardakani_sum = KahanSum(kahan_acc_eucl_ardakani_sum, log(eucl_kNN_distances_ardakani_corrected(0, i)));

      double ardakaniSum = kahan_acc_eucl_ardakani_sum.sum / double(numberOfDraws);
      ardakaniSum *= double(rowIndices.size());
      ardakaniSum += log(pow(pi, double(rowIndices.size()) / 2.) / (tgamma(0.5 * rowIndices.size() + 1)));

      if (func == kNN_FUNCTION::LOMBARDI)
        ardakaniSum += digammal(double(numberOfDraws));
      else if (func == kNN_FUNCTION::GORIA)
        ardakaniSum += log(double(numberOfDraws - 1.));
      else if (func == kNN_FUNCTION::HNIZDO)
        ardakaniSum += log(double(numberOfDraws));
      else
        throw std::runtime_error("Critical Error in NN Entropy.");

      ardakaniSum -= digammal(double(kNN));

      return ardakaniSum;
    }
    else if (norm == kNN_NORM::MAXIMUM && ardakaniCorrection)
    {
      // Maximum Norm ArdakaniSum
      KahanAccumulation<double> kahan_acc_max_ardakani_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_max_ardakani_sum = KahanSum(kahan_acc_max_ardakani_sum, log(maxnorm_kNN_distances_ardakani_corrected(0, i)));

      double maxArdakaniEntropy = kahan_acc_max_ardakani_sum.sum / double(numberOfDraws);
      maxArdakaniEntropy *= double(rowIndices.size());
      maxArdakaniEntropy += log(pow(2., rowIndices.size()));

      if (func == kNN_FUNCTION::LOMBARDI)
        maxArdakaniEntropy += digammal(double(numberOfDraws));
      else if (func == kNN_FUNCTION::GORIA)
        maxArdakaniEntropy += log(double(numberOfDraws - 1.));
      else if (func == kNN_FUNCTION::HNIZDO)
        maxArdakaniEntropy += log(double(numberOfDraws));
      else
        throw std::runtime_error("Critical Error in NN Entropy.");

      maxArdakaniEntropy -= digammal(double(kNN));

      return maxArdakaniEntropy;
    }
    else if (norm == kNN_NORM::EUCLEDEAN && !ardakaniCorrection)
    {
      // ENTROPY according to Hnzido
      KahanAccumulation<double> kahan_acc_eucl_sum;

      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_eucl_sum = KahanSum(kahan_acc_eucl_sum, log(eucl_kNN_distances(0, i)));

      double hnizdoSum = kahan_acc_eucl_sum.sum;
      hnizdoSum /= double(numberOfDraws);
      hnizdoSum *= double(rowIndices.size());

      // Einschub
      // Entropy according to Lombardi (traditional)
      double lombardi_without_correction = hnizdoSum + log(pow(pi, double(rowIndices.size()) / 2.) / (tgamma(0.5 * double(rowIndices.size()) + 1.)));

      if (func == kNN_FUNCTION::LOMBARDI)
        lombardi_without_correction += digammal(double(numberOfDraws));
      else if (func == kNN_FUNCTION::GORIA)
        lombardi_without_correction += log(double(numberOfDraws - 1.));
      else if (func == kNN_FUNCTION::HNIZDO)
        lombardi_without_correction += log(double(numberOfDraws));
      else
        throw std::runtime_error("Critical Error in NN Entropy.");

      lombardi_without_correction -= digammal(double(kNN));

      return lombardi_without_correction; // Lombardi Entropy
    }
    else if (norm == kNN_NORM::MAXIMUM && !ardakaniCorrection)
    {
      KahanAccumulation<double> kahan_acc_max_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_max_sum = KahanSum(kahan_acc_max_sum, log(maxnorm_kNN_distances(0, i)));
      double maxNormSum = kahan_acc_max_sum.sum;

      //
      double lobardi_maximum_norm = maxNormSum / double(numberOfDraws);
      lobardi_maximum_norm *= double(rowIndices.size());
      lobardi_maximum_norm += log(pow(2., rowIndices.size()));

      if (func == kNN_FUNCTION::LOMBARDI)
        lobardi_maximum_norm += digammal(double(numberOfDraws));
      else if (func == kNN_FUNCTION::GORIA)
        lobardi_maximum_norm += log(double(numberOfDraws - 1.));
      else if (func == kNN_FUNCTION::HNIZDO)
        lobardi_maximum_norm += log(double(numberOfDraws));
      else
        throw std::runtime_error("Critical Error in NN Entropy.");

      lobardi_maximum_norm -= digammal(double(kNN));

      return lobardi_maximum_norm;
    }
    else
    {
      throw std::runtime_error("Something went terribly wrong in kNN calculation. Sorry.");
      return 0.;
    }
  }

  // Called recursivly, calculates joint entropies for all possible permutations of the rowIndices up to the full dimensional one
  // and stores those results in the "result" vecor
  void jointEntropiesOfCertainRows(const kNN_NORM norm, const kNN_FUNCTION func, bool const& ardakaniCorrection, std::vector<calcBuffer>& result, std::vector<size_t> rowIndices,
    std::vector<size_t> curRowIndices = std::vector<size_t>{}, const size_t curDim = 0u, size_t startIter = 0u)
  {
    const size_t maxDim = rowIndices.size();
    if (curDim < maxDim)
    {
      for (size_t i = startIter; i < maxDim; i++)
      {
        if (curDim == 0u)
          curRowIndices = std::vector<size_t>{};

        if (curRowIndices.size() == curDim + 1u)
          curRowIndices.at(curDim) = rowIndices.at(i);
        else
          curRowIndices.push_back(rowIndices.at(i));

        size_t nextIter = i + 1u;

        calcBuffer cBuffer;
        cBuffer.dim = curRowIndices.size();
        cBuffer.rowIdent = curRowIndices;

        cBuffer.entropyValue = calculateNNsubentropy(norm, func, ardakaniCorrection, curRowIndices);

        result.push_back(cBuffer);

        jointEntropiesOfCertainRows(norm, func, ardakaniCorrection, result, rowIndices, curRowIndices, curRowIndices.size(), nextIter);
      }
    }
  }

  // Called recursivly, calculates Higher-Order Mutual Information for all possible permutations of the rowIndices up to the full dimensional one
  // and stores those results in the "result" vecor
  void helperFKT(const kNN_NORM norm, const kNN_FUNCTION func, bool const& ardakaniCorrection, const size_t maxDim, std::vector<calcBuffer>& result, std::vector<size_t> rowIndices,
    std::vector<size_t> curRowIndices = std::vector<size_t>{}, const size_t curDim = 0u, size_t startIter = 0u)
  {
    if (curDim < maxDim)
    {
      for (size_t i = startIter; i < rowIndices.size(); i++)
      {
        if (curDim == 0u)
          curRowIndices = std::vector<size_t>{};

        if (curRowIndices.size() == curDim + 1u)
          curRowIndices.at(curDim) = rowIndices.at(i);
        else
          curRowIndices.push_back(rowIndices.at(i));

        size_t nextIter = i + 1u;

        std::vector<calcBuffer> results;
        jointEntropiesOfCertainRows(norm, func, ardakaniCorrection, results, curRowIndices);



        // Debug
        //std::cout << "\n";
        //for (auto&& item : results)
        //{
        //  for (auto&& subitem : item.rowIdent)
        //    std::cout << " " << subitem;
        //  std::cout << "\n";
        //}
        //std::cout << std::endl;

        calcBuffer cBuffer;
        cBuffer.dim = curRowIndices.size();
        cBuffer.rowIdent = curRowIndices;
        cBuffer.entropyValue = 0.;
        // SUMM
        for (size_t j = 0u; j < results.size(); j++)
        {
          cBuffer.entropyValue += std::pow(-1., results.at(j).dim + 1.) * results.at(j).entropyValue;
        }
        result.push_back(cBuffer);

        helperFKT(norm, func, ardakaniCorrection, maxDim, result, rowIndices, curRowIndices, curRowIndices.size(), nextIter);
      }
    }
  }

};

