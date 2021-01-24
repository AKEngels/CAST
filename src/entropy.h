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
#include <algorithm>

#include "constants.h"
#include "histogram.h"
#include "Scon/scon_chrono.h"
#include "TrajectoryMatrixClass.h"
#include "alignment.h"
#include "kahan_summation.h"
#include "PCA.h"

/////////////////
// Some constants
/////////////////
using namespace constants;
using namespace mathFunctions;
using float_type = coords::float_type;

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

  Matrix_Class unmassweightedStdDevFromMWPCAeigenvalues(Matrix_Class const& massVector, Matrix_Class const& pcaEigenvalues, Matrix_Class const& pcaEigenvectors, std::vector<size_t> const& subDims = std::vector<size_t>());
  Matrix_Class calculateReducedMassOfPCAModes(Matrix_Class const& massVector, Matrix_Class const& pca_eigenvalues, Matrix_Class const& pca_eigenvectors, std::vector<size_t> const& subDims = std::vector<size_t>());

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
    // TO-DO: Asser numberOfDraws and Dimension are equal to matrix size!
  }

  entropyobj(TrajectoryMatrixRepresentation const& traj)
  {
    this->drawMatrix = traj.getCoordsMatrix();
    this->dimension = traj.getCoordsMatrix().rows();
    this->numberOfDraws = traj.getCoordsMatrix().cols();
    this->subDims = traj.getSubDims();
    transpose(this->drawMatrix);
  }

  std::vector<std::size_t> const& getSubDims() const
  {
    return this->subDims;
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
  Matrix_Class classicalHarmonicShiftingConstants;

  calculatedentropyobj(size_t k_, entropyobj const& obj) :
    entropyobj(obj),
    kNN(k_),
    mean(std::numeric_limits<double>::quiet_NaN()),
    standardDeviation(std::numeric_limits<double>::quiet_NaN()),
    empiricalNormalDistributionEntropy(std::numeric_limits<double>::quiet_NaN()),
    pcaModes(Matrix_Class(0u, 0u))
  {

  }

  double empiricalGaussianEntropy()
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

      const double gaussentropy = 0.5 * log(cov_matr.determ()) + double(dimensionality) / 2. * std::log(2. * ::constants::pi * ::constants::e);
      empiricalNormalDistributionEntropy = gaussentropy;
    }
    std::cout << "Empirical gaussian entropy (statistical): " << empiricalNormalDistributionEntropy << std::endl;
    return empiricalNormalDistributionEntropy;
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
    DimensionalHistogram<float_type>* histograms_p;
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

  double calculateFulldimensionalNNEntropyOfDraws(const kNN_NORM norm = kNN_NORM::EUCLEDEAN, bool const& ardakaniCorrection = false, const kNN_FUNCTION func = kNN_FUNCTION::HNIZDO) const
  {
    return this->calculateNN(this->drawMatrix,norm,ardakaniCorrection,func);
  }

  double calculateFulldimensionalNNEntropyOfPCAModes(const kNN_NORM norm = kNN_NORM::EUCLEDEAN, bool const& ardakaniCorrection = false, const kNN_FUNCTION func = kNN_FUNCTION::HNIZDO) const
  {
    if (this->pcaModes.rows() == 0u && this->pcaModes.cols() == 0u)
    {
      throw std::logic_error("PCA Modes have not been calculated. You need to call 'pcaTransformDraws()' first. Aborting.");
      return -1.;
    }
    else
    {
      return this->calculateNN(transposed(this->pcaModes), norm, ardakaniCorrection, func);
    }
  }

  // pca Transformation is applied to the Draw-Matrix and the resulting PCA eigenvalues und eigenvectors are stored.
  double pcaTransformDraws(Matrix_Class& eigenvaluesPCA, Matrix_Class& eigenvectorsPCA, Matrix_Class massVector, \
    double temperatureInK = 0.0, bool removeDOF = true )
  {
    //
    if (temperatureInK == 0.0 && Config::get().entropy.entropy_temp != 0.0)
    {
      std::cout << "Assuming a temperature of " << Config::get().entropy.entropy_temp << "K for quasiharmonic treatment.\n";
      temperatureInK = Config::get().entropy.entropy_temp;
    }
    //

    std::cout << "Transforming the input coordinates into their PCA modes.\n";
    std::cout << "This directly yields the marginal Quasi - Harmonic - Approx. according to Knapp et. al. without corrections (Genome Inform. 2007; 18:192 - 205)\n";
    std::cout << "Commencing calculation..." << std::endl;
    std::cout << "DEBUG drawMatrix:\n" << this->drawMatrix << std::endl;
    Matrix_Class input(this->drawMatrix);
    Matrix_Class cov_matr = Matrix_Class{ input };
    cov_matr = cov_matr - Matrix_Class(input.rows(), input.rows(), 1.) * cov_matr / static_cast<float_type>(input.rows());
    cov_matr = transposed(cov_matr) * cov_matr;
    cov_matr *= (1.f / static_cast<float_type>(input.rows()));
    Matrix_Class eigenvalues;
    Matrix_Class eigenvectors;
    float_type cov_determ = 0.;
    int cov_rank = cov_matr.rank();
    std::tie(eigenvalues, eigenvectors) = cov_matr.eigensym(true);
    // Checks
    if (massVector.cols() == 1u)
    {
      if (massVector.rows() != eigenvalues.rows())
      {
        throw std::logic_error("Given Mass Vector has wrong dimensionality, aborting.");
        return -1.0;
      }
    }


    if (removeDOF)
    {
      //Remove Eigenvalues that should be zero if cov_matr is singular
      if ((cov_rank < (int)eigenvalues.rows()) || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90))
      {
        std::cout << "Notice: covariance matrix is singular, attempting to fix by truncation of Eigenvalues.\n";
        std::cout << "Details: rank of covariance matrix is " << cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
        size_t temp = std::max(6, int((cov_matr.rows() - cov_rank)));
        eigenvalues.shed_rows(eigenvalues.rows() - temp, eigenvalues.rows() - 1u);
        eigenvectors.shed_cols(eigenvectors.cols() - temp, eigenvectors.cols() - 1u);

      }
      else
      {
        eigenvectors.shed_cols(0, 5);
        eigenvalues.shed_rows(0, 5);
      }
    }



    eigenvaluesPCA = eigenvalues;
    eigenvectorsPCA = eigenvectors;
    Matrix_Class assocRedMasses = entropy::calculateReducedMassOfPCAModes(massVector, eigenvalues, eigenvectors, this->subDims);
    std::cout << "DEBUG: assoc red masses:\n" << assocRedMasses << std::endl;
    Matrix_Class eigenvectors_t(transposed(eigenvectorsPCA));
    Matrix_Class input2(transposed(this->drawMatrix)); // Root mass weighted cartesian coords, most likely...
    //
    this->pcaModes = Matrix_Class(eigenvectors_t * input2);
    const Matrix_Class covarianceMatrixOfPCAModes = this->pcaModes.covarianceMatrix();
    std::cout << "DEBUG covariance matrix of mass weighted pca modes:\n";
    for (std::size_t i = 0; i < covarianceMatrixOfPCAModes.rows(); i++)
      std::cout << covarianceMatrixOfPCAModes(i, i) << "\n";
    //std::cout << covarianceMatrixOfPCAModes << std::endl;
    std::cout << std::endl;
    this->pcaModes = this->unmassweightPCAModes(assocRedMasses,Matrix_Class(eigenvectors_t * input2));
    const Matrix_Class covarianceMatrixOfUnweightedPCAModes = this->pcaModes.covarianceMatrix();
    std::cout << "DEBUG covariance matrix of unweighted pca modes:\n";
    for (std::size_t i = 0; i < covarianceMatrixOfUnweightedPCAModes.rows(); i++)
      std::cout << covarianceMatrixOfUnweightedPCAModes(i, i) << "\n";
    //std::cout << covarianceMatrixOfPCAModes << std::endl;
    std::cout << std::endl;
    //std::cout << "PCA-Vec_t:\n" << eigenvectors_t << std::endl;
    //std::cout << "PCA-Modes (unweighted):\n" << this->pcaModes << std::endl; // Nrows are 3xDOFs, NCloumns are Nframes

    const Matrix_Class covarianceMatrixOfINPUT = input2.covarianceMatrix();
    std::cout << "DEBUG Covariance Matrix of input DrawMatrix:\n";
    for (std::size_t i = 0; i < covarianceMatrixOfINPUT.rows(); i++)
      std::cout << covarianceMatrixOfINPUT(i, i) << "\n";
    

    //Calculate PCA Frequencies in quasi-harmonic approximation and Entropy in SHO approximation; provides upper limit of entropy
    Matrix_Class pca_frequencies(eigenvalues.rows(), 1u);
    Matrix_Class alpha_i(pca_frequencies.rows(), 1u);
    Matrix_Class quantum_entropy(pca_frequencies.rows(), 1u);
    Matrix_Class classical_entropy(pca_frequencies.rows(), 1u);
    Matrix_Class statistical_entropy(pca_frequencies.rows(), 1u);
    //
    // via https://physics.stackexchange.com/questions/401370/normal-modes-how-to-get-reduced-masses-from-displacement-vectors-atomic-masses
    Matrix_Class constant_C(pca_frequencies.rows(), 1u);
    float_type entropy_qho = 0.;
    float_type entropy_cho = 0.;
    
    

    for (std::size_t i = 0; i < eigenvalues.rows(); i++)
    {
      if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
      {
        //Magic number here, something is wrong with a constant :--(
        pca_frequencies(i, 0u) = sqrt(constants::boltzmann_constant_kb_SI_units * temperatureInK / eigenvalues(i, 0u));
        if (massVector.cols() == 1u && massVector.rows() == eigenvalues.rows())
        {
          std::cout << "....................\n";
          std::cout << "Debug: Mode " << i << std::endl;
          std::cout << "Debug: kB SI " << constants::boltzmann_constant_kb_SI_units << std::endl;
          std::cout << "Debug: eigenvalues " << eigenvalues(i, 0u) << std::endl;
          std::cout << "Debug: pca_frequencies " << pca_frequencies(i, 0u) << std::endl;
          std::cout << "Debug: pca_frequencies cm-1 " << pca_frequencies(i, 0u) / constants::speed_of_light_cm_per_s << std::endl;
          // Assoc red mass of each mode via https://physics.stackexchange.com/questions/401370/normal-modes-how-to-get-reduced-masses-from-displacement-vectors-atomic-masses
          double A___normalizationThisEigenvector = 0.;
          double inv_red_mass = 0.0;
          for (std::size_t j = 0; j < eigenvectors.cols(); j++)
          {
            //Each column is one eigenvector
            const double squaredEigenvecValue = eigenvectors(i, j) * eigenvectors(i, j);
            A___normalizationThisEigenvector += squaredEigenvecValue;
            std::cout << "Debug: A___normalizationThisEigenvector " << A___normalizationThisEigenvector << std::endl;
            const double currentMass = massVector(i, 0u);
            std::cout << "Debug: currentMass " << currentMass << std::endl;
            inv_red_mass += A___normalizationThisEigenvector / currentMass;
            std::cout << "Debug: inv_red_mass currently  " << inv_red_mass << std::endl;
          }
          //
          const double red_mass = 1.0/inv_red_mass;
          std::cout << "Debug: red_mass " << red_mass << std::endl;
          //assocRedMasses(i,0u) = red_mass;
          std::cout << "Debug: Sanity check red_mass: " << assocRedMasses(i, 0u) << std::endl;
          const double squaredStdDev = covarianceMatrixOfPCAModes(i,i);
          std::cout << "Debug: squaredStdDev (convoluted with red mass) " << squaredStdDev << std::endl;
          //
          const double stdDev_ofPCAMode_inSIUnits = std::sqrt(squaredStdDev) / std::sqrt(red_mass);
          std::cout << "SDEBUG: sqrt(" << squaredStdDev << ")/sqrt(" << red_mass << ")= " << stdDev_ofPCAMode_inSIUnits << std::endl;
          const double x_0 = stdDev_ofPCAMode_inSIUnits * std::sqrt(2);
          const double x_0_SI = stdDev_ofPCAMode_inSIUnits * std::sqrt(2);
          const double Sspatial = constants::joules2cal * constants::N_avogadro*(-1.0 * constants::boltzmann_constant_kb_SI_units * (std::log(2 / constants::pi) - std::log(x_0_SI)));
          std::cout << "Debug: StdDev in SI units " << stdDev_ofPCAMode_inSIUnits << std::endl;
          const double gaussianSSpatial = std::log(stdDev_ofPCAMode_inSIUnits * std::sqrt(2*constants::pi*std::exp(1.))); //via https://en.wikipedia.org/wiki/Differential_entropy#:~:text=With%20a%20normal%20distribution%2C%20differential,and%20variance%20is%20the%20Gaussian.
          std::cout << "Debug: Entropy of gaussian with this StdDev " << gaussianSSpatial << std::endl;
          std::cout << "Debug: Sspatial " << Sspatial << std::endl;
          std::cout << "Debug: Sspatial in raw units " << -1.0 * (std::log(2 / constants::pi) - std::log(x_0_SI)) << std::endl;
          std::cout << "Debug: x_0_SI " << x_0_SI << std::endl;
          const double C1 = constants::boltzmann_constant_kb_SI_units * temperatureInK / constants::h_bar_SI_units / pca_frequencies(i, 0u) * 2. / constants::pi / x_0_SI;
          std::cout << "Debug: C1 " << C1 << std::endl;
          const double C2 = std::log(C1) + 1.;
          std::cout << "Debug: C2 " << C2 << std::endl;
          const double C3 = constants::joules2cal * constants::N_avogadro * constants::boltzmann_constant_kb_SI_units * C2;
          std::cout << "Debug: C3 " << C3 << std::endl;
          //C=k*(ln((k*temp)/h_red/freq*2/pi/x_0) + 1)
          //C_dash = C * avogadro
          constant_C(i,0u) = C3;
        }
        alpha_i(i, 0u) = constants::h_bar_SI_units / (sqrt(constants::boltzmann_constant_kb_SI_units * temperatureInK) * sqrt(eigenvalues(i, 0u)));
        std::cout << "Debug: alpha_i " << alpha_i(i, 0u) << std::endl;
        const double sanityCheck = constants::h_bar_SI_units * pca_frequencies(i, 0u) / constants::boltzmann_constant_kb_SI_units / temperatureInK;
        std::cout << "Debug: sanitycheck " << sanityCheck << std::endl;
        //These are in units S/k_B (therefore: not multiplied by k_B)
        quantum_entropy(i, 0u) = ((alpha_i(i, 0u) / (exp(alpha_i(i, 0u)) - 1)) - log(1 - exp(-1 * alpha_i(i, 0u)))) * 1.380648813 * 6.02214129 * 0.239005736;
        std::cout << "Debug: quantum_entropy " << quantum_entropy(i, 0u) << std::endl;
        statistical_entropy(i, 0u) = -1.0 * constants::N_avogadro * constants::boltzmann_constant_kb_SI_units * constants::joules2cal * (log(alpha_i(i, 0u)) -/*this might be plus or minus?!*/ log(sqrt(2. * constants::pi * 2.71828182845904523536)));
        std::cout << "Debug: statistical_entropy " << statistical_entropy(i, 0u) << std::endl;
        classical_entropy(i, 0u) = -1.0 * constants::N_avogadro * constants::boltzmann_constant_kb_SI_units * constants::joules2cal * (log(alpha_i(i, 0u)) - 1.); // should this be +1??? // The formula written HERE NOW is correct, there is a sign error in the original pape rof Knapp/numata
        std::cout << "Debug: classical_entropy " << classical_entropy(i, 0u) << std::endl;
        std::cout << "Debug: classical_entropy in raw units " << classical_entropy(i, 0u) / (constants::N_avogadro * constants::boltzmann_constant_kb_SI_units * constants::joules2cal )<< std::endl;
        //
        //
        if (!std::isnan(quantum_entropy(i, 0u)))
          entropy_qho += quantum_entropy(i, 0u);
        if (!std::isnan(classical_entropy(i, 0u)))
          entropy_cho += classical_entropy(i, 0u);
      }
    }
    if (Config::get().general.verbosity >= 3u)
    {
      std::cout << "----------\nPrinting PCA_Mode Frequencies and quantum QHA-Entropies:\n";
      std::cout << std::setw(20) << "Mode #" << std::setw(20) << "assoc. red. mass" << std::setw(20) << "cm^-1" << std::setw(20) << "entropy" << std::setw(20) << "% contrib." << std::setw(20) << "C[SI]" << "\n";
      //
      for (std::size_t i = 0; i < eigenvalues.rows(); i++)
      {
        if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
        {
          std::cout << std::setw(20);
          std::cout << i + 1 << std::setw(20) << assocRedMasses(i, 0u) << std::setw(20) << pca_frequencies(i, 0u) / constants::speed_of_light_cm_per_s << std::setw(20) << quantum_entropy(i, 0u);
          std::cout << std::setw(20) << std::to_string(std::round(quantum_entropy(i, 0u) / entropy_qho*10000)/100) + " %" ;
          std::cout << std::setw(20) << constant_C(i, 0u);
          std::cout << "\n";
        }
      }
      std::cout << "----------\n";
    }
    std::cout << "Entropy in quantum QH-approximation from PCA-Modes: " << entropy_qho << " cal / (mol * K)" << std::endl;
    std::cout << "Entropy in classical QH-approximation from PCA-Modes: " << entropy_cho << " cal / (mol * K)" << std::endl;
    

    this->classicalHarmonicShiftingConstants = constant_C;
    return entropy_qho;
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
  double numataCorrectionsFromMI(size_t orderOfCorrection, Matrix_Class const& eigenvaluesPCA, Matrix_Class const& stdDevPCAModes,
    const double temperatureInK, const kNN_NORM norm, const kNN_FUNCTION func, const bool removeNegativeMI = true, const float_type anharmonicityCutoff = 0.007)
  {
    std::cout << "\nCommencing entropy calculation:\nHybrid-Approach according to Knapp et. al. with 1st/2nd order MI corrections (Genome Inform. 2007;18:192-205.)" << std::endl;

    //std::cout << "Dimensionality: " << dimensionality << std::endl;
    scon::chrono::high_resolution_timer timer;
    //unmassweightPCAModes
    Matrix_Class const& pca_modes = this->pcaModes;
    //
    Matrix_Class entropy_anharmonic(pca_modes.rows(), 1u, 0.);

    Matrix_Class statistical_entropy(pca_modes.rows(), 1u, 0.);
    Matrix_Class classical_entropy(pca_modes.rows(), 1u, 0.);


    // Modify PCA modes as the PCA eigenvalues have been modified. This is not detailed in the original paper
    // but sensible and reasonable to obtain valid values.
    //scalePCACoordinatesForQuasiHarmonicTreatment(pca_modes, temperatureInK);

    const Matrix_Class storeDrawMatrix = this->drawMatrix;

    const size_t storeDim = this->dimension;
    this->dimension = pca_modes.rows();
    this->drawMatrix = transposed(pca_modes);
    this->calculateNN_MIExpansion(orderOfCorrection, norm, func, false);
    this->drawMatrix = storeDrawMatrix;
    this->dimension = storeDim;

    Matrix_Class pca_frequencies(eigenvaluesPCA.rows(), 1u);
    Matrix_Class alpha_i(pca_frequencies.rows(), 1u);
    Matrix_Class quantum_entropy(pca_frequencies.rows(), 1u);
    float_type entropy_sho = 0;
    for (std::size_t i = 0; i < eigenvaluesPCA.rows(); i++)
    {
      pca_frequencies(i, 0u) = sqrt(constants::boltzmann_constant_kb_SI_units * temperatureInK / eigenvaluesPCA(i, 0u));
      alpha_i(i, 0u) = constants::h_bar_SI_units / (sqrt(constants::boltzmann_constant_kb_SI_units * temperatureInK) * sqrt(eigenvaluesPCA(i, 0u)));
      std::cout << "Debug: alpha_i " << alpha_i(i, 0u) << std::endl;
      const double sanityCheck = constants::h_bar_SI_units * pca_frequencies(i, 0u) / constants::boltzmann_constant_kb_SI_units / temperatureInK;
      std::cout << "Debug: sanitycheck " << sanityCheck << std::endl;
      //These are in units S/k_B (therefore: not multiplied by k_B)
      quantum_entropy(i, 0u) = ((alpha_i(i, 0u) / (exp(alpha_i(i, 0u)) - 1)) - log(1 - exp(-1 * alpha_i(i, 0u))));
      std::cout << "Debug: quantum_entropy " << quantum_entropy(i, 0u) << std::endl;
      statistical_entropy(i, 0u) = -1.0 * (log(alpha_i(i, 0u)) -/*this might be plus or minus?!*/ log(sqrt(2. * constants::pi * 2.71828182845904523536)));
      if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
      {
        if (!std::isnan(quantum_entropy(i, 0u)))
          entropy_sho += quantum_entropy(i, 0u) * constants::joules2cal * constants::N_avogadro * constants::boltzmann_constant_kb_SI_units;
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
            if (pca_frequencies(rowNr, 0u) > (temperatureInK * constants::boltzmann_constant_kb_SI_units / (constants::h_bar_SI_units)))
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
    std::cout.unsetf(std::ios_base::floatfield);

    constexpr double conversionFromRawToCalPerMolPerK = constants::joules2cal * constants::N_avogadro * constants::boltzmann_constant_kb_SI_units;

    for (size_t i = 0; i < entropy_kNN.rows(); i++)
    {
      if (this->subDims == std::vector<size_t>() || std::find(this->subDims.begin(), this->subDims.end(), i) != this->subDims.end())
      {
        //These are in units S/k_B (therefore: not multiplied by k_B)
        statistical_entropy(i, 0u) = -1.0 * (log(alpha_i(i, 0u)) -/*this might be plus or minus?!*/ log(sqrt(2. * 3.14159265358979323846 * 2.71828182845904523536)));
        classical_entropy(i, 0u) = -1.0 * (log(alpha_i(i, 0u)) - 1.); // should this be +1??? // The formula written HERE NOW is correct, there is a sign error in the original pape rof Knapp/numata
        
        const double kNN_in_cal_molK = conversionFromRawToCalPerMolPerK * entropy_kNN(i, 0u);
        const double shifted_kNN_cal_molK = this->classicalHarmonicShiftingConstants(i, 0u) + conversionFromRawToCalPerMolPerK * entropy_kNN(i, 0u);
        const double shifted_kNN_raw = shifted_kNN_cal_molK / (conversionFromRawToCalPerMolPerK);
        //entropy_kNN(i, 0u) *= -1.;
        const double anharmonic_correction = classical_entropy(i, 0u) - shifted_kNN_raw;
        entropy_anharmonic(i, 0u) = anharmonic_correction;

        // Debug output for developers
        if (Config::get().general.verbosity >= 4)
        {
          std::cout << "---------------------" << std::endl;
          std::cout << "--- Units: S/k_B ---" << std::endl;
          std::cout << "Mode " << i << ": entropy kNN/kB in J/K: " << entropy_kNN(i, 0u) << "\n";
          const double kNN_in_cal_molK = conversionFromRawToCalPerMolPerK * entropy_kNN(i, 0u);
          const double shifted_kNN_cal_molK = this->classicalHarmonicShiftingConstants(i, 0u) + conversionFromRawToCalPerMolPerK * entropy_kNN(i, 0u);
          //
          std::cout << "Mode " << i << ": entropy kNN in cal/molK: " << kNN_in_cal_molK << "\n";
          std::cout << "Mode " << i << ": shifted entropy kNN in cal/molK: " << shifted_kNN_cal_molK << "\n";
          std::cout << "Mode " << i << ": entropy anharmonic correction: in cal/molK: " << entropy_anharmonic(i, 0u) * conversionFromRawToCalPerMolPerK << "\n";
          std::cout << "Mode " << i << ": classical entropy in cal/molK: " << classical_entropy(i, 0u) * conversionFromRawToCalPerMolPerK << "\n";
          std::cout << "Mode " << i << ": statistical entropy using alpha in cal/molK: " << statistical_entropy(i, 0u) * conversionFromRawToCalPerMolPerK  << "\n";
          std::cout << "Mode " << i << ": std. dev. in SI units: " << stdDevPCAModes(i, 0u) << "\n";
          //
          Matrix_Class covDebug=transposed(pca_modes).covarianceMatrix();

          std::cout << "Mode " << i << ": statistical entropy using std. dev. in SI units: " << 0.5*std::log(2*constants::pi*std::exp(1.)* stdDevPCAModes(i, 0u)* stdDevPCAModes(i, 0u)) << "\n";
          std::cout << "Mode " << i << ": quantum entropy in cal/molK: " << quantum_entropy(i, 0u) * conversionFromRawToCalPerMolPerK << "\n";
          std::cout << "Mode " << i << ": pca freq cm^-1: " << pca_frequencies(i, 0u) / constants::speed_of_light_cm_per_s << "\n";
          std::cout << "Mode " << i << ": alpha (dimensionless, standard deviation): " << alpha_i(i, 0u) << "\n";
          std::cout << "---------------------" << std::endl;
        }

        if (pca_frequencies(i, 0u) < (temperatureInK * constants::boltzmann_constant_kb_SI_units / (constants::h_bar_SI_units)))
        {
          if (std::abs(entropy_anharmonic(i, 0u) / quantum_entropy(i, 0u)) < anharmonicityCutoff)
          {
            entropy_anharmonic(i, 0u) = 0.0;
            std::cout << "Notice: PCA-Mode " << i << " not corrected for anharmonicity (absolute value too small: S/kB [J/K] " << std::abs(entropy_anharmonic(i, 0u) / quantum_entropy(i, 0u)) << " < " << anharmonicityCutoff << ").\n";
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
          std::cout << (temperatureInK * constants::boltzmann_constant_kb_SI_units / (constants::h_bar_SI_units)) / constants::speed_of_light_cm_per_s << "cm^-1,";
          std::cout << "but is " << pca_frequencies(i, 0u) / constants::speed_of_light_cm_per_s << "cm^-1;";
          std::cout << " equipartition is not a valid assumption).\n";
          entropy_anharmonic(i, 0u) = 0.0;
        }

        // Change dimensionless entropy to cal / K * mol
        entropy_anharmonic(i, 0u) *= constants::N_avogadro * constants::boltzmann_constant_kb_SI_units * constants::joules2cal;
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
      {
        std::cout << "---------------------" << std::endl;
        std::cout << "Correction for entropy (anharmonicity, order 1): " << delta_entropy << " cal / (mol * K)\n";
        std::cout << "---------------------" << std::endl;
      }
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
              higher_order_entropy += element.entropyValue * conversionFromRawToCalPerMolPerK;
            }
            countNegativeSecondOrderMIs++;
            if (Config::get().general.verbosity >= 4)
            {
              std::cout << "Notice: Negative 2nd order MI for modes " << element.rowIdent.at(0u) << " and " << element.rowIdent.at(1u);
              std::cout << ": " << element.entropyValue * conversionFromRawToCalPerMolPerK << " cal / (mol * K)\n";
            }
            sumOfNegativeSecondOrderMIs += element.entropyValue * conversionFromRawToCalPerMolPerK;
          }
          else
          {
            higher_order_entropy += element.entropyValue * conversionFromRawToCalPerMolPerK;
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
    transpose(drawMatrix);
    scon::chrono::high_resolution_timer timer;
    std::cout << "--------------------\nBeginning Mutual Information Expansion Procedure" << std::endl;

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
      std::cout << "--- Order: << " << order_N << " ---" << "\n";
      if (norm == kNN_NORM::EUCLEDEAN)
        std::cout << "--- Distance-Norm: << " << "L2" << " ---" << "\n";
      else
        std::cout << "--- Distance-Norm: << " << "L_inf" << " ---" << "\n";
      if (func == kNN_FUNCTION::GORIA)
        std::cout << "--- Function: << " << "Goria" << " ---" << "\n";
      else if (func == kNN_FUNCTION::HNIZDO)
        std::cout << "--- Distance-Norm: << " << "Hnizdo" << " ---" << "\n";
      else
        std::cout << "--- Distance-Norm: << " << "Lombardi" << " ---" << "\n";
      if (ardakaniCorrection)
        std::cout << "--- " << "Using Rahvar&Ardakani's correction scheme" << " ---" << "\n";
    }
    std::size_t count_negative_2MIE_terms = 0u;
    for (size_t i = 0u; i < buffer.size(); i++)
    {
      if (Config::get().general.verbosity >= 4)
      {
        std::cout << "Contribution " << i << ":  Dimensionality " << buffer.at(i).dim << "\n";
        for (size_t j = 0u; j < buffer.at(i).dim; j++)
        {
          std::cout << "Mode " << buffer.at(i).rowIdent.at(j) << "\n";
        }
        std::cout << "Value [=(-1)^(dim+1) * Contribution]: " << std::pow(-1., buffer.at(i).dim + 1.) * buffer.at(i).entropyValue << std::endl;
      }
      if (buffer.at(i).dim == 2u)
      {
        if (buffer.at(i).entropyValue < 0.0)
        {
          count_negative_2MIE_terms += 1u;
        }
      }
      miEntropy += std::pow(-1., buffer.at(i).dim + 1.) * buffer.at(i).entropyValue;
    }
    if (Config::get().general.verbosity > 1 && count_negative_2MIE_terms > 0u)
    {
      std::cout << "WARNING: Counted " << std::to_string(count_negative_2MIE_terms) << " negative 2nd order Mutual Information contributions.\n";
      std::cout << "THIS IS UNPHYSICAL AND CONCERNING!" << std::endl;
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
      std::cout << "Mutual Information Expansion took " << timer << " ." << std::endl;
      std::cout << "Entropy value of expansion: " << std::fixed << std::setprecision(6u) << miEntropy << " nats." << std::endl;

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

    inline bool operator == (const calcBuffer& b) const
    {
      return this->dim == b.dim && this->rowIdent == b.rowIdent && this->entropyValue == b.entropyValue;
    }

    inline bool operator != (const calcBuffer& b) const
    {
      return !((*this) == b);
    }
  };

  std::vector<calcBuffer> calculatedMIs;

  // Calculates the NN entropy for the specified row indicices
  double calculateNNsubentropy(const kNN_NORM norm, const kNN_FUNCTION func, bool const& ardakaniCorrection, std::vector<size_t> const& rowIndices) const
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

      delete[] buffer;
#ifdef _OPENMP
    }
#endif

    //transpose(drawMatrix);

    // Eucledean ArdakaniSum
    if (norm == kNN_NORM::EUCLEDEAN && ardakaniCorrection)
    {
      KahanAccumulation<double> kahan_acc_eucl_ardakani_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
      {
        double const& current_log_dist = log(eucl_kNN_distances_ardakani_corrected(0, i));
        if (!std::isinf(current_log_dist))
          kahan_acc_eucl_ardakani_sum = KahanSum(kahan_acc_eucl_ardakani_sum, current_log_dist);
        else if (Config::get().general.verbosity >= 5)
        {
          std::cout << "Warning: Detected kNN-distance equal to 0.0. Ignoring this distance.\n";
        }
      }


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
      {
        double const& current_log_dist = log(maxnorm_kNN_distances_ardakani_corrected(0, i));
        if (!std::isinf(current_log_dist))
          kahan_acc_max_ardakani_sum = KahanSum(kahan_acc_max_ardakani_sum, current_log_dist);
        else if (Config::get().general.verbosity >= 5)
        {
          std::cout << "Warning: Detected kNN-distance equal to 0.0. Ignoring this distance.\n";
        }
      }

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
      {
        const double current_log_dist = log(eucl_kNN_distances(0, i));
        if (!std::isinf(current_log_dist))
          kahan_acc_eucl_sum = KahanSum(kahan_acc_eucl_sum, current_log_dist);
        else if (Config::get().general.verbosity >= 5)
        {
          std::cout << "Warning: Detected kNN-distance equal to 0.0. Ignoring this distance.\n";
        }
      }

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
      //for (size_t i = 0u; i < numberOfDraws; i++)
      //  std::cout << eucl_kNN_distances(0, i) << std::endl;
      return lombardi_without_correction; // Lombardi Entropy
    }
    else if (norm == kNN_NORM::MAXIMUM && !ardakaniCorrection)
    {
      KahanAccumulation<double> kahan_acc_max_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
      {
        double const& current_log_dist = log(maxnorm_kNN_distances(0, i));
        if (!std::isinf(current_log_dist))
          kahan_acc_max_sum = KahanSum(kahan_acc_max_sum, current_log_dist);
        else if (Config::get().general.verbosity >= 5)
        {
          std::cout << "Warning: Detected kNN-distance equal to 0.0. Ignoring this distance.\n";
        }
      }
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
    std::vector<size_t> curRowIndices = std::vector<size_t>{}, const size_t curDim = 0u, size_t startIter = 0u) const
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
    std::vector<size_t> curRowIndices = std::vector<size_t>{}, const size_t curDim = 0u, size_t startIter = 0u) const
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

  // Full dimensional nearest neighbor entropy computation (without MI-Expansion)
  double calculateNN(Matrix_Class const& currentData, const kNN_NORM norm, bool const& ardakaniCorrection, const kNN_FUNCTION func = kNN_FUNCTION::HNIZDO) const
  {
    std::cout << "Commencing full-dimensional kNN-Entropy calculation." << std::endl;

    const unsigned int dimensionality = this->subDims != std::vector<size_t>() ? this->subDims.size() : currentData.cols();
    std::cout << "Dimensionality: " << dimensionality << std::endl;
    if (dimensionality != this->dimension)
    {
      std::cout << "Note: NN-Entropy will be evaluated in truncated dimensionality (Original: " << this->dimension << ")." << std::endl;
    }

    
    Matrix_Class dimPurgedDrawMatrix = Matrix_Class(static_cast<uint_type>(currentData.rows()), static_cast<uint_type>(dimensionality));
    if (this->subDims != std::vector<size_t>())
    {
      for (size_t i = 0u; i < this->subDims.size(); i++)
      {
        for (size_t j = 0u; j < numberOfDraws; j++)
        {
          dimPurgedDrawMatrix(j, i) = currentData(j, this->subDims.at(i));
        }
      }
    }
    else
    {
      dimPurgedDrawMatrix = currentData;
    }

    //Neccessarry
    transpose(dimPurgedDrawMatrix);

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
      for (size_t j = 0; j < currentData.cols(); j++)
      {
        for (unsigned int i = 0u; i < dimensionality; i++)
        {
          if (ardakaniCorrection_minimumValueInDataset.at(i) > currentData(i, j))
            ardakaniCorrection_minimumValueInDataset.at(i) = currentData(i, j);

          if (ardakaniCorrection_maximumValueInDataset.at(i) < currentData(i, j))
            ardakaniCorrection_maximumValueInDataset.at(i) = currentData(i, j);
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
      std::cout << std::fixed;
      std::cout << std::setprecision(6u);
      std::cout << ": " << returnValue << " nats [";
      std::cout << std::to_string(returnValue * constants::boltzmann_constant_kb_gaussian_units * constants::eV2kcal_mol * 1000) << " cal/(mol*K)].";
      std::cout << std::endl;
      std::cout << "NN Calculation took " << timer << " ." << std::endl;
    }
    return returnValue;
  }

  

  Matrix_Class unmassweightPCAModes(Matrix_Class const& assocRedMasses, Matrix_Class const& pcaModes) const
  {
    if (assocRedMasses.cols() != 1u || pcaModes.rows() != assocRedMasses.rows())
    {
      throw std::logic_error("Cannot un-massweight PCA Modes: Dimensionality of matrices is wrong. Aborting!");
      return Matrix_Class();
    }
    Matrix_Class unweightedPcaModes(pcaModes);
    for (std::size_t i = 0u; i < assocRedMasses.rows(); ++i)
    {
      for (std::size_t j = 0u; j < pcaModes.cols(); ++j)
      {
        const double toBeDivided = unweightedPcaModes(i, j);
        const double divisor = std::sqrt(assocRedMasses(i, 0u));
        unweightedPcaModes(i,j) = toBeDivided/divisor;
      }
    }
    //std::cout << "Debug unweighted:\n" << unweightedPcaModes << "\n";
    return unweightedPcaModes;
  }
  
};

