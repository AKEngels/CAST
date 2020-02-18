#pragma once
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <omp.h>
#include <limits>
#include <string>
#include <cmath>
#include "entropy.h"
#include <tuple>
#include "kaham_summation.h"
#include "constants.h"
#include "cubature.h"
#include "histogram.h"
#include "scon_chrono.h"

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
float_type ardakaniCorrection1D(float_type const& globMin, float_type const& globMax, float_type const& currentPoint, float_type const& NNdistance)
{
  if (currentPoint - NNdistance * 0.5 < globMin && currentPoint != globMin)
    return currentPoint + NNdistance * 0.5 - globMin;
  else if (currentPoint + NNdistance * 0.5 > globMax && currentPoint != globMax)
    return globMax - (currentPoint - NNdistance * 0.5);
  return NNdistance;
}

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
    if (radiiOfHyperEllipsoid.at(i) < 0.f)
    {
      radiiOfHyperEllipsoid.at(i) *= -1.f;
    }
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
    // due to numerics, this can become zero or below. -NaN results and propagates
    if (radiiOfBox.at(i) < 0.f)
    {
      radiiOfBox.at(i) *= -1.f;
    }
  }

  T volumeOfBox = 1.;
  for (unsigned int i = 0u; i < radiiOfBox.size(); i++)
  {
    volumeOfBox *= radiiOfBox.at(i);
  }

  const T equivalentRadiusOfHyperSquare = std::pow(volumeOfBox, 1. / radiiOfBox.size());
  return equivalentRadiusOfHyperSquare;

}

void scalePCACoordinatesForQuasiHarmonicTreatment(Matrix_Class& modes, float_type const& temperatureInK)
{
  pow(modes, -1.);
  modes = modes * 1.05457172647 * 10e-34 / (sqrt(1.380648813 * 10e-23 * temperatureInK));
  //
}

//For reading in data of a gaussian mixture model from file
//GMM Files can be written using armadillo in the PCA task
struct GMM_data
{
  unsigned int numberOfDimensions = 0u;
  unsigned int numberOfGaussians = 0u;
  std::vector<double> weights;
  std::vector<Matrix_Class> mean;
  std::vector<Matrix_Class> covariance;
  GMM_data() {};
  GMM_data(std::string name)
  {

    std::ifstream gmm_stream(name, std::ios::in);
    std::string line;
    
    //int dimensions = std::stoi(line.substr(13, 2));
    while (std::getline(gmm_stream, line))
    {
      if (line.find("Number of Gaussians") != std::string::npos)
      {
        numberOfGaussians = static_cast<size_t>(stoi(line.substr(21)));
        weights = std::vector<double>(numberOfGaussians, 0.);
        continue;
      }
      if (line.find("Number of Dimensions: ") != std::string::npos)
      {
        numberOfDimensions = static_cast<size_t>(stoi(line.substr(22)));
        mean = std::vector<Matrix_Class>(numberOfGaussians, Matrix_Class(numberOfDimensions, 1u, 0.));
        covariance = std::vector<Matrix_Class>(numberOfGaussians, Matrix_Class(numberOfDimensions, numberOfDimensions, 0u));
        continue;
      }
      if (line.find("Mean of Gaussian ") != std::string::npos)
      {
        std::getline(gmm_stream, line);
        std::istringstream iss(line);
        
        for (unsigned int i = 0u; i < this->numberOfDimensions; i++)
        {
          std::string s;
          std::getline(iss, s, ' ');
          try
          {
            mean.at(0)(i, 0) = std::stod(s);
          }
          catch (std::out_of_range e)
          {
            mean.at(0)(i, 0) = 0;
          }
        }
        for (unsigned int j = 1u; j < this->numberOfGaussians; j++)
        {
          std::getline(gmm_stream, line);
          std::getline(gmm_stream, line);
          iss = std::istringstream(line);
          for (unsigned int i = 0u; i < this->numberOfDimensions; i++)
          {
            std::string s;
            std::getline(iss, s, ' ');
            while (s == std::string())
            {
              std::getline(iss, s, ' ');
            }
            
            try
            {
              mean.at(j)(i, 0) = std::stod(s);
            }
            catch (std::out_of_range e)
            {
              mean.at(j)(i, 0) = 0;
            }
          }
        }
        continue;
      }
      if (line.find("Covariance Matrix of Gaussian ") != std::string::npos)
      {
        for (unsigned int j = 0u; j < covariance.at(0).rows(); j++)
        {
          for (unsigned int k = j; k < covariance.at(0).cols(); k++)
          {
            std::getline(gmm_stream, line);
            line = line.substr(6);
            covariance.at(0)(j, k) = std::stod(line);
            try
            {
              covariance.at(0)(j, k) = std::stod(line);
            }
            catch (std::out_of_range e)
            {
              covariance.at(0)(j, k) = 0;
            }
            covariance.at(0)(k, j) = covariance.at(0)(j, k);
          }
        }

        for (unsigned int i = 1u; i < this->numberOfGaussians; i++)
        {
          std::getline(gmm_stream, line);
          std::getline(gmm_stream, line);
          std::getline(gmm_stream, line);
          for (unsigned int j = 0u; j < covariance.at(i).rows(); j++)
          {
            for (unsigned int k = j; k < covariance.at(i).cols(); k++)
            {
              std::getline(gmm_stream, line);
              line = line.substr(6);
              try
              {
                covariance.at(i)(j, k) = std::stod(line);
              }
              catch (std::out_of_range e)
              {
                covariance.at(i)(j, k) = 0;
              }

              covariance.at(i)(k, j) = covariance.at(i)(j, k);
            }
          }
        }
        continue;
      }
      if (line.find("Weight of Gaussian ") != std::string::npos)
      {
        for (unsigned int i = 0u; i < this->numberOfGaussians; i++)
        {
          if (line == std::string())
            break;

          weights.at(i) = std::stod(line.substr(line.find(":") + 1u));
          try
          {
            weights.at(i) = std::stod(line.substr(line.find(":") + 1u));
          }
          catch (std::out_of_range e)
          {
            weights.at(i) = 0;
          }
          if (!std::getline(gmm_stream, line))
          {
            break;
          }
        }
        continue;
      }
    }
  }
};


// Class for PDFunctions
//
// PDF Function needs to be centered at 0,0
class ProbabilityDensity
{
public:
  GMM_data gmmdatafromfile;

  ProbabilityDensity(int ident_) : identif(ident_), maximumOfPDF(- std::numeric_limits<double>::max())
  {
    if (ident_ == 0)
    {
      this->dimension = 1;
      maximumOfPDF = 1. / sqrt(2. * pi);
      analyticEntropy_ = log(1. * sqrt(2. * pi * e));
      PDF = [&, this](std::vector<double> const& x, std::vector<size_t> const& subdims = std::vector<size_t>())
      {
        if (x.size() != 1)
          throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");
        return (1. / sqrt(2. * pi * 1. * 1.)) * exp(-1.*(x.at(0)*x.at(0)) / double(2. * 1. * 1.));
      };
      //PDFrange = 5.920; // range incorporating function values up to under 10e-7
      this->m_identString = "Normal sig=1";
    }
    else if (ident_ == 1)
    {
      this->dimension = 1;
      // Gaussian with sigma 10
      PDF = [&, this](std::vector<double> const& x, std::vector<size_t> const& subdims = std::vector<size_t>())
      {
        if (x.size() != 1)
          throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");
        return (1. / sqrt(2. * pi * 10. * 10.)) * exp(-1.*(x.at(0)*x.at(0)) / double(2. * 10. * 10.));
      };
      maximumOfPDF = 1. / (10. * sqrt(2. * pi));
      analyticEntropy_ = log(10. * sqrt(2. * pi * e));
      //PDFrange = 55.136; // range incorporating function values up to under 10e-7
      this->m_identString = "Normal sig=10";
    }
    else if (ident_ == 2)
    {
      this->dimension = 1;
      // Parabola normated to integrate to 1 over [0;1]
      PDF = [&, this](std::vector<double> const& x, std::vector<size_t> const& subdims = std::vector<size_t>())
      {
        if (x.size() != 1)
          throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");
        if (x.at(0) < -1. || x.at(0) > 1.)
          return 0.;
        else
          return (3. / 4.) * (-1. * ((x.at(0))*(x.at(0))) + 1.);
      };
      maximumOfPDF = PDF(std::vector<double>{ 0. }, std::vector<size_t>());
      analyticEntropy_ = 0.568054;
      PDFrange = std::make_shared<std::pair<double, double>>(-1., 1.); // Function is only defined in [-1;1]
      this->m_identString = "parabolic";
    }
    else if (ident_ == 3)
    {
      this->dimension = 1;
      // Beta Distribution https://en.wikipedia.org/wiki/Differential_entropy
      PDF = [&, this](std::vector<double> const& x, std::vector<size_t> const& subdims = std::vector<size_t>())
      {
        if (x.size() != 1)
          throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");
        constexpr double alpha = 1.5;
        constexpr double beta = 9.;
        std::vector<double> xtemp(x);
        xtemp.at(0) = abs(x.at(0));
        if (xtemp.at(0) > 1.)
          return 0.;
        else
          return pow(xtemp.at(0), alpha - 1.) * pow(1. - xtemp.at(0), beta - 1.) / (tgamma(alpha)*tgamma(beta) / tgamma(alpha + beta));
      };
      maximumOfPDF = 4.;
      constexpr double alpha = 1.5;
      constexpr double beta = 9.;
      analyticEntropy_ = log(tgamma(alpha)*tgamma(beta) / tgamma(alpha + beta)) - (alpha - 1.) * (digammal(alpha) - digammal(alpha + beta)) - (beta - 1.) * (digammal(beta) - digammal(alpha + beta));
      PDFrange = std::make_shared<std::pair<double, double>>(0., 1.); // Function is only defined in [0;1]
      this->m_identString = "Beta Distr.";
    }
    else if (ident_ == 9)
    {
      // Read GMM from file
      gmmdatafromfile = GMM_data("gmmfile.dat");
      std::cout << "Read gmmfile.dat with " << gmmdatafromfile.numberOfDimensions << " dimensions and " << gmmdatafromfile.numberOfGaussians << " gaussians." << std::endl;
      this->dimension = gmmdatafromfile.numberOfDimensions;
      //
      std::vector<Matrix_Class> covariance_inv(gmmdatafromfile.covariance);
      std::vector<double> covariance_det(covariance_inv.size(),0);
      for (std::size_t i = 0u; i < covariance_inv.size(); ++i)
      {
        covariance_inv[i] = covariance_inv[i].inversed();
        covariance_det[i] = gmmdatafromfile.covariance[i].determ();
      }
      const double renormalizationConstant = Config::get().entropytrails.renormalizationOfPDF;
      if (renormalizationConstant != 1.0)
        std::cout << "Renormalizing PDF as requested by multiplication with " << renormalizationConstant << "." << std::endl;
      // Multivariate gaussian
      PDF = [=](std::vector<double> const& x, std::vector<size_t> const& subdims = std::vector<size_t>())
      {
        if (subdims == std::vector<size_t>())
        {
          if (x.size() != this->dimension)
            throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");

          Matrix_Class input(x.size(), 1u);
          for (std::size_t i = 0u; i < x.size(); i++)
            input(i, 0u) = x.at(i);

          unsigned int const numberOfGaussians = gmmdatafromfile.numberOfGaussians;

          std::vector<double> const weight = gmmdatafromfile.weights;
          std::vector<Matrix_Class> const mean = gmmdatafromfile.mean;

          double returnValue = 0.;
          for (unsigned int i = 0u; i < numberOfGaussians; i++)
          {
            Matrix_Class tempMatr1 = Matrix_Class(covariance_inv[i] * Matrix_Class(Matrix_Class(input) - mean[i]));
            Matrix_Class value = Matrix_Class(transposed(Matrix_Class(input - mean[i])) * tempMatr1);
            const double value_ = value(0u, 0u);
            const double prefactor = 1. / std::sqrt(std::pow(2 * ::constants::pi, x.size()) * covariance_det[i]);

            const double addition = prefactor * std::exp(value_ * -0.5) * weight[i];
            if (std::isnormal(addition))
              returnValue += addition;
            else if (addition != 0.0)
              std::cout << "Non-normal number in drawing from GMM PDF. This should not happen and is cause for concern.\n";
          }

          return returnValue * renormalizationConstant;

        }
        else
        {
          if (x.size() != subdims.size())
            throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");

          Matrix_Class input(x.size(), 1u);
          for (unsigned int i = 0u; i < x.size(); i++)
            input(i, 0u) = x.at(i);

          unsigned int const numberOfGaussians = gmmdatafromfile.numberOfGaussians;

          std::vector<double> weight = gmmdatafromfile.weights;
          std::vector<Matrix_Class> mean;// .gmmdatafromfile.mean;
          std::vector<Matrix_Class> covariance;// = gmmdatafromfile.covariance;


          // Build smaller matrices of mean and covariance for subdimensions
          for (auto&& mean_current : gmmdatafromfile.mean)
          {
            Matrix_Class current_new(subdims.size(), 1u);
            for (unsigned int i = 0u; i < subdims.size(); i++)
            {
              current_new(i, 0u) = mean_current(subdims.at(i), 0u);
            }
            mean.push_back(current_new);
          }

          for (auto&& covariance_current : gmmdatafromfile.covariance)
          {
            Matrix_Class current_new(subdims.size(), subdims.size());
            for (unsigned int i = 0u; i < subdims.size(); i++)
            {
              for (unsigned int j = 0u; j < subdims.size(); j++)
              {
                current_new(i, j) = covariance_current(subdims.at(i), subdims.at(j));
              }
            }
            covariance.push_back(current_new);
          }



          double returnValue = 0.;
          for (unsigned int i = 0u; i < numberOfGaussians; i++)
          {
            Matrix_Class value = Matrix_Class(transposed(Matrix_Class(input - mean[i])) * covariance[i].inversed() * Matrix_Class(input - mean[i]));
            const double cov_determ = covariance[i].determ();
            const double prefactor = 1. / std::sqrt(std::pow(2 * ::constants::pi, x.size()) * cov_determ);

            const double addition = prefactor * std::exp(value(0u, 0u) * -0.5) * weight[i];
            if (addition != addition)
              continue;
            returnValue += addition;
          }

          return returnValue;
        }

      };

      if (Config::get().entropytrails.rangeForGMM.first == Config::get().entropytrails.rangeForGMM.second
        && Config::get().entropytrails.rangeForGMM.first == 0.)
      {
        PDFrange = std::make_shared<std::pair<double, double>>(-1e-11, 4e-11);
      }
      else
      {
        PDFrange = std::make_shared<std::pair<double, double>>(Config::get().entropytrails.rangeForGMM);
      }

      std::cout << "Range for GMM probability density is from " << PDFrange->first << " to " << PDFrange->second << "." << std::endl;


      //
      //find maximum of pdf by drawing many times
      std::random_device rd;
      std::mt19937 gen(rd());
      std::cout << "Finding maximum of PDF by drawing 100k times, taking highest value x10." << std::endl;
      std::uniform_real_distribution<double> unifDistrRange(PDFrange->first, PDFrange->second);
      const double absrange = PDFrange->second - PDFrange->first;

      //if (false)
      //{
        for (unsigned int i = 0u; i < 100000u; i++)
        {
          std::vector<double> drawUnifWithRange;
          for (unsigned int currentDim = 0u; currentDim < this->dimension; currentDim++)
          {
            double randNum = unifDistrRange(gen);
            //std::cout << "RandNum: " << randNum << "\n";
            drawUnifWithRange.push_back(randNum);
          }
          const double pdfvalue = PDF(drawUnifWithRange, std::vector<size_t>());
          //std::cout << pdfvalue << "\n";
          maximumOfPDF = std::max(pdfvalue, maximumOfPDF);
        }
        maximumOfPDF *= 10.;
        std::cout << "Estimated maximum of current PDF is " << maximumOfPDF << "." << std::endl;
        if (!std::isnormal(maximumOfPDF))
        {
          throw std::runtime_error("Estimated maximum of PDF is not a normal number. Critical Error. Aborting.");
        }
      //}
      //
      //
      std::cout << "Additionally searching maximum of PDF by evaluating at means of gaussians.." << std::endl;
      double maximum2 = 0.;
      for (unsigned int i = 0u; i < gmmdatafromfile.numberOfGaussians; i++)
      {
        std::vector<double> currentMeans;
        for (unsigned int j = 0u; j < gmmdatafromfile.mean.at(i).rows(); j++)
        {
          currentMeans.push_back(gmmdatafromfile.mean.at(i)(j, 0));
        }
        const double pdfvalue = PDF(currentMeans, std::vector<size_t>());
        //std::cout << pdfvalue << "\n";
        maximum2 = std::max(pdfvalue, maximum2);
      }
      maximum2 *= 10.;
      std::cout << "Estimated maximum at gaussian means of current PDF is " << maximum2 << "." << std::endl;
      maximumOfPDF = std::max(maximumOfPDF, 2*maximum2);
      std::cout << "Estimated final maximum is " << maximumOfPDF << "." << std::endl;
      if (!std::isnormal(maximumOfPDF))
      {
        throw std::runtime_error("Estimated maximum of PDF is not a normal number. Critical Error. Aborting.");
      }
      analyticEntropy_ = std::numeric_limits<double>::quiet_NaN();
      this->m_identString = "GMMfile";
    }
  };

  std::vector<std::vector<double>> draw(size_t numberOfSamples, std::vector<size_t> const subdims)
  {
    //
    std::cout << "Starting to draw from probability distribution..." << std::endl;
    scon::chrono::high_resolution_timer exec_timer1;
    std::vector<std::vector<double>> draws;
    std::function<double(std::vector<double> const& x, std::vector<size_t> const& subdims)> cPDF = this->PDF;

    //Target PDF in 0.1er Schritten auswerten und schauen wann es kleiner ist
    if (PDFrange == nullptr)
      PDFrange = std::make_shared<std::pair<double, double>>(meaningfulRange());
    auto cPDFrange = PDFrange;

    const size_t currentDimensionality = subdims.size() == 0 ? this->dimension : subdims.size();
    
#ifdef _OPENMP
    std::cout << "Using openMP parallelization for drawing with " + std::to_string(omp_get_max_threads()) + " threads." << std::endl;
    std::cout << "Priting a random number for every openMP thread. Please check that these numbers differ.\n";
#pragma omp parallel firstprivate(currentDimensionality, cPDFrange,cPDF,subdims ) \
    shared(draws)
    {
#endif
      std::random_device rd;
      std::mt19937 gen(rd());
#ifdef _OPENMP
#pragma omp critical
      {
#endif
        gen.seed(rd() * (omp_get_thread_num() + 1));
        std::cout << "OMP Thread Num: " << omp_get_thread_num() + 1 << "\n      Random Number Drawn:" << gen() << std::endl;
#ifdef _OPENMP
      }
#endif
      std::uniform_real_distribution<double> unifDistr(0, 1);
      std::uniform_real_distribution<double> unifDistrRange(cPDFrange->first, cPDFrange->second);
      const double absrange = cPDFrange->second - cPDFrange->first;
      const double unifWithRangeProbability = 1. / absrange;
      const double k = maximumOfPDF * 1.05 / unifWithRangeProbability;
#ifdef _OPENMP
      auto const n_omp = static_cast<std::ptrdiff_t>(numberOfSamples);
#pragma omp for
      for (std::ptrdiff_t n = 0; n < n_omp; ++n)
#else
      for (std::size_t n = 0u; n < numberOfSamples; ++n)
#endif
      {
        bool inRange = false;
        while (!inRange)
        {
          std::vector<double> drawUnifWithRange;
          for (std::size_t currentDim = 0u; currentDim < currentDimensionality; currentDim++)
          {
            const double randNum = unifDistrRange(gen);
//#ifdef _OPENMP
//#pragma omp critical
//            {
//              std::cout << "Randnum: " << randNum << " ThreadID: " << std::to_string(omp_get_thread_num()) << std::endl;
//            }
//#endif
            drawUnifWithRange.push_back(randNum);
          }
          const double drawUnif = unifDistr(gen);
          const double pdfvalue = cPDF(drawUnifWithRange, subdims);
          const double p = pdfvalue / (k * unifWithRangeProbability);
          inRange = drawUnif < p;
          if (inRange)
          {
              draws.push_back(drawUnifWithRange);
              if (draws.size() % 5 == 0)
#ifdef _OPENMP
#pragma omp critical
              {
#endif
                std::cout << "Number of draws currently finished: " << draws.size() << std::endl;
#ifdef _OPENMP
              }
#endif
          }
        }
      }
#ifdef _OPENMP
    }
#endif
    std::cout << "Drawing Samples took " << exec_timer1 << " to complete." << std::endl;
    return draws;
  };

  int ident() const { return this->identif; };

  std::string identString() const
  {
    return this->m_identString;
  }

  double analyticEntropy() const {
    return this->analyticEntropy_;
  };

  auto function() -> std::function<double(std::vector<double> const& x, std::vector<size_t> const& subdims)> const
  {
    return this->PDF;
  };

  size_t getDimension() const
  {
    return this->dimension;
  }

  std::pair<double, double> meaningfulRange() const
  {
    if (PDFrange != nullptr)
      return *PDFrange;
    if (this->dimension == 1)
    {
      if (PDFrange != nullptr)
        return *PDFrange;
      // Target PDF in 0.001er Schritten auswerten und schauen wann es kleiner ist
      // als 10e-7.
      // Ident 0: 5.917
      // Ident 1: 55.1349999999
      // Ident 2: 1.0000000
      // Ident 3:
      double rangeOne = 0.f;
      double newPDFrange = 0.;
      bool run = true;
      newPDFrange = 0.;
      while (run)
      {
        if (PDF(std::vector<double> {newPDFrange}, std::vector<size_t>()) < 0.00000001)
          break;
        else
          newPDFrange += 0.001;
      }
      rangeOne = newPDFrange;
      newPDFrange = 0.;
      while (run)
      {
        if (PDF(std::vector<double> {newPDFrange}, std::vector<size_t>()) < 0.000000001)
          break;
        else
          newPDFrange -= 0.005;
      }

      return std::make_pair(newPDFrange, rangeOne);
    }
    else
    {
      throw std::logic_error("Meaningful range not implemented for multivariate distributions.");
    }
  }

private:
  double analyticEntropy_ = 0.f;
  double maximumOfPDF;
  int identif;
  size_t dimension = 0;
  std::shared_ptr<std::pair<double, double>> PDFrange = nullptr;
  std::string m_identString = "ERROR_NAME";
  std::function<double(std::vector<double> const& x, std::vector<size_t> const& subdims)> PDF;
};


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
  entropyobj(size_t iter_, ProbabilityDensity probdens_, std::vector<size_t> const& subdimsGMM_) :
    numberOfDraws(iter_), dimension(probdens_.getDimension()), drawMatrix(iter_, probdens_.getDimension()),
    subDims(subdimsGMM_)
  {
    // Check if draw exists
    std::ifstream myfile;
    std::string line;
    std::string filename;
    filename = std::string("draw_i" + std::to_string(numberOfDraws) + "_d" + std::to_string(dimension) + "_ident" + std::to_string(probdens_.ident()) + ".txt");
    myfile.open(filename);
    if (myfile.good())
    {
      std::cout << "Reading draw from file." << std::endl;
      int j = 0;
      int currentDim = 0;
      while (std::getline(myfile, line))
      {
        drawMatrix(j, currentDim) = std::stod(line);
        j++;
        if (j == this->numberOfDraws)
        {
          currentDim++;
          j = 0;
        }
      }
    }
    else
    {
      std::cout << "Creating new draw." << std::endl;
      // Draw 
      std::vector<std::vector<double>> draws = probdens_.draw(numberOfDraws, subDims);

      //Sort samples (entirely optional) if dimension == 1
      if (this->dimension == 1)
        std::stable_sort(draws.begin(), draws.end());

      for (unsigned int n = 0; n < numberOfDraws; ++n)
        for (unsigned int currentDim = 0; currentDim < this->dimension; ++currentDim)
          drawMatrix(n, currentDim) = draws[n][currentDim];

      // Write 
      std::ofstream myfile2;
      myfile2.open(std::string("draw_i" + std::to_string(numberOfDraws) + "_d" + std::to_string(dimension) + "_ident" + std::to_string(probdens_.ident()) + ".txt"));
      for (unsigned int currentDim = 0; currentDim < this->dimension; ++currentDim)
      {
        for (size_t w = 0u; w < drawMatrix.rows(); w++)
        {
          myfile2 << std::setw(30) << std::scientific << std::setprecision(15) << drawMatrix(w, currentDim);
          myfile2 << "\n";
        }
      }
      myfile2.close();
    }
    myfile.close();
  }

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

struct informationForCubatureIntegration
{
  ProbabilityDensity& probdens;
  std::vector<size_t>& subdims;
  informationForCubatureIntegration(ProbabilityDensity& probdens_in, std::vector<size_t>& subdims_in)
    : probdens(probdens_in), subdims(subdims_in) {}
};

int cubaturefunctionEntropy(unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval)
{
  informationForCubatureIntegration& info = *((informationForCubatureIntegration *)fdata);

  if (ndim != info.probdens.getDimension())
    throw;
  if (fdim != 1)
    throw;

  for (int j = 0; j < static_cast<int>(npts); ++j) 
  { // evaluate the integrand for npts points
    std::vector<double> input;
    for (unsigned int i = 0u; i < ndim; i++)
    {
      if (info.subdims.size() != 0)
      {
        if (std::find(info.subdims.begin(), info.subdims.end(), i) != info.subdims.end())
          input.push_back(x[j*ndim + i]);
      }
      else
      {
        input.push_back(x[j*ndim + i]);
      }
    }

    double probdens_result = info.probdens.function()(input, info.subdims);
    if (probdens_result != 0.)
    {
      double returner = info.probdens.function()(input, info.subdims) * std::log(info.probdens.function()(input, info.subdims));
      fval[j] = returner;
    }
    else
    {
      fval[j] = 0.;
    }
  }
  return 0;
}

int cubaturefunctionProbDens(unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval)
{
  informationForCubatureIntegration& info = *((informationForCubatureIntegration *)fdata);

  if (ndim != info.probdens.getDimension())
    throw;
  if (fdim != 1)
    throw;

  for (int j = 0; j < static_cast<int>(npts); ++j)
  { // evaluate the integrand for npts points
    std::vector<double> input;
    for (unsigned int i = 0u; i < ndim; i++)
    {
      if (info.subdims.size() != 0)
      {
        if (std::find(info.subdims.begin(), info.subdims.end(), i) != info.subdims.end())
          input.push_back(x[j*ndim + i]);
      }
      else
      {
        input.push_back(x[j*ndim + i]);
      }
    }

    double probdens_result = info.probdens.function()(input, info.subdims);
    fval[j] = probdens_result;
  }
  return 0;
}


// Calculated entropy object calculates estimated entropy
class calculatedentropyobj : public entropyobj
{
public:
  size_t kNN;
  double mean, standardDeviation;
  double analyticalEntropy;
  double empiricalNormalDistributionEntropy;
  Matrix_Class pcaModes;

  calculatedentropyobj(size_t k_, entropyobj const& obj) :
    entropyobj(obj),
    kNN(k_),
    pcaModes(Matrix_Class(0u, 0u)),
    mean(std::numeric_limits<double>::quiet_NaN()),
    standardDeviation(std::numeric_limits<double>::quiet_NaN()),
    analyticalEntropy(std::numeric_limits<double>::quiet_NaN()),
    empiricalNormalDistributionEntropy(std::numeric_limits<double>::quiet_NaN())
  {

  }

  void setAnalyticalEntropy(ProbabilityDensity probdens)
  {
    analyticalEntropy = probdens.analyticEntropy();
    writeToCSV("entropy.csv", "analytic", analyticalEntropy, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::HNIZDO, this->numberOfDraws);
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
      Matrix_Class cov_matr,meanPerDim;
      
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
      const double determinant = cov_matr.determ();

      const double gaussentropy = 0.5 * log(cov_matr.determ()) + double(dimensionality) / 2. * std::log(2. * ::constants::pi * ::constants::e);
      empiricalNormalDistributionEntropy =  gaussentropy;
      std::cout << "Entropy of multivariate Gaussian with sigma of the draws (statistical entropy, no thermodyn. units): " << gaussentropy << std::endl;
      writeToCSV("entropy.csv", "empricial_gaussian_statistical", gaussentropy, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::HNIZDO, this->numberOfDraws);
      //Covariance Matrix


      //Calculate PCA Frequencies in quasi-harmonic approximation and Entropy in SHO approximation; provides upper limit of entropy
      Matrix_Class eigenvalues;
      Matrix_Class eigenvectors;
      cov_matr.eigensym(eigenvalues, eigenvectors);
      Matrix_Class pca_frequencies(eigenvalues.rows(), 1u);
      Matrix_Class alpha_i(pca_frequencies.rows(), 1u);
      Matrix_Class quantum_entropy(pca_frequencies.rows(), 1u);
      float_type entropy_sho = 0;
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "Gaussian Entropy in qQH-approximation is estimated from eigenvalues of draw matrix.\n";
        std::cout << "NOTICE: The draw-matrix will be transformed to obtain thermodynamically valid units." << std::endl;
        std::cout << "Using Temperature: " << Config::get().entropy.entropy_temp << " Kelvin." << std::endl;
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
      writeToCSV("entropy.csv", "qQH_approx", gaussentropy, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::HNIZDO, this->numberOfDraws);
    }
  }

  double cubatureIntegrationEntropy(ProbabilityDensity probdens)
  {
    std::cout << "Commencing cubature integration of entropy." << std::endl;
    const unsigned int dimensionality = this->subDims != std::vector<size_t>() ? this->subDims.size() : this->dimension;
    std::cout << "Dimensionality: " << dimensionality << std::endl;
    std::cout << "Integrating in a hypercube region ranging from " << probdens.meaningfulRange().first;
    std::cout << " to " << probdens.meaningfulRange().second << " in each dimension." << std::endl;

    const double min_ = (probdens.meaningfulRange().first);
    const double max_ = (probdens.meaningfulRange().second);
    double* xmin, *xmax, *val, *err;
    xmin = new double[this->dimension];
    xmax = new double[this->dimension];
    val = new double[this->dimension];
    err = new double[this->dimension];
    for (unsigned int i = 0u; i < this->dimension; i++)
    {
      xmin[i] = min_;
      xmax[i] = max_;
      val[i] = 0.;
      err[i] = 0.;
    }

    informationForCubatureIntegration info(probdens, this->subDims);

    const double range = std::abs(probdens.meaningfulRange().first - probdens.meaningfulRange().second);

    std::cout << "Absolute desired error (L2 norm) for cubature integration is " << std::scientific <<
      std::setw(5) << Config::get().entropytrails.errorThresholdForCubatureIntegration << "." << std::endl;

    if (Config::get().entropytrails.cubatureAlsoProbdens)
    {
      scon::chrono::high_resolution_timer exec_timer1;

      int returncode = hcubature_v(1, cubaturefunctionProbDens, &(info),
        static_cast<unsigned int>(this->dimension), xmin, xmax,
        0u, 0u, Config::get().entropytrails.errorThresholdForCubatureIntegration, ERROR_L2, val, err);

      if (returncode != 0)
      {
        std::cout << "Returncode was not zero. An error occured." << std::endl;
      }

      std::cout << "Integral of the Probability Density is: " << std::defaultfloat << val[0] << " with error " << err[0] << " (this should be close to 1.0)." << std::endl;
      std::cout << "Integration took " << exec_timer1 << " to complete." << std::endl;
    }
    scon::chrono::high_resolution_timer exec_timer2;
    int returncode = hcubature_v(1, cubaturefunctionEntropy, &(info),
      static_cast<unsigned int>(this->dimension), xmin, xmax,
      0u, 0u, Config::get().entropytrails.errorThresholdForCubatureIntegration, ERROR_L2, val, err); // Switched to absolute error

    if (returncode != 0)
    {
      std::cout << "Returncode of cubature integration subroutine was not zero. An error occured." << std::endl;
    }

    const double value = val[0];

    std::cout << "Computed Entropy via cubature integration: " << value << " with error: " << err[0] << std::endl;
    std::cout << "Integration took " << exec_timer2 << " to complete." << std::endl;
    std::cout << "NOTICE: This is a statistical entropy, no thermodynamic transformation has been applied." << std::endl;

    writeToCSV("entropy.csv", "cubature_with_error_" + std::to_string(err[0]), value, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::HNIZDO, this->numberOfDraws);

    //for (std::size_t i = 0u; i < this->dimension; ++i)
    //{
    //  std::cout << "v " << val[i] << "  e " << err[i] << std::endl;
    //}

    delete[] xmin, xmax, err, val;
    return value * -1.;

  }

  void meanNNEntropyFaivishevsky()
  {
    std::cout << "Commencing meanNNEntropyFaivishevsky calculation." << std::endl;


    //http://papers.nips.cc/paper/3500-ica-based-on-a-smooth-estimation-of-the-differential-entropy
    transpose(drawMatrix);
    Matrix_Class drawMatrix_TemporaryCopy = drawMatrix;
    KahanAccumulation<double> kahan_acc_eucl_norm, kahan_acc_max_norm;
    double sum_final_eucl = 0.f;
    double sum_final_max = 0.f;
#ifdef _OPENMP
#pragma omp parallel firstprivate(drawMatrix_TemporaryCopy,kahan_acc_eucl_norm, kahan_acc_max_norm) reduction(+:sum_final_eucl) reduction(+:sum_final_max)
#endif
    {
#ifdef _OPENMP
      auto const n_omp = static_cast<std::ptrdiff_t>(drawMatrix_TemporaryCopy.cols());
#pragma omp for
      for (std::ptrdiff_t i = 0; i < n_omp; ++i)
#else
      for (size_t i = 0u; i < drawMatrix_TemporaryCopy.cols(); i++)
#endif
      {
        for (size_t j = i + 1u; j < drawMatrix_TemporaryCopy.cols(); j++)
        {
          // Eucledean norm
          double eucl_norm = 0.;
          double max_norm = std::numeric_limits<double>::max() * -1.;
          for (size_t currentDim = 0u; currentDim < this->dimension; currentDim++)
          {
            eucl_norm += (drawMatrix_TemporaryCopy(0u, i) - drawMatrix_TemporaryCopy(0u, j))
              * (drawMatrix_TemporaryCopy(0u, i) - drawMatrix_TemporaryCopy(0u, j));

            max_norm = std::max(max_norm, std::abs(drawMatrix_TemporaryCopy(0u, i) - drawMatrix_TemporaryCopy(0u, j)));

          }
          eucl_norm = sqrt(eucl_norm);
          kahan_acc_eucl_norm = KahanSum(kahan_acc_eucl_norm, log(eucl_norm));
          kahan_acc_max_norm = KahanSum(kahan_acc_max_norm, log(max_norm));

        }
      }
      sum_final_max += kahan_acc_max_norm.sum;
      sum_final_eucl += kahan_acc_eucl_norm.sum;
    }
    const double mod_summation_eucl = sum_final_eucl * static_cast<double>(this->dimension) / static_cast<double>(drawMatrix_TemporaryCopy.cols() * (drawMatrix_TemporaryCopy.cols() - 1));
    const double mod_summation_max = sum_final_max * static_cast<double>(this->dimension) / static_cast<double>(drawMatrix_TemporaryCopy.cols() * (drawMatrix_TemporaryCopy.cols() - 1));

    double sum_gamma_k = 0.f;
    const int n_of_samples = static_cast<int>(drawMatrix_TemporaryCopy.cols() - 1u);

#ifdef _OPENMP
#pragma omp parallel firstprivate(n_of_samples) reduction(+:sum_gamma_k)
#endif
    {
#ifdef _OPENMP
      auto const n_omp = static_cast<std::ptrdiff_t>(n_of_samples);
#pragma omp for
      for (std::ptrdiff_t i = 1; i < n_omp; ++i)
#else
      for (size_t i = 1u; i < n_of_samples; i++)
#endif
        sum_gamma_k += digammal(static_cast<long double>(i));
    }
    sum_gamma_k *= -1;
    sum_gamma_k *= static_cast<double>(this->dimension) / static_cast<double>(drawMatrix_TemporaryCopy.cols() * (drawMatrix_TemporaryCopy.cols() - 1));
    sum_gamma_k += digammal(static_cast<long double>(drawMatrix_TemporaryCopy.cols()));


    const double volume_of_unit_ball_for_eucledean_norm = pow(pi, static_cast<double>(this->dimension) / 2.) / tgamma(1. + static_cast<double>(this->dimension) / 2.);
    const double volume_of_unit_ball_for_max_norm = pow(2., this->dimension);


    transpose(drawMatrix);

    writeToCSV("entropy.csv", "faivishevsky", mod_summation_eucl + sum_gamma_k + log(volume_of_unit_ball_for_eucledean_norm), kNN_NORM::EUCLEDEAN, kNN_FUNCTION::HNIZDO, this->numberOfDraws);
    writeToCSV("entropy.csv", "faivishevsky", mod_summation_max + sum_gamma_k + log(volume_of_unit_ball_for_max_norm), kNN_NORM::MAXIMUM, kNN_FUNCTION::HNIZDO, this->numberOfDraws);

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
        temp[i] = drawMatrix(j,dimensionsToBeUsedInHistogramming[i]);
      }
      histograms_p->add_value(temp);
    }

    histograms_p->distribute();
    histograms_p->writeProbabilityDensity("entropytrails_" + filename);
    histograms_p->writeAuxilaryData("auxfile_entropytrails_" + filename);
    delete histograms_p;
  }

  // Full dimensional nearest neighbor entropy computation (without MI-Expansion)
  double calculateNN(const kNN_NORM norm, bool const& ardakaniCorrection , const kNN_FUNCTION func = kNN_FUNCTION::HNIZDO)
  {
    std::cout << "--------------------------\nCommencing Full-Dimensional NNEntropy calculation." << std::endl;

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
      for (size_t currentDraw = 0u; currentDraw < drawMatrix.rows(); currentDraw++)
      {
        for (unsigned int i = 0u; i < dimensionality; i++)
        {
          if (ardakaniCorrection_minimumValueInDataset.at(i) > drawMatrix(currentDraw,i))
            ardakaniCorrection_minimumValueInDataset.at(i) = drawMatrix(currentDraw, i);

          if (ardakaniCorrection_maximumValueInDataset.at(i) < drawMatrix(currentDraw, i))
            ardakaniCorrection_maximumValueInDataset.at(i) = drawMatrix(currentDraw, i);
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



      writeToCSV("entropy.csv", "fullNNentropy_ardakaniCorrected", ardakaniSum_hnizdo, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::HNIZDO, this->numberOfDraws);
      writeToCSV("entropy.csv", "fullNNentropy_ardakaniCorrected", ardakaniSum_goria, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::GORIA, this->numberOfDraws);
      writeToCSV("entropy.csv", "fullNNentropy_ardakaniCorrected", ardakaniSum_lombardi, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::LOMBARDI, this->numberOfDraws);


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



      writeToCSV("entropy.csv", "fullNNentropy_ardakaniCorrected", ardakaniSum_hnizdo, kNN_NORM::MAXIMUM, kNN_FUNCTION::HNIZDO, this->numberOfDraws);
      writeToCSV("entropy.csv", "fullNNentropy_ardakaniCorrected", ardakaniSum_goria, kNN_NORM::MAXIMUM, kNN_FUNCTION::GORIA, this->numberOfDraws);
      writeToCSV("entropy.csv", "fullNNentropy_ardakaniCorrected", ardakaniSum_lombardi, kNN_NORM::MAXIMUM, kNN_FUNCTION::LOMBARDI, this->numberOfDraws);


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



      writeToCSV("entropy.csv", "fullNNentropy", sum_hnizdo, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::HNIZDO, this->numberOfDraws);
      writeToCSV("entropy.csv", "fullNNentropy", sum_goria, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::GORIA, this->numberOfDraws);
      writeToCSV("entropy.csv", "fullNNentropy", sum_lombardi, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::LOMBARDI, this->numberOfDraws);


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



      writeToCSV("entropy.csv", "fullNNentropy", sum_hnizdo, kNN_NORM::MAXIMUM, kNN_FUNCTION::HNIZDO, this->numberOfDraws);
      writeToCSV("entropy.csv", "fullNNentropy", sum_goria, kNN_NORM::MAXIMUM, kNN_FUNCTION::GORIA, this->numberOfDraws);
      writeToCSV("entropy.csv", "fullNNentropy", sum_lombardi, kNN_NORM::MAXIMUM, kNN_FUNCTION::LOMBARDI, this->numberOfDraws);


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
      std::cout << ": " << std::setprecision(6) << returnValue << ".\n";
      std::cout << "NN Calculation took " << timer << " .\n";
      std::cout << "--------------------------\n" << std::endl;
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

    Matrix_Class cov_matr = Matrix_Class{ transposed(input) };
    cov_matr = cov_matr - Matrix_Class(input.cols(), input.cols(), 1.) * cov_matr / static_cast<float_type>(input.cols());
    cov_matr = transposed(cov_matr) * cov_matr;
    cov_matr *= (1.f / static_cast<float_type>(input.cols()));
    Matrix_Class eigenvalues;
    Matrix_Class eigenvectors;
    float_type cov_determ = 0.;
    int *cov_rank = new int;
    cov_matr.eigensym(eigenvalues, eigenvectors, cov_rank);


    //Remove Eigenvalues that should be zero if cov_matr is singular
    if ((*cov_rank < (int)eigenvalues.rows()) || (cov_determ = cov_matr.determ(), abs(cov_determ) < 10e-90))
    {
      std::cout << "Notice: covariance matrix is singular, attempting to fix by truncation of Eigenvalues.\n";
      std::cout << "Details: rank of covariance matrix is " << *cov_rank << ", determinant is " << cov_determ << ", size is " << cov_matr.rows() << ".\n";
      if (removeDOF)
      {
        size_t temp = std::max(6, int((cov_matr.rows() - *cov_rank)));
        eigenvalues.shed_rows(eigenvalues.rows() - temp, eigenvalues.rows() - 1u);
        eigenvectors.shed_cols(eigenvectors.cols() - temp, eigenvectors.cols() - 1u);
      }
      else
      {
        eigenvalues.shed_rows((*cov_rank), eigenvalues.rows() - 1u);
        eigenvectors.shed_cols((*cov_rank), eigenvectors.cols() - 1u);
      }
    }
    else if (removeDOF)
    {
      eigenvectors.shed_cols(0, 5);
      eigenvalues.shed_rows(0, 5);
    }
    delete cov_rank;

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

    writeToCSV("entropy.csv", "numata_no_corrections", entropy_sho, kNN_NORM::EUCLEDEAN, kNN_FUNCTION::HNIZDO, this->numberOfDraws);

    

    //Corrections for anharmonicity and M.I.
    // I. Create PCA-Modes matrix
    Matrix_Class eigenvectors_t(transposed(eigenvectorsPCA));

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
    const double temperatureInK,const kNN_NORM norm, const kNN_FUNCTION func, const bool removeNegativeMI = true, const float_type anharmonicityCutoff = 0.007)
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
      quantum_entropy(i, 0u) = ((alpha_i(i, 0u) / (exp(alpha_i(i, 0u)) - 1)) - log(1 - exp(-1 * alpha_i(i, 0u)))) ; // This was wrong (remove this comment only for github commit)
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
              std::cout << "Notice: PCA-Mode " << i << " not corrected for anharmonicity (estimated correction is < 0.0 : " << entropy_anharmonic(i, 0u)  << ").\n";
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
              std::cout << ": " << element.entropyValue * 1.380648813 * 6.02214129 * 0.239005736 <<  " cal / (mol * K)\n";
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

    writeToCSV("entropy.csv", "numata_without_corrections", entropy_sho, norm, func, this->numberOfDraws);
    writeToCSV("entropy.csv", "numata_with_anharmonicity_corrections", entropy_sho - delta_entropy, norm, func, this->numberOfDraws);

    std::string ident = "numata_with_generalized_corrections_order_" + std::to_string(orderOfCorrection);
    if (removeNegativeMI)
      ident += "secondOrderMIHardZero";
    writeToCSV("entropy.csv", ident, entropy_sho - delta_entropy - higher_order_entropy, norm, func,this->numberOfDraws);

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
    //transpose(drawMatrix); ??
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
      std::cout << "--- Order: " << order_N << "---" << "\n";
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
        std::cout << "Value: " << std::setprecision(6) << std::pow(-1., buffer.at(i).dim + 1.) * buffer.at(i).entropyValue << std::endl;
      }
      miEntropy += std::pow(-1., buffer.at(i).dim + 1.) * buffer.at(i).entropyValue;
    }
    if (Config::get().general.verbosity >= 4)
    {
      std::cout << "--- MI-EXPANSION INFO DONE ---" << "\n---------------------------------\n";
    }

    storeMI mi;
    mi.order = order_N;
    mi.value = miEntropy;
    this->calculatedMIEs.push_back(mi);

    if (Config::get().general.verbosity >= 3)
    {
      std::cout << "NN Calculation took " << timer << " ." << std::endl;
      std::cout << "NN Value: " << std::setprecision(6) << miEntropy << ".\n\n" << std::endl;

    }
    

    writeToCSV("entropy.csv", "MIE_order_" + std::to_string(order_N) + (ardakaniCorrection ? "_ardakanicorrected" : ""), miEntropy, norm, func, this->numberOfDraws);


    //transpose(drawMatrix);
    return miEntropy;

  }

private:

  void writeToCSV(std::string nameToFile, std::string ident, const float_type & value, kNN_NORM norm, kNN_FUNCTION func, size_t numberOfDraws)
  {
    std::ofstream myfile2;
    myfile2.open(std::string(nameToFile), std::ios::app);

    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << kNN << ",";
    if (norm == kNN_NORM::EUCLEDEAN)
      myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << "eucledean,";
    else
      myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << "maxnorm,";

    if (func == kNN_FUNCTION::GORIA)
      myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << "goria,";
    else if (func == kNN_FUNCTION::LOMBARDI)
      myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << "lombardi,";
    else
      myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << "hnizdo,";

    myfile2 << std::setw(25) << std::scientific << std::setprecision(15) << ident << ",";
    myfile2 << std::setw(25) << std::scientific << std::setprecision(15) << value << ",";
    myfile2 << std::setw(25) << std::scientific << std::setprecision(15) << numberOfDraws << "\n";

    myfile2.close();
  }

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

      Matrix_Class copytemp = transposed(drawMatrix); // rows are dimensions
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

      //DEBUG
      //auto vis_this = copytemp.to_std_vector();
      //std::cout << std::endl << vis_this[0][0] << " " << vis_this[0][1] << " " << vis_this[0][2] << std::endl;
      //std::cout << " " << vis_this[1][0] << " " << vis_this[2][0] << std::endl;

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

