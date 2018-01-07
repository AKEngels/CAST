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
          mean.at(0)(i, 0) = std::stod(s);
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
            mean.at(j)(i, 0) = std::stod(s);
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
              covariance.at(i)(j, k) = std::stod(line);
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




// Class for PDF
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
      PDF = [&, this](std::vector<double> const& x, std::vector<unsigned int> const& subdims = std::vector<unsigned int>())
      {
        if (x.size() != 1)
          throw std::runtime_error("Wring dimensionality for chosen Probability Density.");
        return (1. / sqrt(2. * pi * 1. * 1.)) * exp(-1.*(x.at(0)*x.at(0)) / double(2. * 1. * 1.));
      };
      //PDFrange = 5.920; // range incorporating function values up to under 10e-7
      this->m_identString = "Normal sig=1";
    }
    else if (ident_ == 1)
    {
      this->dimension = 1;
      // Gaussian with sigma 10
      PDF = [&, this](std::vector<double> const& x, std::vector<unsigned int> const& subdims = std::vector<unsigned int>())
      {
        if (x.size() != 1)
          throw std::runtime_error("Wring dimensionality for chosen Probability Density.");
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
      PDF = [&, this](std::vector<double> const& x, std::vector<unsigned int> const& subdims = std::vector<unsigned int>())
      {
        if (x.size() != 1)
          throw std::runtime_error("Wring dimensionality for chosen Probability Density.");
        if (x.at(0) < -1. || x.at(0) > 1.)
          return 0.;
        else
          return (3. / 4.) * (-1. * ((x.at(0))*(x.at(0))) + 1.);
      };
      maximumOfPDF = PDF(std::vector<double>{ 0. }, std::vector<unsigned int>());
      analyticEntropy_ = 0.568054;
      PDFrange = std::make_shared<std::pair<double, double>>(-1., 1.); // Function is only defined in [-1;1]
      this->m_identString = "parabolic";
    }
    else if (ident_ == 3)
    {
      this->dimension = 1;
      // Beta Distribution https://en.wikipedia.org/wiki/Differential_entropy
      PDF = [&, this](std::vector<double> const& x, std::vector<unsigned int> const& subdims = std::vector<unsigned int>())
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
    else if (ident_ == 4)
    {
      this->dimension = 3;
      // Multivariate gaussian
      PDF = [&, this](std::vector<double> const& x, std::vector<unsigned int> const& subdims = std::vector<unsigned int>())
      {
        if (x.size() != 3)
          throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");

        Matrix_Class input(x.size(), 1u);
        for (unsigned int i = 0u; i < x.size(); i++)
          input(i, 0u) = x.at(i);

        Matrix_Class mean(x.size(), 1u, 0.);
        Matrix_Class covariance(x.size(), x.size(), 0u);
        covariance(0, 0) = 1.;
        covariance(1, 0) = 0.7;
        covariance(2, 1) = 0.9;
        covariance(0, 1) = 0.7;
        covariance(1, 1) = 1.5;
        covariance(1, 2) = 0.9;
        covariance(2, 2) = 1.5;
        covariance(0, 2) = 0.56;
        covariance(2, 0) = 0.56;

        Matrix_Class value = transposed(Matrix_Class(input - mean)) * covariance.inversed() * Matrix_Class(input - mean);
        const double prefactor = 1. / std::sqrt(std::pow(2 * ::constants::pi, x.size()) * covariance.determ());


        return prefactor * std::exp(value(0u,0u) * -0.5);
      };
      maximumOfPDF = PDF(std::vector<double>{ 0.,0.,0. }, std::vector<unsigned int>());
      Matrix_Class covariance(this->dimension, this->dimension, 0u);
      covariance(0, 0) = 1.;
      covariance(1, 0) = 0.7;
      covariance(2, 1) = 0.9;
      covariance(0, 1) = 0.7;
      covariance(1, 1) = 1.5;
      covariance(1, 2) = 0.9;
      covariance(2, 2) = 1.5;
      covariance(0, 2) = 0.56;
      covariance(2, 0) = 0.56;

      analyticEntropy_ = 0.5 * log(covariance.determ()) + double(this->dimension)/2. * std::log(2. * ::constants::pi * ::constants::e);
      PDFrange = std::make_shared<std::pair<double, double>>(-8, 8.);
      this->m_identString = "3varGauss";
    }
    else if (ident_ == 5)
    {
      this->dimension = 2;
      // Multivariate gaussian
      PDF = [&, this](std::vector<double> const& x, std::vector<unsigned int> const& subdims = std::vector<unsigned int>())
      {
        if (x.size() != 2)
          throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");
        const double coeff = 0.9;

        const double exp_argument = -0.5 * (1. / (1. - coeff * coeff)) * (x.at(0) * x.at(0) - 2. * coeff * x.at(0)* x.at(1) + x.at(1) * x.at(1));
        const double prefactor = 1. / ((2 * ::constants::pi) * sqrt(1 -  coeff *  coeff));


        return prefactor * std::exp(exp_argument);
      };
      maximumOfPDF = PDF(std::vector<double>{ 0., 0. }, std::vector<unsigned int>());
      const double coeff = 0.9;
      Matrix_Class covariance(2, 2, 0u);
      covariance(0, 0) = 1.;
      covariance(1, 0) = coeff;
      covariance(1, 1) = 1.;
      covariance(0, 1) = coeff;
      analyticEntropy_ = 0.5 * log(covariance.determ()) + double(this->dimension) / 2. * std::log(2. * ::constants::pi * ::constants::e);
      PDFrange = std::make_shared<std::pair<double, double>>(-8, 8.);
      this->m_identString = "2varGauss";
    }
    else if (ident_ == 6)
    {
      // 7 GMM fittet to MD of tridecalanine dim 1,2
      this->dimension = 2;
      // Multivariate gaussian
      PDF = [&, this](std::vector<double> const& x, std::vector<unsigned int> const& subdims = std::vector<unsigned int>())
      {
        if (x.size() != 2)
          throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");

        Matrix_Class input(x.size(), 1u);
        for (unsigned int i = 0u; i < x.size(); i++)
          input(i, 0u) = x.at(i);

        std::vector<double> weight(7u,0.);
        std::vector<Matrix_Class> mean(7u, Matrix_Class(x.size(), 1u, 0.));
        std::vector<Matrix_Class> covariance(7u, Matrix_Class(2, 2, 0u));

        weight[0] = 4.59166833e-01;
        weight[1] = 1.81156684e-01;
        weight[2] = 4.40929007e-02;
        weight[3] = 1.11308104e-01;
        weight[4] = 9.41269069e-02;
        weight[5] = 8.58282250e-02;
        weight[6] = 2.43203471e-02;

        mean[0](0, 0) = 8.56364042e-12;  mean[0](1, 0) = 3.34079765e-13;
        mean[1](0, 0) = -2.96791360e-12; mean[1](1, 0) = -2.69003397e-12;
        mean[2](0, 0) = 8.47527692e-12;  mean[2](1, 0) = 9.86326412e-12;
        mean[3](0, 0) = 1.26754874e-11;  mean[3](1, 0) = -4.38271175e-12;
        mean[4](0, 0) = 1.05388512e-12;  mean[4](1, 0) = 4.78499853e-12;
        mean[5](0, 0) = 4.78838598e-12;  mean[5](1, 0) = -6.31405122e-12;
        mean[6](0, 0) = 3.82288827e-12;  mean[6](1, 0) = 1.34958334e-11;

        covariance[0](0, 0) = 1.46117966e-23;
        covariance[0](1, 0) = -9.62323794e-26;
        covariance[0](1, 1) = 5.74205745e-24;
        covariance[0](0, 1) = covariance[0](1, 0);
        covariance[1](0, 0) = 5.56924448e-24;
        covariance[1](1, 0) = 8.28644834e-25;
        covariance[1](1, 1) = 4.53278335e-24;
        covariance[1](0, 1) = covariance[1](1, 0);
        covariance[2](0, 0) = 1.20196030e-23;
        covariance[2](1, 0) = -4.14350002e-24;
        covariance[2](1, 1) = 8.64760265e-24;
        covariance[2](0, 1) = covariance[2](1, 0);
        covariance[3](0, 0) = 4.04029907e-24;
        covariance[3](1, 0) = 4.25933021e-25;
        covariance[3](1, 1) = 1.20733036e-23;
        covariance[3](0, 1) = covariance[3](1, 0);
        covariance[4](0, 0) = 5.61578290e-24;
        covariance[4](1, 0) = -5.66515611e-25;
        covariance[4](1, 1) = 7.56044583e-24;
        covariance[4](0, 1) = covariance[4](1, 0);
        covariance[5](0, 0) = 7.97900136e-24;
        covariance[5](1, 0) = -1.89541174e-24;
        covariance[5](1, 1) = 7.87990277e-24;
        covariance[5](0, 1) = covariance[5](1, 0);
        covariance[6](0, 0) = 4.20157802e-24;
        covariance[6](1, 0) = -8.10342841e-25;
        covariance[6](1, 1) = 1.80154978e-24;
        covariance[6](0, 1) = covariance[6](1, 0);

        double returnValue = 0.;
        for (unsigned int i = 0u; i < 7; i++)
        {
          Matrix_Class value = transposed(Matrix_Class(input - mean[i])) * covariance[i].inversed() * Matrix_Class(input - mean[i]);
          const double prefactor = 1. / std::sqrt(std::pow(2 * ::constants::pi, x.size()) * covariance[i].determ());
          returnValue += prefactor * std::exp(value(0u, 0u) * -0.5) * weight[i];
        }

        return returnValue;
      };
      maximumOfPDF = PDF(std::vector<double>{ 0., 0. }, std::vector<unsigned int>()) + 0.7; //I DONT KNOW IF THIS IS OK, CHECK
      analyticEntropy_ = std::numeric_limits<double>::quiet_NaN();
      PDFrange = std::make_shared<std::pair<double, double>>(-1e-10, 1e-10);
      this->m_identString = "GMM1";
    }
    else if (ident_ == 7)
    {
      // 10 GMM fittet to MD of tridecalanine dim 1,2,3
      this->dimension = 3;
      // Multivariate gaussian
      PDF = [&, this](std::vector<double> const& x, std::vector<unsigned int> const& subdims = std::vector<unsigned int>())
      {
        if (x.size() != 3)
          throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");

        Matrix_Class input(x.size(), 1u);
        for (unsigned int i = 0u; i < x.size(); i++)
          input(i, 0u) = x.at(i);

        std::vector<double> weight(10u, 0.);
        std::vector<Matrix_Class> mean(10u, Matrix_Class(x.size(), 1u, 0.));
        std::vector<Matrix_Class> covariance(10u, Matrix_Class(3, 3, 0u));

        weight[0] = 0.297859;
        weight[1] = 0.0322744;
        weight[2] = 0.147175;
        weight[3] = 0.0242578;
        weight[4] = 0.087544;
        weight[5] = 0.0482769;
        weight[6] = 0.154248;
        weight[7] = 0.0822646;
        weight[8] = 0.105497;
        weight[9] = 0.0206026;

        mean[0](0, 0) = 1.01991e-11; mean[0](1, 0) = 6.2136e-13 ; mean[0](2, 0)   = -2.22955e-12    ;
        mean[1](0, 0) = 2.02376e-12; mean[1](1, 0) = 5.60621e-12; mean[1](2, 0) = -4.73088e-12    ;
        mean[2](0, 0) = 3.67781e-12; mean[2](1, 0) = -1.86976e-12; mean[2](2, 0) = -1.83584e-12   ;
        mean[3](0, 0) = 3.9285e-12; mean[3](1, 0) = 1.25377e-11; mean[3](2, 0) = -8.63069e-12     ;
        mean[4](0, 0) = 8.38703e-12; mean[4](1, 0) = -6.51049e-12; mean[4](2, 0) = -6.30375e-12   ;
        mean[5](0, 0) = 1.0033e-11; mean[5](1, 0) = 8.74454e-12; mean[5](2, 0) = -5.99058e-12     ;
        mean[6](0, 0) = -3.60727e-12; mean[6](1, 0) = -2.87797e-12; mean[6](2, 0) = -4.07751e-12  ;
        mean[7](0, 0) = 1.20153e-11; mean[7](1, 0) = -1.8596e-12; mean[7](2, 0) = 4.33388e-14     ;
        mean[8](0, 0) = 1.37377e-12; mean[8](1, 0) = 3.74631e-12; mean[8](2, 0) = -1.55544e-12    ;
        mean[9](0, 0) = 1.19832e-11; mean[9](1, 0) = -9.49643e-12; mean[9](2, 0) = -3.73347e-13   ;

        covariance[0](0, 0) = 9.03379e-24 ;
        covariance[0](1, 0) = -1.59016e-24;
        covariance[0](2, 0) = -3.71382e-24;
        covariance[0](1, 1) = 5.02031e-24 ;
        covariance[0](2, 1) = -5.47595e-25;
        covariance[0](2, 2) = 1.39e-23    ;
        covariance[0](0, 1) = covariance[0](1, 0);
        covariance[0](0, 2) = covariance[0](2, 0);
        covariance[0](1, 2) = covariance[0](2, 1);

        covariance[1](0, 0) = 4.72767e-24  ;
        covariance[1](1, 0) = 3.64308e-25;
        covariance[1](2, 0) = 4.1357e-25 ;
        covariance[1](1, 1) = 2.70482e-23;
        covariance[1](2, 1) = 6.36459e-24;
        covariance[1](2, 2) = 7.4303e-24 ;
        covariance[1](0, 1) = covariance[1](1, 0);
        covariance[1](0, 2) = covariance[1](2, 0);
        covariance[1](1, 2) = covariance[1](2, 1);


        covariance[2](0, 0) = 4.81709e-24  ;
        covariance[2](1, 0) = 3.40596e-25  ;
        covariance[2](2, 0) = -8.31425e-25 ;
        covariance[2](1, 1) = 6.35518e-24  ;
        covariance[2](2, 1) = 2.17662e-24  ;
        covariance[2](2, 2) = 1.1871e-23   ;
        covariance[2](0, 1) = covariance[2](1, 0);
        covariance[2](0, 2) = covariance[2](2, 0);
        covariance[2](1, 2) = covariance[2](2, 1);

        covariance[3](0, 0) = 2.67491e-24;
        covariance[3](1, 0) = -1.73085e-24;
        covariance[3](2, 0) = -1.32122e-24;
        covariance[3](1, 1) = 4.93419e-24;
        covariance[3](2, 1) = 3.15282e-24;
        covariance[3](2, 2) = 3.7244e-24;
        covariance[3](0, 1) = covariance[3](1, 0);
        covariance[3](0, 2) = covariance[3](2, 0);
        covariance[3](1, 2) = covariance[3](2, 1);

        covariance[4](0, 0) = 1.68382e-23;
        covariance[4](1, 0) = 4.61495e-24;
        covariance[4](2, 0) = -5.11702e-24;
        covariance[4](1, 1) = 6.54348e-24;
        covariance[4](2, 1) = -3.82479e-24;
        covariance[4](2, 2) = 1.5579e-23;
        covariance[4](0, 1) = covariance[4](1, 0);
        covariance[4](0, 2) = covariance[4](2, 0);
        covariance[4](1, 2) = covariance[4](2, 1);

        covariance[5](0, 0) = 9.38482e-24;
        covariance[5](1, 0) = -7.38849e-24;
        covariance[5](2, 0) = -3.31397e-24;
        covariance[5](1, 1) = 1.51979e-23;
        covariance[5](2, 1) = 2.51153e-24;
        covariance[5](2, 2) = 4.87609e-24;
        covariance[5](0, 1) = covariance[5](1, 0);
        covariance[5](0, 2) = covariance[5](2, 0);
        covariance[5](1, 2) = covariance[5](2, 1);

        covariance[6](0, 0) = 3.33945e-24;
        covariance[6](1, 0) = 1.59445e-25;
        covariance[6](2, 0) = -8.61476e-25;
        covariance[6](1, 1) = 3.94627e-24;
        covariance[6](2, 1) = 4.21489e-24;
        covariance[6](2, 2) = 9.25474e-24;
        covariance[6](0, 1) = covariance[6](1, 0);
        covariance[6](0, 2) = covariance[6](2, 0);
        covariance[6](1, 2) = covariance[6](2, 1);

        covariance[7](0, 0) = 3.83925e-24;
        covariance[7](1, 0) = -2.37765e-24;
        covariance[7](2, 0) = -1.75374e-24;
        covariance[7](1, 1) = 5.72715e-24;
        covariance[7](2, 1) = 3.57909e-24;
        covariance[7](2, 2) = 4.99204e-24;
        covariance[7](0, 1) = covariance[7](1, 0);
        covariance[7](0, 2) = covariance[7](2, 0);
        covariance[7](1, 2) = covariance[7](2, 1);

        covariance[8](0, 0) = 6.9919e-24;
        covariance[8](1, 0) = -1.45301e-24;
        covariance[8](2, 0) = -3.2444e-24;
        covariance[8](1, 1) = 7.4847e-24;
        covariance[8](2, 1) = 2.61163e-24;
        covariance[8](2, 2) = 9.62358e-24;
        covariance[8](0, 1) = covariance[8](1, 0);
        covariance[8](0, 2) = covariance[8](2, 0);
        covariance[8](1, 2) = covariance[8](2, 1);

        covariance[9](0, 0) = 5.27517e-24;
        covariance[9](1, 0) = 4.34838e-25;
        covariance[9](2, 0) = -2.87762e-24;
        covariance[9](1, 1) = 2.49859e-24;
        covariance[9](2, 1) = 2.80166e-24;
        covariance[9](2, 2) = 9.32401e-24;
        covariance[9](0, 1) = covariance[9](1, 0);
        covariance[9](0, 2) = covariance[9](2, 0);
        covariance[9](1, 2) = covariance[9](2, 1);


        double returnValue = 0.;
        for (unsigned int i = 0u; i < 10; i++)
        {
          Matrix_Class value = transposed(Matrix_Class(input - mean[i])) * covariance[i].inversed() * Matrix_Class(input - mean[i]);
          const double prefactor = 1. / std::sqrt(std::pow(2 * ::constants::pi, x.size()) * covariance[i].determ());
          returnValue += prefactor * std::exp(value(0u, 0u) * -0.5) * weight[i];
        }

        return returnValue;
      };
      maximumOfPDF = PDF(std::vector<double>{ 0., 0., 0. }, std::vector<unsigned int>()) + 0.7; //I DONT KNOW IF THIS IS OK, CHECK
      analyticEntropy_ = std::numeric_limits<double>::quiet_NaN();
      PDFrange = std::make_shared<std::pair<double, double>>(-4e-11, 4e-11);
      this->m_identString = "GMM1";
    }
    else if (ident_ == 8)
    {
      // 10 GMM fittet to MD of tridecalanine dim 1,2,3
      this->dimension =2;
      // Multivariate gaussian
      PDF = [&, this](std::vector<double> const& x, std::vector<unsigned int> const& subdims = std::vector<unsigned int>())
      {
        if (x.size() !=2)
          throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");

        Matrix_Class input(x.size(), 1u);
        for (unsigned int i = 0u; i < x.size(); i++)
          input(i, 0u) = x.at(i);

        std::vector<double> weight(10u, 0.);
        std::vector<Matrix_Class> mean(10u, Matrix_Class(x.size(), 1u, 0.));
        std::vector<Matrix_Class> covariance(10u, Matrix_Class(2, 2, 0u));

        weight[0] = 0.296898      ;
        weight[1] = 0.11709     ;
        weight[2] = 0.00502877  ;
        weight[3] = 0.186937    ;
        weight[4] = 0.0422829   ;
        weight[5] = 0.0138878   ;
        weight[6] = 0.000274472 ;
        weight[7] = 0.0144979   ;
        weight[8] = 0.286551    ;
        weight[9] = 0.0365531   ;

        mean[0](0, 0) = 4.78981e-12 ; mean[0](1, 0) = -1.36579e-12   ;
        mean[1](0, 0) = 5.51715e-12; mean[1](1, 0) = 4.76378e-12   ;
        mean[2](0, 0) = 4.21048e-12; mean[2](1, 0) = -1.14226e-11  ;
        mean[3](0, 0) = -2.85438e-12; mean[3](1, 0) = -1.61024e-12 ;
        mean[4](0, 0) = 9.96989e-12; mean[4](1, 0) = -8.20098e-12  ;
        mean[5](0, 0) = 3.36488e-12; mean[5](1, 0) = 1.33293e-11   ;
        mean[6](0, 0) = 4.71386e-1;   mean[6](1, 0) = 1.65182e-11  ;
        mean[7](0, 0) = 1.44319e-11; mean[7](1, 0) = -7.62714e-12  ;
        mean[8](0, 0) = 1.15107e-11; mean[8](1, 0) = -4.21604e-13  ;
        mean[9](0, 0) = 6.82143e-12; mean[9](1, 0) = 1.21607e-11   ;

        covariance[0](0, 0) = 1.11731e-23  ;
          covariance[0](1, 0) = 2.00613e-24;
          covariance[0](1, 1) = 1.15219e-23;
        covariance[0](0, 1) = covariance[0](1, 0);

        covariance[1](0, 0) = 2.47926e-23;
          covariance[1](1, 0) = -4.84533e-24 ;
          covariance[1](1, 1) = 8.67598e-24  ;
        covariance[1](0, 1) = covariance[1](1, 0);

        covariance[2](0, 0) = 5.57373e-24     ;
          covariance[2](1, 0) = -1.0595e-24   ;
          covariance[2](1, 1) = 8.4526e-25    ;
        covariance[2](0, 1) = covariance[2](1, 0);

        covariance[3](0, 0) = 5.61383e-24;
        covariance[3](1, 0) = 3.79572e-24;
        covariance[3](1, 1) = 9.50937e-24;
        covariance[3](0, 1) = covariance[3](1, 0);

        covariance[4](0, 0) = 7.01799e-24;
        covariance[4](1, 0) = -1.14537e-25;
        covariance[4](1, 1) = 3.25432e-24;
        covariance[4](0, 1) = covariance[4](1, 0);

        covariance[5](0, 0) = 3.02747e-24;
        covariance[5](1, 0) = -7.55053e-25;
        covariance[5](1, 1) = 2.16706e-24;
        covariance[5](0, 1) = covariance[5](1, 0);

        covariance[6](0, 0) = 5.70972e-25;
        covariance[6](1, 0) = -1.83951e-25;
        covariance[6](1, 1) = 2.22965e-25;
        covariance[6](0, 1) = covariance[6](1, 0);

        covariance[7](0, 0) = 1.90793e-24;
        covariance[7](1, 0) = 1.96845e-24;
        covariance[7](1, 1) = 5.79567e-24;
        covariance[7](0, 1) = covariance[7](1, 0);

        covariance[8](0, 0) = 5.49353e-24;
        covariance[8](1, 0) = -1.62023e-24;
        covariance[8](1, 1) = 6.07475e-24;
        covariance[8](0, 1) = covariance[8](1, 0);

        covariance[9](0, 0) = 1.31776e-23;
        covariance[9](1, 0) = -3.69162e-24;
        covariance[9](1, 1) = 4.42914e-24;
        covariance[9](0, 1) = covariance[9](1, 0);

        double returnValue = 0.;
        for (unsigned int i = 0u; i < 10; i++)
        {
          Matrix_Class value = transposed(Matrix_Class(input - mean[i])) * covariance[i].inversed() * Matrix_Class(input - mean[i]);
          const double prefactor = 1. / std::sqrt(std::pow(2 * ::constants::pi, x.size()) * covariance[i].determ());
          returnValue += prefactor * std::exp(value(0u, 0u) * -0.5) * weight[i];
        }

        return returnValue;
      };
      PDFrange = std::make_shared<std::pair<double, double>>(-9e-11, 9e-11);



      maximumOfPDF = PDF(std::vector<double>{ 0., 0. }, std::vector<unsigned int>()) + 0.7; //I DONT KNOW IF THIS IS OK, CHECK
      analyticEntropy_ = std::numeric_limits<double>::quiet_NaN();

      this->m_identString = "GMM3";
    }
    else if (ident_ == 9)
    {
      // Read GMM from file
      gmmdatafromfile = GMM_data("gmmfile.dat");
      std::cout << "Read gmmfile.dat with " << gmmdatafromfile.numberOfDimensions << " dimensions and " << gmmdatafromfile.numberOfGaussians << " gaussians." << std::endl;
      this->dimension = gmmdatafromfile.numberOfDimensions;
      // Multivariate gaussian
      PDF = [&, this](std::vector<double> const& x, std::vector<unsigned int> const& subdims = std::vector<unsigned int>())
      {
        if (subdims == std::vector<unsigned int>())
        {
          if (x.size() != this->dimension)
            throw std::runtime_error("Wrong dimensionality for chosen Probability Density.");

          Matrix_Class input(x.size(), 1u);
          for (unsigned int i = 0u; i < x.size(); i++)
            input(i, 0u) = x.at(i);

          unsigned int const numberOfGaussians = gmmdatafromfile.numberOfGaussians;

          std::vector<double> weight = gmmdatafromfile.weights;
          std::vector<Matrix_Class> mean = gmmdatafromfile.mean;
          std::vector<Matrix_Class> covariance = gmmdatafromfile.covariance;

          long double returnValue = 0.;
          for (unsigned int i = 0u; i < numberOfGaussians; i++)
          {
            Matrix_Class value = transposed(Matrix_Class(input - mean[i])) * covariance[i].inversed() * Matrix_Class(input - mean[i]);
            const long double cov_determ = covariance[i].determ();
            const long double prefactor = 1. / std::sqrt(std::pow(2 * ::constants::pi, x.size()) * cov_determ);

            const double addition = prefactor * std::exp(value(0u, 0u) * -0.5) * weight[i];
            if (addition != addition)
              continue;
            returnValue += addition;
          }

          return returnValue;

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



          long double returnValue = 0.;
          for (unsigned int i = 0u; i < numberOfGaussians; i++)
          {
            Matrix_Class value = transposed(Matrix_Class(input - mean[i])) * covariance[i].inversed() * Matrix_Class(input - mean[i]);
            const long double cov_determ = covariance[i].determ();
            const long double prefactor = 1. / std::sqrt(std::pow(2 * ::constants::pi, x.size()) * cov_determ);

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
        PDFrange = std::make_shared<std::pair<double, double>>(-2e-11, 4e-11);
      }
      else
      {
        PDFrange = std::make_shared<std::pair<double, double>>(Config::get().entropytrails.rangeForGMM);
      }

      std::cout << "Range for GMM probability density is from " << PDFrange->first << " to " << PDFrange->second << "." << std::endl;


      //
      //find maximum of pdf by drawing many times
      std::random_device rd;
      std::mt19937 gen;
      std::cout << "Finding maximum of PDF by drawing 100k times, taking highest value x10." << std::endl;
      std::uniform_real_distribution<double> unifDistrRange(PDFrange->first, PDFrange->second);
      const double absrange = PDFrange->second - PDFrange->first;

      for (unsigned int i = 0u; i < 100000u; i++)
      {
        std::vector<double> drawUnifWithRange;
        for (unsigned int currentDim = 0u; currentDim < this->dimension; currentDim++)
        {
          drawUnifWithRange.push_back(unifDistrRange(gen));
        }
        const double pdfvalue = PDF(drawUnifWithRange, std::vector<unsigned int>());
        maximumOfPDF = std::max(pdfvalue, maximumOfPDF);
      }
      maximumOfPDF *= 10.;
      std::cout << "Estimated maximum of current PDF is " << maximumOfPDF << "." << std::endl;
      
      //
      //

      //maximumOfPDF = PDF(std::vector<double>(this->dimension, 0.), std::vector<unsigned int>()) * 1.5; //I DONT KNOW IF THIS IS OK, CHECK
      analyticEntropy_ = std::numeric_limits<double>::quiet_NaN();
      this->m_identString = "GMMfile";
    }
  };

  std::vector<std::vector<double>> draw(size_t numberOfSamples, std::vector<unsigned int> const& subdims)
  {
    std::random_device rd;
    std::mt19937 gen;
    //
    std::vector<std::vector<double>> draws;
    std::uniform_real_distribution<double> unifDistr(0, 1);

    //Target PDF in 0.1er Schritten auswerten und schauen wann es kleiner ist
    if (PDFrange == nullptr)
      PDFrange = std::make_shared<std::pair<double, double>>(meaningfulRange());

    const unsigned int currentDimensionality = subdims.size() == 0 ? this->dimension : subdims.size();


    std::uniform_real_distribution<double> unifDistrRange(PDFrange->first, PDFrange->second);
    const long double absrange = PDFrange->second - PDFrange->first;
    const long double unifWithRangeProbability = 1. / absrange;
    const long double k = maximumOfPDF *1.05 / unifWithRangeProbability;


    for (unsigned int n = 0; n < numberOfSamples; ++n)
    {
      std::vector<double> drawUnifWithRange;
      //double drawUnif = unifDistr(gen);
      //double drawUnifWithRange = unifDistrRange(gen);


      for (unsigned int currentDim = 0u; currentDim < currentDimensionality; currentDim++)
      {
        drawUnifWithRange.push_back(unifDistrRange(gen));
      }
      const long double drawUnif = unifDistr(gen);
      const long double pdfvalue = PDF(drawUnifWithRange, subdims);
      const long double p = pdfvalue / (k * unifWithRangeProbability);
      bool inRange = drawUnif < p;

      if (inRange)
      {
        draws.push_back(drawUnifWithRange);

        if (draws.size() % 25 == 0)
          std::cout << "Number of draws: " << draws.size() << std::endl;
      }
      else
      {
        --n;
      }
    }
    return draws;
  };

  int ident() const { return this->identif; };

  std::string identString() const
  {
    return this->m_identString;
  }

  double analyticEntropy() {
    return this->analyticEntropy_;
  };

  auto function() -> std::function<double(std::vector<double> const& x, std::vector<unsigned int> const& subdims)>
  {
    return this->PDF;
  };

  unsigned int getDimension() const
  {
    return this->dimension;
  }

  std::pair<double, double> meaningfulRange()
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
        if (PDF(std::vector<double> {newPDFrange}, std::vector<unsigned int>()) < 0.00000001)
          break;
        else
          newPDFrange += 0.001;
      }
      rangeOne = newPDFrange;
      newPDFrange = 0.;
      while (run)
      {
        if (PDF(std::vector<double> {newPDFrange}, std::vector<unsigned int>()) < 0.000000001)
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
  int dimension = 0;
  std::shared_ptr<std::pair<double, double>> PDFrange = nullptr;
  std::string m_identString = "ERROR_NAME";
  std::function<double(std::vector<double> const& x, std::vector<unsigned int> const& subdims)> PDF;
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
  int identifierOfPDF;
  ProbabilityDensity probdens;
  std::vector<unsigned int> subDimsGMM;
  entropyobj(size_t iter_, ProbabilityDensity probdens_, std::vector<unsigned int> const& subdimsGMM_) :
    numberOfDraws(iter_), dimension(probdens_.getDimension()), drawMatrix(iter_, probdens_.getDimension()),
    probdens(probdens_), identifierOfPDF(probdens_.ident()),
    subDimsGMM(subdimsGMM_)
  {
    // Check if draw exists
    std::ifstream myfile;
    std::string line;
    std::string filename;
    filename = std::string("draw_i" + std::to_string(numberOfDraws) + "_d" + std::to_string(dimension) + "_ident" + std::to_string(identifierOfPDF) + ".txt");
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
      std::vector<std::vector<double>> draws = probdens.draw(numberOfDraws, subDimsGMM);

      //Sort samples (entirely optional) if dimension == 1
      if (this->dimension == 1)
        std::stable_sort(draws.begin(), draws.end());

      for (unsigned int n = 0; n < numberOfDraws; ++n)
        for (unsigned int currentDim = 0; currentDim < this->dimension; ++currentDim)
          drawMatrix(n, currentDim) = draws[n][currentDim];

      // Write 
      std::ofstream myfile2;
      myfile2.open(std::string("draw_i" + std::to_string(numberOfDraws) + "_d" + std::to_string(dimension) + "_ident" + std::to_string(identifierOfPDF) + ".txt"));
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
};

struct informationForCubatureIntegration
{
  ProbabilityDensity& probdens;
  std::vector<unsigned int>& subdims;
  informationForCubatureIntegration(ProbabilityDensity& probdens_in, std::vector<unsigned int>& subdims_in)
    : probdens(probdens_in), subdims(subdims_in) {}
};

int cubaturefunctionEntropy(unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval)
{
  informationForCubatureIntegration& info = *((informationForCubatureIntegration *)fdata);

  if (ndim != info.probdens.getDimension())
    throw;
  if (fdim != 1)
    throw;

  for (int j = 0; j < npts; ++j) 
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

  for (int j = 0; j < npts; ++j)
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

// Calculated entropy object calculates estiamted entropy
class calculatedentropyobj : public entropyobj
{
public:
  size_t k;
  double mean, standardDeviation;
  double calculatedEntropyHnizdo;
  double calculatedEntropyLombardi;
  double calculatedEntropyMeanFaivishevskyMaximum;
  double calculatedEntropyMeanFaivishevskyEucledean;
  double mcintegrationEntropy;
  double mcdrawEntropy;
  double analyticalEntropy;
  double empiricalNormalDistributionEntropy;
  double calculatedEntropyGoria; // http://www.tandfonline.com/doi/abs/10.1080/104852504200026815
  double ardakaniEntropyEucledean;
  double ardakaniEntropyMaximum;
  double calculatedEntropyGoriaMaximum;
  double calculatedEntropyLombardiMaximum;
  double cubatureIntegral;
  calculatedentropyobj(size_t k_, entropyobj const& obj) :
    entropyobj(obj),
    k(k_),
    mean(std::numeric_limits<double>::quiet_NaN()),
    standardDeviation(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyHnizdo(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyLombardi(std::numeric_limits<double>::quiet_NaN()),
    mcintegrationEntropy(std::numeric_limits<double>::quiet_NaN()),
    mcdrawEntropy(std::numeric_limits<double>::quiet_NaN()),
    analyticalEntropy(this->probdens.analyticEntropy()),
    empiricalNormalDistributionEntropy(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyMeanFaivishevskyEucledean(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyMeanFaivishevskyMaximum(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyGoria(std::numeric_limits<double>::quiet_NaN()),
    ardakaniEntropyEucledean(std::numeric_limits<double>::quiet_NaN()),
    ardakaniEntropyMaximum(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyGoriaMaximum(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyLombardiMaximum(std::numeric_limits<double>::quiet_NaN()),
    cubatureIntegral(std::numeric_limits<double>::quiet_NaN())
  {

  }

  double empiricalGaussianEntropy()
  {
    std::cout << "Commencing empirical gaussian entropy calculation." << std::endl;

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

      return log(standardDeviation * sqrt(2. * pi * e));
    }
    else
    {
      Matrix_Class cov_matr;
      cov_matr = Matrix_Class(this->dimension, this->dimension, 0.);
      Matrix_Class meanPerDim(this->dimension, 1u, std::numeric_limits<double>::quiet_NaN());
      for (unsigned int dim = 0u; dim < this->dimension; dim++)
      {
        double meanThisDim = 0u;
        for (unsigned int i = 0u; i < this->numberOfDraws; i++)
        {
          meanThisDim += this->drawMatrix(i, dim);
        }
        meanThisDim /= double(this->numberOfDraws);
        meanPerDim(dim, 0u) = meanThisDim;

        double sum = 0.;
        for (unsigned int i = 0u; i < this->numberOfDraws; i++)
        {
          sum += std::pow(this->drawMatrix(i, dim) - meanThisDim, 2);
        }
        sum /= double(this->numberOfDraws);
        cov_matr(dim, dim) = sum;

        for (unsigned int dim2 = dim + 1; dim2 < this->dimension; dim2++)
        {
          if (meanPerDim(dim2, 0u) != meanPerDim(dim2, 0u))
          {
            double mean = 0u;
            for (unsigned int i = 0u; i < this->numberOfDraws; i++)
            {
              mean += this->drawMatrix(i, dim2);
            }
            mean /= double(this->numberOfDraws);
            meanPerDim(dim2, 0u) = mean;
          }

          double sum2 = 0.;
          for (unsigned int i = 0u; i < this->numberOfDraws; i++)
          {
            sum2 += (this->drawMatrix(i, dim2) - meanPerDim(dim2, 0u)) * (this->drawMatrix(i, dim) - meanPerDim(dim, 0u));
          }
          sum2 /= double(this->numberOfDraws);
          cov_matr(dim, dim2) = sum2;
          cov_matr(dim2, dim) = sum2;
        }
      }

      float_type cov_determ = 0.;
      int *cov_rank = new int;
      Matrix_Class eigenval, eigenvec;
      double determinant = cov_matr.determ();

      const double gaussentropy = 0.5 * log(cov_matr.determ()) + double(this->dimension) / 2. * std::log(2. * ::constants::pi * ::constants::e);
      return gaussentropy;
      //Covariance Matrix
    }
  }

  double cubatureIntegrationEntropy()
  {
    std::cout << "Commencing cubature integration of entropy." << std::endl;
    const double min_ = (this->probdens.meaningfulRange().first);
    const double max_ = (this->probdens.meaningfulRange().second);
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

    informationForCubatureIntegration info(this->probdens, this->subDimsGMM);

    double range = std::abs(this->probdens.meaningfulRange().first - this->probdens.meaningfulRange().second);

    std::cout << "Relative error (L2 norm) for cubature integration is " << std::scientific <<
      std::setw(5) << Config::get().entropytrails.errorThresholdForCubatureIntegration << "." << std::endl;

    int returncode = hcubature_v(1, cubaturefunctionProbDens, &(info),
      this->dimension, xmin, xmax,
      0u, 0, Config::get().entropytrails.errorThresholdForCubatureIntegration, ERROR_L2, val, err);

    if (returncode != 0)
    {
      std::cout << "Returncode was not zero. An error occured." << std::endl;
    }

    std::cout << "Integral of ProbDens is: " << std::defaultfloat << val[0] << " with error " << err[0] << "." << std::endl;
    

    returncode = hcubature_v(1, cubaturefunctionEntropy, &(info),
      this->dimension, xmin, xmax,
      0u, 0, Config::get().entropytrails.errorThresholdForCubatureIntegration, ERROR_L2, val, err);

    if (returncode != 0)
    {
      std::cout << "Returncode was not zero. An error occured." << std::endl;
    }

    const double value = val[0];

    std::cout << "Computed Entropy via cubature integration: " << value << " with error: " << err[0] << std::endl;


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
    const int n_of_samples = drawMatrix_TemporaryCopy.cols() - 1u;

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
        sum_gamma_k += digammal(i);
    }
    sum_gamma_k *= -1;
    sum_gamma_k *= static_cast<double>(this->dimension) / static_cast<double>(drawMatrix_TemporaryCopy.cols() * (drawMatrix_TemporaryCopy.cols() - 1));
    sum_gamma_k += digammal(drawMatrix_TemporaryCopy.cols());


    const double volume_of_unit_ball_for_eucledean_norm = pow(pi, static_cast<double>(this->dimension) / 2.) / tgamma(1. + static_cast<double>(this->dimension) / 2.);
    const double volume_of_unit_ball_for_max_norm = pow(2., this->dimension);


    transpose(drawMatrix);
    this->calculatedEntropyMeanFaivishevskyEucledean = mod_summation_eucl + sum_gamma_k + log(volume_of_unit_ball_for_eucledean_norm);
    this->calculatedEntropyMeanFaivishevskyMaximum = mod_summation_max + sum_gamma_k + log(volume_of_unit_ball_for_max_norm);
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

  void calculate()
  {
    if(Config::get().entropytrails.cubatureIntegration)
      this->cubatureIntegral = cubatureIntegrationEntropy();

    //mcdrawEntropy = this->MCDrawEntropy(drawMatrix);

    //mcintegrationEntropy = this->MCIntegrationEntropy(this->probdens.meaningfulRange(), numberOfDraws);

    empiricalNormalDistributionEntropy = this->empiricalGaussianEntropy();


    if(Config::get().entropytrails.meanNNcalculation)
      meanNNEntropyFaivishevsky();


    // Calculate Hnizdo as well as Lombardi/Pant entropy
    if (Config::get().entropytrails.NNcalculation)
    {
      std::cout << "Commencing NNEntropy calculation." << std::endl;
      // Matrix Layout after calculation:
      // First col: Drawn samples
      // Second col: k / iter * kNNdistance
      // Third col: Value of real distribution
      // Fourth col: kNN Distance

      //Neccessarry
      transpose(drawMatrix);
      Matrix_Class copytemp = drawMatrix;
      Matrix_Class eucl_kNN_distances(1u, numberOfDraws, 0.);
      Matrix_Class maxnorm_kNN_distances(1u, numberOfDraws, 0.);
      Matrix_Class eucl_kNN_distances_ardakani_corrected(1u, numberOfDraws, 0.);
      Matrix_Class maxnorm_kNN_distances_ardakani_corrected(1u, numberOfDraws, 0.);
      //std::function<std::vector<double>(std::vector<double> const& x)> PDFtemporary = this->probdens.function();
      std::vector<float_type> ardakaniCorrection_minimumValueInDataset(this->dimension, std::numeric_limits<float_type>::max());
      std::vector<float_type> ardakaniCorrection_maximumValueInDataset(this->dimension, -std::numeric_limits<float_type>::max());
      for (size_t k = 0; k < drawMatrix.cols(); k++)
      {
        for (unsigned int i = 0u; i < this->dimension; i++)
        {
          if (ardakaniCorrection_minimumValueInDataset.at(i) > drawMatrix(i, k))
            ardakaniCorrection_minimumValueInDataset.at(i) = drawMatrix(i, k);

          if (ardakaniCorrection_maximumValueInDataset.at(i) < drawMatrix(i, k))
            ardakaniCorrection_maximumValueInDataset.at(i) = drawMatrix(i, k);
        }
      }

#ifdef _OPENMP
#pragma omp parallel firstprivate(copytemp, ardakaniCorrection_minimumValueInDataset, ardakaniCorrection_maximumValueInDataset ) \
      shared(eucl_kNN_distances,maxnorm_kNN_distances, eucl_kNN_distances_ardakani_corrected, maxnorm_kNN_distances_ardakani_corrected)
      {
#endif
        float_type* buffer = new float_type[k];
#ifdef _OPENMP
        auto const n_omp = static_cast<std::ptrdiff_t>(numberOfDraws);

#pragma omp for
        for (std::ptrdiff_t i = 0; i < n_omp; ++i)
#else
        for (size_t i = 0u; i < numberOfDraws; i++)
#endif
        {
          std::vector<size_t> rowQueryPts;
          for (unsigned int currentDim = 0u; currentDim < this->dimension; currentDim++)
          {
            rowQueryPts.push_back(currentDim);
          }

          const float_type holdNNdistanceEucl = sqrt(entropy::knn_distance_eucl_squared(copytemp, this->dimension, k, rowQueryPts, i, buffer));
          eucl_kNN_distances(0, i) = holdNNdistanceEucl;

          const float_type holdNNdistanceMax = entropy::maximum_norm_knn_distance(copytemp, this->dimension, k, rowQueryPts, i, buffer);
          maxnorm_kNN_distances(0, i) = holdNNdistanceMax;
          //EvaluateMatrix(2, i) = PDFtemporary(copytemp(0, i));
          //EvaluateMatrix(3, i) = holdNNdistance;
          // Ardakani Correction

          std::vector<double> current;
          for (unsigned int j = 0u; j < this->dimension; j++)
            current.push_back(copytemp(j, i));

          eucl_kNN_distances_ardakani_corrected(0, i) = ardakaniCorrectionGeneralizedEucledeanNorm(ardakaniCorrection_minimumValueInDataset,
            ardakaniCorrection_maximumValueInDataset, current, holdNNdistanceEucl);
          maxnorm_kNN_distances_ardakani_corrected(0, i) = ardakaniCorrectionGeneralizedMaximumNorm(ardakaniCorrection_minimumValueInDataset,
            ardakaniCorrection_maximumValueInDataset, current, holdNNdistanceMax);
          // Lombardi kPN hier, fehlt bisher noch
        }
#ifdef _OPENMP
      }
#endif


      // Eucledean ArdakaniSum
      KahanAccumulation<double> kahan_acc_eucl_ardakani_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_eucl_ardakani_sum = KahanSum(kahan_acc_eucl_ardakani_sum,log(eucl_kNN_distances_ardakani_corrected(0, i)));

      double ardakaniSum = kahan_acc_eucl_ardakani_sum.sum / double(numberOfDraws);
      ardakaniSum *= double(dimension);
      ardakaniSum += log(pow(pi, double(dimension) / 2.) / (tgamma(0.5 * dimension + 1)));
      ardakaniSum += digammal(double(numberOfDraws));
      ardakaniSum -= digammal(double(k));

      // Maximum Norm ArdakaniSum
      KahanAccumulation<double> kahan_acc_max_ardakani_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_max_ardakani_sum = KahanSum(kahan_acc_max_ardakani_sum,log(maxnorm_kNN_distances_ardakani_corrected(0, i)));

      double maxArdakaniEntropy = kahan_acc_max_ardakani_sum.sum / double(numberOfDraws);
      maxArdakaniEntropy *= double(dimension);
      maxArdakaniEntropy += log(pow(2., this->dimension));
      maxArdakaniEntropy += digammal(double(numberOfDraws));
      maxArdakaniEntropy -= digammal(double(k));


      KahanAccumulation<double> kahan_acc_max_sum;
      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_max_sum = KahanSum(kahan_acc_max_sum,log(maxnorm_kNN_distances(0, i)));
      double maxNormSum = kahan_acc_max_sum.sum;

      // ENTROPY according to Hnzido
      KahanAccumulation<double> kahan_acc_eucl_sum;

      for (size_t i = 0u; i < numberOfDraws; i++)
        kahan_acc_eucl_sum = KahanSum(kahan_acc_eucl_sum, log(eucl_kNN_distances(0, i)));

      double hnizdoSum = kahan_acc_eucl_sum.sum;
      hnizdoSum /= double(numberOfDraws);
      hnizdoSum *= double(dimension);

      // Einschub
      // Entropy according to Lombardi (traditional)
      double lombardi_without_correction = hnizdoSum + log(pow(pi, double(dimension) / 2.) / (tgamma(0.5 * dimension + 1)));
      lombardi_without_correction += digammal(double(numberOfDraws));
      lombardi_without_correction -= digammal(double(k));
      //
      double lobardi_maximum_norm = maxNormSum / double(numberOfDraws);
      lobardi_maximum_norm *= double(dimension);
      lobardi_maximum_norm += log(pow(2., this->dimension));
      lobardi_maximum_norm += digammal(double(numberOfDraws));
      lobardi_maximum_norm -= digammal(double(k));

      //Einschub calcualtedEntropyGoria
      double tempsum_knn_goria = hnizdoSum;
      tempsum_knn_goria += log(pow(pi, double(dimension) / 2.) / (tgamma(0.5 * dimension + 1)));
      tempsum_knn_goria += log(double(numberOfDraws - 1.));
      tempsum_knn_goria -= digammal(double(k));

      //MaxNormGoria:
      double maxnorm_knn_goria = maxNormSum;
      maxnorm_knn_goria /= double(numberOfDraws);
      maxnorm_knn_goria *= double(dimension);
      maxnorm_knn_goria += log(pow(2., this->dimension));
      maxnorm_knn_goria += log(double(numberOfDraws - 1.));
      maxnorm_knn_goria -= digammal(double(k));

      hnizdoSum += (log(numberOfDraws * pow(pi, double(dimension) / 2.) / (tgamma(0.5 * dimension + 1))));

      // The following is identical to -digamma(double(k))
      double tempsum2 = 0;
      if (k != 1)
      {
        for (size_t i = 1; i < k; i++)
        {
          tempsum2 += 1.0 / double(i);
        }
      }
      hnizdoSum -= tempsum2;
      hnizdoSum += 0.5772156649015328606065;

      //////////////
      calculatedEntropyHnizdo = hnizdoSum; // Hnizdo Entropy
      calculatedEntropyLombardi = lombardi_without_correction; // Lombardi Entropy
      calculatedEntropyGoria = tempsum_knn_goria; // Goria Entropy
      ardakaniEntropyEucledean = ardakaniSum;
      ardakaniEntropyMaximum = maxArdakaniEntropy;
      calculatedEntropyGoriaMaximum = maxnorm_knn_goria;
      calculatedEntropyLombardiMaximum = lobardi_maximum_norm;
      //Neccessarry
      transpose(drawMatrix);
    }
  }

  void writeToFile()
  {
    //std::ofstream myfile;
    //myfile.open(std::string("out_k" + std::to_string(k) + "_i" + std::to_string(numberOfDraws) + "_d" + std::to_string(dimension) + "_ident" + std::to_string(identifierOfPDF) + ".txt"));
    //for (size_t i = 0u; i < drawAndEvaluateMatrix.rows(); i++)
    //{
    //    for (size_t j = 0u; j < drawAndEvaluateMatrix.cols(); j++)
    //      myfile << std::setw(15) << std::scientific << std::setprecision(5) << drawAndEvaluateMatrix(i, j);
    //
    //    myfile << "\n";
    //  }
    //myfile.close();


    std::ifstream myfile3;
    myfile3.open(std::string("out_entropy.txt"));
    bool writeHeader = false;
    if (!myfile3.good()) 
      writeHeader = true;
    myfile3.close();
    std::ofstream myfile2;
    myfile2.open(std::string("out_entropy.txt"), std::ios::app);
    if (writeHeader)
    {
      // HEADER
      myfile2 << "ENTROPY ESTIMATION\n";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Distribution|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "NumberOfDraws|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "k for NN|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Analytical|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Cubature|";
      //myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "MC-Integral|";
      //myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "MC-Draw|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Hnizdo|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Lombardi(eucl)|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Lombardi(max)|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Ardakani(eucl)|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Ardakani(max)|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Goria(eucl)|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Goria(max)|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Empirical Gauss|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "meanNNEucl|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "meanNNMax|";

      myfile2 << "\n=========================================================================";
      myfile2 << "==================================================================================";
      myfile2 << "=====================================================================================================\n";
    }
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->probdens.identString() << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->numberOfDraws << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->k << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->probdens.analyticEntropy() << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->cubatureIntegral << "|";
    //myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->mcintegrationEntropy << "|";
    //myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->mcdrawEntropy << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyHnizdo << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyLombardi << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyLombardiMaximum << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->ardakaniEntropyEucledean << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->ardakaniEntropyMaximum << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyGoria << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyGoriaMaximum << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->empiricalNormalDistributionEntropy << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyMeanFaivishevskyEucledean << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyMeanFaivishevskyMaximum << "|\n";
    myfile2.close();
  }
};
