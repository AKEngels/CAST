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
    radiiOfHyperEllipsoid.at(i) = std::min(currentPoint.at(i) + NNdistance / 0.5, globMax.at(i)) - std::max(currentPoint.at(i) - NNdistance / 0.5, globMin.at(i));
  }

  // Getting determinant of Mahalanobis distance for calculation of
  // Hyperelipsoid volume. Determinant is product of all eigenvalues. 
  // Eigenvalues of Hyperelipsoid quartic matrix are radius^-2
  T determinant = 1.;
  for (unsigned int i = 0u; i < radiiOfHyperEllipsoid.size(); i++)
  {
    determinant *= std::pow(radiiOfHyperEllipsoid.at(i), -2);
  }

  // For volume of hyperelipsoid https://math.stackexchange.com/questions/332391/volume-of-hyperellipsoid
  // Is radius of hypersphere volume 1? I think, because we have all scaling the eigenvalues
  // But I am not sure...
  const T V_d = 2. / radiiOfHyperEllipsoid.size() * (std::pow(pi, radiiOfHyperEllipsoid.size() / 2.) / tgamma(radiiOfHyperEllipsoid.size() / 2.)) * std::pow(1, radiiOfHyperEllipsoid.size());
  const T volumeOfHyperellipsoid = V_d * std::sqrt(determinant);
  const T equivalentRadiusOfHypersphere = std::pow(volumeOfHyperellipsoid * tgamma(radiiOfHyperEllipsoid.size() / 2. + 1) / std::pow(pi, radiiOfHyperEllipsoid.size() / 2.), 1. / radiiOfHyperEllipsoid.size());
  return equivalentRadiusOfHypersphere;

}

template<typename T>
float_type ardakaniCorrectionGeneralizedMaximumNorm(std::vector<T> const& globMin, std::vector<T> const& globMax, std::vector<T> const& currentPoint, T const& NNdistance)
{
#ifdef _DEBUG
  if (!(globMin.size() == globMax.size() && globMax.size() == currentPoint.size() && NNdistance.size() == currentPoint.size()))
  {
    throw std::runtime_error("Size Mismatch in N Dimensional Eucledean Ardakani Correction. Aborting.");
  }
#endif
  std::vector<T> radiiOfBox(globMin.size());
  for (unsigned int i = 0u; i < radiiOfBox.size(); i++)
  {
    radiiOfBox.at(i) = 
      std::min(currentPoint.at(i) + NNdistance / 0.5, globMax.at(i)) 
      - std::max(currentPoint.at(i) - NNdistance / 0.5, globMin.at(i));
  }

  T volumeOfBox = 1.;
  for (unsigned int i = 0u; i < radiiOfHyperEllipsoid.size(); i++)
  {
    volumeOfBox *= radiiOfHyperEllipsoid.at(i);
  }

  const T equivalentRadiusOfHyperSquare = std::pow(volumeOfBox, 1. / radiiOfHyperEllipsoid.size());
  return equivalentRadiusOfHyperSquare;

}





// Class for PDF
//
// PDF Function needs to be centered at 0,0
class ProbabilityDensity
{
public:
  ProbabilityDensity(int ident_) : identif(ident_)
  {
    if (ident_ == 0)
    {
      this->dimension = 1;
      maximumOfPDF = 1. / sqrt(2. * pi);
      analyticEntropy_ = log(1. * sqrt(2. * pi * e));
      PDF = [&, this](std::vector<double> const& x) 
      {         
        if (x.size() != 1)
          throw std::runtime_error("Wring dimensionality for chosen Probability Density.");
        return (1. / sqrt(2. * pi * 1. * 1.)) * exp(-1.*(x.at(0)*x.at(0)) / double(2. * 1. * 1.)) ;
      };
      //PDFrange = 5.920; // range incorporating function values up to under 10e-7
      this->m_identString = "Normal sig=1";
    }
    else if (ident_ == 1)
    {
      this->dimension = 1;
      // Gaussian with sigma 10
      PDF = [&, this](std::vector<double> const& x)
      { 
        if (x.size() != 1)
          throw std::runtime_error("Wring dimensionality for chosen Probability Density.");
        return (1. / sqrt(2. * pi * 10. * 10.)) * exp(-1.*(x.at(0)*x.at(0)) / double(2. * 10. * 10.)) ;
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
      PDF = [&, this](std::vector<double> const& x) 
      { 
        if (x.size() != 1)
          throw std::runtime_error("Wring dimensionality for chosen Probability Density.");
        if (x.at(0) < -1. || x.at(0) > 1.)
          return 0.;
        else 
          return (3. / 4.) * (-1. * ((x.at(0))*(x.at(0))) + 1.) ;
      };
      maximumOfPDF = PDF(std::vector<double>{ 0. });
      analyticEntropy_ = 0.568054;
      PDFrange = std::make_shared<std::pair<double, double>>(-1.,1.); // Function is only defined in [-1;1]
      this->m_identString = "parabolic";
    }
    else if (ident_ == 3)
    {
      this->dimension = 1;
      // Beta Distribution https://en.wikipedia.org/wiki/Differential_entropy
      PDF = [&, this](std::vector<double> const& x)
      {
        if (x.size() != 1)
          throw std::runtime_error("Wring dimensionality for chosen Probability Density.");
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
      analyticEntropy_ = log(tgamma(alpha)*tgamma(beta) / tgamma(alpha + beta)) - (alpha -1.) * (digammal(alpha) - digammal(alpha+beta)) - (beta -1.) * (digammal(beta) - digammal(alpha + beta));
      PDFrange = std::make_shared<std::pair<double, double>>(0., 1.); // Function is only defined in [0;1]
      this->m_identString = "Beta Distr.";
    }
  };


  std::vector<std::vector<double>> draw(size_t numberOfSamples)
  {
    std::random_device rd;
    std::mt19937 gen;
    //
    std::vector<std::vector<double>> draws;
    //draws.reserve(numberOfSamples);
    std::normal_distribution<double> normDistrSigmaOne(0, 1);
    std::function<double(double x)> normDistrSigmaOnePDF = [&, this](double x) { return (1. / sqrt(2. * pi)) * exp(-1.*(x*x) / double(2.)); };
    double maximumOfNormalDistrWithSigmaOne = 1. / sqrt(2. * pi);
    std::uniform_real_distribution<double> unifDistr(0, 1);

    //Target PDF in 0.1er Schritten auswerten und schauen wann es kleiner ist
    if (PDFrange == nullptr)
      PDFrange = std::make_shared<std::pair<double, double>>(meaningfulRange());



    std::uniform_real_distribution<double> unifDistrRange(PDFrange->first, PDFrange->second);
    const double absrange = PDFrange->second - PDFrange->first;

    for (unsigned int n = 0; n < numberOfSamples; ++n)
    {
      std::vector<double> drawUnifWithRange;
      //double drawUnif = unifDistr(gen);
      //double drawUnifWithRange = unifDistrRange(gen);
      double M = maximumOfPDF / (1. / (2.*absrange));
      bool inRange = true;
      for (unsigned int currentDim = 0u; currentDim < this->dimension; currentDim++)
      {
        drawUnifWithRange.push_back(unifDistrRange(gen));
      }
      double drawUnif = unifDistr(gen);
      inRange = inRange && (drawUnif < PDF(drawUnifWithRange) / (M * 1.5 * 1. / (2.*absrange)));

      if (inRange)
      {
        draws.push_back(drawUnifWithRange);
      }
      else
      {
        --n;
      }
    }
    return draws;
  };

  int ident() const { return this->identif; } ;

  std::string identString() const
  {
    return this->m_identString;
  }

  double analyticEntropy() {
    return this->analyticEntropy_;
  };

  auto function() {
    return this->PDF;
  };

  std::pair<double,double> meaningfulRange()
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
        if (PDF(std::vector<double> {newPDFrange}) < 0.00000001)
          break;
        else
          newPDFrange += 0.001;
      }
      rangeOne = newPDFrange;
      newPDFrange = 0.;
      while (run)
      {
        if (PDF(std::vector<double> {newPDFrange}) < 0.000000001)
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
  std::shared_ptr<std::pair<double,double>> PDFrange = nullptr;
  std::string m_identString = "ERROR_NAME";
  std::function<double(std::vector<double> const& x)> PDF;
};


// entropyobj
// Contains a matrix with draws from a distribution and 
// the number of draws(iter), dimensionality(dimension)
// and standard deviation (sigma)
// 
// No calculations have been done here
//
// CURRENTLY ONLY WORKS WITH NORMAL DISTRIBUTION
// CURRENTLY DIMENSION NEEDS TO BE 1
//
// drawAndEvaluateMatrix:
// Matrix Layout after calculation:
// First col: Drawn samples
// Second col: k / iter * kNNdistance
// Third col: Value of real distribution
// Fourth col: kNN Distance
class entropyobj
{
public:
  size_t numberOfDraws, dimension;
  Matrix_Class drawAndsEvaluateMatrix;
  Matrix_Class drawMatrix;
  int identifierOfPDF;
  ProbabilityDensity probdens;
  double mean, standardDeviation;
  entropyobj(size_t iter_, size_t dimension_, ProbabilityDensity probdens_) :
    drawAndsEvaluateMatrix(iter_, 5),
    numberOfDraws(iter_), dimension(dimension_), drawMatrix(iter_,dimension_),
    probdens(probdens_), identifierOfPDF(probdens_.ident()),
    mean(0.), standardDeviation(0.)
  {
    // Check if draw exists
    std::ifstream myfile;
    std::string line;
    myfile.open(std::string("draw_i" + std::to_string(numberOfDraws) + "_d" + std::to_string(dimension) + "_ident" + std::to_string(identifierOfPDF) + ".txt"));
    if (myfile.good())
    {
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
      // Draw 
      std::vector<std::vector<double>> draws = probdens.draw(numberOfDraws);

      //Sort samples (entirely optional)
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


// Calculated entropy object calculates estiamted entropy
class calculatedentropyobj : public entropyobj
{
public:
  size_t k;
  double calculatedEntropyHnizdo;
  double calculatedEntropyLombardi; 
  double calculatedEntropyMeanFaivishevsky;
  double mcintegrationEntropy;
  double mcdrawEntropy;
  double analyticalEntropy;
  double empiricalNormalDistributionEntropy;
  double calculatedEntropyGoria; // http://www.tandfonline.com/doi/abs/10.1080/104852504200026815
  double ardakaniEntropyEucledean;
  calculatedentropyobj(size_t k_, entropyobj const& obj) : 
    entropyobj(obj), 
    k(k_), 
    calculatedEntropyHnizdo(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyLombardi(std::numeric_limits<double>::quiet_NaN()),
    mcintegrationEntropy(std::numeric_limits<double>::quiet_NaN()),
    mcdrawEntropy(std::numeric_limits<double>::quiet_NaN()),
    analyticalEntropy(this->probdens.analyticEntropy()),
    empiricalNormalDistributionEntropy(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyMeanFaivishevsky(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyGoria(std::numeric_limits<double>::quiet_NaN()),
    ardakaniEntropyEucledean(std::numeric_limits<double>::quiet_NaN())
  {
  }

  double empiricalGaussianEntropy()
  {
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

      return log(standardDeviation * sqrt(2. * pi * e));
    }
    else
    {
      Matrix_Class ones(this->drawMatrix.cols(), this->drawMatrix.cols(), 1.0);
      Matrix_Class cov_matr = Matrix_Class(cov_matr - ones * cov_matr / static_cast<float_type>(this->drawMatrix.cols()));
      cov_matr = Matrix_Class(transposed(cov_matr) * cov_matr);
      cov_matr = cov_matr / static_cast<float_type>(this->drawMatrix.cols());
      float_type cov_determ = 0.;
      int *cov_rank = new int;
      Matrix_Class eigenval, eigenvec;
      double determinant = cov_matr.determ();

      return (log(sqrt(2. * pi)) + 1) -0.5  *double(this->dimension) + 0.5 * log(determinant);
      //Covariance Matrix
    }
  }

  double MCIntegrationEntropy(std::pair<double,double> const& range, size_t numberOfSamples)
  {
    // MC guess with uniform distribution
    // Draw Uniform and do
    // Draw samples
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distr(range.first, range.second);
    std::vector<std::vector<double>> temp;
    //temp.reserve(numberOfSamples);
    for (unsigned int n = 0; n < numberOfSamples; ++n) 
    {
      temp.push_back(std::vector<double>{});
      for (unsigned int currentDim = 0; currentDim < this->dimension; ++currentDim)
      {
        double drawnnum = distr(gen);
        temp.at(currentDim).push_back(drawnnum);
      }
    }

    long double uniformMCvalue = 0, uniformMCvalue2 = 0, c = 0.f;

    // Debug
    //Sort samples
    //std::stable_sort(temp.begin(), temp.end());
    //size_t countHowManyContinuesBecauseOfZeroPDF = 0u;

    for (unsigned int n = 0; n < numberOfSamples; ++n)
    {
      long double gvalue;
      std::vector<double> current;
      for (unsigned int currentDim = 0; currentDim < this->dimension; ++currentDim)
      {
        current.push_back(temp.at(currentDim).at(n));
      }

      gvalue = static_cast<long double>(probdens.function()(current));
      if (gvalue == 0) {
        //++countHowManyContinuesBecauseOfZeroPDF;  
        continue;
      }
      long double y = gvalue * log(gvalue) - c;
      long double t = uniformMCvalue + y;
      c = (t - uniformMCvalue) - y;
      uniformMCvalue = t;
      uniformMCvalue2 += gvalue * log(gvalue);
    }
    uniformMCvalue /= double(numberOfSamples);
    uniformMCvalue2 /= double(numberOfSamples);
    uniformMCvalue *= -1. * (range.second - range.first);
    uniformMCvalue2 *= -1. * (range.second - range.first);
    return uniformMCvalue2;
  };

  double MCDrawEntropy(std::vector<std::vector<double>>& samples)
  {
    // new try MC guess with Kahan Summation
    double drawMCvalue = 0.;
    long double drawMCvalue2 = 0;
    long double c = 0;
    for (size_t i = 0u; i < samples.size(); i++)
    {
      // Value of Distribution at point x
      double temp1 = probdens.function()(samples[i]);
      drawMCvalue += log(temp1);
      long double y = log(temp1) - c;
      long double t = drawMCvalue2 + y;
      c = (t - drawMCvalue2) - y;
      drawMCvalue2 = t;
    }
    drawMCvalue /= -1.* double(samples.size());
    drawMCvalue2 /= -1 * double(samples.size());
    return drawMCvalue2;
  };

  double MCDrawEntropy(size_t& numberOfSamples)
  {
    std::vector<std::vector<double>> samples = this->probdens.draw(numberOfSamples);
    return this->MCDrawEntropy(samples);
  };

  double MCDrawEntropy(Matrix_Class& samples)
  {
    std::vector<std::vector<double>> samples_vec;
    //samples_vec.reserve(samples.rows());

    for (size_t i = 0u; i < samples.rows(); ++i)
    {
      samples_vec.push_back(std::vector<double>{});
      for (size_t j = 0u; j < samples.cols(); ++j)
      {
        samples_vec.at(i).push_back(samples(i, j));
      }
    }

    return this->MCDrawEntropy(samples_vec);
  }

  double meanNNEntropyFaivishevsky()
  {
    //http://papers.nips.cc/paper/3500-ica-based-on-a-smooth-estimation-of-the-differential-entropy
    transpose(drawMatrix);
    Matrix_Class drawMatrix_TemporaryCopy = drawMatrix;
    KahanAccumulation<double> kahan_acc;
    double sum_final = 0.f;

#ifdef _OPENMP
#pragma omp parallel firstprivate(drawMatrix_TemporaryCopy,kahan_acc) reduction(+:sum_final)
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
          double norm = 0.;
          for (size_t currentDim = 0u; currentDim < this->dimension; currentDim++)
          {
            norm += (drawMatrix_TemporaryCopy(0u, i) - drawMatrix_TemporaryCopy(0u, j))
              * (drawMatrix_TemporaryCopy(0u, i) - drawMatrix_TemporaryCopy(0u, j));
          }
          norm = sqrt(norm);
          kahan_acc = KahanSum(kahan_acc, log(norm));
        }
      }
      sum_final += kahan_acc.sum;
    }
    const double mod_summation = sum_final * static_cast<double>(this->dimension) / static_cast<double>(drawMatrix_TemporaryCopy.cols() * (drawMatrix_TemporaryCopy.cols() - 1));
    
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


    transpose(drawMatrix);
    return mod_summation + sum_gamma_k + log(volume_of_unit_ball_for_eucledean_norm);
  }

  void calculate()
  {
    // Calculates entropy inegral using MC 
    // with known PDF
    mcdrawEntropy = this->MCDrawEntropy(drawMatrix);

    mcintegrationEntropy = this->MCIntegrationEntropy(this->probdens.meaningfulRange(), numberOfDraws);

    empiricalNormalDistributionEntropy = this->empiricalGaussianEntropy();

    calculatedEntropyMeanFaivishevsky = meanNNEntropyFaivishevsky();

    // Calculate Hnizdo as well as Lombardi/Pant entropy
    if (Config::get().entropytrails.NNcalculation)
    {
      // Matrix Layout after calculation:
      // First col: Drawn samples
      // Second col: k / iter * kNNdistance
      // Third col: Value of real distribution
      // Fourth col: kNN Distance

      //Neccessarry
      transpose(drawMatrix);
      Matrix_Class copytemp = drawMatrix;
      Matrix_Class EvaluateMatrix(2, numberOfDraws, 0.);
      //std::function<std::vector<double>(std::vector<double> const& x)> PDFtemporary = this->probdens.function();
      std::vector<float_type> ardakaniCorrection_minimumValueInDataset(this->dimension, std::numeric_limits<float_type>::max());
      std::vector<float_type> ardakaniCorrection_maximumValueInDataset(this->dimension, -std::numeric_limits<float_type>::max());
      for (size_t k = 0; k < drawMatrix.cols(); k++)
      {
        std::vector<double> current;
        for (unsigned int i = 0u; i < this->dimension; i++)
          current.push_back(drawMatrix(i, k));
        const double euclDist = eucledeanDistance(ardakaniCorrection_minimumValueInDataset);

        if (eucledeanDistance(ardakaniCorrection_minimumValueInDataset) > euclDist)
          ardakaniCorrection_minimumValueInDataset = current;

        if (eucledeanDistance(ardakaniCorrection_maximumValueInDataset) < euclDist)
          ardakaniCorrection_maximumValueInDataset = current;
      }

#ifdef _OPENMP
#pragma omp parallel firstprivate(copytemp, ardakaniCorrection_minimumValueInDataset, ardakaniCorrection_maximumValueInDataset ) shared(EvaluateMatrix)
     {
#endif
        float_type* buffer = new float_type[k];
#ifdef _OPENMP
        auto const n_omp = static_cast<std::ptrdiff_t>(EvaluateMatrix.cols());

#pragma omp for
        for (std::ptrdiff_t i = 0; i < n_omp; ++i)
#else
        for (size_t i = 0u; i < drawAndEvaluateMatrix_TemporaryCopy.cols(); i++)
#endif
        {
          const float_type holdNNdistance = sqrt(entropy::knn_distance(copytemp, this->dimension, k, 0u, i, buffer));
          EvaluateMatrix(0, i) = holdNNdistance;
          //EvaluateMatrix(2, i) = PDFtemporary(copytemp(0, i));
          //EvaluateMatrix(3, i) = holdNNdistance;
          // Ardakani Correction
          std::vector<double> current;
          for (unsigned int j = 0u; j < this->dimension; j++)
            current.push_back(copytemp(j, i));

          EvaluateMatrix(1, i) = ardakaniCorrectionGeneralizedEucledeanNorm(ardakaniCorrection_minimumValueInDataset,
            ardakaniCorrection_maximumValueInDataset, current, holdNNdistance);
          // Lombardi kPN hier, fehlt bisher noch
        }
#ifdef _OPENMP
      }
#endif

     //drawAndEvaluateMatrix = drawAndEvaluateMatrix_TemporaryCopy;

     double ardakaniSum = 0.f;
     for (size_t i = 0u; i < EvaluateMatrix.cols(); i++)
       ardakaniSum += log(EvaluateMatrix(1, i));
     ardakaniSum /= double(numberOfDraws);
     ardakaniSum *= double(dimension);
     ardakaniSum += (log(pow(pi, double(dimension) / 2.)) / (tgamma(0.5 * dimension + 1)));
     ardakaniSum += digammal(double(numberOfDraws));
     ardakaniSum -= digammal(double(k));

      // ENTROPY according to Hnzido
      double hnizdoSum = 0.;
      for (size_t i = 0u; i < EvaluateMatrix.cols(); i++)
        hnizdoSum += log(EvaluateMatrix(0, i));

      hnizdoSum /= double(numberOfDraws);
      hnizdoSum *= double(dimension);

      // Einschub
      // Entropy according to Lombardi (traditional)
      double lombardi_without_correction = hnizdoSum + (log(pow(pi, double(dimension) / 2.)) / (tgamma(0.5 * dimension + 1)));
      lombardi_without_correction += digammal(double(numberOfDraws));
      lombardi_without_correction -= digammal(double(k));
      //

      //Einschub calcualtedEntropyGoria
      double tempsum_knn_goria = hnizdoSum + (log(pow(pi, double(dimension) / 2.)) / (tgamma(0.5 * dimension + 1)));
      tempsum_knn_goria += log(double(numberOfDraws - 1.));
      tempsum_knn_goria -= digammal(double(k));

      hnizdoSum += (log(numberOfDraws * pow(pi, double(dimension) / 2.)) / (tgamma(0.5 * dimension + 1)));

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
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "MC-Integral|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "MC-Draw|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Hnizdo|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Lombardi(fake)|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Ardakani(eucl)|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Goria|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Empirical Gauss|";
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "meanNN|";

      myfile2 << "\n===============================================================================================================================================================\n";
    }
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->probdens.identString() << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->numberOfDraws << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->k << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->probdens.analyticEntropy() << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->mcintegrationEntropy << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->mcdrawEntropy << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyHnizdo << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyLombardi << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->ardakaniEntropyEucledean << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyGoria << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->empiricalNormalDistributionEntropy << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyMeanFaivishevsky << "|\n";
    myfile2.close();
  }
};
