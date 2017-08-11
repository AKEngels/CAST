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



// Returns the ardakani Corrected NN distance, see PHYSICAL REVIEW E 83, 051121 (2011)
float_type ardakaniCorrection1D(float_type const& globMin, float_type const& globMax, float_type const& currentPoint, float_type const& NNdistance)
{
  if (currentPoint - NNdistance * 0.5 < globMin)
    return currentPoint + NNdistance * 0.5 - globMin;
  else if (currentPoint + NNdistance * 0.5 > globMax)
    return globMax - (currentPoint - NNdistance * 0.5);
  return NNdistance;
}

template<typename T>
float_type ardakaniCorrectionGeneralizedEucledeanNorm(std::vector<T> const& globMin, std::vector<T> const& globMax, std::vector<T> const& currentPoint, T const& NNdistance)
{
#ifdef _DEBUG
  if (!(globMin.size() == globMax.size() && globMax.size() == currentPoint.size() && NNdistance.size() == currentPoint.size()))
  {
    throw std::runtime_error("Size Mismatch in N Dimensional Eucledean Ardakani Correction. Aborting.");
  }
#endif
  std::vector<T> radiiOfHyperEllipsoid(globMin.size());
  for (unsigned int i = 0u; i < radius.size(), i++)
  {
    radiiOfHyperEllipsoid.at(i) = std::min(currentPoint.at(i) + NNdistance / 0.5, globMax.at(i)) - std::max(currentPoint.at(i) - NNdistance / 0.5, globMin.at(i))
  }

  // Getting determinant of Mahalanobis distance for calculation of
  // Hyperelipsoid volume. Determinant is product of all eigenvalues. 
  // Eigenvalues of Hyperelipsoid quartic matrix are radius^-2
  T determinant = 1.;
  for (unsigned int i = 0u; i < radiiOfHyperEllipsoid.size(), i++)
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
      maximumOfPDF = 1. / sqrt(2. * pi);
      analyticEntropy_ = log(1. * sqrt(2. * pi * e));
      PDF = [&, this](double x) { return (1. / sqrt(2. * pi * 1. * 1.)) * exp(-1.*(x*x) / double(2. * 1. * 1.)); };
      //PDFrange = 5.920; // range incorporating function values up to under 10e-7
      this->m_identString = "Normal sig=1";
    }
    else if (ident_ == 1)
    {
      // Gaussian with sigma 10
      PDF = [&, this](double x) { return (1. / sqrt(2. * pi * 10. * 10.)) * exp(-1.*(x*x) / double(2. * 10. * 10.)); };
      maximumOfPDF = 1. / (10. * sqrt(2. * pi));
      analyticEntropy_ = log(10. * sqrt(2. * pi * e));
      //PDFrange = 55.136; // range incorporating function values up to under 10e-7
      this->m_identString = "Normal sig=10";
    }
    else if (ident_ == 2)
    {
      // Parabola normated to integrate to 1 over [0;1]
      PDF = [&, this](double x) { if (x < -1. || x > 1.) return 0.; else return (3. / 4.) * (-1. * ((x)*(x)) + 1.); };
      maximumOfPDF = PDF(0.);
      analyticEntropy_ = 0.568054;
      PDFrange = std::make_shared<std::pair<double, double>>(-1.,1.); // Function is only defined in [-1;1]
      this->m_identString = "parabolic";
    }
    else if (ident_ == 3)
    {
      // Beta Distribution https://en.wikipedia.org/wiki/Differential_entropy
      PDF = [&, this](double x) 
      {
        constexpr double alpha = 1.5;
        constexpr double beta = 9.;
        x = abs(x);
        if (x > 1.) 
          return 0.; 
        else 
          return pow(x,alpha -1.) * pow(1. -x, beta -1.) / (tgamma(alpha)*tgamma(beta) / tgamma(alpha + beta));
      };
      maximumOfPDF = 4.;
      constexpr double alpha = 1.5;
      constexpr double beta = 9.;
      analyticEntropy_ = log(tgamma(alpha)*tgamma(beta) / tgamma(alpha + beta)) - (alpha -1.) * (digammal(alpha) - digammal(alpha+beta)) - (beta -1.) * (digammal(beta) - digammal(alpha + beta));
      PDFrange = std::make_shared<std::pair<double, double>>(0., 1.); // Function is only defined in [0;1]
      this->m_identString = "Beta Distr.";
    }
  };

  double draw()
  {
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> unifDistr(0, 1);

    //Target PDF in 0.1er Schritten auswerten und schauen wann es kleiner ist
    bool run = true;
    double PDFrange = 0.;
    while (run)
    {
      if (PDF(PDFrange) < 0.0000001) break;
      else PDFrange += 0.1;
    }
    std::uniform_real_distribution<double> unifDistrRange(-PDFrange, PDFrange);

    double drawUnif = unifDistr(gen);
    double drawUnifWithRange = unifDistrRange(gen);
    double M = maximumOfPDF / (1. / (2.*PDFrange));
    while (true)
    {
      if (drawUnif < PDF(drawUnifWithRange) / (M * 1.5 * 1. / (2.*PDFrange)))
      {
        return(drawUnifWithRange);
      }
    }
  };

  std::vector<double> draw(size_t numberOfSamples)
  {
    std::random_device rd;
    std::mt19937 gen;
    //
    std::vector<double> draws;
    draws.reserve(numberOfSamples);
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
      double drawUnif = unifDistr(gen);
      double drawUnifWithRange = unifDistrRange(gen);
      double M = maximumOfPDF / (1. / (2.*absrange));
      if (drawUnif < PDF(drawUnifWithRange) / (M * 1.5 * 1. / (2.*absrange)))
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

  std::function<double(double x)> function() {
    return this->PDF;
  };

  std::pair<double,double> meaningfulRange()
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
      if (PDF(newPDFrange) < 0.00000001)
        break;
      else newPDFrange += 0.001;
    }
    rangeOne = newPDFrange;
    newPDFrange = 0.;
    while (run)
    {
      if (PDF(newPDFrange) < 0.000000001)
        break;
      else newPDFrange -= 0.005;
    }
    
    return std::make_pair(newPDFrange, rangeOne); 
  }

private:
  double analyticEntropy_;
  double maximumOfPDF;
  int identif;
  std::shared_ptr<std::pair<double,double>> PDFrange = nullptr;
  std::string m_identString = "ERROR_NAME";
  std::function<double(double x)> PDF;
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
  Matrix_Class drawAndEvaluateMatrix;
  int identifierOfPDF;
  ProbabilityDensity probdens;
  double mean, standardDeviation;
  entropyobj(size_t iter_, size_t dimension_, ProbabilityDensity probdens_) :
    drawAndEvaluateMatrix(iter_, 5),
    numberOfDraws(iter_), dimension(dimension_),
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
      while (std::getline(myfile, line))
      {
        drawAndEvaluateMatrix(j, 0) = std::stod(line);
        j++;
      }
    }
    else
    {
      // Draw 
      std::vector<double> draws = probdens.draw(numberOfDraws);

      //Sort samples (entirely optional)
      std::stable_sort(draws.begin(), draws.end());

      for (unsigned int n = 0; n < numberOfDraws; ++n)
        drawAndEvaluateMatrix(n, 0) = draws[n];

      // Write 
      std::ofstream myfile2;
      myfile2.open(std::string("draw_i" + std::to_string(numberOfDraws) + "_d" + std::to_string(dimension) + "_ident" + std::to_string(identifierOfPDF) + ".txt"));
      for (size_t w = 0u; w < drawAndEvaluateMatrix.rows(); w++)
      {
        myfile2 << std::setw(30) << std::scientific << std::setprecision(15) << drawAndEvaluateMatrix(w, 0);
        myfile2 << "\n";
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
    // Calculate Mean
    for (unsigned int n = 0; n < numberOfDraws; ++n)
      mean += drawAndEvaluateMatrix(n, 0);
    mean /= numberOfDraws;

    // Calculate standard Deviation
    for (unsigned int n = 0; n < numberOfDraws; ++n)
      standardDeviation += (drawAndEvaluateMatrix(n, 0) - mean) * (drawAndEvaluateMatrix(n, 0) - mean);

    standardDeviation /= numberOfDraws;
    standardDeviation = sqrt(standardDeviation);

    return log(standardDeviation * sqrt(2. * pi * e));
  }

  double MCIntegrationEntropy(std::pair<double,double> const& range, size_t numberOfSamples)
  {
    // MC guess with uniform distribution
    // Draw Uniform and do
    // Draw samples
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distr(range.first, range.second);
    std::vector<double> temp;
    temp.reserve(numberOfSamples);
    for (unsigned int n = 0; n < numberOfSamples; ++n) {
      double drawnnum = distr(gen);
      temp.push_back(drawnnum);
    }

    long double uniformMCvalue = 0, uniformMCvalue2 = 0, c = 0.f;

    // Debug
    //Sort samples
    //std::stable_sort(temp.begin(), temp.end());
    //size_t countHowManyContinuesBecauseOfZeroPDF = 0u;

    for (unsigned int n = 0; n < numberOfSamples; ++n)
    {
      long double gvalue;
      gvalue = static_cast<long double>(probdens.function()(temp[n]));
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

  double MCDrawEntropy(std::vector<double>& samples)
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
    std::vector<double> samples = this->probdens.draw(numberOfSamples);
    return this->MCDrawEntropy(samples);
  };

  double MCDrawEntropy(Matrix_Class& samples)
  {
    std::vector<double> samples_vec;
    samples_vec.reserve(samples.rows());

    for (size_t i = 0u; i < samples.rows(); ++i)
      samples_vec.push_back(samples(i, 0u));

    return this->MCDrawEntropy(samples_vec);
  }

  double meanNNEntropyFaivishevsky()
  {
    //http://papers.nips.cc/paper/3500-ica-based-on-a-smooth-estimation-of-the-differential-entropy
    transpose(drawAndEvaluateMatrix);
    Matrix_Class drawAndEvaluateMatrix_TemporaryCopy = drawAndEvaluateMatrix;
    KahanAccumulation<double> kahan_acc;
    double sum_final = 0.f;

#ifdef _OPENMP
#pragma omp parallel firstprivate(drawAndEvaluateMatrix_TemporaryCopy,kahan_acc) reduction(+:sum_final)
#endif
    {
#ifdef _OPENMP
      auto const n_omp = static_cast<std::ptrdiff_t>(drawAndEvaluateMatrix_TemporaryCopy.cols());
#pragma omp for
      for (std::ptrdiff_t i = 0; i < n_omp; ++i)
#else
      for (size_t i = 0u; i < drawAndEvaluateMatrix_TemporaryCopy.cols(); i++)
#endif
      {
        for (size_t j = i + 1u; j < drawAndEvaluateMatrix_TemporaryCopy.cols(); j++)
        {
          kahan_acc = KahanSum(kahan_acc, log(std::abs(drawAndEvaluateMatrix_TemporaryCopy(0u, i) - drawAndEvaluateMatrix_TemporaryCopy(0u, j))));
        }
      }
      sum_final += kahan_acc.sum;
    }
    const double mod_summation = sum_final * static_cast<double>(this->dimension) / static_cast<double>(drawAndEvaluateMatrix.cols() * (drawAndEvaluateMatrix.cols() - 1));
    
    double sum_gamma_k = 0.f;
    const int n_of_samples = drawAndEvaluateMatrix.cols() - 1u;

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
    sum_gamma_k *= static_cast<double>(this->dimension) / static_cast<double>(drawAndEvaluateMatrix.cols() * (drawAndEvaluateMatrix.cols() - 1));
    sum_gamma_k += digammal(drawAndEvaluateMatrix.cols());
    

    const double volume_of_unit_ball_for_eucledean_norm = pow(pi, static_cast<double>(this->dimension) / 2.) / tgamma(1. + static_cast<double>(this->dimension) / 2.);


    transpose(drawAndEvaluateMatrix);
    return mod_summation + sum_gamma_k + log(volume_of_unit_ball_for_eucledean_norm);
  }

  void calculate()
  {
    // Calculates entropy inegral using MC 
    // with known PDF
    mcdrawEntropy = this->MCDrawEntropy(drawAndEvaluateMatrix);

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
      transpose(drawAndEvaluateMatrix);
      Matrix_Class copytemp = drawAndEvaluateMatrix;
      Matrix_Class drawAndEvaluateMatrix_TemporaryCopy = drawAndEvaluateMatrix;
      std::function<double(double x)> PDFtemporary = this->probdens.function();
      float_type ardakaniCorrection_minimumValueInDataset(std::numeric_limits<float_type>::max());
      float_type ardakaniCorrection_maximumValueInDataset(-std::numeric_limits<float_type>::max());
      for (size_t k = 0; k < drawAndEvaluateMatrix.cols(); k++)
      {
        ardakaniCorrection_minimumValueInDataset = std::min(ardakaniCorrection_minimumValueInDataset, drawAndEvaluateMatrix(0, k));
        ardakaniCorrection_maximumValueInDataset = std::max(ardakaniCorrection_maximumValueInDataset, drawAndEvaluateMatrix(0, k));
      }

#ifdef _OPENMP
#pragma omp parallel firstprivate(copytemp, PDFtemporary) shared(drawAndEvaluateMatrix_TemporaryCopy)
     {
#endif
        float_type* buffer = new float_type[k];
#ifdef _OPENMP
        auto const n_omp = static_cast<std::ptrdiff_t>(drawAndEvaluateMatrix_TemporaryCopy.cols());

#pragma omp for
        for (std::ptrdiff_t i = 0; i < n_omp; ++i)
#else
        for (size_t i = 0u; i < drawAndEvaluateMatrix_TemporaryCopy.cols(); i++)
#endif
        {
          const float_type holdNNdistance = sqrt(entropy::knn_distance(copytemp, 1, k, 0u, i, buffer));
          drawAndEvaluateMatrix_TemporaryCopy(1, i) = double(k) / double(numberOfDraws) / holdNNdistance;
          drawAndEvaluateMatrix_TemporaryCopy(2, i) = PDFtemporary(copytemp(0, i));
          drawAndEvaluateMatrix_TemporaryCopy(3, i) = holdNNdistance;
          // Ardakani Correction
          drawAndEvaluateMatrix_TemporaryCopy(4, i) = ardakaniCorrection1D(ardakaniCorrection_minimumValueInDataset, ardakaniCorrection_maximumValueInDataset, drawAndEvaluateMatrix(0, i), holdNNdistance);
         
          // Lombardi kPN hier, fehlt bisher noch
        }
#ifdef _OPENMP
      }
#endif

     drawAndEvaluateMatrix = drawAndEvaluateMatrix_TemporaryCopy;

     double ardakaniSum = 0.f;
     for (size_t i = 0u; i < drawAndEvaluateMatrix.cols(); i++)
       ardakaniSum += log(drawAndEvaluateMatrix(4, i));
     ardakaniSum /= double(numberOfDraws);
     ardakaniSum *= double(dimension);
     ardakaniSum += (log(numberOfDraws * pow(pi, double(dimension) / 2.)) / (tgamma(0.5 * dimension + 1)));
     ardakaniSum += digammal(double(numberOfDraws));
     ardakaniSum -= digammal(double(k));

      // ENTROPY according to Hnzido
      double hnizdoSum = 0.;
      for (size_t i = 0u; i < drawAndEvaluateMatrix.cols(); i++)
        hnizdoSum += log(drawAndEvaluateMatrix(3, i));

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
      transpose(drawAndEvaluateMatrix);
    }
  }

  void writeToFile()
  {
    std::ofstream myfile;
    myfile.open(std::string("out_k" + std::to_string(k) + "_i" + std::to_string(numberOfDraws) + "_d" + std::to_string(dimension) + "_ident" + std::to_string(identifierOfPDF) + ".txt"));
    for (size_t i = 0u; i < drawAndEvaluateMatrix.rows(); i++)
    {
        for (size_t j = 0u; j < drawAndEvaluateMatrix.cols(); j++)
          myfile << std::setw(15) << std::scientific << std::setprecision(5) << drawAndEvaluateMatrix(i, j);

        myfile << "\n";
      }
    myfile.close();


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
