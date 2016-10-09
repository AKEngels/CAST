#pragma once
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <omp.h>
#include <limits>
#include "matop.h"
#include <string>



/////////////////
// Some constants
/////////////////

#ifndef M_PIl
/** The constant Pi in high precision */
#define M_PIl 3.1415926535897932384626433832795029L
#endif
#ifndef M_GAMMAl
/** Euler's constant in high precision */
#define M_GAMMAl 0.5772156649015328606065120900824024L
#endif
#ifndef M_LN2l
/** the natural logarithm of 2 in high precision */
#define M_LN2l 0.6931471805599453094172321214581766L
#endif

const double pi = M_PIl;
const double e = 2.71828182845904523536028747135266249775724709369995;

/////////////////
// Auxiliary Functions
/////////////////

// Gives DiGamma Function value, via Stackoverflow
long double digammal(long double x)
{
  /* force into the interval 1..3 */
  if (x < 0.0L)
    return digammal(1.0L - x) + M_PIl / tanl(M_PIl*(1.0L - x));	/* reflection formula */
  else if (x < 1.0L)
    return digammal(1.0L + x) - 1.0L / x;
  else if (x == 1.0L)
    return -M_GAMMAl;
  else if (x == 2.0L)
    return 1.0L - M_GAMMAl;
  else if (x == 3.0L)
    return 1.5L - M_GAMMAl;
  else if (x > 3.0L)
    /* duplication formula */
    return 0.5L*(digammal(x / 2.0L) + digammal((x + 1.0L) / 2.0L)) + M_LN2l;
  else
  {
    /* Just for your information, the following lines contain
    * the Maple source code to re-generate the table that is
    * eventually becoming the Kncoe[] array below
    * interface(prettyprint=0) :
    * Digits := 63 :
    * r := 0 :
    *
    * for l from 1 to 60 do
    * 	d := binomial(-1/2,l) :
    * 	r := r+d*(-1)^l*(Zeta(2*l+1) -1) ;
    * 	evalf(r) ;
    * 	print(%,evalf(1+Psi(1)-r)) ;
    *o d :
    *
    * for N from 1 to 28 do
    * 	r := 0 :
    * 	n := N-1 :
    *
    *	for l from iquo(n+3,2) to 70 do
    *		d := 0 :
    *		for s from 0 to n+1 do
    *		 d := d+(-1)^s*binomial(n+1,s)*binomial((s-1)/2,l) :
    *		od :
    *		if 2*l-n > 1 then
    *		r := r+d*(-1)^l*(Zeta(2*l-n) -1) :
    *		fi :
    *	od :
    *	print(evalf((-1)^n*2*r)) ;
    *od :
    *quit :
    */
    static long double Kncoe[] = { .30459198558715155634315638246624251L,
      .72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
      .27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
      .17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
      .11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
      .83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
      .59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
      .42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
      .304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
      .21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
      .15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
      .11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
      .80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
      .58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
      .41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L };

    register long double Tn_1 = 1.0L;	/* T_{n-1}(x), started at n=1 */
    register long double Tn = x - 2.0L;	/* T_{n}(x) , started at n=1 */
    register long double resul = Kncoe[0] + Kncoe[1] * Tn;

    x -= 2.0L;

    for (int n = 2; n < sizeof(Kncoe) / sizeof(long double); n++)
    {
      const long double Tn1 = 2.0L * x * Tn - Tn_1;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
      resul += Kncoe[n] * Tn1;
      Tn_1 = Tn;
      Tn = Tn1;
    }
    return resul;
  }
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
    }
    else if (ident_ == 1)
    {
      // Gaussian with sigma 10
      PDF = [&, this](double x) { return (1. / sqrt(2. * pi * 10. * 10.)) * exp(-1.*(x*x) / double(2. * 10. * 10.)); };
      maximumOfPDF = 1. / (10. * sqrt(2. * pi));
      analyticEntropy_ = log(10. * sqrt(2. * pi * e));
    }
    else if (ident_ == 2)
    {
      // Parabola normated to integrate to 1 over [0;1]
      PDF = [&, this](double x) { if (x < -1. || x > 1.) return 0.; else return (3. / 4.) * (-1. * ((x)*(x)) + 1.); };
      maximumOfPDF = PDF(0.);
      analyticEntropy_ = 0.568054;
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
    bool run = true;
    double PDFrange = 0.;
    while (run)
    {
      if (PDF(PDFrange) < 0.000000001) break;
      else PDFrange += 0.005;
    }
    std::uniform_real_distribution<double> unifDistrRange(-PDFrange, PDFrange);


    for (int n = 0; n < numberOfSamples; ++n)
    {
      double drawUnif = unifDistr(gen);
      double drawUnifWithRange = unifDistrRange(gen);
      double M = maximumOfPDF / (1. / (2.*PDFrange));
      if (drawUnif < PDF(drawUnifWithRange) / (M * 1.5 * 1. / (2.*PDFrange)))
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

  int ident() { return this->identif; };

  std::string identString() 
  {
    if (this->ident() == 0) return "Normal sig=1";
    else if (this->ident() == 1) return "Normal sig=10";
    else if (this->ident() == 2) return "parabolic";
    else return "ERROR_NAME";
  }

  double analyticEntropy() {
    return this->analyticEntropy_;
  };

  std::function<double(double x)> function() {
    return this->PDF;
  };

  double MCIntegrationEntropy(double range, size_t numberOfSamples)
  {
    // MC guess with uniform distribution
    // Draw Uniform and do
    // Draw samples
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distr(-1. * range, range);
    std::vector<double> temp;
    temp.reserve(numberOfSamples);
    for (int n = 0; n < numberOfSamples; ++n) {
      double drawnnum = distr(gen);
      temp.push_back(drawnnum);
    }



    long double uniformMCvalue = 0, uniformMCvalue2 = 0, c = 0.f;

    // Debug
    //Sort samples
    //std::stable_sort(temp.begin(), temp.end());
    //size_t countHowManyContinuesBecauseOfZeroPDF = 0u;

    for (int n = 0; n < numberOfSamples; ++n)
    {
      long double gvalue;
      gvalue = PDF(temp[n]);
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
    uniformMCvalue *= -1. * (2*range);
    uniformMCvalue2 *= -1. * (2*range);
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
      double temp1 = PDF(samples[i]);
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
    std::vector<double> samples = this->draw(numberOfSamples);
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

  double meaningfulRange()
  {
    //Target PDF in 0.1er Schritten auswerten und schauen wann es kleiner ist
    bool run = true;
    double PDFrange = 0.;
    while (run)
    {
      if (PDF(PDFrange) < 0.00000001 || PDF(-1.* PDFrange) < 0.00000001) break;
      else PDFrange += 0.001;
    }
    return PDFrange;
  }

private:
  double analyticEntropy_;
  double maximumOfPDF;
  int identif;
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
    drawAndEvaluateMatrix(iter_, 4),
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

      for (int n = 0; n < numberOfDraws; ++n)
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

  double empiricalGaussianEntropy()
  {
    // Calculate Mean
    for (int n = 0; n < numberOfDraws; ++n)
      mean += drawAndEvaluateMatrix(n, 0);
    mean /= numberOfDraws;

    // Calculate standard Deviation
    for (int n = 0; n < numberOfDraws; ++n)
      standardDeviation += (drawAndEvaluateMatrix(n, 0) - mean) * (drawAndEvaluateMatrix(n, 0) - mean);
    standardDeviation /= numberOfDraws;
    standardDeviation = sqrt(standardDeviation);

    return log(standardDeviation * sqrt(2. * pi * e));
  }
};


// Calculated entropy object calculates estiamted entropy
class calculatedentropyobj : public entropyobj
{
public:
  size_t k;
  double calculatedEntropyHnizdo;
  double calculatedEntropyLombardi; 
  double mcintegrationEntropy;
  double mcdrawEntropy;
  double analyticalEntropy;
  double empiricalNormalDistributionEntropy;
  calculatedentropyobj(size_t k_, entropyobj const& obj, double analyticEntropyValue = 0.) : 
    entropyobj(obj), 
    k(k_), 
    calculatedEntropyHnizdo(std::numeric_limits<double>::quiet_NaN()),
    calculatedEntropyLombardi(std::numeric_limits<double>::quiet_NaN()),
    mcintegrationEntropy(std::numeric_limits<double>::quiet_NaN()),
    mcdrawEntropy(std::numeric_limits<double>::quiet_NaN()),
    analyticalEntropy(this->probdens.analyticEntropy()),
    empiricalNormalDistributionEntropy(std::numeric_limits<double>::quiet_NaN())
  {
  }

  void calculate()
  {
    // Calculates entropy inegral using MC 
    // with known PDF
    mcdrawEntropy = this->probdens.MCDrawEntropy(drawAndEvaluateMatrix);

    mcintegrationEntropy = this->probdens.MCIntegrationEntropy(this->probdens.meaningfulRange(), numberOfDraws);

    empiricalNormalDistributionEntropy = this->empiricalGaussianEntropy();

    // Calculate Hnizdo as well as Lombardi/Pant entropy
    if (Config::get().entropytrails.NNcalculation)
    {
      // Matrix Layout after calculation:
      // First col: Drawn samples
      // Second col: k / iter * kNNdistance
      // Third col: Value of real distribution
      // Fourth col: kNN Distance

      //Matrix_Class addition(drawAndEvaluateMatrix.rows(), 3);
      //drawAndEvaluateMatrix.append_right(addition);

      //Neccessarry
      transpose(drawAndEvaluateMatrix);
      Matrix_Class copytemp = drawAndEvaluateMatrix;
      Matrix_Class drawAndEvaluateMatrix_TemporaryCopy = drawAndEvaluateMatrix;
      std::function<double(double x)> PDFtemporary = this->probdens.function();

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
          float_type holdNNdistance = sqrt(matop::entropy::knn_distance(copytemp, 1, k, 0u, i, buffer));
          drawAndEvaluateMatrix_TemporaryCopy(1, i) = double(k) / double(numberOfDraws) / holdNNdistance;
          drawAndEvaluateMatrix_TemporaryCopy(2, i) = PDFtemporary(copytemp(0, i));
          drawAndEvaluateMatrix_TemporaryCopy(3, i) = holdNNdistance;
          // Lombardi kPN
        }
#ifdef _OPENMP
      }
#endif

     drawAndEvaluateMatrix = drawAndEvaluateMatrix_TemporaryCopy;


      // ENTROPY according to Hnzido
      double tempsum = 0.;
      for (size_t i = 0u; i < drawAndEvaluateMatrix.cols(); i++)
        tempsum += log(drawAndEvaluateMatrix(3, i));

      tempsum /= double(numberOfDraws);
      tempsum *= double(dimension);

      // Einschub
      // Entropy according to Lombardi (traditional)
      double tempsum_kpn_alternative = tempsum + (log(pow(pi, double(dimension) / 2.)) / (tgamma(0.5 * dimension + 1)));
      tempsum_kpn_alternative += digammal(numberOfDraws);
      tempsum_kpn_alternative -= digammal(k);

      tempsum += (log(numberOfDraws * pow(pi, double(dimension) / 2.)) / (tgamma(0.5 * dimension + 1)));
      double tempsum2 = 0;
      if (k != 1)
      {
        for (size_t i = 1; i < k; i++)
        {
          tempsum2 += 1.0 / double(i);
        }
      }
      tempsum -= tempsum2;
      tempsum += 0.5772156649015328606065;

      //////////////
      calculatedEntropyHnizdo = tempsum; // Hnizdo Entropy
      calculatedEntropyLombardi = tempsum_kpn_alternative; // Lombardi Entropy

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
    if (!myfile3.good()) writeHeader = true;
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
      myfile2 << std::setw(16) << std::scientific << std::setprecision(5) << "Empirical Gauss|";

      myfile2 << "\n==============================================================================================================================\n";
    }
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->probdens.identString() << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->numberOfDraws << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->k << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->probdens.analyticEntropy() << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->mcintegrationEntropy << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->mcdrawEntropy << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyHnizdo << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->calculatedEntropyLombardi << "|";
    myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << this->empiricalNormalDistributionEntropy << "|\n";
    myfile2.close();
  }
};