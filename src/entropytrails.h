#pragma once
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <omp.h>
#include "matop.h"
#include <string>

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592;
const double e = 2.71828182845904523536028747135266249775724709369995;

class entropyobj
{
public:
  size_t iter, sigma, dimension;
  Matrix_Class matrix;
  entropyobj(Matrix_Class& matrix_, size_t& iter_, size_t& dimension_, double& sigma_) : matrix(matrix_), iter(iter_), dimension(dimension_), sigma(sigma_) {}
};

class calculatedentropyobj : public entropyobj
{
public:
  size_t k;
  double entropy;
  calculatedentropyobj(size_t k_, entropyobj const& obj) : entropyobj(obj), k(k_), entropy(0)
  {
    this->calculate();
  }

  void calculate()
  {
    Matrix_Class addition(matrix.rows(), 3);
    matrix.append_right(addition);
    //Neccessarry
    transpose(matrix);
    Matrix_Class copytemp = matrix;
    Matrix_Class copytemp2 = matrix;
#ifdef _OPENMP
    auto const n_omp = static_cast<std::ptrdiff_t>(matrix.cols());
#pragma omp parallel for firstprivate(copytemp) shared(copytemp2)
    for (std::ptrdiff_t i = 0; i< n_omp; ++i)
#else
    for (size_t i = 0u; i < matrix.cols(); i++)
#endif
    {
      double holdNNdistance = sqrt(matop::entropy::knn_distance(copytemp, 1, k, 0u, i));
      copytemp2(1, i) = k / double(iter) / holdNNdistance;
      copytemp2(2, i) = (1. / sqrt(2 * pi * sigma * sigma)) * exp(-1. * copytemp(0, i) * copytemp(0, i) / (2. * sigma * sigma));
      copytemp2(3, i) = holdNNdistance;
    }
    matrix = copytemp2;

    // ENTROPY
    long double tempsum = 0.;
    for (size_t i = 0u; i < matrix.cols(); i++)
    {
      tempsum += log(matrix(3, i));
    }
    tempsum /= iter;
    tempsum *= double(dimension);
    tempsum += (log(iter * pow(pi, double(dimension) / 2.)) / (tgamma(0.5 * dimension + 1)));
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
    entropy = tempsum;
    //Neccessarry
    transpose(matrix);
  }

  void writeToFile()
  {
    std::ofstream myfile;
    myfile.open(std::string("out_k" + std::to_string(k) + "_i" + std::to_string(iter) + "_s" + std::to_string(sigma) + "_d" + std::to_string(dimension) +".txt"));

    for (size_t i = 0u; i < matrix.rows(); i++)
    {
      for (size_t j = 0u; j < matrix.cols(); j++)
        myfile << std::setw(15) << std::scientific << std::setprecision(5) << matrix(i, j);

      myfile << "\n";
    }
    myfile.close();

    std::ofstream myfile2;
    myfile2.open(std::string("out_entropy_k" + std::to_string(k) + "_i" + std::to_string(iter) + "_s" + std::to_string(sigma) + "_d" + std::to_string(dimension) + ".txt"));

    myfile2 << "Hnizdo: " << entropy << "\n";
    myfile2 << "Gauss analytical: " << log(sigma * sqrt(2 * pi * e)) << "\n";
    myfile2.close();
  }


};

class entropyconfig
{

private:
  // Stacking is
  std::vector<entropyobj> matvec;
  std::vector<size_t> iter;
  std::vector<double> sigma;
  std::vector<size_t> dimension;
  std::vector<size_t> k_values;
  std::vector<calculatedentropyobj> calculatedDistributions;
public:
  entropyconfig()
  {
  }

  void readConfig(int argc, char **argv)
  {
    /*
    for (int i = 0; i < argc; i++) {
      if (std::string(argv[i]).substr(0u, 2) == "k=")
      {
        size_t lastCommaFound = 0;
        while (std::string(argv[i]).substr(lastCommaFound, std::string::npos).find(",") != std::string::npos)
        {
          k_values.push_back(size_t(std::stoi(std::string(argv[i]).substr( std::max(size_t(3u), lastCommaFound), std::string(argv[i]).find(",", lastCommaFound) - std::max(size_t(3u), lastCommaFound)))));
          lastCommaFound = std::string(argv[i]).find(",", lastCommaFound) + 1;
        }
      }
      else if (std::string(argv[i]).substr(0u, 2) == "i=")
      {
          size_t lastCommaFound = 0;
          while (std::string(argv[i]).substr(lastCommaFound, std::string::npos).find(",") != std::string::npos)
          {
            iter.push_back(size_t(std::stoi(std::string(argv[i]).substr(std::max(size_t(3u), lastCommaFound), std::string(argv[i]).find(",", lastCommaFound) - std::max(size_t(3u), lastCommaFound)))));
            lastCommaFound = std::string(argv[i]).find(",", lastCommaFound) + 1;
          }
      }
      else if (std::string(argv[i]).substr(0u, 2) == "d=")
      {
          size_t lastCommaFound = 0;
          while (std::string(argv[i]).substr(lastCommaFound, std::string::npos).find(",") != std::string::npos)
          {
            dimension.push_back(size_t(std::stoi(std::string(argv[i]).substr(std::max(size_t(3u), lastCommaFound), std::string(argv[i]).find(",", lastCommaFound) - std::max(size_t(3u), lastCommaFound)))));
            lastCommaFound = std::string(argv[i]).find(",", lastCommaFound) + 1;
          }
      }
      else if (std::string(argv[i]).substr(0u, 2) == "s=")
      {
          size_t lastCommaFound = 0;
          while (std::string(argv[i]).substr(lastCommaFound, std::string::npos).find(",") != std::string::npos)
          {
            sigma.push_back(size_t(std::stod(std::string(argv[i]).substr(std::max(size_t(3u), lastCommaFound), std::string(argv[i]).find(",", lastCommaFound) - std::max(size_t(3u), lastCommaFound)))));
            lastCommaFound = std::string(argv[i]).find(",", lastCommaFound) + 1;
          }
      }
    }
    */
    k_values = Config::get().entropytrails.k;
    dimension = Config::get().entropytrails.dimension;
    sigma = Config::get().entropytrails.sigma;
    iter = Config::get().entropytrails.iteration;
  }

  void draw()
  {
    for (size_t d = 0u; d < dimension.size(); d++)
    {
      for (size_t s = 0u; s < sigma.size(); s++)
      {
        for (size_t i = 0u; i < iter.size(); i++)
        {
          Matrix_Class currentmat(iter[i]);
          // Check if draw exists
          std::ifstream myfile;
          std::string line;
          myfile.open(std::string("draw_i" + std::to_string(iter[i]) + "_s" + std::to_string(sigma[s]) + "_d" + std::to_string(dimension[d]) + ".txt"));
          if (myfile.good())
          {
            int j = 0;
            while (std::getline(myfile, line))
            {
              currentmat(j, 0) = std::stod(line);
              j++;
            }
            myfile.close();
          }
          else
          {
            myfile.close();
            // Draw samples
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<double> distr(0, sigma[s]);
            std::vector<double> temp;
            temp.reserve(iter[i]);
            for (int n = 0; n < iter[i]; ++n) {
              double drawnnum = distr(gen);
              temp.push_back(drawnnum);
            }
            //Sort samples
            std::stable_sort(temp.begin(), temp.end());
            for (int n = 0; n < iter[i]; ++n)
              currentmat(n, 0) = temp[n];


            // Write 
            std::ofstream myfile2;
            myfile2.open(std::string("draw_i" + std::to_string(iter[i]) + "_s" + std::to_string(sigma[s]) + "_d" + std::to_string(dimension[d]) + ".txt"));

            for (size_t w = 0u; w < currentmat.rows(); w++)
            {
              myfile2 << std::setw(15) << std::scientific << std::setprecision(5) << currentmat(w);
              myfile2 << "\n";
            }
            myfile2.close();
          }
          matvec.push_back(entropyobj(currentmat, iter[i], dimension[d], sigma[s]));
        }
      }
    }
  }

  void calculateall(void)
  {
    for (auto kc : k_values)
    {
      for (auto obj : matvec)
      {
        calculatedentropyobj temp(kc, obj);
        calculatedDistributions.push_back(temp);
        std::cout << "done" << std::endl;
      }
    }
  }

  void writeall(void)
  {
    for (auto obj : calculatedDistributions)
    {
      obj.writeToFile();
    }
  }

};