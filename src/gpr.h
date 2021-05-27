/**
CAST 3
gpr.h
Purpose: Facilities for GPR Interpolation (as alternatives to splines)

@author Christian Schärf
*/

#pragma once

#include <cmath>
#include <functional>
#include <memory>
#include <boost/container/small_vector.hpp>
#include "Scon/scon_mathmatrix.h"

namespace gpr {
  constexpr std::size_t pes_dimensions = 2;
  using PES_Point = boost::container::small_vector<double, pes_dimensions>;

  class KernelFunction;

  class GPR_Interpolator {
  public:
    /**
     * Construct a Gaussian Process
     * @param kf
     * @param training_points
     * @param training_data
     */
    GPR_Interpolator(std::unique_ptr<KernelFunction> kf,
                     std::vector<PES_Point> training_points,
                     std::vector<double> const &training_data);

    double interpolate(PES_Point const& x) const;
    std::vector<double> const& get_weights() const;

  private:
    std::unique_ptr<KernelFunction> kernel_;

    std::vector<PES_Point> training_points_;
    std::vector<double> weights_;
    double y_prior_ = 0;

    /**
     * Train the Gaussian Process
     * @param training_data
     */
    void train_gp(std::vector<double> const& training_data);
  };

  GPR_Interpolator gpr_interpolator_1d(std::unique_ptr<KernelFunction> kf,
                                       std::vector<double> const& training_points,
                                       std::vector<double> const& training_data);
  GPR_Interpolator gpr_interpolator_2d(std::unique_ptr<KernelFunction> kf,
                                       std::vector<std::pair<double, double>> const& training_points,
                                       std::vector<double> const& training_data);

  /**
   * Calculates the squared euclidean norm of @param x - @param y
   */
  inline auto r(PES_Point const& x, PES_Point const& y) {
    return std::inner_product(x.begin(), x.end(), y.begin(), double(0), std::plus{},
                              [](auto a, auto b){return std::pow(a-b, 2);});
  }

  class KernelFunction {
  public:
    /**
     * Evaluates the kernel function between points @param x and @param y
     */
    virtual double evaluate(PES_Point const& x, PES_Point const& y) const = 0;

    /**
     * Evaluates the first derivative of the kernel function wrt. the @param i-th component of @param x
     * at points @param x and @param y
     */
    virtual double first_der_x(PES_Point const& x, PES_Point const& y, std::size_t i) const = 0;

    /**
     * Evaluates the first derivative of the kernel function wrt. the @param i-th component of @param y
     * at points @param x and @param y
     */
    virtual double first_der_y(PES_Point const& x, PES_Point const& y, std::size_t i) const {
      return -first_der_x(x, y, i);
    }

    /**
     * Evaluates the second derivative of the kernel function wrt. the @param i-th component of \param x
     * and the @param j-th component of @param y at points @param x and @param y
     */
    virtual double second_der(PES_Point const& x, PES_Point const& y, std::size_t i, std::size_t j) const = 0;
  };

  /**
   * Squared-Exponential kernel function
   * k(x, y) = exp(-l*|x-y|²/2
   */
  class SqExpKernel: public KernelFunction {
  public:
    SqExpKernel(double l): l_{l}{}

    double evaluate(PES_Point const& x, PES_Point const&y) const final {
      return std::exp(-l_*r(x, y) / 2);
    }

    double first_der_x(PES_Point const& x, PES_Point const& y, std::size_t i) const final {
      return l_ * (y[i]-x[i]) * evaluate(x, y);
    }

    double second_der(PES_Point const& x, PES_Point const& y, std::size_t i, std::size_t j) const final {
      auto delta_ij = i == j;
      return l_ * (delta_ij - l_ * (x[i]-y[i])*(x[j]-y[j])) * evaluate(x, y);
    }

  private:
    double l_;
  };

  /**
   * Matérn kernel function for ν = 5/2
   * k(x, y) = (1 + r + r²/l) * exp(-r) where r = sqrt(5) * |x-y| / l
   */
  class MaternKernel: public KernelFunction{
  public:
    MaternKernel(double l): l_{l}{}

    double evaluate(PES_Point const& x, PES_Point const& y) const final {
      double a = arg(x, y);
      return (1 + a + a*a/3) * std::exp(-a);
    }

    double first_der_x(PES_Point const& x, PES_Point const& y, std::size_t i) const final {
      auto a = arg(x, y);
      return 5/(3*l_*l_) * (y[i] - x[i]) * (1 + a) * std::exp(-a);
    }

    double second_der(PES_Point const& x, PES_Point const& y, std::size_t i, std::size_t j) const final {
      auto a = arg(x, y);
      return 5/(3*std::pow(l_, 4)) * ((i == j) * l_*l_ * (1 + a) - 5 * (x[i]-y[i]) * (x[j]-y[j])) * std::exp(-a);
    }

  private:
    double l_;

    double arg(const PES_Point &x, const PES_Point &y) const {
      return sqrt(5 * r(x, y)) / l_;
    }
  };

  inline auto matern_kernel(double l) {
    return [l](PES_Point const& x, PES_Point const& y) {
      auto arg = std::sqrt(5* r(x, y))/l;
      return (1 + arg + arg*arg/3) * std::exp(-arg);
    };
  }

  void run_gpr_test();

}
