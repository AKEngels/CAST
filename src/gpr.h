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
#include <optional>
#include <boost/container/small_vector.hpp>
#include "Scon/scon_mathmatrix.h"

namespace gpr {
  constexpr std::size_t pes_dimensions = 2;
  using PES_Point = boost::container::small_vector<double, pes_dimensions>;

  class CovarianceFunction;

  class GPR_Interpolator {
  public:
    /**
     * Construct a Gaussian Process
     * @param kf
     * @param training_points
     * @param training_values
     * @param training_gradients
     */
    GPR_Interpolator(std::unique_ptr<CovarianceFunction> cf,
                     std::vector<PES_Point> training_points,
                     std::optional<std::pair<std::vector<double>, double>> const &training_values,
                     std::optional<std::pair<std::vector<PES_Point>, double>> const& training_gradients = std::nullopt);

    double interpolate(PES_Point const& x) const;

    /**
     * Calculates the interpolated derivative wrt. the @param d-th component of @param x at @param x
     */
    double interpolate_derivative(PES_Point const& x, std::size_t d) const;
    std::vector<double> const& get_weights() const;

  private:
    std::unique_ptr<CovarianceFunction> covarianceFunc_;

    std::vector<PES_Point> training_points_;
    std::vector<double> weights_;
    double y_prior_ = 0;
    bool has_values_;
    bool has_derivatives_;

    /**
     * Train the Gaussian Process
     * @param training_values
     * @param training_gradients
     */
    void train_gp(std::optional<std::pair<std::vector<double>, double>> const& training_values,
                  std::optional<std::pair<std::vector<PES_Point>, double>> const& training_gradients);
  };

  GPR_Interpolator gpr_interpolator_1d(std::unique_ptr<CovarianceFunction> kf,
                                       std::vector<double> const& training_points,
                                       std::vector<double> const& training_data,
                                       std::optional<std::vector<double>> const& training_gradients = std::nullopt);
  GPR_Interpolator gpr_interpolator_2d(std::unique_ptr<CovarianceFunction> kf,
                                       std::vector<std::pair<double, double>> const& training_points,
                                       std::vector<double> const& training_data);

  /**
   * Calculates the squared euclidean norm of @param x - @param y
   */
  inline auto r(PES_Point const& x, PES_Point const& y) {
    return std::inner_product(x.begin(), x.end(), y.begin(), double(0), std::plus{},
                              [](auto a, auto b){return std::pow(a-b, 2);});
  }

  class CovarianceFunction {
  public:
    /**
     * Evaluates the covariance function between points @param x and @param y
     */
    virtual double evaluate(PES_Point const& x, PES_Point const& y) const = 0;

    /**
     * Evaluates the first derivative of the covariance function wrt. the @param i-th component of @param x
     * at points @param x and @param y
     */
    virtual double first_der_x(PES_Point const& x, PES_Point const& y, std::size_t i) const = 0;

    /**
     * Evaluates the first derivative of the covariance function wrt. the @param i-th component of @param y
     * at points @param x and @param y
     */
    virtual double first_der_y(PES_Point const& x, PES_Point const& y, std::size_t i) const {
      return -first_der_x(x, y, i);
    }

    /**
     * Evaluates the second derivative of the covariance function wrt. the @param i-th component of \param x
     * and the @param j-th component of @param y at points @param x and @param y
     */
    virtual double second_der(PES_Point const& x, PES_Point const& y, std::size_t i, std::size_t j) const = 0;
  };

  /**
   * Squared-Exponential covariance function
   * k(x, y) = exp(-|x-y|²/(2l²))
   */
  class SqExpCovariance: public CovarianceFunction {
  public:
    SqExpCovariance(double l): l_{l}{}

    double evaluate(PES_Point const& x, PES_Point const&y) const final {
      return std::exp(-r(x, y) / (2*l_*l_));
    }

    double first_der_x(PES_Point const& x, PES_Point const& y, std::size_t i) const final {
      return (y[i]-x[i]) * evaluate(x, y) / (l_*l_);
    }

    double second_der(PES_Point const& x, PES_Point const& y, std::size_t i, std::size_t j) const final {
      auto delta_ij = i == j;
      return (delta_ij * l_*l_ - (x[i]-y[i])*(x[j]-y[j])) * evaluate(x, y) / std::pow(l_, 4);
    }

  private:
    double l_;
  };

  /**
   * Periodic covariance function (currently only for 1-dimensional interpolation)
   * k(x, y) = exp(-2*(sin(c(x-y))*(c*l))²) where c = 2π/period
   */
  class PeriodicCovariance: public CovarianceFunction {
  public:
    PeriodicCovariance(double l, double period): l_{l}, c_(SCON_2PI / period){}

    double evaluate(PES_Point const& x, PES_Point const& y) const final {
      assert(x.size() == y.size());
      auto arg = std::inner_product(x.begin(), x.end(), y.begin(), 0.0, std::plus{}, [c = this->c_](double a, double b) {
        return std::pow(std::sin((a-b)*c/2), 2);
      });
      return std::exp(-2* arg / std::pow(l_*c_, 2));
    }

    double first_der_x(PES_Point const& x, PES_Point const& y, std::size_t i) const final {
      assert(x.size() == y.size());
      auto arg = (x[i] - y[i]) * c_ / 2;
      return -2 * c_ * sin(arg) * cos(arg) * evaluate(x, y) / std::pow(l_ * c_, 2);
    }

    double second_der(PES_Point const& x, PES_Point const& y, std::size_t i, std::size_t j) const final {
      assert(x.size() == y.size());
      auto l_conv = l_ * c_;
      auto factor = [&x, &y, i, j, l_conv, c = this->c_]{
        auto arg = (x[i] - y[i]) * c;
        if (i == j)
          return l_conv * l_conv * cos(arg) - std::pow(sin(arg), 2);
        else
          return (cos(c*(x[i] + x[j] - y[i] - y[j])) - cos(c*(x[i] - x[j] - y[i] + y[j])))/2;
      }();
      return c_ * c_ * factor * evaluate(x, y) / std::pow(l_conv, 4);
    }

  private:
    double l_, c_;
  };

  /**
   * Matérn covariance function for ν = 5/2
   * k(x, y) = (1 + r + r²/(3l)) * exp(-r) where r = sqrt(5) * |x-y| / l
   */
  class MaternCovariance: public CovarianceFunction{
  public:
    MaternCovariance(double l): l_{l}{}

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

  void run_gpr_test();

}
