/**
CAST 3
gpr.h
Purpose: Facilities for GPR Interpolation (as alternatives to splines)

@author Christian Schärf
*/

#pragma once

#include <cmath>
#include <functional>
#include <boost/container/small_vector.hpp>
#include "Scon/scon_mathmatrix.h"

namespace gpr {
  constexpr std::size_t pes_dimensions = 2;
  using PES_Point = boost::container::small_vector<double, pes_dimensions>;

  using KernelFunction = std::function<double(PES_Point, PES_Point)>;

  class GPR_Interpolator {
  public:
    /**
     * Construct a Gaussian Process
     * @param kf
     * @param training_points
     * @param training_data
     */
    GPR_Interpolator(KernelFunction kf,
                     std::vector<PES_Point> training_points,
                     std::vector<double> const &training_data);

    double interpolate(PES_Point const& x) const;
    std::vector<double> const& get_weights() const;

  private:
    KernelFunction kernel_;

    std::vector<PES_Point> training_points_;
    std::vector<double> weights_;
    double y_prior_ = 0;

    /**
     * Train the Gaussian Process
     * @param training_data
     */
    void train_gp(std::vector<double> const& training_data);
  };

  GPR_Interpolator gpr_interpolator_1d(KernelFunction kf,
                                       std::vector<double> const& training_points,
                                       std::vector<double> const& training_data);
  GPR_Interpolator gpr_interpolator_2d(KernelFunction kf,
                                       std::vector<std::pair<double, double>> const& training_points,
                                       std::vector<double> const& training_data);

  inline auto r(PES_Point const& x, PES_Point const& y) {
    return std::inner_product(x.begin(), x.end(), y.begin(), double(0), std::plus{},
                              [](auto a, auto b){return std::pow(a-b, 2);});
  }

  inline auto exponential_kernel(double l) {
    return [l](PES_Point const& x, PES_Point const& y) {
      return std::exp(-r(x, y) / (2*l*l));
    };
  }

  inline auto matern_kernel(double l) {
    return [l](PES_Point const& x, PES_Point const& y) {
      auto arg = std::sqrt(5* r(x, y))/l;
      return (1 + arg + arg*arg/3) * std::exp(-arg);
    };
  }

  void run_gpr_test();

}
