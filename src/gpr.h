/**
CAST 3
gpr.h
Purpose: Facilities for GPR Interpolation (as alternatives to splines)

@author Christian Schärf
*/

#pragma once

#include <cmath>
#include <functional>
#include "Scon/scon_mathmatrix.h"

namespace gpr {
  using PES_Point = double;
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

  inline auto exponential_kernel(double l) {
    return [l](PES_Point x, PES_Point y) {
      return std::exp(-std::pow(x-y, 2) / (2*l*l));
    };
  }

  void run_gpr_test();

}
