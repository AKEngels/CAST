/**
CAST 3
pmf_interpolator.h
Purpose: Interpolator interfaces for PMF-IC

@author Christian Schärf
*/

#pragma once
#include <memory>
#include <utility>

#include "gpr.h"

class PmfInterpolator1DInterface {
public:
  virtual double get_value(double x) const = 0;
  virtual double get_derivative(double x) const = 0;
};

class PmfInterpolator2DInterface {
public:
  virtual double get_value(double x, double y) const = 0;
  virtual std::pair<double, double> get_derivative(double x, double y) const = 0;
};

using PmfInterpolator = std::variant<std::unique_ptr<PmfInterpolator1DInterface>,
                                     std::unique_ptr<PmfInterpolator2DInterface>>;

/** wrapper classes around gpr::GPR_Interpolator for PMF IC */
class GPRInterpolator1D: public PmfInterpolator1DInterface, private gpr::GPR_Interpolator {
public:
  template <typename ...Args>
  GPRInterpolator1D(Args&& ...args):
    gpr::GPR_Interpolator(gpr::gpr_interpolator_1d(std::forward<Args>(args)...))
  {}

  double get_value(double x) const override {
    return gpr::GPR_Interpolator::interpolate({x});
  }

  double get_derivative(double x) const override {
    return gpr::GPR_Interpolator::interpolate_derivative({x}, 0);
  }
};

class GPRInterpolator2D: public PmfInterpolator2DInterface, private gpr::GPR_Interpolator {
public:
  template <typename ...Args>
  GPRInterpolator2D(Args&& ...args):
  gpr::GPR_Interpolator(gpr::gpr_interpolator_2d(std::forward<Args>(args)...))
  {}

  double get_value(double x, double y) const override {
    return gpr::GPR_Interpolator::interpolate({x, y});
  }

  std::pair<double, double> get_derivative(double x, double y) const override {
    return {gpr::GPR_Interpolator::interpolate_derivative({x, y}, 0),
            gpr::GPR_Interpolator::interpolate_derivative({x, y}, 1)};
  }
};
