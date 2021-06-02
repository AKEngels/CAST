/**
CAST 3
interpolator.h
Purpose: Factories for interpolators for PMF-IC

@author Christian Schärf
*/

#pragma once
#include <memory>
#include <utility>

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
