/**
CAST 3
spline.h
Purpose: wrapper for spline fitting (needed for PMF-IC, see https://doi.org/10.1021/jp049633g)
         alglib library is used for spline (normal cubic spline or bicubic spline)

@author Susanne Sauer
@version 1.0
*/

#pragma once
#include<algorithm>
#include<utility>
#include<variant>
#include<src/interpolation.h> // from alglib

#include "pmf_interpolator.h"

/**wrapper for 1D spline from alglib*/
class Spline1D : public PmfInterpolator1DInterface
{
private:

  /**alglib spline*/
  alglib::spline1dinterpolant spline;

public:

  /**create spline from values
  @param x_values: x-values
  @param y_values: y-values*/
  virtual void fill(std::vector<double> const& x_values, std::vector<double> const& y_values);

  /**get value of spline function at given point
  @param x: x-value for which value is returned*/
  virtual double get_value(double x) const;

  /**get value and derivative of spline function at given point
  @param x: x-value for which the derivative is given*/
  virtual double get_derivative(double x) const;
};

/**wrapper for 2D spline from alglib*/
class Spline2D : public PmfInterpolator2DInterface
{
private:

  /**alglib spline*/
  alglib::spline2dinterpolant spline;

public:

  /**create spline from values
  @param x_values: x-values (as pairs <x1,x2>)
  @param y_values: y-values*/
  virtual void fill(std::vector<std::pair<double, double>> const& x_values, std::vector<double> const& y_values);

  /**get value of spline function at given point
  @param x: x-values for which value is returned (as std::pair)*/
  virtual double get_value(double x, double y) const;

  /**get derivatives of spline function at given point
  @param x: x-values for which value is returned (as std::pair)
  returns a pair where first value is derivative in direction x1 and second value is derivative in direction x2*/
  virtual std::pair<double,double> get_derivative(double x, double y) const;
};

/**
 * Mapper to map xi to z
 * for the purpose see https://doi.org/10.1021/jp049633g
 */
class XiToZMapper {
public:
  XiToZMapper(double xi_0, double L): xi_0_{xi_0}, L_{L} {}

  /** map given xi to zeta */
  double map(double xi) const;

  /** calculate derivative dz/dxi for given xi */
  double dz_dxi(double xi) const;

private:
  double xi_0_, L_;
};

/// Wrappers around the respective Spline classes that also perform the xi to z mapping

class Spline1DInterpolator: public Spline1D
{
public:
  Spline1DInterpolator(XiToZMapper const& mapper): mapper_{mapper}{}

  void fill(const std::vector<double> &xi_values, const std::vector<double> &y_values) final {
    std::vector<double> z_values;
    z_values.reserve(xi_values.size());
    std::transform(xi_values.begin(), xi_values.end(), std::back_inserter(z_values),
                   [this](double xi){return mapper_.map(xi);});

    Spline1D::fill(z_values, y_values);
  }

  double get_value(double xi) const final {
    return Spline1D::get_value(mapper_.map(xi));
  }

  double get_derivative(double xi) const final {
    return mapper_.dz_dxi(xi) * Spline1D::get_derivative(mapper_.map(xi));
  }

private:
  XiToZMapper mapper_;
};

class Spline2DInterpolator: public Spline2D
{
public:
  Spline2DInterpolator(XiToZMapper const& mapper1, XiToZMapper const& mapper2):
    mapper1_{mapper1},
    mapper2_{mapper2}
  {}

  void fill(const std::vector<std::pair<double, double>> &xi_values, const std::vector<double> &y_values) override {
    std::vector<std::pair<double, double>> z_values;
    z_values.reserve(xi_values.size());
    std::transform(xi_values.begin(), xi_values.end(), std::back_inserter(z_values),
                   [this](auto xi){return std::make_pair(mapper1_.map(xi.first), mapper2_.map(xi.second));});

    Spline2D::fill(z_values, y_values);
  }

  double get_value(double x, double y) const final {
    return Spline2D::get_value(mapper1_.map(x), mapper2_.map(y));
  }

  std::pair<double, double> get_derivative(double x, double y) const final {
    auto [dz1, dz2] = Spline2D::get_derivative(mapper1_.map(x), mapper2_.map(y));

    return {dz1*mapper1_.dz_dxi(x), dz2*mapper2_.dz_dxi(y)};
  }

private:
  XiToZMapper mapper1_, mapper2_;
};
