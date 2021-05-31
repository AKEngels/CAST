/**
CAST 3
spline.h
Purpose: wrapper for spline fitting (needed for PMF-IC, see https://doi.org/10.1021/jp049633g)
         alglib library is used for spline (normal cubic spline or bicubic spline)

@author Susanne Sauer
@version 1.0
*/

#pragma once
#include<utility>
#include<variant>
#include<src/interpolation.h> // from alglib

/**wrapper for 1D spline from alglib*/
class Spline1D
{
private:

  /**alglib spline*/
  alglib::spline1dinterpolant spline;

public:

  /**create spline from values
  @param x_values: x-values
  @param y_values: y-values*/
  void fill(std::vector<double> const& x_values, std::vector<double> const& y_values);

  /**get value of spline function at given point
  @param x: x-value for which value is returned*/
  double get_value(double const x) const;

  /**get value and derivative of spline function at given point
  @param x: x-value for which the derivative is given*/
  double get_derivative(double const x) const;
};

/**wrapper for 2D spline from alglib*/
class Spline2D
{
private:

  /**alglib spline*/
  alglib::spline2dinterpolant spline;

public:

  /**create spline from values
  @param x_values: x-values (as pairs <x1,x2>)
  @param y_values: y-values*/
  void fill(std::vector<std::pair<double, double>> const& x_values, std::vector<double> const& y_values);

  /**get value of spline function at given point
  @param x: x-values for which value is returned (as std::pair)*/
  double get_value(std::pair<double, double> const& x) const;

  /**get derivatives of spline function at given point
  @param x: x-values for which value is returned (as std::pair)
  returns a pair where first value is derivative in direction x1 and second value is derivative in direction x2*/
  std::pair<double,double> get_derivative(std::pair<double, double> const& x) const;
};

using Spline = std::variant<Spline1D, Spline2D>;

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
