/**
CAST 3
spline.h
Purpose: wrapper for spline fitting (needed for PMF-IC, see https://doi.org/10.1021/jp049633g)
         alglib library is used for spline (normal cubic spline or bicubic spline)

@author Susanne Sauer
@version 1.0
*/

#pragma once
#include<vector>
#include<utility>
#include<src/interpolation.h> // from alglib

/**wrapper for 1D spline from alglib*/
class Spline
{
private:

  /**alglib spline*/
  alglib::spline1dinterpolant spline;

public:

  /**default constructor*/
  Spline() {};

  /**constructor
  @param x_values: x-values
  @param y_values: y-values, corresponding to x values*/
  Spline(std::vector<double> const& x_values, std::vector<double> const& y_values);
  
  /**get value of spline function at given point
  @double x: x-value for which value is returned*/
  double get_value(double const x) const;

  /**get value and derivative of spline function at given point
  @double x: x-value
  returns vector where first element is value and second value is derivative*/
  std::vector<double> get_value_and_derivative(double x) const;
};


/**wrapper for 2D spline from alglib*/
class Spline2D
{
private:

  /**alglib spline*/
  alglib::spline2dinterpolant spline;

public:

  /**constructor
  @param x_values: x-values (as pairs <x1, x2>)
  @param y_values: y-values, corresponding to x values*/
  Spline2D(std::vector<std::pair<double, double>> const& x_values, std::vector<double> const& y_values);

  /**get value of spline function at given point
  @double x1, x2: x-values for which value is returned*/
  double get_value(double const x1, double const x2) const;
};

/**namespace for mapping from xi to z (see https://doi.org/10.1021/jp049633g) */
namespace mapping
{
  /**simple mapping function*/
  double xi_to_z(double const xi, double const xi_0, double const L);
  /**derivative dz/dxi */
  double dz_dxi(double const xi, double const xi_0, double const L);
}