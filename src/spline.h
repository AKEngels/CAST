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


/**wrapper class for both 1D and 2D splines*/
class Spline
{
protected:

  /**alglib spline (only used for 1D)*/
  alglib::spline1dinterpolant spline;

  /**alglib spline (only used for 2D)*/
  alglib::spline2dinterpolant spline2d;

  /**dimension of the spline (can be 1 or 2)*/
  unsigned dimension;

public:

  /**constructor*/
  Spline() {};

  /**get dimension of spline*/
  unsigned get_dimension() { return dimension; };

  // functions for 1D spline

  /**get value of spline function at given point (only used for 1D)
  @double x: x-value for which value is returned*/
  virtual double get_value(double const x) const = 0;

  /**get derivative of spline function at given point (only used for 1D)
  @double x: x-value for which the derivative is given*/
  virtual double get_derivative(double const x) const = 0;

  // functions for 2D spline

  /**get value of spline function at given point (only used for 2D)
  @double x1, x2: x-values for which value is returned*/
  virtual double get_value(double const x1, double const x2) const = 0;

  /**get derivatives of 2D spline function at given point (only used for 2D)
  @double x1, x2: x-values for which the derivative is given
  returns a vector where first value is derivative in direction x1 and second value is derivative in direction x2*/
  virtual std::vector<double> get_derivatives(double const x1, double const x2) const = 0;
};


/**wrapper for 1D spline from alglib*/
class Spline1D : public Spline
{
public:

  /**constructor*/
  Spline1D() { dimension = 1; };

  /**create spline from values
  @param x_values: x-values
  @param y_values: y-values*/
  void fill(std::vector<double> const& x_values, std::vector<double> const& y_values);

  /**get value of spline function at given point
  @double x: x-value for which value is returned*/
  double get_value(double const x) const override;

  /**get value and derivative of spline function at given point
  @double x: x-value for which the derivative is given*/
  double get_derivative(double const x) const override;

  // functions that need to be defined but should never used because they are for 2D spline

  double get_value(double const x1, double const x2) const override {
    throw std::runtime_error("Error: This function should not be used for 1D spline.");
  }

  std::vector<double> get_derivatives(double const x1, double const x2) const override {
    throw std::runtime_error("Error: This function should not be used for 1D spline.");
  }
};


/**wrapper for 2D spline from alglib*/
class Spline2D : public Spline
{
public:

  /**constructor*/
  Spline2D() { dimension = 2; };

  /**create spline from values
  @param x_values: x-values (as pairs <x1,x2>)
  @param y_values: y-values*/
  void fill(std::vector<std::pair<double, double>> const& x_values, std::vector<double> const& y_values);

  /**get value of spline function at given point
  @double x1, x2: x-values for which value is returned*/
  double get_value(double const x1, double const x2) const override;

  /**get derivatives of spline function at given point
  @double x1, x2: x-values for which value is returned
  returns a vector where first value is derivative in direction x1 and second value is derivative in direction x2*/
  std::vector<double> get_derivatives(double const x1, double const x2) const override;

  // functions that need to be defined but should never used because they are for 1D spline

  double get_value(double const x) const override {
    throw std::runtime_error("Error: This function should not be used for 2D spline.");
  }

  double get_derivative(double const x) const override {
    throw std::runtime_error("Error: This function should not be used for 2D spline.");
  }
};


/**namespace for mapping from xi to z (see https://doi.org/10.1021/jp049633g) */
namespace mapping
{
  /**simple mapping function*/
  double xi_to_z(double const xi, double const xi_0, double const L);
  /**derivative dz/dxi */
  double dz_dxi(double const xi, double const xi_0, double const L);
}