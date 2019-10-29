#include<stdexcept>
#include "spline.h"
#include "configuration.h"
#include "helperfunctions.h"
#include "Scon/scon_angle.h"

void Spline1D::fill(std::vector<double> const& x_values, std::vector<double> const& y_values)
{
  if (x_values.size() != y_values.size())
    throw std::runtime_error("Error in creating spline: x-values and y-values don't have the same size.");
  else if (double_element(x_values)) {
    throw std::runtime_error("Error in creating spline: one x-value is there more than once.");
  }
    
  alglib::real_1d_array x, y;
  x.setcontent(x_values.size(), &x_values[0]);
  y.setcontent(y_values.size(), &y_values[0]);

  alglib::spline1dbuildcubic(x, y, spline);                    
}

double Spline1D::get_value(double const x) const
{
  return spline1dcalc(spline, x);
}

double Spline1D::get_derivative(double const x) const
{
  double value, derivative, second;                       
  spline1ddiff(spline, x, value, derivative, second);
  return value;
}

void Spline2D::fill(std::vector<std::pair<double, double>> const& x_values, std::vector<double> const& y_values)
{
  // spline is created in a similar way as here:
  // http://www.alglib.net/translator/man/manual.cpp.html#example_spline2d_fit_blocklls

  // some checks if data is okay
  if (x_values.size() != y_values.size())
    throw std::runtime_error("Error in creating spline: x-values and y-values don't have the same size.");
  else if (double_element(x_values)) {
    throw std::runtime_error("Error in creating spline: one x-value is there more than once.");
  }

  // bring data in correct format (alglib::real_2d_array)
  std::vector <std::vector<double>> xxy_vec;
  for (auto i{ 0u }; i < x_values.size(); ++i) {
    std::vector<double> current{ {x_values[i].first, x_values[i].second, y_values[i]} };
    xxy_vec.emplace_back(current);
  }
  alglib::real_2d_array xxy;
  xxy.setlength(x_values.size(), 3);
  for (int i = 0; i < x_values.size(); i++) {
    for (int j = 0; j < 3; j++) {
      xxy(i, j) = xxy_vec[i][j];
    }
  }

  // build spline
  alglib::spline2dbuilder builder;
  alglib::spline2dbuildercreate(1, builder);
  alglib::spline2dbuildersetpoints(builder, xxy, x_values.size());
  
  alglib::spline2dfitreport rep;
  spline2dfit(builder, spline2d, rep);
}

double Spline2D::get_value(double const x1, double const x2) const
{
  return spline2dcalc(spline2d, x1, x2);
}

std::vector<double> Spline2D::get_derivatives(double const x1, double const x2) const
{
  double value, derivative1, derivative2, second;
  spline2ddiff(spline2d, x1, x2, value, derivative1, derivative2, second);
  return { derivative1, derivative2 };
}

double mapping::xi_to_z(double const xi, double const xi_0, double const L)
{
  auto z = (2.0 / SCON_PI) * atan((xi - xi_0) / L);
  return z;
}

double mapping::dz_dxi(double const xi, double const xi_0, double const L)
{
  auto res = (2 * L) / (SCON_PI * (xi_0 * xi_0 - 2 * xi_0 * xi + L * L + xi * xi));
  return res;
}