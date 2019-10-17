#include<stdexcept>
#include "spline.h"
#include "configuration.h"
#include "helperfunctions.h"
#include "Scon/scon_angle.h"

Spline::Spline(std::vector<double> const& x_values, std::vector<double> const& y_values)
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

double Spline::get_value(double x) const
{
  return spline1dcalc(spline, x);
}

std::vector<double> Spline::get_value_and_derivative(double x) const
{
  double value, derivative, second;                       
  spline1ddiff(spline, x, value, derivative, second);
  return { value, derivative };
}

double mapping::xi_to_z(double xi)
{
  auto const& xi_0 = Config::get().coords.umbrella.pmf_ic.xi0;
  auto const& L = Config::get().coords.umbrella.pmf_ic.L;
  auto z = (2.0 / SCON_PI) * atan((xi - xi_0) / L);
  return z;
}

double mapping::dz_dxi(double xi)
{
  auto const& xi_0 = Config::get().coords.umbrella.pmf_ic.xi0;
  auto const& L = Config::get().coords.umbrella.pmf_ic.L;
  auto res = (2 * L) / (SCON_PI * (xi_0 * xi_0 - 2 * xi_0 * xi + L * L + xi * xi));
  return res;
}