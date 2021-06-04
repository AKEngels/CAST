/**
CAST 3
pmf_interpolator_builder.h
Purpose: Factories for PMF-IC interpolators

@author Christian Schärf
*/

#pragma once

#include "spline.h"

inline std::size_t get_interpolation_dimensionality() {
  return Config::get().coords.umbrella.pmf_ic.indices_xi.size();
}

/* What is up with this return type?
 * It goes like this:
 * 1. Only works if Input is either double or std::pair<double, double>,
 *    will fail to compile otherwise with an error message like "no type in std::enable_if<false, ...>"
 * 2. If Input is double, the type becomes std::unique_ptr<PmfInterpolator1DInterface>
 * 3. Otherwise, i.e. if Input is std::pair<double, double>, the type becomes std::unique_ptr<PmfInterpolator2DInterface>
 */
template<typename Input>
using InterpolatorResult = std::enable_if_t<
        std::is_same_v<Input, double> || std::is_same_v<Input, std::pair<double, double>>,
        std::conditional_t<std::is_same_v<Input, double>, std::unique_ptr<PmfInterpolator1DInterface>, std::unique_ptr<PmfInterpolator2DInterface>>>;

/**
 * Builds an interpolator according to configuration (spline or GPR, parameters etc.)
 * @tparam Input The input of the resulting interpolator. Can be double or std::pair<double, double> for 1d and 2d interpolators, respectively.
 * @param x Training points.
 * @param y Training data, i.e.
 */
template<typename Input>
InterpolatorResult<Input> build_interpolator(std::vector<Input> const& x, std::vector<double> const& y);
