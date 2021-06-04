#include "pmf_interpolator_builder.h"

#include "helperfunctions.h"

template<typename Input>
pmf_ic::InterpolatorResult<Input> pmf_ic::build_interpolator(std::vector<Input> const& x, std::vector<double> const& y) {
  // Compile time stuff
  constexpr bool interpolation_1d = std::is_same_v<Input, double>;
  using GPR_Type = std::conditional_t<interpolation_1d, GPRInterpolator1D, GPRInterpolator2D>;
  using SplineType = std::conditional_t<interpolation_1d, Spline1DInterpolator, Spline2DInterpolator>;

  // Run-time checks
  if (x.size() != y.size())
    throw std::runtime_error("Size mismatch between training inputs and training data");

  auto dimensionality = get_interpolation_dimensionality();
  if constexpr(interpolation_1d)
    assert(dimensionality == 1);
  else
    assert(dimensionality == 2);

  auto interpolation_mode = Config::get().coords.umbrella.pmf_ic.mode;
  using ICModes = config::coords::umbrellas::pmf_ic_conf::ic_mode;
  assert(interpolation_mode != ICModes::OFF);
  if (interpolation_mode == ICModes::SPLINE) {
    auto const& xi0 = Config::get().coords.umbrella.pmf_ic.xi0;
    auto const& L = Config::get().coords.umbrella.pmf_ic.L;
    if (xi0.size() != dimensionality || L.size() != dimensionality)
      throw std::runtime_error("Wrong number of xi0 or L parameters for PMF-IC. Check CAST.txt");

    if constexpr(interpolation_1d)
      return std::make_unique<SplineType>(XiToZMapper(xi0[0], L[0]), x, y);
    else
      return std::make_unique<SplineType>(XiToZMapper(xi0[0], L[0]),
                                          XiToZMapper(xi0[1], L[1]), x, y);
  } else {
    auto kernel = [interpolation_mode]() -> std::unique_ptr<gpr::KernelFunction> {
      auto l = Config::get().coords.umbrella.pmf_ic.gpr_hyperparameter;
      if (interpolation_mode == ICModes::GPR_SQEXP)
        return std::make_unique<gpr::SqExpKernel>(l);
      else
        return std::make_unique<gpr::MaternKernel>(l);
    }();

    return std::make_unique<GPR_Type>(std::move(kernel), x, y);
  }
}

template std::unique_ptr<pmf_ic::Interpolator1DInterface> pmf_ic::build_interpolator(std::vector<double> const& x, std::vector<double> const& y);
template std::unique_ptr<pmf_ic::Interpolator2DInterface> pmf_ic::build_interpolator(std::vector<std::pair<double, double>> const& x, std::vector<double> const& y);
