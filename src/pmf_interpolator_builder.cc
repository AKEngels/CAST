#include "pmf_interpolator_builder.h"

#include "helperfunctions.h"

template<typename Input>
pmf_ic::InterpolatorResult<Input> pmf_ic::build_interpolator(std::vector<Input> const& x, std::vector<double> const& y, std::optional<std::vector<Input>> const& grads) {
  // Compile time stuff
  constexpr bool interpolation_1d = std::is_same_v<Input, double>;
  using GPR_Type = std::conditional_t<interpolation_1d, GPRInterpolator1D, GPRInterpolator2D>;
  using SplineType = std::conditional_t<interpolation_1d, Spline1DInterpolator, Spline2DInterpolator>;

  // Run-time checks
  if (x.size() != y.size())
    throw std::runtime_error("Size mismatch between training inputs and training data");
  using ICModes = config::coords::umbrellas::pmf_ic_conf::ic_mode;
  auto interpolation_mode = Config::get().coords.umbrella.pmf_ic.mode;
  if (interpolation_mode == ICModes::OFF)
    throw std::runtime_error("No PMF-IC mode selected. Check CAST.txt");

  auto dimensionality = get_interpolation_dimensionality();
  if constexpr(interpolation_1d)
    assert(dimensionality == 1);
  else
    assert(dimensionality == 2);

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

    return std::make_unique<GPR_Type>(std::move(kernel), x, y, grads);
  }
}

template std::unique_ptr<pmf_ic::Interpolator1DInterface> pmf_ic::build_interpolator(std::vector<double> const&, std::vector<double> const&, std::optional<std::vector<double>> const&);
template std::unique_ptr<pmf_ic::Interpolator2DInterface> pmf_ic::build_interpolator(std::vector<std::pair<double, double>> const&, std::vector<double> const&, std::optional<std::vector<std::pair<double, double>>> const&);

pmf_ic::Interpolator pmf_ic::load_interpolation() {
  if (!file_exists(Config::get().coords.umbrella.pmf_ic.prepfile_name))
    throw std::runtime_error("File for getting spline function not found.");

  std::ifstream input(Config::get().coords.umbrella.pmf_ic.prepfile_name);
  std::string line;
  std::vector<std::string> linestr;
  std::getline(input, line);        // discard first line

  try {
    if (Config::get().coords.umbrella.pmf_ic.indices_xi.size() == 1)   // one-dimensional
    {
      std::vector<double> xis;
      std::vector<double> deltaEs;

      while (!input.eof())
      {
        std::getline(input, line);
        if (line == "") break;
        linestr = split(line, ',');
        xis.emplace_back(std::stod(linestr.at(0)));
        deltaEs.emplace_back(std::stod(linestr.at(3)));
      }
      return build_interpolator(xis, deltaEs);
    }

    else if (Config::get().coords.umbrella.pmf_ic.indices_xi.size() == 2)          // two-dimensional
    {
      std::vector<std::pair<double, double>> xis;
      std::vector<double> deltaEs;

      while (!input.eof())
      {
        std::getline(input, line);
        if (line == "") break;
        linestr = split(line, ',');
        double xi1 = std::stod(linestr.at(0));
        double xi2 = std::stod(linestr.at(1));
        xis.emplace_back(std::make_pair(xi1, xi2));
        deltaEs.emplace_back(std::stod(linestr.at(4)));
      }
      return build_interpolator(xis, deltaEs);
    }
    else
      throw std::runtime_error("Interpolated corrections enabled but no CVs specified. Check CAST.txt");
  } catch (std::out_of_range&) {
    throw std::runtime_error("Wrong number of columns in PMF-IC preparation file");
  }
}
