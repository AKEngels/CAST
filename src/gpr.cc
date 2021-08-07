#include "gpr.h"

#include <cassert>

#include "configuration.h"

#include "pmf_interpolator_builder.h"

gpr::GPR_Interpolator::GPR_Interpolator(std::unique_ptr <gpr::CovarianceFunction> cf,
                                        std::vector<PES_Point> training_points,
                                        std::optional<std::pair<std::vector<double>, double>> const &training_values,
                                        std::optional<std::pair<std::vector<PES_Point>, double>> const& training_gradients)
  : covarianceFunc_{std::move(cf)}
  , training_points_{std::move(training_points)}
  , has_values_(training_values.has_value())
  , has_derivatives_(training_gradients.has_value())
{
  if (!(training_values.has_value() || training_gradients))
    throw std::runtime_error("Values or gradients required for GP training!");
  if (has_values_)
    assert(training_points_.size() == training_values->first.size());
  if (has_derivatives_)
    assert(training_points_.size() == training_gradients->first.size());
  train_gp(training_values, training_gradients);
}

gpr::GPR_Interpolator gpr::gpr_interpolator_1d(std::unique_ptr<CovarianceFunction> kf,
                                               const std::vector<double> &training_points,
                                               const std::vector<double> &training_data,
                                               std::optional<std::vector<double>> const& training_gradients) {
  std::vector<PES_Point> pes_points;
  pes_points.reserve(training_data.size());
  auto transformer = [](auto x) -> PES_Point {return {x};};
  std::transform(training_points.begin(), training_points.end(), std::back_inserter(pes_points),
                 transformer);

  std::optional<std::pair<std::vector<PES_Point>, double>> gradients;
  if (training_gradients) {
    gradients = std::make_pair(std::vector<PES_Point>(), 0.);
    gradients->first.reserve(training_gradients->size());
    std::transform(training_gradients->begin(), training_gradients->end(), std::back_inserter(gradients->first),
                   transformer);
  }

  return GPR_Interpolator(std::move(kf), std::move(pes_points), std::make_pair(training_data, 0.), gradients);
}

gpr::GPR_Interpolator gpr::gpr_interpolator_2d(std::unique_ptr<CovarianceFunction> kf, const std::vector<std::pair<double, double>> &training_points,
                                               const std::vector<double> &training_data) {
  std::vector<PES_Point> pes_points;
  pes_points.reserve(training_data.size());
  std::transform(training_points.begin(), training_points.end(), std::back_inserter(pes_points),
                 [](auto point) -> PES_Point {return {point.first, point.second};});

  return GPR_Interpolator(std::move(kf), std::move(pes_points), std::make_pair(training_data, 0.));
}

void gpr::GPR_Interpolator::train_gp(std::optional<std::pair<std::vector<double>, double>> const &training_values,
                                     std::optional<std::pair<std::vector<PES_Point>, double>> const& training_gradients) {
  auto n = training_points_.size();
  auto ndim = training_points_.front().size();

  auto cov_size = (has_values_ ? n : 0) + (has_derivatives_ ? n * ndim : 0);

  scon::mathmatrix<double> K(cov_size, cov_size, 0);
  if (has_values_) {
    auto sigma = training_values->second;
    for (std::size_t i=0; i<n; ++i) {
      for (std::size_t j=0; j<n; ++j) {
        auto const& x = training_points_[i];
        auto const& y = training_points_[j];

        // Upper left
        K(i, j) = covarianceFunc_->evaluate(x, y) + (i == j) * sigma * sigma;
      }
    }

    // Why is std::accumulate weired?
    for(auto y_i: training_values->first)
      y_prior_ += y_i;
    y_prior_ /= training_points_.size();
  }

  if (has_derivatives_) {
    for (std::size_t i=0; i<n; ++i) {
      for (std::size_t j=0; j<n; ++j) {
        auto const& x = training_points_[i];
        auto const& y = training_points_[j];
        auto sigma_g = training_gradients->second;

        auto start = has_values_ ? n : 0;

        for (std::size_t d_i = 0; d_i < ndim; ++d_i) {
          // Upper right
          K(i, start + j*ndim + d_i) = covarianceFunc_->first_der_y(x, y, d_i);

          // Lower left
          K(start + i*ndim + d_i, j) = covarianceFunc_->first_der_x(x, y, d_i);

          // Lower right
          for (std::size_t d_j=0; d_j<ndim; ++d_j) {
            K(start + i*ndim + d_i, start + j*ndim + d_j) = covarianceFunc_->second_der(x, y, d_i, d_j) + (i == j && d_i == d_j) * sigma_g * sigma_g;
          }
        }
      }
    }
  }
  if (Config::get().general.verbosity >= 4)
    std::cout << "Covariance matrix:\n" << K << '\n';

  if (Config::get().general.verbosity > 3)
    std::cout << "Prior value: " << y_prior_ << '\n';

  //assert(K.positive_definite_check());
  auto y = [&] {
    auto res = scon::mathmatrix<double>::zero(cov_size, 1);
    if (has_values_) {
      for (int i=0; i<training_points_.size(); ++i)
        res(i, 0) = training_values->first[i] - y_prior_;
    }

    if (has_derivatives_) {
      auto start = has_values_ ? n : 0;
      for (std::size_t i=0; i<n; ++i) {
        for (std::size_t d=0; d<ndim; ++d) {
          res(start + i*ndim + d, 0) = training_gradients->first[i][d];
        }
      }
    }
    return res;
  }();
  auto w_mat = K.solve(y);
  weights_ = scon::mathmatrix<double>(w_mat).col_to_std_vector();

  if (Config::get().general.verbosity >= 4)
    std::cout << "y vector:\n" << y << "\nWeights:\n" << w_mat << '\n';

  assert(weights_.size() == cov_size);
}

double gpr::GPR_Interpolator::interpolate(const gpr::PES_Point &x) const {
  double res = y_prior_;
  auto ndim = training_points_.front().size();

  if (has_values_) {
    for (std::size_t i=0; i<training_points_.size(); ++i) {
      auto const& y = training_points_[i];
      res += weights_[i] * covarianceFunc_->evaluate(x, y);
    }
  }

  if (has_derivatives_) {
    auto start = has_values_ ? training_points_.size() : 0;
    for (std::size_t i=0; i<training_points_.size(); ++i) {
      auto const& y = training_points_[i];
      for (std::size_t d=0; d<ndim; ++d)
        res += weights_[start + i*ndim + d] * covarianceFunc_->first_der_y(x, y, d);
    }
  }
  return res;
}

double gpr::GPR_Interpolator::interpolate_derivative(PES_Point const& x, std::size_t d) const {
  double res = 0;
  auto ndim = training_points_.front().size();

  if (has_values_) {
    for (std::size_t i=0; i<training_points_.size(); ++i) {
      auto const& y = training_points_[i];
      res += weights_[i] * covarianceFunc_->first_der_x(x, y, d);
    }
  }

  if (has_derivatives_) {
    auto start = has_values_ ? training_points_.size() : 0;
    for (std::size_t i=0; i<training_points_.size();  ++i) {
      auto const& y = training_points_[i];
      for (std::size_t d_j=0; d_j<ndim; ++d_j)
        res += weights_[start + i*ndim + d] * covarianceFunc_->second_der(x, y, d, d_j);
    }
  }
  return res;
}

std::vector<double> const& gpr::GPR_Interpolator::get_weights() const {
  return weights_;
}

void gpr::run_gpr_test() {
  auto interpolator = pmf_ic::load_interpolation();

  if (auto interpolator_2d = std::get_if<std::unique_ptr<pmf_ic::Interpolator2DInterface>>(&interpolator)) {
    std::ofstream out("output_PMF_IC_interpolation.txt");
    out.precision(10);
    auto xi1_range = Config::get().coords.umbrella.pmf_ic.ranges[0];
    auto xi2_range = Config::get().coords.umbrella.pmf_ic.ranges[1];
    for (double xi1 = xi1_range.start; xi1 <= xi1_range.stop; xi1 += xi1_range.step) {
      out << ',' << xi1;
    }
    out << ",xi1";
    for (double xi2 = xi2_range.start; xi2 <= xi2_range.stop; xi2 += xi2_range.step) {
      std::cout << xi2 << std::endl;
      out << '\n' << xi2;
      for (double xi1 = xi1_range.start; xi1 <= xi1_range.stop; xi1 += xi1_range.step) {
        out << ',' << (*interpolator_2d)->get_value(xi1, xi2);
      }
    }
    out << "\nxi2";
  }
}
