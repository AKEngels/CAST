#include "gpr.h"

#include <cassert>

#include "configuration.h"

#include "pmf_interpolator_builder.h"

gpr::GPR_Interpolator::GPR_Interpolator(std::unique_ptr <gpr::KernelFunction> kf,
                                        std::vector<PES_Point> training_points,
                                        std::vector<double> const& training_data,
                                        std::optional<std::vector<PES_Point>> const& training_gradients)
  :kernel_{std::move(kf)}
  ,training_points_{std::move(training_points)}
  ,has_derivatives_(training_gradients.has_value())
{
  assert(training_points_.size() == training_data.size());
  train_gp(training_data, training_gradients);
}

gpr::GPR_Interpolator gpr::gpr_interpolator_1d(std::unique_ptr<KernelFunction> kf,
                                               const std::vector<double> &training_points,
                                               const std::vector<double> &training_data,
                                               std::optional<std::vector<double>> const& training_gradients) {
  std::vector<PES_Point> pes_points;
  pes_points.reserve(training_data.size());
  auto transformer = [](auto x) -> PES_Point {return {x};};
  std::transform(training_points.begin(), training_points.end(), std::back_inserter(pes_points),
                 transformer);

  std::optional<std::vector<PES_Point>> gradients;
  if (training_gradients) {
    gradients = std::vector<PES_Point>();
    gradients->reserve(training_gradients->size());
    std::transform(training_gradients->begin(), training_gradients->end(), std::back_inserter(*gradients),
                   transformer);
  }

  return GPR_Interpolator(std::move(kf), std::move(pes_points), training_data, gradients);
}

gpr::GPR_Interpolator gpr::gpr_interpolator_2d(std::unique_ptr<KernelFunction> kf, std::vector<std::pair<double, double>> const& training_points,
                                               std::vector<double> const& training_data,
                                               std::optional<std::vector<std::pair<double, double>>> const& training_gradients) {
  std::vector<PES_Point> pes_points;
  std::optional<std::vector<PES_Point>> pes_gradients;
  pes_points.reserve(training_data.size());
  std::transform(training_points.begin(), training_points.end(), std::back_inserter(pes_points),
                 [](auto point) -> PES_Point {return {point.first, point.second};});

  if (training_gradients) {
    pes_gradients = std::vector<PES_Point>();
    pes_gradients->reserve(training_gradients->size());
    std::transform(training_gradients->begin(), training_gradients->end(), std::back_inserter(*pes_gradients),
                   [](auto point) -> PES_Point {return {point.first, point.second};});
  }

  return GPR_Interpolator(std::move(kf), std::move(pes_points), training_data, std::move(pes_gradients));
}

void gpr::GPR_Interpolator::train_gp(const std::vector<double> &training_data,
                                     std::optional<std::vector<PES_Point>> const& training_gradients) {
  auto n = training_data.size();
  auto ndim = training_points_.front().size();

  auto cov_size = has_derivatives_ ? n * (ndim+1) : training_points_.size();

  scon::mathmatrix<double> K(cov_size, cov_size, 0);
  for (std::size_t i=0; i<n; ++i) {
    for (std::size_t j=0; j<n; ++j) {
      auto const& x = training_points_[i];
      auto const& y = training_points_[j];

      // Upper left
      K(i, j) = kernel_->evaluate(x, y);

      if (has_derivatives_) {
        for (std::size_t d_i = 0; d_i < ndim; ++d_i) {
          // Upper right
          K(i, n + j*ndim + d_i) = kernel_->first_der_y(x, y, d_i);

          // Lower left
          K(n + i*ndim + d_i, j) = kernel_->first_der_x(x, y, d_i);

          // Lower right
          for (std::size_t d_j=0; d_j<ndim; ++d_j) {
            K(n + i*ndim + d_i, n + j*ndim + d_j) = kernel_->second_der(x, y, d_i, d_j);
          }
        }
      }
    }
  }
  //K += scon::mathmatrix<double>::identity(K.rows(), K.cols()) * 0.1;
  if (Config::get().general.verbosity >= 4)
    std::cout << K << '\n';

  // Why is std::accumulate weired?
  for(auto y_i: training_data)
    y_prior_ += y_i;
  y_prior_ /= training_data.size();

  //assert(K.positive_definite_check());
  auto y = [&] {
    if (has_derivatives_) {
      auto res = scon::mathmatrix<double>::zero(n*(ndim+1), 1);
      for (std::size_t i=0; i<n; ++i) {
        res(i, 0) = training_data[i] - y_prior_;

        for (std::size_t d=0; d<ndim; ++d) {
          res(n + i*ndim + d, 0) = (*training_gradients)[i][d];
        }
      }
      return res;
    }
    else
      return scon::mathmatrix<double>::col_from_vec(training_data) - y_prior_;
  }();
  auto w_mat = K.solve(y);
  weights_ = scon::mathmatrix<double>(w_mat).col_to_std_vector();

  assert(weights_.size() == cov_size);
}

double gpr::GPR_Interpolator::interpolate(const gpr::PES_Point &x) const {
  double res = y_prior_;
  auto ndim = training_points_.front().size();
  for (std::size_t i=0; i<training_points_.size(); ++i) {
    auto const& y = training_points_[i];
    res += weights_[i] * kernel_->evaluate(x, y);

    if (has_derivatives_) {
      for (std::size_t d=0; d<ndim; ++d)
        res += weights_[training_points_.size() + i*ndim + d] * kernel_->first_der_y(x, y, d);
    }
  }
  return res;
}

double gpr::GPR_Interpolator::interpolate_derivative(PES_Point const& x, std::size_t d) const {
  double res = 0;
  auto ndim = training_points_.front().size();
  for (std::size_t i=0; i<training_points_.size(); ++i) {
    auto const& y = training_points_[i];
    res += weights_[i] * kernel_->first_der_x(x, y, d);;

    if (has_derivatives_) {
      for (std::size_t d_j=0; d_j<ndim; ++d_j)
        res += weights_[training_points_.size() + i*ndim + d] * kernel_->second_der(x, y, d, d_j);
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
