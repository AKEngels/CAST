#include "gpr.h"

#include <cassert>

#include "configuration.h"

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

gpr::GPR_Interpolator gpr::gpr_interpolator_2d(std::unique_ptr<KernelFunction> kf, const std::vector<std::pair<double, double>> &training_points,
                                               const std::vector<double> &training_data) {
  std::vector<PES_Point> pes_points;
  pes_points.reserve(training_data.size());
  std::transform(training_points.begin(), training_points.end(), std::back_inserter(pes_points),
                 [](auto point) -> PES_Point {return {point.first, point.second};});

  return GPR_Interpolator(std::move(kf), std::move(pes_points), training_data);
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
  std::ifstream in("gpr-input.txt");
  std::vector<PES_Point> x, y_der;
  std::vector<double> y;
  while(!in.eof()) {
    double x1, x2, new_y, y_der1, y_der2;
    in >> x1 >> x2 >> new_y >> y_der1 >> y_der2;
    x.emplace_back(PES_Point{x1, x2});
    y.emplace_back(new_y);
    y_der.emplace_back(PES_Point {y_der1, y_der2});
  }

  // For some reason, the last values are duplicated
  x.erase(x.end()-1);
  y.erase(y.end()-1);
  y_der.erase(y_der.end()-1);

  auto gpr = GPR_Interpolator(std::make_unique<MaternKernel>(10), x, y, y_der);
  std::ofstream out("gpr-output.txt");
  for (double xi=-10; xi<=10; xi+=0.1) {
    out << ',' << xi;
  }
  out << ",xi1";

  for (double xi2=-10; xi2<=10; xi2+=0.1) {
    out << '\n' << xi2;
    for (double xi1=-10; xi1<=10; xi1+=0.1) {
      out << ',' << gpr.interpolate({xi1, xi2});
    }
  }
  out << "\nxi2";

  auto weights = gpr.get_weights();
  std::ofstream weights_file("weights.txt");
  for (auto w: weights)
    weights_file << w << '\n';
}
