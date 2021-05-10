#include "gpr.h"

#include <cassert>


gpr::GPR_Interpolator::GPR_Interpolator(gpr::KernelFunction kf,
                                        std::vector<PES_Point> training_points,
                                        const std::vector<double> &training_data)
  :kernel_{std::move(kf)}
  ,training_points_{std::move(training_points)}
{
  assert(training_points_.size() == training_data.size());
  train_gp(training_data);
}

gpr::GPR_Interpolator gpr::gpr_interpolator_1d(KernelFunction kf, const std::vector<double> &training_points,
                                          const std::vector<double> &training_data) {
  std::vector<PES_Point> pes_points;
  pes_points.reserve(training_data.size());
  std::transform(training_points.begin(), training_points.end(), std::back_inserter(pes_points),
                 [](auto x) -> PES_Point {return {x};});

  return GPR_Interpolator(std::move(kf), std::move(pes_points), training_data);
}

void gpr::GPR_Interpolator::train_gp(const std::vector<double> &training_data) {
  scon::mathmatrix<double> K(training_points_.size(), training_points_.size(), 0);
  for (std::size_t i=0; i<training_points_.size(); ++i) {
    for (std::size_t j=0; j<training_points_.size(); ++j) {
      K(i, j) = kernel_(training_points_[i], training_points_[j]);
    }
  }
  //std::cout << K << '\n';

  // Why is std::accumulate weired?
  for(auto y_i: training_data)
    y_prior_ += y_i;
  y_prior_ /= training_data.size();

  assert(K.positive_definite_check());
  auto y = scon::mathmatrix<double>::col_from_vec(training_data);
  y = y - y_prior_;
  auto w_mat = K.solve(y);
  weights_ = scon::mathmatrix<double>(w_mat).col_to_std_vector();
  assert(weights_.size() == training_points_.size());
}

double gpr::GPR_Interpolator::interpolate(const gpr::PES_Point &x) const {
  double res = y_prior_;
  for (std::size_t i=0; i<training_points_.size(); ++i)
    res += weights_[i] * kernel_(x, training_points_[i]);
  return res;
}

void gpr::run_gpr_test() {
  std::ifstream in("gpr-input.txt");
  std::vector<double> x, y;
  while(!in.eof()) {
    double new_x, new_y;
    in >> new_x >> new_y;
    x.emplace_back(new_x);
    y.emplace_back(new_y);
  }

  // For some reason, the last values are duplicated
  x.erase(x.end()-1);
  y.erase(y.end()-1);

  auto gpr = gpr_interpolator_1d(exponential_kernel(2), x, y);
  std::ofstream out("gpr-output.txt");
  for(double i=x.front(); i<=x.back(); i+=0.02) {
    out << i << ' ' << gpr.interpolate({i}) << '\n';
  }
}
