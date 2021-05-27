#include "gpr.h"

#include <cassert>

gpr::GPR_Interpolator::GPR_Interpolator(std::unique_ptr <gpr::KernelFunction> kf,
                                        std::vector<PES_Point> training_points,
                                        const std::vector<double> &training_data)
  :kernel_{std::move(kf)}
  ,training_points_{std::move(training_points)}
{
  assert(training_points_.size() == training_data.size());
  train_gp(training_data);
}

gpr::GPR_Interpolator gpr::gpr_interpolator_1d(std::unique_ptr<KernelFunction> kf, const std::vector<double> &training_points,
                                          const std::vector<double> &training_data) {
  std::vector<PES_Point> pes_points;
  pes_points.reserve(training_data.size());
  std::transform(training_points.begin(), training_points.end(), std::back_inserter(pes_points),
                 [](auto x) -> PES_Point {return {x};});

  return GPR_Interpolator(std::move(kf), std::move(pes_points), training_data);
}

gpr::GPR_Interpolator gpr::gpr_interpolator_2d(std::unique_ptr<KernelFunction> kf, const std::vector<std::pair<double, double>> &training_points,
                                               const std::vector<double> &training_data) {
  std::vector<PES_Point> pes_points;
  pes_points.reserve(training_data.size());
  std::transform(training_points.begin(), training_points.end(), std::back_inserter(pes_points),
                 [](auto point) -> PES_Point {return {point.first, point.second};});

  return GPR_Interpolator(std::move(kf), std::move(pes_points), training_data);
}

void gpr::GPR_Interpolator::train_gp(const std::vector<double> &training_data) {
  scon::mathmatrix<double> K(training_points_.size(), training_points_.size(), 0);
  for (std::size_t i=0; i<training_points_.size(); ++i) {
    for (std::size_t j=0; j<training_points_.size(); ++j) {
      K(i, j) = kernel_->evaluate(training_points_[i], training_points_[j]);
    }
  }
  //K += scon::mathmatrix<double>::identity(K.rows(), K.cols()) * 0.1;
  //std::cout << K << '\n';

  // Why is std::accumulate weired?
  for(auto y_i: training_data)
    y_prior_ += y_i;
  y_prior_ /= training_data.size();

  //assert(K.positive_definite_check());
  auto y = scon::mathmatrix<double>::col_from_vec(training_data);
  y = y - y_prior_;
  auto w_mat = K.solve(y);
  weights_ = scon::mathmatrix<double>(w_mat).col_to_std_vector();
  assert(weights_.size() == training_points_.size());
}

double gpr::GPR_Interpolator::interpolate(const gpr::PES_Point &x) const {
  double res = y_prior_;
  for (std::size_t i=0; i<training_points_.size(); ++i)
    res += weights_[i] * kernel_->evaluate(x, training_points_[i]);
  return res;
}

double gpr::GPR_Interpolator::interpolate_derivative(PES_Point const& x, std::size_t d) const {
  double res = 0;
  for (std::size_t i=0; i<training_points_.size(); ++i)
    res += weights_[i] * kernel_->first_der_x(x, training_points_[i], d);
  return res;
}

std::vector<double> const& gpr::GPR_Interpolator::get_weights() const {
  return weights_;
}

void gpr::run_gpr_test() {
  std::ifstream in("gpr-input.txt");
  //std::vector<PES_Point> x;
  std::vector<double> x, y, y_der;
  while(!in.eof()) {
    /*std::string ign;
    double x1, x2, new_y;
    in >> x1 >> ign >> x2 >> ign >> ign >> ign >> new_y;
    x.emplace_back(PES_Point{x1, x2});*/
    double new_x, new_y, new_y_der;
    in >> new_x >> new_y >> new_y_der;
    x.emplace_back(new_x);
    y.emplace_back(new_y);
    y_der.emplace_back(new_y_der);
  }

  // For some reason, the last values are duplicated
  x.erase(x.end()-1);
  y.erase(y.end()-1);
  y_der.erase(y_der.end()-1);

  auto gpr = gpr_interpolator_1d(std::make_unique<SqExpKernel>(0.5), x, y);
  std::ofstream out("gpr-output.txt");
  for (double x=-10; x<=10; x+= 0.1)
    out << x << ' ' << gpr.interpolate({x}) << ' ' << gpr.interpolate_derivative({x}, 0) << '\n';
  /*for (double xi=-180; xi<=180; xi+=1) {
    out << ',' << xi;
  }
  out << ",xi1";

  for (double xi2=-180; xi2<=180; xi2+=1) {
    out << '\n' << xi2;
    for (double xi1=-180; xi1<=180; xi1+=1) {
      out << ',' << gpr.interpolate({xi1, xi2});
    }
  }
  out << "\nxi2";*/

  auto weights = gpr.get_weights();
  std::ofstream weights_file("weights.txt");
  for (auto w: weights)
    weights_file << w << '\n';
}
