#ifndef cast_ic_rotation_h_guard
#define cast_ic_rotation_h_guard

#pragma once

#include "ic_util.h"
#include "quaternion.h"
#include "scon_angle.h"
#include "scon_spherical.h"
#include "scon_vect.h"

#include <algorithm>
#include <armadillo>
#include "scon_mathmatrix.h"
#include <array>
#include <utility>
#include <vector>

namespace ic_rotation {

using coords::float_type;

static constexpr auto q_thres{ 1e-6 };

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
arma::Mat<T> correlation_matrix(const coords::Representation_3D& trial,
                                const coords::Representation_3D& target) {
  auto trial_mat = ic_util::Rep3D_to_arma<T>(trial);
  auto target_mat = ic_util::Rep3D_to_arma<T>(target);
  return trial_mat.t() * target_mat;
}

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
arma::Mat<T> F_matrix(const coords::Representation_3D& trial,
                      const coords::Representation_3D& target) {
  using Mat = arma::Mat<T>;

  auto correl_mat = correlation_matrix<T>(trial, target);
  auto R11 = correl_mat(0, 0);
  auto R12 = correl_mat(0, 1);
  auto R13 = correl_mat(0, 2);
  auto R21 = correl_mat(1, 0);
  auto R22 = correl_mat(1, 1);
  auto R23 = correl_mat(1, 2);
  auto R31 = correl_mat(2, 0);
  auto R32 = correl_mat(2, 1);
  auto R33 = correl_mat(2, 2);
  Mat F(4, 4);
  F(0, 0) = R11 + R22 + R33;
  F(0, 1) = R23 - R32;
  F(0, 2) = R31 - R13;
  F(0, 3) = R12 - R21;
  F(1, 0) = R23 - R32;
  F(1, 1) = R11 - R22 - R33;
  F(1, 2) = R12 + R21;
  F(1, 3) = R13 + R31;
  F(2, 0) = R31 - R13;
  F(2, 1) = R12 + R21;
  F(2, 2) = -R11 + R22 - R33;
  F(2, 3) = R23 + R32;
  F(3, 0) = R12 - R21;
  F(3, 1) = R13 + R31;
  F(3, 2) = R23 + R32;
  F(3, 3) = -R11 - R22 + R33;
  return F;
}

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
std::pair<T, ic_util::Quaternion<T>>
quaternion(const coords::Representation_3D& trial,
           const coords::Representation_3D& target) {
  using Col = arma::Col<T>;
  using Mat = arma::Mat<T>;
  using arma::conv_to;
  using arma::eig_sym;
  using arma::size;
  using arma::sort;
  using coords::Cartesian_Point;

  Cartesian_Point mean_trial(0.0, 0.0, 0.0);
  for (auto& i : trial) {
    mean_trial += i;
  }
  mean_trial /= (T)trial.size();
  Cartesian_Point mean_target(0.0, 0.0, 0.0);
  for (auto& i : target) {
    mean_target += i;
  }
  mean_target /= (T)target.size();
  auto trial_shift = trial - mean_trial;
  auto target_shift = target - mean_target;
  auto F_mat = F_matrix<T>(trial_shift, target_shift);
  Col eigval;
  Mat eigvec;
  eig_sym(eigval, eigvec, F_mat);
  auto eigval_sz = eigval.n_rows;
  Col q = eigvec.col(eigval_sz - 1);
  std::vector<T> q_std = conv_to<std::vector<T>>::from(q);
  ic_util::Quaternion<T> res_q(q_std);
  if (res_q.q_.at(0) < 0) {
    res_q = res_q * -1;
  }
  Col eigval_high = sort(eigval, "descend");
  T res_val = eigval_high.at(0);
  return std::make_pair(res_val, res_q);
}

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
std::array<T, 3u> exponential_map(const coords::Representation_3D& trial,
                                  const coords::Representation_3D& target) {
  auto q = quaternion<T>(trial, target);
  auto quat = q.second;
  auto q0 = quat.q_.at(0);
  auto t1 = std::acos(q0);
  auto t2 = std::sqrt(1 - std::pow(q0, 2));
  auto p = 2 * (t1 / t2);
  std::array<T, 3u> result = { p * quat.q_.at(1), p * quat.q_.at(2),
                               p * quat.q_.at(3) };
  return result;
}

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
std::vector<arma::Cube<T>>
correlation_matrix_derivs(const arma::Mat<T>& R,
                          const coords::Representation_3D& target) {
  using Cube = arma::Cube<T>;
  using Mat = arma::Mat<T>;
  using arma::zeros;
  using coords::Cartesian_Point;

  Cartesian_Point mean_target(0.0, 0.0, 0.0);
  for (auto& i : target) {
    mean_target += i;
  }
  mean_target /= (T)target.size();
  auto tshift = target - mean_target;
  auto S = ic_util::Rep3D_to_arma<T>(tshift);
  std::vector<Cube> result;
  for (std::size_t c = 0; c < S.n_rows; ++c) {
    Cube A = zeros<Cube>(3, 3, 3);
    for (std::size_t l = 0; l < 3; ++l) {
      for (std::size_t r = 0; r < 3; ++r) {
        if (l == r) {
          A.slice(l).row(r) = S.row(c);
        }
      }
    }
    result.emplace_back(A);
  }
  return result;
}

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
std::vector<arma::Cube<T>>
F_matrix_derivs(const coords::Representation_3D& trial,
                const coords::Representation_3D& target) {
  using Cube = arma::Cube<T>;
  using Mat = arma::Mat<T>;
  using arma::zeros;
  using coords::Cartesian_Point;

  Cartesian_Point mean_target(0.0, 0.0, 0.0);
  for (auto& i : target) {
    mean_target += i;
  }
  mean_target /= (T)target.size();
  Cartesian_Point mean_trial(0.0, 0.0, 0.0);
  for (auto& i : trial) {
    mean_trial += i;
  }
  mean_trial /= (T)trial.size();
  auto tashift = target - mean_target;
  auto trshift = trial - mean_trial;
  auto R = correlation_matrix<T>(tashift, trshift);
  auto dR = correlation_matrix_derivs(R, tashift);
  std::vector<Cube> result;
  for (auto& S : dR) {
    Cube Q = zeros<Cube>(4, 4, 3);
    for (std::size_t c = 0; c < S.n_slices; ++c) {
      Mat F = zeros<Mat>(4, 4);
      Mat M = S.slice(c);
      auto dR11 = M(0, 0);
      auto dR12 = M(0, 1);
      auto dR13 = M(0, 2);
      auto dR21 = M(1, 0);
      auto dR22 = M(1, 1);
      auto dR23 = M(1, 2);
      auto dR31 = M(2, 0);
      auto dR32 = M(2, 1);
      auto dR33 = M(2, 2);
      F(0, 0) = dR11 + dR22 + dR33;
      F(0, 1) = dR23 - dR32;
      F(0, 2) = dR31 - dR13;
      F(0, 3) = dR12 - dR21;
      F(1, 0) = dR23 - dR32;
      F(1, 1) = dR11 - dR22 - dR33;
      F(1, 2) = dR12 + dR21;
      F(1, 3) = dR13 + dR31;
      F(2, 0) = dR31 - dR13;
      F(2, 1) = dR12 + dR21;
      F(2, 2) = -dR11 + dR22 - dR33;
      F(2, 3) = dR23 + dR32;
      F(3, 0) = dR12 - dR21;
      F(3, 1) = dR13 + dR31;
      F(3, 2) = dR23 + dR32;
      F(3, 3) = -dR11 - dR22 + dR33;
      Q.slice(c) = F;
    }
    result.emplace_back(Q);
  }
  return result;
}

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
std::vector<arma::Mat<T>>
quaternion_derivs(const coords::Representation_3D& trial,
                  const coords::Representation_3D& target) {
  using Mat = arma::Mat<T>;
  using Row = arma::Row<T>;
  using arma::conv_to;
  using arma::eye;
  using arma::pinv;

  auto q_eigval = quaternion<T>(trial, target);
  auto F = F_matrix<T>(trial, target);
  auto F_der = F_matrix_derivs<T>(trial, target);
  auto identity = eye<Mat>(4, 4);
  auto lambda = q_eigval.first;
  auto t1 = lambda * identity - F;
  auto t2 = pinv(t1);
  std::vector<Mat> result;
  auto q = q_eigval.second;
  auto q_in_vec = ic_util::arr_to_vec(q.q_);
  Row qrow = conv_to<Row>::from(q_in_vec);
  for (std::size_t i = 0; i < target.size(); ++i) {
    auto F_dtemp = F_der.at(i);
    Mat qtemp(3, 4);
    for (std::size_t c = 0; c < 3; ++c) {
      auto F_sl = F_dtemp.slice(c);
      auto temp_res = t2 * F_sl * qrow.t();
      qtemp.row(c) = temp_res.t();
    }
    result.emplace_back(qtemp);
  }
  return result;
}

template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
std::vector<arma::Mat<T>>
exponential_derivs(const coords::Representation_3D& trial,
                   const coords::Representation_3D& target) {
  using Mat = arma::Mat<T>;

  auto q_val = quaternion<T>(trial, target);
  auto q = q_val.second;
  auto q0 = std::get<0>(q.q_);
  T d{ 0.0 };
  if (std::abs(q0 - 1) < q_thres) {
    d = -(T)2 / 3;
  } else {
    auto t1 = std::acos(q0);
    auto t2 = 1 - std::pow(q0, 2);
    auto t3 = std::pow(t2, 3);
    auto term1 = (2 * q0 * t1) / std::sqrt(t3);
    auto term2 = 2 / t2;
    d = term1 - term2;
  }
  auto temp1 = std::acos(q0);
  auto temp2 = std::sqrt(1 - std::pow(q0, 2));
  auto p = 2 * (temp1 / temp2);
  Mat dv_mat = { { d * q.q_.at(1), d * q.q_.at(2), d * q.q_.at(3) },
                 { p, 0, 0 },
                 { 0, p, 0 },
                 { 0, 0, p } };
  auto dv_trans = dv_mat.t();
  auto dq_mat = quaternion_derivs<T>(trial, target);
  std::vector<Mat> result;
  for (std::size_t i = 0; i < target.size(); ++i) {
    Mat M = dq_mat.at(i);
    Mat R = dv_trans * M.t();
    result.emplace_back(R);
  }
  return result;
}
}

#endif // cast_ic_rotation_h_guard