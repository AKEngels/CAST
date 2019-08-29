#ifndef cast_ic_rotation_h_guard
#define cast_ic_rotation_h_guard

#include "InternalCoordinateUtilities.h"
#include "quaternion.h"

#include <algorithm>
#include <array>
#include <utility>
#include <vector>
#include <functional>
#include <numeric>



namespace ic_rotation {

using coords::float_type;

auto constexpr q_thres{ 1e-6 };


template <template<typename> class MatrixType, typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
typename std::enable_if<std::is_arithmetic<T>::value, MatrixType<T>>::type
correlation_matrix(ContainerType<CoordType<T>, ContainerArgs...> const& old_xyz,
                   ContainerType<CoordType<T>, ContainerArgs...> const& new_xyz) {

  auto new_xyz_mat = ic_util::Rep3D_to_Mat<MatrixType>(new_xyz - ic_util::get_mean(new_xyz));
  auto old_xyz_mat = ic_util::Rep3D_to_Mat<MatrixType>(old_xyz - ic_util::get_mean(old_xyz));
  return new_xyz_mat.t() * old_xyz_mat;
}

template <template<typename> class MatrixType, typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
typename std::enable_if<std::is_arithmetic<T>::value, MatrixType<T>>::type
F_matrix(ContainerType<CoordType<T>, ContainerArgs...> const& old_xyz,
         ContainerType<CoordType<T>, ContainerArgs...> const& new_xyz) {

  using Mat = MatrixType<T>;

  auto correl_mat = correlation_matrix<MatrixType>(old_xyz, new_xyz);
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

template <template<typename> class MatrixType, typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
typename std::enable_if<std::is_arithmetic<T>::value, std::pair<T, ic_util::Quaternion<T>>>::type
quaternion(ContainerType<CoordType<T>, ContainerArgs...> const& old_xyz,
           ContainerType<CoordType<T>, ContainerArgs...> const& new_xyz) {
  using Mat = MatrixType<T>;

  auto F_mat = F_matrix<MatrixType>(old_xyz, new_xyz);

  Mat eigvec, eigval;

  std::tie(eigval, eigvec) = F_mat.eigensym();
  //shouldn't that be col_to_vector(0)? See rotate.py line 272 (in get_quat) //Seems to be right like described here. They are sorting it the other way round
  auto q_std = eigvec.col_to_std_vector(eigvec.cols()-1);
  ic_util::Quaternion<T> res_q(q_std);
  if (res_q.q_.at(0) < 0.) {
    res_q = res_q * -1.;
  }

  return std::make_pair(eigval(eigval.rows()-1,0), res_q);
}

template<template<typename> class MatrixType, typename T>
MatrixType<T> form_rot(ic_util::Quaternion<T> const& q){
  auto al = ic_util::al(q);
  auto ar = ic_util::ar(ic_util::conj(q));
  auto ret = al*ar;
  ret.shed_cols(0);
  ret.shed_rows(0);
  return ret;
}

template <template<typename> class MatrixType, typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
typename std::enable_if<std::is_arithmetic<T>::value, std::array<T, 3u>>::type
exponential_map(ContainerType<CoordType<T>, ContainerArgs...> const& old_xyz,
                ContainerType<CoordType<T>, ContainerArgs...> const& new_xyz) {

  auto get_fac = [](auto const & q0) {
    auto qm1 = q0 - 1.;
    if (std::fabs(qm1) < 1.e-8) {
      return 2. - 2.*qm1/3.;
    }
    else {
      auto t1 = std::acos(q0);
      auto t2 = std::sqrt(1 - std::pow(q0, 2));
      return 2. * (t1 / t2);
    }
  };

  auto q = quaternion<MatrixType>(old_xyz, new_xyz);
  auto quat = q.second;
  auto p = get_fac(quat.q_.at(0));
  return { p * quat.q_.at(1), p * quat.q_.at(2), p * quat.q_.at(3) };
}

template <template<typename> class MatrixType, typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
  typename std::enable_if<std::is_arithmetic<T>::value, std::vector<std::vector<MatrixType<T>> >>::type
correlation_matrix_derivs(ContainerType<CoordType<T>, ContainerArgs...> const& new_xyz) {
  using Mat = MatrixType<T>;

  //auto const & add_cp = std::plus<coords::Cartesian_Point>();
  auto S = ic_util::Rep3D_to_Mat<MatrixType>(new_xyz);
  std::vector<std::vector<Mat> > result;
  for (auto c = 0u; c < S.rows(); ++c) {
    std::vector<Mat> A(3, Mat::zero(3, 3));
    for (std::size_t l = 0; l < A.size(); ++l) {
      A.at(l).set_row(l, S.row(c));
    }
    result.emplace_back(A);
  }
  return result;
}

template <template<typename> class MatrixType, typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
typename std::enable_if<std::is_arithmetic<T>::value, std::vector<std::vector<MatrixType<T>> >>::type
F_matrix_derivs(ContainerType<CoordType<T>, ContainerArgs...> const& new_xyz) {
  using Mat = MatrixType<T>;

  auto new_shift = new_xyz - ic_util::get_mean(new_xyz);
  auto dR = correlation_matrix_derivs<MatrixType>(new_shift);
  std::vector<std::vector<Mat> > result;
  for (auto&& S : dR) {
    std::vector<Mat> Q;
    for (auto&& M : S) {
      auto F = MatrixType<float_type>::zero(4, 4);
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
      Q.emplace_back(std::move(F));
    }
    result.emplace_back(Q);
  }
  return result;
}

template <template<typename> class MatrixType, typename T, template <typename> class CoordType,
          template <typename, typename...> class ContainerType,
          typename... ContainerArgs>
typename std::enable_if<std::is_arithmetic<T>::value,
                        std::vector<MatrixType<T>>>::type
quaternion_derivs(
    ContainerType<CoordType<T>, ContainerArgs...> const& old_xyz,
    ContainerType<CoordType<T>, ContainerArgs...> const& new_xyz) {
  using Mat = MatrixType<T>;

  auto q_eigval = quaternion<scon::mathmatrix>(old_xyz, new_xyz);
  auto F = F_matrix<scon::mathmatrix>(old_xyz, new_xyz);
  auto F_der = F_matrix_derivs<scon::mathmatrix>(old_xyz); // is it really the old xyz???

  auto t = (Mat::fill_diag(4, 4, q_eigval.first) - F).pinv();

  std::vector<Mat> result;
  auto q = q_eigval.second;
  auto q_in_vec = ic_util::arr_to_vec(q.q_);
  auto qrow = Mat::col_from_vec(q_in_vec);
  for (auto&& Fd : F_der) {
    Mat qtemp(3, 4);
    for (std::size_t c = 0; c < Fd.size(); ++c) {
      qtemp.set_row(c, (t * Fd.at(c) * qrow).t());
    }
    result.emplace_back(qtemp);
  }
  return result;
}

template <template<typename> class MatrixType, typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
typename std::enable_if<std::is_arithmetic<T>::value, std::vector<MatrixType<T>>>::type
exponential_derivs(ContainerType<CoordType<T>, ContainerArgs...> const& old_xyz,
                   ContainerType<CoordType<T>, ContainerArgs...> const& new_xyz) {
  using Mat = MatrixType<T>;

  auto const fac_and_dfac = [](auto const & q0) {
    if (std::abs(q0 - 1.) < q_thres) {
      return std::make_pair(2. - 2. * (q0 - 1.) / 3., -2. / 3.);
    }
    auto acosq0 = std::acos(q0);
    auto q0_sq = 1. - q0 * q0;
    return std::make_pair(2.*acosq0 / std::sqrt(q0_sq), 2.*q0*acosq0 / std::pow(q0_sq, 1.5) - 2./q0_sq /*<-So at least in their paper (Lee Ping)*/);
  };

  auto q_val = quaternion<MatrixType>(old_xyz, new_xyz);
  auto q = q_val.second;
  auto q0 = std::get<0>(q.q_);
  //std::cout << q << "    " << q0 << "\n\n";
  T p{ 0.0 }, d{ 0.0 };
  std::tie(p, d) = fac_and_dfac(q0);
  //std::cout << p << " " << d << "\n\n";
  auto dv_mat = Mat({ { d * q.q_.at(1), d * q.q_.at(2), d * q.q_.at(3) },
                 { p, 0, 0 },
                 { 0, p, 0 },
                 { 0, 0, p } });
  //std::cout << "dv_mat:\n" << dv_mat << "\n\n";
  auto dq_mat = quaternion_derivs<MatrixType>(old_xyz, new_xyz);
  std::vector<Mat> result;
  for (auto i{ 0u }; i < old_xyz.size(); ++i) {
    Mat R(3,3);
    auto const& dqdx_i = dq_mat.at(i);
    for(auto j{ 0u }; j < dqdx_i.rows(); ++j){
      auto const& dqdx_ij = dqdx_i.row(j);
      auto new_row = Mat::zero(1,3);
      for(auto k{ 0u }; k<dqdx_i.cols(); ++k){
        new_row += dv_mat.row(k) * dqdx_ij(0, k);
      }
      R.set_row(j,new_row);
    }
    result.emplace_back(R);
  }
  return result;
}

template<template<typename> class MatrixType, typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
typename std::enable_if<std::is_arithmetic<T>::value, std::pair<T, T>>::type
 displacementRmsValAndMaxTwoStructures(ContainerType<CoordType<T>, ContainerArgs...> const& lhs, ContainerType<CoordType<T>, ContainerArgs...> const& rhs) {
	auto q = quaternion<MatrixType>(lhs, rhs);
	auto U = form_rot<MatrixType>(q.second);

	auto new_xyz_mat = ic_util::Rep3D_to_Mat<MatrixType>(rhs - ic_util::get_mean(rhs));
	auto old_xyz_mat = ic_util::Rep3D_to_Mat<MatrixType>(lhs - ic_util::get_mean(lhs));
	auto rot = (U*new_xyz_mat.t()).t();

	old_xyz_mat -= rot;
	old_xyz_mat *= -energy::bohr2ang;
	auto norms = ic_util::atomsNorm(old_xyz_mat);

	return { norms.rmsd(), norms.max() };
}


}
/*template <typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
typename std::enable_if<std::is_arithmetic<T>::value, std::vector<scon::mathmatrix<T>>>::type
ic_rotation::quaternion_derivs(ContainerType<CoordType<T>, ContainerArgs...> const& old_xyz,
                               ContainerType<CoordType<T>, ContainerArgs...> const& new_xyz) {
  using Mat = scon::mathmatrix<T>;

  auto q_eigval = quaternion(old_xyz, new_xyz);
  auto F = F_matrix(old_xyz, new_xyz);
  auto F_der = F_matrix_derivs(new_xyz);

  auto t = (Mat::fill_diag(4, 4, q_eigval.first) -F).pinv();

  std::vector<Mat> result;
  auto q = q_eigval.second;
  auto q_in_vec = ic_util::arr_to_vec(q.q_);
  auto qrow = Mat::col_from_vec(q_in_vec);
  for (auto&& Fd : F_der) {
    Mat qtemp(3, 4);
    for (std::size_t c = 0; c < Fd.size(); ++c) {
      qtemp.set_row(c, (t * Fd.at(c) * qrow).t());
    }
    result.emplace_back(qtemp);
  }
  return result;
}*/
#endif // cast_ic_rotation_h_guard
