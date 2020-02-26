#include"BestFitRotation.h"

#include<vector>

namespace{
Eigen::MatrixXd matrixToEigenMatrixXd(AbstractMatrix const& matrix){
    Eigen::MatrixXd matrixXd(matrix.rows(), matrix.cols());
    for(std::size_t i = 0u; i < matrix->rows(); ++i){
        for(std::size_t j = 0u; j < matrix->cols(); ++j){
            matrixXd(i,j) = matrix(i, j);
        }
    }
    return matrixXd;
}

Eigen::Vector3d addAllRows(Eigen::MatrixXd const& matrix){
    auto sum = Eigen::Vector3d::Zero();

    for(std::size_t row = 0u; row < matrix.rows(); ++row){
        sum += matrixXd.row(row);
    }

    return sum;
}

}

namespace internals{
Eigen::MatrixXd toBaryCenter(AbstractMatrix const& matrix){
    Eigen::MatrixXd matrixXd = matrixToEigenMatrixXd(matrix);
    auto mean = addAllRows(matrixXd) / static_cast<float_type>(matrixXd.rows());

    for(std::size_t row = 0u; row < matrixXd.rows(); ++row){
        matrixXd.row(row) -= mean;
    }

    return matrixXd;
}

Eigen::MatrixXd getCorrelationMatrix(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz) {
  auto new_xyz_mat = toBaryCenter(new_xyz);
  auto old_xyz_mat = toBaryCenter(old_xyz);
  return new_xyz_mat.transpose() * old_xyz_mat;
}

Eigen::MatrixXd getFmatrix(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz) {

  auto correl_mat = getCorrelationMatrix(old_xyz, new_xyz);
  auto R11 = correl_mat(0, 0);
  auto R12 = correl_mat(0, 1);
  auto R13 = correl_mat(0, 2);
  auto R21 = correl_mat(1, 0);
  auto R22 = correl_mat(1, 1);
  auto R23 = correl_mat(1, 2);
  auto R31 = correl_mat(2, 0);
  auto R32 = correl_mat(2, 1);
  auto R33 = correl_mat(2, 2);
  Eigen::MatrixXd F(4, 4);
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

std::pair<float_type, Eigen::Vector4d> getQuaternion(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz) {

  auto fMatrix = getFmatrix(old_xyz, new_xyz);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(fMatrix);

  Eigen::Vector4d quaternion = es.eigenvalues().col(3);
  if (quaternion(0) < 0.) {
    quaternion *= -1.;
  }

  return std::make_pair(es.eigenvalues()(3), quaternion);
}

Eigen::MatrixXd getRotationMatrix(Eigen::Vector4d const& q){
  Eigen::MatrixXd rotationMatrix(3, 3);
  rotationMatrix <<
    q(0)*q(0) + q(1)*q(1) - q(2)*q(2) - q(3)*q(3), 2.*(q(1)*q(2) - q(0)*q(3))                   , 2.*(q(1)*q(3) + q(0)*q(2)),
    2.*(q(1)*q(2) + q(0)*q(3))                   , q(0)*q(0) - q(1)*q(1) + q(2)*q(2) - q(3)*q(3), 2.*(q(2)*q(3) - q(0)*q(1)),
    2.*(q(1)*q(3) - q(0)*q(2))                   , 2.*(q(2)*q(3) + q(0)*q(1))                   , q(0)*q(0) - q(1)*q(1) - q(2)*q(2) + q(3)*q(3)
  ;
  return rotationMatrix;
}

Eigen::Vector3d exponential_map(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz); {

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

  auto q = getQuaternion(old_xyz, new_xyz);
  auto quat = q.second;
  auto p = get_fac(quat(0));
  return { p * quat(1), p * quat(2), p * quat(3) };
}

std::vector<std::vector<Eigen::MatrixXd>> getCorrelationMatrixDerivatives(Eigen::MatrixXd const& new_xyz); {

  std::vector<std::vector<Eigen::MatrixXd> > result;
  for (auto c = 0u; c < new_xyz.rows(); ++c) {
    std::vector<Eigen::MatrixXd> A(3, Eigen::MatrixXd::Zero(3, 3));
    for (std::size_t l = 0; l < A.size(); ++l) {
      A.at(l).row(l) = new_xyz.row(c);
    }
    result.emplace_back(A);
  }
  return result;
}

std::vector<std::vector<Eigen::MatrixXd>> getFmatrixDerivatives(AbstractMatrix const& new_xyz) {

  auto new_shift = toBaryCenter(new_xyz);
  auto dR = getCorrelationMatrixDerivatives(new_shift);
  std::vector<std::vector<Eigen::MatrixXd> > result;
  for (auto const& S : dR) {
    std::vector<Eigen::MatrixXd> Q;
    for (auto const& M : S) {
      Eigen::MatrixXd F(4, 4);
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
    result.emplace_back(std::move(Q));
  }
  return result;
}

Eigen::MatrixXd pinvJacobi(Eigen::MatrixXd const& matrix){
    Eigen::JacobiSVD<Eigen::MatrixXd> svd = matrix.jacobiSvd();

  typename Eigen::MatrixXd::Scalar tolerance = 1.e-7 * std::max(matrix.cols(), matrix.rows()) * svd.singularValues().array().abs().maxCoeff();

  return svd.matrixV() * Eigen::MatrixXd(Eigen::MatrixXd( (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0) ).diagonal()) * svd.matrixU().adjoint();
}

std::vector<Eigen::MatrixXd> quaternion_derivs(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz) {

  auto q_eigval = getQuaternion(old_xyz, new_xyz);
  auto F = getFmatrix(old_xyz, new_xyz);
  auto F_der = getFmatrixDerivatives(old_xyz); // is it really the old xyz???

  auto t = pinvJacobi((Eigen::MatrixXd::Identity(4, 4) * q_eigval.first) - F);

  std::vector<Eigen::MatrixXd> result;
  auto q = q_eigval.second;
  for (auto const& Fd : F_der) {
    Eigen::MatrixXd qtemp(3, 4);
    for (std::size_t c = 0; c < Fd.size(); ++c) {
      qtemp.row(c) = (t * Fd.at(c) * q);
    }
    result.emplace_back(std::move(qtemp));
  }
  return result;
}

std::vector<Eigen::MatrixXd> exponentialDerivatives(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz) {
  auto const fac_and_dfac = [](auto const & q0) {
    if (std::abs(q0 - 1.) < q_thres) {
      return std::make_pair(2. - 2. * (q0 - 1.) / 3., -2. / 3.);
    }
    auto acosq0 = std::acos(q0);
    auto q0_sq = 1. - q0 * q0;
    return std::make_pair(2.*acosq0 / std::sqrt(q0_sq), 2.*q0*acosq0 / std::pow(q0_sq, 1.5) - 2./q0_sq /*<-So at least in their paper (Lee Ping)*/);
  };

  auto q_val = getQuaternion(old_xyz, new_xyz);
  auto q = q_val.second;
  //std::cout << q << "    " << q0 << "\n\n";
  float_type p{ 0.0 }, d{ 0.0 };
  std::tie(p, d) = fac_and_dfac(q(0));
  //std::cout << p << " " << d << "\n\n";
  Eigen::MatrixXd dv_mat(4,3);
  dv_mat << d * q(1), d * q(2), d * q(3),
            p       , 0.      , 0.      ,
            0.      , p       , 0.      ,
            0.      , 0.      , p       ;
  //std::cout << "dv_mat:\n" << dv_mat << "\n\n";
  auto dq_mat = getQuaternionDerivatives(old_xyz, new_xyz);
  std::vector<Eigen::MatrixXd> result;
  for (auto i{ 0u }; i < old_xyz.rows(); ++i) {
    Eigen::MatrixXd R(3,3);
    auto const& dqdx_i = dq_mat.at(i);
    for(auto j{ 0u }; j < dqdx_i.rows(); ++j){
      auto const& dqdx_ij = dqdx_i.row(j);
      Eigen::MatrixXd new_row = Eigen::MatrixXd::Zero(1,3);
      for(auto k{ 0u }; k<dqdx_i.cols(); ++k){
        new_row += dv_mat.row(k) * dqdx_ij(0, k);
      }
      R.row(j) = new_row;
    }
    result.emplace_back(std::move(R));
  }
  return result;
}

}