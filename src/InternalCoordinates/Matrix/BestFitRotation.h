#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_MATRIX_BESTFITROTATION_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_MATRIX_BESTFITROTATION_H_

#ifndef CAST_INTERNALCOORDINATES_BESTFITROTATION_H_
#define CAST_INTERNALCOORDINATES_BESTFITROTATION_H_

#include "AbstractMatrix.h"

#include<Eigen/Dense>


namespace internals {

float_type constexpr q_thres{ 1e-6 };

Eigen::MatrixXd toBaryCenter(AbstractMatrix const&);
Eigen::MatrixXd getCorrelationMatrix(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz);
Eigen::MatrixXd getFmatrix(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz);

std::pair<float_type, Eigen::Vector4d> getQuaternion(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz);

Eigen::MatrixXd getRotationMatrix(Eigen::Vector4d const& q);

Eigen::Vector3d exponential_map(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz);

std::vector<std::vector<Eigen::MatrixXd>> getCorrelationMatrixDerivatives(Eigen::MatrixXd const& new_xyz);

std::vector<std::vector<Eigen::MatrixXd>> getFmatrixDerivatives(AbstractMatrix const& new_xyz);

std::vector<Eigen::MatrixXd> getQuaternionDerivatives(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz);

std::vector<Eigen::MatrixXd> exponentialDerivatives(AbstractMatrix const& old_xyz, AbstractMatrix const& new_xyz);
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
#endif // CAST_INTERNALCOORDINATES_BESTFITROTATION_H_


#endif