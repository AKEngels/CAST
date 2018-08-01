#include "TranslationRotationInternalCoordinates.h"

namespace internals {

  scon::mathmatrix<coords::float_type>& TRIC::delocalize_ic_system() {
    using Mat = scon::mathmatrix<coords::float_type>;

    Mat eigval, eigvec;
    std::tie(eigval, eigvec) = system::Gmat().eigensym(false);

    auto row_index_vec = eigval.sort_idx();
    auto col_index_vec = eigval.find_idx([](coords::float_type const & a) {
      return std::abs(a) > 1e-6;
    });

    del_mat = eigvec.submat(row_index_vec, col_index_vec);
    new_B_matrix = new_G_matrix = true; //B and G got to be calculated for the new ic_system
    return del_mat;
  }


  scon::mathmatrix<coords::float_type>& TRIC::Bmat() {
    if (!new_B_matrix) {
      return B_matrix;
    }
    B_matrix = del_mat.t()*system::Bmat();
    new_B_matrix = false;
    return B_matrix;
  }

  scon::mathmatrix<coords::float_type>& TRIC::Gmat() {
    if (!new_G_matrix) {
      return G_matrix;
    }
    Bmat();
    G_matrix = B_matrix * B_matrix.t();
    new_G_matrix = false;
    return G_matrix;
  }

  scon::mathmatrix<coords::float_type>& TRIC::guess_hessian() {
    hessian = del_mat.t() * system::guess_hessian() * del_mat;
    return hessian;
  }

  scon::mathmatrix<coords::float_type> TRIC::calc(coords::Representation_3D const& xyz) const {
    auto prims = system::calc(xyz);
    return (prims * del_mat).t();
  }

  scon::mathmatrix<coords::float_type> TRIC::calc_diff(coords::Representation_3D const& lhs, coords::Representation_3D const& rhs) const {
    auto diff = system::calc_diff(lhs, rhs);
    return (diff * del_mat).t();
  }
}