#include "ic_core.h"
#include "ic_util.h"

using coords::float_type;

float_type ic_core::distance::dist() {
  return ic_util::euclid_dist<double>(a_, b_);
}

std::vector<scon::c3<float_type>> ic_core::distance::bond_der() {
  using std::sqrt;
  using scon::dot;

  std::vector<scon::c3<float_type>> bond_der;
  auto bond = a_ - b_;
  auto bond_dot = dot(bond, bond);
  auto length = sqrt(bond_dot);
  auto der = bond / length;
  bond_der.emplace_back(der);
  bond_der.emplace_back(-der);
  return bond_der;
}

std::vector<float_type>
ic_core::distance::bond_der_vec(const std::size_t& sys_size) {
  using scon::c3;

  auto firstder = ic_core::distance::bond_der();
  c3<float_type> temp(0.0, 0.0, 0.0);
  std::vector<c3<float_type>> der_vec(sys_size, temp);
  der_vec.at(index_a_ - 1) = firstder.front();
  der_vec.at(index_b_ - 1) = firstder.back();
  auto result = ic_util::flatten_c3_vec(der_vec);
  return result;
}

coords::float_type ic_core::angle::ang() {
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto u = a_ - b_;
  auto v = c_ - b_;
  auto uXv = cross(u, v);
  auto uDv = dot(u, v);
  auto l = len(uXv);
  return std::atan2(l, uDv) * SCON_180PI;
}

std::vector<scon::c3<float_type>> ic_core::angle::angle_der() {
  using coords::Cartesian_Point;
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto u_p = a_ - b_;
  auto v_p = c_ - b_;
  auto u = ic_util::normalize(u_p);
  auto v = ic_util::normalize(v_p);
  Cartesian_Point cp1(1.0, -1.0, 1.0);
  Cartesian_Point cp2(-1.0, 1.0, 1.0);
  Cartesian_Point w_p(0.0, 0.0, 0.0);
  auto epsilon{ 0.1 };
  if (dot(u, v) / (len(u) * len(v)) > (1 - epsilon) ||
      dot(u, v) / (len(u) * len(v)) < (-1 + epsilon)) {
    if (dot(u, cp1) / (len(u) * len(cp1)) > (1 - epsilon) ||
        dot(u, cp1) / (len(u) * len(cp1)) < (-1 + epsilon)) {
      w_p = cross(u, cp2);
    } else {
      w_p = cross(u, cp1);
    }
  } else {
    w_p = cross(u, v);
  }
  auto w = ic_util::normalize(w_p);
  auto ad0 = cross(u, w) / len(u_p);
  auto ad1 = cross(w, v) / len(v_p);
  std::vector<scon::c3<float_type>> angle_der;
  angle_der.emplace_back(ad0);
  angle_der.emplace_back(ad1);
  angle_der.emplace_back(-ad0 - ad1);
  return angle_der;
}

std::vector<float_type>
ic_core::angle::angle_der_vec(const std::size_t& sys_size) {
  using scon::c3;

  auto firstder = ic_core::angle::angle_der();
  c3<float_type> temp(0.0, 0.0, 0.0);
  std::vector<c3<float_type>> der_vec(sys_size, temp);
  der_vec.at(index_a_ - 1) = firstder.at(0);
  der_vec.at(index_b_ - 1) = firstder.at(2);
  der_vec.at(index_c_ - 1) = firstder.at(1);
  auto result = ic_util::flatten_c3_vec(der_vec);
  return result;
}

coords::float_type ic_core::dihedral::dihed() {
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto b1 = b_ - a_;
  auto b2 = c_ - b_;
  auto b3 = d_ - c_;
  auto b1Xb2 = cross(b1, b2);
  auto lb1Xb2 = len(b1Xb2);
  auto n1 = b1Xb2 / lb1Xb2;
  auto n2 = cross(b2, b3) / len(cross(b2, b3));
  auto m1 = cross(n1, ic_util::normalize(b2));
  auto x = dot(n1, n2);
  auto y = dot(m1, n2);
  return std::atan2(y, x) * SCON_180PI;
}

std::vector<scon::c3<float_type>> ic_core::dihedral::dihed_der() {
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto u_p = a_ - b_;
  auto w_p = c_ - b_;
  auto v_p = d_ - c_;
  auto u = ic_util::normalize(u_p);
  auto w = ic_util::normalize(w_p);
  auto v = ic_util::normalize(v_p);
  auto t1 = cross(u, w) / (len(u_p) * (1 - std::pow(dot(u, w), 2)));
  auto t2 = cross(v, w) / (len(v_p) * (1 - std::pow(dot(v, w), 2)));
  auto cuwd = cross(u, w) * dot(u, w);
  auto t3 = cuwd / (len(w_p) * (1 - std::pow(dot(u, w), 2)));
  auto cvwd = cross(v, w) * dot(v, w);
  auto t4 = cvwd / (len(w_p) * (1 - std::pow(dot(u, w), 2)));
  std::vector<scon::c3<float_type>> dih_der;
  dih_der.emplace_back(t1);
  dih_der.emplace_back(-t2);
  dih_der.emplace_back(-t1 + t2 - t4);
  dih_der.emplace_back(t2 - t3 + t4);
  return dih_der;
}

std::vector<float_type>
ic_core::dihedral::dihed_der_vec(const std::size_t& sys_size) {
  using scon::c3;

  auto firstder = ic_core::dihedral::dihed_der();
  c3<float_type> temp(0.0, 0.0, 0.0);
  std::vector<c3<float_type>> der_vec(sys_size, temp);
  der_vec.at(index_a_ - 1) = firstder.at(0);
  der_vec.at(index_b_ - 1) = firstder.at(2);
  der_vec.at(index_c_ - 1) = firstder.at(3);
  der_vec.at(index_d_ - 1) = firstder.at(1);
  auto result = ic_util::flatten_c3_vec(der_vec);
  return result;
}

coords::float_type ic_core::out_of_plane::oop() {
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto b1 = b_ - a_;
  auto b2 = c_ - b_;
  auto b3 = d_ - c_;
  auto b1Xb2 = cross(b1, b2);
  auto lb1Xb2 = len(b1Xb2);
  auto n1 = b1Xb2 / lb1Xb2;
  auto n2 = cross(b2, b3) / len(cross(b2, b3));
  auto m1 = cross(n1, ic_util::normalize(b2));
  auto x = dot(n1, n2);
  auto y = dot(m1, n2);
  return std::atan2(y, x) * SCON_180PI;
}

std::vector<scon::c3<float_type>> ic_core::out_of_plane::oop_der() {
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto u_p = a_ - b_;
  auto w_p = c_ - b_;
  auto v_p = d_ - c_;
  auto u = ic_util::normalize(u_p);
  auto w = ic_util::normalize(w_p);
  auto v = ic_util::normalize(v_p);
  auto t1 = cross(u, w) / (len(u_p) * (1 - std::pow(dot(u, w), 2)));
  auto t2 = cross(v, w) / (len(v_p) * (1 - std::pow(dot(v, w), 2)));
  auto cuwd = cross(u, w) * dot(u, w);
  auto t3 = cuwd / (len(w_p) * (1 - std::pow(dot(u, w), 2)));
  auto cvwd = cross(v, w) * dot(v, w);
  auto t4 = cvwd / (len(w_p) * (1 - std::pow(dot(u, w), 2)));
  std::vector<scon::c3<float_type>> oop_der;
  oop_der.emplace_back(t1);
  oop_der.emplace_back(-t2);
  oop_der.emplace_back(-t1 + t2 - t4);
  oop_der.emplace_back(t2 - t3 + t4);
  return oop_der;
}

std::vector<float_type>
ic_core::out_of_plane::oop_der_vec(const std::size_t& sys_size) {
  using scon::c3;

  auto firstder = ic_core::out_of_plane::oop_der();
  c3<float_type> temp(0.0, 0.0, 0.0);
  std::vector<c3<float_type>> der_vec(sys_size, temp);
  der_vec.at(index_a_ - 1) = firstder.at(0);
  der_vec.at(index_b_ - 1) = firstder.at(2);
  der_vec.at(index_c_ - 1) = firstder.at(3);
  der_vec.at(index_d_ - 1) = firstder.at(1);
  auto result = ic_util::flatten_c3_vec(der_vec);
  return result;
}

std::vector<float_type>
ic_core::trans_x::trans_x_der_vec() {
  using cp = coords::Cartesian_Point;

  return ic_util::flatten_c3_vec(trans_der([](auto const & s) {
    return cp(1. / static_cast<float_type>(s), 0., 0.);
  }));
}

std::vector<float_type>
ic_core::trans_y::trans_y_der_vec() {
  
  using cp = coords::Cartesian_Point;

  return ic_util::flatten_c3_vec(trans_der([](auto const & s) {
    return cp(0., 1. / static_cast<float_type>(s), 0.);
  }));
}

std::vector<float_type>
ic_core::trans_z::trans_z_der_vec() {
  using cp = coords::Cartesian_Point;

  return ic_util::flatten_c3_vec(trans_der([](auto const & s) {
    return cp(0., 0., 1. / static_cast<float_type>(s));
  }));
}

std::array<float_type, 3u>
ic_core::rotation::rot_val(const coords::Representation_3D& trial) {
  auto result = ic_rotation::exponential_map(trial, reference_);
  return result;
}

std::vector<scon::mathmatrix<float_type>>
ic_core::rotation::rot_der(const coords::Representation_3D& trial) {
  return ic_rotation::exponential_derivs(trial, reference_);
}

scon::mathmatrix<float_type>
ic_core::rotation::rot_der_mat(const std::size_t& sys_size,
                               const coords::Representation_3D& trial) {
  using Mat = scon::mathmatrix<float_type>;
  auto const & zero = scon::mathmatrix<float_type>::zero;

  Mat X = zero(sys_size, 3);
  Mat Y = zero(sys_size, 3);
  Mat Z = zero(sys_size, 3);
  auto first_ders = ic_core::rotation::rot_der(trial);
  std::size_t index{ 1 };
  for (auto const & i : first_ders) {
    auto val_index = indices_.at(index - 1) - 1;
    X.row(val_index) = i.row(0);
    Y.row(val_index) = i.row(1);
    Z.row(val_index) = i.row(2);
    ++index;
  }
  Mat result(3, sys_size * 3);
  result.row(0) = X.vectorise_row();
  result.row(1) = Y.vectorise_row();
  result.row(2) = Z.vectorise_row();
  return result;
}

coords::float_type
ic_core::rotation::radius_gyration(const coords::Representation_3D& struc) {
  coords::float_type result = ic_util::rad_gyr(struc);
  return result;
}

std::vector<ic_core::trans_x> ic_core::system::create_trans_x(
    const std::vector<coords::Representation_3D>& rep_vec,
    const std::vector<std::vector<std::size_t>>& index_vec) {
  std::vector<ic_core::trans_x> result;
  auto it1 = begin(rep_vec);
  auto end_it1 = end(rep_vec);
  auto it2 = begin(index_vec);
  for (; it1 != end_it1; ++it1, ++it2) {
    ic_core::trans_x temp(*it1, *it2);
    result.emplace_back(temp);
  }
  return result;
}

std::vector<ic_core::trans_y> ic_core::system::create_trans_y(
    const std::vector<coords::Representation_3D>& rep_vec,
    const std::vector<std::vector<std::size_t>>& index_vec) {
  std::vector<ic_core::trans_y> result;
  auto it1 = begin(rep_vec);
  auto end_it1 = end(rep_vec);
  auto it2 = begin(index_vec);
  for (; it1 != end_it1; ++it1, ++it2) {
    ic_core::trans_y temp(*it1, *it2);
    result.emplace_back(temp);
  }
  return result;
}

std::vector<ic_core::trans_z> ic_core::system::create_trans_z(
    const std::vector<coords::Representation_3D>& rep_vec,
    const std::vector<std::vector<std::size_t>>& index_vec) {
  std::vector<ic_core::trans_z> result;
  auto it1 = begin(rep_vec);
  auto end_it1 = end(rep_vec);
  auto it2 = begin(index_vec);
  for (; it1 != end_it1; ++it1, ++it2) {
    ic_core::trans_z temp(*it1, *it2);
    result.emplace_back(temp);
  }
  return result;
}

std::vector<ic_core::rotation> ic_core::system::create_rotations(
    const std::vector<coords::Representation_3D>& rep_vec,
    const std::vector<std::vector<std::size_t>>& index_vec) {
  std::vector<ic_core::rotation> result;
  auto it1 = begin(rep_vec);
  auto end_it1 = end(rep_vec);
  auto it2 = begin(index_vec);
  for (; it1 != end_it1; ++it1, ++it2) {
    ic_core::rotation temp(*it1, *it2);
    result.emplace_back(temp);
  }
  return result;
}

std::pair<scon::mathmatrix<float_type>, scon::mathmatrix<float_type>>
ic_core::system::delocalize_ic_system(const coords::Representation_3D& trial) {

  using Mat = scon::mathmatrix<float_type>;

  auto const & sys_size = trial.size();

  std::vector<std::vector<float_type>> result;
  for (auto& i : trans_x_vec_) {
    result.emplace_back(i.trans_x_der_vec());
  }
  for (auto& i : trans_y_vec_) {
    result.emplace_back(i.trans_y_der_vec());
  }
  for (auto& i : trans_z_vec_) {
    result.emplace_back(i.trans_z_der_vec());
  }
  for (auto& i : distance_vec_) {
    result.emplace_back(i.bond_der_vec(sys_size));
  }
  for (auto& i : angle_vec_) {
    result.emplace_back(i.angle_der_vec(sys_size));
  }
  for (auto& i : dihed_vec_) {
    result.emplace_back(i.dihed_der_vec(sys_size));
  }
  for (auto& i : oop_vec_) {
    result.emplace_back(i.oop_der_vec(sys_size));
  }
  for (auto& i : rotation_vec_) {
    auto temp = i.rot_der_mat(sys_size, trial);
    auto vec_X = temp.col_to_std_vector(0);
    auto vec_Y = temp.col_to_std_vector(1);
    auto vec_Z = temp.col_to_std_vector(2);
    result.emplace_back(vec_X);
    result.emplace_back(vec_Y);
    result.emplace_back(vec_Z);
  }
  std::size_t n_cols = result.size();
  std::size_t n_rows = result.at(0).size();
  Mat B_matrix_t(n_rows, n_cols);
  for (std::size_t h = 0; h < n_cols; ++h) {
    Mat col_vec = Mat::col_from_vec(result.at(h));
    B_matrix_t.col(h) = col_vec;
  }
  auto B_matrix = B_matrix_t.t();
  auto G_matrix = B_matrix * B_matrix.t();
  Mat eigval, eigvec;
  std::tie(eigval, eigvec) = G_matrix.eigensym();
  auto row_index_vec = eigval.sort_idx();
  auto col_index_vec = eigval.find_idx([](float_type const & a) {
    return a > 1e-6;
  });
  Mat del_mat = eigvec.submat(row_index_vec, col_index_vec);
  return std::make_pair(del_mat, G_matrix);
}

scon::mathmatrix<float_type> ic_core::system::initial_hessian() {
  using ic_atom::period_one;
  using ic_atom::period_two;
  using ic_atom::period_three;
  using Mat = scon::mathmatrix<float_type>;
  using scon::cross;
  using scon::dot;
  using scon::len;

  std::vector<float_type> values;
  for (auto& i : distance_vec_) {
    bool a_period_one = std::find(period_one.begin(), period_one.end(),
                                  i.elem_a_) != period_one.end();
    bool a_period_two = std::find(period_two.begin(), period_two.end(),
                                  i.elem_a_) != period_two.end();
    bool a_period_three = std::find(period_three.begin(), period_three.end(),
                                    i.elem_a_) != period_three.end();
    bool b_period_one = std::find(period_one.begin(), period_one.end(),
                                  i.elem_b_) != period_one.end();
    bool b_period_two = std::find(period_two.begin(), period_two.end(),
                                  i.elem_b_) != period_two.end();
    bool b_period_three = std::find(period_three.begin(), period_three.end(),
                                    i.elem_b_) != period_three.end();
    auto B_val{ 0.0 };
    if (a_period_one && b_period_one) {
      B_val = -0.244;
    } else if (a_period_one && b_period_two || b_period_one && a_period_two) {
      B_val = 0.352;
    } else if (a_period_two && b_period_two) {
      B_val = 1.085;
    } else if (a_period_one && b_period_three ||
               b_period_one && a_period_three) {
      B_val = 0.660;
    } else if (a_period_two && b_period_three ||
               b_period_two && a_period_three) {
      B_val = 1.522;
    } else if (a_period_three && b_period_three) {
      B_val = 2.068;
    }
    auto A_val{ 1.734 };
    auto dist_bohr = i.dist() / ic_atom::bohr;
    auto temp = std::pow(dist_bohr - B_val, 3);
    values.emplace_back(A_val / temp);
  }
  for (auto& i : angle_vec_) {
    if (i.elem_a_ == "H" || i.elem_c_ == "H") {
      values.emplace_back(0.160);
    } else {
      values.emplace_back(0.250);
    }
  }
  values.insert(values.end(), dihed_vec_.size(), 0.023);
  for (auto& i : oop_vec_) {
    auto a_bohr = i.a_ / ic_atom::bohr;
    auto b_bohr = i.b_ / ic_atom::bohr;
    auto c_bohr = i.c_ / ic_atom::bohr;
    auto d_bohr = i.d_ / ic_atom::bohr;
    auto r1 = b_bohr - a_bohr;
    auto r2 = b_bohr - c_bohr;
    auto r3 = b_bohr - d_bohr;
    auto r2Xr3 = cross(r2, r3);
    auto rd = dot(r1, r2Xr3);
    auto t2 = rd / (len(r1) * len(r2) * len(r3));
    auto d = 1 - t2;
    auto d_pow = std::pow(d, 4);
    values.emplace_back(0.045 * d_pow);
  }
  values.insert(values.end(), trans_x_vec_.size() * 3, 0.05);
  values.insert(values.end(), rotation_vec_.size() * 3, 0.05);

  auto val_col = Mat::col_from_vec(values);
  //val_col.print(); //test output
  Mat hessian = val_col.diagmat();
  return hessian;
}

scon::mathmatrix<float_type>
ic_core::system::delocalize_hessian(const scon::mathmatrix<float_type>& del_mat,
                                    const scon::mathmatrix<float_type>& hessian) {
  using Mat = scon::mathmatrix<float_type>;

  Mat del_hessian = del_mat.t() * hessian * del_mat;

  return del_hessian;
}

scon::mathmatrix<float_type>
ic_core::system::G_mat_inversion(const scon::mathmatrix<float_type>& G_matrix) {

  using Mat = scon::mathmatrix<float_type>;

  Mat U, s, V;
  G_matrix.singular_value_decomposition(U, s, V);
  //remplace values smaller than 1e-6 with 0.0
  s.replace_idx_with(s.find_idx([](float_type const & a) {
    return a <= 1e-6;
  }), 0.0);
  auto index = s.find_idx([](float_type const & a) {
    return a > 1e-6;
  });
  for (auto && el : s.elem(index)) {
    el.get() = 1. / el.get();
  }
  auto G_mat_inv = U * s.diagmat() * V.t();
  return G_mat_inv;
}
