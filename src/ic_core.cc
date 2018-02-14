#include "ic_core.h"
#include "ic_util.h"
#include<iterator>

using coords::float_type;

float_type ic_core::distance::dist() {
  return ic_util::euclid_dist<double>(a_, b_);
}

std::pair<scon::c3<float_type>, scon::c3<float_type>> ic_core::distance::bond_der() {
  auto bond = ic_util::normalize(a_ - b_);
  return std::make_pair(bond,-bond);
}

std::vector<float_type>
ic_core::distance::bond_der_vec(const std::size_t& sys_size) {
  using scon::c3;

  auto firstder = bond_der();
  std::vector<c3<float_type>> der_vec(sys_size, c3<float_type>(0.0, 0.0, 0.0));
  der_vec.at(index_a_ - 1) = firstder.first;
  der_vec.at(index_b_ - 1) = firstder.second;
  return ic_util::flatten_c3_vec(der_vec);
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

std::tuple<scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>> 
ic_core::angle::angle_der() const {
  using coords::Cartesian_Point;
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto u = a_ - b_;
  auto v = c_ - b_;
  auto lu = len(u);
  auto lv = len(v);
  u = ic_util::normalize(u);
  v = ic_util::normalize(v);
  auto cp1 = ic_util::normalize(Cartesian_Point(1.0, -1.0, 1.0));
  auto cp2 = ic_util::normalize(Cartesian_Point(-1.0, 1.0, 1.0));
  Cartesian_Point w_p(0.0, 0.0, 0.0);
  auto epsilon{ 0.1 };
  if (dot(u, v) > (1. - epsilon) ||
      dot(u, v) < (-1. + epsilon)) {
    if (dot(u, cp1) > (1. - epsilon) ||
        dot(u, cp1) < (-1. + epsilon)) {
      w_p = cross(u, cp2);
    } else {
      w_p = cross(u, cp1);
    }
  } else {
    w_p = cross(u, v);
  }
  w_p = ic_util::normalize(w_p);
  auto ad0 = cross(u, w_p) / lu;
  auto ad1 = cross(w_p, v) / lv;
  return std::make_tuple(ad0, ad1, -ad0 - ad1);
}

std::vector<float_type>
ic_core::angle::angle_der_vec(const std::size_t& sys_size) {
  using scon::c3;

  auto firstder = ic_core::angle::angle_der();
  std::vector<c3<float_type>> der_vec(sys_size, c3<float_type>(0.0, 0.0, 0.0));
  der_vec.at(index_a_ - 1) = std::get<0>(firstder);
  der_vec.at(index_b_ - 1) = std::get<2>(firstder);
  der_vec.at(index_c_ - 1) = std::get<1>(firstder);
  return ic_util::flatten_c3_vec(der_vec);
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

std::tuple<scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>> 
ic_core::dihedral::dihed_der() const{
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
  return std::make_tuple(t1, -t2, -t1 + t2 - t4, t2 - t3 + t4);
}

std::vector<float_type>
ic_core::dihedral::dihed_der_vec(const std::size_t& sys_size) {
  using scon::c3;

  auto firstder = ic_core::dihedral::dihed_der();
  c3<float_type> temp(0.0, 0.0, 0.0);
  std::vector<c3<float_type>> der_vec(sys_size, temp);
  der_vec.at(index_a_ - 1) = std::get<0>(firstder);
  der_vec.at(index_b_ - 1) = std::get<2>(firstder);
  der_vec.at(index_c_ - 1) = std::get<3>(firstder);
  der_vec.at(index_d_ - 1) = std::get<1>(firstder);
  return ic_util::flatten_c3_vec(der_vec);
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
ic_core::trans_x::trans_der_vec(std::size_t const & system_size) const{
  using cp = coords::Cartesian_Point;

  return ic_util::flatten_c3_vec(trans_der(system_size, [](auto const & s) {
    return cp(1. / static_cast<float_type>(s), 0., 0.);
  }));
}

std::vector<float_type>
ic_core::trans_y::trans_der_vec(std::size_t const & system_size) const{
  
  using cp = coords::Cartesian_Point;

  return ic_util::flatten_c3_vec(trans_der(system_size, [](auto const & s) {
    return cp(0., 1. / static_cast<float_type>(s), 0.);
  }));
}

std::vector<float_type>
ic_core::trans_z::trans_der_vec(std::size_t const & system_size) const{
  using cp = coords::Cartesian_Point;

  return ic_util::flatten_c3_vec(trans_der(system_size, [](auto const & s) {
    return cp(0., 0., 1. / static_cast<float_type>(s));
  }));
}

std::array<float_type, 3u>
ic_core::rotation::rot_val(const coords::Representation_3D& trial) {
  coords::Representation_3D curr_xyz;
  curr_xyz.reserve(indices_.size());
  for (auto const & i : indices_) {
    curr_xyz.emplace_back(trial.at(i - 1));
  }
  auto result = ic_rotation::exponential_map(curr_xyz, reference_);
  return result;
}

std::vector<scon::mathmatrix<float_type>>
ic_core::rotation::rot_der(const coords::Representation_3D& trial) {
  return ic_rotation::exponential_derivs(trial, reference_);
}

scon::mathmatrix<float_type>
ic_core::rotation::rot_der_mat(std::size_t const & sys_size, const coords::Representation_3D& trial) {
  using Mat = scon::mathmatrix<float_type>;
  auto const & zero = scon::mathmatrix<float_type>::zero;

  Mat X = zero(sys_size, 3);
  Mat Y = zero(sys_size, 3);
  Mat Z = zero(sys_size, 3);
  coords::Representation_3D curr_xyz;
  curr_xyz.reserve(indices_.size());
  for (auto const & i : indices_) {
    curr_xyz.emplace_back(trial.at(i-1));
  }
  auto first_ders = rot_der(curr_xyz);
  std::size_t index{ 1 };
  for (auto const & i : first_ders) {
    auto val_index = indices_.at(index - 1) - 1;
    X.set_row(val_index, i.row(0));
    Y.set_row(val_index, i.row(1));
    Z.set_row(val_index, i.row(2));
    ++index;
  }
  Mat result(sys_size * 3, 3);
  auto xv = X.vectorise_row();
  std::cout << xv << std::endl;
  result.set_col(0, X.vectorise_row());
  result.set_col(1, Y.vectorise_row());
  result.set_col(2, Z.vectorise_row());
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

  /*auto const & trans_derivs = [](auto const & trans_vec) {
    std::vector<float_type> temp;
    for (auto const & trans : trans_vec) {
      auto trans_temp = trans.trans_der_vec();
      temp.insert(std::begin(temp), std::make_move_iterator(std::begin(trans_temp)), std::make_move_iterator(std::end(trans_temp)));
    }
    return temp;
  };*/

  auto const & sys_size = trial.size();

  std::vector<std::vector<float_type>> result;
  
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
  /*result.emplace_back(trans_derivs(trans_x_vec_));
  result.emplace_back(trans_derivs(trans_y_vec_));
  result.emplace_back(trans_derivs(trans_z_vec_));*/
  for (auto& i : trans_x_vec_) {
    result.emplace_back(i.trans_der_vec(sys_size));
  }
  for (auto& i : trans_y_vec_) {
    result.emplace_back(i.trans_der_vec(sys_size));
  }
  for (auto& i : trans_z_vec_) {
    result.emplace_back(i.trans_der_vec(sys_size));
  }
  for (auto& i : rotation_vec_) {
    auto temp = i.rot_der_mat(sys_size, trial);
    //look for th cols and rows!
    result.emplace_back(temp.col_to_std_vector(0));
    result.emplace_back(temp.col_to_std_vector(1));
    result.emplace_back(temp.col_to_std_vector(2));
  }
  std::size_t n_cols = result.size();
  std::size_t n_rows = result.at(0).size();
  Mat B_matrix_t(n_rows, n_cols);
  for (std::size_t h = 0; h < n_cols; ++h) {
    std::cout << h << std::endl;
    B_matrix_t.set_col(h, Mat::col_from_vec(result.at(h)));
  }
  auto B_mat = B_matrix_t.t();
  std::ofstream o_file("B_mat_t.dat");
  for (auto i = 0; i < B_mat.rows(); ++i) {
    for (auto j = 0; j < B_mat.cols(); ++j) {
      o_file << std::setw(10) << std::setprecision(5) << std::fixed << B_mat(i, j);
    }
    o_file << "\n";
  }
  std::cout << B_matrix_t << "\n";
  auto G_matrix = B_matrix_t.t() * B_matrix_t;
  Mat eigval, eigvec;
  std::tie(eigval, eigvec) = G_matrix.eigensym(true);
  auto eigval_vec = eigval.col_to_std_vector();
  auto row_index_vec = eigval.sort_idx();
  auto col_index_vec = eigval.find_idx([](float_type const & a) {
    return std::abs(a) > 1e-6;
  });
  Mat del_mat = eigvec.submat(row_index_vec, col_index_vec);
  std::cout << del_mat << "\n";
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
