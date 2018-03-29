#include "ic_core.h"
#include "ic_util.h"
#include "energy.h"
#include<iterator>

using coords::float_type;

coords::Representation_3D ic_core::grads_to_bohr(coords::Representation_3D const& grads){
  coords::Representation_3D bohr_grads;
  bohr_grads.reserve(grads.size());
  for(auto const& g: grads){
    bohr_grads.emplace_back(g / energy::Hartree_Bohr2Kcal_MolAng);
  }
  return bohr_grads;
}

float_type ic_core::distance::val(coords::Representation_3D const& xyz) const {

  auto const& a = xyz.at(index_a_ - 1u);
  auto const& b = xyz.at(index_b_ - 1u);

  return ic_util::euclid_dist<double>(a, b);
}

std::pair<scon::c3<float_type>, scon::c3<float_type>> ic_core::distance::der(coords::Representation_3D const& xyz) const {

  auto const& a = xyz.at(index_a_ - 1u);
  auto const& b = xyz.at(index_b_ - 1u);

  auto bond = ic_util::normalize(a - b);
  return std::make_pair(bond,-bond);
}

std::vector<float_type>
ic_core::distance::der_vec(coords::Representation_3D const& xyz) const {
  using scon::c3;

  auto firstder = der(xyz);
  std::vector<c3<float_type>> der_vec(xyz.size(), c3<float_type>(0.0, 0.0, 0.0));
  der_vec.at(index_a_ - 1) = firstder.first;
  der_vec.at(index_b_ - 1) = firstder.second;
  return ic_util::flatten_c3_vec(der_vec);
}

float_type ic_core::distance::hessian_guess(coords::Representation_3D const& xyz) const{
  using ic_atom::element_period;
  using ic_atom::period;

  auto el_a = element_period(elem_a_);
  auto el_b = element_period(elem_b_);

  auto B_val{ 0.0 };
  if (el_a == period::one &&
      el_b == period::two) {
    B_val = -0.244;
  }
  else if ((el_a == period::one && el_b == period::two) ||
           (el_a == period::two && el_b == period::one)) {
    B_val = 0.352;
  }
  else if (el_a == period::two && el_b == period::two) {
    B_val = 1.085;
  }
  else if ((el_a == period::one && el_b == period::three) ||
           (el_a == period::three && el_b == period::one)) {
    B_val = 0.660;
  }
  else if ((el_a == period::two && el_b == period::three) ||
           (el_a == period::three && el_b == period::two)) {
    B_val = 1.522;
  }
  else if (el_a == period::three && el_b == period::three) {
    B_val = 2.068;
  }
  auto A_val{ 1.734 };
  auto temp = std::pow(val(xyz) - B_val, 3);
  return A_val / temp;
}

std::string ic_core::distance::info(coords::Representation_3D const & xyz) const
{
  std::ostringstream oss;
  oss << val(xyz) << "||" << index_a_ << "||" << index_b_ << "\n";
  return oss.str();
}

coords::float_type ic_core::angle::val(coords::Representation_3D const& xyz) const{
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto const& a = xyz.at(index_a_ - 1u);
  auto const& b = xyz.at(index_b_ - 1u);
  auto const& c = xyz.at(index_c_ - 1u);

  auto u = a - b;
  auto v = c - b;
  auto uXv = cross(u, v);
  auto uDv = dot(u, v);
  auto l = len(uXv);
  return std::atan2(l, uDv);
}

std::tuple<scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>>
ic_core::angle::der(coords::Representation_3D const& xyz) const {
  using coords::Cartesian_Point;
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto const& a = xyz.at(index_a_ - 1u);
  auto const& b = xyz.at(index_b_ - 1u);
  auto const& c = xyz.at(index_c_ - 1u);

  auto u = a - b;
  auto v = c - b;
  auto lu = len(u);
  auto lv = len(v);
  u = ic_util::normalize(u);
  v = ic_util::normalize(v);
  auto cp1 = ic_util::normalize(Cartesian_Point(1.0, -1.0, 1.0));
  auto cp2 = ic_util::normalize(Cartesian_Point(-1.0, 1.0, 1.0));
  Cartesian_Point w_p(0.0, 0.0, 0.0);
  auto constexpr epsilon{ 0.01 };
  if (std::fabs(dot(u, v)) > (1. - epsilon)) {
    if (std::fabs(dot(u, cp1)) > (1. - epsilon)) {
      w_p = cross(u, cp2);
    } else {
      w_p = cross(u, cp1);
    }
  } else {
    w_p = cross(u, v);
  }
  auto w = ic_util::normalize(w_p);
  auto ad0 = cross(u, w) / lu;
  auto ad1 = cross(w, v) / lv;
  //                      a     b           c
  //                      m     o           n
  return std::make_tuple(ad0, -ad0 - ad1, ad1);
}

std::vector<float_type>
ic_core::angle::der_vec(coords::Representation_3D const& xyz) const {
  using scon::c3;

  auto firstder = ic_core::angle::der(xyz);
  std::vector<c3<float_type>> der_vec(xyz.size(), c3<float_type>(0.0, 0.0, 0.0));
  der_vec.at(index_a_ - 1) = std::get<0>(firstder);
  der_vec.at(index_b_ - 1) = std::get<1>(firstder);
  der_vec.at(index_c_ - 1) = std::get<2>(firstder);
  return ic_util::flatten_c3_vec(der_vec);
}

float_type ic_core::angle::hessian_guess(coords::Representation_3D const& /*xyz*/) const
{
  if (elem_a_ == "H" || elem_c_ == "H") {
    return 0.160;
  }
  else {
    return 0.250;
  }
}

std::string ic_core::angle::info(coords::Representation_3D const & xyz) const
{
  std::ostringstream oss;
  oss << val(xyz) * SCON_180PI << "||" << index_a_ << "||" << index_b_ << "||" << index_c_ << "\n";
  return oss.str();
}

coords::float_type ic_core::dihedral::val(coords::Representation_3D const& xyz) const {
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto const & a = xyz.at(index_a_ - 1u);
  auto const & b = xyz.at(index_b_ - 1u);
  auto const & c = xyz.at(index_c_ - 1u);
  auto const & d = xyz.at(index_d_ - 1u);

  auto b1 = b - a;
  auto b2 = c - b;
  auto b3 = d - c;
  auto b1Xb2 = cross(b1, b2);
  auto lb1Xb2 = len(b1Xb2);
  auto n1 = b1Xb2 / lb1Xb2;
  auto n2 = cross(b2, b3) / len(cross(b2, b3));
  auto m1 = cross(n1, ic_util::normalize(b2));
  auto x = dot(n1, n2);
  auto y = dot(m1, n2);
  return std::atan2(y, x);
}

std::tuple<scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>>
ic_core::dihedral::der(coords::Representation_3D const& xyz) const {
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto const & a = xyz.at(index_a_ - 1u);
  auto const & b = xyz.at(index_b_ - 1u);
  auto const & c = xyz.at(index_c_ - 1u);
  auto const & d = xyz.at(index_d_ - 1u);

  auto u_p = a - b;
  auto w_p = c - b;
  auto v_p = d - c;
  auto u = ic_util::normalize(u_p);
  auto w = ic_util::normalize(w_p);
  auto v = ic_util::normalize(v_p);
  auto sin2_u = 1. - std::pow(dot(u, w), 2.);
  auto sin2_v = 1. - std::pow(dot(v, w), 2.);

  auto t1 = cross(u, w) / (len(u_p) * sin2_u);
  auto t2 = cross(v, w) / (len(v_p) * sin2_v);
  auto cuwd = cross(u, w) * dot(u, w);
  auto t3 = cuwd / (len(w_p) * sin2_u);
  //changed v's sign according to J. Chem. Ohys Vol 117 No. 20, 22 2002 p. 9160-9174
  auto cvwd = cross(v, w) * dot(-v, w);
  auto t4 = cvwd / (len(w_p) * sin2_v);
  //exchanged +t2 with +t3 in o's derivative
  //                      a   b               c             d
  //                      m   o               p             n
  return std::make_tuple(t1, -t1 + t3 - t4, t2 - t3 + t4, -t2);
}

std::vector<float_type>
ic_core::dihedral::der_vec(coords::Representation_3D const& xyz) const {
  using scon::c3;

  auto firstder = der(xyz);
  c3<float_type> temp(0.0, 0.0, 0.0);
  std::vector<c3<float_type>> der_vec(xyz.size(), temp);
  der_vec.at(index_a_ - 1) = std::get<0>(firstder);
  der_vec.at(index_b_ - 1) = std::get<1>(firstder);
  der_vec.at(index_c_ - 1) = std::get<2>(firstder);
  der_vec.at(index_d_ - 1) = std::get<3>(firstder);
  return ic_util::flatten_c3_vec(der_vec);
}

float_type ic_core::dihedral::hessian_guess(coords::Representation_3D const & /*xyz*/) const
{
  return 0.023;
}

std::string ic_core::dihedral::info(coords::Representation_3D const & xyz) const
{
  std::ostringstream oss;
  oss << val(xyz) * SCON_180PI << "||" << index_a_ << "||" << index_b_ << "||" << index_c_ << "||" << index_d_ << "\n";
  return oss.str();
}

coords::float_type ic_core::out_of_plane::val(coords::Representation_3D const& xyz) const {
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto const & a = xyz.at(index_a_ - 1u);
  auto const & b = xyz.at(index_b_ - 1u);
  auto const & c = xyz.at(index_c_ - 1u);
  auto const & d = xyz.at(index_d_ - 1u);

  auto b1 = b - a;
  auto b2 = c - b;
  auto b3 = d - c;
  auto b1Xb2 = cross(b1, b2);
  auto lb1Xb2 = len(b1Xb2);
  auto n1 = b1Xb2 / lb1Xb2;
  auto n2 = cross(b2, b3) / len(cross(b2, b3));
  auto m1 = cross(n1, ic_util::normalize(b2));
  auto x = dot(n1, n2);
  auto y = dot(m1, n2);
  return std::atan2(y, x);
}

std::vector<scon::c3<float_type>> ic_core::out_of_plane::der(coords::Representation_3D const& xyz) const{
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto const & a = xyz.at(index_a_ - 1u);
  auto const & b = xyz.at(index_b_ - 1u);
  auto const & c = xyz.at(index_c_ - 1u);
  auto const & d = xyz.at(index_d_ - 1u);

  auto u_p = a - b;
  auto w_p = c - b;
  auto v_p = d - c;
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
ic_core::out_of_plane::der_vec(coords::Representation_3D const& xyz) const {
  using scon::c3;

  auto firstder = ic_core::out_of_plane::der(xyz);
  c3<float_type> temp(0.0, 0.0, 0.0);
  std::vector<c3<float_type>> der_vec(xyz.size(), temp);
  der_vec.at(index_a_ - 1) = firstder.at(0);
  der_vec.at(index_b_ - 1) = firstder.at(2);
  der_vec.at(index_c_ - 1) = firstder.at(3);
  der_vec.at(index_d_ - 1) = firstder.at(1);
  auto result = ic_util::flatten_c3_vec(der_vec);
  return result;
}

float_type ic_core::out_of_plane::hessian_guess(coords::Representation_3D const & xyz) const
{
  auto const& a = xyz.at(index_a_ - 1u);
  auto const& b = xyz.at(index_b_ - 1u);
  auto const& c = xyz.at(index_c_ - 1u);
  auto const& d = xyz.at(index_d_ - 1u);
  auto r1 = b - a;
  auto r2 = b - c;
  auto r3 = b - d;
  auto r2Xr3 = cross(r2, r3);
  auto rd = dot(r1, r2Xr3);
  auto t2 = rd / (len(r1) * len(r2) * len(r3));
  auto dd = 1. - t2;
  auto d_pow = std::pow(dd, 4);
  return 0.045 * d_pow;
}

std::string ic_core::out_of_plane::info(coords::Representation_3D const & xyz) const
{
  std::ostringstream oss;
  oss << val(xyz) << "||" << index_a_ << "||" << index_b_ << "||" << index_c_ << "||" << index_d_ << "\n";
  return oss.str();
}

std::vector<float_type>
ic_core::trans_x::der_vec(coords::Representation_3D const& xyz) const{
  using cp = coords::Cartesian_Point;

  return ic_util::flatten_c3_vec(der(xyz.size(), [](auto const & s) {
    return cp(1. / static_cast<float_type>(s), 0., 0.);
  }));
}

std::string ic_core::trans_x::info(coords::Representation_3D const & xyz) const
{
  std::ostringstream oss;
  oss << "Trans X: " << val(xyz) << "\n";
  return oss.str();
}

std::vector<float_type>
ic_core::trans_y::der_vec(coords::Representation_3D const& xyz) const{

  using cp = coords::Cartesian_Point;

  return ic_util::flatten_c3_vec(der(xyz.size(), [](auto const & s) {
    return cp(0., 1. / static_cast<float_type>(s), 0.);
  }));
}

std::string ic_core::trans_y::info(coords::Representation_3D const & xyz) const
{
  std::ostringstream oss;
  oss << "Trans Y: " << val(xyz) << "\n";
  return oss.str();
}

std::vector<float_type>
ic_core::trans_z::der_vec(coords::Representation_3D const& xyz) const{
  using cp = coords::Cartesian_Point;

  return ic_util::flatten_c3_vec(der(xyz.size(), [](auto const & s) {
    return cp(0., 0., 1. / static_cast<float_type>(s));
  }));
}

std::string ic_core::trans_z::info(coords::Representation_3D const & xyz) const
{
  std::ostringstream oss;
  oss << "Trans Z: " << val(xyz) << "\n";
  return oss.str();
}

//coords::Representation_3D ic_core::rotation::xyz0;

ic_core::rotation ic_core::build_rotation(coords::Representation_3D const& target,
  std::vector<std::size_t> const& index_vec){
    coords::Representation_3D reference;
  for (auto const & ind : index_vec) {
    reference.emplace_back(target.at(ind - 1));
  }
  return rotation(std::move(reference), index_vec);
}

std::array<float_type, 3u>
ic_core::rotation::rot_val(const coords::Representation_3D& new_xyz) const {
  coords::Representation_3D curr_xyz_;
  curr_xyz_.reserve(indices_.size());
  for (auto const & i : indices_) {
    curr_xyz_.emplace_back(new_xyz.at(i - 1));
  }
  auto result = ic_rotation::exponential_map(curr_xyz_, reference_);
  return result;
}

std::vector<scon::mathmatrix<float_type>>
ic_core::rotation::rot_der(const coords::Representation_3D& new_xyz) const{
  coords::Representation_3D new_xyz_;
  for (auto const & indi : indices_) {
    new_xyz_.emplace_back(new_xyz.at(indi - 1));
  }

  return ic_rotation::exponential_derivs(new_xyz_, reference_);
}

scon::mathmatrix<float_type>
ic_core::rotation::rot_der_mat(std::size_t const & sys_size, const coords::Representation_3D& new_xyz)const {
  using Mat = scon::mathmatrix<float_type>;
  auto const & zero = scon::mathmatrix<float_type>::zero;

  Mat X = zero(sys_size, 3);
  Mat Y = zero(sys_size, 3);
  Mat Z = zero(sys_size, 3);

  auto first_ders = rot_der(new_xyz);
  for(auto i{0u}; i<first_ders.size();++i){
    auto const& ind = indices_.at(i);
    auto const& der = first_ders.at(i);
    X.set_row(ind - 1, der.col(0).t());
    Y.set_row(ind - 1, der.col(1).t());
    Z.set_row(ind - 1, der.col(2).t());
  }
  X *= rad_gyr_;
  Y *= rad_gyr_;
  Z *= rad_gyr_;
  Mat result(sys_size * 3, 3);
  result.set_col(0, X.vectorise_row());
  result.set_col(1, Y.vectorise_row());
  result.set_col(2, Z.vectorise_row());
  return result;
}

coords::float_type
ic_core::rotation::radius_gyration(const coords::Representation_3D& struc) {
  return ic_util::rad_gyr(struc);
}

std::vector<std::unique_ptr<ic_core::internal_coord>> ic_core::system::create_trans_x(
    const std::vector<std::vector<std::size_t>>& index_vec) const {

  std::vector<std::unique_ptr<internal_coord>> result;
  for (auto const & indices : index_vec) {
    result.emplace_back(std::make_unique<trans_x>(indices));
  }
  return result;
}

std::vector<std::unique_ptr<ic_core::internal_coord>> ic_core::system::create_trans_y(
    const std::vector<std::vector<std::size_t>>& index_vec) const {

  std::vector<std::unique_ptr<internal_coord>> result;
  for (auto const & indices : index_vec) {
    result.emplace_back(std::make_unique<trans_y>(indices));
  }
  return result;
}

std::vector<std::unique_ptr<ic_core::internal_coord>> ic_core::system::create_trans_z(
    const std::vector<std::vector<std::size_t>>& index_vec) const {

  std::vector<std::unique_ptr<internal_coord>> result;
  for (auto const & indices : index_vec) {
    result.emplace_back(std::make_unique<trans_z>(indices));
  }
  return result;
}

std::vector<ic_core::rotation> ic_core::system::create_rotations(
    const coords::Representation_3D& rep,
    const std::vector<std::vector<std::size_t>>& index_vec) {
  std::vector<ic_core::rotation> result;
  for(auto const& iv : index_vec){
    result.emplace_back(ic_core::build_rotation(rep, iv));
  }
  return result;
}

std::vector<std::vector<float_type>> ic_core::system::deriv_vec(){
  std::vector<std::vector<float_type>> result;

  for (auto const& pic : primitive_internals) {
    result.emplace_back(pic->der_vec(xyz_));
  }
  //Refactor the next few lines!
  std::vector<std::vector<float_type>> X_rot;
  std::vector<std::vector<float_type>> Y_rot;
  std::vector<std::vector<float_type>> Z_rot;
  for (auto const& i : rotation_vec_) {
    auto temp = i.rot_der_mat(xyz_.size(), xyz_);
    //look for th cols and rows!
    X_rot.emplace_back(temp.col_to_std_vector(0));
    Y_rot.emplace_back(temp.col_to_std_vector(1));
    Z_rot.emplace_back(temp.col_to_std_vector(2));
  }
  auto move_vecs = [&result](auto&& vec){
    for(auto&&v : vec){
      result.emplace_back(std::move(v));
    }
  };
  move_vecs(X_rot);
  move_vecs(Y_rot);
  move_vecs(Z_rot);

  return result;
}

scon::mathmatrix<float_type>& ic_core::system::ic_Bmat(){
  if(!new_B_matrix){
    return B_matrix;
  }
  B_matrix = del_mat.t()*Bmat();
  new_B_matrix = false;
  return B_matrix;
}

scon::mathmatrix<float_type>& ic_core::system::Bmat() {
  if(!new_B_matrix){
    return B_matrix;
  }
  using Mat = scon::mathmatrix<float_type>;

  auto ders = deriv_vec();

  std::size_t n_rows = ders.size(), n_cols = ders.at(0).size();
  B_matrix = Mat(n_rows, n_cols);
  for (std::size_t i{ 0 }; i < n_rows; ++i) {
    B_matrix.set_row(i, Mat::row_from_vec(ders.at(i)));
  }
  new_B_matrix = false;
  return B_matrix;
}

scon::mathmatrix<float_type>& ic_core::system::Gmat(){
  if(!new_G_matrix){
    return G_matrix;
  }
  Bmat();
  G_matrix = B_matrix * B_matrix.t();
  new_G_matrix=false;
  return G_matrix;
}

scon::mathmatrix<float_type>& ic_core::system::ic_Gmat(){
  if(!new_G_matrix){
    return G_matrix;
  }
  ic_Bmat();
  G_matrix = B_matrix * B_matrix.t();
  new_G_matrix=false;
  return G_matrix;
}

scon::mathmatrix<float_type> ic_core::system::guess_hessian() {
  using Mat = scon::mathmatrix<float_type>;
  using scon::cross;
  using scon::dot;
  using scon::len;

  std::vector<float_type> values;
  for (auto const & pic : primitive_internals) {
    values.emplace_back(pic->hessian_guess(xyz_));
  }
  values.insert(values.end(), rotation_vec_.size() * 3, 0.05);
  return Mat::col_from_vec(values).diagmat();
}

/*scon::mathmatrix<float_type>
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
}*/

scon::mathmatrix<float_type>& ic_core::system::delocalize_ic_system() {
  using Mat = scon::mathmatrix<float_type>;
  Gmat();
  Mat eigval, eigvec;
  std::tie(eigval, eigvec) = G_matrix.eigensym(true);
  auto row_index_vec = eigval.sort_idx();
  auto col_index_vec = eigval.find_idx([](float_type const & a) {
    return std::abs(a) > 1e-6;
  });
  del_mat = eigvec.submat(row_index_vec, col_index_vec);
  new_B_matrix = new_G_matrix = true; //B and G got to be calculated fot the new ic_system
  return del_mat;
}

scon::mathmatrix<float_type> ic_core::system::calc() const
{
  std::vector<float_type> primitives;
  primitives.reserve(primitive_internals.size() + rotation_vec_.size()*3);

  for(auto const & pic : primitive_internals){
    primitives.emplace_back(pic->val(xyz_));
  }
  for (auto const & rot : rotation_vec_) {
    auto rv = rot.rot_val(xyz_);
    primitives.emplace_back(rv.at(0));
    primitives.emplace_back(rv.at(1));
    primitives.emplace_back(rv.at(2));
  }
  return scon::mathmatrix<float_type>::row_from_vec(primitives) * del_mat;
}

scon::mathmatrix<float_type> ic_core::system::calculate_internal_grads(scon::mathmatrix<float_type> const& g) {
  auto&& gmat = ic_Gmat();
  return gmat.pinv() * B_matrix * g;
}

scon::mathmatrix<float_type>& ic_core::system::initial_delocalized_hessian(){
  hessian = del_mat.t() * guess_hessian() * del_mat;
  return hessian;
}
