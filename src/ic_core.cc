#include "ic_core.h"
#include "ic_util.h"
#include "energy.h"
#include<iterator>

using coords::float_type;

coords::Representation_3D ic_core::rep3d_bohr_to_ang(coords::Representation_3D const& bohr){
  coords::Representation_3D ang;
  ang.reserve(bohr.size());
  for(auto const& b: bohr){
    ang.emplace_back(b * energy::bohr2ang);
  }
  return ang;
}

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
  if(el_a == period::one && el_b == period::one){
    B_val = -0.244;
  }
  else if ((el_a == period::one && el_b == period::two) ||
           (el_a == period::two && el_b == period::one)) {
    B_val = 0.352;
  }
  else if ((el_a == period::one && el_b == period::three) ||
           (el_a == period::three && el_b == period::one)) {
    B_val = 0.660;
  }
  else if (el_a == period::two && el_b == period::two) {
    B_val = 1.085;
  }
  else if ((el_a == period::two && el_b == period::three) ||
           (el_a == period::three && el_b == period::two)) {
    B_val = 1.522;
  }
  else {
    B_val = 2.068;
  }
  return 1.734 / std::pow(val(xyz) - B_val, 3);
}

std::string ic_core::distance::info(coords::Representation_3D const & xyz) const
{
  std::ostringstream oss;
  oss << "Bond: " << val(xyz) << "||" << index_a_ << "||" << index_b_ << "\n";
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
  using ic_atom::element_period;
  using ic_atom::period;

  auto el_a = element_period(elem_a_);
  auto el_b = element_period(elem_b_);
  auto el_c = element_period(elem_c_);
  if (el_a == period::one ||
      el_b == period::one ||
      el_c == period::one) {
    return 0.160;
  }
  else {
    return 0.250;
  }
}

std::string ic_core::angle::info(coords::Representation_3D const & xyz) const
{
  std::ostringstream oss;
  oss << "Angle: " << val(xyz) * SCON_180PI << "||" << index_a_ << "||" << index_b_ << "||" << index_c_ << "\n";
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
  return std::atan2(x,y);
}
/*Implementation using armas routines instread of CASTs
std::tuple<scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>>
ic_core::dihedral::der(coords::Representation_3D const& xyz) const {
  using scon::cross;
  using scon::dot;
  using scon::len;

  auto const & a_ = xyz.at(index_a_ - 1u);
  auto const & b_ = xyz.at(index_b_ - 1u);
  auto const & c_ = xyz.at(index_c_ - 1u);
  auto const & d_ = xyz.at(index_d_ - 1u);

  Eigen::Vector3d a(a_.x(),a_.y(),a_.z());
  Eigen::Vector3d b(b_.x(),b_.y(),b_.z());
  Eigen::Vector3d c(c_.x(),c_.y(),c_.z());
  Eigen::Vector3d d(d_.x(),d_.y(),d_.z());

  auto u_p = a - b;
  auto w_p = c - b;
  auto v_p = d - c;
  auto u = u_p/u_p.norm();
  auto w = w_p/w_p.norm();
  auto v = v_p/v_p.norm();
  //auto u = ic_util::normalize(u_p);
  //auto w = ic_util::normalize(w_p);
  //auto v = ic_util::normalize(v_p);
  auto dot_uw = u.dot(w);
  auto dot_vw = v.dot(w);
  auto sin2_u = 1. - dot_uw*dot_uw;
  auto sin2_v = 1. - dot_vw*dot_vw;

  auto t1 = u.cross(w) / (u_p.norm() * sin2_u);
  auto t2 = v.cross(w) / (v_p.norm() * sin2_v);
  auto cuwd = u.cross(w) * u.dot(w);
  auto t3 = cuwd / (w_p.norm() * sin2_u);
  //rechanged it: Lee-Ping pointed out that his numeric evaluated derivatives match with analytics without this sign ---> changed v's sign according to J. Chem. Ohys Vol 117 No. 20, 22 2002 p. 9160-9174
  auto cvwd = v.cross(w) * v.dot(w);//auto cvwd = cross(v, w) * dot(-v, w);
  auto t4 = cvwd / (w_p.norm() * sin2_v);
  auto make_vec = [](auto vec){
    return scon::c3<float_type>(vec(0), vec(1), vec(2));
  };
*/
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
  //rechanged it: Lee-Ping pointed out that his numeric evaluated derivatives match with analytics without this sign ---> changed v's sign according to J. Chem. Ohys Vol 117 No. 20, 22 2002 p. 9160-9174
  auto cvwd = cross(v, w) * dot(v, w);//auto cvwd = cross(v, w) * dot(-v, w);
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
  oss << "Dihedral: " << val(xyz) * SCON_180PI << "||" << index_a_ << "||" << index_b_ << "||" << index_c_ << "||" << index_d_ << "\n";
  return oss.str();
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
  oss << "Out of plane: " << val(xyz) * SCON_180PI << "||" << index_a_ << "||" << index_b_ << "||" << index_c_ << "||" << index_d_ << "\n";
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
  /*auto ref = reference_;
  ref.at(0) = coords::Cartesian_Point(-5.380,   2.994,  -0.001);
  ref.at(1) = coords::Cartesian_Point(-3.982,   2.994,   0.013);
  ref.at(2) = coords::Cartesian_Point(-5.771,   3.923,   0.468);
  ref.at(3) = coords::Cartesian_Point(-5.771,   2.093,   0.520);
  ref.at(4) = coords::Cartesian_Point(-5.729,   2.964,  -1.054);
  ref.at(5) = coords::Cartesian_Point(-3.718,   3.021,   0.969);
  curr_xyz_.at(0) = coords::Cartesian_Point(-6.40029,       -0.54123,       -0.00037);
  curr_xyz_.at(1) = coords::Cartesian_Point(-5.00210,       -0.52808,        0.00408);
  curr_xyz_.at(2) = coords::Cartesian_Point(-6.79116,        0.22071,       -0.70912);
  curr_xyz_.at(3) = coords::Cartesian_Point(-6.79116,       -0.36516,        1.02525);
  curr_xyz_.at(4) = coords::Cartesian_Point(-6.74902,       -1.53902,       -0.33742);
  curr_xyz_.at(5) = coords::Cartesian_Point(-4.73863,        0.37832,        0.31026);
  Test Ergebnis: 2.7526415805955393, -1.0064325510255895e-05, -0.00046769875255223972
  */

  auto result = ic_rotation::exponential_map(reference_, curr_xyz_);
  for(auto & r : result){
    r*=rad_gyr_;
  }
  return result;
}

std::vector<scon::mathmatrix<float_type>>
ic_core::rotation::rot_der(const coords::Representation_3D& new_xyz) const{
  coords::Representation_3D new_xyz_;
  for (auto const & indi : indices_) {
    new_xyz_.emplace_back(new_xyz.at(indi - 1));
  }

  return ic_rotation::exponential_derivs(reference_, new_xyz_);
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

  Mat eigval, eigvec;
  std::tie(eigval, eigvec) = Gmat().eigensym(false);

  auto row_index_vec = eigval.sort_idx();
  auto col_index_vec = eigval.find_idx([](float_type const & a) {
    return std::abs(a) > 1e-6;
  });

  del_mat = eigvec.submat(row_index_vec, col_index_vec);
  //The Del Matrix of Lee-Pings Code. The results are the same
  /*del_mat <<    0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,   -0.4214840792,   -0.0000000000,    0.0000475867,    0.3389630690,    0.0497763428,    0.0505143352,   -0.0045016042,    0.4745152490,    0.0029898210,    0.0000860452,    0.6830371535,   -0.1033859621,    0.0018659334,
   0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0511919632,   -0.0000000000,   -0.1817409988,   -0.0937377047,   -0.4054358696,   -0.2276238256,   -0.4250044002,    0.4150877982,   -0.0599637319,    0.2592102085,   -0.1872555758,    0.2878904914,    0.4446506462,
   0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0511801633,    0.0000000000,    0.1816405260,   -0.0939221680,   -0.4056226532,   -0.1485893819,    0.4585773957,    0.4149381079,   -0.0598878041,   -0.2594215153,   -0.1874384803,    0.2688736211,   -0.4562655526,
  -0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,    0.0255435857,    0.0000000000,    0.0000739972,    0.1927626986,   -0.4060928570,    0.5924710557,   -0.0531501073,    0.1702081287,    0.2588064044,    0.0000449700,   -0.2914197967,   -0.5122878637,    0.0107633769,
  -0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,   -0.0762221428,   -0.0000000000,    0.0000853751,    0.3173917019,   -0.1124076281,    0.0614945216,   -0.0055314262,   -0.1167473148,   -0.9077997692,   -0.0001900132,   -0.1374663761,   -0.1427175964,    0.0030541927,
  -0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0000000000,   -0.3161070027,   -0.0000000000,   -0.6257921778,   -0.1619263503,    0.2830974701,    0.0084894540,    0.0329583501,    0.1553676849,   -0.0282203458,    0.0672036210,   -0.2601976595,   -0.1147797145,   -0.2556116589,
  -0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0000000000,   -0.3163425769,    0.0000000000,    0.6257035744,   -0.1624155935,    0.2829772915,    0.0026805744,   -0.0338212678,    0.1553205560,   -0.0282199859,   -0.0672701296,   -0.2598944293,   -0.1038483376,    0.2603526595,
  -0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,    0.2703433008,    0.0000000000,    0.0004411247,    0.7180447068,    0.1590767703,   -0.0577875667,    0.0053054355,    0.0576593076,    0.1778204328,   -0.0000339527,   -0.1910397760,    0.3077148269,   -0.0062963156,
  -0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,    0.1860771607,    0.0000000000,    0.0000081932,    0.0953737406,   -0.3701851429,   -0.4811164471,    0.0429740721,   -0.1892665578,    0.0806574291,    0.0000911068,    0.2066818785,   -0.3618163482,    0.0075845874,
  -0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,   -0.0000000000,    0.0982229742,    0.0000000000,    0.1947756921,   -0.2387637810,   -0.1824119767,    0.2318623574,   -0.4766694183,   -0.0926091783,   -0.1010983241,   -0.1506250706,    0.2591632576,    0.1420216091,   -0.2848402869,
  -0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,    0.0982582793,    0.0000000000,   -0.1950346458,   -0.2385108856,   -0.1825049779,    0.3128822433,    0.4279112735,   -0.0925088808,   -0.1011087309,    0.1506782364,    0.2592074278,    0.1538520434,    0.2785129255,
  -0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,   -0.6817439310,   -0.0000000000,   -0.0001217631,    0.1829037993,   -0.3137345402,    0.0169918950,   -0.0015724351,   -0.5318067826,    0.1988472897,   -0.0000118470,   -0.0800120395,    0.2725025813,   -0.0057058833,
   0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,   -0.0971335394,   -0.0000000000,   -0.1007735682,   -0.0498768406,   -0.0418351301,   -0.2865836670,    0.1800170492,   -0.0141429826,    0.0319073894,   -0.4675111069,   -0.0457237789,   -0.2950757887,    0.3861589291,
   0.0000000000,   -0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,    0.0971093153,   -0.0000000000,   -0.1009314627,    0.0498047071,    0.0417480374,    0.3140666391,    0.1262856509,    0.0141482470,   -0.0317334066,   -0.4675588027,    0.0461703735,    0.3108022171,    0.3732652239,
  -0.0000000000,    0.0000000000,    0.0000000000,    0.0000000000,    0.0000000000,    0.0000933452,   -0.0000000000,   -0.2335323014,    0.0002089553,    0.0001804930,   -0.0339034093,   -0.3787066397,    0.0000143569,    0.0001494796,   -0.6114729603,   -0.0000140309,   -0.0024262540,   -0.1165300369,
  -0.1963503436,    0.8533055572,   -0.4830281241,   -0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,
  -0.8523809969,   -0.3920312477,   -0.3460608863,    0.0000000000,   -0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0000000000,    0.0000000000,
  -0.4846577956,    0.3437748200,    0.8043169116,   -0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,
  -0.0000000000,    0.0000000000,   -0.0000000000,   -0.0008611594,   -0.4314072789,    0.0000000000,   -0.9021568700,    0.0000000000,    0.0000000000,    0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,   -0.0000000000,
   0.0000000000,   -0.0000000000,   -0.0000000000,   -0.9995650157,    0.0269664556,    0.0000000000,   -0.0119410944,    0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,   -0.0000000000,    0.0000000000,
  -0.0000000000,   -0.0000000000,   -0.0000000000,    0.0294794482,    0.9017541627,   -0.0000000000,   -0.4312428459,    0.0000000000,   -0.0000000000,    0.0000000000,    0.0000000000,    0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,   -0.0000000000,    0.0000000000,   -0.0000000000;*/
  new_B_matrix = new_G_matrix = true; //B and G got to be calculated for the new ic_system
  return del_mat;
}

scon::mathmatrix<float_type> ic_core::system::calculate_internal_grads(scon::mathmatrix<float_type> const& g) {
  return ic_Gmat().pinv() * ic_Bmat() * g;
}

scon::mathmatrix<float_type>& ic_core::system::initial_delocalized_hessian(){
  hessian = del_mat.t() * guess_hessian() * del_mat;
  return hessian;
}

void ic_core::system::optimize(coords::DL_Coordinates & coords){
  coords.set_xyz(ic_core::rep3d_bohr_to_ang(xyz_));
  initial_delocalized_hessian();
  coords::output::formats::tinker output(coords);
  output.to_stream(std::cout);

  auto old_E = coords.g();
  auto old_Q = calc(xyz_);
  auto gxyz = scon::mathmatrix<coords::float_type>::col_from_vec(ic_util::flatten_c3_vec(
    ic_core::grads_to_bohr(coords.g_xyz())
  ));
  auto old_xyz = xyz_;
  auto old_gq = calculate_internal_grads(gxyz);
  for(auto i = 0; i< 10; ++i){

    auto dq_step = get_internal_step(old_gq);
    apply_internal_change(dq_step);

    coords.set_xyz(ic_core::rep3d_bohr_to_ang(xyz_));
    auto new_E = coords.g();
    auto gxyz = scon::mathmatrix<coords::float_type>::col_from_vec(ic_util::flatten_c3_vec(
      ic_core::grads_to_bohr(coords.g_xyz())
    ));

    auto new_gq = calculate_internal_grads(gxyz);
    auto new_Q = calc(xyz_);

    auto d_gq = old_gq - new_gq;
    auto dq = old_Q - new_Q;
    //std::cout << "Delta Grads: \n" << d_gq << "\n\n";
    //std::cout << "Delta Internals: \n" << dq << "\n\n";
    auto term1 = (d_gq*d_gq.t())/(d_gq.t()*dq)(0,0);
    auto term2 = ((hessian*dq)*(dq.t()*hessian))/(dq.t()*hessian*dq)(0,0);
    hessian += term1 - term2;
    //std::cout << "New Hessian:\n" << hessian << "\n\n";
    old_gq = std::move(new_gq);
    old_Q = std::move(new_Q);
    std::cout << "Energy now: " << new_E << " Energy diff: " << new_E-old_E <<"\n";
    old_E = new_E;
    std::cout << "RMSD Cartesian: " << ic_util::Rep3D_to_Mat(old_xyz - xyz_).rmsd() << "\n";
    old_xyz = xyz_;

    output.to_stream(std::cout);
  }
}
