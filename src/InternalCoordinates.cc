#include "InternalCoordinates.h"

#include "ic_util.h"
#include "ic_rotation.h"
#include "helperfunctions.h"

namespace InternalCoordinates {

  coords::float_type BondDistance::val(coords::Representation_3D const& cartesians) const {
    auto const& a = cartesians.at(index_a_);
    auto const& b = cartesians.at(index_b_);

    return scon::len(a - b);
  }

  std::pair<coords::r3, coords::r3> BondDistance::der(coords::Representation_3D const& cartesians) const {

    auto const& a = cartesians.at(index_a_);
    auto const& b = cartesians.at(index_b_);

    auto bond = ic_util::normalize(a - b);
    return std::make_pair(bond, -bond);
  }

  std::vector<coords::float_type>
    BondDistance::der_vec(coords::Representation_3D const& cartesians) const {
    using scon::c3;

    auto firstder = der(cartesians);
    std::vector<coords::r3> der_vec(cartesians.size(), coords::r3(0.0, 0.0, 0.0));
    der_vec.at(index_a_) = firstder.first;
    der_vec.at(index_b_) = firstder.second;
    return ic_util::flatten_c3_vec(der_vec);
  }

  //enums hold an integral thus I pass them by value
  bool BondDistance::bothElementsInPeriodOne(ic_atom::period const atomA, ic_atom::period const atomB) const {
    return atomA == ic_atom::period::one && atomB == ic_atom::period::one;
  }

  bool BondDistance::oneElementInPeriodOneTheOtherInPeriodTwo(ic_atom::period const atomA, ic_atom::period const atomB) const {
    return (atomA == ic_atom::period::one && atomB == ic_atom::period::two) 
		|| (atomA == ic_atom::period::two && atomB == ic_atom::period::one);
  }

  bool BondDistance::oneElementInPeriodOneTheOtherInPeriodThree(ic_atom::period const atomA, ic_atom::period const atomB) const {
    return (atomA == ic_atom::period::one && atomB == ic_atom::period::three) 
		|| (atomA == ic_atom::period::three && atomB == ic_atom::period::one);
  }

  bool BondDistance::bothElementsInPeriodTwo(ic_atom::period const atomA, ic_atom::period const atomB) const {
    return atomA == ic_atom::period::two && atomB == ic_atom::period::two;
  }

  bool BondDistance::oneElementInPeriodTwoTheOtherInPeriodThree(ic_atom::period const atomA, ic_atom::period const atomB) const {
    return (atomA == ic_atom::period::two && atomB == ic_atom::period::three)
	   	|| (atomA == ic_atom::period::three && atomB == ic_atom::period::two);
  }

  coords::float_type BondDistance::hessian_guess(coords::Representation_3D const& cartesians) const {
    using ic_atom::element_period;
    using ic_atom::period;

    auto el_a = element_period(elem_a_);
    auto el_b = element_period(elem_b_);

    auto B_val{ 0.0 };
    if (bothElementsInPeriodOne(el_a, el_b)) {
      B_val = -0.244;
    }
    else if (oneElementInPeriodOneTheOtherInPeriodTwo(el_a, el_b)) {
      B_val = 0.352;
    }
    else if (oneElementInPeriodOneTheOtherInPeriodThree(el_a, el_b)) {
      B_val = 0.660;
    }
    else if (bothElementsInPeriodTwo(el_a, el_b)) {
      B_val = 1.085;
    }
    else if (oneElementInPeriodTwoTheOtherInPeriodThree(el_a, el_b)) {
      B_val = 1.522;
    }
    else {
      B_val = 2.068;
    }
    return 1.734 / std::pow(val(cartesians) - B_val, 3);
  }

  std::string BondDistance::info(coords::Representation_3D const & cartesians) const {
    std::ostringstream oss;
    auto l_bohr = val(cartesians);
    oss << "Bond: " << l_bohr << " (a. u.) " << energy::bohr2ang*l_bohr << " (Angstrom) || " << index_a_ + 1u << " || " << index_b_ + 1u << " || " << "Constrained: " << std::boolalpha << is_constrained();
    return oss.str();
  }

  bool BondDistance::operator==(BondDistance const & other) const {
	  return index_a_ == other.index_a_ && index_b_ == other.index_b_ && elem_a_ == other.elem_a_ && elem_b_ == other.elem_b_;
  }

  /*bool operator==(BondDistance const & lhs, BondDistance const & rhs) {
    return lhs.operator==(rhs);
  }*/

  coords::float_type BondAngle::val(coords::Representation_3D const& cartesians) const {
    using scon::cross;
    using scon::dot;
    using scon::len;

    auto const& a = cartesians.at(index_a_);
    auto const& b = cartesians.at(index_b_);
    auto const& c = cartesians.at(index_c_);

    auto u = a - b;
    auto v = c - b;
    auto uXv = cross(u, v);
    auto uDv = dot(u, v);
    auto l = len(uXv);
    return std::atan2(l, uDv);
  }

  std::tuple<coords::r3, coords::r3, coords::r3>
    BondAngle::der(coords::Representation_3D const& cartesians) const {
    using coords::Cartesian_Point;
    using scon::cross;
    using scon::dot;
    using scon::len;

    auto const& a = cartesians.at(index_a_);
    auto const& b = cartesians.at(index_b_);
    auto const& c = cartesians.at(index_c_);

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
      }
      else {
        w_p = cross(u, cp1);
      }
    }
    else {
      w_p = cross(u, v);
    }
    auto w = ic_util::normalize(w_p);
    auto ad0 = cross(u, w) / lu;
    auto ad1 = cross(w, v) / lv;
    //                      a     b           c
    //                      m     o           n
    return std::make_tuple(ad0, -ad0 - ad1, ad1);
  }

  std::vector<coords::float_type>
    BondAngle::der_vec(coords::Representation_3D const& cartesians) const {
    using scon::c3;

    auto firstder = BondAngle::der(cartesians);
    std::vector<coords::r3> der_vec(cartesians.size(), coords::r3(0.0, 0.0, 0.0));
    der_vec.at(index_a_) = std::get<0>(firstder);
    der_vec.at(index_b_) = std::get<1>(firstder);
    der_vec.at(index_c_) = std::get<2>(firstder);
    return ic_util::flatten_c3_vec(der_vec);
  }

  coords::float_type BondAngle::hessian_guess(coords::Representation_3D const& /*cartesians*/) const {
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

  std::string BondAngle::info(coords::Representation_3D const & cartesians) const {
    std::ostringstream oss;
    oss << "Angle: " << val(cartesians) * SCON_180PI << " || " << index_a_ + 1 << " || " << index_b_ + 1 << " || " << index_c_ + 1 << " || " << "Constrained: " << std::boolalpha << is_constrained();
    return oss.str();
  }

  bool BondAngle::operator==(BondAngle const & other) const {
    return index_a_ == other.index_a_ && index_b_ == other.index_b_ && index_c_ == other.index_c_ 
      && elem_a_ == other.elem_a_ && elem_b_ == other.elem_b_ && elem_c_ == other.elem_c_;
  }

  coords::float_type DihedralAngle::val(coords::Representation_3D const& cartesians) const {
    using scon::cross;
    using scon::dot;
    using scon::len;


    //Eigentlich sollte das hier richtig sein:
    auto const & a = cartesians.at(index_a_);
    auto const & b = cartesians.at(index_b_);
    auto const & c = cartesians.at(index_c_);
    auto const & d = cartesians.at(index_d_);

    auto b1 = b - a;
    b1 /= len(b1);
    auto b2 = c - b;
    b2 /= len(b2);
    auto b3 = d - c;
    b3 /= len(b3);
    auto n1 = cross(b1, b2);
    auto n2 = cross(b2, b3);
    return std::atan2(dot(b1, n2), dot(n1, n2));

  }
 
  std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
    DihedralAngle::der(coords::Representation_3D const& cartesians) const {
    using scon::cross;
    using scon::dot;
    using scon::len;

    auto const & a = cartesians.at(index_a_);
    auto const & b = cartesians.at(index_b_);
    auto const & c = cartesians.at(index_c_);
    auto const & d = cartesians.at(index_d_);

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

  std::vector<coords::float_type>
    DihedralAngle::der_vec(coords::Representation_3D const& cartesians) const {
    using scon::c3;

    auto firstder = der(cartesians);
    coords::r3 temp(0.0, 0.0, 0.0);
    std::vector<coords::r3> der_vec(cartesians.size(), temp);
    der_vec.at(index_a_) = std::get<0>(firstder);
    der_vec.at(index_b_) = std::get<1>(firstder);
    der_vec.at(index_c_) = std::get<2>(firstder);
    der_vec.at(index_d_) = std::get<3>(firstder);
    return ic_util::flatten_c3_vec(der_vec);
  }

  coords::float_type DihedralAngle::hessian_guess(coords::Representation_3D const & /*cartesians*/) const{
    return 0.023;
  }

  std::string DihedralAngle::info(coords::Representation_3D const & cartesians) const{
    std::ostringstream oss;
    oss << "Dihedral: " << val(cartesians) * SCON_180PI << " || " << index_a_ + 1u << " || " << index_b_ + 1u << " || " << index_c_ + 1u << " || " << index_d_ + 1u << " || " << "Constrained: " << std::boolalpha << is_constrained();
    return oss.str();
  }

  bool DihedralAngle::operator==(DihedralAngle const & other) const{
    return index_a_ == other.index_a_ && index_b_ == other.index_b_ && index_c_ == other.index_c_ && index_d_ == other.index_d_;
  }

  coords::float_type OutOfPlane::hessian_guess(coords::Representation_3D const & cartesians) const
  {
    auto const& a = cartesians.at(index_a_);
    auto const& b = cartesians.at(index_b_);
    auto const& c = cartesians.at(index_c_);
    auto const& d = cartesians.at(index_d_);
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

  std::string OutOfPlane::info(coords::Representation_3D const & cartesians) const
  {
    std::ostringstream oss;
    oss << "Out of plane: " << val(cartesians) * SCON_180PI << "||" << index_a_ + 1u << "||" << index_b_ + 1u << "||" << index_c_ + 1u << "||" << index_d_ + 1u << " || " << "Constrained: " << std::boolalpha << is_constrained();
    return oss.str();
  }

  coords::float_type Translations::val(coords::Representation_3D const &cartesians) const {
    auto coord_sum{0.0};
    for (auto &i : indices_) {
      coord_sum += coordinate_value(cartesians.at(i));
    }
    return coord_sum / indices_.size();
  }

  std::string Translations::info(coords::Representation_3D const & cartesians) const
  {
    std::ostringstream oss;
    oss << "Trans " << coordinate_letter() << ": " << val(cartesians) << " | Constrained: " << std::boolalpha << is_constrained();
    return oss.str();
  }

  std::vector<coords::float_type>
  Translations::der_vec(coords::Representation_3D const& cartesians) const {
    return ic_util::flatten_c3_vec(der(cartesians.size(), [this](std::size_t s) {
      return size_reciprocal(s);
    }));
  }

  void RotatorObserver::setNewRotator(std::shared_ptr<AbstractRotatorListener> const rotator) { this->rotator = rotator; }

  void RotatorObserver::update(){
    notify();
  }

  void RotatorObserver::notify(){
    rotator->setAllFlag();
  }

  std::array<coords::float_type, 3u> const&
    Rotator::valueOfInternalCoordinate(const coords::Representation_3D& new_xyz) {
    //TODO This needs to be activated again
    //if (!updateStoredValues) {
    //  return storedValuesForRotations;
    //}
    updateStoredValues = false;
    coords::Representation_3D curr_xyz_;
    curr_xyz_.reserve(indices_.size());
    for (auto const & i : indices_) {
      curr_xyz_.emplace_back(new_xyz.at(i));
    }

    storedValuesForRotations = ic_rotation::exponential_map(reference_, curr_xyz_);
    for (auto & r : storedValuesForRotations) {
      r *= rad_gyr_;
    }
    return storedValuesForRotations;
  }

  std::vector<scon::mathmatrix<coords::float_type>>
    InternalCoordinates::Rotator::rot_der(coords::Representation_3D const& new_xyz) const {
    
    coords::Representation_3D new_xyz_;
    for (auto const & indi : indices_) {
      new_xyz_.emplace_back(new_xyz.at(indi));
    }

    return ic_rotation::exponential_derivs(reference_, new_xyz_);
  }

  scon::mathmatrix<coords::float_type> const&
    InternalCoordinates::Rotator::rot_der_mat(coords::Representation_3D const& new_xyz) {
    using Mat = scon::mathmatrix<coords::float_type>;
    //TODO This needs to be activated again
    //if (!updateStoredDerivatives) {
    //  return storedDerivativesForRotations;
    //}
    updateStoredDerivatives = false;

    auto const & zero = scon::mathmatrix<coords::float_type>::zero;

    auto first_ders = rot_der(new_xyz);

    Mat X = zero(new_xyz.size(), 3);
    Mat Y = zero(new_xyz.size(), 3);
    Mat Z = zero(new_xyz.size(), 3);

    for (auto i{ 0u }; i<first_ders.size(); ++i) {
      auto const& ind = indices_.at(i);
      auto const& der = first_ders.at(i);
      X.set_row(ind, der.col(0).t());
      Y.set_row(ind, der.col(1).t());
      Z.set_row(ind, der.col(2).t());
    }

    X *= rad_gyr_;
    Y *= rad_gyr_;
    Z *= rad_gyr_;

    //TODO initialize storedDerivativesForRotations in ctor
    storedDerivativesForRotations = Mat(new_xyz.size() * 3, 3);
    storedDerivativesForRotations.set_col(0, X.vectorise_row());
    storedDerivativesForRotations.set_col(1, Y.vectorise_row());
    storedDerivativesForRotations.set_col(2, Z.vectorise_row());
    return storedDerivativesForRotations;
  }

  Rotations Rotator::makeRotations() {
    return Rotations{
      shared_from_this(),
      std::make_unique<RotationA>(shared_from_this()),
      std::make_unique<RotationB>(shared_from_this()),
      std::make_unique<RotationC>(shared_from_this())
    };
  }

  bool Rotator::operator==(Rotator const & other) const {
    if (indices_.size() != other.indices_.size()) return false;
    for (auto i = 0u; i < indices_.size(); ++i) {
      if (indices_.at(i) != other.indices_.at(i)) return false;
    }
    return true;
  }

  coords::float_type
    InternalCoordinates::Rotator::radiusOfGyration(const coords::Representation_3D& struc) {
    return ic_util::rad_gyr(struc);
  }

  void Rotator::registerCartesians(InternalCoordinates::CartesiansForInternalCoordinates & cartesianCoordinates){
    auto observer = std::make_shared<InternalCoordinates::RotatorObserver>();
    observer->setNewRotator(shared_from_this());
    cartesianCoordinates.registerObserver(observer);
  }

  bool Translations::operator==(Translations const & other) const{
    if(indices_.size() != other.indices_.size()) return false;
    for (auto i = 0u; i < indices_.size(); ++i) {
      if (indices_.at(i) != other.indices_.at(i)) return false;
    }
    return true;
  }
  
  
  
}
