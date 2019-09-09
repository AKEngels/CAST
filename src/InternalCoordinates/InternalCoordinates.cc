#include "InternalCoordinates.h"

#include "../Scon/scon_mathmatrix.h"

#include "InternalCoordinateUtilities.h"
#include "BestFitRotation.h"
#include "../helperfunctions.h"
#include "Init/ConstraintManager.h"
#include "BondGraph/BondGraph.h"

namespace InternalCoordinates {

	CartesiansForInternalCoordinates::CartesiansForInternalCoordinates(CartesiansForInternalCoordinates && cartesians)
		: observerList(std::move(cartesians.observerList)), coordinates(std::move(cartesians.coordinates)) {}

	CartesiansForInternalCoordinates::CartesiansForInternalCoordinates(CartesiansForInternalCoordinates const& cartesians)
		: observerList(cartesians.observerList), coordinates(cartesians.coordinates) {}

	CartesiansForInternalCoordinates::CartesiansForInternalCoordinates(coords::Representation_3D && cartesians)
		: observerList(), coordinates(std::move(cartesians)) {}

	CartesiansForInternalCoordinates::CartesiansForInternalCoordinates(coords::Representation_3D const& cartesians)
		: observerList(), coordinates(cartesians) {}

	CartesiansForInternalCoordinates & CartesiansForInternalCoordinates::operator=(CartesiansForInternalCoordinates const& cartesians) {
		observerList = cartesians.observerList;
		setCartesianCoordnatesIntern(cartesians.coordinates);
		return *this;
	}

	CartesiansForInternalCoordinates & CartesiansForInternalCoordinates::operator=(CartesiansForInternalCoordinates && cartesians) {
		observerList = std::move(cartesians.observerList);
		setCartesianCoordnatesIntern(std::move(cartesians.coordinates));
		return *this;
	}

  std::pair<coords::float_type, coords::float_type> CartesiansForInternalCoordinates::displacementRmsValAndMaxTwoStructures(coords::Representation_3D const & other) const {
	return ic_rotation::displacementRmsValAndMaxTwoStructures<scon::mathmatrix>(coordinates, other);
  }

  std::pair<coords::float_type, coords::float_type> CartesiansForInternalCoordinates::displacementRmsValAndMaxTwoStructures(CartesiansForInternalCoordinates const& other) const {
	  return displacementRmsValAndMaxTwoStructures(other.coordinates);
  }

	coords::Cartesian_Point const& CartesiansForInternalCoordinates::at(std::size_t const i) const { return coordinates.at(i); }
	coords::Cartesian_Point & CartesiansForInternalCoordinates::at(std::size_t const i) { return coordinates.at(i); }

	coords::float_type CartesiansForInternalCoordinates::getInternalValue(InternalCoordinate const& in) const { return in.val(coordinates); }

	std::vector<coords::float_type> CartesiansForInternalCoordinates::getInternalDerivativeVector(InternalCoordinate const& in) const { return in.der_vec(coordinates); }

	coords::float_type CartesiansForInternalCoordinates::getInternalHessianGuess(InternalCoordinate const& in) const { return in.hessian_guess(coordinates); }

	coords::float_type CartesiansForInternalCoordinates::getInternalDifference(CartesiansForInternalCoordinates const& other, InternalCoordinate const& in) const {
		return in.difference(coordinates, other.coordinates);
	}

	coords::Representation_3D CartesiansForInternalCoordinates::toAngstrom() const { return ic_util::rep3d_bohr_to_ang(coordinates); }

	void CartesiansForInternalCoordinates::registerObserver(std::shared_ptr<RotatorObserver> const observer) {
		observerList.emplace_back(observer);
	}

	void CartesiansForInternalCoordinates::setCartesianCoordnates(coords::Representation_3D const& newCartesianCoordinates) {
		setCartesianCoordnatesIntern(newCartesianCoordinates);
	}
	void CartesiansForInternalCoordinates::setCartesianCoordnates(coords::Representation_3D&& newCartesianCoordinates) {
		setCartesianCoordnatesIntern(std::move(newCartesianCoordinates));
	}

	void CartesiansForInternalCoordinates::reset() { notify(); }

	void CartesiansForInternalCoordinates::notify() {
		for (auto const& observer : observerList) {
			observer->update();
		}
	}

	coords::Representation_3D operator+(CartesiansForInternalCoordinates const& lhs, coords::Representation_3D const& rhs) {
		return lhs.coordinates + rhs;
	}

  coords::float_type BondDistance::val(coords::Representation_3D const& cartesians) const {
    auto const& a = cartesians.at(index_a_);
    auto const& b = cartesians.at(index_b_);

    return scon::len(a - b);
  }

  coords::float_type BondDistance::difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const {
	  return val(newCoordinates) - val(oldCoordinates);
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
  bool BondDistance::bothElementsInPeriodOne(period const atomA, period const atomB) const {
    return atomA == period::one && atomB == period::one;
  }

  bool BondDistance::oneElementInPeriodOneTheOtherInPeriodTwo(period const atomA, period const atomB) const {
    return (atomA == period::one && atomB == period::two) 
		|| (atomA == period::two && atomB == period::one);
  }

  bool BondDistance::oneElementInPeriodOneTheOtherInPeriodThree(period const atomA, period const atomB) const {
    return (atomA == period::one && atomB == period::three) 
		|| (atomA == period::three && atomB == period::one);
  }

  bool BondDistance::bothElementsInPeriodTwo(period const atomA, period const atomB) const {
    return atomA == period::two && atomB == period::two;
  }

  bool BondDistance::oneElementInPeriodTwoTheOtherInPeriodThree(period const atomA, period const atomB) const {
    return (atomA == period::two && atomB == period::three)
	   	|| (atomA == period::three && atomB == period::two);
  }

  coords::float_type BondDistance::hessian_guess(coords::Representation_3D const& cartesians) const {
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

	bool BondDistance::hasIndices(std::vector<std::size_t> const& indices) const {
		return indices.size() == 2u && ((index_a_ == indices[0u] && index_b_ == indices[1u]) || (index_a_ == indices[1u] && index_b_ == indices[0u]));
	}

	std::vector<std::size_t> BondDistance::getIndices() const {
		return { index_a_, index_b_ };
	}

	void BondDistance::makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) { 
		auto constraint = manager->checkForBonds(*this);
		if (constraint) {
			if (constraint->isFrozen()) makeConstrained();
			else releaseConstraint();
		}
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

  coords::float_type BondAngle::difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const {
	  return val(newCoordinates) - val(oldCoordinates);
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

	bool BondAngle::hasIndices(std::vector<std::size_t> const& indices) const {
		return indices.size() == 3u && ((index_a_ == indices[0u] && index_b_ == indices[1u] && index_c_ == indices[2u]) || (index_a_ == indices[2u] && index_b_ == indices[1u] && index_c_ == indices[0u]));
	}

	std::vector<std::size_t> BondAngle::getIndices() const{
		return { index_a_, index_b_, index_c_ };
	}

	void BondAngle::makeConstrained(std::shared_ptr<AbstractConstraintManager> manager){
		auto constraint = manager->checkForAngles(*this);
		if (constraint) {
			if (constraint->isFrozen()) makeConstrained();
			else releaseConstraint();
		}
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

  coords::float_type DihedralAngle::difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const {
	  auto diff = val(newCoordinates) - val(oldCoordinates);
	  if (std::fabs(diff) > SCON_PI) {
		  if (diff < 0.0) {
			  diff += 2. * SCON_PI;
		  }
		  else {
			  diff -= 2. * SCON_PI;
		  }
	  }

	  return diff;
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

	bool DihedralAngle::hasIndices(std::vector<std::size_t> const& indices) const {
		return indices.size() == 3u && ((index_a_ == indices[0u] && index_b_ == indices[1u] && index_c_ == indices[2u] && index_d_ == indices[3u]) || (index_a_ == indices[3u] && index_b_ == indices[2u] && index_c_ == indices[1u] && index_d_ == indices[0u]));
	}

	std::vector<std::size_t> DihedralAngle::getIndices() const {
		return { index_a_, index_b_, index_c_, index_d_ };
	}

	void DihedralAngle::makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) {
		auto constraint = manager->checkForDihedrals(*this);
		if (constraint) {
			if (constraint->isFrozen()) makeConstrained();
			else releaseConstraint();
		}

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
      coord_sum += coord_func(cartesians.at(i));
    }
    return coord_sum / indices_.size();
  }

  coords::float_type Translations::difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const {
	  return val(newCoordinates) - val(oldCoordinates);
  }

  std::string Translations::info(coords::Representation_3D const & cartesians) const
  {
    std::ostringstream oss;
    oss << "Trans " << coordinate_letter() << ": " << val(cartesians) << " | Constrained: " << std::boolalpha << is_constrained();
    return oss.str();
  }

  std::vector<coords::float_type>
  Translations::der_vec(coords::Representation_3D const& cartesians) const {
    coords::Representation_3D result(cartesians.size(), coords::Cartesian_Point (0., 0., 0.));

    for (auto const& i: indices_){
      result.at(i) = size_reciprocal();
    }

    return ic_util::flatten_c3_vec(result);
  }

	bool isSameSet(std::vector<std::size_t> lhs, std::vector<std::size_t> rhs) {
		if (lhs.size() != rhs.size()) return false;
		std::sort(lhs.begin(), lhs.end());
		std::sort(rhs.begin(), rhs.end());
		return std::includes(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
	}

	bool Translations::hasIndices(std::vector<std::size_t> const& indices) const {
		return indices.size() == indices_.size() && isSameSet(indices, indices_);
	}

	std::vector<std::size_t> Translations::getIndices() const {
		return indices_;
	}

	void Translations::makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) {
		auto constraint = manager->checkForTranslation(*this);
		if (constraint) {
			if (constraint->isFrozen()) makeConstrained();
			else releaseConstraint();
		}
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

	storedValuesForRotations = calculateValueOfInternalCoordinate(new_xyz);
    return storedValuesForRotations;
  }

  std::array<coords::float_type, 3u>
	  Rotator::calculateValueOfInternalCoordinate(coords::Representation_3D const& newXyz) const {
	  coords::Representation_3D curr_xyz_;
	  curr_xyz_.reserve(indices_.size());
	  for (auto const & i : indices_) {
		  curr_xyz_.emplace_back(newXyz.at(i));
	  }

	  auto ret = ic_rotation::exponential_map<scon::mathmatrix>(reference_, curr_xyz_);
	  for (auto & r : ret) {
		  r *= rad_gyr_;
	  }

	  return ret;
  }

  std::vector<scon::mathmatrix<coords::float_type>>
    InternalCoordinates::Rotator::rot_der(coords::Representation_3D const& new_xyz) const {
    
    coords::Representation_3D new_xyz_;
    for (auto const & indi : indices_) {
      new_xyz_.emplace_back(new_xyz.at(indi));
    }

    return ic_rotation::exponential_derivs<scon::mathmatrix>(reference_, new_xyz_);
  }

  scon::mathmatrix<coords::float_type> const&
    InternalCoordinates::Rotator::rot_der_mat(coords::Representation_3D const& new_xyz) {
    using Mat = scon::mathmatrix<coords::float_type>;
    //TODO This needs to be activated again
    //if (!updateStoredDerivatives) {
    //  return *storedDerivativesForRotations;
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
    *storedDerivativesForRotations = Mat(new_xyz.size() * 3, 3);
    storedDerivativesForRotations->set_col(0, X.vectorise_row());
    storedDerivativesForRotations->set_col(1, Y.vectorise_row());
    storedDerivativesForRotations->set_col(2, Z.vectorise_row());
    return *storedDerivativesForRotations;
  }

  Rotator::Rotator(coords::Representation_3D const& reference, std::vector<std::size_t> const& index_vec) :
	storedDerivativesForRotations{std::make_unique<scon::mathmatrix<coords::float_type>>()},
	  updateStoredValues{ true }, updateStoredDerivatives{ true }, reference_{ reference }, rad_gyr_{ radiusOfGyration(reference_) }{
	  for (auto index : index_vec) {
		  indices_.emplace_back(index - 1u);
	  }
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

  Rotator::~Rotator() = default;

  bool Translations::operator==(Translations const & other) const{
    if(indices_.size() != other.indices_.size()) return false;
    for (auto i = 0u; i < indices_.size(); ++i) {
      if (indices_.at(i) != other.indices_.at(i)) return false;
    }
    return true;
  }
  
  std::vector<coords::float_type> Rotation::der_vec(coords::Representation_3D const& cartesians) const {
	  auto const& derivativeMatrix = rotator->rot_der_mat(cartesians);
	  return derivativeMatrix.col_to_std_vector(index());
  }

	bool Rotation::hasIndices(std::vector<std::size_t> const& indices) const {
		return indices.size() == rotator->indices_.size() && isSameSet(indices, rotator->indices_);
	}

	std::vector<std::size_t> Rotation::getIndices() const {
		return rotator->indices_;
	}

	void Rotation::makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) {

	}

	InternalCoordinatesBuilder::InternalCoordinatesBuilder(ic_util::BondGraph const& graph, InternalCoordinates::CartesiansForInternalCoordinates & coordinates, std::vector<std::shared_ptr<InternalCoordinates::Rotator>> & rotators) :
		bondGraph{ graph }, cartesians{ coordinates }, rotatorContainer{ rotators }, numberOfAtoms{ graph.getNumberOfAtoms() } {}

	std::unique_ptr<InternalCoordinate> InternalCoordinatesBuilder::buildBondDistance(std::size_t const a, std::size_t const b) const	{
		if (a >= numberOfAtoms || b >= numberOfAtoms) {
			throw std::runtime_error("Cannot create bond length coordinate: Atom index out of range");
		}
		return std::make_unique<BondDistance>(bondGraph.getAtom(a), bondGraph.getAtom(b));
	}

	std::unique_ptr<InternalCoordinate> InternalCoordinatesBuilder::buildBondAngle(std::size_t const a, std::size_t const b, std::size_t const c) const {
		if (a >= numberOfAtoms || b >= numberOfAtoms || c >= numberOfAtoms) {
			throw std::runtime_error("Cannot create bond angle coordinate: Atom index out of range");
		}
		return std::make_unique<BondAngle>(bondGraph.getAtom(a), bondGraph.getAtom(b), bondGraph.getAtom(c));
	}

	std::unique_ptr<InternalCoordinate> InternalCoordinatesBuilder::buildDihedralAngle(std::size_t const a, std::size_t const b, std::size_t const c, std::size_t const d) const {
		if (a >= numberOfAtoms || b >= numberOfAtoms || c >= numberOfAtoms || d >= numberOfAtoms) {
			throw std::runtime_error("Cannot create dihedral angle coordinate: Atom index out of range");
		}
		return std::make_unique<BondAngle>(bondGraph.getAtom(a), bondGraph.getAtom(b), bondGraph.getAtom(c), bondGraph.getAtom(d));
	}

	std::unique_ptr<InternalCoordinate> InternalCoordinatesBuilder::buildTranslationX(std::vector<std::size_t> const& indices) const {
		if (std::any_of(indices.cbegin(), indices.cend(), [this] (auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create translation in x-direction: Some atom index is out of range");
		}
		return std::make_unique<TranslationX>(indices);
	}

	std::unique_ptr<InternalCoordinate> InternalCoordinatesBuilder::buildTranslationY(std::vector<std::size_t> const& indices) const {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create translation in y-direction: Some atom index is of range");
		}
		return std::make_unique<TranslationY>(indices);
	}

	std::unique_ptr<InternalCoordinate> InternalCoordinatesBuilder::buildTranslationZ(std::vector<std::size_t> const& indices) const {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create translation in z-direction: Some atom index is out of range");
		}
		return std::make_unique<TranslationZ>(indices);
	}

	InternalCoordinatesBuilder::TwoInternals InternalCoordinatesBuilder::buildTranslationXY(std::vector<std::size_t> const& indices) const {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create translation in x- or y-direction: Some atom index is out of range");
		}
		return { std::make_unique<TranslationX>(indices) , std::make_unique<TranslationY>(indices) };
	}

	InternalCoordinatesBuilder::TwoInternals InternalCoordinatesBuilder::buildTranslationXZ(std::vector<std::size_t> const& indices) const {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create translation in x- or z-direction: Some atom index is out of range");
		}
		return { std::make_unique<TranslationX>(indices) , std::make_unique<TranslationZ>(indices) };
	}

	InternalCoordinatesBuilder::TwoInternals InternalCoordinatesBuilder::buildTranslationYZ(std::vector<std::size_t> const& indices) const {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create translation in y- or z-direction: Some atom index is out of range");
		}
		return { std::make_unique<TranslationY>(indices) , std::make_unique<TranslationZ>(indices) };
	}

	InternalCoordinatesBuilder::ThreeInternals InternalCoordinatesBuilder::buildTranslationXYZ(std::vector<std::size_t> const& indices) const {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create translation in x-, y-, or z-direction: Some atom index is out of range");
		}
		return { std::make_unique<TranslationX>(indices) , std::make_unique<TranslationY>(indices), std::make_unique<TranslationZ>(indices) };
	}

	std::unique_ptr<InternalCoordinate> InternalCoordinatesBuilder::buildRotationA(std::vector<std::size_t> const& indices) {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create rotation in plane a: Atom index out of range");
		}
		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, indices)->makeRotations();

		rotatorContainer.emplace_back(curr_rotations.rotator);

		return std::move(curr_rotations.rotationA);
	}

	std::unique_ptr<InternalCoordinate> InternalCoordinatesBuilder::buildRotationB(std::vector<std::size_t> const& indices) {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create rotation in plane b: Atom index out of range");
		}
		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, indices)->makeRotations();

		rotatorContainer.emplace_back(curr_rotations.rotator);

		return std::move(curr_rotations.rotationB);
	}

	std::unique_ptr<InternalCoordinate> InternalCoordinatesBuilder::buildRotationC(std::vector<std::size_t> const& indices) {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create rotation in plane a: Atom index out of range");
		}
		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, indices)->makeRotations();

		rotatorContainer.emplace_back(curr_rotations.rotator);

		return std::move(curr_rotations.rotationC);
	}

	InternalCoordinatesBuilder::TwoInternals InternalCoordinatesBuilder::buildRotationAB(std::vector<std::size_t> const& indices) {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create rotation in plane a and b: Atom index out of range");
		}
		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, indices)->makeRotations();

		rotatorContainer.emplace_back(curr_rotations.rotator);

		return { std::move(curr_rotations.rotationA), std::move(curr_rotations.rotationB) };
	}

	InternalCoordinatesBuilder::TwoInternals InternalCoordinatesBuilder::buildRotationAC(std::vector<std::size_t> const& indices) {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create rotation in plane a and c: Atom index out of range");
		}
		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, indices)->makeRotations();

		rotatorContainer.emplace_back(curr_rotations.rotator);

		return { std::move(curr_rotations.rotationA), std::move(curr_rotations.rotationC) };
	}

	InternalCoordinatesBuilder::TwoInternals InternalCoordinatesBuilder::buildRotationBC(std::vector<std::size_t> const& indices) {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create rotation in plane b and c: Atom index out of range");
		}
		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, indices)->makeRotations();

		rotatorContainer.emplace_back(curr_rotations.rotator);

		return { std::move(curr_rotations.rotationB), std::move(curr_rotations.rotationC) };
	}

	InternalCoordinatesBuilder::ThreeInternals InternalCoordinatesBuilder::buildRotationABC(std::vector<std::size_t> const& indices) {
		if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
			throw std::runtime_error("Cannot create rotation in plane a, b, and c: Atom index out of range");
		}
		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, indices)->makeRotations();

		rotatorContainer.emplace_back(curr_rotations.rotator);

		return { std::move(curr_rotations.rotationA), std::move(curr_rotations.rotationB), std::move(curr_rotations.rotationC) };
	}

}


