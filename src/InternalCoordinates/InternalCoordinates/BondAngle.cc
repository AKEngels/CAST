#include "BondAngle.h"

#include "../InternalCoordinateUtilities.h"
#include "../BondGraph/ElementInformations.h"

namespace internals{


float_type BondAngle::value(Eigen::MatrixXd const& cartesians) const {
	auto a = cartesians.row(index_a_);
	auto b = cartesians.row(index_b_);
	auto c = cartesians.row(index_c_);

	Eigen::Vector3d u = a - b;
	Eigen::Vector3d v = c - b;
	auto uXv = u.cross(v);
	auto uDv = u.dot(v);
	auto l = euclideanLength(uXv);
	return std::atan2(l, uDv);
}

float_type BondAngle::difference(Eigen::MatrixXd const& newCoordinates, Eigen::MatrixXd const& oldCoordinates) const {
	return value(newCoordinates) - value(oldCoordinates);
}

std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>
BondAngle::derivatives(Eigen::MatrixXd const& cartesians) const {
	auto constexpr normOfVec111 = 1.7320508075688771932;

	auto a = cartesians.row(index_a_);
	auto b = cartesians.row(index_b_);
	auto c = cartesians.row(index_c_);

	Eigen::Vector3d u = a - b;
	Eigen::Vector3d v = c - b;
	auto lu = u.norm();
	auto lv = v.norm();
	u /= lu;
	v /= lv;
	auto cp1 = Eigen::Vector3d{ 1.0, -1.0, 1.0 } / normOfVec111;
	auto cp2 = Eigen::Vector3d{ -1.0, 1.0, 1.0 } / normOfVec111;
	Eigen::Vector3d w{ 0.0, 0.0, 0.0 };
	auto constexpr epsilon{ 0.01 };
	if (std::fabs(u.dot(v)) > (1. - epsilon)) {
		if (std::fabs(u.dot(cp1)) > (1. - epsilon)) {
			w = u.cross(cp2);
		}
		else {
			w = u.cross(cp1);
		}
	}
	else {
		w = u.cross(v);
	}
	w /= w.norm();
	auto ad0 = u.cross(w) / lu;
	auto ad1 = w.cross(v) / lv;
	//                      a     b           c
	//                      m     o           n
	return std::make_tuple(ad0, -ad0 - ad1, ad1);
}

Eigen::VectorXd
BondAngle::derivativeVector(scon::mathmatrix<float_type> const& cartesians) const {
	auto firstDerivatives = derivatives(cartesians);

	Eigen::VectorXd result = Eigen::VectorXd::Zero(cartesians.size());

	Eigen::Map<RowMajorMatrixXd> mappedResult(result.data(), cartesians.rows(), cartesians.cols());

	mappedResult.row(index_a_) = std::get<0u>(firstder);
	mappedResult.row(index_b_) = std::get<1u>(firstder);
	mappedResult.row(index_c_) = std::get<2u>(firstder);

	return result;
}

float_type BondAngle::hessianGuess(Eigen::MatrixXd const& /*cartesians*/) const {

	auto el_a = ic_util::element_period(elem_a_);
	auto el_b = ic_util::element_period(elem_b_);
	auto el_c = ic_util::element_period(elem_c_);
	if (el_a == ic_util::period::one ||
		el_b == ic_util::period::one ||
		el_c == ic_util::period::one) {
		return 0.160;
	}
	else {
		return 0.250;
	}
}

std::string BondAngle::info(Eigen::MatrixXd const & cartesians) const {
	std::ostringstream oss;
	oss << "Angle: " << value(cartesians) * SCON_180PI << " || " << index_a_ + 1 << " || " << index_b_ + 1 << " || " << index_c_ + 1 << " || " << "Constrained: " << std::boolalpha << is_constrained();
	return oss.str();
}

bool BondAngle::hasIndices(std::vector<std::size_t> const& indices) const {
	return indices.size() == 3u && ((index_a_ == indices[0u] && index_b_ == indices[1u] && index_c_ == indices[2u]) || (index_a_ == indices[2u] && index_b_ == indices[1u] && index_c_ == indices[0u]));
}

std::vector<std::size_t> BondAngle::getIndices() const {
	return { index_a_, index_b_, index_c_ };
}

void BondAngle::makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) {
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

}