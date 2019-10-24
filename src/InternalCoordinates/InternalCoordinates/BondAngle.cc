#include "BondAngle.h"

#include "../InternalCoordinateUtilities.h"
#include "../BondGraph/ElementInformations.h"
#include "../../Scon/scon_mathmatrix.h"

namespace internals{

float_type BondAngle::val(scon::mathmatrix<float_type> const& cartesians) const {
	auto a = getAtom(cartesians, index_a_);
	auto b = getAtom(cartesians, index_b_);
	auto c = getAtom(cartesians, index_c_);

	auto u = a - b;
	auto v = c - b;
	auto uXv = crossProduct(u, v);
	auto uDv = dotProduct(u, v);
	auto l = euclideanLength(uXv);
	return std::atan2(l, uDv);
}

float_type BondAngle::difference(scon::mathmatrix<float_type> const& newCoordinates, scon::mathmatrix<float_type> const& oldCoordinates) const {
	return val(newCoordinates) - val(oldCoordinates);
}

std::tuple<CartesianPoint, CartesianPoint, CartesianPoint>
BondAngle::der(scon::mathmatrix<float_type> const& cartesians) const {
	auto a = getAtom(cartesians, index_a_);
	auto b = getAtom(cartesians, index_b_);
	auto c = getAtom(cartesians, index_c_);

	auto u = a - b;
	auto v = c - b;
	auto lu = euclideanLength(u);
	auto lv = euclideanLength(v);
	u = normalize(u);
	v = normalize(v);
	auto cp1 = normalize(CartesianPoint{ 1.0, -1.0, 1.0 });
	auto cp2 = normalize(CartesianPoint{ -1.0, 1.0, 1.0 });
	CartesianPoint w_p{ 0.0, 0.0, 0.0 };
	auto constexpr epsilon{ 0.01 };
	if (std::fabs(dotProduct(u, v)) > (1. - epsilon)) {
		if (std::fabs(dotProduct(u, cp1)) > (1. - epsilon)) {
			w_p = crossProduct(u, cp2);
		}
		else {
			w_p = crossProduct(u, cp1);
		}
	}
	else {
		w_p = crossProduct(u, v);
	}
	auto w = normalize(w_p);
	auto ad0 = crossProduct(u, w) / lu;
	auto ad1 = crossProduct(w, v) / lv;
	//                      a     b           c
	//                      m     o           n
	return std::make_tuple(ad0, -ad0 - ad1, ad1);
}

scon::mathmatrix<float_type>
BondAngle::der_vec(scon::mathmatrix<float_type> const& cartesians) const {
	using scon::c3;


	auto firstder = der(cartesians);

	BmatrixRowCreator rowCreator(cartesians.cols()*cartesians.rows());
	rowCreator.insertAtomDerivative(std::get<0u>(firstder), index_a_);
	rowCreator.insertAtomDerivative(std::get<1u>(firstder), index_b_);
	rowCreator.insertAtomDerivative(std::get<2u>(firstder), index_b_);

	// TODO: test if this line is mandatory to trigger return optimization
	scon::mathmatrix<float_type> result = rowCreator.getRow();

	return result;
}

float_type BondAngle::hessian_guess(scon::mathmatrix<float_type> const& /*cartesians*/) const {

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

std::string BondAngle::info(scon::mathmatrix<float_type> const & cartesians) const {
	std::ostringstream oss;
	oss << "Angle: " << val(cartesians) * SCON_180PI << " || " << index_a_ + 1 << " || " << index_b_ + 1 << " || " << index_c_ + 1 << " || " << "Constrained: " << std::boolalpha << is_constrained();
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