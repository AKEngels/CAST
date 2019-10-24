#include"BondDistance.h"

#include"../InternalCoordinateUtilities.h"
#include"../BondGraph/ElementInformations.h"
#include"../../Scon/scon_mathmatrix.h"

namespace internals{

float_type BondDistance::val(scon::mathmatrix<float_type> const& cartesians) const {
	auto a = getAtom(cartesians, index_a_);
	auto b = getAtom(cartesians, index_b_);

	return euclideanLength(a - b);
}

float_type BondDistance::difference(scon::mathmatrix<float_type> const& newCoordinates, scon::mathmatrix<float_type> const& oldCoordinates) const {
	return val(newCoordinates) - val(oldCoordinates);
}

std::pair<CartesianPoint, CartesianPoint> BondDistance::der(scon::mathmatrix<float_type> const&  cartesians) const {

	auto const& a = getAtom(cartesians, index_a_);
	auto const& b = getAtom(cartesians, index_b_);

	auto bond = normalize(a - b);
	return std::make_pair(bond, -bond);
}

scon::mathmatrix<float_type>
BondDistance::der_vec(scon::mathmatrix<float_type> const& cartesians) const {

	auto firstder = der(cartesians);
	
	BmatrixRowCreator rowCreator(cartesians.cols()*cartesians.rows());
	rowCreator.insertAtomDerivative(firstder.first, index_a_);
	rowCreator.insertAtomDerivative(firstder.second, index_b_);

	// TODO: test if this line is mandatory to trigger return optimization
	scon::mathmatrix<float_type> result = rowCreator.getRow();

	return result;
}

//enums hold an integral thus I pass them by value
bool BondDistance::bothElementsInPeriodOne(ic_util::period const atomA, ic_util::period const atomB) const {
	return atomA == ic_util::period::one && atomB == ic_util::period::one;
}

bool BondDistance::oneElementInPeriodOneTheOtherInPeriodTwo(ic_util::period const atomA, ic_util::period const atomB) const {
	return (atomA == ic_util::period::one && atomB == ic_util::period::two)
		|| (atomA == ic_util::period::two && atomB == ic_util::period::one);
}

bool BondDistance::oneElementInPeriodOneTheOtherInPeriodThree(ic_util::period const atomA, ic_util::period const atomB) const {
	return (atomA == ic_util::period::one && atomB == ic_util::period::three)
		|| (atomA == ic_util::period::three && atomB == ic_util::period::one);
}

bool BondDistance::bothElementsInPeriodTwo(ic_util::period const atomA, ic_util::period const atomB) const {
	return atomA == ic_util::period::two && atomB == ic_util::period::two;
}

bool BondDistance::oneElementInPeriodTwoTheOtherInPeriodThree(ic_util::period const atomA, ic_util::period const atomB) const {
	return (atomA == ic_util::period::two && atomB == ic_util::period::three)
		|| (atomA == ic_util::period::three && atomB == ic_util::period::two);
}

coords::float_type BondDistance::hessian_guess(scon::mathmatrix<float_type> const& cartesians) const {
	auto el_a = ic_util::element_period(elem_a_);
	auto el_b = ic_util::element_period(elem_b_);

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

std::string BondDistance::info(scon::mathmatrix<float_type> const & cartesians) const {
	std::ostringstream oss;
	auto l_bohr = val(cartesians);
	oss << "Bond: " << l_bohr << " (a. u.) " << energy::bohr2Ang()*l_bohr << " (Angstrom) || " << index_a_ + 1u << " || " << index_b_ + 1u << " || " << "Constrained: " << std::boolalpha << is_constrained();
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

}