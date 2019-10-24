#include"CartesiansForInternals.h"
#include"CartesiansForTRIC.h"

#include"../../Scon/scon_mathmatrix.h"

namespace internals {

CartesiansForInternalCoordinates::CartesiansForInternalCoordinates(std::unique_ptr<CoordinatesContainer> && cartesians)
	: coordinates{ std::move(cartesians) } {}

CartesiansForInternalCoordinates::CartesiansForInternalCoordinates(scon::mathmatrix<float_type> && cartesians)
	: coordinates{ std::make_unique<CartesianCoordinates>(std::move(cartesians)) } {}

CartesiansForInternalCoordinates::CartesiansForInternalCoordinates(scon::mathmatrix<float_type> const& cartesians)
	: coordinates{ std::make_unique<CartesianCoordinates>(cartesians) } {}

CartesiansForInternalCoordinates & CartesiansForInternalCoordinates::operator=(CartesiansForInternalCoordinates const& cartesians) {
	cartesians.copyTo(*this);
	return *this;
}

CartesiansForInternalCoordinates & CartesiansForInternalCoordinates::operator=(CartesiansForInternalCoordinates && cartesians) {
	std::move(cartesians).copyTo(*this);
	return *this;
}

void CartesiansForInternalCoordinates::copyTo(CartesiansForInternalCoordinates & cartesians) const&{
	cartesians.coordinates->setCoordinates(coordinates->getCartesians());
}

void CartesiansForInternalCoordinates::copyTo(CartesiansForInternalCoordinates & cartesians) && {
	cartesians.coordinates.swap(coordinates);
}


void CartesiansForInternalCoordinates::copyTo(CartesiansForTRIC & cartesians) const&{
	cartesians.coordinates->setCoordinates(coordinates->getCartesians());
}
void CartesiansForInternalCoordinates::copyTo(CartesiansForTRIC & cartesians) && {
	cartesians.coordinates.swap(coordinates);
}

std::pair<coords::float_type, coords::float_type> CartesiansForInternalCoordinates::displacementRmsValAndMaxTwoStructures(scon::mathmatrix<float_type> const& other) const {
	return ic_rotation::displacementRmsValAndMaxTwoStructures<scon::mathmatrix>(coordinates, other);
}

std::pair<coords::float_type, coords::float_type> CartesiansForInternalCoordinates::displacementRmsValAndMaxTwoStructures(CartesiansForInternalCoordinates const& other) const {
	return displacementRmsValAndMaxTwoStructures(other.coordinates->getCartesians());
}

scon::mathmatrix<float_type> const & CartesiansForInternalCoordinates::getCartesians() const
{
	// TODO: hier Rückgabeanweisung eingeben
}

scon::mathmatrix<float_type> CartesiansForInternalCoordinates::inBohr() const
{
	return scon::mathmatrix<float_type>();
}

scon::mathmatrix<float_type> CartesiansForInternalCoordinates::inAngstrom() const
{
	return scon::mathmatrix<float_type>();
}

//TODO check if this is working
scon::mathmatrix<float_type> CartesiansForInternalCoordinates::at(std::size_t const i) const {
	scon::mathmatrix<float_type> result(1,3);
	result(0,0) = coordinates->getCartesians()(i,0);
	result(0, 1) = coordinates->getCartesians()(i, 1);
	result(0, 2) = coordinates->getCartesians()(i, 2);
	return result;
}

coords::float_type CartesiansForInternalCoordinates::getInternalValue(InternalCoordinate const& in) const { return in.val(coordinates->getCartesians()); }

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

}