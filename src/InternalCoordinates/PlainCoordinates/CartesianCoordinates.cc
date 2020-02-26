#include "CartesianCoordinates.h"

namespace internals{

CartesianCoordinates::CartesianCoordinates(std::unique_ptr<AbstractMatrix> && cartesians) : cartesians{ std::move(cartesians) } {}

AbstractMatrix const& CartesianCoordinates::getCartesians() const {
	return *cartesians;
}

AbstractMatrix & CartesianCoordinates::getCartesians() {
	return *cartesians;
}

void CartesianCoordinates::setCoordinates(std::unique_ptr<AbstractMatrix> && coords) {

}

}