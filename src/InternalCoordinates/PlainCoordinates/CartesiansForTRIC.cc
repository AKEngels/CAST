#include "CartesiansForTRIC.h"

namespace internals{

CartesiansForTRIC & CartesiansForTRIC::operator=(CartesiansForInternalCoordinates const& cartesians) {
	cartesians.copyTo(*this);
	return *this;
}

CartesiansForTRIC & CartesiansForTRIC::operator=(CartesiansForInternalCoordinates && cartesians) {
	std::move(cartesians).copyTo(*this);
	return *this;
}

void CartesiansForTRIC::copyTo(CartesiansForInternalCoordinates & cartesians) const& {
	cartesians.coordinates->setCoordinates(coordinates->getCartesians());
}

void CartesiansForTRIC::copyTo(CartesiansForInternalCoordinates & cartesians) && {
	cartesians.coordinates.swap(coordinates);
}

void CartesiansForTRIC::copyTo(CartesiansForTRIC & cartesians) const& {
	cartesians.coordinates->setCoordinates(coordinates->getCartesians());
	cartesians.observerList = observerList;
}

void CartesiansForTRIC::copyTo(CartesiansForTRIC & cartesians) && {
	cartesians.coordinates.swap(coordinates);
	cartesians.observerList = std::move(observerList);
}

}