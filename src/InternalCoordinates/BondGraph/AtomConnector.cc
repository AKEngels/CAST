#include "AtomConnector.h"
#include "ElementInformations.h"

std::vector<std::pair<std::size_t, std::size_t>> AtomConnector::operator()() {
	std::vector<std::pair<std::size_t, std::size_t>> connectedAtoms;
	findTheFirstAtom(connectedAtoms);
	return connectedAtoms;
}

void AtomConnector::findTheFirstAtom(returnType & connectedAtoms) {
	for (firstAtomIndex = 0u; firstAtomIndex < sequenceOfSymbols.size(); ++firstAtomIndex) {
		findAtomWithIndexHigherThanFirst(connectedAtoms);
	}
}

void AtomConnector::findAtomWithIndexHigherThanFirst(returnType & connectedAtoms) {
	for (secondAtomIndex = firstAtomIndex + 1u; secondAtomIndex < sequenceOfSymbols.size(); ++secondAtomIndex) {
		connectIfCloseEnough(connectedAtoms);
	}
}

void AtomConnector::connectIfCloseEnough(returnType & connectedAtoms) {
	if (areTheyCloseEnough()) {
		connectedAtoms.emplace_back(firstAtomIndex, secondAtomIndex);
	}
}

bool AtomConnector::areTheyCloseEnough() {
	double const threshold = getThresholdForBeingNotConnected(sequenceOfSymbols.at(firstAtomIndex), sequenceOfSymbols.at(secondAtomIndex));
	double const actualDistance = scon::len(cartesianRepresentation.at(firstAtomIndex) - cartesianRepresentation.at(secondAtomIndex));
	return actualDistance < threshold;
}



double AtomConnector::getThresholdForBeingNotConnected(std::string const& oneAtom, std::string const& otherAtom) {
	return 1.2 * (element_radius(oneAtom) + element_radius(otherAtom));
}
