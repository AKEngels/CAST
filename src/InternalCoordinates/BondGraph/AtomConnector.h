#ifndef CAST_INTERNALCOORDINATES_BONDGRAPH_ATOMCONNECTOR_H_
#define CAST_INTERNALCOORDINATES_BONDGRAPH_ATOMCONNECTOR_H_

#include "../../coords.h"

#include "../InternalCoordinatesAliases.h"

#include<vector>

struct AtomConnector {
	using returnType = std::vector<std::pair<std::size_t, std::size_t>>;

	AtomConnector(std::vector<std::string> const& elem_vec, coords::Representation_3D const& cp_vec) : sequenceOfSymbols{ elem_vec }, cartesianRepresentation{ cp_vec }, firstAtomIndex{ 0u }, secondAtomIndex{ 0u } {}
	returnType operator()();

private:
	double getThresholdForBeingNotConnected(std::string const& oneAtom, std::string const& otherAtom);
	void findTheFirstAtom(returnType & connectedAtoms);
	void findAtomWithIndexHigherThanFirst(returnType & connectedAtoms);
	void connectIfCloseEnough(returnType & connectedAtoms);
	bool areTheyCloseEnough();

	std::vector<std::string> const& sequenceOfSymbols;
	coords::Representation_3D const& cartesianRepresentation;

	std::size_t firstAtomIndex;
	std::size_t secondAtomIndex;
};

#endif CAST_INTERNALCOORDINATES_BONDGRAPH_ATOMCONNECTOR_H_