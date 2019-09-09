#include "Constraints.h"

#include"../InternalCoordinates.h"
#include "../BondGraph/BondGraph.h"

std::unique_ptr<InternalCoordinates::InternalCoordinate> BondConstraint::makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const{
	return builder.buildBondDistance(atomIndices[0u] - 1u, atomIndices[1u] - 1u);
}

std::unique_ptr<InternalCoordinates::InternalCoordinate> AngleConstraint::makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const{
	return builder.buildBondAngle(atomIndices[0u] - 1u, atomIndices[1u] - 1u, atomIndices[2u] - 1u);
}

std::unique_ptr<InternalCoordinates::InternalCoordinate> DihedralConstraint::makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const{
	return builder.buildDihedralAngle(atomIndices[0u] - 1u, atomIndices[1u] - 1u, atomIndices[2u] - 1u, atomIndices[3u] - 1u);
}

std::unique_ptr<InternalCoordinates::InternalCoordinate> TranslationXConstraint::makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const{
	return builder.buildTranslationX(atomIndices);
}

std::unique_ptr<InternalCoordinates::InternalCoordinate> TranslationYConstraint::makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const{
	return builder.buildTranslationY(atomIndices);
}

std::unique_ptr<InternalCoordinates::InternalCoordinate> TranslationZConstraint::makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const{
	return builder.buildTranslationZ(atomIndices);
}

// Thiere should be rotational constraints wich pair up. E. g.: rot AB or ABC or something like that
std::unique_ptr<InternalCoordinates::InternalCoordinate> RotationAConstraint::makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const{
	return builder.buildRotationA(atomIndices);
}

std::unique_ptr<InternalCoordinates::InternalCoordinate> RotationBConstraint::makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const{
	return builder.buildRotationB(atomIndices);
}

std::unique_ptr<InternalCoordinates::InternalCoordinate> RotationCConstraint::makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const {
	return builder.buildRotationC(atomIndices);
}
