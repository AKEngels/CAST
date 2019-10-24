/**
CAST 3
ConstrainedInternalCoordinate.h
Purpose: Internal Coordinate System for constrained optimization


@author Christian Schärf, Julian Erdmannsdörfer
@version 3.0
*/

#ifndef CAST_INTERNALCOORDINATES_CONSTRAINEDINTERNALCOORDINATES_H_
#define CAST_INTERNALCOORDINATES_CONSTRAINEDINTERNALCOORDINATES_H_

#include "PrimitiveInternalCoordinates.h"

namespace internals {
	class ConstrainedInternalCoordinates : public PrimitiveInternalCoordinates {
	public:
		using PrimitiveInternalCoordinates::PrimitiveInternalCoordinates;

		virtual scon::mathmatrix<coords::float_type> projectorMatrix(CartesianType const& cartesian) override;
		virtual scon::mathmatrix<coords::float_type> constraintMatrix() const;
	};
}

#endif // CONSTRAINED_INTERNAL_COORDINATES