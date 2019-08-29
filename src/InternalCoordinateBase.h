/**
CAST 3
InternalCoordinateBase.h
Purpose: base (intreface) class for PrimitiveInternalCoordinates and its decorators


@author Julian Erdmannsdörfer, Christian Schärf
@version 3.0
*/

#ifndef INTERNAL_COORDINATES_BASE
#define INTERNAL_COORDINATES_BASE

#include<vector>

#include "configuration.h"

namespace internals {
	using float_type = double;
}

namespace scon {
	template<typename T>
	class c3;
}

namespace coords {
	using r3 = scon::c3<internals::float_type>;
	using Cartesian_Point = r3;
}


namespace InternalCoordinates {
	class InternalCoordinate;
	class Rotator;
	class temporaryCartesian;
	class CartesiansForInternalCoordinates;
}



namespace internals {
	class ICDecoratorBase;
	class ICDecorator;
	class AppropriateStepFinder;
	class InternalToCartesianConverter;
	using CartesianType = InternalCoordinates::CartesiansForInternalCoordinates;
	using InternalVec = std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>;
	using IndexVec = std::vector<std::vector<std::size_t>>;

	class AbstractConstraintManager {
	public:
		using ConstrainVec = std::vector<std::shared_ptr<config::AbstractConstraint>>;
		virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintPrimitive(std::vector<std::size_t> const&) = 0;
		virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintTrans(std::vector<std::size_t> const&, config::Constraint const) = 0;
		virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintRot(std::vector<std::size_t> const&, config::Constraint const) = 0;
		virtual ConstrainVec getConstraintsOfType(config::Constraint const) = 0;
	};
}

#endif // INTERNAL_COORDIANTES_BASE