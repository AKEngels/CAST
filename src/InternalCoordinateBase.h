/**
CAST 3
InternalCoordinateBase.h
Purpose: base (intreface) class for PrimitiveInternalCoordinates and its decorators


@author Julian Erdmannsdörfer, Christian Schärf
@version 3.0
*/

#ifndef INTERNAL_COORDINATES_BASE
#define INTERNAL_COORDINATES_BASE

#include "InternalCoordinates.h"
#include "BondGraph.h"

namespace internals {

	class AbstractConstraintManager {
	public:
		using ConstrainVec = std::vector<std::shared_ptr<config::AbstractConstraint>>;
		virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintPrimitive(std::vector<std::size_t> const&) = 0;
		virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintTrans(std::vector<std::size_t> const&, config::Constraint const) = 0;
		virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintRot(std::vector<std::size_t> const&, config::Constraint const) = 0;
		virtual ConstrainVec getConstraintsOfType(config::Constraint const) = 0;
		//virtual ConstrainVec && getAllConstraints() = 0;
	};

	using CartesianType = InternalCoordinates::CartesiansForInternalCoordinates;
	using BondGraph = ic_util::Graph<ic_util::Node>;
	using IndexVec = std::vector<std::vector<std::size_t>>;
	using InternalVec = std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>;
}

#endif // INTERNAL_COORDIANTES_BASE