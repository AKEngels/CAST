#ifndef CAST_INTERNALCOORDINATES_COORDINATESYSTEM_H_
#define CAST_INTERNALCOORDINATES_COORDINATESYSTEM_H_

#include"InternalCoordinatesAliases.h"

namespace internals {

class CoordinateSystem {
public:
	virtual ~CoordinateSystem() = 0;

	virtual scon::mathmatrix<float_type> initialHessian() const = 0;
	virtual scon::mathmatrix<float_type> getGradients(scon::mathmatrix<float_type> const&) const = 0;
	virtual scon::mathmatrix<float_type> getValues() const = 0;
	virtual scon::mathmatrix<float_type> const& getCartesians() const = 0;
	virtual void reset() = 0;
	virtual StepFinder constructStepFinder() const = 0;
};

CoordinateSystem::~CoordinateSystem() = default;

}

#endif // CAST_INTERNALCOORDINATES_COORDINATESYSTEM_H_