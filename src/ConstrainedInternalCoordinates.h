/**
CAST 3
ConstrainedInternalCoordinate.h
Purpose: Internal Coordinate System for constrained optimization


@author Christian Schärf, Julian Erdmannsdörfer
@version 3.0
*/

#ifndef CONSTRAINED_INTERNAL_COORDINATES
#define CONSTRAINED_INTERNAL_COORDINATES

#include "PrimitiveInternalCoordinates.h"

namespace internals{
  class ConstrainedInternalCoordinates : public PrimitiveInternalCoordinates{
  public:
    
    virtual scon::mathmatrix<coords::float_type> projectorMatrix(CartesianType const& cartesian) override;
    virtual scon::mathmatrix<coords::float_type> constraintMatrix() const;
    
    virtual std::unique_ptr<AppropriateStepFinder> constructStepFinder(
      InternalToCartesianConverter const& converter,
      scon::mathmatrix<coords::float_type> const& gradients,
      scon::mathmatrix<coords::float_type> const& hessian,
      CartesianType const& cartesians
    ) override;
  };
}

#endif // CONSTRAINED_INTERNAL_COORDINATES