#include "ConstrainedInternalCoordinates.h"

#include "InternalCoordinates/InternalCoordinate.h"
#include "../Scon/scon_mathmatrix.h"

namespace internals {
	/*std::unique_ptr<AppropriateStepFinder> ConstrainedInternalCoordinates::constructStepFinder(
		InternalToCartesianConverter const& converter,
		scon::mathmatrix<coords::float_type> const& gradients,
		scon::mathmatrix<coords::float_type> const& hessianMatrix,
		CartesianType const& cartesians
	) {
		auto pmat = projectorMatrix(cartesians);
		auto imat = scon::mathmatrix<coords::float_type>::identity(pmat.rows(), pmat.cols());
		auto projectedHessian = pmat * hessianMatrix * pmat + (imat - pmat) * static_cast<coords::float_type>(1000.0);
		auto projectedGradient = pmat * gradients;

		return std::make_unique<AppropriateStepFinder>(converter, projectedGradient, projectedHessian);
	}*/

	scon::mathmatrix<coords::float_type> ConstrainedInternalCoordinates::projectorMatrix(CartesianType const& cartesian) {
		auto P = PrimitiveInternalCoordinates::projectorMatrix(cartesian);
		auto C = constraintMatrix();
		auto CPC = C * P * C;
		return P - P * C * CPC.pinv() * C * P;
	}

	scon::mathmatrix<coords::float_type> ConstrainedInternalCoordinates::constraintMatrix() const {
		auto s = primitive_internals.size();
		auto ret = scon::mathmatrix<coords::float_type>::zero(s, s);
		for (std::size_t i = 0; i < s; ++i) {
			ret(i, i) = primitive_internals.at(i)->is_constrained() ? 1.0 : 0.0;
		}
		return ret;
	}
}
