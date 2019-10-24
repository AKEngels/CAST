#ifndef CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANCOORDINATES_H_
#define CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANCOORDINATES_H_

#include "CoordinatesContainer.h"

#include<memory>

namespace internals{

class CartesianCoordinates : public CoordinatesContainer {
	public:
		CartesianCoordinates(scon::mathmatrix<float_type> && cartesians);
		CartesianCoordinates(scon::mathmatrix<float_type> const& cartesians);
		virtual scon::mathmatrix<float_type> const& getCartesians() const override;
	protected:
	std::unique_ptr<scon::mathmatrix<float_type>> cartesians;
};

}

#endif // CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANCOORDINATES_H_