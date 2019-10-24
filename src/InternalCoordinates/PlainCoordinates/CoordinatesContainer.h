#ifndef CAST_INTERNALCOORDINATES_PLAINCOORDINATES_COORDINATESCONTAINER_H_
#define CAST_INTERNALCOORDINATES_PLAINCOORDINATES_COORDINATESCONTAINER_H_

#include "../InternalCoordinatesAliases.h"

namespace internals{

class CoordinatesContainer {
public:
	virtual scon::mathmatrix<float_type> const& getCartesians() const = 0;
	virtual scon::mathmatrix<float_type> & getCartesians() = 0;
	virtual void setCoordinates(scon::mathmatrix<float_type> && coords) = 0;
	virtual void setCoordinates(scon::mathmatrix<float_type> const& coords) = 0;
	virtual ~CoordinatesContainer() = 0;
};

CoordinatesContainer::~CoordinatesContainer() = default;

}

#endif // CAST_INTERNALCOORDINATES_PLAINCOORDINATES_COORDINATESCONTAINER_H_