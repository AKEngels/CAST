#ifndef CAST_INTERNALCOORDINATES_PLAINCOORDINATES_COORDINATESCONTAINER_H_
#define CAST_INTERNALCOORDINATES_PLAINCOORDINATES_COORDINATESCONTAINER_H_

#include "../Matrix/AbstractMatrix.h"
#include "../InternalCoordinatesAliases.h"

namespace internals{

class CoordinatesContainer {
public:
	virtual AbstractMatrix const& getCartesians() const = 0;
	virtual AbstractMatrix & getCartesians() = 0;
	virtual void setCoordinates(std::unique_ptr<AbstractMatrix> && coords) = 0;
	virtual ~CoordinatesContainer() = 0;
};

CoordinatesContainer::~CoordinatesContainer() = default;

}

#endif // CAST_INTERNALCOORDINATES_PLAINCOORDINATES_COORDINATESCONTAINER_H_