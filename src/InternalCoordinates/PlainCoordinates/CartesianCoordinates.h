#ifndef CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANCOORDINATES_H_
#define CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANCOORDINATES_H_

#include "CoordinatesContainer.h"

#include<memory>

namespace internals{

class CartesianCoordinates : public CoordinatesContainer {
	public:
		CartesianCoordinates(std::unique_ptr<AbstractMatrix> && cartesians);

		virtual AbstractMatrix const& getCartesians() const override;
		virtual AbstractMatrix & getCartesians() override;
		virtual void setCoordinates(std::unique_ptr<AbstractMatrix> && coords) override;

		virtual ~CartesianCoordinates() = default;
	protected:
		std::unique_ptr<AbstractMatrix> cartesians;
};

}

#endif // CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANCOORDINATES_H_