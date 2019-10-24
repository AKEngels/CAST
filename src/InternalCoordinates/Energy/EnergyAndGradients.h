#ifndef CAST_INTERNALCOORDINATES_ENERGY_ENERGYANDGRADIENTS_H_
#define CAST_INTERNALCOORDINATES_ENERGY_ENERGYANDGRADIENTS_H_

#include"../InternalCoordinatesAliases.h"

namespace internals{

class EnergyAndGradients {
	public:
	virtual ~EnergyAndGradients() = 0;
	virtual void evaluate() = 0;
	virtual float_type getEnergy() const = 0;
	virtual scon::mathmatrix<float_type> const& getGradients() const = 0;
};

EnergyAndGradients::~EnergyAndGradients() = default;

}

#endif // CAST_INTERNALCOORDINATES_ENERGY_ENERGYANDGRADIENTS_H_