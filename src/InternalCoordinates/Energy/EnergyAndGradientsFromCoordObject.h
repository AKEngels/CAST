#ifndef CAST_INTERNALCOORDINATES_ENERGY_ENERGYANDGRADIENTSFROMCOORDSOBJECT_H_
#define CAST_INTERNALCOORDINATES_ENERGY_ENERGYANDGRADIENTSFROMCOORDSOBJECT_H_

#include<memory>

#include"EnergyAndGradients.h"

namespace internals {

	class EnergyAndGradientsFromCoordsObject : public EnergyAndGradients {
	public:
		EnergyAndGradientsFromCoordsObject(coords::Coordinates &);
		~EnergyAndGradientsFromCoordsObject() override;
		void evaluate() override;
		float_type getEnergy() const override;
		scon::mathmatrix<float_type> const& getGradients() const override;
	private:
		coords::Coordinates & coords;
		std::unique_ptr<scon::mathmatrix<float_type>> gradients;
		float_type energy;
	};

	EnergyAndGradients::~EnergyAndGradients() = default;

}

#endif // CAST_INTERNALCOORDINATES_ENERGY_ENERGYANDGRADIENTSFROMCOORDSOBJECT_H_