#ifndef CAST_INTERNALCOORDINATES_ENERGY_ENERGYANDGRADIENTSFROMCOORDSOBJECT_H_
#define CAST_INTERNALCOORDINATES_ENERGY_ENERGYANDGRADIENTSFROMCOORDSOBJECT_H_

#include<memory>

#include "../../coords.h"
#include"EnergyAndGradients.h"

namespace internals {

	class EnergyAndGradientsFromCoordsObject : public EnergyAndGradients {
	public:
		virtual ~EnergyAndGradientsFromCoordsObject();

		EnergyAndGradientsFromCoordsObject(coords::Coordinates &);
		void evaluate() override;
		float_type getEnergy() const override;
		scon::mathmatrix<float_type> const& getGradients() const override;
	private:
		coords::Coordinates & coords;
		std::unique_ptr<scon::mathmatrix<float_type>> gradients;
		float_type energy;
	};

}

#endif // CAST_INTERNALCOORDINATES_ENERGY_ENERGYANDGRADIENTSFROMCOORDSOBJECT_H_