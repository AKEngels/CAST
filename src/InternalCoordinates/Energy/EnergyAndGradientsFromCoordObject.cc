#include "EnergyAndGradientsFromCoordObject.h"

#include "../../coords.h"
#include "../../Scon/scon_mathmatrix.h"
#include"../InternalCoordinateUtilities.h"

internals::EnergyAndGradientsFromCoordsObject::EnergyAndGradientsFromCoordsObject(coords::Coordinates & coordinatesObject) : coords{ coordinatesObject }, gradients{ std::make_unique<scon::mathmatrix<float_type>>() }, energy{ 0.0 } {
}

internals::EnergyAndGradientsFromCoordsObject::~EnergyAndGradientsFromCoordsObject() = default;

void internals::EnergyAndGradientsFromCoordsObject::evaluate(){

	energy = coords.g() / energy::au2kcal_mol;
	*gradients = scon::mathmatrix<internals::float_type>::col_from_vec(ic_util::flatten_c3_vec(
		ic_util::grads_to_bohr(coords.g_xyz())
	));
}

internals::float_type internals::EnergyAndGradientsFromCoordsObject::getEnergy() const{
	return energy;
}

scon::mathmatrix<internals::float_type> const & internals::EnergyAndGradientsFromCoordsObject::getGradients() const{
	return *gradients;
}
