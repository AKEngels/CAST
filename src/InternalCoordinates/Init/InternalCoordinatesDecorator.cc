#include "InternalCoordinatesDecorator.h"

#include "../InternalCoordinates.h"
#include "../InternalCoordinateUtilities.h"
#include "../../Scon/scon_vect.h" // for scon::dot

namespace internals {

	void ConstraintAdder::applyConstraints(InternalVector & primitives) {
		for (auto & primitive : primitives) {
			primitive->makeConstrained(manager);
		}
		auto rest = manager->restConstraintsToInternalCoordinates(coordinatesBuilder);
		primitives.insert(primitives.end(), std::make_move_iterator(rest.begin()), std::make_move_iterator(rest.end()));
	}

	InternalVector BondCreator::findInternalCoordinates() {
		auto primitives = graph.getBonds();
		InternalVector result;
		for (auto const& primitive : primitives) {
			result.emplace_back(std::make_unique<InternalCoordinates::BondDistance>(graph.getAtom(std::get<0u>(primitive) - 1u), graph.getAtom(std::get<1u>(primitive) - 1u)));
		}
		return result;
	}

	InternalVector AngleCreator::findInternalCoordinates() {
		auto primitives = graph.getAngles();
		InternalVector result;
		for (auto const& primitive : primitives) {
			result.emplace_back(std::make_unique<InternalCoordinates::BondAngle>(graph.getAtom(std::get<0u>(primitive) - 1u), graph.getAtom(std::get<1u>(primitive) - 1u), graph.getAtom(std::get<2u>(primitive) - 1u)));
		}
		return result;
	}

	InternalVector DihedralCreator::findInternalCoordinates() {
		auto primitives = graph.getDihedrals();
		InternalVector result;
		for (auto const& primitive : primitives) {
			result.emplace_back(std::make_unique<InternalCoordinates::DihedralAngle>(graph.getAtom(std::get<0u>(primitive) - 1u), graph.getAtom(std::get<1u>(primitive) - 1u), graph.getAtom(std::get<2u>(primitive) - 1u), graph.getAtom(std::get<3u>(primitive) - 1u)));
		}
		return result;
	}

	InternalVector TranslationCreator::findInternalCoordinates() {
		auto const& molecules = graph.getMolecules();
		InternalVector result;
		for (auto const& molecule : molecules) {
			result.emplace_back(std::make_unique<InternalCoordinates::TranslationX>(molecule));
			result.emplace_back(std::make_unique<InternalCoordinates::TranslationY>(molecule));
			result.emplace_back(std::make_unique<InternalCoordinates::TranslationZ>(molecule));
		}
		return result;
	}

	InternalVector RotationCereator::findInternalCoordinates() {
		auto const& molecules = graph.getMolecules();
		InternalVector result;
		for (auto const& molecule : molecules) {
			auto curr_rotations = InternalCoordinates::Rotator::buildRotator(coordinates, molecule)->makeRotations();
			result.emplace_back(std::move(curr_rotations.rotationA));
			result.emplace_back(std::move(curr_rotations.rotationB));
			result.emplace_back(std::move(curr_rotations.rotationC));

			rotators.emplace_back(curr_rotations.rotator);
		}
		return result;
	}

}
