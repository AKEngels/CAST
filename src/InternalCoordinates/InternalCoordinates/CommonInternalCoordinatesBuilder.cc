#include "CommonInternalCoordinatesBuilder.h"

#include"AllInternalCoordinates.h"
#include "../BondGraph/BondGraph.h"

namespace internals {

CommonInternalCoordinatesBuilder::CommonInternalCoordinatesBuilder(ic_util::BondGraph const& graph, InternalCoordinates::CartesiansForInternalCoordinates & coordinates, std::vector<std::shared_ptr<InternalCoordinates::Rotator>> & rotators) :
	bondGraph{ graph }, cartesians{ coordinates }, rotatorContainer{ rotators }, numberOfAtoms{ graph.getNumberOfAtoms() } {}

std::unique_ptr<InternalCoordinate> CommonInternalCoordinatesBuilder::buildBondDistance(std::size_t const a, std::size_t const b) const {
	if (a >= numberOfAtoms || b >= numberOfAtoms) {
		throw std::runtime_error("Cannot create bond length coordinate: Atom index out of range");
	}
	return std::make_unique<BondDistance>(bondGraph.getAtom(a), bondGraph.getAtom(b));
}

std::unique_ptr<InternalCoordinate> CommonInternalCoordinatesBuilder::buildBondAngle(std::size_t const a, std::size_t const b, std::size_t const c) const {
	if (a >= numberOfAtoms || b >= numberOfAtoms || c >= numberOfAtoms) {
		throw std::runtime_error("Cannot create bond angle coordinate: Atom index out of range");
	}
	return std::make_unique<BondAngle>(bondGraph.getAtom(a), bondGraph.getAtom(b), bondGraph.getAtom(c));
}

std::unique_ptr<InternalCoordinate> CommonInternalCoordinatesBuilder::buildDihedralAngle(std::size_t const a, std::size_t const b, std::size_t const c, std::size_t const d) const {
	if (a >= numberOfAtoms || b >= numberOfAtoms || c >= numberOfAtoms || d >= numberOfAtoms) {
		throw std::runtime_error("Cannot create dihedral angle coordinate: Atom index out of range");
	}
	return std::make_unique<BondAngle>(bondGraph.getAtom(a), bondGraph.getAtom(b), bondGraph.getAtom(c), bondGraph.getAtom(d));
}

std::unique_ptr<InternalCoordinate> CommonInternalCoordinatesBuilder::buildTranslationX(std::vector<std::size_t> const& indices) const {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create translation in x-direction: Some atom index is out of range");
	}
	return std::make_unique<TranslationX>(indices);
}

std::unique_ptr<InternalCoordinate> CommonInternalCoordinatesBuilder::buildTranslationY(std::vector<std::size_t> const& indices) const {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create translation in y-direction: Some atom index is of range");
	}
	return std::make_unique<TranslationY>(indices);
}

std::unique_ptr<InternalCoordinate> CommonInternalCoordinatesBuilder::buildTranslationZ(std::vector<std::size_t> const& indices) const {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create translation in z-direction: Some atom index is out of range");
	}
	return std::make_unique<TranslationZ>(indices);
}

CommonInternalCoordinatesBuilder::TwoInternals CommonInternalCoordinatesBuilder::buildTranslationXY(std::vector<std::size_t> const& indices) const {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create translation in x- or y-direction: Some atom index is out of range");
	}
	return { std::make_unique<TranslationX>(indices) , std::make_unique<TranslationY>(indices) };
}

CommonInternalCoordinatesBuilder::TwoInternals CommonInternalCoordinatesBuilder::buildTranslationXZ(std::vector<std::size_t> const& indices) const {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create translation in x- or z-direction: Some atom index is out of range");
	}
	return { std::make_unique<TranslationX>(indices) , std::make_unique<TranslationZ>(indices) };
}

CommonInternalCoordinatesBuilder::TwoInternals CommonInternalCoordinatesBuilder::buildTranslationYZ(std::vector<std::size_t> const& indices) const {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create translation in y- or z-direction: Some atom index is out of range");
	}
	return { std::make_unique<TranslationY>(indices) , std::make_unique<TranslationZ>(indices) };
}

CommonInternalCoordinatesBuilder::ThreeInternals CommonInternalCoordinatesBuilder::buildTranslationXYZ(std::vector<std::size_t> const& indices) const {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create translation in x-, y-, or z-direction: Some atom index is out of range");
	}
	return { std::make_unique<TranslationX>(indices) , std::make_unique<TranslationY>(indices), std::make_unique<TranslationZ>(indices) };
}

std::unique_ptr<InternalCoordinate> CommonInternalCoordinatesBuilder::buildRotationA(std::vector<std::size_t> const& indices) {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create rotation in plane a: Atom index out of range");
	}
	auto curr_rotations = Rotator::buildRotator(cartesians, indices)->makeRotations();

	rotatorContainer.emplace_back(curr_rotations.rotator);

	return std::move(curr_rotations.rotationA);
}

std::unique_ptr<InternalCoordinate> CommonInternalCoordinatesBuilder::buildRotationB(std::vector<std::size_t> const& indices) {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create rotation in plane b: Atom index out of range");
	}
	auto curr_rotations = Rotator::buildRotator(cartesians, indices)->makeRotations();

	rotatorContainer.emplace_back(curr_rotations.rotator);

	return std::move(curr_rotations.rotationB);
}

std::unique_ptr<InternalCoordinate> CommonInternalCoordinatesBuilder::buildRotationC(std::vector<std::size_t> const& indices) {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create rotation in plane a: Atom index out of range");
	}
	auto curr_rotations = Rotator::buildRotator(cartesians, indices)->makeRotations();

	rotatorContainer.emplace_back(curr_rotations.rotator);

	return std::move(curr_rotations.rotationC);
}

CommonInternalCoordinatesBuilder::TwoInternals CommonInternalCoordinatesBuilder::buildRotationAB(std::vector<std::size_t> const& indices) {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create rotation in plane a and b: Atom index out of range");
	}
	auto curr_rotations = Rotator::buildRotator(cartesians, indices)->makeRotations();

	rotatorContainer.emplace_back(curr_rotations.rotator);

	return { std::move(curr_rotations.rotationA), std::move(curr_rotations.rotationB) };
}

CommonInternalCoordinatesBuilder::TwoInternals CommonInternalCoordinatesBuilder::buildRotationAC(std::vector<std::size_t> const& indices) {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create rotation in plane a and c: Atom index out of range");
	}
	auto curr_rotations = Rotator::buildRotator(cartesians, indices)->makeRotations();

	rotatorContainer.emplace_back(curr_rotations.rotator);

	return { std::move(curr_rotations.rotationA), std::move(curr_rotations.rotationC) };
}

CommonInternalCoordinatesBuilder::TwoInternals CommonInternalCoordinatesBuilder::buildRotationBC(std::vector<std::size_t> const& indices) {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create rotation in plane b and c: Atom index out of range");
	}
	auto curr_rotations = Rotator::buildRotator(cartesians, indices)->makeRotations();

	rotatorContainer.emplace_back(curr_rotations.rotator);

	return { std::move(curr_rotations.rotationB), std::move(curr_rotations.rotationC) };
}

CommonInternalCoordinatesBuilder::ThreeInternals CommonInternalCoordinatesBuilder::buildRotationABC(std::vector<std::size_t> const& indices) {
	if (std::any_of(indices.cbegin(), indices.cend(), [this](auto const i) { return i >= this->numberOfAtoms; })) {
		throw std::runtime_error("Cannot create rotation in plane a, b, and c: Atom index out of range");
	}
	auto curr_rotations = Rotator::buildRotator(cartesians, indices)->makeRotations();

	rotatorContainer.emplace_back(curr_rotations.rotator);

	return { std::move(curr_rotations.rotationA), std::move(curr_rotations.rotationB), std::move(curr_rotations.rotationC) };
}

}