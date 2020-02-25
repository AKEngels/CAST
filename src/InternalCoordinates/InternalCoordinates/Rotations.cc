#include"Rotations.h"

#include"../InternalCoordinateUtilities.h"

#include "../BestFitRotation.h"
#include "../../Scon/scon_mathmatrix.h"

#include"RotatorObserver.h"

#include<sstream>

namespace internals{


	Rotator::Rotator(coords::Representation_3D const& reference, std::vector<std::size_t> const& index_vec) :
		storedDerivativesForRotations{ std::make_unique<scon::mathmatrix<coords::float_type>>() },
		updateStoredValues{ true }, updateStoredDerivatives{ true }, reference_{ std::make_unique<coords::Representation_3D>(reference) }, rad_gyr_{ radiusOfGyration(*reference_) }{
		for (auto index : index_vec) {
			indices_.emplace_back(index - 1u);
		}
	}

std::array<float_type, 3u> const&
Rotator::valueOfInternalCoordinate(const coords::Representation_3D& new_xyz) {
	//TODO This needs to be activated again
	//if (!updateStoredValues) {
	//  return storedValuesForRotations;
	//}
	updateStoredValues = false;

	storedValuesForRotations = calculateValueOfInternalCoordinate(new_xyz);
	return storedValuesForRotations;
}

std::array<float_type, 3u>
Rotator::calculateValueOfInternalCoordinate(coords::Representation_3D const& newXyz) const {
	coords::Representation_3D curr_xyz_;
	curr_xyz_.reserve(indices_.size());
	for (auto const & i : indices_) {
		curr_xyz_.emplace_back(newXyz.at(i));
	}

	auto ret = ic_rotation::exponential_map<scon::mathmatrix>(*reference_, curr_xyz_);
	for (auto & r : ret) {
		r *= rad_gyr_;
	}

	return ret;
}

std::vector<scon::mathmatrix<float_type>>
Rotator::rot_der(coords::Representation_3D const& new_xyz) const {

	coords::Representation_3D new_xyz_;
	for (auto const & indi : indices_) {
		new_xyz_.emplace_back(new_xyz.at(indi));
	}

	return ic_rotation::exponential_derivs<scon::mathmatrix>(*reference_, new_xyz_);
}

scon::mathmatrix<float_type> const&
Rotator::rot_der_mat(coords::Representation_3D const& new_xyz) {
	using Mat = scon::mathmatrix<coords::float_type>;
	//TODO This needs to be activated again
	//if (!updateStoredDerivatives) {
	//  return *storedDerivativesForRotations;
	//}
	updateStoredDerivatives = false;

	auto const & zero = scon::mathmatrix<coords::float_type>::zero;

	auto first_ders = rot_der(new_xyz);

	Mat X = zero(new_xyz.size(), 3);
	Mat Y = zero(new_xyz.size(), 3);
	Mat Z = zero(new_xyz.size(), 3);

	for (auto i{ 0u }; i < first_ders.size(); ++i) {
		auto const& ind = indices_.at(i);
		auto const& der = first_ders.at(i);
		X.set_row(ind, der.col(0).t());
		Y.set_row(ind, der.col(1).t());
		Z.set_row(ind, der.col(2).t());
	}

	X *= rad_gyr_;
	Y *= rad_gyr_;
	Z *= rad_gyr_;

	//TODO initialize storedDerivativesForRotations in ctor
	*storedDerivativesForRotations = Mat(new_xyz.size() * 3, 3);
	storedDerivativesForRotations->set_col(0, X.vectorise_row());
	storedDerivativesForRotations->set_col(1, Y.vectorise_row());
	storedDerivativesForRotations->set_col(2, Z.vectorise_row());
	return *storedDerivativesForRotations;
}


Rotations Rotator::makeRotations() {
	return Rotations{
		shared_from_this(),
		std::make_unique<RotationA>(shared_from_this()),
		std::make_unique<RotationB>(shared_from_this()),
		std::make_unique<RotationC>(shared_from_this())
	};
}

bool Rotator::operator==(Rotator const & other) const {
	if (indices_.size() != other.indices_.size()) return false;
	for (auto i = 0u; i < indices_.size(); ++i) {
		if (indices_.at(i) != other.indices_.at(i)) return false;
	}
	return true;
}

coords::float_type
Rotator::radiusOfGyration(const coords::Representation_3D& struc) {
	return ic_util::rad_gyr(struc);
}

void Rotator::registerCartesians(InternalCoordinates::CartesiansForInternalCoordinates & cartesianCoordinates) {
	auto observer = std::make_shared<RotatorObserver>();
	observer->setNewRotator(shared_from_this());
	cartesianCoordinates.registerObserver(observer);
}

Rotator::~Rotator() = default;

float_type Rotation::val(scon::mathmatrix<float_type> const& cartesians) const {
	auto const& returnValues = rotator->valueOfInternalCoordinate(cartesians);
	return returnValues.at(index());
}

float_type Rotation::difference(scon::mathmatrix<float_type> const& newCoordinates, scon::mathmatrix<float_type> const& oldCoordinates) const {
	auto const previousVals = rotator->calculateValueOfInternalCoordinate(oldCoordinates);
	return val(newCoordinates) - previousVals.at(index());
}

scon::mathmatrix<float_type> Rotation::der_vec(scon::mathmatrix<float_type> const& cartesians) const {
	auto const& derivativeMatrix = rotator->rot_der_mat(cartesians);
	return derivativeMatrix.col_to_std_vector(index());
}

std::string Rotation::info(scon::mathmatrix<float_type> const & cartesians) const {
	std::ostringstream oss;
	oss << "Rotation " << name() << ": " << val(cartesians) << " | Constrained: " << std::boolalpha << is_constrained();
	return oss.str();
}

bool Rotation::hasIndices(std::vector<std::size_t> const& indices) const {
	return indices.size() == rotator->indices_.size() && ic_util::isSameSet(indices, rotator->indices_);
}

std::vector<std::size_t> Rotation::getIndices() const {
	return rotator->indices_;
}

void Rotation::makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) {
	auto constraint = manager->checkForRotation(*this);
	if (constraint) {
		if (constraint->isFrozen()) makeConstrained();
		else releaseConstraint();
	}
}


}