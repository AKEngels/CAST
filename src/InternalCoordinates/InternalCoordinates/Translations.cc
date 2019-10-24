#include"Translations.h"

#include "../InternalCoordinateUtilities.h"

#include "../../coords_rep.h"

namespace internals{

float_type Translations::val(coords::Representation_3D const &cartesians) const {
	auto coord_sum{ 0.0 };
	for (auto &i : indices_) {
		coord_sum += coord_func(cartesians.at(i));
	}
	return coord_sum / indices_.size();
}

float_type Translations::difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const {
	return val(newCoordinates) - val(oldCoordinates);
}

std::string Translations::info(coords::Representation_3D const & cartesians) const
{
	std::ostringstream oss;
	oss << "Trans " << coordinate_letter() << ": " << val(cartesians) << " | Constrained: " << std::boolalpha << is_constrained();
	return oss.str();
}

std::vector<float_type>
Translations::der_vec(coords::Representation_3D const& cartesians) const {
	coords::Representation_3D result(cartesians.size(), coords::Cartesian_Point(0., 0., 0.));

	for (auto const& i : indices_) {
		result.at(i) = size_reciprocal();
	}

	return ic_util::flatten_c3_vec(result);
}

bool Translations::hasIndices(std::vector<std::size_t> const& indices) const {
	return indices.size() == indices_.size() && ic_util::isSameSet(indices, indices_);
}

std::vector<std::size_t> Translations::getIndices() const {
	return indices_;
}

void Translations::makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) {
	auto constraint = manager->checkForTranslation(*this);
	if (constraint) {
		if (constraint->isFrozen()) makeConstrained();
		else releaseConstraint();
	}
}

coords::Cartesian_Point TranslationX::size_reciprocal() const
{
	return coords::Cartesian_Point{ 1. / float_type(indices_.size()), 0., 0. };
}

coords::Cartesian_Point TranslationY::size_reciprocal() const
{
	return coords::Cartesian_Point{ 0., 1. / coords::float_type(indices_.size()), 0. };
}

coords::Cartesian_Point TranslationZ::size_reciprocal() const
{
	return coords::Cartesian_Point{ 0., 0., 1. / coords::float_type(indices_.size()) };
}

bool Translations::operator==(Translations const & other) const {
	if (indices_.size() != other.indices_.size()) return false;
	for (auto i = 0u; i < indices_.size(); ++i) {
		if (indices_.at(i) != other.indices_.at(i)) return false;
	}
	return true;
}

}