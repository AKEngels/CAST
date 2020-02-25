#include"Translations.h"

#include "../InternalCoordinateUtilities.h"
#include "../../Scon/scon_mathmatrix.h"

namespace internals{

float_type Translations::val(scon::mathmatrix<float_type> const &cartesians) const {
	auto coord_sum{ 0.0 };
	for (auto &i : indices_) {
		coord_sum += coord_func(getAtom(cartesians, i));
	}
	return coord_sum / indices_.size();
}

float_type Translations::difference(scon::mathmatrix<float_type> const& newCoordinates, scon::mathmatrix<float_type> const& oldCoordinates) const {
	return val(newCoordinates) - val(oldCoordinates);
}

std::string Translations::info(scon::mathmatrix<float_type> const & cartesians) const
{
	std::ostringstream oss;
	oss << "Trans " << coordinate_letter() << ": " << val(cartesians) << " | Constrained: " << std::boolalpha << is_constrained();
	return oss.str();
}

scon::mathmatrix<float_type>
Translations::der_vec(scon::mathmatrix<float_type> const& cartesians) const {
	auto firstder = sizeReciprocal();

	BmatrixRowCreator rowCreator(cartesians.cols()*cartesians.rows());

	for (auto const& i : indices_) {
		rowCreator.insertAtomDerivative(firstder, i);
	}

	// TODO: test if this line is mandatory to trigger return optimization
	scon::mathmatrix<float_type> result = rowCreator.getRow();

	return result;
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

CartesianPoint TranslationX::sizeReciprocal() const
{
	return CartesianPoint{ 1. / float_type(indices_.size()), 0., 0. };
}

CartesianPoint TranslationY::sizeReciprocal() const
{
	return CartesianPoint{ 0., 1. / float_type(indices_.size()), 0. };
}

CartesianPoint TranslationZ::sizeReciprocal() const
{
	return CartesianPoint{ 0., 0., 1. / float_type(indices_.size()) };
}

bool Translations::operator==(Translations const & other) const {
	if (indices_.size() != other.indices_.size()) return false;
	for (auto i = 0u; i < indices_.size(); ++i) {
		if (indices_.at(i) != other.indices_.at(i)) return false;
	}
	return true;
}

}