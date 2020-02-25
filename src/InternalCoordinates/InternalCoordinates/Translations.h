#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_TRANSLATIONS_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_TRANSLATIONS_H_

#include "InternalCoordinate.h"
#include "Helpers.h"

namespace internals{

struct Translations : public InternalCoordinate {

	enum class Direction : int { X, Y, Z };


	virtual ~Translations() = default;

	virtual float_type val(scon::mathmatrix<float_type> const& cartesians) const override;
	virtual float_type difference(scon::mathmatrix<float_type> const& newCoordinates, scon::mathmatrix<float_type> const& oldCoordinates) const override;
	virtual std::string info(scon::mathmatrix<float_type> const& cartesians) const override;
	virtual scon::mathmatrix<float_type> der_vec(scon::mathmatrix<float_type> const& cartesians) const override;

	std::vector<std::size_t> indices_;

	virtual float_type hessian_guess(scon::mathmatrix<float_type> const& /*cartesians*/) const override {
		return 0.05;
	}

	bool operator==(Translations const&) const;

	virtual bool hasIndices(std::vector<std::size_t> const& indices) const override;
	virtual std::vector<std::size_t> getIndices() const override;
	virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) override;
	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

	bool constrained_;
	virtual bool is_constrained() const override { return constrained_; }

	virtual Direction getDirection() const = 0;

protected:
	Translations(std::vector<std::size_t> const& index_vec) :
		constrained_{ false }
	{
		for (auto index : index_vec) {
			indices_.emplace_back(index - 1u);
		}
	}
	virtual float_type coord_func(CartesianPoint const& cp) const = 0;
	virtual char coordinate_letter() const = 0;
	virtual CartesianPoint sizeReciprocal() const = 0;
};

struct TranslationX : Translations {
	explicit TranslationX(std::vector<std::size_t> const& index_vec) :
		Translations(index_vec)
	{}

protected:
	float_type coord_func(CartesianPoint const& cp) const override {
		return std::get<0u>(cp);
	}

	Direction getDirection() const override { return Direction::X; }

	char coordinate_letter() const override { return 'X'; }

	virtual CartesianPoint sizeReciprocal() const override;
};

struct TranslationY : Translations {
	explicit TranslationY(std::vector<std::size_t> const& index_vec) :
		Translations(index_vec)
	{}

protected:
	coords::float_type coord_func(CartesianPoint const& cp) const override {
		return std::get<1u>(cp);
	}

	Direction getDirection() const override { return Direction::Y; }

	char coordinate_letter() const override { return 'Y'; }

	CartesianPoint sizeReciprocal() const override;
};

struct TranslationZ : Translations {
	explicit TranslationZ(std::vector<std::size_t> const& index_vec) :
		Translations(index_vec)
	{}
protected:
	coords::float_type coord_func(CartesianPoint const& cp) const override {
		return std::get<2u>(cp);
	}

	Direction getDirection() const override { return Direction::Z; }

	char coordinate_letter() const override { return 'Z'; }

	virtual CartesianPoint sizeReciprocal() const override;
};

}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_TRANSLATIONS_H_