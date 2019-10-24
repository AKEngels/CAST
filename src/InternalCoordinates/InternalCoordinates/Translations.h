#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_TRANSLATIONS_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_TRANSLATIONS_H_

#include "InternalCoordinate.h"

namespace internals{

struct Translations : public InternalCoordinate {

	enum class Direction : int { X, Y, Z };


	virtual ~Translations() = default;

	virtual float_type val(coords::Representation_3D const& cartesians) const override;
	float_type difference(coords::Representation_3D const& newCoordinates, coords::Representation_3D const& oldCoordinates) const override;
	virtual std::string info(coords::Representation_3D const& cartesians) const override;
	virtual std::vector<float_type> der_vec(coords::Representation_3D const& cartesians) const override;

	std::vector<std::size_t> indices_;

	float_type hessian_guess(coords::Representation_3D const& /*cartesians*/) const override {
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
	virtual float_type coord_func(coords::Cartesian_Point const& cp) const = 0;
	virtual char coordinate_letter() const = 0;
	virtual coords::Cartesian_Point size_reciprocal() const = 0;
};

struct TranslationX : Translations {
	explicit TranslationX(std::vector<std::size_t> const& index_vec) :
		Translations(index_vec)
	{}

protected:
	float_type coord_func(coords::Cartesian_Point const& cp) const override {
		return cp.x();
	}

	Direction getDirection() const override { return Direction::X; }

	char coordinate_letter() const override { return 'X'; }

	coords::Cartesian_Point size_reciprocal() const override;
};

struct TranslationY : Translations {
	explicit TranslationY(std::vector<std::size_t> const& index_vec) :
		Translations(index_vec)
	{}

protected:
	coords::float_type coord_func(coords::Cartesian_Point const& cp) const override {
		return cp.y();
	}

	Direction getDirection() const override { return Direction::Y; }

	char coordinate_letter() const override { return 'Y'; }

	coords::Cartesian_Point size_reciprocal() const override;
};

struct TranslationZ : Translations {
	explicit TranslationZ(std::vector<std::size_t> const& index_vec) :
		Translations(index_vec)
	{}
protected:
	coords::float_type coord_func(coords::Cartesian_Point const& cp) const override {
		return cp.z();
	}

	Direction getDirection() const override { return Direction::Z; }

	char coordinate_letter() const override { return 'Z'; }

	coords::Cartesian_Point size_reciprocal() const override;
};

}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_TRANSLATIONS_H_