#ifndef CAST_INTERNALCOORDINATES_INIT_CONSTRAINTS_H_
#define CAST_INTERNALCOORDINATES_INIT_CONSTRAINTS_H_

#include<vector>
#include<memory>

#include "../InternalCoordinatesAliases.h"

/* Constraints on internal coordinates
	 */

	 // Proposal: give all Constraints a shared pointer to the Internal Coorindates in order, to controll the constraint coordinates.
	 //class ::InternalCoordinates::InternalCoordinate;

enum class Constraint : int { NONE = 0, BONDS, ANGLE, DIHEDRAL, TRANSLATION_X, TRANSLATION_Y, TRANSLATION_Z, ROTATION_A, ROTATION_B, ROTATION_C };

struct AbstractConstraint {
	virtual Constraint getType() const = 0;
	virtual bool isFrozen() const = 0;
	virtual std::vector<std::size_t> const& getAtomIndices() const = 0;
	virtual std::unique_ptr<InternalCoordinates::InternalCoordinate> makeInternal(InternalCoordinates::InternalCoordinatesBuilder &) const = 0;
};

struct NormalConstraint : AbstractConstraint {
	NormalConstraint(std::vector<std::size_t>&& indices, bool freeze) : atomIndices(std::move(indices)), isFreeze(freeze) {}
	virtual std::vector<std::size_t> const& getAtomIndices() const override {
		return atomIndices;
	}
	virtual bool isFrozen() const override { return isFreeze; }
	std::vector<std::size_t> atomIndices;
	bool isFreeze;
};

struct BondConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::BONDS;
	}
	std::unique_ptr<InternalCoordinates::InternalCoordinate> makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const override;
};

struct AngleConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::ANGLE;
	}
	std::unique_ptr<InternalCoordinates::InternalCoordinate> makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const override;
};

struct DihedralConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::DIHEDRAL;
	}
	std::unique_ptr<InternalCoordinates::InternalCoordinate> makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const override;
};

struct TranslationXConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::TRANSLATION_X;
	}
	std::unique_ptr<InternalCoordinates::InternalCoordinate> makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const override;
};

struct TranslationYConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::TRANSLATION_Y;
	}
	std::unique_ptr<InternalCoordinates::InternalCoordinate> makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const override;
};

struct TranslationZConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::TRANSLATION_Z;
	}
	std::unique_ptr<InternalCoordinates::InternalCoordinate> makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const override;
};

struct RotationAConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::ROTATION_A;
	}
	std::unique_ptr<InternalCoordinates::InternalCoordinate> makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const override;
};

struct RotationBConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::ROTATION_B;
	}
	std::unique_ptr<InternalCoordinates::InternalCoordinate> makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const override;
};

struct RotationCConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::ROTATION_C;
	}
	std::unique_ptr<InternalCoordinates::InternalCoordinate> makeInternal(InternalCoordinates::InternalCoordinatesBuilder & builder) const override;
};

#endif // CAST_INTERNALCOORDINATES_INIT_CONSTRAINTS_H_