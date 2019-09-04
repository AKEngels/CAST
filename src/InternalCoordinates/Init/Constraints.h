#ifndef CAST_INTERNALCOORDINATES_INIT_CONSTRAINTS_H_
#define CAST_INTERNALCOORDINATES_INIT_CONSTRAINTS_H_

#include<vector>

/* Constraints on internal coordinates
	 */

	 // Proposal: give all Constraints a shared pointer to the Internal Coorindates in order, to controll the constraint coordinates.
	 //class ::InternalCoordinates::InternalCoordinate;

enum class Constraint : int { NONE = 0, BONDS, ANGLE, DIHEDRAL, TRANSLATION_X, TRANSLATION_Y, TRANSLATION_Z, ROTATION_A, ROTATION_B, ROTATION_C };

struct AbstractConstraint {
	virtual Constraint getType() const = 0;
	virtual bool isFrozen() const = 0;
	virtual std::vector<std::size_t> const& getAtomIndices() const = 0;
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
};

struct AngleConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::ANGLE;
	}
};

struct DihedralConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::DIHEDRAL;
	}
};

struct TranslationXConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::TRANSLATION_X;
	}
};

struct TranslationYConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::TRANSLATION_Y;
	}
};

struct TranslationZConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::TRANSLATION_Z;
	}
};

struct RotationAConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::ROTATION_A;
	}
};

struct RotationBConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::ROTATION_B;
	}
};

struct RotationCConstraint : NormalConstraint {
	using NormalConstraint::NormalConstraint;
	virtual Constraint getType() const override {
		return Constraint::ROTATION_C;
	}
};

#endif // CAST_INTERNALCOORDINATES_INIT_CONSTRAINTS_H_