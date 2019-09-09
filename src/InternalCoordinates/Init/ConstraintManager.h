#ifndef CAST_INTERNALCOORDINATES_INIT_CONSTRAINTMANAGER_H_
#define CAST_INTERNALCOORDINATES_INIT_CONSTRAINTMANAGER_H_

#include<memory>

#include "Constraints.h"

namespace InternalCoordinates {
	class InternalCoordinate;
	class Translations;
	class Rotation;
}

class AbstractConstraintManager {
public:
	using ConstrainVec = std::vector<std::shared_ptr<AbstractConstraint>>;
	virtual std::shared_ptr<AbstractConstraint> checkForBonds(InternalCoordinates::InternalCoordinate const&) = 0;
	virtual std::shared_ptr<AbstractConstraint> checkForAngles(InternalCoordinates::InternalCoordinate const&) = 0;
	virtual std::shared_ptr<AbstractConstraint> checkForDihedrals(InternalCoordinates::InternalCoordinate const&) = 0;
	virtual std::shared_ptr<AbstractConstraint> checkForTranslation(InternalCoordinates::Translations const&) = 0;
	virtual std::shared_ptr<AbstractConstraint> checkForRotation(InternalCoordinates::Rotation const&) = 0;
	virtual std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> restConstraintsToInternalCoordinates(InternalCoordinates::InternalCoordinatesBuilder &) const = 0;
	virtual ~AbstractConstraintManager() = 0;
};

inline AbstractConstraintManager::~AbstractConstraintManager() = default;

class NoConstraintManager : public AbstractConstraintManager {
public:
	virtual std::shared_ptr<AbstractConstraint> checkForBonds(InternalCoordinates::InternalCoordinate const&) override;
	virtual std::shared_ptr<AbstractConstraint> checkForAngles(InternalCoordinates::InternalCoordinate const&) override;
	virtual std::shared_ptr<AbstractConstraint> checkForDihedrals(InternalCoordinates::InternalCoordinate const&) override;
	virtual std::shared_ptr<AbstractConstraint> checkForTranslation(InternalCoordinates::Translations const&) override;
	virtual std::shared_ptr<AbstractConstraint> checkForRotation(InternalCoordinates::Rotation const&) override;
	virtual std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> restConstraintsToInternalCoordinates(InternalCoordinates::InternalCoordinatesBuilder &) const override;
	//virtual ConstrainVec && getAllConstraints() override;
};

class ConstraintManager : public AbstractConstraintManager {
public:

	virtual std::shared_ptr<AbstractConstraint> checkForBonds(InternalCoordinates::InternalCoordinate const&) override;
	virtual std::shared_ptr<AbstractConstraint> checkForAngles(InternalCoordinates::InternalCoordinate const&) override;
	virtual std::shared_ptr<AbstractConstraint> checkForDihedrals(InternalCoordinates::InternalCoordinate const&) override;
	virtual std::shared_ptr<AbstractConstraint> checkForTranslation(InternalCoordinates::Translations const&) override;
	virtual std::shared_ptr<AbstractConstraint> checkForRotation(InternalCoordinates::Rotation const&) override;
	virtual std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> restConstraintsToInternalCoordinates(InternalCoordinates::InternalCoordinatesBuilder &) const override;
	//virtual ConstrainVec && getAllConstraints() override;

	ConstraintManager(std::shared_ptr<ConstrainVec> const& constraints) : masterConstraints(constraints), constrainDistances{ false }, constrainAngles{ false }, constrainDihedrals{ false }, constrainOOP{ false }, constrainTranslations{ false }, constrainRotations{ false }, copiedConstraints{ getSharedConstraints() }{};

	ConstraintManager& constrainAllDistances(bool const constrain) { constrainDistances = constrain; return *this; }
	ConstraintManager& constrainAllAngles(bool const constrain) { constrainAngles = constrain; return *this; }
	ConstraintManager& constrainAllDihedrals(bool const constrain) { constrainDihedrals = constrain; return *this; }
	ConstraintManager& constrainAllOOPs(bool const constrain) { constrainOOP = constrain; return *this; }
	ConstraintManager& constrainAllTranslations(bool const constrain) { constrainTranslations = constrain; return *this; }
	ConstraintManager& constrainAllRotations(bool const constrain) { constrainRotations = constrain; return *this; }

	ConstrainVec getSharedConstraints() const { return *masterConstraints; }

private:
	std::shared_ptr<ConstrainVec> masterConstraints;
	bool constrainDistances, constrainAngles, constrainDihedrals, constrainOOP, constrainTranslations, constrainRotations;

	std::shared_ptr<AbstractConstraint> findAndEraseInternalConstraint(InternalCoordinates::InternalCoordinate const&);
	std::shared_ptr<AbstractConstraint> checkForInternalCoordinate(InternalCoordinates::InternalCoordinate const& internalCoordinate, bool cosntrainAnyway);

	ConstrainVec copiedConstraints;
};

#endif // CAST_INTERNALCOORDINATES_INIT_CONSTRAINTMANAGER_H_