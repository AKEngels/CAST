#ifndef CAST_INTERNALCOORDINATES_INIT_CONSTRAINTMANAGER_H_
#define CAST_INTERNALCOORDINATES_INIT_CONSTRAINTMANAGER_H_

#include<memory>

#include "Constraints.h"

class AbstractConstraintManager {
public:
	using ConstrainVec = std::vector<std::shared_ptr<AbstractConstraint>>;
	virtual std::shared_ptr<AbstractConstraint> checkIfConstraintPrimitive(std::vector<std::size_t> const&) = 0;
	virtual std::shared_ptr<AbstractConstraint> checkIfConstraintTrans(std::vector<std::size_t> const&, Constraint const) = 0;
	virtual std::shared_ptr<AbstractConstraint> checkIfConstraintRot(std::vector<std::size_t> const&, Constraint const) = 0;
	virtual ConstrainVec getConstraintsOfType(Constraint const) = 0;
	virtual ~AbstractConstraintManager() = 0;
};

inline AbstractConstraintManager::~AbstractConstraintManager() = default;

class NoConstraintManager : public AbstractConstraintManager {
public:
	virtual std::shared_ptr<AbstractConstraint> checkIfConstraintPrimitive(std::vector<std::size_t> const&) override;
	virtual std::shared_ptr<AbstractConstraint> checkIfConstraintTrans(std::vector<std::size_t> const&, Constraint const) override;
	virtual std::shared_ptr<AbstractConstraint> checkIfConstraintRot(std::vector<std::size_t> const&, Constraint const) override;
	virtual ConstrainVec getConstraintsOfType(Constraint const) override;
	//virtual ConstrainVec && getAllConstraints() override;
};

class ConstraintManager : public AbstractConstraintManager {
public:

	virtual std::shared_ptr<AbstractConstraint> checkIfConstraintPrimitive(std::vector<std::size_t> const&) override;
	virtual std::shared_ptr<AbstractConstraint> checkIfConstraintTrans(std::vector<std::size_t> const&, Constraint const) override;
	virtual std::shared_ptr<AbstractConstraint> checkIfConstraintRot(std::vector<std::size_t> const&, Constraint const) override;
	virtual ConstrainVec getConstraintsOfType(Constraint const) override;
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

	std::shared_ptr<AbstractConstraint> checkForBonds(std::vector<std::size_t> const&);
	std::shared_ptr<AbstractConstraint> checkForAngles(std::vector<std::size_t> const&);
	std::shared_ptr<AbstractConstraint> checkForDihedrals(std::vector<std::size_t> const&);

	ConstrainVec copiedConstraints;
};

#endif // CAST_INTERNALCOORDINATES_INIT_CONSTRAINTMANAGER_H_