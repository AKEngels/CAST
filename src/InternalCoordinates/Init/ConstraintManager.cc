#include<algorithm>

#include "ConstraintManager.h"
#include "../InternalCoordinates.h"

std::shared_ptr<AbstractConstraint> ConstraintManager::findAndEraseInternalConstraint(InternalCoordinates::InternalCoordinate const& internalCoordinate) {
	for (auto it = copiedConstraints.begin(); it != copiedConstraints.end(); ++it) {
		auto const& constraint = (*it);
		if (internalCoordinate.hasIndices(constraint->getAtomIndices())) {
			auto ret = *it;
			copiedConstraints.erase(it);
			return ret;
		}
	}
	return std::shared_ptr<AbstractConstraint>();
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForInternalCoordinate(InternalCoordinates::InternalCoordinate const& internalCoordinate, bool cosntrainAnyway) {
	auto constraint = findAndEraseInternalConstraint(internalCoordinate);
	if (constraint) return constraint;

	if (cosntrainAnyway) {
		masterConstraints->emplace_back(std::make_shared<BondConstraint>(internalCoordinate.getIndices(), true));
		return masterConstraints->back();
	}

	return constraint;
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForBonds(InternalCoordinates::InternalCoordinate const& internalCoordinate) {
	return checkForInternalCoordinate(internalCoordinate, constrainDistances);
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForAngles(InternalCoordinates::InternalCoordinate const& internalCoordinate) {
	return checkForInternalCoordinate(internalCoordinate, constrainAngles);
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForDihedrals(InternalCoordinates::InternalCoordinate const& internalCoordinate) {
	return checkForInternalCoordinate(internalCoordinate, constrainDihedrals);
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForTranslation(InternalCoordinates::Translations const& internalCoordinate) {
	using Direction = InternalCoordinates::Translations::Direction;

	auto constraint = findAndEraseInternalConstraint(internalCoordinate);
	if (constraint) return constraint;


	if (constrainTranslations) {
		if (internalCoordinate.getDirection() == Direction::X) {
			masterConstraints->emplace_back(std::make_shared<TranslationXConstraint>(internalCoordinate.getIndices(), true));
		}
		else if (internalCoordinate.getDirection() == Direction::Y) {
			masterConstraints->emplace_back(std::make_shared<TranslationYConstraint>(internalCoordinate.getIndices(), true));
		}
		else if (internalCoordinate.getDirection() == Direction::Z) {
			masterConstraints->emplace_back(std::make_shared<TranslationZConstraint>(internalCoordinate.getIndices(), true));
		}

		return masterConstraints->back();
	}


	return constraint;
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForRotation(InternalCoordinates::Rotation const& internalCoordinate) {

	using Direction = InternalCoordinates::Rotation::Direction;

	auto constraint = findAndEraseInternalConstraint(internalCoordinate);
	if (constraint) return constraint;

	if (constrainTranslations) {
		if (internalCoordinate.getDirection() == Direction::A) {
			masterConstraints->emplace_back(std::make_shared<RotationAConstraint>(internalCoordinate.getIndices(), true));
		}
		else if (internalCoordinate.getDirection() == Direction::B) {
			masterConstraints->emplace_back(std::make_shared<RotationBConstraint>(internalCoordinate.getIndices(), true));
		}
		else if (internalCoordinate.getDirection() == Direction::C) {
			masterConstraints->emplace_back(std::make_shared<RotationCConstraint>(internalCoordinate.getIndices(), true));
		}

		return masterConstraints->back();
	}


	return constraint;

}

std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> ConstraintManager::restConstraintsToInternalCoordinates(InternalCoordinates::InternalCoordinatesBuilder & builder) const{
	std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
	for (auto const& restConstraints : copiedConstraints) {
		result.emplace_back(restConstraints->makeInternal(builder));
		if (restConstraints->isFrozen()) result.back()->makeConstrained();
		else result.back()->releaseConstraint();
	}
	return result;
}

std::shared_ptr<AbstractConstraint> NoConstraintManager::checkForBonds(InternalCoordinates::InternalCoordinate const& /*internalCoordinate*/) { return std::shared_ptr<AbstractConstraint>{}; }
std::shared_ptr<AbstractConstraint> NoConstraintManager::checkForAngles(InternalCoordinates::InternalCoordinate const& /*internalCoordinate*/) { return std::shared_ptr<AbstractConstraint>{}; }
std::shared_ptr<AbstractConstraint> NoConstraintManager::checkForDihedrals(InternalCoordinates::InternalCoordinate const& /*internalCoordinate*/) { return std::shared_ptr<AbstractConstraint>{}; }
std::shared_ptr<AbstractConstraint> NoConstraintManager::checkForTranslation(InternalCoordinates::Translations const& /*internalCoordinate*/) { return std::shared_ptr<AbstractConstraint>{}; }
std::shared_ptr<AbstractConstraint> NoConstraintManager::checkForRotation(InternalCoordinates::Rotation const& /*internalCoordinate*/) { return std::shared_ptr<AbstractConstraint>{}; }
std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> NoConstraintManager::restConstraintsToInternalCoordinates(InternalCoordinates::InternalCoordinatesBuilder & /*builder*/) const {	return std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>(); }