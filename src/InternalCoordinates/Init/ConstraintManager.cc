#include<algorithm>

#include "ConstraintManager.h"
#include "../InternalCoordinates/InternalCoordinate.h"
#include "../InternalCoordinates/Translations.h"

namespace internals {

std::shared_ptr<AbstractConstraint> ConstraintManager::findAndEraseInternalConstraint(internals::InternalCoordinate const& internalCoordinate) {
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

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForInternalCoordinate(internals::InternalCoordinate const& internalCoordinate, bool cosntrainAnyway) {
	auto constraint = findAndEraseInternalConstraint(internalCoordinate);
	if (constraint) return constraint;

	if (cosntrainAnyway) {
		masterConstraints->emplace_back(std::make_shared<BondConstraint>(internalCoordinate.getIndices(), true));
		return masterConstraints->back();
	}

	return constraint;
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForBonds(internals::InternalCoordinate const& internalCoordinate) {
	return checkForInternalCoordinate(internalCoordinate, constrainDistances);
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForAngles(internals::InternalCoordinate const& internalCoordinate) {
	return checkForInternalCoordinate(internalCoordinate, constrainAngles);
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForDihedrals(internals::InternalCoordinate const& internalCoordinate) {
	return checkForInternalCoordinate(internalCoordinate, constrainDihedrals);
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForTranslation(internals::Translations const& internalCoordinate) {
	using Direction = internals::Translations::Direction;

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

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForRotation(internals::Rotation const& internalCoordinate) {

	using Direction = internals::Rotation::Direction;

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

std::vector<std::unique_ptr< internals::InternalCoordinate>> ConstraintManager::restConstraintsToInternalCoordinates(internals::InternalCoordinatesBuilder & builder) const{
	std::vector<std::unique_ptr< internals::InternalCoordinate>> result;
	for (auto const& restConstraints : copiedConstraints) {
		result.emplace_back(restConstraints->makeInternal(builder));
		if (restConstraints->isFrozen()) result.back()->makeConstrained();
		else result.back()->releaseConstraint();
	}
	return result;
}

std::shared_ptr<AbstractConstraint> NoConstraintManager::checkForBonds(internals::InternalCoordinate const& /*internalCoordinate*/) { return std::shared_ptr<AbstractConstraint>{}; }
std::shared_ptr<AbstractConstraint> NoConstraintManager::checkForAngles(internals::InternalCoordinate const& /*internalCoordinate*/) { return std::shared_ptr<AbstractConstraint>{}; }
std::shared_ptr<AbstractConstraint> NoConstraintManager::checkForDihedrals(internals::InternalCoordinate const& /*internalCoordinate*/) { return std::shared_ptr<AbstractConstraint>{}; }
std::shared_ptr<AbstractConstraint> NoConstraintManager::checkForTranslation(internals::Translations const& /*internalCoordinate*/) { return std::shared_ptr<AbstractConstraint>{}; }
std::shared_ptr<AbstractConstraint> NoConstraintManager::checkForRotation(internals::Rotation const& /*internalCoordinate*/) { return std::shared_ptr<AbstractConstraint>{}; }
std::vector<std::unique_ptr<internals::InternalCoordinate>> NoConstraintManager::restConstraintsToInternalCoordinates(internals::InternalCoordinatesBuilder & /*builder*/) const {	return std::vector<std::unique_ptr<internals::InternalCoordinate>>(); }
	
} // namespace internals