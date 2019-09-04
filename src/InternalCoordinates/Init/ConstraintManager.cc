#include<algorithm>

#include "ConstraintManager.h"

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForBonds(std::vector<std::size_t> const& constraints) {
	for (auto it = copiedConstraints.begin(); it != copiedConstraints.end(); ++it) {
		auto const& constraint = (*it);
		if (constraint->getAtomIndices().size() != constraints.size()) continue;
		if ((constraint->getAtomIndices()[0] == constraints[0] + 1 && constraint->getAtomIndices()[1] == constraints[1] + 1) ||
			(constraint->getAtomIndices()[0] == constraints[1] + 1 && constraint->getAtomIndices()[1] == constraints[0] + 1)) {
			auto ret = *it;
			copiedConstraints.erase(it);
			return ret;
		}
	}
	if (constrainDistances) {
		masterConstraints->emplace_back(std::make_shared<BondConstraint>(std::vector<std::size_t>(constraints), true));
		return masterConstraints->back();
	}

	return std::shared_ptr<AbstractConstraint>();
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForAngles(std::vector<std::size_t> const& constraints) {
	for (auto it = copiedConstraints.begin(); it != copiedConstraints.end(); ++it) {
		auto const& constraint = (*it);
		if (constraint->getAtomIndices().size() != constraints.size()) continue;
		if ((constraint->getAtomIndices()[0] == constraints[0] + 1 && constraint->getAtomIndices()[1] == constraints[1] + 1 && constraint->getAtomIndices()[2] == constraints[2] + 1) ||
			(constraint->getAtomIndices()[0] == constraints[2] + 1 && constraint->getAtomIndices()[1] == constraints[1] + 1 && constraint->getAtomIndices()[2] == constraints[0] + 1)) {
			auto ret = *it;
			copiedConstraints.erase(it);
			return ret;
		}
	}

	if (constrainAngles) {
		masterConstraints->emplace_back(std::make_shared<AngleConstraint>(std::vector<std::size_t>(constraints), true));
		return masterConstraints->back();
	}

	return std::shared_ptr<AbstractConstraint>();
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkForDihedrals(std::vector<std::size_t> const& constraints) {
	for (auto it = copiedConstraints.begin(); it != copiedConstraints.end(); ++it) {
		auto const& constraint = (*it);
		if (constraint->getAtomIndices().size() != constraints.size()) continue;
		if ((constraint->getAtomIndices()[0] == constraints[0] + 1 && constraint->getAtomIndices()[1] == constraints[1] + 1 && constraint->getAtomIndices()[2] == constraints[2] + 1 && constraint->getAtomIndices()[3] == constraints[3] + 1) ||
			(constraint->getAtomIndices()[3] == constraints[0] + 1 && constraint->getAtomIndices()[2] == constraints[1] + 1 && constraint->getAtomIndices()[1] == constraints[2] + 1 && constraint->getAtomIndices()[0] == constraints[3] + 1)) {
			auto ret = *it;
			copiedConstraints.erase(it);
			return ret;
		}
	}

	if (constrainDihedrals) {
		masterConstraints->emplace_back(std::make_shared<DihedralConstraint>(std::vector<std::size_t>(constraints), true));
		return masterConstraints->back();
	}

	return std::shared_ptr<AbstractConstraint>();
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkIfConstraintPrimitive(std::vector<std::size_t> const& constraints) {
	if (constraints.size() == 2) {
		return checkForBonds(constraints);
	}
	else if (constraints.size() == 3) {
		return checkForAngles(constraints);
	}
	else if (constraints.size() == 4) {
		return checkForDihedrals(constraints);
	}
	else throw std::runtime_error("Unexpected number of constraints!");
}


bool isSameSet(std::vector<std::size_t> lhs, std::vector<std::size_t> rhs) {
	if (lhs.size() != rhs.size()) return false;
	std::sort(lhs.begin(), lhs.end());
	std::sort(rhs.begin(), rhs.end());
	return std::includes(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkIfConstraintTrans(std::vector<std::size_t> const& constraints, Constraint type) {
	for (auto it = copiedConstraints.begin(); it != copiedConstraints.end(); ++it) {
		auto const& constraint = *it;
		if (constraint->getType() == type && isSameSet(constraint->getAtomIndices(), constraints)) {
			auto ret = *it;
			copiedConstraints.erase(it);
			return ret;
		}
	}

	if (constrainTranslations) {
		if (type == Constraint::TRANSLATION_X) {
			masterConstraints->emplace_back(std::make_shared<TranslationXConstraint>(std::vector<std::size_t>(constraints), true));
		}
		else if (type == Constraint::TRANSLATION_Y) {
			masterConstraints->emplace_back(std::make_shared<TranslationYConstraint>(std::vector<std::size_t>(constraints), true));
		}
		else if (type == Constraint::TRANSLATION_Z) {
			masterConstraints->emplace_back(std::make_shared<TranslationZConstraint>(std::vector<std::size_t>(constraints), true));
		}

		return masterConstraints->back();
	}

	return std::shared_ptr<AbstractConstraint>();
}

std::shared_ptr<AbstractConstraint> ConstraintManager::checkIfConstraintRot(std::vector<std::size_t> const& constraints, Constraint type) {

	std::vector<std::size_t> newConstraints;
	newConstraints.reserve(constraints.size());
	std::transform(constraints.cbegin(), constraints.cend(), std::back_inserter(newConstraints), [](std::size_t s) {return ++s; });

	for (auto it = copiedConstraints.begin(); it != copiedConstraints.end(); ++it) {
		auto const& constraint = *it;
		if (constraint->getType() == type && isSameSet(constraint->getAtomIndices(), newConstraints)) {
			auto ret = *it;
			copiedConstraints.erase(it);
			return ret;
		}
	}

	if (constrainTranslations) {
		if (type == Constraint::ROTATION_A) {
			masterConstraints->emplace_back(std::make_shared<RotationAConstraint>(std::vector<std::size_t>(constraints), true));
		}
		else if (type == Constraint::ROTATION_B) {
			masterConstraints->emplace_back(std::make_shared<RotationBConstraint>(std::vector<std::size_t>(constraints), true));
		}
		else if (type == Constraint::ROTATION_C) {
			masterConstraints->emplace_back(std::make_shared<RotationCConstraint>(std::vector<std::size_t>(constraints), true));
		}

		return masterConstraints->back();
	}

	return std::shared_ptr<AbstractConstraint>();

}

ConstraintManager::ConstrainVec ConstraintManager::getConstraintsOfType(Constraint const type) {
	auto it = std::stable_partition(copiedConstraints.begin(), copiedConstraints.end(), [type](auto const& c) {
		return c->getType() != type;
	});
	ConstrainVec y(std::make_move_iterator(it), std::make_move_iterator(copiedConstraints.end()));
	copiedConstraints.erase(it, copiedConstraints.end());

	return y;
}


/*ConstraintManager::ConstrainVec && ConstraintManager::getAllConstraints() {
	return std::move(copiedConstraints);
}*/

std::shared_ptr<AbstractConstraint> NoConstraintManager::checkIfConstraintPrimitive(std::vector<std::size_t> const&) { return std::shared_ptr<AbstractConstraint>{}; }
std::shared_ptr<AbstractConstraint> NoConstraintManager::checkIfConstraintTrans(std::vector<std::size_t> const&, Constraint const) { return std::shared_ptr<AbstractConstraint>{}; }
std::shared_ptr<AbstractConstraint> NoConstraintManager::checkIfConstraintRot(std::vector<std::size_t> const&, Constraint const) { return std::shared_ptr<AbstractConstraint>{}; }
NoConstraintManager::ConstrainVec NoConstraintManager::getConstraintsOfType(Constraint const) { return ConstrainVec{}; }