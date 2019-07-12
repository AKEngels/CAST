#include "InternalCoordinateDecorator.h"

#include "PrimitiveInternalCoordinates.h"

namespace internals{

  ICAbstractDecorator::ICAbstractDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    parent_{ parent }
  {}
  
  void ICAbstractDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager){
    parent_->buildCoordinates(cartesians, graph, indexVec, manager);
  }
  
  void ICAbstractDecorator::appendCoordinates(std::shared_ptr<InternalCoordinateAppenderInterface> appender){
    parent_->appendCoordinates(appender);
  }
  
  ICAbstractDecorator::InternalCoordinatesCreator::InternalCoordinatesCreator(BondGraph const& graph):
    bondGraph(graph)
  {}
  
  ICGeneralAppender::ICGeneralAppender(InternalVec && internal_coords):
    internal_coords_{ std::move (internal_coords) }
  {}
  
  void ICGeneralAppender::append(std::shared_ptr<PrimitiveInternalCoordinates> primitives){
    primitives->appendPrimitives(std::move(internal_coords_));
  }
  
  ICBondDecorator::ICBondDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    ICAbstractDecorator{ parent }
  {}
    
  void ICBondDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager){
    DistanceCreator dc(graph);
    appendCoordinates(std::make_shared<ICGeneralAppender>(dc.getInternals(manager)));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec, manager);
  }
  
  ICAbstractDecorator::DistanceCreator::DistanceCreator(BondGraph const& graph):
     InternalCoordinatesCreator{ graph },
     source{ 0u },
     target{ 0u },
     edgeIterators{ boost::edges(bondGraph) },
	  pointerToResult{ nullptr }
  {}
  
  InternalVec ICAbstractDecorator::DistanceCreator::getInternals(AbstractConstraintManager& manager) {
    InternalVec result;
	pointerToResult = &result;
    while (nextEdgeDistances()) {
		auto bond = std::make_unique<InternalCoordinates::BondDistance>(bondGraph[source], bondGraph[target]);

		auto constrtaint = manager.checkIfConstraintPrimitive({ bond->index_a_, bond->index_b_ });
		if (constrtaint) {
			if (constrtaint->isFrozen()) bond->makeConstrained();
			else bond->releaseConstraint();
		}
      result.emplace_back(std::move(bond));
    }
    // Add new bonds if we still have constraints left
	auto num_atoms = boost::num_vertices(bondGraph);
	auto rest = manager.getConstraintsOfType(config::Constraint::BONDS);
    for (auto const& curr_constraint : rest){
        auto index1 = curr_constraint->getAtomIndices()[0]-1, index2 = curr_constraint->getAtomIndices()[1]-1;
        if (index1 >= num_atoms || index2 >= num_atoms){
          throw std::runtime_error("Cannot create constrained bond length coordinate: Atom index out of range");
        }
        result.emplace_back(std::make_unique<InternalCoordinates::BondDistance>(bondGraph[index1], bondGraph[index2]));
		if(curr_constraint->isFrozen()) result.back()->makeConstrained();
		else result.back()->releaseConstraint();
    }
    return result;
  }

  bool ICAbstractDecorator::DistanceCreator::nextEdgeDistances() {
    if (edgeIterators.first == edgeIterators.second) return false;
    source = boost::source(*edgeIterators.first, bondGraph);
    target = boost::target(*edgeIterators.first, bondGraph);
    ++edgeIterators.first;
    return true;
  }
  
  ICAngleDecorator::ICAngleDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    ICAbstractDecorator{ parent }
  {}
    
  void ICAngleDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager){
    AngleCreator ac(graph);
    appendCoordinates(std::make_shared<ICGeneralAppender>(ac.getInternals(manager)));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec, manager);
  }
  
  ICAbstractDecorator::AngleCreator::AngleCreator(BondGraph const& graph):
    InternalCoordinatesCreator{ graph },
    leftAtom{ 0u },
    middleAtom{ 0u },
    rightAtom{ 0u },
    vertexIterators{ boost::vertices(graph) },
	  pointerToResult{ nullptr },
	  pointerToManager{ nullptr }
  {}
  
  InternalVec ICAbstractDecorator::AngleCreator::getInternals(AbstractConstraintManager& manager) {
    InternalVec result;
	pointerToManager = &manager;
    pointerToResult = &result;
    while (nextVertex()) {
      addAngleForAllNeighbors();
    }
    // Add new angles if we still have constraints left
	auto num_atoms = boost::num_vertices(bondGraph);
    for (auto const& curr_constraint : manager.getConstraintsOfType(config::Constraint::ANGLE)){
		auto index1 = curr_constraint->getAtomIndices()[0] - 1, index2 = curr_constraint->getAtomIndices()[1] - 1, index3 = curr_constraint->getAtomIndices()[2] - 1;
        if (index1 >= num_atoms || index2 >= num_atoms || index3 >= num_atoms){
          throw std::runtime_error("Cannot create constrained bond angle coordinate: Atom index out of range");
        }
		pointerToResult->emplace_back(std::make_unique<InternalCoordinates::BondAngle>(bondGraph[index1], bondGraph[index2], bondGraph[index3]));
		if (curr_constraint->isFrozen()) pointerToResult->back()->makeConstrained();
		else pointerToResult->back()->releaseConstraint();
      
    }
	pointerToResult = nullptr;
	pointerToManager = nullptr;
    return result;
  }

  bool ICAbstractDecorator::AngleCreator::nextVertex() {
    if (vertexIterators.first == vertexIterators.second) return false;
    middleAtom = *vertexIterators.first;
    vertexIterators.first++;
    return true;
  }

  void ICAbstractDecorator::AngleCreator::addAngleForAllNeighbors() {
    auto allNeighbors = boost::adjacent_vertices(middleAtom, bondGraph);
    spanLeftAndRightNeighborsForAngle(allNeighbors);
  }

  void ICAbstractDecorator::AngleCreator::spanLeftAndRightNeighborsForAngle(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& neighbors) {
    while (findLeftAtom(neighbors)) {
      auto copyOfLeftNeighbors = neighbors;
      while (findRightAtom(copyOfLeftNeighbors)) {
		  auto angle = std::make_unique<InternalCoordinates::BondAngle>(
			  bondGraph[leftAtom], bondGraph[middleAtom], bondGraph[rightAtom]);

		  auto constrtaint = pointerToManager->checkIfConstraintPrimitive({ angle->index_a_, angle->index_b_, angle->index_c_ });
		  if (constrtaint) {
			  if (constrtaint->isFrozen()) angle->makeConstrained();
			  else angle->releaseConstraint();
		  }
		  pointerToResult->emplace_back(std::move(angle));
      }
    }
  }

  bool ICAbstractDecorator::AngleCreator::findLeftAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& neighbors) {
    if (neighbors.first == neighbors.second) return false;
    leftAtom = *neighbors.first;
    ++neighbors.first;
    return true;
  }

  bool ICAbstractDecorator::AngleCreator::findRightAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& neighborsLeft) {
    if (neighborsLeft.first == neighborsLeft.second) return false;
    rightAtom = *neighborsLeft.first;
    ++neighborsLeft.first;
    return true;
  }
  
  ICDihedralDecorator::ICDihedralDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    ICAbstractDecorator(parent)
  {}
  
  void ICDihedralDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager){
    DihedralCreator dc(graph);
    appendCoordinates(std::make_shared<ICGeneralAppender>(dc.getInternals(manager)));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec, manager);
  }
  
  ICAbstractDecorator::DihedralCreator::DihedralCreator(BondGraph const& graph):
    DistanceCreator{ graph },
    outerLeft{ 0u },
    outerRight{ 0u },
	  pointerToResult{nullptr},
	  pointerToManager{nullptr}
  {}
  
  InternalVec ICAbstractDecorator::DihedralCreator::getInternals(AbstractConstraintManager& manager) {
    InternalVec result;
    pointerToResult = &result;
	pointerToManager = &manager;
    while (nextEdgeDistances()) {
      findLeftAndRightAtoms();
    }
    // Add new dihedrals if we still have constraints left
	auto num_atoms = boost::num_vertices(bondGraph);
    for (auto const& curr_constraint : manager.getConstraintsOfType(config::Constraint::DIHEDRAL)){
        auto index1 = curr_constraint->getAtomIndices()[0]-1, index2 = curr_constraint->getAtomIndices()[1]-1, index3 = curr_constraint->getAtomIndices()[2]-1, index4 = curr_constraint->getAtomIndices()[3]-1;
        if (index1 >= num_atoms || index2 >= num_atoms || index3 >= num_atoms || index4 >= num_atoms){
          throw std::runtime_error("Cannot create constrained dihedral coordinate: Atom index out of range");
        }
        result.emplace_back(std::make_unique<InternalCoordinates::DihedralAngle>(bondGraph[index1], bondGraph[index2], bondGraph[index3], bondGraph[index4]));
		if (curr_constraint->isFrozen()) pointerToResult->back()->makeConstrained();
		else pointerToResult->back()->releaseConstraint();
      
    }
	pointerToResult = nullptr;
	pointerToManager = nullptr;

    return result;
  }

  void ICAbstractDecorator::DihedralCreator::findLeftAndRightAtoms() {
    auto leftVertices = boost::adjacent_vertices(source, bondGraph);
    while (findLeftAtoms(leftVertices)) {
      auto rightVertices = boost::adjacent_vertices(target, bondGraph);
      while (findRightAtoms(rightVertices)) {
		  auto dihedral = std::make_unique<InternalCoordinates::DihedralAngle>(
			  bondGraph[outerLeft], bondGraph[source], bondGraph[target], bondGraph[outerRight]);
		  auto constrtaint = pointerToManager->checkIfConstraintPrimitive({ dihedral->index_a_, dihedral->index_b_, dihedral->index_c_, dihedral->index_d_ });
		  if (constrtaint) { 
			  if(constrtaint->isFrozen()) dihedral->makeConstrained();
			  else dihedral->releaseConstraint();
		  }
		  pointerToResult->emplace_back(std::move(dihedral));
      }
    }
  }

  bool ICAbstractDecorator::DihedralCreator::findLeftAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& sourceNeighbors) {
    if (sourceNeighbors.first == sourceNeighbors.second) return false;
    outerLeft = *sourceNeighbors.first;
    ++sourceNeighbors.first;
    if (outerLeft == target) return findLeftAtoms(sourceNeighbors);
    return true;
  }

  bool ICAbstractDecorator::DihedralCreator::findRightAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& targetNeighbors) {
    if (targetNeighbors.first == targetNeighbors.second) return false;
    outerRight = *targetNeighbors.first;
    ++targetNeighbors.first;
    if (outerRight == source) return findRightAtoms(targetNeighbors);
    return true;
  }

  ICTranslationDecorator::ICTranslationDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    ICAbstractDecorator(parent)
  {}
  
  void ICTranslationDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager){
    InternalVec result;
	pointerToResult = &result;
	pointerToManager = &manager;
    for (auto const& indices : indexVec){
		auto transX = std::make_unique<InternalCoordinates::TranslationX>(indices);
		auto transY = std::make_unique<InternalCoordinates::TranslationY>(indices);
		auto transZ = std::make_unique<InternalCoordinates::TranslationZ>(indices);

		std::shared_ptr<config::AbstractConstraint> constraint;
		if(constraint = pointerToManager->checkIfConstraintTrans(indices, config::Constraint::TRANSLATION_X)) {
			if (constraint->isFrozen()) transX->makeConstrained();
			else transX->releaseConstraint();
		}
		if (constraint = pointerToManager->checkIfConstraintTrans(indices, config::Constraint::TRANSLATION_Y)) {
				if (constraint->isFrozen()) transY->makeConstrained();
				else transY->releaseConstraint();
		}
		if (constraint = pointerToManager->checkIfConstraintTrans(indices, config::Constraint::TRANSLATION_Z)) {
			if (constraint->isFrozen()) transZ->makeConstrained();
			else transZ->releaseConstraint();
		}
		pointerToResult->emplace_back(std::move(transX));
		pointerToResult->emplace_back(std::move(transY));
		pointerToResult->emplace_back(std::move(transZ));
    }

	auto num_atoms = boost::num_vertices(graph);
	auto rest = manager.getConstraintsOfType(config::Constraint::TRANSLATION_X);
	for (auto const& curr_constraint : rest) {
		auto index1 = curr_constraint->getAtomIndices()[0] - 1, index2 = curr_constraint->getAtomIndices()[1] - 1;
		if (index1 >= num_atoms || index2 >= num_atoms) {
			throw std::runtime_error("Cannot create constrained bond length coordinate: Atom index out of range");
		}
		result.emplace_back(std::make_unique<InternalCoordinates::TranslationX>(curr_constraint->getAtomIndices()));
		if (curr_constraint->isFrozen()) result.back()->makeConstrained();
		else result.back()->releaseConstraint();
	}

	rest = manager.getConstraintsOfType(config::Constraint::TRANSLATION_Y);
	for (auto const& curr_constraint : rest) {
		auto index1 = curr_constraint->getAtomIndices()[0] - 1, index2 = curr_constraint->getAtomIndices()[1] - 1;
		if (index1 >= num_atoms || index2 >= num_atoms) {
			throw std::runtime_error("Cannot create constrained bond length coordinate: Atom index out of range");
		}
		result.emplace_back(std::make_unique<InternalCoordinates::TranslationY>(curr_constraint->getAtomIndices()));
		if (curr_constraint->isFrozen()) result.back()->makeConstrained();
		else result.back()->releaseConstraint();
	}

	rest = manager.getConstraintsOfType(config::Constraint::TRANSLATION_Z);
	for (auto const& curr_constraint : rest) {
		auto index1 = curr_constraint->getAtomIndices()[0] - 1, index2 = curr_constraint->getAtomIndices()[1] - 1;
		if (index1 >= num_atoms || index2 >= num_atoms) {
			throw std::runtime_error("Cannot create constrained bond length coordinate: Atom index out of range");
		}
		result.emplace_back(std::make_unique<InternalCoordinates::TranslationZ>(curr_constraint->getAtomIndices()));
		if (curr_constraint->isFrozen()) result.back()->makeConstrained();
		else result.back()->releaseConstraint();
	}

    appendCoordinates(std::make_shared<ICGeneralAppender>(std::move(result)));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec, manager);
  }
  
  ICRotationAppender::ICRotationAppender(InternalVec && internal_coords, std::vector<std::shared_ptr<InternalCoordinates::Rotator>> && rotators):
    ICGeneralAppender{ std::move(internal_coords) },
    rotators_{ std::move(rotators) }
  {}
  
  void ICRotationAppender::append(std::shared_ptr<PrimitiveInternalCoordinates> primitives){
    primitives->appendPrimitives(std::move(internal_coords_));
    primitives->appendRotators(rotators_);
  }
  
  ICRotationDecorator::ICRotationDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    ICAbstractDecorator(parent)
  {}
  
  void ICRotationDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager){
    InternalVec result;
	pointerToResult = &result;
	pointerToManager = &manager;

	std::vector<std::shared_ptr<InternalCoordinates::Rotator>> rotators;
	for (auto const& curr_indices : indexVec) {

		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, curr_indices)->makeRotations();

		std::shared_ptr<config::AbstractConstraint> constraint;
		if (constraint = pointerToManager->checkIfConstraintRot(curr_indices, config::Constraint::ROTATION_A)) {
			if (constraint->isFrozen()) curr_rotations.rotationA->makeConstrained();
			else curr_rotations.rotationA->releaseConstraint();
		}
		if(constraint = pointerToManager->checkIfConstraintRot(curr_indices, config::Constraint::ROTATION_B)){
			if (constraint->isFrozen()) curr_rotations.rotationB->makeConstrained();
			else curr_rotations.rotationB->releaseConstraint();
		}
		if (constraint = pointerToManager->checkIfConstraintRot(curr_indices, config::Constraint::ROTATION_C)) {
			if (constraint->isFrozen()) curr_rotations.rotationC->makeConstrained();
			else curr_rotations.rotationC->releaseConstraint();
		}

		pointerToResult->emplace_back(std::move(curr_rotations.rotationA));
		pointerToResult->emplace_back(std::move(curr_rotations.rotationB));
		pointerToResult->emplace_back(std::move(curr_rotations.rotationC));

		rotators.emplace_back(curr_rotations.rotator);
	}

	auto num_atoms = boost::num_vertices(graph);
	auto rest = manager.getConstraintsOfType(config::Constraint::ROTATION_A);
	for (auto const& curr_constraint : rest) {
		auto index1 = curr_constraint->getAtomIndices()[0] - 1, index2 = curr_constraint->getAtomIndices()[1] - 1;
		if (index1 >= num_atoms || index2 >= num_atoms) {
			throw std::runtime_error("Cannot create constrained bond length coordinate: Atom index out of range");
		}
		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, curr_constraint->getAtomIndices())->makeRotations();

		if (curr_constraint->isFrozen()) curr_rotations.rotationA->makeConstrained();
		else curr_rotations.rotationA->releaseConstraint();

		result.emplace_back(std::move(curr_rotations.rotationA));
	}

	rest = manager.getConstraintsOfType(config::Constraint::ROTATION_B);
	for (auto const& curr_constraint : rest) {
		auto index1 = curr_constraint->getAtomIndices()[0] - 1, index2 = curr_constraint->getAtomIndices()[1] - 1;
		if (index1 >= num_atoms || index2 >= num_atoms) {
			throw std::runtime_error("Cannot create constrained bond length coordinate: Atom index out of range");
		}
		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, curr_constraint->getAtomIndices())->makeRotations();

		if (curr_constraint->isFrozen()) curr_rotations.rotationB->makeConstrained();
		else curr_rotations.rotationB->releaseConstraint();

		result.emplace_back(std::move(curr_rotations.rotationB));
	}

	rest = manager.getConstraintsOfType(config::Constraint::ROTATION_C);
	for (auto const& curr_constraint : rest) {
		auto index1 = curr_constraint->getAtomIndices()[0] - 1, index2 = curr_constraint->getAtomIndices()[1] - 1;
		if (index1 >= num_atoms || index2 >= num_atoms) {
			throw std::runtime_error("Cannot create constrained bond length coordinate: Atom index out of range");
		}
		auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, curr_constraint->getAtomIndices())->makeRotations();

		if (curr_constraint->isFrozen()) curr_rotations.rotationC->makeConstrained();
		else curr_rotations.rotationC->releaseConstraint();

		result.emplace_back(std::move(curr_rotations.rotationC));
	}

    appendCoordinates(std::make_shared<ICRotationAppender>(std::move(result), std::move(rotators)));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec, manager);
  }
  
  ICOutOfPlaneDecorator::ICOutOfPlaneDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    ICAbstractDecorator(parent)
  {}
  
  void ICOutOfPlaneDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager){
    appendCoordinates(std::make_shared<ICGeneralAppender>(create_oops(cartesians, graph)));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec, manager);
  }
  
  //This function surely does not work.
  inline InternalVec ICOutOfPlaneDecorator::create_oops(const coords::Representation_3D& coords, const BondGraph& g) const {
    using boost::adjacent_vertices;
    using boost::vertices;
    using scon::dot;

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    auto vert = vertices(g);
    for (auto it = vert.first; it != vert.second; ++it) {
      auto core = g[*it].atom_serial;
      auto core_cp = coords.at(core - 1u);
      auto vert_i = adjacent_vertices(*it, g);
      auto sym = g[*it].element;
      if ((vert_i.second - vert_i.first) >= 3) {
        auto permutations = possible_sets_of_3(vert_i.first, vert_i.second);
        for (auto & combination : permutations) {
          std::sort(combination.begin(), combination.end());


          //auto permutation_vec = ic_util::permutation_from_vec(neighbours);
          //for (auto& permutation : permutation_vec) {
          auto u_cp = coords.at(combination.at(0));
          auto v_cp = coords.at(combination.at(1));
          auto w_cp = coords.at(combination.at(2));
          auto n_vec1 = ic_util::normal_unit_vector(u_cp, v_cp, w_cp);
          auto n_vec2 = ic_util::normal_unit_vector(core_cp, u_cp, v_cp);
          auto dot_n_vecs = dot(n_vec1, n_vec2);
          if (0.95 < std::fabs(dot_n_vecs)) {
            result.emplace_back(std::make_unique<InternalCoordinates::OutOfPlane>(g[*it], g[combination.at(0)], g[combination.at(1)], g[combination.at(2)]));
          }
        }
        //}
      }
    }
    return result;
  }
  
  std::vector<std::vector<std::size_t>> ICOutOfPlaneDecorator::possible_sets_of_3(BondGraph::adjacency_iterator const vbegin, BondGraph::adjacency_iterator const vend) {
    std::vector<std::vector<std::size_t>> result;
    for (auto first = vbegin; first < vend - 2; ++first) {
      for (auto second = first + 1; second < vend - 1; ++second) {
        for (auto third = first + 2; third < vend; ++third) {
          result.emplace_back(std::vector<std::size_t>{*first, *second, *third});
        }
      }
    }
    return result;
  }

  std::shared_ptr<config::AbstractConstraint> ConstraintManager::checkForBonds(std::vector<std::size_t> const& constraints) {
	  for (auto it = copiedConstraints.begin(); it != copiedConstraints.end(); ++it) {
		  auto const& constraint = (*it);
		  if (constraint->getAtomIndices().size() != constraints.size()) continue;
		  if ((constraint->getAtomIndices()[0] == constraints[0]+1 && constraint->getAtomIndices()[1] == constraints[1]+1) ||
			  (constraint->getAtomIndices()[0] == constraints[1]+1 && constraint->getAtomIndices()[1] == constraints[0]+1)) {
			  auto ret = *it;
			  copiedConstraints.erase(it);
			  return ret;
		  }
	  }
	  if (constrainDistances) {
		  masterConstraints->emplace_back(std::make_shared < config::BondConstraint>(std::vector<std::size_t>(constraints), true));
		  return masterConstraints->back();
	  }

	  return std::shared_ptr<config::AbstractConstraint>();
  }

  std::shared_ptr<config::AbstractConstraint> ConstraintManager::checkForAngles(std::vector<std::size_t> const& constraints) {
	  for (auto it = copiedConstraints.begin(); it != copiedConstraints.end(); ++it) {
		  auto const& constraint = (*it);
		  if (constraint->getAtomIndices().size() != constraints.size()) continue;
		  if ((constraint->getAtomIndices()[0] == constraints[0]+1 && constraint->getAtomIndices()[1] == constraints[1]+1 && constraint->getAtomIndices()[2] == constraints[2]+1) ||
			  (constraint->getAtomIndices()[0] == constraints[2]+1 && constraint->getAtomIndices()[1] == constraints[1]+1 && constraint->getAtomIndices()[2] == constraints[0]+1)) {
			  auto ret = *it;
			  copiedConstraints.erase(it);
			  return ret;
		  }
	  }

	  if (constrainAngles) {
		  masterConstraints->emplace_back(std::make_shared < config::AngleConstraint>(std::vector<std::size_t>(constraints), true));
		  return masterConstraints->back();
	  }

	  return std::shared_ptr<config::AbstractConstraint>();
  }

  std::shared_ptr<config::AbstractConstraint> ConstraintManager::checkForDihedrals(std::vector<std::size_t> const& constraints) {
	  for (auto it = copiedConstraints.begin(); it != copiedConstraints.end(); ++it) {
		  auto const& constraint = (*it);
		  if (constraint->getAtomIndices().size() != constraints.size()) continue;
		  if ((constraint->getAtomIndices()[0] == constraints[0]+1 && constraint->getAtomIndices()[1] == constraints[1]+1 && constraint->getAtomIndices()[2] == constraints[2]+1 && constraint->getAtomIndices()[3] == constraints[3]+1) ||
			  (constraint->getAtomIndices()[3] == constraints[0]+1 && constraint->getAtomIndices()[2] == constraints[1]+1 && constraint->getAtomIndices()[1] == constraints[2]+1 && constraint->getAtomIndices()[0] == constraints[3]+1)) {
			  auto ret = *it;
			  copiedConstraints.erase(it);
			  return ret;
		  }
	  }

	  if (constrainDihedrals) {
		  masterConstraints->emplace_back(std::make_shared < config::DihedralConstraint>(std::vector<std::size_t>(constraints), true));
		  return masterConstraints->back();
	  }

	  return std::shared_ptr<config::AbstractConstraint>();
  }

  std::shared_ptr<config::AbstractConstraint> ConstraintManager::checkIfConstraintPrimitive(std::vector<std::size_t> const& constraints) {
	  if (constraints.size() == 2) {
		  return checkForBonds(constraints);
	  }
	  else if (constraints.size() == 3) {
		  return checkForAngles(constraints);
	  }
	  else if (constraints.size() == 4) {
		  return checkForDihedrals(constraints);
	  }
  }


  bool isSameSet(std::vector<std::size_t> lhs, std::vector<std::size_t> rhs) {
	  if (lhs.size() != rhs.size()) return false;
	  std::sort(lhs.begin(), lhs.end());
	  std::sort(rhs.begin(), rhs.end());
	  return std::includes(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
  }

  std::shared_ptr<config::AbstractConstraint> ConstraintManager::checkIfConstraintTrans(std::vector<std::size_t> const& constraints, config::Constraint type) {
	  /*std::vector<std::size_t> newConstraints;
	  newConstraints.reserve(constraints.size());
	  std::transform(constraints.cbegin(), constraints.cend(), std::back_inserter(newConstraints), [](std::size_t s) {return ++s; });*/
	  for (auto it = copiedConstraints.begin(); it != copiedConstraints.end(); ++it) {
		  auto const& constraint = *it;
		  if (constraint->getType() == type && isSameSet(constraint->getAtomIndices(), constraints)) {
			  auto ret = *it;
			  copiedConstraints.erase(it);
			  return ret;
		  }
	  }

	  if (constrainTranslations) {
		  if (type == config::Constraint::TRANSLATION_X) {
			  masterConstraints->emplace_back(std::make_shared < config::TranslationXConstraint>(std::vector<std::size_t>(constraints), true));
		  }
		  else if (type == config::Constraint::TRANSLATION_Y) {
			  masterConstraints->emplace_back(std::make_shared < config::TranslationYConstraint>(std::vector<std::size_t>(constraints), true));
		  }
		  else if (type == config::Constraint::TRANSLATION_Z) {
			  masterConstraints->emplace_back(std::make_shared < config::TranslationZConstraint>(std::vector<std::size_t>(constraints), true));
		  }
		  
		  return masterConstraints->back();
	  }

	  return std::shared_ptr<config::AbstractConstraint>();
  }

  std::shared_ptr<config::AbstractConstraint> ConstraintManager::checkIfConstraintRot(std::vector<std::size_t> const& constraints, config::Constraint type) {

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
		  if (type == config::Constraint::ROTATION_A) {
			  masterConstraints->emplace_back(std::make_shared < config::RotationAConstraint>(std::vector<std::size_t>(constraints), true));
		  }
		  else if (type == config::Constraint::ROTATION_B) {
			  masterConstraints->emplace_back(std::make_shared < config::RotationBConstraint>(std::vector<std::size_t>(constraints), true));
		  }
		  else if (type == config::Constraint::ROTATION_C) {
			  masterConstraints->emplace_back(std::make_shared < config::RotationCConstraint>(std::vector<std::size_t>(constraints), true));
		  }

		  return masterConstraints->back();
	  }

	  return std::shared_ptr<config::AbstractConstraint>();

  }

  ConstraintManager::ConstrainVec ConstraintManager::getConstraintsOfType(config::Constraint const type) {
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

  std::shared_ptr<config::AbstractConstraint> NoConstraintManager::checkIfConstraintPrimitive(std::vector<std::size_t> const&) { return std::shared_ptr<config::AbstractConstraint>{}; }
  std::shared_ptr<config::AbstractConstraint> NoConstraintManager::checkIfConstraintTrans(std::vector<std::size_t> const&, config::Constraint const) { return std::shared_ptr<config::AbstractConstraint>{}; }
  std::shared_ptr<config::AbstractConstraint> NoConstraintManager::checkIfConstraintRot(std::vector<std::size_t> const&, config::Constraint const) { return std::shared_ptr<config::AbstractConstraint>{}; }
  NoConstraintManager::ConstrainVec NoConstraintManager::getConstraintsOfType(config::Constraint const) { return ConstrainVec{}; }
}
