#include "InternalCoordinateDecorator.h"

#include "PrimitiveInternalCoordinates.h"

namespace internals{

  ICAbstractDecorator::ICAbstractDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    parent_{ parent }
  {}
  
  void ICAbstractDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec){
    parent_->buildCoordinates(cartesians, graph, indexVec);
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
    
  void ICBondDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec){
    DistanceCreator dc(graph);
    appendCoordinates(std::make_shared<ICGeneralAppender>(dc.getInternals()));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec);
  }
  
  ICAbstractDecorator::DistanceCreator::DistanceCreator(BondGraph const& graph):
     InternalCoordinatesCreator{ graph },
     source{ 0u },
     target{ 0u },
     edgeIterators{ boost::edges(bondGraph) },
     constraintManager_{ Config::get().constrained_internals.constrained_bond_lengths }
  {}
  
  InternalVec ICAbstractDecorator::DistanceCreator::getInternals() {
    InternalVec result;
    while (nextEdgeDistances()) {
      result.emplace_back(std::make_unique<InternalCoordinates::BondDistance>(bondGraph[source], bondGraph[target], constraintManager_.pop_constraint({source+1, target+1})));
    }
    // Add new bonds if we still have constraints left
    for (auto const& curr_constraint : constraintManager_.get_constraints()){
      if(curr_constraint.second != Config::get().constrained_internals.constrain_bond_lengths){
        auto const& atom_indices = curr_constraint.first;
        auto index1 = atom_indices[0]-1, index2 = atom_indices[1]-1;
        auto num_atoms = boost::num_vertices(bondGraph);
        if (index1 >= num_atoms || index2 >= num_atoms){
          throw std::runtime_error("Cannot create constrained bond length coordinate: Atom index out of range");
        }
        result.emplace_back(std::make_unique<InternalCoordinates::BondDistance>(bondGraph[index1], bondGraph[index2], curr_constraint.second));
      }
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
    
  void ICAngleDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec){
    AngleCreator ac(graph);
    appendCoordinates(std::make_shared<ICGeneralAppender>(ac.getInternals()));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec);
  }
  
  ICAbstractDecorator::AngleCreator::AngleCreator(BondGraph const& graph):
    InternalCoordinatesCreator{ graph },
    leftAtom{ 0u },
    middleAtom{ 0u },
    rightAtom{ 0u },
    vertexIterators{ boost::vertices(graph) },
    constraintManager_{ Config::get().constrained_internals.constrained_bond_angles }
  {}
  
  InternalVec ICAbstractDecorator::AngleCreator::getInternals() {
    InternalVec result;
    pointerToResult = &result;
    while (nextVertex()) {
      addAngleForAllNeighbors();
    }
    // Add new angles if we still have constraints left
    for (auto const& curr_constraint : constraintManager_.get_constraints()){
      if(curr_constraint.second != Config::get().constrained_internals.constrain_bond_angles){
        auto const& atom_indices = curr_constraint.first;
        auto index1 = atom_indices[0]-1, index2 = atom_indices[1]-1, index3 = atom_indices[2]-1;
        auto num_atoms = boost::num_vertices(bondGraph);
        if (index1 >= num_atoms || index2 >= num_atoms || index3 >= num_atoms){
          throw std::runtime_error("Cannot create constrained bond angle coordinate: Atom index out of range");
        }
        result.emplace_back(std::make_unique<InternalCoordinates::BondAngle>(bondGraph[index1], bondGraph[index2], bondGraph[index3], curr_constraint.second));
      }
    }
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
        pointerToResult->emplace_back(std::make_unique<InternalCoordinates::BondAngle>(
          bondGraph[leftAtom], bondGraph[middleAtom], bondGraph[rightAtom],
          constraintManager_.pop_constraint({leftAtom+1, middleAtom+1, rightAtom+1})));
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
  
  void ICDihedralDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec){
    DihedralCreator dc(graph);
    appendCoordinates(std::make_shared<ICGeneralAppender>(dc.getInternals()));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec);
  }
  
  ICAbstractDecorator::DihedralCreator::DihedralCreator(BondGraph const& graph):
    DistanceCreator{ graph },
    outerLeft{ 0u },
    outerRight{ 0u },
    constraintManager_{ Config::get().constrained_internals.constrained_dihedrals }
  {}
  
  InternalVec ICAbstractDecorator::DihedralCreator::getInternals() {
    InternalVec result;
    pointerToResult = &result;
    while (nextEdgeDistances()) {
      findLeftAndRightAtoms();
    }
    // Add new dihedrals if we still have constraints left
    for (auto const& curr_constraint : constraintManager_.get_constraints()){
      if(curr_constraint.second != Config::get().constrained_internals.constrain_bond_angles){
        auto const& atom_indices = curr_constraint.first;
        auto index1 = atom_indices[0]-1, index2 = atom_indices[1]-1, index3 = atom_indices[2]-1, index4 = atom_indices[3]-1;
        auto num_atoms = boost::num_vertices(bondGraph);
        if (index1 >= num_atoms || index2 >= num_atoms || index3 >= num_atoms || index4 >= num_atoms){
          throw std::runtime_error("Cannot create constrained dihedral coordinate: Atom index out of range");
        }
        result.emplace_back(std::make_unique<InternalCoordinates::DihedralAngle>(bondGraph[index1], bondGraph[index2], bondGraph[index3], bondGraph[index4], curr_constraint.second));
      }
    }
    return result;
  }

  void ICAbstractDecorator::DihedralCreator::findLeftAndRightAtoms() {
    auto leftVertices = boost::adjacent_vertices(source, bondGraph);
    while (findLeftAtoms(leftVertices)) {
      auto rightVertices = boost::adjacent_vertices(target, bondGraph);
      while (findRightAtoms(rightVertices)) {
        pointerToResult->emplace_back(std::make_unique<InternalCoordinates::DihedralAngle>(
          bondGraph[outerLeft], bondGraph[source], bondGraph[target], bondGraph[outerRight],
          constraintManager_.pop_constraint({outerLeft+1, source+1, target+1, outerRight+1})));
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
  
  void ICTranslationDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec){
    InternalVec result;
    for (auto const& indices : indexVec){
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationX>(indices));
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationY>(indices));
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationZ>(indices));
    }
    appendCoordinates(std::make_shared<ICGeneralAppender>(std::move(result)));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec);
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
  
  void ICRotationDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec){
    InternalVec result;
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> rotators;
    for (auto const& curr_indices : indexVec){
      auto curr_rotations = InternalCoordinates::Rotator::buildRotator(cartesians, curr_indices)->makeRotations();
      result.emplace_back(std::move(curr_rotations.rotationA));
      result.emplace_back(std::move(curr_rotations.rotationB));
      result.emplace_back(std::move(curr_rotations.rotationC));
      rotators.emplace_back(curr_rotations.rotator);
    }

    appendCoordinates(std::make_shared<ICRotationAppender>(std::move(result), std::move(rotators)));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec);
  }
  
  ICOutOfPlaneDecorator::ICOutOfPlaneDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    ICAbstractDecorator(parent)
  {}
  
  void ICOutOfPlaneDecorator::buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec){
    appendCoordinates(std::make_shared<ICGeneralAppender>(create_oops(cartesians, graph)));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec);
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
}
