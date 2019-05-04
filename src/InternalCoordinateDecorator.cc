#include "InternalCoordinateDecorator.h"

namespace internals{

  ICAbstractDecorator::ICAbstractDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    parent_{ parent }
  {}
  
  void ICAbstractDecorator::buildCoordinates(CartesianType const& cartesians, BondGraph const& graph, IndexVec const& indexVec){
    parent_->buildCoordinates(cartesians, graph, indexVec);
  }
  
  void ICAbstractDecorator::appendCoordinates(InternalVec && newCoordinates){
    parent_->appendCoordinates(std::move(newCoordinates));
  }
  
  ICAbstractDecorator::InternalCoordinatesCreator::InternalCoordinatesCreator(BondGraph const& graph):
    bondGraph(graph)
  {}
  
  ICBondDecorator::ICBondDecorator(std::shared_ptr<InternalCoordinatesBase> parent):
    ICAbstractDecorator{ parent }
  {}
    
  void ICBondDecorator::buildCoordinates(CartesianType const& cartesians, BondGraph const& graph, IndexVec const& indexVec){
    DistanceCreator dc(graph);
    appendCoordinates(dc.getInternals());
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec);
  }
  
  ICAbstractDecorator::DistanceCreator::DistanceCreator(BondGraph const& graph):
     InternalCoordinatesCreator{ graph },
     source{ 0u },
     target{ 0u },
     edgeIterators{ boost::edges(bondGraph) } 
  {}
  
  InternalCoordinatesBase::InternalVec ICAbstractDecorator::DistanceCreator::getInternals() {
    InternalVec result;
    while (nextEdgeDistances()) {
      result.emplace_back(std::make_unique<InternalCoordinates::BondDistance>(bondGraph[source], bondGraph[target]));
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
    
  void ICAngleDecorator::buildCoordinates(CartesianType const& cartesians, BondGraph const& graph, IndexVec const& indexVec){
    AngleCreator ac(graph);
    appendCoordinates(ac.getInternals());
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec);
  }
  
  ICAbstractDecorator::AngleCreator::AngleCreator(BondGraph const& graph):
    InternalCoordinatesCreator{ graph },
    leftAtom{ 0u },
    middleAtom{ 0u },
    rightAtom{ 0u },
    vertexIterators{ boost::vertices(graph) }
  {}
  
  InternalCoordinatesBase::InternalVec ICAbstractDecorator::AngleCreator::getInternals() {
    InternalVec result;
    pointerToResult = &result;
    while (nextVertex()) {
      addAngleForAllNeighbors();
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
          bondGraph[leftAtom], bondGraph[middleAtom], bondGraph[rightAtom]));
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
  
  void ICDihedralDecorator::buildCoordinates(CartesianType const& cartesians, BondGraph const& graph, IndexVec const& indexVec){
    DihedralCreator dc(graph);
    appendCoordinates(dc.getInternals());
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec);
  }
  
  ICAbstractDecorator::DihedralCreator::DihedralCreator(BondGraph const& graph):
    DistanceCreator{ graph },
    outerLeft{ 0u },
    outerRight{ 0u }
  {}
  
  InternalCoordinatesBase::InternalVec ICAbstractDecorator::DihedralCreator::getInternals() {
    InternalVec result;
    pointerToResult = &result;
    while (nextEdgeDistances()) {
      findLeftAndRightAtoms();
    }
    return result;
  }

  void ICAbstractDecorator::DihedralCreator::findLeftAndRightAtoms() {
    auto leftVertices = boost::adjacent_vertices(source, bondGraph);
    while (findLeftAtoms(leftVertices)) {
      auto rightVertices = boost::adjacent_vertices(target, bondGraph);
      while (findRightAtoms(rightVertices)) {
        pointerToResult->emplace_back(std::make_unique<InternalCoordinates::DihedralAngle>(
          bondGraph[outerLeft], bondGraph[source], bondGraph[target], bondGraph[outerRight]));
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
  
  void ICTranslationDecorator::buildCoordinates(CartesianType const& cartesians, BondGraph const& graph, IndexVec const& indexVec){
    InternalVec result;
    for (auto const& indices : indexVec){
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationX>(indices));
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationY>(indices));
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationZ>(indices));
    }
    appendCoordinates(std::move(result));
    ICAbstractDecorator::buildCoordinates(cartesians, graph, indexVec);
  }
}

