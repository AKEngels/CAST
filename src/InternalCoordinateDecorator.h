/**
CAST 3
InternalCoordinateBase.h
Purpose: Decorators for InternalCoordinateBase


@author Julian Erdmannsdörfer, Christian Schärf
@version 3.0
*/


#ifndef INTERNAL_COORDINATE_DECORATOR
#define INTERNAL_COORDIANTE_DECORATOR

#include "InternalCoordinateBase.h"

namespace internals{
  class ICAbstractDecorator : public InternalCoordinatesBase{

  public:
    ICAbstractDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType const& cartesians, BondGraph const& graph, IndexVec const& indexVec) override;
    virtual void appendCoordinates(InternalVec && newCoordinates) override;
  
  protected:
    std::shared_ptr<InternalCoordinatesBase> parent_;
    
    class InternalCoordinatesCreator {
    //protected:
    //  using BondGraph = ic_util::Graph<ic_util::Node>;
      
    public:
      InternalCoordinatesCreator(BondGraph const& graph);
      virtual ~InternalCoordinatesCreator() = default;
      virtual InternalVec getInternals() = 0;
    protected:
      BondGraph const& bondGraph;
    };
    
    class DistanceCreator : public InternalCoordinatesCreator {
    public:
      DistanceCreator(BondGraph const& graph);
      virtual ~DistanceCreator() = default;

      virtual InternalVec getInternals() override;

    protected:
      bool nextEdgeDistances();
      std::size_t source, target;
      std::pair<ic_util::Graph<ic_util::Node>::edge_iterator, ic_util::Graph<ic_util::Node>::edge_iterator> edgeIterators;
    };
    
    class AngleCreator : public InternalCoordinatesCreator {
    public:
      AngleCreator(BondGraph const& graph);
      virtual ~AngleCreator() = default;

      virtual InternalVec getInternals() override;
    protected:
      bool nextVertex();
      void addAngleForAllNeighbors();
      void spanLeftAndRightNeighborsForAngle(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & neighbors);
      bool findLeftAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & neighbors);

      bool findRightAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & neighborsLeft);

      std::size_t leftAtom, middleAtom, rightAtom;
      std::pair<ic_util::Graph<ic_util::Node>::vertex_iterator, ic_util::Graph<ic_util::Node>::vertex_iterator> vertexIterators;
      InternalVec * pointerToResult;
    };
    
    class DihedralCreator : public DistanceCreator {
    public:
      DihedralCreator(BondGraph const& graph);
      virtual ~DihedralCreator() = default;

      virtual InternalVec getInternals() override;
    protected:
      void findLeftAndRightAtoms();
      bool findLeftAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & sourceNeighbors);
      bool findRightAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & targetNeighbors);

      std::size_t outerLeft, outerRight;

      InternalVec * pointerToResult;
    };
  };
  
  class ICBondDecorator : public ICAbstractDecorator{
  public:
    ICBondDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType const& cartesians, BondGraph const& graph, IndexVec const& indexVec) override;    
  };
  
  class ICAngleDecorator : public ICAbstractDecorator{
  public:
    ICAngleDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType const& cartesians, BondGraph const& graph, IndexVec const& indexVec) override;    
  };
  
  class ICDihedralDecorator : public ICAbstractDecorator{
  public:
    ICDihedralDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType const& cartesians, BondGraph const& graph, IndexVec const& indexVec) override;    
  };

  class ICTranslationDecorator : public ICAbstractDecorator{
  public:
    ICTranslationDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType const&cartesians, BondGraph const& graph, IndexVec const& indexVec) override;
  };
}

#endif // INTERNAL_COORDINATE_DECORATOR
