/**
CAST 3
InternalCoordinateBase.h
Purpose: Decorators for InternalCoordinateBase


@author Julian Erdmannsdörfer, Christian Schärf
@version 3.0
*/


#ifndef INTERNAL_COORDINATE_DECORATOR
#define INTERNAL_COORDINATE_DECORATOR

#include "InternalCoordinateBase.h"

namespace internals{
  class ICAbstractDecorator : public InternalCoordinatesBase{

  public:
    ICAbstractDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager) override;
    virtual void appendCoordinates(std::shared_ptr<InternalCoordinateAppenderInterface> appender) override;
  
  protected:
    std::shared_ptr<InternalCoordinatesBase> parent_;
    
    class InternalCoordinatesCreator {
    //protected:
    //  using BondGraph = ic_util::Graph<ic_util::Node>;
      
    public:
      InternalCoordinatesCreator(BondGraph const& graph);
      virtual ~InternalCoordinatesCreator() = default;
      virtual InternalVec getInternals(AbstractConstraintManager& manager) = 0;
    protected:
      BondGraph const& bondGraph;
    };
    
    class DistanceCreator : public InternalCoordinatesCreator {
    public:
      DistanceCreator(BondGraph const& graph);
      virtual ~DistanceCreator() = default;

      virtual InternalVec getInternals(AbstractConstraintManager& manager) override;

    protected:
      bool nextEdgeDistances();
      std::size_t source, target;
      std::pair<ic_util::Graph<ic_util::Node>::edge_iterator, ic_util::Graph<ic_util::Node>::edge_iterator> edgeIterators;
	  InternalVec * pointerToResult;
    
    };
    
    class AngleCreator : public InternalCoordinatesCreator {
    public:
      AngleCreator(BondGraph const& graph);
      virtual ~AngleCreator() = default;

      virtual InternalVec getInternals(AbstractConstraintManager& manager) override;
    protected:
      bool nextVertex();
      void addAngleForAllNeighbors();
      void spanLeftAndRightNeighborsForAngle(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & neighbors);
      bool findLeftAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & neighbors);

      bool findRightAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & neighborsLeft);

      std::size_t leftAtom, middleAtom, rightAtom;
      std::pair<ic_util::Graph<ic_util::Node>::vertex_iterator, ic_util::Graph<ic_util::Node>::vertex_iterator> vertexIterators;
      InternalVec * pointerToResult;
	  AbstractConstraintManager * pointerToManager;
    
    };
    
    class DihedralCreator : public DistanceCreator {
    public:
      DihedralCreator(BondGraph const& graph);
      virtual ~DihedralCreator() = default;

      virtual InternalVec getInternals(AbstractConstraintManager& manager) override;
    protected:
      void findLeftAndRightAtoms();
      bool findLeftAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & sourceNeighbors);
      bool findRightAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & targetNeighbors);

      std::size_t outerLeft, outerRight;

      InternalVec * pointerToResult;
	  AbstractConstraintManager * pointerToManager;
    
    };
  };
  
  class ICGeneralAppender : public InternalCoordinateAppenderInterface{
  public:
    ICGeneralAppender(InternalVec && internal_coords);
  
    virtual void append(std::shared_ptr<PrimitiveInternalCoordinates> primitives) override;
    
  protected:
    InternalVec internal_coords_;
  };
  
  class ICBondDecorator : public ICAbstractDecorator{
  public:
    ICBondDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager) override;
  };
  
  class ICAngleDecorator : public ICAbstractDecorator{
  public:
    ICAngleDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager) override;
  };
  
  class ICDihedralDecorator : public ICAbstractDecorator{
  public:
    ICDihedralDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager) override;
  };

  class ICTranslationDecorator : public ICAbstractDecorator{
  public:
    ICTranslationDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager) override;

  protected:
	  InternalVec * pointerToResult;
	  AbstractConstraintManager * pointerToManager;
  };
  
  class ICRotationAppender : public ICGeneralAppender{
  public:
    ICRotationAppender(InternalVec && internal_coords, std::vector<std::shared_ptr<InternalCoordinates::Rotator>> && rotators);
    
    void append(std::shared_ptr<PrimitiveInternalCoordinates> primitives) override;
    
  protected:
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> rotators_;
  };
  
  class ICRotationDecorator : public ICAbstractDecorator{
  public:
    ICRotationDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager) override;

  protected:
	  InternalVec * pointerToResult;
	  AbstractConstraintManager * pointerToManager;
  };
  
  class ICOutOfPlaneDecorator : public ICAbstractDecorator{
  public:
    ICOutOfPlaneDecorator(std::shared_ptr<InternalCoordinatesBase> parent);
    
    virtual void buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec, AbstractConstraintManager& manager) override;
    
  protected:
    InternalVec create_oops(const coords::Representation_3D& coords, const BondGraph& g) const;
    
    static std::vector<std::vector<std::size_t>> possible_sets_of_3(BondGraph::adjacency_iterator const vbegin, BondGraph::adjacency_iterator const vend);
  };

  class NoConstraintManager : public AbstractConstraintManager {
  public:
	  virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintPrimitive(std::vector<std::size_t> const&) override;
	  virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintTrans(std::vector<std::size_t> const&, config::Constraint const) override;
	  virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintRot(std::vector<std::size_t> const&, config::Constraint const) override;
	  virtual ConstrainVec getConstraintsOfType(config::Constraint const) override;
	  //virtual ConstrainVec && getAllConstraints() override;
  };

  class ConstraintManager : public AbstractConstraintManager {
  public:

	  virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintPrimitive(std::vector<std::size_t> const&) override;
	  virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintTrans(std::vector<std::size_t> const&, config::Constraint const) override;
	  virtual std::shared_ptr<config::AbstractConstraint> checkIfConstraintRot(std::vector<std::size_t> const&, config::Constraint const) override;
	  virtual ConstrainVec getConstraintsOfType(config::Constraint const) override;
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

	  std::shared_ptr<config::AbstractConstraint> checkForBonds(std::vector<std::size_t> const&);
	  std::shared_ptr<config::AbstractConstraint> checkForAngles(std::vector<std::size_t> const&);
	  std::shared_ptr<config::AbstractConstraint> checkForDihedrals(std::vector<std::size_t> const&);

	  ConstrainVec copiedConstraints;
  };
}

#endif // INTERNAL_COORDINATE_DECORATOR
