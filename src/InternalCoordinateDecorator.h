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
  class PrimitiveInternalCoordinates;

  class ICAbstractDecorator{

  public:
    ICAbstractDecorator(std::unique_ptr<ICAbstractDecorator> parent);
    
    virtual void buildCoordinates(CartesianType & cartesians,
                                  BondGraph const& graph,
                                  IndexVec const& indexVec,
                                  AbstractConstraintManager& manager);

    virtual void appendCoordinates(PrimitiveInternalCoordinates & primitiveInternals);

  protected:
    std::unique_ptr<ICAbstractDecorator> parent_;

    void storeInternals(InternalVec && new_internals);

  private:
    InternalVec created_internals_;

  protected:
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
  
  class ICBondDecorator : public ICAbstractDecorator{
  public:
    using ICAbstractDecorator::ICAbstractDecorator;

    virtual void buildCoordinates(CartesianType & cartesians,
                                  BondGraph const& graph,
                                  IndexVec const& indexVec,
                                  AbstractConstraintManager& manager);
  };
  
  class ICAngleDecorator : public ICAbstractDecorator{
  public:
    using ICAbstractDecorator::ICAbstractDecorator;

    virtual void buildCoordinates(CartesianType& cartesians,
                                  BondGraph const& graph,
                                  IndexVec const& indexVec,
                                  AbstractConstraintManager& manager);
  };
  
  class ICDihedralDecorator : public ICAbstractDecorator{
  public:
    using ICAbstractDecorator::ICAbstractDecorator;

    virtual void buildCoordinates(CartesianType& cartesians,
                                  BondGraph const& graph,
                                  IndexVec const& indexVec,
                                  AbstractConstraintManager& manager);
  };

  class ICTranslationDecorator : public ICAbstractDecorator{
  public:
    using ICAbstractDecorator::ICAbstractDecorator;
    
    virtual void buildCoordinates(CartesianType& cartesians,
                                  BondGraph const& graph,
                                  IndexVec const& indexVec,
                                  AbstractConstraintManager& manager);

  protected:
	  InternalVec * pointerToResult;
	  AbstractConstraintManager * pointerToManager;
  };
  
  class ICRotationDecorator : public ICAbstractDecorator{
  public:
    using ICAbstractDecorator::ICAbstractDecorator;
    
    virtual void buildCoordinates(CartesianType& cartesians,
                                  BondGraph const& graph,
                                  IndexVec const& indexVec,
                                  AbstractConstraintManager& manager);

    virtual void appendCoordinates(PrimitiveInternalCoordinates & primitiveInternals) override;

  protected:
    InternalVec * pointerToResult;
    AbstractConstraintManager * pointerToManager;
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> createdRotators_;
  };
  
  class ICOutOfPlaneDecorator : public ICAbstractDecorator{
  public:
    using ICAbstractDecorator::ICAbstractDecorator;

    virtual void buildCoordinates(CartesianType& cartesians,
                                  BondGraph const& graph,
                                  IndexVec const& indexVec,
                                  AbstractConstraintManager& manager);
    
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
