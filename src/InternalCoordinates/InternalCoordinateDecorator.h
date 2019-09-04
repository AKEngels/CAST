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
#include "BondGraph.h"

namespace internals {

	using BondGraph = ic_util::Graph<ic_util::Node>;

	

	class PrimitiveInternalCoordinates;

	class ICDecoratorBase {

	public:
		ICDecoratorBase(std::unique_ptr<ICDecoratorBase> parent);

		virtual void buildCoordinates(CartesianType& cartesians,
			BondGraph const& graph,
			IndexVec const& indexVec,
			AbstractConstraintManager& manager);

		virtual void appendCoordinates(PrimitiveInternalCoordinates& primitiveInternals);

		virtual ~ICDecoratorBase() = default;

	protected:
		std::unique_ptr<ICDecoratorBase> parent_;

		void storeInternals(InternalVec&& new_internals);

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
			InternalVec* pointerToResult;

		};

		class AngleCreator : public InternalCoordinatesCreator {
		public:
			AngleCreator(BondGraph const& graph);
			virtual ~AngleCreator() = default;

			virtual InternalVec getInternals(AbstractConstraintManager& manager) override;
		protected:
			bool nextVertex();
			void addAngleForAllNeighbors();
			void spanLeftAndRightNeighborsForAngle(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& neighbors);
			bool findLeftAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& neighbors);

			bool findRightAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& neighborsLeft);

			std::size_t leftAtom, middleAtom, rightAtom;
			std::pair<ic_util::Graph<ic_util::Node>::vertex_iterator, ic_util::Graph<ic_util::Node>::vertex_iterator> vertexIterators;
			InternalVec* pointerToResult;
			AbstractConstraintManager* pointerToManager;

		};

		class DihedralCreator : public DistanceCreator {
		public:
			DihedralCreator(BondGraph const& graph);
			virtual ~DihedralCreator() = default;

			virtual InternalVec getInternals(AbstractConstraintManager& manager) override;
		protected:
			void findLeftAndRightAtoms();
			bool findLeftAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& sourceNeighbors);
			bool findRightAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& targetNeighbors);

			std::size_t outerLeft, outerRight;

			InternalVec* pointerToResult;
			AbstractConstraintManager* pointerToManager;

		};
	};

	class ICBondDecorator : public ICDecoratorBase {
	public:
		using ICDecoratorBase::ICDecoratorBase;

		virtual void buildCoordinates(CartesianType& cartesians,
			BondGraph const& graph,
			IndexVec const& indexVec,
			AbstractConstraintManager& manager);
	};

	class ICAngleDecorator : public ICDecoratorBase {
	public:
		using ICDecoratorBase::ICDecoratorBase;

		virtual void buildCoordinates(CartesianType& cartesians,
			BondGraph const& graph,
			IndexVec const& indexVec,
			AbstractConstraintManager& manager);
	};

	class ICDihedralDecorator : public ICDecoratorBase {
	public:
		using ICDecoratorBase::ICDecoratorBase;

		virtual void buildCoordinates(CartesianType& cartesians,
			BondGraph const& graph,
			IndexVec const& indexVec,
			AbstractConstraintManager& manager);
	};

	class ICTranslationDecorator : public ICDecoratorBase {
	public:
		using ICDecoratorBase::ICDecoratorBase;

		virtual void buildCoordinates(CartesianType& cartesians,
			BondGraph const& graph,
			IndexVec const& indexVec,
			AbstractConstraintManager& manager);

	protected:
		InternalVec* pointerToResult;
		AbstractConstraintManager* pointerToManager;
	};

	class ICRotationDecorator : public ICDecoratorBase {
	public:
		using ICDecoratorBase::ICDecoratorBase;

		virtual void buildCoordinates(CartesianType& cartesians,
			BondGraph const& graph,
			IndexVec const& indexVec,
			AbstractConstraintManager& manager);

		virtual void appendCoordinates(PrimitiveInternalCoordinates& primitiveInternals) override;

	protected:
		InternalVec* pointerToResult;
		AbstractConstraintManager* pointerToManager;
		std::vector<std::shared_ptr<InternalCoordinates::Rotator>> createdRotators_;
	};

	class ICOutOfPlaneDecorator : public ICDecoratorBase {
	public:
		using ICDecoratorBase::ICDecoratorBase;

		virtual void buildCoordinates(CartesianType& cartesians,
			BondGraph const& graph,
			IndexVec const& indexVec,
			AbstractConstraintManager& manager);

	protected:
		InternalVec create_oops(const CartesianType& coords, const BondGraph& g) const;

		static std::vector<std::vector<std::size_t>> possible_sets_of_3(BondGraph::adjacency_iterator const vbegin, BondGraph::adjacency_iterator const vend);
	};

	

	
}

#endif // INTERNAL_COORDINATE_DECORATOR
