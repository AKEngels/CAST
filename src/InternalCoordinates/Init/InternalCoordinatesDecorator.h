/**
CAST 3
InternalCoordinateBase.h
Purpose: Decorators for InternalCoordinateBase


@author Julian Erdmannsdörfer, Christian Schärf
@version 3.0
*/


#ifndef CAST_INTERNALCOORDINATES_INIT_INTERNALCOORDINATESDECORATOR_H_
#define CAST_INTERNALCOORDINATES_INIT_INTERNALCOORDINATESDECORATOR_H_

#include "../InternalCoordinates.h"
#include "../BondGraph/BondGraph.h"

#include "ConstraintManager.h"

#include<memory>

namespace internals {

	using InternalVector = std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>;

	using BondGraph = ic_util::BondGraph;

	

	class PrimitiveInternalCoordinates;

	class PrimitiveInternalsCreator {

	public:
		virtual void appendInternalCoordinates(InternalVector & vec) = 0;
		virtual ~PrimitiveInternalsCreator() = 0;
	};

	PrimitiveInternalsCreator::~PrimitiveInternalsCreator() = default;
	
	class ConstraintAdder : public PrimitiveInternalsCreator {
	public:
		ConstraintAdder(std::shared_ptr<ConstraintManager> m, InternalCoordinates::InternalCoordinatesBuilder & builder) : manager{ std::move(m) }, coordinatesBuilder{ builder } {}
		virtual void appendInternalCoordinates(InternalVector & vec) { applyConstraints(vec); }

	private:
		void applyConstraints(InternalVector & vec);
		std::shared_ptr<ConstraintManager> manager;

		InternalCoordinates::InternalCoordinatesBuilder & coordinatesBuilder;
	};

	class BaseCreator : public PrimitiveInternalsCreator {
	public:
		virtual void appendInternalCoordinates(InternalVector & vec) override {
			auto internals = findInternalCoordinates();
			vec.insert(vec.end(), std::make_move_iterator(internals.begin()), std::make_move_iterator(internals.end()));
		}
	protected:
		virtual InternalVector findInternalCoordinates() = 0;
	};

	class BondCreator : public BaseCreator {
	public:
		BondCreator(BondGraph const& g) : graph{ g } {}

	protected:
		virtual InternalVector findInternalCoordinates() override;
	private:
		BondGraph const& graph;
	};

	class AngleCreator : public BaseCreator {
	public:
		AngleCreator(BondGraph const& g) : graph{ g } {}
	protected:
		virtual InternalVector findInternalCoordinates() override;
	private:
		BondGraph const& graph;
	};

	class DihedralCreator : public BaseCreator {
	public:
		DihedralCreator(BondGraph const& g) : graph{ g } {}
	protected:
		virtual InternalVector findInternalCoordinates() override;
	private:
		BondGraph const& graph;
	};

	class TranslationCreator : public BaseCreator {
	public:
		TranslationCreator(BondGraph & g) : graph{ g } {}
	protected:
		virtual InternalVector findInternalCoordinates() override;
	private:
		BondGraph & graph;
	};

	class RotationCereator : public BaseCreator {
	public:
		RotationCereator(BondGraph & g, InternalCoordinates::CartesiansForInternalCoordinates & coords, std::vector<std::shared_ptr<InternalCoordinates::Rotator>> & listOfRotators) : graph{ g }, coordinates{ coords }, rotators{ listOfRotators } {}
	protected:
		virtual InternalVector findInternalCoordinates() override;
	private:
		BondGraph & graph;
		InternalCoordinates::CartesiansForInternalCoordinates & coordinates;
		std::vector<std::shared_ptr<InternalCoordinates::Rotator>> & rotators;
	};

	class OutOfPlaneCreator : public BaseCreator {
	public:
		OutOfPlaneCreator(BondGraph const& g) : graph{ g } {}
	protected:
		virtual InternalVector findInternalCoordinates() override;
	private:
		BondGraph const& graph;
	};

	

	
}

#endif // CAST_INTERNALCOORDINATES_INIT_INTERNALCOORDINATESDECORATOR_H_
