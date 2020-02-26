#ifndef CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANSFORINTERNALS_H_
#define CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANSFORINTERNALS_H_

#include "CartesianCoordinates.h"
#include "../InternalCoordinates/InternalCoordinate.h"

#include <vector>

namespace internals {
	class CartesiansForTRIC;

	// TODO: Rename as soon as it works
	class CartesiansForInternalCoordinates {
	public:


		CartesiansForInternalCoordinates(std::unique_ptr<CoordinatesContainer> && cartesians);

		CartesiansForInternalCoordinates(std::unique_ptr<AbstractMatrix> && cartesians);
		CartesiansForInternalCoordinates(std::unique_ptr<AbstractMatrix> const& cartesians);

		virtual CartesiansForInternalCoordinates & operator=(CartesiansForInternalCoordinates const& cartesians);
		virtual CartesiansForInternalCoordinates & operator=(CartesiansForInternalCoordinates && cartesians);

		virtual void copyTo(CartesiansForInternalCoordinates & cartesians) const&;
		// be sure to really call these functions with a std::move(X).copyTo(Y);
		virtual void copyTo(CartesiansForInternalCoordinates & cartesians) &&;

		virtual void copyTo(CartesiansForTRIC & cartesians) const&;
		virtual void copyTo(CartesiansForTRIC & cartesians) &&;

		CartesiansForInternalCoordinates() = default;
		virtual ~CartesiansForInternalCoordinates() = default;

		friend scon::mathmatrix<float_type> operator+(CartesiansForInternalCoordinates const& lhs, scon::mathmatrix<float_type> const& rhs);
		friend scon::mathmatrix<float_type> operator+(scon::mathmatrix<float_type> const& rhs, CartesiansForInternalCoordinates const& lhs);

		// This is unsatisfactory I'd rather have a return which is reference like. Though this is hard using an extern algebra library
		scon::mathmatrix<float_type> at(std::size_t const i) const;

		float_type getInternalValue(InternalCoordinate const& in) const;
		std::vector<float_type> getInternalDerivativeVector(InternalCoordinate const& in) const;
		float_type getInternalHessianGuess(InternalCoordinate const& in) const;
		float_type getInternalDifference(CartesiansForInternalCoordinates const& other, InternalCoordinate const& in) const;

		std::pair<float_type, float_type> displacementRmsValAndMaxTwoStructures(AbstractMatrix const& other) const;

		std::pair<float_type, float_type> displacementRmsValAndMaxTwoStructures(CartesiansForInternalCoordinates const& other) const;

		scon::mathmatrix<float_type> const& getCartesians() const;
		scon::mathmatrix<float_type> inBohr() const;
		scon::mathmatrix<float_type> inAngstrom() const;

		virtual void setCartesianCoordnates(scon::mathmatrix<float_type> const& newCartesianCoordinates);
		virtual void setCartesianCoordnates(scon::mathmatrix<float_type> && newCartesianCoordinates);

		void reset();

	protected:
		friend class CartesiansForTRIC;
		std::unique_ptr<CoordinatesContainer> coordinates;
	};

	/*class temporaryCartesian {
	public:
		temporaryCartesian(CartesiansForInternalCoordinates & cartesians) : coordinates{ cartesians.coordinates }, coordinatesPtr{ &cartesians } {}
		coords::Container<coords::Cartesian_Point> coordinates;
		void stolenNotify() {
			coordinatesPtr->notify();
		}
	private:
		CartesiansForInternalCoordinates * coordinatesPtr;
	};*/
}

#endif // CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANSFORINTERNALS_H_