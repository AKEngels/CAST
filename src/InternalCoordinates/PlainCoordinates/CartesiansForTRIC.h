#ifndef CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANSFORTRIC_H_
#define CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANSFORTRIC_H_


#include "CartesiansForInternals.h"

#include"GeometryObserver.h"
#include "../InternalCoordinates/RotatorObserver.h"

namespace internals {

	class CartesiansForTRIC : public CartesiansForInternalCoordinates {
	public:



		CartesiansForTRIC(CartesiansForTRIC && cartesians);
		CartesiansForTRIC(CartesiansForTRIC const& cartesians);

		CartesiansForTRIC(scon::mathmatrix<float_type> && cartesians);
		CartesiansForTRIC(scon::mathmatrix<float_type> const& cartesians);

		virtual CartesiansForTRIC & operator=(CartesiansForInternalCoordinates const& cartesians) override;
		virtual CartesiansForTRIC & operator=(CartesiansForInternalCoordinates && cartesians) override;


		virtual void copyTo(CartesiansForInternalCoordinates & cartesians) const& override;
		// be sure to really call these functions with a std::move(X).copyTo(Y);
		virtual void copyTo(CartesiansForInternalCoordinates & cartesians) && override;

		virtual void copyTo(CartesiansForTRIC & cartesians) const& override;
		virtual void copyTo(CartesiansForTRIC & cartesians) && override;

		void registerObserver(std::shared_ptr<RotatorObserver> const observer);
	private:
		friend class CartesiansForInternalCoordinates;
		std::vector<std::shared_ptr<GeometryObserver>> observerList;
		void notify();
	};

}

#endif // CAST_INTERNALCOORDINATES_PLAINCOORDINATES_CARTESIANSFORTRIC_H_