#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_INTERNALCOORDINATE_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_INTERNALCOORDINATE_H_

#include"../InternalCoordinatesAliases.h"
#include"../Init/ConstraintManager.h"

#include<vector>
#include<memory>

namespace internals{

struct InternalCoordinate {
	virtual float_type val(scon::mathmatrix<float_type> const& cartesians) const = 0;
	virtual float_type difference(scon::mathmatrix<float_type> const& newCoordinates, scon::mathmatrix<float_type> const&  oldCoordinates) const = 0;
	virtual scon::mathmatrix<float_type> der_vec(scon::mathmatrix<float_type> const& cartesians) const = 0;
	virtual float_type hessian_guess(scon::mathmatrix<float_type> const& cartesians) const = 0;
	virtual std::string info(scon::mathmatrix<float_type> const & cartesians) const = 0;
	virtual bool hasIndices(std::vector<std::size_t> const& indices) const = 0;
	virtual std::vector<std::size_t> getIndices() const = 0;
	virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) = 0;
	virtual void makeConstrained() = 0;
	virtual void releaseConstraint() = 0;
	virtual bool is_constrained() const = 0;
	virtual ~InternalCoordinate() = 0;
};

InternalCoordinate::~InternalCoordinate() = default;

}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_INTERNALCOORDINATE_H_