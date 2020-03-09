#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_INTERNALCOORDINATE_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_INTERNALCOORDINATE_H_

#include"../InternalCoordinatesAliases.h"
#include"../Init/ConstraintManager.h"
#include"../Matrix/Matrix.h"

#include<vector>
#include<memory>

namespace internals{

struct InternalCoordinate {
	virtual float_type value(Eigen::MatrixXd const& cartesians) const = 0;
	virtual float_type difference(Eigen::MatrixXd const& newCoordinates, Eigen::MatrixXd const&  oldCoordinates) const = 0;
	virtual Eigen::VectorXd derivativeVector(Eigen::MatrixXd const& cartesians) const = 0;
	virtual float_type hessianGuess(Eigen::MatrixXd const& cartesians) const = 0;
	virtual std::string info(Eigen::MatrixXd const & cartesians) const = 0;
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