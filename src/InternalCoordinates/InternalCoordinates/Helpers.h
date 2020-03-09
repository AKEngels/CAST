#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_HELPERS_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_HELPERS_H_

#include<Eigen/Dense>

#include "../InternalCoordinatesAliases.h"

namespace internals {
	// struct CartesianPointView {
	// 	CartesianPointView(double & x, double & y, double & z);
	// 	double & x;
	// 	double & y;
	// 	double & z;
	// };
	// using CartesianPoint = std::tuple<float_type, float_type, float_type>;
	// CartesianPoint getAtom(scon::mathmatrix<float_type> const& coordinates, std::size_t const index);
	// CartesianPoint operator-(CartesianPoint const&);
	// CartesianPoint operator-(CartesianPoint const&, CartesianPoint const&);
	// CartesianPoint operator+(CartesianPoint const&, CartesianPoint const&);
	// CartesianPoint & operator/=(CartesianPoint &, float_type const);
	// CartesianPoint operator/(CartesianPoint const&, float_type const);
	// CartesianPoint & operator*=(CartesianPoint &, float_type const);
	// CartesianPoint operator*(CartesianPoint const&, float_type const);
	// float_type euclideanLength(CartesianPoint const&);
	// float_type dotProduct(CartesianPoint const&);
	// float_type dotProduct(CartesianPoint const&, CartesianPoint const&);
	// CartesianPoint crossProduct(CartesianPoint const&, CartesianPoint const&);
	// CartesianPoint normalize(CartesianPoint const& a);

	// class BmatrixRowCreator {
	// public:
	// 	BmatrixRowCreator(std::size_t const);
	// 	void insertAtomDerivative(CartesianPoint const&, std::size_t const index);
	// 	scon::mathmatrix<float_type> && getRow();
	// private:
	// 	std::unique_ptr<scon::mathmatrix<float_type>> row;
	// };
}

#endif