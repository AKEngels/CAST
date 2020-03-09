#include "Helpers.h"

#include <cmath>

namespace internals {
	// CartesianPoint getAtom(scon::mathmatrix<float_type> const& coordinates, std::size_t const index) {
	// 	return CartesianPoint{ coordinates(index, 0), coordinates(index, 1), coordinates(index, 2) };
	// }

	// CartesianPoint operator-(CartesianPoint const& point) {
	// 	return CartesianPoint{
	// 		-std::get<0u>(point),
	// 		-std::get<1u>(point),
	// 		-std::get<2u>(point)
	// 	};
	// }

	// CartesianPoint operator-(CartesianPoint const& lhs, CartesianPoint const& rhs) {
	// 	return CartesianPoint{
	// 		std::get<0u>(lhs) - std::get<0u>(rhs),
	// 		std::get<1u>(lhs) - std::get<1u>(rhs),
	// 		std::get<2u>(lhs) - std::get<2u>(rhs)
	// 	};
	// }

	// CartesianPoint operator+(CartesianPoint const& lhs, CartesianPoint const& rhs) {
	// 	return CartesianPoint{
	// 		std::get<0u>(lhs) + std::get<0u>(rhs),
	// 		std::get<1u>(lhs) + std::get<1u>(rhs),
	// 		std::get<2u>(lhs) + std::get<2u>(rhs)
	// 	};
	// }

	// CartesianPoint & operator/=(CartesianPoint & point, float_type const scalar) {
	// 	std::get<0u>(point) /= scalar;
	// 	std::get<1u>(point) /= scalar;
	// 	std::get<2u>(point) /= scalar;
	// 	return point;
	// }

	// CartesianPoint operator/(CartesianPoint const& point, float_type const scalar) {
	// 	CartesianPoint newPoint = point;
	// 	return newPoint /= scalar;
	// }

	// CartesianPoint & operator*=(CartesianPoint & point, float_type const scalar) {
	// 	std::get<0u>(point) *= scalar;
	// 	std::get<1u>(point) *= scalar;
	// 	std::get<2u>(point) *= scalar;
	// 	return point;
	// }

	// CartesianPoint operator*(CartesianPoint const& point, float_type const scalar) {
	// 	CartesianPoint newPoint = point;
	// 	return newPoint *= scalar;
	// }

	// float_type euclideanLength(CartesianPoint const& point) {
	// 	return std::sqrt(dotProduct(point));
	// }

	// float_type dotProduct(CartesianPoint const& point) {
	// 	return std::get<0u>(point)*std::get<0u>(point) + std::get<1u>(point)*std::get<1u>(point) + std::get<2u>(point)*std::get<2u>(point);
	// }

	// float_type dotProduct(CartesianPoint const& lhs, CartesianPoint const& rhs) {
	// 	return std::get<0u>(lhs)*std::get<0u>(rhs) + std::get<1u>(lhs)*std::get<1u>(rhs) + std::get<2u>(lhs)*std::get<2u>(rhs);
	// }

	// CartesianPoint crossProduct(CartesianPoint const& lhs, CartesianPoint const& rhs){
	// 	return CartesianPoint {
	// 		std::get<1u>(lhs)*std::get<2u>(rhs) - std::get<2u>(lhs)*std::get<1u>(rhs),
	// 		std::get<2u>(lhs)*std::get<0u>(rhs) - std::get<0u>(lhs)*std::get<2u>(rhs),
	// 		std::get<0u>(lhs)*std::get<1u>(rhs) - std::get<1u>(lhs)*std::get<0u>(rhs)
	// 	};
	// }

	// CartesianPoint normalize(CartesianPoint const& a) {
	// 	return a / euclideanLength(a);
	// }

	// void insertDerivatives(scon::mathmatrix<float_type> & row, CartesianPoint const& der, std::size_t const index) {
	// }

	// BmatrixRowCreator::BmatrixRowCreator(std::size_t const size) : row{ std::make_unique<scon::mathmatrix<float_type>>(1, size) } {}

	// void BmatrixRowCreator::insertAtomDerivative(CartesianPoint const & der, std::size_t const index) {
	// 	auto begin = index * 3u;
	// 	(*row)(0, index) = std::get<0u>(der);
	// 	(*row)(0, index + 1u) = std::get<1u>(der);
	// 	(*row)(0, index + 2u) = std::get<2u>(der);
	// }
	// scon::mathmatrix<float_type>&& BmatrixRowCreator::getRow(){
	// 	return std::move(*row);
	// }
}