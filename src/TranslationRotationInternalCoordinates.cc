#include "TranslationRotationInternalCoordinates.h"
#include "InternalCoordinateDecorator.h"

#include "Scon/scon_mathmatrix.h"

namespace internals {

	TRIC::TRIC(ICDecoratorBase& decorator, const CartesianType& cartesians) :
		PrimitiveInternalCoordinates{ decorator }, del_mat{ std::make_unique<scon::mathmatrix<coords::float_type>>() } {
		delocalize_ic_system(cartesians);
	}


	scon::mathmatrix<coords::float_type>& TRIC::delocalize_ic_system(CartesianType const& cartesians) {
		using Mat = scon::mathmatrix<coords::float_type>;

		Mat eigval, eigvec;
		std::tie(eigval, eigvec) = PrimitiveInternalCoordinates::Gmat(cartesians).eigensym(false);

		auto row_index_vec = eigval.sort_idx();
		auto col_index_vec = eigval.find_idx([](coords::float_type const& a) {
			return std::abs(a) > 1e-6;
		});

		*del_mat = eigvec.submat(row_index_vec, col_index_vec);
		new_B_matrix = new_G_matrix = true; //B and G got to be calculated for the new ic_system
		return *del_mat;
	}


	scon::mathmatrix<coords::float_type>& TRIC::Bmat(CartesianType const& cartesians) {
		//TODO activate it again!!!!!
		/*if (!new_B_matrix) {
			return B_matrix;
		}*/
		*B_matrix = del_mat->t() * PrimitiveInternalCoordinates::Bmat(cartesians);
		new_B_matrix = false;
		return *B_matrix;
	}

	scon::mathmatrix<coords::float_type> TRIC::transposeOfBmat(CartesianType const& cartesian) {
		return Bmat(cartesian).t();
	}

	scon::mathmatrix<coords::float_type> TRIC::pseudoInverseOfGmat(CartesianType const& cartesian) {
		return Gmat(cartesian).pinv();
	}

	scon::mathmatrix<coords::float_type>& TRIC::Gmat(CartesianType const& cartesians) {
		//TODO activate it again!!!!!
		/*if (!new_G_matrix) {
			return G_matrix;
		}*/
		Bmat(cartesians);
		*G_matrix = *B_matrix * B_matrix->t();
		new_G_matrix = false;
		return *G_matrix;
	}

	scon::mathmatrix<coords::float_type> TRIC::guess_hessian(CartesianType const& cartesians) const {
		return del_mat->t() * PrimitiveInternalCoordinates::guess_hessian(cartesians) * (*del_mat);
	}

	scon::mathmatrix<coords::float_type> TRIC::calc(CartesianType const& xyz) const {
		auto prims = PrimitiveInternalCoordinates::calc(xyz);
		return (prims * (*del_mat)).t();
	}

	scon::mathmatrix<coords::float_type> TRIC::calc_diff(CartesianType const& lhs, CartesianType const& rhs) const {
		auto diff = PrimitiveInternalCoordinates::calc_diff(lhs, rhs);
		return diff * (*del_mat);
	}

	scon::mathmatrix<coords::float_type> TRIC::projectorMatrix(CartesianType const& /*cartesian*/) {
		auto s = getDelMat().cols();
		return scon::mathmatrix<coords::float_type>::identity(s, s);
	}

	scon::mathmatrix<coords::float_type> const& TRIC::getDelMat() const { return *del_mat; }
}
