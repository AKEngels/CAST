#include "DihedralAngle.h"

#include "../InternalCoordinateUtilities.h"
#include "../../Scon/scon_mathmatrix.h"

namespace internals{

float_type DihedralAngle::val(scon::mathmatrix<float_type> const& cartesians) const {
	//Eigentlich sollte das hier richtig sein:
	auto a = getAtom(cartesians, index_a_);
	auto b = getAtom(cartesians, index_b_);
	auto c = getAtom(cartesians, index_c_);
	auto d = getAtom(cartesians, index_d_);

	auto b1 = b - a;
	b1 /= euclideanLength(b1);
	auto b2 = c - b;
	b2 /= euclideanLength(b2);
	auto b3 = d - c;
	b3 /= euclideanLength(b3);
	auto n1 = crossProduct(b1, b2);
	auto n2 = crossProduct(b2, b3);
	return std::atan2(dotProduct(b1, n2), dotProduct(n1, n2));

}

float_type DihedralAngle::difference(scon::mathmatrix<float_type> const& newCoordinates, scon::mathmatrix<float_type> const& oldCoordinates) const {
	auto diff = val(newCoordinates) - val(oldCoordinates);
	if (std::fabs(diff) > SCON_PI) {
		if (diff < 0.0) {
			diff += 2. * SCON_PI;
		}
		else {
			diff -= 2. * SCON_PI;
		}
	}

	return diff;
}

std::tuple<CartesianPoint, CartesianPoint, CartesianPoint, CartesianPoint>
DihedralAngle::der(scon::mathmatrix<float_type> const& cartesians) const {
	using scon::cross;
	using scon::dot;
	using scon::len;

	auto const & a = getAtom(cartesians, index_a_);
	auto const & b = getAtom(cartesians, index_b_);
	auto const & c = getAtom(cartesians, index_c_);
	auto const & d = getAtom(cartesians, index_d_);

	auto u_p = a - b;
	auto w_p = c - b;
	auto v_p = d - c;
	auto u = normalize(u_p);
	auto w = normalize(w_p);
	auto v = normalize(v_p);
	auto sin2_u = 1. - std::pow(dotProduct(u, w), 2.);
	auto sin2_v = 1. - std::pow(dotProduct(v, w), 2.);

	auto t1 = crossProduct(u, w) / (euclideanLength(u_p) * sin2_u);
	auto t2 = crossProduct(v, w) / (euclideanLength(v_p) * sin2_v);
	auto cuwd = crossProduct(u, w) * dotProduct(u, w);
	auto t3 = cuwd / (euclideanLength(w_p) * sin2_u);
	//rechanged it: Lee-Ping pointed out that his numeric evaluated derivatives match with analytics without this sign ---> changed v's sign according to J. Chem. Ohys Vol 117 No. 20, 22 2002 p. 9160-9174
	auto cvwd = crossProduct(v, w) * dotProduct(v, w);//auto cvwd = cross(v, w) * dot(-v, w);
	auto t4 = cvwd / (euclideanLength(w_p) * sin2_v);
	//exchanged +t2 with +t3 in o's derivative
	//                      a   b               c             d
	//                      m   o               p             n
	return std::make_tuple(t1, -t1 + t3 - t4, t2 - t3 + t4, -t2);
}

scon::mathmatrix<float_type>
DihedralAngle::der_vec(scon::mathmatrix<float_type> const& cartesians) const {
	using scon::c3;


	auto firstder = der(cartesians);

	BmatrixRowCreator rowCreator(cartesians.cols()*cartesians.rows());
	rowCreator.insertAtomDerivative(std::get<0u>(firstder), index_a_);
	rowCreator.insertAtomDerivative(std::get<1u>(firstder), index_b_);
	rowCreator.insertAtomDerivative(std::get<2u>(firstder), index_c_);
	rowCreator.insertAtomDerivative(std::get<2u>(firstder), index_d_);

	// TODO: test if this line is mandatory to trigger return optimization
	scon::mathmatrix<float_type> result = rowCreator.getRow();

	return result;
}

float_type DihedralAngle::hessian_guess(scon::mathmatrix<float_type> const & /*cartesians*/) const {
	return 0.023;
}

std::string DihedralAngle::info(scon::mathmatrix<float_type> const & cartesians) const {
	std::ostringstream oss;
	oss << "Dihedral: " << val(cartesians) * SCON_180PI << " || " << index_a_ + 1u << " || " << index_b_ + 1u << " || " << index_c_ + 1u << " || " << index_d_ + 1u << " || " << "Constrained: " << std::boolalpha << is_constrained();
	return oss.str();
}

bool DihedralAngle::hasIndices(std::vector<std::size_t> const& indices) const {
	return indices.size() == 3u && ((index_a_ == indices[0u] && index_b_ == indices[1u] && index_c_ == indices[2u] && index_d_ == indices[3u]) || (index_a_ == indices[3u] && index_b_ == indices[2u] && index_c_ == indices[1u] && index_d_ == indices[0u]));
}

std::vector<std::size_t> DihedralAngle::getIndices() const {
	return { index_a_, index_b_, index_c_, index_d_ };
}

void DihedralAngle::makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) {
	auto constraint = manager->checkForDihedrals(*this);
	if (constraint) {
		if (constraint->isFrozen()) makeConstrained();
		else releaseConstraint();
	}

}

bool DihedralAngle::operator==(DihedralAngle const & other) const {
	return index_a_ == other.index_a_ && index_b_ == other.index_b_ && index_c_ == other.index_c_ && index_d_ == other.index_d_;
}

}