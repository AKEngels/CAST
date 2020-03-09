#include "DihedralAngle.h"

#include "../InternalCoordinateUtilities.h"
#include "../../Scon/scon_mathmatrix.h"

namespace internals{

float_type DihedralAngle::value(Eigen::MatrixXd const& cartesians) const {
	//Eigentlich sollte das hier richtig sein:
	auto a = cartesians.row(index_a_);
	auto b = cartesians.row(index_b_);
	auto c = cartesians.row(index_c_);
	auto d = cartesians.row(index_d_);

	Eigen::Vector3d b1 = b - a;
	b1 /= b1.norm();
	Eigen::Vector3d b2 = c - b;
	b2 /= b2.norm();
	Eigen::Vector3d b3 = d - c;
	b3 /= b3.norm();
	auto n1 = b1.cross(b2);
	auto n2 = b2.cross(b3);
	return std::atan2(b1.dot(n2), n1.dot(n2));

}

float_type DihedralAngle::difference(Eigen::MatrixXd const& newCoordinates, Eigen::MatrixXd const& oldCoordinates) const {
	auto diff = value(newCoordinates) - value(oldCoordinates);
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

std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>
DihedralAngle::der(Eigen::MatrixXd const& cartesians) const {

	auto a = cartesians.row(index_a_);
	auto b = cartesians.row(index_b_);
	auto c = cartesians.row(index_c_);
	auto d = cartesians.row(index_d_);

	Eigen::Vector3d u = a - b;
	Eigen::Vector3d w = c - b;
	Eigen::Vector3d v = d - c;

	auto ul = u.norm();
	auto wl = w.nomr();
	auto vl = v.norm();

	auto u /= ul;
	auto w /= wl;
	auto v /= vl;
	auto sin2_u = 1. - std::pow(u.dot(w), 2.);
	auto sin2_v = 1. - std::pow(v.dot(w), 2.);

	auto t1 = u.cross(w) / (ul * sin2_u);
	auto t2 = v.cross(w) / (vl * sin2_v);
	auto cuwd = u.cross(w) * u.dot(w);
	auto t3 = cuwd / (wl * sin2_u);
	//rechanged it: Lee-Ping pointed out that his numeric evaluated derivatives match with analytics without this sign ---> changed v's sign according to J. Chem. Ohys Vol 117 No. 20, 22 2002 p. 9160-9174
	auto cvwd = v.cross(w) * v.dot(w);//auto cvwd = cross(v, w) * dot(-v, w);
	auto t4 = cvwd / (wl * sin2_v);
	//exchanged +t2 with +t3 in o's derivative
	//                      a   b               c             d
	//                      m   o               p             n
	return std::make_tuple(t1, -t1 + t3 - t4, t2 - t3 + t4, -t2);
}

Eigen::VectorXd
DihedralAngle::der_vec(Eigen::MatrixXd const& cartesians) const {

	auto firstDerivatives = derivatives(cartesians);

	Eigen::VectorXd result = Eigen::VectorXd::Zero(cartesians.size());

	Eigen::Map<RowMajorMatrixXd> mappedResult(result.data(), cartesians.rows(), cartesians.cols());

	mappedResult.row(index_a_) = std::get<0u>(firstder);
	mappedResult.row(index_b_) = std::get<1u>(firstder);
	mappedResult.row(index_c_) = std::get<2u>(firstder);
	mappedResult.row(index_d_) = std::get<3u>(firstder);

	return result;
}

float_type DihedralAngle::hessian_guess(Eigen::MatrixXd const & /*cartesians*/) const {
	return 0.023;
}

std::string DihedralAngle::info(Eigen::MatrixXd const & cartesians) const {
	std::ostringstream oss;
	oss << "Dihedral: " << value(cartesians) * SCON_180PI << " || " << index_a_ + 1u << " || " << index_b_ + 1u << " || " << index_c_ + 1u << " || " << index_d_ + 1u << " || " << "Constrained: " << std::boolalpha << is_constrained();
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