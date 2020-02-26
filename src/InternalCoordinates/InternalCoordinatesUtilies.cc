#include"InternalCoordinateUtilities.h"
#include"Matrix/BestFitRotation.h"


namespace internals {

std::pair<float_type, float_type> displacementRmsValAndMaxTwoStructures(AbstractMatrix const& lhs, AbstractMatrix const& rhs) {
	auto q = getQuaternion(lhs, rhs);
	auto U = getRotationMatrix(q.second);

	auto new_xyz_mat = toBaryCenter(rhs);
	auto old_xyz_mat = toBaryCenter(lhs);
	auto rot = (U*new_xyz_mat.transpose()).transpose();

	old_xyz_mat -= rot;
	old_xyz_mat *= -energy::bohr2ang;
	auto norms = ic_util::atomsNorm(old_xyz_mat);

	return { norms.rmsd(), norms.max() };
}

}