#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_DIHEDRALANGLE_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_DIHEDRALANGLE_H_

#include "InternalCoordinate.h"
#include "Helpers.h"

namespace internals{

struct DihedralAngle : public InternalCoordinate {

	template<typename Atom>
	DihedralAngle(Atom const& outerLeftAtom, Atom const& leftAtom,
		Atom const& rightAtom, Atom const& outerRightAtom)
		: index_a_{ outerLeftAtom.atom_serial - 1u },
		index_b_{ leftAtom.atom_serial - 1u }, index_c_{ rightAtom.atom_serial - 1u }, index_d_{ outerRightAtom.atom_serial - 1u },
		constrained_{ false }{}
	virtual ~DihedralAngle() = default;

	float_type val(scon::mathmatrix<float_type> const& cartesians) const override;
	float_type difference(scon::mathmatrix<float_type> const& newCoordinates, scon::mathmatrix<float_type> const& oldCoordinates) const override;
	std::tuple<CartesianPoint, CartesianPoint, CartesianPoint , CartesianPoint>
		der(scon::mathmatrix<float_type> const& cartesians) const;
	scon::mathmatrix<float_type> der_vec(scon::mathmatrix<float_type> const& cartesians) const override;
	float_type hessian_guess(scon::mathmatrix<float_type> const& cartesians) const override;
	std::string info(scon::mathmatrix<float_type> const& cartesians) const override;

	virtual bool hasIndices(std::vector<std::size_t> const& indices) const override;
	virtual std::vector<std::size_t> getIndices() const override;
	virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) override;
	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

	std::size_t index_a_;
	std::size_t index_b_;
	std::size_t index_c_;
	std::size_t index_d_;

	bool operator==(DihedralAngle const&) const;

	bool constrained_;
	virtual bool is_constrained() const override { return constrained_; }
};

}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_DIHEDRALANGLE_H_