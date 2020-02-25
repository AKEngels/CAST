#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_BONDANGLE_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_BONDANGLE_H_

#include"InternalCoordinate.h"
#include"Helpers.h"

namespace internals{

struct BondAngle : InternalCoordinate {
	template<typename Atom>
	BondAngle(Atom const& leftAtom, Atom const& middleAtom, Atom const& rightAtom)
		: index_a_{ leftAtom.atom_serial - 1u }, index_b_{ middleAtom.atom_serial - 1u },
		index_c_{ rightAtom.atom_serial - 1u }, elem_a_{ leftAtom.element }, elem_b_{ middleAtom.element },
		elem_c_{ rightAtom.element },
		constrained_{ false }
	{}

	std::size_t index_a_;
	std::size_t index_b_;
	std::size_t index_c_;
	std::string elem_a_;
	std::string elem_b_;
	std::string elem_c_;

	virtual float_type val(scon::mathmatrix<float_type> const& cartesians) const override;
	virtual float_type difference(scon::mathmatrix<float_type> const& newCoordinates, scon::mathmatrix<float_type> const& oldCoordinates) const override;
	std::tuple<CartesianPoint, CartesianPoint, CartesianPoint> der(scon::mathmatrix<float_type> const& cartesians) const;
	virtual scon::mathmatrix<float_type> der_vec(scon::mathmatrix<float_type> const& cartesians) const override;
	virtual float_type hessian_guess(scon::mathmatrix<float_type> const& cartesians) const override;
	virtual std::string info(scon::mathmatrix<float_type> const& cartesians) const override;

	virtual bool hasIndices(std::vector<std::size_t> const& indices) const override;
	virtual std::vector<std::size_t> getIndices() const override;
	virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) override;
	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

	bool operator==(BondAngle const&) const;

	bool constrained_;
	virtual bool is_constrained() const override { return constrained_; }
};

}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_BONDANGLE_H_