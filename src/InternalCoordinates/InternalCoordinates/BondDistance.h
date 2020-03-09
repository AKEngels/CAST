#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_BONDDISTANCE_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_BONDDISTANCE_H_

#include "InternalCoordinate.h"
#include"Helpers.h"

namespace internals{

struct BondDistance : public InternalCoordinate {
	template<typename Atom>
	BondDistance(Atom const& atomOne, Atom const& atomTwo)
		: index_a_{ atomOne.atom_serial - 1u }, index_b_{ atomTwo.atom_serial - 1u },
		elem_a_{ atomOne.element }, elem_b_{ atomTwo.element },
		constrained_{ false /*constraint? *constraint : Config::get().constrained_internals.constrain_bond_lengths*/ }
	{}

	std::size_t index_a_;
	std::size_t index_b_;
	std::string elem_a_;
	std::string elem_b_;

	float_type value(Eigen::MatrixXd const& cartesians) const override;
	float_type difference(Eigen::MatrixXd const& newCoordinates, Eigen::MatrixXd const&  oldCoordinates) const override;
	std::pair<Eigen::Vector3d, Eigen::Vector3d> derivatives(Eigen::MatrixXd const& cartesians) const;
	Eigen::VectorXd derivativeVector(Eigen::MatrixXd const& cartesians) const override;
	float_type hessianGuess(Eigen::MatrixXd const& cartesians) const override;
	std::string info(Eigen::MatrixXd const& cartesians) const override;

	virtual bool hasIndices(std::vector<std::size_t> const& indices) const override;
	virtual std::vector<std::size_t> getIndices() const override;
	virtual void makeConstrained(std::shared_ptr<AbstractConstraintManager> manager) override;
	virtual void makeConstrained() override { constrained_ = true; }
	virtual void releaseConstraint() override { constrained_ = false; }

	bool operator==(BondDistance const&) const;

	bool constrained_;
	virtual bool is_constrained() const override { return constrained_; }

private:
	bool bothElementsInPeriodOne(ic_util::period const atomA, ic_util::period const atomB)const;
	bool oneElementInPeriodOneTheOtherInPeriodTwo(ic_util::period const atomA, ic_util::period const atomB)const;
	bool oneElementInPeriodOneTheOtherInPeriodThree(ic_util::period const atomA, ic_util::period const atomB)const;
	bool bothElementsInPeriodTwo(ic_util::period const atomA, ic_util::period const atomB)const;
	bool oneElementInPeriodTwoTheOtherInPeriodThree(ic_util::period const atomA, ic_util::period const atomB)const;
};

}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_BONDDISTANCE_H_