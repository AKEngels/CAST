#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_COMMONINTERNALCOORDINATESBUILDER_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_COMMONINTERNALCOORDINATESBUILDER_H_

#include"InternalCoordinatesBuilder.h"

namespace internals{

class CommonInternalCoordinatesBuilder : public InternalCoordinatesBuilder {
public:
	CommonInternalCoordinatesBuilder(ic_util::BondGraph const& graph, InternalCoordinates::CartesiansForInternalCoordinates & coordinates, std::vector<std::shared_ptr<InternalCoordinates::Rotator>> & rotators);
	virtual std::unique_ptr<InternalCoordinate> buildBondDistance(std::size_t const, std::size_t const) const override;
	virtual std::unique_ptr<InternalCoordinate> buildBondAngle(std::size_t const, std::size_t const, std::size_t const) const override;
	virtual std::unique_ptr<InternalCoordinate> buildDihedralAngle(std::size_t const, std::size_t const, std::size_t const, std::size_t const) const override;
	virtual std::unique_ptr<InternalCoordinate> buildTranslationX(std::vector<std::size_t> const&) const override;
	virtual std::unique_ptr<InternalCoordinate> buildTranslationY(std::vector<std::size_t> const&) const override;
	virtual std::unique_ptr<InternalCoordinate> buildTranslationZ(std::vector<std::size_t> const&) const override;
	virtual TwoInternals buildTranslationXY(std::vector<std::size_t> const&) const override;
	virtual TwoInternals buildTranslationXZ(std::vector<std::size_t> const&) const override;
	virtual TwoInternals buildTranslationYZ(std::vector<std::size_t> const&) const override;
	virtual ThreeInternals buildTranslationXYZ(std::vector<std::size_t> const&) const override;
	virtual std::unique_ptr<InternalCoordinate> buildRotationA(std::vector<std::size_t> const&) override;
	virtual std::unique_ptr<InternalCoordinate> buildRotationB(std::vector<std::size_t> const&) override;
	virtual std::unique_ptr<InternalCoordinate> buildRotationC(std::vector<std::size_t> const&) override;
	virtual TwoInternals buildRotationAB(std::vector<std::size_t> const&) override;
	virtual TwoInternals buildRotationAC(std::vector<std::size_t> const&) override;
	virtual TwoInternals buildRotationBC(std::vector<std::size_t> const&) override;
	virtual ThreeInternals buildRotationABC(std::vector<std::size_t> const&) override;
private:
	ic_util::BondGraph const& bondGraph;
	InternalCoordinates::CartesiansForInternalCoordinates & cartesians;
	std::vector<std::shared_ptr<InternalCoordinates::Rotator>> & rotatorContainer;
	std::size_t numberOfAtoms;
};

}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_COMMONINTERNALCOORDINATESBUILDER_H_