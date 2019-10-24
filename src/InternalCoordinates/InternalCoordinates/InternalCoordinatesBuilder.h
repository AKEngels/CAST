#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_INTERNALCOORDINATESBUILDER_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_INTERNALCOORDINATESBUILDER_H_

#include<memory>
#include<tuple>
#include<vector>

#include"InternalCoordinate.h"

namespace internals{

class InternalCoordinatesBuilder {
public:
	using TwoInternals = std::pair<std::unique_ptr<InternalCoordinate>, std::unique_ptr<InternalCoordinate>>;
	using ThreeInternals = std::tuple<std::unique_ptr<InternalCoordinate>, std::unique_ptr<InternalCoordinate>, std::unique_ptr<InternalCoordinate>>;
	virtual std::unique_ptr<InternalCoordinate> buildBondDistance(std::size_t const, std::size_t const) const = 0;
	virtual std::unique_ptr<InternalCoordinate> buildBondAngle(std::size_t const, std::size_t const, std::size_t const) const = 0;
	virtual std::unique_ptr<InternalCoordinate> buildDihedralAngle(std::size_t const, std::size_t const, std::size_t const, std::size_t const) const = 0;
	virtual std::unique_ptr<InternalCoordinate> buildTranslationX(std::vector<std::size_t> const&) const = 0;
	virtual std::unique_ptr<InternalCoordinate> buildTranslationY(std::vector<std::size_t> const&) const = 0;
	virtual std::unique_ptr<InternalCoordinate> buildTranslationZ(std::vector<std::size_t> const&) const = 0;
	virtual TwoInternals buildTranslationXY(std::vector<std::size_t> const&) const = 0;
	virtual TwoInternals buildTranslationXZ(std::vector<std::size_t> const&) const = 0;
	virtual TwoInternals buildTranslationYZ(std::vector<std::size_t> const&) const = 0;
	virtual ThreeInternals buildTranslationXYZ(std::vector<std::size_t> const&) const = 0;
	virtual std::unique_ptr<InternalCoordinate> buildRotationA(std::vector<std::size_t> const&) = 0;
	virtual std::unique_ptr<InternalCoordinate> buildRotationB(std::vector<std::size_t> const&) = 0;
	virtual std::unique_ptr<InternalCoordinate> buildRotationC(std::vector<std::size_t> const&) = 0;
	virtual TwoInternals buildRotationAB(std::vector<std::size_t> const&) = 0;
	virtual TwoInternals buildRotationAC(std::vector<std::size_t> const&) = 0;
	virtual TwoInternals buildRotationBC(std::vector<std::size_t> const&) = 0;
	virtual ThreeInternals buildRotationABC(std::vector<std::size_t> const&) = 0;
	virtual ~InternalCoordinatesBuilder() = 0;
};

InternalCoordinatesBuilder::~InternalCoordinatesBuilder() = default;

}
#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_INTERNALCOORDINATESBUILDER_H_