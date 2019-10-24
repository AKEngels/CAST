#ifndef CAST_INTERNALCOORDINATES_GEOMETRYOPTIMIZATION_H_
#define CAST_INTERNALCOORDINATES_GEOMETRYOPTIMIZATION_H_

#include"InternalCoordinatesAliases.h"

class GeometryOptimization {
public:
	virtual 
};

class PrimitiveInternalCoordinatesOptimization : public GeometryOptimization {

};

class AbstractGeometryOptimizationFactory {
public:
	virtual ~AbstractGeometryOptimizationFactory() = 0;

	virtual CoordinateSystem
};

AbstractGeometryOptimizationFactory::~AbstractGeometryOptimizationFactory() = default;

class GeometryOptimizaionFromCoordsObj : public AbstractGeometryOptimizationFactory {
public:
	GeometryOptimizaionFromCoordsObj(coords::Coordinates const&);


};

#endif // CAST_INTERNALCOORDINATES_GEOMETRYOPTIMIZATION_H_