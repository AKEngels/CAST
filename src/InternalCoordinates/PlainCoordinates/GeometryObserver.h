#ifndef CAST_INTERNALCOORDINATES_PLAINCOORDINATES_GEOMETRYOBSERVER_H_
#define CAST_INTERNALCOORDINATES_PLAINCOORDINATES_GEOMETRYOBSERVER_H_

class GeometryObserver {
public:
	virtual ~GeometryObserver() = 0;
	virtual void notify() = 0;
};

GeometryObserver::~GeometryObserver() = default;

#endif // CAST_INTERNALCOORDINATES_PLAINCOORDINATES_GEOMETRYOBSERVER_H_