#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_ROTATOROBSERVER_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_ROTATOROBSERVER_H_

#include<memory>

#include "../PlainCoordinates/GeometryObserver.h"
#include "RotatorListener.h"

namespace internals{

class RotatorObserver : public GeometryObserver {
public:
	virtual ~RotatorObserver() = default;
	void setNewRotator(std::shared_ptr<RotatorListener> const rotator);
	void notify() override;
private:
	std::shared_ptr<RotatorListener> rotator;
};


void RotatorObserver::setNewRotator(std::shared_ptr<RotatorListener> const rotator) { this->rotator = rotator; }

void RotatorObserver::notify() {
	rotator->setAllFlag();
}

}

#endif // CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_ROTATOROBSERVER_H_