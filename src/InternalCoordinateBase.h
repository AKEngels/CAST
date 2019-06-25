/**
CAST 3
InternalCoordinateBase.h
Purpose: base (intreface) class for PrimitiveInternalCoordinates and its decorators


@author Julian Erdmannsdörfer, Christian Schärf
@version 3.0
*/

#ifndef INTERNAL_COORDINATES_BASE
#define INTERNAL_COORDINATES_BASE

#include "InternalCoordinates.h"
#include "graph.h"

namespace internals {
  // Must be here so it can be used by the Appender classes
  using InternalVec = std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>;
  
  class PrimitiveInternalCoordinates;
  
  class InternalCoordinateAppenderInterface{
  public:
    virtual void append(std::shared_ptr<PrimitiveInternalCoordinates> primitives) = 0;
  };

  class InternalCoordinatesBase{
  public: // Must be public to make them accessible for the implementation outside of class definition
    using CartesianType = InternalCoordinates::CartesiansForInternalCoordinates;
    using BondGraph = ic_util::Graph<ic_util::Node>;
    using IndexVec = std::vector<std::vector<std::size_t>>;
  
  public:
    // cartesians cannot be a const& because of InternalCoordinates::Rotator::buildRotator
    // (which is called by ICRotationDecorator::buildCoordinates which implements this method)
    virtual void buildCoordinates(CartesianType & cartesians, BondGraph const& graph, IndexVec const& indexVec) = 0;
    virtual void appendCoordinates(std::shared_ptr<InternalCoordinateAppenderInterface> appender) = 0;
    
  };


}

#endif // INTERNAL_COORDIANTES_BASE