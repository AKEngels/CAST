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

  class InternalCoordinatesBase{
  public: // Must be public to make them accessible for the implementation outside of class definition
    using InternalVec = std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>;
    using CartesianType = InternalCoordinates::CartesiansForInternalCoordinates;
    using BondGraph = ic_util::Graph<ic_util::Node>;
    using IndexVec = std::vector<std::vector<std::size_t>>;
  
  public:
    virtual void buildCoordinates(CartesianType const& cartesians, BondGraph const& graph, IndexVec const& indexVec) = 0;
    virtual void appendCoordinates(InternalVec && newCoordinates) = 0;
    
  };


}

#endif // INTERNAL_COORDIANTES_BASE