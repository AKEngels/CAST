/**
CAST 3
PrimitiveInternalsTransRot.h
Purpose: Extension of primitive Internal Coordinate System by translational and rotational coordinates


@author Julian Erdmannsdörfer, Michael Prem, Christian Schärf
@version 3.0
*/

#ifndef PRIMITIVE_INTERNALs_TRANS_ROT_H
#define PRIMITIVE_INTERNALs_TRANS_ROT_H

#include "PrimitiveInternalCoordinates.h"

namespace internals{
  class PrimitiveInternalsTransRot : public PrimitiveInternalCoordinates{
  public:
    PrimitiveInternalsTransRot (
      const std::vector<coords::Representation_3D>& res_init,
      const std::vector<std::vector<std::size_t>>& res_index,
      CartesianType & xyz_init,
      BondGraph const& graph
    );
    
    void create_translations_rotations(CartesianType & cartesians);
  
    virtual void prepare_rotations() const override;
    
  protected:
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> registeredRotators;
    
    InternalVec create_trans_x() const;
    InternalVec create_trans_y() const;
    InternalVec create_trans_z() const;
    
    //std::tuple<InternalVec, InternalVec, InternalVec>
    //  createRotationABC(std::vector<InternalCoordinates::Rotations> & rotations);
    std::tuple<InternalVec, InternalVec, InternalVec>
      create_rotations(CartesianType & cartesians);
      
    std::shared_ptr<InternalCoordinates::Rotator> build_rotation (
      InternalCoordinates::CartesiansForInternalCoordinates & target,
      std::vector<std::size_t> const& index_vec
    );
      
  };

}


#endif // PRIMITIVE_INTERNALs_TRANS_ROT_H
