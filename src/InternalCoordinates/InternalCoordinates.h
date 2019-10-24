#ifndef H_INTERNAL_COORDINATES
#define H_INTERNAL_COORDINATES

#include<array>

#include"InternalCoordinatesAliases.h"
#include "BondGraph/ElementInformations.h"
#include "../coords.h"

class AbstractConstraintManager;

namespace ic_util {
	enum class period;
}

namespace InternalCoordinates {

  

  

  
  
  

  inline coords::Representation_3D sliceCartesianCoordinates(CartesiansForInternalCoordinates const& cartesians, std::vector<std::size_t> const& indexVector) {
    coords::Representation_3D slicedCoordinates;
    for (auto const& index : indexVector) {
      slicedCoordinates.emplace_back(cartesians.at(index - 1));
    }
    return slicedCoordinates;
  }



  

  

  struct OutOfPlane : public DihedralAngle {
    template <typename Atom>
    OutOfPlane(Atom const& outerLeftAtom, Atom const& leftAtom,
      Atom const& rightAtom, Atom const& outerRightAtom)
    : DihedralAngle{ outerLeftAtom, leftAtom, rightAtom, outerRightAtom }
    {}
    
    using DihedralAngle::DihedralAngle;

    coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
    std::string info(coords::Representation_3D const& cartesians) const override;

    virtual bool is_constrained() const override {return constrained_;}
  };

  
  



	

	
}

#endif
