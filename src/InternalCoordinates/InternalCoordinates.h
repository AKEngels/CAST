#ifndef H_INTERNAL_COORDINATES
#define H_INTERNAL_COORDINATES

#include<array>

#include"InternalCoordinatesAliases.h"
#include "BondGraph/ElementInformations.h"
#include "../coords.h"

namespace InternalCoordinates {

  // struct OutOfPlane : public DihedralAngle {
  //   template <typename Atom>
  //   OutOfPlane(Atom const& outerLeftAtom, Atom const& leftAtom,
  //     Atom const& rightAtom, Atom const& outerRightAtom)
  //   : DihedralAngle{ outerLeftAtom, leftAtom, rightAtom, outerRightAtom }
  //   {}
    
  //   using DihedralAngle::DihedralAngle;

  //   coords::float_type hessian_guess(coords::Representation_3D const& cartesians) const override;
  //   std::string info(coords::Representation_3D const& cartesians) const override;

  //   virtual bool is_constrained() const override {return constrained_;}
  // };
	
}

#endif
