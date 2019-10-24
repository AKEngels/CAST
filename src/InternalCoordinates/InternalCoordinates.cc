#include "InternalCoordinates.h"

#include "../Scon/scon_mathmatrix.h"

#include "InternalCoordinateUtilities.h"
#include "BestFitRotation.h"
#include "../helperfunctions.h"
#include "Init/ConstraintManager.h"
#include "BondGraph/BondGraph.h"

namespace InternalCoordinates {

	



  /*bool operator==(BondDistance const & lhs, BondDistance const & rhs) {
    return lhs.operator==(rhs);
  }*/

  

  

  coords::float_type OutOfPlane::hessian_guess(coords::Representation_3D const & cartesians) const
  {
    auto const& a = cartesians.at(index_a_);
    auto const& b = cartesians.at(index_b_);
    auto const& c = cartesians.at(index_c_);
    auto const& d = cartesians.at(index_d_);
    auto r1 = b - a;
    auto r2 = b - c;
    auto r3 = b - d;
    auto r2Xr3 = cross(r2, r3);
    auto rd = dot(r1, r2Xr3);
    auto t2 = rd / (len(r1) * len(r2) * len(r3));
    auto dd = 1. - t2;
    auto d_pow = std::pow(dd, 4);
    return 0.045 * d_pow;
  }

  std::string OutOfPlane::info(coords::Representation_3D const & cartesians) const
  {
    std::ostringstream oss;
    oss << "Out of plane: " << val(cartesians) * SCON_180PI << "||" << index_a_ + 1u << "||" << index_b_ + 1u << "||" << index_c_ + 1u << "||" << index_d_ + 1u << " || " << "Constrained: " << std::boolalpha << is_constrained();
    return oss.str();
  }

  




  
  


	

}


