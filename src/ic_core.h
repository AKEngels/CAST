#ifndef cast_ic_core_h_guard
#define cast_ic_core_h_guard

#include "coords.h"
#include "coords_rep.h"
#include "ic_atom.h"
#include "ic_rotation.h"
#include "pdb.h"
#include "scon_angle.h"
#include "scon_spherical.h"
#include "scon_vect.h"

#include <algorithm>
#include "scon_mathmatrix.h"
#include <array>
#include <boost/graph/adjacency_list.hpp>
#include <cmath>
#include <iterator>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <memory>

#include"InternalCoordinates.h"

#include "coords_io_pdb.h"

namespace coords {
  template<typename CoordKind>
  class DL_Coordinates : public Coordinates {
  public:
    std::shared_ptr<typename CoordKind::helper::template Parser<double>> parser;
    DL_Coordinates(Coordinates const& coords, std::unique_ptr<CoordKind> format) : Coordinates(coords) {
      if (!format) {
        throw std::runtime_error("You need to pass a format with fragments i. e. a pdb format.\n");
      }
      parser = format->parser;
    }
  };
}
namespace ic_core {

using coords::float_type;

inline coords::Representation_3D grads_to_bohr(coords::Representation_3D const& grads) {
  coords::Representation_3D bohr_grads;
  bohr_grads.reserve(grads.size());
  for (auto const& g : grads) {
    bohr_grads.emplace_back(g / energy::Hartree_Bohr2Kcal_MolAng);
  }
  return bohr_grads;
}
inline coords::Representation_3D rep3d_bohr_to_ang(coords::Representation_3D const& bohr) {
  coords::Representation_3D ang;
  ang.reserve(bohr.size());
  for (auto const& b : bohr) {
    ang.emplace_back(b * energy::bohr2ang);
  }
  return ang;
}









}
#endif // cast_ic_core_h_guard
