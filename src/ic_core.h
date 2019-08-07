#ifndef cast_ic_core_h_guard
#define cast_ic_core_h_guard

#include "coords.h"

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
