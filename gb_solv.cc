#include "gb_solv.h"

//coords::float_type std_solv_rad(size_t const atomic_number, size_t const bonds, size_t const first_bound_atomic)
//{
//  switch (atomic_number)
//  {
//  case 1:
//  {
//    if (bonds == 1)
//    {
//      if (first_bound_atomic == 7) return coords::float_type(1.15);
//      else if (first_bound_atomic == 8) return coords::float_type(1.05);
//    }
//  }
//  case 6:
//  {
//    if (bonds == 3) return coords::float_type(1.875);
//    else if (bonds == 2) return coords::float_type(1.825);
//  }
//  case 7:
//  {
//    if (bonds == 4) return coords::float_type(1.625);
//    else if (bonds == 1) return coords::float_type(1.6);
//  }
//  case 8:
//  {
//    if (bonds == 1) return coords::float_type(1.48);
//  }
//  };
//  if (atomic_number < ATOMIC_H) return atomic::radiusMap[atomic_number];
//  else return 2.0;
//}
//
//coords::float_type gbsa::gborn::solv_rad(size_t const atom_index) const
//{
//  if (Config::get().gbsa.radius_type == config::gbsa_conf::radius_types::STD)
//  {
//    return std_solv_rad(coord_ptr->atoms(atom_index).number(),
//      coord_ptr->atoms(atom_index).bonds().size(),
//      coord_ptr->atoms(atom_index).bonds().empty() ? 0u : coord_ptr->atoms(atom_index).bonds(0u));
//  }
//  else
//  {
//    return param_ptr->vdws()[param_ptr->type(atom_index, tinker::VDW)].r;
//  }
//}
//
//void gbsa::gborn::solv_radii()
//{
//  size_t const n = coord_ptr ? 0u : coord_ptr->size();
//  solvation_radius.resize(n);
//  for (size_t i = 0; i < n; ++i)
//  {
//    solvation_radius[i] = solv_rad(i);
//  }
//}
//
//
//
//void gbsa::gborn::sa_tinker()
//{
//  size_t const n = coord_ptr ? 0u : coord_ptr->size();
//  solvent_access.resize(n);
//  for (size_t i = 0; i < n; ++i)
//  {
//    coords::float_type p(solvation_radius[i] + probe);
//    p *= p;
//    coords::float_type q(solvation_radius[i] / born_radius[i]);
//    q = q*q*q; // q^3
//    q *= q; // q^6
//    solvent_access[i] = coords::float_type(4.)*GBSA_PI*p*q;
//  }
//}

/*

  Still helper functions

*/

//coords::float_type pi(coords::Coordinates const * const , size_t )
//{
//  //size_t const atomic_i = coord_ptr->atoms(i).number();
//  //if (atomic_i == 1)
//  //{
//  //  coords::float_type ret = 0.9822;
//  //  size_t const bound = coord_ptr->atoms(i).bonds(0);
//  //  if (coord_ptr->atoms(bound).number() == 7)
//  //  {
//  //    size_t hn = 0;
//  //    for (size_t j = 0; j < coord_ptr->atoms(i).bonds().size(); ++j)
//  //    {
//  //      if (coord_ptr->atoms(coord_ptr->atoms(i).bonds(j)).number() == 1) ++hn;
//  //      if (coord_ptr->atoms(coord_ptr->atoms(i).bonds(j)).number() == 8)
//  //      {
//  //        return 0.9808;
//  //      }
//  //    }
//  //    if (hn == 3) return 1.1170;
//  //    return 1.2245;
//  //  }
//  //  return 0.9564;
//  //}
//}
//
//coords::float_type pij(coords::Coordinates const * const cp, size_t i, size_t j)
//{
//  size_t const bi = cp->atoms(i).bonds().size();
//  size_t const bj = cp->atoms(i).bonds().size();
//  for (size_t a = 0; a < bi; ++a)
//  {
//    if (cp->atoms(i).bonds(a) == j) return coords::float_type(0.95);
//    for (size_t b = 0; b < bj; ++b)
//    {
//      if (cp->atoms(i).bonds(a) == cp->atoms(j).bonds(b)) return coords::float_type(0.20);
//    }
//  }
//  return coords::float_type(0.35);
//}
//
//coords::float_type pijk(coords::Coordinates const * const cp, size_t i, size_t j, coords::float_type const d_ij2)
//{
//  size_t const bi = cp->atoms(i).bonds().size();
//  coords::float_type ret = coords::float_type();
//  for (size_t a = 0; a < bi; ++a)
//  {
//    size_t const k = cp->atoms(i).bonds(a);
//    if (k != j)
//    {
//      coords::Cartesian_Point const jk = cp->xyz(j) - cp->xyz(k);
//      coords::Cartesian_Point const ij = cp->xyz(i) - cp->xyz(j);
//      ret += jk.scalar(jk) / d_ij2;
//    }
//  }
//  return ret;
//}
//
//void gbsa::gborn::sa_still()
//{
//  size_t const n = coord_ptr ? 0u : coord_ptr->size();
//  solvent_access.resize(n);
//  std::vector<coords::float_type> const & r = solvation_radius;
//  for (size_t i = 0; i < n; ++i)
//  {
//    solvent_access[i] = coords::float_type(1.);
//    coords::float_type const rip = r[i] + probe;
//    coords::float_type const Si = coords::float_type(4.)*GBSA_PI*rip*rip;
//    for (size_t j = 0; j < n; ++j)
//    {
//      coords::Cartesian_Point const ij = coord_ptr->xyz(i) - coord_ptr->xyz(j);
//      coords::float_type const rij2 = ij.scalar(ij);
//      coords::float_type const rij = std::sqrt(rij2);
//      coords::float_type const or = rip + r[j] + probe;
//      if (rij < or)
//      {
//        coords::float_type const bij = GBSA_PI * rip *(or - rij)*((coords::float_type(1.) + r[j] - r[i]) / rij);
//        solvent_access[i] *= (coords::float_type(1.) - pi(coord_ptr, i)*pij(coord_ptr, i, j)*pijk(coord_ptr, i, j, rij2)*bij) / Si;
//      }
//    }
//  }
//}