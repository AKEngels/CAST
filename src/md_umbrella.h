#pragma once

#include "md.h"
namespace md
{
  class CoordinatesUBIAS {
  public:
    CoordinatesUBIAS(coords::Coordinates* const coordinates) : coordinates(coordinates), broken_bonds() {}
    //umbrella
    void ubias(std::vector<double>& uout)
    {
      if (!coordinates->m_potentials.uempty())
        coordinates->m_potentials.umbrellaapply(coordinates->m_representation.structure.cartesian,
          coordinates->m_representation.gradient.cartesian,
          uout);
    }
  
    // very ugly approach to get the correct return type. No other way seen yet
    decltype(std::declval<coords::Coordinates>().atoms(std::declval<int>()))
      atoms(int const i) const {
      return coordinates->atoms(i);
    }
    decltype(std::declval<coords::Coordinates>().size())
      size() const {
      return coordinates->size();
    }
    decltype(std::declval<coords::Coordinates>().g())
      g() {
      return coordinates->g();
    }
    decltype(std::declval<coords::Coordinates>().xyz(std::declval<int>()))
      xyz(int const i) const {
      return coordinates->xyz(i);
    }
    decltype(std::declval<coords::Coordinates>().xyz())
      xyz() const {
      return coordinates->xyz();
    }
    bool validate_bonds();
    decltype(std::declval<coords::Coordinates>().getFep())
      getFep() const {
      return coordinates->getFep();
    }
    std::vector<std::vector<float>> const& getBrokenBonds() const { return broken_bonds; }
    auto center_of_mass() const {
      return coordinates->center_of_mass();
    }
    auto center_of_geometry() const {
      return coordinates->center_of_geometry();
    }
    template<typename ... Args>
    auto move_all_by(Args ... args) {
      return coordinates->move_all_by(args...);
    }
    template<typename ... Args>
    void move_atom_by(Args ... args) {
      coordinates->move_atom_by(args...);
    }
    template<typename ... Args>
    void move_atom_to(Args ... args) {
      coordinates->move_atom_to(args...);
    }
    /** add gradients to an atom (spherical boundaries)
    @param index: atom index
    @param g: gradients that should be added to the gradients of index*/
    void add_sp_gradients(std::size_t const index, coords::Cartesian_Point const& g)
    {
      coordinates->m_representation.gradient.cartesian[index] += g;
    }
    /**scale the coordinates of an atom (used for pressure control)
    @param index: index of atom that is to be moved
    @param p: factor by which the coordinates should be scaled
    @param force_move: if set to true also move fixed atoms*/
    void scale_atom_by(std::size_t const index, double& p, bool const force_move = false)
    {
      if (!atoms(index).fixed() || force_move)
      {
        coordinates->m_representation.structure.cartesian[index] *= p;
        coordinates->energy_valid = false;
      }
      coordinates->m_stereo.update(xyz());
    }
    /** set new atom coordinates if anisotropic pressure control is enabled*/
    void set_atom_aniso(std::size_t const index, double tvec[3][3], bool const force_move = false) {
      if (!atoms(index).fixed() || force_move)
      {
        double tcor;
        tcor = tvec[0][0] * coordinates->m_representation.structure.cartesian[index].x() +
          tvec[0][1] * coordinates->m_representation.structure.cartesian[index].y() + tvec[0][2] *
          coordinates->m_representation.structure.cartesian[index].z();
        coordinates->m_representation.structure.cartesian[index].x() = tcor;
        tcor = tvec[1][0] * coordinates->m_representation.structure.cartesian[index].x() +
          tvec[1][1] * coordinates->m_representation.structure.cartesian[index].y() + tvec[1][2] *
          coordinates->m_representation.structure.cartesian[index].z();
        coordinates->m_representation.structure.cartesian[index].y() = tcor;
        tcor = tvec[2][0] * coordinates->m_representation.structure.cartesian[index].x() +
          tvec[2][1] * coordinates->m_representation.structure.cartesian[index].y() + tvec[2][2] *
          coordinates->m_representation.structure.cartesian[index].z();
        coordinates->m_representation.structure.cartesian[index].z() = tcor;
        coordinates->energy_valid = false;
      }
      coordinates->m_stereo.update(xyz());
    }
    template<typename ... Args>
    void set_xyz(Args ... args) const {
      return coordinates->set_xyz(args...);
    }
    template<typename ... Args>
    auto g_xyz(Args ... args) const {
      return coordinates->g_xyz(args...);
    }
  
    /**updates the topology*/
    void energy_update(bool const skip_topology = false) { coordinates->energy_update(skip_topology); }
    /**returns matrix with interactions between subsystems*/
    coords::sub_ia_matrix_t& interactions()
    {
      return coordinates->interactions();
    }
    /**returns matrix with interactions between subsystems*/
    coords::sub_ia_matrix_t const& interactions() const
    {
      return coordinates->interactions();
    }
    /**adds V to virial coefficients
    @param V: values that should be added*/
    void add_to_virial(std::array<std::array<double, 3>, 3> & V)
    {
      coordinates->m_virial[0][0] += V[0][0];
      coordinates->m_virial[1][0] += V[1][0];
      coordinates->m_virial[2][0] += V[2][0];
      coordinates->m_virial[0][1] += V[0][1];
      coordinates->m_virial[1][1] += V[1][1];
      coordinates->m_virial[2][1] += V[2][1];
      coordinates->m_virial[0][2] += V[0][2];
      coordinates->m_virial[1][2] += V[1][2];
      coordinates->m_virial[2][2] += V[2][2];
    };
    /**returns the virial coefficients*/
    coords::virial_t const& virial() const { return coordinates->m_virial; }
    decltype(std::declval<coords::Coordinates>().pes())
      pes() const {
      return coordinates->pes();
    }
  protected:
    coords::Coordinates* const coordinates;
  
    /**vector of broken bonds (determined by validate_bonds())
    i.e. bondlength either too short or too long
    each element of the vector is a vector which contains the numbers of the two atoms that form the bond and the bond length*/
    std::vector<std::vector<float>> broken_bonds;
  };
}