#pragma once

#include "md.h"

namespace md
{
  /**special Coordinates class for MD simulations*/
  class CoordinatesUBIAS : public coords::Coordinates 
  {
  public:
    /**constructor (from "normal" Coordinates object*/
    CoordinatesUBIAS(coords::Coordinates* const coordinates) : coords::Coordinates(*coordinates), broken_bonds() {}
    
    /**add biased potential for umbrella sampling*/
    void ubias(std::vector<double>& uout, std::optional<Spline> const& s)
    {
      if (!m_potentials.uempty())
        m_potentials.umbrellaapply(m_representation.structure.cartesian,
          m_representation.gradient.cartesian,
          uout, s);
    }
  
    /**check if the lenths of bonds is still reasonable
    i. e. between 0.3 and 5 angstrom*/
    bool validate_bonds();

    /**return information about broken bonds*/
    std::vector<std::vector<float>> const& getBrokenBonds() const { return broken_bonds; }

    /**scale the coordinates of an atom (used for pressure control)
    @param index: index of atom that is to be moved
    @param p: factor by which the coordinates should be scaled
    @param force_move: if set to true also move fixed atoms*/
    void scale_atom_by(std::size_t const index, double& p, bool const force_move = false)
    {
      if (!atoms(index).fixed() || force_move)
      {
        m_representation.structure.cartesian[index] *= p;
        energy_valid = false;
      }
      m_stereo.update(xyz());
    }
    /** set new atom coordinates if anisotropic pressure control is enabled*/
    void set_atom_aniso(std::size_t const index, double tvec[3][3], bool const force_move = false) {
      if (!atoms(index).fixed() || force_move)
      {
        double tcor;
        tcor = tvec[0][0] * m_representation.structure.cartesian[index].x() +
          tvec[0][1] * m_representation.structure.cartesian[index].y() + tvec[0][2] *
          m_representation.structure.cartesian[index].z();
        m_representation.structure.cartesian[index].x() = tcor;
        tcor = tvec[1][0] * m_representation.structure.cartesian[index].x() +
          tvec[1][1] * m_representation.structure.cartesian[index].y() + tvec[1][2] *
          m_representation.structure.cartesian[index].z();
        m_representation.structure.cartesian[index].y() = tcor;
        tcor = tvec[2][0] * m_representation.structure.cartesian[index].x() +
          tvec[2][1] * m_representation.structure.cartesian[index].y() + tvec[2][2] *
          m_representation.structure.cartesian[index].z();
        m_representation.structure.cartesian[index].z() = tcor;
        energy_valid = false;
      }
      m_stereo.update(xyz());
    }

    /**adds V to virial coefficients
    @param V: values that should be added*/
    void add_to_virial(std::array<std::array<double, 3>, 3> & V)
    {
      m_virial[0][0] += V[0][0];
      m_virial[1][0] += V[1][0];
      m_virial[2][0] += V[2][0];
      m_virial[0][1] += V[0][1];
      m_virial[1][1] += V[1][1];
      m_virial[2][1] += V[2][1];
      m_virial[0][2] += V[0][2];
      m_virial[1][2] += V[1][2];
      m_virial[2][2] += V[2][2];
    };
    /**returns the virial coefficients*/
    coords::virial_t const& virial() const { return m_virial; }

  protected:
  
    /**vector of broken bonds (determined by validate_bonds())
    i.e. bondlength either too short or too long
    each element of the vector is a vector which contains the numbers of the two atoms that form the bond and the bond length*/
    std::vector<std::vector<float>> broken_bonds;
  };
}