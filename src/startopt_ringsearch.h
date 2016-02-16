#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <array>
#include <iterator>

#include "coords.h"
#include "histogram.h"
#include "scon_vect.h"
#include "scon_matrix.h"
#include "startopt.h"

namespace startopt
{

  namespace ringsearch
  {

    static const coords::angle_type RING_DIH[3][4] =
    { // dihedral values for hydrogen bond mediated rings
      {
        coords::angle_type::from_deg(0.0),
        coords::angle_type::from_deg(35.0),
        coords::angle_type::from_deg(0.0),
        coords::angle_type::from_deg(0.0)
      }, // 5 [first two]
      {
        coords::angle_type::from_deg(60.0),
        coords::angle_type::from_deg(-60.0),
        coords::angle_type::from_deg(60.0),
        coords::angle_type::from_deg(0.0)
      }, // 6 [frist three]
      {
        coords::angle_type::from_deg(60.0),
        coords::angle_type::from_deg(-60.0),
        coords::angle_type::from_deg(95.0),
        coords::angle_type::from_deg(-35.0)
      } // 7 [all]
    };

    class ring
    {
    public:
      const std::size_t size;
      std::vector<std::size_t> atoms, dihedrals;
      std::vector<double> dihedral_direction;
      ring(std::size_t _s);
      ring(std::size_t _s, std::vector<std::size_t> _v);
      ring(ring const &r);
      ring& operator= (const ring &r);
    };

    class Search
    {
      coords::Representation_3D const m_init_xyz;
      coords::Ensemble_PES     m_final_ensemble;
      coords::Coordinates &    m_coord;
      std::vector<std::size_t> m_donors;
      std::vector<bool>        m_is_acceptor;
      std::vector<ring>        m_rings;
      scon::matrix<bool, true> m_overlap;
      std::vector<std::size_t> m_ringcontainer;
    public:
      Search(coords::Coordinates & coord_object);
      Search& operator= (Search const &) = delete;
      ~Search();
      // Ensemble Access
      coords::Ensemble_PES & ensemble() { return m_final_ensemble; }
      // Coord access
      coords::Coordinates const & get_coords() const { return m_coord; }
      // rings access
      std::vector<ring> const & rings() const { return m_rings; }
      // Get energy for specific closure combination
      coords::float_type energy(std::vector<bool> const &close_it, coords::float_type const force = 0.1);
    private:
      // Find starting points and connect them via find_route
      std::size_t find_rings();
      // Find route from donor to acceptor
      void find_route(std::size_t const i, std::size_t const size, std::vector<bool> tabulist);
      // Find internal main torsion
      std::size_t find_torsion(double &direction, std::size_t const index_b, std::size_t const index_c) const;
      // Look for main torsion for every ring dihedral via find_torsion
      void find_torsions();
      // Set overlap matrix
      void find_overlap();
      // Apply bias for specific ring
      void bias_ring(std::size_t const index, coords::float_type const force = 0.1);
      // Apply bias for all rings that have a true value in the vector
      std::vector<std::size_t> bias_rings(std::vector<bool> const &close_it, coords::float_type const force = 0.1);

    };

    coords::Ensemble_PES genetic_propagation(coords::Coordinates & coord_object,
      std::size_t const population, std::size_t const generations);

  }


  namespace preoptimizers
  {

    class R_evolution
      : public Preoptimizer
    {
    public:

      // construct preoptimizer
      R_evolution(coords::Coordinates const & init_coords)
        : Preoptimizer(init_coords) {}
      // disallow assignment
      R_evolution& operator= (R_evolution const &) = delete;
      // generate 'm' preoptimized structures for each element of 'init_ensemble'
      void generate(coords::Ensemble_PES const & init_ensemble, std::size_t const m);

    };

  }


}
