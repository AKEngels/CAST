#pragma once

#include <cstring>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <unordered_map>

#include "global.h"
#include "configuration.h"
#include "scon_linkedcell.h"
#include "scon_vect.h"
#include "startopt.h"
#include "representation.h"


namespace startopt
{

  namespace solvadd
  {

    struct site
    {
      // v is site direction
      // p is site origin
      coords::Cartesian_Point v, p;
      // atom is origin atom index
      std::size_t atom;
      // donor determines whether atom === H
      bool tabu, donor;
      site() : v{}, p{}, atom{}, tabu{ false }, donor{ false } {}
    };

    struct water
    {
      coords::Cartesian_Point o, h[2];
      bool check_geometry() const;
    };

    // build sites
    std::vector<site> build_hb_sites(coords::Representation_3D const &xyz, 
      coords::Atoms const &atms, std::vector<bool> const &tabu);

    // remove non-accessible sites
    std::vector<site> accessible_sites(std::vector<site> const &sites, 
      coords::Representation_3D const &xyz);

    // fit water into site and return identified position
    // bool return part tells about success
    std::pair<bool, water> fit_water_into_site(site const &site, 
      std::vector<std::size_t> const &atoms_around, 
      coords::Representation_3D const &xyz,
      coords::Atoms const &atms);

    struct wp_optimum
    {
    private:
      std::size_t m_atoms;
    public:
      wp_optimum (std::size_t atoms) : m_atoms(atoms) { }
      std::size_t ia(std::size_t const num_particles) const { return (num_particles*num_particles - num_particles) / 2; }
      std::size_t npp(void) const { return ia(m_atoms); }
      std::size_t nww (std::size_t const num_water) const { return ia(3*num_water); }
      std::size_t nwp (std::size_t const num_water) const { return m_atoms*num_water*3; }
      double ratio (std::size_t const num_water) const
      {
        std::size_t const n_pp((m_atoms*m_atoms - m_atoms) / 2), 
          n_ww((9 * num_water*num_water - 3 * num_water) / 2), 
          n_wp(m_atoms*num_water * 3);
        return (static_cast<double>(n_wp-n_ww)/static_cast<double>(n_pp+n_ww+n_wp));
      }
      operator std::size_t (void) const
      {
        std::size_t w(0u);
        while (ratio(w+1) > ratio(w)) ++w;
        return w;
      }
    };

  }

  namespace preoptimizers
  {

    class Solvadd 
      : public Preoptimizer
    {

    public:

      using cells_type = scon::linked::Cells < coords::float_type,
        coords::Cartesian_Point, coords::Representation_3D > ;

      Solvadd (coords::Coordinates const & init_coords, double const bound = 10.0)
        : Preoptimizer(init_coords), 
        m_init_center(init_coords.center_of_mass()),
        m_box_max(m_init_center + bound /2.0),
        m_box_min(m_init_center - bound /2.0),
        solvated_coords(init_coords), 
        m_solvated_positions(init_coords.xyz()), m_solvated_atoms(), 
        m_solvated_cells(m_solvated_positions, 3.0, 
          Config::get().energy.periodic, 
          Config::get().energy.periodic ? Config::get().energy.pb_box : coords::Cartesian_Point(),
          Config::get().energy.periodic ? 1.0 : bound + 5.0),
        m_init_cells(init_coords.xyz(), std::max(bound/2.0, 4.0),
          Config::get().energy.periodic, 
          Config::get().energy.periodic ? Config::get().energy.pb_box : coords::Cartesian_Point(), 
          Config::get().energy.periodic ? 1.0 : bound + 2.0),
        m_boundary(bound)
      { }

      void generate (coords::Ensemble_PES const & init_ensemble, std::size_t const multiplier);

    private:

      Solvadd& operator= (Solvadd const &);

      // Solvation sites
      std::vector<solvadd::site> m_sites;
      // Center initial
      coords::Cartesian_Point    m_init_center, m_box_max, m_box_min;
      // Tabu atoms
      std::vector<bool>          tabu_atoms;
       // coordinates object, saving the solvated coords
      coords::Coordinates        solvated_coords;
      // Solvated Positions
      coords::Representation_3D  m_solvated_positions;
      // Solvated Atoms
      coords::Atoms              m_solvated_atoms;
      // linked cells obj
      cells_type                 m_solvated_cells, m_init_cells;
      // links
      std::unordered_map<std::size_t, std::size_t> m_interlinks;
      // boundary                
      double                     m_boundary;

      //! Creates all sites from current solvated_coords
      void build_sites (void);
      //! builds a single extension using internal coords (index, bond, angle, dihedral)
      void build_site(std::size_t const index[3], coords::angle_type, coords::angle_type, bool);
      //! call buildExtension multiple times 
      void build_multisite(std::size_t const index[3], std::size_t const, coords::angle_type const, coords::angle_type const, coords::float_type const, bool const);
      //! Functions specific for atomic numbers
      void build_site_group_1  (std::size_t atom);
      void build_site_group_15 (std::size_t atom);
      void build_site_group_16 (std::size_t atom);
      void build_site_group_17 (std::size_t atom);
      //! populate single site
      bool populate_site (solvadd::site const &);
      // check whether site is applicable
      bool check_sterics (solvadd::water const &w, std::vector<std::size_t> const &atoms_around) const;
      bool check_out_of_boundary (coords::Cartesian_Point const &p) const;
      void push_boundary (void);
      void add_water (solvadd::water const &w);
      void populate_coords (std::size_t const added=0U);
      std::size_t purge_coords (void);

    };

    class GOSol
    {
      // deleted copy constructor
      GOSol& operator= (GOSol const &);
      // Data
      coords::Ensemble_PES m_pes;
      coords::Coordinates m_coords;
      std::size_t const m_num_init_atoms;
      // Functions
      std::size_t f_solvate(std::size_t const max_num_w);
      void f_optimize (std::string const &);
    public:
      GOSol(coords::Coordinates & init_coords, coords::Ensemble_PES const & init_pes)
        : m_pes(init_pes), m_coords(init_coords), m_num_init_atoms(init_coords.size())
      {}
      void run(std::size_t const num_water);
    };

  }

}
