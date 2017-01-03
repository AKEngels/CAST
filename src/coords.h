#ifndef coords_h_guard_
#define coords_h_guard_  

#include <vector>
#include <string>
#include <utility> 
#include <array>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cstddef>
#include <memory>

#include "configuration.h"
#include "energy.h"
#include "scon_vect.h"
#include "scon_utility.h"
#include "scon_matrix.h"
#include "coords_rep.h"
#include "scon_log.h"
#include "coords_atoms.h"

namespace coords
{


  /* ######################################################


   ######  ######## ######## ########  ########  #######
  ##    ##    ##    ##       ##     ## ##       ##     ##
  ##          ##    ##       ##     ## ##       ##     ##
   ######     ##    ######   ########  ######   ##     ##
        ##    ##    ##       ##   ##   ##       ##     ##
  ##    ##    ##    ##       ##    ##  ##       ##     ##
   ######     ##    ######## ##     ## ########  #######


  ###################################################### */


  class Stereo
  {
  public:

    struct pair
    {
      std::size_t m_atom;
      std::size_t bonds[4];
      bool m_dir;
      pair() : m_atom(0U), bonds(), m_dir(false)
      { }
      pair(std::size_t atom, bool direction)
        : m_atom(atom), bonds(), m_dir(direction)
      { }
    };

    Stereo() : m_centers() { }
    Stereo(Atoms const &, coords::Representation_3D const &);

    std::vector<pair> const & centers() const { return m_centers; }
    void update(coords::Representation_3D const &);
    void swap(Stereo &r) { m_centers.swap(r.m_centers); }
    void clear() { m_centers.clear(); }

    friend std::ostream& operator<< (std::ostream &stream, Stereo const & stereo);

  private:

    std::vector< pair > m_centers;

  };

  inline bool operator== (Stereo::pair const &a, Stereo::pair const &b)
  {
    return a.m_atom == b.m_atom && a.m_dir == b.m_dir;
  }
  inline bool operator!= (Stereo::pair const &a, Stereo::pair const &b)
  {
    return !(a == b);
  }

  inline bool operator== (Stereo const &a, Stereo const &b)
  {
    std::size_t const N(a.centers().size());
    if (N != b.centers().size()) return false;
    for (std::size_t i(0); i < N; ++i)
    {
      if (a.centers()[i] != b.centers()[i]) return false;
    }
    return true;
  }
  inline bool operator!= (Stereo const &a, Stereo const &b)
  {
    return !(a == b);
  }

  inline void swap(Stereo &a, Stereo &b)
  {
    a.swap(b);
  }



  /* ######################################################


  ########  ####    ###     ######
  ##     ##  ##    ## ##   ##    ##
  ##     ##  ##   ##   ##  ##
  ########   ##  ##     ##  ######
  ##     ##  ##  #########       ##
  ##     ##  ##  ##     ## ##    ##
  ########  #### ##     ##  ######


  ###################################################### */


  namespace bias
  {
    struct Potentials
    {

      Potentials();

      bool empty() const;
      bool uempty() const { return m_utors.empty() && m_udist.empty(); }
      double energy() const { return b + a + d + s + c; }

      void clear()
      {
        b = a = d = s = c = 0.0;
        scon::clear(m_dihedrals, m_angles, m_distances,
          m_spherical, m_cubic, m_utors, m_udist);
      }

      double e_dist() const { return b; }
      double e_angle() const { return a; }
      double e_dihedral() const { return d; }
      double e_spherical() const { return s; }
      double e_cubic() const { return c; }

      void add(config::biases::dihedral const &new_d) { m_dihedrals.push_back(new_d); }
      void add(config::biases::angle const &new_a) { m_angles.push_back(new_a); }
      void add(config::biases::distance const &new_d) { m_distances.push_back(new_d); }
      void add(config::biases::spherical const &new_d) { m_spherical.push_back(new_d); }
      void add(config::biases::cubic const &new_d) { m_cubic.push_back(new_d); }

      std::vector<config::biases::dihedral> const & dihedrals() const { return m_dihedrals; }
      std::vector<config::biases::angle> const & angles() const { return m_angles; }
      std::vector<config::biases::distance> const & distances() const { return m_distances; }
      std::vector<config::biases::spherical> const & sphericals() const { return m_spherical; }
      std::vector<config::biases::cubic> const & cubic() const { return m_cubic; }

      double apply(Representation_3D const & xyz, Representation_3D & g_xyz,
        Cartesian_Point const & center = Cartesian_Point());
      void umbrellaapply(Representation_3D const & xyz,
        Representation_3D & g_xyz, std::vector<double> &uout);

      void append_config();

      void swap(Potentials &rhs);

    private:

      double b, a, d, s, c;
      std::vector<config::biases::dihedral>  m_dihedrals;
      std::vector<config::biases::angle>     m_angles;
      std::vector<config::biases::distance>  m_distances;
      std::vector<config::biases::spherical> m_spherical;
      std::vector<config::biases::cubic>     m_cubic;
      std::vector<config::coords::umbrellas::umbrella_tor> m_utors;
      std::vector<config::coords::umbrellas::umbrella_dist> m_udist;

      double dih(Representation_3D const & xyz, Gradients_3D & g_xyz);
      double dist(Representation_3D const & xyz, Gradients_3D & g_xyz);
      double ang(Representation_3D const & xyz, Gradients_3D & g_xyz);
      double spherical(Representation_3D const & xyz, Gradients_3D & g_xyz,
        Cartesian_Point const & center = Cartesian_Point());
      double cubic(Representation_3D const & xyz, Gradients_3D & g_xyz,
        Cartesian_Point const & center = Cartesian_Point());
      void umbrelladih(Representation_3D const & xyz, Gradients_3D & g_xyz, std::vector<double> &uout)  const;
      void umbrelladist(Representation_3D const & xyz, Gradients_3D & g_xyz, std::vector<double> &uout)  const;
    };
  }

  struct fep_data
  {
    energy::fepvect feptemp;
    std::vector<energy::fepvect> fepdata;
    std::vector<energy::fepvar> window;
  };

  /* struct pme_data
   {
     energy::pmevar pmetemp;
   };*/

   /* ##################################################################################################



     ######   #######   #######  ########  ########  #### ##    ##    ###    ######## ########  ######
    ##    ## ##     ## ##     ## ##     ## ##     ##  ##  ###   ##   ## ##      ##    ##       ##    ##
    ##       ##     ## ##     ## ##     ## ##     ##  ##  ####  ##  ##   ##     ##    ##       ##
    ##       ##     ## ##     ## ########  ##     ##  ##  ## ## ## ##     ##    ##    ######    ######
    ##       ##     ## ##     ## ##   ##   ##     ##  ##  ##  #### #########    ##    ##             ##
    ##    ## ##     ## ##     ## ##    ##  ##     ##  ##  ##   ### ##     ##    ##    ##       ##    ##
     ######   #######   #######  ##     ## ########  #### ##    ## ##     ##    ##    ########  ######


    ################################################################################################## */

  using Tensor = std::array<std::array<float_type, 3u>, 3u>;

  typedef Tensor virial_t;


  static inline coords::virial_t empty_virial()
  {
    coords::virial_t v;
    for (std::size_t i = 0; i < 3; ++i)
    {
      for (std::size_t j = 0; j < 3; ++j) v[i][j] = 0.0;
    }
    return v;
  }


  class Coordinates
  {

    Atoms                     m_atoms;
    PES_Point                 m_representation;
    Stereo                    m_stereo;
    bias::Potentials          m_potentials;
    virial_t                  m_virial;
    energy::interface_base    *m_interface;
    energy::interface_base    *m_preinterface;
    bool                      energy_valid;

    coords::float_type prelbfgs();
    coords::float_type lbfgs();

    coords::float_type m_e(energy::interface_base * const p)
    {
      if (p) 
      {
        energy_valid = true;
        if (Config::get().energy.periodic) periodic_boxjump();
        m_representation.energy = p->e();
        m_representation.integrity = p->intact();
        apply_bias();
        zero_fixed_g();
        return m_representation.energy;
      }
      return coords::float_type();
    }

    coords::float_type m_g(energy::interface_base * const p)
    {
      if (p)
      {
        energy_valid = true;
        if (Config::get().energy.periodic) periodic_boxjump();
        m_representation.energy = p->g();
        m_representation.integrity = p->intact();
        this->apply_bias();
        zero_fixed_g();
        return m_representation.energy;
      }
      return coords::float_type();
    }

  public:

    fep_data              fep;

    /*pme_data              pme;*/

    bool                  NEB_control, PathOpt_control;

    std::size_t                mult_struc_counter;
    //aditional variables for perpendicular g() and o() and NEB

    bool									orthogonalize, hessian_def, use_fep;
    //void topo(std::size_t const i, std::size_t const j)
    //{ 
    //  scon::sorted::insert_unique(m_topology[j], i); 
    //}    
    //
    typedef PES_Point::size_type size_type;

    //fep_data    fep;
    //bool        NEB_control, PathOpt_control;
    //size_type   mult_struc_counter;
    //// aditional variables for perpendicular g() and o() and NEB
    //Ensemble_3d R_o, taux, NEB_taux, image_inix;
    //bool				orthogonalize, hessian_def, use_fep;

    Coordinates();
    Coordinates(Coordinates const &r);
    Coordinates(Coordinates &&r);
    Coordinates& operator= (Coordinates&& rhs);
    Coordinates& operator= (Coordinates const & rhs);
    ~Coordinates();

    void set_fix(size_t const atom, bool const fix_it = true);

    void clear()
    {
      *this = Coordinates{};
    }

    void apply_bias()
    {
      if (!m_potentials.empty())
      {
        m_representation.energy += m_potentials.apply(
          m_representation.structure.cartesian,
          m_representation.gradient.cartesian,
          Cartesian_Point());
      }
    }

    void ubias(std::vector<double> &uout)
    {
      if (!m_potentials.uempty())
        m_potentials.umbrellaapply(m_representation.structure.cartesian,
          m_representation.gradient.cartesian,
          uout);
    }
    //prelim pme function
    //void pme_stuff(int);
    //void pme_dftmodulus(std::vector<double> &);
    //void roughgrid();
    //void getfepinfo();

    coords::float_type pe()
    { // energy+gradients
      return m_e(m_preinterface);
    }
    coords::float_type e()
    { // energy
      return m_e(m_interface);
    }

    coords::float_type pg()
    { // energy+gradients
      return m_g(m_preinterface);
    }
    coords::float_type g()
    { // energy+gradients
      return m_g(m_interface);
    }

    coords::float_type po()
    {
      if (m_preinterface)
      {
        if (Config::get().general.verbosity >= 4)
          std::cout << "Preotimization will be performed.\n";
        energy_valid = true;
        if (m_preinterface->has_optimizer()
          && m_potentials.empty()
          && !Config::get().energy.periodic)
        {
          m_representation.energy = m_preinterface->o();
        }
        else m_representation.energy = prelbfgs();
        m_representation.integrity = m_preinterface->intact();
        m_stereo.update(xyz());
        zero_fixed_g();
        return m_representation.energy;
      }
      return 0.;
    }

    coords::float_type o()
    {
      if (preoptimize()) po();
      energy_valid = true;
      if (m_interface->has_optimizer()
        && m_potentials.empty()
        && !Config::get().energy.periodic)
      {
        m_representation.energy = m_interface->o();
      }
      else
      {
        m_representation.energy = lbfgs();
      }
      m_representation.integrity = m_interface->intact();
      m_stereo.update(xyz());
      zero_fixed_g();
      return m_representation.energy;
    }

    bool preoptimize() const { return m_preinterface ? true : false; }

    coords::Gradients_Main dimermethod_dihedral(std::vector<coords::Gradients_Main> const &tabu_direction);
    coords::Gradients_Main dimermethod_dihedral()
    {
      return dimermethod_dihedral(std::vector<coords::Gradients_Main>());
    }

    energy::interface_base const * energyinterface() const { return m_interface; }
    energy::interface_base const * preinterface() const { return m_preinterface; }

    // helpers
    bool               integrity() const { return pes().integrity; }
    size_type          size() const { return m_atoms.size(); }
    void               swap(Coordinates &rhs); // object swap

    void set_virial(virial_t const & V) { m_virial = V; }

    void add_to_virial(std::array<std::array<double, 3>, 3>& V)
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

    virial_t const &   virial() const { return m_virial; }
    Cartesian_Point    center_of_mass() const; // center of mass
    Cartesian_Point    center_of_mass_mol(size_type const index) const;
    Cartesian_Point    center_of_geometry() const; // center of geometry
    coords::float_type weight() const; // total mass of system

	/**vector of broken bonds (determined by validate_bonds())
	i.e. bondlength either too short or too long
	each element of the vector is a vector which contains the numbers of the two atoms that form the bond and the bond length*/
	std::vector<std::vector<float>> broken_bonds;

    void init_swap_in(Atoms &a, PES_Point &p, bool const energy_update = true); // initial data swap in
    void init_in(Atoms a, PES_Point p, bool const energy_update = true); // new data input
    void energy_update(bool const skip_topology = false) { m_interface->update(skip_topology); }
    void fix_all(bool const fix_it = true) { m_atoms.fix_all(fix_it); }
    void fix(size_type const index, bool const fix = true) { m_atoms.fix(index, fix); }

    void e_head_tostream_short(std::ostream &strm, energy::interface_base const * const ep = nullptr) const;
    void e_tostream_short(std::ostream &strm, energy::interface_base const * const ep = nullptr) const;

    // PES_Point
    PES_Point const & pes() const { return m_representation; }
    PES_Point & pes() { return m_representation; }

    // xyz representation, const
    Representation_3D const & xyz() const
    {
      return m_representation.structure.cartesian;
    }
    // internal representation, const
    Representation_Internal const & intern() const
    {
      return m_representation.structure.intern;
    }
    // Main internal representation, const
    Representation_Main const & main() const
    {
      return m_representation.structure.main;
    }

    // single position, const
    cartesian_type const & xyz(size_type index) const
    {
      return m_representation.structure.cartesian[index];
    }
    // single internal, const
    internal_type const & intern(size_type index) const
    {
      return m_representation.structure.intern[index];
    }
    // single main dihedral, const
    main_type const & main(size_type index) const
    {
      return m_representation.structure.main[index];
    }

    // gradient representations
    Gradients_3D const & g_xyz() const
    {
      return m_representation.gradient.cartesian;
    }
    Gradients_Internal const & g_intern() const
    {
      return m_representation.gradient.intern;
    }
    Gradients_Main const & g_main() const
    {
      return m_representation.gradient.main;
    }

    cartesian_gradient_type const & g_xyz(size_type const index) const
    {
      return m_representation.gradient.cartesian[index];
    }
    internal_gradient_type const & g_intern(size_type const index) const
    {
      return m_representation.gradient.cartesian[index];
    }
    main_gradient_type const & g_main(size_type const index) const
    {
      return m_representation.gradient.main[index];
    }

    size_2d const & subsystems() const { return m_atoms.subsystems(); }
    size_1d const & subsystems(size_type index) const { return m_atoms.subsystems(index); }
    size_2d const & molecules() const { return m_atoms.molecules(); }
    size_1d const & molecules(size_type index) const { return m_atoms.molecules(index); }

    std::vector< Stereo::pair > const & stereos() const { return m_stereo.centers(); }

	/**looks if all bonds are okay (reasonable bond length)
	and saves the ones that aren't into the vector broken_bonds*/
    bool validate_bonds();

    // Setters
    void periodic_boxjump();

    // move atom
    void move_atom_by(size_type const index, cartesian_type const &p, bool const force_move = false)
    {
      if (!atoms(index).fixed() || force_move)
      {
        m_representation.structure.cartesian[index] += p;
        energy_valid = false;
        m_stereo.update(xyz());
      }
    }
    //scale atom
    void scale_atom_by(size_type const index, double &p, bool const force_move = false)
    {
      if (!atoms(index).fixed() || force_move)
      {
        m_representation.structure.cartesian[index] *= p;
        energy_valid = false;
      }
      m_stereo.update(xyz());
    }
    // set new atom coordinates if anisotropic pressure control is enabled
    void set_atom_aniso(size_type const index, double tvec[3][3], bool const force_move = false) {
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
    // add gradients form spherical boundaries
    void add_sp_gradients(size_type const index, Cartesian_Point const &g)
    {
      m_representation.gradient.cartesian[index] += g;
    }
    void move_atom_to(size_type const index, Cartesian_Point const &p, bool const force_move = false)
    {
      if (!atoms(index).fixed() || force_move)
      {
        m_representation.structure.cartesian[index] = p;
        energy_valid = false;
        m_stereo.update(xyz());
      }
    }

    // move all atoms along specified vector
    void move_all_by(Cartesian_Point const &p, bool const force_move = false)
    {
      size_type const N(size());
      for (size_type i(0U); i < N; ++i)
      {
        if (force_move || !atoms(i).fixed()) m_representation.structure.cartesian[i] += p;
      }
      energy_valid = false;
    }

    // Set Internals
    void set_internal(Representation_Internal const & ri, bool const force_move = false)
    {
      size_type const N = m_representation.structure.intern.size();
      if (N == ri.size() && force_move) m_representation.structure.intern = ri;
      else
      {
        for (size_type i(0u); i < N; ++i)
        {
          if (!atoms(i).ifix()) m_representation.structure.intern[i] = ri[i];
        }
      }

    }

    // rotate dihedrals
    void set_dih(size_type const int_index, coords::angle_type const target_angle,
      bool const move_dependants_along = true, bool const move_fixed_dih = false);
    void rotate_dih(size_type const int_index, coords::angle_type const rot_angle,
      bool const move_dependants_along = true, bool const move_fixed_dih = false);
    void rotate_main(size_type const main_index, coords::angle_type const rot_angle,
      bool const move_dependants_along = true, bool const move_fixed_dih = false);
    void set_all_main(Representation_Main const & new_values, bool const aplly_to_xyz = true,
      bool const move_dependants_along = true, bool const move_fixed_dih = false);

    // update gradients
    void update_g_xyz(size_type const index, Cartesian_Point const &p)
    {
      if (!atoms(index).fixed()) m_representation.gradient.cartesian[index] = p;
    }

    // sum gradients
    void sum_g_xyz(Gradients_3D const & rep)
    {
      m_representation.gradient.cartesian += rep;
      zero_fixed_g();
    }

    // zero fixed
    void zero_fixed_g()
    {
      size_type const N(size());
      for (size_type i(0U); i < N; ++i)
      {
        if (atoms(i).fixed()) m_representation.gradient.cartesian[i] = cartesian_gradient_type(0);
      }
    }

    // put in pes point
    void set_pes(PES_Point const & pes_point, bool const overwrite_fixed = false);
    void set_pes(PES_Point && pes_point, bool const overwrite_fixed = false);

    void set_xyz(Representation_3D const & new_xyz, bool const overwrite_fixed = false)
    {
      size_type const N(size());
      if (new_xyz.size() != N) throw std::logic_error("Wrong sized coordinates in set_xyz.");
      if (!overwrite_fixed)
      {
        for (size_type i(0U); i < N; ++i)
        {
          if (!atoms(i).fixed()) m_representation.structure.cartesian[i] = new_xyz[i];
        }
      }
      else
        m_representation.structure.cartesian = Representation_3D(new_xyz.begin(), new_xyz.end());
      m_stereo.update(xyz());
    }

#if defined(SCON_CC11_RVALUE_REF) && defined(SCON_CC11_MOVE)
    void set_xyz(Representation_3D && new_xyz, bool const overwrite_fixed = false)
    {
      size_type const N(size());
      if (new_xyz.size() != N) throw std::logic_error("Wrong sized coordinates in set_xyz.");
      m_representation.structure.cartesian.swap(new_xyz);
      if (!overwrite_fixed)
      {
        for (size_type i(0U); i < N; ++i)
        {
          if (atoms(i).fixed()) m_representation.structure.cartesian[i] = std::move(new_xyz[i]);
        }
      }
      else
        m_representation.structure.cartesian = Representation_3D(new_xyz.begin(), new_xyz.end());
      m_stereo.update(xyz());
    }
#endif

    // swap gradients in/out
    void swap_g_xyz(Gradients_3D & rhs, bool const overwrite_fixed = false)
    {
      size_type const N(m_representation.gradient.cartesian.size());
      if (rhs.size() != N) throw std::logic_error("Wrong sized gradients in swap_g_xyz.");
      rhs.swap(m_representation.gradient.cartesian); // swap xyz gradients in
      if (!overwrite_fixed)
      {
        for (size_type i(0U); i < N; ++i)
        {
          if (atoms(i).fixed()) m_representation.gradient.cartesian[i] = rhs[i];
        }
      }
      m_atoms.c_to_i(m_representation); // update internals
    }

    void get_g_xyz(Gradients_3D & out_g_xyz) const
    {
      out_g_xyz = m_representation.gradient.cartesian;
    }
    void clear_g_xyz()
    {
      m_representation.gradient.cartesian.assign(size(), Cartesian_Point());
    }

    Atom const & atoms(size_type const index) const
    {
      return m_atoms.atom(index);
    }
    Atoms const & atoms() const
    {
      return m_atoms;
    }

    bias::Potentials & potentials()
    {
      return m_potentials;
    }
    bias::Potentials const & potentials() const
    {
      return m_potentials;
    }

    // Subsystem interactions
    sub_ia_matrix_t & interactions()
    {
      return m_representation.ia_matrix;
    }
    sub_ia_matrix_t const & interactions() const
    {
      return m_representation.ia_matrix;
    }
    sub_ia & interactions(size_type const index)
    {
      return m_representation.ia_matrix(index);
    }
    sub_ia const & interactions(size_type const index) const
    {
      return m_representation.ia_matrix(index);
    }
    sub_ia & interactions(size_type const x, size_type const y)
    {
      return m_representation.ia_matrix(x, y);
    }
    sub_ia const & interactions(size_type const x, size_type const y) const
    {
      return m_representation.ia_matrix(x, y);
    }

    void to_internal() { m_atoms.c_to_i(m_representation); }
    void to_internal_light() { m_atoms.c_to_i_light(m_representation); }

    void to_xyz()
    {
      m_atoms.i_to_c(m_representation);
      if (Config::get().energy.periodic) { periodic_boxjump(); }
      m_stereo.update(xyz());
    }

    bool check_superposition_xyz(Representation_3D const &a,
      Representation_3D const &b, double const x = 0.35) const;

    bool equal_structure(coords::PES_Point const &a, coords::PES_Point const &b, 
      coords::main_type const md = coords::main_type::from_deg(8.0), 
      coords::internal_type const &id = coords::internal_type{
        0.2, 
        coords::angle_type::from_deg(1.0), 
        coords::angle_type::from_deg(8.0)}, 
      coords::Cartesian_Point const &cd = coords::Cartesian_Point{0.1, 0.1, 0.1}) const;

  };

  std::ostream& operator<< (std::ostream &stream, Coordinates const & coord);

  struct internal_float_callback
  {
    coords::Coordinates * cp;
    internal_float_callback(coords::Coordinates & coordinates_object)
      : cp(&coordinates_object)
    { }
    float operator() (scon::vector<float> const & v,
      scon::vector<float> & g, std::size_t const S, bool & go_on)
    {
      std::size_t i = 0;
      coords::Representation_Internal rin(v.size() / 3);
      for (auto & e : rin)
      {
        e.radius() = v[i++];
        e.inclination() = coords::angle_type(v[i++]);
        e.azimuth() = coords::angle_type(v[i++]);
      }
      cp->set_internal(rin);
      cp->to_xyz();
      float E = float(cp->g());
      cp->to_internal();
      go_on = cp->integrity();
      g.resize(v.size());
      i = 0;
      for (auto const & e : cp->g_intern())
      {
        g[i++] = static_cast<float>(e.x());
        g[i++] = static_cast<float>(e.y());
        g[i++] = static_cast<float>(e.z());
      }
      if (Config::get().general.verbosity >= 4)
      {
        std::cout << "Optimization: Energy of step " << S;
        std::cout << " is " << E << " integrity " << go_on << '\n';
      }
      return E;
    }
    scon::vector<float> from_rep(coords::Representation_Internal const & v)
    {
      scon::vector<float> r(v.size() * 3);
      std::size_t i = 0;
      for (auto e : v)
      {
        r[i++] = static_cast<float>(e.radius());
        r[i++] = static_cast<float>(e.inclination().radians());
        r[i++] = static_cast<float>(e.azimuth().radians());
      }
      return r;
    }
    scon::vector<float> from_grad(coords::Gradients_Internal const & v)
    {
      scon::vector<float> r(v.size() * 3);
      std::size_t i = 0;
      for (auto e : v)
      {
        r[i++] = static_cast<float>(e.x());
        r[i++] = static_cast<float>(e.y());
        r[i++] = static_cast<float>(e.z());
      }
      return r;
    }
    coords::Representation_Internal to_rep(scon::vector<float> const & v)
    {
      coords::Representation_Internal r(v.size() / 3);
      std::size_t i = 0;
      for (auto & e : r)
      {
        auto radius = v[i++];
        auto inclination = v[i++];
        auto azimuth = v[i++];
        e = coords::internal_type(radius, coords::angle_type(inclination), coords::angle_type(azimuth));
      }
      return r;
    }
    coords::Gradients_Internal to_grad(scon::vector<float> const & v)
    {
      coords::Gradients_Internal r(v.size() / 3);
      std::size_t i = 0;
      for (auto & e : r)
      {
        auto radius = v[i++];
        auto inclination = v[i++];
        auto azimuth = v[i++];
        e = coords::internal_gradient_type(radius, inclination, azimuth);
      }
      return r;
    }
  };

  struct cartesian_float_callback
  {
    coords::Coordinates * cp;
    cartesian_float_callback(coords::Coordinates & coordinates_object)
      : cp(&coordinates_object)
    { }
    float operator() (scon::vector<float> const & v,
      scon::vector<float> & g, std::size_t const S, bool & go_on)
    {
      coords::Representation_3D x;
      x.reserve(v.size() / 3);
      std::size_t n = 0;
      for (auto & e : x)
      {
        e.x() = v[n++];
        e.y() = v[n++];
        e.z() = v[n++];
      }
      cp->set_xyz(std::move(x), true);
      float E = float(cp->g());
      go_on = cp->integrity();
      g.resize(v.size());
      n = 0;
      for (auto const & e : cp->g_xyz())
      {
        g[n++] = static_cast<float>(e.x());
        g[n++] = static_cast<float>(e.y());
        g[n++] = static_cast<float>(e.z());
      }
      if (Config::get().general.verbosity >= 4)
      {
        std::cout << "Optimization: Energy of step " << S;
        std::cout << " is " << E << " integrity " << go_on << '\n';
        std::cout << "|g| = " << scon::geometric_length(g) << "\n";
      }
      return E;
    }
  };

  struct cartesian_float_stream_log_drain
  {
    coords::Coordinates * cp;
    std::ostream * strm;
    cartesian_float_stream_log_drain(coords::Coordinates & coords,
      std::ostream & stream) : cp(&coords), strm(&stream) { }
    void operator() (optimization::Point<scon::vector<float>,
      scon::vector<float>, float> & v)
    {
      if (strm && cp)
      {
        auto tmp = cp->xyz();
        coords::Representation_3D x;
        x.reserve(v.x.size() / 3);
        std::size_t n = 0;
        for (auto & e : x)
        {
          e.x() = v.x[n++];
          e.y() = v.x[n++];
          e.z() = v.x[n++];
        }
        cp->set_xyz(x);
        *strm << *cp;
        cp->set_xyz(std::move(tmp), true);
      }
    }
  };

  struct internal_float_stream_log_drain
  {
    coords::Coordinates * cp;
    std::ostream * strm;
    internal_float_stream_log_drain(coords::Coordinates & coords,
      std::ostream & stream) : cp(&coords), strm(&stream)
    { }
    void operator() (optimization::Point<scon::vector<float>,
      scon::vector<float>, float> & v)
    {
      if (strm && cp)
      {
        auto p = cp->pes();
        std::size_t i = 0;
        coords::Representation_Internal rin(v.x.size() / 3);
        for (auto & e : rin)
        {
          e.radius() = v.x[i++];
          e.inclination() = coords::angle_type(v.x[i++]);
          e.azimuth() = coords::angle_type(v.x[i++]);
        }
        cp->set_internal(rin);
        cp->to_xyz();
        *strm << *cp;
        cp->set_pes(std::move(p), true);
      }
    }
  };

  using internal_log_drain = internal_float_stream_log_drain;
  using internal_log = scon::vector_offset_buffered_callable<
    optimization::Point<scon::vector<float>,
    scon::vector<float>, float>, internal_log_drain >;

  struct Coords_3d_float_callback
  {
    coords::Coordinates * cp;
    //std::unique_ptr<std::ofstream> ls;
    Coords_3d_float_callback(coords::Coordinates & coordpointer) :
      cp(&coordpointer)/*, ls(new std::ofstream("st.arc"))*/ { }
    float operator() (scon::vector<scon::c3<float>> const & v,
      scon::vector<scon::c3<float>> & g, std::size_t const S, bool & go_on);
    scon::vector<scon::c3<float>> from(coords::Gradients_3D const & g);
    coords::Representation_3D to(scon::vector<scon::c3<float>> const & v);
  };

  struct Coords_3d_float_pre_callback
  {
    coords::Coordinates * cp;
    Coords_3d_float_pre_callback(coords::Coordinates & coordpointer) : cp(&coordpointer) { }
    float operator() (scon::vector<scon::c3<float>> const & v,
      scon::vector<scon::c3<float>> & g, std::size_t const S, bool & go_on);
  };

  struct Internal_Callback
  {
    coords::Coordinates * cp;
    Internal_Callback(coords::Coordinates & coordpointer)
      : cp(&coordpointer)
    { }
    coords::float_type operator() (scon::vector<coords::float_type> const & v,
      scon::vector<coords::float_type> & g, std::size_t const S, bool & go_on);
    coords::Gradients_Internal from(coords::Representation_Internal const & p);
    coords::Representation_Internal to(coords::Gradients_Internal const & p);
  };

  struct Internal_Log
  {

    void operator() (optimization::Point<coords::Gradients_Internal,
      coords::Gradients_Internal, coords::float_type> &);
  };

  struct Main_Callback
  {
    coords::Coordinates * cp;
    Main_Callback(coords::Coordinates & coordpointer)
      : cp(&coordpointer)
    { }
    coords::float_type operator() (coords::Gradients_Main const & v,
      coords::Gradients_Main & g, std::size_t const S, bool & go_on);
  };

  class cartesian_logfile_drain
  {
    coords::Coordinates * cp;
    std::unique_ptr<std::ofstream> strm;
    bool opt;
  public:
    cartesian_logfile_drain() : cp(), strm(), opt() {}
    cartesian_logfile_drain(coords::Coordinates &c, char const * const filename, bool optimize = false) :
      cp(&c), strm(new std::ofstream(filename, std::ios::out)), opt(optimize)
    {}
    void operator() (coords::Representation_3D && xyz);
  };

  using offset_buffered_cartesian_logfile =
    scon::vector_offset_buffered_callable<coords::Representation_3D, cartesian_logfile_drain>;

  offset_buffered_cartesian_logfile make_buffered_cartesian_log(Coordinates &c,
    std::string file_suffix, std::size_t buffer_size,
    std::size_t log_offset, bool optimize = false);

  inline void swap(Coordinates &a, Coordinates &b)
  {
    a.swap(b);
  }
}

#endif // coords_h_guard_
