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
        b = a = d = s = c = thr = 0.0;
        scon::clear(m_dihedrals, m_angles, m_distances,
          m_spherical, m_cubic, m_utors, m_udist, m_thresh);
      }

      double e_dist() const { return b; }
      double e_angle() const { return a; }
      double e_dihedral() const { return d; }
      double e_spherical() const { return s; }
      double e_cubic() const { return c; }
      double e_thresh() const {return thr;}

      void add(config::biases::dihedral const &new_d) { m_dihedrals.push_back(new_d); }
      void add(config::biases::angle const &new_a) { m_angles.push_back(new_a); }
      void add(config::biases::distance const &new_d) { m_distances.push_back(new_d); }
      void add(config::biases::spherical const &new_d) { m_spherical.push_back(new_d); }
      void add(config::biases::cubic const &new_d) { m_cubic.push_back(new_d); }
      void add(config::biases::thresholdstr const &new_thr) {m_thresh.push_back(new_thr); }

      std::vector<config::biases::dihedral> const & dihedrals() const { return m_dihedrals; }
      std::vector<config::biases::angle> const & angles() const { return m_angles; }
      std::vector<config::biases::distance> const & distances() const { return m_distances; }
      std::vector<config::biases::spherical> const & sphericals() const { return m_spherical; }
      std::vector<config::biases::cubic> const & cubic() const { return m_cubic; }
      std::vector<config::biases::thresholdstr> const & thresholds() const { return m_thresh; }

      double apply(Representation_3D const & xyz, Representation_3D & g_xyz,
        Cartesian_Point maxPos, Cartesian_Point minPos, Cartesian_Point const & center = Cartesian_Point());
      void umbrellaapply(Representation_3D const & xyz,
        Representation_3D & g_xyz, std::vector<double> &uout);

      void append_config();

      void swap(Potentials &rhs);

    private:

      double b, a, d, s, c, thr, thrB;
      std::vector<config::biases::dihedral>  m_dihedrals;
      std::vector<config::biases::angle>     m_angles;
      std::vector<config::biases::distance>  m_distances;
      std::vector<config::biases::spherical> m_spherical;
      std::vector<config::biases::cubic>     m_cubic;
      std::vector<config::biases::thresholdstr>  m_thresh;
      std::vector<config::biases::thresholdstr>  m_threshBottom;
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
      double thresh(Representation_3D const & xyz, Gradients_3D & g_xyz, Cartesian_Point maxPos);
      double thresh_bottom(Representation_3D const & xyz, Gradients_3D & g_xyz, Cartesian_Point minPos);
    };
  }

  struct fep_data
  {
    energy::fepvect feptemp;
    /**object in which data for one window is saved
    which is relevant for FEP calculation
    (every element of vector is one conformation)*/
    std::vector<energy::fepvect> fepdata;
    /**fep parameters for all windows (every element of vector corresponds to one window)*/
    std::vector<energy::fepvar> window;
  };


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
        if (Config::get().periodics.periodic)
          periodic_boxjump();
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
        if (Config::get().periodics.periodic)
          periodic_boxjump();
        m_representation.energy = p->g();
        m_representation.integrity = p->intact();
        this->apply_bias();
        zero_fixed_g();
        return m_representation.energy;
      }
      return coords::float_type();
    }

  public:

    energy::interface_base   *catch_interface = m_interface;

    void get_catch_interface()
    {
      energy::interface_base   *catch_interface = m_interface;
    }

    fep_data              fep;

    bool                  NEB_control, PathOpt_control;

    std::size_t           mult_struc_counter;
    //aditional variables for perpendicular g() and o() and NEB

    // Orthogonalize is used in Dimer method
    bool									orthogonalize, use_fep;

    typedef PES_Point::size_type size_type;

    Coordinates();
    Coordinates(Coordinates const &r);
    Coordinates(Coordinates &&r);
    Coordinates& operator= (Coordinates&& rhs);
    Coordinates& operator= (Coordinates const & rhs);
    ~Coordinates();

    /**fix an atom, i.e. this atom can't be moved
    @param atom: atom index*/
    void set_fix(size_t const atom, bool const fix_it = true);

    /**delete everything in the Coordinates object -> empty object*/
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
          max_valuePosfix(), min_valuePosfix(),
          Cartesian_Point());
      }
    }

    //umbrella
    void ubias(std::vector<double> &uout)
    {
      if (!m_potentials.uempty())
        m_potentials.umbrellaapply(m_representation.structure.cartesian,
          m_representation.gradient.cartesian,
          uout);
    }

    /**calculates energy with preinterface*/
    coords::float_type pe()
    { // energy+gradients
      return m_e(m_preinterface);
    }
    /**calculates energy*/
    coords::float_type e()
    { // energy
      return m_e(m_interface);
    }
    /**calculates energy+gradients with preinterface*/
    coords::float_type pg()
    { // energy+gradients
      return m_g(m_preinterface);
    }
    /**calculates energy+gradients*/
    coords::float_type g()
    { // energy+gradients
      return m_g(m_interface);
    }
    /**performs an optimisation by steepest gradient method with preinterface*/
    coords::float_type po()
    {
      if (m_preinterface)
      {
        if (Config::get().general.verbosity >= 4)
          std::cout << "Preoptimization will be performed.\n";
        energy_valid = true;
        if (m_preinterface->has_optimizer()
          && m_potentials.empty()
          && !Config::get().periodics.periodic)
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

    /**performs an optimisation by steepest gradient method*/
    coords::float_type o()
    {
      if (preoptimize()) po();
      energy_valid = true;
      if (m_interface->has_optimizer()
        && m_potentials.empty() //bias
        && !Config::get().periodics.periodic)
      {
        m_representation.energy = m_interface->o();
      }
      else
      {
        m_representation.energy = lbfgs();
      }
      m_representation.integrity = m_interface->intact();
      m_stereo.update(xyz());
      zero_fixed_g(); //nullt gradienten alelr fixed atrome
      return m_representation.energy;
    }

    /**calculate hessian matrix*/
    coords::float_type h()
    {
     return m_representation.energy = m_interface->h();
    }

    bool preoptimize() const { return m_preinterface ? true : false; }

    coords::Gradients_Main dimermethod_dihedral(std::vector<coords::Gradients_Main> const &tabu_direction);

    coords::Gradients_Main dimermethod_dihedral()
    {
      return dimermethod_dihedral(std::vector<coords::Gradients_Main>());
    }
    /**returns the energy interface*/
    energy::interface_base const * energyinterface() const { return m_interface; }
    /**returns the energy preinterface*/
    energy::interface_base const * preinterface() const { return m_preinterface; }

    /**determines if structure is intact,
    i.e. bond lengths do not differ from ideal bond length by more than 5 angstrom,
    angles do not differ from ideal angle by more than 20 degrees
    and similar criteria for dihedrals and ureys
    details see in energy_int_..._pot.cc*/
    bool               integrity() const { return pes().integrity; }
    /**returns the size of the coordinates object, i.e. the total number of atoms*/
    size_type          size() const { return m_atoms.size(); }
    /**object swap*/
    void               swap(Coordinates &rhs);

    /**saves virial coefficients into coordinates object
    @param V: virials coefficients*/
    void set_virial(virial_t const & V) { m_virial = V; }

    /**adds V to virial coefficients
    @param V: values that should be added*/
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
    /**returns the virial coefficients*/
    virial_t const &   virial() const { return m_virial; }
    /**calculates the center of mass*/
    Cartesian_Point    center_of_mass() const;
    /**calculates the center of mass for one molecule
    @param index: index of the molecule*/
    Cartesian_Point    center_of_mass_mol(size_type const index) const;
    /**calculates the geometrical center*/
    Cartesian_Point    center_of_geometry() const;
    /**calculates the total mass of the system*/
    coords::float_type weight() const;

    /**vector of broken bonds (determined by validate_bonds())
    i.e. bondlength either too short or too long
    each element of the vector is a vector which contains the numbers of the two atoms that form the bond and the bond length*/
    std::vector<std::vector<float>> broken_bonds;

    /**fills the coordinates object with data
    @param a: atoms object
    @param p: PES_point*/
    void init_swap_in(Atoms &a, PES_Point &p, bool const energy_update = true);
    /**does the same as init_swap_in
    @param a: atoms object
    @param p: PES_point*/
    void init_in(Atoms a, PES_Point p, bool const energy_update = true);
    /**updates the topology*/
    void energy_update(bool const skip_topology = false) { m_interface->update(skip_topology); }
    /**fixes all atoms, i.e. nothing can move anymore*/
    void fix_all(bool const fix_it = true) { m_atoms.fix_all(fix_it); }
    /**fixes a given atom
    @param index: index of atom that is to be fixed*/
    void fix(size_type const index, bool const fix = true) { m_atoms.fix(index, fix); }

    void e_head_tostream_short(std::ostream &strm, energy::interface_base const * const ep = nullptr) const;
    void e_tostream_short(std::ostream &strm, energy::interface_base const * const ep = nullptr) const;
    /**writes hessian matrix
    @param strm: can be std::cout or ofstream file*/
    void h_tostream(std::ostream &strm, energy::interface_base const * const ep = nullptr) const;

    /**returns the PES point*/
    PES_Point const & pes() const { return m_representation; }
    /**also returns the PES point*/
    PES_Point & pes() { return m_representation; }

    /** returns the xyz representation*/
    Representation_3D const & xyz() const
    {
      return m_representation.structure.cartesian;
    }
    /** returns the internal representation*/
    Representation_Internal const & intern() const
    {
      return m_representation.structure.intern;
    }
    /** returns the representation of main dihedrals*/
    Representation_Main const & main() const
    {
      return m_representation.structure.main;
    }

    /**returns the xyz coordinates of a given atom
  @param index: atom index*/
    cartesian_type const & xyz(size_type index) const
    {
      return m_representation.structure.cartesian[index];
    }
    /**returns the internal coordinates of a given atom
    @param index: atom index*/
    internal_type const & intern(size_type index) const
    {
      return m_representation.structure.intern[index];
    }
    /**returns the main dihedral coordinates of a given atom
    @param index: atom index*/
    main_type const & main(size_type index) const
    {
      return m_representation.structure.main[index];
    }

    /**returns the gradients in cartesian space*/
    Gradients_3D const & g_xyz() const
    {
      return m_representation.gradient.cartesian;
    }
    /**returns the gradients in internal space*/
    Gradients_Internal const & g_intern() const
    {
      return m_representation.gradient.intern;
    }
    /**returns the gradients in main dihedral space*/
    Gradients_Main const & g_main() const
    {
      return m_representation.gradient.main;
    }

    /**returns the gradients in cartesian space for a given atom
    @param index: atom index*/
    cartesian_gradient_type const & g_xyz(size_type const index) const
    {
      return m_representation.gradient.cartesian[index];
    }
    /**returns the gradients in internal space for a given atom
    @param index: atom index*/
    internal_gradient_type const & g_intern(size_type const index) const
    {
      return m_representation.gradient.cartesian[index];
    }
    /**returns the gradients in main dihedral space for a given atom
    @param index: atom index*/
    main_gradient_type const & g_main(size_type const index) const
    {
      return m_representation.gradient.main[index];
    }

    /**sets hessian matrix
    @param hess: vector of vectors of doubles (e.g. matrix of doubles) that contains values for hessian matrix*/
    void set_hessian(std::vector<std::vector<double>> hess)
    {
      m_representation.hessian = hess;
    }

    /**returns the hessian matrix*/
    std::vector<std::vector<double>> get_hessian()
    {
      return m_representation.hessian;
    }

    /**returns all subsystems*/
    size_2d const & subsystems() const { return m_atoms.subsystems(); }
    /**returns a given subsystem
    @param index: index of subsystem*/
    size_1d const & subsystems(size_type index) const { return m_atoms.subsystems(index); }
    /**returns all molecules*/
    size_2d const & molecules() const { return m_atoms.molecules(); }
    /**returns a given molecule
    @param index: index of molecule*/
    size_1d const & molecule(size_type index) const { return m_atoms.molecule(index); }

    /**returns stereo centers?*/
    std::vector< Stereo::pair > const & stereos() const { return m_stereo.centers(); }

    /**looks if all bonds are okay (reasonable bond length)
    and saves the ones that aren't into the vector broken_bonds*/
    bool validate_bonds();

    /**if periodic boundaries are activated:
  moves atoms that are outside of the box into the box*/
    void periodic_boxjump();

    /**move atom
  @param index: index of atom that is to be moved
  @param p: "space vector" by which it should be moved
  @param force_move: if set to true also move fixed atoms*/
    void move_atom_by(size_type const index, cartesian_type const &p, bool const force_move = false)
    {
      if (!atoms(index).fixed() || force_move)
      {
        m_representation.structure.cartesian[index] += p;
        energy_valid = false;
        m_stereo.update(xyz());
      }
    }
    /**scale the coordinates of an atom (used for pressure control)
    @param index: index of atom that is to be moved
    @param p: factor by which the coordinates should be scaled
    @param force_move: if set to true also move fixed atoms*/
    void scale_atom_by(size_type const index, double &p, bool const force_move = false)
    {
      if (!atoms(index).fixed() || force_move)
      {
        m_representation.structure.cartesian[index] *= p;
        energy_valid = false;
      }
      m_stereo.update(xyz());
    }
    /** set new atom coordinates if anisotropic pressure control is enabled*/
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
    /** add gradients to an atom (spherical boundaries)
  @param index: atom index
  @param g: gradients that should be added to the gradients of index*/
    void add_sp_gradients(size_type const index, Cartesian_Point const &g)
    {
      m_representation.gradient.cartesian[index] += g;
    }
    /**move given atom to given coordinates
    @param index: index of atom that is to be moved
    @param p: coordinates where the atom should be moved to
    @param force_move: if set to true also move fixed atoms*/
    void move_atom_to(size_type const index, Cartesian_Point const &p, bool const force_move = false)
    {
      if (!atoms(index).fixed() || force_move)
      {
        m_representation.structure.cartesian[index] = p;
        energy_valid = false;
        m_stereo.update(xyz());
      }
    }

    /** move all atoms along specified vector
  @param p: vector along which the system should be moved
  @param force_move: if set to true also move fixed atoms*/
    void move_all_by(Cartesian_Point const &p, bool const force_move = false)
    {
      size_type const N(size());
      for (size_type i(0U); i < N; ++i)
      {
        if (force_move || !atoms(i).fixed()) m_representation.structure.cartesian[i] += p;
      }
      energy_valid = false;
    }

    /** Set Internals
  @param ri: internals to which the internal coordinates of the coordintes object should be set
  @param force_move: if set to true also move fixed atoms*/
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

    void set_dih(size_type const int_index, coords::angle_type const target_angle,
      bool const move_dependants_along = true, bool const move_fixed_dih = false);
    void rotate_dih(size_type const int_index, coords::angle_type const rot_angle,
      bool const move_dependants_along = true, bool const move_fixed_dih = false);
    void rotate_main(size_type const main_index, coords::angle_type const rot_angle,
      bool const move_dependants_along = true, bool const move_fixed_dih = false);
    void set_all_main(Representation_Main const & new_values, bool const aplly_to_xyz = true,
      bool const move_dependants_along = true, bool const move_fixed_dih = false);

    /** set gradients of a given atom to a given value
  @param index: atom index
  @param p: new gradients*/
    void update_g_xyz(size_type const index, Cartesian_Point const &p)
    {
      if (!atoms(index).fixed()) m_representation.gradient.cartesian[index] = p;
    }

    /** add additional gradients to all atoms
  @param rep: gradient object that is to be added*/
    void sum_g_xyz(Gradients_3D const & rep)
    {
      m_representation.gradient.cartesian += rep;
      zero_fixed_g();
    }

    /** set gradients of all fixed atoms to zero*/
    void zero_fixed_g()
    {
      size_type const N(size());
      for (size_type i(0U); i < N; ++i)
      {
        if (atoms(i).fixed()) m_representation.gradient.cartesian[i] = cartesian_gradient_type(0);
      }
    }

    /** put in pes point*/
    void set_pes(PES_Point const & pes_point, bool const overwrite_fixed = false);
    /** put in pes point*/
    void set_pes(PES_Point && pes_point, bool const overwrite_fixed = false);

    /**set new cartisian coordinates
    @param new_xyz: new cartesian coordinates
    @param overwrite_fixed: if true also change coordinates of fixed atoms*/
    void set_xyz(Representation_3D const & new_xyz, bool const overwrite_fixed = false)
    {
      size_type const N(size());
      if (new_xyz.size() != N)
      {
        throw std::logic_error("Wrong sized coordinates in set_xyz.");
      }
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
    //Irgendwer der das perfect forwarden will? Der move macht naemlich probleme
    //ifdef kann weg
#if defined(SCON_CC11_RVALUE_REF) && defined(SCON_CC11_MOVE)
  /**set new cartisian coordinates
  @param new_xyz: new cartesian coordinates
  @param overwrite_fixed: if true also change coordinates of fixed atoms*/
    void set_xyz(Representation_3D && new_xyz, bool const overwrite_fixed = false)
    {
      size_type const N(size());
      if (new_xyz.size() != N) throw std::logic_error("Wrong sized coordinates in set_xyz.");
      if (!overwrite_fixed)
      {
        for (size_type i(0U); i < N; ++i)
        {
          if (!atoms(i).fixed()) m_representation.structure.cartesian[i] = std::move(new_xyz[i]);
        }
      }
      else
        m_representation.structure.cartesian = std::move(new_xyz);
      m_stereo.update(xyz());
    }
#endif

    /** set new gradients and update internal coordinates
  @param rhs: new gradients
  @param overwrite_fixed: if true also change gradients of fixed atoms*/
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
    /**sets another Gradients_3D object to the cartesian gradients of coordinates object
    @param out_g_xyz: name of object that should be set to the gradients*/
    void get_g_xyz(Gradients_3D & out_g_xyz) const
    {
      out_g_xyz = m_representation.gradient.cartesian;
    }
    /**delete gradients -> empty object*/
    void clear_g_xyz()
    {
      m_representation.gradient.cartesian.assign(size(), Cartesian_Point());
    }

    /**returns a given atom object
    @param index: atom index*/
    Atom const & atoms(size_type const index) const
    {
      return m_atoms.atom(index);
    }
    /**returns all atoms*/
    Atoms const & atoms() const
    {
      return m_atoms;
    }

    /**returns biased potentials*/
    bias::Potentials & potentials()
    {
      return m_potentials;
    }
    /**returns biased potentials*/
    bias::Potentials const & potentials() const
    {
      return m_potentials;
    }

    /**returns matrix with interactions between subsystems*/
    sub_ia_matrix_t & interactions()
    {
      return m_representation.ia_matrix;
    }
    /**returns matrix with interactions between subsystems*/
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

    /**converts cartesian to internal coordinates*/
    void to_internal() { m_atoms.c_to_i(m_representation); }
    /**converts cartesian to internal coordinates (light???)*/
    void to_internal_light() { m_atoms.c_to_i_light(m_representation); }

    /**converts internal to cartesian coordinates*/
    void to_xyz()
    {
      m_atoms.i_to_c(m_representation);
      if (Config::get().periodics.periodic)
      {
        periodic_boxjump();
      }
      m_stereo.update(xyz());
    }

    bool check_superposition_xyz(Representation_3D const &a,
      Representation_3D const &b, double const x = 0.35) const;

    bool is_equal_structure(coords::PES_Point const &a, coords::PES_Point const &b) const;
    //returns if the atom is terminal for every atom
    std::vector<bool> const terminal();
    //returns 1 for terminal atoms, 2 for atoms that are terminal when ignoring actually terminal atoms, etc.
    std::vector<size_t> const terminal_enum();
    //returns a Coordinates object reduced by all atoms with bool "false"
    coords::Coordinates get_red_replic(std::vector<bool> criterion);
    //for changing atoms list. Only to be used with replics of importent coords::Coordinates objects
    Atoms & atoms_changable()
    {
      return m_atoms;
    }
    Atom & atoms_changeable(size_type const index)
    {
      return m_atoms.atom(index);
    }
    //adapts indexation of coords object to initially read structure
    void adapt_indexation(size_t no_dist, size_t no_angle, size_t no_dihedral,
      std::vector<std::vector<std::pair<std::vector<size_t>, double>>> const &reference,
      coords::Coordinates const *cPtr);

    //returns maximal found values of cartesian coordiantes as a Cartesian_Point for fixed atoms
    Cartesian_Point max_valuePosfix()
    {
      Cartesian_Point maxV;
      bool check_fix = false;

      maxV = m_representation.structure.cartesian[0];

      for (std::size_t i=1u;i < m_atoms.size();i++)
      {
        if (m_atoms.check_fix(i) == true)
        {
          check_fix=true;
          if (m_representation.structure.cartesian[i].x() > maxV.x()) { maxV.x() = m_representation.structure.cartesian[i].x(); }
          if (m_representation.structure.cartesian[i].y() > maxV.y()) { maxV.y() = m_representation.structure.cartesian[i].y(); }
          if (m_representation.structure.cartesian[i].z() > maxV.z()) { maxV.z() = m_representation.structure.cartesian[i].z(); }
        }
      }
      if(check_fix == false){maxV.x()=0.0; maxV.y() = 0.0; maxV.z() = 0.0;}

      return maxV;
    }

    Cartesian_Point min_valuePosfix()
    {
      Cartesian_Point minV;
      bool check_fix = false;

      minV = m_representation.structure.cartesian[0];

      for (std::size_t i = 1u; i < m_atoms.size(); i++)
      {
        if (m_atoms.check_fix(i) == true)
        {
          check_fix = true;
          if (m_representation.structure.cartesian[i].x() < minV.x()) { minV.x() = m_representation.structure.cartesian[i].x(); }
          if (m_representation.structure.cartesian[i].y() < minV.y()) { minV.y() = m_representation.structure.cartesian[i].y(); }
          if (m_representation.structure.cartesian[i].z() < minV.z()) { minV.z() = m_representation.structure.cartesian[i].z(); }
        }
      }
      if (check_fix == false) { minV.x() = 0.0; minV.y() = 0.0; minV.z() = 0.0;}

      return minV;
    }
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

  /**
   * Purpose: Creates a callback for optimizer,
   * since they don't care.
   *
   *
   */
  struct Coords_3d_float_callback
  {
    coords::Coordinates * cp;
    Coords_3d_float_callback(coords::Coordinates & coordpointer) :
      cp(&coordpointer)
    { }

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
