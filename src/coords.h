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
#include "Scon/scon_vect.h"
#include "Scon/scon_utility.h"
#include "Scon/scon_matrix.h"
#include "coords_rep.h"
#include "Scon/scon_log.h"
#include "coords_atoms.h"

// forward declaration of different Coordinate Classes

namespace md {
  class CoordinatesUBIAS;
}

namespace optimization {
  namespace global {
    class CoordsOptimizationTS;
  }
}

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
    Stereo(Atoms const&, coords::Representation_3D const&);

    std::vector<pair> const& centers() const { return m_centers; }
    void update(coords::Representation_3D const&);
    void swap(Stereo& r) { m_centers.swap(r.m_centers); }
    void clear() { m_centers.clear(); }

    friend std::ostream& operator<< (std::ostream& stream, Stereo const& stereo);

  private:

    std::vector< pair > m_centers;

  };

  inline bool operator== (Stereo::pair const& a, Stereo::pair const& b)
  {
    return a.m_atom == b.m_atom && a.m_dir == b.m_dir;
  }
  inline bool operator!= (Stereo::pair const& a, Stereo::pair const& b)
  {
    return !(a == b);
  }

  inline bool operator== (Stereo const& a, Stereo const& b)
  {
    std::size_t const N(a.centers().size());
    if (N != b.centers().size()) return false;
    for (std::size_t i(0); i < N; ++i)
    {
      if (a.centers()[i] != b.centers()[i]) return false;
    }
    return true;
  }
  inline bool operator!= (Stereo const& a, Stereo const& b)
  {
    return !(a == b);
  }

  inline void swap(Stereo& a, Stereo& b)
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

  /**namespace that contains only bias potentials*/
  namespace bias
  {
    /**struct for every kind of bias potentials*/
    struct Potentials
    {
      /**constructor, fetches biases from config object*/
      Potentials();

      /**are there any biases?*/
      bool empty() const;
      /**are there any umbrella restraints?*/
      bool uempty() const { return m_utors.empty() && m_udist.empty() && m_ucombs.empty() && m_uangles.empty(); }
      /**returns bias energy*/
      double energy() const { return b + a + d + s + c; }

      /**clear all biases*/
      void clear()
      {
        b = a = d = s = c = thr = u = 0.0;
        scon::clear(m_dihedrals, m_angles, m_distances,
          m_spherical, m_cubic, m_utors, m_udist, m_ucombs, m_thresh, m_uangles);
      }

      /**energy for distance biases*/
      double e_dist() const { return b; }
      /**energy for angle biases*/
      double e_angle() const { return a; }
      /**energy for dihedral biases*/
      double e_dihedral() const { return d; }
      /**energy for spherical biases*/
      double e_spherical() const { return s; }
      /**energy for cubic biases*/
      double e_cubic() const { return c; }
      /**energy for threshold biases*/
      double e_thresh() const { return thr; }
      /**energy for umbrella combination biases*/
      double e_ucomb() const { return u; }

      /**add a new dihedral bias*/
      void add(config::biases::dihedral const& new_d) { m_dihedrals.push_back(new_d); }
      /**add a new angle bias*/
      void add(config::biases::angle const& new_a) { m_angles.push_back(new_a); }
      /**add a new distance bias*/
      void add(config::biases::distance const& new_d) { m_distances.push_back(new_d); }
      /**add a new spherical bias*/
      void add(config::biases::spherical const& new_d) { m_spherical.push_back(new_d); }
      /**add a new cubic bias*/
      void add(config::biases::cubic const& new_d) { m_cubic.push_back(new_d); }
      /**add a new threshold bias*/
      void add(config::biases::thresholdstr const& new_thr) { m_thresh.push_back(new_thr); }

      /**return all dihedral biases*/
      std::vector<config::biases::dihedral> const& dihedrals() const { return m_dihedrals; }
      /**return all angle biases*/
      std::vector<config::biases::angle> const& angles() const { return m_angles; }
      /**return all distance biases*/
      std::vector<config::biases::distance> const& distances() const { return m_distances; }
      /**return all spherical biases*/
      std::vector<config::biases::spherical> const& sphericals() const { return m_spherical; }
      /**return all cubic biases*/
      std::vector<config::biases::cubic> const& cubic() const { return m_cubic; }
      /**return all threshold biases*/
      std::vector<config::biases::thresholdstr> const& thresholds() const { return m_thresh; }
      /**return all umbrella combination biases*/
      std::vector<config::coords::umbrellas::umbrella_comb> const& ucombs() const { return m_ucombs; }
      /**function to change umbrella combinations
      necessary for raising force constant during equilibration*/
      std::vector<config::coords::umbrellas::umbrella_comb>& set_ucombs() { return m_ucombs; }

      /**apply "normal" biases (dihedrals, angles, distances, spherical, cubic, threshold, if desired umbrella_combs)
      returns bias energy*/
      double apply(Representation_3D const& xyz, Representation_3D& g_xyz,
        Cartesian_Point maxPos, Cartesian_Point minPos, Cartesian_Point const& center = Cartesian_Point());
      /**apply umbrella biases, i.e. apply bias gradients and fill values for reaction coordinate into uout
      @param xyz: cartesian coordinates of molecule
      @param g_xyz: cartesian gradients of molecule (are changed according to bias)
      @param iout: vector of values for umbrella reaction coordinate that are later written into 'umbrella.txt'*/
      void umbrellaapply(Representation_3D const& xyz,
        Representation_3D& g_xyz, std::vector<double>& uout);

      /**get biases from config*/
      void append_config();

      void swap(Potentials& rhs);

    private:

      //energies of biases
      /**energy of distance bias*/
      double b;
      /**energy of angle bias*/
      double a;
      /**energy of dihedral bias*/
      double d;
      /**energy of spherical bias*/
      double s;
      /**energy of cubic bias*/
      double c;
      /**energy of threshold*/
      double thr;
      /**energy of bottom-threshold*/
      double thrB;
      /**energy of umbrella combination bias*/
      double u;

      //biases

      /**distance biases*/
      std::vector<config::biases::distance>  m_distances;
      /**angle biases*/
      std::vector<config::biases::angle>     m_angles;
      /**dihedral biases*/
      std::vector<config::biases::dihedral>  m_dihedrals;
      /**spherical biases*/
      std::vector<config::biases::spherical> m_spherical;
      /**cubic biases*/
      std::vector<config::biases::cubic>     m_cubic;
      /**threshold biases*/
      std::vector<config::biases::thresholdstr>  m_thresh;
      std::vector<config::biases::thresholdstr>  m_threshBottom;
      // umbrella biases
      /**distance biases for umbrella*/
      std::vector<config::coords::umbrellas::umbrella_dist> m_udist;
      /**angle biases for umbrella*/
      std::vector<config::coords::umbrellas::umbrella_angle> m_uangles;
      /**dihedral biases for umbrella*/
      std::vector<config::coords::umbrellas::umbrella_tor> m_utors;
      /**linear combinations biases for umbrella (can also be used as a "normal" bias)*/
      std::vector<config::coords::umbrellas::umbrella_comb> m_ucombs;

      // applying biases

      /**function to apply bias potential on dihedral*/
      double dih(Representation_3D const& xyz, Gradients_3D& g_xyz);
      /**function to apply bias potential on distance*/
      double dist(Representation_3D const& xyz, Gradients_3D& g_xyz);
      /**function to apply bias potential on angle*/
      double ang(Representation_3D const& xyz, Gradients_3D& g_xyz);
      /**function to apply spherical potential*/
      double spherical(Representation_3D const& xyz, Gradients_3D& g_xyz,
        Cartesian_Point const& center = Cartesian_Point());
      /**function to apply cubic potential*/
      double cubic(Representation_3D const& xyz, Gradients_3D& g_xyz,
        Cartesian_Point const& center = Cartesian_Point());
      /**function to apply a restraint on an umbrella dihedral and saves the values for the 'umbrella.txt' file
      @param xyz: coordinates of system
      @param g_xyz: cartesian gradients of system
      @param uout: vector with the real values for the restraint coordinate*/
      void umbrelladih(Representation_3D const& xyz, Gradients_3D& g_xyz, std::vector<double>& uout) const;
      /**function to apply a restraint on an umbrella angle and saves the values for the 'umbrella.txt' file
      @param xyz: coordinates of system
      @param g_xyz: cartesian gradients of system
      @param uout: vector with the real values for the restraint coordinate*/
      void umbrellaang(Representation_3D const& xyz, Gradients_3D& g_xyz, std::vector<double>& uout) const;
      /**function to apply a restraint on an umbrella distance and saves the values for the 'umbrella.txt' file
      @param xyz: coordinates of system
      @param g_xyz: cartesian gradients of system
      @param uout: vector with the real values for the restraint coordinate*/
      void umbrelladist(Representation_3D const& xyz, Gradients_3D& g_xyz, std::vector<double>& uout)  const;
      /**function to apply a restraint on an umbrella combination of distances and saves the values for the 'umbrella.txt' file
      @param xyz: coordinates of system
      @param g_xyz: cartesian gradients of system
      @param uout: vector with the real values for the restraint coordinate*/
      void umbrellacomb(Representation_3D const& xyz, Gradients_3D& g_xyz, std::vector<double>& uout);
      /**function to apply a umbrella combination biases on energy and gradients outside of umbrella sampling
      @param xyz: coordinates of system
      @param g_xyz: cartesian gradients of system
      returns the additional bias energy*/
      double umbrellacomb(Representation_3D const& xyz, Gradients_3D& g_xyz);
      /**function to apply threshold potential*/
      double thresh(Representation_3D const& xyz, Gradients_3D& g_xyz, Cartesian_Point maxPos);
      double thresh_bottom(Representation_3D const& xyz, Gradients_3D& g_xyz, Cartesian_Point minPos);
    };
  }

  /**all important information collected during FEP run*/
  struct fep_data
  {
    /**object in which data for one conformation is saved
    which is relevant for FEP calculation*/
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

  /**set viral tensor to all zeroes*/
  static inline coords::virial_t empty_virial()
  {
    coords::virial_t v;
    for (std::size_t i = 0; i < 3; ++i)
    {
      for (std::size_t j = 0; j < 3; ++j) v[i][j] = 0.0;
    }
    return v;
  }
  /**coordinates object
  most important part of CAST program as nearly every task works with one (or more) objects of this class*/
  class Coordinates
  {
    friend class md::CoordinatesUBIAS;
    friend class optimization::global::CoordsOptimizationTS;
    /**atoms (without xyz coordinates)*/
    Atoms                     m_atoms;

    /**other information about system (i.e. xyz coordinates)*/
  protected:
    PES_Point                 m_representation;
  private:
    /**stereo information*/
    Stereo                    m_stereo;
    /**external potentials*/
  protected:
    bias::Potentials          m_potentials;
  private:

    /**virial (needed for pressure)*/
    virial_t                  m_virial;
    /**pointer to energy interface*/
    energy::interface_base* m_interface;
    /**pointer to energy interface that is preinterface*/
    energy::interface_base* m_preinterface;
    /**???*/
    bool                      energy_valid;
    /**number of iterations needed for last optimization*/
    std::size_t               m_iter{ 0u };

    /**lbfgs optimizer with preinterface, returns energy of optimized structure*/
    std::pair<coords::float_type, std::size_t> prelbfgs();
    /**lbfgs optimizer, returns energy of optimized structure and number of iterations*/
    std::pair<coords::float_type, std::size_t> lbfgs();

    /**function to calculate energy*/
    coords::float_type m_e(energy::interface_base* const p)
    {
      if (p)
      {
        energy_valid = true;
        if (Config::get().periodics.periodic)
          periodic_boxjump_prep();
        m_representation.energy = p->e();
        m_representation.integrity = p->intact();
        apply_bias();
        zero_fixed_g();
        return m_representation.energy;
      }
      return coords::float_type();
    }

    /**function to calculate energy and gradients*/
    coords::float_type m_g(energy::interface_base* const p)
    {
      if (p)
      {
        energy_valid = true;
        if (Config::get().periodics.periodic)
          periodic_boxjump_prep();
        m_representation.energy = p->g();
        m_representation.integrity = p->intact();
        this->apply_bias();
        zero_fixed_g();
        if (Config::get().general.verbosity > 4u)
        {
          coords::float_type sumGrad = 0.;
          for (auto&& ele : this->g_xyz())
          {
            coords::float_type thisAtomValueAbs = std::sqrt(
              ele.x() * ele.x()
              + ele.y() * ele.y()
              + ele.z() * ele.z());
            sumGrad += thisAtomValueAbs;
          }
          std::cout << "Obtained gradient. Summed gradient norm of all atoms is: " << std::to_string(sumGrad) << " [kcal/(mol*Angst)].\n";
        }
        return m_representation.energy;
      }
      return coords::float_type();
    }

  public:

    energy::interface_base* catch_interface = m_interface;

    /*void get_catch_interface()
    {
      energy::interface_base   *catch_interface = m_interface;
    }*/

    /**contains all important information collected during FEP run*/
    fep_data const& grtFep() const { return fep; } // only used in MD and read in Coordinates for ACO and Amoeba
    fep_data& getFep() { return fep; }
    fep_data              fep;


    bool                  NEB_control, PathOpt_control;

    std::size_t           mult_struc_counter;
    //aditional variables for perpendicular g() and o() and NEB

    // Orthogonalize is used in Dimer method
    bool									orthogonalize;

    typedef PES_Point::size_type size_type;

    /**default constructor*/
    Coordinates();
    /**copy constructor*/
    Coordinates(Coordinates const& r);
    /**move constructor)*/
    Coordinates(Coordinates&& r);
    /**move assign operator*/
    Coordinates& operator= (Coordinates&& rhs);
    /**copy assign operator*/
    Coordinates& operator= (Coordinates const& rhs);
    /**destructor*/
    virtual ~Coordinates();

    void apply_bias() {
      if (!m_potentials.empty()) {
        m_representation.energy +=
          m_potentials.apply(m_representation.structure.cartesian,
            m_representation.gradient.cartesian,
            max_valuePosfix(), Cartesian_Point());
      }
    }

    /**move atom
    @param index: index of atom that is to be moved
    @param p: "space vector" by which it should be moved
    @param force_move: if set to true also move fixed atoms*/

    /**fix an atom, i.e. this atom can't be moved
    @param atom: atom index*/
    void set_fix(size_t const atom, bool const fix_it = true);

    /**set all atoms to their original fixation state, given in inputfile*/
    void reset_fixation();

    /**delete everything in the Coordinates object -> empty object*/
    void clear()
    {
      *this = Coordinates{};
    }

    void move_atom_by(std::size_t const index, coords::cartesian_type const& p, bool const force_move = false)
    {
      if (!atoms(index).fixed() || force_move)
      {
        m_representation.structure.cartesian[index] += p;
        energy_valid = false;
        m_stereo.update(xyz());
      }
    }


    /**function to get bias potentials*/
    bias::Potentials& get_biases() { return m_potentials; }

    /**function that applies umbrella potentials and adds values to uout (later written in 'umbrella.txt')*/
    void ubias(std::vector<double>& uout)
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
        else {
          auto lbfgs_result = prelbfgs();
          m_representation.energy = lbfgs_result.first;  // energy
          m_iter = lbfgs_result.second;                  // number of optmization steps
        }
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
        && m_potentials.empty()                // no bias
        && !Config::get().periodics.periodic)  // no periodic boundaries
      {
        m_representation.energy = m_interface->o();
      }
      else
      {
        auto lbfgs_result = lbfgs();
        m_representation.energy = lbfgs_result.first;  // energy
        m_iter = lbfgs_result.second;                  // number of optmization steps
      }
      m_representation.integrity = m_interface->intact();
      m_stereo.update(xyz());
      zero_fixed_g();     // sets gradients of all fixed atoms to zero
      return m_representation.energy;
    }

    /**calculate hessian matrix*/
    coords::float_type h()
    {
      return m_representation.energy = m_interface->h();
    }

    bool preoptimize() const { return m_preinterface ? true : false; }

    /**returns the energy interface*/
    energy::interface_base* energyinterface() const { return m_interface; }
    /**returns the energy preinterface*/
    energy::interface_base* preinterface() const { return m_preinterface; }

    /**determines if structure is intact,
    i.e. bond lengths do not differ from ideal bond length by more than 5 angstrom,
    angles do not differ from ideal angle by more than 30 degrees
    and similar criteria for dihedrals and ureys
    details see in energy_int_..._pot.cc
    for other energy interfaces integrity is destroyed when external program can't calculate an energy*/
    bool               integrity() const { return pes().integrity; }
    /**returns the size of the coordinates object, i.e. the total number of atoms*/
    size_type          size() const { return m_atoms.size(); }
    /**object swap*/
    void               swap(Coordinates& rhs);


    /**checks if structure is complete, i.e. no coordinates are NaN
        coordinates become NaN sometimes in TS (dimer method)*/
    bool check_structure();
    /**checks if non-bound atoms in structure are crashing,
    i.e. their distance is less than 1.2 times sum of covalent radii though there is no bond between them
    true means everything is okay, false means there are crashes*/
    bool check_for_crashes() const;
    /**checks if all bonds are still intact (bond length smaller than 1.2 sum of covalent radii but bigger than 0.3 Angstrom)*/
    bool check_bond_preservation() const;
    /**looks if currently a valid z-matrix exists*/
    bool has_valid_internals() { return m_atoms.z_matrix_valid(); }

    /**saves virial coefficients into coordinates object
    @param V: virials coefficients*/
    void set_virial(virial_t const& V) { m_virial = V; }

    /**calculates the center of mass*/
    Cartesian_Point    center_of_mass() const;
    /**calculates the center of mass for one molecule
    @param index: index of the molecule*/
    Cartesian_Point    center_of_mass_mol(size_type const index) const;
    /**calculates the geometrical center*/
    Cartesian_Point    center_of_geometry() const;
    /**calculates the total mass of the system*/
    coords::float_type weight() const;


    /**function that rebinds atoms after a geometry change
    due to a distance criterion (d < 1.2 * sum of covalent radii)*/
    void rebind();

    /**fills the coordinates object with data
    @param a: atoms object
    @param p: PES_point*/
    void init_swap_in(Atoms& a, PES_Point& p, bool const energy_update = true);
    /**does the same as init_swap_in
    @param a: atoms object
    @param p: PES_point*/
    template<typename T, typename U> // Without templates some calls might end with an unwanted copy ctor call.
    void init_in(T&& a, U&& p, bool const energy_update = true) {
      init_swap_in(a, p, energy_update);
    }
    /**updates the topology*/
    void energy_update(bool const skip_topology = false) { m_interface->update(skip_topology); }
    /**fixes all atoms, i.e. nothing can move anymore*/
    void fix_all(bool const fix_it = true) { m_atoms.fix_all(fix_it); }
    /**fixes a given atom
    @param index: index of atom that is to be fixed*/
    void fix(size_type const index, bool const fix = true) { m_atoms.fix(index, fix); }

    void e_head_tostream_short(std::ostream& strm, energy::interface_base const* const ep = nullptr) const;
    void e_tostream_short(std::ostream& strm, energy::interface_base const* const ep = nullptr) const;
    /**writes hessian matrix
    @param strm: can be std::cout or ofstream file*/
    void h_tostream(std::ostream& strm) const;

    /**returns the PES point*/
    PES_Point const& pes() const { return m_representation; }
    /**also returns the PES point*/
    PES_Point& pes() { return m_representation; }

    /** returns the xyz representation*/
    Representation_3D const& xyz() const
    {
      return m_representation.structure.cartesian;
    }
    /** returns the internal representation*/
    Representation_Internal const& intern() const
    {
      return m_representation.structure.intern;
    }
    /** returns the representation of main dihedrals*/
    Representation_Main const& main() const
    {
      return m_representation.structure.main;
    }

    /**returns the xyz coordinates of a given atom
  @param index: atom index*/
    cartesian_type const& xyz(size_type index) const
    {
      return m_representation.structure.cartesian[index];
    }
    /**returns the internal coordinates of a given atom
    @param index: atom index*/
    internal_type const& intern(size_type index) const
    {
      return m_representation.structure.intern[index];
    }
    /**returns the main dihedral coordinates of a given atom
    @param index: atom index*/
    main_type const& main(size_type index) const
    {
      return m_representation.structure.main[index];
    }

    /**returns the gradients in cartesian space*/
    Gradients_3D const& g_xyz() const
    {
      return m_representation.gradient.cartesian;
    }
    /**returns the gradients in internal space*/
    Gradients_Internal const& g_intern() const
    {
      return m_representation.gradient.intern;
    }
    /**returns the gradients in main dihedral space*/
    Gradients_Main const& g_main() const
    {
      return m_representation.gradient.main;
    }

    /**returns the gradients in cartesian space for a given atom
    @param index: atom index*/
    cartesian_gradient_type const& g_xyz(size_type const index) const
    {
      return m_representation.gradient.cartesian[index];
    }
    /**returns the gradients in internal space for a given atom
    @param index: atom index*/
    internal_gradient_type const& g_intern(size_type const index) const
    {
      return m_representation.gradient.cartesian[index];
    }
    /**returns the gradients in main dihedral space for a given atom
    @param index: atom index*/
    main_gradient_type const& g_main(size_type const index) const
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
    size_2d const& subsystems() const { return m_atoms.subsystems(); }
    /**returns a given subsystem
    @param index: index of subsystem*/
    size_1d const& subsystems(size_type index) const { return m_atoms.subsystems(index); }
    /**returns all molecules*/
    size_2d const& molecules() const { return m_atoms.molecules(); }
    /**returns a given molecule
    @param index: index of molecule*/
    size_1d const& molecule(size_type index) const { return m_atoms.molecule(index); }

    /**returns stereo centers?*/
    std::vector< Stereo::pair > const& stereos() const { return m_stereo.centers(); }

    /**returns name of the given molecule (at the moment only water and sodium)*/
    std::string molecule_name(Container<std::size_t> const& molecule) const;

    /**if periodic boundaries are activated:
    preperation for moving molecules that are outside of the box into the box
    (for QM/MM interfaces: whole QM part consists of one molecule)*/
    void periodic_boxjump_prep();

    /**if periodic boundaries are activated:
    move molecules that are outside of the box into the box*/
    void periodic_boxjump(std::vector<std::vector<std::size_t>> const& molecules);

    /**move given atom to given coordinates
    @param index: index of atom that is to be moved
    @param p: coordinates where the atom should be moved to
    @param force_move: if set to true also move fixed atoms*/
    void move_atom_to(size_type const index, Cartesian_Point const& p, bool const force_move = false)
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
    void move_all_by(Cartesian_Point const& p, bool const force_move = false)
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
    void set_internal(Representation_Internal const& ri, bool const force_move = false)
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
    void set_all_main(Representation_Main const& new_values, bool const aplly_to_xyz = true,
      bool const move_dependants_along = true, bool const move_fixed_dih = false);

    /** set gradients of a given atom to a given value
  @param index: atom index
  @param p: new gradients*/
    void update_g_xyz(size_type const index, Cartesian_Point const& p)
    {
      if (!atoms(index).fixed()) m_representation.gradient.cartesian[index] = p;
    }

    /** add additional gradients to all atoms
  @param rep: gradient object that is to be added*/
    void sum_g_xyz(Gradients_3D const& rep)
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
    void set_pes(PES_Point const& pes_point, bool const overwrite_fixed = false);
    /** put in pes point*/
    void set_pes(PES_Point&& pes_point, bool const overwrite_fixed = false);

    /**set new cartisian coordinates
    @param new_xyz: new cartesian coordinates
    @param overwrite_fixed: if true also change coordinates of fixed atoms*/
    void set_xyz(Representation_3D const& new_xyz, bool const overwrite_fixed = false)
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
    void set_xyz(Representation_3D&& new_xyz, bool const overwrite_fixed = false)
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
    void swap_g_xyz(Gradients_3D& rhs, bool const overwrite_fixed = false)
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

    template<typename Grad3D>
    void set_g_xyz(Grad3D&& new_grads, bool const overwrite_fixed = false) {
      size_type const N(m_representation.gradient.cartesian.size());
      if (new_grads.size() != N) throw std::logic_error("Wrong sized gradients in set_g_xyz.");
      auto old_grads = std::move(m_representation.gradient.cartesian);
      m_representation.gradient.cartesian = std::forward<Grad3D>(new_grads);
      if (!overwrite_fixed)
      {
        for (size_type i(0U); i < N; ++i)
        {
          if (atoms(i).fixed()) m_representation.gradient.cartesian[i] = old_grads[i];
        }
      }
      m_atoms.c_to_i(m_representation); // update internals
    }

    /**sets another Gradients_3D object to the cartesian gradients of coordinates object
    @param out_g_xyz: name of object that should be set to the gradients*/
    void get_g_xyz(Gradients_3D& out_g_xyz) const
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
    Atom const& atoms(size_type const index) const
    {
      return m_atoms.atom(index);
    }

    /**returns all atoms*/
    Atoms const& atoms() const
    {
      return m_atoms;
    }

    /**returns biased potentials*/
    bias::Potentials& potentials()
    {
      return m_potentials;
    }
    /**returns biased potentials*/
    bias::Potentials const& potentials() const
    {
      return m_potentials;
    }

    /**returns matrix with interactions between subsystems*/
    coords::sub_ia_matrix_t& interactions()
    {
      return m_representation.ia_matrix;
    }
    /**returns matrix with interactions between subsystems*/
    coords::sub_ia_matrix_t const& interactions() const
    {
      return m_representation.ia_matrix;
    }
    sub_ia& interactions(size_type const index)
    {
      return m_representation.ia_matrix(index);
    }
    sub_ia const& interactions(size_type const index) const
    {
      return m_representation.ia_matrix(index);
    }
    sub_ia& interactions(size_type const x, size_type const y)
    {
      return m_representation.ia_matrix(x, y);
    }
    sub_ia const& interactions(size_type const x, size_type const y) const
    {
      return m_representation.ia_matrix(x, y);
    }
    /**get number of optimization steps*/
    std::size_t get_opt_steps() const { return m_iter; }

    /**converts coordinates first to internal coordinates then back to cartesian ones
    if conversion to internals fails, no backwards conversion is performed*/
    void to_internal_to_xyz() {
      to_internal();
      if (has_valid_internals()) to_xyz();
      else std::cout << "No back and forth conversion of coordinates possible!\n";
    };

    /**converts cartesian to internal coordinates*/
    void to_internal() { return m_atoms.c_to_i(m_representation); }
    /**converts cartesian to internal coordinates but without gradients*/
    void to_internal_light() { return m_atoms.c_to_i_light(m_representation); }

    /**converts internal to cartesian coordinates*/
    void to_xyz()
    {
      if (!has_valid_internals())
        std::cout << "WARNING! Internal coordinates are broken. The cartesian coordinates resulting of this conversion might be crap.\n";
      m_atoms.i_to_c(m_representation);
      if (Config::get().periodics.periodic)
      {
        periodic_boxjump_prep();
      }
      m_stereo.update(xyz());
    }

    bool check_superposition_xyz(Representation_3D const& a,
      Representation_3D const& b, double const x = 0.35) const;

    /**checks if two structures are equal*/
    bool is_equal_structure(coords::PES_Point const& a, coords::PES_Point const& b) const;
    //returns if the atom is terminal for every atom
    std::vector<bool> const terminal();
    //returns 1 for terminal atoms, 2 for atoms that are terminal when ignoring actually terminal atoms, etc.
    std::vector<size_t> const terminal_enum();
    //returns a Coordinates object reduced by all atoms with bool "false"
    coords::Coordinates get_red_replic(std::vector<bool> criterion);
    //for changing atoms list. Only to be used with replics of importent coords::Coordinates objects
    Atoms& atoms_changable()
    {
      return m_atoms;
    }
    Atom& atoms_changeable(size_type const index)
    {
      return m_atoms.atom(index);
    }
    //adapts indexation of coords object to initially read structure
    void adapt_indexation(std::vector<std::vector<std::pair<std::vector<size_t>, double>>> const& reference,
      coords::Coordinates const* cPtr);

    //returns maximal found values of cartesian coordiantes as a Cartesian_Point for fixed atoms used for thresh potential

    Cartesian_Point max_valuePosfix()
    {
      Cartesian_Point maxV;

      maxV = m_representation.structure.cartesian[0];//maxV must contain values so the compare in the loop works
      bool check_fix = false; //the test for fixed atoms is done with check_fix this way so the information from an fixed atom is not overwritten by a non fixed atom with a higher index in the m_atoms object.
      for (std::size_t i = 1u; i < m_atoms.size(); i++)
      {

        if (m_atoms.check_fix(i) == true)
        {
          check_fix = true;
          if (m_representation.structure.cartesian[i].x() > maxV.x()) { maxV.x() = m_representation.structure.cartesian[i].x(); }
          if (m_representation.structure.cartesian[i].y() > maxV.y()) { maxV.y() = m_representation.structure.cartesian[i].y(); }
          if (m_representation.structure.cartesian[i].z() > maxV.z()) { maxV.z() = m_representation.structure.cartesian[i].z(); }
        }
      }
      if (check_fix == false) { maxV.x() = 0.0; maxV.y() = 0.0; maxV.z() = 0.0; }

      std::cout << "Reference position for threshold potential is: " << maxV.x() << "  " << maxV.y() << "  " << maxV.z() << '\n';
      return maxV;
    }

    //returns minimal found values of cartesian coordiantes as a Cartesian_Point for fixed atoms used for thresh_bottom potential
    Cartesian_Point min_valuePosfix()
    {
      Cartesian_Point minV;

      minV = m_representation.structure.cartesian[0];
      bool check_fix = false;
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
      if (check_fix == false) { minV.x() = 0.0; minV.y() = 0.0; minV.z() = 0.0; }

      std::cout << "Reference position for threshold potential bottom is: " << minV.x() << "  " << minV.y() << "  " << minV.z() << '\n';

      return minV;
    }
  };

  std::ostream& operator<< (std::ostream& stream, Coordinates const& coord);

  struct internal_float_callback
  {
    coords::Coordinates* cp;
    internal_float_callback(coords::Coordinates& coordinates_object)
      : cp(&coordinates_object)
    { }
    float operator() (scon::vector<float> const& v,
      scon::vector<float>& g, std::size_t const S, bool& go_on)
    {
      std::size_t i = 0;
      coords::Representation_Internal rin(v.size() / 3);
      for (auto& e : rin)
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
      for (auto const& e : cp->g_intern())
      {
        g[i++] = static_cast<float>(e.x());
        g[i++] = static_cast<float>(e.y());
        g[i++] = static_cast<float>(e.z());
      }
      if (Config::get().general.verbosity >= 4)
      {
        std::cout << "Optimization: Energy of step " << S;
        std::cout << " is " << E;
        if (go_on)
          std::cout << " with intact integrity.\n";
        else
          std::cout << " with BROKEN INTEGRITY!\n";
      }
      return E;
    }
    scon::vector<float> from_rep(coords::Representation_Internal const& v)
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
    scon::vector<float> from_grad(coords::Gradients_Internal const& v)
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
    coords::Representation_Internal to_rep(scon::vector<float> const& v)
    {
      coords::Representation_Internal r(v.size() / 3);
      std::size_t i = 0;
      for (auto& e : r)
      {
        auto radius = v[i++];
        auto inclination = v[i++];
        auto azimuth = v[i++];
        e = coords::internal_type(radius, coords::angle_type(inclination), coords::angle_type(azimuth));
      }
      return r;
    }
    coords::Gradients_Internal to_grad(scon::vector<float> const& v)
    {
      coords::Gradients_Internal r(v.size() / 3);
      std::size_t i = 0;
      for (auto& e : r)
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
    coords::Coordinates* cp;
    cartesian_float_callback(coords::Coordinates& coordinates_object)
      : cp(&coordinates_object)
    { }
    float operator() (scon::vector<float> const& v,
      scon::vector<float>& g, std::size_t const S, bool& go_on)
    {
      coords::Representation_3D x;
      x.reserve(v.size() / 3);
      std::size_t n = 0;
      for (auto& e : x)
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
      for (auto const& e : cp->g_xyz())
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
    coords::Coordinates* cp;
    std::ostream* strm;
    cartesian_float_stream_log_drain(coords::Coordinates& coords,
      std::ostream& stream) : cp(&coords), strm(&stream) { }
    void operator() (optimization::Point<scon::vector<float>,
      scon::vector<float>, float>& v)
    {
      if (strm && cp)
      {
        auto tmp = cp->xyz();
        coords::Representation_3D x;
        x.reserve(v.x.size() / 3);
        std::size_t n = 0;
        for (auto& e : x)
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
    coords::Coordinates* cp;
    std::ostream* strm;
    internal_float_stream_log_drain(coords::Coordinates& coords,
      std::ostream& stream) : cp(&coords), strm(&stream)
    { }
    void operator() (optimization::Point<scon::vector<float>,
      scon::vector<float>, float>& v)
    {
      if (strm && cp)
      {
        auto p = cp->pes();
        std::size_t i = 0;
        coords::Representation_Internal rin(v.x.size() / 3);
        for (auto& e : rin)
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
    coords::Coordinates* cp;
    Coords_3d_float_callback(coords::Coordinates& coordpointer) :
      cp(&coordpointer)
    { }

    float operator() (scon::vector<scon::c3<float>> const& v,
      scon::vector<scon::c3<float>>& g, std::size_t const S, bool& go_on);

    scon::vector<scon::c3<float>> from(coords::Gradients_3D const& g);

    coords::Representation_3D to(scon::vector<scon::c3<float>> const& v);
  };

  struct Coords_3d_float_pre_callback
  {
    coords::Coordinates* cp;
    Coords_3d_float_pre_callback(coords::Coordinates& coordpointer) : cp(&coordpointer) { }
    float operator() (scon::vector<scon::c3<float>> const& v,
      scon::vector<scon::c3<float>>& g, std::size_t const S, bool& go_on);
  };

  struct Internal_Callback
  {
    coords::Coordinates* cp;
    Internal_Callback(coords::Coordinates& coordpointer)
      : cp(&coordpointer)
    { }
    coords::float_type operator() (scon::vector<coords::float_type> const& v,
      scon::vector<coords::float_type>& g, std::size_t const S, bool& go_on);
    coords::Gradients_Internal from(coords::Representation_Internal const& p);
    coords::Representation_Internal to(coords::Gradients_Internal const& p);
  };

  struct Main_Callback
  {
    coords::Coordinates* cp;
    Main_Callback(coords::Coordinates& coordpointer)
      : cp(&coordpointer)
    { }
    coords::float_type operator() (coords::Gradients_Main const& v,
      coords::Gradients_Main& g, std::size_t const S, bool& go_on);
  };

  /**class for writing the snapshot file during MD*/
  class cartesian_logfile_drain
  {
    /**pointer to coordinates object*/
    coords::Coordinates* cp;
    /**unique pointer to std::ofstream object*/
    std::unique_ptr<std::ofstream> strm;
    /**should structure be optimized before written out?*/
    bool opt;

  public:
    /**default constructor*/
    cartesian_logfile_drain() : cp(), strm(), opt() {}
    /**another constructor*/
    cartesian_logfile_drain(coords::Coordinates& c, char const* const filename, bool optimize = false) :
      cp(&c), strm(std::make_unique<std::ofstream>(std::ofstream(filename, std::ios::out))), opt(optimize) {}
    /**callback operator: writes structure (with given coordinates xyz) into ofstream*/
    void operator() (coords::Representation_3D&& xyz);
  };

  using offset_buffered_cartesian_logfile =
    scon::vector_offset_buffered_callable<coords::Representation_3D, cartesian_logfile_drain>;

  /**I guess this is a function to write the snapshots only after buffer is full*/
  offset_buffered_cartesian_logfile make_buffered_cartesian_log(Coordinates& c,
    std::string file_suffix, std::size_t buffer_size,
    std::size_t log_offset, bool optimize = false);

  /**swap two coordinates objects*/
  inline void swap(Coordinates& a, Coordinates& b)
  {
    a.swap(b);
  }
}

#endif // coords_h_guard_
