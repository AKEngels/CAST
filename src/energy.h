#pragma once

#if defined _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <memory>
#include <string>
#include "coords_rep.h"

namespace energy
{
  static coords::float_type constexpr au2kcal_mol{ 627.5096080306 }; //1 au = 627.5095 kcal/mol
  static coords::float_type constexpr eV2kcal_mol{ 23.061078 };
  static coords::float_type constexpr bohr2ang{ 0.52917721067 };
  static coords::float_type constexpr Hartree_Bohr2Kcal_MolAng{ au2kcal_mol / bohr2ang };
  static coords::float_type constexpr Hartree_Bohr2Kcal_MolAngSquare{ Hartree_Bohr2Kcal_MolAng / bohr2ang };

  /**struct for a point charge (used for QM/MM)*/
  struct PointCharge
  {
    /**position*/
    double x, y, z;
    /**charge (scaling already applied)*/
    double scaled_charge;
    /**charge (without scaling)*/
    double original_charge;

    /**function to set position*/
    void set_xyz(double m_x, double m_y, double m_z)
    {
      x = m_x;
      y = m_y;
      z = m_z;
    }
  };
  /**overloaded output operator for PointCharge*/
  inline std::ostream& operator<< (std::ostream& stream, const PointCharge& c)
  {
    stream << c.x << ", " << c.y << ", " << c.z << ", charge: " << c.scaled_charge << ", original: " << c.original_charge;
    return stream;
  }

  /**object where fep parameters for one window are saved*/
  struct fepvar
  {
    /**lambda_el of former window for appearing atoms*/
    double mein;
    /**lambda_el of former window for disappearing atoms*/
    double meout;
    /**lambda_vdw of former window for appearing atoms*/
    double mvin;
    /**lambda_vdw of former window for disappearing atoms*/
    double mvout;
    /**current lambda_el for appearing atoms*/
    coords::float_type ein;
    /**current lambda_el for disappearing atoms*/
    coords::float_type eout;
    /**current lambda_vdw for appearing atoms*/
    coords::float_type vin;
    /**current lambda_vdw for disappearing atoms*/
    coords::float_type vout;
    /**(lambda + dlambda)_el for appearing atoms*/
    coords::float_type dein;
    /**(lambda + dlambda)_el for disappearing atoms*/
    coords::float_type deout;
    /**(lambda + dlambda)_vdw for appearing atoms*/
    coords::float_type dvin;
    /**(lambda + dlambda)_vdw for disappearing atoms*/
    coords::float_type dvout;
    /**number of current window*/
    int step;
    fepvar(void)
      : ein(1.0), eout(1.0), vin(1.0), vout(1.0),
      dein(1.0), deout(1.0), dvin(1.0), dvout(1.0),
      step(0)
    { }
  };

  /**object in which data for one conformation is saved
  which is relevant for FEP calculation*/
  struct fepvect
  {
    /**coulomb-energy for this conformation with lambda-dlambda*/
    coords::float_type e_c_l0;
    /**vdw-energy for this conformation with lambda-dlambda*/
    coords::float_type e_vdw_l0;
    /**coulomb-energy for this conformation with lambda*/
    coords::float_type e_c_l1;
    /**coulomb-energy for this conformation with lambda+dlambda*/
    coords::float_type e_c_l2;
    /**vdw-energy for this conformation with lambda*/
    coords::float_type e_vdw_l1;
    /**vdw-energy for this conformation with lambda+dlambda*/
    coords::float_type e_vdw_l2;
    /**difference in potential energy for this conformation (calculated in energy_int_aco_pot.cc line 1234)
    = (e_c_l2 + e_vdw_l2) - (e_c_l1 + e_vdw_l1)*/
    coords::float_type dE;
    /**difference in potential energy for backwards transformation
    = (e_c_l1 + e_vdw_l1) - (e_c_l0 + e_vdw_l0)*/
    coords::float_type dE_back;
    /**difference in free energy calculated for all conformations
    in this window until current conformation*/
    coords::float_type dG;
    /**difference in free energy calculated for all conformations
    in this window until current conformation in backwards transformation*/
    coords::float_type dG_back;
    /**exp((-1 / (k_B*T))*dE ) for this conformation*/
    coords::float_type de_ens;
    /**exp((1 / (k_B*T))*dE_back ) for this conformation*/
    coords::float_type de_ens_back;
    /**temperature*/
    coords::float_type T;
    fepvect(void) :
      e_c_l0{ 0.0 }, e_vdw_l0{ 0.0 }, e_c_l1{ 0.0 }, e_c_l2{ 0.0 }, e_vdw_l1{ 0.0 },
      e_vdw_l2{ 0.0 }, dE{ 0.0 }, dE_back{ 0.0 }, dG{ 0.0 }, dG_back{ 0.0 }, de_ens{ 0.0 }, de_ens_back{ 0.0 }, T{ 0.0 }
    { }
  };


  /** Abstract  base class for interfaces,
  * parent class for all inrterface classes used
  * by CAST for example FF, MOPAC, terachem , gaussian etc.
  * Output should be in kcal/mol!
  */

  class interface_base
  {
  private:
    /**vector of external charges (used in QM/MM methods)*/
    static std::vector<PointCharge> external_charges;

  protected:

    /**pointer to coord object*/
    coords::Coordinates* const coords;
    bool periodic;
    /**is system intact?*/
    bool integrity;
    /**has energy interface its own optimizer?*/
    bool optimizer;
    bool interactions, internal_optimizer;
    /**function that returns external charges*/
    static std::vector<PointCharge> const& get_external_charges() {
      return external_charges;
    }
    /**function to set external charges*/
    static void set_external_charges(std::vector<PointCharge> const& new_ext_charges) {
      external_charges = new_ext_charges;
    }
    /**function that deletes all external charges*/
    static void clear_external_charges() {
      external_charges.clear();
    }
    
  public:

    // some static functions (can be called independant of interface object)

    static std::string create_random_file_name(std::string const& output) {
      std::stringstream ss;
      std::srand(std::time(0));
      ss << (std::size_t(std::rand()) | (std::size_t(std::rand()) << 15));
      return output + "_tmp_" + ss.str();
    }
    /**add an external point charge*/
    static void add_external_charge(PointCharge const& p) {
      external_charges.emplace_back(p);
    }

    /**total energy, in dftbaby interface this is called e_tot*/
    coords::float_type energy;
    /**name of input and output files in some interfaces (GAUSSIAN, PSI4, MOPAC)*/
    std::string id;
    /**total charge of system (needed for QM interfaces)*/
    int charge;
    

    interface_base(coords::Coordinates* coord_pointer) :
      coords(coord_pointer), periodic(false), integrity(true),
      optimizer(false), interactions(false), energy(0.0)
    {
      if (!coord_pointer)
      {
        throw std::runtime_error(
          "Interface without valid coordinates prohibited.");
      }
    }

    interface_base& operator= (interface_base const& other)
    {
      energy = other.energy;
      periodic = other.periodic;
      integrity = other.integrity;
      optimizer = other.optimizer;
      internal_optimizer = other.internal_optimizer;
      interactions = other.interactions;
      return *this;
    }

    virtual void swap(interface_base& other) = 0;

    /** create an copy-instance of derived via new and return pointer*/
    virtual interface_base* clone(coords::Coordinates* coord_object) const = 0;
    // create new instance of derived and move in data
    virtual interface_base* move(coords::Coordinates* coord_object) = 0;

    /** update interface information from coordinates*/
    virtual void update(bool const skip_topology = false) = 0;
    /** delete interface*/
    virtual ~interface_base(void) { }

    /** Energy function*/
    virtual coords::float_type e(void) = 0;

    /** Energy+Gradient function*/
    virtual coords::float_type g(void) = 0;

    /** Energy+Hessian function*/
    virtual coords::float_type h(void) = 0;

    /** Optimization in the intface or interfaced program*/
    virtual coords::float_type o(void) = 0;

    /** Return charges */
    virtual std::vector<coords::float_type> charges() const = 0;
    /**returns the coulomb gradients on external charges (used for QM/MM methods)*/
    virtual coords::Gradients_3D get_g_ext_chg() const = 0;

    /**This function is called in the non-bonding part of energy calculation with periodic boundaries.
    Before calling it the vector between two atoms whose interactions should be calculated is determined
    and given as input to this function.
    The function then modifies the vector: If any of the components (x, y or z) is longer than half the
    box size, the box size is subtracted from this component.
    In this way the resulting vector represents the shortest distance between these two atoms
    in any unit cells next to each other.
    @param r: input vector*/
    void boundary(coords::Cartesian_Point& r);

    /** are periodic boundaries applied? */
    bool has_periodics() const { return periodic; }
    /**does energy interface has its own optimizer*/
    bool has_optimizer() const { return optimizer; }
    bool has_interactions() const { return interactions; }
    /**is system intact*/
    bool intact() const { return integrity; }

    // Output functions
    /**print total energy*/
    virtual void print_E(std::ostream&) const = 0;
    /**print headline for energies*/
    virtual void print_E_head(std::ostream&, bool const endline = true) const = 0;
    /**print partial energies*/
    virtual void print_E_short(std::ostream&, bool const endline = true) const = 0;
    /**print gradients*/
    virtual void print_G_tinkerlike(std::ostream& S, bool const endline = true) const;

    virtual void to_stream(std::ostream&) const = 0;

    coords::Coordinates* cop() const
    {
      return coords;
    }

    //Functions to fetch special outputdata from gaussian interface
    std::vector <double> get_occMO() { return occMO; }
    std::vector <double> get_virtMO() { return virtMO; }
    std::vector <double> get_excitE() { return excitE; }
    std::vector <int> get_state_i() { return state_i; }
    std::vector <int> get_state_j() { return state_j; }
    std::vector <int> get_gz_i_state() { return gz_i_state; }
    coords::Representation_3D  get_ex_ex_trans() { return ex_ex_trans; }
    coords::Representation_3D  get_gz_ex_trans() { return gz_ex_trans; }

    //MO, excitation energies and dipolemoments
    std::vector <double> occMO, virtMO, excitE;
    coords::Representation_3D  ex_ex_trans, gz_ex_trans;
    std::vector <int> state_i, state_j, gz_i_state;

  };

  interface_base* new_interface(coords::Coordinates*);
  interface_base* pre_interface(coords::Coordinates*);
}
