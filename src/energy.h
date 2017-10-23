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
    fepvar (void) 
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
    fepvect (void) :
      e_c_l1(0.0), e_c_l2(0.0), e_vdw_l1(0.0), e_c_l0(0.0), e_vdw_l0(0.0),
      e_vdw_l2(0.0), dE(0.0), dG(0.0), de_ens(0.0), T(0.0), dE_back(0.0), dG_back(0.0)
    { }
  };
  

  /** Abstract  base class for interfaces, 
  * parent class for all inrterface classes used 
  * by CAST for example FF, MOPAC, terachem , gaussian etc.
  * Output should be in kcal/mol! 
  */

  class interface_base
  {
  protected:

    /**pointer to coord object*/
    coords::Coordinates * const coords;
    bool periodic;
    /**is system intact?*/
    bool integrity;
    /**has energy interface its own optimizer?*/
    bool optimizer;
    bool interactions, internal_optimizer;

  public:

    /**total energy, in dftb interface this is called e_tot*/
    coords::float_type energy;
    coords::Cartesian_Point pb_max, pb_min, pb_dim;

    interface_base (coords::Coordinates *coord_pointer) : 
      coords(coord_pointer), periodic(false), integrity(true), 
      optimizer(false), interactions(false), energy(0.0)
    { 
      if(!coord_pointer) throw std::runtime_error("Interface without valid coordinates prohibited."); 
    }

    interface_base(); 

    interface_base& operator= (interface_base const &other)
    {
      energy = other.energy;
      pb_max = other.pb_max;
      pb_min = other.pb_min;
      pb_dim = other.pb_dim;
      periodic = other.periodic;
      integrity = other.integrity;
      optimizer = other.optimizer;
      internal_optimizer = other.internal_optimizer;
      interactions = other.interactions;
      return *this;
    }

    virtual void swap (interface_base &other) = 0;
    
    /** create an copy-instance of derived via new and return pointer*/
    virtual interface_base * clone (coords::Coordinates * coord_object) const = 0;
    // create new instance of derived and move in data
    virtual interface_base * move (coords::Coordinates * coord_object) = 0;

    /** update interface information from coordinates*/
    virtual void update (bool const skip_topology = false) = 0;
    /** delete interface*/
    virtual ~interface_base (void) { }

    /** Energy function*/
    virtual coords::float_type e (void) = 0;

    /** Energy+Gradient function*/
    virtual coords::float_type g (void) = 0;

    /** Energy+Gradient function*/
    //virtual coords::float_type gi(void) = 0;

    /** Energy+Hessian function*/
    virtual coords::float_type h (void) = 0;

    /** Optimization in the intface or interfaced program*/
    virtual coords::float_type o (void) = 0;

    // Feature getter
    bool has_periodics() const { return periodic; }
    /**does energy interface has its own optimizer*/
    bool has_optimizer() const { return optimizer; }
    bool has_internal_optimizer() const { return optimizer; }
    bool has_interactions() const { return interactions; }
    /**is system intact*/
    bool intact() const { return integrity; }

    // Output functions
    /**print total energy*/
    virtual void print_E (std::ostream&) const = 0;
    /**print headline for energies*/
    virtual void print_E_head (std::ostream&, bool const endline = true) const = 0;
    /**print partial energies*/
    virtual void print_E_short (std::ostream&, bool const endline = true) const = 0;
    /**print gradients*/
    virtual void print_G_tinkerlike (std::ostream&, bool const aggregate = false) const = 0;
    virtual void to_stream (std::ostream&) const = 0;


    //Functions to fetch special outputdata from gaussian interface
    std::vector <double> get_occMO()   {return occMO;}
    std::vector <double> get_virtMO()  {return virtMO;}
    std::vector <double> get_excitE()  {return excitE;}
    std::vector <int> get_state_i()    {return state_i;}
    std::vector <int> get_state_j()    {return state_j;}
    std::vector <int> get_gz_i_state() {return gz_i_state;}
    coords::Representation_3D  get_ex_ex_trans() {return ex_ex_trans;}
    coords::Representation_3D  get_gz_ex_trans() {return gz_ex_trans;}
     
    //MO, excitation energies and dipolemoments
    std::vector <double> occMO, virtMO, excitE;
    coords::Representation_3D  ex_ex_trans, gz_ex_trans;
    std::vector <int> state_i, state_j, gz_i_state;

  };

  interface_base* new_interface (coords::Coordinates *);
  interface_base* pre_interface (coords::Coordinates *);
}
