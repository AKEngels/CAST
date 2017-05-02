#pragma once 

#if defined _OPENMP
  #include <omp.h>
#endif

#include <iostream>
#include <memory>
#include <string>
//#include <fftw3.h>
#include "coords_rep.h"
//#include "Array.h"
//#include "Complex.h"

namespace energy
{
 /* struct pmevar
  {
    struct feinfo{
      int atom, flag;
    };
    fftw_complex *in;
    fftw_plan forward;
    fftw_plan backward;
    std::vector<int> fepi;
    std::vector<int> fepo;
    std::vector<int> fepa;
    std::vector<feinfo> feinf;
    coords::float_type ewaldcoeff;
    std::vector<double>atomcharges;
    std::vector<int> gridnumbers;
    std::vector<double> moduli1;
    std::vector<double> moduli2;
    std::vector<double> moduli3;
    Array::array3<Complex> charges;
    Array::array3<double> fractionalcharges;
    Array::array3<double> bscx;
    Array::array3<double> bscy;
    Array::array3<double> bscz;
    Array::array2<double> recivectors;
    Array::array2<double> initgrid;
    Array::array2<double> initingrid;
    Array::array2<double> initoutgrid;
    Array::array2<double> initallgrid;
    Array::array2<int>parallelpme;
    Array::array2<int>parallelinpme;
    Array::array2<int>paralleloutpme;
    Array::array2<int>parallelallpme;
    int natoms;
    int rgrid1, rgrid2, rgrid3, rgridtotal, nrough1, nrough2, nrough3, roughoffset;
    int roughleft, roughright, bsoffset;
    int forwardplan, backwardplan;
    int dim1, dim2, dim3;
    int nxpoints, nypoints, nzpoints;
    std::string wisdom;
    int fftgridmax;
    int bsplineorder;
    coords::float_type treshold, elecfac;
    pmevar(void)
      : ewaldcoeff(1.0), dim1(0), dim2(0), dim3(0),
      nxpoints(0), nypoints(0), nzpoints(0),
      fftgridmax(200), bsplineorder(5), treshold(1e-8)
    {}
  };*/
  struct fepvar
  {
	  coords::float_type ein, eout, vin, vout;
	  coords::float_type dein, deout, dvin, dvout;
    int step;
    fepvar (void) 
      : ein(1.0), eout(1.0), vin(1.0), vout(1.0),
      dein(1.0), deout(1.0), dvin(1.0), dvout(1.0),
      step(0)
    { }
  };

  struct fepvect
  {
	  coords::float_type e_c_l1, e_c_l2;
	  coords::float_type e_vdw_l1, e_vdw_l2;
	  coords::float_type dE, dG, de_ens;
    coords::float_type T;
    fepvect (void) :
      e_c_l1(0.0), e_c_l2(0.0), e_vdw_l1(0.0), 
      e_vdw_l2(0.0), dE(0.0), dG(0.0), de_ens(0.0), T(0.0)
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

    coords::Coordinates * const coords;
    bool periodic, integrity, optimizer, interactions, internal_optimizer;

  public:

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

    /** Energy+Gradient+Hessian function*/
    virtual coords::float_type h (void) = 0;

    /** Optimization in the intface or interfaced program*/
    virtual coords::float_type o (void) = 0;

    // Feature getter
    bool has_periodics() const { return periodic; }
    bool has_optimizer() const { return optimizer; }
    bool has_internal_optimizer() const { return optimizer; }
    bool has_interactions() const { return interactions; }

    bool intact() const { return integrity; }

    // Output functions
    virtual void print_E (std::ostream&) const = 0;
    virtual void print_E_head (std::ostream&, bool const endline = true) const = 0;
    virtual void print_E_short (std::ostream&, bool const endline = true) const = 0;
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
    


  private:
    //MO, excitation energies and dipolemoments
    std::vector <double> occMO, virtMO, excitE;
    coords::Representation_3D  ex_ex_trans, gz_ex_trans;
    std::vector <int> state_i, state_j, gz_i_state;

  };

  interface_base* new_interface (coords::Coordinates *);
  interface_base* pre_interface (coords::Coordinates *);


}
