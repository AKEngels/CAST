#ifdef USE_OPTPP

#ifndef OPTIMIZATION_OPTPP
#define OPTIMIZATION_OPTPP

#include "coords.h"

// from OPT++
#include "NLF.h"
#include "OptQNIPS.h"
#include "OptFDNIPS.h"
#include "NonLinearEquation.h"
#include "CompoundConstraint.h"

/**namespace for OPT++ optimization that can't be defined inside OptppObj class*/
namespace optpp
{
  /**struct that contains necessary information for a constraint bond function*/
  struct constraint_bond
  {
    /**indices of the constrained coordinates in the ColumnVector for OPT++*/
    std::size_t x1, y1, z1, x2, y2, z2;
    /**distance to which the bond should be constrained*/
    double dist;
  };

  // these are the functions that are put into the OPT++ optimizer
  // for more information about them see: https://software.sandia.gov/opt++/opt++2.4_doc/html/SetUp.html#problem

  /**initializing function for optimization problem
  sets starting values for function to 'inital_values'
  @param ndim: dimension of problem
  @param x: reference to vector that will be filled with inital values*/
  void init_function(int ndim, NEWMAT::ColumnVector& x);
  /**function of optimization problem (more or less energy function of coordinates object in 'coordptr')
  @param mode: integer encoding calculation mode (function value and/or gradient)
  @param ndim: dimension of problem
  @param x: values for which the function should be calculated
  @param fx: will be filled with result of function
  @param gx: will be filled with gradients of function
  @param result: integer encoding which evaluations are available after having run function (value and/or gradients)
  @param vptr: void pointer to object of class OptppObj (is needed mainly to get coordinates object
                  which contains the information to perform an energy/gradient calculation)*/
  void opt_function(int mode, int ndim, const NEWMAT::ColumnVector& x, double& fx, 
                    NEWMAT::ColumnVector& gx, int& result, void *vptr);
  /**constraint function for first bond (i. e. constraint_bonds[0])
  @param mode: integer encoding calculation mode (function value and/or gradient)
  @param ndim: dimension of problem
  @param x: values for which the function should be calculated
  @param fx: will be filled with result of function (ColumnVector because several constraints should be possible)
  @param gx: will be filled with gradients of function (Matrix because several constraints should be possible)
  @param result: integer encoding which evaluations are available after having run function (value and/or gradients)*/
  void first_constraint_bond_function(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& cx, 
                                      NEWMAT::Matrix& cgx, int& result);
  /**constraint function for second bond (i. e. constraint_bonds[1])
  @param mode: integer encoding calculation mode (function value and/or gradient)
  @param ndim: dimension of problem
  @param x: values for which the function should be calculated
  @param fx: will be filled with result of function (ColumnVector because several constraints should be possible)
  @param gx: will be filled with gradients of function (Matrix because several constraints should be possible)
  @param result: integer encoding which evaluations are available after having run function (value and/or gradients)*/
  void second_constraint_bond_function(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& cx, 
                                       NEWMAT::Matrix& cgx, int& result);
  /**helperfunction that is called for first and second constrained bond and performs the stuff with the correct constraint
  @param mode: integer encoding calculation mode (function value and/or gradient)
  @param ndim: dimension of problem
  @param x: values for which the function should be calculated
  @param fx: will be filled with result of function (ColumnVector because several constraints should be possible)
  @param gx: will be filled with gradients of function (Matrix because several constraints should be possible)
  @param result: integer encoding which evaluations are available after having run function (value and/or gradients)
  @param cb: constrained bond that is currently treated*/
  void constraint_bond_helperfunction(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& cx, 
                                       NEWMAT::Matrix& cgx, int& result, optpp::constraint_bond const& cb);
}

/**OPT++ optimization class*/
class OptppObj
{
public:
  /**constructor
  @param c: coordinates object
            a reference of this object is saved here an used for optimization*/
  OptppObj(coords::Coordinates& c);
  /**This is the only function (besides the constructor) to be called from the outside. It performs the optimization.
  @param c: reference to coordinates object. The optimized positions will be filled into this. 
  It returns the energy of the optimized coordinates object.*/
  double perform_optimization();

  // some functions to get class members from the outside 
  // only used in functions of namespace optpp
  
  /**returns dimension*/
  unsigned int get_dimension() const {return dimension;}
  /**returns coordintes object*/
  coords::Coordinates get_coordobj() const {return coordobj;}
  /**also returns coordinates object but in a way so that it can be changed*/
  coords::Coordinates& get_changeable_coordobj() {return coordobj;}
  
  // some functions for converting stuff from CAST to OPT++ format and back

  /**convert gradients to ColumnVector in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  NEWMAT::ColumnVector convert_gradients_to_columnvector(coords::Gradients_3D const& g) const;
  /**convert gradients of atoms that are not fixed to ColumnVector*/
  NEWMAT::ColumnVector convert_nonfixed_gradients_to_columnvector(coords::Gradients_3D const& g) const;
  /**fill atomic coordinates of the object behind 'coordptr' with the content of a ColumnVector 
  if there are fixed atoms the elements of ColumnVector are filled into those coordinates that are not fixed
  where they are stored in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  void fill_columnvector_into_coords(NEWMAT::ColumnVector const& x);

private:
  /**coordinates object that should be optimized
  (as reference because changes on this should be applied on the coordinates 
  that are given to this class)*/
  coords::Coordinates& coordobj;
  /**dimensions of optimization problem (= three times number of movable atoms)*/
  unsigned int dimension;
  
  // functions to convert the atom coordinates of coordobj to OPT++ format
  
  /**convert coordinates to ColumnVector in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  NEWMAT::ColumnVector convert_coords_to_columnvector() const;
  /**convert coordinates of atoms that are not fixed to ColumnVector*/
  NEWMAT::ColumnVector convert_nonfixed_coords_to_columnvector() const;

  // functions for performing the optimization
  
  /**creates the global variable 'constraint_bonds'
  called during constructor if constraints should be set*/
  void prepare_constraints();
  /**choosing OPT++ optimizer and setting it to configuration options
  for explanation of options see: https://software.sandia.gov/opt++/opt++2.4_doc/html/ControlParameters.html*/
  void setting_up_optimizer(std::unique_ptr<OPTPP::OptNIPSLike> & opt_ptr, OPTPP::NLF1 & nlf) const;
};

#endif
#endif 