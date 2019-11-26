#ifndef OPTIMIZATION_OPTPP
#define OPTIMIZATION_OPTPP

#include "coords.h"

// from OPT++
#include "NLF.h"
#include "OptQNIPS.h"
#include "OptFDNIPS.h"
#include "NonLinearEquation.h"
#include "CompoundConstraint.h"

/**namespace for performing optimizations with OPT++*/
namespace optpp
{
  /**This is the only function to be called from the outside. It performs the optimization.
  @param c: reference to coordinates object. The optimized positions will be filled into this. 
  It returns the energy of the optimized coordinates object.*/
  double perform_optimization(coords::Coordinates& c);
  
  // INTERNAL STUFF 
  
  /**struct that contains necessary information for a constraint bond function*/
  struct constraint_bond
  {
    /**indices of the constrained coordinates in the ColumnVector for OPT++*/
    std::size_t x1, y1, z1, x2, y2, z2;
    /**distance to which the bond should be constrained*/
    double dist;
  };

  // some global variables that are needed in the functions for the OPT++ optimizer
  
  /**pointer to coordinates object that should be optimized*/
  extern std::unique_ptr<coords::Coordinates> coordptr;
  /**dimensions of optimization problem (= 3*number_of_atoms)*/
  extern unsigned int dimension;
  /**initial values for the atomic coordinates that should be optimized
  as CoulmnVector in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  extern NEWMAT::ColumnVector initial_values;
  /**vector that contains information for all constrained bonds*/
  extern std::vector<optpp::constraint_bond> constraint_bonds;

  // some functions for converting stuff from CAST to OPT++ format and back
  
  /**convert coordinates to ColumnVector in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  NEWMAT::ColumnVector convert_coords_to_columnvector(coords::Coordinates const& c);
  /**convert coordinates of atoms that are not fixed to ColumnVector*/
  NEWMAT::ColumnVector convert_nonfixed_coords_to_columnvector(coords::Coordinates const& c);
  /**convert gradients to ColumnVector in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  NEWMAT::ColumnVector convert_gradients_to_columnvector(coords::Gradients_3D const& g);
  /**convert gradients of atoms that are not fixed to ColumnVector*/
  NEWMAT::ColumnVector convert_nonfixed_gradients_to_columnvector(coords::Gradients_3D const& g);
  /**fill atomic coordinates of the object behind 'coordptr' with the content of a ColumnVector 
  if there are fixed atoms the elements of ColumnVector are filled into those coordinates that are not fixed
  where they are stored in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  void fill_columnvector_into_coords(NEWMAT::ColumnVector const& x);

  // functions for performing the optimization
  
  /**preparation function
  sets the global variables that are needed in the functions for the OPT++ optimizer*/
  void prepare(coords::Coordinates& c);
  /**creates the global variable 'constraint_bonds'
  called during prepare() if constraints should be set*/
  void prepare_constraints();
  /**choosing OPT++ optimizer and setting it to configuration options
  for explanation of options see: https://software.sandia.gov/opt++/opt++2.4_doc/html/ControlParameters.html*/
  void setting_up_optimizer(std::unique_ptr<OPTPP::OptNIPSLike> & opt_ptr, OPTPP::NLF1 & nlf);
  /**performs optimization, fills coords into 'coordptr' and returns energy of optimized structure*/
  double optimize();

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
  @param result: integer encoding which evaluations are available after having run function (value and/or gradients)*/
  void function_to_be_optimized(int mode, int ndim, const NEWMAT::ColumnVector& x, double& fx, NEWMAT::ColumnVector& gx, int& result);
  /**constraint function for first bond (i. e. constraint_bonds[0])
  @param mode: integer encoding calculation mode (function value and/or gradient)
  @param ndim: dimension of problem
  @param x: values for which the function should be calculated
  @param fx: will be filled with result of function (for some reason this must be a ColumnVector)
  @param gx: will be filled with gradients of function (for some reason this must be a Matrix)
  @param result: integer encoding which evaluations are available after having run function (value and/or gradients)*/
  void first_constraint_bond_function(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& fx, NEWMAT::Matrix& gx, int& result);
  /**constraint function for second bond (i. e. constraint_bonds[1])
  @param mode: integer encoding calculation mode (function value and/or gradient)
  @param ndim: dimension of problem
  @param x: values for which the function should be calculated
  @param fx: will be filled with result of function (for some reason this must be a ColumnVector)
  @param gx: will be filled with gradients of function (for some reason this must be a Matrix)
  @param result: integer encoding which evaluations are available after having run function (value and/or gradients)*/
  void second_constraint_bond_function(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& fx, NEWMAT::Matrix& gx, int& result);
}

#endif 
