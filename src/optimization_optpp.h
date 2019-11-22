#ifndef OPTIMIZATION_OPTPP
#define OPTIMIZATION_OPTPP

#include "coords.h"

// from OPT++
#include "NLF.h"
#include "OptQNIPS.h"
#include "OptFDNIPS.h"

/**namespace for performing optimizations with OPT++*/
namespace optpp
{
  /**This is the only function to be called from the outside. It performs the optimization.
  @param c: reference to coordinates object. The optimized positions will be filled into this. 
  It returns the energy of the optimized coordinates object.*/
  double perform_optimization(coords::Coordinates& c);
  
  // INTERNAL STUFF 

  // some global variables that are needed in the functions for the OPT++ optimizer
  
  /**pointer to coordinates object that should be optimized*/
  extern std::unique_ptr<coords::Coordinates> coordptr;
  /**dimensions of optimization problem (= 3*number_of_atoms)*/
  extern unsigned int dimension;
  /**initial values for the atomic coordinates that should be optimized
  as CoulmnVector in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  extern NEWMAT::ColumnVector initial_values;

  // some functions for converting stuff from CAST to OPT++ format and back
  
  /**convert coordinates to ColumnVector in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  NEWMAT::ColumnVector convert_coords_to_columnvector(coords::Coordinates const& c);
  /**convert gradients to ColumnVector in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  NEWMAT::ColumnVector convert_gradients_to_columnvector(coords::Gradients_3D const& g);
  /**fill atomic coordinates of the object behind 'coordptr' with the content of a ColumnVector 
  where they are stored in the order <x1, y1, z1, x2, y2, z2, x3, ...>*/
  void fill_columnvector_into_coords(NEWMAT::ColumnVector const& x);

  // functions for performing the optimization
  
  /**preparation function
  sets the global variables that are needed in the functions for the OPT++ optimizer
  and creates a non-linear problem out of these functions (functions see below)*/
  OPTPP::NLF1 prepare(coords::Coordinates& c);
  /**setting options of the OPT++ optimizer to configuration options
  for explanation of options see: https://software.sandia.gov/opt++/opt++2.4_doc/html/ControlParameters.html*/
  void setting_up_optimizer(std::unique_ptr<OPTPP::OptNIPSLike> & opt_ptr);
  /**performs optimization with chosen optimizer
  @param nlf: reference to non-linear problem
              after running this function nlf is at the local minimum*/
  void optimize(OPTPP::NLF1 & nlf);

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
}

#endif 
