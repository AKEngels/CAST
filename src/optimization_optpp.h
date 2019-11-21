#ifndef OPTIMIZATION_OPTPP
#define OPTIMIZATION_OPTPP

#include "coords.h"

// from OPT++
#include "NLF.h"

namespace optpp
{
  extern coords::Coordinates* coordptr;
  extern unsigned int dimension;
  extern NEWMAT::ColumnVector initial_values;

  NEWMAT::ColumnVector convert_coords_to_columnvector(coords::Coordinates const& c);

  // void init_function(int ndim, NEWMAT::ColumnVector& x);
  // void function_to_be_optimized(int mode, int ndim, const NEWMAT::ColumnVector& x, double& fx, NEWMAT::ColumnVector& gx, int& result);
  
  void prepare(coords::Coordinates& c);
  void perform_optimization(coords::Coordinates& c);
}

#endif 
