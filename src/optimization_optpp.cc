#include "optimization_optpp.h"

coords::Coordinates* optpp::coordptr;
unsigned int optpp::dimension;
NEWMAT::ColumnVector optpp::initial_values;

void optpp::perform_optimization(coords::Coordinates& c)
{
  optpp::prepare(c);
}

NEWMAT::ColumnVector optpp::convert_coords_to_columnvector(coords::Coordinates const& c)
{
  NEWMAT::ColumnVector res(dimension);
  for (auto i = 0u; i < c.size(); ++i)
  {
    auto pos = c.xyz(i);
    res(3*i+1) = pos.x();
    res(3*i+2) = pos.y();
    res(3*i+3) = pos.z();
  }
  return res;
}

void optpp::prepare(coords::Coordinates& c)
{
  optpp::coordptr = &c;
  optpp::dimension = 3*c.size();
  optpp::initial_values = optpp::convert_coords_to_columnvector(c);
  std::cout<<"Initial values:\n";
  for (auto i=0u; i<optpp::dimension; ++i) std::cout<<optpp::initial_values(i+1)<<"  ";
  std::cout<<"\n";
  // OPTPP::NLF1 nlf(dimension, function_to_be_optimized, init_function);
  // std::cout<<"created nlf\n";
}

// void optpp::init_function(int ndim, NEWMAT::ColumnVector& x)
// {
//   if (ndim != 3*optpp::coordobj.size()) throw std::runtime_error("Something went wrong");
//   x = initial_values;
// }

// void optpp::function_to_be_optimized(int mode, int ndim, const NEWMAT::ColumnVector& x, double& fx, NEWMAT::ColumnVector& gx, int& result)
// {
//   if (ndim != 3*optpp::coordobj.size()) throw std::runtime_error("Something went wrong");

//   if (mode & OPTPP::NLPFunction)   
//   {
//     fx = 5;   
//     result = OPTPP::NLPFunction;    
//   }
//   if (mode & OPTPP::NLPGradient)   
//   {
//     for (auto i=0u; i < dimension; ++i) gx(i+1) = 5;
//     result = OPTPP::NLPGradient;
//   }
// }