#include "optimization_optpp.h"

std::unique_ptr<coords::Coordinates> optpp::coordptr;
unsigned int optpp::dimension;
NEWMAT::ColumnVector optpp::initial_values;

double optpp::perform_optimization(coords::Coordinates& c)
{
  auto nlf = optpp::prepare(c);                      // preparation
  optpp::optimize(nlf);                              // optimization
  optpp::fill_columnvector_into_coords(nlf.getXc()); // fill optimized values into coordptr
  c = *coordptr;                                     // set c to optimized coordobj
  return nlf.getF();                                 // return energy
}

NEWMAT::ColumnVector optpp::convert_coords_to_columnvector(coords::Coordinates const& c)
{
  NEWMAT::ColumnVector res(dimension);  // Attention! counting in ColumnVector starts with 1!
  for (auto i = 0u; i < c.size(); ++i)
  {
    auto pos = c.xyz(i);
    res(3*i+1) = pos.x();
    res(3*i+2) = pos.y();
    res(3*i+3) = pos.z();
  }
  return res;
}

NEWMAT::ColumnVector optpp::convert_gradients_to_columnvector(coords::Gradients_3D const& g)
{
  NEWMAT::ColumnVector res(dimension);  // Attention! counting in ColumnVector starts with 1!
  for (auto i = 0u; i < g.size(); ++i)
  {
    auto grad = g[i];
    res(3*i+1) = grad.x();
    res(3*i+2) = grad.y();
    res(3*i+3) = grad.z();
  }
  return res;
}

void optpp::fill_columnvector_into_coords(NEWMAT::ColumnVector const& x)
{
  coords::Representation_3D xyz_tmp;
  for (auto i = 0u; i < coordptr->size(); ++i)
  {
    coords::Cartesian_Point xyz(x(3*i+1), x(3*i+2), x(3*i+3)); // Attention! counting in ColumnVector starts with 1!
    xyz_tmp.push_back(xyz);
  }
  coordptr->set_xyz(std::move(xyz_tmp));
}

OPTPP::NLF1 optpp::prepare(coords::Coordinates& c)
{
  // set global variables
  optpp::coordptr = std::make_unique<coords::Coordinates>(c);
  optpp::dimension = 3*c.size();
  optpp::initial_values = optpp::convert_coords_to_columnvector(c);
  // create non-linear problem
  OPTPP::NLF1 nlf(dimension, function_to_be_optimized, init_function);
  // create constraints

  return nlf;
}

void optpp::setting_up_optimizer(std::unique_ptr<OPTPP::OptNIPSLike> const& optptr)  
{ 
  optptr->setFcnTol(Config::get().optimization.local.optpp_conf.fcnTol);
  optptr->setGradTol(Config::get().optimization.local.optpp_conf.gradTol); 
  optptr->setStepTol(Config::get().optimization.local.optpp_conf.stepTol);
  optptr->setMaxIter(Config::get().optimization.local.optpp_conf.maxIter);              
  optptr->setMaxFeval(Config::get().optimization.local.optpp_conf.maxFeval);
  optptr->setMaxBacktrackIter(Config::get().optimization.local.optpp_conf.maxBacktrackIter);     
  optptr->setMinStep(Config::get().optimization.local.optpp_conf.minStep);
}

void optpp::optimize(OPTPP::NLF1 & nlf)
{
  // setting up the desired optimizer
  std::unique_ptr<OPTPP::OptNIPSLike> optptr;
  if (Config::get().optimization.local.optpp_conf.optimizer == 0) optptr = std::make_unique<OPTPP::OptQNIPS>(&nlf);   
  else if (Config::get().optimization.local.optpp_conf.optimizer == 1) optptr = std::make_unique<OPTPP::OptFDNIPS>(&nlf); 
  else throw std::runtime_error("Unknown optimizer chosen for OPT++!");
  optpp::setting_up_optimizer(optptr);
  // perform optimization and print status into file 'OPT_DEFAULT.out'
  optptr->optimize();
  optptr->printStatus("STATUS");                                                   
  optptr->cleanup();                                                   
}

void optpp::init_function(int ndim, NEWMAT::ColumnVector& x)
{
  if (ndim != optpp::dimension) throw std::runtime_error("Something went wrong");
  x = initial_values;
}

void optpp::function_to_be_optimized(int mode, int ndim, const NEWMAT::ColumnVector& x, double& fx, NEWMAT::ColumnVector& gx, int& result)
{
  if (ndim != optpp::dimension) throw std::runtime_error("Something went wrong");
  
  optpp::fill_columnvector_into_coords(x);
  double energy = coordptr->g();

  if (mode & OPTPP::NLPFunction)   // function evaluation
  {
    fx = energy;   
    result = OPTPP::NLPFunction;    
  }
  if (mode & OPTPP::NLPGradient)   // gradient evaluation
  {
    gx = optpp::convert_gradients_to_columnvector(coordptr->g_xyz());
    result = OPTPP::NLPGradient;
  }
}