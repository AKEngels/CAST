#include "optimization_optpp.h"

std::unique_ptr<coords::Coordinates> optpp::coordptr;
unsigned int optpp::dimension;
NEWMAT::ColumnVector optpp::initial_values;
std::vector<optpp::constraint_bond> optpp::constraint_bonds;
OPTPP::CompoundConstraint optpp::compound_constraint;

// void constraint(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& fx, NEWMAT::Matrix& gx, int& result)
// {
//   std::cout<<"constraint\n";
//   if (ndim != optpp::dimension) throw std::runtime_error("Something went wrong");
//   auto cb = optpp::constraint_bonds[0];  // first constraint bond

//   if (mode & OPTPP::NLPFunction)   // function evaluation
//   {
//     fx = std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2)) - cb.dist;
//     result = OPTPP::NLPFunction;
//   }
//   if (mode & OPTPP::NLPGradient)   // gradient evaluation
//   {
//     for (auto i = 0u; i < optpp::dimension; ++i) gx(i+1, 1) = 0.0;  // set all gradients to 0 first
//     gx(cb.x1,1) = (cb.x1 - cb.x2) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
//     gx(cb.y1,1) = (cb.y1 - cb.y2) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
//     gx(cb.z1,1) = (cb.z1 - cb.z2) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
//     gx(cb.x2,1) = (cb.x2 - cb.x1) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
//     gx(cb.y2,1) = (cb.y2 - cb.y1) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
//     gx(cb.z2,1) = (cb.z2 - cb.z1) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
//     result = OPTPP::NLPGradient;
//   }
// }

// double optpp::perform_optimization(coords::Coordinates& c)
// {
//   std::cout<<"start performing optmization\n";
//   auto nlf = optpp::prepare(c);                            // preparation
//   optpp::coordptr = std::make_unique<coords::Coordinates>(c);
//   optpp::dimension = 3*c.size();
//   optpp::initial_values = optpp::convert_coords_to_columnvector(c);
//   std::cout<<"preparation finished\n";
//   optpp::optimize(nlf);                                    // optimization
//   std::cout<<"optmization finished\n";
//   optpp::fill_columnvector_into_coords(nlf.getXc());       // fill optimized values into coordptr
//   c = *coordptr;                                           // set c to optimized coordobj
//   return nlf.getF();                                       // return energy
// }

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

// OPTPP::NLF1 optpp::prepare(coords::Coordinates& c)
// {
//   // set global variables
//   optpp::coordptr = std::make_unique<coords::Coordinates>(c);
//   optpp::dimension = 3*c.size();
//   optpp::initial_values = optpp::convert_coords_to_columnvector(c);
//   // create non-linear problem without constraints
//   if (Config::get().optimization.local.optpp_conf.constraints.size() == 0)
//   {
//     std::cout<<"no constraints\n";
//     OPTPP::NLF1 nlf(dimension, function_to_be_optimized, init_function);
//     return nlf;
//   }
//   // create non-linear problem with constraints
//   else 
//   {
//     // std::cout<<"constraints\n";
//     // //optpp::compound_constraint = optpp::create_constraint();
//     // std::cout<<"created compound constraint\n";
//     // OPTPP::NLF1 nlf(dimension, function_to_be_optimized, init_function, &compound_constraint); 
//     // std::cout<<"created NLF\n";
//     // return nlf;
//   }
// }

// OPTPP::CompoundConstraint optpp::create_constraint()
// {
//   optpp::prepare_constraints();
//   std::cout<<"preparation of constraints finished\n";
//   OPTPP::NLF1 nlf_constr(optpp::dimension,1,first_constraint_bond_function,init_function); 
//   std::cout<<"1\n";            
//   OPTPP::NLP* nlp_constr = new OPTPP::NLP(&nlf_constr);  
//   std::cout<<"2\n";              
//   OPTPP::Constraint constr = new OPTPP::NonLinearEquation(nlp_constr); 
//   std::cout<<"created constraint\n";
//   if (optpp::constraint_bonds.size() == 1) {
//     OPTPP::CompoundConstraint comp_constr(constr);  
//     std::cout<<"return compound constraint\n";
//     return comp_constr;
//   }
//   else std::cout<<"WARNING! Too many constraints. Ignoring all except the first!\n";
// }

void optpp::prepare_constraints()
{
  for (auto const&c : Config::get().optimization.local.optpp_conf.constraints)
  {
    optpp::constraint_bond cb;
    cb.x1 = 3*c.index1 +1;
    cb.y1 = 3*c.index1 +2;
    cb.z1 = 3*c.index1 +3;
    cb.x2 = 3*c.index2 +1;
    cb.y2 = 3*c.index2 +2;
    cb.z2 = 3*c.index2 +3;
    cb.dist = c.distance;
    optpp::constraint_bonds.emplace_back(cb);
  }
}

// void optpp::setting_up_optimizer(std::unique_ptr<OPTPP::OptNIPSLike> const& optptr)  
// { 
//   optptr->setFcnTol(Config::get().optimization.local.optpp_conf.fcnTol);
//   optptr->setGradTol(Config::get().optimization.local.optpp_conf.gradTol); 
//   optptr->setStepTol(Config::get().optimization.local.optpp_conf.stepTol);
//   optptr->setMaxIter(Config::get().optimization.local.optpp_conf.maxIter);              
//   optptr->setMaxFeval(Config::get().optimization.local.optpp_conf.maxFeval);
//   optptr->setMaxBacktrackIter(Config::get().optimization.local.optpp_conf.maxBacktrackIter);     
//   optptr->setMinStep(Config::get().optimization.local.optpp_conf.minStep);
// }

// void optpp::optimize(OPTPP::NLF1 & nlf2)
// {
//   optpp::prepare_constraints();
//   std::cout<<"preparation of constraints finished\n";
//   OPTPP::NLF1 nlf_constr(optpp::dimension,1,constraint,init_function); 
//   std::cout<<"1\n";            
//   OPTPP::NLP* nlp_constr = new OPTPP::NLP(&nlf_constr);  
//   std::cout<<"2\n";              
//   OPTPP::Constraint constr = new OPTPP::NonLinearEquation(nlp_constr); 
//   OPTPP::CompoundConstraint comp_constr(constr); 
//   OPTPP::NLF1 nlf(dimension, function_to_be_optimized, init_function, &compound_constraint);  
//   std::cout<<"in optimization\n";
//   // setting up the desired optimizer
//   std::unique_ptr<OPTPP::OptNIPSLike> optptr;
//   if (Config::get().optimization.local.optpp_conf.optimizer == 0) {
//     optptr = std::make_unique<OPTPP::OptQNIPS>(&nlf);   
//     std::cout<<"QNIPS\n";
//   }
//   else if (Config::get().optimization.local.optpp_conf.optimizer == 1){
//     optptr = std::make_unique<OPTPP::OptFDNIPS>(&nlf); 
//     std::cout<<"FDNIPS\n";
//   } 
//   else throw std::runtime_error("Unknown optimizer chosen for OPT++!");
//   optpp::setting_up_optimizer(optptr);
//   std::cout<<"finshed optimizer setup\n";
//   // perform optimization and print status into file 'OPT_DEFAULT.out'
//   optptr->optimize();
//   optptr->printStatus("STATUS");                                                   
//   optptr->cleanup();                                                   
// }

void optpp::init_function(int ndim, NEWMAT::ColumnVector& x)
{
  if (ndim != optpp::dimension) {
    std::cout<<ndim<<" , "<<optpp::dimension<<"\n";
    throw std::runtime_error("Something went wrong in init-function.");
  }
  x = initial_values;
}

void optpp::function_to_be_optimized(int mode, int ndim, const NEWMAT::ColumnVector& x, double& fx, NEWMAT::ColumnVector& gx, int& result)
{
  std::cout<<"function\n";
  if (ndim != optpp::dimension) {
    std::cout<<ndim<<" , "<<optpp::dimension<<"\n";
    throw std::runtime_error("Something went wrong in function.");
  }
  
  optpp::fill_columnvector_into_coords(x);
  double energy = coordptr->g();

  if (mode & OPTPP::NLPFunction)   // function evaluation
  {
    std::cout<<"energy\n";
    fx = energy;   
    result = OPTPP::NLPFunction;    
  }
  if (mode & OPTPP::NLPGradient)   // gradient evaluation
  {
    std::cout<<"gradient\n";
    gx = optpp::convert_gradients_to_columnvector(coordptr->g_xyz());
    std::cout<<gx.nrows()<<" , "<<gx.ncols()<<"\n";
    result = OPTPP::NLPGradient;
  }
  std::cout<<"leaving\n";
}

void first_constraint_bond_function(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& fx, NEWMAT::Matrix& gx, int& result)
{
  std::cout<<"constraint\n";
  if (ndim != optpp::dimension) {
    std::cout<<ndim<<" , "<<optpp::dimension<<"\n";
    throw std::runtime_error("Something went wrong in constraint-function.");
  }
  auto cb = optpp::constraint_bonds[0];  // first constraint bond

  if (mode & OPTPP::NLPFunction)   // function evaluation
  {
    fx = std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2)) - cb.dist;
    result = OPTPP::NLPFunction;
  }
  if (mode & OPTPP::NLPGradient)   // gradient evaluation
  {
    for (auto i = 0u; i < optpp::dimension; ++i) gx(i+1, 1) = 0.0;  // set all gradients to 0 first
    gx(cb.x1,1) = (cb.x1 - cb.x2) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
    gx(cb.y1,1) = (cb.y1 - cb.y2) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
    gx(cb.z1,1) = (cb.z1 - cb.z2) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
    gx(cb.x2,1) = (cb.x2 - cb.x1) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
    gx(cb.y2,1) = (cb.y2 - cb.y1) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
    gx(cb.z2,1) = (cb.z2 - cb.z1) / std::sqrt( (cb.x1-cb.x2)*(cb.x1-cb.x2) +  (cb.y1-cb.y2)*(cb.y1-cb.y2) + (cb.y1-cb.y2)*(cb.y1-cb.y2));
    result = OPTPP::NLPGradient;
  }
}