#include "optimization_optpp.h"

std::unique_ptr<coords::Coordinates> optpp::coordptr;
unsigned int optpp::dimension;
NEWMAT::ColumnVector optpp::initial_values;
std::vector<optpp::constraint_bond> optpp::constraint_bonds;

double optpp::perform_optimization(coords::Coordinates& c)
{
  optpp::prepare(c);            // preparation
  auto E = optpp::optimize();   // optimization
  c = *coordptr;                // set c to optimized coordobj
  return E;                     // return energy
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

void optpp::prepare(coords::Coordinates& c)
{
  optpp::coordptr = std::make_unique<coords::Coordinates>(c);
  optpp::dimension = 3*c.size();
  optpp::initial_values = optpp::convert_coords_to_columnvector(c);
  if (Config::get().optimization.local.optpp_conf.constraints.size() != 0) optpp::prepare_constraints();
}

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

void optpp::setting_up_optimizer(std::unique_ptr<OPTPP::OptNIPSLike> & optptr, OPTPP::NLF1 & nlf)  
{ 
  // choose optimizer
  if (Config::get().optimization.local.optpp_conf.optimizer == 0) optptr = std::make_unique<OPTPP::OptQNIPS>(&nlf);   
  else if (Config::get().optimization.local.optpp_conf.optimizer == 1) optptr = std::make_unique<OPTPP::OptFDNIPS>(&nlf); 
  else throw std::runtime_error("Unknown optimizer chosen for OPT++!");
  // apply settings
  optptr->setFcnTol(Config::get().optimization.local.optpp_conf.fcnTol);
  optptr->setGradTol(Config::get().optimization.local.optpp_conf.gradTol); 
  optptr->setStepTol(Config::get().optimization.local.optpp_conf.stepTol);
  optptr->setMaxIter(Config::get().optimization.local.optpp_conf.maxIter);              
  optptr->setMaxFeval(Config::get().optimization.local.optpp_conf.maxFeval);
  optptr->setMaxBacktrackIter(Config::get().optimization.local.optpp_conf.maxBacktrackIter);     
  optptr->setMinStep(Config::get().optimization.local.optpp_conf.minStep);
}

double optpp::optimize()
{
  // declaring some variables we need (most of them only for constraints)
  OPTPP::NLF1 nlf;  
  OPTPP::NLF1 nlf_constr;
  OPTPP::NLP* nlp_constr;
  OPTPP::Constraint constr;
  OPTPP::CompoundConstraint comp_constr;

  // set up non-linear problem
  if (optpp::constraint_bonds.size() == 0)  // without constraints
  {
    nlf = OPTPP::NLF1(optpp::dimension, optpp::function_to_be_optimized, optpp::init_function);
  }
  else if (optpp::constraint_bonds.size() == 1)  // with one constraint
  {
    nlf_constr = OPTPP::NLF1(optpp::dimension,1,optpp::first_constraint_bond_function,optpp::init_function);            
    nlp_constr = new OPTPP::NLP(&nlf_constr);          // 'new' because of conversion from NLF1 to NLP              
    constr = new OPTPP::NonLinearEquation(nlp_constr); // 'new' because of conversion from *NonLinearEquation to Constraint
    comp_constr = OPTPP::CompoundConstraint(constr);                      
    nlf = OPTPP::NLF1(optpp::dimension, optpp::function_to_be_optimized, optpp::init_function, &comp_constr);  
  }
  else throw std::runtime_error("Too many constraints!");

  // set up the desired optimizer
  std::unique_ptr<OPTPP::OptNIPSLike> optptr;
  optpp::setting_up_optimizer(optptr, nlf);
  // perform optimization and print status into file 'OPT_DEFAULT.out'
  optptr->optimize();
  optptr->printStatus("STATUS");                                                   
  optptr->cleanup();  
  // save coordinates and energy
  optpp::fill_columnvector_into_coords(nlf.getXc());       
  return nlf.getF();                                            
}

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
  if (ndim != optpp::dimension) {
    throw std::runtime_error("Something went wrong in function.");
  }
  
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

void optpp::first_constraint_bond_function(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& fx, NEWMAT::Matrix& gx, int& result)
{
  if (ndim != optpp::dimension) {
    throw std::runtime_error("Something went wrong in constraint-function.");
  }
  auto cb = optpp::constraint_bonds[0];  // first constraint bond
  double x1 = x(cb.x1);
  double y1 = x(cb.y1);
  double z1 = x(cb.z1);
  double x2 = x(cb.x2);
  double y2 = x(cb.y2);
  double z2 = x(cb.z2);

  if (mode & OPTPP::NLPFunction)   // function evaluation
  {
    fx = std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) - cb.dist;
    result = OPTPP::NLPFunction;
  }
  if (mode & OPTPP::NLPGradient)   // gradient evaluation
  {
    for (auto i = 0u; i < optpp::dimension; ++i) gx(i+1, 1) = 0.0;  // set all gradients to 0 first
    gx(cb.x1,1) = (x1 - x2) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    gx(cb.y1,1) = (y1 - y2) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    gx(cb.z1,1) = (z1 - z2) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    gx(cb.x2,1) = (x2 - x1) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    gx(cb.y2,1) = (y2 - y1) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    gx(cb.z2,1) = (z2 - z1) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    result = OPTPP::NLPGradient;
  }
}