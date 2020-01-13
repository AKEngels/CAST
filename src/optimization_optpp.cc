#ifdef USE_OPTPP

#include "optimization_optpp.h"
#include "helperfunctions.h"
#include "coords_io.h"

/**some global variables that are needed for init-function and constraints
because of 'static' keyword they can only be used within this file*/
namespace optpp
{
  /**initial values for the atomic coordinates that should be optimized
  as CoulmnVector in the order <x1, y1, z1, x2, y2, z2, x3, ...>
  (needed in init function)*/
  static NEWMAT::ColumnVector initial_values;
  /**vector that contains information for all constrained bonds
  (needed in constraint functions)*/
  static std::vector<optpp::constraint_bond> constraint_bonds;
}

OptppObj::OptppObj(coords::Coordinates& c) : coordobj(c)  // set coordobj
{
  // set dimension
  if (Config::get().coords.fixed.size() == 0) dimension = 3*coordobj.size();
  else dimension = 3* (coordobj.size() - Config::get().coords.fixed.size());
  // set initial values (global variable)
  if (Config::get().coords.fixed.size() == 0) optpp::initial_values = convert_coords_to_columnvector();
  else optpp::initial_values = convert_nonfixed_coords_to_columnvector();
  // set constraints (global variable)
  if (Config::get().optimization.local.optpp_conf.constraints.size() != 0) prepare_constraints();
}

double OptppObj::perform_optimization()
{
  // declaring some variables we need for constraints
  OPTPP::NLF1 nlf_constr, nlf_constr_2;
  OPTPP::NLP *nlp_constr, *nlp_constr_2;
  OPTPP::Constraint constr, constr_2;
  OPTPP::CompoundConstraint comp_constr;
  
  // convert this-pointer to void pointer in order to give it to function
  void* this_as_void = this;

  // set up non-linear problem
  OPTPP::NLF1 nlf;  
  if (optpp::constraint_bonds.size() == 0)  // without constraints
  {
    nlf = OPTPP::NLF1(dimension, optpp::function_to_be_optimized, optpp::init_function, this_as_void);
  }
  else if (optpp::constraint_bonds.size() == 1)  // with one constraint
  {
    nlf_constr = OPTPP::NLF1(dimension,1,optpp::first_constraint_bond_function,optpp::init_function);            
    nlp_constr = new OPTPP::NLP(&nlf_constr);          // 'new' because of conversion from NLF1 to NLP              
    constr = new OPTPP::NonLinearEquation(nlp_constr); // 'new' because of conversion from *NonLinearEquation to Constraint
    comp_constr = OPTPP::CompoundConstraint(constr);                      
    nlf = OPTPP::NLF1(dimension, optpp::function_to_be_optimized, optpp::init_function, &comp_constr, this_as_void);  
  }
  else if (optpp::constraint_bonds.size() == 2)  // with two constraints
  {
    // create first constraint
    nlf_constr = OPTPP::NLF1(dimension,1,optpp::first_constraint_bond_function,optpp::init_function);            
    nlp_constr = new OPTPP::NLP(&nlf_constr);                     
    constr = new OPTPP::NonLinearEquation(nlp_constr); 
    // create second constraint
    nlf_constr_2 = OPTPP::NLF1(dimension,1,optpp::second_constraint_bond_function,optpp::init_function);            
    nlp_constr_2 = new OPTPP::NLP(&nlf_constr_2);                       
    constr_2 = new OPTPP::NonLinearEquation(nlp_constr_2); 
    // create compound constraint and non-linear problem
    comp_constr = OPTPP::CompoundConstraint(constr, constr_2);                      
    nlf = OPTPP::NLF1(dimension, optpp::function_to_be_optimized, optpp::init_function, &comp_constr, this_as_void);  
  }
  else throw std::runtime_error("Too many constraints!");

  // set up the desired optimizer
  std::unique_ptr<OPTPP::OptNIPSLike> optptr;
  setting_up_optimizer(optptr, nlf);
  // perform optimization and print status into file 'OPT_DEFAULT.out'
  optptr->optimize();
  optptr->printStatus((char*)"STATUS");                                                   
  optptr->cleanup();  
  // save coordinates and energy
  fill_columnvector_into_coords(nlf.getXc());       
  return nlf.getF();  
}

NEWMAT::ColumnVector OptppObj::convert_coords_to_columnvector() const
{
  NEWMAT::ColumnVector res(dimension);  // Attention! counting in ColumnVector starts with 1!
  for (auto i = 0u; i < coordobj.size(); ++i)
  {
    auto pos = coordobj.xyz(i);
    res(3*i+1) = pos.x();
    res(3*i+2) = pos.y();
    res(3*i+3) = pos.z();
  }
  return res;
}

NEWMAT::ColumnVector OptppObj::convert_nonfixed_coords_to_columnvector() const
{
  NEWMAT::ColumnVector res(dimension);
  auto counter {0u};
  for (auto i = 0u; i < coordobj.size(); ++i)
  {
    if (!is_in(i, Config::get().coords.fixed))  // only non-fixed
    {
      auto pos = coordobj.xyz(i);
      res(3*counter+1) = pos.x();
      res(3*counter+2) = pos.y();
      res(3*counter+3) = pos.z();
      counter+=1;
    }
  }
  return res;
}

NEWMAT::ColumnVector OptppObj::convert_gradients_to_columnvector(coords::Gradients_3D const& g) const
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

NEWMAT::ColumnVector OptppObj::convert_nonfixed_gradients_to_columnvector(coords::Gradients_3D const& g) const
{
  NEWMAT::ColumnVector res(dimension);
  auto counter {0u};
  for (auto i = 0u; i < g.size(); ++i)
  {
    if (!is_in(i, Config::get().coords.fixed))  // only non-fixed
    {
      auto grad = g[i];
      res(3*counter+1) = grad.x();
      res(3*counter+2) = grad.y();
      res(3*counter+3) = grad.z();
      counter+=1;
    }
  }
  return res;
}

void OptppObj::fill_columnvector_into_coords(NEWMAT::ColumnVector const& x)
{
  coords::Representation_3D xyz_tmp;
  if (Config::get().coords.fixed.size() == 0)   // if there are no fixed atoms
  {
    for (auto i = 0u; i < coordobj.size(); ++i)
    {
      coords::Cartesian_Point xyz(x(3*i+1), x(3*i+2), x(3*i+3)); // Attention! counting in ColumnVector starts with 1!
      xyz_tmp.emplace_back(xyz);
    }
  }
  else                                          // if there are fixed atoms
  {
    auto counter{0u};
    for (auto i = 0u; i < coordobj.size(); ++i)
    {
      if (!is_in(i, Config::get().coords.fixed))   // atom not fixed -> get coordinates from ColumnVector
      {
        coords::Cartesian_Point xyz(x(3*counter+1), x(3*counter+2), x(3*counter+3));
        xyz_tmp.emplace_back(xyz);
        counter+=1;
      }
      else xyz_tmp.emplace_back(coordobj.xyz(i)); // atom fixed -> just take coordinates from before
    }
  }
  coordobj.set_xyz(std::move(xyz_tmp));
}

void OptppObj::prepare_constraints()
{
  optpp::constraint_bonds.clear();
  if (Config::get().coords.fixed.size() == 0)                         // if there are no fixed atoms
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
  else                                                                // if there are fixed atoms
  {
    std::vector<std::size_t> nonfixed_indices;  // indices of atoms that are not fixed
    for (auto i{0u}; i < coordobj.size(); ++i) {
      if (!is_in(i, Config::get().coords.fixed)) {
        nonfixed_indices.emplace_back(i);
      }
    }
    for (auto const&c : Config::get().optimization.local.optpp_conf.constraints)
    {
      optpp::constraint_bond cb;
      auto new_index1 = find_index(c.index1, nonfixed_indices);
      auto new_index2 = find_index(c.index2, nonfixed_indices);
      cb.x1 = 3*new_index1 +1;
      cb.y1 = 3*new_index1 +2;
      cb.z1 = 3*new_index1 +3;
      cb.x2 = 3*new_index2 +1;
      cb.y2 = 3*new_index2 +2;
      cb.z2 = 3*new_index2 +3;
      cb.dist = c.distance;
      optpp::constraint_bonds.emplace_back(cb);
    }
  }
}

void OptppObj::setting_up_optimizer(std::unique_ptr<OPTPP::OptNIPSLike> & optptr, OPTPP::NLF1 & nlf) const
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
  optptr->setLineSearchTol(Config::get().optimization.local.optpp_conf.lsTol);
  optptr->setMaxBacktrackIter(Config::get().optimization.local.optpp_conf.maxBacktrackIter);     
  optptr->setMinStep(Config::get().optimization.local.optpp_conf.minStep);
  optptr->setMaxStep(Config::get().optimization.local.optpp_conf.maxStep);
  optptr->setConTol(Config::get().optimization.local.optpp_conf.conTol);
}

void optpp::init_function(int ndim, NEWMAT::ColumnVector& x)
{
  if (ndim != optpp::initial_values.size()) throw std::runtime_error("Wrong dimension of init function.");
  x = optpp::initial_values;
}

void optpp::function_to_be_optimized(int mode, int ndim, const NEWMAT::ColumnVector& x, double& fx, 
                                     NEWMAT::ColumnVector& gx, int& result, void *vptr)
{
  // convert void-pointer back to pointer to OptppObj, otherwise it can't be dereferenced
  OptppObj* optpp_ptr = static_cast<OptppObj*>(vptr);  

  if (ndim != (int)optpp_ptr->get_dimension()) throw std::runtime_error("Wrong dimension of function to be optimized.");
  
  optpp_ptr->fill_columnvector_into_coords(x);
  double energy = optpp_ptr->get_changeable_coordobj().g();  // gradients are filled into coordobj, so it must be changable

  if (mode & OPTPP::NLPFunction)   // function evaluation
  {
    fx = energy;   
    result = OPTPP::NLPFunction;    
  }
  if (mode & OPTPP::NLPGradient)   // gradient evaluation
  {
    if (Config::get().coords.fixed.size() == 0) gx = optpp_ptr->convert_gradients_to_columnvector(optpp_ptr->get_coordobj().g_xyz());
    else gx = optpp_ptr->convert_nonfixed_gradients_to_columnvector(optpp_ptr->get_coordobj().g_xyz());
    result = OPTPP::NLPGradient;

    if (Config::get().optimization.local.trace)  // write trace (done here because gradients are only calculated once per iteration)
    {
      std::ofstream trace("trace.arc", std::ios_base::app);
      trace << coords::output::formats::tinker(optpp_ptr->get_coordobj()) << std::flush;
    }
  }
}

void optpp::first_constraint_bond_function(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& cx, 
                                          NEWMAT::Matrix& cgx, int& result)
{
  int dimension = optpp::initial_values.size();
  if (ndim != dimension) throw std::runtime_error("Wrong dimension of first constraint function.");
  
  auto cb = optpp::constraint_bonds[0];  // first constraint bond
  constraint_bond_helperfunction(mode, ndim, x, cx, cgx, result, cb);
}

void optpp::second_constraint_bond_function(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& cx, 
                                            NEWMAT::Matrix& cgx, int& result)
{
  int dimension = optpp::initial_values.size();
  if (ndim != dimension) throw std::runtime_error("Wrong dimension of second constraint function.");

  auto cb = optpp::constraint_bonds[1];  // second constraint bond
  constraint_bond_helperfunction(mode, ndim, x, cx, cgx, result, cb);
}

void optpp::constraint_bond_helperfunction(int mode, int ndim, const NEWMAT::ColumnVector& x, NEWMAT::ColumnVector& cx, 
                                            NEWMAT::Matrix& cgx, int& result, optpp::constraint_bond const& cb)
{
  int dimension = optpp::initial_values.size();
  if (ndim != dimension) throw std::runtime_error("Wrong dimension of constraint helperfunction.");

  double x1 = x(cb.x1);
  double y1 = x(cb.y1);
  double z1 = x(cb.z1);
  double x2 = x(cb.x2);
  double y2 = x(cb.y2);
  double z2 = x(cb.z2);

  if (mode & OPTPP::NLPFunction)   // function evaluation
  {
    cx = std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) - cb.dist;
    result = OPTPP::NLPFunction;
  }
  if (mode & OPTPP::NLPGradient)   // gradient evaluation
  {
    for (auto i = 0; i < dimension; ++i) cgx(i+1, 1) = 0.0;  // set all gradients to 0 first
    cgx(cb.x1,1) = (x1 - x2) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    cgx(cb.y1,1) = (y1 - y2) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    cgx(cb.z1,1) = (z1 - z2) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    cgx(cb.x2,1) = (x2 - x1) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    cgx(cb.y2,1) = (y2 - y1) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    cgx(cb.z2,1) = (z2 - z1) / std::sqrt( (x1-x2)*(x1-x2) +  (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    result = OPTPP::NLPGradient;
  }
}

#endif