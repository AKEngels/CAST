#ifdef USE_PYTHON
#include "energy_int_dftb.h"


/*
dftb sysCall functions
*/

std::string energy::interfaces::dftb::get_python_modulepath(std::string modulename)
{
    std::string find = "import "+modulename+"\nwith open('tmpfile.txt','w') as fn:\n    fn.write("+modulename+".__file__)";
    const char *c_find = find.c_str();
    PyRun_SimpleString(c_find);  //call a python programme to find the modulepath and write it to tmpfile
    std::ifstream file("tmpfile.txt");  //open tmpfile and read content
    std::string content;
    file >> content;
    remove("tmpfile.txt");  //delete tmpfile
    return content.substr(0,content.size()-14-modulename.size());  //give back path without filename __init__.pyc and modulename
}

void energy::interfaces::dftb::create_dftbaby_configfile()
{
    std::ofstream file("dftbaby.cfg");
    file << "[DFTBaby]\n\n";
    file << "gradient_file = "+Config::get().energy.dftb.gradfile+"\n\n"; // has to be set in any case cause otherwise gradients are not calculated
    file << "gradient_state = "+std::to_string(Config::get().energy.dftb.gradstate)+"\n\n";  // set because standard not ground state
    file << "verbose = "+std::to_string(Config::get().energy.dftb.verbose)+"\n\n";
    if (Config::get().energy.dftb.longrange == false)
    {    //in dftbaby long range correction is standard
      file << "long_range_correction = 0\n\n";
    }
    if (Config::get().energy.dftb.cutoff != 0)
        file << "distance_cutoff = "+std::to_string(Config::get().energy.dftb.cutoff)+"\n\n";
    if (Config::get().energy.dftb.lr_dist != 0)
        file << "long_range_radius = "+std::to_string(Config::get().energy.dftb.lr_dist)+"\n\n";
    if (Config::get().energy.dftb.maxiter != 0)
        file << "maxiter = "+std::to_string(Config::get().energy.dftb.maxiter)+"\n\n";
    if (Config::get().energy.dftb.conv_threshold != "0")
        file << "scf_conv = "+Config::get().energy.dftb.conv_threshold+"\n\n";
    if (Config::get().energy.dftb.states != 0)
        file << "nstates = "+std::to_string(Config::get().energy.dftb.states)+"\n\n";
    if (Config::get().energy.dftb.orb_occ != 0)
        file << "nr_active_occ = "+std::to_string(Config::get().energy.dftb.orb_occ)+"\n\n";
    if (Config::get().energy.dftb.orb_virt != 0)
        file << "nr_active_virt = "+std::to_string(Config::get().energy.dftb.orb_virt)+"\n\n";
    if (Config::get().energy.dftb.diag_conv != "0")
        file << "diag_conv = "+Config::get().energy.dftb.diag_conv+"\n\n";
    if (Config::get().energy.dftb.diag_maxiter != 0)
        file << "diag_maxiter = "+std::to_string(Config::get().energy.dftb.diag_maxiter)+"\n\n";
    file.close();
}

std::string energy::interfaces::dftb::create_pythonpath(std::string numpath, std::string scipath)
{
    std::string pythonpaths_str = Py_GetPath();
    std::string path;
    std::vector<std::string> pythonpaths = split(pythonpaths_str,':');
    path = "import sys\n";
    for (auto p : pythonpaths)
    {
      path += "sys.path.append('"+p+"')\n";
    }
    path += "sys.path.append('"+Config::get().energy.dftb.path+"')\n";
    path += "sys.path.append('"+numpath+"')\n";
    path += "sys.path.append('"+scipath+"')\n";
    return path;
}

energy::interfaces::dftb::sysCallInterface::sysCallInterface(coords::Coordinates * cp) :
  energy::interface_base(cp),
  e_bs(0.0), e_coul(0.0), e_rep(0.0), e_tot(0.0)
{
    std::string numpath = get_python_modulepath("numpy");
    std::string scipath = get_python_modulepath("scipy");
    add_path = create_pythonpath(numpath, scipath);
    create_dftbaby_configfile();
    optimizer = Config::get().energy.dftb.opt;
}

energy::interfaces::dftb::sysCallInterface::sysCallInterface(sysCallInterface const & rhs, coords::Coordinates *cobj) :
  interface_base(cobj),
  e_bs(rhs.e_bs), e_coul(rhs.e_coul), e_rep(rhs.e_rep), e_tot(rhs.e_tot)
{
  interface_base::operator=(rhs);
}

energy::interface_base * energy::interfaces::dftb::sysCallInterface::clone(coords::Coordinates * coord_object) const
{
  sysCallInterface * tmp = new sysCallInterface(*this, coord_object);
  return tmp;
}

energy::interface_base * energy::interfaces::dftb::sysCallInterface::move(coords::Coordinates * coord_object)
{
  sysCallInterface * tmp = new sysCallInterface(*this, coord_object);
  return tmp;
}

void energy::interfaces::dftb::sysCallInterface::swap(interface_base &rhs)
{
  swap(dynamic_cast<sysCallInterface&>(rhs));
}

void energy::interfaces::dftb::sysCallInterface::swap(sysCallInterface &rhs)
{
  interface_base::swap(rhs);
}

energy::interfaces::dftb::sysCallInterface::~sysCallInterface(void)
{

}

/*
Energy class functions that need to be overloaded
*/

// Energy function
double energy::interfaces::dftb::sysCallInterface::e(void)
{
  integrity = true;
  grad_var = false;

    //write inputstructure
    std::ofstream file("tmp_struc.xyz");
    file << coords::output::formats::xyz_dftb(*this->coords);
    file.close();
    
    //call programme
    std::string result_str; 
    PyObject *modul, *funk, *prm, *ret;
    
    PySys_SetPath("."); //set path
    const char *c = add_path.c_str();  //add paths from variable add_path
    PyRun_SimpleString(c);

    modul = PyImport_ImportModule("dftbaby_interface"); //import module 

    if(modul) 
        { 
        funk = PyObject_GetAttrString(modul, "calc_energies"); //create function
        prm = Py_BuildValue("(ss)", "tmp_struc.xyz", "dftbaby.cfg"); //give parameters
        ret = PyObject_CallObject(funk, prm);  //call function with parameters

        result_str = PyString_AsString(ret); //read function return (has to be a string)
        result_str = result_str.substr(1,result_str.size()-2);  //process return
        std::vector<std::string> result_vec = split(result_str, ',');

        //read energies and convert them to kcal/mol
        e_bs = std::stod(result_vec[0])*627.503; 
        e_coul = std::stod(result_vec[1])*627.503;
        e_rep = std::stod(result_vec[3])*627.503;
        e_tot = std::stod(result_vec[4])*627.503;
        if (result_vec.size() == 6) e_lr = std::stod(result_vec[5])*627.503;
        
        //delete PyObjects
        Py_DECREF(prm); 
        Py_DECREF(ret); 
        Py_DECREF(funk); 
        Py_DECREF(modul); 
        } 
    else 
    {
        printf("ERROR: module dftbaby_interface not found\n"); 
        std::exit(0);
    }
    std::remove("tmp_struc.xyz"); // delete file
  return e_tot;
}

// Energy+Gradient function
double energy::interfaces::dftb::sysCallInterface::g(void)
{
  integrity = true;
  grad_var = true;

  // write inputstructure
  std::ofstream file("tmp_struc.xyz");
  file << coords::output::formats::xyz_dftb(*this->coords);
  file.close();
    
  //call programme
  std::string result_str; 
  PyObject *modul, *funk, *prm, *ret;
    
  PySys_SetPath("."); //set path
  const char *c = add_path.c_str();  //add paths from variable add_path
  PyRun_SimpleString(c);

  modul = PyImport_ImportModule("dftbaby_interface"); //import module 

  if(modul) 
    { 
      funk = PyObject_GetAttrString(modul, "calc_gradients"); //create function
      prm = Py_BuildValue("(ss)", "tmp_struc.xyz", "dftbaby.cfg"); //give parameters
      ret = PyObject_CallObject(funk, prm);  //call function with parameters

      result_str = PyString_AsString(ret); //read function return (has to be a string)
      result_str = result_str.substr(1,result_str.size()-2);  //process return
      std::vector<std::string> result_vec = split(result_str, ',');

      //read energies and convert them to kcal/mol
        e_bs = std::stod(result_vec[0])*627.503; 
        e_coul = std::stod(result_vec[1])*627.503;
        e_rep = std::stod(result_vec[3])*627.503;
        e_tot = std::stod(result_vec[4])*627.503;
        if (result_vec.size() == 6) e_lr = std::stod(result_vec[5])*627.503;
        
      //delete PyObjects
      Py_DECREF(prm); 
      Py_DECREF(ret); 
      Py_DECREF(funk); 
      Py_DECREF(modul); 
    } 
    else 
    {
      printf("ERROR: module dftbaby_interface not found\n"); 
      std::exit(0);
    }
    
    double CONVERSION_FACTOR = 627.503 / 0.5291172107;  // hartree/bohr -> kcal/(mol*A)
    //read gradients
    std::string line;
    coords::Representation_3D g_tmp;
    std::ifstream infile(Config::get().energy.dftb.gradfile);
    std::getline(infile, line);  //discard fist two lines
    std::getline(infile, line);
    std::string element;
    double x,y,z;
    while (infile >> element >> x >> y >> z)  //read gradients and convert them to kcal/mol
    {
        coords::Cartesian_Point g(x*CONVERSION_FACTOR,y*CONVERSION_FACTOR,z*CONVERSION_FACTOR);
        g_tmp.push_back(g);
    }
    infile.close();
    const char *gradfile = Config::get().energy.dftb.gradfile.c_str();
    std::remove(gradfile); // delete file
    std::remove("tmp_struc.xyz"); // delete file
    coords->swap_g_xyz(g_tmp); //give gradients to coordobject

  return e_tot;
}

// Hessian function
double energy::interfaces::dftb::sysCallInterface::h(void)
{
  integrity = true;
  grad_var = false;

  std::cout<<"calc hessian\n";
      //write inputstructure
      std::ofstream file("tmp_struc.xyz");
      file << coords::output::formats::xyz_dftb(*this->coords);
      file.close();
      
      //call programme
      std::string result_str; 
      PyObject *modul, *funk, *prm, *ret;
      
      PySys_SetPath("."); //set path
      const char *c = add_path.c_str();  //add paths from variable add_path
      PyRun_SimpleString(c);
  
      modul = PyImport_ImportModule("dftbaby_interface"); //import module 
  
      if(modul) 
          { 
          funk = PyObject_GetAttrString(modul, "hessian"); //create function
          prm = Py_BuildValue("(ss)", "tmp_struc.xyz", "dftbaby.cfg"); //give parameters
          ret = PyObject_CallObject(funk, prm);  //call function with parameters
  
          result_str = PyString_AsString(ret); //read function return (has to be a string)
          result_str = result_str.substr(1,result_str.size()-2);  //process return
          std::vector<std::string> result_vec = split(result_str, ',');
  
          //read energies and convert them to kcal/mol
          e_bs = std::stod(result_vec[0])*627.503; 
          e_coul = std::stod(result_vec[1])*627.503;
          e_rep = std::stod(result_vec[3])*627.503;
          e_tot = std::stod(result_vec[4])*627.503;
          if (result_vec.size() == 6) e_lr = std::stod(result_vec[5])*627.503;
          
          //delete PyObjects
          Py_DECREF(prm); 
          Py_DECREF(ret); 
          Py_DECREF(funk); 
          Py_DECREF(modul); 
          } 
      else 
      {
          printf("ERROR: module dftbaby_interface not found\n"); 
          std::exit(0);
      }
      std::remove("tmp_struc.xyz"); // delete file

  return energy;
}

// Optimization
double energy::interfaces::dftb::sysCallInterface::o(void)
{
      //write inputstructure
    std::ofstream file("tmp_struc.xyz");
    file << coords::output::formats::xyz_dftb(*this->coords);
    file.close();
    
    //call programme
    std::string result_str; 
    PyObject *modul, *funk, *prm, *ret;
    
    PySys_SetPath("."); //set path
    const char *c = add_path.c_str();  //add paths from variable add_path
    PyRun_SimpleString(c);

    modul = PyImport_ImportModule("dftbaby_interface"); //import module 

    if(modul) 
        { 
        funk = PyObject_GetAttrString(modul, "opt"); //create function
        prm = Py_BuildValue("(ss)", "tmp_struc.xyz", "dftbaby.cfg"); //give parameters
        ret = PyObject_CallObject(funk, prm);  //call function with parameters

        result_str = PyString_AsString(ret); //read function return (has to be a string)
        result_str = result_str.substr(1,result_str.size()-2);  //process return
        std::vector<std::string> result_vec = split(result_str, ',');

        //read energies and convert them to kcal/mol
        e_bs = std::stod(result_vec[0])*627.503; 
        e_coul = std::stod(result_vec[1])*627.503;
        e_rep = std::stod(result_vec[3])*627.503;
        e_tot = std::stod(result_vec[4])*627.503;
        if (result_vec.size() == 6) e_lr = std::stod(result_vec[5])*627.503;
        
        //delete PyObjects
        Py_DECREF(prm); 
        Py_DECREF(ret); 
        Py_DECREF(funk); 
        Py_DECREF(modul); 
        } 
    else 
    {
        printf("ERROR: module dftbaby_interface not found\n"); 
        std::exit(0);
    }
    
    //read new geometry
    std::string line;
    coords::Representation_3D xyz_tmp;
    std::ifstream infile("tmp_struc_opt.xyz");
    std::getline(infile, line);  //discard fist two lines
    std::getline(infile, line);
    std::string element;
    double x,y,z;
    while (infile >> element >> x >> y >> z)  //read gradients and convert them to kcal/mol
    {
        coords::Cartesian_Point xyz(x*627.503,y*627.503,z*627.503);
        xyz_tmp.push_back(xyz);
    }
    infile.close();
    coords->set_xyz(std::move(xyz_tmp));

    std::remove("tmp_struc_opt.xyz"); // delete file
    std::remove("tmp_struc.xyz"); // delete file
  return e_tot;
}

// Output functions
void energy::interfaces::dftb::sysCallInterface::print_E(std::ostream &S) const
{
  S << "Total Energy:      ";
  S << std::right << std::setw(16) << std::fixed << std::setprecision(8) << e_tot;
}

void energy::interfaces::dftb::sysCallInterface::print_E_head(std::ostream &S, bool const endline) const
{
  S << "Energies\n";
  S << std::right << std::setw(24) << "E_bs";
  S << std::right << std::setw(24) << "E_coul";
  S << std::right << std::setw(24) << "E_lr";
  S << std::right << std::setw(24) << "E_rep";
  S << std::right << std::setw(24) << "SUM\n\n";
}

void energy::interfaces::dftb::sysCallInterface::print_E_short(std::ostream &S, bool const endline) const
{
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << e_bs;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << e_coul;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << e_lr;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << e_rep;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << e_tot << '\n';
  S << "\n";
}

void energy::interfaces::dftb::sysCallInterface::print_G_tinkerlike(std::ostream &S, bool const) const 
{ 
  S << " Cartesian Gradient Breakdown over Individual Atoms :" << std::endl << std::endl;
  S << "  Type      Atom              dE/dX       dE/dY       dE/dZ          Norm" << std::endl << std::endl;
  for(std::size_t k=0; k < coords->size(); ++k)
  {
    S << " Anlyt";
    S << std::right << std::setw(10) << k+1U;
    S << "       ";
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).x();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).y();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).z();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4);
    S << std::sqrt(
      coords->g_xyz(k).x() * coords->g_xyz(k).x()
    + coords->g_xyz(k).y() * coords->g_xyz(k).y()
    + coords->g_xyz(k).z() * coords->g_xyz(k).z()) << std::endl;
  }
}

void energy::interfaces::dftb::sysCallInterface::to_stream(std::ostream&) const { }

bool energy::interfaces::dftb::sysCallInterface::check_bond_preservation(void) const
{
  std::size_t const N(coords->size());
  for (std::size_t i(0U); i < N; ++i)
  { // cycle over all atoms i
    if (!coords->atoms(i).bonds().empty())
    {
      std::size_t const M(coords->atoms(i).bonds().size());
      for (std::size_t j(0U); j < M && coords->atoms(i).bonds(j) < i; ++j)
      { // cycle over all atoms bound to i
        double const L(geometric_length(coords->xyz(i) - coords->xyz(coords->atoms(i).bonds(j))));
        if (L > 2.2) return false;
      }
    }
  }
  return true;
}
#endif