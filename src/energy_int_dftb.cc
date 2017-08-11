
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

energy::interfaces::dftb::sysCallInterface::sysCallInterface(coords::Coordinates * cp) :
  energy::interface_base(cp),
  e_bs(0.0), e_coul(0.0), e_rep(0.0), e_tot(0.0), id(Config::get().general.outputFilename), failcounter(0u)
{
    //find paths to numpy and scipy
    std::string numpath = get_python_modulepath("numpy");
    std::string scipath = get_python_modulepath("scipy");

    //create pythonpath
    std::string pythonpaths_str = Py_GetPath();
    std::vector<std::string> pythonpaths = split(pythonpaths_str,':');
    add_path = "import sys\n";
    for (auto p : pythonpaths)
    {
      add_path += "sys.path.append('"+p+"')\n";
    }
    add_path += "sys.path.append('"+Config::get().energy.dftb.path+"')\n";
    add_path += "sys.path.append('"+numpath+"')\n";
    add_path += "sys.path.append('"+scipath+"')\n";
}

energy::interfaces::dftb::sysCallInterface::sysCallInterface(sysCallInterface const & rhs, coords::Coordinates *cobj) :
  interface_base(cobj),
  e_bs(rhs.e_bs), e_coul(rhs.e_coul), e_rep(rhs.e_rep), e_tot(rhs.e_tot), id(rhs.id), failcounter(rhs.failcounter)
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
  std::swap(failcounter, rhs.failcounter);
}

energy::interfaces::dftb::sysCallInterface::~sysCallInterface(void)
{

}





namespace
{
  int dftb_system_call(std::string const & command_line)
  {

  }
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
    
    std::string modulepath = Config::get().energy.dftb.path+"/DFTB";  //path to python module
    PySys_SetPath(const_cast<char*>(modulepath.c_str())); //set path
    const char *c = add_path.c_str();  //add paths from variable add_path
    PyRun_SimpleString(c);

    modul = PyImport_ImportModule("DFTB2_cast"); //import module 

    if(modul) 
        { 
        funk = PyObject_GetAttrString(modul, "main"); //create function
        prm = Py_BuildValue("(ss)", "tmp_struc.xyz", "dftbaby.cfg"); //give parameters
        ret = PyObject_CallObject(funk, prm);  //call function with parameters

        result_str = PyString_AsString(ret); //read function return (has to be a string)
        result_str = result_str.substr(1,result_str.size()-2);  //process return
        std::vector<std::string> result_vec = split(result_str, ',');

        //read energies and convert them to kcal/mol
        e_bs = std::stod(result_vec[0])*627.503; 
        e_coul = std::stod(result_vec[1])*627.503;
        e_rep = std::stod(result_vec[3])*627.503;
        e_lr = std::stod(result_vec[4])*627.503;
        e_tot = std::stod(result_vec[5])*627.503;
        
        //delete PyObjects
        Py_DECREF(prm); 
        Py_DECREF(ret); 
        Py_DECREF(funk); 
        Py_DECREF(modul); 
        } 
    else 
    {
        printf("Fehler: Modul DFTB2_cast nicht gefunden\n"); 
        std::exit(0);
    }
  return e_tot;
}

// Energy+Gradient function
double energy::interfaces::dftb::sysCallInterface::g(void)
{
  integrity = true;
  grad_var = true;
  std::cout<<"no gradients yet\n";
  return energy;
}

// Energy+Gradient+Hessian function
double energy::interfaces::dftb::sysCallInterface::h(void)
{
  integrity = true;
  grad_var = false;
  std::cout<<"no hessian yet\n";
  return energy;
}

// Optimization
double energy::interfaces::dftb::sysCallInterface::o(void)
{
  integrity = true;
  grad_var = false;
  std::cout<<"no optimizer yet\n";
  return energy;
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

void energy::interfaces::dftb::sysCallInterface::print_G_tinkerlike(std::ostream &, bool const) const { }

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