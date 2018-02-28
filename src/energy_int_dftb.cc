#include "energy_int_dftb.h"

energy::interfaces::dftb::sysCallInterface::sysCallInterface(coords::Coordinates * cp) :
  energy::interface_base(cp),energy(0.0)
{
  if (Config::get().energy.dftb.opt > 0) optimizer = true;
  else optimizer = false;
}

energy::interfaces::dftb::sysCallInterface::sysCallInterface(sysCallInterface const & rhs, coords::Coordinates *cobj) :
  interface_base(cobj), energy(rhs.energy)
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

/**checks if all atom coordinates are numbers*/
bool energy::interfaces::dftb::sysCallInterface::check_structure()
{
  bool structure = true;
  double x, y, z;
  for (auto i : (*this->coords).xyz())
  {
    double x = i.x();
    double y = i.y();
    double z = i.z();

    if (std::isnan(x) || std::isnan(y) || std::isnan(z))
    {
      std::cout << "Atom coordinates are not a number. Treating structure as broken.\n";
      structure = false;
      break;
    }
  }
  return structure;
}

void energy::interfaces::dftb::sysCallInterface::write_inputfile(int t)
{
  // create a vector with all element symbols that are found in the structure
  // are needed for writing angular momenta into inputfile
  std::vector<std::string> elements;  
  for (auto a : (*this->coords).atoms())
  {
    if (is_in(a.symbol(), elements) == false)
    {
      elements.push_back(a.symbol());
    }
  }

  // create inputfile
  std::ofstream file("dftb_in.hsd");

  // write geometry
  file << "Geometry = GenFormat {\n";
  file << coords::output::formats::xyz_gen(*this->coords);
  file << "}\n\n";

  if (t == 2)  // hessian matrix is to be calculated
  {
    file << "Driver = SecondDerivatives {}\n\n";
  }
  else if (t == 3) // optimization is to be performed
  {
    std::string driver;
    if (Config::get().energy.dftb.opt == 1) driver = "SteepestDescent";
    else if (Config::get().energy.dftb.opt == 2) driver = "ConjugateGradient";
    else throw std::runtime_error("Cannot write correct inputfile for DFTB+. Unknown driver\n");

    file << "Driver = " << driver << " {\n";
    file << "  MaxSteps = " << Config::get().energy.dftb.max_steps_opt << "\n";
    file << "}\n\n";
  }

  // write information that is needed for SCC calculation
  file << "Hamiltonian = DFTB {\n";
  file << "  SCC = Yes\n";
  file << "  SCCTolerance = " <<std::scientific << Config::get().energy.dftb.scctol << "\n";
  file << "  MaxSCCIterations = " << Config::get().energy.dftb.max_steps << "\n";
  file << "  Charge = " << Config::get().energy.dftb.charge << "\n";
  file << "  SlaterKosterFiles = Type2FileNames {\n";
  file << "    Prefix = '"<<Config::get().energy.dftb.sk_files<<"'\n";
  file << "    Separator = '-'\n";
  file << "    Suffix = '.skf'\n";
  file << "  }\n";
  if (Config::get().energy.qmmm.use == true)  // if QM/MM calculation: tell DFTB+ to read external charges
  {                          
    file << "  ElectricField = {\n";
    file << "    PointCharges = {\n";
    file << "      CoordsAndCharges [Angstrom] = DirectRead {\n";
    file << "        Records = " << Config::get().energy.qmmm.mm_atoms_number << "\n";
    file << "        File = 'charges.dat'\n";
    file << "      }\n";
    file << "    }\n";
    file << "  }\n";
  }
  file << "  MaxAngularMomentum {\n";
  for (auto s : elements)
  {
    char angular_momentum = atomic::angular_momentum_by_symbol(s);
    if (angular_momentum == 'e')
    {
      std::cout << "Angular momentum for element " << s << " not defined. \n";
      std::cout << "Please go to file 'atomic.h' and define an angular momentum (s, p, d or f) in the array 'angular_momentum'.\n";
      std::cout << "Talk to a CAST developer if this is not possible for you.";
      std::exit(0);
    }
    file << "    " << s << " = " << angular_momentum << "\n";
  }
  file << "  }\n";
  file << "}\n\n";

  // which information will be saved after calculation?
  file << "Options {\n";
  file << "  WriteResultsTag = Yes\n";
  if (t < 2) file << "  RestartFrequency = 0\n";   // does not work together with driver
  if (Config::get().energy.dftb.verbosity < 2) file << "  WriteDetailedOut = No\n";
  file << "}\n\n";

  // additional analysis that should be performed
  file << "Analysis = {\n";
  file << "  WriteBandOut = No\n";
  if (t == 1) file << "  CalculateForces = Yes\n";
  file << "}\n\n";

  // parser version (recommended so it is possible to use newer DFTB+ versions without adapting inputfile)
  file << "ParserOptions {\n";
  file << "  ParserVersion = 5\n";
  file << "}";
}

double energy::interfaces::dftb::sysCallInterface::read_output(int t)
{
  if (file_exists("results.tag") == false) // if SCC does not converge DFTB+ doesn't produce this file
  {
    std::cout << "DFTB+ did not produce an output file. Treating structure as broken.\n";
    integrity = false;
  }

  // successfull SCC -> read output
  else
  {
    int N = (*this->coords).size();

    std::ifstream in_file("results.tag", std::ios_base::in);
    std::string line;

    while (!in_file.eof())
    {
      std::getline(in_file, line);
      if (line == "total_energy        :real:0:")  // read energy
      {
        std::getline(in_file, line);
        energy = std::stod(line)*627.503; // convert hartree to kcal/mol
      }

      else if (line.substr(0, 29) == "forces              :real:2:3" && t == 1)  // read gradients
      {
        double x, y, z;
        coords::Representation_3D g_tmp;

        for (int i = 0; i < N; i++)
        {
          std::getline(in_file, line);
          std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
          x = (-x) * (627.503 / 0.5291172107);  // hartree/bohr -> kcal/(mol*A)
          y = (-y) * (627.503 / 0.5291172107);
          z = (-z) * (627.503 / 0.5291172107);
          coords::Cartesian_Point g(x, y, z);
          g_tmp.push_back(g);
        }
        coords->swap_g_xyz(g_tmp);  // set gradients
      }

      else if (Config::get().energy.qmmm.use == true)
      {    // in case of QM/MM calculation: read forces on external charges
        if (line.substr(0, 29) == "forces_ext_charges  :real:2:3")
        {
          std::vector<coords::Cartesian_Point> grad_tmp;
          double x, y, z;

          for (int i = 0; i < Config::get().energy.qmmm.mm_atoms_number; i++)
          {
            std::getline(in_file, line);
            std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
            x = (-x) * (627.503 / 0.5291172107);  // hartree/bohr -> kcal/(mol*A)
            y = (-y) * (627.503 / 0.5291172107);
            z = (-z) * (627.503 / 0.5291172107);
            coords::Cartesian_Point g(x, y, z);
            grad_tmp.push_back(g);
          }
          grad_ext_charges = grad_tmp;
        }
      }

      else if (line.substr(0, 27) == "hessian_numerical   :real:2" && t == 2)  // read hessian
      {
        double CONVERSION_FACTOR = 627.503 / (0.5291172107*0.5291172107); // hartree/bohr^2 -> kcal/(mol*A^2)

        // read all values into one vector (tmp)
        std::vector<double> tmp;
        double x, y, z;
        for (int i=0; i < (3*N*3*N)/3; i++)  
        {
          std::getline(in_file, line);
          std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
          tmp.push_back(x*CONVERSION_FACTOR);
          tmp.push_back(y*CONVERSION_FACTOR);
          tmp.push_back(z*CONVERSION_FACTOR);
        }

        // format values in vector tmp into matrix hess
        std::vector<std::vector<double>> hess;
        std::vector<double> linevec;
        for (int i = 0; i < 3*N; i++)
        {
          linevec.resize(0);
          for (int j = 0; j < 3*N; j++)
          {
            linevec.push_back(tmp[(i*3*N)+j]);
          }
          hess.push_back(linevec);
        }

        coords->set_hessian(hess);  //set hessian
        std::remove("hessian.out"); // delete file
      }
    }

    if (t == 3)  // optimization
    {
      if (file_exists("geo_end.gen") == false)  // sometimes happens when moving randomly (tasks MC and TS)
      {
        std::cout << "DFTB+ did not produce a geometry file. Treating structure as broken.\n";
        integrity = false;
      }
      else  // if optimized geometry present -> read geometry from gen-file
      {
        std::ifstream geom_file("geo_end.gen", std::ios_base::in);

        std::getline(geom_file, line);  
        std::getline(geom_file, line);

        coords::Representation_3D xyz_tmp;
        std::string number, type;
        double x, y, z;

        while (geom_file >> number >> type >> x >> y >> z) 
        {
          coords::Cartesian_Point xyz(x, y, z);
          xyz_tmp.push_back(xyz);
        }
        geom_file.close();
        coords->set_xyz(std::move(xyz_tmp));  // set new coordinates

        std::remove("geo_end.gen"); // delete file
        std::remove("geo_end.xyz"); // delete file

        std::ifstream file("output_dftb.txt");  // check for geometry convergence
        line = last_line(file);
        if (line == "Geometry converged") {}
        else
        {
          std::cout << "Geometry did not converge. Treating structure as broken.\n";
          std::cout << "If this is a problem for you use a bigger value for 'DFTB+max_steps_opt'!\n";
          integrity = false;
        }
      }
    }
  }

  // check if geometry is still intact
  if (check_bond_preservation() == false) integrity = false;
  else if (check_atom_dist() == false) integrity = false;
  
  // remove files
  if (t > 1) std::remove("charges.bin");
  if (Config::get().energy.dftb.verbosity < 2)
  {
    std::remove("dftb_pin.hsd");
    if (Config::get().energy.dftb.verbosity == 0)
    {
      std::remove("dftb_in.hsd");
      std::remove("output_dftb.txt");
    }
  }
   
  return energy;
}

/*
Energy class functions that need to be overloaded
*/

// Energy function
double energy::interfaces::dftb::sysCallInterface::e(void)
{
  integrity = check_structure();
  if (integrity == true)
  {
    write_inputfile(0);
    scon::system_call(Config::get().energy.dftb.path + " > output_dftb.txt");
    energy = read_output(0);
    return energy;
  }
}

// Energy+Gradient function
double energy::interfaces::dftb::sysCallInterface::g(void)
{
  integrity = check_structure();
  if (integrity == true)
  {
    write_inputfile(1);
    scon::system_call(Config::get().energy.dftb.path + " > output_dftb.txt");
    energy = read_output(1);
    return energy;
  }
}

// Hessian function
double energy::interfaces::dftb::sysCallInterface::h(void)
{
  integrity = check_structure();
  if (integrity == true)
  {
    write_inputfile(2);
    scon::system_call(Config::get().energy.dftb.path + " > output_dftb.txt");
    energy = read_output(2);
    return energy;
  } 
}

// Optimization
double energy::interfaces::dftb::sysCallInterface::o(void)
{
  integrity = check_structure();
  if (integrity == true)
  {
    write_inputfile(3);
    scon::system_call(Config::get().energy.dftb.path + " > output_dftb.txt");
    energy = read_output(3);
    return energy;
  }
}

// Output functions
void energy::interfaces::dftb::sysCallInterface::print_E(std::ostream &S) const
{
  S << "Total Energy:      ";
  S << std::right << std::setw(16) << std::fixed << std::setprecision(8) << energy;
}

void energy::interfaces::dftb::sysCallInterface::print_E_head(std::ostream &S, bool const endline) const
{
  S << "Energies\n";
  S << std::right << std::setw(24) << "E_bs";
  S << std::right << std::setw(24) << "E_coul";
  S << std::right << std::setw(24) << "E_rep";
  S << std::right << std::setw(24) << "SUM\n\n";
}

void energy::interfaces::dftb::sysCallInterface::print_E_short(std::ostream &S, bool const endline) const
{
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << energy << '\n';
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

bool energy::interfaces::dftb::sysCallInterface::check_atom_dist(void) const
{
  std::size_t const N(coords->size());
  for (std::size_t i(0U); i < N; ++i)
  {
    for (std::size_t j(0U); j < i; j++)
    {
      if (dist(coords->xyz(i), coords->xyz(j)) < 0.3)
      {
        return false;
      }
    }
  }
  return true;
}

std::vector<coords::float_type>
energy::interfaces::dftb::sysCallInterface::charges() const
{
  std::vector<coords::float_type> charges;

  auto in_string = "results.tag";
  if (file_exists(in_string) == false)
  {
    throw std::runtime_error("DFTB+ logfile results.tag not found.");
  }

  else
  {
    std::ifstream in_file("results.tag", std::ios_base::in);
    std::string line;
    double q;

    while (!in_file.eof())
    {
      std::getline(in_file, line);
      if (line.substr(0, 27) == "net_atomic_charges  :real:1")  // read energy
      {
        for (int i = 0; i < coords->size(); i++)
        {
          in_file >> q;
          charges.push_back(q);
        }
      }
    }
  }
  return charges;
}

std::vector<coords::Cartesian_Point>
energy::interfaces::dftb::sysCallInterface::get_g_coul_mm() const
{
  return grad_ext_charges;
}