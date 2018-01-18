#include "energy_int_dftb.h"



energy::interfaces::dftb::sysCallInterface::sysCallInterface(coords::Coordinates * cp) :
  energy::interface_base(cp),energy(0.0)
{
    
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

void energy::interfaces::dftb::sysCallInterface::write_inputfile(int t)
{
  std::vector<std::string> elements;
  for (auto a : (*this->coords).atoms())
  {
    if (is_in(a.symbol(), elements) == false)
    {
      elements.push_back(a.symbol());
    }
  }

  std::ofstream file("dftb_in.hsd");
  file << "Geometry = GenFormat {\n";
  file << coords::output::formats::xyz_gen(*this->coords);
  file << "}\n\n";

  file << "Hamiltonian = DFTB {\n";
  file << "  SCC = Yes\n";
  file << "  SlaterKosterFiles = Type2FileNames {\n";
  file << "    Prefix = '"<<Config::get().energy.dftb.sk_files<<"'\n";
  file << "    Separator = '-'\n";
  file << "    Suffix = '.skf'\n";
  file << "  }\n";
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

  file << "Options {\n";
  file << "  WriteResultsTag = Yes\n";
  file << "  RestartFrequency = 0\n";
  if (Config::get().energy.dftb.verbosity < 2) file << "  WriteDetailedOut = No\n";
  file << "}\n\n";

  file << "Analysis = {\n";
  file << "  WriteBandOut = No\n";
  if (t == 1) file << "  CalculateForces = Yes\n";
  file << "}\n\n";

  file << "ParserOptions {\n";
  file << "  ParserVersion = 5\n";
  file << "}";
}

double energy::interfaces::dftb::sysCallInterface::read_output(int t)
{
  std::ifstream in_file("results.tag", std::ios_base::in);
  std::string line;
  
  while (!in_file.eof())
  {
    std::getline(in_file, line);
    if (line == "total_energy        :real:0:")
    {
      std::getline(in_file, line);
      energy = std::stod(line)*627.503;
    }

    else if (line.substr(0, 29) == "forces              :real:2:3" && t == 1)
    {
      double x, y, z;
      coords::Representation_3D g_tmp;
      
      for (int i = 0; i < (*this->coords).size(); i++)
      {
        std::getline(in_file, line);
        std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
        x = (-x) * (627.503 / 0.5291172107);
        y = (-y) * (627.503 / 0.5291172107);
        z = (-z) * (627.503 / 0.5291172107);
        coords::Cartesian_Point g(x,y,z);
        g_tmp.push_back(g);
      }
      coords->swap_g_xyz(g_tmp);
    }
  }

  if (Config::get().energy.dftb.verbosity < 2)
  {
    std::remove("dftb_pin.hsd");
    if (Config::get().energy.dftb.verbosity == 0)
    {
      std::remove("dftb_in.hsd");
      std::remove("output_dftb.hsd");
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
  write_inputfile(0);
  scon::system_call(Config::get().energy.dftb.path +" > output_dftb.txt");
  energy = read_output(0);

  return energy;
}

// Energy+Gradient function
double energy::interfaces::dftb::sysCallInterface::g(void)
{
  write_inputfile(1);
  scon::system_call(Config::get().energy.dftb.path + " > output_dftb.txt");
  energy = read_output(1);
  return energy;
}

// Hessian function
double energy::interfaces::dftb::sysCallInterface::h(void)
{
  throw std::runtime_error("Hessian not implemented in DFTBplus interface.");
}

// Optimization
double energy::interfaces::dftb::sysCallInterface::o(void)
{
  throw std::runtime_error("DFTBplus has no own optimizer.");
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

std::vector<coords::float_type>
energy::interfaces::dftb::sysCallInterface::charges() const
{
  std::vector<coords::float_type> charges;
  return charges;
}