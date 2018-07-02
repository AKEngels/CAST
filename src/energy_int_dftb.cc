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
    x = i.x();
    y = i.y();
    z = i.z();

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

  // create a chargefile for external charges if desired (needed for QM/MM methods)
  if (Config::get().energy.qmmm.mm_charges.size() != 0)
  {
	  std::vector<PointCharge> charge_vector = Config::get().energy.qmmm.mm_charges;
	  std::ofstream chargefile("charges.dat");
	  for (int j = 0; j < charge_vector.size(); j++)
	  {
		  chargefile << charge_vector[j].x << " " << charge_vector[j].y << " " << charge_vector[j].z << "  " << charge_vector[j].charge << "\n";
	  }
	  chargefile.close();
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

  if (Config::get().energy.qmmm.mm_charges.size() != 0)  // include chargefile
  {
	  file << "  ElectricField = {\n";
	  file << "    PointCharges = {\n";
	  file << "      CoordsAndCharges [Angstrom] = DirectRead {\n";
	  file << "        Records = " << Config::get().energy.qmmm.mm_charges.size() << "\n";  
	  file << "        File = 'charges.dat'\n";
	  file << "      }\n";
	  file << "    }\n";
	  file << "  }\n";
  }

  file << "  MaxAngularMomentum {\n";
  for (auto s : elements)
  {
    char angular_momentum = angular_momentum_by_symbol(s);
    file << "    " << s << " = " << angular_momentum << "\n";
  }
  file << "  }\n";
  if (Config::get().energy.dftb.dftb3 == true)
  {
    file << "  ThirdOrderFull = Yes\n";
    file << "  DampXH = Yes\n";
    file << "  DampXHExponent = " << get_zeta() << "\n";
    file << "  HubbardDerivs {\n";
    for (auto s : elements)
    {
      double hubbard_deriv = hubbard_deriv_by_symbol(s);
      file << "    " << s << " = " << hubbard_deriv << "\n";
    }
    file << "  }\n";
  }
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
  std::string res_filename{ "results.tag" };
  if (file_is_empty(res_filename)) // if SCC does not converge this file is empty
  {
    std::cout << "DFTB+ produced an empty output file. Treating structure as broken.\n";
    integrity = false;
		return 0.0;      // return zero-energy because no energy was calculated
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
        energy = std::stod(line)*energy::au2kcal_mol; // convert hartree to kcal/mol
      }

      else if (line.substr(0, 29) == "forces              :real:2:3" && t == 1)  // read gradients
      {
        int link_atom_number = std::stoi(line.substr(30)) - N;
        double x, y, z;
        coords::Representation_3D g_tmp;

        for (int i = 0; i < N; i++)  // read "normal" gradients
        {
          std::getline(in_file, line);
          std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
          x *= -energy::Hartree_Bohr2Kcal_MolAng;
          y *= -energy::Hartree_Bohr2Kcal_MolAng;
          z *= -energy::Hartree_Bohr2Kcal_MolAng;
          coords::Cartesian_Point g(x, y, z);
          g_tmp.push_back(g);
        }
        coords->swap_g_xyz(g_tmp);  // set gradients

        for (int i=0; i < link_atom_number; i++)  // read gradients of link atoms
        {
          std::getline(in_file, line);
          std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
          x *= -energy::Hartree_Bohr2Kcal_MolAng;  // hartree/bohr -> kcal/(mol*A)
          y *= -energy::Hartree_Bohr2Kcal_MolAng;
          z *= -energy::Hartree_Bohr2Kcal_MolAng;
          coords::Cartesian_Point g(x, y, z);
          link_atom_grad.push_back(g);
        }
      }

      else if (line.substr(0, 27) == "hessian_numerical   :real:2" && t == 2)  // read hessian
      {
        double CONVERSION_FACTOR = energy::Hartree_Bohr2Kcal_MolAngSquare; // hartree/bohr^2 -> kcal/(mol*A^2)

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

  if (Config::get().energy.qmmm.mm_charges.size() != 0)
  {
    double ext_chg_energy = calc_self_interaction_of_external_charges();  // calculates self interaction energy of the external charges
    energy -= ext_chg_energy;  // subtract self-interaction because it's already in MM calculation
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
  else return 0;  // energy = 0 if structure contains NaN
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
  else return 0;  // energy = 0 if structure contains NaN
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
  else return 0;  // energy = 0 if structure contains NaN
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
  else return 0;  // energy = 0 if structure contains NaN
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
  S << std::right << std::setw(24) << "SUM\n";
  if (endline) S << '\n';
}

void energy::interfaces::dftb::sysCallInterface::print_E_short(std::ostream &S, bool const endline) const
{
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << energy << '\n';
  if (endline) S << '\n';
}

void energy::interfaces::dftb::sysCallInterface::to_stream(std::ostream&) const { }

double energy::interfaces::dftb::sysCallInterface::calc_self_interaction_of_external_charges()
{
  double energy{ 0.0 };
  for (int i=0; i < Config::get().energy.qmmm.mm_charges.size(); ++i)
  {
    auto c1 = Config::get().energy.qmmm.mm_charges[i];
    for (int j = 0; j < i; ++j)
    {
      auto c2 = Config::get().energy.qmmm.mm_charges[j];

      double dist_x = c1.x - c2.x;
      double dist_y = c1.y - c2.y;
      double dist_z = c1.z - c2.z;
      double dist = std::sqrt(dist_x*dist_x + dist_y * dist_y + dist_z * dist_z);  // distance in angstrom

      energy += 332.0 * c1.charge * c2.charge / dist;  // energy in kcal/mol
    }
  }
  return energy;
}

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
        double const max = 1.2 * (coords->atoms(i).cov_radius() + coords->atoms(coords->atoms(i).bonds(j)).cov_radius());
        if (L > max) return false;
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
      if (line.substr(0, 27) == "gross_atomic_charges:real:1")
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
energy::interfaces::dftb::sysCallInterface::get_g_ext_chg() const
{
  auto elec_factor = 332.0;  // factor for conversion of charge product into amber units
  auto atom_charges = charges();

  std::vector<coords::Cartesian_Point> grad_ext_charges;
  grad_ext_charges.resize(Config::get().energy.qmmm.mm_charges.size());

  for (int i = 0; i < coords->size(); ++i)  // for every atom
  {
    double charge_i = atom_charges[i];

    for (int j = 0; j < Config::get().energy.qmmm.mm_charges.size(); ++j)  // for every external charge
    {
      auto current_charge = Config::get().energy.qmmm.mm_charges[j];
      double charge_j = current_charge.charge;

      auto dx = current_charge.x - coords->xyz(i).x();
      auto dy = current_charge.y - coords->xyz(i).y();
      auto dz = current_charge.z - coords->xyz(i).z();
      auto r_ij = coords::r3{ dx, dy, dz };   // vector between atom and charge
      coords::float_type d = len(r_ij);     // distance between atom and charge

      coords::float_type db = -elec_factor * (charge_i*charge_j) / (d*d);  // derivative of coulomb energy (only number)
      auto c_gradient_ij = (r_ij / d) * db;                                // now gradient gets a direction
      grad_ext_charges[j] += c_gradient_ij;   // add gradient 
    }
  }
  return grad_ext_charges;
}

