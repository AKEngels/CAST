#include "energy_int_dftb.h"

energy::interfaces::dftb::sysCallInterface::sysCallInterface(coords::Coordinates* cp) :
  energy::interface_base(cp), energy(0.0)
{
  if (Config::get().energy.dftb.opt > 0) optimizer = true;
  else optimizer = false;
  charge = Config::get().energy.dftb.charge;
}

energy::interfaces::dftb::sysCallInterface::sysCallInterface(sysCallInterface const& rhs, coords::Coordinates* cobj) :
  interface_base(cobj), energy(rhs.energy)
{
  optimizer = rhs.optimizer;
  charge = rhs.charge;
  interface_base::operator=(rhs);
}

energy::interface_base* energy::interfaces::dftb::sysCallInterface::clone(coords::Coordinates* coord_object) const
{
  sysCallInterface* tmp = new sysCallInterface(*this, coord_object);
  return tmp;
}

energy::interface_base* energy::interfaces::dftb::sysCallInterface::move(coords::Coordinates* coord_object)
{
  sysCallInterface* tmp = new sysCallInterface(*this, coord_object);
  return tmp;
}

void energy::interfaces::dftb::sysCallInterface::swap(interface_base& rhs)
{
  swap(dynamic_cast<sysCallInterface&>(rhs));
}

void energy::interfaces::dftb::sysCallInterface::swap(sysCallInterface& rhs)
{
  interface_base::swap(rhs);
}

energy::interfaces::dftb::sysCallInterface::~sysCallInterface(void) {}

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
    for (auto j = 0u; j < charge_vector.size(); j++)
    {
      chargefile << charge_vector[j].x << " " << charge_vector[j].y << " " << charge_vector[j].z << "  " << charge_vector[j].scaled_charge << "\n";
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
  file << "  SCCTolerance = " << std::scientific << Config::get().energy.dftb.scctol << "\n";
  file << "  MaxSCCIterations = " << Config::get().energy.dftb.max_steps << "\n";
  file << "  Charge = " << charge << "\n";
  file << "  SlaterKosterFiles = Type2FileNames {\n";
  file << "    Prefix = '" << Config::get().energy.dftb.sk_files << "'\n";
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
  if (Config::get().periodics.periodic)   // if periodic boundaries are used, k points must be set
  {
    file << "  KPointsAndWeights = {\n";
    file << "    0.0 0.0 0.0  1.0\n";     // we only set gamma point
    file << "  }\n";
  }
  if (Config::get().energy.dftb.d3)   // D3 correction
  {
    std::array<double, 4> params;  // which params to use?
    if (Config::get().energy.dftb.d3param > 3) {
      throw std::runtime_error("Unvalid parameters for D3 correction in DFTB+ interface!");
    }
    else params = D3PARAMS[Config::get().energy.dftb.d3param];
    file << "  Dispersion = DftD3 {\n";
    file << "    Damping = BeckeJohnson {\n";
    file << "      a1 = " << params[0] << "\n";
    file << "      a2 = " << params[1] << "\n";
    file << "    }\n";
    file << "    s6 = " << params[2] << "\n";
    file << "    s8 = " << params[3] << "\n";
    file << "  }\n";
  }
  if (Config::get().energy.dftb.dftb3)   // DFTB3
  {
    file << "  ThirdOrderFull = Yes\n";
    file << "  HubbardDerivs {\n";
    for (auto s : elements)
    {
      double hubbard_deriv = hubbard_deriv_by_symbol(s);
      file << "    " << s << " = " << hubbard_deriv << "\n";
    }
    file << "  }\n";
    file << "  HCorrection = Damping {\n";
    file << "    Exponent = " << get_zeta() << "\n";
    file << "  }\n";
  }
  if (Config::get().energy.dftb.range_sep)   // range separation
  {
    file << "  RangeSeparated = LC {\n";
    file << "    Screening = Thresholded { }\n";
    file << "  }\n";
  }
  if (Config::get().energy.dftb.fermi_temp > 0)  // Fermi filling
  {
    file << "  Filling = Fermi {\n";
    file << "    Temperature [K] = " << Config::get().energy.dftb.fermi_temp << "\n";
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
  file << "  ParserVersion = 7\n";
  file << "}";
  file.close();
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
        energy = std::stod(line) * energy::au2kcal_mol; // convert hartree to kcal/mol
      }

      else if (line.substr(0, 29) == "forces              :real:2:3" && t == 1)  // read gradients
      {
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
      }

      else if (Config::get().energy.qmmm.mm_charges.size() != 0)
      {    // in case of QM/MM calculation: read forces on external charges
        if (line.substr(0, 29) == "forces_ext_charges  :real:2:3")
        {
          int ext_charge_number = std::stoi(line.substr(30));

          std::vector<coords::Cartesian_Point> grad_tmp;
          double x, y, z;

          for (int i = 0; i < ext_charge_number; i++)
          {
            std::getline(in_file, line);
            std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
            x *= -energy::Hartree_Bohr2Kcal_MolAng;  // hartree/bohr -> kcal/(mol*A)
            y *= -energy::Hartree_Bohr2Kcal_MolAng;
            z *= -energy::Hartree_Bohr2Kcal_MolAng;
            coords::Cartesian_Point g(x, y, z);
            grad_tmp.push_back(g);
          }
          grad_ext_charges = grad_tmp;
        }
      }

      if (line.substr(0, 27) == "gross_atomic_charges:real:1")     // read atomic charges
      {
        partial_charges.clear();
        double q;
        for (auto i = 0u; i < coords->size(); i++)
        {
          in_file >> q;
          partial_charges.push_back(q);
        }
      }

      else if (line.substr(0, 27) == "hessian_numerical   :real:2" && t == 2)  // read hessian
      {
        double CONVERSION_FACTOR = energy::Hartree_Bohr2Kcal_MolAngSquare; // hartree/bohr^2 -> kcal/(mol*A^2)

                                                                           // read all values into one vector (tmp)
        std::vector<double> tmp;
        double x, y, z;
        for (int i = 0; i < (3 * N * 3 * N) / 3; i++)
        {
          std::getline(in_file, line);
          std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
          tmp.push_back(x * CONVERSION_FACTOR);
          tmp.push_back(y * CONVERSION_FACTOR);
          tmp.push_back(z * CONVERSION_FACTOR);
        }

        // format values in vector tmp into matrix hess
        std::vector<std::vector<double>> hess;
        std::vector<double> linevec;
        for (int i = 0; i < 3 * N; i++)
        {
          linevec.resize(0);
          for (int j = 0; j < 3 * N; j++)
          {
            linevec.push_back(tmp[(i * 3 * N) + j]);
          }
          hess.push_back(linevec);
        }

        coords->set_hessian(hess);  //set hessian
        std::remove("hessian.out"); // delete file
      }
    }
    in_file.close();

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
        bool converged = false;
        while (!file.eof())
        {
          std::getline(file, line);
          if (line == "Geometry converged") {
            converged = true;
            break;
          }
        }
        if (converged == false)
        {
          std::cout << "Geometry did not converge. Treating structure as broken.\n";
          std::cout << "If this is a problem for you use a bigger value for 'DFTB+max_steps_opt'!\n";
          integrity = false;
        }
      }
    }
  }

  // check if geometry is still intact
  if (coords->check_bond_preservation() == false) integrity = false;
  else if (coords->check_for_crashes() == false) integrity = false;

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
  integrity = coords->check_structure();
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
  integrity = coords->check_structure();
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
  integrity = coords->check_structure();
  throw std::runtime_error("Hessian for DFTB-Interface not implemented in CAST as yet.");
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
  integrity = coords->check_structure();
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
void energy::interfaces::dftb::sysCallInterface::print_E(std::ostream& S) const
{
  S << "Total Energy:      ";
  S << std::right << std::setw(16) << std::fixed << std::setprecision(8) << energy;
}

void energy::interfaces::dftb::sysCallInterface::print_E_head(std::ostream& S, bool const endline) const
{
  S << "Energies\n";
  S << std::right << std::setw(24) << "E_bs";
  S << std::right << std::setw(24) << "E_coul";
  S << std::right << std::setw(24) << "E_rep";
  S << std::right << std::setw(24) << "SUM\n";
  if (endline) S << '\n';
}

void energy::interfaces::dftb::sysCallInterface::print_E_short(std::ostream& S, bool const endline) const
{
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << 0;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << energy << '\n';
  if (endline) S << '\n';
}

void energy::interfaces::dftb::sysCallInterface::to_stream(std::ostream&) const { }


