#include "energy_int_orca.h"

energy::interfaces::orca::sysCallInterface::sysCallInterface(coords::Coordinates * cp) :
  energy::interface_base(cp), energy(0.0), nuc_rep(0.0), elec_en(0.0), one_elec(0.0), two_elec(0.0)
{
  if (Config::get().energy.orca.opt > 0) optimizer = true;
  else optimizer = false;
	charge = Config::get().energy.orca.charge;
}

energy::interfaces::orca::sysCallInterface::sysCallInterface(sysCallInterface const & rhs, coords::Coordinates *cobj) :
  interface_base(cobj), energy(rhs.energy), nuc_rep(rhs.nuc_rep), elec_en(rhs.elec_en), one_elec(rhs.one_elec), two_elec(rhs.two_elec)
{
	optimizer = rhs.optimizer;
	charge = rhs.charge;
  interface_base::operator=(rhs);
}

energy::interface_base * energy::interfaces::orca::sysCallInterface::clone(coords::Coordinates * coord_object) const
{
  sysCallInterface * tmp = new sysCallInterface(*this, coord_object);
  return tmp;
}

energy::interface_base * energy::interfaces::orca::sysCallInterface::move(coords::Coordinates * coord_object)
{
  sysCallInterface * tmp = new sysCallInterface(*this, coord_object);
  return tmp;
}

void energy::interfaces::orca::sysCallInterface::swap(interface_base &rhs)
{
  swap(dynamic_cast<sysCallInterface&>(rhs));
}

void energy::interfaces::orca::sysCallInterface::swap(sysCallInterface &rhs)
{
  interface_base::swap(rhs);
}

energy::interfaces::orca::sysCallInterface::~sysCallInterface(void) {}

/**checks if all atom coordinates are numbers*/
bool energy::interfaces::orca::sysCallInterface::check_structure()
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

void energy::interfaces::orca::sysCallInterface::write_external_pointcharges(std::string const &filename)
{
  std::vector<PointCharge> chargevec = Config::get().energy.qmmm.mm_charges;

  std::ofstream pointcharges;
  pointcharges.open(filename);

  pointcharges << chargevec.size() <<"\n";
  for (auto &charge : chargevec)
  {
    pointcharges << charge.charge << "  " << charge.x << "  " << charge.y << "  " << charge.z << "\n";
  }
  pointcharges.close();
}
  

void energy::interfaces::orca::sysCallInterface::write_inputfile(int t)
{
  if (Config::get().energy.qmmm.mm_charges.size() != 0) write_external_pointcharges("pointcharges.pc");  // write external pointcharges into file

	std::ofstream inp;
	inp.open("orca.inp");

	inp << "! " << Config::get().energy.orca.method << " " << Config::get().energy.orca.basisset << "\n";  // method and basisset
	if (t == 1) inp << "! EnGrad\n";                                                                       // request gradients
	if (t == 2) inp << "! Freq\n";                                                                         // request hessian
	if (t == 3) inp << "! Opt\n";                                                                          // request optimization

  if (Config::get().energy.qmmm.mm_charges.size() != 0) inp << "\n% pointcharges \"pointcharges.pc\"\n";     // tell orca that there are pointcharges in this file

	inp << "\n";  // empty line
	inp << "*xyz " << charge << " " << Config::get().energy.orca.multiplicity << "\n";  // headline for geometry input
	for (auto i{ 0u }; i < coords->atoms().size(); ++i)  // coordinates definition for every atom
	{
		inp << coords->atoms(i).symbol() << "  " << std::setw(12) << std::setprecision(6) << coords->xyz(i).x()
			<< "  " << std::setw(12) << std::setprecision(6) << coords->xyz(i).y()
			<< "  " << std::setw(12) << std::setprecision(6) << coords->xyz(i).z() << "\n";
	}
	inp << "*\n";  // end of coordinates definition
  inp.close();
}

double energy::interfaces::orca::sysCallInterface::read_output(int t)
{
  if (file_exists("output_orca.txt") == false)   // if orca produced no outputfile
  {
    std::cout << "ORCA output file not present. Integrity is broken\n";
    integrity = false;
    return 0;
  }

  int N = (*this->coords).size();  // number of atoms

	std::ifstream out;
	out.open("output_orca.txt");

	std::string line;
	std::vector<std::string> linevec;
	while (out.eof() == false)          // reading output file
	{
		std::getline(out, line);  // reading

		if (line.substr(0, 25) == "FINAL SINGLE POINT ENERGY")  // get single point energy
		{
			linevec = split(line, ' ', true);
			energy = std::stod(linevec[4]) * energy::au2kcal_mol;
		}

		else if (line.substr(0, 17) == "Nuclear Repulsion")    // get partial energies
		{
			linevec = split(line, ' ', true);
			nuc_rep = std::stod(linevec[3]) * energy::au2kcal_mol;
		}
		else if (line.substr(0, 17) == "Electronic Energy")   
		{
			linevec = split(line, ' ', true);
			elec_en = std::stod(linevec[3]) * energy::au2kcal_mol;
		}
		else if (line.substr(0, 19) == "One Electron Energy")   
		{
			linevec = split(line, ' ', true);
			one_elec = std::stod(linevec[3]) * energy::au2kcal_mol;
		}
		else if (line.substr(0, 19) == "Two Electron Energy")    
		{
			linevec = split(line, ' ', true);
			two_elec = std::stod(linevec[3]) * energy::au2kcal_mol;
		}

		if (t == 1)   // if gradients requested
		{
			if (line.substr(0, 18) == "CARTESIAN GRADIENT")  // get gradients
			{
				std::getline(out, line);    // ------------------
				std::getline(out, line);    // empty line

				double x, y, z;
				coords::Representation_3D g_tmp;

				for (int i = 0; i < N; ++i)     // for every atom
				{
					std::getline(out, line);
					linevec = split(line, ' ', true);
					x = std::stod(linevec[3]) * energy::Hartree_Bohr2Kcal_MolAng;
					y = std::stod(linevec[4]) * energy::Hartree_Bohr2Kcal_MolAng;
					z = std::stod(linevec[5]) * energy::Hartree_Bohr2Kcal_MolAng;

					coords::Cartesian_Point g(x, y, z);
					g_tmp.push_back(g);
				}
				coords->swap_g_xyz(g_tmp);  // set gradients
			}
		}

    if (line.substr(0, 23) == "MULLIKEN ATOMIC CHARGES")  // reading charges
    {
      mulliken_charges.clear(); // delete whatever is in there before
      std::getline(out, line);  // '---------------'

      std::string buffer;
      double charge;
      for (int i = 0; i < N; ++i)
      {
        out >> buffer >> buffer >> buffer >> charge;
        mulliken_charges.emplace_back(charge);
      }
    }
	}
  out.close();

	if (t == 2)
	{
		if (file_exists("orca.hess") == true) read_hessian_from_file("orca.hess");
    else std::cout<<"ORCA hessian file not present\n";
	}

  if (t == 3)        // if optimization requested
  {
    if (file_exists("orca.xyz") == false)
    {
      std::cout << "Optimization produced no output file.\n";
      integrity = false;
    }
    else
    {
      std::ifstream geom_file;
      geom_file.open("orca.xyz");

      std::getline(geom_file, line);        // first line: number of atoms
      int number_of_atoms = std::stoi(line);
      if (N != number_of_atoms) throw std::runtime_error("wrong number of atoms in structure");
      std::getline(geom_file, line);        // second line: stuff

      std::string element;
      double x, y, z;
      coords::Representation_3D xyz_tmp;

      while (geom_file >> element >> x >> y >> z)
      {
        coords::Cartesian_Point xyz(x, y, z);
        xyz_tmp.push_back(xyz);
      }
      coords->set_xyz(std::move(xyz_tmp));  // set new coordinates
      geom_file.close();
    }
  }

  if (t == 1 && Config::get().energy.qmmm.mm_charges.size() != 0)      // read gradients on external point charges
  {
    grad_ext_charges.clear();  // delete former gradients

    if (file_exists("orca.pcgrad") == false) throw std::runtime_error("can't read gradients on external point charges from file 'orca.pcgrad'");

    std::ifstream pcgrad;
    pcgrad.open("orca.pcgrad");

    unsigned number_of_pointcharges;         // read number of point charges
    pcgrad >> number_of_pointcharges;
    if (number_of_pointcharges != Config::get().energy.qmmm.mm_charges.size()) throw std::runtime_error("wrong number of gradients on external point charges");

    double x, y, z;                                  // read gradients
    for (auto i = 0u; i < number_of_pointcharges; ++i)
    {
      pcgrad >> x >> y >> z;                    // read
      x = x * energy::Hartree_Bohr2Kcal_MolAng; // convert to correct units
      y = y * energy::Hartree_Bohr2Kcal_MolAng;
      z = z * energy::Hartree_Bohr2Kcal_MolAng;
      coords::Cartesian_Point grad{ x,y,z };    // create gradient on one atom
      grad_ext_charges.emplace_back(grad);      // push gradient into vector
    }

    pcgrad.close();
  }

  // check if geometry is still intact
  if (check_bond_preservation() == false) integrity = false;
  else if (check_atom_dist() == false) integrity = false;

  // deleting files
  if (Config::get().energy.orca.verbose < 4)
  {
    std::remove("orca.ges");
    std::remove("orca.prop");
  }
  if (Config::get().energy.orca.verbose < 3)
  {
    std::remove("orca.opt");
    std::remove("orca.trj");
    std::remove("orca.engrad");
    std::remove("orca_property.txt");
  }
  if (Config::get().energy.orca.verbose < 2)
  {
    std::remove("orca.xyz");
    std::remove("orca.hess");
  }
  if (Config::get().energy.orca.verbose < 1)
  {
    std::remove("output_orca.txt");
    std::remove("orca.inp");
  }

  return energy;
}

void energy::interfaces::orca::sysCallInterface::read_hessian_from_file(std::string const &filename)
{
	std::ifstream hess;
	hess.open(filename);

	std::string line;
	std::vector<std::string> linevec;
	while (hess.eof() == false)
	{
		std::getline(hess, line);

		if (line.substr(0, 8) == "$hessian")
		{
			std::getline(hess, line);  // get size if hessian (= 3x number of atoms)
			unsigned size = stoi(line) ;
			if (size != 3 * coords->size()) throw std::runtime_error("wrong size of hessian: " + size);

			else   // if correct size
			{
				// initialize hessian matrix
				std::vector<double> hessian_line;
				hessian_line.resize(size);
				std::vector<std::vector<double>> hessian;
				for (auto i = 0u; i < size; ++i) hessian.push_back(hessian_line);

				// start reading
				std::getline(hess, line);         // headline of first block
				linevec = split(line, ' ', true);  

				// get a bit of information how much is there to read
				int columns = linevec.size();
				int number_of_reading_blocks = size / columns;
				int number_of_columns_in_last_block = size % columns;

				// really reading numbers
				for (int block = 0; block < number_of_reading_blocks; ++block)
				{
					for (auto j = 0u; j < size; ++j)  
					{
						std::getline(hess, line);         // "real" line
						linevec = split(line, ' ', true);
						for (int k = 1; k <= columns; ++k)
						{
							hessian[j][k-1+block*columns] = std::stod(linevec[k]) * energy::Hartree_Bohr2Kcal_MolAngSquare;
						}
					}
					std::getline(hess, line); // headline of next block
				}

				// last block
				for (auto j = 0u; j < size; ++j)
				{
					std::getline(hess, line);         // "real" line
					linevec = split(line, ' ', true);
					for (int k = 1; k <= number_of_columns_in_last_block; ++k)
					{
						hessian[j][k-1+number_of_reading_blocks * columns] = std::stod(linevec[k]) * energy::Hartree_Bohr2Kcal_MolAngSquare;
					}
				}

				coords->set_hessian(hessian);  //set hessian
				break;
			}
		}
	}
  hess.close();
}

/*
Energy class functions that need to be overloaded
*/

// Energy function
double energy::interfaces::orca::sysCallInterface::e(void)
{
  integrity = check_structure();
  if (integrity == true)
  {
    write_inputfile(0);
    int res = scon::system_call(Config::get().energy.orca.path + " orca.inp > output_orca.txt");
		if (res == 0) energy = read_output(0);
    else
    {
      if (Config::get().general.verbosity >= 2)
      {
        std::cout << "Orca call (e) return value was not 0. Treating structure as broken.\n";
      }
      integrity = false;
    }
    return energy;
  }
  else return 0;  // energy = 0 if structure contains NaN
}

// Energy+Gradient function
double energy::interfaces::orca::sysCallInterface::g(void)
{
  integrity = check_structure();
  if (integrity == true)
  {
    write_inputfile(1);
		int res = scon::system_call(Config::get().energy.orca.path + " orca.inp > output_orca.txt");
    if (res == 0) energy = read_output(1);
    else
    {
      if (Config::get().general.verbosity >= 2)
      {
        std::cout << "Orca call (g) return value was not 0. Treating structure as broken.\n";
      }
      integrity = false;
    }
    return energy;
  }
  else return 0;  // energy = 0 if structure contains NaN
}

// Hessian function
double energy::interfaces::orca::sysCallInterface::h(void)
{
  integrity = check_structure();
  if (integrity == true)
  {
    write_inputfile(2);
		int res = scon::system_call(Config::get().energy.orca.path + " orca.inp > output_orca.txt");
    if (res == 0) energy = read_output(2);
    else
    {
      if (Config::get().general.verbosity >= 2)
      {
        std::cout << "Orca call (h) return value was not 0. Treating structure as broken.\n";
      }
      integrity = false;
    }
    return energy;
  }
  else return 0;  // energy = 0 if structure contains NaN
}

// Optimization
double energy::interfaces::orca::sysCallInterface::o(void)
{
	integrity = check_structure();
  if (integrity == true)
  {
    write_inputfile(3);
		int res = scon::system_call(Config::get().energy.orca.path + " orca.inp > output_orca.txt");
    if (res == 0) energy = read_output(3);
    else
    {
      if (Config::get().general.verbosity >= 2)
      {
        std::cout << "Orca call (o) return value was not 0. Treating structure as broken.\n";
      }
      integrity = false;
    }
    return energy;
  }
  else return 0;  // energy = 0 if structure contains NaN
}

// Output functions
void energy::interfaces::orca::sysCallInterface::print_E(std::ostream &S) const
{
  S << "Total Energy:      ";
  S << std::right << std::setw(16) << std::fixed << std::setprecision(8) << energy;
}

void energy::interfaces::orca::sysCallInterface::print_E_head(std::ostream &S, bool const endline) const
{
  S << "Partial energies (only for SCF):\n";
  S << std::right << std::setw(24) << "E_nuc";
  S << std::right << std::setw(24) << "E_elec";
  S << std::right << std::setw(24) << "E_one_elec";
	S << std::right << std::setw(24) << "E_two_elec";
  S << std::right << std::setw(24) << "SUM\n";
  if (endline) S << '\n';
}

void energy::interfaces::orca::sysCallInterface::print_E_short(std::ostream &S, bool const endline) const
{
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << nuc_rep;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << elec_en;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << one_elec;
	S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << two_elec;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << energy << '\n';
  if (endline) S << '\n';
}

void energy::interfaces::orca::sysCallInterface::to_stream(std::ostream&) const { }

bool energy::interfaces::orca::sysCallInterface::check_bond_preservation(void) const
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

bool energy::interfaces::orca::sysCallInterface::check_atom_dist(void) const
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
energy::interfaces::orca::sysCallInterface::charges() const
{
  return mulliken_charges;
}

std::vector<coords::Cartesian_Point>
energy::interfaces::orca::sysCallInterface::get_g_ext_chg() const
{
  return grad_ext_charges;
}

