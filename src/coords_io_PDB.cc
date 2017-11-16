/**
CAST 3
Purpose: Reading structures from PDB-files
bullshit atom types are assigned
bonds are created by distance criterion (1.2 times sum of covalent radii)

@author Susanne Sauer
@version 1.0
*/
#pragma once
#include "coords_io.h"
#include "helperfunctions.h"

std::vector<std::string> RESIDUE_NAMES = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "CYX", "CYM", "HID", "HIE", "HIP" };

/**struct that contains information about a residue*/
struct residue
{
  /**name of the residue*/
  std::string res_name;
  /**atoms in the residue*/
  std::vector<coords::Atom> atoms;
  /**is an amino acid terminal or not
  possible values: no (not terminal), C (C terminal), N (N terminal)*/
  std::string terminal;
};

int find_at_sidechain(std::string atom_name, std::string res_name)
{
  if (res_name == "ALA")
  {
    if (atom_name.substr(0, 1) == "C") return 80;
    else if (atom_name.substr(0, 1) == "H") return 85;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "PRO")
  {
    if (atom_name.substr(0, 1) == "H") return 85;
    else if (atom_name == "CD") return 187;
    else if (atom_name.substr(0, 1) == "C") return 81;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "VAL")
  {
    if (atom_name.substr(0, 1) == "H") return 85;
    else if (atom_name == "CB") return 82;
    else if (atom_name.substr(0, 1) == "C") return 80;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "ASP")
  {
    if (atom_name == "HA" || atom_name.substr(0, 2) == "HB") return 85;
    else if (atom_name == "CB") return 81;
    else if (atom_name.substr(0, 1) == "C") return 213;
    else if (atom_name.substr(0, 1) == "O") return 214;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "TRP")
  {
    if (atom_name == "CD1") return 455;
    else if (atom_name.substr(0, 2) == "CD") return 442;
    else if (atom_name == "CB") return 81;
    else if (atom_name.substr(0, 2) == "CG") return 441;
    else if (atom_name.substr(0, 2) == "CE") return 443;
    else if (atom_name.substr(0, 1) == "N") return 444;
    else if (atom_name.substr(0, 2) == "HE") return 445;
    else if (atom_name.substr(0, 2) == "HB") return 85;
    else if (atom_name.substr(0, 1) == "H") return 91;
    else if (atom_name.substr(0, 1) == "C") return 90;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "ARG")
  {
    if (atom_name.substr(0, 2) == "CD") return 250;
    else if (atom_name.substr(0, 2) == "CG") return 251;
    else if (atom_name.substr(0, 2) == "CZ") return 245;
    else if (atom_name.substr(0, 2) == "CB") return 81;
    else if (atom_name.substr(0, 2) == "NH") return 243;
    else if (atom_name.substr(0, 2) == "NE") return 246;
    else if (atom_name.substr(0, 2) == "HE") return 247;
    else if (atom_name.substr(0, 2) == "HH") return 244;
    else if (atom_name.substr(0, 1) == "H") return 85;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "GLU")
  {
    if (atom_name.substr(0, 2) == "OE") return 214;
    else if (atom_name.substr(0, 2) == "CD") return 213;
    else if (atom_name.substr(0, 1) == "C") return 81;
    else if (atom_name.substr(0, 2) == "HB" || atom_name.substr(0, 2) == "HG" || atom_name.substr(0, 2) == "HA") return 85;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "LYS")
  {
    if (atom_name.substr(0, 1) == "N") return 230;
    else if (atom_name.substr(0, 2) == "HZ") return 233;
    else if (atom_name.substr(0, 2) == "CE") return 236;
    else if (atom_name.substr(0, 1) == "H") return 85;
    else if (atom_name.substr(0, 1) == "C") return 81;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "GLY")
  {
    if (atom_name.substr(0, 1) == "H") return 85;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "THR")
  {
    if (atom_name.substr(0, 1) == "O") return 96;
    else if (atom_name.substr(0, 2) == "HG") return 97;
    else if (atom_name.substr(0, 2) == "CB") return 99;
    else if (atom_name.substr(0, 2) == "CG") return 80;
    else if (atom_name.substr(0, 1) == "H") return 85;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "GLN")
  {
    if (atom_name.substr(0, 1) == "O") return 178;
    else if (atom_name.substr(0, 1) == "N") return 179;
    else if (atom_name.substr(0, 2) == "HE") return 182;
    else if (atom_name.substr(0, 2) == "CD") return 177;
    else if (atom_name.substr(0, 1) == "C") return 81;
    else if (atom_name.substr(0, 1) == "H") return 85;
    else
    {
      std::cout << "Strange atom in residue " << res_name << ": " << atom_name << "\nNo atom type assigned.\n";
      return 0;
    }
  }
  else if (res_name == "CYX")
  {
    // TODO: hier geht's weiter
    return 0;
  }
  else if (res_name == "MET")
  {
    if (atom_name == "S") return 144;
    else return 0;
  }
  else return 0;
}

/**function that assigns atom types (oplsaa) to atoms of protein backbone
(they are not suitable for force field calucations)*/
int find_energy_type(std::string atom_name, std::string res_name, std::string terminal)
{
  if (is_in(res_name, RESIDUE_NAMES))
  {
    if (terminal == "no")
    {
      if (atom_name == "N" && res_name != "PRO") return 180;  // amid N 
      else if (atom_name == "N" && res_name == "PRO") return 181;  // amid N 
      else if (atom_name == "H") return 183;  // amid H
      else if (atom_name == "C") return 177;  // amid C
      else if (atom_name == "O") return 178;  // amid O 
      else if (atom_name == "CA" && res_name != "GLY") return 166; // alpha C atom 
      else if (atom_name == "CA" && res_name == "GLY") return 165; // alpha C atom
      else return find_at_sidechain(atom_name,res_name);
    }
    else if (terminal == "C")
    {
      if (atom_name == "N" && res_name != "PRO") return 180;  // amid N 
      else if (atom_name == "N" && res_name == "PRO") return 181;  // amid N 
      else if (atom_name == "H") return 183; // amid H
      else if (atom_name == "C") return 213;  // carbonyl C
      else if (atom_name == "O" || atom_name == "OXT") return 214;  // C-terminal O
      else if (atom_name == "CA" && res_name != "GLY") return 166; // alpha C atom 
      else if (atom_name == "CA" && res_name == "GLY") return 165; // alpha C atom 
      else return find_at_sidechain(atom_name, res_name);
    }
    else if (terminal == "N")
    {
      if (atom_name == "C") return 177;  // amid C
      else if (atom_name == "O") return 178;  // amid O 
      else if (atom_name == "CA" && res_name != "GLY" && res_name != "PRO") return 236; // alpha C atom 
      else if (atom_name == "CA" && res_name == "GLY") return 235; // alpha C atom
      else if (atom_name == "CA" && res_name == "PRO") return 237; // alpha C atom
      else if (atom_name == "N" && res_name != "PRO") return 230; // terminal N
      else if (atom_name == "N" && res_name == "PRO") return 252; // terminal N
      else if (atom_name.substr(0, 1) == "H" && isdigit(atom_name.substr(1, 1)) && res_name == "PRO")
      {
        std::cout << "N terminal Prolin not implemented\n";
        std::cout << "no atom type is assigned for H atom\n";
        return 0;
      }
      else if (atom_name.substr(0, 1) == "H" && isdigit(atom_name.substr(1, 1))) return 233; // terminal H(N)      
      else return find_at_sidechain(atom_name, res_name);
    }
  }
  else return 0;
  
}

/**finds element symbol and energy type
@param atom_name: atom name from pdb file
@param res_name: residue name from pdb file
returns element symbol*/
std::string find_element_symbol(std::string atom_name, std::string res_name)
{
  std::string element;
  if (is_in(res_name, RESIDUE_NAMES))  // protein
  {
    element = atom_name.substr(0, 1);
  }
  else if (atom_name.substr(atom_name.size() - 1, 1) == "+" || atom_name.substr(atom_name.size() - 1, 1) == "-")  // ion
  {
    element = atom_name.substr(0, atom_name.size() - 1);
  }
  else if (res_name == "WAT" || res_name == "LIG")  // water or ligand
  {
    element = "";
    for (auto s : atom_name) // every sign in atom name
    {
      if (!isdigit(s)) element += s;
    }
  }
  else  // unknown residue
  {
    std::cout << "WARNING: unknown residue: " << res_name << ". Guessing element symbol.\n";
    element = "";
    for (auto s : atom_name) // every sign in atom name
    {
      if (!isdigit(s)) element += s;
    }
  }
  return element;
}

/**function that reads the structure
@ param file: name of the pdb-file
@ return: Coordinates object that is created out of file*/
coords::Coordinates coords::input::formats::pdb::read(std::string file)
{
	if ((Config::get().general.energy_interface == config::interface_types::T::AMBER) ||
		(Config::get().general.energy_interface == config::interface_types::T::AMOEBA) ||
		(Config::get().general.energy_interface == config::interface_types::T::CHARMM22) ||
		(Config::get().general.energy_interface == config::interface_types::T::OPLSAA))
	{
		std::cout << "ERROR: It is not possible to use PDB files with a forcefield interface because no atom types are assigned!\n";
		if (Config::get().general.task == config::tasks::WRITE_TINKER || Config::get().general.task == config::tasks::CUT_RESIDUES)
		{
			std::cout << "Yes, I know you just want to write a tinkerstructure and you don't need any energies. But it doesn't work like this. So just use GAUSSIAN or MOPAC as energy interface and all will be fine (even if you don't have access to any of these programmes).\n";
		}
		std::exit(0);
	}

	Coordinates coord_object;
	std::ifstream config_file_stream(file.c_str(), std::ios_base::in);  // read file to ifstream

    std::string line, element;  // some variables
    Representation_3D positions;
    std::vector<residue> residues;

    int N = 0; // number of atoms
	while (std::getline(config_file_stream, line))
	{
		if (line.substr(0, 4) == "ATOM")
		{
			std::string atom_name = line.substr(12, 4);  // read atom name and remove spaces
            atom_name.erase(remove_if(atom_name.begin(), atom_name.end(), isspace), atom_name.end());

			std::string res_name = line.substr(17, 3);  // read residue name
            std::string res_number = line.substr(22, 4);  // read residue id

            // find element symbol
            element = find_element_symbol(atom_name, res_name);

            // create atom
            Atom current(element);
            current.set_residue(res_name);  
            current.set_res_id(std::stoi(res_number));
            current.set_pdb_atom_name(atom_name);
            atoms.add(current);

            // create residues
            if (std::stoi(res_number) > residues.size())
            {
              residue new_res;
              new_res.res_name = res_name;
              new_res.atoms.push_back(current);
              residues.push_back(new_res);
            }
            else residues[std::stoi(res_number) - 1].atoms.push_back(current);

            // read and create positions
            std::string x = line.substr(30, 8);
            std::string y = line.substr(38, 8);
            std::string z = line.substr(46, 8);

            position.x() = std::stof(x);
            position.y() = std::stof(y);
            position.z() = std::stof(z);
            positions.push_back(position);

            N += 1; // count atoms
		}
	}

	input_ensemble.push_back(positions);  // fill the positions into PES_Point
	coords::PES_Point pes(input_ensemble[0u]);

	// loop over all atompairs and bind them if they fulfill distance criterion 
	// i.e. the distance is smaller than 1.2 * sum of covalent radiuses
	for (unsigned i = 0; i<N; i++)
	{
		for (unsigned j = 0; j<i; j++)
		{
			double d = dist(positions[i], positions[j]);
			double d_max = 1.2*(atoms.atom(i).cov_radius() + atoms.atom(j).cov_radius());
			if (d < d_max)
			{
              // do not bond ions
              std::string res1 = atoms.atom(i).get_residue();
              std::string res2 = atoms.atom(j).get_residue();
              if (res1.substr(res1.size() - 1, 1) == "+" || res1.substr(res1.size() - 1, 1) == "-") {}
              else if (res2.substr(res2.size() - 1, 1) == "+" || res2.substr(res2.size() - 1, 1) == "-") {}

              //save bonds
              else
              {   
                atoms.atom(i).bind_to(j);
                atoms.atom(j).bind_to(i);
              }
			}
		}
	}

    // determine if protein residues are terminal or not
    for (auto &r: residues)
    {
      r.terminal = "no";
      if (is_in(r.res_name, RESIDUE_NAMES))
      {
        for (auto a : r.atoms)
        {
          if (a.get_pdb_atom_name() == "OXT") r.terminal = "C";
          else if (a.get_pdb_atom_name() == "H1") r.terminal = "N";
        }
      }
    }

    // set energy type of every atom
    for (auto &a : atoms)
    {
      if (Config::get().general.verbosity > 3)
      {
        std::cout << "Atom " << a.symbol() << " belongs to residue " << a.get_residue() << " that is terminal: " << residues[a.get_res_id() - 1].terminal << "\n";
      }
      a.set_energy_type(find_energy_type(a.get_pdb_atom_name(), a.get_residue(), residues[a.get_res_id() - 1].terminal));
    }

	coord_object.init_swap_in(atoms, pes);  // fill atoms and positions into coord_object

	for (auto & p : input_ensemble)  // do some important stuff (see coords_io_AMBER.cc)
	{
		p.gradient.cartesian.resize(p.structure.cartesian.size());
		coord_object.set_xyz(p.structure.cartesian);
		coord_object.to_internal_light();
		p = coord_object.pes();
	}

    return coord_object;
}