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

/**function that assigns atom types (oplsaa) to atoms of protein backbone
(they are not suitable for force field calucations)*/
int find_energy_type(std::string atom_name)
{
  if (atom_name == "N") return 180;  // amid N 
  else if (atom_name == "H") return 182;  // amid H 
  else if (atom_name == "C") return 177;  // amid C
  else if (atom_name == "O") return 178;  // amid O (also one of the C-terminal Os for now)
  else if (atom_name == "CA") return 166; // alpha C atom (not differentiated if terminal or not for now)
  else if (atom_name == "OXT") return 214; // C-terminal O
  else if (atom_name.substr(0, 1) == "H" && isdigit(atom_name.substr(1, 1))) return 233; // terminal N
  else if (atom_name.substr(0, 1) == "N" && isdigit(atom_name.substr(1, 1))) return 230; // N-terminal H
  else return 0;
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
		std::cout << "ERROR: It is not possible to use XYZ files with a forcefield interface because no atom types are assigned!\n";
		if (Config::get().general.task == config::tasks::WRITE_TINKER)
		{
			std::cout << "Yes, I know you just want to write a tinkerstructure and you don't need any energies. But it doesn't work like this. So just use GAUSSIAN or MOPAC as energy interface and all will be fine (even if you don't have access to any of these programmes).\n";
		}
		std::exit(0);
	}

	Coordinates coord_object;
	std::ifstream config_file_stream(file.c_str(), std::ios_base::in);  // read file to ifstream

    std::string line, element;  // some variables
    Representation_3D positions;

    int N = 0; // number of atoms
	while (std::getline(config_file_stream, line))
	{
		if (line.substr(0, 4) == "ATOM")
		{
			std::string atom_name = line.substr(12, 4);  // read atom name and remove spaces
            atom_name.erase(remove_if(atom_name.begin(), atom_name.end(), isspace), atom_name.end());

			std::string res_name = line.substr(17, 3);  // read residue name
            std::string res_number = line.substr(22, 4);  // read residue id

            int et = 0;
            
            // find element symbol
            if (is_in(res_name, RESIDUE_NAMES))
            {
              element = atom_name.substr(0, 1);
              et = find_energy_type(atom_name);
            }
            else if (atom_name.substr(atom_name.size() - 1,1) == "+" || atom_name.substr(atom_name.size() - 1,1) == "-")
            {
              element = atom_name.substr(0, atom_name.size() - 1);
            }
            else
            {
              element = "";
              for (auto s : atom_name) // every sign in atom name
              {
                if (!isdigit(s)) element += s; 
              }
            }

            // create atom
            Atom current(element);
            current.set_energy_type(et);  // because of this no forcefield interfaces are available with PDB inputfile
            current.set_residue(res_name);  
            current.set_res_id(std::stoi(res_number));
            atoms.add(current);

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