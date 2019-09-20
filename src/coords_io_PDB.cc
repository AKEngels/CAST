/**
CAST 3
Purpose: Reading structures from PDB-files and writing structures into PDB files

PDB reading:
  no atom types are assigned
  bonds are created by distance criterion (1.2 times sum of covalent radii)

@author Susanne Sauer
@version 1.0
*/
#include "coords_io.h"
#include "helperfunctions.h"

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
    std::cout << "ERROR: It is not possible to use PDB files with that interface because wrong atom types are assigned!\n";
    if (Config::get().general.task == config::tasks::WRITE_TINKER || Config::get().general.task == config::tasks::WRITE_GAUSSVIEW)
    {
      std::cout << "Yes, I know you just want to write a tinkerstructure and you don't need any energies. But it doesn't work like this. So just use GAUSSIAN or MOPAC as energy interface and all will be fine (even if you don't have access to any of these programmes).\n";
    }
    std::exit(-1);
  }

  Coordinates coord_object;
  std::ifstream config_file_stream(file.c_str(), std::ios_base::in);  // read file to ifstream

  std::string line, element;  // some variables
  Representation_3D positions;

  std::size_t N = 0; // number of atoms
  while (std::getline(config_file_stream, line))
  {
    if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM")
    {
      std::string alternate_location = line.substr(16, 1);                   // if there are several alternate locations
      if (alternate_location != " " && alternate_location != "A") continue;  // always choose the first one

      try {
        element = remove_spaces(line.substr(76, 2)); // read element directly from PDB file
      }
      catch (...) {
        element = remove_spaces(line.substr(12, 4));  // take atom name as element symbol (like in gaussian output)
      }

      // create atom
      Atom current(element);
      current.set_energy_type(0);
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
  // ions are not bonded to anything
  for (unsigned i = 0; i < N; i++)
  {
    for (unsigned j = 0; j < i; j++)
    {
      double d = dist(positions[i], positions[j]);
      double d_max = 1.2 * (atoms.atom(i).cov_radius() + atoms.atom(j).cov_radius());
      if (d < d_max)
      {
        if (atoms.atom(i).symbol() == "Na" || atoms.atom(j).symbol() == "Na")
        {                           // Na ions often have a small distance to their neighbors but no bonds
          std::cout << "creating no bond between atoms " << i + 1 << " and " << j + 1 << " because one of the atoms is a Na\n";
        }
        else
        {
          atoms.atom(i).bind_to(j);
          atoms.atom(j).bind_to(i);
        }
      }
    }
  }

  if (!Config::get().coords.fixed.empty())    // fix atoms
  {
    for (auto fix : Config::get().coords.fixed)
    {
      if (fix < atoms.size()) atoms.atom(fix).fix(true);
    }
  }
  if (Config::get().coords.fix_sphere.use)  // fix everything outside of a given sphere
  {
    for (auto i{ 0u }; i < atoms.size(); ++i)
    {
      double d = dist(positions[i], positions[Config::get().coords.fix_sphere.central_atom]);
      if (d > Config::get().coords.fix_sphere.radius) atoms.atom(i).fix(true);
    }
  }

  coord_object.init_swap_in(atoms, pes);  // fill atoms and positions into coord_object

  for (auto& p : input_ensemble)  // do some important stuff (see coords_io_AMBER.cc)
  {
    p.gradient.cartesian.resize(p.structure.cartesian.size());
    coord_object.set_xyz(p.structure.cartesian);
    coord_object.to_internal_light();
    p = coord_object.pes();
  }

  return coord_object;
}

void coords::output::formats::pdb::preparation()
{
  auto atoms = ref.atoms();
  AtomtypeFinder atf(atoms);
  auto amino_acids = atf.get_aminoacids();   // find aminoacids

  auto res_counter{ 0u }; // count residues

  for (auto& as : amino_acids)
  {
    as.determine_aminoacid(atoms);       // determine which aminoacid

    if (as.get_res_name() != "XXX")
    {
      res_counter++;
      as.correct_residue_names(atoms);    // determine_aminoacid() not always assigns correct three-letter code

      for (auto j{ 0u }; j < as.get_indices().size(); ++j)   // create PDB atoms from aminoacids
      {
        auto i = as.get_indices()[j];
        pdb_atoms[i].record_name = "ATOM";
        if (j == 2) pdb_atoms[i].atom_name = "CA";   // C_alpha atoms need to have this name, otherwise protein chain will not be recognized in Pymol
        else pdb_atoms[i].atom_name = atoms.atom(i).symbol();
        pdb_atoms[i].residue_name = as.get_res_name();
        pdb_atoms[i].residue_number = res_counter;
        pdb_atoms[i].x = ref.xyz(i).x();
        pdb_atoms[i].y = ref.xyz(i).y();
        pdb_atoms[i].z = ref.xyz(i).z();
        pdb_atoms[i].symbol = atoms.atom(i).symbol();
      }
    }
  }

  for (auto m : ref.molecules()) {
    for (auto a : m) {
      if (atf.recognized_atom(a) == false)   // for every molecule that contains atoms that are not recognized yet
      {
        res_counter++;
        set_pdb_atoms_of_molecule(m, res_counter);
      }
    }
  }
}

void coords::output::formats::pdb::set_pdb_atoms_of_molecule(Container<std::size_t> const& molecule, int residue_counter)
{
  for (auto a : molecule) {
    if (pdb_atoms[a].symbol == "X")   // for every atom is molecule that has not been recognized yet
    {
      pdb_atoms[a].record_name = "HETATM";
      pdb_atoms[a].atom_name = ref.atoms().atom(a).symbol();
      pdb_atoms[a].residue_name = ref.molecule_name(molecule);
      pdb_atoms[a].residue_number = residue_counter;
      pdb_atoms[a].x = ref.xyz(a).x();
      pdb_atoms[a].y = ref.xyz(a).y();
      pdb_atoms[a].z = ref.xyz(a).z();
      pdb_atoms[a].symbol = ref.atoms().atom(a).symbol();
    }
  }
}

void coords::output::formats::pdb::to_stream(std::ostream& os) const
{
  auto atom_counter{ 0u };
  for (auto a : pdb_atoms)  // see: ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf
  {
    atom_counter++;
    os << std::left << std::setw(6) << a.record_name <<
      std::right << std::setw(5) << atom_counter << " " <<
      std::left << std::setw(4) << a.atom_name <<
      std::right << std::setw(4) << a.residue_name <<
      std::setw(6) << a.residue_number << "    " <<
      std::setw(8) << std::fixed << std::setprecision(3) << a.x << std::setw(8) << std::setprecision(3) << a.y << std::setw(8) << std::setprecision(3) << a.z <<
      std::setw(24) << std::right << a.symbol << "\n";
  }
}