/**
CAST 3
Purpose: Reading structures from XYZ-files
if desired oplsaa atom types are assigned to proteins, water and other easy stuff
bonds are created by distance criterion (1.2 times sum of covalent radii)

@author Susanne Sauer
@version 1.0
*/
#include "coords_io.h"
#include "helperfunctions.h"
#include "find_as.h"

/**function that reads the structure
@ param file: name of the xyz-file
@ return: Coordinates object that is created out of file*/
coords::Coordinates coords::input::formats::xyz::read(std::string file)
{
  if ((Config::get().general.energy_interface == config::interface_types::T::AMOEBA) ||
    (Config::get().general.energy_interface == config::interface_types::T::FORCEFIELD))
  {
    std::cout << "ERROR: It is not possible to use XYZ files with a forcefield interface because no atom types are assigned!\n";
    if (Config::get().general.task == config::tasks::WRITE_TINKER)
    {
      std::cout << "Yes, I know you just want to write a tinkerstructure and you don't need any energies. But it doesn't work like this. So just use GAUSSIAN or MOPAC as energy interface and all will be fine (even if you don't have access to any of these programmes).\n";
    }
    std::exit(-1);
  }

  Coordinates coord_object;
  std::ifstream config_file_stream(file.c_str(), std::ios_base::in);  // read file to ifstream

  std::string line, element;  // a few variables
  double x, y, z;
  Representation_3D positions;

  std::getline(config_file_stream, line);  // first line: number of atoms
  const std::size_t N = std::stoi(line);

  std::getline(config_file_stream, line);  // discard second line (comment)

  positions = Representation_3D();
  while (config_file_stream.good())  // for every line
  {
    // read line into element and coordinates
    std::getline(config_file_stream, line);
    std::stringstream linestream(line);
    if (check_if_integer(line) || only_whitespace(line)) break;  // either beginning of new structure or empty line
    linestream >> element >> x >> y >> z;

    // create atom
    Atom current(element);
    current.set_energy_type(0);  // because of this no forcefield interfaces are available with XYZ inputfile
    atoms.add(current);

    // create position
    position.x() = x;
    position.y() = y;
    position.z() = z;
    positions.push_back(position);
  }
  input_ensemble.push_back(positions);
  positions.clear();

  while (config_file_stream.good())     // if there is still file left there are probably some more structures in it
  {
    if (check_if_integer(line))   // first line = number of atoms
    {
      std::size_t number_of_atoms = std::stoi(line);
      if (N != number_of_atoms) throw std::runtime_error("There is a structure in your inputfile that has not the same number of atoms as the first.");
      std::getline(config_file_stream, line);    // throwing away second (empty) line of structure

      while (config_file_stream.good())  // for every line
      {
        // read line into element (thrown away) and coordinates
        std::getline(config_file_stream, line);
        std::stringstream linestream(line);
        if (check_if_integer(line) || only_whitespace(line)) break;  // either beginning of new structure or empty line
        linestream >> element >> x >> y >> z;

        // create position
        position.x() = x;
        position.y() = y;
        position.z() = z;
        positions.push_back(position);
      }
      input_ensemble.push_back(positions);
      positions.clear();
    }
    else std::getline(config_file_stream, line);
  }

  positions = input_ensemble.at(0u).structure.cartesian;
  coords::PES_Point pes(input_ensemble[0u]);

  // loop over all atompairs and bind them if they fulfill distance criterion 
  // i.e. the distance is smaller than 1.2 * sum of covalent radiuses
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

  if (Config::get().stuff.xyz_atomtypes)   // try to create atomtypes
  {
    AtomtypeFinder energytype_creator(atoms);
    energytype_creator.find_energy_types();
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