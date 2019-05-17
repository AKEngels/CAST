/**
CAST 3
Purpose: Reading structures from XYZ-files
no atom types are assigned
bonds are created by distance criterion (1.2 times sum of covalent radii)

@author Susanne Sauer
@version 1.0
*/
#include "coords_io.h"
#include "helperfunctions.h"


std::vector<coords::input::formats::xyz::AminoAcid> coords::input::formats::xyz::ETfromAA::get_aminoacids()
{
  std::vector<AminoAcid> amino_acids;

  for (auto i{ 0u }; i < atoms.size(); ++i)     // for all atoms
  {
    auto a = atoms.atom(i);
    if (a.symbol() == "O" && got_it[i] == false)      // O atom that is not in an amino acid yet
    {
      if (a.bonds().size() == 1 && atoms.atom(a.bonds()[0]).symbol() == "C")  // if only bound to C (=> carbonyle)
      {
        auto j = a.bonds()[0];
        auto b = atoms.atom(j);
        if (b.bonds().size() == 3 && got_it[j] == false)                      // if carbonyle C is not in an amino acid yet (could be for a terminal amino acid)
        {
          std::vector<std::string> symbolvec_b = get_bonding_symbols(b, atoms);
          auto index = find_index("C", symbolvec_b);
          if (index < std::numeric_limits<int>::max())                        // if other C atom is bound to carbonyle C (=> C alpha)
          {
            auto k = b.bonds()[index];
            auto c = atoms.atom(k);
            if (c.bonds().size() == 4)
            {
              auto symbolvec_c = get_bonding_symbols(c, atoms);
              index = find_index("N", symbolvec_c);
              if (index < std::numeric_limits<int>::max())                    // if this C is bound to N 
              {
                auto l = c.bonds()[index];
                auto d = atoms.atom(l);

                // set all 4 atoms to got_it
                got_it[i] = true;      
                got_it[j] = true;
                got_it[k] = true;
                got_it[l] = true;

                // determine terminal state
                auto terminal = terminalState::no;
                if (count_element("O", symbolvec_b) == 2) terminal = terminalState::C;
                if (count_element("H", get_bonding_symbols(d, atoms)) > 1) terminalState::N;

                // create amino acid and add it to vector
                AminoAcid as({ i, j, k, l }, terminal);
                amino_acids.emplace_back(as);
              }
            }
          }
        }
      }
    }
  }
  return amino_acids;
}

void coords::input::formats::xyz::ETfromAA::complete_atoms_of_aminoacids(std::vector<AminoAcid>& amino_acids)
{
  // TODO
}

void coords::input::formats::xyz::ETfromAA::find_energy_types()
{
  auto amino_acids = get_aminoacids();
  complete_atoms_of_aminoacids(amino_acids);
}

/**function that reads the structure
@ param file: name of the xyz-file
@ return: Coordinates object that is created out of file*/
coords::Coordinates coords::input::formats::xyz::read(std::string file)
{
  if ((Config::get().general.energy_interface == config::interface_types::T::AMBER) ||
    (Config::get().general.energy_interface == config::interface_types::T::AMOEBA) ||
    (Config::get().general.energy_interface == config::interface_types::T::CHARMM22) ||
    (Config::get().general.energy_interface == config::interface_types::T::OPLSAA))
  {
    std::cout<<"ERROR: It is not possible to use XYZ files with a forcefield interface because no atom types are assigned!\n";
    if (Config::get().general.task == config::tasks::WRITE_TINKER)
    {
      std::cout<<"Yes, I know you just want to write a tinkerstructure and you don't need any energies. But it doesn't work like this. So just use GAUSSIAN or MOPAC as energy interface and all will be fine (even if you don't have access to any of these programmes).\n";
    }
    std::exit(0);
  }

  Coordinates coord_object;
  std::ifstream config_file_stream(file.c_str(), std::ios_base::in);  // read file to ifstream
    
  std::string line, element;  // a few variables
  double x,y,z;
  Representation_3D positions;

  std::getline(config_file_stream, line);  // first line: number of atoms
  double N = std::stoi(line); 
    
  std::getline(config_file_stream, line);  // discard second line (comment)
    
  while (config_file_stream >> element >> x >> y >> z)  // for every line
  {
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

  input_ensemble.push_back(positions);  // fill the positions into PES_Point
  coords::PES_Point pes(input_ensemble[0u]);

  // loop over all atompairs and bind them if they fulfill distance criterion 
  // i.e. the distance is smaller than 1.2 * sum of covalent radiuses
  for (unsigned i=0; i<N; i++)
  {
    for(unsigned j=0; j<i; j++)
    {
      double d = dist(positions[i], positions[j]);
      double d_max = 1.2*(atoms.atom(i).cov_radius() + atoms.atom(j).cov_radius());
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

  if (Config::get().stuff.et_from_aa)   // try to create atomtypes
  {
    ETfromAA energytype_creator(atoms);
    energytype_creator.find_energy_types();
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