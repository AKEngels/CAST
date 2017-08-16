#pragma once
#include "coords.h"
#include "scon_utility.h"

namespace periodicsHelperfunctions
{

  bool isInsideCenteredBox(scon::c3<coords::float_type> const& point, scon::c3<coords::float_type> const& halfbox)
  {
    if (point.x() > halfbox.x())
    {
      return false;
    }
    else if (point.x() < -halfbox.x())
    {
      return false;
    }
    else if (point.y() > halfbox.y())
    {
      return false;
    }
    else if (point.y() < -halfbox.y())
    {
      return false;
    }
    else if (point.z() > halfbox.z())
    {
      return false;
    }
    else if (point.z() < -halfbox.z())
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  coords::Coordinates periodicCutout(coords::Coordinates const& inputStructure)
  {
    coords::Cartesian_Point const halfbox(Config::get().periodics.pb_box / 2.0);
    coords::Atoms const& atoms = inputStructure.atoms();
    coords::Atoms truncatedAtoms;
    coords::Representation_3D positions;
    std::vector<std::size_t> new_index_of_atom(inputStructure.atoms().size(), 0u);
    /**
     * Check for atom criterion
     */
    if (Config::get().periodics.criterion == 0u)
    {
      for (unsigned int i = 0u; i < inputStructure.atoms().molecules().size(); i++)
      {
        bool checkIfAllPartsOfMoleculeAreInsideBox = true;
        for (std::vector<coords::Atom>::size_type j = 0u; j < atoms.molecule(i).size(); ++j)
        {

          coords::cartesian_type const& currentAtom = inputStructure.xyz(atoms.atomOfMolecule(i, j));
          if (!isInsideCenteredBox(currentAtom, halfbox - scon::c3<coords::float_type>(Config::get().periodics.cutout_distance_to_box)))
          {
            checkIfAllPartsOfMoleculeAreInsideBox = false;
            break;
          }
        }

        if (checkIfAllPartsOfMoleculeAreInsideBox)
        {
          for (std::vector<coords::Atom>::size_type j = 0u; j < atoms.molecule(i).size(); ++j)
          {
            truncatedAtoms.add(inputStructure.atoms(atoms.atomOfMolecule(i, j)));
            positions.push_back(inputStructure.xyz(atoms.atomOfMolecule(i, j)));
            new_index_of_atom.at(atoms.atomOfMolecule(i, j)) = truncatedAtoms.size() - 1u;
          }
          unsigned int helperIterator = 0u;
          for (std::vector<coords::Atom>::size_type j = truncatedAtoms.size() - atoms.molecule(i).size(); j < truncatedAtoms.size(); ++j, helperIterator++)
          {
            for (auto& bonding_partner : inputStructure.atoms(atoms.atomOfMolecule(i, helperIterator)).bonds())
            {
              truncatedAtoms.atom(j).detach_from(bonding_partner);
              truncatedAtoms.atom(j).bind_to(new_index_of_atom[bonding_partner]);
            }
          }
        }
      }
    }
    /**
     * Check for center of mass criterion
     */
    else if (Config::get().periodics.criterion == 1u)
    {
      for (unsigned int i = 0u; i < inputStructure.atoms().molecules().size(); i++)
      {
        coords::Cartesian_Point COM = coords::Cartesian_Point();
        coords::float_type M = coords::float_type();
        for (std::vector<coords::Atom>::size_type j = 0u; j < atoms.molecule(i).size(); ++j)
        {
          coords::float_type mass(inputStructure.atoms(atoms.atomOfMolecule(i, j)).mass());
          M += mass;
          COM += inputStructure.xyz(atoms.atomOfMolecule(i, j)) * mass;
        }
        COM /= M;

        if (isInsideCenteredBox(COM, halfbox - scon::c3<coords::float_type>(Config::get().periodics.cutout_distance_to_box)))
        {
          for (std::vector<coords::Atom>::size_type j = 0u; j < atoms.molecule(i).size(); ++j)
          {
            truncatedAtoms.add(inputStructure.atoms(atoms.atomOfMolecule(i, j)));
            positions.push_back(inputStructure.xyz(atoms.atomOfMolecule(i, j)));
            new_index_of_atom.at(atoms.atomOfMolecule(i, j)) = truncatedAtoms.size() - 1u;
          }
          unsigned int helperIterator = 0u;
          for (std::vector<coords::Atom>::size_type j = truncatedAtoms.size() - atoms.molecule(i).size(); j < truncatedAtoms.size(); ++j, helperIterator++)
          {
            for (auto& bonding_partner : inputStructure.atoms(atoms.atomOfMolecule(i, helperIterator)).bonds())
            {
              truncatedAtoms.atom(j).detach_from(bonding_partner);
              truncatedAtoms.atom(j).bind_to(new_index_of_atom[bonding_partner]);
            }
          }
        }
      }
    }
    else
      throw std::logic_error("Unknown criterion in periodics cutout. Aborting.");

    coords::Coordinates newCoords;
    coords::PES_Point x(positions);

    newCoords.init_in(truncatedAtoms, x, true);
    std::stringstream temporaryStringstream;
    return newCoords;
  }

  /**delete molecules from structure
  @param inputStructure: coords object containing structure
  @param indices: vector containing indices of molecules to delete*/
  coords::Coordinates delete_molecules(coords::Coordinates const& inputStructure, std::vector<std::size_t> indices)
  {
    coords::Atoms const& atoms = inputStructure.atoms();
    coords::Atoms truncatedAtoms;
    coords::Representation_3D positions;
    std::vector<std::size_t> new_index_of_atom(inputStructure.atoms().size(), 0u);

    for (std::size_t i = 0u; i < inputStructure.atoms().molecules().size(); i++)
    {
      if (std::find(indices.begin(), indices.end(), i) == indices.end())
      {
        for (std::vector<coords::Atom>::size_type k = 0u; k < atoms.molecule(i).size(); ++k)
        {
          truncatedAtoms.add(inputStructure.atoms(atoms.atomOfMolecule(i, k)));
          positions.push_back(inputStructure.xyz(atoms.atomOfMolecule(i, k)));
          new_index_of_atom.at(atoms.atomOfMolecule(i, k)) = truncatedAtoms.size() - 1u;
        }//k
        unsigned int helperIterator = 0u;
        for (std::vector<coords::Atom>::size_type k = truncatedAtoms.size() - atoms.molecule(i).size(); k < truncatedAtoms.size(); ++k, helperIterator++)
        {
          for (auto& bonding_partner : inputStructure.atoms(atoms.atomOfMolecule(i, helperIterator)).bonds())
          {
            truncatedAtoms.atom(k).detach_from(bonding_partner);
            truncatedAtoms.atom(k).bind_to(new_index_of_atom[bonding_partner]);
          }
        }//k
      }
    }//i
    coords::Coordinates newCoords;
    coords::PES_Point x(positions);

    newCoords.init_in(truncatedAtoms, x, true);
    return newCoords;
  }

  coords::Coordinates delete_random_molecules(coords::Coordinates const& coords, std::size_t del_amount)
  {
    coords::Coordinates new_coords;
    std::size_t N_mol = coords.molecules().size();
    std::vector<std::size_t> del_indices;                //Vector for the indices of the molecules to delete                 

    std::default_random_engine rn_generator;
    std::uniform_int_distribution<int> uni_int_distr(0, N_mol - 1);
    std::size_t tmp_rn;

    while (del_indices.size() < del_amount) //Loop to get secified amount of different random numbers
    {
      tmp_rn = uni_int_distr(rn_generator);
      if (del_indices.empty())                             //Check if vector for indicesis empty
      {
        del_indices.push_back(tmp_rn);
      }
      else
      {
        if (std::find(del_indices.begin(), del_indices.end(), tmp_rn) == del_indices.end()) //Chck if generated random number is already part of del_indices
        {
          del_indices.push_back(tmp_rn);
        }
      }
    }//end while(del_indices)

    new_coords = periodicsHelperfunctions::delete_molecules(coords, del_indices);

    return new_coords;
  }

}