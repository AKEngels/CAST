#pragma once
#include <vector>
#include "coords.h"
#include "matop.h"

namespace monomerManipulation
{

  /** create coords object for a single molecule from a coords object
  @param inputStructure: coords object containing structure
  @param index: index of molecule a coords_object is to be created for*/
  coords::Coordinates monomer_structure(coords::Coordinates const& inputStructure, std::size_t index)
  {
    coords::Atoms truncatedAtoms;
    coords::Representation_3D positions;
    coords::Atoms const& atoms = inputStructure.atoms();
    std::vector<std::size_t> molecule(atoms.molecule(index));
    std::vector<std::size_t> new_index_of_atom(inputStructure.atoms().size(), 0u);


    if (index >= atoms.molecules().size())//check if a valid molecule is requested
    {
      throw std::logic_error("Index too big, strucutre contains not so many monomers.");
    }
    else
    {
      for (std::vector<coords::Atom>::size_type j = 0u; j < atoms.molecule(index).size(); ++j)
      {
        truncatedAtoms.add(inputStructure.atoms(atoms.atomOfMolecule(index, j)));
        positions.push_back(inputStructure.xyz(atoms.atomOfMolecule(index, j)));
        new_index_of_atom.at(molecule.at(j)) = truncatedAtoms.size() - 1u;
      }//j
      unsigned int helperIterator = 0u;
      for (std::vector<coords::Atom>::size_type k = truncatedAtoms.size() - atoms.molecule(index).size(); k < truncatedAtoms.size(); ++k, helperIterator++)//ashures only atoms of molecule considered at the moment are used
      {
        for (auto& bonding_partner : inputStructure.atoms(atoms.atomOfMolecule(index, helperIterator)).bonds())
        {
          truncatedAtoms.atom(k).detach_from(bonding_partner);
          truncatedAtoms.atom(k).bind_to(new_index_of_atom[bonding_partner]);
        }
      }//k

      coords::Coordinates newCoords;
      coords::PES_Point x(positions);

      newCoords.init_in(truncatedAtoms, x, true);
      return newCoords;
    }
  }

  /** Replace monomers in structure with another molecule structure
  @param inputStructure: coords object containing structure
  @param reference: structure of monomer for replacement*/
  coords::Coordinates replaceMonomers(coords::Coordinates const& inputCoords, coords::Coordinates const& reference1, coords::Coordinates const& reference2, std::size_t threshold_index)
  {
    std::size_t N = inputCoords.molecules().size(); //Number of Molecules in system
    std::vector<coords::Cartesian_Point> com;

    for (std::size_t i = 0u; i < N; i++)//i iterates over all molecules
    {
      com.push_back(inputCoords.center_of_mass_mol(i));
    }

    coords::Coordinates monomer;
    coords::Coordinates reference1_cog_aligned = align::centerOfGeometryAligned(reference1);
    coords::Coordinates reference2_cog_aligned = align::centerOfGeometryAligned(reference2);
    coords::Coordinates newCoords;
    coords::Atoms truncatedAtoms;
    coords::Representation_3D positions;
    std::vector<std::size_t> new_index_of_atom(inputCoords.atoms().size(), 0u);
    std::size_t new_index_of_atom_offset = 0;

    for (std::size_t i = 0u; i < N; i++)//i iterates over all molecules
    {
      monomer = align::centerOfGeometryAligned(monomer_structure(inputCoords, i));

      if (i < threshold_index)
      {
        newCoords = align::kabschAligned(reference1_cog_aligned, monomer);
      }
      else
      {
        newCoords = align::kabschAligned(reference2_cog_aligned, monomer);
      }

      for (std::vector<coords::Atom>::size_type j = 0u; j < newCoords.size(); ++j)
      {
        truncatedAtoms.add(newCoords.atoms(j));
        positions.push_back(newCoords.xyz(j) + com[i]);//positions of alinged molecule atoms are shifted by the position of the com of the original molecule
        new_index_of_atom.at(j + new_index_of_atom_offset) = truncatedAtoms.size() - 1u;
      }//j

      std::size_t helperIterator = 0u;

      for (std::vector<coords::Atom>::size_type j = truncatedAtoms.size() - newCoords.size(); j < truncatedAtoms.size(); ++j, helperIterator++)
      {
        std::vector<std::size_t> new_bond_indices;
        for (std::size_t k = 0u; 0 < truncatedAtoms.atom(j).bonds().size();)
        {
          auto detach = truncatedAtoms.atom(j).bonds(k);
          new_bond_indices.push_back(detach);
          truncatedAtoms.atom(j).detach_from(detach);
        }//k
        for (std::size_t l = 0; l < new_bond_indices.size(); l++)
        {
          truncatedAtoms.atom(j).bind_to(new_index_of_atom[new_bond_indices[l] + new_index_of_atom_offset]);
        }//l
      }//j

      new_index_of_atom_offset += newCoords.size();

    }//i

    coords::Coordinates returnCoords;
    coords::PES_Point x(positions);

    returnCoords.init_in(truncatedAtoms, x, true);
    return returnCoords;
  }

}