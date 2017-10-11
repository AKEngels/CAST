#include "Center.h"
#include<iostream>
#include<iomanip>
#include "coords.h"





void center(coords::Coordinates coords)
{
  std::size_t N = coords.molecules().size(); //Number of Molecules in system
  std::vector <coords::Cartesian_Point> com; //vector to save all center of masses
  std::ofstream masscenters;

  for (std::size_t i = 0u; i < N; i++)//i iterates over all molecules
  {
    com.push_back(coords.center_of_mass_mol(i));
  }

  masscenters.open("CenterofMasses.out");

  masscenters << N << '\n' << '\n';

  for (std::size_t i = 0u; i < N; i++)
  {
    masscenters << std::right << std::fixed << std::setprecision(7) << std::setw(5) << i + 1 << std::setw(13) << com[i].x() << std::setw(13) << com[i].y()
      << std::setw(13) << com[i].z() << '\n';
  }

  if (Config::set().center.dimer == true)
  {
    double abstkrit = Config::get().center.distance;
    double com_dist(0u);
    std::ofstream dimerstrukt;
    std::size_t size_i(0u), size_j(0u), size_dimer(0u);

    for (std::size_t i = 0u; i < N; i++) //Monomer i
    {
      for (std::size_t j = 0u; j < N ; j++) //Monomer j
      {
        if(j > i)
        {
        //Calculate distance between momomers
        com_dist = sqrt((com[i].x() - com[j].x())*(com[i].x() - com[j].x()) + (com[i].y() - com[j].y())*(com[i].y() - com[j].y()) + (com[i].z() - com[j].z())*(com[i].z() - com[j].z()));

        if (com_dist <= abstkrit)
        {
          size_i = coords.molecule(i).size();//Number of atoms in monomer i
          size_j = coords.molecule(j).size();;//Number of atoms in monomer j
          size_dimer = size_i + size_j;
          int korr_bond_i(0);
          int korr_bond_j(0);

          if (i > 0)//calculate coretion for bondingpartnerindices
          {
            for (int ki = i-1; ki >= 0; ki--)
            korr_bond_i += coords.molecule(ki).size();
          }

          if (j > 1)//calculate coretion for bondingpartnerindices
          {
            korr_bond_j = korr_bond_i;
            for (int kj = j-1; kj > i; kj--)
              korr_bond_j += coords.molecule(kj).size();
          }
      
          coords::Atoms const& atoms = coords.atoms();
          coords::Atoms truncatedAtoms;
          coords::Representation_3D positions;
          std::vector<std::size_t> new_index_of_atom(coords.size(), 0u);

          for (std::size_t k = 0u; k < size_i; k++) //loop for first monomer
          {
            truncatedAtoms.add(coords.atoms(atoms.atomOfMolecule(i,k)));
            positions.push_back(coords.xyz(atoms.atomOfMolecule(i, k)));
            new_index_of_atom.at(atoms.atomOfMolecule(i, k)) = truncatedAtoms.size() - 1u;
          }

          unsigned int helperIterator = 0u;
          for (std::vector<coords::Atom>::size_type l = truncatedAtoms.size() - atoms.molecule(i).size(); l < truncatedAtoms.size(); ++l, helperIterator++)//ashures only atoms of molecule considered at the moment are used
          {
            std::vector<std::size_t>  old_bonds = coords.atoms(atoms.atomOfMolecule(i, helperIterator)).bonds();

            for (auto& bonding_partner : coords.atoms(atoms.atomOfMolecule(i, helperIterator)).bonds())//range based loop for removal of old bondingpartners
            {
              truncatedAtoms.atom(l).detach_from(bonding_partner);
            }

            for (auto& bonding_partner : old_bonds) //range based loop to add new bonding partners | splitted removal and adding necessary to prevent deletion of false bonding index
            {
              truncatedAtoms.atom(l).bind_to(new_index_of_atom[bonding_partner]);
            }

          }//end bindingindices first monomer

          for (std::size_t k = 0u; k < size_j; k++) //loop for second monomer
          {
            truncatedAtoms.add(coords.atoms(atoms.atomOfMolecule(j, k)));
            positions.push_back(coords.xyz(atoms.atomOfMolecule(j, k)));
            new_index_of_atom.at(atoms.atomOfMolecule(j, k)) = truncatedAtoms.size() - 1u;
          }

          helperIterator = 0u;
          for (std::vector<coords::Atom>::size_type l = truncatedAtoms.size() - atoms.molecule(j).size(); l < truncatedAtoms.size(); ++l, helperIterator++)//ashures only atoms of molecule considered at the moment are used
          {
            std::vector<std::size_t>  old_bonds = coords.atoms(atoms.atomOfMolecule(j, helperIterator)).bonds();

            for (auto& bonding_partner : coords.atoms(atoms.atomOfMolecule(j, helperIterator)).bonds())//range based loop for removal of old bondingpartners
            {
              truncatedAtoms.atom(l).detach_from(bonding_partner);
            }

            for (auto& bonding_partner : old_bonds) //range based loop to add new bonding partners | splitted removal and adding necessary to prevent deletion of false bonding index
            {
              truncatedAtoms.atom(l).bind_to(new_index_of_atom[bonding_partner]);
            }

          }//end bindingindices second monomer

          coords::Coordinates newCoords;
          coords::PES_Point x(positions);

          newCoords.init_in(truncatedAtoms, x, true);

          std::stringstream oname;
          oname << "Dimerstrukt_" << i + 1 << "_" << j + 1 << ".xyz";
          dimerstrukt.open(oname.str());
          dimerstrukt << newCoords;
          dimerstrukt.close();

        }//dist check end
       }//if j>i end
      }//j
    }//i
  }//dimer check end
}