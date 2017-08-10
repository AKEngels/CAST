#include "Center.h"
#include<iostream>
#include<iomanip>





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

          std::stringstream oname;
          oname << "Dimerstrukt_" << i + 1 << "_" << j + 1 << ".xyz";
          dimerstrukt.open(oname.str());

          dimerstrukt << size_dimer << '\n';

          for (std::size_t k = 0u; k < size_i; k++) //loop for writing first monomer
          {
            std::size_t atom_index = coords.atoms().atomOfMolecule(i, k);
            dimerstrukt << std::right << std::fixed << std::setprecision(7) << std::setw(4) << k + 1 << std::setw(6) << coords.atoms(atom_index).symbol()
              << std::setw(13) << coords.xyz(atom_index).x() << std::setw(13) << coords.xyz(atom_index).y() << std::setw(13) << coords.xyz(atom_index).z()
              << std::setw(8) << coords.atoms(atom_index).energy_type();

            for (std::size_t l = 0u; l < coords.atoms(atom_index).bonds().size(); l++) //loop for bondpartners
            {
              dimerstrukt << std::right << std::setw(7) << coords.atoms(atom_index).bonds(l) + 1 - korr_bond_i;
            }
            dimerstrukt << '\n';
          }

          for (std::size_t k = 0u; k < size_j; k++) //loop for writing second monomer
          {
            std::size_t atom_index = coords.atoms().atomOfMolecule(j, k);
            dimerstrukt << std::right << std::fixed << std::setprecision(7) << std::setw(4) << k + size_i + 1 << std::setw(6) << coords.atoms(atom_index).symbol()
              << std::setw(13) << coords.xyz(atom_index).x() << std::setw(13) << coords.xyz(atom_index).y() << std::setw(13) << coords.xyz(atom_index).z()
              << std::setw(8) << coords.atoms(atom_index).energy_type();

            for (std::size_t l = 0u; l < coords.atoms(atom_index).bonds().size(); l++)//loop for bondpartners
            {
              dimerstrukt << std::right << std::setw(7) << coords.atoms(atom_index).bonds(l) + 1 - korr_bond_j;
            }
            dimerstrukt << '\n';
          }
          dimerstrukt.close();
        }
       }
      }

    }
  }
}