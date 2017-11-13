#include "cut_residues.h"

void write_tinker(coords::Coordinates ref, std::vector<int> indizes, std::ostream & stream)
{
  std::size_t const N(indizes.size());
  stream << N << '\n';
  int counter = 1;
  for (auto i : indizes)
  {
    stream << std::right << std::setw(6) << counter << "  ";
    stream << std::left << std::setw(3) <<    ref.atoms(i).symbol().substr(0U, 2U);
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z();
    stream << std::right << std::setw(6) << ref.atoms(i).energy_type();
    std::size_t const bSize(ref.atoms(i).bonds().size());
    for (std::size_t j(0U); j < bSize; ++j)
    {
      stream << std::right << std::setw(6) << ref.atoms(i).bonds()[j] + 1U;
    }
    if (ref.atoms(i).sub_type() == coords::Atom::sub_types::ST_IN) stream << " IN";
    else if (ref.atoms(i).sub_type() == coords::Atom::sub_types::ST_OUT) stream << " OUT";
    stream << '\n';
    counter += 1;
  }
}

//void rebind(coords::Coordinates coordobj, std::vector<int> indizes)
//{
//  int counter = 0;
//  std::vector<int> bonds;
//  for (auto i : indizes)
//  {
//    std::vector<
//  }
//}


std::string cut_residues(coords::Coordinates coordobj, std::ostream & stream)
{
  coords::Cartesian_Point sum; 
  for (auto atom : Config::get().cut.react_atoms)
  {
    sum += coordobj.xyz(atom-1);
  }
  coords::Cartesian_Point react_site = sum / double(Config::get().cut.react_atoms.size());

  int counter = 0;
  std::vector<int> remaining_atoms;
  std::vector<int> remaining_resids;
  for (int i=0; i<coordobj.size();i++)
  {
    double distance = dist(react_site, coordobj.xyz(i));
    if (distance < Config::get().cut.distance)
    {
      counter += 1;
      int current_resid = coordobj.atoms(i).get_res_id();
      if (!is_in(current_resid, remaining_resids))
      {
        remaining_resids.push_back(current_resid);
      }
    }
  }
  if (Config::get().general.verbosity > 2)
  {
    std::cout << counter << " atoms in cutout radius\n";
    std::cout << remaining_resids.size() << " residues remaining\n";
  }

  for (int i = 0; i < coordobj.size(); i++)
  {
    if (is_in(coordobj.atoms(i).get_res_id(), remaining_resids))
    {
      remaining_atoms.push_back(i);
    }
  }
  if (Config::get().general.verbosity > 2)
  {
    std::cout << remaining_atoms.size() << " atoms remaining\n";
  }
  
  //rebind(coordobj,remaining_atoms);

  write_tinker(coordobj, remaining_atoms, stream);
};