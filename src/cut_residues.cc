#include "cut_residues.h"

void write_tinker(coords::Coordinates ref, std::vector<std::vector<int>> bondvector, std::ostream & stream)
{
  std::size_t const N(bondvector.size());
  stream << N << '\n';
  int counter = 1;
  for (auto i : bondvector)
  {
    stream << std::right << std::setw(6) << counter << "  ";
    stream << std::left << std::setw(3) <<    ref.atoms(i[0]).symbol().substr(0U, 2U);
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i[0]).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i[0]).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i[0]).z();
    stream << std::right << std::setw(6) << ref.atoms(i[0]).energy_type();
    std::size_t const bSize(i.size());
    for (std::size_t j(1U); j < bSize; ++j)
    {
      stream << std::right << std::setw(6) << i[j] + 1U;
    }
    if (ref.atoms(i[0]).sub_type() == coords::Atom::sub_types::ST_IN) stream << " IN";
    else if (ref.atoms(i[0]).sub_type() == coords::Atom::sub_types::ST_OUT) stream << " OUT";
    stream << '\n';
    counter += 1;
  }
}

std::vector<std::vector<int>> rebind(coords::Coordinates coordobj, std::vector<int> indizes)
{
  std::vector<std::vector<int>> bondvector; // fist element is index (old indexation), other elements are atoms to which it is bound (new indexation)
  for (auto i : indizes)
  {
    std::vector<int> bond;
    bond.push_back(i);
    for (auto b : coordobj.atoms(i).bonds())
    {
      bond.push_back(find_index(b,indizes));
    }
    bondvector.push_back(bond);
  }
  return bondvector;
}


void cut_residues(coords::Coordinates coordobj, std::ostream & stream)
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
  
  std::vector<std::vector<int>> bondvector = rebind(coordobj,remaining_atoms);

  write_tinker(coordobj, bondvector, stream);
};