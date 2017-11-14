#include "cut_residues.h"

/**writes the new tinkerstructure
@param ref: (old) coordobject
@param bondvector: the vector that is created in function rebind (see there for detailed description)
@param new_atoms: vector of newly created atoms
@param new_positions: vector of positions of newly created atoms
@param stream: stream where the tinkerstructure should be written to*/
void write_tinker(coords::Coordinates ref, std::vector<std::vector<int>> bondvector, std::vector<coords::Atom> new_atoms, std::vector<coords::cartesian_type> new_positions, std::ostream & stream)
{
  std::size_t const N(bondvector.size()+new_atoms.size());
  stream << N << '\n';
  int counter = 1;
  for (auto i : bondvector)  // remaining atoms
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
  for (int i = 0; i < new_atoms.size(); i++)  // new atoms
  {
    stream << std::right << std::setw(6) << counter << "  ";
    stream << std::left << std::setw(3) << new_atoms[i].symbol().substr(0U, 2U);
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << new_positions[i].x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << new_positions[i].y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << new_positions[i].z();
    stream << std::right << std::setw(6) << new_atoms[i].energy_type();
    std::size_t const bSize(new_atoms[i].bonds().size());
    for (std::size_t j(0U); j < bSize; ++j)
    {
      stream << std::right << std::setw(6) << new_atoms[i].bonds()[j] + 1U;
    }
    stream << "\n";
    counter += 1;
  }
}

/**find positions for H-atoms that fill the amino group where a peptide bond is broken
@param coordobj: (old) coordobject
@param i: index of the N-atom which is to be saturated (old indexation)*/
std::vector<coords::cartesian_type> find_position_Hs(coords::Coordinates coordobj, int i)
{
  std::vector<size_t> bonded_atoms_indizes = coordobj.atoms(i).bonds();
  coords::cartesian_type n; // normal vector
  coords::cartesian_type p; // point on plane (H-atom)
  for (auto b : bonded_atoms_indizes)
  {
    if (coordobj.atoms(b).symbol() == "C" && coordobj.atoms(b).get_res_id() == coordobj.atoms(i).get_res_id())
    {
      n = coordobj.xyz(b) - coordobj.xyz(i);
    }
    else if (coordobj.atoms(b).symbol() == "H")
    {
      p = coordobj.xyz(b);
    }
  }
  double d = n.x()*p.x() + n.y()*p.y() + n.z()*p.z(); // equation for plane: ax+by+cz=d with (a,b,c) is normal vector

  // calculate intersection of normal vector with plane
  double lambda = (d - coordobj.xyz(i).x()*n.x() - coordobj.xyz(i).y()*n.y() - coordobj.xyz(i).z()*n.z()) / (n.x()*n.x() + n.y()*n.y() + n.z()*n.z());

  double sx = coordobj.xyz(i).x() + lambda * n.x();
  double sy = coordobj.xyz(i).y() + lambda * n.y();
  double sz = coordobj.xyz(i).z() + lambda * n.z();
  coords::cartesian_type s(sx, sy, sz);  // intersection of normal vector with plane

  coords::cartesian_type v = p - s; // vector on plane
  coords::cartesian_type v2 = scon::cross(v, n); 
  double desired_length = sqrt(v.x()*v.x() + v.y()*v.y() + v.z()*v.z());
  double current_length = sqrt(v2.x()*v2.x() + v2.y()*v2.y() + v2.z()*v2.z());
  v2 = v2 * (desired_length / current_length); // second vector on plane with the same length (perpendicular to first)

  // calculate position of the H atoms (= S - 0.5*v +- sqrt(3)/2 * v2)
  double vx = 0.5*v.x();
  double vy = 0.5*v.y();
  double vz = 0.5*v.z();
  coords::cartesian_type v05(vx, vy, vz);

  vx = sqrt(3)/2*v2.x();
  vy = sqrt(3) / 2 *v2.y();
  vz = sqrt(3) / 2 *v2.z();
  coords::cartesian_type v32(vx, vy, vz);

  coords::cartesian_type pos_H1 = s - v05 + v32;
  coords::cartesian_type pos_H2 = s - v05 - v32;

  std::vector<coords::cartesian_type> positions;
  positions.push_back(pos_H1);
  positions.push_back(pos_H2);

  return positions;
}

/**find positions for O-atom that fill the carboxyl group where a peptide bond is broken
@param coordobj: (old) coordobject
@param i: index of the C-atom which is to be saturated (old indexation)*/
coords::cartesian_type find_position_O(coords::Coordinates coordobj, int i)
{
  std::vector<size_t> bonded_atoms_indizes = coordobj.atoms(i).bonds();
  coords::cartesian_type n; // normal vector
  coords::cartesian_type p; // point on plane (O-atom)
  for (auto b : bonded_atoms_indizes)
  {
    if (coordobj.atoms(b).symbol() == "C")
    {
      n = coordobj.xyz(b) - coordobj.xyz(i);
    }
    else if (coordobj.atoms(b).symbol() == "O")
    {
      p = coordobj.xyz(b);
    }
  }
  double d = n.x()*p.x() + n.y()*p.y() + n.z()*p.z(); // equation for plane: ax+by+cz=d with (a,b,c) is normal vector

  // calculate intersection of normal vector with plane
  double lambda = (d - coordobj.xyz(i).x()*n.x() - coordobj.xyz(i).y()*n.y() - coordobj.xyz(i).z()*n.z()) / (n.x()*n.x() + n.y()*n.y() + n.z()*n.z());

  double sx = coordobj.xyz(i).x() + lambda * n.x();
  double sy = coordobj.xyz(i).y() + lambda * n.y();
  double sz = coordobj.xyz(i).z() + lambda * n.z();
  coords::cartesian_type s(sx, sy, sz);  // intersection of normal vector with plane

  // calculate position of O atom
  double vx = p.x() - sx;
  double vy = p.y() - sy;
  double vz = p.z() - sz;
  coords::cartesian_type v(vx, vy, vz);

  coords::cartesian_type pos_O = s - v;
  return pos_O;
}

/**creates the new bonds in the following form:
every element of the result is a vector that corresponds to one atom
the first element of the vector is the the index of the atom (old indexation)
the other elements are the atoms to which the current atom is bound (new indexation)
the index of the current atom in new indexation is identical to the index in the vector
@param coordobj: coordinates object (with all atoms)
@param indizes: vector of indizes that remain in the new structure*/
std::vector<std::vector<int>> rebind(coords::Coordinates coordobj, std::vector<int> indizes, std::vector<coords::Atom>&new_atoms, std::vector<coords::cartesian_type> &new_positions)
{
  std::vector<std::vector<int>> bondvector; 
  int counter = 0;
  for (auto i : indizes)
  {
    std::vector<int> bond;
    bond.push_back(i);
    for (auto b : coordobj.atoms(i).bonds())
    {
      int new_bond = find_index(b, indizes);
      if (new_bond == 99998)
      {
        if (coordobj.atoms(i).energy_type() == 180) // amide N
        {
          // add bond to 2 H atoms
          bond.push_back(indizes.size()+counter);  
          counter += 1;
          bond.push_back(indizes.size()+counter);
          counter += 1;
          
          // add 2 H atoms
          coords::Atom current_atom("H");
          current_atom.set_energy_type(233);
          current_atom.bind_to(find_index(i, indizes));
          new_atoms.push_back(current_atom);  // first atom
          new_atoms.push_back(current_atom);  // second atom

          // add positions for 2 H atoms
          std::vector<coords::cartesian_type> positions = find_position_Hs(coordobj, i);
          new_positions.push_back(positions[0]);
          new_positions.push_back(positions[1]);
        }
        else if (coordobj.atoms(i).energy_type() == 177) // amide C
        {
          // add bond to an O atom
          bond.push_back(indizes.size() + counter);
          counter += 1;

          // add a C-terminal O atom
          coords::Atom current_atom("O");
          current_atom.set_energy_type(214);
          current_atom.bind_to(find_index(i, indizes));
          new_atoms.push_back(current_atom);  

          // add position for the O atom
          coords::cartesian_type position = find_position_O(coordobj, i);
          new_positions.push_back(position);
        }
        else std::cout << "Strange things are happening.\n";
      }
      else bond.push_back(new_bond);
    }
    bondvector.push_back(bond);
  }
  return bondvector;
}


void cut_residues(coords::Coordinates coordobj, std::ostream & stream)
{
  // calculate geometrical center of reaction site
  coords::Cartesian_Point sum; 
  for (auto atom : Config::get().cut.react_atoms)
  {
    sum += coordobj.xyz(atom-1);
  }
  coords::Cartesian_Point react_site = sum / double(Config::get().cut.react_atoms.size());

  // determine which residues to keep
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

  // determine which atoms to keep
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
  
  // find new bonds
  std::vector<coords::Atom> new_atoms;
  std::vector<coords::cartesian_type> new_positions;
  std::vector<std::vector<int>> bondvector = rebind(coordobj,remaining_atoms,new_atoms,new_positions);

  if (Config::get().general.verbosity > 2)
  {
    std::cout << new_atoms.size() << " added for saturation\n";
    std::cout << remaining_atoms.size()+new_atoms.size() << " atoms in total\n";
  }

  // write tinkerstucture
  write_tinker(coordobj, bondvector, new_atoms, new_positions, stream);
};