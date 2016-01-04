#include <cmath>
#include <stdexcept>
//#include <fftw3.h>
#include "atomic.h"
#include "coords.h"
#include "configuration.h"
#include "coords_io.h"
#include "lbfgs.h"
#include "Array.h"
#define SUPERPI 3.141592653589793238

#include "optimization_dimer.h"
std::ostream& coords::operator<< (std::ostream &stream, internal_relations const & inter)
{
  stream << "[Rel: B,A,D: " << inter.rel_b << ", " << inter.rel_a << ", " << inter.rel_d;
  stream << "]:[Atom: " << inter.atom_of_inter_index << ", Inter: " << inter.inter_of_atom_index << "]:[Follow: ";
  for (std::size_t i = 0; i < inter.follow_ups.size(); ++i)
  {
    stream << (i != 0?", ":"") << inter.follow_ups[i];
  }
  stream << "]:[Dep: ";
  for (std::size_t i = 0; i < inter.dependant_torsions.size(); ++i)
  {
    stream << (i != 0?", ":"") << inter.dependant_torsions[i];
  }
  stream << "] Main: " << inter.is_main_torsion << "\n";
  return stream;
}

std::ostream& coords::operator<< (std::ostream &stream, Atom const & atom)
{
  stream << atom.m_symbol << " (" << atom.m_number << ", " << atom.m_mass << "):[";
  for (std::size_t i = 0; i < atom.m_bonds.size(); ++i)
  {
    stream << (i != 0?", ":"") << atom.m_bonds[i];
  }
  stream << "] Sys: " << atom.m_system << ", Type: " << atom.m_etype << ", Sub: " << atom.m_sub_id << '\n';
  stream << "[B,A,D: " << atom.rel_b << ", " << atom.rel_a << ", " << atom.rel_d;
  stream << "]:[Atom: " << atom.atom_of_inter_index << ", Inter: " << atom.inter_of_atom_index << "]:[Follow: ";
  for (std::size_t i = 0; i < atom.follow_ups.size(); ++i)
  {
    stream << (i != 0?", ":"") << atom.follow_ups[i];
  }
  stream << "]:[Dep: ";
  for (std::size_t i = 0; i < atom.dependant_torsions.size(); ++i)
  {
    stream << (i != 0?", ":"") << atom.dependant_torsions[i];
  }
  stream << "] Main: " << atom.is_main_torsion << '\n';
  return stream;
}

std::ostream& coords::operator<< (std::ostream &stream, Atoms const & atoms)
{
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    stream << atoms.atom(i);
  }
  return stream;
}

std::ostream& coords::operator<< (std::ostream &stream, Stereo const & stereo)
{
  for (auto const & center : stereo.centers())
  {
    stream << "Stereo, Atom: " << center.m_atom << ", Dir: " << center.m_dir << '\n';
  }
  return stream;
}

/* #############################################################################

  class PES_Point and helper static function implementations

  ########  ########  ######          ########   #######  #### ##    ## ######## 
  ##     ## ##       ##    ##         ##     ## ##     ##  ##  ###   ##    ##    
  ##     ## ##       ##               ##     ## ##     ##  ##  ####  ##    ##    
  ########  ######    ######          ########  ##     ##  ##  ## ## ##    ##    
  ##        ##             ##         ##        ##     ##  ##  ##  ####    ##    
  ##        ##       ##    ##         ##        ##     ##  ##  ##   ###    ##    
  ##        ########  ######  ####### ##         #######  #### ##    ##    ##  

############################################################################# */




static inline bool in_3d_delta (coords::Representation_3D const &r, 
  coords::Cartesian_Point const &d)
{
  std::size_t const N(r.size());
  for (std::size_t i=0; i<N; ++i)
  {
    if (std::abs(r[i].x()) > d.x() || 
      std::abs(r[i].y()) > d.y() || 
      std::abs(r[i].z()) > d.z())
    {
      return false;
    }
  }
  return true;
}

static inline bool out_of_angle_delta(coords::angle_type const & a, 
  coords::angle_type const & l)
{
  return (a > l || a < -l);
}

static inline bool in_internal_delta(coords::Representation_Internal const &r, 
  coords::internal_type const &d)
{
  std::size_t const N(r.size());
  for (std::size_t i = 0; i<N; ++i)
  {
    if (std::abs(r[i].radius()) > d.radius() || 
      out_of_angle_delta(r[i].inclination(), d.inclination()) || 
      out_of_angle_delta(r[i].azimuth(), d.azimuth()))
    {
      return false;
    }
  }
  return true;
}


static inline bool in_main_delta(coords::Representation_Main const &r, 
  coords::main_type const f)
{
  std::size_t const N(r.size());
  for (std::size_t i(0u); i<N; ++i)
  {
    if (r[i] > f || r[i] < -f)
    {
      //std::cout << "Structure not the same cause main delta ";
      //std::cout << i << " is " << r[i] << " which is bigger than " << f << "\n";
      return false;
    }
  }
  return true;
}


bool coords::PES_Point::equal_compare (PES_Point const &rhs) const
{
  using scon::operator-;
  //std::cout << "Checking delta: ";
  //std::cout << "M: " << structure.main << "\n";
  //std::cout << "RM: " << rhs.structure.main << "\n";
  auto main_delta = structure.main - rhs.structure.main;
  if (in_main_delta(main_delta, Config::get().coords.equals.main))
  {
    return true;
  }
  else if (in_internal_delta(structure.intern - rhs.structure.intern,
    Config::get().coords.equals.intern))
  {
    return true;
  }
  else if (in_3d_delta(structure.cartesian - rhs.structure.cartesian,
    Config::get().coords.equals.xyz))
  {
    return true;
  }
  //std::cout << "MAINDELTA: " << scon::vector_delimeter(' ') << main_delta << "\n";
  return false;
}




/* ############################################

  class Atom 
  non-inline implementation

        ###    ########  #######  ##     ## 
       ## ##      ##    ##     ## ###   ### 
      ##   ##     ##    ##     ## #### #### 
     ##     ##    ##    ##     ## ## ### ## 
     #########    ##    ##     ## ##     ## 
     ##     ##    ##    ##     ## ##     ## 
     ##     ##    ##     #######  ##     ## 

############################################ */




coords::Atom::Atom (std::string s)
  : m_symbol(s), m_number(atomic::atomic_number_by_symbol(s)), 
  m_mass(atomic::massMap[m_number]), m_system(0U), m_etype(0U), 
  m_sub_id(ST_DEFAULT), m_fix(false), m_intern_root(false)
{ }


coords::Atom::Atom (std::size_t n)
  : m_symbol(atomic::symbolMap[n]), m_number(n),
  m_mass(atomic::massMap[m_number]), m_system(0U), m_etype(0U), 
  m_sub_id(ST_DEFAULT), m_fix(false), m_intern_root(false)
{ }


coords::Atom::Atom (double m)
  : m_symbol(atomic::symbolMap[atomic::atomic_number_by_mass(m)]), 
    m_number(atomic::atomic_number_by_mass(m)), m_mass(m), m_system(0U), 
    m_etype(0U), m_sub_id(ST_DEFAULT), m_fix(false), m_intern_root(false)
{ }


void coords::Atom::swap (Atom &r)
{
  using std::swap;
  swap(rel_b, r.rel_b);
  swap(rel_a, r.rel_a);
  swap(rel_d, r.rel_d);
  swap(atom_of_inter_index, r.atom_of_inter_index);
  swap(inter_of_atom_index, r.inter_of_atom_index);
  swap(main_torsion_index, r.main_torsion_index);
  follow_ups.swap(r.follow_ups);
  dependant_torsions.swap(r.dependant_torsions);
  m_successors.swap(r.m_successors);
  swap(is_main_torsion, r.is_main_torsion);
  m_symbol.swap(r.m_symbol);
  swap(m_number, r.m_number);
  swap(m_mass, r.m_mass);
  m_bonds.swap(r.m_bonds);
  swap(m_system, r.m_system);
  swap(m_fix, r.m_fix);
  swap(m_intern_root, r.m_intern_root);
}




/* ##############################################

  class Atoms 
  non-inline implementation

     ###    ########  #######  ##     ##  ######  
    ## ##      ##    ##     ## ###   ### ##    ## 
   ##   ##     ##    ##     ## #### #### ##       
  ##     ##    ##    ##     ## ## ### ##  ######  
  #########    ##    ##     ## ##     ##       ## 
  ##     ##    ##    ##     ## ##     ## ##    ## 
  ##     ##    ##     #######  ##     ##  ######  

############################################## */




void coords::Atoms::refine ()
{
  refine_internals();
  //refine_followups();
  refine_subsystems();
  //refine_torsions();
}


bool coords::Atoms::sub_io_transition (std::size_t a, std::size_t b) const
{
  return ((m_atoms[a].sub_type() == Atom::ST_IN && 
	   m_atoms[b].sub_type() == Atom::ST_OUT) ||
      (m_atoms[b].sub_type() == Atom::ST_IN && 
	    m_atoms[a].sub_type() == Atom::ST_OUT));
}

// Helpers for external info for mains
namespace
{
  inline bool main_list_allowed(std::pair<std::size_t, std::size_t> const & p)
  {
    if (Config::get().coords.internal.main_whitelist.empty() &&
      !scon::sorted::exists(Config::get().coords.internal.main_blacklist, p) &&
      !scon::sorted::exists(Config::get().coords.internal.main_blacklist, std::make_pair(p.second, p.first)))
    {
      return true;
    }
    else if (scon::sorted::exists(Config::get().coords.internal.main_whitelist, p) ||
      scon::sorted::exists(Config::get().coords.internal.main_whitelist, std::make_pair(p.second, p.first)))
    {
      return true;
    }
    return false;
  }

  template<class ATOMS>
  inline std::vector<std::size_t> find_ring(std::vector<std::size_t> const & indices,
    ATOMS const &atoms)
  {
    using sv = std::vector<std::size_t>;
    std::size_t const IS = indices.size();
    if (IS < 1 || indices.size() > 6)
    {
      return sv();
    }

    for (std::size_t const b : atoms[indices.back()].bonds())
    {
      if (IS > 2)
      {
        bool was_in = false;
        for (std::size_t k = IS - 2; k > 0; --k)
        {
          if (b == indices[k])
          {
            was_in = true;
            break;
          }
        }
        if (was_in) continue;

        if (b == indices[0]) return indices;
      }
      sv ni(indices);
      ni.push_back(b);
      ni = find_ring(ni, atoms);
      if (!ni.empty()) return ni;
    }
    return sv();
  }

  void fix_rotation(std::vector<coords::Atom> &atoms)
  {
    auto const N = atoms.size();
    std::vector<bool> checked(atoms.size(), false);
    for (auto & a : scon::reverse_range(atoms))
    {
      auto const aa = a.i_to_a();
      if (checked[aa]) continue;
      if (atoms[aa].fixed())
      { // if we have cartesian fix we fix internals
        a.ifix(true);
        checked[aa] = true;
        auto ba = a.ibond();
        while (ba < N)
        {
          atoms[ba].ifix(true);
          checked[atoms[ba].i_to_a()] = true;
          ba = atoms[ba].ibond();
        }
      }
    }
  }

}


struct Part_of_Ring
{
  std::vector<bool> atom_is_part_of_ring;
  Part_of_Ring(std::vector<coords::Atom> const &atms) 
    : atom_is_part_of_ring(atms.size(), false)
  {
    std::size_t const N = atms.size();
    for (std::size_t i = 0; i < N; ++i)
    {
      if (!atom_is_part_of_ring[i])
      {
        auto const ring_indices =
          find_ring(std::vector<std::size_t>(1, i), atms);
        for (auto ri : ring_indices)
        {
          atom_is_part_of_ring[ri] = true;
        }
      }
    }
  }
  bool operator() (std::size_t const i) 
  { return atom_is_part_of_ring[i]; }
};

//#define PRINT_MAIN_AXIS

void coords::Atoms::refine_mains()
{
  fix_rotation(m_atoms);
  std::size_t const N = m_atoms.size();
  Part_of_Ring ringpart(m_atoms);
  // vector saving the atoms which have main torsions attached
  std::vector<std::size_t> atom_has_main_torsion_attached(N+1U, N+1U);
  for (auto const & mti : main_torsion_indices)
  {
    if (atom(mti).ibond() <= N)
      atom_has_main_torsion_attached[atom(mti).ibond()] = mti;
  }
  std::cout << std::boolalpha;

  bool const decoupled_molecules = m_molecules.size() > 1 && 
    Config::get().coords.decouple_internals && 
    Config::get().coords.internal.main_whitelist.empty();
  
#ifdef PRINT_MAIN_AXIS
  std::vector<std::size_t> main_axis_atoms;
#endif

  for (std::size_t i = 0; i < N; ++i)
  {
    std::size_t const ib = atom(i).ibond(), // bond internal
      ia = atom(i).iangle(), id = atom(i).idihedral(); // angle internal
    bool saturated_bond = false;
    bool list_allowed = false;
    bool molecule_rotation = ib >= N || ia >= N || 
      (ib < N && ia < N && id >= N && atom(atom(ia).i_to_a()).bonds().size() == 1 
        && atom(atom(ia).i_to_a()).bonds()[0u] == atom(ib).i_to_a());
    bool rotates_no_ringbond = true;

    if (ib < N && ia < N)
    {
      std::size_t const b = atom(ib).i_to_a(), // bond atom
        a = atom(ia).i_to_a(); // angle atom
      // saturated ? (num bonds = saturated num bonds)
      bool const a_saturated = atomic::saturated(atom(a).number(), atom(a).bonds().size()),
        b_saturated = atomic::saturated(atom(b).number(), atom(b).bonds().size());
      // a and b bound to each other
      bool const axis_is_bond = scon::sorted::exists(atom(b).bonds(), a) &&
        scon::sorted::exists(atom(a).bonds(), b);
      // rotation axis = bond between unsaturated atoms or no bond?
      saturated_bond = a_saturated || b_saturated || !axis_is_bond;
      // is the current main allowed with respect to black- or whitelisting?
      list_allowed = main_list_allowed({ a,b });
      // if one is not part of a ring
      rotates_no_ringbond = !ringpart(b) || !ringpart(a);
    }
    // No further internals attached to current
    bool const rotates_single_atom = atom(i).bound_internals().size() < 1u;
    // Is current internal rotating hydrogen?
    bool const not_only_hydrogen = !rotates_single_atom || 
      atom(atom(i).i_to_a()).number() > 1u;
    // is hydrogen only rotation to be removed?
    bool const non_removed_rotation = not_only_hydrogen || 
      !Config::get().coords.remove_hydrogen_rot;
    //! is this atom fixed?
    bool const not_fixed = !m_atoms[i].ifix();
    // does a main torsion exist for the current axis?
    bool main_does_not_yet_exist = true;
    // Add dependant dihedral if main exists
    if (ib <= N && atom_has_main_torsion_attached[ib] < N)
    {
      main_does_not_yet_exist = false;
      atom(atom_has_main_torsion_attached[ib]).append_dependant_idihedral(i);
    }
    // Add main if possible
    if ((molecule_rotation && decoupled_molecules) || 
      (!molecule_rotation && saturated_bond && non_removed_rotation && list_allowed && 
        main_does_not_yet_exist && rotates_no_ringbond && not_fixed))
    {
#ifdef PRINT_MAIN_AXIS
      std::cout << "Atom " << atom(i).i_to_a() << " is rotated as main with axis ";
      std::cout << (ib < N ? atom(ib).i_to_a() : ib) << "--";
      std::cout << (ia < N ? atom(ia).i_to_a() : ia) << "\n";
      if (ib < N) scon::sorted::insert_unique(main_axis_atoms, atom(ib).i_to_a());
      if (ia < N) scon::sorted::insert_unique(main_axis_atoms, atom(ia).i_to_a());
#endif
      atom_has_main_torsion_attached[ib] = i;
      scon::sorted::insert_unique(main_torsion_indices, i);
    }

  }
  if (Config::get().general.verbosity > 4)
  {
    std::cout << "Identified " << main_torsion_indices.size() << " main torsions.\n";
  }
  
#ifdef PRINT_MAIN_AXIS
  for (std::size_t i = 0; i < N; ++i)
  {
    if (!scon::sorted::exists(main_axis_atoms, i))
    {
      std::cout << i << " ";
    }
  }
  std::cout << "\n";

  for (auto a : main_axis_atoms)
  {
    std::cout << a << " ";
  }
  std::cout << "\n";
#endif
  

}


void coords::Atoms::refine_internals ()
{
  // First of all we assume that we do not have any mains
  main_bond_indices.clear();
  main_angle_indices.clear();
  main_torsion_indices.clear();
  // Build internals
  std::size_t const N(size());
  std::vector<bool> done(N, false);
  std::size_t current_internal(0u);
  // cycle rest
  for (std::size_t i(0U); i<N; ++i)
  {
    if (done[i]) continue;
    auto connect_it = Config::get().coords.internal.connect.find(i);
    std::size_t j = 0;
    //std::cout << "Atom " << i << " begins new molecule with internal " << current_internal << lineend;
    //// If current internal is not done yet we have a new molecule
    //// Internals of new molcules are always "main" coordinates
    //if (current_internal != 0 || Config::get().coords.decouple_internals)
    //{
    //  scon::sorted::insert_unique(main_bond_indices, current_internal);
    //  scon::sorted::insert_unique(main_angle_indices, current_internal);
    //  scon::sorted::insert_unique(main_torsion_indices, current_internal);
    //}
    if (connect_it != Config::get().coords.internal.connect.end() &&
      (*connect_it).second < N && done[(*connect_it).second])
    { // ...  if we find a required connection we add it
      j = atom((*connect_it).second).a_to_i();
    }
    else if (Config::get().coords.decouple_internals || current_internal == 0)
    { // else if we decouple internals we root every new molecule ...
      atom(current_internal).root();
      j = N;
    }
    // ...  otherwise we attach to the last known internal
    else
    { // append new molecule to last internal
      j = current_internal - 1u;
    }
    get_relatives(current_internal, j);
    atom(current_internal).set_a_of_i(i);
    atom(i).set_i_of_a(current_internal);
    size_1d this_molecule(1, i);
    done[i] = true;
    append_atoms(0U, i, this_molecule, ++current_internal, done);
    m_molecules.push_back(this_molecule);
  }
  refine_mains();
}


void coords::Atoms::append_atoms(std::size_t const lvl, std::size_t const A, 
  size_1d &molecule, std::size_t &index_size, std::vector<bool> &done)
{
  std::size_t const nBound = m_atoms[A].bonds().size();
  //std::cout << "Appending " << nBound << " atoms bound to " << A << "\n";
  for (std::size_t i{ 0U }; i < nBound; ++i)
  {
    std::size_t index = m_atoms[A].bonds()[i]; 
    //std::cout << "Bound " << i << " is " << index << "\n";
    if (!done[index])
    {
      //std::cout << "Not done yet: " << index << " binding to " << m_atoms[A].a_to_i() <<  "\n";
      get_relatives(index_size, m_atoms[A].a_to_i());
      done[index] = true;
      molecule.push_back(index);
      m_atoms[index_size].set_a_of_i(index);
      m_atoms[index].set_i_of_a(index_size);
      append_atoms(lvl+1, index, molecule, ++index_size, done);
    }
  }
}


void coords::Atoms::get_relatives (std::size_t const i, const std::size_t b)
{
  std::size_t const S = size();
  // if i is root, relatives are S / S+1 and S+2
  if (atom(i).is_root())
  {
    atom(i).set_ibond(S);
    atom(i).set_iangle(S + 1);
    atom(i).set_idihedral(S + 2);
  }
  // otherwise i binds to b
  else
  {
    // bind i to b
    atom(i).set_ibond(b);
    // attach b to i
    atom(b).attach(i);
    // b has valid bond relative?
    if (atom(b).ibond() < S)
    {
      // set angle relative to bond relative of b
      atom(i).set_iangle(atom(b).ibond());
      // b has valid angle relative?
      if (atom(b).iangle() < S)
      {
        // set dihedral relative to angle relative of b
        atom(i).set_idihedral(atom(b).iangle());
      }
      else
      {
        std::size_t const BBS = atom(atom(b).ibond()).bound_internals().size();
        if (BBS < 2u)
        {
          // if no other internals than b are bound to 
          // bond relative of b we utilize the angle relative anyways
          // iangle should be S or S+1
          if (atom(b).iangle() > (S + 1))
          {
            throw std::logic_error("Invalid angle partner.");
          }
          atom(i).set_idihedral(atom(b).iangle());
        }
        else
        {
          // otherwise we use the firste relative that is not b
          atom(i).set_idihedral(
            atom(atom(b).ibond()).bound_internals(0u) == b ?
              atom(atom(b).ibond()).bound_internals(1u) 
              : atom(atom(b).ibond()).bound_internals(0u)
            );
        }
      }
    }
    else
    {
      switch (atom(b).bound_internals().size())
      {
        case 0u:
        {
          throw std::logic_error("No atom bound to selection. At least current atom should be.");
        }
        // 1 internal bound is just i
        case 1u:
        {
          atom(i).set_iangle(atom(b).ibond());
          atom(i).set_idihedral(atom(b).iangle());
          break;
        }
        // besides i we have one other internal for angle relative
        case 2u:
        {
          atom(i).set_iangle(
            atom(b).bound_internals(0u) == i ? 
              atom(b).bound_internals(1u) : atom(b).bound_internals(0u)
            );
          atom(i).set_idihedral(atom(b).ibond());
          break;
        }
        // size > 2 means we have 2 other internals bound to b to set relation of i
        default:
        {
          if (atom(b).bound_internals(0u) == i)
          {
            atom(i).set_iangle(atom(b).bound_internals(1u));
            atom(i).set_idihedral(atom(b).bound_internals(2u));
          }
          else if (atom(b).bound_internals(1u) == i)
          {
            atom(i).set_iangle(atom(b).bound_internals(0u));
            atom(i).set_idihedral(atom(b).bound_internals(2u));
          }
          else
          {
            atom(i).set_iangle(atom(b).bound_internals(0u));
            atom(i).set_idihedral(atom(b).bound_internals(1u));
          }
          break;
        }
      }
    }
  }
}

bool coords::Atoms::common_torsion_axis(std::size_t a, std::size_t b, bool & direction) const
{
  if ((m_atoms[a].ibond() == m_atoms[b].ibond() && 
    m_atoms[a].iangle() == m_atoms[b].iangle()))
  {
    direction = true;
    return true;
  }
  else if (m_atoms[a].iangle() == m_atoms[b].ibond() && 
    m_atoms[a].ibond() == m_atoms[b].iangle())
  {
    direction = false;
    return true;
  }
  return false;
}


bool coords::Atoms::common_bond(std::size_t a, std::size_t b) const
{
  return (a < size() && b < size() && m_atoms[a].ibond() == m_atoms[b].ibond());
}


void coords::Atoms::refine_subsystems ()
{
  m_sub_systems.clear();
  std::size_t const N = m_atoms.size();
  std::size_t max_sys = 0;
  for (std::size_t i = 0; i < N; ++i)
  {
    max_sys = atom(i).system() > max_sys ? atom(i).system() : max_sys;
  }
  ++max_sys;
  std::vector<bool> found(max_sys, false);
  std::vector<std::size_t> map(max_sys);
  std::size_t sys;
  for (std::size_t i=0; i<N; ++i)
  {
    if (!found[atom(i).system()])
    {
      sys = m_sub_systems.size();
      map[atom(i).system()] = sys;
      m_sub_systems.push_back(size_1d());
      found[atom(i).system()] = true;
    }
    else
    {
      sys = map[atom(i).system()];
    }
    m_sub_systems[sys].push_back(i);
    m_atoms[i].assign_to_system(sys);
    if (m_atoms[i].sub_type() == Atom::ST_IN) 
    {
      m_in_exists = true;
      m_sub_in_index = sys;
    }
    else if (m_atoms[i].sub_type() == Atom::ST_OUT) 
    {
      m_out_exists = true;
      m_sub_out_index = sys;
    }
  }
  if (m_in_exists && m_out_exists) m_sub_io = true;
}

coords::Cartesian_Point coords::Atoms::rel_xyz(std::size_t const index, 
  coords::Representation_3D const & xyz) const
{
  std::size_t const N = xyz.size();
  if (index < N) return xyz[atom(index).i_to_a()];
  else if (index == N) return coords::Cartesian_Point(1, 0, 0);
  else if (index == N + 1) return coords::Cartesian_Point(0, 1, 0);
  else if (index == N + 2) return coords::Cartesian_Point(0, 0, 1);
  else throw std::logic_error("Wrong relative position requested. i > N + 2");
}

void coords::Atoms::c_to_i_light(PES_Point &p) const
{
  using scon::dot;
  using scon::spherical;
  using scon::geometric_length;
  std::size_t const N(m_atoms.size());
  std::size_t const M = main_torsion_indices.size();
  // Check size
  if (N != p.structure.cartesian.size())
    throw std::logic_error("ERR_COORD_internal_INDEXATION");
  // internal coords & gradients
  Representation_Internal & intern(p.structure.intern);
  Gradients_Internal & gintern(p.gradient.intern);
  // cartesian coords & gradients
  Representation_3D const & xyz(p.structure.cartesian);
  Gradients_3D const & gxyz(p.gradient.cartesian);
  // set correct sizes for internal data vectors
  intern.resize(N);
  gintern.resize(N);
  // Mains
  p.structure.main.resize(M);
  p.gradient.main.resize(M);
  for (std::size_t i(0u); i < N; ++i)
  {
    Representation_3D::size_type const ind_i = atom(i).i_to_a();
    coords::Cartesian_Point const & rel_bond = rel_xyz(atom(i).ibond(), xyz);
    coords::Cartesian_Point const & rel_angle = rel_xyz(atom(i).iangle(), xyz);
    coords::Cartesian_Point const & rel_dihedral = rel_xyz(atom(i).idihedral(), xyz);
    intern[i] = spherical(xyz[ind_i], rel_bond, rel_angle - rel_bond, rel_dihedral - rel_bond);
  }
  for (std::size_t j = 0; j < M; ++j)
  {
    auto const mti = main_torsion_indices[j];
    //std::cout << "Main Torsion " << j << " which is internal " << mti << " is " << intern[mti].azimuth() << lineend;
    p.structure.main[j] = intern[mti].azimuth();
  }
}

void coords::Atoms::c_to_i(PES_Point &p) const
{
  using scon::dot; 
  using scon::spherical;
  using scon::geometric_length;
  std::size_t const N(m_atoms.size());
  std::size_t const M = main_torsion_indices.size();
  // Check size
  if (N != p.structure.cartesian.size())
    throw std::logic_error("ERR_COORD_internal_INDEXATION");
  // internal coords & gradients
  Representation_Internal & intern(p.structure.intern);
  Gradients_Internal & gintern(p.gradient.intern);
  // cartesian coords & gradients
  Representation_3D const & xyz(p.structure.cartesian);
  Gradients_3D const & gxyz(p.gradient.cartesian);
  // set correct sizes for internal data vectors
  intern.resize(N);
  gintern.resize(N);
  gintern.assign(N, internal_gradient_type());

  // Mains
  p.structure.main.resize(M);
  p.gradient.main.resize(M);
  // zero main gradient vectors
  p.structure.main.assign(M, main_type());
  p.gradient.main.assign(M, main_gradient_type());

  //auto const sector_arc_length = 

  // Calculation
  for (std::size_t i(0u); i < N; ++i)
  {

    Representation_3D::size_type const ind_i = atom(i).i_to_a();
    coords::Cartesian_Point const & rel_bond = rel_xyz(atom(i).ibond(), xyz);
    coords::Cartesian_Point const & rel_angle = rel_xyz(atom(i).iangle(), xyz);
    coords::Cartesian_Point const & rel_dihedral = rel_xyz(atom(i).idihedral(), xyz);

    intern[i] = spherical(xyz[ind_i], rel_bond, rel_angle - rel_bond, rel_dihedral - rel_bond);

    //intern[i] = xyz[ind_i].spherical(rel_bond, rel_angle - rel_bond, rel_dihedral - rel_bond);

    auto j(i);
    while (j < N)
    { 

      auto const ind_j = atom(j).i_to_a();
      
      auto const rel_bond_j = rel_xyz(atom(j).ibond(), xyz);
      auto const rel_angle_j = rel_xyz(atom(j).iangle(), xyz);

      auto const Zj = normalized((rel_angle_j - rel_bond_j));
      auto const Dj = xyz[ind_j] - rel_bond_j;
      auto const Di = xyz[ind_i] - rel_bond_j;
      auto const PD = normalized(Dj);
      auto const RA = normalized(cross(Zj, PD));
      auto const PA = normalized(cross(RA, Di));
      auto const PT = normalized(cross(Zj, Di));
      auto const dA = i == j ? geometric_length(Di) : geometric_length((Di - RA*dot(Di, RA)));
      auto const dT = geometric_length(Di - Zj*dot(Di, Zj));

      coords::internal_gradient_type ig(
        dot(gxyz[ind_i], PD),
        dot(gxyz[ind_i], PA)*dA,
        dot(gxyz[ind_i], PT)*dT
        );

      gintern[j].x() += dot(gxyz[ind_i], PD);    // dot(gxyz[ind_i], PD)
      gintern[j].y() += dot(gxyz[ind_i], PA)*dA; // dot(gxyz[ind_i], PA)*dA
      gintern[j].z() += dot(gxyz[ind_i], PT)*dT; // dot(gxyz[ind_i], PT)*dT

      //if (i == j)
      //{
      //  coords::Cartesian_Point p;
      //  p += PD*ig.x();
      //  p += PA*(ig.y()/dA);
      //  p += PT*(ig.z()/dT);
      //  auto delt = gxyz[ind_i] - p;
      //  if (abs(geometric_length(delt)) > 1.e-10) std::cout << delt << "(" << gxyz[ind_i] << " <-> " << p << ")" << "\n";
      //  //std::cout << "Dots: " << dot(PD, PA) << ", " << dot(PD, PT) << ", " << dot(PA, PT) << "\n";
      //  //std::cout << "i was " << gxyz[ind_i] << ", is " << p << "\n";
      //}

      j = atom(j).ibond();

    }

  }

  for (std::size_t j = 0; j < M; ++j)
  {
    auto const mti = main_torsion_indices[j];
    //std::cout << "Main Torsion " << j << " which is internal " << mti << " is " << intern[mti].azimuth() << lineend;
    p.structure.main[j] = intern[mti].azimuth();
    if (atom(mti).ibond() < N)
    {
      for (auto const rotating : atom(atom(mti).ibond()).bound_internals())
      {
        //std::cout << "Rotates: " << rotating << " which is atom " << atom(rotating).i_to_a() << lineend;
        p.gradient.main[j] += gintern[rotating].z();
      }
      
    }
    
  }

}

void coords::Atoms::i_to_c(PES_Point &p) const
{
  using scon::append_spherical_NERF;
  std::size_t const N(m_atoms.size());
  // Check size
  if (N != p.structure.intern.size())
    throw std::logic_error("ERR_COORD_internal_INDEXATION");
  // internal coords & gradients
  Representation_Internal const & intern(p.structure.intern);
  Gradients_Internal gintern(p.gradient.intern);
  // cartesian coords & gradients
  Representation_3D & xyz(p.structure.cartesian);
  Gradients_3D & gxyz(p.gradient.cartesian);
  xyz.resize(N);
  // zero gradient vectors
  gxyz.assign(N, cartesian_gradient_type());

  for (std::size_t i(0U); i < N; ++i)
  {

    coords::Cartesian_Point const & rel_bond = rel_xyz(atom(i).ibond(), xyz);
    coords::Cartesian_Point const & rel_angle = rel_xyz(atom(i).iangle(), xyz);
    coords::Cartesian_Point const & rel_dihedral = rel_xyz(atom(i).idihedral(), xyz);
    coords::Cartesian_Point const & zenith = rel_angle - rel_bond;
    coords::Cartesian_Point const & azimuth_reference = rel_dihedral - rel_angle;
    xyz[atom(i).i_to_a()] = append_spherical_NERF(rel_bond, zenith, azimuth_reference, intern[i]);
    //xyz[atom(i).i_to_a()] = rel_bond.append_spherical_NERF(zenith, azimuth_reference, intern[i]);

  }

  // Unwind internal gradients

  for (std::size_t k(1U); k <= N; ++k)
  {

    std::size_t const i(N - k);
    std::size_t const ind_i(atom(i).i_to_a());

    coords::Cartesian_Point const & rel_bond = rel_xyz(atom(i).ibond(), xyz);
    coords::Cartesian_Point const & rel_angle = rel_xyz(atom(i).iangle(), xyz);
    coords::Cartesian_Point const bond_direction(normalized((xyz[ind_i] - rel_bond)));
    coords::Cartesian_Point const & zenith = rel_angle - rel_bond;
    coords::Cartesian_Point const dihedral_direction(normalized(cross(zenith, bond_direction)));
    coords::Cartesian_Point const angle_direction(normalized(cross(dihedral_direction, bond_direction)));

    gxyz[ind_i] = bond_direction*gintern[i].x() + angle_direction*gintern[i].y() + dihedral_direction*gintern[i].z();

    // iterate down the tree and remove gradient of i from lower internal gradients
    std::size_t j(i);
    while (j < N)
    {

      Representation_3D::size_type const ind_j = atom(j).i_to_a();

      coords::Cartesian_Point const & rel_bond_j = rel_xyz(atom(j).ibond(), xyz);
      coords::Cartesian_Point const & rel_angle_j = rel_xyz(atom(j).iangle(), xyz);

      auto const Zj = normalized((rel_angle_j - rel_bond_j));
      auto const Dj = xyz[ind_j] - rel_bond_j;
      auto const Di = xyz[ind_i] - rel_bond_j;
      auto const PD = normalized(Dj);
      //Cartesian_Point const RA = normalized(cross(Zj, Dj));
      auto const RA = normalized(cross(Zj, Dj));
      //Cartesian_Point const PA = normalized(cross(RA, Di));
      auto const PA = normalized(cross(RA, Di));
      //float_type const dA = len(Di - RA*dot(Di, RA));
      auto const dA = geometric_length(Di - RA*dot(Di, RA));
      //Cartesian_Point const PT = normalized(Zj.crossd(Di));
      auto const PT = normalized(cross(Zj, Di));
      //float_type const dT = len(Di - Zj*dot(Di, Zj));
      auto const dT = geometric_length(Di - Zj*dot(Di, Zj));

      gintern[j].x() -= dot(gxyz[ind_i], PD);
      gintern[j].y() -= dot(gxyz[ind_i], PA)*dA;
      gintern[j].z() -= dot(gxyz[ind_i], PT)*dT;

      j = atom(j).ibond();

    }

  }

}

//void coords::Atoms::internal_to_cartesian (PES_Point &p) const
//{
//  Representation_Internal const & intern(p.structure.intern);
//  Representation_3D & xyz(p.structure.cartesian);
//  const std::size_t N = m_atoms.size();
//  if (N != intern.size()) throw std::logic_error("ERR_COORD_internal_WRONG_COUNT");
//  xyz.resize(N);
//  if (N > 0U)
//  { 
//    xyz[m_atoms[0U].i_to_a()].init(0.0, 0.0, 0.0);
//    if (N > 1U)
//    {
//      xyz[m_atoms[1U].i_to_a()].init(intern[1U].radius(), 0.0, 0.0);
//      if (N > 2U)
//      {
//        std::size_t const atom_bound_to2(m_atoms[m_atoms[2U].ibond()].i_to_a());
//        float_type const cosine = (m_atoms[2U].ibond() == 1U ? -intern[2U].inclination().cos() : intern[2U].inclination().cos());
//        xyz[m_atoms[2U].i_to_a()].init(intern[2U].radius()*cosine + xyz[atom_bound_to2].x(), intern[2U].radius()*intern[2U].inclination().sin(), 0.0);
//        if (N > 3U)
//        {
//          for (std::size_t i(3U); i<N; ++i)
//          {
//            std::size_t const atom_bound_to_i(m_atoms[m_atoms[i].ibond()].i_to_a());
//            std::size_t const atom_angle_to_i(m_atoms[m_atoms[i].iangle()].i_to_a());
//            std::size_t const atom_dihed_to_i(m_atoms[m_atoms[i].idihedral()].i_to_a());
//              xyz[m_atoms[i].i_to_a()] = xyz[atom_bound_to_i].appendNERF(
//                 xyz[atom_angle_to_i]-xyz[atom_dihed_to_i], 
//                 xyz[atom_angle_to_i]-xyz[atom_bound_to_i], 
//                 intern[i].radius(), intern[i].inclination(), intern[i].azimuth()
//              );
//          }
//        }
//      }
//    }
//  }
//}



//
//void coords::Atoms::cartesian_to_internal (PES_Point &p) const
//{
//  // if the internals size does not match the number of atoms we need to recreate the indices
//  std::size_t const N(m_atoms.size());
//  if (N != p.structure.cartesian.size()) 
//    throw std::logic_error("ERR_COORD_internal_INDEXATION");
//  // internal coords & gradients
//  Representation_Internal & intern(p.structure.intern);
//  Gradients_Internal & gintern(p.gradient.intern);
//  // cartesian coords & gradients
//  Representation_3D const & xyz(p.structure.cartesian);
//  Gradients_3D const & gxyz(p.gradient.cartesian);
//  // set correct sizes for internal data vectors
//  intern.resize(N);
//  gintern.resize(N);
//  gintern.assign(N, internal_gradient_type());
//  // Mains
//  p.structure.main.resize(main_torsion_indices.size());
//  p.gradient.main.resize(main_torsion_indices.size());
//  // zero main gradient vectors
//  p.structure.main.assign(main_torsion_indices.size(), main_type());
//  p.gradient.main.assign(main_torsion_indices.size(), main_gradient_type());
//  if (N > 1U) 
//  {
//    coords::Cartesian_Point bond_axis(xyz[atom_by_intern(1U)] - xyz[atom_by_intern(m_atoms[1u].ibond())]);
//    intern[1U] = internal_type(len(bond_axis), angle_type(0.0), angle_type(0.0));
//    gintern[1U].x() = dot(gxyz[atom_by_intern(1U)], bond_axis.norm());
//    for (auto follower : m_atoms[1U].followers())
//    {
//      gintern[1U].x() += dot(gxyz[atom_by_intern(follower)], bond_axis);
//    }
//  }
//  if (N > 2U) 
//  {
//    std::size_t const a2(atom_by_intern(2u)), 
//      rb2(atom_by_intern(m_atoms[2u].ibond())),
//      ra2(atom_by_intern(m_atoms[2u].iangle()));
//    coords::Cartesian_Point bond_axis(xyz[a2] -  xyz[rb2]),
//      progenitor_bond_axis_inv(xyz[ra2] - xyz[rb2]),
//      bonds_cross(progenitor_bond_axis_inv.crossd(bond_axis));
//    intern[2U].init(len(bond_axis), progenitor_bond_axis_inv.angle(bond_axis), angle_type(0.0));
//    gintern[2U].x() = dot(gxyz[a2], bond_axis.norm());
//    gintern[2U].y() = dot(gxyz[a2], bonds_cross.crossd(bond_axis).norm());
//    for (auto follower : m_atoms[2U].followers())
//    {
//      std::size_t const following_atom(atom_by_intern(follower));
//      gintern[2U].x() += dot(gxyz[following_atom], bond_axis);
//      gintern[2U].y() += dot(gxyz[following_atom], bonds_cross.crossd(xyz[following_atom]-xyz[rb2]).norm());
//    }
//  }
//  if (N > 3U)
//  {
//    for (std::size_t i(3U); i<N; ++i)
//    {
//
//      std::size_t const a(atom_by_intern(i)), 
//        rbi(atom_by_intern(m_atoms[i].ibond())), // atom that is "bond-related"
//        rai(atom_by_intern(m_atoms[i].iangle())), // atom that is "angle-related"
//        rdi(atom_by_intern(m_atoms[i].idihedral())); // atom that is "dihedral-related"
//
//      coords::Cartesian_Point bond_axis(xyz[a] - xyz[rbi]),
//        prequel_bond_axis_inv(xyz[rai] - xyz[rbi]),
//        bonds_cross(prequel_bond_axis_inv.crossd(bond_axis)),
//        torsional_relative(xyz[rdi] - xyz[rai]),
//        torsion_cross(prequel_bond_axis_inv.crossd(torsional_relative));
//
//      intern[i].init(len(bond_axis), 
//        prequel_bond_axis_inv.angle(bond_axis), 
//        torsion_cross.angle(bonds_cross));
//
//      coords::Cartesian_Point const cross_cross(bonds_cross.crossd(torsion_cross));
//      double const norm = len(prequel_bond_axis_inv)*len(cross_cross);
//
//      // reverse direction if necessary
//      if (std::fabs(norm) > 0.0 && dot(prequel_bond_axis_inv, cross_cross)/norm < 0.0) 
//        intern[i].azimuth() = -intern[i].azimuth();
//
//      // apply to main dihedral coordinate as well
//      if (m_atoms[i].is_main_idihedral()) p.structure.main[m_atoms[i].main_idihedral_index()] = intern[i].azimuth();
//
//      gintern[i].x() = dot(gxyz[a], bond_axis.norm());
//      gintern[i].y() = dot(gxyz[a], bonds_cross.crossd(bond_axis).norm());
//      gintern[i].z() = dot(gxyz[a], bonds_cross.norm().invertd());
//
//      // sum follower gradients
//      for (auto const follower : m_atoms[i].followers())
//      {
//        std::size_t const following_atom(atom_by_intern(follower));
//        coords::Cartesian_Point const delta(xyz[following_atom] - xyz[rbi]);
//        gintern[i].x() += dot(gxyz[following_atom], bond_axis);
//        gintern[i].y() += dot(gxyz[following_atom], normalized(bonds_cross.crossd(delta)));
//        gintern[i].z() += dot(gxyz[following_atom], -normalized(cross(prequel_bond_axis_inv, delta)));
//      }
//
//      // sum main torsional gradients
//      std::size_t const num_maintors(main_torsion_indices.size());
//      for (std::size_t j(0U); j<num_maintors; ++j)
//      {
//        bool direction(true);
//        if (common_torsion_axis(main_torsion_indices[j], i, direction)) 
//        {
//          if (direction) p.gradient.main[j] += gintern[i].z();
//          else p.gradient.main[j] -= gintern[i].z();
//        }
//      }
//    }
//  }
//}


bool coords::Atoms::res_is_equal(std::size_t const a, std::size_t const b, 
  std::size_t const from_a, std::size_t const from_b, std::size_t depth) const
{
  // Different atoms (different atomic number) .. do not tend to be equal ;)
  if (m_atoms[a].number() != m_atoms[b].number()) return false;
  // Do we go in a circle? Then we are probably identical.
  if (a == b || (a == from_b && b == from_a)) return true;
  size_1d::size_type const TA(m_atoms[a].bonds().size());
  // If we did not find any differences in 100 layers we stop
  // and claim both residues to be equal
  if (depth > 100) return true;
  // not the same number of bound atoms? -> not an equal residue -> return false
  if (TA != m_atoms[b].bonds().size()) return false;
  size_1d atomicA, atomicB, otherA, otherB;
  // append all atoms bound to a / b that are not source atoms from_a / from_b
  for (size_1d::size_type i=0; i<TA; ++i)
  {
    std::size_t const bai = m_atoms[a].bonds()[i]; // bond i of a
    if (bai != from_a) 
    {
      scon::sorted::insert(atomicA, m_atoms[bai].number());
      otherA.push_back(bai);
    }
    std::size_t const bbi(m_atoms[b].bonds()[i]); // bond i of b
    if (bbi != from_b)
    {
      scon::sorted::insert(atomicB, m_atoms[bbi].number());
      otherB.push_back(bbi);
    }
  }
  size_1d::size_type const O = otherA.size();
  std::vector<bool> done_b(O, false);
  // cycle all residues i on A
  for (size_1d::size_type i(0U); i<O; ++i)
  {
    // atomicA and B are sorted, so for the same i the should have the same numbers if equal
    if (atomicA[i] != atomicB[i]) return false;
    // we search for a residue j at B that is equal to residue i at A
    // j has to be marked as "not found yet" (done_b[j] = false)
    bool found = false;
    for (std::vector<bool>::size_type j=0; j<O; ++j)
    {
      if (!done_b[j] && res_is_equal(otherA[i], otherB[j], a, b, depth + 1))
      {
        done_b[j] = true;
        found = true;
        break;
      }
    }
    // if no equivalent j has been found for i 
    // the residues a and b are different
    if (!found) return false;
  }
  return true;
}

static void fix_internal(std::vector<coords::Atom> & atoms, std::size_t intdex)
{
  while (intdex < atoms.size() && intdex > 0)
  {
    if (atoms[intdex].ifix()) break;
    atoms[intdex].ifix(true);
    intdex = atoms[intdex].ibond();
  }
  atoms[0].ifix(true);
}

static void update_fixations(std::vector<coords::Atom> & atoms)
{
  std::size_t const N = atoms.size();
  // unfix all internal fixations
  for (std::size_t i = 0; i < N; ++i) atoms[i].ifix(false);
  // refix internals
  for (std::size_t i = 0; i < N; ++i)
  {
    if (atoms[i].fixed()) fix_internal(atoms, atoms[i].a_to_i());
  }
}

void coords::Atoms::fix_all(bool const fix_it)
{
  for (auto & atom : m_atoms) atom.fix(fix_it);
  update_fixations(m_atoms);
}

void coords::Atoms::fix(std::size_t const atom, bool const fix_it)
{
  if (atom >= m_atoms.size() || m_atoms.empty()) return;
  if (!fix_it)
  {
    m_atoms[atom].fix(false);
    update_fixations(m_atoms);
  }
  else
  {
    m_atoms[atom].fix();
    fix_internal(m_atoms, m_atoms[atom].a_to_i());
  }
}

std::size_t coords::Atoms::intern_of_dihedral(std::size_t a,
  std::size_t b, std::size_t c, std::size_t d) const
{
  std::size_t const n = size();
  std::size_t const di = m_atoms[d].a_to_i();
  std::size_t const diba = m_atoms[m_atoms[di].ibond()].i_to_a();
  if (diba == c)
  {
    std::size_t const diaa = m_atoms[m_atoms[di].iangle()].i_to_a();
    if (diaa == b)
    {
      std::size_t const dida = m_atoms[m_atoms[di].idihedral()].i_to_a();
      if (dida == a) return di;
    }
  }
  return n;
}


void coords::Atoms::swap (Atoms &r)
{
  m_atoms.swap(r.m_atoms);
  m_sub_systems.swap(r.m_sub_systems);
  m_molecules.swap(r.m_molecules);
  main_bond_indices.swap(r.main_bond_indices);
  main_angle_indices.swap(r.main_angle_indices);
  main_torsion_indices.swap(r.main_torsion_indices);
  std::swap(m_sub_in_index, r.m_sub_in_index);
  std::swap(m_sub_out_index, r.m_sub_out_index);
  std::swap(m_sub_io, r.m_sub_io);
  std::swap(m_in_exists, r.m_in_exists);
  std::swap(m_out_exists, r.m_out_exists);
}


/* ######################################################

  class Stereo implementation

   ######  ######## ######## ########  ########  #######  
  ##    ##    ##    ##       ##     ## ##       ##     ## 
  ##          ##    ##       ##     ## ##       ##     ## 
   ######     ##    ######   ########  ######   ##     ## 
        ##    ##    ##       ##   ##   ##       ##     ## 
  ##    ##    ##    ##       ##    ##  ##       ##     ## 
   ######     ##    ######## ##     ## ########  #######  

###################################################### */




coords::Stereo::Stereo (Atoms const &atoms, coords::Representation_3D const &xyz)
  : m_centers()
{
  using scon::cross;
  std::vector<Atom>::size_type const N = atoms.size();
  m_centers.clear();
  for (std::size_t i(0U); i<N; ++i)
  {
    if (atoms.atom(i).bonds().size() == 4)
    {
      pair stereo_pair;
      stereo_pair.m_atom = i;
      // Bond indices
      stereo_pair.bonds[0] = atoms.atom(i).bonds()[0];
      stereo_pair.bonds[1] = atoms.atom(i).bonds()[1];
      stereo_pair.bonds[2] = atoms.atom(i).bonds()[2];
      stereo_pair.bonds[3] = atoms.atom(i).bonds()[3];
      // Two hydrogens? 
      //if (atoms.atom(stereo_pair.bonds[0]).number() == 1 )
      // All of the residues need to differ if we have a stereo center
      if (
        !atoms.res_is_equal(stereo_pair.bonds[0], stereo_pair.bonds[1], i, i, 0) &&
        !atoms.res_is_equal(stereo_pair.bonds[0], stereo_pair.bonds[2], i, i, 0) &&
        !atoms.res_is_equal(stereo_pair.bonds[0], stereo_pair.bonds[3], i, i, 0) &&
        !atoms.res_is_equal(stereo_pair.bonds[1], stereo_pair.bonds[2], i, i, 0) &&
        !atoms.res_is_equal(stereo_pair.bonds[1], stereo_pair.bonds[3], i, i, 0) &&
        !atoms.res_is_equal(stereo_pair.bonds[2], stereo_pair.bonds[3], i, i, 0)
        )
      {
        coords::Cartesian_Point crossv(xyz[stereo_pair.bonds[1]] - xyz[stereo_pair.bonds[2]]);
        crossv = scon::cross(crossv, xyz[stereo_pair.bonds[3]] - xyz[stereo_pair.bonds[2]]);
        double const scalar(dot(xyz[stereo_pair.bonds[0]] - xyz[stereo_pair.bonds[2]], crossv));
        stereo_pair.m_dir = !(scalar < 0.0);
        m_centers.push_back(stereo_pair);
      }
    }
  }
}

void coords::Stereo::update(coords::Representation_3D const &xyz)
{
  for (auto & stereo_pair : m_centers)
  {
    coords::Cartesian_Point crossv(
      cross(xyz[stereo_pair.bonds[1]] - xyz[stereo_pair.bonds[2]],
            xyz[stereo_pair.bonds[3]] - xyz[stereo_pair.bonds[2]])
      );
    float_type const scalar(dot(xyz[stereo_pair.bonds[0]] - 
      xyz[stereo_pair.bonds[2]], crossv));
    stereo_pair.m_dir = !(scalar < 0.0);
  }
}




/* ######################################################


  ########  ####    ###     ######  
  ##     ##  ##    ## ##   ##    ## 
  ##     ##  ##   ##   ##  ##       
  ########   ##  ##     ##  ######  
  ##     ##  ##  #########       ## 
  ##     ##  ##  ##     ## ##    ## 
  ########  #### ##     ##  ######  


###################################################### */




void coords::bias::Potentials::append_config ()
{
  m_distances.insert(m_distances.end(), Config::get().coords.bias.distance.begin(), Config::get().coords.bias.distance.end());
  m_angles.insert(m_angles.end(), Config::get().coords.bias.angle.begin(), Config::get().coords.bias.angle.end());
  m_dihedrals.insert(m_dihedrals.end(), Config::get().coords.bias.dihedral.begin(), Config::get().coords.bias.dihedral.end());
  m_utors.insert(m_utors.end(), Config::get().coords.bias.utors.begin(), Config::get().coords.bias.utors.end());
}

void coords::bias::Potentials::swap(Potentials & rhs)
{
  std::swap(b, rhs.b);
  std::swap(a, rhs.a);
  std::swap(d, rhs.d);
  m_dihedrals.swap(rhs.m_dihedrals);
  m_angles.swap(rhs.m_angles);
  m_distances.swap(rhs.m_distances);
  m_spherical.swap(rhs.m_spherical);
  m_cubic.swap(rhs.m_cubic);
  m_utors.swap(rhs.m_utors);
  m_udist.swap(rhs.m_udist);
}


coords::bias::Potentials::Potentials()
  : b(), a(), d(), s(), c(),
  m_dihedrals(Config::get().coords.bias.dihedral),
  m_angles(Config::get().coords.bias.angle),
  m_distances(Config::get().coords.bias.distance),
  m_spherical(Config::get().coords.bias.spherical),
  m_cubic(Config::get().coords.bias.cubic),
  m_utors(Config::get().coords.bias.utors),
  m_udist(Config::get().coords.bias.udist)
{ }

bool coords::bias::Potentials::empty() const
{
  return scon::empty(m_dihedrals, m_angles, m_distances,
    m_spherical, m_cubic, m_utors, m_udist);
}

double coords::bias::Potentials::apply(Representation_3D const & xyz,
  Gradients_3D & g_xyz, Cartesian_Point const & center)
{
  if (!m_dihedrals.empty()) d = dih(xyz, g_xyz);
  if (!m_angles.empty()) a = ang(xyz, g_xyz);
  if (!m_distances.empty()) b = dist(xyz, g_xyz);
  if (!m_spherical.empty()) s = spherical(xyz, g_xyz, center);
  if (!m_cubic.empty()) c = cubic(xyz, g_xyz, center);
  return b+a+d+s+c;
}


void coords::bias::Potentials::umbrellaapply(Representation_3D const & xyz, 
  Gradients_3D & g_xyz, std::vector<double> &uout){
  if (!m_utors.empty()) umbrelladih(xyz, g_xyz, uout);
  if (!m_udist.empty()) umbrelladist(xyz, g_xyz, uout);
}

static inline double check_angle (double angle)
{
  return angle > 180.0 ? -(360.0 - angle) : (angle < -180.0 ? 360.0 + angle : angle);
}


double coords::bias::Potentials::dih (Representation_3D const &positions, Gradients_3D & gradients)
{
  float_type E(0.0);
  using std::abs;
  using std::atan2;
  for (auto &dih : m_dihedrals)
  {
    std::cout << "Current value: " << dih.value << std::endl;

    coords::r3 const b01(positions[dih.b] - positions[dih.a]);
    coords::r3 const b12(positions[dih.c] - positions[dih.b]);
    coords::r3 const b23(positions[dih.d] - positions[dih.c]);
    coords::r3 const b02(positions[dih.c] - positions[dih.a]);
    coords::r3 const b13(positions[dih.d] - positions[dih.b]);
    coords::r3 const t(cross(b01, b12));
    coords::r3 const u(cross(b12, b23));
    coords::r3 const tu(cross(t, u));

    coords::float_type const tl2(dot(t,t));
    coords::float_type const ul2(dot(u, u));
    coords::float_type const r12(geometric_length(b12));
    //coords::float_type const norm = r12*geometric_length(tu);

    auto ta = coords::angle_type::from_rad(atan2(geometric_length(tu), dot(t, u)));
    
    dih.value = ta;

    auto const torsion((ta - dih.ideal).degrees());

    //auto const tor_rad = torsion*SCON_PI180;

    auto const dE(2.0*dih.force*torsion);
    auto const tE(dih.force*torsion*torsion);


    E += tE;

    std::cout << "Target: " << dih.ideal << ", Is: " << ta << ", delta = " << torsion << ", E = " << tE << "\n";

    // Derivatives
    auto const dt(cross(t, b12) * (dE / (tl2*r12)));
    auto const du(cross(u, b12) * (-dE / (ul2*r12)));

    auto const g1 = cross(dt, b12);
    auto const g2 = cross(b02, dt) + cross(du, b23);
    auto const g3 = cross(dt, b01) + cross(b13, du);
    auto const g4 = cross(du, b12);

    std::cout << "  10 XX " << scon::c3_delimeter(' ') << (positions[dih.a] + g1) << " 0\n";
    std::cout << "  11 XX " << scon::c3_delimeter(' ') << (positions[dih.b] + g2) << " 0\n";
    std::cout << "  12 XX " << scon::c3_delimeter(' ') << (positions[dih.c] + g3) << " 0\n";
    std::cout << "  13 XX " << scon::c3_delimeter(' ') << (positions[dih.d] + g4) << " 0\n";
      
    gradients[dih.a] += cross(dt, b12);
    gradients[dih.b] += cross(b02, dt) + cross(du, b23);
    gradients[dih.c] += cross(dt, b01) + cross(b13, du);
    gradients[dih.d] += cross(du, b12);

    //std::cout << "dihg " << scon::vector_delimeter('\n') << gradients << "\n";

  }
  return E;
}

void coords::bias::Potentials::umbrelladih(Representation_3D const &positions, 
  Gradients_3D & gradients, std::vector<double> &uout) const
{
  using std::abs;
  for (auto const &dih : m_utors)
  {
    float_type dE(0.0), diff(0.0);
    Cartesian_Point const b01(positions[dih.index[1]] - positions[dih.index[0]]);
    Cartesian_Point const b12(positions[dih.index[2]] - positions[dih.index[1]]);
    Cartesian_Point const b23(positions[dih.index[3]] - positions[dih.index[2]]);
    Cartesian_Point const b02(positions[dih.index[2]] - positions[dih.index[0]]);
    Cartesian_Point const b13(positions[dih.index[3]] - positions[dih.index[1]]);
    Cartesian_Point t(cross(b01, b12));
    Cartesian_Point u(cross(b12, b23));
    float_type const tl2(dot(t, t));
    float_type const ul2(dot(u, u));
    float_type const r12(geometric_length(b12));
    Cartesian_Point const tu(cross(t, u));
    float_type torsion = angle(t,u).degrees();
    float_type const norm = r12*geometric_length(tu);
    torsion = (abs(norm) > float_type(0) && (dot(b12, tu) / norm) < 0.0) ? -torsion : torsion;
    // Apply half harmonic bias potential according to torsion value
    if (dih.angle > 0){
      if (torsion > 0){
        uout.push_back(torsion);
        diff = torsion - dih.angle;
        dE = dih.force * diff * RATIO180PI;
      }
      else{
        if (((360 + torsion) > 180) && dih.angle > 90)
        {
          uout.push_back(torsion + 360);
        }
        else uout.push_back(torsion);
        diff = torsion - dih.angle;
        if (diff < -180){
          diff = 360 + diff;
          dE = dih.force * diff * RATIO180PI;
        }
        else{
          dE = dih.force * diff * RATIO180PI;
        }
      }
    }// end of angle > 0
    else if (dih.angle < 0){
      if (torsion < 0){
        uout.push_back(torsion);
        diff = torsion - dih.angle;
        dE = dih.force * diff * RATIO180PI;
      }
      else{
        if (((-360 + torsion) < -180) && dih.angle < -90)
        {
          uout.push_back(-360 + torsion);
        }
        else  uout.push_back(torsion);
        diff = torsion - dih.angle;
        if (diff > 180){
          diff = 360 - diff;
          dE = -dih.force * diff * RATIO180PI;
        }
        else dE = dih.force * diff * RATIO180PI;
      }
    }// end of angle < 0
    else if (dih.angle == 0){
      diff = torsion;
      uout.push_back(torsion);
      dE = dih.force * diff * RATIO180PI;
    }
    // Derivatives
    t = cross(t, b12);
    t *= dE / (tl2*r12);
    u = cross(u, b12);
    u *= -dE / (ul2*r12);
    //std::cout << "U W: " << torsion << ", AW: " << diff << ", SOLL: " << dih.angle;
    //std::cout << "; FORCE: " << dih.force << ", DE: " << dE << ", PW: " << uout.back() << std::endl;
    //std::cout << dih.index[0] << "   " << dih.index[1] << "   " << dih.index[2] << "   " << dih.index[3] << std::endl;
    gradients[dih.index[0]] += cross(t, b12);
    gradients[dih.index[1]] += cross(b02, t) + cross(u, b23);
    gradients[dih.index[2]] += cross(t, b01) + cross(b13, u);
    gradients[dih.index[3]] += cross(u, b12);
    //std::cout << "U GRAD:  " << t.crossd(b12) << "   " << (b02.crossd(t) + u.crossd(b23));
    //std::cout << "   " << (t.crossd(b01) + b13.crossd(u)) << "  " << u.crossd(b12) << std::endl;
  }

}

void coords::bias::Potentials::umbrelladist(Representation_3D const &positions, 
  Gradients_3D & gradients, std::vector<double> &uout) const
{
  for (auto const &dist : m_udist)
  {
    coords::Cartesian_Point bv(positions[dist.index[0]] - positions[dist.index[1]]);
    float_type md = geometric_length(bv);
    uout.push_back(md);
    float_type diff(md - dist.dist);
    //apply half harmonic potential
    float_type dE = dist.force * diff / md;
    Cartesian_Point gv = bv*dE;
    gradients[dist.index[0]] += gv; 
    gradients[dist.index[1]] -= gv;
  }
}

double coords::bias::Potentials::dist (Representation_3D const &positions, Gradients_3D &gradients )
{
  double E(0.0);
  for (auto &distance : m_distances){
    coords::Cartesian_Point bv(positions[distance.a] - positions[distance.b]);
    double md = geometric_length(bv);
    distance.value = md;
    double diff(md - distance.ideal);
    //apply half harmonic potential
    E = 0.5 * diff * distance.force;
    double dE = distance.force * diff / md;
    Cartesian_Point gv = bv*dE;
    gradients[distance.a] += gv; 
    gradients[distance.b] -= gv;
  }
  return E;
}


double coords::bias::Potentials::ang (Representation_3D const &, Gradients_3D & )
{
  return 0.0;
}

double coords::bias::Potentials::spherical(Representation_3D const &positions, 
  Gradients_3D & gradients, Cartesian_Point const & center)
{
  using std::max;
  using std::pow;
  double E(0.0);
  std::size_t const N(positions.size());
  for (auto const & orb : m_spherical)
  {
    if (Config::get().general.verbosity > 29)
    {
      std::cout << "Applying spherical boundary with radius " << orb.radius;
      std::cout << " around " << center << " (exponent = " << orb.exponent << ")\n";
    }
    for (std::size_t i = 0; i < N; ++i)
    {
      Cartesian_Point db = positions[i] - center;
      double const l = geometric_length(db);
      if (l < orb.radius)
      {
        continue;
      }
      double const md = l - orb.radius;
      db *= md / l;
      double const e = max(1.0, orb.exponent);
      double const ene = orb.force * pow(md, e);
      E += ene;
      double dE = orb.force * e * pow(md, e - 1.0);
      db *= dE;
      if (Config::get().general.verbosity > 49)
      {
        std::cout << "Spherical boundary gradient of atom " << i + 1 << " at ";
        std::cout << positions[i] << " (which is " << l << " from " << center;
        std::cout << " making delta = " << md << ") is " << db << " with energy " << ene << lineend;
      }
      gradients[i] += db;
    }
    if (Config::get().general.verbosity > 29)
    {
      std::cout << "E = " << E << lineend;
    }
  }
  return E;
}

double coords::bias::Potentials::cubic(Representation_3D const &positions, 
  Gradients_3D & gradients, Cartesian_Point const & center)
{
  double E(0.0);
  using scon::abs;
  std::size_t const N(positions.size());
  for (auto const & box : m_cubic)
  {
    Cartesian_Point const halfdim(box.dim / 2.);
    if (Config::get().general.verbosity > 29)
    {
      std::cout << "Applying cubic boundary with side length " << box.dim;
      std::cout << " and halfdim " << halfdim << " around " << center;
      std::cout << " (exponent = " << box.exponent << ")\n";
    }
    for (std::size_t i = 0; i < N; ++i)
    {
      Cartesian_Point const db = positions[i] - center;
      Cartesian_Point const delta = scon::abs(db) - scon::abs(halfdim);
      double const expo = std::max(1.0, box.exponent);
      double const fexpo = box.force * expo;
      double energy = double();
      Cartesian_Point de;
      if (delta.x() > 0.)
      {
        energy += box.force * std::pow(delta.x(), expo);
        de.x() = (db.x() < 0. ? -fexpo : fexpo) * std::pow(delta.x(), expo - 1.);
        if (Config::get().general.verbosity > 49)
        {
          std::cout << positions[i].x() << " x out of " << center.x();
          std::cout << " +/- " << halfdim.x() << " by " << delta.x() << " (" << db.x() << ") ";
        }
      }
      if (delta.y() > 0.)
      {
        energy += box.force * std::pow(delta.y(), expo);
        de.y() = (db.y() < 0.?-fexpo:fexpo) * std::pow(delta.y(), expo - 1.);
        if (Config::get().general.verbosity > 49)
        {
          std::cout << positions[i].y() << " y out of " << center.y();
          std::cout << " +/- " << halfdim.y() << " by " << delta.y() << " (" << db.y() << ") ";
        }
      }
      if (delta.z() > 0.)
      {
        energy += box.force * std::pow(delta.z(), expo);
        de.z() = (db.z() < 0.?-fexpo:fexpo) * std::pow(delta.z(), expo - 1.);
        if (Config::get().general.verbosity > 49)
        {
          std::cout << positions[i].z() << " z out of " << center.z();
          std::cout << " +/- " << halfdim.z() << " by " << delta.z() << " (" << db.z() << ") ";
        }
      }
      Cartesian_Point const g = delta*de;
      E += energy;
      if (Config::get().general.verbosity > 49 && energy > 0.)
      {
        std::cout << " E = " << energy << " and grad " << g << lineend;
      }
      gradients[i] += g;
    }
  }
  return E;
}



/* ################################################################################################

 ######   #######   #######  ########  ########  #### ##    ##    ###    ######## ########  ######  
##    ## ##     ## ##     ## ##     ## ##     ##  ##  ###   ##   ## ##      ##    ##       ##    ## 
##       ##     ## ##     ## ##     ## ##     ##  ##  ####  ##  ##   ##     ##    ##       ##       
##       ##     ## ##     ## ########  ##     ##  ##  ## ## ## ##     ##    ##    ######    ######  
##       ##     ## ##     ## ##   ##   ##     ##  ##  ##  #### #########    ##    ##             ## 
##    ## ##     ## ##     ## ##    ##  ##     ##  ##  ##   ### ##     ##    ##    ##       ##    ## 
 ######   #######   #######  ##     ## ########  #### ##    ## ##     ##    ##    ########  ######  

################################################################################################ */


coords::Coordinates::Coordinates() :
  m_atoms(), m_representation(), m_stereo(),
  m_potentials(), m_virial(empty_virial()),
  m_interface(energy::new_interface(this)),
  m_preinterface(energy::pre_interface(this)),
  energy_valid(false),
  NEB_control(false),
  PathOpt_control(false),
  mult_struc_counter(0)
  //, m_topology(0)
{
}

coords::Coordinates::Coordinates(Coordinates &&r) :
  m_atoms(std::move(r.m_atoms)),
  m_representation(std::move(r.m_representation)),
  m_stereo(std::move(r.m_stereo)),
  m_potentials(std::move(r.m_potentials)),
  m_virial(std::move(r.m_virial)),
  m_interface(r.m_interface->move(this)),
  m_preinterface(r.m_preinterface ? r.m_preinterface->move(this) : nullptr),
  energy_valid(r.energy_valid),
  NEB_control(r.NEB_control), PathOpt_control(r.PathOpt_control),
  mult_struc_counter(r.mult_struc_counter)
{
}

coords::Coordinates::Coordinates(Coordinates const &r) :
  m_atoms(r.m_atoms),
  m_representation(r.m_representation),
  m_stereo(r.m_stereo),
  m_potentials(r.m_potentials),
  m_virial(r.m_virial),
  m_interface(r.m_interface->clone(this)),
  m_preinterface(r.m_preinterface ? r.m_preinterface->clone(this) : nullptr),
  energy_valid(false),
  fep(r.fep),
  NEB_control(r.NEB_control),
  PathOpt_control(r.PathOpt_control),
  mult_struc_counter(r.mult_struc_counter)
{
}

coords::Coordinates& coords::Coordinates::operator=(Coordinates const & rhs)
{
  if (this != &rhs)
  {
    Coordinates tmp{ rhs };
    this->swap(tmp);
  }
  return *this;
}

coords::Coordinates& coords::Coordinates::operator=(Coordinates&& rhs)
{
  if (this != &rhs)
  {
    Coordinates tmp{ std::move(rhs) };
    this->swap(tmp);
  }
  return *this;
}

coords::Coordinates::~Coordinates()
{
  if (m_interface) delete m_interface;
  if (m_preinterface) delete m_preinterface;
}



void coords::Coordinates::set_fix(size_t const atom, bool const fix_it)
{
  fix(atom, fix_it);
}


void coords::Coordinates::init_swap_in (Atoms &a, PES_Point &p, bool const update) 
{ 
  if (a.size() != p.size())
  {
    throw std::logic_error("Initialization of Coordinates failed. Invalid sizes.");
  }
  energy_valid = false;
  a.swap(m_atoms); 
  m_atoms.refine();
  p.swap(m_representation); 
  Stereo(m_atoms, m_representation.structure.cartesian).swap(m_stereo);
  m_representation.ia_matrix.resize(subsystems().size());
  for (auto & interact : m_representation.ia_matrix)
  {
    interact.grad.resize(size());
    interact.energy = 0.0;
  }
  if (update) 
  {
    std::size_t const N = m_atoms.size(), M = m_atoms.mains().size();
    m_representation.structure.cartesian.resize(N);
    m_representation.gradient.cartesian.resize(N);
    m_representation.structure.intern.resize(N);
    m_representation.gradient.intern.resize(N);
    m_representation.structure.main.resize(M);
    m_representation.gradient.main.resize(M);
    //m_atoms.c_to_i(m_representation); 
    m_interface->update(false);
    if (m_preinterface) m_preinterface->update(false);
  }
}


void coords::Coordinates::init_in (Atoms a, PES_Point p, bool const update) 
{ 
  init_swap_in(a, p, update);
}

struct OCB
{
  coords::Coordinates * cp;
  OCB(coords::Coordinates & coordpointer) : cp(&coordpointer) {}
  float operator() (scon::vector<scon::c3<float>> const & v, 
    scon::vector<scon::c3<float>> & g, std::size_t const S, bool & go_on)
  {
    cp->set_xyz(coords::Representation_3D(v.begin(), v.end()), false);
    float E = float(cp->g());
    go_on = cp->integrity();
    g.resize(cp->g_xyz().size());
    scon::explicit_transform(cp->g_xyz(), g);
    if (Config::get().general.verbosity > 19)
    {
      std::cout << "Optimization: Energy of step " << S;
      std::cout << " is " << E << " integrity " << go_on << '\n';
    }
    return E;
  }
};

coords::float_type coords::Coordinates::lbfgs ()
{
  using namespace  optimization::local;
  typedef coords::Container<scon::c3<float>> nc3_type;
  // Create optimizer
  auto optimizer = make_lbfgs(
    make_more_thuente(Coords_3d_float_callback(*this))
  );
  optimizer.ls.config.ignore_callback_stop = true;
  // Create Point
  using op_type = decltype(optimizer);
  op_type::point_type x(nc3_type(xyz().begin(), xyz().end()));
  // Optimize point
  optimizer.config.max_iterations = 
    Config::get().optimization.local.bfgs.maxstep;
  optimizer.config.epsilon = 
    (float)Config::get().optimization.local.bfgs.grad;
  optimizer(x);
  // Get structure, gradients and energy into coords
  m_representation.energy = optimizer.p().f;
  m_representation.structure.cartesian = 
    coords::Representation_3D(optimizer.p().x.begin(), optimizer.p().x.end());
  //g();
  //std::cout << "Ene = " << m_representation.energy << "\n";
  m_representation.gradient.cartesian = 
    coords::Gradients_3D(optimizer.p().g.begin(), optimizer.p().g.end());
  // Output
  if (Config::get().general.verbosity > 19 || 
    (optimizer.state() < 0 && Config::get().general.verbosity > 14))
  {
    std::cout << "Optimization done (status " << optimizer.state() << 
      "). Evaluations:" << optimizer.iter() << lineend;
  }
  if (Config::get().general.verbosity > 19 && integrity())
  {
    std::cout << "Energy after optimization: " << lineend;
    e_head_tostream_short(std::cout, energyinterface());
    e_tostream_short(std::cout, energyinterface());
  }
  // Return floating point
  return optimizer.p().f;
}

double coords::Coordinates::prelbfgs()
{
  using namespace  optimization::local;
  typedef coords::Container<scon::c3<float>> nc3_type;
  // Create optimizer
  auto optimizer = make_lbfgs(
    make_more_thuente(Coords_3d_float_pre_callback(*this))
  );
  // Create Point
  using op_type = decltype(optimizer);
  op_type::point_type x(nc3_type(xyz().begin(), xyz().end()));
  // Optimize point
  optimizer(x);
  m_representation.energy = x.f;
  m_representation.structure.cartesian = 
    coords::Representation_3D(x.x.begin(), x.x.end());
  m_representation.gradient.cartesian = 
    coords::Gradients_3D(x.g.begin(), x.g.end());
  if (Config::get().general.verbosity > 19)
  {
    std::cout << "Optimization done (status " << optimizer.state() << 
      "). Evaluations:" << optimizer.iter() << lineend;
  }
    
  if (Config::get().general.verbosity > 19 && m_interface->intact())
  {
    std::cout << "Energy after optimization: " << lineend;
    e_head_tostream_short(std::cout, m_interface);
    e_tostream_short(std::cout, m_interface);
  }
  return x.f;
}

coords::Gradients_Main coords::Coordinates::dimermethod_dihedral
  (std::vector<coords::Gradients_Main> const &D)
{
  if (this->atoms().mains().empty())
  {
    throw std::runtime_error("System does not contain any main dihedrals. Dimermethod cannot be applied.");
  }
  using Move_T = optimization::CG_DimerMover < Main_Callback > ;
  Main_Callback C(*this);
  Move_T mover(0.1, 1000u);
  coords::Gradients_Main structure;
  std::size_t const N = m_representation.structure.main.size();
  structure.reserve(N);
  for (std::size_t i(0u); i < N; ++i)
  {
    structure.emplace_back(m_representation.structure.main[i].radians());
  }
  Move_T::minimum_type minimum(structure);
  minimum.directions = D;
  C = mover(minimum, C);
  return minimum.directions.back();
}


void coords::Coordinates::swap(Coordinates &rhs) // object swap
{
  m_atoms.swap(rhs.m_atoms);
  m_representation.swap(rhs.m_representation);
  m_stereo.swap(rhs.m_stereo);
  m_potentials.swap(rhs.m_potentials);
  m_virial.swap(rhs.m_virial);
  if (m_interface && rhs.m_interface) m_interface->swap(*rhs.m_interface);
  if (m_preinterface && rhs.m_preinterface) m_preinterface->swap(*rhs.m_preinterface);
  using std::swap;
  //m_sub_interaction.swap(rhs.m_sub_interaction);
  swap(energy_valid, rhs.energy_valid);
  swap(this->fep, rhs.fep);
  swap(this->hessian_def, rhs.hessian_def);
  swap(this->use_fep, rhs.use_fep);
  swap(this->mult_struc_counter, rhs.mult_struc_counter);
  swap(this->NEB_control, rhs.NEB_control);
  swap(this->orthogonalize, rhs.orthogonalize);
  swap(this->PathOpt_control, rhs.PathOpt_control);
}

coords::Cartesian_Point coords::Coordinates::center_of_mass() const
{
  coords::Cartesian_Point COM;
  std::vector<Atom>::size_type const N(m_atoms.size());
  coords::float_type M(0.0);
  for (std::vector<Atom>::size_type i(0U); i<N; ++i)
  {
    coords::float_type const mass(m_atoms.atom(i).mass());
    M += mass;
    COM += xyz(i)*mass;
  }
  COM /= M;
  return COM;
}


coords::Cartesian_Point coords::Coordinates::center_of_geometry () const
{
  std::size_t const N = xyz().size();
  coords::Cartesian_Point p(0);
  for (std::size_t i(0u); i < N; ++i)
  {
    p += xyz()[i];
  }
  p /= float_type(N);
  return p;
}


double coords::Coordinates::weight () const
{
  std::vector<Atom>::size_type const N(m_atoms.size());
  coords::float_type M(0.0);
  for (std::vector<Atom>::size_type i(0U); i<N; ++i)
  {
    M += m_atoms.atom(i).mass();
  }
  return M;
}


void coords::Coordinates::e_head_tostream_short(std::ostream &strm, 
  energy::interface_base const * const ep) const
{
  if (ep) ep->print_E_head(strm);
  else m_interface->print_E_head(strm);
  if (!m_potentials.empty())
  {
    strm << "Bias Potentials: " << lineend;
    strm << std::setw(24) << "DIH";
    strm << std::setw(24) << "ANG";
    strm << std::setw(24) << "DIST";
    strm << std::setw(24) << "SPHERICAL";
    strm << std::setw(24) << "CUBIC" << lineend;
    strm << std::setw(24) << m_potentials.dihedrals().size();
    strm << std::setw(24) << m_potentials.angles().size();
    strm << std::setw(24) << m_potentials.distances().size();
    strm << std::setw(24) << m_potentials.sphericals().size();
    strm << std::setw(24) << m_potentials.cubic().size() << lineend;
  }
}


void coords::Coordinates::e_tostream_short(std::ostream &strm, 
  energy::interface_base const * const ep) const
{
  if (ep) ep->print_E_short(strm);
  else m_interface->print_E_short(strm);
  if (!m_potentials.empty())
  {
    strm << "Bias Energies: " << lineend;
    strm << std::setw(24) << std::fixed << std::setprecision(8) << m_potentials.e_dihedral();
    strm << std::setw(24) << std::fixed << std::setprecision(8) << m_potentials.e_angle();
    strm << std::setw(24) << std::fixed << std::setprecision(8) << m_potentials.e_dist();
    strm << std::setw(24) << std::fixed << std::setprecision(8) << m_potentials.e_spherical();
    strm << std::setw(24) << std::fixed << std::setprecision(8) << m_potentials.e_cubic() << lineend;
  }
  strm << lineend;
}


coords::Cartesian_Point coords::Coordinates::center_of_mass_mol (std::size_t index) const
{
  if (index >= m_atoms.molecules().size())
  {
    throw std::logic_error("Wrong molecule index for center_of_mass_mol() given.");
  }
  coords::Cartesian_Point COM = coords::Cartesian_Point();
  size_2d::size_type const N(m_atoms.molecules()[index].size());
  coords::float_type M = coords::float_type();
  for (std::vector<Atom>::size_type i(0U); i<N; ++i)
  {
    coords::float_type mass(m_atoms.atom(m_atoms.molecules(index, i)).mass());
    M += mass;
    COM += xyz(m_atoms.molecules(index, i))*mass;
  }
  return COM/M;
}


void coords::Coordinates::set_dih(size_type const int_index, 
  coords::angle_type const target_angle, 
  bool const move_dependants_along, bool const move_fixed_dih)
{
  if (int_index < size() && (!m_atoms.atom(int_index).ifix() || move_fixed_dih))
  {
    angle_type rot = target_angle - intern(int_index).azimuth();
    if (move_dependants_along)
    {
      std::size_t int_bond = m_atoms.atom(int_index).ibond();
      size_type const N = m_atoms.atom(int_bond).bound_internals().size();
      for (size_type i(0u); i<N; ++i)
      {
        m_representation.structure.intern[m_atoms.atom(int_bond).bound_internals()[i]].azimuth() += rot;
      }
    }
    else
    {
      m_representation.structure.intern[int_index].azimuth() = target_angle;
    }
  }
}


void coords::Coordinates::rotate_dih(size_type const int_index, 
  coords::angle_type const rot_angle, 
  bool const move_dependants_along, bool const move_fixed_dih)
{
  if (int_index < size() && (!m_atoms.atom(int_index).ifix() || move_fixed_dih))
  {
    if (move_dependants_along)
    {
      // If we rotate all we need to rotate every dihedral of every atom attached to the
      // atom to which the selected one is attached
      std::size_t int_bond = m_atoms.atom(int_index).ibond();
      size_type const N = m_atoms.atom(int_bond).bound_internals().size();
      for (size_type i(0u); i<N; ++i)
      {
        m_representation.structure.intern[m_atoms.atom(int_bond).bound_internals()[i]].azimuth() += rot_angle;
      }
    }
    else
    {
      m_representation.structure.intern[int_index].azimuth() += rot_angle;
    }
  }
}


void coords::Coordinates::rotate_main(size_type const main_index, 
  coords::angle_type const rot_angle, 
  bool const move_dependants_along, bool const move_fixed_dih)
{
  rotate_dih(m_atoms.intern_of_main_idihedral(main_index), rot_angle, move_dependants_along, move_fixed_dih);
}


void coords::Coordinates::set_all_main(coords::Representation_Main const & new_values, 
  bool const apply_to_xyz, bool const move_dependants_along, bool const move_fixed_dih)
{
  size_type const N(main().size());
  if (new_values.size() != N)
  {
    throw std::logic_error("set_all_main called with wrong-sized vector");
  }
  for (size_type i(0u); i<N; ++i)
  {
    set_dih(m_atoms.intern_of_main_idihedral(i), new_values[i], 
      move_dependants_along, move_fixed_dih);
    m_representation.structure.main[i] = new_values[i];
  }
  if (apply_to_xyz) to_xyz();
}

void coords::Coordinates::periodic_boxjump()
{
  std::size_t const N(molecules().size());
  Cartesian_Point const halfbox(Config::get().energy.pb_box/2.0);
  for (std::size_t i = 0; i < N; ++i)
  {
    Cartesian_Point tmp_com(-center_of_mass_mol(i));
    //tmp_com /= halfbox;
    tmp_com.x() = (std::abs(tmp_com.x()) > halfbox.x()) ? tmp_com.x() 
      / Config::get().energy.pb_box.x() : float_type(0.);
    tmp_com.y() = (std::abs(tmp_com.y()) > halfbox.y()) ? tmp_com.y() 
      / Config::get().energy.pb_box.y() : float_type(0.);
    tmp_com.z() = (std::abs(tmp_com.z()) > halfbox.z()) ? tmp_com.z() 
      / Config::get().energy.pb_box.z() : float_type(0.);
    round(tmp_com); 
    tmp_com *= Config::get().energy.pb_box;
    for (auto const atom : molecules(i)) move_atom_by(atom, tmp_com, true);
	}
}


bool coords::Coordinates::validate_bonds() const
{
  std::size_t const N(m_atoms.size());
  for (std::size_t i = 0; i < N; ++i)
  {
    for (auto const & bound : m_atoms.atom(i).bonds())
    {
      double const L(geometric_length(xyz(i) - xyz(bound)));
      if (L < 0.3 || L > 5.0) return false;
    }
  }
  return true;
}


void coords::Coordinates::set_pes(PES_Point const & pes_point, 
  bool const overwrite_fixed /*= false*/)
{
  if (!m_representation.same_size(pes_point))
  {
    std::cout << "CS: " << (m_representation.structure.cartesian.size() == 
      pes_point.structure.cartesian.size()) << "\n";
    std::cout << "CG: " << (m_representation.gradient.cartesian.size() == 
      pes_point.gradient.cartesian.size()) << "\n";
    std::cout << "IS: " << (m_representation.structure.intern.size() == 
      pes_point.structure.intern.size()) << "\n";
    std::cout << "IG: " << (m_representation.gradient.intern.size() == 
      pes_point.gradient.intern.size()) << "\n";
    std::cout << "MS: " << (m_representation.structure.main.size() == 
      pes_point.structure.main.size()) << "\n";
    std::cout << "MG: " << (m_representation.gradient.main.size() == 
      pes_point.gradient.main.size()) << "\n";
    throw std::logic_error("Invalid PES size for coordinates.");
  }
  if (overwrite_fixed)
  {
    m_representation = pes_point;
    m_stereo.update(xyz());
  }
  else
  {
    std::size_t const N(size());
    m_representation.gradient.cartesian.resize(pes_point.gradient.cartesian.size());
    m_representation.structure.cartesian.resize(pes_point.structure.cartesian.size());
    for (std::size_t i(0U); i < N; ++i)
    {
      if (atoms(i).fixed())
      {
        m_representation.gradient.cartesian[i] = coords::cartesian_gradient_type();
      }
      else
      {
        m_representation.structure.cartesian[i] = pes_point.structure.cartesian[i];
        m_representation.gradient.cartesian[i] = pes_point.gradient.cartesian[i];
      }
      if (atoms(i).ifix())
      {
        m_representation.gradient.intern[i] = coords::internal_gradient_type();
      }
      else
      {
        m_representation.gradient.intern[i] = pes_point.gradient.intern[i];
        m_representation.structure.intern[i] = pes_point.structure.intern[i];
      }
    }
    std::size_t const MG(pes_point.gradient.main.size()), 
      MS(pes_point.structure.main.size());
    m_representation.gradient.main.resize(MG);
    m_representation.structure.main.resize(MS);
    for (std::size_t i(0U); i < MG && i < MS; ++i)
    {
      std::size_t const ii = atoms().intern_of_main_idihedral(i);
      if (atoms(ii).ifix())
      {
        m_representation.gradient.main[i] = coords::main_gradient_type();
      }
      else
      {
        m_representation.gradient.main[i] = pes_point.gradient.main[i];
        m_representation.structure.main[i] = pes_point.structure.main[i];
      }
    }
    m_stereo.update(xyz());
  }
}


bool coords::Coordinates::check_superposition_xyz(Representation_3D const &a, Representation_3D const &b, double const x /*= 0.35*/) const
{
  std::size_t const N(a.size());
  if (size() != N || b.size() != N)
    throw std::logic_error("Evaluating the structural overlap between different sized structures.");
  for (std::size_t i(0u); i < N; ++i)
  {
    bool is_superposed(false);
    for (std::size_t j(0u); j < N; ++j)
    {
      if (j != i && m_atoms.atom(i).number() == m_atoms.atom(j).number())
      {
        if (scon::geometric_length(a[i] - b[j]) < x) is_superposed = true;
      }
    }
    if (!is_superposed) return false;
  }
  return true;
}

void coords::Coordinates::set_pes(PES_Point && pes_point,
  bool const overwrite_fixed /*= false*/)
{
  bool const valid_size = scon::equal_size_ranges(pes_point.structure.cartesian,
    pes_point.gradient.cartesian,
    pes_point.structure.intern,
    pes_point.gradient.intern,
    m_representation.structure.cartesian,
    m_representation.gradient.cartesian,
    m_representation.structure.intern,
    m_representation.gradient.intern) && 
    scon::equal_size_ranges(pes_point.structure.main, 
    m_representation.structure.main,
    pes_point.gradient.main,
    m_representation.gradient.main);

  if (!valid_size)
  {
    throw std::logic_error("Invalid moved PES size for coordinates.");
  }
  m_representation.swap(pes_point);
  m_representation.resize(pes_point.structure.cartesian.size(), 
    pes_point.structure.main.size());
  if (!overwrite_fixed)
  {
    std::size_t const N(m_representation.structure.cartesian.size());
    for (std::size_t i(0U); i < N; ++i)
    {
      if (atoms(i).fixed())
      {
        m_representation.structure.cartesian[i] = pes_point.structure.cartesian[i];
        m_representation.gradient.cartesian[i] = coords::cartesian_gradient_type();
      }
      if (atoms(i).ifix())
      {
        m_representation.structure.intern[i] = pes_point.structure.intern[i];
        m_representation.gradient.intern[i] = coords::internal_gradient_type();
      }
    }
    std::size_t const M(pes_point.structure.main.size());
    for (std::size_t i(0U); i < M; ++i)
    {
      std::size_t const ii = atoms().intern_of_main_idihedral(i);
      if (atoms(ii).ifix())
      {
        m_representation.structure.main[i] = pes_point.structure.main[i];
        m_representation.gradient.main[i] = coords::main_gradient_type();
      }
    }
  }
  m_stereo.update(xyz());
}

coords::float_type coords::Internal_Callback::operator() 
  (scon::vector<coords::float_type> const & v, 
    scon::vector<coords::float_type>& g, std::size_t const S, bool & go_on)
{
  std::size_t i = 0;
  coords::Representation_Internal rin(v.size() / 3);
  for (auto & e : rin)
  {
    e.radius() = v[i++];
    e.inclination() = coords::angle_type(v[i++]);
    e.azimuth() = coords::angle_type(v[i++]);
  }
  cp->set_internal(rin);
  cp->to_xyz();
  float E = float(cp->g());
  cp->to_internal();
  go_on = cp->integrity();
  g.resize(v.size());
  i = 0;
  for (auto const & e : cp->g_intern())
  {
    g[i++] = e.x();
    g[i++] = e.y();
    g[i++] = e.z();
  }
  if (Config::get().general.verbosity > 19)
  {
    std::cout << "Optimization: Energy of step " << S;
    std::cout << " is " << E << " integrity " << go_on << lineend;
  }
  return E;
}

coords::Gradients_Internal coords::Internal_Callback::from(
  coords::Representation_Internal const & v)
{
  std::size_t const N = v.size();
  coords::Gradients_Internal r;
  r.reserve(N);
  for (std::size_t i(0u); i < N; ++i)
  {
    r.emplace_back(v[i].radius(), 
      v[i].inclination().radians(),
      v[i].azimuth().radians());
  }
  return r;
}
coords::Representation_Internal coords::Internal_Callback::to(
  scon::vector<scon::c3<coords::float_type>> const & v)
{
  std::size_t const N = v.size();
  coords::Representation_Internal r;
  r.reserve(N);
  for (std::size_t i(0u); i < N; ++i)
  {
    r.emplace_back(v[i].x(),
      scon::ang<coords::float_type>::from_rad(v[i].y()),
      scon::ang<coords::float_type>::from_rad(v[i].z()));
  }
  return r;
}

coords::float_type coords::Main_Callback::operator() (coords::Gradients_Main const & v,
  coords::Gradients_Main & g, std::size_t const S, bool & go_on)
{
  std::size_t const N = v.size();
  coords::Representation_Main r;
  r.reserve(N);
  for (std::size_t i(0u); i < N; ++i)
  {
    r.emplace_back(coords::main_type::from_rad(v[i]));
  }
  cp->set_all_main(r);
  auto E = cp->g();
  cp->to_internal();
  go_on = cp->integrity();
  g = cp->g_main();
  if (Config::get().general.verbosity > 19)
  {
    std::cout << "Optimization: Energy of step " << S;
    std::cout << " is " << E << " integrity " << go_on << lineend;
  }
  return E;
}

void coords::cartesian_logfile_drain::operator() (coords::Representation_3D && xyz)
{
  if (!cp || !strm) return;
  // Save current state
  auto const tmp = (*cp).xyz();
  // plug snapshot into coords
  (*cp).set_xyz(std::move(xyz));
  // optimize snapshot
  if (opt) { (*cp).o(); }
  // Print to stream
  *strm << *cp;
  // reset state
  (*cp).set_xyz(std::move(tmp));
}

coords::offset_buffered_cartesian_logfile coords::make_buffered_cartesian_log(Coordinates &c,
  std::string file_suffix, std::size_t buffer_size,
  std::size_t log_offset, bool optimize)
{
  return scon::offset_call_buffer<coords::Representation_3D>(buffer_size, log_offset,
    cartesian_logfile_drain{ c, output::filename(file_suffix).c_str(), optimize });
}

float coords::Coords_3d_float_pre_callback::operator() (scon::vector<scon::c3<float>> const & v,
  scon::vector<scon::c3<float>> & g, std::size_t const S, bool & go_on)
{
  cp->set_xyz(coords::Representation_3D(v.begin(), v.end()));
  float E = float(cp->pg());
  go_on = cp->integrity();
  g = scon::vector<scon::c3<float>>(cp->g_xyz().begin(), cp->g_xyz().end());
  if (Config::get().general.verbosity > 19)
    std::cout << "Optimization: Energy of step " << S << " is " << E << " integrity " << go_on << lineend;
  return E;
}

coords::Representation_3D coords::Coords_3d_float_callback::to(scon::vector<scon::c3<float>> const & v)
{
  return coords::Representation_3D(v.begin(), v.end());
}

scon::vector<scon::c3<float>> coords::Coords_3d_float_callback::from(coords::Gradients_3D const & g)
{
  scon::vector<scon::c3<float>> r(g.size());
  std::transform(g.begin(), g.end(), r.begin(),
    [](coords::Cartesian_Point const &p) -> scon::c3<float>
  {
    return scon::c3<float>(static_cast<float>(p.x()),
      static_cast<float>(p.y()), static_cast<float>(p.z()));
  });
  return r;
}


float coords::Coords_3d_float_callback::operator() (scon::vector<scon::c3<float>> const & v,
  scon::vector<scon::c3<float>> & g, std::size_t const S, bool & go_on)
{
  cp->set_xyz(to(v), false);
  float E = float(cp->g());
  go_on = cp->integrity();
  g = from(cp->g_xyz());
  //*ls << *cp;
  if (Config::get().general.verbosity > 19)
  {
    std::cout << "Optimization: Energy of step " << S;
    std::cout << " is " << E << " integrity " << go_on << lineend;
    //std::cout << "totg " << scon::vector_delimeter('\n') << g << "\n";
  }
  return E;
}


//##################################################################################
//##################################################################################
//                    SMOOTH PARTICLE MESH EWALD
//##################################################################################
//##################################################################################

//void coords::Coordinates::pme_stuff(int atoms)
//{
//  pme.pmetemp.natoms = atoms;
//  double x, nenner;
//  std::vector<double> stupidarray(250);
//  std::vector<double> brecurs(6);
//  double eps = 1e-8;
//  double tempcoeff;
//  double rate, j, l, blow, bhigh;
//  //Define PME grid size. Must be even number with prime factores 2, 3 and 5
//  // get even numbers factorized by 2, 3 and 5
//  for (int i = 2; i < 600; i++)
//  {
//    int n = i;
//    // check if numnber is even
//    if (n % 2 != 0) continue;
//    // divide by 2 as long as possible
//    while (n % 2 == 0)
//    {
//      n = n / 2;
//    }
//    // get prime factors 3 and 5
//    for (int j = 3; j <= 5; j = j + 2)
//    {
//      // divide by 3 as long as possible, then by 5 as long as possible
//      while (n%j == 0)
//      {
//        n = n / j;
//      }
//    }
//    // if resulting number is 1 push number i into vector
//    if (n == 1) pme.pmetemp.gridnumbers.push_back(i);
//  }
//  //Get initial Ewald coefficient based on cutoff radius
//  rate = eps + 1.0;
//  pme.pmetemp.ewaldcoeff = 0.5;
//  j = 0;
//  // trial and error function for first initial guess
//  while (rate > eps){
//    j += 1;
//    pme.pmetemp.ewaldcoeff *= 2.0;
//    l = pme.pmetemp.ewaldcoeff * Config::get().energy.cutoff;
//    rate = erfc(l) / Config::get().energy.cutoff;
//  }
//  //run binary search according to Sander program to refine coefficient
//  blow = 0;
//  bhigh = pme.pmetemp.ewaldcoeff;
//  for (int i = 1; i <= 50; i++)
//  {
//    pme.pmetemp.ewaldcoeff = (blow + bhigh) / 2;
//    if (erfc(pme.pmetemp.ewaldcoeff*Config::get().energy.cutoff) / Config::get().energy.cutoff > eps)
//    {
//      blow = pme.pmetemp.ewaldcoeff;
//    }
//    else bhigh = pme.pmetemp.ewaldcoeff;
//  }
//  tempcoeff = pme.pmetemp.ewaldcoeff;
//  //Calculate the grid size from the periodic box dimensions
//  pme.pmetemp.treshold = 1e-8;
//  pme.pmetemp.dim1 = (Config::get().energy.pb_box.x() * 1.2 - pme.pmetemp.treshold) + 1;
//  pme.pmetemp.dim2 = (Config::get().energy.pb_box.y() * 1.2 - pme.pmetemp.treshold) + 1;
//  pme.pmetemp.dim3 = (Config::get().energy.pb_box.z() * 1.2 - pme.pmetemp.treshold) + 1;
//  // check that grid size is even-numbered for effective fftw
//  pme.pmetemp.nxpoints = pme.pmetemp.nypoints = pme.pmetemp.nzpoints = pme.pmetemp.fftgridmax;
//  // get grid size in each dimension
//  for (int i = pme.pmetemp.gridnumbers.size(); i >= 1; i--){
//    int temp;
//    temp = pme.pmetemp.gridnumbers[i];
//    if (temp < 250)
//    {
//      if (temp >= pme.pmetemp.dim1) pme.pmetemp.nxpoints = temp;
//      if (temp >= pme.pmetemp.dim2) pme.pmetemp.nypoints = temp;
//      if (temp >= pme.pmetemp.dim3) pme.pmetemp.nzpoints = temp;
//    }
//  }
//  // compute the moduli for the inverse discrete fourier transformation
//  // ########
//  // B-Spline Stuff
//  // ########
//  x = 0.0;
//  brecurs[1] = (1.0 - x);
//  brecurs[2] = x;
//  // recursive loop over the b-spline order 1 to x (here 5) to get spline coefficients
//  for (int k = 3; k <= pme.pmetemp.bsplineorder; k++)
//  {
//    nenner = 1.0 / (k - 1);
//    brecurs[k] = x * brecurs[k - 1] * nenner;
//    for (int i = 1; i <= (k - 2); i++)
//    {
//      brecurs[k - i] = ((x + i)*brecurs[k - i - 1] + (k - i - x)*brecurs[k - i]) * nenner;
//    }
//    brecurs[1] = (1.0 - x) * brecurs[1] * nenner;
//  }
//  // ########
//  // Fill spline array
//  // #######
//  for (int i = 0; i < 250; i++)
//  {
//    stupidarray[i] = 0.0;
//  }
//  // i = 2 till b-spline order + 1
//  for (int i = 2; i <= pme.pmetemp.bsplineorder + 1; i++)
//  {
//    stupidarray[i] = brecurs[i - 1];
//  }
//  // Calculate the modulus for the spline arrays
//  coords::Coordinates::pme_dftmodulus(stupidarray);
//  // Set chunks for parallel execution
//#ifdef _OPENMP
//  coords::Coordinates::roughgrid();
//#endif
//  // Allocate memory for needed arrays
//  pme.pmetemp.bscx.Allocate(4, pme.pmetemp.bsplineorder, pme.pmetemp.natoms);
//  pme.pmetemp.bscy.Allocate(4, pme.pmetemp.bsplineorder, pme.pmetemp.natoms);
//  pme.pmetemp.bscz.Allocate(4, pme.pmetemp.bsplineorder, pme.pmetemp.natoms);
//  size_t align = sizeof(Complex);
//  pme.pmetemp.charges.Allocate(pme.pmetemp.nxpoints, pme.pmetemp.nypoints, pme.pmetemp.nzpoints, align);
//  pme.pmetemp.fractionalcharges.Allocate(pme.pmetemp.nxpoints, pme.pmetemp.nypoints, pme.pmetemp.nzpoints);
//  pme.pmetemp.recivectors.Allocate(3, 3);
//  pme.pmetemp.initgrid.Allocate(3, pme.pmetemp.natoms);
//  pme.pmetemp.parallelpme.Allocate(pme.pmetemp.natoms, pme.pmetemp.rgridtotal);
//  pme.pmetemp.feinf.resize(pme.pmetemp.natoms);
//
//  // set up FFTW variables if OPENMP is used
//#ifdef _OPENMP
//  int threads;
//#pragma omp parallel
//  {
//    threads = omp_get_num_threads();
//  }
//  fftw_init_threads();
//  fftw_plan_with_nthreads(threads);
//#endif
//  // Assign memory for charge grid
//  pme.pmetemp.in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* (pme.pmetemp.nxpoints*pme.pmetemp.nypoints*pme.pmetemp.nzpoints));
//  // Set up FFTW plans
//  pme.pmetemp.forward = fftw_plan_dft_3d(pme.pmetemp.nxpoints, pme.pmetemp.nypoints, pme.pmetemp.nzpoints, pme.pmetemp.in, pme.pmetemp.in, FFTW_FORWARD, FFTW_PATIENT);
//  pme.pmetemp.backward = fftw_plan_dft_3d(pme.pmetemp.nxpoints, pme.pmetemp.nypoints, pme.pmetemp.nzpoints, pme.pmetemp.in, pme.pmetemp.in, FFTW_BACKWARD, FFTW_PATIENT);
//  // Calculate reciprocal lattice vectors
//  double boxvolume = Config::get().energy.pb_box.x() * Config::get().energy.pb_box.y() * Config::get().energy.pb_box.z();
//  double x1, x2, x3, y1, y2, y3, z1, z2, z3;
//  x1 = Config::get().energy.pb_box.x();
//  x2 = x3 = 0.0;
//  y1 = y3 = 0.0;
//  y2 = Config::get().energy.pb_box.y();
//  z1 = z2 = 0.0;
//  z3 = Config::get().energy.pb_box.z();
//  pme.pmetemp.recivectors(0, 0) = (y2*z3 - z2*y3) / boxvolume;
//  pme.pmetemp.recivectors(1, 0) = (y3*z1 - z3*y1) / boxvolume;
//  pme.pmetemp.recivectors(2, 0) = (y1*z2 - z1*y1) / boxvolume;
//  pme.pmetemp.recivectors(0, 1) = (z2*x3 - x2*z3) / boxvolume;
//  pme.pmetemp.recivectors(1, 1) = (z3*x1 - z1*x3) / boxvolume;
//  pme.pmetemp.recivectors(2, 1) = (z1*x2 - z2*x1) / boxvolume;
//  pme.pmetemp.recivectors(0, 2) = (x2*y3 - y2*x3) / boxvolume;
//  pme.pmetemp.recivectors(1, 2) = (x3*y1 - y3*x1) / boxvolume;
//  pme.pmetemp.recivectors(2, 2) = (x1*y2 - y1*x2) / boxvolume;
//  // set up charge array
//  pme.pmetemp.atomcharges.resize(pme.pmetemp.natoms);
//  std::ifstream ifs2;
//  ifs2.open(Config::get().general.paramFilename.c_str());
//  int iterator = 0;
//  char buffer[150];
//  std::string bufferc, tempname;
//  std::vector<std::string> tokens;
//  std::string bla;
//  bla = Config::get().general.paramFilename.substr(0, 5);
//  pme.pmetemp.elecfac = 1.0;
//  if (bla.substr(0, 3) == "cha") pme.pmetemp.elecfac = 332.0716;
//  else if (bla.substr(0, 3) == "opl") pme.pmetemp.elecfac = 332.06;
//  //else if (bla.substr(0, 3) == "amb") pme.pmetemp.elecfac = 1.0;
//  //else if (bla.substr(0, 3) == "amo") pme.pmetemp.elecfac
//}
//
//// get FEP data for the three different grids (IN, OUT, ALL) and set flags for total charge vector
//void coords::Coordinates::getfepinfo()
//{
//  std::ifstream ifs;
//  int iterator = 0;
//  char buffer[150];
//  std::string bufferc;
//  std::vector<std::string> tokens;
//
//  ifs.open(Config::get().general.inputFilename.c_str());
//  while (!ifs.eof())
//  {
//    ifs.getline(buffer, 150);
//    bufferc = buffer;
//    tokens.clear();
//    std::istringstream iss(bufferc);
//    std::copy(std::istream_iterator <std::string>(iss), std::istream_iterator <std::string>(), std::back_inserter <std::vector <std::string >>(tokens));
//    int g = tokens.size();
//    if (tokens.size() > 2)
//    {
//
//      if (tokens[g - 1] == "IN")
//      {
//        pme.pmetemp.fepi.push_back(atoi(tokens[0].c_str()));
//        pme.pmetemp.feinf[iterator].flag = 1;
//      }
//      else if (tokens[g - 1] == "OUT")
//      {
//        pme.pmetemp.fepo.push_back(atoi(tokens[0].c_str()));
//        pme.pmetemp.feinf[iterator].flag = 0;
//      }
//      else
//      {
//        pme.pmetemp.fepi.push_back(atoi(tokens[0].c_str()));
//        pme.pmetemp.fepo.push_back(atoi(tokens[0].c_str()));
//        pme.pmetemp.fepa.push_back(atoi(tokens[0].c_str()));
//        pme.pmetemp.feinf[iterator].flag = 2;
//      }
//      iterator += 1;
//    }
//
//  }
//}
//
//// Precalculate the DFT moduli
//void coords::Coordinates::pme_dftmodulus(std::vector<double> & stupidarray)
//{
//  double eps = 1e-7;
//  pme.pmetemp.moduli1.resize(250);
//  pme.pmetemp.moduli2.resize(250);
//  pme.pmetemp.moduli3.resize(250);
//  double sinsum, cossum, fact, sumarg, numbersum, zeta, tempsum1, tempsum2;
//  int acutoff, doubleorder, indexi, indexk, indexj;
//  fact = 2.0 * SUPERPI / pme.pmetemp.nxpoints;
//  for (int i = 0; i < pme.pmetemp.nxpoints; i++)
//  {
//    sinsum = cossum = 0.0;
//    for (int k = 0; k < pme.pmetemp.nxpoints; k++)
//    {
//      numbersum = double(i*k);
//      sumarg = fact * numbersum;
//      cossum += stupidarray[k] * cos(sumarg);
//      sinsum += stupidarray[k] * sin(sumarg);
//    }
//    pme.pmetemp.moduli1[i] = sinsum*sinsum + cossum*cossum;
//  }
//  // Coorection for the euler interpolation
//  double tresh = 1e-7;
//  if (pme.pmetemp.moduli1[0] < eps) pme.pmetemp.moduli1[0] = 0.5 * pme.pmetemp.moduli1[1];
//  for (int i = 1; i < pme.pmetemp.nxpoints - 1; i++)
//  {
//    if (pme.pmetemp.moduli1[i] < eps) pme.pmetemp.moduli1[i] = 0.5 * (pme.pmetemp.moduli1[i - 1] + pme.pmetemp.moduli1[i + 1]);
//  }
//  if (pme.pmetemp.moduli1[pme.pmetemp.nxpoints] < eps)  pme.pmetemp.moduli1[pme.pmetemp.nxpoints] = 0.5 * pme.pmetemp.moduli1[pme.pmetemp.nxpoints - 1];
//  // calculate some stupid factor someone called zeta...whatever this thing does....
//  acutoff = 50;
//  doubleorder = 2 * 5;
//  for (int i = 0; i < pme.pmetemp.nxpoints; i++)
//  {
//    indexk = i;
//    if (i > pme.pmetemp.nxpoints / 2) indexk -= pme.pmetemp.nxpoints;
//    if (indexk == 0) zeta = 1.0;
//    else
//    {
//      tempsum1 = tempsum2 = 1.0;
//      fact = SUPERPI * double(indexk) / double(pme.pmetemp.nxpoints);
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact + SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact - SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      zeta = tempsum2 / tempsum1;
//    }
//    pme.pmetemp.moduli1[i] *= (zeta*zeta);
//  }
//  // ###################
//  // SECOND MODULUS ARRAY
//  //####################
//  fact = 2.0 * SUPERPI / pme.pmetemp.nypoints;
//  for (int i = 0; i < pme.pmetemp.nypoints; i++)
//  {
//    sinsum = cossum = 0.0;
//    for (int k = 0; k < pme.pmetemp.nypoints; k++)
//    {
//      numbersum = double(i*k);
//      sumarg = fact * numbersum;
//      cossum += stupidarray[k] * cos(sumarg);
//      sinsum += stupidarray[k] * sin(sumarg);
//    }
//    pme.pmetemp.moduli2[i] = sinsum*sinsum + cossum*cossum;
//  }
//  // Coorection for the euler interpolation
//  if (pme.pmetemp.moduli2[0] < eps) pme.pmetemp.moduli2[0] = 0.5 * pme.pmetemp.moduli2[1];
//  for (int i = 1; i < pme.pmetemp.nypoints - 1; i++)
//  {
//    if (pme.pmetemp.moduli2[i] < eps) pme.pmetemp.moduli2[i] = 0.5 * (pme.pmetemp.moduli2[i - 1] + pme.pmetemp.moduli2[i + 1]);
//  }
//  if (pme.pmetemp.moduli2[pme.pmetemp.nypoints] < eps)  pme.pmetemp.moduli2[pme.pmetemp.nypoints] = 0.5 * pme.pmetemp.moduli2[pme.pmetemp.nypoints - 1];
//  // calculate some stupid factor someone called zeta...whatever this thing does....
//  acutoff = 50;
//  doubleorder = 2 * 5;
//  for (int i = 0; i < pme.pmetemp.nypoints; i++)
//  {
//    indexk = i;
//    if (i > pme.pmetemp.nypoints / 2) indexk -= pme.pmetemp.nypoints;
//    if (indexk == 0) zeta = 1.0;
//    else
//    {
//      tempsum1 = tempsum2 = 1.0;
//      fact = SUPERPI * double(indexk) / double(pme.pmetemp.nypoints);
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact + SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact - SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      zeta = tempsum2 / tempsum1;
//    }
//    pme.pmetemp.moduli2[i] *= (zeta*zeta);
//  }
//  // ###################
//  // THIRD MODULUS ARRAY
//  //####################
//  fact = 2.0 * SUPERPI / pme.pmetemp.nzpoints;
//  for (int i = 0; i < pme.pmetemp.nzpoints; i++)
//  {
//    sinsum = cossum = 0.0;
//    for (int k = 0; k < pme.pmetemp.nzpoints; k++)
//    {
//      numbersum = double(i*k);
//      sumarg = fact * numbersum;
//      cossum += stupidarray[k] * cos(sumarg);
//      sinsum += stupidarray[k] * sin(sumarg);
//    }
//    pme.pmetemp.moduli3[i] = sinsum*sinsum + cossum*cossum;
//  }
//  // Coorection for the euler interpolation
//  if (pme.pmetemp.moduli3[0] < eps) pme.pmetemp.moduli3[0] = 0.5 * pme.pmetemp.moduli3[1];
//  for (int i = 1; i < pme.pmetemp.nzpoints - 1; i++)
//  {
//    if (pme.pmetemp.moduli3[i] < eps) pme.pmetemp.moduli3[i] = 0.5 * (pme.pmetemp.moduli3[i - 1] + pme.pmetemp.moduli3[i + 1]);
//  }
//  if (pme.pmetemp.moduli3[pme.pmetemp.nzpoints] < eps)  pme.pmetemp.moduli3[pme.pmetemp.nzpoints] = 0.5 * pme.pmetemp.moduli3[pme.pmetemp.nzpoints - 1];
//  // calculate some stupid factor someone called zeta...whatever this thing does....
//  acutoff = 50;
//  doubleorder = 2 * 5;
//  for (int i = 0; i < pme.pmetemp.nzpoints; i++)
//  {
//    indexk = i;
//    if (i > pme.pmetemp.nzpoints / 2) indexk -= pme.pmetemp.nzpoints;
//    if (indexk == 0) zeta = 1.0;
//    else
//    {
//      tempsum1 = tempsum2 = 1.0;
//      fact = SUPERPI * double(indexk) / double(pme.pmetemp.nzpoints);
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact + SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact - SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      zeta = tempsum2 / tempsum1;
//    }
//    pme.pmetemp.moduli3[i] *= (zeta*zeta);
//  }
//}
//
//#ifdef _OPENMP
//// calculate spatial sites in the box for partial parallelization
//void coords::Coordinates::roughgrid()
//{
//
//  pme.pmetemp.nrough1 = pme.pmetemp.nrough2 = pme.pmetemp.nrough3 = pme.pmetemp.rgridtotal = 1;
//  int threads;
//#pragma omp parallel
//  {
//    threads = omp_get_num_threads();
//  }
//  for (int i = 2; i <= 6; i++)
//  {
//    if (threads > pme.pmetemp.rgridtotal && (pme.pmetemp.nxpoints%i == 0))
//    {
//      pme.pmetemp.nrough1 = i;
//      pme.pmetemp.rgridtotal = pme.pmetemp.nrough1 * pme.pmetemp.nrough2 * pme.pmetemp.nrough3;
//    }
//    if (threads > pme.pmetemp.rgridtotal && (pme.pmetemp.nypoints%i == 0))
//    {
//      pme.pmetemp.nrough2 = i;
//      pme.pmetemp.rgridtotal = pme.pmetemp.nrough1 * pme.pmetemp.nrough2 * pme.pmetemp.nrough3;
//    }
//    if (threads > pme.pmetemp.rgridtotal && (pme.pmetemp.nzpoints%i == 0))
//    {
//      pme.pmetemp.nrough3 = i;
//      pme.pmetemp.rgridtotal = pme.pmetemp.nrough1 * pme.pmetemp.nrough2 * pme.pmetemp.nrough3;
//    }
//  }
//  // x,y,z-axis number of points per chunk
//  pme.pmetemp.rgrid1 = pme.pmetemp.nxpoints / pme.pmetemp.nrough1;
//  pme.pmetemp.rgrid2 = pme.pmetemp.nypoints / pme.pmetemp.nrough2;
//  pme.pmetemp.rgrid3 = pme.pmetemp.nzpoints / pme.pmetemp.nrough3;
//  // offset for bsplines and left and right points of central point
//  pme.pmetemp.roughleft = (pme.pmetemp.bsplineorder - 1) / 2;
//  pme.pmetemp.roughright = pme.pmetemp.bsplineorder - pme.pmetemp.roughleft - 1;
//  pme.pmetemp.bsoffset = (pme.pmetemp.bsplineorder + 1) / 2 + 1;
//
//}
//#endif

