#include "atomic.h"
#include "configuration.h"
#include "coords_atoms.h"

#include <string>
#include <cstddef>
#include <utility>

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




coords::Atom::Atom(std::string s)
  : m_symbol(s), m_number(atomic::atomic_number_by_symbol(s)),
  m_mass(atomic::massMap[m_number]), m_system(0U), m_etype(0U),
  m_sub_id(ST_DEFAULT), m_fix(false), m_intern_root(false)
{ }


coords::Atom::Atom(std::size_t n)
  : m_symbol(atomic::symbolMap[n]), m_number(n),
  m_mass(atomic::massMap[m_number]), m_system(0U), m_etype(0U),
  m_sub_id(ST_DEFAULT), m_fix(false), m_intern_root(false)
{ }


coords::Atom::Atom(double m)
  : m_symbol(atomic::symbolMap[atomic::atomic_number_by_mass(m)]),
  m_number(atomic::atomic_number_by_mass(m)), m_mass(m), m_system(0U),
  m_etype(0U), m_sub_id(ST_DEFAULT), m_fix(false), m_intern_root(false)
{ }


void coords::Atom::swap(Atom &r)
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




void coords::Atoms::refine()
{
  refine_internals();
  //refine_followups();
  refine_subsystems();
  //refine_torsions();
}


bool coords::Atoms::sub_io_transition(std::size_t a, std::size_t b) const
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
  {
    return atom_is_part_of_ring[i];
  }
};

//#define PRINT_MAIN_AXIS

void coords::Atoms::refine_mains()
{
  fix_rotation(m_atoms);
  std::size_t const N = m_atoms.size();
  Part_of_Ring ringpart(m_atoms);
  // vector saving the atoms which have main torsions attached
  std::vector<std::size_t> atom_has_main_torsion_attached(N + 1U, N + 1U);
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
      m_atoms[i].make_main_idihedral(
        scon::sorted::find(main_torsion_indices, i));
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


void coords::Atoms::refine_internals()
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
  for (std::size_t i(0U); i < N; ++i)
  {
    if (done[i]) continue;
    auto connect_it = Config::get().coords.internal.connect.find(i);
    std::size_t j = 0;
    //std::cout << "Atom " << i << " begins new molecule with internal " << current_internal << '\n';
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
      append_atoms(lvl + 1, index, molecule, ++index_size, done);
    }
  }
}


void coords::Atoms::get_relatives(std::size_t const i, const std::size_t b)
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
          if (atom(b).iangle() >(S + 1))
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


void coords::Atoms::refine_subsystems()
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
  for (std::size_t i = 0; i < N; ++i)
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
  //Gradients_3D const & gxyz(p.gradient.cartesian);
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
    //std::cout << "Main Torsion " << j << " which is internal " << mti << " is " << intern[mti].azimuth() << '\n';
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
    //std::cout << "Main Torsion " << j << " which is internal " << mti << " is " << intern[mti].azimuth() << '\n';
    p.structure.main[j] = intern[mti].azimuth();
    if (atom(mti).ibond() < N)
    {
      for (auto const rotating : atom(atom(mti).ibond()).bound_internals())
      {
        //std::cout << "Rotates: " << rotating << " which is atom " << atom(rotating).i_to_a() << '\n';
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
  for (size_1d::size_type i = 0; i < TA; ++i)
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
  for (size_1d::size_type i(0U); i < O; ++i)
  {
    // atomicA and B are sorted, so for the same i the should have the same numbers if equal
    if (atomicA[i] != atomicB[i]) return false;
    // we search for a residue j at B that is equal to residue i at A
    // j has to be marked as "not found yet" (done_b[j] = false)
    bool found = false;
    for (std::vector<bool>::size_type j = 0; j < O; ++j)
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

void coords::Atoms::swap(Atoms &r)
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

size_t coords::Atoms::getNumberOfAtomsWithAtomicNumber(size_t searchedNumber) const
{
  // Find number of atoms
  size_t counter = 0u;
  for (size_t i = 0u; i < this->size(); i++)
  {
    if (this->m_atoms[i].number() == searchedNumber) counter++;
  }
  return counter;
}