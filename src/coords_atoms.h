#ifndef coords_atoms_h_3336f4e7_81a8_4d3f_994f_e4104ee90926
#define coords_atoms_h_3336f4e7_81a8_4d3f_994f_e4104ee90926

#pragma once

#include <cstddef>

#include "coords_rep.h"
#include "scon_matrix.h"

namespace coords
{

  /* ##############################################################################


  #### ##    ## ######## ######## ########  ##    ##    ###    ##
   ##  ###   ##    ##    ##       ##     ## ###   ##   ## ##   ##
   ##  ####  ##    ##    ##       ##     ## ####  ##  ##   ##  ##
   ##  ## ## ##    ##    ######   ########  ## ## ## ##     ## ##
   ##  ##  ####    ##    ##       ##   ##   ##  #### ######### ##
   ##  ##   ###    ##    ##       ##    ##  ##   ### ##     ## ##
  #### ##    ##    ##    ######## ##     ## ##    ## ##     ## ########

  ########  ######## ##          ###    ######## ####  #######  ##    ##  ######
  ##     ## ##       ##         ## ##      ##     ##  ##     ## ###   ## ##    ##
  ##     ## ##       ##        ##   ##     ##     ##  ##     ## ####  ## ##
  ########  ######   ##       ##     ##    ##     ##  ##     ## ## ## ##  ######
  ##   ##   ##       ##       #########    ##     ##  ##     ## ##  ####       ##
  ##    ##  ##       ##       ##     ##    ##     ##  ##     ## ##   ### ##    ##
  ##     ## ######## ######## ##     ##    ##    ####  #######  ##    ##  ######


  ############################################################################## */


  class internal_relations // internal relation
  {
    /*
    rel_b : internal index in bonding relation to this
    rel_a : internal index in 1-3(angle) relation to this
    rel_d : internal index in 1-4(dihedral) relation to this
    atom_of_inter_index : the atom index of the current intern index
    i.e. relations[i].atom_of_inter_index is the atom that corresponds to
    internal index i
    inter_of_atom_index : the internal index of the current atom index
    vice versa
    main_torsion_index : valid if this is a main_torsion
    follow_ups : internal indices that follow the current one (moved if this moves)
    dependant_torsions : torsions that share the same rel_b and rel_a
    is_main_torsion : is this a main torsion?
    */
  protected:

    std::size_t rel_b, rel_a, rel_d, atom_of_inter_index, inter_of_atom_index, main_torsion_index;
    size_1d follow_ups, dependant_torsions, m_successors;
    bool is_main_torsion, m_ifix;

  public:

    internal_relations()
      : rel_b(0U), rel_a(0U), rel_d(0U), atom_of_inter_index(rel_b),
      inter_of_atom_index(rel_b), main_torsion_index(0U), is_main_torsion(false), m_ifix(false)
    { }
    internal_relations(std::size_t atom, std::size_t intern, std::size_t rb, std::size_t ra, std::size_t rd)
      : rel_b(rb), rel_a(ra), rel_d(rd), atom_of_inter_index(atom),
      inter_of_atom_index(intern), main_torsion_index(0U), is_main_torsion(false), m_ifix(false)
    { }

    // get/set "bond-partner" of index
    std::size_t ibond() const { return rel_b; }
    void set_ibond(std::size_t index) { rel_b = index; }
    // get/set "angle-partner" of index
    std::size_t iangle() const { return rel_a; }
    void set_iangle(std::size_t index) { rel_a = index; }
    // get/set "dihedral-partner" of index
    std::size_t idihedral() const { return rel_d; }
    void set_idihedral(std::size_t index) { rel_d = index; }
    // append a new follower
    void append_follower(std::size_t index) { scon::sorted::insert_unique(follow_ups, index); }
    // get followers
    size_1d const & followers() const { return follow_ups; }
    // get single follower
    std::size_t follower(std::size_t index) const { return follow_ups[index]; }
    // get successors
    size_1d const & successors() const { return m_successors; }
    // get successor
    std::size_t successors(std::size_t const index) const { return m_successors[index]; }
    // append a new dependant torsion
    void append_dependant_idihedral(std::size_t index) { scon::sorted::insert_unique(dependant_torsions, index); }
    // get deptorsions
    size_1d const & dependant_idihedrals() const { return dependant_torsions; }
    std::size_t const & dependant_idihedrals(std::size_t index) const { return dependant_torsions[index]; }
    // is this torsional relation a main torsion?
    bool is_main_idihedral() const { return is_main_torsion; }
    // make this torsional relation a main torsion
    void make_main_idihedral(std::size_t index)
    {
      main_torsion_index = index;
      is_main_torsion = true;
    }
    // main_torsion index of this
    std::size_t main_idihedral_index() const { return main_torsion_index; }
    // atom index of this internal
    std::size_t i_to_a() const { return atom_of_inter_index; }
    void set_a_of_i(std::size_t index) { atom_of_inter_index = index; }
    // internal index of this atom
    std::size_t a_to_i() const { return inter_of_atom_index; }
    void set_i_of_a(std::size_t index) { inter_of_atom_index = index; }

    bool ifix() const { return m_ifix; }
    void ifix(bool const val) { m_ifix = val; }

    friend std::ostream& operator<< (std::ostream &stream, internal_relations const & inter);

  };


  /* ############################################


     ###    ########  #######  ##     ##
    ## ##      ##    ##     ## ###   ###
   ##   ##     ##    ##     ## #### ####
  ##     ##    ##    ##     ## ## ### ##
  #########    ##    ##     ## ##     ##
  ##     ##    ##    ##     ## ##     ##
  ##     ##    ##     #######  ##     ##


  ############################################ */


  class Atom
    : public internal_relations
  {

  public:

    enum sub_types { ST_DEFAULT, ST_IN, ST_OUT };

  private:

    // types
    typedef std::vector<std::size_t> bondvector_t;

    // data
    std::string m_symbol;
    std::size_t m_number;
    double m_mass;
    size_1d m_bonds, m_ibound;
    std::size_t m_system, m_etype;
    sub_types m_sub_id;
    bool m_fix, m_intern_root;

  public:

    Atom()
      : m_number(0U), m_mass(0.0), m_system(0U), m_etype(0U),
      m_sub_id(sub_types::ST_DEFAULT), m_fix(false), m_intern_root(false)
    { }

    // constructor from symbol
    explicit Atom(std::string);
    // constructor from number
    explicit Atom(std::size_t);
    // constructor from mass
    explicit Atom(double);

    // return mass
    double mass() const { return m_mass; }
    // return atomic number
    std::size_t number() const { return m_number; }
    // return symbol
    std::string const & symbol() const { return m_symbol; }
    // return bonds
    size_1d const & bonds() const { return m_bonds; }
    std::size_t const & bonds(std::size_t const i) const { return m_bonds[i]; }
    // add bond to i
    void bind_to(std::size_t const i) { scon::sorted::insert_unique(m_bonds, i); }
    // remove bond to i
    void detach_from(std::size_t const i) { scon::sorted::remove_equal(m_bonds, i); }
    // check if bound to i
    bool is_bound_to(std::size_t const i) const { return scon::sorted::exists(m_bonds, i); }
    // get internal bonds
    size_1d const & bound_internals() const { return m_ibound; }
    std::size_t const & bound_internals(std::size_t const i) const { return m_ibound[i]; }
    // add internally bound atom
    void attach(std::size_t const i) { scon::sorted::insert_unique(m_ibound, i); }
    // remove internally bound atom
    void detach(std::size_t const i) { scon::sorted::remove_equal(m_ibound, i); }
    // is this atom fixed?
    bool fixed() const { return m_fix; }
    // change fixation
    void fix(bool const fix_it = true) { m_fix = fix_it; }
    // is this atom internal root?
    bool is_root() const { return m_intern_root; }
    // change fixation
    void root(bool const root_it = true) { m_intern_root = root_it; }
    // system number
    std::size_t system() const { return m_system; }
    // assign new system
    void assign_to_system(std::size_t const system_index) { m_system = system_index; }
    // subsystem type
    sub_types sub_type() const { return m_sub_id; }
    // set type
    void set_sub_type(sub_types const id) { m_sub_id = id; }
    // subsystem type
    std::size_t energy_type() const { return m_etype; }
    // set type
    void set_energy_type(std::size_t const id) { m_etype = id; }

    void set_relation(internal_relations const &r) { internal_relations::operator=(r); }
    // swap contents
    void swap(Atom &r);

    friend std::ostream& operator<< (std::ostream &stream, Atom const & atom);

  };

  inline void swap(Atom &a, Atom &b)
  {
    a.swap(b);
  }


  /* ##############################################


     ###    ########  #######  ##     ##  ######
    ## ##      ##    ##     ## ###   ### ##    ##
   ##   ##     ##    ##     ## #### #### ##
  ##     ##    ##    ##     ## ## ### ##  ######
  #########    ##    ##     ## ##     ##       ##
  ##     ##    ##    ##     ## ##     ## ##    ##
  ##     ##    ##     #######  ##     ##  ######


  ############################################## */


  class Atoms
  {

  private:

    // the atoms themselves
    std::vector<Atom> m_atoms;
    // the systems the atoms are grouped into
    size_2d m_sub_systems, m_molecules;
    // indices of the internal bonds/angles/torsions that are "main"
    size_1d main_bond_indices, main_angle_indices, main_torsion_indices;
    // subsystem index "incoming" and "outgoing"
    std::size_t m_sub_in_index, m_sub_out_index;
    // do we have in AND out?
    bool m_sub_io, m_in_exists, m_out_exists;

    // subsystem helper stuff
    void refine_subsystems();
    // internal helper stuff
    void refine_internals();
    void get_relatives(std::size_t const i, std::size_t const b);
    void append_atoms(std::size_t const lvl, std::size_t const A, size_1d &molecule, 
      std::size_t &index_size, std::vector<bool> &done);
    //void refine_followups();
    void refine_torsions();
    bool atom_fixed;

    // New stuff
    void refine_mains();

    bool common_torsion_axis(std::size_t a, std::size_t b, bool & direction) const;
    std::size_t atom_by_intern(std::size_t index) const { return m_atoms[index].i_to_a(); }
    bool common_bond(std::size_t a, std::size_t b) const;
    void add_main_idihedral(std::size_t index);

    Cartesian_Point rel_xyz(std::size_t const index, Representation_3D const & xyz) const;

  public:

    Atoms()
      : m_sub_in_index(0U), m_sub_out_index(0U),
      m_sub_io(false), m_in_exists(false), m_out_exists(false)
    { }
    std::vector<Atom>::size_type size() const { return m_atoms.size(); }
    void swap(Atoms&);

    std::vector<Atom>::iterator begin() { return m_atoms.begin(); }
    std::vector<Atom>::iterator end() { return m_atoms.end(); }
    std::vector<Atom>::const_iterator begin() const { return m_atoms.begin(); }
    std::vector<Atom>::const_iterator end() const { return m_atoms.end(); }

    void clear() { Atoms().swap(*this); }

    // general stuff
    void add(Atom const &atom) { m_atoms.push_back(atom); }
    void pop_back() { m_atoms.pop_back(); }
    Atom & atom(std::size_t index) { return m_atoms[index]; }
    Atom const & atom(std::size_t index) const { return m_atoms[index]; }
    size_2d const & molecules() const { return m_molecules; }
    size_1d const & molecules(std::size_t index) const { return m_molecules[index]; }
    std::size_t const & molecules(std::size_t molecule, std::size_t atom) const { return m_molecules[molecule][atom]; }
    void refine();
    // subsystem stuff
    size_2d const & subsystems() const { return m_sub_systems; }
    size_1d const & subsystems(std::size_t index) const { return m_sub_systems[index]; }
    bool in_exists() const { return m_in_exists; }
    bool out_exists() const { return m_out_exists; }
    std::size_t sub_in() const { return m_sub_in_index; }
    std::size_t sub_out() const { return m_sub_out_index; }
    std::size_t sub_ia_index() const { return scon::triangularIndex<true>(m_sub_in_index, m_sub_out_index); }
    bool sub_io() const { return m_sub_io; }
    bool sub_io_transition(std::size_t a, std::size_t b) const;
    //bool valid_sub_ia (std::size_t a, std::size_t b) const { return }
    // internal stuff
    std::size_t intern_of_main_idihedral(std::size_t index) const { return main_torsion_indices[index]; }
    size_1d const & mains() const { return main_torsion_indices; }
    // stereo stuff
    bool res_is_equal(std::size_t a, std::size_t b, std::size_t from_a, std::size_t from_b, std::size_t deepth) const;
    // internal to cartesian et vice versa
    //void internal_to_cartesian(PES_Point&) const;
    //void cartesian_to_internal(PES_Point&) const;

    // New Translation
    void i_to_c(PES_Point&) const;
    void c_to_i(PES_Point&) const;
    void c_to_i_light(PES_Point&) const;

    // Helper

    //void update_int_from_xyz_movement(std::size_t index, PES_Point&);
    void fix_all(bool const fix_it = true);
    void fix(std::size_t const atom, bool const fix_it = true);
    std::size_t intern_of_dihedral(std::size_t a, std::size_t b, std::size_t c, std::size_t d) const;
    friend std::ostream& operator<< (std::ostream &stream, Atoms const & atoms);
  };

  inline void swap(Atoms &a, Atoms &b)
  {
    a.swap(b);
  }

}

#endif // coords_atoms_h_3336f4e7_81a8_4d3f_994f_e4104ee90926