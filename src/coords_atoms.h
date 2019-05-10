#ifndef coords_atoms_h_3336f4e7_81a8_4d3f_994f_e4104ee90926
#define coords_atoms_h_3336f4e7_81a8_4d3f_994f_e4104ee90926

#pragma once

#include <cstddef>
#include <string>
#include <cstddef>
#include <utility>
#include <deque>
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


	/**class for relations between cartesian and internal coordinates
	class Atom is a child class of the class so the functions are available for Atom also*/
  class internal_relations
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

		/**default constructor*/
    internal_relations()
      : rel_b(0U), rel_a(0U), rel_d(0U), atom_of_inter_index(rel_b),
      inter_of_atom_index(rel_b), main_torsion_index(0U), is_main_torsion(false), m_ifix(false)
    { }
		/**constructor*/
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
    std::size_t dependant_idihedrals(std::size_t index) const { return dependant_torsions[index]; }
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
    /**get internal atom index */
    std::size_t i_to_a() const { return atom_of_inter_index; }
		/**set internal atom index*/
    void set_a_of_i(std::size_t index) { atom_of_inter_index = index; }
		/**get cartesian atom index */
    std::size_t a_to_i() const { return inter_of_atom_index; }
		/**set cartesian atom index*/
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


	/**class for an atom*/
  class Atom
    : public internal_relations
  {

  public:

		/**enum for subsystem types (used for FEP)*/
    enum sub_types { ST_DEFAULT, ST_IN, ST_OUT };

  private:

    // types
    typedef std::vector<std::size_t> bondvector_t;

    /**element symbol*/
    std::string m_symbol;
		/**element number in PES*/
    std::size_t m_number;
		/**atomic mass*/
		double m_mass;
		/**covalent radius*/
		double m_cov_rad;
		/**bonding partners (cartesian indices)*/
		size_1d m_bonds;
		/**bonding partners (internal indices)*/
		size_1d m_ibound;
		/**subsystem type*/
		std::size_t m_system;
		/**forcefield energy type*/
		std::size_t m_etype;
		/**subsystem id*/
    sub_types m_sub_id;
		/**is this atom fixed?*/
		bool m_fix;
		/**is this atom root atom of internal coordinates 
		(i.e. its internal coordinates are not defined by other atoms
		but by absolute position in space, relative to virtual atoms)*/
		bool m_intern_root;
    /**name of the residue (from pdb)*/
    std::string residue;
    /**unique residue id (from pdb)*/
    int(res_id);
    /**atom name from pdb*/
    std::string pdb_atom_name;

  public:

		/**default constructor*/
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
    // return covalent radius
    double cov_radius() const {return m_cov_rad;}
    // return bonds
    size_1d const & bonds() const { return m_bonds; }
    std::size_t bonds(std::size_t const i) const { return m_bonds[i]; }
    // add bond to i
    void bind_to(std::size_t const i) { scon::sorted::insert_unique(m_bonds, i); }
    // remove bond to i
    void detach_from(std::size_t const i) { scon::sorted::remove_equal(m_bonds, i); }
    // check if bound to i
    bool is_bound_to(std::size_t const i) const { return scon::sorted::exists(m_bonds, i); }
    // get internal bonds
    size_1d const & bound_internals() const { return m_ibound; }
    std::size_t bound_internals(std::size_t const i) const { return m_ibound[i]; }
		/**delete all bound internals*/
		void clear_bound_internals() { m_ibound.clear(); }
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
    // change if atom is root
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
    /**set residue name*/
    void set_residue(std::string s) { residue = s; }
    /**get residue name*/
    std::string get_residue() const { return residue; }
    /**set res_id*/
    void set_res_id(int i) { res_id = i; }
    /**get red_id*/
    int get_res_id() const { return res_id; }
    /**set pdb atom name*/
    void set_pdb_atom_name(std::string i) { pdb_atom_name = i; }
    /**get pdb atom name*/
    std::string get_pdb_atom_name() const { return pdb_atom_name; }

    // Note (DK): This function does not seem to be used currently
    void set_relation(internal_relations const &r) { internal_relations::operator=(r); }

    // swap contents
    void swap(Atom &r);

    friend std::ostream& operator<< (std::ostream &stream, Atom const & atom);
    friend bool operator== (Atom const &lhs, Atom const &rhs)
    {
      if (lhs.m_mass == rhs.m_mass)
      {
        return true;
      }
      else
      {
        return false;
      }
    }
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
    // internal indices of the internal torsions that are "main"
    size_1d main_bond_indices, main_angle_indices, main_torsion_indices;
    // subsystem index "incoming" and "outgoing"
    std::size_t m_sub_in_index, m_sub_out_index;
    // do we have in AND out?
    bool m_sub_io, m_in_exists, m_out_exists;

    // subsystem helper stuff
    void refine_subsystems();

    // internal helper stuff

		/**get atoms which are used for creating internal coordinates (relative atoms) for all atoms*/
    void refine_internals();
		/** function to find relatives (atoms that define internal coordinates) of an atom
		@param i: cartesian atom index of which we want to find the atoms that define the internal coordinates
		@param b: cartesian atom index of an atom which is bound to atom i*/
    void get_relatives(std::size_t const i, std::size_t const b);
    void append_atoms(std::size_t const lvl, std::size_t const A, size_1d &molecule, 
      std::size_t &index_size, std::deque<bool> &done);
    //void refine_followups();

    // New stuff
    void refine_mains();

    bool common_torsion_axis(std::size_t a, std::size_t b, bool & direction) const;
    std::size_t atom_by_intern(std::size_t index) const { return m_atoms[index].i_to_a(); }
    bool common_bond(std::size_t a, std::size_t b) const;

		/**gets cartesian coordinates of the relative atom 'index'
		in most cases these are just the 'normal' Cartesian coordinates of this atom
		but if index is bigger than number of atoms in molecule it's a virtual atom
		then there are previously defined coordinates created in order to be able to define whole molecule
		@param index: cartesian atom index for which the position should be returned
		@param xyz: cartesian coordinates of all atoms*/
    Cartesian_Point rel_xyz(std::size_t const index, Representation_3D const & xyz) const;

  public:

		/**default constructor*/
    Atoms()
      : m_sub_in_index(0U), m_sub_out_index(0U),
      m_sub_io(false), m_in_exists(false), m_out_exists(false)
    { }

    /* Returns number of atoms in coords object */
    std::vector<Atom>::size_type size() const 
    { 
      return m_atoms.size(); 
    }
		/**swap in*/
    void swap(Atoms&);

    std::vector<Atom>::iterator begin() { return m_atoms.begin(); }
    std::vector<Atom>::iterator end() { return m_atoms.end(); }
    std::vector<Atom>::const_iterator begin() const { return m_atoms.begin(); }
    std::vector<Atom>::const_iterator end() const { return m_atoms.end(); }

    void clear() { Atoms().swap(*this); }

    // general stuff

		/**add an atom*/
    void add(Atom const &atom) { m_atoms.push_back(atom); }
		/**delete last atom*/
    void pop_back() { m_atoms.pop_back(); }
		/**return atom with index 'index' to change it*/
    Atom & atom(std::size_t index) { return m_atoms[index]; }
		/**return atom with index 'index'*/
    Atom const & atom(std::size_t index) const { return m_atoms[index]; }
		/**return all molecules*/
    size_2d const & molecules() const { return m_molecules; }
		/**return one molecule
		@param index: index of molecule*/
    size_1d const & molecule(std::size_t index) const { return m_molecules[index]; }
		/**returns a given atom of a given molecule
		@param molecule: index of molecule
		@param atom: index of atom in molecule*/
    std::size_t atomOfMolecule(std::size_t molecule, std::size_t atom) const { return m_molecules[molecule][atom]; }
		/**refine, e.g. built up of z-matrix and a lot of other stuff*/
    void refine();
    /**get all subsystems*/
    size_2d const & subsystems() const { return m_sub_systems; }
		/**get subsystem with index 'index'*/
    size_1d const & subsystems(std::size_t index) const { return m_sub_systems[index]; }
		/**does subsystem IN exist?*/
    bool in_exists() const { return m_in_exists; }
		/**does subsystem OUT exist?*/
    bool out_exists() const { return m_out_exists; }
		/**index of IN*/
    std::size_t sub_in() const { return m_sub_in_index; }
		/**index of OUT*/
    std::size_t sub_out() const { return m_sub_out_index; }
    // Matrix mit indices geh�rig zu PESpoint::sub_ia_matrix_t interaction matrix
    std::size_t sub_ia_index() const { return scon::triangularIndex<true>(m_sub_in_index, m_sub_out_index); }
		/**???*/
    bool sub_io() const { return m_sub_io; }
    // Sagt mir ob a out und b in ist oder vice versa
    bool sub_io_transition(std::size_t a, std::size_t b) const;
    //bool valid_sub_ia (std::size_t a, std::size_t b) const { return }

    // internal stuff, index is main, what is the internal index?
    std::size_t intern_of_main_idihedral(std::size_t index) const { return main_torsion_indices[index]; }
		/**returns all mains*/
    size_1d const & mains() const { return main_torsion_indices; }

    // stereo stuff
    bool res_is_equal(std::size_t a, std::size_t b, std::size_t from_a, std::size_t from_b, std::size_t deepth) const;
    
    //check fixation of an atom
    bool check_fix(std::size_t atom){ return m_atoms[atom].fixed(); }

    /**convert internal to cartesian coordinates*/
    void i_to_c(PES_Point&) const;
		/**convert cartesian to internal coordinates*/
    void c_to_i(PES_Point&) const;
		/**convert cartesian to internal coordinates, but without gradients*/
    void c_to_i_light(PES_Point&) const;

    // Helper
    size_t getNumberOfAtomsWithAtomicNumber(size_t searchedNumber) const;
		/**fix or unfix all atoms
		@param fix_it: true if atom should be fixed, false if it should be unfixed*/
    void fix_all(bool const fix_it = true);
		/**fix or unfix one atom
		@param atom: index of atom that should be fixed
		@param fix_it: true if atom should be fixed, false if it should be unfixed*/
    void fix(std::size_t const atom, bool const fix_it = true);
    // @todo document return vlaue
    std::size_t intern_of_dihedral(std::size_t a, std::size_t b, std::size_t c, std::size_t d) const;


    friend std::ostream& operator<< (std::ostream &stream, Atoms const & atoms);
    friend bool operator== (Atoms const &lhs, Atoms const &rhs)
    {
      bool ret = true;
      for (size_t i = 0; i < lhs.size(); ++i)
      {
        ret = ret && lhs.atom(i) == rhs.atom(i);
      }
      return ret;
    }

    friend bool operator!= (Atoms const &lhs, Atoms const &rhs)
    {
      return !operator== (lhs, rhs);
    }

  };

  inline void swap(Atoms &a, Atoms &b)
  {
    a.swap(b);
  }

}

#endif // coords_atoms_h_3336f4e7_81a8_4d3f_994f_e4104ee90926