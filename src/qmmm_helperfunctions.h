/**
CAST 3
qmmm_helperfunctions.h
Purpose: some stuff that is used by all QM/MM interfaces, i. e. additive and subtractive QM/MM as well as three-layer

@author Susanne Sauer
@version 1.0
*/

#ifndef QMMM_HELPERFUNCTIONS_H
#define QMMM_HELPERFUNCTIONS_H

#include<vector>
#include"atomic.h"
#include"configuration.h"
#include"coords.h"
#include"coords_io.h"
#include"helperfunctions.h"
#include"tinker_parameters.h"

/**struct with all relevant information about a link atom*/
struct LinkAtom
{
  /**position*/
  coords::Cartesian_Point position;
  /**equilibrium distance to QM atom*/
  double deq_L_QM{ 0.0 };
  /**index of QM atom*/
  unsigned int qm;
  /**index of MM atom*/
  unsigned int mm;
  /**force field atom type of link atom*/
  size_t energy_type;

  /**default constructor*/
  LinkAtom() { }

  /**another constructor
  @param b: index of QM atom
  @param a: index of MM atom
  @param atomtype: energytype of link atom in forcefield
  @param coords: pointer to coordinates object
  @param tp: tinker parameter object*/
  LinkAtom(unsigned int b, unsigned int a, int atomtype, coords::Coordinates* coords, tinker::parameter::parameters const& tp);

  /**function to calculate position of link atom (see: doi 10.1002/jcc.20857)
  @param cp: pointer to coordinates object*/
  void calc_position(coords::Coordinates* cp);
};

namespace energy
{
  namespace interfaces
  {
    /**namespace for QMMM interfaces*/
    namespace qmmm
    {
      /**function that returns a vector of link atoms
      @param qm_indices: indizes of QM atoms
      @param coords: pointer to coordinates object
      @param tp: tinker parameter object
      @param linkatomtypes: forcefield types of link atoms (can be empty if they don't matter)*/
      std::vector<LinkAtom> create_link_atoms(std::vector<size_t> const& qm_indices, coords::Coordinates* coords,
        tinker::parameter::parameters const& tp, std::vector<int> const& linkatomtypes);

      /**function that creates several link atom sets (one for each QM system)
      @param qm_indices: indizes of QM atoms (one vector for each QM system)
      @param coords: pointer to coordinates object
      @param tp: tinker parameter object
      @param linkatomtypes: forcefield types of link atoms (one vector for each QM system)*/
      std::vector<std::vector<LinkAtom>> create_several_linkatomsets(std::vector < std::vector<size_t>> const& qm_indices, coords::Coordinates* coords,
        tinker::parameter::parameters const& tp, std::vector < std::vector<int>> const& linkatomtypes);

      /**calculate gradients on QM and MM atom from link atom (see DOI 10.1002/(SICI)1096-987X(199703)18:4<463::AID-JCC2>3.0.CO;2-R)
      @param l: link atom
      @param G_L: gradient on link atom
      @param coords: coordinates object of whole system
      @param G_QM: reference to QM gradient (is calculated)
      @param G_MM: reference to MM gradient (is calculated)*/
      void calc_link_atom_grad(LinkAtom const& l, coords::r3 const& G_L, coords::Coordinates* coords, coords::r3& G_QM, coords::r3& G_MM);

      /**creates a vector with the indizes of all MM atoms
      @param num_atoms: number of atoms in whole system*/
      std::vector<std::size_t> get_mm_atoms(std::size_t const num_atoms);

      /**creates a vector new_indices
      for description of the vector see energy_int_qmmm.h
      @param indices: vector with original indizes
      @param num_atoms: number of atoms in whole system*/
      std::vector<std::size_t> make_new_indices(std::vector<std::size_t> const& indices, coords::Coordinates::size_type const num_atoms);

      /**creates several vectors of new indices new_indices, one for each QM system
      for description of the vector see energy_int_qmmm.h
      @param indices: vectors with original indizes (one for each QM system)
      @param num_atoms: number of atoms in whole system*/
      std::vector<std::vector<std::size_t>> make_several_new_indices(std::vector<std::vector<std::size_t>> const& indices, coords::Coordinates::size_type const num_atoms);

      /**creates a coordobject for partial systems (if QM system link atoms are added, too)
      @param cp: pointer to original coordobject
      @param indices: indizes of atoms that should be in new coordobject
      @param new_indices: vector of length total number of atoms
                          only those elements are filled whose position corresponds to atoms of new coordobject
                          they are filled with successive numbers starting from 0
                          purpose: faciliate mapping between total coordinates object and subsystems
      @param energy_interface: energy interface of new coordobject
      @param system_information: which system is currenty created? (string to be printed if verbosity > 3)
      @param write_into_file: if true writes subsystem into tinkerfile
      @param link_atoms: vector with link atoms (if present)
      @param filename: name of the tinkerfile (only when write_into_file = true)*/
      coords::Coordinates make_small_coords(coords::Coordinates const* cp,
        std::vector<std::size_t> const& indices, std::vector<std::size_t> const& new_indices, config::interface_types::T energy_interface, std::string const& system_information,
        bool const write_into_file = false, std::vector<LinkAtom> const& link_atoms = std::vector<LinkAtom>(), std::string const& filename = "qm_system.arc");

      /**creates several coordobjects for several partial systems
      @param cp: pointer to original coordobject
      @param indices: indizes of atoms that should be in new coordobject (one vector for each coordobject that should be created)
      @param new_indices: vector of length total number of atoms (one vector for each coordobject that should be created)
                          only those elements are filled whose position corresponds to atoms of new coordobject
                          they are filled with successive numbers starting from 0
                          purpose: faciliate mapping between total coordinates object and subsystems
      @param energy_interface: energy interface of new coordobject
      @param write_into_file: if true writes subsystem into tinkerfile
      @param link_atoms: vector with link atoms (if present), one vector for each coordobject that should be created*/
      std::vector<coords::Coordinates> make_several_small_coords(coords::Coordinates const* cp,
        std::vector<std::vector<std::size_t>> const& indices, std::vector < std::vector<std::size_t>> const& new_indices, config::interface_types::T energy_interface, bool const write_into_file = false,
        std::vector < std::vector<LinkAtom>> const& link_atoms = std::vector < std::vector<LinkAtom>>());

      /**selects only those charges from atom_charges vector of Coordinates object cp which correspond to the indices
      all other charges are removed*/
      std::vector<coords::float_type> select_from_atomcharges(std::vector<std::size_t> const& indices, coords::Coordinates const* cp);

      /**This function modifies the coordinates of the current charge in case of periodic boundaries:
      The distance between the center of the QM system and the position of the charge is determined.
      If any of the components (x, y or z) of the connecting vector is longer than half the box size,
      the charge is moved in such a way that it is closer by the QM region by adding or subtracting the box size.
      In this way the distance between the charge and the QM region is minimized.
      @param current_coords: reference to the position of the current charge, will by modified
      @param center_of_QM: geometrical center of QM region*/
      void move_periodics(coords::Cartesian_Point& current_coords, coords::Cartesian_Point const& center_of_QM);

      /**function to determine the atom index of the center of the QM region
      @param default_index: if this atom is part of QM region it will be directly returned
      @param qm_indizes: indices of the atoms that define the QM region
      @param coords: pointer to original coordobject*/
      std::size_t get_index_of_QM_center(std::size_t const default_index, std::vector<size_t> const& qm_indizes, coords::Coordinates* coords);

      std::vector<std::size_t> get_indices_of_several_QMcenters(std::vector<std::size_t> const default_indices, std::vector<std::vector<std::size_t>> const& qm_indices,
        coords::Coordinates* coords);

      /**adds external charges to the following calculations
      @param ignore_indizes: indizes of atoms that should be ignored
      @param charges: vector of charge values that might be added to the calculation
      @param indizes_of_charges: indizes of the charges in the overall coordinates object
      @param link_atoms: vector of link atoms for the current "QM system"
      @param charge_indizes: reference to a vector where the indizes of the atoms whose charges are taken into account are added
      @param coords: pointer to original coordobject
      @param QMcenter: index of atom that defines center of QM region*/
      void add_external_charges(std::vector<size_t> const& ignore_indizes, std::vector<double> const& charges, std::vector<size_t> const& indizes_of_charges,
        std::vector<LinkAtom> const& link_atoms, std::vector<int>& charge_indizes, coords::Coordinates* coords, std::size_t const QMcenter);

      /**renames outputfiles for calculations with external energyinterfaces to prevent them from being overwritten
      @param interface: energy interface for which files should be renamed (can be DFTB, MOPAC, ORCA, GAUSSIAN or PSI4)
      @param id: id from which filesnames in that interface are created (should be member of energy interface)
      @param systemname: string which is inserted in filenames*/
      void save_outputfiles(config::interface_types::T const& interface, std::string const& id, std::string const& systemname);
    }
  }
}

#endif
