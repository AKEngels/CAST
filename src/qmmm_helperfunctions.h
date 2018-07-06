#ifndef QMMM_HELPERFUNCTIONS_H
#define QMMM_HELPERFUNCTIONS_H

#include<vector>
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
  double deq_L_QM;
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
  LinkAtom(unsigned int b, unsigned int a, int atomtype, coords::Coordinates *coords, tinker::parameter::parameters const &tp) : qm(b), mm(a), energy_type(atomtype)
  {
    // determine equilibrium distance between link atom and QM atom from force field
    deq_L_QM = 0.0;
    auto b_type_qm = tp.type(coords->atoms().atom(b).energy_type(), tinker::potential_keys::BOND); // bonding energy type for QM atom
    auto b_type = tp.type(energy_type, tinker::potential_keys::BOND);                              // bonding energy type for link atom
    for (auto b_param : tp.bonds())
    {
      if (b_param.index[0] == b_type_qm && b_param.index[1] == b_type)  deq_L_QM = b_param.ideal;
      else if (b_param.index[0] == b_type && b_param.index[1] == b_type_qm) deq_L_QM = b_param.ideal;
    }
    if (deq_L_QM == 0.0)  throw std::runtime_error("Determining position of link atom is not possible.\n");

    // calculate position of link atom
    calc_position(coords);
  }

  /**function to calculate position of link atom (see: doi 10.1002/jcc.20857)
  @param cp: pointer to coordinates object*/
  void calc_position(coords::Coordinates *cp)
  {
    coords::cartesian_type r_MM = cp->xyz(mm);
    coords::cartesian_type r_QM = cp->xyz(qm);
    double d_MM_QM = dist(r_MM, r_QM);

    position = r_QM + ((r_MM - r_QM) / d_MM_QM) * deq_L_QM;
  }
};

/**namespace for diverse functions needed in QM/MM calculations*/
namespace qmmm_helpers
{
  /**function that returns a vector of link atoms
  @param coords: pointer to coordinates object
  @param qm_indices: indizes of QM atoms
  @param mm_indices: indizes of MM atoms
  @param tp: tinker parameter object*/
  std::vector<LinkAtom> create_link_atoms(coords::Coordinates* coords, std::vector<size_t> &qm_indices, std::vector<size_t> &mm_indices, tinker::parameter::parameters const &tp);

  /**calculate gradients on QM and MM atom from link atom (see DOI 10.1002/(SICI)1096-987X(199703)18:4<463::AID-JCC2>3.0.CO;2-R)
  @param l: link atom
  @param G_L: gradient on link atom
  @param coords: coordinates object of whole system
  @param G_QM: reference to QM gradient (is calculated)
  @param G_MM: reference to MM gradient (is calculated)*/
  void calc_link_atom_grad(LinkAtom l, coords::r3 const &G_L, coords::Coordinates* coords, coords::r3 &G_QM, coords::r3 &G_MM);

  /**creates a vector with the indizes of all MM atoms
  @param num_atoms: number of atoms in whole system*/
  std::vector<std::size_t> get_mm_atoms(std::size_t const num_atoms);

  /**creates a vector new_indices_qm
  for description of the vector see energy_int_qmmm.h
  @param num_atoms: number of atoms in whole system*/
  std::vector<std::size_t> make_new_indices_qm(std::size_t const num_atoms);

  /**creates a vector new_indices_mm
  for description of the vector see energy_int_qmmm.h
  @param num_atoms: number of atoms in whole system
  @param mmi: vector with indizes of MM atoms*/
  std::vector<std::size_t> make_new_indices_mm(std::size_t const num_atoms, std::vector<std::size_t> const& mmi);

  /**creates coordobject for MM interface (maybe also be replaced by make_small_coords?)
  @param cp: coordobj for whole system (QM + MM)
  @param indices: indizes of MM atoms
  @param new_indices: new indizes (see new_indices_mm)*/
  coords::Coordinates make_aco_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices);

  /**creates a coordobject out of the whole system with MM interface (used for ONIOM)
  @param cp: pointer to original coordobject*/
  coords::Coordinates make_mmbig_coords(coords::Coordinates const * cp);

  /**creates a coordobject from the QM atoms and link atoms with either QM or MM interface
  @param cp: pointer to original coordobject
  @param indices: indizes of QM atoms
  @param new_indices: vector new_indices_qm
  @param link_atoms: vector with link atom
  @param energy_interface: energy interface*/
  coords::Coordinates make_small_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices, 
    std::vector<LinkAtom> link_atoms, config::interface_types::T energy_interface);

  /**selects only those charges from amber_charges vector which correspond to the indices
  all other charges are removed*/
  void select_from_ambercharges(std::vector<std::size_t> const & indices);
}

#endif
