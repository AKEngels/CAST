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
  /**force field atom type of link atom*/
  size_t energy_type;
  /**equilibrium distance to QM atom*/
  double deq_L_QM;
  /**index of QM atom*/
  int qm;
  /**index of MM atom*/
  int mm;

  /**default constructor*/
  LinkAtom() { }

  /**another constructor
  @param b: index of QM atom
  @param a: index of MM atom
  @param coords: pointer to coordinates object
  @param tp: tinker parameter object*/
  LinkAtom(int b, int a, coords::Coordinates *coords, tinker::parameter::parameters const &tp) : qm(b), mm(a)
  {
    // determine energy type for link atom
    if (Config::get().energy.qmmm.mminterface == config::interface_types::T::OPLSAA)
    {
      energy_type = 85;
    }
    else if (Config::get().energy.qmmm.mminterface == config::interface_types::T::AMBER)
    {
      energy_type = 3024;
    }
    else throw("Something went wrong. Invalid MM interface for QM/MM.");

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
    position = calc_position(coords);
  }

  /**function to calculate position of link atom (see: doi 10.1002/jcc.20857)
  @param cp: pointer to coordinates object*/
  coords::cartesian_type calc_position(coords::Coordinates *cp)
  {
    coords::cartesian_type r_MM = cp->xyz(mm);
    coords::cartesian_type r_QM = cp->xyz(qm);
    double d_MM_QM = dist(r_MM, r_QM);

    coords::cartesian_type pos = r_QM + ((r_MM - r_QM) / d_MM_QM) * deq_L_QM;
    return pos;
  }
};

/**namespace for diverse functions needed in QM/MM calculations*/
namespace qmmm_helpers
{
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

  /**creates coordobject for QM interface
  @param cp: coordobj for whole system (QM + MM)
  @param indices: indizes of QM atoms
  @param new_indices: new indizes (see new_indices_qm)*/
  coords::Coordinates make_qm_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices);

  /**creates coordobject for MM interface
  @param cp: coordobj for whole system (QM + MM)
  @param indices: indizes of MM atoms
  @param new_indices: new indizes (see new_indices_mm)*/
  coords::Coordinates make_aco_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices);

  /**creates a coordobject out of the whole system with MM interface (used for ONIOM)
  @param cp: pointer to original coordobject*/
  coords::Coordinates make_mmbig_coords(coords::Coordinates const * cp);

  /**creates a coordobject from the QM atoms and link atoms but with MM interface
  @param cp: pointer to original coordobject
  @param indices: indizes of QM atoms
  @param new_indices: vector new_indices_qm
  @param link_atoms: vector with link atoms*/
  coords::Coordinates make_mmsmall_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices, std::vector<LinkAtom> link_atoms);

}

#endif
