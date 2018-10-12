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

		if (file_exists(Config::get().get().general.paramFilename))  // if parameterfile -> equilibrium distance from forcefield
		{
			auto b_type_qm = tp.type(coords->atoms().atom(b).energy_type(), tinker::potential_keys::BOND); // bonding energy type for QM atom
			auto b_type = tp.type(energy_type, tinker::potential_keys::BOND);                              // bonding energy type for link atom
			for (auto b_param : tp.bonds())
			{
				if (b_param.index[0] == b_type_qm && b_param.index[1] == b_type)  deq_L_QM = b_param.ideal;
				else if (b_param.index[0] == b_type && b_param.index[1] == b_type_qm) deq_L_QM = b_param.ideal;
			}
		}
    
    if (deq_L_QM == 0.0)   // equilibrium distance cannot be determined by forcefield -> sum of covalent radii
    {
			if (Config::get().general.verbosity > 3)
			{
				std::cout << "determining link atom position from sum of covalent radii because forcefield parameter not available\n";
			}

      size_t atomic_number_QM = atomic::atomic_number_by_symbol(coords->atoms().atom(b).symbol());
      size_t atomic_number_LA = atomic::atomic_number_by_symbol("H");  // only H-atoms

      deq_L_QM = atomic::cov_radiusMap[atomic_number_QM] + atomic::cov_radiusMap[atomic_number_LA];
    }

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
  @param tp: tinker parameter object*/
  std::vector<LinkAtom> create_link_atoms(coords::Coordinates* coords, std::vector<size_t> &qm_indices, tinker::parameter::parameters const &tp);

  /**calculate gradients on QM and MM atom from link atom (see DOI 10.1002/(SICI)1096-987X(199703)18:4<463::AID-JCC2>3.0.CO;2-R)
  @param l: link atom
  @param G_L: gradient on link atom
  @param coords: coordinates object of whole system
  @param G_QM: reference to QM gradient (is calculated)
  @param G_MM: reference to MM gradient (is calculated)*/
  void calc_link_atom_grad(LinkAtom &l, coords::r3 const &G_L, coords::Coordinates* coords, coords::r3 &G_QM, coords::r3 &G_MM);

  /**creates a vector with the indizes of all MM atoms
  @param num_atoms: number of atoms in whole system*/
  std::vector<std::size_t> get_mm_atoms(std::size_t const num_atoms);

  /**creates a vector new_indices
  for description of the vector see energy_int_qmmm.h
  @param num_atoms: number of atoms in whole system
  @param mmi: vector with original indizes*/
  std::vector<std::size_t> make_new_indices(std::size_t const num_atoms, std::vector<std::size_t> const &indices);

  /**creates a coordobject for partial systems (if QM system link atoms are added, too)
  @param cp: pointer to original coordobject
  @param indices: indizes of atoms that should be in new coordobject
  @param new_indices: vector of length total number of atoms
                      only those elements are filled whose position corresponds to atoms of new coordobject
                      they are filled with successive numbers starting from 0
                      purpose: faciliate mapping between total coordinates object and subsystems
  @param link_atoms: vector with link atoms (if present)
  @param energy_interface: energy interface of new coordobject
	@param write_into_file: if true writes subsystem into tinkerfile
	@param filename: name of the tinkerfile (only when write_into_file = true)*/
  coords::Coordinates make_small_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices, config::interface_types::T energy_interface, bool const write_into_file = false,
    std::vector<LinkAtom> const &link_atoms = std::vector<LinkAtom>(), std::string const& filename = "qm_system.arc");

  /**selects only those charges from amber_charges vector which correspond to the indices
  all other charges are removed*/
  void select_from_ambercharges(std::vector<std::size_t> const & indices);

	/**adds external charges to the following calculations
	@param qm_indizes: indizes of current "QM system", can be the "real" QM system or the intermediate system of QM and SE atoms
	@param ignore_indizes: indizes of atoms that should be ignored
	@param charges: vector of charge values that might be added to the calculation
	@param indizes_of_charges: indizes of the charges in the overall coordinates object
	@param link_atoms: vector of link atoms for the current "QM system"
	@param charge_indizes: reference to a vector where the indizes of the atoms whose charges are taken into account are added
	@param coords: pointer to original coordobject*/
	void add_external_charges(std::vector<size_t> &qm_indizes, std::vector<size_t> &ignore_indizes, std::vector<double> &charges, std::vector<size_t> &indizes_of_charges, 
		std::vector<LinkAtom> &link_atoms, std::vector<int> &charge_indizes, coords::Coordinates *coords);
}

#endif
