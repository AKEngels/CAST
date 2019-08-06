#include <cstddef>

#include "energy_int_qmmm.h"
#include "Scon/scon_utility.h"

::tinker::parameter::parameters energy::interfaces::qmmm::QMMM::tp;

energy::interfaces::qmmm::QMMM::QMMM(coords::Coordinates * cp) :
  interface_base(cp),
  qm_indices(Config::get().energy.qmmm.qm_systems[0]),
  mm_indices(qmmm_helpers::get_mm_atoms(cp->size())),
  new_indices_qm(qmmm_helpers::make_new_indices(qm_indices, cp->size())),
  new_indices_mm(qmmm_helpers::make_new_indices(mm_indices, cp->size())),
	link_atoms(qmmm_helpers::create_link_atoms(qm_indices, cp, tp, Config::get().energy.qmmm.linkatom_sets[0])),
	qmc(qmmm_helpers::make_small_coords(cp, qm_indices, new_indices_qm, Config::get().energy.qmmm.qminterface, "QM system: ", Config::get().energy.qmmm.qm_to_file, link_atoms)),
  mmc(qmmm_helpers::make_small_coords(cp, mm_indices, new_indices_mm, Config::get().energy.qmmm.mminterface, "MM system: ")),
	index_of_QM_center(qmmm_helpers::get_index_of_QM_center(Config::get().energy.qmmm.centers[0], qm_indices, coords)),
  qm_energy(0.0), mm_energy(0.0), vdw_energy(0.0), bonded_energy(0.0), coulomb_energy(0.0)
{
  if (!tp.valid())
  {
    tp.from_file(Config::get().get().general.paramFilename);
  }
	

	if (Config::get().periodics.periodic)
	{
		if (Config::get().energy.qmmm.mminterface == config::interface_types::T::OPLSAA || Config::get().energy.qmmm.mminterface == config::interface_types::T::AMBER)
		{
			double const min_cut = std::min({ Config::get().periodics.pb_box.x(), Config::get().periodics.pb_box.y(), Config::get().periodics.pb_box.z() }) / 2.0;
			if (Config::get().energy.cutoff > min_cut)
			{
				std::cout << "\n!!! WARNING! Forcefield cutoff too big! Your cutoff should be smaller than " << min_cut << "! !!!\n\n";
			}
		}

		double const min_cut = std::min({ Config::get().periodics.pb_box.x(), Config::get().periodics.pb_box.y(), Config::get().periodics.pb_box.z() }) / 2.0;
		if (Config::get().energy.qmmm.cutoff > min_cut)
		{
			std::cout << "\n!!! WARNING! QM/MM cutoff too big! Your cutoff should be smaller than " << min_cut << "! !!!\n\n";
		}
	}

  std::vector<std::size_t> types;
  for (auto atom : (*cp).atoms())
  {
    scon::sorted::insert_unique(types, atom.energy_type());
  }
  cparams = tp.contract(types);
  torsionunit = cparams.torsionunit();
  prepare_bonded_qmmm();
}

energy::interfaces::qmmm::QMMM::QMMM(QMMM const & rhs,
  coords::Coordinates *cobj) : interface_base(cobj),
  cparams(rhs.cparams), qm_indices(rhs.qm_indices), mm_indices(rhs.mm_indices),
  new_indices_qm(rhs.new_indices_qm), new_indices_mm(rhs.new_indices_mm), link_atoms(rhs.link_atoms),
  qmc(rhs.qmc), mmc(rhs.mmc), qm_energy(rhs.qm_energy), mm_energy(rhs.mm_energy), 
  vdw_energy(rhs.vdw_energy), bonded_energy(rhs.bonded_energy), coulomb_energy(rhs.coulomb_energy),
  c_gradient(rhs.c_gradient), vdw_gradient(rhs.vdw_gradient), bonded_gradient(rhs.bonded_gradient)
{
  interface_base::operator=(rhs);
}

energy::interfaces::qmmm::QMMM::QMMM(QMMM&& rhs, coords::Coordinates *cobj)
  : interface_base(cobj),
  cparams(std::move(rhs.cparams)), 
  qm_indices(std::move(rhs.qm_indices)), mm_indices(std::move(rhs.mm_indices)),
  new_indices_qm(std::move(rhs.new_indices_qm)),
  new_indices_mm(std::move(rhs.new_indices_mm)),
  link_atoms(std::move(rhs.link_atoms)),
  qmc(std::move(rhs.qmc)), mmc(std::move(rhs.mmc)),
  qm_energy(std::move(rhs.qm_energy)), mm_energy(std::move(rhs.mm_energy)), vdw_energy(std::move(rhs.vdw_energy)),
  bonded_energy(std::move(rhs.bonded_energy)), coulomb_energy(std::move(rhs.coulomb_energy)),
  c_gradient(std::move(rhs.c_gradient)), 
  vdw_gradient(std::move(rhs.vdw_gradient)), bonded_gradient(std::move(rhs.bonded_gradient))
{
  interface_base::operator=(rhs);
}

/**function to find force field parameters for bonds, angles and so on between QM and MM system*/
void energy::interfaces::qmmm::QMMM::find_parameters()
{
  for (auto &b : qmmm_bonds)  // find bond parameters
  {
	  bool found = false;
    auto b_type_a = cparams.type(coords->atoms().atom(b.a).energy_type(), tinker::potential_keys::BOND);
    auto b_type_b = cparams.type(coords->atoms().atom(b.b).energy_type(), tinker::potential_keys::BOND);
    for (auto b_param : cparams.bonds())
    {
      if (b_param.index[0] == b_type_a && b_param.index[1] == b_type_b)
      {
        b.ideal = b_param.ideal;
        b.force = b_param.f;
		    found = true;
      }
      else if (b_param.index[0] == b_type_b && b_param.index[1] == b_type_a)
      {
        b.ideal = b_param.ideal;
        b.force = b_param.f;
		    found = true;
      }
    }
	  if (found == false) std::cout << "Parameters for bond between " << b.a+1 << " and " << b.b+1 << " not found.\n";
  }

  for (auto &a : qmmm_angles)  // find angle parameters
  {
	  bool found = false;
    auto a_type_a = cparams.type(coords->atoms().atom(a.a).energy_type(), tinker::potential_keys::ANGLE);
    auto a_type_b = cparams.type(coords->atoms().atom(a.b).energy_type(), tinker::potential_keys::ANGLE);
    auto a_type_c = cparams.type(coords->atoms().atom(a.c).energy_type(), tinker::potential_keys::ANGLE);
    for (auto a_param : cparams.angles())
    {
      if (a_param.index[1] == a_type_c)
      {
        if (a_param.index[0] == a_type_a && a_param.index[2] == a_type_b)
        {
          a.ideal = a_param.ideal;
          a.force = a_param.f;
		      found = true;
        }
        else if (a_param.index[0] == a_type_b && a_param.index[2] == a_type_a)
        {
          a.ideal = a_param.ideal;
          a.force = a_param.f;
		      found = true;
        }
      }
    }
	  if (found == false) std::cout << "Parameters for angle made of " << a.a+1 <<", "<<a.c+1 << " and " << a.b+1 << " not found.\n";
  }

  for (auto &d : qmmm_dihedrals)  // find parameters for dihedrals
  {
	  bool found = false;
    auto d_type_a = cparams.type(coords->atoms().atom(d.a).energy_type(), tinker::potential_keys::TORSION);
    auto d_type_b = cparams.type(coords->atoms().atom(d.b).energy_type(), tinker::potential_keys::TORSION);
    auto d_type_c1 = cparams.type(coords->atoms().atom(d.c1).energy_type(), tinker::potential_keys::TORSION);
    auto d_type_c2 = cparams.type(coords->atoms().atom(d.c2).energy_type(), tinker::potential_keys::TORSION);
    for (auto d_param : cparams.torsions())  // find parameters only with real atom types (no 0)
    {
      if (d_param.index[0] == d_type_a && d_param.index[1] == d_type_c1 && d_param.index[2] == d_type_c2 && d_param.index[3] == d_type_b)
      {
        d.max_order = d_param.max_order;
        d.number = d_param.number;
        d.orders = d_param.order;
        d.forces = d_param.force;
        d.ideals = d_param.ideal;
		    found = true;
      }
      else if (d_param.index[0] == d_type_b && d_param.index[1] == d_type_c2 && d_param.index[2] == d_type_c1 && d_param.index[3] == d_type_a)
      {
        d.max_order = d_param.max_order;
        d.number = d_param.number;
        d.orders = d_param.order;
        d.forces = d_param.force;
        d.ideals = d_param.ideal;
		    found = true;
      }
    }
	  if (found == false)   // find parameters where one of the outer atoms is a 0 instead of the "real" atom type
	  {
		  for (auto d_param : cparams.torsions())
		  {
			  if (d_param.index[0] == 0 && d_param.index[1] == d_type_c1 && d_param.index[2] == d_type_c2 && d_param.index[3] == d_type_b)
			  {
				  d.max_order = d_param.max_order;
				  d.number = d_param.number;
				  d.orders = d_param.order;
				  d.forces = d_param.force;
				  d.ideals = d_param.ideal;
				  found = true;
			  }
			  else if (d_param.index[0] == d_type_a && d_param.index[1] == d_type_c1 && d_param.index[2] == d_type_c2 && d_param.index[3] == 0)
			  {
				  d.max_order = d_param.max_order;
				  d.number = d_param.number;
				  d.orders = d_param.order;
				  d.forces = d_param.force;
				  d.ideals = d_param.ideal;
				  found = true;
			  }
			  else if (d_param.index[0] == 0 && d_param.index[1] == d_type_c2 && d_param.index[2] == d_type_c1 && d_param.index[3] == d_type_a)
			  {
				  d.max_order = d_param.max_order;
				  d.number = d_param.number;
				  d.orders = d_param.order;
				  d.forces = d_param.force;
				  d.ideals = d_param.ideal;
				  found = true;
			  }
			  else if (d_param.index[0] == d_type_b && d_param.index[1] == d_type_c2 && d_param.index[2] == d_type_c1 && d_param.index[3] == 0)
			  {
				  d.max_order = d_param.max_order;
				  d.number = d_param.number;
				  d.orders = d_param.order;
				  d.forces = d_param.force;
				  d.ideals = d_param.ideal;
				  found = true;
			  }
		  }
	  }
    if (found == false)
    {
      for (auto d_param : cparams.torsions())  // find parameters where both outer atoms are 0
      {
        if (d_param.index[0] == 0 && d_param.index[1] == d_type_c1 && d_param.index[2] == d_type_c2 && d_param.index[3] == 0)
        {
          d.max_order = d_param.max_order;
          d.number = d_param.number;
          d.orders = d_param.order;
          d.forces = d_param.force;
          d.ideals = d_param.ideal;
          found = true;
        }
        else if (d_param.index[0] == 0 && d_param.index[1] == d_type_c2 && d_param.index[2] == d_type_c1 && d_param.index[3] == 0)
        {
          d.max_order = d_param.max_order;
          d.number = d_param.number;
          d.orders = d_param.order;
          d.forces = d_param.force;
          d.ideals = d_param.ideal;
          found = true;
        }
      }
    }
	  if (found == false)
	  {
		  std::cout << "Parameters for dihedral made of " << d.a+1 << ", " << d.c1+1 << ", " << d.c2+1 << " and " << d.b+1 << " not found.\n";
	  }
  }
}

/**function to find bonds, angles and so on between QM and MM system*/
void energy::interfaces::qmmm::QMMM::find_bonds_etc()
{
  // find bonds between QM and MM region
  for (auto mma : mm_indices)
  {
    for (auto b : coords->atoms().atom(mma).bonds())
    {
      if (scon::sorted::exists(qm_indices, b))
      {
        bonded::Bond bond(mma, b);
        qmmm_bonds.push_back(bond);
      }
    }
  }

  // find angles between QM and MM region
  for (auto b : qmmm_bonds) // for every bond
  {
    for (auto p : coords->atoms().atom(b.a).bonds()) // atom a as angle center
    {
      if (b.b != p)
      {
        bonded::Angle angle(b.b, p, b.a);
        if (!bonded::is_in(angle, qmmm_angles))
        {
          qmmm_angles.push_back(angle);
        }
      }
    }
    for (auto p : coords->atoms().atom(b.b).bonds()) // atom b as angle center
    {
      if (b.a != p)
      {
        bonded::Angle angle(b.a, p, b.b);
        if (!bonded::is_in(angle, qmmm_angles))
        {
          qmmm_angles.push_back(angle);
        }
      }
    }
  }

  // find dihedrals between QM and MM region
  for (auto a : qmmm_angles)
  {
    for (auto p : coords->atoms().atom(a.a).bonds())   // expand angle at atom a
    {
      if (a.c != p)
      {
        bonded::Dihedral dihed(p, a.b, a.a, a.c);
        if (!bonded::is_in(dihed, qmmm_dihedrals))
        {
          qmmm_dihedrals.push_back(dihed);
        }
      }
    }
    for (auto p : coords->atoms().atom(a.b).bonds())   // expand angle at atom b
    {
      if (a.c != p)
      {
        bonded::Dihedral dihed(p, a.a, a.b, a.c);
        if (!bonded::is_in(dihed, qmmm_dihedrals))
        {
          qmmm_dihedrals.push_back(dihed);
        }
      }
    }
  }
}

void energy::interfaces::qmmm::QMMM::prepare_bonded_qmmm()
{
  find_bonds_etc();    // find bonds, angles and dihedrals between QM and MM region
  find_parameters();   // find force field parameters for energy calculation

  if (Config::get().general.verbosity > 3)  // Output
  {
	  std::cout << "QM/MM-Bonds\n";
	  for (auto b : qmmm_bonds)
	  {
		  std::cout << b.info() << "\n";
	  }
	  std::cout << "QM/MM-Angles\n";
	  for (auto a : qmmm_angles)
	  {
		  std::cout << a.info() << "\n";
	  }
	  std::cout << "QM/MM-Dihedrals\n";
	  for (auto a : qmmm_dihedrals)
	  {
		  std::cout << a.info() << "\n";
	  }
  }
}

void energy::interfaces::qmmm::QMMM::update_representation()
{
  std::size_t qi = 0u;   // update positions of QM system
  for (auto i : qm_indices)
  {
    qmc.move_atom_to(qi, coords->xyz()[i], true);
    ++qi;
  }

	for (auto &l : link_atoms) l.calc_position(coords); // update positions of link atoms in QM system
	for (auto i = 0u; i < link_atoms.size(); ++i)
	{
		int index = qm_indices.size() + i;
		coords::cartesian_type &new_pos = link_atoms[i].position;
		qmc.move_atom_to(index, new_pos, true);
	}

  std::size_t mi = 0u;       // update positions of MM system
  for (auto j : mm_indices)
  {
    mmc.move_atom_to(mi, coords->xyz()[j], true);
    ++mi;
  }

}

/**calculates energies and gradients
@param if_gradient: true if gradients should be calculated, false if not*/
coords::float_type energy::interfaces::qmmm::QMMM::qmmm_calc(bool if_gradient)
{
  integrity = true;
	if (link_atoms.size() != Config::get().energy.qmmm.linkatom_sets[0].size())  // test if correct number of link atom types is given
	{                                                                          // can't be done in constructor because interface is first constructed without atoms 
		std::cout << "Wrong number of link atom types given. You have " << link_atoms.size() << " in the following order:\n";
		for (auto &l : link_atoms)
		{
			std::cout << "QM atom: " << l.qm + 1 << ", MM atom: " << l.mm + 1 << "\n";
		}
		throw std::runtime_error("wrong number of link atom types");
	}
  
  // ############ UPDATE STUFF ##############################

  update_representation();  // update positions of QM and MM subsystem to those of coordinates object
  for (auto &l : link_atoms) l.calc_position(coords); // update positions of link atoms

	// ############### CREATE MM CHARGES ######################

	if (Config::get().coords.amber_charges.size() > mm_indices.size())
	{
		total_amber_charges = Config::get().coords.amber_charges;  // save all amber charges (for mechanical embedding)
		qmmm_helpers::select_from_ambercharges(mm_indices);
	}

  if (Config::get().energy.qmmm.zerocharge_bonds != 0)
  {
    std::vector<double> mm_charge_vector = mmc.energyinterface()->charges();

    charge_indices.clear();
    qmmm_helpers::add_external_charges(qm_indices, mm_charge_vector, mm_indices, link_atoms, charge_indices, coords, index_of_QM_center);
  }
	
  // ################### DO CALCULATION ###########################################

  if (Config::get().energy.qmmm.qminterface == config::interface_types::T::MOPAC) {                 // these QM interfaces are allowed
    Config::set().energy.mopac.link_atoms = Config::get().energy.qmmm.linkatom_sets[0].size(); // set number of link atoms for MOPAC
  }   
  else if (Config::get().energy.qmmm.qminterface == config::interface_types::T::GAUSSIAN) {} 
	else if (Config::get().energy.qmmm.qminterface == config::interface_types::T::DFTB) {}
  else if (Config::get().energy.qmmm.qminterface == config::interface_types::T::PSI4) {}
	else if (Config::get().energy.qmmm.qminterface == config::interface_types::T::ORCA) {}
  else throw std::runtime_error("Chosen QM interface not implemented for QM/MM!");

	bool periodic = Config::get().periodics.periodic;   // switch off periodic boundaries for QM calculation
	Config::set().periodics.periodic = false;

  try {
    if (if_gradient)
    {
      qm_energy = qmc.g();  // get energy for QM part and save gradients for QM part
      if (Config::get().energy.qmmm.zerocharge_bonds != 0) {
        g_coul_mm = qmc.energyinterface()->get_g_ext_chg();  // coulomb gradients on MM atoms
      }
    }
    else qm_energy = qmc.e();  // get only energy for QM part
  }
  catch(...)
  {
	  std::cout << "QM programme failed. Treating structure as broken.\n";
    integrity = false;  // if QM programme fails: integrity is destroyed
  }
  Config::set().energy.qmmm.mm_charges.clear();  // no external charges for MM calculation
	Config::set().periodics.periodic = periodic;   // reset periodic boundaries

  ww_calc(if_gradient);  // calculate interactions between QM and MM part

  if (integrity == true)
  {
    if (if_gradient)  // if gradients should be calculated
    {
      mm_energy = mmc.g(); // get energy for MM part

      // get gradients: QM + MM + vdW + Coulomb + bonded
	    auto new_grad = vdw_gradient + c_gradient + bonded_gradient;  // vdW + Coulomb + bonded
      auto g_qm = qmc.g_xyz(); // QM
      auto g_mm = mmc.g_xyz(); // MM

      // calculate gradients from link atoms
      for (auto j=0u; j<link_atoms.size(); j++)
      {
        LinkAtom l = link_atoms[j];
				auto link_atom_grad = g_qm[qm_indices.size() + j];
        coords::r3 G_QM, G_MM;

        qmmm_helpers::calc_link_atom_grad(l, link_atom_grad, coords, G_QM, G_MM);

        new_grad[l.qm] += G_QM;
        new_grad[l.mm] += G_MM;

        if (Config::get().general.verbosity > 4)
        {
		      std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
		      std::cout<< "It causes a gradient on QM atom " << G_QM << " and on MM atom " << G_MM << ".\n";
        }
      }

      // get gradients from pure QM and MM system
			for (auto&& qmi : qm_indices)
			{
				new_grad[qmi] += g_qm[new_indices_qm[qmi]];
			}
			for (auto && mmi : mm_indices)
			{
				new_grad[mmi] += g_mm[new_indices_mm[mmi]];
			}
      coords->swap_g_xyz(new_grad);
    }
    else  // only energy
    {
      mm_energy = mmc.e();  // get energy for MM part
    }

    // energy = QM + MM + vdW + bonded (+ coulomb)
    energy = qm_energy + mm_energy + vdw_energy + bonded_energy + coulomb_energy;
		if (coords->check_bond_preservation() == false) integrity = false;
		else if (coords->check_for_crashes() == false) integrity = false;
  }
  return energy;
}

/**calculate bonded energy and gradients*/
double energy::interfaces::qmmm::QMMM::calc_bonded(bool if_gradient)
{
  double E(0.0);
  if (if_gradient == false)  // only energy calculation
  {   // bonded gradient is given to functions because the function needs a parameter, is not used
    for (auto b : qmmm_bonds) E += b.calc_energy(coords, bonded_gradient);
    for (auto a : qmmm_angles) E += a.calc_energy(coords, bonded_gradient);
    for (auto d : qmmm_dihedrals) E += d.calc_energy(coords, torsionunit, bonded_gradient);
  }
  else   // gradient calculation
  {
    for (auto b : qmmm_bonds) E += b.calc_energy(coords, bonded_gradient, true);
    for (auto a : qmmm_angles) E += a.calc_energy(coords, bonded_gradient, true);
    for (auto d : qmmm_dihedrals) E += d.calc_energy(coords, torsionunit, bonded_gradient, true);
  }
  return E;
}

/**determines if a van der waals interaction between a QM and a MM atom should be calculated
@param qm: index of QM atom
@param mm: index of MM atom
returns 0 if no vdW is calculated (1 or 2 bonds between the atoms), 1 if vdW is calculated normally
and 2 if vdW is scaled down by 1/2 (3 bonds between the atoms)*/
int energy::interfaces::qmmm::QMMM::calc_vdw(unsigned qm, unsigned mm)
{
  for (auto b : qmmm_bonds)
  {
    if (qm == b.b && mm == b.a) return 0;
  }
  for (auto a : qmmm_angles)
  {
    if (qm == a.a && mm == a.b) return 0;
    else if (qm == a.b && mm == a.a) return 0;
  }
  for (auto d : qmmm_dihedrals)
  {
    if (qm == d.a && mm == d.b) return 2;
    else if (qm == d.b && mm == d.a) return 2;
  }
  return 1;
}

/**calculates interaction between QM and MM part
energy is only vdW interactions
for MOPAC gradients are coulomb and vdW
for all the other QM interfaces gradients are vdW and coulomb on MM atoms
@param if_gradient: true if gradients should be calculated, false if not*/
void energy::interfaces::qmmm::QMMM::ww_calc(bool if_gradient)
{
  // bonded interactions
  bonded_gradient.assign(coords->size(), coords::r3{});
  bonded_energy = calc_bonded(if_gradient);

  // preparation for calculation of non-bonded interactions
  std::vector<double> qm_charge_vector;                                     // vector with all charges of QM atoms
  std::vector<double> mm_charge_vector = mmc.energyinterface()->charges();  // vector with all charges of MM atoms
  try {
    qm_charge_vector = qmc.energyinterface()->charges(); // still link atoms in it
    for (auto i = 0u; i < link_atoms.size(); ++i)
    {                                    // remove charges from link atoms
      qm_charge_vector.pop_back();
    }
  }
  catch (...) { integrity = false; }

  if (integrity == true)
  {
    //########## reset vdw gradients and energies #################################

    c_gradient.assign(coords->size(), coords::r3{});     
    vdw_gradient.assign(coords->size(), coords::r3{});
    vdw_energy = 0.0;
		coulomb_energy = 0.0;

		// ########## calculate vdw interactions ##########################################

    auto vdw_params = cparams.vdws();  // get vdw parameters

		std::size_t i2 = 0u;
    for (auto i : qm_indices)  // for every QM atom
    {
      auto z = qmc.atoms(i2).energy_type();  // get atom type
      std::size_t i_vdw = cparams.type(z, tinker::potential_keys::VDW);  // get index to find this atom type in vdw parameters
      auto vparams_i = vdw_params[i_vdw - 1];  // get vdw parameters for QM atom
      std::size_t j2 = 0u;
      for (auto j : mm_indices)  // for every MM atom
      {
        auto e_type = mmc.atoms(j2).energy_type();  // get atom type
        std::size_t j_vdw = cparams.type(e_type, tinker::potential_keys::VDW); // get index to find this atom type in vdw parameters
        auto vparams_j = vdw_params[j_vdw - 1];  // get vdw parameters for MM atom

        auto r_ij = coords->xyz(j) - coords->xyz(i);            // distance between QM and MM atom
        if (Config::get().periodics.periodic) boundary(r_ij);   // if periodic boundaries: take shortest distance
        coords::float_type d = len(r_ij);

				double scaling = 1.0;     // determine scaling factor, default: no scaling
				if (Config::get().energy.cutoff < std::numeric_limits<double>::max())   
				{
					double const& c = Config::get().energy.cutoff;       // cutoff distance
					double const& s = Config::get().energy.switchdist;   // distance where cutoff starts to kick in (only vdW)

					if (d > c) scaling = 0.0;
					else if (d > s) scaling = ( (c*c - d*d)* (c * c - d * d) * (c*c + 2*d*d - 3*s*s) ) / ((c * c - s * s) * (c * c - s * s) * (c * c - s * s));
				}

        double R_0;  // r_min or sigma
        if (cparams.general().radiustype.value ==
          ::tinker::parameter::radius_types::T::SIGMA)
        {
          R_0 = sqrt(vparams_i.r*vparams_j.r);  // sigma
        }
        else if (cparams.general().radiustype.value ==
          ::tinker::parameter::radius_types::T::R_MIN)
        {
          R_0 = vparams_i.r + vparams_j.r;  // r_min
        }
        else throw std::runtime_error("no valid radius_type");

        double epsilon = sqrt(vparams_i.e * vparams_j.e);  // epsilon
        auto R_r = std::pow(R_0 / d, 6);

        int calc_modus = calc_vdw(i, j);  // will the vdw interaction be calculated? if yes, will it be scaled down?
        if (Config::get().general.verbosity > 4)
        {
          std::cout << "VdW calc_modus between atoms " << i + 1 << " and " << j + 1 << " is " << calc_modus << ".\n";
        }

        double vdw;  // vdw energy for current atom pair
        if (calc_modus != 0)  // calculate vdW interaction
        {
          if (cparams.general().radiustype.value ==
            ::tinker::parameter::radius_types::T::SIGMA)
          {
            vdw = 4 * R_r * epsilon*(R_r - 1.0) * scaling;
            if (calc_modus == 2) vdw = vdw / 2;
          }
          else if (cparams.general().radiustype.value ==
            ::tinker::parameter::radius_types::T::R_MIN)
          {
            vdw = R_r * epsilon*(R_r - 2.0) * scaling;
            if (calc_modus == 2) vdw = vdw / 2;
          }
          else
          {
            throw std::runtime_error("no valid radius_type");
          }
          vdw_energy += vdw;
        }

        if (if_gradient)  // gradients
        {
          if (calc_modus != 0)  // gradients of vdW interaction
          {
            if (cparams.general().radiustype.value
              == ::tinker::parameter::radius_types::T::SIGMA)
            {
              coords::float_type const V = 4 * epsilon * R_r;
              auto vdw_r_grad_sigma = (V / d)*(6.0 - 12.0 * R_r);
              auto vdw_gradient_ij_sigma = ((r_ij*vdw_r_grad_sigma) / d) * scaling;

							if (d > Config::get().energy.switchdist && d <= Config::get().energy.cutoff)  // additional gradient as charge changes with distance
							{
								double const& c = Config::get().energy.cutoff;       
								double const& s = Config::get().energy.switchdist;  

								auto E_unscaled = 4 * R_r * epsilon * (R_r - 1.0);
								auto deriv_S = (-12*d* (c*c - d*d) * (d*d-s*s)) / ( (c*c - s*s)* (c * c - s * s)* (c * c - s * s));
								auto abs_grad = E_unscaled * deriv_S;

								vdw_gradient_ij_sigma += (r_ij / d) * abs_grad;  // give additional gradient a direction
							}
							
              if (calc_modus == 2) vdw_gradient_ij_sigma = vdw_gradient_ij_sigma / 2;
              vdw_gradient[i] -= vdw_gradient_ij_sigma;
              vdw_gradient[j] += vdw_gradient_ij_sigma;
            }
            else
            {
              coords::float_type const V = epsilon * R_r;
              auto vdw_r_grad_R_MIN = (V / d) * 12 * (1.0 - R_r);
              auto vdw_gradient_ij_R_MIN = ((r_ij*vdw_r_grad_R_MIN) / d) * scaling;

							if (d > Config::get().energy.switchdist && d <= Config::get().energy.cutoff)  // additional gradient as charge changes with distance
							{
								double const& c = Config::get().energy.cutoff;
								double const& s = Config::get().energy.switchdist;

								auto E_unscaled = R_r * epsilon * (R_r - 2.0);
								auto deriv_S = (-12 * d * (c * c - d * d) * (d * d - s * s)) / ((c * c - s * s) * (c * c - s * s) * (c * c - s * s));
								auto abs_grad = E_unscaled * deriv_S;

								vdw_gradient_ij_R_MIN += (r_ij / d) * abs_grad;  // give additional gradient a direction
							}

              if (calc_modus == 2) vdw_gradient_ij_R_MIN = vdw_gradient_ij_R_MIN / 2;
              vdw_gradient[i] -= vdw_gradient_ij_R_MIN;
              vdw_gradient[j] += vdw_gradient_ij_R_MIN;
            }
          }
        }
        ++j2;
      }
      ++i2;
    }

		// ########## calculate coulomb interactions ##########################################

    if (if_gradient == true && Config::get().energy.qmmm.zerocharge_bonds != 0)  // electrostatic embedding -> coulomb gradients on MM atoms
    {
      for (auto i = 0u; i < charge_indices.size(); ++i)
      {
        int mma = charge_indices[i];
				auto grad = g_coul_mm[i];

				coords::r3 deriv_Q{ 0.0, 0.0, 0.0 };  // additional gradients because charge also changes with position
				if (Config::get().energy.qmmm.cutoff != std::numeric_limits<double>::max())
				{
					double constexpr elec_factor = 332.06;
					double const& c = Config::get().energy.qmmm.cutoff;
					double const& ext_chg = Config::get().energy.qmmm.mm_charges[i].original_charge;
					double const& scaling = Config::get().energy.qmmm.mm_charges[i].scaled_charge / Config::get().energy.qmmm.mm_charges[i].original_charge;
					auto chargesQM = qmc.energyinterface()->charges();

					double sum_of_QM_interactions{ 0.0 };   // calculate sum(Q_qm * Q_ext / r)
					for (auto j{ 0u }; j < qm_indices.size(); ++j)
					{
						double const QMcharge = chargesQM[j];
						coords::r3 MMpos{ Config::get().energy.qmmm.mm_charges[i].x,  Config::get().energy.qmmm.mm_charges[i].y,  Config::get().energy.qmmm.mm_charges[i].z };
						double const dist = len(MMpos - coords->xyz(qm_indices[j]));
						sum_of_QM_interactions += (QMcharge * ext_chg * elec_factor) / dist;
					}

					// calculate additional gradient
					deriv_Q.x() = sum_of_QM_interactions * 4 * (coords->xyz(index_of_QM_center).x() - Config::get().energy.qmmm.mm_charges[i].x) * std::sqrt(scaling) / (c * c);
					deriv_Q.y() = sum_of_QM_interactions * 4 * (coords->xyz(index_of_QM_center).y() - Config::get().energy.qmmm.mm_charges[i].y) * std::sqrt(scaling) / (c * c);
					deriv_Q.z() = sum_of_QM_interactions * 4 * (coords->xyz(index_of_QM_center).z() - Config::get().energy.qmmm.mm_charges[i].z) * std::sqrt(scaling) / (c * c);

					grad += deriv_Q;                             // gradient on external charge
					c_gradient[index_of_QM_center] -= deriv_Q;   // gradient on center of QM region
				}

				c_gradient[mma] += grad;
      }
    }

    else if (Config::get().energy.qmmm.zerocharge_bonds == 0)  // mechanical embedding -> coulomb energy and gradients from MM interface
    {
      auto c_params = cparams.charges();     // get forcefield parameters (charge)
      double current_coul_energy, current_coul_grad, qm_charge;  // some variables

      auto i2{ 0u };
      for (auto i : qm_indices)  // for every QM atom
      {
        if (Config::get().general.input == config::input_types::AMBER || Config::get().general.chargefile)  // amber charges
        {
          qm_charge = total_amber_charges[i] / 18.2223;
        }
        else   // normally (i.e. OPLSAA)
        {
          auto z = qmc.atoms(i2).energy_type();  // get atom type
          std::size_t i_coul = cparams.type(z, tinker::potential_keys::CHARGE);  // get index to find this atom type in charge parameters
          qm_charge = c_params[i_coul - 1].c;  // get charge parameter for QM atom
        }
        
        auto j2{ 0u };
        for (auto j : mm_indices)
        {
          double mm_charge = mm_charge_vector[j2];    // MM charges were filled before
          int calc_modus = calc_vdw(i, j);            // calc modus for coulomb interactions is the same as for vdw
          
          current_coul_energy = 0.0;
          current_coul_grad = 0.0;
          if (calc_modus != 0)
          {
            auto r_ij = coords->xyz(j) - coords->xyz(i);            // distance between QM and MM atom
            if (Config::get().periodics.periodic) boundary(r_ij);   // if periodic boundaries: take shortest distance
            coords::float_type d = len(r_ij);

						double scaling = 1.0;     // determine scaling factor, default: no scaling
						if (Config::get().energy.cutoff < std::numeric_limits<double>::max())   
						{
							double const& c = Config::get().energy.cutoff;       // cutoff distance

							if (d > c) scaling = 0.0;
							else scaling = (1 - (d / c) * (d / c)) * (1 - (d / c) * (d / c));
						}

            current_coul_energy = ((qm_charge*mm_charge*cparams.general().electric) / d);   // calculate energy (unscaled)
            if (calc_modus == 2)
            {
              if (Config::get().energy.qmmm.mminterface == config::interface_types::T::OPLSAA) {
                current_coul_energy = current_coul_energy / 2;
              }
              else if(Config::get().energy.qmmm.mminterface == config::interface_types::T::AMBER) {
                current_coul_energy = current_coul_energy / 1.2;
              }
            }
						auto scaled_energy = current_coul_energy * scaling;
            coulomb_energy += scaled_energy;

            if (if_gradient == true)   // calculate gradient
            {
              current_coul_grad = scaled_energy / d;      // gradient without direction (scaling already included)
							if (d < Config::get().energy.cutoff)        // additional gradient because charge changes with distance
							{
								double const& c = Config::get().energy.cutoff;
								current_coul_grad -= current_coul_energy * 4 * d * (d * d - c * c) / (c * c * c * c);  // minus because the vector r_ij is the other way round
							}
              auto new_grad = (r_ij / d) * current_coul_grad;   // gradient gets a direction
              c_gradient[i] += new_grad;
              c_gradient[j] -= new_grad;
            }
          }
          j2++;
        }
        i2++;
      }
    }
  }
}


energy::interface_base * energy::interfaces::qmmm::QMMM::clone(coords::Coordinates * c) const
{
  QMMM * tmp = new QMMM(*this, c);
  return tmp;
}

energy::interface_base * energy::interfaces::qmmm::QMMM::move(coords::Coordinates * c)
{
  QMMM * tmp = new QMMM(std::move(*this), c);
  return tmp;
}


void energy::interfaces::qmmm::QMMM::swap(interface_base& rhs)
{
  swap(dynamic_cast<QMMM&>(rhs));
}

void energy::interfaces::qmmm::QMMM::swap(QMMM& rhs)
{
  interface_base::swap(rhs);
  std::swap(cparams, rhs.cparams);
  qm_indices.swap(rhs.qm_indices);
  mm_indices.swap(rhs.mm_indices);
  new_indices_mm.swap(rhs.new_indices_mm);
  new_indices_qm.swap(rhs.new_indices_qm);
	link_atoms.swap(rhs.link_atoms);
  qmc.swap(rhs.qmc);
  mmc.swap(rhs.mmc);
  std::swap(vdw_energy, rhs.vdw_energy);
  std::swap(qm_energy, rhs.qm_energy);
  std::swap(mm_energy, rhs.mm_energy);
  c_gradient.swap(rhs.c_gradient);
  vdw_gradient.swap(rhs.vdw_gradient);
}

// update structure (account for topology or rep change)
void energy::interfaces::qmmm::QMMM::update(bool const skip_topology)
{
  if (!skip_topology)
  {
    *this = QMMM(this->coords);
  }
  else
  {
    update_representation();
    qmc.energy_update(true);
    mmc.energy_update(true);
  }
}

coords::float_type energy::interfaces::qmmm::QMMM::g()
{
	integrity = coords->check_structure();
	if (integrity == true) return qmmm_calc(true);
	else return 0;
}

coords::float_type energy::interfaces::qmmm::QMMM::e()
{
	integrity = coords->check_structure();
	if (integrity == true) return qmmm_calc(false);
	else return 0;
}

coords::float_type energy::interfaces::qmmm::QMMM::h()
{
  throw std::runtime_error("no QMMM-function yet");
}

coords::float_type energy::interfaces::qmmm::QMMM::o()
{
  throw std::runtime_error("QMMM-cannot optimize");
}

std::vector<coords::float_type> energy::interfaces::qmmm::QMMM::charges() const
{
	std::vector<double> qm_charge_vector;                                   // vector with all charges of QM atoms
	std::vector<double> mm_charge_vector = mmc.energyinterface()->charges();  // vector with all charges of MM atoms
	qm_charge_vector = qmc.energyinterface()->charges(); // still link atoms in it
	for (auto i = 0u; i < link_atoms.size(); ++i)
	{                                    // remove charges from link atoms
		qm_charge_vector.pop_back();
	}

  std::vector<coords::float_type> v;
  v.assign(coords->size(), 0);
  std::size_t i = 0;
  std::size_t j = 0;
  if (qm_charge_vector.size() + mm_charge_vector.size() == coords->size())
  {
    for (auto && qmi : qm_indices)
    {
      v[qmi] = qm_charge_vector[i];
      i++;
    }

    for (auto && mmi : mm_indices)
    {
      v[mmi] = mm_charge_vector[j];
      j++;
    }

    if (v.size() != coords->size())
    {
      throw std::logic_error("Found " + std::to_string(v.size()) +
        " charges instead of " + std::to_string(coords->size()) + " charges.");
    }
  }
  else
  {
    throw std::logic_error("QM/MM-charges not yet calculated");
  }
  return v;
}

void energy::interfaces::qmmm::QMMM::print_E(std::ostream &) const
{
  throw std::runtime_error("no QMMM-function yet");
}

void energy::interfaces::qmmm::QMMM::print_E_head(std::ostream &S, bool const endline) const
{
  S << "QM-atoms: " << qm_indices.size() << '\n';
  S << "MM-atoms: " << mm_indices.size() << '\n';
  S << "Potentials\n";
  S << std::right << std::setw(24) << "QM";
  S << std::right << std::setw(24) << "MM";
  S << std::right << std::setw(24) << "VDW";
  S << std::right << std::setw(24) << "BONDED";
  S << std::right << std::setw(24) << "COUL";
  S << std::right << std::setw(24) << "TOTAL";
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::QMMM::print_E_short(std::ostream &S, bool const endline) const
{
  S << '\n';
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << qm_energy;
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << mm_energy;
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << vdw_energy;
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << bonded_energy;
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << coulomb_energy;
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << energy;
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::QMMM::to_stream(std::ostream &S) const
{
  S << '\n';
  interface_base::to_stream(S);
  throw std::runtime_error("no QMMM-function yet");
}