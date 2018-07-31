#include <cstddef>

#include "energy_int_3layer.h"
#include "scon_utility.h"

::tinker::parameter::parameters energy::interfaces::three_layer::THREE_LAYER::tp;

energy::interfaces::three_layer::THREE_LAYER::THREE_LAYER(coords::Coordinates *cp):
  energy::interface_base(cp), qm_indices(Config::get().energy.qmmm.qmatoms),
  qm_se_indices(add_vectors(qm_indices, Config::get().energy.qmmm.seatoms)),
  new_indices_qm(qmmm_helpers::make_new_indices(cp->size(), qm_indices)),
  new_indices_middle(qmmm_helpers::make_new_indices(cp->size(), qm_se_indices)),
  link_atoms_small(qmmm_helpers::create_link_atoms(cp, qm_indices, tp)),
  link_atoms_middle(qmmm_helpers::create_link_atoms(cp, qm_se_indices, tp)),
  qmc(qmmm_helpers::make_small_coords(cp, qm_indices, new_indices_qm, link_atoms_small, Config::get().energy.qmmm.qminterface, "small_system.arc")),
  sec_small(qmmm_helpers::make_small_coords(cp, qm_indices, new_indices_qm, link_atoms_small, Config::get().energy.qmmm.seinterface, "small_system.arc")),
  sec_middle(qmmm_helpers::make_small_coords(cp, qm_se_indices, new_indices_middle, link_atoms_middle, Config::get().energy.qmmm.seinterface, "intermediate_system.arc")),
  mmc_middle(qmmm_helpers::make_small_coords(cp, qm_se_indices, new_indices_middle, link_atoms_middle, Config::get().energy.qmmm.mminterface, "intermediate_system.arc")),
	mmc_big(qmmm_helpers::make_mmbig_coords(cp)),
  qm_energy(0.0), se_energy_small(0.0), se_energy_middle(0.0), mm_energy_middle(0.0), mm_energy_big(0.0)
{
	if ((Config::get().energy.qmmm.qminterface != config::interface_types::T::OPLSAA && Config::get().energy.qmmm.qminterface != config::interface_types::T::AMBER &&
		Config::get().energy.qmmm.qminterface != config::interface_types::T::DFTB && Config::get().energy.qmmm.qminterface != config::interface_types::T::GAUSSIAN
    && Config::get().energy.qmmm.qminterface != config::interface_types::T::PSI4 && Config::get().energy.qmmm.qminterface != config::interface_types::T::MOPAC)
		||
    (Config::get().energy.qmmm.seinterface != config::interface_types::T::OPLSAA && Config::get().energy.qmmm.seinterface != config::interface_types::T::AMBER &&
      Config::get().energy.qmmm.seinterface != config::interface_types::T::DFTB && Config::get().energy.qmmm.seinterface != config::interface_types::T::GAUSSIAN
      && Config::get().energy.qmmm.seinterface != config::interface_types::T::PSI4 && Config::get().energy.qmmm.seinterface != config::interface_types::T::MOPAC)
    ||
		(Config::get().energy.qmmm.mminterface != config::interface_types::T::OPLSAA && Config::get().energy.qmmm.mminterface != config::interface_types::T::AMBER &&
			Config::get().energy.qmmm.mminterface != config::interface_types::T::DFTB && Config::get().energy.qmmm.mminterface != config::interface_types::T::GAUSSIAN
      && Config::get().energy.qmmm.mminterface != config::interface_types::T::PSI4 && Config::get().energy.qmmm.mminterface != config::interface_types::T::MOPAC))
	{
		throw std::runtime_error("One of your chosen interfaces is not suitable for THREE_LAYER.");
	}
  if (!file_exists(Config::get().get().general.paramFilename) &&    // if forcefield is desired but no parameterfile is given -> throw error
    (Config::get().energy.qmmm.qminterface == config::interface_types::T::OPLSAA || Config::get().energy.qmmm.qminterface == config::interface_types::T::AMBER
      || Config::get().energy.qmmm.seinterface == config::interface_types::T::OPLSAA || Config::get().energy.qmmm.seinterface == config::interface_types::T::AMBER
      || Config::get().energy.qmmm.mminterface == config::interface_types::T::OPLSAA || Config::get().energy.qmmm.qminterface == config::interface_types::T::AMBER))
  {
    throw std::runtime_error("You need a tinker-like parameterfile for your chosen forcefield.");
  }

  if (!tp.valid() && file_exists(Config::get().get().general.paramFilename))
  {
    tp.from_file(Config::get().get().general.paramFilename);
  }
}

energy::interfaces::three_layer::THREE_LAYER::THREE_LAYER(THREE_LAYER const & rhs,
  coords::Coordinates *cobj) : interface_base(cobj),
  qm_indices(rhs.qm_indices), qm_se_indices(rhs.qm_se_indices), 
  new_indices_qm(rhs.new_indices_qm), new_indices_middle(rhs.new_indices_middle), link_atoms_small(rhs.link_atoms_small),
  link_atoms_middle(rhs.link_atoms_middle), qmc(rhs.qmc), sec_small(rhs.sec_small),
  sec_middle(rhs.sec_middle), mmc_middle(rhs.mmc_middle), mmc_big(rhs.mmc_big), 
  qm_energy(rhs.qm_energy), se_energy_small(rhs.se_energy_small), se_energy_middle(rhs.se_energy_middle), 
  mm_energy_middle(rhs.mm_energy_middle), mm_energy_big(rhs.mm_energy_big)
{
  interface_base::operator=(rhs);
}

energy::interfaces::three_layer::THREE_LAYER::THREE_LAYER(THREE_LAYER&& rhs, coords::Coordinates *cobj)
  : interface_base(cobj),
  qm_indices(std::move(rhs.qm_indices)), qm_se_indices(std::move(rhs.qm_se_indices)), 
  new_indices_qm(std::move(rhs.new_indices_qm)), new_indices_middle(std::move(rhs.new_indices_middle)), 
  link_atoms_small(std::move(rhs.link_atoms_small)), link_atoms_middle(std::move(rhs.link_atoms_middle)),
  qmc(std::move(rhs.qmc)), sec_small(std::move(rhs.sec_small)), sec_middle(std::move(rhs.sec_middle)), 
  mmc_middle(std::move(rhs.mmc_middle)), mmc_big(std::move(rhs.mmc_big)), qm_energy(std::move(rhs.qm_energy)),
  se_energy_small(std::move(rhs.se_energy_small)), se_energy_middle(std::move(rhs.se_energy_middle)), 
  mm_energy_middle(std::move(rhs.mm_energy_middle)), mm_energy_big(std::move(rhs.mm_energy_big))
{
  interface_base::operator=(rhs);
}


energy::interface_base * energy::interfaces::three_layer::THREE_LAYER::clone(coords::Coordinates * c) const
{
  THREE_LAYER * tmp = new THREE_LAYER(*this, c);
  return tmp;
}

energy::interface_base * energy::interfaces::three_layer::THREE_LAYER::move(coords::Coordinates * c)
{
  THREE_LAYER * tmp = new THREE_LAYER(std::move(*this), c);
  return tmp;
}


void energy::interfaces::three_layer::THREE_LAYER::swap(interface_base& rhs)
{
  swap(dynamic_cast<THREE_LAYER&>(rhs));
}

void energy::interfaces::three_layer::THREE_LAYER::swap(THREE_LAYER& rhs)
{
  interface_base::swap(rhs);
  qm_indices.swap(rhs.qm_indices);
  qm_se_indices.swap(rhs.qm_se_indices);
  link_atoms_small.swap(rhs.link_atoms_small);
  link_atoms_middle.swap(rhs.link_atoms_middle);
  new_indices_qm.swap(rhs.new_indices_qm);
  new_indices_middle.swap(rhs.new_indices_middle);
  qmc.swap(rhs.qmc);
  mmc_middle.swap(rhs.mmc_middle);
  mmc_big.swap(rhs.mmc_big);
  sec_small.swap(rhs.sec_small);
  sec_middle.swap(rhs.sec_middle);
  std::swap(qm_energy, rhs.qm_energy);
  std::swap(se_energy_small, rhs.se_energy_small);
  std::swap(se_energy_middle, rhs.se_energy_middle);
  std::swap(mm_energy_middle, rhs.mm_energy_middle);
  std::swap(mm_energy_big, rhs.mm_energy_big);
}

// update structure (account for topology or rep change)
void energy::interfaces::three_layer::THREE_LAYER::update(bool const skip_topology)
{
  if (!skip_topology)
  {
    *this = THREE_LAYER(this->coords);
  }
  else
  {
    update_representation();
    qmc.energy_update(true);
    sec_small.energy_update(true);
    sec_middle.energy_update(true);
    mmc_middle.energy_update(true);
    mmc_big.energy_update(true);
  }
}

void energy::interfaces::three_layer::THREE_LAYER::update_representation()
{
  std::size_t qi = 0u;    // update position of QM atoms in small systems
  for (auto i : qm_indices)
  {
    qmc.move_atom_to(qi, coords->xyz(i), true);
    sec_small.move_atom_to(qi, coords->xyz(i), true);
    ++qi;
  }
  for (auto &l : link_atoms_small) l.calc_position(coords); // update positions of link atoms in small systems
  for (auto i = 0u; i < link_atoms_small.size(); ++i)
  {
	  int index = qm_indices.size() + i;
	  coords::cartesian_type &new_pos = link_atoms_small[i].position;
	  qmc.move_atom_to(index, new_pos, true);
	  sec_small.move_atom_to(index, new_pos, true);
  }

  std::size_t si = 0u;    // update position of atoms in middle systems
  for (auto i : qm_se_indices)
  {
    sec_middle.move_atom_to(si, coords->xyz(i), true);
    mmc_middle.move_atom_to(si, coords->xyz(i), true);
    ++si;
  }
  for (auto &l : link_atoms_middle) l.calc_position(coords); // update positions of link atoms in middle systems
  for (auto i = 0u; i < link_atoms_middle.size(); ++i)
  {
    int index = qm_se_indices.size() + i;
    coords::cartesian_type &new_pos = link_atoms_middle[i].position;
    sec_middle.move_atom_to(index, new_pos, true);
    mmc_middle.move_atom_to(index, new_pos, true);
  }
     
  for (std::size_t mi = 0u; mi<coords->size(); mi++)   // update positions of all atoms in big system
  {
    mmc_big.move_atom_to(mi, coords->xyz()[mi], true);
  }

}

coords::float_type energy::interfaces::three_layer::THREE_LAYER::qmmm_calc(bool if_gradient)
{
	if (link_atoms_middle.size() != Config::get().energy.qmmm.linkatom_types.size())  // test if correct number of link atom types is given
	{                                                                          // can't be done in constructor because interface is first constructed without atoms 
		std::cout << "Wrong number of link atom types given. You have " << link_atoms_middle.size() << " in the following order:\n";
		for (auto &l : link_atoms_middle)
		{
			std::cout << "QM atom: " << l.qm + 1 << ", MM atom: " << l.mm + 1 << "\n";
	  }
		std::cout << "This is assuming you are using a forcefield for your big system. \nIf you want to use one for the intermediate system talk to a CAST developer!\n";
		throw std::runtime_error("wrong number of link atom types");
	}

  update_representation(); // update positions of QM and MM subsystems to those of coordinates object
  
	mm_energy_big = 0.0;     // set energies to zero
	mm_energy_middle = 0.0;
	se_energy_middle = 0.0;
	se_energy_small = 0.0;
	qm_energy = 0.0;
  coords::Gradients_3D new_grads;  // save gradients in case of gradient calculation

  // ############### MM ENERGY AND GRADIENTS FOR WHOLE SYSTEM ######################

	try {
		if (!if_gradient)
		{
			mm_energy_big = mmc_big.e();   // calculate MM energy of whole system
		}
		else   // gradient calculation
		{
			mm_energy_big = mmc_big.g();
			new_grads = mmc_big.g_xyz();
		}
		if (mmc_big.integrity() == false) integrity = false;
		if (Config::get().general.verbosity > 4)
		{
			std::cout << "MM energy of big system: \n";
			mmc_big.e_head_tostream_short(std::cout);
			mmc_big.e_tostream_short(std::cout);
		}
	}
	catch (...)
	{
		std::cout << "MM programme (for big system) failed. Treating structure as broken.\n";
		integrity = false;  // if MM programme fails: integrity is destroyed
	}

	// if program didn't calculate an energy: return zero-energy (otherwise CAST will break because it doesn't find charges)
	if (integrity == false) return 0.0;   

  // ############### CREATE EXTERNAL CHARGES FOR MIDDLE SYSTEM ######################

  std::vector<double> charge_vector = mmc_big.energyinterface()->charges();

  std::vector<int> charge_indices;  // indizes of all atoms that are in charge_vector

  bool use_charge;
  charge_indices.clear();
  for (auto i = 0u; i < coords->size(); ++i) // go through all atoms
  {
	  use_charge = true;
	  for (auto &l : link_atoms_middle) // ignore those atoms that are connected to a QM or SE atom...
	  {
		  if (l.mm == i) use_charge = false;
	  }
	  for (auto &qs : qm_se_indices) // ...and the QM and SE atoms
	  {
		  if (qs == i) use_charge = false;
	  }

	  if (use_charge)  // for the other 
	  {
			if (Config::get().energy.qmmm.cutoff != 0.0)  // if cutoff given: test if one QM or SE atom is nearer than cutoff
			{
				use_charge = false;
				for (auto qs : qm_se_indices)
				{
					auto dist = len(coords->xyz(i) - coords->xyz(qs));
					if (dist < Config::get().energy.qmmm.cutoff)
					{
						use_charge = true;
						break;
					}
				}
				
			}

			if (use_charge)  // if yes create a PointCharge and add it to vector
			{
				PointCharge new_charge;
				new_charge.charge = charge_vector[i];
				new_charge.set_xyz(coords->xyz(i).x(), coords->xyz(i).y(), coords->xyz(i).z());
				Config::set().energy.qmmm.mm_charges.push_back(new_charge);

				charge_indices.push_back(i);  // add index to charge_indices
			}
	  }
  }

	// ############### SE ENERGY AND GRADIENTS FOR MIDDLE SYSTEM ######################
	try {
		if (!if_gradient)
		{
			se_energy_middle = sec_middle.e();  // get se energy for intermediate part 
		}
		else  // gradient calculation
		{
			se_energy_middle = sec_middle.g();    // get energy and calculate gradients
			auto g_se_middle = sec_middle.g_xyz();        // get gradients
			for (auto&& qsi : qm_se_indices)
			{
				new_grads[qsi] += g_se_middle[new_indices_middle[qsi]];
			}

			for (auto i = 0u; i < link_atoms_middle.size(); ++i)   // take into account link atoms
			{
				LinkAtom l = link_atoms_middle[i];

				coords::r3 g_qm, g_mm;        // divide link atom gradient to QM and MM atom
				auto link_atom_grad = g_se_middle[qm_se_indices.size() + i];
				qmmm_helpers::calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);
				new_grads[l.qm] += g_qm;
				new_grads[l.mm] += g_mm;

				if (Config::get().general.verbosity > 4)
				{
					std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
					std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
				}
			}
		}
		if (sec_middle.integrity() == false) integrity = false;

		if (Config::get().general.verbosity > 4)
		{
			std::cout << "SE energy of intermediate system: \n";
			sec_middle.e_head_tostream_short(std::cout);
			sec_middle.e_tostream_short(std::cout);
		}
	}
	catch (...)
	{
		std::cout << "SE programme (for intermediate system) failed. Treating structure as broken.\n";
		integrity = false;  // if SE programme fails: integrity is destroyed
	}

 // // ############### ONLY AMBER: PREPARATION OF CHARGES FOR INTERMEDIATE SYSTEM ################

	//// temporarily: only QM charges and those of link atoms in amber_charges
	//std::vector<double> old_amber_charges;
	//if (Config::get().general.input == config::input_types::AMBER || Config::get().general.chargefile)
	//{
	//	old_amber_charges = Config::get().coords.amber_charges;                       // save old amber_charges
	//	qmmm_helpers::select_from_ambercharges(qm_indices);                           // only QM charges in amber_charges
	//	for (auto i = 0u; i < link_atoms.size(); ++i)                                   // add charges of link atoms
	//	{
	//		double la_charge = qmc.energyinterface()->charges()[qm_indices.size() + i]; // get charge
	//		Config::set().coords.amber_charges.push_back(la_charge*18.2223);            // convert it to AMBER units and add it to vector
	//	}
	//}

  // ############### MM ENERGY AND GRADIENTS FOR MIDDLE SYSTEM ######################

	try {
		if (!if_gradient)
		{
			mm_energy_middle = mmc_middle.e();  // calculate energy of intermediate MM system
		}
		else  // gradient calculation
		{
			mm_energy_middle = mmc_middle.g();     // get energy and calculate gradients
			auto g_mm_middle = mmc_middle.g_xyz(); // get gradients
			for (auto&& qsi : qm_se_indices)
			{
				new_grads[qsi] -= g_mm_middle[new_indices_qm[qsi]];
			}

			for (auto i = 0u; i < link_atoms_middle.size(); ++i)  // take into account link atoms
			{
				LinkAtom l = link_atoms_middle[i];

				coords::r3 g_qm, g_mm;             // divide link atom gradient to QM and MM atom
				auto link_atom_grad = g_mm_middle[qm_se_indices.size() + i];
				qmmm_helpers::calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);
				new_grads[l.qm] -= g_qm;
				new_grads[l.mm] -= g_mm;
				if (Config::get().general.verbosity > 4)
				{
					std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
					std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
				}
			}
		}
    if (mmc_middle.integrity() == false) integrity = false;

		if (Config::get().general.verbosity > 4)
		{
			std::cout << "MM energy of intermediate system: \n";
			mmc_middle.e_head_tostream_short(std::cout);
			mmc_middle.e_tostream_short(std::cout);
		}
	}
	catch (...)
	{
		std::cout << "MM programme (for intermediate system) failed. Treating structure as broken.\n";
		integrity = false;  // if MM programme fails: integrity is destroyed
	}

  // ############### GRADIENTS ON MM ATOMS DUE TO COULOMB INTERACTION WITH MIDDLE REGION ###

  if (if_gradient && integrity == true)
  {
		auto sec_middle_g_ext_charges = sec_middle.energyinterface()->get_g_ext_chg();
    auto mmc_middle_g_ext_charges = mmc_middle.energyinterface()->get_g_ext_chg();

    for (auto i=0u; i<charge_indices.size(); ++i)
    {
      int mma = charge_indices[i];
      new_grads[mma] += sec_middle_g_ext_charges[i];
      new_grads[mma] -= mmc_middle_g_ext_charges[i];
    }
  }

	// ############### CREATE EXTERNAL CHARGES FOR SMALL SYSTEM ######################

	Config::set().energy.qmmm.mm_charges.clear();
	charge_vector = sec_middle.energyinterface()->charges();
	charge_indices.clear();
	for (auto i = qm_indices.size(); i < qm_se_indices.size(); ++i)  // go through all se atoms
	{
		use_charge = true;
		
		for (auto b : sec_middle.atoms().atom(i).bonds())
		{
			if (b < qm_indices.size())  // and ignore those which are bonded to a QM atom
			{
				use_charge = false;
				break;
			}
		}

		if (use_charge)  // for the other 
		{
			if (Config::get().energy.qmmm.cutoff != 0.0)  // if cutoff given: test if one QM atom is nearer than cutoff
			{
				use_charge = false;
				for (auto qm : qm_indices)
				{
					auto dist = len(sec_middle.xyz(i) - coords->xyz(qm));
					if (dist < Config::get().energy.qmmm.cutoff)
					{
						use_charge = true;
						break;
					}
				}

			}

			if (use_charge)  // if yes create a PointCharge and add it to vector
			{
				PointCharge new_charge;
				new_charge.charge = charge_vector[i];
				new_charge.set_xyz(sec_middle.xyz(i).x(), sec_middle.xyz(i).y(), sec_middle.xyz(i).z());
				Config::set().energy.qmmm.mm_charges.push_back(new_charge);

				charge_indices.push_back(i);  // add index to charge_indices
			}
		}
	}

	charge_vector = mmc_big.energyinterface()->charges();
	for (auto i = 0u; i < coords->size(); ++i)    // go through all atoms
	{
		use_charge = true;
		for (auto &l : link_atoms_small) // ignore those atoms that are connected to a QM atom...
		{
			if (l.mm == i) use_charge = false;
		}
		for (auto &qs : qm_se_indices) // ...and the QM and SE atoms
		{
			if (qs == i) use_charge = false;
		}

		if (use_charge)  // for the other 
		{
			if (Config::get().energy.qmmm.cutoff != 0.0)  // if cutoff given: test if one QM atom is nearer than cutoff
			{
				use_charge = false;
				for (auto qi : qm_indices)
				{
					auto dist = len(coords->xyz(i) - coords->xyz(qi));
					if (dist < Config::get().energy.qmmm.cutoff)
					{
						use_charge = true;
						break;
					}
				}

			}

			if (use_charge)  // if yes create a PointCharge and add it to vector
			{
				PointCharge new_charge;
				new_charge.charge = charge_vector[i];
				new_charge.set_xyz(coords->xyz(i).x(), coords->xyz(i).y(), coords->xyz(i).z());
				Config::set().energy.qmmm.mm_charges.push_back(new_charge);

				charge_indices.push_back(i);  // add index to charge_indices
			}
		}
	}

	// ############### QM ENERGY AND GRADIENTS FOR SMALL SYSTEM ######################
	try {
		if (!if_gradient)
		{
			qm_energy = qmc.e();  // get qm energy  
		}
		else  // gradient calculation
		{
			qm_energy = qmc.g();    // get energy and calculate gradients
			auto g_qm_small = qmc.g_xyz();        // get gradients
			for (auto&& qmi : qm_indices)
			{
				new_grads[qmi] += g_qm_small[new_indices_middle[qmi]];
			}

			for (auto i = 0u; i < link_atoms_small.size(); ++i)   // take into account link atoms
			{
				LinkAtom l = link_atoms_small[i];

				coords::r3 g_qm, g_mm;        // divide link atom gradient to QM and MM atom
				auto link_atom_grad = g_qm_small[qm_indices.size() + i];
				qmmm_helpers::calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);
				new_grads[l.qm] += g_qm;
				new_grads[l.mm] += g_mm;

				if (Config::get().general.verbosity > 4)
				{
					std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
					std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
				}
			}
		}
		if (qmc.integrity() == false) integrity = false;

		if (Config::get().general.verbosity > 4)
		{
			std::cout << "QM energy of small system: \n";
			qmc.e_head_tostream_short(std::cout);
			qmc.e_tostream_short(std::cout);
		}
	}
	catch (...)
	{
		std::cout << "QM programme (for small system) failed. Treating structure as broken.\n";
		integrity = false;  // if QM programme fails: integrity is destroyed
	}

	// // ############### ONLY AMBER: PREPARATION OF CHARGES FOR SMALL SYSTEM ################

	//// temporarily: only QM charges and those of link atoms in amber_charges
	//std::vector<double> old_amber_charges;
	//if (Config::get().general.input == config::input_types::AMBER || Config::get().general.chargefile)
	//{
	//	old_amber_charges = Config::get().coords.amber_charges;                       // save old amber_charges
	//	qmmm_helpers::select_from_ambercharges(qm_indices);                           // only QM charges in amber_charges
	//	for (auto i = 0u; i < link_atoms.size(); ++i)                                   // add charges of link atoms
	//	{
	//		double la_charge = qmc.energyinterface()->charges()[qm_indices.size() + i]; // get charge
	//		Config::set().coords.amber_charges.push_back(la_charge*18.2223);            // convert it to AMBER units and add it to vector
	//	}
	//}

	// ############### SE ENERGY AND GRADIENTS FOR SMALL SYSTEM ######################

	try {
		if (!if_gradient)
		{
			se_energy_small = sec_small.e();  // calculate energy of small SE system
		}
		else  // gradient calculation
		{
			se_energy_small = sec_small.g();     // get energy and calculate gradients
			auto g_se_small = sec_small.g_xyz(); // get gradients
			for (auto&& qmi : qm_indices)
			{
				new_grads[qmi] -= g_se_small[new_indices_qm[qmi]];
			}

			for (auto i = 0u; i < link_atoms_middle.size(); ++i)  // take into account link atoms
			{
				LinkAtom l = link_atoms_middle[i];

				coords::r3 g_qm, g_mm;             // divide link atom gradient to QM and MM atom
				auto link_atom_grad = g_se_small[qm_se_indices.size() + i];
				qmmm_helpers::calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);
				new_grads[l.qm] -= g_qm;
				new_grads[l.mm] -= g_mm;
				if (Config::get().general.verbosity > 4)
				{
					std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
					std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
				}
			}
		}
		if (sec_small.integrity() == false) integrity = false;

		if (Config::get().general.verbosity > 4)
		{
			std::cout << "SE energy of small system: \n";
			sec_small.e_head_tostream_short(std::cout);
			sec_small.e_tostream_short(std::cout);
		}
	}
	catch (...)
	{
		std::cout << "SE programme (for small system) failed. Treating structure as broken.\n";
		integrity = false;  // if SE programme fails: integrity is destroyed
	}

	// ############### GRADIENTS ON MM ATOMS DUE TO COULOMB INTERACTION WITH MIDDLE REGION ###

	if (if_gradient && integrity == true)
	{
		auto qmc_g_ext_charges = qmc.energyinterface()->get_g_ext_chg();
		auto sec_small_g_ext_charges = sec_small.energyinterface()->get_g_ext_chg();

		for (auto i = 0u; i<charge_indices.size(); ++i)
		{
			int mma = charge_indices[i];
			new_grads[mma] += qmc_g_ext_charges[i];
			new_grads[mma] -= sec_small_g_ext_charges[i];
		}
	}

  // ############### STUFF TO DO AT THE END OF CALCULATION ######################

  Config::set().energy.qmmm.mm_charges.clear();  // clear vector -> no point charges in calculation of mmc_big
	//Config::set().coords.amber_charges = old_amber_charges;  // set AMBER charges back to total AMBER charges

  if (check_bond_preservation() == false) integrity = false;
  else if (check_atom_dist() == false) integrity = false;
  
  if (if_gradient) coords->swap_g_xyz(new_grads);     // swap gradients into coordobj
 return mm_energy_big + se_energy_middle - mm_energy_middle + qm_energy - se_energy_small; // return total energy
}

coords::float_type energy::interfaces::three_layer::THREE_LAYER::g()
{
  integrity = true;
  energy = qmmm_calc(true);
  return energy;
}

coords::float_type energy::interfaces::three_layer::THREE_LAYER::e()
{
  integrity = true;
  energy = qmmm_calc(false);
  return energy;
}

coords::float_type energy::interfaces::three_layer::THREE_LAYER::h()
{
  throw std::runtime_error("no hessian function implemented for this interface");
}

coords::float_type energy::interfaces::three_layer::THREE_LAYER::o()
{
  throw std::runtime_error("this interface doesn't have an own optimizer");
}

void energy::interfaces::three_layer::THREE_LAYER::print_E(std::ostream &) const
{
  throw std::runtime_error("function not implemented");
}

void energy::interfaces::three_layer::THREE_LAYER::print_E_head(std::ostream &S, bool const endline) const
{
  S << "QM-atoms: " << qm_indices.size() << '\n';
  //S << "MM-atoms: " << mm_indices.size() << '\n';
  S << "Potentials\n";
  S << std::right << std::setw(24) << "MM_big";
  S << std::right << std::setw(24) << "MM_middle";
  S << std::right << std::setw(24) << "SE_middle";
  S << std::right << std::setw(24) << "SE_small";
  S << std::right << std::setw(24) << "QM";
  S << std::right << std::setw(24) << "TOTAL";
  if (endline) S << '\n';
}

void energy::interfaces::three_layer::THREE_LAYER::print_E_short(std::ostream &S, bool const endline) const
{
  S << '\n';
  S << std::right << std::setw(24) << mm_energy_big;
  S << std::right << std::setw(24) << mm_energy_middle;
  S << std::right << std::setw(24) << se_energy_middle;
  S << std::right << std::setw(24) << se_energy_small;
  S << std::right << std::setw(24) << qm_energy;
  S << std::right << std::setw(24) << energy;
  if (endline) S << '\n';
}

bool energy::interfaces::three_layer::THREE_LAYER::check_bond_preservation(void) const
{
  std::size_t const N(coords->size());
  for (std::size_t i(0U); i < N; ++i)
  { // cycle over all atoms i
    if (!coords->atoms(i).bonds().empty())
    {
      std::size_t const M(coords->atoms(i).bonds().size());
      for (std::size_t j(0U); j < M && coords->atoms(i).bonds(j) < i; ++j)
      { // cycle over all atoms bound to i
        double const L(geometric_length(coords->xyz(i) - coords->xyz(coords->atoms(i).bonds(j))));
        double const max = 1.2 * (coords->atoms(i).cov_radius() + coords->atoms(coords->atoms(i).bonds(j)).cov_radius());
        if (L > max) return false;
      }
    }
  }
  return true;
}

bool energy::interfaces::three_layer::THREE_LAYER::check_atom_dist(void) const
{
  std::size_t const N(coords->size());
  for (std::size_t i(0U); i < N; ++i)
  {
    for (std::size_t j(0U); j < i; j++)
    {
      if (dist(coords->xyz(i), coords->xyz(j)) < 0.3)
      {
        return false;
      }
    }
  }
  return true;
}


