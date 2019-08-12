#include <cstddef>

#include "energy_int_3layer.h"
#include "Scon/scon_utility.h"


::tinker::parameter::parameters energy::interfaces::three_layer::THREE_LAYER::tp;

energy::interfaces::three_layer::THREE_LAYER::THREE_LAYER(coords::Coordinates* cp) :
	energy::interface_base(cp), qm_indices(Config::get().energy.qmmm.qm_systems[0]),
	qm_se_indices(add_vectors(qm_indices, Config::get().energy.qmmm.seatoms, true)),
	new_indices_qm(qmmm_helpers::make_new_indices(qm_indices, cp->size())),
	new_indices_middle(qmmm_helpers::make_new_indices(qm_se_indices, cp->size())),
	link_atoms_small(qmmm_helpers::create_link_atoms(qm_indices, cp, tp, std::vector<int>())),
	link_atoms_middle(qmmm_helpers::create_link_atoms(qm_se_indices, cp, tp, Config::get().energy.qmmm.linkatom_sets[0])),
	qmc(qmmm_helpers::make_small_coords(cp, qm_indices, new_indices_qm, Config::get().energy.qmmm.qminterface, "Small system: ",
		Config::get().energy.qmmm.qm_to_file, link_atoms_small, "small_system.arc")),
	sec_small(qmmm_helpers::make_small_coords(cp, qm_indices, new_indices_qm, Config::get().energy.qmmm.seinterface, "Small system: ",
		Config::get().energy.qmmm.qm_to_file, link_atoms_small, "small_system.arc")),
	sec_middle(qmmm_helpers::make_small_coords(cp, qm_se_indices, new_indices_middle, Config::get().energy.qmmm.seinterface, "Intermediate system: ",
		Config::get().energy.qmmm.qm_to_file, link_atoms_middle, "intermediate_system.arc")),
	mmc_middle(qmmm_helpers::make_small_coords(cp, qm_se_indices, new_indices_middle, Config::get().energy.qmmm.mminterface, "Intermediate system: ",
		Config::get().energy.qmmm.qm_to_file, link_atoms_middle, "intermediate_system.arc")),
	mmc_big(qmmm_helpers::make_small_coords(cp, range(cp->size()), range(cp->size()), Config::get().energy.qmmm.mminterface, "Big system: ")),
	index_of_middle_center(qmmm_helpers::get_index_of_QM_center(Config::get().energy.qmmm.centers[0], qm_se_indices, coords)),
	index_of_small_center(qmmm_helpers::get_index_of_QM_center(Config::get().energy.qmmm.small_center, qm_indices, coords)),
	qm_energy(0.0), se_energy_small(0.0), se_energy_middle(0.0), mm_energy_middle(0.0), mm_energy_big(0.0)
{
	sec_small.energyinterface()->charge = qmc.energyinterface()->charge;          // set correct charges for small and intermediate system
	mmc_middle.energyinterface()->charge = sec_middle.energyinterface()->charge;

	if (Config::get().energy.qmmm.qminterface == config::interface_types::T::MOPAC) Config::set().energy.mopac.link_atoms = link_atoms_small.size();   // set number of link atoms for MOPAC
	if (Config::get().energy.qmmm.seinterface == config::interface_types::T::MOPAC) Config::set().energy.mopac.link_atoms = link_atoms_middle.size();

	if ((Config::get().energy.qmmm.qminterface != config::interface_types::T::DFTB && Config::get().energy.qmmm.qminterface != config::interface_types::T::GAUSSIAN
		&& Config::get().energy.qmmm.qminterface != config::interface_types::T::PSI4 && Config::get().energy.qmmm.qminterface != config::interface_types::T::MOPAC
		&& Config::get().energy.qmmm.qminterface != config::interface_types::T::ORCA)
		||
		(Config::get().energy.qmmm.seinterface != config::interface_types::T::DFTB && Config::get().energy.qmmm.seinterface != config::interface_types::T::GAUSSIAN
			&& Config::get().energy.qmmm.seinterface != config::interface_types::T::PSI4 && Config::get().energy.qmmm.seinterface != config::interface_types::T::MOPAC
			&& Config::get().energy.qmmm.seinterface != config::interface_types::T::ORCA)
		||
		(Config::get().energy.qmmm.mminterface != config::interface_types::T::OPLSAA && Config::get().energy.qmmm.mminterface != config::interface_types::T::AMBER &&
			Config::get().energy.qmmm.mminterface != config::interface_types::T::DFTB && Config::get().energy.qmmm.mminterface != config::interface_types::T::GAUSSIAN
			&& Config::get().energy.qmmm.mminterface != config::interface_types::T::PSI4 && Config::get().energy.qmmm.mminterface != config::interface_types::T::MOPAC
			&& Config::get().energy.qmmm.mminterface != config::interface_types::T::ORCA))
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
}

energy::interfaces::three_layer::THREE_LAYER::THREE_LAYER(THREE_LAYER const& rhs,
	coords::Coordinates* cobj) : interface_base(cobj),
	qm_indices(rhs.qm_indices), qm_se_indices(rhs.qm_se_indices),
	new_indices_qm(rhs.new_indices_qm), new_indices_middle(rhs.new_indices_middle), link_atoms_small(rhs.link_atoms_small),
	link_atoms_middle(rhs.link_atoms_middle), qmc(rhs.qmc), sec_small(rhs.sec_small),
	sec_middle(rhs.sec_middle), mmc_middle(rhs.mmc_middle), mmc_big(rhs.mmc_big),
	qm_energy(rhs.qm_energy), se_energy_small(rhs.se_energy_small), se_energy_middle(rhs.se_energy_middle),
	mm_energy_middle(rhs.mm_energy_middle), mm_energy_big(rhs.mm_energy_big)
{
	interface_base::operator=(rhs);
}

energy::interfaces::three_layer::THREE_LAYER::THREE_LAYER(THREE_LAYER&& rhs, coords::Coordinates* cobj)
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


energy::interface_base* energy::interfaces::three_layer::THREE_LAYER::clone(coords::Coordinates* c) const
{
	THREE_LAYER* tmp = new THREE_LAYER(*this, c);
	return tmp;
}

energy::interface_base* energy::interfaces::three_layer::THREE_LAYER::move(coords::Coordinates* c)
{
	THREE_LAYER* tmp = new THREE_LAYER(std::move(*this), c);
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
	for (auto& l : link_atoms_small) l.calc_position(coords); // update positions of link atoms in small systems
	for (auto i = 0u; i < link_atoms_small.size(); ++i)
	{
		int index = qm_indices.size() + i;
		coords::cartesian_type& new_pos = link_atoms_small[i].position;
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
	for (auto& l : link_atoms_middle) l.calc_position(coords); // update positions of link atoms in middle systems
	for (auto i = 0u; i < link_atoms_middle.size(); ++i)
	{
		int index = qm_se_indices.size() + i;
		coords::cartesian_type& new_pos = link_atoms_middle[i].position;
		sec_middle.move_atom_to(index, new_pos, true);
		mmc_middle.move_atom_to(index, new_pos, true);
	}

	for (std::size_t mi = 0u; mi < coords->size(); mi++)   // update positions of all atoms in big system
	{
		mmc_big.move_atom_to(mi, coords->xyz()[mi], true);
	}
}

coords::float_type energy::interfaces::three_layer::THREE_LAYER::qmmm_calc(bool if_gradient)
{
	if (link_atoms_middle.size() != Config::get().energy.qmmm.linkatom_sets[0].size())  // test if correct number of link atom types is given
	{                                                                                 // can't be done in constructor because interface is first constructed without atoms 
		std::cout << "Wrong number of link atom types given. You have " << link_atoms_middle.size() << " in the following order:\n";
		for (auto& l : link_atoms_middle)
		{
			std::cout << "QM atom: " << l.qm + 1 << ", MM atom: " << l.mm + 1 << "\n";
		}
		std::cout << "This is assuming you are using a forcefield for your big system. \nIf you want to use one for the intermediate system talk to a CAST developer!\n";
		throw std::runtime_error("wrong number of link atom types");
	}

	// test if no atom is double in intermediate system
	if (double_element(qm_se_indices) == true) throw std::runtime_error("ERROR! You have at least one atom in QM as well as in SE atoms.");

	update_representation(); // update positions of QM and MM subsystems to those of coordinates object

	mm_energy_big = 0.0;     // set energies to zero
	mm_energy_middle = 0.0;
	se_energy_middle = 0.0;
	se_energy_small = 0.0;
	qm_energy = 0.0;
	coords::Gradients_3D new_grads;  // save gradients in case of gradient calculation
	bool periodic = Config::get().periodics.periodic;

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
		if (mm_energy_big == 0) integrity = false;
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

	std::vector<int> charge_indices;  // indizes of all atoms that are in charge_vector
	charge_indices.clear();

	if (Config::get().energy.qmmm.zerocharge_bonds != 0)
	{
		auto mmc_big_charges = mmc_big.energyinterface()->charges();
		if (mmc_big_charges.size() == 0) throw std::runtime_error("no charges found in MM interface");
		auto all_indices = range(coords->size());
		qmmm_helpers::add_external_charges(qm_se_indices, mmc_big_charges, all_indices, link_atoms_middle, charge_indices, coords, index_of_middle_center);
	}

	Config::set().periodics.periodic = false;

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
		if (se_energy_middle == 0) integrity = false;

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

	// ############### ONLY AMBER: PREPARATION OF CHARGES FOR INTERMEDIATE SYSTEM ################

	// temporarily: only QM charges, SE charges and those of link atoms in amber_charges
	std::vector<double> old_amber_charges;
	if (Config::get().general.input == config::input_types::AMBER || Config::get().general.chargefile)
	{
		old_amber_charges = Config::get().coords.amber_charges;                       // save old amber_charges
		qmmm_helpers::select_from_ambercharges(qm_se_indices);                        // only QM and SE charges in amber_charges
		for (auto i = 0u; i < link_atoms_middle.size(); ++i)                                    // add charges of link atoms
		{
			double la_charge = sec_middle.energyinterface()->charges()[qm_se_indices.size() + i]; // get charge
			Config::set().coords.amber_charges.push_back(la_charge * 18.2223);                      // convert it to AMBER units and add it to vector
		}
	}

	// ################ SAVE OUTPUT FOR BIG MM SYSTEM ########################################################

	qmmm_helpers::save_outputfiles(Config::get().energy.qmmm.mminterface, mmc_big.energyinterface()->id, "big");

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
				new_grads[qsi] -= g_mm_middle[new_indices_middle[qsi]];
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
		if (mm_energy_middle == 0) integrity = false;

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

	if (if_gradient && integrity == true && Config::get().energy.qmmm.zerocharge_bonds != 0)
	{
		auto sec_middle_g_ext_charges = sec_middle.energyinterface()->get_g_ext_chg();
		auto mmc_middle_g_ext_charges = mmc_middle.energyinterface()->get_g_ext_chg();

		for (auto i = 0u; i < charge_indices.size(); ++i)
		{
			int mma = charge_indices[i];
			auto grad_sec = sec_middle_g_ext_charges[i];
			auto grad_mmc = mmc_middle_g_ext_charges[i];

			coords::r3 derivQ_sec{ 0.0, 0.0, 0.0 }, derivQ_mmc{ 0.0, 0.0, 0.0 };   // additional gradients because charge also changes with position
			if (Config::get().energy.qmmm.cutoff != std::numeric_limits<double>::max())
			{
				double constexpr elec_factor = 332.06;
				double const& c = Config::get().energy.qmmm.cutoff;
				double const& ext_chg = Config::get().energy.qmmm.mm_charges[i].original_charge;
				double const& scaling = Config::get().energy.qmmm.mm_charges[i].scaled_charge / Config::get().energy.qmmm.mm_charges[i].original_charge;
				auto chargesSE = sec_middle.energyinterface()->charges();
				auto chargesMM = mmc_middle.energyinterface()->charges();

				// calculate sum(Q_qm * Q_ext / r) 
				double sum_of_QM_interactions_sec{ 0.0 };
				double sum_of_QM_interactions_mmc{ 0.0 };
				for (auto j{ 0u }; j < qm_se_indices.size(); ++j)
				{
					double const QMcharge_sec = chargesSE[j];
					double const QMcharge_mmc = chargesMM[j];
					coords::r3 MMpos{ Config::get().energy.qmmm.mm_charges[i].x,  Config::get().energy.qmmm.mm_charges[i].y,  Config::get().energy.qmmm.mm_charges[i].z };
					double const dist = len(MMpos - coords->xyz(qm_se_indices[j]));
					sum_of_QM_interactions_sec += (QMcharge_sec * ext_chg * elec_factor) / dist;
					sum_of_QM_interactions_mmc += (QMcharge_mmc * ext_chg * elec_factor) / dist;
				}

				// additional gradient on external charge due to interaction with SEC_MIDDLE
				derivQ_sec.x() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_middle_center).x() - Config::get().energy.qmmm.mm_charges[i].x) * std::sqrt(scaling) / (c * c);
				derivQ_sec.y() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_middle_center).y() - Config::get().energy.qmmm.mm_charges[i].y) * std::sqrt(scaling) / (c * c);
				derivQ_sec.z() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_middle_center).z() - Config::get().energy.qmmm.mm_charges[i].z) * std::sqrt(scaling) / (c * c);
				grad_sec += derivQ_sec;

				// additional gradient on external charge due to interaction with MMC_MIDDLE
				derivQ_mmc.x() = sum_of_QM_interactions_mmc * 4 * (coords->xyz(index_of_middle_center).x() - Config::get().energy.qmmm.mm_charges[i].x) * std::sqrt(scaling) / (c * c);
				derivQ_mmc.y() = sum_of_QM_interactions_mmc * 4 * (coords->xyz(index_of_middle_center).y() - Config::get().energy.qmmm.mm_charges[i].y) * std::sqrt(scaling) / (c * c);
				derivQ_mmc.z() = sum_of_QM_interactions_mmc * 4 * (coords->xyz(index_of_middle_center).z() - Config::get().energy.qmmm.mm_charges[i].z) * std::sqrt(scaling) / (c * c);
				grad_mmc += derivQ_mmc;

				// additional gradient on QM atom that defines distance
				new_grads[index_of_middle_center] += (derivQ_mmc - derivQ_sec);
			}

			new_grads[mma] += grad_sec;
			new_grads[mma] -= grad_mmc;
		}
	}

	// ############### EXTERNAL CHARGES FOR SMALL SYSTEM ######################

	Config::set().coords.amber_charges = old_amber_charges;  // set AMBER charges back to total AMBER charges
	Config::set().periodics.periodic = periodic;

	if (Config::get().energy.qmmm.zerocharge_bonds != 0)
	{
		if (Config::get().energy.qmmm.emb_small == 0)   // if EEx: no external charges for small system
		{
			Config::set().energy.qmmm.mm_charges.clear();
		}

		else if (Config::get().energy.qmmm.emb_small > 1)  // EE+ and MMSE
		{
			Config::set().energy.qmmm.mm_charges.clear();
			charge_indices.clear();

			if (Config::get().energy.qmmm.emb_small == 3)  // only for MMSE
			{
				auto sec_middle_charges = sec_middle.energyinterface()->charges();
				if (sec_middle.size() == 0) throw std::runtime_error("no charges found in SE interface");
				qmmm_helpers::add_external_charges(qm_indices, sec_middle_charges, qm_se_indices, link_atoms_small, charge_indices, coords, index_of_small_center);   // add charges from SE atoms
			}

			auto mmc_big_charges = mmc_big.energyinterface()->charges();
			if (mmc_big_charges.size() == 0) throw std::runtime_error("no charges found in MM interface");
			auto all_indices = range(coords->size());
			qmmm_helpers::add_external_charges(qm_se_indices, mmc_big_charges, all_indices, link_atoms_small, charge_indices, coords, index_of_small_center);     // add charges from MM atoms
		}
	}

	Config::set().periodics.periodic = false;

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
				new_grads[qmi] += g_qm_small[new_indices_qm[qmi]];
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
		if (qm_energy == 0) integrity = false;

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

	// ################ SAVE OUTPUT FOR INTERMEDIATE SE SYSTEM ########################################################

	qmmm_helpers::save_outputfiles(Config::get().energy.qmmm.seinterface, sec_middle.energyinterface()->id, "intermediate");

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

			for (auto i = 0u; i < link_atoms_small.size(); ++i)  // take into account link atoms
			{
				LinkAtom l = link_atoms_small[i];

				coords::r3 g_qm, g_mm;             // divide link atom gradient to QM and MM atom
				auto link_atom_grad = g_se_small[qm_indices.size() + i];
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
		if (se_energy_small == 0) integrity = false;

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

	// ############### GRADIENTS ON MM ATOMS DUE TO COULOMB INTERACTION WITH SMALL REGION ###

	if (Config::get().energy.qmmm.emb_small != 0 && if_gradient && integrity == true && Config::get().energy.qmmm.zerocharge_bonds != 0)
	{
		auto qmc_g_ext_charges = qmc.energyinterface()->get_g_ext_chg();
		auto sec_small_g_ext_charges = sec_small.energyinterface()->get_g_ext_chg();

		for (auto i = 0u; i < charge_indices.size(); ++i)
		{
			int mma = charge_indices[i];
			auto grad_qmc = qmc_g_ext_charges[i];
			auto grad_sec = sec_small_g_ext_charges[i];

			coords::r3 derivQ_qmc{ 0.0, 0.0, 0.0 }, derivQ_sec{ 0.0, 0.0, 0.0 };   // additional gradients because charge also changes with position
			if (Config::get().energy.qmmm.cutoff != std::numeric_limits<double>::max())
			{
				double constexpr elec_factor = 332.06;
				double const& c = Config::get().energy.qmmm.cutoff;
				double const& ext_chg = Config::get().energy.qmmm.mm_charges[i].original_charge;
				double const& scaling = Config::get().energy.qmmm.mm_charges[i].scaled_charge / Config::get().energy.qmmm.mm_charges[i].original_charge;
				auto chargesQM = qmc.energyinterface()->charges();
				auto chargesSE = sec_small.energyinterface()->charges();

				// calculate sum(Q_qm * Q_ext / r) 
				double sum_of_QM_interactions_qmc{ 0.0 };
				double sum_of_QM_interactions_sec{ 0.0 };
				for (auto j{ 0u }; j < qm_indices.size(); ++j)
				{
					double const QMcharge_qmc = chargesQM[j];
					double const QMcharge_sec = chargesSE[j];
					coords::r3 MMpos{ Config::get().energy.qmmm.mm_charges[i].x,  Config::get().energy.qmmm.mm_charges[i].y,  Config::get().energy.qmmm.mm_charges[i].z };
					double const dist = len(MMpos - coords->xyz(qm_indices[j]));
					sum_of_QM_interactions_qmc += (QMcharge_qmc * ext_chg * elec_factor) / dist;
					sum_of_QM_interactions_sec += (QMcharge_sec * ext_chg * elec_factor) / dist;
				}

				// additional gradient on external charge due to interaction with QMC
				derivQ_qmc.x() = sum_of_QM_interactions_qmc * 4 * (coords->xyz(index_of_small_center).x() - Config::get().energy.qmmm.mm_charges[i].x) * std::sqrt(scaling) / (c * c);
				derivQ_qmc.y() = sum_of_QM_interactions_qmc * 4 * (coords->xyz(index_of_small_center).y() - Config::get().energy.qmmm.mm_charges[i].y) * std::sqrt(scaling) / (c * c);
				derivQ_qmc.z() = sum_of_QM_interactions_qmc * 4 * (coords->xyz(index_of_small_center).z() - Config::get().energy.qmmm.mm_charges[i].z) * std::sqrt(scaling) / (c * c);
				grad_qmc += derivQ_qmc;

				// additional gradient on external charge due to interaction with SEC_SMALL
				derivQ_sec.x() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_small_center).x() - Config::get().energy.qmmm.mm_charges[i].x) * std::sqrt(scaling) / (c * c);
				derivQ_sec.y() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_small_center).y() - Config::get().energy.qmmm.mm_charges[i].y) * std::sqrt(scaling) / (c * c);
				derivQ_sec.z() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_small_center).z() - Config::get().energy.qmmm.mm_charges[i].z) * std::sqrt(scaling) / (c * c);
				grad_sec += derivQ_sec;

				// additional gradient on QM atom that defines distance
				new_grads[index_of_small_center] += (derivQ_sec - derivQ_qmc);
			}

			new_grads[mma] += grad_qmc;
			new_grads[mma] -= grad_sec;
		}
	}

	// ############### STUFF TO DO AT THE END OF CALCULATION ######################

	Config::set().energy.qmmm.mm_charges.clear();          // clear vector -> no point charges in calculation of mmc_big
	Config::set().periodics.periodic = periodic;
	if (file_exists("orca.gbw")) std::remove("orca.gbw");  // delete orca MOs for small system, otherwise orca will try to use them for middle system and fail

	if (coords->check_bond_preservation() == false) integrity = false;
	else if (coords->check_for_crashes() == false) integrity = false;

	if (if_gradient) coords->swap_g_xyz(new_grads);     // swap gradients into coordobj
	energy = mm_energy_big + se_energy_middle - mm_energy_middle + qm_energy - se_energy_small;
	return energy; // return total energy
}

coords::float_type energy::interfaces::three_layer::THREE_LAYER::g()
{
	integrity = coords->check_structure();
	if (integrity == true) return qmmm_calc(true);
	else return 0;
}

coords::float_type energy::interfaces::three_layer::THREE_LAYER::e()
{
	integrity = coords->check_structure();
	if (integrity == true) return qmmm_calc(false);
	else return 0;
}

coords::float_type energy::interfaces::three_layer::THREE_LAYER::h()
{
	throw std::runtime_error("no hessian function implemented for this interface");
}

coords::float_type energy::interfaces::three_layer::THREE_LAYER::o()
{
	throw std::runtime_error("this interface doesn't have an own optimizer");
}

void energy::interfaces::three_layer::THREE_LAYER::print_E(std::ostream&) const
{
	throw std::runtime_error("function not implemented");
}

void energy::interfaces::three_layer::THREE_LAYER::print_E_head(std::ostream& S, bool const endline) const
{
	S << "QM-atoms: " << qm_indices.size() << '\n';
	S << "SE-atoms: " << qm_se_indices.size() - qm_indices.size() << '\n';
	S << "MM-atoms: " << coords->size() - qm_se_indices.size() << '\n';
	S << "Potentials\n";
	S << std::right << std::setw(24) << "MM_big";
	S << std::right << std::setw(24) << "MM_middle";
	S << std::right << std::setw(24) << "SE_middle";
	S << std::right << std::setw(24) << "SE_small";
	S << std::right << std::setw(24) << "QM";
	S << std::right << std::setw(24) << "TOTAL";
	if (endline) S << '\n';
}

void energy::interfaces::three_layer::THREE_LAYER::print_E_short(std::ostream& S, bool const endline) const
{
	S << '\n';
	S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << mm_energy_big;
	S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << mm_energy_middle;
	S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << se_energy_middle;
	S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << se_energy_small;
	S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << qm_energy;
	S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << energy;
	if (endline) S << '\n';
}


