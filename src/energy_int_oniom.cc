#include <cstddef>

#include "energy_int_oniom.h"
#include "scon_utility.h"

::tinker::parameter::parameters energy::interfaces::oniom::ONIOM::tp;

energy::interfaces::oniom::ONIOM::ONIOM(coords::Coordinates *cp):
  energy::interface_base(cp), qm_indices(Config::get().energy.qmmm.qmatoms),
  mm_indices(qmmm_helpers::get_mm_atoms(cp->size())),
  new_indices_qm(qmmm_helpers::make_new_indices_qm(cp->size())),
  link_atoms(qmmm_helpers::create_link_atoms(cp, qm_indices, mm_indices, tp)),
  qmc(qmmm_helpers::make_small_coords(cp, qm_indices, new_indices_qm, link_atoms, Config::get().energy.qmmm.qminterface)),
  mmc_big(qmmm_helpers::make_mmbig_coords(cp)),
  mmc_small(qmmm_helpers::make_small_coords(coords, qm_indices, new_indices_qm, link_atoms, Config::get().energy.qmmm.mminterface)),
  qm_energy(0.0), mm_energy_small(0.0), mm_energy_big(0.0)
{
  if (!tp.valid())
  {
    tp.from_file(Config::get().get().general.paramFilename);
  }
}

energy::interfaces::oniom::ONIOM::ONIOM(ONIOM const & rhs,
  coords::Coordinates *cobj) : interface_base(cobj),
  qm_indices(rhs.qm_indices), mm_indices(rhs.mm_indices),
  new_indices_qm(rhs.new_indices_qm), link_atoms(rhs.link_atoms),
  qmc(rhs.qmc), mmc_big(rhs.mmc_big), mmc_small(rhs.mmc_small), 
  qm_energy(rhs.qm_energy), mm_energy_big(rhs.mm_energy_big), mm_energy_small(rhs.mm_energy_small)
{
  interface_base::operator=(rhs);
}

energy::interfaces::oniom::ONIOM::ONIOM(ONIOM&& rhs, coords::Coordinates *cobj)
  : interface_base(cobj),
  qm_indices(std::move(rhs.qm_indices)), mm_indices(std::move(rhs.mm_indices)),
  new_indices_qm(std::move(rhs.new_indices_qm)), link_atoms(std::move(rhs.link_atoms)),
  qmc(std::move(rhs.qmc)), mmc_big(std::move(rhs.mmc_big)), mmc_small(std::move(rhs.mmc_small)),
  qm_energy(std::move(rhs.qm_energy)), mm_energy_big(std::move(rhs.mm_energy_big)), mm_energy_small(std::move(rhs.mm_energy_small))
{
  interface_base::operator=(rhs);
}


energy::interface_base * energy::interfaces::oniom::ONIOM::clone(coords::Coordinates * c) const
{
  ONIOM * tmp = new ONIOM(*this, c);
  return tmp;
}

energy::interface_base * energy::interfaces::oniom::ONIOM::move(coords::Coordinates * c)
{
  ONIOM * tmp = new ONIOM(std::move(*this), c);
  return tmp;
}


void energy::interfaces::oniom::ONIOM::swap(interface_base& rhs)
{
  swap(dynamic_cast<ONIOM&>(rhs));
}

void energy::interfaces::oniom::ONIOM::swap(ONIOM& rhs)
{
  interface_base::swap(rhs);
  qm_indices.swap(rhs.qm_indices);
  mm_indices.swap(rhs.mm_indices);
  link_atoms.swap(rhs.link_atoms);
  new_indices_qm.swap(rhs.new_indices_qm);
  qmc.swap(rhs.qmc);
  mmc_big.swap(rhs.mmc_big);
  mmc_small.swap(rhs.mmc_small);
  std::swap(qm_energy, rhs.qm_energy);
  std::swap(mm_energy_big, rhs.mm_energy_big);
  std::swap(mm_energy_small, rhs.mm_energy_small);
}

// update structure (account for topology or rep change)
void energy::interfaces::oniom::ONIOM::update(bool const skip_topology)
{
  if (!skip_topology)
  {
    *this = ONIOM(this->coords);
  }
  else
  {
    update_representation();
    qmc.energy_update(true);
    mmc_small.energy_update(true);
    mmc_big.energy_update(true);
  }
}

void energy::interfaces::oniom::ONIOM::update_representation()
{
  std::size_t qi = 0u;    // update position of QM atoms in small systems
  for (auto i : qm_indices)
  {
    qmc.move_atom_to(qi, coords->xyz(i), true);
    mmc_small.move_atom_to(qi, coords->xyz(i), true);
    ++qi;
  }

  for (auto &l : link_atoms) l.calc_position(coords); // update positions of link atoms in small systems
  for (int i = 0; i < link_atoms.size(); ++i)
  {
	  int index = qm_indices.size() + i;
	  coords::cartesian_type &new_pos = link_atoms[i].position;
	  qmc.move_atom_to(index, new_pos, true);
	  mmc_small.move_atom_to(index, new_pos, true);
  }

  std::size_t mi = 0u;       // update positions of all atoms in big system
  for (std::size_t mi = 0u; mi<coords->size(); mi++)
  {
    mmc_big.move_atom_to(mi, coords->xyz()[mi], true);
  }

}

coords::float_type energy::interfaces::oniom::ONIOM::qmmm_calc(bool if_gradient)
{
	if (link_atoms.size() != Config::get().energy.qmmm.linkatom_types.size())  // test if correct number of link atom types is given
	{                                                                          // can't be done in constructor because interface is first constructed without atoms 
		std::cout << "Wrong number of link atom types given. You have " << link_atoms.size() << " in the following order:\n";
		for (auto &l : link_atoms)
		{
			std::cout << "MM atom: " << l.mm + 1 << ", QM atom: " << l.qm + 1 << "\n";
	  }
		throw std::runtime_error("wrong number of link atom types");
	}

  update_representation(); // update positions of QM and MM subsystems to those of coordinates object
  
  coords::Gradients_3D new_grads;  // save gradients in case of gradient calculation

  // ############### MM ENERGY AND GRADIENTS FOR WHOLE SYSTEM ######################
  if (!if_gradient)
  {
    mm_energy_big = mmc_big.e();   // calculate MM energy of whole system
  }
  else   // gradient calculation
  {
    mm_energy_big = mmc_big.g();
    new_grads = mmc_big.g_xyz();
  }
  if (Config::get().general.verbosity > 4)
  {
    std::cout << "Energy of big MM system: \n";
    mmc_big.e_head_tostream_short(std::cout);
    mmc_big.e_tostream_short(std::cout);
  }

  // ############### CREATE MM CHARGES ######################

  //std::vector<double> charge_vector = mmc_big.energyinterface()->charges();

  //bool use_charge;
  //for (int i{ 0u }; i < coords->size(); ++i) // go through all atoms
  //{
	 // use_charge = true;
	 // for (auto &l : link_atoms) // ignore those atoms that are connected to a QM atom...
	 // {
		//  if (l.mm == i) use_charge = false;
	 // }
	 // for (auto &qm : qm_indices) // ... and the QM atoms themselves
	 // {
		//  if (qm == i) use_charge = false;
	 // }
	 // if (use_charge)  // for the other create a PointCharge and add it to vector
	 // {
		//  PointCharge new_charge;
		//  new_charge.charge = charge_vector[i];
		//  new_charge.set_xyz(coords->xyz(i).x(), coords->xyz(i).y(), coords->xyz(i).z());
		//  Config::set().energy.qmmm.mm_charges.push_back(new_charge);
	 // }
  //}

  // ############### MM ENERGY AND GRADIENTS FOR SMALL MM SYSTEM ######################
  if (!if_gradient)
  {
    mm_energy_small = mmc_small.e();  // calculate energy of small MM system
  }
  else  // gradient calculation
  {
    mm_energy_small = mmc_small.g();     // get energy and calculate gradients
    auto g_mm_small = mmc_small.g_xyz(); // get gradients
    for (auto&& qmi : qm_indices)
    {
      new_grads[qmi] -= g_mm_small[new_indices_qm[qmi]];
    }

    for (int i = 0; i < link_atoms.size(); ++i)  // take into account link atoms
    {
      LinkAtom l = link_atoms[link_atoms.size()-1-i];

      coords::r3 g_qm, g_mm;             // divide link atom gradient to QM and MM atom
      auto link_atom_grad = g_mm_small[qm_indices.size() + link_atoms.size() - 1 - i];
      qmmm_helpers::calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);
      new_grads[l.qm] -= g_qm;
      new_grads[l.mm] -= g_mm;
      if (Config::get().general.verbosity > 4)
      {
        std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
        std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
      }

      g_mm_small.pop_back();  // delete LinkAtom from gradients
    }
  }
  if (Config::get().general.verbosity > 4)
  {
    std::cout << "Energy of small MM system: \n";
    mmc_small.e_head_tostream_short(std::cout);
    mmc_small.e_tostream_short(std::cout);
  }

  // ############### QM ENERGY AND GRADIENTS FOR QM SYSTEM ######################
  try {
    if (!if_gradient)
    {
      qm_energy = qmc.e();  // get energy for QM part 
    }
    else  // gradient calculation
    {
      qm_energy = qmc.g();            // get energy and calculate gradients
      auto g_qm_small = qmc.g_xyz();  // get gradients
      for (auto&& qmi : qm_indices)
      {
        new_grads[qmi] += g_qm_small[new_indices_qm[qmi]];
      }

      for (int i = 0; i < link_atoms.size(); ++i)   // take into account link atoms
      {
        LinkAtom l = link_atoms[link_atoms.size() - 1 - i];

        coords::r3 g_qm, g_mm;        // divide link atom gradient to QM and MM atom
        auto link_atom_grad = g_qm_small[qm_indices.size() + link_atoms.size() - 1 - i];
        qmmm_helpers::calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);
        new_grads[l.qm] += g_qm;
        new_grads[l.mm] += g_mm;

        if (Config::get().general.verbosity > 4)
        {
          std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
          std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
        }

        g_qm_small.pop_back();  // delete LinkAtom from gradients
      }
    }
    if (Config::get().general.verbosity > 4)
    {
      std::cout << "Energy of QM system: \n";
      qmc.e_head_tostream_short(std::cout);
      qmc.e_tostream_short(std::cout);
    }
  }
  catch (...)
  {
    std::cout << "QM programme failed. Treating structure as broken.\n";
    integrity = false;  // if QM programme fails: integrity is destroyed
  }

  // ############### GRADIENTS ON MM ATOMS DUE TO COULOMB INTERACTION WITH QM REGION ###

  //if (if_gradient)
  //{
  //  for (int i=0; i<mm_indices.size(); ++i)
  //  {
  //    int mma = mm_indices[i];
  //    new_grads[mma] += qmc.energyinterface()->get_g_coul_mm()[i];
  //    new_grads[mma] -= mmc_small.energyinterface()->get_g_coul_mm()[i];
  //  }
  //}

  // ############### STUFF TO DO AT THE END OF CALCULATION ######################

  Config::set().energy.qmmm.mm_charges.clear();  // clear vector -> no point charges in calculation of mmc_big

  if (check_bond_preservation() == false) integrity = false;
  else if (check_atom_dist() == false) integrity = false;
  
  if (if_gradient) coords->swap_g_xyz(new_grads);     // swap gradients into coordobj
  return mm_energy_big - mm_energy_small + qm_energy; // return total energy
}

coords::float_type energy::interfaces::oniom::ONIOM::g()
{
  integrity = true;
  energy = qmmm_calc(true);
  return energy;
}

coords::float_type energy::interfaces::oniom::ONIOM::e()
{
  integrity = true;
  energy = qmmm_calc(false);
  return energy;
}

coords::float_type energy::interfaces::oniom::ONIOM::h()
{
  throw std::runtime_error("no QMMM-function yet");
}

coords::float_type energy::interfaces::oniom::ONIOM::o()
{
  throw std::runtime_error("QMMM-cannot optimize");
}

void energy::interfaces::oniom::ONIOM::print_E(std::ostream &) const
{
  throw std::runtime_error("no QMMM-function yet");
}

void energy::interfaces::oniom::ONIOM::print_E_head(std::ostream &S, bool const endline) const
{
  S << "QM-atoms: " << qm_indices.size() << '\n';
  S << "MM-atoms: " << mm_indices.size() << '\n';
  S << "Potentials\n";
  S << std::right << std::setw(24) << "MM_big";
  S << std::right << std::setw(24) << "MM_small";
  S << std::right << std::setw(24) << "QM";
  S << std::right << std::setw(24) << "TOTAL";
  if (endline) S << '\n';
}

void energy::interfaces::oniom::ONIOM::print_E_short(std::ostream &S, bool const endline) const
{
  S << '\n';
  S << std::right << std::setw(24) << mm_energy_big;
  S << std::right << std::setw(24) << mm_energy_small;
  S << std::right << std::setw(24) << qm_energy;
  S << std::right << std::setw(24) << energy;
  if (endline) S << '\n';
}

void energy::interfaces::oniom::ONIOM::to_stream(std::ostream &S) const
{
  throw std::runtime_error("no QMMM-function yet");
}

bool energy::interfaces::oniom::ONIOM::check_bond_preservation(void) const
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

bool energy::interfaces::oniom::ONIOM::check_atom_dist(void) const
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


