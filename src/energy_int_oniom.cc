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
  //cparams(rhs.cparams),
  qm_indices(rhs.qm_indices), mm_indices(rhs.mm_indices),
  new_indices_qm(rhs.new_indices_qm), //new_indices_mm(rhs.new_indices_mm),
  //qmc(rhs.qmc), mmc_big(rhs.mmc_big), mmc_small(rhs.mmc_small), //qm_charge_vector(rhs.qm_charge_vector),
  //mm_charge_vector(rhs.mm_charge_vector), vdw_energy(rhs.vdw_energy),
  qm_energy(rhs.qm_energy), mm_energy_big(rhs.mm_energy_big), mm_energy_small(rhs.mm_energy_small)
{
  interface_base::operator=(rhs);
}

energy::interfaces::oniom::ONIOM::ONIOM(ONIOM&& rhs, coords::Coordinates *cobj)
  : interface_base(cobj),
  //cparams(std::move(rhs.cparams)),
  qm_indices(std::move(rhs.qm_indices)), mm_indices(std::move(rhs.mm_indices)),
  new_indices_qm(std::move(rhs.new_indices_qm)),
  //new_indices_mm(std::move(rhs.new_indices_mm)),
  //qmc(std::move(rhs.qmc)), mmc_big(std::move(rhs.mmc_big)), mmc_small(std::move(rhs.mmc_small)),
  /*qm_charge_vector(std::move(rhs.qm_charge_vector)),
  mm_charge_vector(std::move(rhs.mm_charge_vector)),
  vdw_energy(std::move(rhs.vdw_energy)),*/
  qm_energy(std::move(rhs.qm_energy)), mm_energy_big(std::move(rhs.mm_energy_big)), mm_energy_small(std::move(rhs.mm_energy_small))
  /*c_gradient(std::move(rhs.c_gradient)), vdw_gradient(std::move(rhs.vdw_gradient)),
  bonded_energy(std::move(rhs.bonded_energy)), bonded_gradient(std::move(rhs.bonded_gradient))*/
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
  //std::swap(cparams, rhs.cparams);
  qm_indices.swap(rhs.qm_indices);
  mm_indices.swap(rhs.mm_indices);
  //new_indices_mm.swap(rhs.new_indices_mm);
  new_indices_qm.swap(rhs.new_indices_qm);
  /*qmc.swap(rhs.qmc);
  mmc_big.swap(rhs.mmc_big);
  mmc_small.swap(rhs.mmc_small);*/
  /*qm_charge_vector.swap(rhs.qm_charge_vector);
  mm_charge_vector.swap(rhs.mm_charge_vector);*/
  std::swap(qm_energy, rhs.qm_energy);
  std::swap(mm_energy_big, rhs.mm_energy_big);
  //c_gradient.swap(rhs.c_gradient);
  //vdw_gradient.swap(rhs.vdw_gradient);
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
  std::size_t qi = 0u;
  for (auto i : qm_indices)
  {
    qmc.move_atom_to(qi, coords->xyz(i), true);
    mmc_small.move_atom_to(qi, coords->xyz(i), true);
    ++qi;
  }
  std::size_t mi = 0u;
  for (std::size_t mi = 0u; mi<coords->size(); mi++)
  {
    mmc_big.move_atom_to(mi, coords->xyz()[mi], true);
  }

}

coords::float_type energy::interfaces::oniom::ONIOM::qmmm_calc(bool if_gradient)
{
  update_representation(); // update positions of QM and MM subsystem to those of coordinates object
  for (auto &l : link_atoms) l.calc_position(coords); // update positions of link atoms

  coords::Gradients_3D new_grads;  // save gradients in case of gradient calculation

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

  if (!if_gradient)
  {
    mm_energy_small = mmc_small.e();  // calculate energy of small MM system
  }
  else  // gradient calculation
  {
    mm_energy_small = mmc_small.g();
    auto g_mm_small = mmc_small.g_xyz();
    for (int i = 0; i < link_atoms.size(); ++i)
    {
      auto link_atom_grad = g_mm_small[qm_indices.size() + link_atoms.size() -1 -i];
      LinkAtom l = link_atoms[i];

      coords::r3 g_qm, g_mm;
      qmmm_helpers::calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);

      if (Config::get().general.verbosity > 4)
      {
        std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
        std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
      }

      new_grads[l.qm] -= g_qm;
      new_grads[l.mm] -= g_mm;

      g_mm_small.pop_back();  // delete LinkAtom from gradients
    }
    for (auto&& qmi : qm_indices)
    {
      new_grads[qmi] -= g_mm_small[new_indices_qm[qmi]];
    }
  }
  if (Config::get().general.verbosity > 4)
  {
    std::cout << "Energy of small MM system: \n";
    mmc_small.e_head_tostream_short(std::cout);
    mmc_small.e_tostream_short(std::cout);
  }

  try {
    if (!if_gradient)
    {
      qm_energy = qmc.e();  // get energy for QM part 
    }
    else  // gradient calculation
    {
      qm_energy = qmc.g();
      auto g_qm_small = qmc.g_xyz();
      for (int i = 0; i < link_atoms.size(); ++i)
      {
        auto link_atom_grad = g_qm_small[qm_indices.size() + link_atoms.size() - 1 - i];
        LinkAtom l = link_atoms[i];

        coords::r3 g_qm, g_mm;
        qmmm_helpers::calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);

        if (Config::get().general.verbosity > 4)
        {
          std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
          std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
        }

        new_grads[l.qm] += g_qm;
        new_grads[l.mm] += g_mm;

        g_qm_small.pop_back();
      }
      for (auto&& qmi : qm_indices)
      {
        new_grads[qmi] += g_qm_small[new_indices_qm[qmi]];
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


