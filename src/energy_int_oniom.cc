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
    //update_representation();
    /*qmc.energy_update(true);
    mmc.energy_update(true);*/
  }
}

coords::float_type energy::interfaces::oniom::ONIOM::qmmm_calc(bool if_gradient)
{
  mm_energy_big = mmc_big.e();   // calculate MM energy of whole system
  if (Config::get().general.verbosity > 3)
  {
    std::cout << "energy of big MM system: \n";
    mmc_big.e_head_tostream_short(std::cout);
    mmc_big.e_tostream_short(std::cout);
  }

  mm_energy_small = mmc_small.e();  // calculate energy of small MM system
  if (Config::get().general.verbosity > 3)
  {
    std::cout << "energy of small MM system: \n";
    mmc_small.e_head_tostream_short(std::cout);
    mmc_small.e_tostream_short(std::cout);
  }

  try {
    qm_energy = qmc.e();  // get energy for QM part 
    if (Config::get().general.verbosity > 3)
    {
      std::cout << "energy of QM system: \n";
      mmc_small.e_head_tostream_short(std::cout);
      mmc_small.e_tostream_short(std::cout);
    }
  }
  catch (...)
  {
    std::cout << "QM programme failed. Treating structure as broken.\n";
    integrity = false;  // if QM programme fails: integrity is destroyed
  }
  
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

