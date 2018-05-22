#include <cstddef>

#include "energy_int_oniom.h"
#include "scon_utility.h"

energy::interfaces::oniom::ONIOM::ONIOM(coords::Coordinates * cp) :
	energy::interface_base(cp)
{
}

energy::interfaces::oniom::ONIOM::ONIOM(ONIOM const & rhs,
  coords::Coordinates *cobj) : interface_base(cobj)
{
  interface_base::operator=(rhs);
}

energy::interfaces::oniom::ONIOM::ONIOM(ONIOM&& rhs, coords::Coordinates *cobj)
  : interface_base(cobj)
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
  //interface_base::swap(rhs);
  //std::swap(cparams, rhs.cparams);
  //qm_indices.swap(rhs.qm_indices);
  //mm_indices.swap(rhs.mm_indices);
  //new_indices_mm.swap(rhs.new_indices_mm);
  //new_indices_qm.swap(rhs.new_indices_qm);
  //qmc.swap(rhs.qmc);
  //mmc.swap(rhs.mmc);
  //qm_charge_vector.swap(rhs.qm_charge_vector);
  //mm_charge_vector.swap(rhs.mm_charge_vector);
  //std::swap(vdw_energy, rhs.vdw_energy);
  //std::swap(qm_energy, rhs.qm_energy);
  //std::swap(mm_energy, rhs.mm_energy);
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

coords::float_type energy::interfaces::oniom::ONIOM::g()
{
  std::cout << "Oniom gradient calcultation\n";
  return 0.0;
}

coords::float_type energy::interfaces::oniom::ONIOM::e()
{
	std::cout << "Oniom energy calcultation\n";
	return 0.0;
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
  //S << "QM-atoms: " << qm_indices.size() << '\n';
  //S << "MM-atoms: " << mm_indices.size() << '\n';
  //S << "Potentials\n";
  //S << std::right << std::setw(24) << "QM";
  //S << std::right << std::setw(24) << "MM";
  //S << std::right << std::setw(24) << "VDW";
  //S << std::right << std::setw(24) << "BONDED";
  //S << std::right << std::setw(24) << "TOTAL";
  //if (endline) S << '\n';
}

void energy::interfaces::oniom::ONIOM::print_E_short(std::ostream &S, bool const endline) const
{
  //S << '\n';
  //S << std::right << std::setw(24) << qm_energy;
  //S << std::right << std::setw(24) << mm_energy;
  //S << std::right << std::setw(24) << vdw_energy;
  ////S << std::right << std::setw(24) << bonded_energy;
  //S << std::right << std::setw(24) << energy;
  //if (endline) S << '\n';
}

void energy::interfaces::oniom::ONIOM::to_stream(std::ostream &S) const
{
  //S << '\n';
  //interface_base::to_stream(S);
  //throw std::runtime_error("no QMMM-function yet");
}

