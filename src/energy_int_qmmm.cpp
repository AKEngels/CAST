#include <cstddef>

#include "energy_int_qmmm.h"
#include "scon_utility.h"


namespace
{
  std::vector<std::size_t> get_mm_atoms(std::size_t const num_atoms)
  {
    std::vector<std::size_t> mm_atoms;
    auto qm_size = Config::get().energy.qmmm.qmatoms.size();
    if (num_atoms > qm_size)
    {
      auto mm_size = num_atoms - qm_size;
      mm_atoms.reserve(mm_size);
      for (unsigned i = 1; i <= num_atoms; i++)
      {
        auto i2 = scon::sorted::exists(Config::get().energy.qmmm.qmatoms, i);
        if (!i2) { mm_atoms.emplace_back(i); }
      }
    }
    return mm_atoms;
  }

  std::vector<std::size_t> make_new_indices_qm(std::size_t const num_atoms)
  {
    std::vector<std::size_t> new_indices_qm;
    if (num_atoms > 0)
    {
      new_indices_qm.resize(num_atoms);
      std::size_t current_index = 0u;
      for (auto && a : Config::get().energy.qmmm.qmatoms)
      {
        new_indices_qm.at(a - 1) = current_index++;
      }
    }
    return new_indices_qm;
  }

  std::vector<std::size_t> make_new_indices_mm(std::size_t const num_atoms)
  {
    std::vector<std::size_t> new_indices_mm;
    new_indices_mm.resize(num_atoms);
    if (num_atoms > 0u)
    {
      std::size_t current_index = 0u;

      for (auto && a : get_mm_atoms(num_atoms))
      {
        new_indices_mm.at(a - 1) = current_index++;
      }
    }
    return new_indices_mm;
  }

  coords::Coordinates make_mopac_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices)
  {
    auto tmp_i = Config::get().general.energy_interface;
    Config::set().general.energy_interface = config::interface_types::T::MOPAC;
    std::cout << "Want MOPAC!\n";
    coords::Coordinates new_mopac_coords;
    if (cp->size() >= indices.size())
    {
      coords::Atoms new_mopac_atoms;
      coords::PES_Point pes;
      pes.structure.cartesian.reserve(indices.size());
      for (auto && a : indices)
      {
        auto && ref_at = (*cp).atoms().atom(a - 1);
        coords::Atom at{ (*cp).atoms().atom(a - 1).number() };
        at.set_energy_type(ref_at.energy_type());
        auto bonds = ref_at.bonds();
        for (auto && b : bonds)
        {
          at.detach_from(b);
        }
        for (auto && b : bonds)
        {
          at.bind_to(new_indices.at(b));
        }
        new_mopac_atoms.add(at);
        pes.structure.cartesian.push_back(cp->xyz(a - 1));
      }
      new_mopac_coords.init_swap_in(new_mopac_atoms, pes);
    }
    Config::set().general.energy_interface = tmp_i;
    return new_mopac_coords;
  }

  coords::Coordinates make_aco_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices)
  {
    auto tmp_i = Config::get().general.energy_interface;
    Config::set().general.energy_interface = Config::get().energy.qmmm.mminterface;
    coords::Coordinates new_aco_coords;
    if (cp->size() >= indices.size())
    {
      coords::Atoms new_aco_atoms;
      coords::PES_Point pes;
      pes.structure.cartesian.reserve(indices.size());
      for (auto && a : indices)
      {
        auto && ref_at = (*cp).atoms().atom(a - 1);
        coords::Atom at{ (*cp).atoms().atom(a - 1).number() };
        at.set_energy_type(ref_at.energy_type());
        auto bonds = ref_at.bonds();
        for (auto && b : bonds)
        {
          at.detach_from(b);
        }
        for (auto && b : bonds)
        {
          at.bind_to(new_indices[b]);
        }
        new_aco_atoms.add(at);
        pes.structure.cartesian.push_back(cp->xyz(a - 1));
      }
      new_aco_coords.init_swap_in(new_aco_atoms, pes);
    }
    Config::set().general.energy_interface = tmp_i;
    return new_aco_coords;
  }
}


energy::interfaces::qmmm::QMMM::QMMM(QMMM const & rhs,
  coords::Coordinates *cobj) : interface_base(cobj), 
  qm_indices(rhs.qm_indices), mm_indices(rhs.mm_indices),
  qmc(rhs.qmc), mmc(rhs.mmc), cparams(rhs.cparams)
{
  interface_base::operator=(rhs);
}

::tinker::parameter::parameters energy::interfaces::qmmm::QMMM::tp;

energy::interfaces::qmmm::QMMM::QMMM(coords::Coordinates * cp) :
  interface_base(cp),
  qm_indices(Config::get().energy.qmmm.qmatoms), 
  mm_indices(get_mm_atoms(cp->size())),
  new_indices_qm(make_new_indices_qm(cp->size())),
  new_indices_mm(make_new_indices_mm(cp->size())),
  qmc(make_mopac_coords(cp, qm_indices, new_indices_qm)),
  mmc(make_aco_coords(cp, mm_indices, new_indices_mm))
{
  if (!tp.valid())
  {
    tp.from_file(Config::get().get().general.paramFilename);
  }
  std::vector<std::size_t> types;
  for (auto atom : (*cp).atoms())
  {
    scon::sorted::insert_unique(types, atom.energy_type());
  }
  for (auto mma : mm_indices)
  {
    for (auto b : cp->atoms().atom(mma).bonds())
    {
      if (scon::sorted::exists(qm_indices, b))
      {
        throw std::runtime_error("Cannot break bonds.\n");
      }
    }
  }
  cparams = tp.contract(types);
  /*if (cp->size() != 0)
  {
    charges();
  }*/
 
  /*qmc.g();
  qmc.energyinterface()->print_E_head(std::cout);
  qmc.energyinterface()->print_E_short(std::cout);
    
  mmc.g();
  mmc.energyinterface()->print_E_head(std::cout);
  mmc.energyinterface()->print_E_short(std::cout);*/
}

coords::float_type energy::interfaces::qmmm::QMMM::qmmm_calc(bool x)
{
  auto elec_factor = 332.0;
  mm_charge_vector = mmc.energyinterface()->charges();
  auto aco_p = dynamic_cast<energy::interfaces::aco::aco_ff const*>(mmc.energyinterface());
  if (aco_p)
  {
    elec_factor = aco_p->params().general().electric;
  }

  // mol.in schreiben
  std::cout << std::setprecision(6);

  std::ofstream molstream{ "mol.in" };
  if (molstream)
  {
    auto const n_qm = qm_indices.size();
    molstream << '\n';
    molstream << n_qm << " 0\n";
    for (std::size_t i = 0; i < n_qm; ++i)
    {
      double qi{};
      for (std::size_t j = 0; j < mm_charge_vector.size(); ++j)
      {
        auto d = len(coords->xyz(qm_indices[i] - 1) - coords->xyz(mm_indices[j] - 1));
        qi += mm_charge_vector[j] / d;
      }
      qi *= elec_factor;
      molstream << "0 0 0 0 " << qi << "\n";
    }
    molstream.close();
  }
  else
  {
    throw std::runtime_error("Cannot write mol.in file.");
  }
  
  qm_energy = qmc.g();

  if (Config::get().energy.mopac.delete_input) std::remove("mol.in");
  
  ww_calc(x);

  if (x)
  {
    mm_energy = mmc.g();
    auto new_grad = c_gradient + vdw_gradient;
    for (auto g : new_grad)
    {
      std::cout << "WW-Gradienten: " << g << '\n';
    }

    auto g_qm = qmc.g_xyz();
    auto g_mm = mmc.g_xyz();

    for (auto && qmi : qm_indices)
    {
      new_grad[qmi - 1] += g_qm[new_indices_qm[qmi - 1]];
      std::cout << "QM-Gradient von " << qmi << "=" << g_qm[new_indices_qm[qmi - 1]] << '\n';
    }

    for (auto && mmi : mm_indices)
    {
      new_grad[mmi - 1] += g_mm[new_indices_mm[mmi - 1]];
      std::cout << "MM-Gradient von " << mmi << "=" << g_mm[new_indices_mm[mmi - 1]] << '\n';
    }
    coords->swap_g_xyz(new_grad);
  }
  else
  {
    mm_energy = mmc.e();
  }

  this->energy = qm_energy + mm_energy + c_energy + vdw_energy;
 
  return energy;
}

void energy::interfaces::qmmm::QMMM::ww_calc(bool x)
{
  auto elec_factor = 332.0;
  auto aco_p = dynamic_cast<energy::interfaces::aco::aco_ff const*>(mmc.energyinterface());
  if (aco_p)
  {
    elec_factor = aco_p->params().general().electric;
  }
  qm_charge_vector = qmc.energyinterface()->charges();
  c_gradient.assign(coords->size(), coords::r3{});
  vdw_gradient.assign(coords->size(), coords::r3{});

  std::size_t i2 = 0u;
  auto const & p = cparams.vdwc_matrices();
  auto const & p_vdw = p[5];
  std::size_t i3 = 1;
  vdw_energy = 0;
  c_energy = 0;
  for (auto i : qm_indices)
  {
    std::cout << "I: " << i << '\n';
    auto z = qmc.atoms(i3-1).energy_type();
    std::cout << "TI: " << z << '\n';
    std::size_t i_vdw = cparams.type(z, tinker::potential_keys::VDW);
    std::cout << "CTI: " << i_vdw << '\n';
    coords::float_type charge_i = qm_charge_vector[i2];
    std::cout << "QM_CHARGE: " << charge_i << '\n';
    std::size_t j2 = 0u;
    std::size_t j3 = 1;
    for (auto j : mm_indices)
    {
      std::cout << "J: " << j << '\n';
      auto e_type = mmc.atoms(j3-1).energy_type();
      std::cout << "TJ: " << e_type << '\n';
      std::size_t j_vdw = cparams.type(e_type, tinker::potential_keys::VDW);
      std::cout << "CTJ: " << j_vdw << '\n';
      auto const & p_ij = p_vdw(i_vdw, j_vdw);
      std::cout <<"Parameter: "<< p_ij << '\n';
      
      coords::float_type charge_j = mm_charge_vector[j2];
      std::cout << "MM_CHARGE: " << charge_j << '\n';
      auto r_ij = (coords->xyz(j - 1) - coords->xyz(i - 1));
      std::cout << "Verbindungsvektor: " << r_ij << '\n';
      coords::float_type d = len(r_ij);
      set_distance(d);
      std::cout << "Länge: " << d << '\n';
      coords::float_type b = (charge_i*charge_j) / d * elec_factor;
      c_energy += b;
      //EVDW
      auto R_r = std::pow(p_ij.R / d, 6);
      
      if (cparams.general().radiustype.value == ::tinker::parameter::radius_types::T::SIGMA)
      {
        vdw_energy += R_r*p_ij.E*(R_r - 1.0);
      }
      else if (cparams.general().radiustype.value == ::tinker::parameter::radius_types::T::R_MIN)
      {
        vdw_energy += R_r*p_ij.E*(R_r - 2.0);
      }
      else
      {
        throw std::runtime_error("no valid radius_type");
      }
      if (x)
      {
        coords::float_type db = b / d;
        auto c_gradient_ij = r_ij * db / d;
        std::cout << "Db: " << db << '\n';
        c_gradient[i - 1] += c_gradient_ij;
        c_gradient[j - 1] -= c_gradient_ij;

        coords::float_type const V = p_ij.E*R_r;

        if (cparams.general().radiustype.value == ::tinker::parameter::radius_types::T::SIGMA)
        {
          auto vdw_r_grad = (V / d)*(6.0 - 12.0 * R_r);
          auto vdw_gradient_ij = r_ij*vdw_r_grad;
          std::cout << "vdw_g: " << vdw_gradient_ij;
          vdw_gradient[i - 1] += vdw_gradient_ij;
          vdw_gradient[j - 1] -= vdw_gradient_ij;
        }
        else 
        {
          auto vdw_gradient_ij = r_ij*12 * (V / d)*(1.0 - R_r);
          vdw_gradient[i - 1] += vdw_gradient_ij;
          vdw_gradient[j - 1] -= vdw_gradient_ij;
        }
      }
      ++j2;
      ++j3;
    }
    ++i2;
    ++i3;
  }
}


energy::interface_base * energy::interfaces::qmmm::QMMM::clone(coords::Coordinates * c) const
{
  QMMM * tmp = new QMMM (*this, c);
  return tmp;
}

energy::interface_base * energy::interfaces::qmmm::QMMM::move(coords::Coordinates * c)
{
  QMMM * tmp = new QMMM(std::move(*this), c);
  return tmp;
}


void energy::interfaces::qmmm::QMMM::swap(interface_base &rhs)
{
  swap(dynamic_cast<QMMM&>(rhs));
}

void energy::interfaces::qmmm::QMMM::swap(QMMM &rhs)
{
  interface_base::swap(rhs);
  qm_indices.swap(rhs.qm_indices);
  mm_indices.swap(rhs.mm_indices);
  
  new_indices_mm.swap(rhs.new_indices_mm);
  new_indices_qm.swap(rhs.new_indices_qm);
  qmc.swap(rhs.qmc);
  mmc.swap(rhs.mmc);
  cparams.swap(rhs.cparams);
  qm_charge_vector.swap(rhs.qm_charge_vector);
  mm_charge_vector.swap(rhs.mm_charge_vector);
}

// update structure (account for topology or rep change)
void energy::interfaces::qmmm::QMMM::update(bool const skip_topology)
{
  if (!skip_topology)
  {

    *this = QMMM(this->coords);
  }
  //
}

coords::float_type energy::interfaces::qmmm::QMMM::g()
{
  return qmmm_calc(true);
}

coords::float_type energy::interfaces::qmmm::QMMM::e()
{
  return qmmm_calc(false);
}

coords::float_type energy::interfaces::qmmm::QMMM::h()
{
  throw std::runtime_error("no QMMM-function yet");
}

coords::float_type energy::interfaces::qmmm::QMMM::o()
{
  throw std::runtime_error("no QMMM-function yet");
}

std::vector<coords::float_type> energy::interfaces::qmmm::QMMM::charges() const
{
  std::vector<coords::float_type> v;
  v.assign(coords->size(), 0);
  std::size_t i = 0;
  std::size_t j = 0;
  if (qm_charge_vector.size() + mm_charge_vector.size() == coords->size())
  {
    for (auto && qmi : qm_indices)
    {
      v[qmi - 1] = qm_charge_vector[i];
      i++;
    }

    for (auto && mmi : mm_indices)
    {
      v[mmi - 1] = mm_charge_vector[j];
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
  S << std::right << std::setw(24) << "C";
  S << std::right << std::setw(24) << "VDW";
  S << std::right << std::setw(24) << "TOTAL";
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::QMMM::print_E_short(std::ostream &S, bool const endline) const
{
  S << '\n';
  S << std::right << std::setw(24) << qm_energy;
  S << std::right << std::setw(24) << mm_energy;
  S << std::right << std::setw(24) << c_energy;
  S << std::right << std::setw(24) << vdw_energy;
  S << std::right << std::setw(24) << energy;
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::QMMM::print_gnuplot(std::ostream &S, bool const endline) const
{
  
  S << std::right << std::setw(24) << distance;
  S << std::right << std::setw(24) << c_energy;
  if (endline) S << '\n';
}
void energy::interfaces::qmmm::QMMM::print_G_tinkerlike(std::ostream &S, bool const endline) const
{
  if (endline) S << '\n';
  throw std::runtime_error("no QMMM-function yet");
  
}
void energy::interfaces::qmmm::QMMM::to_stream(std::ostream &S) const
{
  S << '\n';
  throw std::runtime_error("no QMMM-function yet");
}