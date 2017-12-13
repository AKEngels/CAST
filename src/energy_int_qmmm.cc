#include <cstddef>

#include "energy_int_qmmm.h"
#include "scon_utility.h"

::tinker::parameter::parameters energy::interfaces::qmmm::QMMM::tp;

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
      for (unsigned i = 0; i < num_atoms; i++)
      {
        if (!scon::sorted::exists(Config::get().energy.qmmm.qmatoms, i)) 
        { 
          mm_atoms.emplace_back(i); 
        }
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
      for (auto&& a : Config::get().energy.qmmm.qmatoms)
      {
        new_indices_qm.at(a) = current_index++;
      }
    }
    return new_indices_qm;
  }

  std::vector<std::size_t> make_new_indices_mm(std::size_t const num_atoms, 
    std::vector<std::size_t> const& mmi)
  {
    std::vector<std::size_t> new_indices_mm;
    new_indices_mm.resize(num_atoms);
    if (num_atoms > 0u)
    {
      std::size_t current_index = 0u;

      for (auto && a : mmi)
      {
        new_indices_mm.at(a) = current_index++;
      }
    }
    return new_indices_mm;
  }

  coords::Coordinates make_qm_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices)
  {
    auto tmp_i = Config::get().general.energy_interface;
    Config::set().general.energy_interface = Config::get().energy.qmmm.qminterface;
    coords::Coordinates new_qm_coords;
    if (cp->size() >= indices.size())
    {
      coords::Atoms new_qm_atoms;
      coords::PES_Point pes;
      pes.structure.cartesian.reserve(indices.size());
      for (auto && a : indices)
      {
        auto && ref_at = (*cp).atoms().atom(a);
        coords::Atom at{ (*cp).atoms().atom(a).number() };
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
        new_qm_atoms.add(at);
        pes.structure.cartesian.push_back(cp->xyz(a));
      }
      new_qm_coords.init_swap_in(new_qm_atoms, pes);
    }
    Config::set().general.energy_interface = tmp_i;
    return new_qm_coords;
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
        auto && ref_at = (*cp).atoms().atom(a);
        coords::Atom at{ (*cp).atoms().atom(a).number() };
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
        pes.structure.cartesian.push_back(cp->xyz(a));
      }
      new_aco_coords.init_swap_in(new_aco_atoms, pes);
    }
    Config::set().general.energy_interface = tmp_i;
    return new_aco_coords;
  }
}

energy::interfaces::qmmm::QMMM::QMMM(coords::Coordinates * cp) :
  interface_base(cp),
  qm_indices(Config::get().energy.qmmm.qmatoms),
  mm_indices(get_mm_atoms(cp->size())),
  new_indices_qm(make_new_indices_qm(cp->size())),
  new_indices_mm(make_new_indices_mm(cp->size(), mm_indices)),
  qmc(make_qm_coords(cp, qm_indices, new_indices_qm)),
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
}

energy::interfaces::qmmm::QMMM::QMMM(QMMM const & rhs, 
  coords::Coordinates *cobj) : interface_base(cobj), 
  cparams(rhs.cparams), distance(rhs.distance), 
  qm_indices(rhs.qm_indices), mm_indices(rhs.mm_indices),
  new_indices_qm(rhs.new_indices_qm), new_indices_mm(rhs.new_indices_mm),
  qmc(rhs.qmc), mmc(rhs.mmc), qm_charge_vector(rhs.qm_charge_vector), 
  mm_charge_vector(rhs.mm_charge_vector), vdw_energy(rhs.vdw_energy),  
  qm_energy(rhs.qm_energy), mm_energy(rhs.mm_energy), vdw_gradient(rhs.vdw_gradient),
  c_gradient(rhs.c_gradient)
{
  interface_base::operator=(rhs);
}

energy::interfaces::qmmm::QMMM::QMMM(QMMM&& rhs, coords::Coordinates *cobj) 
  : interface_base(cobj),
  cparams(std::move(rhs.cparams)), distance(std::move(rhs.distance)),
  qm_indices(std::move(rhs.qm_indices)), mm_indices(std::move(rhs.mm_indices)),
  new_indices_qm(std::move(rhs.new_indices_qm)), 
  new_indices_mm(std::move(rhs.new_indices_mm)),
  qmc(std::move(rhs.qmc)), mmc(std::move(rhs.mmc)), 
  qm_charge_vector(std::move(rhs.qm_charge_vector)),
  mm_charge_vector(std::move(rhs.mm_charge_vector)),
  vdw_energy(std::move(rhs.vdw_energy)),
  qm_energy(std::move(rhs.qm_energy)), mm_energy(std::move(rhs.mm_energy)),
  c_gradient(std::move(rhs.c_gradient)), 
  vdw_gradient(std::move(rhs.vdw_gradient))
{
  interface_base::operator=(rhs);
}

void energy::interfaces::qmmm::QMMM::update_representation()
{
  std::size_t qi = 0u;
  for (auto i : qm_indices)
  {
    qmc.move_atom_to(qi, coords->xyz()[i], true);
    ++qi;
  }
  std::size_t mi = 0u;
  for (auto j : mm_indices)
  {
    mmc.move_atom_to(mi, coords->xyz()[j], true);
    ++mi;
  }

}

// mol.in schreiben (for MOPAC, see http://openmopac.net/manual/QMMM.html)
void energy::interfaces::qmmm::QMMM::write_mol_in()
{
  auto elec_factor = 332.0;
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
        auto d = len(coords->xyz(qm_indices[i]) - coords->xyz(mm_indices[j]));
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
}

// write gaussian inputfile
void energy::interfaces::qmmm::QMMM::write_gaussian_in(char calc_type)
{
  std::string id = qmc.energyinterface()->get_id();
  std::string outstring(id);
  outstring.append(".gjf");

  std::ofstream out_file(outstring.c_str(), std::ios_base::out);

  if (out_file)
  {
    if (Config::get().energy.gaussian.link.length() != 0) { //if no link commands are issued the line wil be skipped

      std::istringstream iss(Config::get().energy.gaussian.link);

      std::vector<std::string> splitted_str(
        std::istream_iterator<std::string>{iss},
        std::istream_iterator<std::string>{}
      );

      for (std::size_t i = 0; i < splitted_str.size(); i++)
      {
        out_file << '%' << splitted_str[i] << '\n';
      }

    }
    out_file << "# " << Config::get().energy.gaussian.method << " " << Config::get().energy.gaussian.basisset << " " << Config::get().energy.gaussian.spec << " " << "Charge ";

    switch (calc_type) {// to ensure the needed gaussian keywords are used in gausian inputfile for the specified calculation
    case 'g':
      out_file << " Force";
      break;
    }

    out_file << '\n';
    out_file << '\n';
    out_file << Config::get().general.outputFilename;
    out_file << '\n';
    out_file << '\n';
    out_file << Config::get().energy.gaussian.charge << " ";
    out_file << Config::get().energy.gaussian.multipl;
    out_file << '\n';
    out_file << coords::output::formats::xyz(qmc);
    out_file << '\n';
    out_file << '\n';
    for (std::size_t j = 0; j < mm_charge_vector.size(); ++j)  // writing additional point charges (from MM atoms)
    {
      out_file << coords->xyz(mm_indices[j]).x() << " " << coords->xyz(mm_indices[j]).y() << " " << coords->xyz(mm_indices[j]).z() << " " << mm_charge_vector[j] << "\n";
    }
    out_file.close();
  }
  else std::runtime_error("Writing Gaussian Inputfile failed.");
}

/**calculates energies and gradients
@paran if_gradient: true if gradients should be calculated, false if not*/
coords::float_type energy::interfaces::qmmm::QMMM::qmmm_calc(bool if_gradient)
{
  integrity = true;
  auto elec_factor = 332.0;
  
  mm_charge_vector = mmc.energyinterface()->charges();
  auto aco_p = dynamic_cast<energy::interfaces::aco::aco_ff const*>(mmc.energyinterface());
  if (aco_p)
  {
    elec_factor = aco_p->params().general().electric;
  }

  if (Config::get().energy.qmmm.qminterface == config::interface_types::T::MOPAC)
  {
    write_mol_in();
  }
  else if (Config::get().energy.qmmm.qminterface == config::interface_types::T::GAUSSIAN)
  {
    char calc_type = 'e';
    if (if_gradient == true) calc_type = 'g';
    write_gaussian_in(calc_type);
  }
  else throw std::runtime_error("Chosen QM interface not implemented for QM/MM!");
  
  update_representation();

  try {
    qm_energy = qmc.g();  // get energy for QM part and save gradients for QM part
  }
  catch(...)
  {
    integrity = false;  // if QM programme fails: integrity is destroyed
  }

  if (Config::get().energy.qmmm.qminterface == config::interface_types::T::MOPAC && Config::get().energy.mopac.delete_input) std::remove("mol.in");
  
  ww_calc(if_gradient);  // calculate interactions between QM and MM part

  if (integrity == true)
  {
    if (if_gradient)  // if gradients should be calculated
    {
      mm_energy = mmc.g(); // get energy for MM part

                           // get gradients are QM + MM + vdW + Coulomb
      auto new_grad = vdw_gradient + c_gradient;  // vdW + Coulomb
      auto g_qm = qmc.g_xyz(); // QM
      auto g_mm = mmc.g_xyz(); // MM

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

    this->energy = qm_energy + mm_energy + vdw_energy;
    if (check_bond_preservation() == false) integrity = false;
    else if (check_atom_dist() == false) integrity = false;
  }
  return energy;
}

/**calculates interaction between QM and MM part
energy is only vdW interactions, gradients are coulomb and vdW
@param if_gradient: true if gradients should be calculated, false if not*/
void energy::interfaces::qmmm::QMMM::ww_calc(bool if_gradient)
{
  // preparation for calculation
  auto elec_factor = 332.0;
  auto aco_p = dynamic_cast<energy::interfaces::aco::aco_ff const*>(mmc.energyinterface());
  if (aco_p)
  {
    elec_factor = aco_p->params().general().electric;
  }
  try { qm_charge_vector = qmc.energyinterface()->charges(); }
  catch (...) { integrity = false; }
  if (integrity == true)
  {
    c_gradient.assign(coords->size(), coords::r3{});
    vdw_gradient.assign(coords->size(), coords::r3{});

    std::size_t i2 = 0u;
    auto const & p = cparams.vdwc_matrices();
    auto const & p_vdw = p[5];
    vdw_energy = 0;

    for (auto i : qm_indices)  // for every QM atom
    {
      auto z = qmc.atoms(i2).energy_type();
      std::size_t i_vdw = cparams.type(z, tinker::potential_keys::VDW);
      coords::float_type charge_i = qm_charge_vector[i2];
      std::size_t j2 = 0u;
      for (auto j : mm_indices)  // for every MM atom
      {
        // preparation (parameters and so on)
        auto e_type = mmc.atoms(j2).energy_type();
        std::size_t j_vdw = cparams.type(e_type, tinker::potential_keys::VDW);
        auto const & p_ij = p_vdw(i_vdw, j_vdw);
        coords::float_type charge_j = mm_charge_vector[j2];
        auto r_ij = coords->xyz(j) - coords->xyz(i);
        coords::float_type d = len(r_ij);
        set_distance(d);
        coords::float_type b = (charge_i*charge_j) / d * elec_factor;

        // calculate vdW interaction
        auto R_r = std::pow(p_ij.R / d, 6);
        if (cparams.general().radiustype.value ==
          ::tinker::parameter::radius_types::T::SIGMA)
        {
          vdw_energy += R_r*p_ij.E*(R_r - 1.0);
        }
        else if (cparams.general().radiustype.value ==
          ::tinker::parameter::radius_types::T::R_MIN)
        {
          vdw_energy += R_r*p_ij.E*(R_r - 2.0);
        }
        else
        {
          throw std::runtime_error("no valid radius_type");
        }

        if (if_gradient)  // gradients
        {
          // gradients of coulomb interaction
          coords::float_type db = b / d;
          auto c_gradient_ij = r_ij * db / d;
          c_gradient[i] += c_gradient_ij;
          c_gradient[j] -= c_gradient_ij;

          // gradients of vdW interaction
          coords::float_type const V = p_ij.E*R_r;

          if (cparams.general().radiustype.value
            == ::tinker::parameter::radius_types::T::SIGMA)
          {
            auto vdw_r_grad_sigma = (V / d)*(6.0 - 12.0 * R_r);
            auto vdw_gradient_ij_sigma = (r_ij*vdw_r_grad_sigma) / d;
            vdw_gradient[i] -= vdw_gradient_ij_sigma;
            vdw_gradient[j] += vdw_gradient_ij_sigma;
          }
          else
          {
            auto vdw_r_grad_R_MIN = (V / d) * 12 * (1.0 - R_r);
            auto vdw_gradient_ij_R_MIN = (r_ij*vdw_r_grad_R_MIN) / d;
            vdw_gradient[i] -= vdw_gradient_ij_R_MIN;
            vdw_gradient[j] += vdw_gradient_ij_R_MIN;
          }
        }
        ++j2;
      }
      ++i2;
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
  qmc.swap(rhs.qmc);
  mmc.swap(rhs.mmc);
  qm_charge_vector.swap(rhs.qm_charge_vector);
  mm_charge_vector.swap(rhs.mm_charge_vector);
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
  throw std::runtime_error("QMMM-cannot optimize");
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
  //S << std::right << std::setw(24) << "C";
  S << std::right << std::setw(24) << "VDW";
  S << std::right << std::setw(24) << "TOTAL";
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::QMMM::print_E_short(std::ostream &S, bool const endline) const
{
  S << '\n';
  S << std::right << std::setw(24) << qm_energy;
  S << std::right << std::setw(24) << mm_energy;
  //S << std::right << std::setw(24) << c_energy;
  S << std::right << std::setw(24) << vdw_energy;
  S << std::right << std::setw(24) << energy;
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::QMMM::print_gnuplot(std::ostream &S, bool const endline) const
{
  
  S << std::right << std::setw(24) << distance;
  S << std::right << std::setw(24) << vdw_gradient + c_gradient;
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::QMMM::print_G_tinkerlike(std::ostream &S, bool const endline) const
{
  S << " Cartesian Gradient Breakdown over Individual Atoms :" << std::endl << std::endl;
  S << "  Type      Atom              dE/dX       dE/dY       dE/dZ          Norm" << std::endl << std::endl;
  for (std::size_t k = 0; k < coords->size(); ++k)
  {
    S << " Anlyt";
    S << std::right << std::setw(10) << k + 1U;
    S << "       ";
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).x();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).y();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).z();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4);
    S << std::sqrt(
      coords->g_xyz(k).x() * coords->g_xyz(k).x()
      + coords->g_xyz(k).y() * coords->g_xyz(k).y()
      + coords->g_xyz(k).z() * coords->g_xyz(k).z()) << std::endl;
  }
  
}

void energy::interfaces::qmmm::QMMM::to_stream(std::ostream &S) const
{
  S << '\n';
  interface_base::to_stream(S);
  throw std::runtime_error("no QMMM-function yet");
}

bool energy::interfaces::qmmm::QMMM::check_bond_preservation(void) const
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
        if (L > 2.2) return false;
      }
    }
  }
  return true;
}

bool energy::interfaces::qmmm::QMMM::check_atom_dist(void) const
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
