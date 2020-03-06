#include <cstddef>

#include "energy_int_qmmm_a.h"
#include "Scon/scon_utility.h"

::tinker::parameter::parameters energy::interfaces::qmmm::QMMM_A::tp;

energy::interfaces::qmmm::QMMM_A::QMMM_A(coords::Coordinates* cp) :
  interface_base(cp),
  qm_indices(Config::get().energy.qmmm.qm_systems[0]),
  mm_indices(get_mm_atoms(cp->size())),
  new_indices_qm(make_new_indices(qm_indices, cp->size())),
  new_indices_mm(make_new_indices(mm_indices, cp->size())),
  link_atoms(create_link_atoms(qm_indices, cp, tp, Config::get().energy.qmmm.linkatom_sets[0])),
  qmc(make_small_coords(cp, qm_indices, new_indices_qm, Config::get().energy.qmmm.qminterface, "QM system: ", Config::get().energy.qmmm.qm_to_file, link_atoms)),
  mmc(make_small_coords(cp, mm_indices, new_indices_mm, Config::get().energy.qmmm.mminterface, "MM system: ")),
  index_of_QM_center(get_index_of_QM_center(Config::get().energy.qmmm.centers[0], qm_indices, coords)),
  qm_energy(0.0), mm_energy(0.0), vdw_energy(0.0), bonded_energy(0.0), coulomb_energy(0.0)
{
  // read force field parameter file if necessary
  if (!tp.valid()) tp.from_file(Config::get().get().general.paramFilename);

  if (coords->size() != 0)     // only do "real" initialisation if there are coordinates (interface is first created without)
  {
    // get force field parameters
    std::vector<std::size_t> types;
    for (auto atom : (*cp).atoms()) scon::sorted::insert_unique(types, atom.energy_type());
    cparams = tp.contract(types);
    torsionunit = cparams.torsionunit();

    // prepare bonded QM/MM
    prepare_bonded_qmmm();

    // check if cutoff is okay for periodics
    if (Config::get().periodics.periodic)  
    {
      double const min_cut = std::min({ Config::get().periodics.pb_box.x(), Config::get().periodics.pb_box.y(), Config::get().periodics.pb_box.z() }) / 2.0;
      if (Config::get().energy.qmmm.mminterface == config::interface_types::T::OPLSAA || Config::get().energy.qmmm.mminterface == config::interface_types::T::AMBER)
      {
        if (Config::get().energy.cutoff > min_cut) {
          std::cout << "\n!!! WARNING! Forcefield cutoff too big! Your cutoff should be smaller than " << min_cut << "! !!!\n\n";
        }
      }
      if (Config::get().energy.qmmm.cutoff > min_cut) {
        std::cout << "\n!!! WARNING! QM/MM cutoff too big! Your cutoff should be smaller than " << min_cut << "! !!!\n\n";
      }
    }
    // check if correct number of link atom types is given
    if (link_atoms.size() != Config::get().energy.qmmm.linkatom_sets[0].size())  // 
    {
      std::cout << "Wrong number of link atom types given. You have " << link_atoms.size() << " in the following order:\n";
      for (auto& l : link_atoms) std::cout << "QM atom: " << l.qm + 1 << ", MM atom: " << l.mm + 1 << "\n";
      throw std::runtime_error("wrong number of link atom types");
    }
  }
}

energy::interfaces::qmmm::QMMM_A::QMMM_A(QMMM_A const& rhs,
  coords::Coordinates* cobj) : interface_base(cobj),
  cparams(rhs.cparams), qm_indices(rhs.qm_indices), mm_indices(rhs.mm_indices), charge_indices(rhs.charge_indices),
  new_indices_qm(rhs.new_indices_qm), new_indices_mm(rhs.new_indices_mm), link_atoms(rhs.link_atoms),
  qmc(rhs.qmc), mmc(rhs.mmc), qmmm_bonds(rhs.qmmm_bonds), qmmm_angles(rhs.qmmm_angles), qmmm_dihedrals(rhs.qmmm_dihedrals),
	torsionunit(rhs.torsionunit), index_of_QM_center(rhs.index_of_QM_center), qm_energy(rhs.qm_energy), mm_energy(rhs.mm_energy),
  vdw_energy(rhs.vdw_energy), bonded_energy(rhs.bonded_energy), coulomb_energy(rhs.coulomb_energy),
  coulomb_gradient(rhs.coulomb_gradient), vdw_gradient(rhs.vdw_gradient), bonded_gradient(rhs.bonded_gradient),
	g_coul_mm(rhs.g_coul_mm), total_atom_charges(rhs.total_atom_charges)
{
  interface_base::operator=(rhs);
}

energy::interfaces::qmmm::QMMM_A::QMMM_A(QMMM_A&& rhs, coords::Coordinates* cobj)
  : interface_base(cobj),
  cparams(std::move(rhs.cparams)),
  qm_indices(std::move(rhs.qm_indices)), mm_indices(std::move(rhs.mm_indices)), charge_indices(std::move(rhs.charge_indices)),
  new_indices_qm(std::move(rhs.new_indices_qm)), new_indices_mm(std::move(rhs.new_indices_mm)),
  link_atoms(std::move(rhs.link_atoms)), qmc(std::move(rhs.qmc)), mmc(std::move(rhs.mmc)),
	qmmm_bonds(std::move(rhs.qmmm_bonds)), qmmm_angles(std::move(rhs.qmmm_angles)), qmmm_dihedrals(std::move(rhs.qmmm_dihedrals)),
	torsionunit(std::move(rhs.torsionunit)), index_of_QM_center(std::move(rhs.index_of_QM_center)),
  qm_energy(std::move(rhs.qm_energy)), mm_energy(std::move(rhs.mm_energy)), vdw_energy(std::move(rhs.vdw_energy)),
  bonded_energy(std::move(rhs.bonded_energy)), coulomb_energy(std::move(rhs.coulomb_energy)),
  coulomb_gradient(std::move(rhs.coulomb_gradient)), vdw_gradient(std::move(rhs.vdw_gradient)), bonded_gradient(std::move(rhs.bonded_gradient)),
	g_coul_mm(std::move(rhs.g_coul_mm)), total_atom_charges(std::move(rhs.total_atom_charges))
{
  interface_base::operator=(rhs);
}

/**function to find force field parameters for bonds, angles and so on between QM and MM system*/
void energy::interfaces::qmmm::QMMM_A::find_parameters()
{
  for (auto& b : qmmm_bonds)  // find bond parameters
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
    if (found == false) std::cout << "Parameters for bond between " << b.a + 1 << " and " << b.b + 1 << " not found.\n";
  }

  for (auto& a : qmmm_angles)  // find angle parameters
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
    if (found == false) std::cout << "Parameters for angle made of " << a.a + 1 << ", " << a.c + 1 << " and " << a.b + 1 << " not found.\n";
  }

  for (auto& d : qmmm_dihedrals)  // find parameters for dihedrals
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
      std::cout << "Parameters for dihedral made of " << d.a + 1 << ", " << d.c1 + 1 << ", " << d.c2 + 1 << " and " << d.b + 1 << " not found.\n";
    }
  }
}

/**function to find bonds, angles and so on between QM and MM system*/
void energy::interfaces::qmmm::QMMM_A::find_bonds_etc()
{
  // find bonds between QM and MM region
  for (auto qma : qm_indices)
  {
    for (auto b : coords->atoms().atom(qma).bonds())
    {
      if (scon::sorted::exists(mm_indices, b))
      {
        bonded::Bond bond(qma, b);
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
        if (!is_in(angle, qmmm_angles))
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
        if (!is_in(angle, qmmm_angles))
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
        if (!is_in(dihed, qmmm_dihedrals))
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
        if (!is_in(dihed, qmmm_dihedrals))
        {
          qmmm_dihedrals.push_back(dihed);
        }
      }
    }
  }
}

void energy::interfaces::qmmm::QMMM_A::prepare_bonded_qmmm()
{
  find_bonds_etc();    // find bonds, angles and dihedrals between QM and MM region
  find_parameters();   // find force field parameters for energy calculation

  if (Config::get().general.verbosity > 3)  // Output
  {
    std::cout << "QM/MM-Bonds\n";
    for (auto const& b : qmmm_bonds) std::cout << b << "\n";
    std::cout << "QM/MM-Angles\n";
    for (auto const& a : qmmm_angles) std::cout << a << "\n";
    std::cout << "QM/MM-Dihedrals\n";
    for (auto const& d : qmmm_dihedrals) std::cout << d << "\n";
  }
}

void energy::interfaces::qmmm::QMMM_A::update_representation()
{
  std::size_t qi = 0u;   // update positions of QM system
  for (auto i : qm_indices)
  {
    qmc.move_atom_to(qi, coords->xyz()[i], true);
    ++qi;
  }

  for (auto& l : link_atoms) l.calc_position(coords); // update positions of link atoms in QM system
  for (auto i = 0u; i < link_atoms.size(); ++i)
  {
    int index = qm_indices.size() + i;
    coords::cartesian_type& new_pos = link_atoms[i].position;
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
coords::float_type energy::interfaces::qmmm::QMMM_A::qmmm_calc(bool const if_gradient)
{
  // ############ UPDATE STUFF ##############################

  integrity = true;
  update_representation();  // update positions of QM and MM subsystem to those of coordinates object

  // ############### CREATE MM CHARGES ######################

  if (Config::get().coords.atom_charges.size() > mm_indices.size())
  {
    total_atom_charges = Config::get().coords.atom_charges;  // save all amber charges (for mechanical embedding)
    select_from_atomcharges(mm_indices);
  }

  if (Config::get().energy.qmmm.zerocharge_bonds != 0)
  {
    std::vector<double> mm_charge_vector = mmc.energyinterface()->charges();

    charge_indices.clear();
    add_external_charges(qm_indices, mm_charge_vector, mm_indices, link_atoms, charge_indices, coords, index_of_QM_center);
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
  catch (...)
  {
    std::cout << "QM programme failed. Treating structure as broken.\n";
    integrity = false;  // if QM programme fails: integrity is destroyed
  }
  clear_external_charges();  // no external charges for MM calculation
  Config::set().periodics.periodic = periodic;   // reset periodic boundaries

  ww_calc(if_gradient);  // calculate interactions between QM and MM part

  if (integrity == true)
  {
    if (if_gradient)  // if gradients should be calculated
    {
      mm_energy = mmc.g(); // get energy for MM part

      // get gradients: QM + MM + vdW + Coulomb + bonded
      auto new_grad = vdw_gradient + coulomb_gradient + bonded_gradient;  // vdW + Coulomb + bonded
      auto g_qm = qmc.g_xyz(); // QM
      auto g_mm = mmc.g_xyz(); // MM

      // calculate gradients from link atoms
      for (auto j = 0u; j < link_atoms.size(); j++)
      {
        LinkAtom l = link_atoms[j];
        auto link_atom_grad = g_qm[qm_indices.size() + j];
        coords::r3 G_QM, G_MM;

        calc_link_atom_grad(l, link_atom_grad, coords, G_QM, G_MM);

        new_grad[l.qm] += G_QM;
        new_grad[l.mm] += G_MM;

        if (Config::get().general.verbosity > 4)
        {
          std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
          std::cout << "It causes a gradient on QM atom " << G_QM << " and on MM atom " << G_MM << ".\n";
        }
      }

      // get gradients from pure QM and MM system
      for (auto&& qmi : qm_indices)
      {
        new_grad[qmi] += g_qm[new_indices_qm[qmi]];
      }
      for (auto&& mmi : mm_indices)
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
double energy::interfaces::qmmm::QMMM_A::calc_bonded(bool const if_gradient)
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
int energy::interfaces::qmmm::QMMM_A::calc_vdw(unsigned const qm, unsigned const mm) const
{
  for (auto b : qmmm_bonds)
  {
    if (qm == b.a && mm == b.b) return 0;
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
@param if_gradient: true if gradients should be calculated, false if not*/
void energy::interfaces::qmmm::QMMM_A::ww_calc(bool const if_gradient)
{
  // bonded interactions
  bonded_gradient.assign(coords->size(), coords::r3{});
  bonded_energy = calc_bonded(if_gradient);

  // preparation for calculation of non-bonded interactions
  std::vector<double> qm_charge_vector;                                     // vector with all charges of QM atoms
  std::vector<double> mm_charge_vector = mmc.energyinterface()->charges();  // vector with all charges of MM atoms
  try {
    qm_charge_vector = qmc.energyinterface()->charges(); // still link atoms in it
    for (auto i = 0u; i < link_atoms.size(); ++i) {      // remove charges from link atoms
      qm_charge_vector.pop_back();
    }
  }
  catch (...) { integrity = false; }

  if (integrity == true)
  {
    //########## reset vdw gradients and energies #################################

    coulomb_gradient.assign(coords->size(), coords::r3{});
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
          else if (d > s) scaling = ((c * c - d * d) * (c * c - d * d) * (c * c + 2 * d * d - 3 * s * s)) / ((c * c - s * s) * (c * c - s * s) * (c * c - s * s));
        }

        double R_0;  // r_min or sigma
        if (cparams.general().radiustype.value ==
          ::tinker::parameter::radius_types::T::SIGMA)
        {
          R_0 = sqrt(vparams_i.r * vparams_j.r);  // sigma
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
            vdw = 4 * R_r * epsilon * (R_r - 1.0) * scaling;
            if (calc_modus == 2) vdw = vdw / cparams.general().vdw_scale.factor(3);
          }
          else if (cparams.general().radiustype.value ==
            ::tinker::parameter::radius_types::T::R_MIN)
          {
            vdw = R_r * epsilon * (R_r - 2.0) * scaling;
            if (calc_modus == 2) vdw = vdw / cparams.general().vdw_scale.factor(3);
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
              auto vdw_r_grad_sigma = (V / d) * (6.0 - 12.0 * R_r);
              auto vdw_gradient_ij_sigma = ((r_ij * vdw_r_grad_sigma) / d) * scaling;

              if (d > Config::get().energy.switchdist && d <= Config::get().energy.cutoff)  // additional gradient as charge changes with distance
              {
                double const& c = Config::get().energy.cutoff;
                double const& s = Config::get().energy.switchdist;

                auto E_unscaled = 4 * R_r * epsilon * (R_r - 1.0);
                auto deriv_S = (-12 * d * (c * c - d * d) * (d * d - s * s)) / ((c * c - s * s) * (c * c - s * s) * (c * c - s * s));
                auto abs_grad = E_unscaled * deriv_S;

                vdw_gradient_ij_sigma += (r_ij / d) * abs_grad;  // give additional gradient a direction
              }

              if (calc_modus == 2) vdw_gradient_ij_sigma = vdw_gradient_ij_sigma / cparams.general().vdw_scale.factor(3);
              vdw_gradient[i] -= vdw_gradient_ij_sigma;
              vdw_gradient[j] += vdw_gradient_ij_sigma;
            }
            else
            {
              coords::float_type const V = epsilon * R_r;
              auto vdw_r_grad_R_MIN = (V / d) * 12 * (1.0 - R_r);
              auto vdw_gradient_ij_R_MIN = ((r_ij * vdw_r_grad_R_MIN) / d) * scaling;

              if (d > Config::get().energy.switchdist && d <= Config::get().energy.cutoff)  // additional gradient as charge changes with distance
              {
                double const& c = Config::get().energy.cutoff;
                double const& s = Config::get().energy.switchdist;

                auto E_unscaled = R_r * epsilon * (R_r - 2.0);
                auto deriv_S = (-12 * d * (c * c - d * d) * (d * d - s * s)) / ((c * c - s * s) * (c * c - s * s) * (c * c - s * s));
                auto abs_grad = E_unscaled * deriv_S;

                vdw_gradient_ij_R_MIN += (r_ij / d) * abs_grad;  // give additional gradient a direction
              }

              if (calc_modus == 2) vdw_gradient_ij_R_MIN = vdw_gradient_ij_R_MIN / cparams.general().vdw_scale.factor(3);
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
          double const& ext_chg = get_external_charges()[i].original_charge;
          double const& scaling = get_external_charges()[i].scaled_charge / get_external_charges()[i].original_charge;
          auto chargesQM = qmc.energyinterface()->charges();

          double sum_of_QM_interactions{ 0.0 };   // calculate sum(Q_qm * Q_ext / r)
          for (auto j{ 0u }; j < qm_indices.size(); ++j)
          {
            double const QMcharge = chargesQM[j];
            coords::r3 MMpos{ get_external_charges()[i].x,  get_external_charges()[i].y,  get_external_charges()[i].z };
            double const dist = len(MMpos - coords->xyz(qm_indices[j]));
            sum_of_QM_interactions += (QMcharge * ext_chg * elec_factor) / dist;
          }

          // calculate additional gradient
          deriv_Q.x() = sum_of_QM_interactions * 4 * (coords->xyz(index_of_QM_center).x() - get_external_charges()[i].x) * std::sqrt(scaling) / (c * c);
          deriv_Q.y() = sum_of_QM_interactions * 4 * (coords->xyz(index_of_QM_center).y() - get_external_charges()[i].y) * std::sqrt(scaling) / (c * c);
          deriv_Q.z() = sum_of_QM_interactions * 4 * (coords->xyz(index_of_QM_center).z() - get_external_charges()[i].z) * std::sqrt(scaling) / (c * c);

          grad += deriv_Q;                             // gradient on external charge
          coulomb_gradient[index_of_QM_center] -= deriv_Q;   // gradient on center of QM region
        }

        coulomb_gradient[mma] += grad;
      }
    }

    else if (Config::get().energy.qmmm.zerocharge_bonds == 0)  // mechanical embedding -> coulomb energy and gradients from MM interface
    {
      auto c_params = cparams.charges();     // get forcefield parameters (charge)
      double current_coul_energy, current_coul_grad, qm_charge;  // some variables

      auto i2{ 0u };
      for (auto i : qm_indices)  // for every QM atom
      {
        if (Config::get().general.single_charges)  // amber charges
        {
          qm_charge = total_atom_charges[i];
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

            current_coul_energy = ((qm_charge * mm_charge * cparams.general().electric) / d);   // calculate energy (unscaled)
            if (calc_modus == 2) {
              current_coul_energy = current_coul_energy / cparams.general().chg_scale.factor(3);
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
              coulomb_gradient[i] += new_grad;
              coulomb_gradient[j] -= new_grad;
            }
          }
          j2++;
        }
        i2++;
      }
    }
  }
}


energy::interface_base* energy::interfaces::qmmm::QMMM_A::clone(coords::Coordinates* c) const
{
  QMMM_A* tmp = new QMMM_A(*this, c);
  return tmp;
}

energy::interface_base* energy::interfaces::qmmm::QMMM_A::move(coords::Coordinates* c)
{
  QMMM_A* tmp = new QMMM_A(std::move(*this), c);
  return tmp;
}


void energy::interfaces::qmmm::QMMM_A::swap(interface_base& rhs)
{
  swap(dynamic_cast<QMMM_A&>(rhs));
}

void energy::interfaces::qmmm::QMMM_A::swap(QMMM_A& rhs)
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
  coulomb_gradient.swap(rhs.coulomb_gradient);
  vdw_gradient.swap(rhs.vdw_gradient);
}

// update structure (account for topology or rep change)
void energy::interfaces::qmmm::QMMM_A::update(bool const skip_topology)
{
  if (!skip_topology)
  {
    *this = QMMM_A(this->coords);
  }
  else
  {
    update_representation();
    qmc.energy_update(true);
    mmc.energy_update(true);
  }
}

coords::float_type energy::interfaces::qmmm::QMMM_A::g()
{
  integrity = coords->check_structure();
  if (integrity == true) return qmmm_calc(true);
  else return 0;
}

coords::float_type energy::interfaces::qmmm::QMMM_A::e()
{
  integrity = coords->check_structure();
  if (integrity == true) return qmmm_calc(false);
  else return 0;
}

coords::float_type energy::interfaces::qmmm::QMMM_A::h()
{
  throw std::runtime_error("no QMMM-function yet");
}

coords::float_type energy::interfaces::qmmm::QMMM_A::o()
{
  throw std::runtime_error("QMMM-cannot optimize");
}

std::vector<coords::float_type> energy::interfaces::qmmm::QMMM_A::charges() const
{
  std::vector<double> qm_charge_vector;                                     // vector with all charges of QM atoms
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
    for (auto&& qmi : qm_indices)
    {
      v[qmi] = qm_charge_vector[i];
      i++;
    }

    for (auto&& mmi : mm_indices)
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

void energy::interfaces::qmmm::QMMM_A::print_E(std::ostream&) const
{
  throw std::runtime_error("no QMMM-function yet");
}

void energy::interfaces::qmmm::QMMM_A::print_E_head(std::ostream& S, bool const endline) const
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

void energy::interfaces::qmmm::QMMM_A::print_E_short(std::ostream& S, bool const endline) const
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

void energy::interfaces::qmmm::QMMM_A::to_stream(std::ostream& S) const
{
  S << '\n';
  interface_base::to_stream(S);
  throw std::runtime_error("no QMMM-function yet");
}

//// STUFF FOR BONDED STRUCTS

/**function to calculate force field energy and gradients (copied from energy_int_aco.cc)
@param cp: pointer to coordinates object
@param gradients: gradients that are filled during gradient calculation (necessary but not used also in only-energy calculation)
@param grad: true if gradients shold be calculated*/
double energy::interfaces::qmmm::bonded::Bond::calc_energy(coords::Coordinates* cp, coords::Gradients_3D& gradients, bool const grad)
{
  double E(0.0);
  auto const bv(cp->xyz(a) - cp->xyz(b)); // r_ij (i=1, j=2)
  auto const d = len(bv);
  auto const r = d - ideal;
  auto dE = force * r;
  E += dE * r;  // kcal/mol
  if (grad == true)
  {
    dE *= 2;  // kcal/(mol*Angstrom)  gradient without direction
    dE /= d;  // kcal/(mol*A^2)   gradient divided by distance because later it is multiplied with it again
    auto const gv = bv * dE;   // "force" on atom i due to atom j (kcal/(mol*A)), gradient with direction
    gradients[a] += gv;   //(kcal/(mol*A))
    gradients[b] -= gv;
  }
  return E;
}

/** overloaded output operator for Bond */
std::ostream& energy::interfaces::qmmm::bonded::operator<<(std::ostream& os, const energy::interfaces::qmmm::bonded::Bond& b)
{
  os << b.a + 1 << " , " + b.b + 1 << " dist: " << b.ideal << ", force constant: " << b.force;
  return os;
}

/**function to calculate force field energy and gradients (copied from energy_int_aco.cc)
@param cp: pointer to coordinates object
@param gradients: gradients that are filled during gradient calculation (necessary but not used also in only-energy calculation)
@param grad: true if gradients shold be calculated*/
double energy::interfaces::qmmm::bonded::Angle::calc_energy(coords::Coordinates* cp, coords::Gradients_3D& gradients, bool grad)
{
  double E(0.0);
  auto av1(cp->xyz(a) - cp->xyz(c));
  auto av2(cp->xyz(b) - cp->xyz(c));
  auto const d(scon::angle(av1, av2).degrees() - ideal);
  auto const r(d * SCON_PI180);
  E += force * r * r;
  if (grad == true)
  {
    auto dE = force * r;
    coords::Cartesian_Point const cv(cross(av1, av2));
    coords::float_type const cvl(len(cv));
    dE *= 2.0 / cvl;
    coords::Cartesian_Point const gv1(cross(av1, cv) * (dE / dot(av1, av1)));
    coords::Cartesian_Point const gv2(cross(cv, av2) * (dE / dot(av2, av2)));
    gradients[a] += gv1;
    gradients[b] += gv2;
    gradients[c] += -(gv1 + gv2);
  }
  return E;
}

/** overloaded == operator for Angle */
bool energy::interfaces::qmmm::bonded::operator==(energy::interfaces::qmmm::bonded::Angle const& lhs, energy::interfaces::qmmm::bonded::Angle const& rhs)
{
  if (lhs.c == rhs.c)  // central atom has to be the same
  { // outer atoms can be switched
    if (lhs.a == rhs.a && rhs.b == lhs.b) return true;
    else if (lhs.a == rhs.b && lhs.b == rhs.a) return true;
  }
  return false;
}

/** overloaded output operator for Angle */
std::ostream& energy::interfaces::qmmm::bonded::operator<<(std::ostream& os, const energy::interfaces::qmmm::bonded::Angle& a)
{
  os << a.a + 1 << " , " << a.c + 1 << " , " << a.b + 1 << " angle: " << a.ideal << ", force constant: " << a.force;
  return os;
}

/**function to calculate force field energy and gradients (copied from energy_int_aco.cc)
@param cp: pointer to coordinates object
@param gradients: gradients that are filled during gradient calculation (necessary but not used also in only-energy calculation)
@param grad: true if gradients shold be calculated*/
double energy::interfaces::qmmm::bonded::Dihedral::calc_energy(coords::Coordinates* cp, double const torsionunit, coords::Gradients_3D& gradients, bool const grad)
{
  double E(0.0);
  // Get bonding vectors
  coords::Cartesian_Point const b01 = cp->xyz(c1) - cp->xyz(a);
  coords::Cartesian_Point const b12 = cp->xyz(c2) - cp->xyz(c1);
  coords::Cartesian_Point const b23 = cp->xyz(b) - cp->xyz(c2);
  // Cross terms
  coords::Cartesian_Point const t = cross(b01, b12);
  coords::Cartesian_Point const u = cross(b12, b23);
  // Get length and variations
  coords::float_type const tl2 = dot(t, t);
  coords::float_type const ul2 = dot(u, u);
  // ...
  coords::float_type const tlul = sqrt(tl2 * ul2);
  coords::float_type const r12 = len(b12);
  // cross of cross
  coords::Cartesian_Point const tu = cross(t, u);
  // scalar and length variations
  coords::float_type const cos_scalar0 = dot(t, u);
  coords::float_type const cos_scalar1 = tlul;
  coords::float_type const sin_scalar0 = dot(b12, tu);
  coords::float_type const sin_scalar1 = r12 * tlul;
  // Get multiple sine and cosine values
  coords::float_type cos[7], sin[7];
  cos[1] = cos_scalar0 / cos_scalar1;
  sin[1] = sin_scalar0 / sin_scalar1;
  for (auto j{ 2 }; j <= max_order; ++j)
  {
    std::size_t const k = j - 1;
    sin[j] = sin[k] * cos[1] + cos[k] * sin[1];
    cos[j] = cos[k] * cos[1] - sin[k] * sin[1];
  }

  coords::float_type tE(0.0), dE(0.0);
  for (auto j{ 0 }; j < number; ++j)
  {
    coords::float_type const F = forces[j] * torsionunit;
    std::size_t const k = orders[j];
    coords::float_type const l = std::abs(ideals[j]) > 0.0 ? -1.0 : 1.0;
    tE += F * (1.0 + cos[k] * l);
    dE += -static_cast<coords::float_type>(k)* F* sin[k] * l;
  }
  E += tE;
  if (grad == true)
  {
    coords::Cartesian_Point const b02 = cp->xyz(c2) - cp->xyz(a);
    coords::Cartesian_Point const b13 = cp->xyz(b) - cp->xyz(c1);

    coords::Cartesian_Point const dt(cross(t, b12) * (dE / (tl2 * r12)));
    coords::Cartesian_Point const du(cross(u, b12) * (-dE / (ul2 * r12)));

    coords::Cartesian_Point const vir1 = cross(dt, b12);
    coords::Cartesian_Point const vir2 = cross(b02, dt) + cross(du, b23);
    coords::Cartesian_Point const vir3 = cross(dt, b01) + cross(b13, du);
    coords::Cartesian_Point const vir4 = cross(du, b12);

    gradients[a] += vir1;
    gradients[c1] += vir2;
    gradients[c2] += vir3;
    gradients[b] += vir4;
  }
  return E;
}

// overloaded == operator for Dihedral
bool energy::interfaces::qmmm::bonded::operator==(energy::interfaces::qmmm::bonded::Dihedral const& lhs, energy::interfaces::qmmm::bonded::Dihedral const& rhs)
{
  if (lhs.c1 == rhs.c1 && lhs.c2 == rhs.c2 && lhs.a == rhs.a && lhs.b == rhs.b) return true;
  else if (lhs.c1 == rhs.c2 && lhs.c2 == rhs.c1 && lhs.a == rhs.b && lhs.b == rhs.a) return true;
  return false;
}

/** overloaded output operator for Dihedral */
std::ostream& energy::interfaces::qmmm::bonded::operator<<(std::ostream& os, const energy::interfaces::qmmm::bonded::Dihedral& d)
{
  os << "Atoms: " << d.a + 1 << " , " << d.c1 + 1 << " , " << d.c2 + 1 << " , " << d.b + 1 << "\n";
  os << "  max order: " << d.max_order << ", number: " << d.number << "\n";
  os << "  orders: ";
  for (auto const& o : d.orders) os << o << "  ";
  os << "\n  force constants: ";
  for (auto const& f : d.forces) os << f << "  ";
  os << "\n  ideals: ";
  for (auto const& i : d.ideals) os << i << "  ";
  return os;
}