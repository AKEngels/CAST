#include <cstddef>

#include "energy_int_3layer.h"
#include "Scon/scon_utility.h"
#include "alignment.h"


::tinker::parameter::parameters energy::interfaces::qmmm::THREE_LAYER::tp;

energy::interfaces::qmmm::THREE_LAYER::THREE_LAYER(coords::Coordinates* cp) :
  energy::interface_base(cp), qm_indices(Config::get().energy.qmmm.qm_systems[0]),
  qmse_indices(add_vectors(qm_indices, Config::get().energy.qmmm.seatoms, true)),
  new_indices_qm(make_new_indices(qm_indices, cp->size())),
  new_indices_qmse(make_new_indices(qmse_indices, cp->size())),
  link_atoms_small(create_link_atoms(qm_indices, cp, tp, std::vector<int>())),
  link_atoms_medium(create_link_atoms(qmse_indices, cp, tp, Config::get().energy.qmmm.linkatom_sets[0])),
  qmc_small(make_small_coords(cp, qm_indices, new_indices_qm, Config::get().energy.qmmm.qminterface, "Small system: ",
    Config::get().energy.qmmm.qm_to_file, link_atoms_small, "small_system.arc")),
  sec_small(make_small_coords(cp, qm_indices, new_indices_qm, Config::get().energy.qmmm.seinterface, "Small system: ",
    Config::get().energy.qmmm.qm_to_file, link_atoms_small, "small_system.arc")),
  sec_medium(make_small_coords(cp, qmse_indices, new_indices_qmse, Config::get().energy.qmmm.seinterface, "Intermediate system: ",
    Config::get().energy.qmmm.qm_to_file, link_atoms_medium, "intermediate_system.arc")),
  mmc_medium(make_small_coords(cp, qmse_indices, new_indices_qmse, Config::get().energy.qmmm.mminterface, "Intermediate system: ",
    Config::get().energy.qmmm.qm_to_file, link_atoms_medium, "intermediate_system.arc")),
  mmc_big(make_small_coords(cp, range(cp->size()), range(cp->size()), Config::get().energy.qmmm.mminterface, "Big system: ")),
  index_of_medium_center(get_index_of_QM_center(Config::get().energy.qmmm.centers[0], qmse_indices, coords)),
  index_of_small_center(get_index_of_QM_center(Config::get().energy.qmmm.small_center, qm_indices, coords)),
  qm_energy_small(0.0), se_energy_small(0.0), se_energy_medium(0.0), mm_energy_medium(0.0), mm_energy_big(0.0)
{
  // should own optimizer be used?
  if (Config::get().energy.qmmm.opt) optimizer = true;
  else optimizer = false;

  // set correct total charges for sec_small and mmc_medium
  sec_small.energyinterface()->charge = qmc_small.energyinterface()->charge;         
  mmc_medium.energyinterface()->charge = sec_medium.energyinterface()->charge;

  // set number of link atoms for MOPAC
  if (Config::get().energy.qmmm.qminterface == config::interface_types::T::MOPAC) Config::set().energy.mopac.link_atoms = link_atoms_small.size();  
  if (Config::get().energy.qmmm.seinterface == config::interface_types::T::MOPAC) Config::set().energy.mopac.link_atoms = link_atoms_medium.size();

  // check if applied interfaces are valid for THREE_LAYER
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

  // if forcefield is desired either read parameters from file or throw error (no parameter file present)
  if (Config::get().energy.qmmm.mminterface == config::interface_types::T::OPLSAA || Config::get().energy.qmmm.qminterface == config::interface_types::T::AMBER)
  {
    if (!file_exists(Config::get().general.paramFilename)) throw std::runtime_error("You need a tinker-like parameterfile for your chosen forcefield.");
    else if (!tp.valid()) tp.from_file(Config::get().general.paramFilename);
  }

  if (coords->size() != 0)     // only do "real" initialisation if there are coordinates (interface is first created without)
  {
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

    // test if correct number of link atom types is given
    if (link_atoms_medium.size() != Config::get().energy.qmmm.linkatom_sets[0].size())  
    {                                                                              
      std::cout << "Wrong number of link atom types given. You have " << link_atoms_medium.size() << " in the following order:\n";
      for (auto& l : link_atoms_medium) std::cout << "QM atom: " << l.qm + 1 << ", MM atom: " << l.mm + 1 << "\n";
      std::cout << "This is assuming you are using a forcefield for your big system. \nIf you want to use one for the intermediate system talk to a CAST developer!\n";
      throw std::runtime_error("wrong number of link atom types");
    }

    // test if no atom is double in intermediate system (i. e. given both in QM and SE atoms)
    if (double_element(qmse_indices) == true) throw std::runtime_error("ERROR! You have at least one atom in QM as well as in SE atoms.");
  }
}

energy::interfaces::qmmm::THREE_LAYER::THREE_LAYER(THREE_LAYER const& rhs,
  coords::Coordinates* cobj) : interface_base(cobj),
  qm_indices(rhs.qm_indices), qmse_indices(rhs.qmse_indices),
  new_indices_qm(rhs.new_indices_qm), new_indices_qmse(rhs.new_indices_qmse), link_atoms_small(rhs.link_atoms_small),
  link_atoms_medium(rhs.link_atoms_medium), qmc_small(rhs.qmc_small), sec_small(rhs.sec_small),
  sec_medium(rhs.sec_medium), mmc_medium(rhs.mmc_medium), mmc_big(rhs.mmc_big),
  index_of_medium_center(rhs.index_of_medium_center), index_of_small_center(rhs.index_of_small_center),
  qm_energy_small(rhs.qm_energy_small), se_energy_small(rhs.se_energy_small), se_energy_medium(rhs.se_energy_medium),
  mm_energy_medium(rhs.mm_energy_medium), mm_energy_big(rhs.mm_energy_big)
{
  interface_base::operator=(rhs);
}

energy::interfaces::qmmm::THREE_LAYER::THREE_LAYER(THREE_LAYER&& rhs, coords::Coordinates* cobj)
  : interface_base(cobj),
  qm_indices(std::move(rhs.qm_indices)), qmse_indices(std::move(rhs.qmse_indices)),
  new_indices_qm(std::move(rhs.new_indices_qm)), new_indices_qmse(std::move(rhs.new_indices_qmse)),
  link_atoms_small(std::move(rhs.link_atoms_small)), link_atoms_medium(std::move(rhs.link_atoms_medium)),
  qmc_small(std::move(rhs.qmc_small)), sec_small(std::move(rhs.sec_small)), sec_medium(std::move(rhs.sec_medium)),
  mmc_medium(std::move(rhs.mmc_medium)), mmc_big(std::move(rhs.mmc_big)), 
  index_of_medium_center(std::move(rhs.index_of_medium_center)), 
  index_of_small_center(std::move(rhs.index_of_small_center)),
  qm_energy_small(std::move(rhs.qm_energy_small)), se_energy_small(std::move(rhs.se_energy_small)),
  se_energy_medium(std::move(rhs.se_energy_medium)), mm_energy_medium(std::move(rhs.mm_energy_medium)), 
  mm_energy_big(std::move(rhs.mm_energy_big))
{
  interface_base::operator=(rhs);
}


energy::interface_base* energy::interfaces::qmmm::THREE_LAYER::clone(coords::Coordinates* c) const
{
  THREE_LAYER* tmp = new THREE_LAYER(*this, c);
  return tmp;
}

energy::interface_base* energy::interfaces::qmmm::THREE_LAYER::move(coords::Coordinates* c)
{
  THREE_LAYER* tmp = new THREE_LAYER(std::move(*this), c);
  return tmp;
}


void energy::interfaces::qmmm::THREE_LAYER::swap(interface_base& rhs)
{
  swap(dynamic_cast<THREE_LAYER&>(rhs));
}

void energy::interfaces::qmmm::THREE_LAYER::swap(THREE_LAYER& rhs)
{
  interface_base::swap(rhs);
  qm_indices.swap(rhs.qm_indices);
  qmse_indices.swap(rhs.qmse_indices);
  new_indices_qm.swap(rhs.new_indices_qm);
  new_indices_qmse.swap(rhs.new_indices_qmse);
  link_atoms_small.swap(rhs.link_atoms_small);
  link_atoms_medium.swap(rhs.link_atoms_medium);
  qmc_small.swap(rhs.qmc_small);
  sec_small.swap(rhs.sec_small);
  sec_medium.swap(rhs.sec_medium);
  mmc_medium.swap(rhs.mmc_medium);
  mmc_big.swap(rhs.mmc_big);
  std::swap(index_of_medium_center, rhs.index_of_medium_center);
  std::swap(index_of_small_center, rhs.index_of_small_center);
  std::swap(qm_energy_small, rhs.qm_energy_small);
  std::swap(se_energy_small, rhs.se_energy_small);
  std::swap(se_energy_medium, rhs.se_energy_medium);
  std::swap(mm_energy_medium, rhs.mm_energy_medium);
  std::swap(mm_energy_big, rhs.mm_energy_big);
}

// update structure (account for topology or rep change)
void energy::interfaces::qmmm::THREE_LAYER::update(bool const skip_topology)
{
  if (!skip_topology)
  {
    *this = THREE_LAYER(this->coords);
  }
  else
  {
    update_representation();
    qmc_small.energy_update(true);
    sec_small.energy_update(true);
    sec_medium.energy_update(true);
    mmc_medium.energy_update(true);
    mmc_big.energy_update(true);
  }
}

void energy::interfaces::qmmm::THREE_LAYER::update_representation()
{
  std::size_t qi = 0u;    // update position of QM atoms in small systems
  for (auto i : qm_indices)
  {
    qmc_small.move_atom_to(qi, coords->xyz(i), true);
    sec_small.move_atom_to(qi, coords->xyz(i), true);
    ++qi;
  }
  for (auto& l : link_atoms_small) l.calc_position(coords); // update positions of link atoms in small systems
  for (auto i = 0u; i < link_atoms_small.size(); ++i)
  {
    int index = qm_indices.size() + i;
    coords::cartesian_type& new_pos = link_atoms_small[i].position;
    qmc_small.move_atom_to(index, new_pos, true);
    sec_small.move_atom_to(index, new_pos, true);
  }

  std::size_t si = 0u;    // update position of atoms in medium systems
  for (auto i : qmse_indices)
  {
    sec_medium.move_atom_to(si, coords->xyz(i), true);
    mmc_medium.move_atom_to(si, coords->xyz(i), true);
    ++si;
  }
  for (auto& l : link_atoms_medium) l.calc_position(coords); // update positions of link atoms in medium systems
  for (auto i = 0u; i < link_atoms_medium.size(); ++i)
  {
    int index = qmse_indices.size() + i;
    coords::cartesian_type& new_pos = link_atoms_medium[i].position;
    sec_medium.move_atom_to(index, new_pos, true);
    mmc_medium.move_atom_to(index, new_pos, true);
  }

  for (std::size_t mi = 0u; mi < coords->size(); mi++)   // update positions of all atoms in big system
  {
    mmc_big.move_atom_to(mi, coords->xyz()[mi], true);
  }
}

coords::float_type energy::interfaces::qmmm::THREE_LAYER::qmmm_calc(bool if_gradient)
{
  // ############ INITIALISATION AND UPDATE ##############################

  update_representation(); // update positions of QM and MM subsystems to those of coordinates object

  mm_energy_big = 0.0;     // set energies to zero
  mm_energy_medium = 0.0;
  se_energy_medium = 0.0;
  se_energy_small = 0.0;
  qm_energy_small = 0.0;
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

  // ############### CREATE EXTERNAL CHARGES FOR MEDIUM SYSTEM ######################

  std::vector<int> charge_indices;  // indizes of all atoms that are in charge_vector
  charge_indices.clear();

  if (Config::get().energy.qmmm.zerocharge_bonds != 0)
  {
    auto mmc_big_charges = mmc_big.energyinterface()->charges();
    if (mmc_big_charges.size() == 0) throw std::runtime_error("no charges found in MM interface");
    auto all_indices = range(coords->size());
    add_external_charges(qmse_indices, mmc_big_charges, all_indices, link_atoms_medium, charge_indices, coords, index_of_medium_center);
  }

  Config::set().periodics.periodic = false;

  // ############### SE ENERGY AND GRADIENTS FOR MEDIUM SYSTEM ######################
  try {
    if (!if_gradient)
    {
      se_energy_medium = sec_medium.e();  // get se energy for intermediate part 
    }
    else  // gradient calculation
    {
      se_energy_medium = sec_medium.g();    // get energy and calculate gradients
      auto g_se_medium = sec_medium.g_xyz();        // get gradients
      for (auto&& qsi : qmse_indices)
      {
        new_grads[qsi] += g_se_medium[new_indices_qmse[qsi]];
      }

      for (auto i = 0u; i < link_atoms_medium.size(); ++i)   // take into account link atoms
      {
        LinkAtom l = link_atoms_medium[i];

        coords::r3 g_qm, g_mm;        // divide link atom gradient to QM and MM atom
        auto link_atom_grad = g_se_medium[qmse_indices.size() + i];
        calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);
        new_grads[l.qm] += g_qm;
        new_grads[l.mm] += g_mm;

        if (Config::get().general.verbosity > 4)
        {
          std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
          std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
        }
      }
    }
    if (se_energy_medium == 0) integrity = false;

    if (Config::get().general.verbosity > 4)
    {
      std::cout << "SE energy of intermediate system: \n";
      sec_medium.e_head_tostream_short(std::cout);
      sec_medium.e_tostream_short(std::cout);
    }
  }
  catch (...)
  {
    std::cout << "SE programme (for intermediate system) failed. Treating structure as broken.\n";
    integrity = false;  // if SE programme fails: integrity is destroyed
  }

  // ############### ONLY AMBER: PREPARATION OF CHARGES FOR MEDIUM SYSTEM ################

  // temporarily: only QM charges, SE charges and those of link atoms in amber_charges
  std::vector<double> old_amber_charges;
  if (Config::get().general.single_charges)
  {
    old_amber_charges = Config::get().coords.atom_charges;                       // save old amber_charges
    select_from_atomcharges(qmse_indices);                        // only QM and SE charges in amber_charges
    for (auto i = 0u; i < link_atoms_medium.size(); ++i)                         // add charges of link atoms
    {
      double la_charge = sec_medium.energyinterface()->charges()[qmse_indices.size() + i]; // get charge
      Config::set().coords.atom_charges.push_back(la_charge);                      // add it to vector
    }
  }

  // ################ SAVE OUTPUT FOR BIG MM SYSTEM ########################################################

  save_outputfiles(Config::get().energy.qmmm.mminterface, mmc_big.energyinterface()->id, "big");

  // ############### MM ENERGY AND GRADIENTS FOR MEDIUM SYSTEM ######################

  try {
    if (!if_gradient)
    {
      mm_energy_medium = mmc_medium.e();  // calculate energy of intermediate MM system
    }
    else  // gradient calculation
    {
      mm_energy_medium = mmc_medium.g();     // get energy and calculate gradients
      auto g_mm_medium = mmc_medium.g_xyz(); // get gradients
      for (auto&& qsi : qmse_indices)
      {
        new_grads[qsi] -= g_mm_medium[new_indices_qmse[qsi]];
      }

      for (auto i = 0u; i < link_atoms_medium.size(); ++i)  // take into account link atoms
      {
        LinkAtom l = link_atoms_medium[i];

        coords::r3 g_qm, g_mm;             // divide link atom gradient to QM and MM atom
        auto link_atom_grad = g_mm_medium[qmse_indices.size() + i];
        calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);
        new_grads[l.qm] -= g_qm;
        new_grads[l.mm] -= g_mm;
        if (Config::get().general.verbosity > 4)
        {
          std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
          std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
        }
      }
    }
    if (mm_energy_medium == 0) integrity = false;

    if (Config::get().general.verbosity > 4)
    {
      std::cout << "MM energy of intermediate system: \n";
      mmc_medium.e_head_tostream_short(std::cout);
      mmc_medium.e_tostream_short(std::cout);
    }
  }
  catch (...)
  {
    std::cout << "MM programme (for intermediate system) failed. Treating structure as broken.\n";
    integrity = false;  // if MM programme fails: integrity is destroyed
  }

  // ############### GRADIENTS ON MM ATOMS DUE TO COULOMB INTERACTION WITH MEDIUM REGION ###

  if (if_gradient && integrity == true && Config::get().energy.qmmm.zerocharge_bonds != 0)
  {
    auto sec_medium_g_ext_charges = sec_medium.energyinterface()->get_g_ext_chg();
    auto mmc_medium_g_ext_charges = mmc_medium.energyinterface()->get_g_ext_chg();

    for (auto i = 0u; i < charge_indices.size(); ++i)
    {
      int mma = charge_indices[i];
      auto grad_sec = sec_medium_g_ext_charges[i];
      auto grad_mmc = mmc_medium_g_ext_charges[i];

      // additional gradients because charge also changes with position (only if cutoff is applied)
      coords::r3 derivQ_sec{ 0.0, 0.0, 0.0 }, derivQ_mmc{ 0.0, 0.0, 0.0 };  
      if (Config::get().energy.qmmm.cutoff != std::numeric_limits<double>::max())
      {
        double constexpr elec_factor = 332.06;
        double const& c = Config::get().energy.qmmm.cutoff;
        double const& ext_chg = get_external_charges()[i].original_charge;
        double const& scaling = get_external_charges()[i].scaled_charge / get_external_charges()[i].original_charge;
        auto chargesSE = sec_medium.energyinterface()->charges();
        auto chargesMM = mmc_medium.energyinterface()->charges();

        // calculate sum(Q_qm * Q_ext / r) 
        double sum_of_QM_interactions_sec{ 0.0 };
        double sum_of_QM_interactions_mmc{ 0.0 };
        for (auto j{ 0u }; j < qmse_indices.size(); ++j)
        {
          double const QMcharge_sec = chargesSE[j];
          double const QMcharge_mmc = chargesMM[j];
          coords::r3 MMpos{ get_external_charges()[i].x,  get_external_charges()[i].y, get_external_charges()[i].z };
          double const dist = len(MMpos - coords->xyz(qmse_indices[j]));
          sum_of_QM_interactions_sec += (QMcharge_sec * ext_chg * elec_factor) / dist;
          sum_of_QM_interactions_mmc += (QMcharge_mmc * ext_chg * elec_factor) / dist;
        }

        // additional gradient on external charge due to interaction with SEC_medium
        derivQ_sec.x() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_medium_center).x() - get_external_charges()[i].x) * std::sqrt(scaling) / (c * c);
        derivQ_sec.y() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_medium_center).y() - get_external_charges()[i].y) * std::sqrt(scaling) / (c * c);
        derivQ_sec.z() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_medium_center).z() - get_external_charges()[i].z) * std::sqrt(scaling) / (c * c);
        grad_sec += derivQ_sec;

        // additional gradient on external charge due to interaction with MMC_medium
        derivQ_mmc.x() = sum_of_QM_interactions_mmc * 4 * (coords->xyz(index_of_medium_center).x() - get_external_charges()[i].x) * std::sqrt(scaling) / (c * c);
        derivQ_mmc.y() = sum_of_QM_interactions_mmc * 4 * (coords->xyz(index_of_medium_center).y() - get_external_charges()[i].y) * std::sqrt(scaling) / (c * c);
        derivQ_mmc.z() = sum_of_QM_interactions_mmc * 4 * (coords->xyz(index_of_medium_center).z() - get_external_charges()[i].z) * std::sqrt(scaling) / (c * c);
        grad_mmc += derivQ_mmc;

        // additional gradient on QM atom that defines distance
        new_grads[index_of_medium_center] += (derivQ_mmc - derivQ_sec);
      }

      new_grads[mma] += grad_sec;
      new_grads[mma] -= grad_mmc;
    }
  }

  // ############### EXTERNAL CHARGES FOR SMALL SYSTEM ######################

  Config::set().coords.atom_charges = old_amber_charges;  // set AMBER charges back to total AMBER charges
  Config::set().periodics.periodic = periodic;

  if (Config::get().energy.qmmm.zerocharge_bonds != 0)
  {
    if (Config::get().energy.qmmm.emb_small == 0)   // if EEx: no external charges for small system
    {
      clear_external_charges();
    }

    else if (Config::get().energy.qmmm.emb_small > 1)  // EE+ and EE+X
    {
      clear_external_charges();
      charge_indices.clear();

      if (Config::get().energy.qmmm.emb_small == 3)    // only for EE+X
      {
        auto sec_medium_charges = sec_medium.energyinterface()->charges();
        if (sec_medium.size() == 0) throw std::runtime_error("no charges found in SE interface");
        add_external_charges(qm_indices, sec_medium_charges, qmse_indices, link_atoms_small, charge_indices, coords, index_of_small_center);   // add charges from SE atoms
      }

      auto mmc_big_charges = mmc_big.energyinterface()->charges();
      if (mmc_big_charges.size() == 0) throw std::runtime_error("no charges found in MM interface");
      auto all_indices = range(coords->size());
      add_external_charges(qmse_indices, mmc_big_charges, all_indices, link_atoms_small, charge_indices, coords, index_of_small_center);     // add charges from MM atoms
    }
  }

  Config::set().periodics.periodic = false;

  // ############### QM ENERGY AND GRADIENTS FOR SMALL SYSTEM ######################
  try {
    if (!if_gradient)
    {
      qm_energy_small = qmc_small.e();  // get qm energy  
    }
    else  // gradient calculation
    {
      qm_energy_small = qmc_small.g();    // get energy and calculate gradients
      auto g_qm_small = qmc_small.g_xyz();        // get gradients
      for (auto&& qmi : qm_indices)
      {
        new_grads[qmi] += g_qm_small[new_indices_qm[qmi]];
      }

      for (auto i = 0u; i < link_atoms_small.size(); ++i)   // take into account link atoms
      {
        LinkAtom l = link_atoms_small[i];

        coords::r3 g_qm, g_mm;        // divide link atom gradient to QM and MM atom
        auto link_atom_grad = g_qm_small[qm_indices.size() + i];
        calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);
        new_grads[l.qm] += g_qm;
        new_grads[l.mm] += g_mm;

        if (Config::get().general.verbosity > 4)
        {
          std::cout << "Link atom between " << l.qm + 1 << " and " << l.mm + 1 << " has a gradient " << link_atom_grad << ".\n";
          std::cout << "It causes a gradient on QM atom " << g_qm << " and on MM atom " << g_mm << ".\n";
        }
      }
    }
    if (qm_energy_small == 0) integrity = false;

    if (Config::get().general.verbosity > 4)
    {
      std::cout << "QM energy of small system: \n";
      qmc_small.e_head_tostream_short(std::cout);
      qmc_small.e_tostream_short(std::cout);
    }
  }
  catch (...)
  {
    std::cout << "QM programme (for small system) failed. Treating structure as broken.\n";
    integrity = false;  // if QM programme fails: integrity is destroyed
  }

  // ################ SAVE OUTPUT FOR MEDIUM SE SYSTEM ########################################################

  save_outputfiles(Config::get().energy.qmmm.seinterface, sec_medium.energyinterface()->id, "intermediate");

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
        calc_link_atom_grad(l, link_atom_grad, coords, g_qm, g_mm);
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
    auto qmc_g_ext_charges = qmc_small.energyinterface()->get_g_ext_chg();
    auto sec_small_g_ext_charges = sec_small.energyinterface()->get_g_ext_chg();

    for (auto i = 0u; i < charge_indices.size(); ++i)
    {
      int mma = charge_indices[i];
      auto grad_qmc = qmc_g_ext_charges[i];
      auto grad_sec = sec_small_g_ext_charges[i];

      // additional gradients because charge also changes with position (only if cutoff is applied)
      coords::r3 derivQ_qmc{ 0.0, 0.0, 0.0 }, derivQ_sec{ 0.0, 0.0, 0.0 };   
      if (Config::get().energy.qmmm.cutoff != std::numeric_limits<double>::max())
      {
        double constexpr elec_factor = 332.06;
        double const& c = Config::get().energy.qmmm.cutoff;
        double const& ext_chg = get_external_charges()[i].original_charge;
        double const& scaling = get_external_charges()[i].scaled_charge / get_external_charges()[i].original_charge;
        auto chargesQM = qmc_small.energyinterface()->charges();
        auto chargesSE = sec_small.energyinterface()->charges();

        // calculate sum(Q_qm * Q_ext / r) 
        double sum_of_QM_interactions_qmc{ 0.0 };
        double sum_of_QM_interactions_sec{ 0.0 };
        for (auto j{ 0u }; j < qm_indices.size(); ++j)
        {
          double const QMcharge_qmc = chargesQM[j];
          double const QMcharge_sec = chargesSE[j];
          coords::r3 MMpos{ get_external_charges()[i].x,  get_external_charges()[i].y,  get_external_charges()[i].z };
          double const dist = len(MMpos - coords->xyz(qm_indices[j]));
          sum_of_QM_interactions_qmc += (QMcharge_qmc * ext_chg * elec_factor) / dist;
          sum_of_QM_interactions_sec += (QMcharge_sec * ext_chg * elec_factor) / dist;
        }

        // additional gradient on external charge due to interaction with QMC
        derivQ_qmc.x() = sum_of_QM_interactions_qmc * 4 * (coords->xyz(index_of_small_center).x() - get_external_charges()[i].x) * std::sqrt(scaling) / (c * c);
        derivQ_qmc.y() = sum_of_QM_interactions_qmc * 4 * (coords->xyz(index_of_small_center).y() - get_external_charges()[i].y) * std::sqrt(scaling) / (c * c);
        derivQ_qmc.z() = sum_of_QM_interactions_qmc * 4 * (coords->xyz(index_of_small_center).z() - get_external_charges()[i].z) * std::sqrt(scaling) / (c * c);
        grad_qmc += derivQ_qmc;

        // additional gradient on external charge due to interaction with SEC_SMALL
        derivQ_sec.x() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_small_center).x() - get_external_charges()[i].x) * std::sqrt(scaling) / (c * c);
        derivQ_sec.y() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_small_center).y() - get_external_charges()[i].y) * std::sqrt(scaling) / (c * c);
        derivQ_sec.z() = sum_of_QM_interactions_sec * 4 * (coords->xyz(index_of_small_center).z() - get_external_charges()[i].z) * std::sqrt(scaling) / (c * c);
        grad_sec += derivQ_sec;

        // additional gradient on QM atom that defines distance
        new_grads[index_of_small_center] += (derivQ_sec - derivQ_qmc);
      }

      new_grads[mma] += grad_qmc;
      new_grads[mma] -= grad_sec;
    }
  }

  // ############### STUFF TO DO AT THE END OF CALCULATION ######################

  clear_external_charges();                              // clear vector -> no point charges in calculation of mmc_big
  Config::set().periodics.periodic = periodic;           // set back periodics
  if (file_exists("orca.gbw")) std::remove("orca.gbw");  // delete orca MOs for small system, otherwise orca will try to use them for medium system and fail

  if (coords->check_bond_preservation() == false) integrity = false;
  else if (coords->check_for_crashes() == false) integrity = false;

  if (if_gradient) coords->swap_g_xyz(new_grads);     // swap gradients into coordobj
  energy = mm_energy_big + se_energy_medium - mm_energy_medium + qm_energy_small - se_energy_small;
  return energy; // return total energy
}

void energy::interfaces::qmmm::THREE_LAYER::fix_qmse_atoms(coords::Coordinates& coordobj)
{
  for (std::size_t i = 0u; i < coordobj.size(); ++i) {
    if (is_in(i, qmse_indices) == true) coordobj.set_fix(i, true);             // fix all QMSE atoms
  }
  for (auto const& link : link_atoms_medium) coordobj.set_fix(link.mm, true);  // fix all M1 atoms
}

void energy::interfaces::qmmm::THREE_LAYER::fix_mm_atoms(coords::Coordinates& coordobj)
{
  for (std::size_t i = 0u; i < coordobj.size(); ++i) {
    if (is_in(i, qmse_indices) == false) coordobj.set_fix(i, true);             // fix all MM atoms
  }
  for (auto const& link : link_atoms_medium) coordobj.set_fix(link.mm, false);  // unfix M1 atoms
}

coords::float_type energy::interfaces::qmmm::THREE_LAYER::g()
{
  integrity = coords->check_structure();
  if (integrity == true) return qmmm_calc(true);
  else return 0;
}

coords::float_type energy::interfaces::qmmm::THREE_LAYER::e()
{
  integrity = coords->check_structure();
  if (integrity == true) return qmmm_calc(false);
  else return 0;
}

coords::float_type energy::interfaces::qmmm::THREE_LAYER::h()
{
  throw std::runtime_error("no hessian function implemented for this interface");
}

coords::float_type energy::interfaces::qmmm::THREE_LAYER::o()
{
  if (Config::get().optimization.local.method == config::optimization_conf::lo_types::INTERNAL) {
    throw std::runtime_error("Microiterations do not work with INTERNAL optimizer!");
  }

  // set optimizer to false in order to go into general o() function of coordinates object
  optimizer = false;

  // some variables we need
  double rmsd{ 0.0 };                       // current RMSD value
  double energy_old{ 0.0 };                 // energy before microiteration
  double dEnergy{ 0.0 };                    // current energy difference 
  coords::Coordinates oldC;                 // coordinates before microiteration
  std::vector<std::size_t> mm_iterations;   // number of MM optimization steps for each microiteration 
  std::vector<std::size_t> qm_iterations;   // number of Three-layer optimization steps for each microiteration 
  std::vector<double> energies;             // energy after each microiteration
  std::vector<double> rmsds;                // RMSD value for each microiteration
  std::size_t total_mm_iterations{ 0u };    // total number of MM optimization steps
  std::size_t total_qm_iterations{ 0u };    // total number of Three-layer optimization steps

  // some configuration stuff that needs to be saved if adaption of coulomb interactions is switched on
  bool original_single_charges = Config::get().general.single_charges;
  std::vector<double> original_atom_charges = Config::get().coords.atom_charges;

  // file for writing trace if desired
  std::ofstream trace("trace_microiterations.arc");

  do {    // microiterations

    // save coordinates from before microiteration
    oldC = *coords;
    energy_old = energy;

    if (Config::get().energy.qmmm.coulomb_adjust)
    {
      Config::set().coords.atom_charges = charges();
      if (Config::get().energy.qmmm.emb_small == 0) {  // if embedding scheme = EEx: charges of small system are the MM charges
        for (auto i{ 0u }; i < coords->size(); ++i) {  //                            as interaction between QM and MM system is calculated only by SE interface
          if (is_in(i, qm_indices)) {
            Config::set().coords.atom_charges[i] = sec_medium.energyinterface()->charges()[new_indices_qmse[i]];
          }
        }
      }
      Config::set().general.single_charges = true;
    }

    // optimize MM atoms with MM interface
    mmc_big.set_xyz(coords->xyz());
    fix_qmse_atoms(mmc_big);
    mmc_big.o();
    mmc_big.reset_fixation();
    mm_iterations.emplace_back(mmc_big.get_opt_steps());
    total_mm_iterations += mmc_big.get_opt_steps();

    if (Config::get().energy.qmmm.coulomb_adjust)
    {
      Config::set().coords.atom_charges = original_atom_charges;
      Config::set().general.single_charges = original_single_charges;
    }

    // optimize QM and SE atoms with Three-layer interface
    coords->set_xyz(mmc_big.xyz());
    fix_mm_atoms(*coords);
    energies.emplace_back(coords->o());
    coords->reset_fixation();
    qm_iterations.emplace_back(coords->get_opt_steps());
    total_qm_iterations += coords->get_opt_steps();

    // write structure into tracefile
    if (Config::get().energy.qmmm.write_opt) trace << coords::output::formats::tinker(*coords);

    // determine if convergence is reached
    rmsd = align::rmsd_aligned(oldC, *coords);
    rmsds.emplace_back(rmsd);
    dEnergy = energy_old - energy;
    if (Config::get().general.verbosity > 2) {
      std::cout << "RMSD of microiteration is " << std::setprecision(3) << rmsd << " and energy difference is " << dEnergy << " kcal/mol\n";
    }
  } while (rmsd > 0.01 || dEnergy > 0.1);

  // writing information into microiterations.csv
  std::ofstream out("microiterations.csv");
  out << "It.,MM,Three-layer,Energy,RMSD\n";
  for (auto i{ 0u }; i < rmsds.size(); ++i)
  {
    out << i + 1 << "," << mm_iterations[i] << "," << qm_iterations[i] << "," << energies[i] << "," << rmsds[i] << "\n";
  }
  out << "TOTAL," << total_mm_iterations << "," << total_qm_iterations << "," << energy << ",";
  out.close();

  // close file, set optimizer to true again and return energy
  trace.close();
  optimizer = true;
  return energy;
}

std::vector<coords::float_type> energy::interfaces::qmmm::THREE_LAYER::charges() const
{
  std::vector<coords::float_type> charges;
  for (std::size_t i{ 0u }; i < coords->size(); ++i)
  {
    if (is_in(i, qmse_indices))   // either QM or SE atom
    {
      if (is_in(i, qm_indices)) {  // QM atom
        charges.emplace_back(qmc_small.energyinterface()->charges()[new_indices_qm[i]]);
      }
      else {                       // SE atom
        charges.emplace_back(sec_medium.energyinterface()->charges()[new_indices_qmse[i]]);
      }
    }
    else {    // MM atom
      charges.emplace_back(mmc_big.energyinterface()->charges()[i]);
    }
  }
  return charges;
}

void energy::interfaces::qmmm::THREE_LAYER::print_E(std::ostream&) const
{
  throw std::runtime_error("function not implemented");
}

void energy::interfaces::qmmm::THREE_LAYER::print_E_head(std::ostream& S, bool const endline) const
{
  S << "QM-atoms: " << qm_indices.size() << '\n';
  S << "SE-atoms: " << qmse_indices.size() - qm_indices.size() << '\n';
  S << "MM-atoms: " << coords->size() - qmse_indices.size() << '\n';
  S << "Potentials\n";
  S << std::right << std::setw(24) << "MM_big";
  S << std::right << std::setw(24) << "MM_medium";
  S << std::right << std::setw(24) << "SE_medium";
  S << std::right << std::setw(24) << "SE_small";
  S << std::right << std::setw(24) << "QM";
  S << std::right << std::setw(24) << "TOTAL";
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::THREE_LAYER::print_E_short(std::ostream& S, bool const endline) const
{
  S << '\n';
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << mm_energy_big;
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << mm_energy_medium;
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << se_energy_medium;
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << se_energy_small;
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << qm_energy_small;
  S << std::fixed << std::setprecision(1) << std::right << std::setw(24) << energy;
  if (endline) S << '\n';
}


