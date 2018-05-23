#include"qmmm_helperfunctions.h"


std::vector<std::size_t> qmmm_helpers::get_mm_atoms(std::size_t const num_atoms)
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

  std::vector<std::size_t> qmmm_helpers::make_new_indices_qm(std::size_t const num_atoms)
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

  std::vector<std::size_t> qmmm_helpers::make_new_indices_mm(std::size_t const num_atoms, std::vector<std::size_t> const& mmi)
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

  /**creates coordobject for QM interface
  @param cp: coordobj for whole system (QM + MM)
  @param indices: indizes of QM atoms
  @param new_indices: new indizes (see new_indices_qm)*/
  coords::Coordinates qmmm_helpers::make_qm_coords(coords::Coordinates const * cp,
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
          if (is_in(b, indices)) at.bind_to(new_indices.at(b)); // only bind if bonding partner is also in QM coords
        }
        new_qm_atoms.add(at);
        pes.structure.cartesian.push_back(cp->xyz(a));
      }
      new_qm_coords.init_swap_in(new_qm_atoms, pes);
    }
    Config::set().general.energy_interface = tmp_i;
    return new_qm_coords;
  }

  /**creates coordobject for MM interface
  @param cp: coordobj for whole system (QM + MM)
  @param indices: indizes of MM atoms
  @param new_indices: new indizes (see new_indices_mm)*/
  coords::Coordinates qmmm_helpers::make_aco_coords(coords::Coordinates const * cp,
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
          if (is_in(b, indices)) at.bind_to(new_indices[b]);  // only bind if bonding partner is also in ACO coords
        }
        new_aco_atoms.add(at);
        pes.structure.cartesian.push_back(cp->xyz(a));
      }
      new_aco_coords.init_swap_in(new_aco_atoms, pes);
    }
    Config::set().general.energy_interface = tmp_i;
    return new_aco_coords;
  }

  /***/
  coords::Coordinates qmmm_helpers::make_mmbig_coords(coords::Coordinates const * cp)
  {
    auto tmp_i = Config::get().general.energy_interface;
    Config::set().general.energy_interface = Config::get().energy.qmmm.mminterface;
    coords::Coordinates new_aco_coords;
    coords::Atoms new_aco_atoms;
    coords::PES_Point pes;

    pes.structure.cartesian.reserve(cp->atoms().size());
    for (int a{ 0u }; a < cp->atoms().size(); ++a)
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
        at.bind_to(b);
      }
      new_aco_atoms.add(at);
      pes.structure.cartesian.push_back(cp->xyz(a));
    }
    new_aco_coords.init_swap_in(new_aco_atoms, pes);
    Config::set().general.energy_interface = tmp_i;
    return new_aco_coords;
  }

  coords::Coordinates qmmm_helpers::make_mmsmall_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices, std::vector<LinkAtom> link_atoms)
  {
    auto tmp_i = Config::get().general.energy_interface;
    Config::set().general.energy_interface = Config::get().energy.qmmm.mminterface;
    coords::Coordinates new_qm_coords;
    if (cp->size() >= indices.size())
    {
      coords::Atoms new_qm_atoms;
      coords::Atoms tmp_link_atoms;
      coords::PES_Point pes;
      pes.structure.cartesian.reserve(indices.size()+link_atoms.size());

      for (auto && a : indices)      // add and bind "normal" atoms
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
          if (is_in(b, indices)) at.bind_to(new_indices.at(b)); // only bind if bonding partner is also in QM coords
        }
        
        new_qm_atoms.add(at);
        pes.structure.cartesian.push_back(cp->xyz(a));
      }
      
      int counter = 0;
      for (auto a : link_atoms)  // create temporary vector with link atoms
      {
        coords::Atom current("H");
        current.set_energy_type(a.energy_type);
        current.bind_to(new_indices.at(a.mm));
        tmp_link_atoms.add(current);
        pes.structure.cartesian.push_back(a.position);

        for (int i{ 0u }; i < new_qm_atoms.size(); ++i)  // bind bonding partner in "normal" atoms to link atom
        {
          if (current.is_bound_to(i))
          {
            coords::Atom& current2 = new_qm_atoms.atom(i);
            current2.bind_to(new_qm_atoms.size() + counter);
			break; // link atom is hydrogen so only one bonding partner
          }
        }
        counter++;
      }
      for (auto a : tmp_link_atoms) new_qm_atoms.add(a);  // add link atoms
      
      new_qm_coords.init_swap_in(new_qm_atoms, pes);
    }

    if (Config::set().energy.qmmm.qm_to_file)  // if desired: write QM region into file
    {
      std::ofstream output("qm_region.arc");
      output << coords::output::formats::tinker(new_qm_coords);
    }
    
    Config::set().general.energy_interface = tmp_i;
    return new_qm_coords;
  }

