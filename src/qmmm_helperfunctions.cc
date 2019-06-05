#include"qmmm_helperfunctions.h"


std::vector<LinkAtom> qmmm_helpers::create_link_atoms(coords::Coordinates* coords, std::vector<size_t> const &qm_indices, tinker::parameter::parameters const &tp)
{
  std::vector<LinkAtom> links;
	unsigned type, counter = 0;

	for (auto q : qm_indices)
	{
		if (coords->atoms().size() == 0) break;  // this is necessary because this function will be called once in the beginning when coordinates object is not created yet
		for (auto b : coords->atoms().atom(q).bonds())
		{
			if (!is_in(b, qm_indices))
			{
				if (counter < Config::get().energy.qmmm.linkatom_types.size())
				{
					type = Config::get().energy.qmmm.linkatom_types[counter];
				}
				else type = 85;    // if atomtype not found -> 85 (should mostly be correct for OPLSAA force field)
				LinkAtom link(q, b, type, coords, tp);
				links.push_back(link);
				counter += 1;

				if (Config::get().general.verbosity > 3)
				{
					std::cout << "created link atom between QM atom " << q + 1 << " and MM atom " << b + 1 << " with atom type " << link.energy_type << ", position: " << link.position << "\n";
				}
			}
		}
	}
  return links;
}

void qmmm_helpers::calc_link_atom_grad(LinkAtom const &l, coords::r3 const &G_L, coords::Coordinates* coords, coords::r3 &G_QM, coords::r3 &G_MM)
{
  double x, y, z;

  double g = l.deq_L_QM / dist(coords->xyz(l.mm), coords->xyz(l.qm));
  coords::r3 n = (coords->xyz(l.mm) - coords->xyz(l.qm)) / dist(coords->xyz(l.mm), coords->xyz(l.qm));

  x = g * scon::dot(G_L, n) *n.x() + (1 - g)*G_L.x();
  y = g * scon::dot(G_L, n) *n.y() + (1 - g)*G_L.y();
  z = g * scon::dot(G_L, n) *n.z() + (1 - g)*G_L.z();
  coords::r3 new_grad(x, y, z);
  G_QM = new_grad;

  x = g * G_L.x() - g * scon::dot(G_L, n) * n.x();
  y = g * G_L.y() - g * scon::dot(G_L, n) * n.y();
  z = g * G_L.z() - g * scon::dot(G_L, n) * n.z(); 
  coords::r3 new_grad2(x, y, z);
  G_MM = new_grad2;
}


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

  std::vector<std::size_t> qmmm_helpers::make_new_indices(std::size_t const num_atoms, std::vector<std::size_t> const& indices)
  {
    std::vector<std::size_t> new_indices;
    new_indices.resize(num_atoms);
    if (num_atoms > 0u)
    {
      std::size_t current_index = 0u;

      for (auto && a : indices)
      {
				new_indices.at(a) = current_index++;
      }
    }
    return new_indices;
  }

  coords::Coordinates qmmm_helpers::make_small_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices, config::interface_types::T energy_interface, bool const write_into_file, std::vector<LinkAtom> const &link_atoms,
		std::string const& filename)
  {
    auto tmp_i = Config::get().general.energy_interface;
    Config::set().general.energy_interface = energy_interface;
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
          if (is_in(b, indices)) at.bind_to(new_indices.at(b)); // only bind if bonding partner is also in new subsystem
        }
        
        new_qm_atoms.add(at);
        pes.structure.cartesian.push_back(cp->xyz(a));
      }
      
      int counter = 0;
      for (auto a : link_atoms)  // create temporary vector with link atoms
      {
        coords::Atom current("H");
        current.set_energy_type(a.energy_type);
        current.bind_to(new_indices.at(a.qm));
        tmp_link_atoms.add(current);
        pes.structure.cartesian.push_back(a.position);

        for (auto i{ 0u }; i < new_qm_atoms.size(); ++i)  // bind bonding partner in "normal" atoms to link atom
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

    if (write_into_file)  // if desired: write QM region into file
    {
      std::ofstream output(filename);
      output << coords::output::formats::tinker(new_qm_coords);
    }
    
    Config::set().general.energy_interface = tmp_i;
    return new_qm_coords;
  }

	void qmmm_helpers::select_from_ambercharges(std::vector<std::size_t> const & indices)
	{
		std::vector<coords::float_type> c = Config::get().coords.amber_charges;  // get AMBER charges
		std::vector<coords::float_type> charges_temp;
		for (auto i = 0u; i<c.size(); i++)
		{
			if (is_in(i, indices))  // find atom charges for indizes
			{
				charges_temp.push_back(c[i]);   // add those charges to new vector
			}
		}
		Config::set().coords.amber_charges = charges_temp; // set new AMBER charges
	}

  void qmmm_helpers::move_periodics(coords::Cartesian_Point &current_coords, coords::Cartesian_Point const& center_of_QM)
  {
    // determine vector to QM system
    auto vec_to_QMcenter = current_coords - center_of_QM;

    // move if distance is bigger than half the box size
    static coords::Cartesian_Point const halfbox(Config::get().periodics.pb_box / 2.0);
    if (vec_to_QMcenter.x() > halfbox.x())
    {
      current_coords.x() -= Config::get().periodics.pb_box.x();
    }
    else if (vec_to_QMcenter.x() < -halfbox.x())
    {
      current_coords.x() += Config::get().periodics.pb_box.x();
    }
    if (vec_to_QMcenter.y() > halfbox.y())
    {
      current_coords.y() -= Config::get().periodics.pb_box.y();
    }
    else if (vec_to_QMcenter.y() < -halfbox.y())
    {
      current_coords.y() += Config::get().periodics.pb_box.y();
    }
    if (vec_to_QMcenter.z() > halfbox.z())
    {
      current_coords.z() -= Config::get().periodics.pb_box.z();
    }
    else if (vec_to_QMcenter.z() < -halfbox.z())
    {
      current_coords.z() += Config::get().periodics.pb_box.z();
    }
  }

	void qmmm_helpers::add_external_charges(std::vector<size_t> const &qm_indizes, std::vector<size_t> const &ignore_indizes, 
    std::vector<double> const &charges, std::vector<size_t> const &indizes_of_charges,
		std::vector<LinkAtom> const &link_atoms, std::vector<int> &charge_indizes, coords::Coordinates *coords)
	{
    coords::Cartesian_Point center_of_QM{ 0.0, 0.0, 0.0 };
    if (Config::get().periodics.periodic)    // calculate center of QM system
    {
      for (auto q : qm_indizes) center_of_QM += coords->xyz(q);
      center_of_QM = center_of_QM / qm_indizes.size();
    }

		for (auto i : indizes_of_charges)  // go through all atoms from which charges are looked at
		{
			bool use_charge = true;

			for (auto &l : link_atoms)      // look at every link atom
			{
				if (l.mm == i) use_charge = false;   // ignore those atoms that are connected to an atom to the "QM system"...
				else
				{
					if (Config::get().energy.qmmm.zerocharge_bonds > 1)   // if desired: also ignore atoms that are two bonds away from "QM system"
					{
						for (auto b : coords->atoms(i).bonds())
						{
              if (b == l.mm) use_charge = false;

              else
              {
                if (Config::get().energy.qmmm.zerocharge_bonds > 2)  // if desired: also ignore atoms that are three bonds away from "QM system"
                {
                  for (auto b2 : coords->atoms(b).bonds())
                  {
                    if (b2 == l.mm) use_charge = false;
                  }
                }
              }
						}
					}
				}
			}
			for (auto &qs : ignore_indizes) // ...and those which are destined to be ignored
			{
				if (qs == i) use_charge = false;
			}

			if (use_charge)  // for the other 
			{
				if (Config::get().energy.qmmm.cutoff != 0.0)  // if cutoff given: test if one "QM atom" is nearer than cutoff
				{
					use_charge = false;
          auto current_coords = coords->xyz(i);

          if (Config::get().periodics.periodic)    // if periodic boundaries -> move current_coords next to QM 
          {
            move_periodics(current_coords, center_of_QM);
          }

					for (auto qs : qm_indizes)
					{
						auto dist = len(current_coords - coords->xyz(qs));
						if (dist < Config::get().energy.qmmm.cutoff)
						{
							use_charge = true;
							break;
						}
					}
				}

				if (use_charge)  // if yes create a PointCharge and add it to vector
				{
					PointCharge new_charge;
					new_charge.charge = charges[find_index(i, indizes_of_charges)];
          auto current_xyz = coords->xyz(i);
          if (Config::get().periodics.periodic) move_periodics(current_xyz, center_of_QM);  // if periodics: move charge next to QM
					new_charge.set_xyz(current_xyz.x(), current_xyz.y(), current_xyz.z());
					Config::set().energy.qmmm.mm_charges.push_back(new_charge);

					charge_indizes.push_back(i);  // add index to charge_indices
				}
			}
		}
	}

	void qmmm_helpers::save_outputfiles(config::interface_types::T const &interface, std::string const &id, std::string const &systemname)
	{
		if (interface == config::interface_types::T::DFTB && Config::get().energy.dftb.verbosity > 0)
		{
			if (file_exists("dftb_in.hsd")) rename("dftb_in.hsd", ("dftb_in_"+systemname+".hsd").c_str());
			if (file_exists("output_dftb.txt")) rename("output_dftb.txt", ("output_dftb_"+systemname+".txt").c_str());
			if (file_exists("charges.dat")) rename("charges.dat", ("charges_"+systemname+".dat").c_str());
			if (file_exists("results.tag")) rename("results.tag", ("results_"+systemname+".tag").c_str());
		}
		if (interface == config::interface_types::T::MOPAC && Config::get().energy.mopac.delete_input == false)
		{
			if (file_exists(id + ".xyz")) rename((id + ".xyz").c_str(), (id + "_"+systemname+".xyz").c_str());
			if (file_exists(id + ".out")) rename((id + ".out").c_str(), (id + "_"+systemname+".out").c_str());
			if (file_exists(id + ".arc")) rename((id + ".arc").c_str(), (id + "_"+systemname+".arc").c_str());
			if (file_exists(id + "_sys.out")) rename((id + "_sys.out").c_str(), (id + "_"+systemname+"_sys.out").c_str());
			if (file_exists(id + ".xyz.out")) rename((id + ".xyz.out").c_str(), (id + "_"+systemname+".xyz.out").c_str());
			if (file_exists(id + ".xyz.aux")) rename((id + ".xyz.aux").c_str(), (id + "_"+systemname+".xyz.aux").c_str());
			if (file_exists("mol.in")) rename("mol.in", ("mol_"+systemname+".in").c_str());
		}
		if (interface == config::interface_types::T::GAUSSIAN && Config::get().energy.gaussian.delete_input == false)
		{
			if (file_exists(id + ".gjf")) rename((id + ".gjf").c_str(), (id + "_"+systemname+".gjf").c_str());
			if (file_exists(id + ".log")) rename((id + ".log").c_str(), (id + "_"+systemname+".log").c_str());
			if (file_exists(id + "_G_.gjf")) rename((id + "_G_.gjf").c_str(), (id + "_G_"+systemname+".gjf").c_str());
			if (file_exists(id + "_G_.log")) rename((id + "_G_.log").c_str(), (id + "_G_"+systemname+".log").c_str());
		}
		if (interface == config::interface_types::T::PSI4)
		{
			if (file_exists(id + "_inp.dat")) rename((id + "_inp.dat").c_str(), (id + "_"+systemname+"_inp.dat").c_str());
			if (file_exists(id + "_out.dat")) rename((id + "_out.dat").c_str(), (id + "_"+systemname+"_out.dat").c_str());
			if (file_exists("grid.dat")) rename("grid.dat", ("grid_"+systemname+".dat").c_str());
			if (file_exists("grid_field.dat")) rename("grid_field.dat", ("grid_field_"+systemname+".dat").c_str());
		}
    if (interface == config::interface_types::T::ORCA)
    {
      if (file_exists("orca.ges")) rename("orca.ges", ("orca_" + systemname + ".ges").c_str());
      if (file_exists("orca.prop")) rename("orca.prop", ("orca_" + systemname + ".prop").c_str());
      if (file_exists("orca.opt")) rename("orca.opt", ("orca_" + systemname + ".opt").c_str());
      if (file_exists("orca.trj")) rename("orca.trj", ("orca_" + systemname + ".trj").c_str());
      if (file_exists("orca.engrad")) rename("orca.engrad", ("orca_" + systemname + ".engrad").c_str());
      if (file_exists("orca_property.txt")) rename("orca_property.txt", ("orca_property_" + systemname + ".txt").c_str());
      if (file_exists("orca.xyz")) rename("orca.xyz", ("orca_" + systemname + ".xyz").c_str());
      if (file_exists("orca.hess")) rename("orca.hess", ("orca_" + systemname + ".hess").c_str());
      if (file_exists("orca.inp")) rename("orca.inp", ("orca_" + systemname + ".inp").c_str());
      if (file_exists("pointcharges.pc")) rename("pointcharges.pc", ("pointcharges_" + systemname + ".pc").c_str());
      if (file_exists("output_orca.txt")) rename("output_orca.txt", ("output_orca_" + systemname + ".txt").c_str());
      if (file_exists("orca.gbw")) std::remove("orca.gbw");       // this is important because otherwise orca will try to read MOs from other system
    }
	}