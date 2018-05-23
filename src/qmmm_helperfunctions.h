#ifndef QMMM_HELPERFUNCTIONS_H
#define QMMM_HELPERFUNCTIONS_H

#include<vector>
#include"configuration.h"
#include"coords.h"
#include"helperfunctions.h"

namespace qmmm_helpers
{

	inline std::vector<std::size_t> get_mm_atoms(std::size_t const num_atoms)
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

	inline std::vector<std::size_t> make_new_indices_qm(std::size_t const num_atoms)
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

	inline std::vector<std::size_t> make_new_indices_mm(std::size_t const num_atoms,
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

	/**creates coordobject for QM interface
	@param cp: coordobj for whole system (QM + MM)
	@param indices: indizes of QM atoms
	@param new_indices: new indizes (see new_indices_qm)*/
	inline coords::Coordinates make_qm_coords(coords::Coordinates const * cp,
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
	inline coords::Coordinates make_aco_coords(coords::Coordinates const * cp,
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

	/**creates coordobject for MM interface
	@param cp: coordobj for whole system (QM + MM)
	@param indices: indizes of MM atoms
	@param new_indices: new indizes (see new_indices_mm)*/
	inline coords::Coordinates make_mmbig_coords(coords::Coordinates const * cp)
	{
		auto tmp_i = Config::get().general.energy_interface;
		Config::set().general.energy_interface = Config::get().energy.qmmm.mminterface;
		coords::Coordinates new_aco_coords;
		coords::Atoms new_aco_atoms;
		coords::PES_Point pes;

			pes.structure.cartesian.reserve(cp->atoms().size());
			for (int a{ 0u }; a<cp->atoms().size();++a)
			{
				auto && ref_at = (*cp).atoms().atom(a);
				coords::Atom at{ (*cp).atoms().atom(a).number() };
				at.set_energy_type(ref_at.energy_type());
				new_aco_atoms.add(at);
				pes.structure.cartesian.push_back(cp->xyz(a));
			}
			new_aco_coords.init_swap_in(new_aco_atoms, pes);
		Config::set().general.energy_interface = tmp_i;
		return new_aco_coords;
	}
}

#endif
