/**
CAST 3
energy_int_3layer.h
Purpose: three-layer-interface (analogous to ONIOM)

available interfaces: OPLSAA, AMBER, MOPAC, DFTB+, PSI4, GAUSSIAN, ORCA


@author Susanne Sauer
@version 1.0
*/

#pragma once
#include "coords.h"
#include "coords_io.h"
#include <vector>
#include "coords_atoms.h"
#include "energy_int_aco.h"
#include "energy_int_mopac.h"
#include "tinker_parameters.h"
#include "helperfunctions.h"
#include "modify_sk.h"
#include "qmmm_helperfunctions.h"

namespace energy
{
  namespace interfaces
  {
    /**namespace for ONIOM interface*/
    namespace three_layer
    {
      /**ONIOM interface class*/
      class THREE_LAYER
        : public interface_base
      {

        /**uncontracted forcefield parameters*/
        static ::tinker::parameter::parameters tp;

      public:

        /**Constructor*/
        THREE_LAYER(coords::Coordinates*);
        /**overloaded Constructor*/
        THREE_LAYER(THREE_LAYER const&, coords::Coordinates*);
        /**another overload of Constructor*/
        THREE_LAYER(THREE_LAYER&&, coords::Coordinates*);

        /*
        Energy class functions that need to be overloaded (for documentation see also energy.h)
        */

        interface_base * clone(coords::Coordinates*) const;
        interface_base * move(coords::Coordinates*);

        void swap(interface_base &);
        void swap(THREE_LAYER &);

        /** update structure (account for topology or rep change)*/
        void update(bool const skip_topology = false);

        /**sets the atom coordinates of the subsystems (QM and MM) to those of the whole coordobject*/
        void update_representation();

        /** Energy function*/
        coords::float_type e() override;
        /** Energy+Gradient function */
        coords::float_type g() override;
        /** Energy+Hessian function (not existent for this interface)*/
        coords::float_type h() override;
        /** Optimization in the interface or interfaced program (not existent for this interface)*/
        coords::float_type o() override;

        /** Return charges (for QM und MM atoms) */
        std::vector<coords::float_type> charges() const override
		    {
			    throw std::runtime_error("function not implemented\n");
		    }
        /**overwritten function, should not be called*/
        std::vector<coords::Cartesian_Point> get_g_ext_chg() const override
        {
          throw std::runtime_error("function not implemented\n");
        }

        /**prints total energy (not implemented)*/
        void print_E(std::ostream&) const  final override;
        /**prints 'headline' for energies*/
        void print_E_head(std::ostream&, bool const endline = true) const  final override;
        /**prints partial energies*/
        void print_E_short(std::ostream&, bool const endline = true) const  final override;
        /**function not implemented*/
        void to_stream(std::ostream&) const final override {};

      private:        

        /**calculates energies and gradients
        @param if_gradient: true if gradients should be calculated, false if not*/
        coords::float_type qmmm_calc(bool if_gradient);

		    /**indizes of QM atoms*/
        std::vector<size_t> qm_indices;
        /**indizes of QM + SE atoms*/
        std::vector<size_t> qm_se_indices;
        /**vector of length total number of atoms
        only those elements are filled whose position corresponds to QM atoms
        they are filled with successive numbers starting from 0
        purpose: faciliate mapping between total coordinates object and subsystems*/
        std::vector<size_t> new_indices_qm;
        /**vector of length total number of atoms
        only those elements are filled whose position corresponds to atoms of the middle system (i.e. QM + SE)
        they are filled with successive numbers starting from 0
        purpose: faciliate mapping between total coordinates object and subsystems*/
        std::vector<size_t> new_indices_middle;

        /**vector with link atoms for small system*/
        std::vector<LinkAtom> link_atoms_small;
        /**vector with link atoms for middle system*/
        std::vector<LinkAtom> link_atoms_middle;

        /**coordinates object for QM part*/
        coords::Coordinates qmc;
        /**SE coordinates object for QM part*/
        coords::Coordinates sec_small;
        /**SE coordinates object for middle system*/
        coords::Coordinates sec_middle;
        /**MM coordinates object for middle system*/
        coords::Coordinates mmc_middle;
        /**coordinates object for whole system*/
        coords::Coordinates mmc_big;

				/**atom index that determines center of middle region*/
				std::size_t index_of_middle_center;
				/**atom index that determines center of small region*/
				std::size_t index_of_small_center;
        
        /**energy of QM system*/
        coords::float_type qm_energy;
        /**energy of small SE system*/
        coords::float_type se_energy_small;
        /**energy of middle SE system*/
        coords::float_type se_energy_middle;
        /**energy of middle MM system*/
        coords::float_type mm_energy_middle;
        /**energy of big MM system*/
        coords::float_type mm_energy_big;

        /**checks if all bonds are still intact (bond length smaller than 1.2 sum of covalent radii)*/
        bool check_bond_preservation(void) const;

        /**checks if there is a minimum atom distance (0.3 Angstrom) between atoms*/
        bool check_atom_dist(void) const;

      };
    }
  }
}
