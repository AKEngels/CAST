/**
CAST 3
energy_int_oniom.h
Purpose: subtractive QM/MM interface

This is a subtractive QM/MM interface between the interfaces OPLSAA, AMBER, MOPAC, DFTB+, PSI4, GAUSSIAN and ORCA.vel

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
    namespace oniom
    {
      /**ONIOM interface class*/
      class ONIOM
        : public interface_base
      {

        /**uncontracted forcefield parameters*/
        static ::tinker::parameter::parameters tp;

      public:

        /**Constructor*/
        ONIOM(coords::Coordinates*);
        /**overloaded Constructor*/
        ONIOM(ONIOM const&, coords::Coordinates*);
        /**another overload of Constructor*/
        ONIOM(ONIOM&&, coords::Coordinates*);

        /*
        Energy class functions that need to be overloaded (for documentation see also energy.h)
        */

        interface_base* clone(coords::Coordinates*) const;
        interface_base* move(coords::Coordinates*);

        void swap(interface_base&);
        void swap(ONIOM&);

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
        std::vector<coords::float_type> charges() const override;
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

        /**fix all QM atoms
        @coordobj: coordinates object where atoms should be fixed*/
        void fix_qm_atoms(coords::Coordinates& coordobj);

        /**fix all MM atoms
        @coordobj: coordinates object where atoms should be fixed*/
        void fix_mm_atoms(coords::Coordinates& coordobj);

        /**indizes of QM atoms (one element of outer vector for every QM system)*/
        std::vector<std::vector<size_t>> qm_indices;
        /**for each QM system: vector of length total number of atoms
        only those elements are filled whose position corresponds to QM atoms
        they are filled with successive numbers starting from 0
        purpose: faciliate mapping between total coordinates object and subsystems*/
        std::vector < std::vector<size_t> >new_indices_qm;

        /**vector with link atoms (one set of link atoms for every QM system)*/
        std::vector < std::vector<LinkAtom>> link_atoms;

        /**coordinates objects for QM parts*/
        std::vector < coords::Coordinates> qmc;
        /**MM coordinates objects for QM parts*/
        std::vector < coords::Coordinates> mmc_small;
        /**coordinates object for whole system*/
        coords::Coordinates mmc_big;

        /**atom index that determines center of QM region*/
        std::vector<std::size_t> QMcenter_indices;

        /**energy of QM system*/
        coords::float_type qm_energy;
        /**energy of small MM system*/
        coords::float_type mm_energy_small;
        /**energy of big MM system*/
        coords::float_type mm_energy_big;

        /**total number of QM systems*/
        std::size_t number_of_qm_systems;
      };
    }
  }
}
