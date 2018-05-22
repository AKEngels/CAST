/**
CAST 3
bla
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

namespace energy
{
  namespace interfaces
  {
    /**namespace for QM/MM interface*/
    namespace oniom
    {
      /**QMMM interface class*/
      class ONIOM
        : public interface_base
      {

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

        interface_base * clone(coords::Coordinates*) const;
        interface_base * move(coords::Coordinates*);

        void swap(interface_base &);
        void swap(ONIOM &);

        /** update structure (account for topology or rep change)*/
        void update(bool const skip_topology = false);

        /** Energy function*/
        coords::float_type e() override;
        /** Energy+Gradient function */
        coords::float_type g() override;
        /** Energy+Hessian function*/
        coords::float_type h() override;
        /** Optimization in the interface or interfaced program (not existent for this interface)*/
        coords::float_type o() override;

        /** Return charges (for QM und MM atoms) */
        std::vector<coords::float_type> charges() const override
		{
			throw std::runtime_error("function not implemented\n");
		}
        /**overwritten function, should not be called*/
        std::vector<coords::Cartesian_Point> get_g_coul_mm() const override
        {
          throw std::runtime_error("function not implemented\n");
        }
        /**overwritten function, should not be called*/
        coords::Gradients_3D get_link_atom_grad() const override
        {
          throw std::runtime_error("function not implemented\n");
        }
        /**overwritten function, should not be called*/
        std::string get_id() const override
        {
          throw std::runtime_error("function not implemented\n");
        }
        /**prints total energy (not implemented)*/
        void print_E(std::ostream&) const  final override;
        /**prints 'headline' for energies*/
        void print_E_head(std::ostream&,
          bool const endline = true) const  final override;
        /**prints partial energies*/
        void print_E_short(std::ostream&,
          bool const endline = true) const  final override;
        /**function not implemented*/
        void to_stream(std::ostream&) const;

      private:

		///**indizes of QM atoms*/
  //      std::vector<size_t> qm_indices;
  //      /**indizes of MM atoms*/
  //      std::vector<size_t> mm_indices;

  //      /**vector of length total number of atoms
  //      only those elements are filled whose position corresponds to QM atoms
  //      they are filled with successive numbers starting from 0
  //      purpose: faciliate mapping between total coordinates object and subsystems*/
  //      std::vector<size_t> new_indices_qm;
  //      /**vector of length total number of atoms
  //      only those elements are filled whose position corresponds to MM atoms
  //      they are filled with successive numbers starting from 0
  //      purpose: faciliate mapping between total coordinates object and subsystems*/
  //      std::vector<size_t> new_indices_mm;

  //      /**coordinates object for QM part*/
  //      coords::Coordinates qmc;
  //      /**coordinates object for MM part*/
  //      coords::Coordinates mmc;

  //      /**atom charges of QM atoms*/
  //      std::vector<double> qm_charge_vector;
  //      /**atom charges of MM atoms*/
  //      std::vector<double> mm_charge_vector;

  //      /**van der Waals interaction energy between QM and MM atoms*/
  //      coords::float_type vdw_energy;
  //      /**energy of only QM system*/
  //      coords::float_type qm_energy;
  //      /**energy of only MM system*/
  //      coords::float_type mm_energy;

  //      /**gradients of electrostatic interaction between QM and MM atoms
  //      for MOPAC gradients on QM as well as on MM atoms
  //      for GAUSSIAN and DFTB+ only gradients on MM atoms (those on QM atoms are calculated by QM program)*/
  //      coords::Gradients_3D c_gradient;
  //      /**gradients of van der waals interaction energy between QM and MM atoms*/
  //      coords::Gradients_3D vdw_gradient;

  //      /**information needed to calculate coulomb gradients on MM atoms
  //      for GAUSSIAN: electric field from gaussian calculation for QM and MM atoms (first QM, then MM)
  //      only those of the MM atoms are used to calculate the gradients of the electrostatic interaction
  //      between QM and MM atoms on the MM atoms
  //      for DFTB+: coulomb gradients on MM atoms due to QM atoms*/
  //      std::vector<coords::Cartesian_Point> g_coul_mm;

  //      /**checks if all bonds are still intact (bond length smaller than 1.2 sum of covalent radii)*/
  //      bool check_bond_preservation(void) const;

  //      /**checks if there is a minimum atom distance (0.3 Angstrom) between atoms*/
  //      bool check_atom_dist(void) const;

      };
    }
  }
}
