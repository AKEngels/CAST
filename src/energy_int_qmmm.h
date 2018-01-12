/**
CAST 3
energy_int_qmmm.h
Purpose: QM/MM interface

This is a QM/MM interface between one of the forcefields OPLSAA, AMBER and CHARM with MOPAC or GAUSSIAN.
Interactions between QM and MM part are done by electrostatic embedding (see Gerrit Groenhof "Introduction to QM/MM Simulations", figure 4).
Only non-bonded interactions are implemented so there must not be bonds between the QM and the MM part.

MOPAC: Gradients of coulomb interactions between QM and MM part are calculated by CAST using the derived charge distribution for QM atoms from MOPAC.
(see http://openmopac.net/manual/QMMM.html)

GAUSSIAN: Gradients of coulomb interactions between QM and MM part on QM atoms are calculated by GAUSSIAN,
on MM atoms they are calculated by CAST using the electric field from GAUSSIAN. 
(see T. Okamoto et. al., A minimal implementation of the AMBER-GAUSSIAN interface for Ab Initio QM/MM-MD Simulation, DOI 10.1002/jcc.21678)

Attention: Problems occur if charged atoms are in MM part near QM part! This situation has to be avoided!!!

@author Sara Wirsing, Susanne Sauer
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

namespace energy
{
  namespace interfaces
  {
    namespace qmmm
    {
      class QMMM
        : public interface_base
      {

        static ::tinker::parameter::parameters tp;
        ::tinker::parameter::parameters cparams;

      public:

        /**Constructor*/
        QMMM(coords::Coordinates*);
        /**overloaded Constructor*/
        QMMM(QMMM const&, coords::Coordinates*);
        /**another overload of Constructor*/
        QMMM(QMMM&&, coords::Coordinates*);

        /**variable to save a distance*/
        double distance;

        /*
        Energy class functions that need to be overloaded (for documentation see also energy.h)
        */

        interface_base * clone(coords::Coordinates*) const;
        interface_base * move(coords::Coordinates*);

        void swap(interface_base &);
        void swap(QMMM &);

        /** update structure (account for topology or rep change)*/
        void update(bool const skip_topology = false);
        
        /** Energy function*/
        coords::float_type e() override;
        /** Energy+Gradient function */
        coords::float_type g() override;
        /** Energy+Hessian function*/
        coords::float_type h() override;
        /** Optimization in the intface or interfaced program (not existent for this interface)*/
        coords::float_type o() override;

        /** Return charges (for QM und MM atoms) */
        std::vector<coords::float_type> charges() const override;
        /**overwritten function, should not be called*/
        std::vector<coords::Cartesian_Point> get_el_field() const override
        {
          throw std::runtime_error("TODO: Implement electric field.\n");
        }
        /**overwritten function*/
        std::string get_id() const override { return "bullshit"; }
        /**prints total energy (not implemented)*/
        void print_E(std::ostream&) const override;
        /**prints 'headline' for energies*/
        void print_E_head(std::ostream&, 
          bool const endline = true) const override;
        /**???*/
        void print_gnuplot(std::ostream&, 
          bool const endline = true) const;
        /**prints partial energies*/
        void print_E_short(std::ostream&, 
          bool const endline = true) const override;
        /**prints gradients*/
        void print_G_tinkerlike(std::ostream&, 
          bool const aggregate = false) const override;
        /**function not implemented*/
        void to_stream(std::ostream&) const;

        /**sets distance
        @param d: new distance*/
        void set_distance(double d)
        {
          distance = d;
        }

        /**sets the atom coordinates of the subsystems (QM and MM) to those of the whole coordobject*/
        void update_representation();

      private:

        /**
        writes inputfile for MOPAC calculation (see http://openmopac.net/manual/QMMM.html)
        */
        void write_mol_in();
        /**writes inputfile for gaussian calculation*/
        void write_gaussian_in(char);

        /**calculates interaction between QM and MM part
        energy is only vdW interactions, gradients are coulomb and vdW
        @param if_gradient: true if gradients should be calculated, false if not*/
        void ww_calc(bool);
        /**calculates energies and gradients
        @paran if_gradient: true if gradients should be calculated, false if not*/
        coords::float_type qmmm_calc(bool);

        /**indizes of QM atoms*/
        std::vector<std::size_t> qm_indices;
        /**indizes of MM atoms*/
        std::vector<std::size_t> mm_indices;
  
        std::vector<std::size_t> new_indices_qm;
        std::vector<std::size_t> new_indices_mm;
       
        /**coordinates object for QM part*/
        coords::Coordinates qmc;
        /**coordinates object for MM part*/
        coords::Coordinates mmc;
        
        /**atom charges of QM atoms*/
        std::vector<double> qm_charge_vector;
        /**atom charges of MM atoms*/
        std::vector<double> mm_charge_vector;

        /**van der Waals interaction energy between QM and MM atoms*/
        coords::float_type vdw_energy;
        /**energy of only QM system*/
        coords::float_type qm_energy;
        /**energy of only MM system*/
        coords::float_type mm_energy;

        /**gradients of electrostatic interaction between QM and MM atoms
        for MOPAC gradients on QM as well as on MM atoms
        for GAUSSIAN only gradients on MM atoms (those on QM atoms are calculated by GAUSSIAN)*/
        coords::Gradients_3D c_gradient;
        /**gradients of van der waals interaction energy between QM and MM atoms*/
        coords::Gradients_3D vdw_gradient;

        /**electric field from gaussian calculation for QM and MM atoms(first QM, then MM)
        only those of the MM atoms are used to calculate the gradients of the electrostatic interaction
        between QM and MM atoms on the MM atoms*/
        std::vector<coords::Cartesian_Point> qm_electric_field;
               
        /**checks if all bonds are still intact (bond length smaller than 2.2 Angstrom)*/
        bool check_bond_preservation(void) const;

        /**checks if there is a minimum atom distance (0.3 Angstrom) between atoms*/
        bool check_atom_dist(void) const;

      };
    }
  }
}
