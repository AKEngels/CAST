/**
CAST 3
energy_int_qmmm.h
Purpose: QM/MM interface

This is a QM/MM interface between one of the forcefields OPLSAA, AMBER and CHARM with MOPAC.
Interactions between QM and MM part are done by electrostatic embedding (see Gerrit Groenhof "Introduction to QM/MM Simulations", figure 4).
Only non-bonded interactions are implemented so there must not be bonds between the QM and the MM part.

Gradients of coulomb interactions between QM and MM part are calculated by CAST using the derived charge distribution for QM atoms from MOPAC.
(see http://openmopac.net/manual/QMMM.html)

@version 1.0
*/

#pragma once
#include "coords.h"
#include <vector>
#include "coords_atoms.h"
#include "energy_int_aco.h"
#include "energy_int_mopac.h"
#include "tinker_parameters.h"

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

        /**???*/
        void update_representation();

      private:

        /**calculates interaction between QM and MM part
        energy is only vdW interactions, gradients are coulomb and vdW
        @param if_gradient: true if gradients should be calculated, false if not*/
        void ww_calc(bool);
        /**calculates energies and gradients
        @paran if_gradient: true if gradients should be calculated, false if not*/
        coords::float_type qmmm_calc(bool);

        std::vector<std::size_t> qm_indices;
        std::vector<std::size_t> mm_indices;
  
        std::vector<std::size_t> new_indices_qm;
        std::vector<std::size_t> new_indices_mm;
       
        coords::Coordinates qmc;
        coords::Coordinates mmc;
        
        std::vector<double> qm_charge_vector;
        std::vector<double> mm_charge_vector;

        coords::float_type vdw_energy;
        coords::float_type qm_energy;
        coords::float_type mm_energy;

        coords::Gradients_3D c_gradient;
        coords::Gradients_3D vdw_gradient;
               
        
        
        

        
      };
    }
  }
}
