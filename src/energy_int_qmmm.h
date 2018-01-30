/**
CAST 3
energy_int_qmmm.h
Purpose: QM/MM interface

This is a QM/MM interface between one of the forcefields OPLSAA, AMBER and CHARM with MOPAC, GAUSSIAN or DFTB+.
Interactions between QM and MM part are done by electrostatic embedding (see Gerrit Groenhof "Introduction to QM/MM Simulations", figure 4).
Only non-bonded interactions are implemented so there must not be bonds between the QM and the MM part.

MOPAC: Gradients of coulomb interactions between QM and MM part are calculated by CAST using the derived charge distribution for QM atoms from MOPAC.
(see http://openmopac.net/manual/QMMM.html)

GAUSSIAN: Gradients of coulomb interactions between QM and MM part on QM atoms are calculated by GAUSSIAN,
on MM atoms they are calculated by CAST using the electric field from GAUSSIAN. 
(see T. Okamoto et. al., A minimal implementation of the AMBER-GAUSSIAN interface for Ab Initio QM/MM-MD Simulation, DOI 10.1002/jcc.21678)

DFTB+: Gradients of coulomb interactions between QM and MM part are read from DFTB+ outputfile "results.tag".

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
    /**namespaces for QM/MM interface*/
    namespace qmmm
    {
      /**namespace for bonded QM/MM*/
      namespace bonded
      {
        /**struct with all relevant information about a bond*/
        struct Bond
        {
          /**index of MM atom (starting with 0)*/
          int a;
          /**index of QM atom (starting with 0)*/
          int b;
          /**ideal bond length (from force field)*/
          double ideal;
          /**force constant*/
          double force;
          /**constructur
          @param p1: index of one binding partner
          @param p2: index of the other binding partner*/
          Bond(int p1, int p2) { a = p1; b = p2; }
          std::string info()
          {
            return std::to_string(a) + " , "+ std::to_string(b) + " dist: " + std::to_string(ideal) + ", force constant: " + std::to_string(force);
          }
          /**function to calculate force field energy*/
          double calc_energy(coords::Coordinates *cp)
          {
            double E;
            auto const bv(cp->xyz(a) - cp->xyz(b)); // r_ij (i=1, j=2)
            auto const d = len(bv);
            auto const r = d - ideal;
            auto dE = force * r;
            E += dE * r;  // kcal/mol
            return E;
          }
        };

        /**struct with all relevant information about an angle*/
        struct Angle
        {
          /**index of one of the outer atoms (starting with 0)*/
          int a;
          /**index of the other outer atom (starting with 0)*/
          int b;
          /**index of the central atom (starting with 0)*/
          int c;
          /**constructor
          @param p1: index of one of the outer atoms
          @param p2: index of the other outer atom
          @param center: index of the central atom*/
          Angle(int p1, int p2, int center) { a = p1; b = p2; c = center; }
          /**returns all relevant information as a string*/
          std::string info()
          {
            return std::to_string(a) +" , "+ std::to_string(c) + " , " + std::to_string(b);
          }
          /**looks if angle a2 is identical to Angle itself*/
          bool is_equal(Angle a2)
          {
            if (c == a2.c)  // central atom has to be the same
            { // outer atoms can be switched
              if (a == a2.a && b == a2.b) return true;
              else if (a = a2.b && b == a2.a) return true;
            }
            return false;
          }
        };

        /**struct with all relevant information about a dihedral*/
        struct Dihedral
        {
          /**index of one of the outer atoms (starting with 0)*/
          int a;
          /**index of the other outer atom (starting with 0)*/
          int b;
          /**index of the central atom bound to a (starting with 0)*/
          int c1;
          /**index of the central atom bound to b (starting with 0)*/
          int c2;
          /**constructor
          @param p1: index of one of the outer atoms
          @param p2: index of the other outer atom
          @param center1: index of the central atom bound to a
          @param center2: index of the central atom bound to b*/
          Dihedral(int p1, int p2, int center1, int center2)
          {
            a = p1; b = p2; c1 = center1, c2 = center2;
          }
          /**returns all relevant information as a string*/
          std::string info()
          {
            return std::to_string(a) + " , " + std::to_string(c1) + " , "+ std::to_string(c2) + " , " + std::to_string(b);
          }
          /**looks if dihedral d2 is identical to Dihedral itself*/
          bool is_equal(Dihedral d2)
          {
            if (c1 == d2.c1 && c2 == d2.c2 && a == d2.a && b == d2.b) return true;
            else if (c1 == d2.c2 && c2 == d2.c1 && a == d2.b && b == d2.a) return true;
            return false;
          }
        };

        /**looks if vector v contains Angle x
        returns true if yes and false if no*/
        inline bool is_in(Angle x, std::vector<Angle> v)
        {
          for (auto ang : v)
          {
            if (x.is_equal(ang)) return true;
          }
          return false;
        }

        /**looks if vector v contains Dihedral x
        returns true if yes and false if no*/
        inline bool is_in(Dihedral x, std::vector<Dihedral> v)
        {
          for (auto dih : v)
          {
            if (x.is_equal(dih)) return true;
          }
          return false;
        }
      }

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
        std::vector<coords::Cartesian_Point> get_g_coul_mm() const override
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
        /**writes charges inputfile for DFTB+ calculation*/
        void write_dftb_in();
        
        /**function to find bonds, angles and so on between QM and MM system*/
        void find_bonds_etc();
        /**function to find force field parameters for bonds, angles and so on between QM and MM system*/
        void find_parameters();

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

        /**bonds between QM and MM part*/
        std::vector<bonded::Bond> qmmm_bonds;
        /**angles between QM and MM part*/
        std::vector<bonded::Angle> qmmm_angles;
        /**dihedrals between QM and MM part*/
        std::vector<bonded::Dihedral> qmmm_dihedrals;
        
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

        /**information needed to calculate coulomb gradients on MM atoms
        for GAUSSIAN: electric field from gaussian calculation for QM and MM atoms (first QM, then MM)
        only those of the MM atoms are used to calculate the gradients of the electrostatic interaction
        between QM and MM atoms on the MM atoms
        for DFTB+: coulomb gradients on MM atoms due to QM atoms*/
        std::vector<coords::Cartesian_Point> g_coul_mm;
               
        /**checks if all bonds are still intact (bond length smaller than 2.2 Angstrom)*/
        bool check_bond_preservation(void) const;

        /**checks if there is a minimum atom distance (0.3 Angstrom) between atoms*/
        bool check_atom_dist(void) const;

      }; 
    }
  }
}
