/**
CAST 3
energy_int_qmmm_a.h
Purpose: QM/MM interface


This is a QM/MM interface between one of the forcefields OPLSAA and AMBER with MOPAC, GAUSSIAN, DFTB+, Psi4 or ORCA.
Interactions between QM and MM part are done by electrostatic embedding (see Gerrit Groenhof "Introduction to QM/MM Simulations", figure 4).

MOPAC: Gradients of coulomb interactions between QM and MM part are calculated by CAST using the derived charge distribution for QM atoms from MOPAC.
(see http://openmopac.net/manual/QMMM.html)

GAUSSIAN and PSI4: Gradients of coulomb interactions between QM and MM part on QM atoms are calculated by GAUSSIAN,
on MM atoms they are calculated by CAST using the electric field from GAUSSIAN.
(see T. Okamoto et. al., A minimal implementation of the AMBER-GAUSSIAN interface for Ab Initio QM/MM-MD Simulation, DOI 10.1002/jcc.21678)

DFTB+ and ORCA: Gradients of coulomb interactions between QM and MM part are directly read from outputfiles of QM program.

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
#include "qmmm_helperfunctions.h"

namespace energy
{
  namespace interfaces
  {
    /**namespace for QM/MM interfaces*/
    namespace qmmm
    {
      /**namespace for bonded QM/MM*/
      namespace bonded
      {
        /**struct with all relevant information about a bond*/
        struct Bond
        {
          /**index of QM atom (starting with 0)*/
          unsigned int a;
          /**index of MM atom (starting with 0)*/
          unsigned int b;
          /**ideal bond length (from force field)*/
          double ideal;
          /**force constant*/
          double force;
          /**constructur
          @param p1: index of one binding partner
          @param p2: index of the other binding partner*/
          Bond(int const p1, int const p2) { a = p1; b = p2; }
          std::string info() const
          {
            return std::to_string(a + 1) + " , " + std::to_string(b + 1) + " dist: " + std::to_string(ideal) + ", force constant: " + std::to_string(force);
          }
          /**function to calculate force field energy and gradients (copied from energy_int_aco.cc)
          @param cp: pointer to coordinates object
          @param gradients: gradients that are filled during gradient calculation (necessary but not used also in only-energy calculation)
          @param grad: true if gradients shold be calculated*/
          double calc_energy(coords::Coordinates* cp, coords::Gradients_3D& gradients, bool const grad = false)
          {
            double E(0.0);
            auto const bv(cp->xyz(a) - cp->xyz(b)); // r_ij (i=1, j=2)
            auto const d = len(bv);
            auto const r = d - ideal;
            auto dE = force * r;
            E += dE * r;  // kcal/mol
            if (grad == true)
            {
              dE *= 2;  // kcal/(mol*Angstrom)  gradient without direction
              dE /= d;  // kcal/(mol*A^2)   gradient divided by distance because later it is multiplied with it again
              auto const gv = bv * dE;   // "force" on atom i due to atom j (kcal/(mol*A)), gradient with direction
              gradients[a] += gv;   //(kcal/(mol*A))
              gradients[b] -= gv;
            }
            return E;
          }
        };

        /**struct with all relevant information about an angle*/
        struct Angle
        {
          /**index of one of the outer atoms (starting with 0)*/
          unsigned int a;
          /**index of the other outer atom (starting with 0)*/
          unsigned int b;
          /**index of the central atom (starting with 0)*/
          unsigned int c;
          /**ideal angle (from force field)*/
          double ideal;
          /**force constant*/
          double force;
          /**constructor
          @param p1: index of one of the outer atoms
          @param p2: index of the other outer atom
          @param center: index of the central atom*/
          Angle(int const p1, int const p2, int const center) { a = p1; b = p2; c = center; }
          /**returns all relevant information as a string*/
          std::string info() const
          {
            return std::to_string(a + 1) + " , " + std::to_string(c + 1) + " , " + std::to_string(b + 1) + " angle: " + std::to_string(ideal) + ", force constant: " + std::to_string(force);
          }
          /**function to calculate force field energy and gradients (copied from energy_int_aco.cc)
          @param cp: pointer to coordinates object
          @param gradients: gradients that are filled during gradient calculation (necessary but not used also in only-energy calculation)
          @param grad: true if gradients shold be calculated*/
          double calc_energy(coords::Coordinates* cp, coords::Gradients_3D& gradients, bool grad = false)
          {
            double E(0.0);
            auto av1(cp->xyz(a) - cp->xyz(c));
            auto av2(cp->xyz(b) - cp->xyz(c));
            auto const d(scon::angle(av1, av2).degrees() - ideal);
            auto const r(d * SCON_PI180);
            E += force * r * r;
            if (grad == true)
            {
              auto dE = force * r;
              coords::Cartesian_Point const cv(cross(av1, av2));
              coords::float_type const cvl(len(cv));
              dE *= 2.0 / cvl;
              coords::Cartesian_Point const gv1(cross(av1, cv) * (dE / dot(av1, av1)));
              coords::Cartesian_Point const gv2(cross(cv, av2) * (dE / dot(av2, av2)));
              gradients[a] += gv1;
              gradients[b] += gv2;
              gradients[c] += -(gv1 + gv2);
            }
            return E;
          }
        };
        // overloaded == operator for Angle
        inline bool operator==(Angle const& lhs, Angle const& rhs)
        {
          if (lhs.c == rhs.c)  // central atom has to be the same
          { // outer atoms can be switched
            if (lhs.a == rhs.a && rhs.b == lhs.b) return true;
            else if (lhs.a == rhs.b && lhs.b == rhs.a) return true;
          }
          return false;
        }

        /**struct with all relevant information about a dihedral*/
        struct Dihedral
        {
          /**index of one of the outer atoms (starting with 0)*/
          unsigned int a;
          /**index of the other outer atom (starting with 0)*/
          unsigned int b;
          /**index of the central atom bound to a (starting with 0)*/
          unsigned int c1;
          /**index of the central atom bound to b (starting with 0)*/
          unsigned int c2;
          /**parameter for force field*/
          int max_order;
          /**parameter for force field*/
          int number;
          /**parameter for force field*/
          std::array<double, 4> forces;
          /**parameter for force field*/
          std::array<double, 4> ideals;
          /**parameter for force field*/
          std::array<size_t, 4> orders;
          /**constructor
          @param p1: index of one of the outer atoms
          @param p2: index of the other outer atom
          @param center1: index of the central atom bound to a
          @param center2: index of the central atom bound to b*/
          Dihedral(int const p1, int const p2, int const center1, int const center2)
          {
            a = p1; b = p2; c1 = center1, c2 = center2;
          }
          /**returns all relevant information as a string*/
          std::string info() const
          {
            std::string result = "Atoms: " + std::to_string(a + 1) + " , " + std::to_string(c1 + 1) + " , " + std::to_string(c2 + 1) + " , " + std::to_string(b + 1) + "\n";
            result += "  max order: " + std::to_string(max_order) + ", number: " + std::to_string(number) + "\n";
            result += "  orders: ";
            for (auto o : orders) result += std::to_string(o) + "  ";
            result += "\n  force constants: ";
            for (auto f : forces) result += std::to_string(f) + "  ";
            result += "\n  ideals: ";
            for (auto i : ideals) result += std::to_string(i) + "  ";
            return result;
          }
          
          /**function to calculate force field energy and gradients (copied from energy_int_aco.cc)
          @param cp: pointer to coordinates object
          @param gradients: gradients that are filled during gradient calculation (necessary but not used also in only-energy calculation)
          @param grad: true if gradients shold be calculated*/
          double calc_energy(coords::Coordinates* cp, double const torsionunit, coords::Gradients_3D& gradients, bool const grad = false)
          {
            double E(0.0);
            // Get bonding vectors
            coords::Cartesian_Point const b01 = cp->xyz(c1) - cp->xyz(a);
            coords::Cartesian_Point const b12 = cp->xyz(c2) - cp->xyz(c1);
            coords::Cartesian_Point const b23 = cp->xyz(b) - cp->xyz(c2);
            // Cross terms
            coords::Cartesian_Point const t = cross(b01, b12);
            coords::Cartesian_Point const u = cross(b12, b23);
            // Get length and variations
            coords::float_type const tl2 = dot(t, t);
            coords::float_type const ul2 = dot(u, u);
            // ...
            coords::float_type const tlul = sqrt(tl2 * ul2);
            coords::float_type const r12 = len(b12);
            // cross of cross
            coords::Cartesian_Point const tu = cross(t, u);
            // scalar and length variations
            coords::float_type const cos_scalar0 = dot(t, u);
            coords::float_type const cos_scalar1 = tlul;
            coords::float_type const sin_scalar0 = dot(b12, tu);
            coords::float_type const sin_scalar1 = r12 * tlul;
            // Get multiple sine and cosine values
            coords::float_type cos[7], sin[7];
            cos[1] = cos_scalar0 / cos_scalar1;
            sin[1] = sin_scalar0 / sin_scalar1;
            for (auto j{ 2 }; j <= max_order; ++j)
            {
              std::size_t const k = j - 1;
              sin[j] = sin[k] * cos[1] + cos[k] * sin[1];
              cos[j] = cos[k] * cos[1] - sin[k] * sin[1];
            }

            coords::float_type tE(0.0), dE(0.0);
            for (auto j{ 0 }; j < number; ++j)
            {
              coords::float_type const F = forces[j] * torsionunit;
              std::size_t const k = orders[j];
              coords::float_type const l = std::abs(ideals[j]) > 0.0 ? -1.0 : 1.0;
              tE += F * (1.0 + cos[k] * l);
              dE += -static_cast<coords::float_type>(k) * F * sin[k] * l;
            }
            E += tE;
            if (grad == true)
            {
              coords::Cartesian_Point const b02 = cp->xyz(c2) - cp->xyz(a);
              coords::Cartesian_Point const b13 = cp->xyz(b) - cp->xyz(c1);

              coords::Cartesian_Point const dt(cross(t, b12) * (dE / (tl2 * r12)));
              coords::Cartesian_Point const du(cross(u, b12) * (-dE / (ul2 * r12)));

              coords::Cartesian_Point const vir1 = cross(dt, b12);
              coords::Cartesian_Point const vir2 = cross(b02, dt) + cross(du, b23);
              coords::Cartesian_Point const vir3 = cross(dt, b01) + cross(b13, du);
              coords::Cartesian_Point const vir4 = cross(du, b12);

              gradients[a] += vir1;
              gradients[c1] += vir2;
              gradients[c2] += vir3;
              gradients[b] += vir4;
            }
            return E;
          }
        };
        // overloaded == operator for Dihedral
        inline bool operator==(Dihedral const& lhs, Dihedral const& rhs)
        {
          if (lhs.c1 == rhs.c1 && lhs.c2 == rhs.c2 && lhs.a == rhs.a && lhs.b == rhs.b) return true;
          else if (lhs.c1 == rhs.c2 && lhs.c2 == rhs.c1 && lhs.a == rhs.b && lhs.b == rhs.a) return true;
          return false;
        }
      }

      /**QMMM interface class*/
      class QMMM_A
        : public interface_base
      {

        /**uncontracted forcefield parameters*/
        static ::tinker::parameter::parameters tp;
        /**contracted forcefield parameters*/
        ::tinker::parameter::parameters cparams;

      public:

        /**Constructor*/
        QMMM_A(coords::Coordinates*);
        /**overloaded Constructor*/
        QMMM_A(QMMM_A const&, coords::Coordinates*);
        /**another overload of Constructor*/
        QMMM_A(QMMM_A&&, coords::Coordinates*);

        /*
        Energy class functions that need to be overloaded (for documentation see also energy.h)
        */

        interface_base* clone(coords::Coordinates*) const;
        interface_base* move(coords::Coordinates*);

        void swap(interface_base&);
        void swap(QMMM_A&);

        /** update structure (account for topology or rep change)*/
        void update(bool const skip_topology = false);

        /**sets the atom coordinates of the subsystems (QM and MM) to those of the whole coordobject*/
        void update_representation();

        /** Energy function*/
        coords::float_type e() override;
        /** Energy+Gradient function */
        coords::float_type g() override;
        /** Energy+Hessian function*/
        coords::float_type h() override;
        /** Optimization in the interface or interfaced program (not existent for this interface)*/
        coords::float_type o() override;

        /** Return charges (for QM und MM atoms) */
        std::vector<coords::float_type> charges() const override;
        /**overwritten function, should not be called*/
        coords::Gradients_3D get_g_ext_chg() const override {
          throw std::runtime_error("function not implemented\n");
        }

        /**prints total energy (not implemented)*/
        void print_E(std::ostream&) const final override;
        /**prints 'headline' for energies*/
        void print_E_head(std::ostream&,
          bool const endline = true) const final override;
        /**prints partial energies*/
        void print_E_short(std::ostream&,
          bool const endline = true) const final override;
        /**function not implemented*/
        void to_stream(std::ostream&) const override;


      private:

        /**function where QM/MM calculation is prepared*/
        void prepare_bonded_qmmm();
        /**function to find bonds, angles and so on between QM and MM system*/
        void find_bonds_etc();
        /**function to find force field parameters for bonds, angles and so on between QM and MM system*/
        void find_parameters();
        /**determines if a van der waals interaction between a QM and a MM atom should be calculated
        @param qm: index of QM atom
        @param mm: index of MM atom
        returns 0 if no vdW is calculated (1 or 2 bonds between the atoms), 1 if vdW is calculated normally
        and 2 if vdW is scaled down by 1/2 (3 bonds between the atoms)*/
        int calc_vdw(unsigned const qm, unsigned const mm) const;

        /**calculates interaction between QM and MM part
        energy is only vdW interactions, gradients are coulomb and vdW
        @param if_gradient: true if gradients should be calculated, false if not*/
        void ww_calc(bool const if_gradient);
        /**calculates energies and gradients
        @param if_gradient: true if gradients should be calculated, false if not*/
        coords::float_type qmmm_calc(bool const if_gradient);
        /**calculates bonded energy and gradients
        @param if_gradient: true if gradients should be calculated, false if not*/
        double calc_bonded(bool const if_gradient);

        /**indizes of QM atoms*/
        std::vector<size_t> qm_indices;
        /**indizes of MM atoms*/
        std::vector<size_t> mm_indices;
        /**indizes of MM atoms that are taken into acoount for electrostatic interaction with QM region*/
        std::vector<int> charge_indices;

        /**vector of length total number of atoms
        only those elements are filled whose position corresponds to QM atoms
        they are filled with successive numbers starting from 0
        purpose: faciliate mapping between total coordinates object and subsystems*/
        std::vector<size_t> new_indices_qm;
        /**vector of length total number of atoms
        only those elements are filled whose position corresponds to MM atoms
        they are filled with successive numbers starting from 0
        purpose: faciliate mapping between total coordinates object and subsystems*/
        std::vector<size_t> new_indices_mm;

        /**link atoms*/
        std::vector<LinkAtom> link_atoms;

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
        /**some parameter needed to calculate dihedral energy*/
        double torsionunit;

        /**atom index that determines center of QM region*/
        std::size_t index_of_QM_center;

        /**energy of only QM system*/
        coords::float_type qm_energy;
        /**energy of only MM system*/
        coords::float_type mm_energy;
        /**van der Waals interaction energy between QM and MM atoms*/
        coords::float_type vdw_energy;
        /**energy of bonded interactions between QM and MM atoms (bonds, angles and dihedrals)*/
        coords::float_type bonded_energy;
        /**coulomb energy between QM and MM system (only used for mechanical embedding)*/
        coords::float_type coulomb_energy;

        /**gradients of electrostatic interaction between QM and MM atoms
        electrostatic embedding: only gradients on MM atoms (those on QM atoms are calculated by QM interface)
        mechanical embedding: gradients on QM and on MM atoms*/
        coords::Gradients_3D coulomb_gradient;
        /**gradients of van der waals interaction energy between QM and MM atoms*/
        coords::Gradients_3D vdw_gradient;
        /**gradients of bonded interactions energy between QM and MM atoms*/
        coords::Gradients_3D bonded_gradient;

        /**gradients on external charges due to QM atoms 
        (only used for electrostatic embedding)
        from this the variable coulomb_gradient will be filled*/
        coords::Gradients_3D g_coul_mm;

        /**total charges (only used for mechanical embedding)*/
        std::vector<double> total_atom_charges;
      };
    }
  }
}
