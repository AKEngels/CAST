/*
* GAUSSIAN Interface
*/

#pragma once

#include <vector>
#include <string>
#include "energy.h"
#include "helperfunctions.h"
#include "modify_sk.h"


namespace energy
{
  namespace interfaces
  {
    namespace gaussian
    {

      class sysCallInterfaceGauss
        : public energy::interface_base
      {

      public:

        sysCallInterfaceGauss(coords::Coordinates*);
        ~sysCallInterfaceGauss(void);

        /*
        Override of virtual voids derrived from parent class
        */

        void swap(interface_base&);
        void swap(sysCallInterfaceGauss&);

        interface_base * clone(coords::Coordinates * coord_object) const;
        interface_base * move(coords::Coordinates * coord_object);

        void update(bool const) { }

        // Energy function
        double e(void);
        // Energy+Gradient function
        double g(void);
        // Energy+Gradient+Hessian function
        double h(void);
        // Optimization in the interface(d program)
        double o(void);

        // Output functions
        void print_E(std::ostream&) const;
        void print_E_head(std::ostream&, bool const endline = true) const;
        void print_E_short(std::ostream&, bool const endline = true) const;
        void to_stream(std::ostream&) const;

        /**return charges*/
        std::vector<coords::float_type> charges() const override;
        /**returns gradients on external charges (calculated by electric field from GAUSSIAN)
        (is used for QM/MM)*/
        std::vector<coords::Cartesian_Point> get_g_ext_chg() const override;

      private:

        ////MO, excitation energies and dipolemoments
        //std::vector <double> occMO, virtMO, excitE;
        //coords::Representation_3D  ex_ex_trans, gz_ex_trans;
        //std::vector <int> state_i, state_j, gz_i_state;

        //constructor for clone and move functions
        sysCallInterfaceGauss(sysCallInterfaceGauss const & rhs, coords::Coordinates *cobj);

        // heat of formation
				double hof_kcal_mol, hof_kj_mol;

				// energies
				double e_total, e_electron, e_core;

				std::string id;

				// FAILCOUNTER
				size_t failcounter;

        /**mulliken charges*/
        std::vector<double> atom_charges;
        /**electric field that is given by GAUSSIAN
				contains the atoms and the external charges
        (is used for QM/MM)*/
        std::vector<coords::Cartesian_Point> electric_field;

        /**gradients of link atoms*/
        coords::Gradients_3D link_atom_grad;

        /*
        Gaussian sysCall funcntions
        */

        int callGaussian(void);
        void print_gaussianInput(char);
        void read_gaussianOutput(bool const grad = true, bool const opt = true, bool const qmmm=false);
        void removeTempFiles(void);

        /**checks if all bonds are still intact (bond length smaller than 1.2 sum of covalent radii)*/
        bool check_bond_preservation(void) const;

        /*! Pointer to the single instance of the gaussian interface class
        *
        * A pointer to it is contained here.
        * If no object exists (yet), this will be a nullpointer.
        */
        /*static sysCallInterfaceGauss *m_interface;*/

      };

    }
  }
}
