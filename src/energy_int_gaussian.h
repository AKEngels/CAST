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
        void print_G_tinkerlike(std::ostream&, bool const aggregate = false) const;
        void to_stream(std::ostream&) const;

        /**return charges*/
        std::vector<coords::float_type> charges() const override;
        /**returns electric field that is given by GAUSSIAN
        (is used for QM/MM)*/
        std::vector<coords::Cartesian_Point> get_g_coul_mm() const override;
        /**return id for gaussian call*/
        std::string get_id() const override { return id;};
        /**overwritten function, should not be called*/
        coords::Gradients_3D get_link_atom_grad() const override
        {
          throw std::runtime_error("function not implemented\n");
        }

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

        /**electric field that is given by GAUSSIAN
        (is used for QM/MM)*/
        std::vector<coords::Cartesian_Point> electric_field;

        /*
        Gaussian sysCall funcntions
        */

        int callGaussian(void);
        void print_gaussianInput(char);
        void read_gaussianOutput(bool const grad = true, bool const opt = true, bool const qmmm=false);
        void removeTempFiles(void);

        /*
        Helper functions
        */

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