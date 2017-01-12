/*
* GAUSSIAN Interface
*/

#pragma once 

#include <vector>
#include <string>

#include "energy.h"


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

      private:

        //constructor for clone and move functions
        sysCallInterfaceGauss(sysCallInterfaceGauss const & rhs, coords::Coordinates *cobj);

        // heat of formation
				double hof_kcal_mol, hof_kj_mol;
				// energies
				double e_total, e_electron, e_core;
				std::string id;

				// FAILCOUNTER
				size_t failcounter;

        /*
        Gaussian sysCall funcntions
        */

        int callGaussian(void);
        void print_gaussianInput();
        void read_gaussianOutput();
        void removeTempFiles(void);

        /*
        Helper functions
        */

        bool check_bond_preservation(void) const;

      };

    }
  }
}