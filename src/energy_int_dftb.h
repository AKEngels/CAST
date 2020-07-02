/**
CAST 3
energy_int_dftb.h
Purpose: interface to DFTB+

@author Susanne Sauer
@version 1.0
*/

#pragma once

#include <vector>
#include <string>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <cstdlib>
#include <utility>
#include "helperfunctions.h"
#include "atomic.h"
#include "energy.h"
#include "configuration.h"
#include "coords.h"
#include "coords_io.h"
#include "modify_sk.h"

#if defined (_MSC_VER)
#include "win_inc.h"
#endif

#ifdef _MSC_VER
#pragma warning (disable: 4996)
#endif

namespace energy
{
  namespace interfaces
  {
    namespace dftb
    {
      /**D3 parameters for 3OB, 2OB(base), 2OB(shift) and 2OB(split) in this order
      each element of array has the form {a1, a2, s6, s8}
      values are taken from DFTB+ manual, appendix F*/
      static constexpr std::array<std::array<double, 4>, 4> D3PARAMS = { {
        {0.746, 4.191, 1.0, 3.209},
        {0.717, 2.565, 1.0, 0.011},
        {0.816, 2.057, 1.0, 0.010},
        {0.497, 3.622, 1.0, 0.010},
      } };

      /**interface for DFTB+*/
      class sysCallInterface
        : public energy::interface_base
      {

      public:
        /**constructor: sets optimizer and charge*/
        sysCallInterface(coords::Coordinates*);
        /**delete interface?*/
        ~sysCallInterface(void);

        /*
        Energy class functions that need to be overloaded (for documentation see also energy.h)
        */

        interface_base* clone(coords::Coordinates* coord_object) const;
        interface_base* move(coords::Coordinates* coord_object);

        void swap(interface_base&);
        void swap(sysCallInterface&);

        /** Energy function*/
        double e(void);
        /** Energy+Gradient function*/
        double g(void);
        /** Energy+Hessian function*/
        double h(void);
        /** Optimization in the interface(d program)*/
        double o(void);
        /**prints total energy*/
        void print_E(std::ostream&) const;
        /**prints 'headline' for energies*/
        void print_E_head(std::ostream&, bool const endline = true) const;
        /**prints partial energies (not much sense in it because no partial energies are read)*/
        void print_E_short(std::ostream&, bool const endline = true) const;
        /**does nothing*/
        void to_stream(std::ostream&) const;
        // "update" function
        void update(bool const) { }
        /**returns partial atomic charges*/
        std::vector<coords::float_type> charges() const override { return partial_charges; };
        /**returns gradients on external charges due to the molecular system (used for QM/MM)*/
        coords::Gradients_3D get_g_ext_chg() const override { return grad_ext_charges; };

      private:

        /**constructor for clone and move functions*/
        sysCallInterface(sysCallInterface const& rhs, coords::Coordinates* cobj);

        /**writes dftb+ inputfile
        @param t: type of calculation (0 = energy, 1 = gradient, 2 = hessian, 3 = optimize)*/
        void write_inputfile(int t);

        /**reads dftb+ outputfile (results.tag)
        @param t: type of calculation (0 = energy, 1 = gradient, 2 = hessian, 3 = optimize)*/
        double read_output(int t);

        /**total energy*/
        double energy;

        /**partial atomic charges*/
        std::vector<coords::float_type> partial_charges;

        /**gradients of external charges*/
        coords::Gradients_3D grad_ext_charges;
      };
    }
  }
}
