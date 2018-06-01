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


			class sysCallInterface
				: public energy::interface_base
			{

			public:
        /**constructor: sets optimizer*/
				sysCallInterface(coords::Coordinates*);
         /**delete interface?*/
				~sysCallInterface(void);

				/*
				Energy class functions that need to be overloaded (for documentation see also energy.h)
				*/

				interface_base * clone(coords::Coordinates * coord_object) const;
				interface_base * move(coords::Coordinates * coord_object);

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
				/**prints partial energies (not much sense in it because not partial energies are read)*/
				void print_E_short(std::ostream&, bool const endline = true) const;
				/**does nothing*/
				void to_stream(std::ostream&) const;
				// "update" function
				void update(bool const) { }
        /**returns partial atomic charges*/
        std::vector<coords::float_type> charges() const override;
        /**returns gradients on external charges due to the molecular system (used for QM/MM)*/
        std::vector<coords::Cartesian_Point> get_g_ext_chg() const override;
        /**overwritten function, should not be called*/
        std::string get_id() const override
        {
          throw std::runtime_error("Function not implemented.\n");
        }
        /**get gradients on link atoms (for QM/MM)*/
        coords::Gradients_3D get_link_atom_grad() const override;

			private:

				/**constructor for clone and move functions*/
				sysCallInterface(sysCallInterface const & rhs, coords::Coordinates *cobj);

        /**checks if structure is complete, i.e. no coordinates are NaN
        coordinates become NaN sometimes in TS (dimer method)*/
        bool check_structure();

        /**writes dftb+ inputfile
        @param t: type of calculation (0 = energy, 1 = gradient, 2 = hessian, 3 = optimize)*/
        void write_inputfile(int t);

        /**reads dftb+ outputfile (results.tag)
        @param t: type of calculation (0 = energy, 1 = gradient, 2 = hessian, 3 = optimize)*/
        double read_output(int t);

        /**total energy*/
				double energy;

        /**gradients of external charges*/
        std::vector<coords::Cartesian_Point> grad_ext_charges;

        /**gradients of link atoms*/
        coords::Gradients_3D link_atom_grad;

        /**checks if all bonds are still intact (bond length smaller than 1.2 sum of covalent radii)*/
        bool check_bond_preservation(void) const;

        /**checks if there is a minimum atom distance (0.3 Angstrom) between atoms*/
        bool check_atom_dist(void) const;

			};

		}
	}
}
