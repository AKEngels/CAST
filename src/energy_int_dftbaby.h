#ifdef USE_PYTHON

#pragma once

#include <Python.h>
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
		namespace dftbaby
		{
			/*
			function that creates a string that can be run as a python programme
			this programme sets the correct pythonpath
			@param numpath: path to module numpy
			@param scipath: path to param scipy
			*/
			std::string create_pythonpath(std::string numpath, std::string scipath);

			/*
			function that creates a configfile for dftbaby (dftbaby.cfg)
			out of configuration options in CAST inputfile
			*/
			void create_dftbaby_configfile();

			class sysCallInterface
				: public energy::interface_base
			{

			public:
                /**constructor
				when called first, sets all partial energies to 0,
				fills the string add_path for adding stuff to pythonpath
				and creates configuration file for dftbaby*/
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
				// Output functions
				/**prints total energy*/
				void print_E(std::ostream&) const;
				/**prints 'headline' for energies*/
				void print_E_head(std::ostream&, bool const endline = true) const;
				/**prints partial energies*/
				void print_E_short(std::ostream&, bool const endline = true) const;
				/**does nothing*/
				void to_stream(std::ostream&) const;
				// "update" function
				void update(bool const) { }
        /**returns partial atomic charges*/
        std::vector<coords::float_type> charges() const override;
        /**overwritten function, should not be called*/
        std::vector<coords::Cartesian_Point> get_g_coul_mm() const override
        {
          throw std::runtime_error("TODO: Implement electric field.\n");
        }
        /**overwritten function, should not be called*/
        coords::Gradients_3D get_link_atom_grad() const override
        {
          throw std::runtime_error("function not implemented\n");
        }
        /**overwritten function*/
        std::string get_id() const override { return "bullshit"; }

			private:

				// constructor for clone and move functions
				sysCallInterface(sysCallInterface const & rhs, coords::Coordinates *cobj);

				// energies
				/**band structure energy*/
				double e_bs;
				/**coulomb energy*/
				double e_coul;
				/**long range correction*/
				double e_lr;
				/**repulsion energy (nuclear)*/
				double e_rep;
				/**total energy*/
				double e_tot;

				/**string that contains pythonprogramme to
				add all necessary paths to pythonpath*/
				std::string add_path;

        /**checks if all bonds are still intact (bond length smaller than 1.2 sum of covalent radii)*/
				bool check_bond_preservation(void) const;

			};

		}
	}
}
#endif
