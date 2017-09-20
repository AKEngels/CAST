#ifdef USE_PYTHON

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
#include <Python.h>

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
			/*
			function that returns the path to a pythonmodule
			(path has to be appended to pythonpath if you want to call this module)
			@param modulename: name of the module
			*/
			std::string get_python_modulepath(std::string modulename);
            
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
				~sysCallInterface(void);

				/*
				Energy class functions that need to be overloaded
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
				/**prints gradients*/
				void print_G_tinkerlike(std::ostream&, bool const aggregate = false) const;
				/**does nothing*/
				void to_stream(std::ostream&) const;
				// "update" function
				void update(bool const) { }

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

				/*
				Helper functions
				*/
				bool check_bond_preservation(void) const;

			};

		}
	}
}
#endif