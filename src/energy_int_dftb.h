/*DW                            *
MOPAC Interface
*                            DW*/

#pragma once 

#include <vector>
#include <string>
//#include "helperfunctions.h"
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
//#include "python2.7/Python.h"

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

			class sysCallInterface
				: public energy::interface_base
			{

			public:

				sysCallInterface(coords::Coordinates*);
				~sysCallInterface(void);

				/*
				Energy class functions that need to be overloaded
				*/

				interface_base * clone(coords::Coordinates * coord_object) const;
				interface_base * move(coords::Coordinates * coord_object);

				void swap(interface_base&);
				void swap(sysCallInterface&);

				// Energy function
				double e(void);
				// Energy+Gradient function
				double g(void);
				// Energy+Gradient+Hessian function
				double h(void);
				// Optimization in the interface(d program)
				double o(void);
				//MOPAC7_HB VAR
				bool grad_var;
				// Output functions
				void print_E(std::ostream&) const;
				void print_E_head(std::ostream&, bool const endline = true) const;
				void print_E_short(std::ostream&, bool const endline = true) const;
				void print_G_tinkerlike(std::ostream&, bool const aggregate = false) const;
				void to_stream(std::ostream&) const;
				// "update" function
				void update(bool const) { }

			private:

				// constructor for clone and move functions
				sysCallInterface(sysCallInterface const & rhs, coords::Coordinates *cobj);

				// energies
				double e_bs, e_coul, e_rep, e_lr, e_tot;
				
				std::string id;

				// FAILCOUNTER
				size_t failcounter;

				/*
				Helper functions
				*/
				bool check_bond_preservation(void) const;

			};

		}
	}
}
