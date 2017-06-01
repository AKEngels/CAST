#ifndef ENERGY_INT_CHEMSHELL_H
#define ENERGY_INT_CHEMSHELL_H

#include"energy.h"
#include"coords.h"
#include<memory>
#include<cstdio>
#include<string>
#include<fstream>
#include<sstream>
#include<istream>
#include<iostream>

namespace energy {
	namespace interfaces {
		namespace chemshell {

			class sysCallInterface : public interface_base{
			public:
				sysCallInterface(coords::Coordinates * coord_ptr)
					: interface_base(coord_ptr){
					tmp_file_name = std::tmpnam(nullptr);
				}
				~sysCallInterface() override {};
				sysCallInterface(sysCallInterface const & other)
					: interface_base(other.coords) {
					interface_base::operator=(other);
				}

			private:
				std::string tmp_file_name;

				std::string create_pdb();
				void write_input(std::ostream & os);

			};
		}
	}
}

#endif