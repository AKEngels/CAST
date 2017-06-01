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

			class sysCallInterface final: public interface_base{
			public:
				sysCallInterface(coords::Coordinates * coord_ptr)
					: interface_base(coord_ptr){
					tmp_file_name = std::tmpnam(nullptr);
				}
				~sysCallInterface() final {};
				sysCallInterface(sysCallInterface const & other)
					: interface_base(other.coords) {
					interface_base::operator=(other);
				}

				void swap(interface_base & other) final;
				interface_base * clone(coords::Coordinates * coord_object) const final;
				interface_base * move(coords::Coordinates * coord_object) final;

				void update(bool const skip_topology = false) final;

				coords::float_type e(void) final;
				coords::float_type g(void) final;
				coords::float_type h(void) final;
				coords::float_type o(void) final;

				void print_E(std::ostream&) const final;
				void print_E_head(std::ostream&, bool const endline = true) const final;
				void print_E_short(std::ostream&, bool const endline = true) const final;
				void print_G_tinkerlike(std::ostream&, bool const aggregate = false) const final;
				void to_stream(std::ostream&) const final;

			private:
				std::string tmp_file_name;

				std::string create_pdb();
				void write_input(std::ostream & os);



			};
		}
	}
}

#endif