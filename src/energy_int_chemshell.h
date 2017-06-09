#ifndef ENERGY_INT_CHEMSHELL_H
#define ENERGY_INT_CHEMSHELL_H

#include<memory>
#include<stdlib.h>
#include<string>
#include<fstream>
#include<sstream>
#include<istream>
#include<iostream>
#include"coords_io.h"
#include"energy.h"
#include"coords.h"

namespace energy {
	namespace interfaces {
		namespace chemshell {

			class sysCallInterface final: public interface_base{
			public:
				sysCallInterface(coords::Coordinates * coord_ptr)
					: interface_base(coord_ptr), tmp_file_name(Config::get().general.outputFilename){
					std::stringstream ss;
					std::srand(std::time(0));
					ss << (std::size_t(std::rand()) | (std::size_t(std::rand()) << 15));
					tmp_file_name.append("_tmp_").append(ss.str());
				}
				~sysCallInterface() final {
					if (Config::get().energy.gaussian.delete_input)
					{
						std::string rem_xyz(tmp_file_name);
						std::string rem_pdb(tmp_file_name);

						rem_xyz.append(".xyz");
						rem_pdb.append(".pdb");

						remove(rem_xyz.c_str());
						remove(rem_pdb.c_str());
					}
				};
				sysCallInterface(sysCallInterface const & other) = default;
				sysCallInterface(sysCallInterface const & other, coords::Coordinates * coord)
					: interface_base(coord) {
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

				void create_pdb() const;
				void write_xyz(std::string const & os) const;
				void write_input() const;
				void write_chemshell_file(std::string const & o_file) const;
				void call_tleap()const;
				void make_tleap_input(std::string const & o_file)const;
				void call_chemshell() const;
				void actual_call()const;



			};
		}
	}
}

#endif