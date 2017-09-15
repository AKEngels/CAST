#ifndef ENERGY_INT_CHEMSHELL_H
#define ENERGY_INT_CHEMSHELL_H

#include<memory>
#include<stdlib.h>
#include<string>
#include<fstream>
#include<sstream>
#include<istream>
#include<iostream>
#include<unordered_set>
#include"coords_io.h"
#include"energy.h"
#include"coords.h"

namespace energy {
	namespace interfaces {
		namespace chemshell {

			class sysCallInterface final: public interface_base{
			public:
				/*
				template<typename T>
				using remove_cr = typename std::remove_const<typename std::remove_reference<T>::type>::type;

				template<typename T>
				struct get_spherical_types {
				};

				template<typename T, typename U>
				struct get_spherical_types<scon::sphericals<T, U>> {
					using length_type = T;
					using angle_type = U;
				};

				using length_type = typename get_spherical_types<remove_cr<decltype(coords->intern(0))>>::length_type;
				using angle_type = typename get_spherical_types<remove_cr<decltype(coords->intern(0))>>::angle_type;
				*/
				sysCallInterface(coords::Coordinates * coord_ptr)
					: interface_base(coord_ptr), 
					tmp_file_name(Config::get().general.outputFilename)
				{
					optimizer = true;
					std::stringstream ss;
					std::srand(std::time(0));
					ss << (std::size_t(std::rand()) | (std::size_t(std::rand()) << 15));
					tmp_file_name.append("_tmp_").append(ss.str());
					first_call = true;
					
				}
				~sysCallInterface() final {
					if (Config::get().energy.chemshell.delete_input)
					{
						get_rid_of_dump_files();
					}
				};
				/*sysCallInterface(sysCallInterface const & other) {
					optimizer = true;
					std::stringstream ss;
					std::srand(std::time(0));
					ss << (std::size_t(std::rand()) | (std::size_t(std::rand()) << 15));
					tmp_file_name.append("_tmp_").append(ss.str());
					first_call = other.first_call;

					change_input_file_names(other.tmp_file_name);
				};

				sysCallInterface(sysCallInterface && other) = default;*/

				sysCallInterface(sysCallInterface const & other, coords::Coordinates * coord)
					: interface_base(coord),
					tmp_file_name(Config::get().general.outputFilename) {
					interface_base::operator=(other);
					optimizer = true;
					std::stringstream ss;
					std::srand(std::time(0));
					ss << (std::size_t(std::rand()) | (std::size_t(std::rand()) << 15));
					tmp_file_name.append("_tmp_").append(ss.str());
					first_call = other.first_call;

					change_input_file_names(other.tmp_file_name);
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

				//TODO: Get better value from better source
				static auto constexpr au_to_kcalmol = 627.503;

			private:
				std::string tmp_file_name;
				bool first_call;

				void create_pdb() const;
				void write_xyz(std::string const & os) const;
				void write_chemshell_coords()const;
				void write_chemshell_file(bool const & sp=true) const;
				void call_tleap()const;
				void make_tleap_input(std::string const & o_file)const;
				void call_chemshell(bool single_point = true) const;
				void actual_call()const;
				std::pair<std::string, std::string> find_active_and_inactive_atoms(std::string const & qm_atoms)const;
				void read_gradients();

				static std::string trim_space_and_tabs(std::string const & str);

				bool check_if_number(std::string const & number)const;

				coords::float_type read_energy()const;
				void read_coords();

				bool check_if_line_is_coord(std::vector<std::string> const & coords)const;
				coords::Cartesian_Point make_coords(std::vector<std::string> const & line)const;
				void change_name_of_energy_and_grad()const;
				coords::Representation_3D extract_gradients(std::vector<coords::float_type> const & grads) const;

				void change_input_file_names(std::string const & filename, std::string const & copy_or_move = "cp")const;

				void make_sp()const;
				void make_opti()const;
				void initialize_before_first_use()const;
				inline void check_for_first_call() {
					if (first_call) {
						first_call = false;
						initialize_before_first_use();
					}
				}

				void get_rid_of_dump_files()const {

					std::vector<std::string> dump_files;

					//Just all the files which are to dump
					dump_files.emplace_back(tmp_file_name + ".frcmod");
					dump_files.emplace_back(tmp_file_name + ".in");
					dump_files.emplace_back(tmp_file_name + ".out");
					dump_files.emplace_back(tmp_file_name + ".inpcrd");
					dump_files.emplace_back(tmp_file_name + ".mol2");
					dump_files.emplace_back(tmp_file_name + ".pdb");
					dump_files.emplace_back(tmp_file_name + ".xyz");
					dump_files.emplace_back(tmp_file_name + ".prmtop");
					dump_files.emplace_back(tmp_file_name + ".chm");
					dump_files.emplace_back("leap.log");
					dump_files.emplace_back("ANTECHAMBER_AC.AC");
					dump_files.emplace_back("ANTECHAMBER_AC.AC0");
					dump_files.emplace_back("ANTECHAMBER_BOND_TYPE.AC");
					dump_files.emplace_back("ANTECHAMBER_BOND_TYPE.AC0");
					dump_files.emplace_back("ATOMTYPE.INF");

					for (auto const & dump : dump_files) {
						remove(dump.c_str());
					}
				}

				void make_sp_inp(std::ofstream & ofs)const;
				void make_opt_inp(std::ofstream & ofs)const;


			};
		}
	}
}

#endif
