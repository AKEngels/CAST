#include "energy_int_chemshell.h"

std::string energy::interfaces::chemshell::sysCallInterface::create_pdb() {
	auto inp_file_name = tmp_file_name;
	auto xyz_file_name = inp_file_name;
	xyz_file_name.append(".xyz");
	
	std::ofstream xyz_file(xyz_file_name);
	xyz_file << coords->xyz();

	std::stringstream ss;
	ss << "babel -ixyz " << xyz_file_name << " -opdb " << inp_file_name << ".pdb";

	scon::system_call(ss.str());
}





void energy::interfaces::chemshell::sysCallInterface::swap(interface_base & other){}
energy::interface_base * energy::interfaces::chemshell::sysCallInterface::clone(coords::Coordinates * coord_object) const { return new sysCallInterface(*this); }
energy::interface_base * energy::interfaces::chemshell::sysCallInterface::move(coords::Coordinates * coord_object) { return new sysCallInterface(coord_object); }

void energy::interfaces::chemshell::sysCallInterface::update(bool const skip_topology){}

coords::float_type energy::interfaces::chemshell::sysCallInterface::e(void) { return 0.0; }
coords::float_type energy::interfaces::chemshell::sysCallInterface::g(void) { return 0.0; }
coords::float_type energy::interfaces::chemshell::sysCallInterface::h(void) { return 0.0; }
coords::float_type energy::interfaces::chemshell::sysCallInterface::o(void) { return 0.0; }

void energy::interfaces::chemshell::sysCallInterface::print_E(std::ostream&) const{}
void energy::interfaces::chemshell::sysCallInterface::print_E_head(std::ostream&, bool const endline) const {}
void energy::interfaces::chemshell::sysCallInterface::print_E_short(std::ostream&, bool const endline) const {}
void energy::interfaces::chemshell::sysCallInterface::print_G_tinkerlike(std::ostream&, bool const aggregate) const {}
void energy::interfaces::chemshell::sysCallInterface::to_stream(std::ostream&) const {}