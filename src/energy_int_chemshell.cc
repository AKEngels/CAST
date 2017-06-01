#include "energy_int_chemshell.h"

template<typename T, typename U>
auto zip(T && a, U && b) {
	std::vector<decltype(std::make_pair(*a.begin(), *b.begin()))> ret_vec;

	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(ret_vec),
		[](auto && a_i, auto && b_i) {
		return std::make_pair(a_i, b_i);
	});
	return std::move(ret_vec);
}

std::string energy::interfaces::chemshell::sysCallInterface::create_pdb() {
	auto inp_file_name = tmp_file_name;
	auto xyz_file_name = inp_file_name;
	xyz_file_name.append(".xyz");
	
	write_xyz(xyz_file_name);

	std::stringstream ss;
	ss << "babel -ixyz " << xyz_file_name << " -opdb " << inp_file_name << ".pdb";

	scon::system_call(ss.str());
}

void energy::interfaces::chemshell::sysCallInterface::write_xyz(std::string const & o_file) {
	std::ofstream xyz_file(o_file);
	xyz_file << coords->size() << "\n";
	for (auto && el : zip(coords->atoms() ,coords->xyz())) {
		coords::Atom atom;
		coords::Coordinates coo;
		std::tie(atom, coo) = el;
		xyz_file << atom.symbol() << " " << coo;
	}
}



void energy::interfaces::chemshell::sysCallInterface::swap(interface_base & other){}
energy::interface_base * energy::interfaces::chemshell::sysCallInterface::clone(coords::Coordinates * coord_object) const { return new sysCallInterface(*this, coord_object); }
energy::interface_base * energy::interfaces::chemshell::sysCallInterface::move(coords::Coordinates * coord_object) { return new sysCallInterface(*this, coord_object); }

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