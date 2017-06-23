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

void energy::interfaces::chemshell::sysCallInterface::create_pdb() const {
	
	write_xyz(tmp_file_name + ".xyz");

	std::stringstream ss;
	ss << "babel -ixyz " << tmp_file_name << ".xyz -opdb " << tmp_file_name << ".pdb";

	scon::system_call(ss.str());

}

void energy::interfaces::chemshell::sysCallInterface::write_xyz(std::string const & o_file) const {
	std::ofstream xyz_file(o_file);
	xyz_file << coords->size() << "\n\n";
	xyz_file << coords::output::formats::xyz(*coords);
	xyz_file.close();
}

void energy::interfaces::chemshell::sysCallInterface::write_input() const {

	call_tleap();
	write_chemshell_file(tmp_file_name + ".chm");
	call_chemshell();

}

void energy::interfaces::chemshell::sysCallInterface::call_tleap()const {
	
	make_tleap_input(tmp_file_name);
	std::stringstream ss;

	ss << "tleap -s -f " << tmp_file_name << ".in > " << tmp_file_name << ".out";

	scon::system_call(ss.str());

}

void energy::interfaces::chemshell::sysCallInterface::make_tleap_input(std::string const & o_file)const {

	std::stringstream ss;

	ss << "antechamber -i " << o_file << ".pdb -fi pdb -o " << o_file << ".mol2 -fo mol2";

	scon::system_call(ss.str());

	// To empty ss
	std::stringstream().swap(ss);

	ss << "parmchk -i " << o_file << ".mol2 -f mol2 -o " << o_file << ".frcmod";

	scon::system_call(ss.str());

	std::stringstream().swap(ss);

	std::ofstream tleap_input(o_file + ".in");

	tleap_input <<
		"source leaprc.gaff\n"
		"LIG = loadmol2 " << o_file << ".mol2\n"
		"check LIG\n"
		"saveof LIG " << o_file << ".lib\n"
		"saveamberparm LIG " << o_file << ".prmtop " << o_file << ".inpcrd\n"
		"savepdb LIG " << o_file << ".pdb\n"
		"quit"
		;

	tleap_input.close();

}

void energy::interfaces::chemshell::sysCallInterface::write_chemshell_file(std::string const & o_file) const {
	
	//auto qm_atoms = parse_qm_atoms();

	constexpr auto maxcycle = 1000;
	constexpr auto maxcyc = 2000;
	constexpr auto tolerance = 0.00045;
	constexpr auto mxlist = 45000;
	constexpr auto cutoff = 1000;

	std::ofstream chem_shell_input_stream(o_file);

	std::string active_atoms = find_active_atoms();

	chem_shell_input_stream <<
		"global sys_name_id\n"
		"global qm_theory\n"
		"global ftupd\n"
		"\n"
		"set dir .\n"
		"set sys_name_id " << tmp_file_name << "\n"
		"\n"
		"set amber_prmtop " << tmp_file_name << ".prmtop\n"
		"set amber_inpcrd " << tmp_file_name << ".inpcrd\n"
		"\n"
		"set control_input_settings [ open control_input.${sys_name_id}  a ]\n"
		"\n"
		"read_pdb file=${dir}/${sys_name_id}.pdb coords=${dir}/${sys_name_id}.c\n"
		"\n"
		"load_amber_coords inpcrd=$amber_inpcrd prmtop=$amber_prmtop coords=${sys_name_id}.c\n"
		"\n"
		"set embedding_scheme " << Config::get().energy.chemshell.scheme << "\n"
		"\n"
		"set qm_theory " << Config::get().energy.chemshell.qm_theory << "\n"
		"puts $control_input_settings \" QM method: $qm_theory \"\n"
		"\n"
		"set qm_ham " << Config::get().energy.chemshell.qm_ham << "\n"
		"\n"
		"set qm_basis " << Config::get().energy.chemshell.qm_basis << "\n"
		"\n"
		"set qm_ch " << Config::get().energy.chemshell.qm_charge << "\n"
		"\n"
		"set qm_atoms  { " << Config::get().energy.chemshell.qm_atoms << " }\n"
		"\n"
		"set residues [pdb_to_res \"${sys_name_id}.pdb\"]\n"
		"\n"
		"flush $control_input_settings\n"
		"\n"
		"dl-find coords = ${dir}/${sys_name_id}.c \\\n"
		"    coordinates=hdlc \\\n"
		"    maxcycle=" << maxcycle << " \\\n"
		"    tolerance=" << tolerance << " \\\n"
		"    active_atoms= { " << active_atoms << "} \\\n"
		"    residues= $residues \\\n"
		"    theory=hybrid : [ list \\\n"
		"        coupling= $embedding_scheme \\\n"
		"        qm_theory= $qm_theory : [ list hamiltonian = $qm_ham \\\n"
		"            basis= $qm_basis \\\n"
		"            maxcyc= " << maxcyc << " \\\n"
		"            dispersion_correction= $qm_ham \\\n"
		"            charge= $qm_ch ] \\\n"
		"    qm_region = $qm_atoms \\\n"
		"    debug=no \\\n"
		"    mm_theory= dl_poly : [ list \\\n"
		"        list_option=none \\\n"
		"        conn= ${sys_name_id}.c \\\n"
		"        mm_defs=$amber_prmtop \\\n"
		"        exact_srf=yes \\\n"
		"        mxlist=" << mxlist << " \\\n"
		"        cutoff=" << cutoff << " \\\n"
		"        scale14 = {1.2 2.0}\\\n"
		"        amber_prmtop_file=$amber_prmtop ] ] \\\n"
		"\n"
		"\n"
		"write_xyz file = ${ sys_name_id }_opt.xyz coords = ${ sys_name_id }_opt.c\n"
		"read_pdb  file = ${ sys_name_id }.pdb  coords = dummy.coords\n"
		"write_pdb file = ${ sys_name_id }_opt.pdb coords = ${ sys_name_id }_opt.c\n"
		"write_xyz file = ${ sys_name_id }_qm_region_opt.xyz coords = hybrid.${ qm_theory }.coords\n"
		"delete_object hybrid.${ qm_theory }.coords\n"
		"catch {file delete dummy.coords}\n"
		"\n"
		"close $control_input_settings\n"
		;

	chem_shell_input_stream.close();


}

std::string energy::interfaces::chemshell::sysCallInterface::find_active_atoms() const {
	
	std::vector<int> indices(coords->size());
	std::iota(indices.begin(), indices.end(), 1);
	std::vector<int> final_vec;

	std::transform(coords->atoms().begin(), coords->atoms().end(), indices.begin(), std::back_inserter(final_vec), 
		[](auto const & a, auto const & b) {
			if (a.fixed()) {
				return 0;
			}
			else {
				return b;
			}
	});

	std::string final_atoms = "";
	for (auto const & i : final_vec) {
		if (i != 0) {
			final_atoms += std::to_string(i) + " ";
		}
	}

	return final_atoms;
}

/*std::vector<std::string> energy::interfaces::chemshell::sysCallInterface::parse_qm_atoms() const {

}*/

void energy::interfaces::chemshell::sysCallInterface::call_chemshell() const {
	
	create_pdb();
	write_input();
	//actual_call();

}

void energy::interfaces::chemshell::sysCallInterface::actual_call()const {

	std::stringstream chemshell_stream;
	chemshell_stream << Config::get().energy.chemshell.path << " < " << tmp_file_name << ".chm > " << tmpfile << ".out";

	auto failcount = 0;

	for (auto failcount = 1; failcount >= 10; ++failcount) {
		auto ret = scon::system_call(chemshell_stream.str());
		if (ret == 0) {
			break;
		}
		else {
			std::cout << "I am failing to call Chemshell! Are you sure you passed the right chemshell path?\n";
			std::cout << "The path passed is: \"" << Config::get().energy.chemshell.path << "\"\n";
		}
	}

	if (failcount == 10) {
		throw std::runtime_error("10 Chemshell calls failed!");
	}
}

void energy::interfaces::chemshell::sysCallInterface::swap(interface_base & other){}
energy::interface_base * energy::interfaces::chemshell::sysCallInterface::clone(coords::Coordinates * coord_object) const { return new sysCallInterface(*this, coord_object); }
energy::interface_base * energy::interfaces::chemshell::sysCallInterface::move(coords::Coordinates * coord_object) { return new sysCallInterface(*this, coord_object); }

void energy::interfaces::chemshell::sysCallInterface::update(bool const skip_topology){}

coords::float_type energy::interfaces::chemshell::sysCallInterface::e(void) { 
	call_chemshell();
	return 123.; 
}
coords::float_type energy::interfaces::chemshell::sysCallInterface::g(void) { return 0.0; }
coords::float_type energy::interfaces::chemshell::sysCallInterface::h(void) { return 0.0; }
coords::float_type energy::interfaces::chemshell::sysCallInterface::o(void) { return 0.0; }

void energy::interfaces::chemshell::sysCallInterface::print_E(std::ostream&) const{}
void energy::interfaces::chemshell::sysCallInterface::print_E_head(std::ostream&, bool const endline) const {}
void energy::interfaces::chemshell::sysCallInterface::print_E_short(std::ostream&, bool const endline) const {}
void energy::interfaces::chemshell::sysCallInterface::print_G_tinkerlike(std::ostream&, bool const aggregate) const {}
void energy::interfaces::chemshell::sysCallInterface::to_stream(std::ostream&) const {}