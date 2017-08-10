#include "energy_int_chemshell.h"

template<typename T, typename U>
auto zip(T && a, U && b) {
	std::vector<decltype(std::make_pair(*a.begin(), *b.begin()))> ret_vec;

	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(ret_vec),
		[](auto && a_i, auto && b_i) {
		return std::make_pair(a_i, b_i);
	});
	return ret_vec;
}

void energy::interfaces::chemshell::sysCallInterface::initialize_before_first_use()const {
	//Fails if coordinates are not known until this point
	if (Config::get().energy.chemshell.extra_pdb == "") {
		create_pdb();
	}
	//How to deal with coordinates...
	else {
		std::stringstream ss;
		ss << "cp " << Config::get().energy.chemshell.extra_pdb << " " << tmp_file_name << ".pdb";

		auto ret = scon::system_call(ss.str());

		if (ret) {
			throw std::runtime_error("Failed to copy the given PDB file!");
		}
	}
	//TODO: Hopefully needs to be called only once so try! <- Antechamber fails the second time it is called ...
	call_tleap();
}

void energy::interfaces::chemshell::sysCallInterface::create_pdb() const {
	
	write_xyz(tmp_file_name + ".xyz");

	std::stringstream ss;
	ss << Config::get().energy.chemshell.babel_path << " -ixyz " << tmp_file_name << ".xyz -opdb " << tmp_file_name << ".pdb";

	auto ret = scon::system_call(ss.str());

	if (ret) {
		throw std::runtime_error("Failed to call babel!");
	}

}

void energy::interfaces::chemshell::sysCallInterface::write_xyz(std::string const & o_file) const {
	std::ofstream xyz_file(o_file);
	xyz_file << coords->size() << "\n\n";
	xyz_file << coords::output::formats::xyz(*coords);
	xyz_file.close();
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

	auto ret = scon::system_call(ss.str());

	if (ret) {
		throw std::runtime_error("Failed to call antechamber!");
	}

	// To empty ss
	std::stringstream().swap(ss);

	ss << "parmchk -i " << o_file << ".mol2 -f mol2 -o " << o_file << ".frcmod";

	ret = scon::system_call(ss.str());

	if (ret) {
		throw std::runtime_error("Failed to call parmchk!");
	}

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

void energy::interfaces::chemshell::sysCallInterface::make_sp_inp(std::ofstream & ofs) const {

	auto const & mxlist = Config::get().energy.chemshell.mxlist;
	auto const & cutoff = Config::get().energy.chemshell.cutoff;
	auto const & embedding_sheme = Config::get().energy.chemshell.scheme;
	auto const & qm_ham = Config::get().energy.chemshell.qm_ham;
	auto const & qm_theory = Config::get().energy.chemshell.qm_theory;
	auto const & qm_region = Config::get().energy.chemshell.qm_atoms;

	ofs << "eandg coords = ${dir}/${sys_name_id}.c \\\n"
		"    theory=hybrid : [ list \\\n";
	if (embedding_sheme != "") {
		ofs << "        coupling= " << embedding_sheme << " \\\n";
	}
	if (qm_theory != "") {
		ofs << "        qm_theory= " << qm_theory << " : [ list \\\n";
		if (qm_ham != "") {
			ofs << "            hamiltonian= " << qm_ham << " \\\n";
		}
		//This can't be right ask Daniel!
		//		"            basis=$qm_basis\\\n"
		ofs << "            ] \\\n";
	}
	if (qm_region != "") {
		ofs << "    qm_region= { " << qm_region << " } \\\n";
	}
	ofs << "    debug=no \\\n"
		"    mm_theory= dl_poly : [ list \\\n"
		"        list_option=none \\\n"
		"        conn= ${sys_name_id}.c \\\n"
		"        mm_defs=$amber_prmtop \\\n"
		"        exact_srf=yes \\\n";
	if (mxlist != "") {
		ofs << "        mxlist=" << mxlist << " \\\n";
	}
	if (cutoff != "") {
		ofs << "        cutoff=" << cutoff << " \\\n";
	}
	ofs << "        scale14 = {1.2 2.0}\\\n"
		"        amber_prmtop_file=$amber_prmtop ] ] \n"
		"    energy=" << tmp_file_name << ".energy\\\n"
		"    gradient=" << tmp_file_name << ".gradient\n"
		"\n\n\n";
//		"\n"
//		"close $control_input_settings\n";
}

void energy::interfaces::chemshell::sysCallInterface::make_opt_inp(std::ofstream & ofs) const {

	auto const & maxcycle = Config::get().energy.chemshell.maxcycle;
	auto const & maxcyc = Config::get().energy.chemshell.maxcyc;
	auto const & tolerance = Config::get().energy.chemshell.tolerance;
	auto const & mxlist = Config::get().energy.chemshell.mxlist;
	auto const & cutoff = Config::get().energy.chemshell.cutoff;
	auto const & qm_basis = Config::get().energy.chemshell.qm_basis;
	auto const & embedding_sheme = Config::get().energy.chemshell.scheme;
	auto const & qm_ham = Config::get().energy.chemshell.qm_ham;
	auto const & qm_ch = Config::get().energy.chemshell.qm_charge;
	auto const & qm_theory = Config::get().energy.chemshell.qm_theory;
	auto const & qm_region = Config::get().energy.chemshell.qm_atoms;

	std::string active_atoms = find_active_atoms();


	ofs << "dl-find coords = ${dir}/${sys_name_id}.c \\\n"
		"    coordinates=hdlc \\\n"
		"    result=" << tmp_file_name << ".c \\\n";
	if (maxcycle != "") {
		ofs << "    maxcycle=" << maxcycle << " \\\n";
	}
	if (tolerance != "") {
		ofs << "    tolerance=" << tolerance << " \\\n";
	}
	ofs << "    active_atoms= { " << active_atoms << "} \\\n"
		"    residues= $residues \\\n"
		"    theory=hybrid : [ list \\\n";
	if (embedding_sheme != "") {
		ofs << "        coupling=" << embedding_sheme << " \\\n";
	}
	if (qm_theory != "") {
		ofs << "        qm_theory= " << qm_theory << " : [ list \\\n";
		if (qm_ham != "") {
			ofs << "            hamiltonian = " << qm_ham << " \\\n";
		}
		if (qm_basis != "") {
			ofs << "            basis= " << qm_basis << " \\\n";
		}
		if (maxcyc != "") {
			ofs << "            maxcyc= " << maxcyc << " \\\n";
		}
		if (Config::get().energy.chemshell.dispersion) {
			ofs << "            dispersion_correction= " << qm_ham << " \\\n";
		}
		if (qm_ch != "") {
			ofs << "            charge= " << qm_ch << " ] \\\n";
		}
	}
	if (qm_region != "") {
		ofs << "    qm_region = { " << qm_region << " } \\\n";
	}
	ofs << "    debug=no \\\n"
		"    mm_theory= dl_poly : [ list \\\n"
		"        list_option=none \\\n"
		"        conn=" << tmp_file_name << ".c \\\n"
		"        mm_defs=$amber_prmtop \\\n"
		"        exact_srf=yes \\\n";
	if (mxlist != "") {
		ofs << "        mxlist=" << mxlist << " \\\n";
	}
	if (cutoff != "") {
		ofs << "        cutoff=" << cutoff << " \\\n";
	}
	ofs << "        scale14 = {1.2 2.0} \\\n"
		"        amber_prmtop_file=$amber_prmtop ] ] \n"
		"\n"
		"write_xyz file=" << tmp_file_name << ".xyz coords=${sys_name_id}_opt.c\n"
		"read_pdb  file=" << tmp_file_name << ".pdb  coords=dummy.coords\n"
		"write_pdb file=" << tmp_file_name << ".pdb coords=${sys_name_id}_opt.c\n\n\n";
/*		"read_pdb  file=${ sys_name_id }.pdb  coords=dummy.coords\n"
		"write_pdb file=${ sys_name_id }_opt.pdb coords=${ sys_name_id }_opt.c\n"
		"write_xyz file=${ sys_name_id }_qm_region_opt.xyz coords=hybrid.${ qm_theory }.coords\n"
		"delete_object hybrid.${ qm_theory }.coords\n"
		"catch {file delete dummy.coords}\n"
		"\n"
		"close $control_input_settings\n";*/
}

void energy::interfaces::chemshell::sysCallInterface::write_chemshell_coords()const {

	auto o_file = tmp_file_name + ".chm";
	write_xyz(tmp_file_name + ".xyz");

	std::ofstream chemshell_file_to_prepare_coords(o_file);

	chemshell_file_to_prepare_coords <<
		"set sys_name_id " << tmp_file_name << "\n"
		"read_xyz file=./${sys_name_id}.xyz coords=./${sys_name_id}.c";

	chemshell_file_to_prepare_coords.close();

	actual_call();

}


/*
void energy::interfaces::chemshell::sysCallInterface::write_chemshell_coords()const{
	auto o_file = tmp_file_name + ".chm";

	std::ofstream chemshell_file_to_prepare_coords(o_file);

	chemshell_file_to_prepare_coords << 
		"set sys_name_id " << tmp_file_name << "\n"
		"read_pdb file=./${sys_name_id}.pdb coords=./${sys_name_id}.c";

	chemshell_file_to_prepare_coords.close();

	actual_call();

}
*/
void energy::interfaces::chemshell::sysCallInterface::write_chemshell_file(bool const & sp) const {
	
	auto o_file = tmp_file_name + ".chm";

	std::ofstream chem_shell_input_stream(o_file);

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
		"load_amber_coords inpcrd=$amber_inpcrd prmtop=$amber_prmtop coords=${sys_name_id}.c\n"
		"\n"
		"set residues [pdb_to_res \"${sys_name_id}.pdb\"]\n";
	//Refactoring NEEDED!!!!!!
	auto const & cov_res = Config::get().energy.chemshell.cov_residues;
	
	if (cov_res != "") {
		chem_shell_input_stream << "set residues [ inlist function= combine residues= $residues sets= {" << cov_res << "} target= MOX ]\n\n";
	}
	else {
		chem_shell_input_stream << "\n";
	}
	//	"flush $control_input_settings\n"
	//	"\n"
		;
	if (sp) {
		make_sp_inp(chem_shell_input_stream);
	}
	else {
		make_opt_inp(chem_shell_input_stream);
	}
		
	chem_shell_input_stream.close();

}

void energy::interfaces::chemshell::sysCallInterface::make_opti() const {
	call_chemshell(false);
}

void energy::interfaces::chemshell::sysCallInterface::make_sp()const {
	call_chemshell();
}

std::string energy::interfaces::chemshell::sysCallInterface::find_active_atoms() const {
	
	/*std::vector<int> indices(coords->size());
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

	
	for (auto const & i : final_vec) {
		if (i != 0) {
			final_atoms += std::to_string(i) + " ";
		}
	}*/
	std::string final_atoms = "";
	auto const & atoms = coords->atoms();
	for (auto i = 0; i < coords->size(); ++i) {
		auto const & atom = coords->atoms().atom(i);
		if (!atom.fixed()) {
			final_atoms += std::to_string(i + 1) + " ";
		}
	}
	std::cout << final_atoms << std::endl;
	return final_atoms;
}


void energy::interfaces::chemshell::sysCallInterface::call_chemshell(bool singlepoint) const {
	
	//create_pdb();
	write_chemshell_file(singlepoint);
	actual_call();

}

void energy::interfaces::chemshell::sysCallInterface::actual_call()const {

	std::stringstream chemshell_stream;
	chemshell_stream << Config::get().energy.chemshell.path << " " << tmp_file_name << ".chm";

	auto failcount = 0;

	for (; failcount <= 10; ++failcount) {
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

bool energy::interfaces::chemshell::sysCallInterface::check_if_number(std::string const & number) const {

	return !number.empty() && std::find_if(number.cbegin(), number.cend(), [](char n) {
		return n != 'E' && n != 'e' && n != '-' && n != '+' && n != '.' && !std::isdigit(n); //check if the line contains digits, a minus or a dot to determine if its a floating point number
	}) == number.end();

}

coords::float_type energy::interfaces::chemshell::sysCallInterface::read_energy()const {
	std::ifstream ifile(tmp_file_name + ".energy");

	std::string line;
	while (getline(ifile, line)) {
		std::istringstream iss(line);
		std::vector<std::string> words{
			std::istream_iterator<std::string>{iss},
			std::istream_iterator<std::string>{}
		};
		if (words.size()==0) {
			continue;
		}
		if(check_if_number(words.at(0))){
			return std::stod(words.at(0));
		}
	}

	return 0.0;

}

coords::Representation_3D energy::interfaces::chemshell::sysCallInterface::extract_gradients(std::vector<coords::float_type> const & grads) const {
	coords::Representation_3D new_grads;
	for (auto b = grads.cbegin(); b < grads.cend(); b += 3) {
		new_grads.emplace_back(coords::Cartesian_Point(
			*b,
			*(b+1),
			*(b+2)
		));
	}
	return new_grads;
}

void energy::interfaces::chemshell::sysCallInterface::read_gradients() {
	std::ifstream ifile(tmp_file_name + ".gradient");

	std::string line;

	std::vector<coords::float_type> gradients;

	while (getline(ifile, line)) {
		std::istringstream iss(line);
		std::vector<std::string> words{
			std::istream_iterator<std::string>{iss},
			std::istream_iterator<std::string>{}
		};
		if (words.size() == 0) {
			continue;
		}
		if (check_if_number(words.at(0))) {
			gradients.emplace_back(std::stod(words.at(0)));
		}
	}

	auto new_gradients = extract_gradients(gradients);

	coords->swap_g_xyz(new_gradients);

}

bool energy::interfaces::chemshell::sysCallInterface::check_if_line_is_coord(std::vector<std::string> const & coords)const {
	return 
		coords.size() == 4 &&
		check_if_number(coords.at(1)) && 
		check_if_number(coords.at(2)) && 
		check_if_number(coords.at(3));
}

coords::Cartesian_Point energy::interfaces::chemshell::sysCallInterface::make_coords(std::vector<std::string> const & line) const {
	std::vector<std::string> coord_words(line.cbegin() + 1, line.cend());
	coords::Cartesian_Point cp(
		std::stod(coord_words.at(0)),
		std::stod(coord_words.at(1)),
		std::stod(coord_words.at(2))
	);
	return cp;
}

void energy::interfaces::chemshell::sysCallInterface::change_name_of_energy_and_grad()const {
	std::stringstream ss;

	ss << "mv dl-find.energy " << tmp_file_name << ".energy";

	scon::system_call(ss.str());

	std::stringstream().swap(ss);

	ss << "mv dl-find.gradient " << tmp_file_name << ".gradient";

	scon::system_call(ss.str());

}

void energy::interfaces::chemshell::sysCallInterface::make_optimized_coords_to_actual_coords(coords::Representation_3D const & xyz) {
	coords->set_xyz(xyz);

	/*std::stringstream ss;
	ss << "mv " << tmp_file_name << "_opt.c " << tmp_file_name << ".c";

	auto ret = scon::system_call(ss.str());
	/*
	if (ret) {
		throw std::runtime_error("Failed to replace unoptimzed .c file with the optimized one.");
	}
	*/
}

void energy::interfaces::chemshell::sysCallInterface::change_input_file_names(std::string const & filename, std::string const & copy_or_move) const {

	std::string what_happens = copy_or_move == "cp" ? "copy" : "move";

	std::stringstream ss;
	ss << copy_or_move << " " << filename << ".c " << tmp_file_name << ".c";

	auto ret = scon::system_call(ss.str());

	if (ret) {
		throw std::runtime_error("Failed to " + what_happens + " the old .c file.");
	}

	std::stringstream().swap(ss);

	ss << copy_or_move << " " << filename << ".prmtop " << tmp_file_name << ".prmtop";

	ret = scon::system_call(ss.str());

	if (ret) {
		throw std::runtime_error("Failed to " + what_happens + " the old .prmtop file.");
	}

	std::stringstream().swap(ss);

	ss << copy_or_move << " " << filename << ".inpcrd " << tmp_file_name << ".inpcrd";

	ret = scon::system_call(ss.str());

	if (ret) {
		throw std::runtime_error("Failed to " + what_happens + " the old .inpcrd file.");
	}

	std::stringstream().swap(ss);

	ss << copy_or_move << " " << filename << ".pdb " << tmp_file_name << ".pdb";

	ret = scon::system_call(ss.str());

	if (ret) {
		throw std::runtime_error("Failed to " + what_happens + " the old .pdb file.");
	}
}

void energy::interfaces::chemshell::sysCallInterface::read_coords() {
	std::ifstream ifile(tmp_file_name+".xyz");

	std::string line;
	coords::Representation_3D xyz;

	while (getline(ifile, line)) {
		std::istringstream iss(line);
		std::vector<std::string> words{
			std::istream_iterator<std::string>{iss},
			std::istream_iterator<std::string>{}
		};
		if(words.size()==0){
			continue;
		}
		
		if (check_if_line_is_coord(words)) {
			xyz.emplace_back(make_coords(words));
		}
	}

	make_optimized_coords_to_actual_coords(xyz);

	ifile.close();
}

void energy::interfaces::chemshell::sysCallInterface::swap(interface_base & other){}
energy::interface_base * energy::interfaces::chemshell::sysCallInterface::clone(coords::Coordinates * coord_object) const { return new sysCallInterface(*this, coord_object); }
energy::interface_base * energy::interfaces::chemshell::sysCallInterface::move(coords::Coordinates * coord_object) { return new sysCallInterface(*this, coord_object); }

void energy::interfaces::chemshell::sysCallInterface::update(bool const skip_topology){}

coords::float_type energy::interfaces::chemshell::sysCallInterface::e(void) { 
	check_for_first_call();
	write_chemshell_coords();
	make_sp();
	return read_energy();
}
coords::float_type energy::interfaces::chemshell::sysCallInterface::g(void) {
	check_for_first_call();
	write_chemshell_coords();
	make_sp();
	read_gradients();
	return read_energy();
}
coords::float_type energy::interfaces::chemshell::sysCallInterface::h(void) {
	check_for_first_call();
	write_chemshell_coords();
	return 0.0; 
}
coords::float_type energy::interfaces::chemshell::sysCallInterface::o(void) {
	check_for_first_call();
	write_chemshell_coords();
	make_opti();
	change_name_of_energy_and_grad();
	read_gradients();
	read_coords();
	return read_energy(); 
}

void energy::interfaces::chemshell::sysCallInterface::print_E(std::ostream&) const{}

void energy::interfaces::chemshell::sysCallInterface::print_E_head(std::ostream&, bool const endline) const {}

void energy::interfaces::chemshell::sysCallInterface::print_E_short(std::ostream&, bool const endline) const {}

void energy::interfaces::chemshell::sysCallInterface::print_G_tinkerlike(std::ostream & S, bool const aggregate) const {

	coords::Representation_3D gradients;

	coords->get_g_xyz(gradients);

	for (auto const & grad : gradients) {

		S << std::right << std::setw(16) << std::scientific << std::setprecision(5) << grad << "\n";

	}

}

void energy::interfaces::chemshell::sysCallInterface::to_stream(std::ostream&) const {}
