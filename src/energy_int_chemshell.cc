#include "energy_int_chemshell.h"
#include "helperfunctions.h"

template<typename T, typename U>
auto zip(T&& a, U&& b) {
  std::vector<decltype(std::make_pair(*a.begin(), *b.begin()))> ret_vec;

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(ret_vec),
    [](auto&& a_i, auto&& b_i) {
    return std::make_pair(a_i, b_i);
  });
  return ret_vec;
}

void energy::interfaces::chemshell::sysCallInterface::initialize_before_first_use()const {

  if (Config::get().energy.chemshell.extra_pdb == "") {
    create_pdb();
  }
  else {
    std::stringstream ss;
    ss << "cp " << Config::get().energy.chemshell.extra_pdb << " " << tmp_file_name << ".pdb";

    auto ret = scon::system_call(ss.str());

    if (ret) {
      throw std::runtime_error("Failed to copy the given PDB file!");
    }
  }
  if (Config::get().energy.chemshell.optional_inpcrd == "" || Config::get().energy.chemshell.optional_prmtop == "") {
    call_tleap();
  }
  else {
    std::stringstream ss;
    ss << "cp " << Config::get().energy.chemshell.optional_inpcrd << " " << tmp_file_name << ".inpcrd";

    auto ret = scon::system_call(ss.str());

    if (ret) {
      throw std::runtime_error("Failed to copy the given frcmod file!");
    }

    std::stringstream().swap(ss);

    ss << "cp " << Config::get().energy.chemshell.optional_prmtop << " " << tmp_file_name << ".prmtop";

    ret = scon::system_call(ss.str());

    if (ret) {
      throw std::runtime_error("Failed to copy the given prmtop file!");
    }
  }
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

void energy::interfaces::chemshell::sysCallInterface::write_xyz(std::string const& o_file) const {
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

void energy::interfaces::chemshell::sysCallInterface::make_tleap_input(std::string const& o_file)const {

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

}

void energy::interfaces::chemshell::sysCallInterface::make_sp_inp(std::ofstream& ofs) const {

  auto const& mxlist = Config::get().energy.chemshell.mxlist;
  auto const& cutoff = Config::get().energy.chemshell.cutoff;
  auto const& embedding_sheme = Config::get().energy.chemshell.scheme;
  auto const& qm_ham = Config::get().energy.chemshell.qm_ham;
  auto const& qm_theory = Config::get().energy.chemshell.qm_theory;
  auto const& qm_region = Config::get().energy.chemshell.qm_atoms;
  auto const& qm_basis = Config::get().energy.chemshell.qm_basis;

  ofs << "eandg coords = ./" << tmp_file_name << ".c \\\n"
    "    theory=hybrid : [ list \\\n";
  if (embedding_sheme != "") {
    ofs << "        coupling=" << embedding_sheme << " \\\n";
  }
  if (qm_theory != "") {
    ofs << "        qm_theory=" << qm_theory << " : [ list \\\n";
    if (qm_ham != "") {
      ofs << "            hamiltonian=" << qm_ham << " \\\n";
    }
    if (qm_basis != "") {
      ofs << "            basis=" << qm_basis << " ] \\\n";
    }
  }
  if (qm_region != "") {
    ofs << "    qm_region= { " << qm_region << " } \\\n";
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
  ofs << "        scale14 = {1.2 2.0}\\\n"
    "        amber_prmtop_file=$amber_prmtop ] ]\\\n"
    "    energy=" << tmp_file_name << ".energy\\\n"
    "    gradient=" << tmp_file_name << ".gradient\n"
    "\n\n\n";
  //		"\n"
  //		"close $control_input_settings\n";
}

void energy::interfaces::chemshell::sysCallInterface::make_opt_inp(std::ofstream& ofs) const {

  auto const& coordinates = Config::get().energy.chemshell.coords;
  auto const& maxcycle = Config::get().energy.chemshell.maxcycle;
  auto const& maxcyc = Config::get().energy.chemshell.maxcyc;
  auto const& tolerance = Config::get().energy.chemshell.tolerance;
  auto const& mxlist = Config::get().energy.chemshell.mxlist;
  auto const& cutoff = Config::get().energy.chemshell.cutoff;
  auto const& qm_basis = Config::get().energy.chemshell.qm_basis;
  auto const& embedding_sheme = Config::get().energy.chemshell.scheme;
  auto const& qm_ham = Config::get().energy.chemshell.qm_ham;
  auto const& qm_ch = Config::get().energy.chemshell.qm_charge;
  auto const& qm_theory = Config::get().energy.chemshell.qm_theory;
  auto const& qm_region = Config::get().energy.chemshell.qm_atoms;
  auto const& scale14 = Config::get().energy.chemshell.scale14;

  std::string active_atoms, inactive_atoms;
  std::tie(active_atoms, inactive_atoms) = find_active_and_inactive_atoms(qm_region);

  auto make_constraints = [&](constraints const& c) {
    std::stringstream ret;
    ret << "{ " << c.kind;
    for (auto const& a : c.atoms) {
      ret << " " << std::to_string(a);
    }
    ret << " } ";
    return ret.str();
  };

  ofs << "dl-find coords = ./" << tmp_file_name << ".c \\\n"
    "    coordinates=";
  if (coordinates != "") {
    ofs << coordinates;
  }
  else {
    ofs << "hdlc";
  }
  ofs << " \\\n";
  if (!cons.empty()) {
    std::string constraint_str = "";
    for (auto const& c : cons) {
      constraint_str += make_constraints(c);
    }
    ofs << "    constraints= { " << constraint_str << "} \\\n";
  }
  ofs << "    result=" << tmp_file_name << "_opt.c \\\n";
  if (maxcycle != "") {
    ofs << "    maxcycle=" << maxcycle << " \\\n";
  }
  if (tolerance != "") {
    ofs << "    tolerance=" << tolerance << " \\\n";
  }
  ofs << "    active_atoms= {" << active_atoms << "} \\\n";
  /*if (inactive_atoms != "") {
    ofs << "    frozen= {" << inactive_atoms << "} \\\n";
  }*/
  ofs << "    residues= $residues \\\n"
    "    theory=hybrid : [ list \\\n";
  if (embedding_sheme != "") {
    ofs << "        coupling=" << embedding_sheme << " \\\n";
  }
  if (qm_theory != "") {
    ofs << "        qm_theory= " << qm_theory << " : [ list \\\n";
    if (qm_ham != "") {
      ofs << "            hamiltonian= " << qm_ham << " \\\n";
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
  if (scale14 != "") {
    ofs << "        scale14 = {" << scale14 << "} \\\n";
  }
  else {
    ofs << "        scale14 = {1.2 2.0} \\\n";
  }
  ofs << "        amber_prmtop_file=$amber_prmtop ] ] \n"
    "\n"
    /*Consider qm_theory.energy and qm_theory.gradient as a way to produce output with the right name in the first place*/
    "write_xyz file=" << tmp_file_name << ".xyz coords=" << tmp_file_name << "_opt.c\n"
    "read_pdb file=" << tmp_file_name << ".pdb coords=dummy.coords\n"
    "write_pdb file=" << tmp_file_name << ".pdb coords=" << tmp_file_name << "_opt.c\n\n\n";
  /*		"read_pdb  file=${ sys_name_id }.pdb  coords=dummy.coords\n"
      "write_pdb file=${ sys_name_id }_opt.pdb coords=${ sys_name_id }_opt.c\n"
      "write_xyz file=${ sys_name_id }_qm_region_opt.xyz coords=hybrid.${ qm_theory }.coords\n"
      "delete_object hybrid.${ qm_theory }.coords\n"
      "catch {file delete dummy.coords}\n"
      "\n"
      "close $control_input_settings\n";*/
}

void energy::interfaces::chemshell::sysCallInterface::eval_constraints()
{
  auto which_bond = [&](std::vector<std::string> const& strs) {
    constraints c;
    c.kind = "bond";
    c.atoms.emplace_back(std::stoi(strs[1]));
    c.atoms.emplace_back(std::stoi(strs[2]));
    this->cons.emplace_back(c);
  };
  auto which_angle = [&](std::vector<std::string> const& strs) {
    constraints c;
    c.kind = "angle";
    c.atoms.emplace_back(std::stoi(strs[1]));
    c.atoms.emplace_back(std::stoi(strs[2]));
    c.atoms.emplace_back(std::stoi(strs[3]));
    this->cons.emplace_back(c);
  };
  auto which_dihedral = [&](std::vector<std::string> const& strs) {
    constraints c;
    c.kind = "torsion";
    c.atoms.emplace_back(std::stoi(strs[1]));
    c.atoms.emplace_back(std::stoi(strs[2]));
    c.atoms.emplace_back(std::stoi(strs[3]));
    c.atoms.emplace_back(std::stoi(strs[4]));
    this->cons.emplace_back(c);
  };

  auto which_const = [&](std::vector<std::string> const& strs) {
    if (strs[0] == "bond") {
      which_bond(strs);
    }
    else if (strs[0] == "angle") {
      which_angle(strs);
    }
    else if (strs[0] == "dihedral") {
      which_dihedral(strs);
    }
  };

  if (!Config::get().scan2d.constraints) return;
  if (cons.empty()) {
    for (auto const& el : Config::get().scan2d.AXES) {
      std::istringstream iss(el);
      std::vector<std::string> split{ std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{} };
      which_const(split);
    }
  }
}

void energy::interfaces::chemshell::sysCallInterface::write_chemshell_coords()const {

  auto o_file = tmp_file_name + ".chm";
  write_xyz(tmp_file_name + ".xyz");
  write_xyz(tmp_file_name + std::to_string(x) + ".xyz");

  std::ofstream chemshell_file_to_prepare_coords(o_file);

  chemshell_file_to_prepare_coords <<
    "read_xyz file=./" << tmp_file_name << ".xyz coords=./" << tmp_file_name << ".c";

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
void energy::interfaces::chemshell::sysCallInterface::write_chemshell_file(bool const& sp) const {

  auto o_file = tmp_file_name + ".chm";

  std::ofstream chem_shell_input_stream(o_file);




  write_xyz(tmp_file_name + ".xyz");
  write_xyz(tmp_file_name + std::to_string(x) + ".xyz");

  chem_shell_input_stream <<
    "global qm_theory\n"
    //		"global ftupd\n"
    "\n"
    "set amber_prmtop " << tmp_file_name << ".prmtop\n"
    "set amber_inpcrd " << tmp_file_name << ".inpcrd\n"
    "\n"
    //		"set control_input_settings [ open control_input.${sys_name_id}  a ]\n"
    //		"\n"
    "load_amber_coords inpcrd=$amber_inpcrd prmtop=$amber_prmtop coords=" << tmp_file_name << ".c\n"
    "\n"
    "set residues [pdb_to_res \"" << tmp_file_name << ".pdb\"]\n"
    "\n"
    "read_xyz file=" << tmp_file_name << ".xyz coords=" << tmp_file_name << ".c\n";
  //Refactoring NEEDED!!!!!!
  auto const& com_res = Config::get().energy.chemshell.com_residues;

  if (com_res != "") {
    chem_shell_input_stream << "set residues [ inlist function= combine residues= $residues sets= {" << com_res << "} target= MOX ]\n\n";
  }
  else {
    chem_shell_input_stream << "\n";
  }
  //	"flush $control_input_settings\n"
  //	"\n"
  //	;
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

std::pair<std::string, std::string> energy::interfaces::chemshell::sysCallInterface::find_active_and_inactive_atoms(std::string const& qm_atoms) const {

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

  auto make_string_to_numbers = [&](std::string const& str) {

    std::istringstream iss(str);
    std::vector<std::string> str_num{ std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{} };

    std::vector<std::size_t> ret;

    for (auto i = 0u; i < str_num.size(); ++i) {
      auto const& str_i = str_num[i];
      if (i == str_num.size() - 1) {
        ret.emplace_back(std::stoi(str_i));
        break;
      }

      if (str_num[i + 1u] != "-" && check_if_number(str_i)) {
        ret.emplace_back(std::stoi(str_i));
      }
      else if (str_num[i + 1u] == "-") {
        auto start = std::stoi(str_i);
        auto end = std::stoi(str_num[i + 2u]);
        std::vector<std::size_t> append_vec(end - start + 1u);
        std::iota(append_vec.begin(), append_vec.end(), start);
        ret.insert(ret.end(), append_vec.begin(), append_vec.end());
        i += 2u;
      }

    }

    return ret;
  };

  auto const border = std::stod(Config::get().energy.chemshell.active_radius);

  std::string active_atoms = "";

  auto qm_list = make_string_to_numbers(qm_atoms);

  std::set<std::size_t> active_atoms_set;

  for (auto const& qm_atom : qm_list) {
    active_atoms_set.insert(qm_atom);
    auto const& xyz = coords->xyz();
    auto const& center = xyz[qm_atom - 1u];
    for (auto i = 0u; i < xyz.size(); ++i) {
      if (scon::geometric_length(center - xyz[i]) < border) {
        active_atoms_set.insert(i + 1u);
      }
    }
  }
  /*for(auto const & qm_atom : qm_list){
      active_atoms_set.erase(qm_atom);
  }*/
  for (auto const& c : cons) {
    for (auto const& a : c.atoms) {
      active_atoms_set.insert(a);
    }
  }

  for (auto const& aa : active_atoms_set) {
    active_atoms += std::to_string(aa) + " ";
  }

  std::string fixed_atoms = "";
  for (auto i = 0u; i < coords->size(); ++i) {
    auto const& atom = coords->atoms().atom(i);
    if (atom.fixed()) {
      fixed_atoms += std::to_string(i + 1) + " ";
    }
  }
  return std::make_pair(trim_space_and_tabs(active_atoms), trim_space_and_tabs(fixed_atoms));
}

std::string energy::interfaces::chemshell::sysCallInterface::trim_space_and_tabs(std::string const& str) {

  if (str == "") {
    return str;
  }
  auto const first = str.find_first_not_of(" \t");
  auto const last = str.find_last_not_of(" \t");

  auto const length = last - first;

  return str.substr(first, length + 1);
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

coords::float_type energy::interfaces::chemshell::sysCallInterface::read_energy()const {
  std::ifstream ifile(tmp_file_name + ".energy");

  std::string line;
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
      return std::stod(words.at(0));
    }
  }

  return 0.0;

}

coords::Representation_3D energy::interfaces::chemshell::sysCallInterface::extract_gradients(std::vector<coords::float_type> const& grads) const {
  coords::Representation_3D new_grads;
  for (auto b = grads.cbegin(); b < grads.cend(); b += 3) {
    new_grads.emplace_back(coords::Cartesian_Point(
      *b,
      *(b + 1),
      *(b + 2)
    ));
  }
  return new_grads;
}

coords::Representation_3D energy::interfaces::chemshell::sysCallInterface::read_gradients() const {
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

  return extract_gradients(gradients);
}

bool energy::interfaces::chemshell::sysCallInterface::check_if_line_is_coord(std::vector<std::string> const& coordobj)const {
  return
    coordobj.size() == 4 &&
    check_if_number(coordobj.at(1)) &&
    check_if_number(coordobj.at(2)) &&
    check_if_number(coordobj.at(3));
}

coords::Cartesian_Point energy::interfaces::chemshell::sysCallInterface::make_coords(std::vector<std::string> const& line) const {
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

  std::stringstream().swap(ss);

  ss << "rm " << tmp_file_name << "_opt.c " << tmp_file_name << ".c";

  scon::system_call(ss.str());


}

void energy::interfaces::chemshell::sysCallInterface::change_input_file_names(std::string const& filename, std::string const& copy_or_move) const {

  std::string what_happens = copy_or_move == "cp" ? "copy" : "move";

  std::stringstream ss;
  /*ss << copy_or_move << " " << filename << ".c " << tmp_file_name << ".c";

  auto ret = scon::system_call(ss.str());

  if (ret) {
    throw std::runtime_error("Failed to " + what_happens + " the old .c file.");
  }

  std::stringstream().swap(ss);*/

  ss << copy_or_move << " " << filename << ".prmtop " << tmp_file_name << ".prmtop";

  auto ret = scon::system_call(ss.str());

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
  std::ifstream ifile(tmp_file_name + ".xyz");

  std::string line;
  coords::Representation_3D xyz;

  while (getline(ifile, line)) {
    std::istringstream iss(line);
    std::vector<std::string> words{
      std::istream_iterator<std::string>{iss},
      std::istream_iterator<std::string>{}
    };
    if (words.size() == 0) {
      continue;
    }

    if (check_if_line_is_coord(words)) {
      xyz.emplace_back(make_coords(words));
    }
  }

  coords->set_xyz(xyz, true);

  ifile.close();

  std::stringstream ss;

  ss << "rm " << tmp_file_name << ".xyz";

  auto ret = scon::system_call(ss.str());

  if (ret) {
    throw std::runtime_error("Failed to remove the red out .xyz file.");
  }
}

void energy::interfaces::chemshell::sysCallInterface::swap(interface_base& other)
{
  swap(dynamic_cast<sysCallInterface&>(other));
}
energy::interface_base* energy::interfaces::chemshell::sysCallInterface::clone(coords::Coordinates* coord_object) const { return new sysCallInterface(*this, coord_object); }
energy::interface_base* energy::interfaces::chemshell::sysCallInterface::move(coords::Coordinates* coord_object) { return new sysCallInterface(*this, coord_object); }

coords::float_type energy::interfaces::chemshell::sysCallInterface::e(void) {
  check_for_first_call();
  //write_chemshell_coords();
  make_sp();
  return read_energy() * au2kcal_mol;
}
coords::float_type energy::interfaces::chemshell::sysCallInterface::g(void) {
  check_for_first_call();
  //write_chemshell_coords();
  make_sp();
  set_gradients(read_gradients());
  return read_energy() * au2kcal_mol;
}
coords::float_type energy::interfaces::chemshell::sysCallInterface::h(void) {
  check_for_first_call();
  //write_chemshell_coords();
  return 0.0;
}
coords::float_type energy::interfaces::chemshell::sysCallInterface::o(void) {
  check_for_first_call();
  //write_chemshell_coords();
  make_opti();
  ++x;
  change_name_of_energy_and_grad();
  set_gradients(read_gradients());
  read_coords();
  return read_energy() * au2kcal_mol;
}
