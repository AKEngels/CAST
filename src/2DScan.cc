#include "2DScan.h"

using length_type = Scan2D::length_type;

Scan2D::Scan2D(coords::Coordinates& coords)
  : _coords(coords), original_verbosity(Config::get().general.verbosity), change_from_atom_to_atom(Config::get().scan2d.change_from_atom_to_atom),
  max_change_rotation(Config::get().scan2d.max_change_to_rotate_whole_molecule)
{
  logfile.open(structures_file);
  energies.open(energie_file);
  before.open("Before_Opti.arc");
  Config::set().scan2d.constraints = true;
  if (Config::get().scan2d.verbose_off) Config::set().general.verbosity = 1;
}

void Scan2D::execute_scan() {
  auto& both_whats = Config::get().scan2d.AXES;

  if (both_whats.size() > 2) {
    throw std::runtime_error("You can't pass more than two axis!");
  }

  auto x_input_parser = parse_input(both_whats.front());
  auto y_input_parser = parse_input(both_whats.back());

  x_input_parser->set_coords(_coords.xyz());
  y_input_parser->set_coords(_coords.xyz());

  auto x_changes = x_input_parser->make_axis();
  auto y_changes = y_input_parser->make_axis();

  parser = std::make_unique<XY_Parser>(std::move(x_input_parser), std::move(y_input_parser));

  axis = std::make_unique<XY_Steps>(x_changes, y_changes);

  make_scan();

  Config::set().general.verbosity = original_verbosity;
}

Scan2D::~Scan2D() {
  //Seems like it is empty ;)
}

void Scan2D::Normal_Input::fill_what(std::vector<std::string>& splitted_vals, coords::Representation_3D const& xyz) {

  what = std::make_unique<Scan2D::what>();

  what->what_kind = splitted_vals.front();

  std::transform(splitted_vals.begin() + 1, splitted_vals.end() - 3, std::back_inserter(what->atoms),
    [](std::string const& atom) { return std::stoi(atom); }
  );

  set_coords(xyz);

  auto const position_of_begin_val = splitted_vals.size() - 3;
  auto const position_of_end_val = splitted_vals.size() - 2;

  auto do_prepare_position = [&]() {
    //      this->what->prepare_position = true;
    this->what->from_position = std::stod(splitted_vals[position_of_begin_val]);
  };

  auto do_not_prepare_position = [&]() {
    //      this->what->prepare_position = false;
    this->what->from_position = say_val();
  };

  splitted_vals[position_of_begin_val] != "current" ?
    do_prepare_position() : do_not_prepare_position();
  what->to_position = std::stod(splitted_vals[position_of_end_val]);
  what->scans = std::stoi(splitted_vals.back());

}

std::unique_ptr<Scan2D::Input_types> Scan2D::parse_input(std::string const& input) {

  std::istringstream iss_axis(input);

  std::vector<std::string> splitted_vals{
    std::istream_iterator<std::string>{iss_axis},
    std::istream_iterator<std::string>{}
  };

  std::unique_ptr<Input_types> unique_parser;

  try {
    unique_parser = Input_Factory(splitted_vals.size(), splitted_vals.front());
    unique_parser->fill_what(splitted_vals, _coords.xyz());
  }
  catch (std::runtime_error const& e) {
    std::cout << e.what() << std::endl;
    exit(1);
  }

  return unique_parser;//std::move(x_filler);

}

std::unique_ptr<Scan2D::Input_types> Scan2D::Input_Factory(std::size_t size, std::string kind) {

  if (kind.compare("bond") == 0) {
    if (size == 6) {
      return std::make_unique<Normal_Bond_Input>(shared_from_this());
    }
  }
  else if (kind.compare("angle") == 0) {
    if (size == 7) {
      return std::make_unique<Normal_Angle_Input>(shared_from_this());
    }
  }
  else if (kind.compare("dihedral") == 0) {
    if (size == 8) {
      return std::make_unique<Normal_Dihedral_Input>(shared_from_this());
    }
  }
  throw std::runtime_error("Your Input is wrong. Excpected 'bond', 'angle' or 'dihedral' by specifieing your 2D Scan!");

}

Scan2D::length_type Scan2D::get_length(cbond const& ab) {

  auto const a_to_b = ab.a - ab.b;
  return scon::geometric_length(a_to_b);
}

Scan2D::angle_type Scan2D::get_angle(cangle const& abc) {

  auto const b_to_a = abc.b - abc.a;
  auto const b_to_c = abc.b - abc.c;
  auto const AB_x_BC = scon::cross(b_to_c, b_to_a);
  auto const AB_O_BC = scon::dot(b_to_c, b_to_a);

  return Scan2D::angle_type::from_rad(atan2(scon::geometric_length(AB_x_BC), AB_O_BC));
}

Scan2D::angle_type Scan2D::get_dihedral(cdihedral const& abcd) {

  auto const b_to_a = abcd.a - abcd.b;
  auto const b_to_c = abcd.b - abcd.c;//Shouldn't this be c-b?
  auto const b_to_d = abcd.d - abcd.b;

  auto const normal_vec_to_plain_ABC = scon::cross(b_to_c, b_to_a);
  auto const normal_vec_to_plain_BCD = scon::cross(b_to_c, b_to_d);

  auto const normal_to_normals = scon::cross(normal_vec_to_plain_BCD, normal_vec_to_plain_ABC);

  auto const orthonormal_reference = scon::normalized(b_to_c);

  auto in_degrees = Scan2D::angle_type::from_rad(
    atan2(scon::geometric_length(normal_to_normals), scon::dot(normal_vec_to_plain_BCD, normal_vec_to_plain_ABC))
  );

  if (scon::dot(orthonormal_reference, normal_to_normals) > 0.0) {
    return in_degrees;
  }
  else {
    return -in_degrees;
  }
}

coords::Cartesian_Point Scan2D::change_length_of_bond(cbond const& ab, length_type const& new_length) {

  auto direction = scon::normalized(ab.a - ab.b);

  return ab.b + direction * new_length;

}

coords::Cartesian_Point Scan2D::rotate_a_to_new_angle(cangle const& abc, Scan2D::angle_type const& new_angle) {

  cbond ab(abc.a, abc.b);

  auto radius = Scan2D::get_length(ab);
  auto sin_inclination = sin(new_angle);
  auto cos_inclination = cos(new_angle);

  auto new_x = radius * cos_inclination;
  auto new_y = radius * sin_inclination;

  auto CB = scon::normalized(abc.c - abc.b);
  auto ZxA = scon::normalized(scon::cross(CB, abc.a - abc.b));
  auto ZxAxZ = scon::normalized(scon::cross(ZxA, CB));

  auto ret = abc.b + CB * new_x + ZxAxZ * new_y;

  return ret;

}

coords::Cartesian_Point Scan2D::rotate_a_to_new_dihedral(cdihedral const& abcd, Scan2D::angle_type const& new_angle) {

  cbond ab(abcd.a, abcd.b);
  cangle abc(abcd.a, abcd.b, abcd.c);

  auto radius = Scan2D::get_length(ab);
  auto inclination = Scan2D::get_angle(abc);
  auto sin_inclination = sin(inclination);
  auto cos_inclination = cos(inclination);
  auto sin_azimuth = -sin(new_angle);
  auto cos_azimuth = cos(new_angle);

  coords::Cartesian_Point helper_point(
    radius * sin_inclination * cos_azimuth,
    radius * sin_inclination * sin_azimuth,
    radius * cos_inclination
  );

  auto BC = scon::normalized(abcd.c - abcd.b);

  auto ZxA = scon::normalized(scon::cross(BC, abcd.d - abcd.b));
  auto ZxAxZ = scon::normalized(scon::cross(ZxA, BC));

  return abcd.b + ZxAxZ * helper_point.x() + ZxA * helper_point.y() + BC * helper_point.z();

}

coords::Representation_3D Scan2D::Move_Handler::rotate_molecule_behind_a_dih(Scan2D::length_type const& deg) const {
  using std::placeholders::_1;

  auto xyz = _coords.xyz();
  auto tmp_axis = xyz[atoms[1] - 1] - xyz[atoms[2] - 1];
  scon::RotationMatrix::Vector ax{ tmp_axis.x(),tmp_axis.y(),tmp_axis.z() };
  scon::RotationMatrix::Vector center{ xyz[atoms[1] - 1].x(), xyz[atoms[1] - 1].y(), xyz[atoms[1] - 1].z() };
  //coroutine_type::pull_type source{ std::bind(&Scan2D::go_along_backbone, this, _1, abcd[0],abcd[1]) };
  auto source = go_along_backbone(atoms);

  //    auto const posssss = parser->y_parser->say_val();

  auto const max_change_rotation_rad = angle_type::from_deg(parent->max_change_rotation).radians();

  for (auto const& dd : source) {
    std::size_t atom_number, bond_count;
    std::tie(atom_number, bond_count) = dd;

    auto const change = angle_type::from_deg(deg < 0.0 ?
      deg + parent->change_from_atom_to_atom * static_cast<double>(bond_count) :
      deg - parent->change_from_atom_to_atom * static_cast<double>(bond_count)
    ).radians();
    if (fabs(change) <= max_change_rotation_rad || (deg * change) < 0.) {
      continue;
    }

    //	std::cout << atom_number << " " << bond_count << " " << deg << " " << change << " " << " " << max_change_rotation_rad << " " << parent->change_from_atom_to_atom << std::endl;

    auto rot = scon::RotationMatrix::rotate_around_axis_with_center(change, ax, center);

    auto&& coord = xyz[atom_number];
    scon::RotationMatrix::Vector tmp_coord{ coord.x(),coord.y(),coord.z() };
    tmp_coord = rot * tmp_coord;
    coord = coords::r3(tmp_coord[0], tmp_coord[1], tmp_coord[2]);
  }
  return xyz;
}


coords::Representation_3D Scan2D::Move_Handler::rotate_molecule_behind_a_ang(Scan2D::length_type const& deg) const {
  using std::placeholders::_1;

  auto xyz = _coords.xyz();

  auto ba = xyz[atoms[0u] - 1u] - xyz[atoms[1u] - 1u];
  auto bc = xyz[atoms[2u] - 1u] - xyz[atoms[1u] - 1u];
  auto const& b = xyz[atoms[1u] - 1u];

  scon::RotationMatrix::Vector center{ b.x(),b.y(),b.z() };

  auto ax = scon::RotationMatrix::Vector{ ba.x(), ba.y(), ba.z() }.cross(scon::RotationMatrix::Vector{ bc.x(),bc.y(),bc.z() }).normalized();

  //coroutine_type::pull_type source{ std::bind(&Scan2D::go_along_backbone, this, _1, abc[0], abc[1]) };
  auto source = go_along_backbone(atoms);

  auto const max_change_rotation_rad = angle_type::from_deg(parent->max_change_rotation).radians();

  for (auto const& aa : source) {
    std::size_t atom_number, bond_count;
    std::tie(atom_number, bond_count) = aa;

    length_type const change = angle_type::from_deg(deg - parent->change_from_atom_to_atom * static_cast<double>(bond_count)).radians();

    if (fabs(change) >= max_change_rotation_rad && (deg * change) < 0.) {
      continue;
    }

    auto rot = scon::RotationMatrix::rotate_around_axis_with_center(change, -ax, center);

    auto&& coord = xyz[atom_number];
    scon::RotationMatrix::Vector tmp_coord{ coord.x(), coord.y(), coord.z() };
    tmp_coord = rot * tmp_coord;
    coord = coords::r3(tmp_coord[0], tmp_coord[1], tmp_coord[2]);
  }
  return xyz;

}

coords::Representation_3D Scan2D::Move_Handler::transform_molecule_behind_a_bond(length_type const& length) const {
  using std::placeholders::_1;
  auto xyz = _coords.xyz();

  auto ba = xyz[atoms[0u] - 1u] - xyz[atoms[1u] - 1u];

  auto ax = scon::RotationMatrix::Vector{ ba.x(),ba.y(),ba.z() };

  //coroutine_type::pull_type source{std::bind(&Scan2D::go_along_backbone, this, _1, ab[0], ab[1])};
  auto source = go_along_backbone(atoms);

  for (auto const& bb : source) {
    std::size_t atom_number, bond_count;
    std::tie(atom_number, bond_count) = bb;

    auto const change = length - parent->change_from_atom_to_atom * static_cast<double>(bond_count);

    if (fabs(change) <= parent->max_change_rotation) {
      continue;
    }

    auto vec = ax.normalized() * change;

    scon::RotationMatrix::Translation trans(vec);

    auto&& coord = xyz[atom_number];
    scon::RotationMatrix::Vector tmp_coord{ coord.x(),coord.y(),coord.z() };
    tmp_coord = trans * tmp_coord;
    coord = coords::r3(tmp_coord[0], tmp_coord[1], tmp_coord[2]);
  }
  return xyz;
}

//void Scan2D::go_along_backbone(coroutine_type::push_type & sink, std::size_t const & atom, std::size_t const & border) {
//    auto const & atoms = _coords.atoms();
//    std::unordered_set<int> check_list;
//    std::size_t recursion_count = 1;
//
//    check_list.insert(atom);
//
//    std::function<void(std::vector<std::size_t>)> parse_neighbors = [&](std::vector<std::size_t> const & neigh) -> void {
//        std::vector<std::size_t> unchecked_bonds;
//        for (auto const & n : neigh) {
//            if (n == border) continue;
//            if (check_list.insert(n).second) {
//                sink(std::make_pair(n,recursion_count));
//                unchecked_bonds.emplace_back(n);
//            }
//        }
//        ++recursion_count;
//        parse_neighbors(unchecked_bonds);
//    };
//
//  parse_neighbors(atoms.atom(atom-1).bonds());
//
//}

Scan2D::bond_set Scan2D::Move_Handler::go_along_backbone(std::vector<std::size_t> const& kind) const {
  auto const& atoms = _coords.atoms();
  bond_set ret;
  std::size_t recursion_count = 0;

  std::size_t atom = 0u;
  std::size_t border = 0u;

  std::tie(atom, border) = kind.size() == 2 ?
    std::make_pair(kind[0u] - 1u, kind[1u] - 1u) :
    std::make_pair(kind[1u] - 1u, kind[2u] - 1u);

  ret.insert(std::make_pair(atom, recursion_count));

  std::function<void(std::vector<std::size_t>)> parse_neighbors = [&](std::vector<std::size_t> const& neigh) -> void {
    if (neigh.size() <= 1) return;
    std::vector<std::size_t> unchecked_bonds;

    auto insert_element = [&](std::size_t const n, std::size_t const m) {
      if (ret.insert(std::make_pair(n, m)).second) {
        auto const& atoms_ref = atoms.atom(n).bonds();
        unchecked_bonds.insert(unchecked_bonds.cend(), atoms_ref.cbegin(), atoms_ref.cend());
      }
    };

    for (auto const& n : neigh) {
      if (n == border) continue;
      /*if(this->_coords.atoms().atom(n).fixed() && n != atom){
        fixed_depth = fixed_depth == 0 ? recursion_count : fixed_depth;
        insert_element(n, fixed_depth);
      }*/
      insert_element(n, recursion_count);
    }
    ++recursion_count;
    parse_neighbors(unchecked_bonds);
    --recursion_count;
  };

  parse_neighbors(atoms.atom(atom).bonds());

  return ret;
}

void Scan2D::make_scan() {

  prepare_scan();

  coords::output::formats::tinker output(_coords);
  //std::cout << parser->x_parser->say_val() << " " << parser->y_parser->say_val() << std::endl;
  //parser->fix_atoms(_coords);
  //write_energy_entry(optimize(_coords));
  //parser->x_parser->set_coords(_coords.xyz());

  //output.to_stream(logfile);

  //go_along_y_axis(_coords);

  for (auto const& x_step : axis->x_steps) {    // for each step on x axis
    
    if (Config::get().general.verbosity > 0) {
      std::cout << "X: " << x_step <<", Y: "<<axis->y_steps[0]<< std::endl;
    }

    auto const& x_atoms = parser->x_parser->what->atoms;
    auto const& y_atoms = parser->y_parser->what->atoms;

    // set current x value (this is important!)
    Move_Handler mh_x(_coords, x_atoms, shared_from_this());
    mh_x.set_new_pos(x_step);

    parser->x_parser->set_coords(_coords.xyz());

    _coords.set_xyz(
      parser->x_parser->make_move(mh_x),  
      true
    );
    parser->x_parser->set_coords(_coords.xyz());

    // set y back to first value (this is only for small changes during optimization)
    Move_Handler mh_y(_coords, y_atoms, shared_from_this());
    mh_y.set_new_pos(axis->y_steps[0]);

    parser->y_parser->set_coords(_coords.xyz());

    _coords.set_xyz(
      parser->y_parser->make_move(mh_y),
      true
    );
    parser->y_parser->set_coords(_coords.xyz());

    // print output before optimization
    before << output << std::flush;

    // set constraints
    if (Config::get().optimization.local.method == config::optimization_conf::lo_types::LBFGS) {
      parser->fix_atoms(_coords);
    } 
    else if (Config::get().optimization.local.method == config::optimization_conf::lo_types::OPTPP) {
      parser->set_constraints(x_step, axis->y_steps[0]);   // set distance constraints
    }

    // perform optimization for first step on y axis
    write_energy_entry(optimize(_coords));  

    // print output after optimization
    logfile << output << std::flush;

    // now do the rest of the steps on y axis
    go_along_y_axis(_coords, x_step);
  }

}

length_type Scan2D::optimize(coords::Coordinates& c) {
  auto E_o = 0.0;
  if (Config::get().scan2d.fixed_scan) E_o = c.e();
  else E_o = c.o();          //<- SUPPPPPPER WICHTIG!!!!!!!!!!!
  parser->x_parser->set_coords(c.xyz());
  parser->y_parser->set_coords(c.xyz());

  // save output files of OPT++
  if (Config::get().optimization.local.method == config::optimization_conf::lo_types::OPTPP) {
    if (Config::get().general.verbosity > 3) {
      auto dist_x = Config::get().optimization.local.optpp_conf.constraints[0].distance;
      auto dist_y = Config::get().optimization.local.optpp_conf.constraints[1].distance;
      std::string new_name = "OPT_" + convert_number_to_string_without_dots(dist_x) + "_" + convert_number_to_string_without_dots(dist_y) + ".out";
      std::rename("OPT_DEFAULT.out", new_name.c_str());
      if (file_exists("trace.arc")) {
        new_name = "trace_" + convert_number_to_string_without_dots(dist_x) + "_" + convert_number_to_string_without_dots(dist_y) + ".arc";
        std::rename("trace.arc", new_name.c_str());
      }
    }
  }
  return E_o;
}

void Scan2D::prepare_scan() {
  auto xyz = _coords.xyz();

  auto const x_move = parser->x_parser->what->from_position;
  auto const y_move = parser->y_parser->what->from_position;
  auto const& x_atoms = parser->x_parser->what->atoms;
  auto const& y_atoms = parser->y_parser->what->atoms;

  parser->x_parser->set_coords(_coords.xyz());

  Move_Handler xmh(_coords, x_atoms, shared_from_this());
  xmh.set_new_pos(x_move);

  _coords.set_xyz(
    parser->x_parser->make_move(xmh),
    true
  );

  parser->y_parser->set_coords(_coords.xyz());

  Move_Handler ymh(_coords, y_atoms, shared_from_this());
  ymh.set_new_pos(y_move);

  _coords.set_xyz(
    parser->y_parser->make_move(std::move(ymh)),
    true
  );
  parser->x_parser->set_coords(_coords.xyz());
  parser->y_parser->set_coords(_coords.xyz());
}

void Scan2D::go_along_y_axis(coords::Coordinates coords, double const x_step) {

  coords::output::formats::tinker output(coords);
  parser->y_parser->set_coords(coords.xyz());

  std::for_each(axis->y_steps.cbegin() + 1, axis->y_steps.cend(), [&](auto&& y_step) {  // for each step on x axis
    
    if (Config::get().general.verbosity > 0) {
      std::cout <<"X: "<<x_step<< ", Y: " << y_step << std::endl;
    }

    auto const& y_atoms = parser->y_parser->what->atoms;
    auto const& x_atoms = parser->x_parser->what->atoms;

    // move y to current value (this is important!)
    Move_Handler mh_y(coords, y_atoms, this->shared_from_this());
    mh_y.set_new_pos(y_step);

    this->parser->y_parser->set_coords(coords.xyz());

    coords.set_xyz(
      parser->y_parser->make_move(mh_y),
      true
    );

    // move x back to current value (this is only for small changes during optimization)
    Move_Handler mh_x(coords, x_atoms, this->shared_from_this());
    mh_x.set_new_pos(x_step);

    this->parser->x_parser->set_coords(coords.xyz());

    coords.set_xyz(
      parser->x_parser->make_move(mh_x),
      true
    );

    // set constraints
    if (Config::get().optimization.local.method == config::optimization_conf::lo_types::LBFGS) {
      parser->fix_atoms(coords);
    }
    else if (Config::get().optimization.local.method == config::optimization_conf::lo_types::OPTPP) {
      parser->set_constraints(x_step, y_step);  // set distance constraints
    }

    // print output before optimization
    before << output << std::flush;

    // perform optimization
    this->write_energy_entry(this->optimize(coords));  

    // save optimized coordinates
    parser->y_parser->set_coords(coords.xyz());
    parser->x_parser->set_coords(coords.xyz());

    // print output after optimization
    logfile << output << std::flush;
  });
  energies << "\n";

}

void Scan2D::Normal_Bond_Input::set_coords(coords::Representation_3D const& xyz) {
  auto& atoms = what->atoms;
  bond = std::make_unique<Scan2D::cbond>(xyz[atoms[0] - 1u], xyz[atoms[1] - 1u]);
}

void Scan2D::Normal_Angle_Input::set_coords(coords::Representation_3D const& xyz) {
  auto& atoms = what->atoms;
  angle = std::make_unique<Scan2D::cangle>(xyz[atoms[0] - 1u], xyz[atoms[1] - 1u], xyz[atoms[2] - 1u]);
}

void Scan2D::Normal_Dihedral_Input::set_coords(coords::Representation_3D const& xyz) {
  auto& atoms = what->atoms;
  dihedral = std::make_unique<Scan2D::cdihedral>(xyz[atoms[0] - 1u], xyz[atoms[1] - 1u], xyz[atoms[2] - 1u], xyz[atoms[3] - 1u]);
}

std::vector<Scan2D::length_type> Scan2D::Normal_Bond_Input::make_axis() {

  auto current_position{ what->from_position };
  auto step_width{ (what->to_position - current_position) / static_cast<length_type>(what->scans - 1u) };
  auto made_vec = std::vector<length_type>(what->scans);

  for (auto&& el : made_vec) {
    el = current_position;
    current_position += step_width;
  }

  return made_vec;
}

std::vector<Scan2D::length_type> Scan2D::Normal_Angle_Input::make_axis() {

  auto current_inclination{ what->from_position };
  auto step_width{ (what->to_position - current_inclination) / static_cast<length_type>(what->scans - 1u) };
  auto made_vec = std::vector<length_type>(what->scans);

  for (auto&& el : made_vec) {
    el = current_inclination;
    current_inclination += step_width;
  }

  return made_vec;

}

std::vector<Scan2D::length_type> Scan2D::Normal_Dihedral_Input::make_axis() {

  auto current_azimuth{ what->from_position };
  auto step_width{ (what->to_position - current_azimuth) / static_cast<length_type>(what->scans - 1u) };
  auto made_vec = std::vector<length_type>(what->scans);

  for (auto&& el : made_vec) {
    el = current_azimuth;
    current_azimuth += step_width;
  }

  return made_vec;

}

//coords::Representation_3D Scan2D::Normal_Bond_Input::make_move(Scan2D::Move_Handler const & mh) {
//    auto p = parent.lock();
//
//      auto new_molecule = p->_coords.xyz();
//      new_molecule[mh.atoms.at(0) - 1u] = change_length_of_bond(*bond, mh.new_pos);
//      return new_molecule;
//}

//This needs to be implemented rather than the other make_move methods
coords::Representation_3D Scan2D::Normal_Bond_Input::make_move(Scan2D::Move_Handler const& mh) {
  auto p = parent.lock();

  auto const change = mh.new_pos - say_val();

  if (fabs(change) > p->max_change_rotation) {
    return mh.transform_molecule_behind_a_bond(change);
  }
  else {
    auto new_molecule = mh._coords.xyz();
    new_molecule[mh.atoms.at(0) - 1u] = change_length_of_bond(*bond, mh.new_pos);
    return new_molecule;
  }
}

//coords::Representation_3D Scan2D::Normal_Angle_Input::make_move(Scan2D::Move_Handler const & mh) {
//    auto p = parent.lock();
//
//    auto new_molecule = mh._coords.xyz();
//    new_molecule[mh.atoms.at(0) - 1u] = rotate_a_to_new_angle(*angle, angle_type::from_deg(mh.new_pos));
//    return new_molecule;
//}

coords::Representation_3D Scan2D::Normal_Angle_Input::make_move(Scan2D::Move_Handler const& mh) {
  auto p = parent.lock();
  auto const change = mh.new_pos - say_val();

  if (fabs(change) > p->max_change_rotation) {
    return mh.rotate_molecule_behind_a_ang(change);
  }
  else {
    auto new_molecule = mh._coords.xyz();
    new_molecule[mh.atoms.at(0) - 1u] = rotate_a_to_new_angle(*angle, angle_type::from_deg(mh.new_pos));
    return new_molecule;
  }
}

//coords::Representation_3D Scan2D::Normal_Dihedral_Input::make_move(Scan2D::Move_Handler const & mh) {
//    auto p = parent.lock();
//
//    auto new_molecule = mh._coords.xyz();
//    new_molecule[mh.atoms.at(0) - 1u] = rotate_a_to_new_dihedral(*dihedral, angle_type::from_deg(mh.new_pos));
//    return new_molecule;
//}

coords::Representation_3D Scan2D::Normal_Dihedral_Input::make_move(Scan2D::Move_Handler const& mh) {
  auto p = parent.lock();
  auto change = mh.new_pos - say_val();
  if (change < 0.) change = 360. + change;  // change should always be a positive value between 0 and 360
  if (fabs(change) > p->max_change_rotation) {
    return mh.rotate_molecule_behind_a_dih(change);
  }
  else {
    auto new_molecule = mh._coords.xyz();
    new_molecule[mh.atoms.at(0) - 1u] = rotate_a_to_new_dihedral(*dihedral, angle_type::from_deg(mh.new_pos));
    return new_molecule;
  }
}

Scan2D::length_type Scan2D::Normal_Bond_Input::say_val() {
  return Scan2D::get_length(*bond);
}

Scan2D::length_type Scan2D::Normal_Angle_Input::say_val() {
  return Scan2D::get_angle(*angle).degrees();
}

Scan2D::length_type Scan2D::Normal_Dihedral_Input::say_val() {
  return Scan2D::get_dihedral(*dihedral).degrees();
}

void Scan2D::XY_Parser::fix_atoms(coords::Coordinates& coords) const 
{
  for (auto&& atom : x_parser->what->atoms) {
    coords.fix(atom - 1);
  }
  for (auto&& atom : y_parser->what->atoms) {
    coords.fix(atom - 1);
  }
}

void Scan2D::XY_Parser::set_constraints(coords::float_type const x_constr, coords::float_type const y_constr) const 
{
  if(x_parser->what->atoms.size() != 2 || y_parser->what->atoms.size() != 2){
    throw std::runtime_error("It is not possible to set constraints except distances (i. e. 2 atoms)");
  }
  
  // remove old constraints
  Config::set().optimization.local.optpp_conf.constraints.clear();

  // set constraint in x-direction
  std::size_t a = x_parser->what->atoms[0];
  std::size_t b = x_parser->what->atoms[1];
  config::optimization_conf::constraint_bond x_constraint{a-1,b-1,x_constr};
  Config::set().optimization.local.optpp_conf.constraints.emplace_back(x_constraint);
  // set constraint in y-direction
  std::size_t c = y_parser->what->atoms[0];
  std::size_t d = y_parser->what->atoms[1];
  config::optimization_conf::constraint_bond y_constraint{c-1,d-1,y_constr};
  Config::set().optimization.local.optpp_conf.constraints.emplace_back(y_constraint);
  
  // console output
  if (Config::get().general.verbosity > 3)
  {
    std::cout<<"Setting the following constraints:\n";
    std::cout<<a-1<<" , "<<b-1<<": "<<x_constr<<"\n";
    std::cout<<c-1<<" , "<<d-1<<": "<<y_constr<<"\n";
  }
}
