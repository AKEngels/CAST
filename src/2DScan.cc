#include "2DScan.h"

Scan2D::Scan2D(coords::Coordinates & coords) 
  : _coords(coords), change_from_atom_to_atom(Config::get().scan2d.change_from_atom_to_atom), max_change_rotation(Config::get().scan2d.max_change_to_rotate_whole_molecule)
{
	logfile.open(structures_file);
	energies.open(energie_file);
}

void Scan2D::execute_scan(){
	auto & both_whats = Config::get().scan2d.AXES;

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

        axis = std::make_unique<XY_steps>(x_changes, y_changes);

        make_scan();
}

Scan2D::~Scan2D() {
	//Seems like it is empty ;)
}

void Scan2D::Normal_Input::fill_what(std::vector<std::string> & splitted_vals, coords::Representation_3D const & xyz) {

	what = std::make_unique<Scan2D::what>();

	what->what_kind = splitted_vals.front();

	std::transform(splitted_vals.begin() + 1, splitted_vals.end() - 3, std::back_inserter(what->atoms),
		[](std::string const & atom) { return std::stoi(atom); }
	);

	set_coords(xyz);

	auto position_of_begin_val = splitted_vals.size() - 3;
	auto position_of_end_val = splitted_vals.size() - 2;

	splitted_vals[position_of_begin_val] != "current" ?
		what->from_position = std::stod(splitted_vals[position_of_begin_val]) : what->from_position = say_val();
	what->to_position = std::stod(splitted_vals[position_of_end_val]);
	what->scans = std::stoi(splitted_vals.back());

}

std::unique_ptr<Scan2D::Input_types> Scan2D::parse_input(std::string const & input) {
	
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
	catch (std::runtime_error e) {
		std::cout << e.what() << std::endl;
		exit(1);
	}

	return unique_parser;//std::move(x_filler);

}

std::unique_ptr<Scan2D::Input_types> Scan2D::Input_Factory(std::size_t size, std::string kind) {

	if(kind.compare("bond") == 0) {
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

Scan2D::length_type Scan2D::get_length(cbond const & ab) {

	auto const a_to_b = ab.a - ab.b;
	return scon::geometric_length(a_to_b);
}

Scan2D::angle_type Scan2D::get_angle(cangle const & abc) {
	
	auto const b_to_a = abc.b - abc.a;
	auto const b_to_c = abc.b - abc.c;
	auto const AB_x_BC = scon::cross(b_to_c, b_to_a);
	auto const AB_O_BC = scon::dot(b_to_c, b_to_a);

	return Scan2D::angle_type::from_rad(atan2(scon::geometric_length(AB_x_BC), AB_O_BC));
}

Scan2D::angle_type Scan2D::get_dihedral(cdihedral const & abcd) {
	
	auto const b_to_a = abcd.a - abcd.b;
	auto const b_to_c = abcd.b - abcd.c;
	auto const b_to_d = abcd.d - abcd.b;

	auto const normal_vec_to_plain_ABC = scon::cross(b_to_c, b_to_a);
	auto const normal_vec_to_plain_BCD = scon::cross(b_to_c, b_to_d);

	auto const notmal_to_normals = scon::cross(normal_vec_to_plain_BCD, normal_vec_to_plain_ABC);

	auto const orthonormal_reference = scon::normalized(b_to_c);

	auto in_degrees = Scan2D::angle_type::from_rad(
		atan2(scon::geometric_length(notmal_to_normals), scon::dot(normal_vec_to_plain_BCD, normal_vec_to_plain_ABC))
	);

	if (scon::dot(orthonormal_reference, notmal_to_normals) > 0.0) {
		return in_degrees;
	}
	else {
		return -in_degrees;
	}
}

coords::Cartesian_Point Scan2D::change_length_of_bond(cbond const & ab, length_type const & new_length) {

	auto direction = scon::normalized(ab.a - ab.b);

	return ab.b + direction * new_length;

}

coords::Cartesian_Point Scan2D::rotate_a_to_new_angle(cangle const & abc, Scan2D::angle_type const & new_angle) {

	cbond ab(abc.a, abc.b);
	
	auto radius = Scan2D::get_length(ab);
	auto sin_inclination = sin(new_angle);
	auto cos_inclination = cos(new_angle);
	
	auto new_x = radius * cos_inclination;
	auto new_y = radius * sin_inclination;
	
	auto CB = scon::normalized(abc.c - abc.b);
	auto ZxA = scon::normalized(scon::cross(CB, abc.a - abc.b));
	auto ZxAxZ = scon::normalized(scon::cross(ZxA, CB));

	return abc.b + CB * new_x + ZxAxZ * new_y;

}

coords::Cartesian_Point Scan2D::rotate_a_to_new_dihedral(cdihedral const & abcd, Scan2D::angle_type const & new_angle) {

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

coords::Representation_3D Scan2D::rotate_molecule_behind_a_dih(std::vector<std::size_t> const & abcd, Scan2D::length_type const & deg) {
  using std::placeholders::_1;
    auto xyz = _coords.xyz(); 
    auto tmp_axis = xyz[abcd[2] - 1] - xyz[abcd[1] - 1];
    RotationMatrix::Vector axis{ tmp_axis.x(),tmp_axis.y(),tmp_axis.z() };
    RotationMatrix::Vector center{ xyz[abcd[1] - 1].x(), xyz[abcd[1] - 1].y(), xyz[abcd[1] - 1].z() };
    //coroutine_type::pull_type source{ std::bind(&Scan2D::go_along_backbone, this, _1, abcd[0],abcd[1]) };
    auto source = go_along_backbone(abcd[0], abcd[1]);

    for (auto const & dd : source) {
        std::size_t atom_number, bond_count;
        std::tie(atom_number, bond_count) = dd;

        length_type const rad = angle_type::from_deg(deg - change_from_atom_to_atom*static_cast<double>(bond_count)).radians();

        if (rad <= 0.0) {
            break;
        }

        auto rot = RotationMatrix::rotate_around_axis_with_center(rad, axis, center);

        auto && coord = xyz[atom_number - 1];
        RotationMatrix::Vector tmp_coord{ coord.x(),coord.y(),coord.z() };
        tmp_coord = rot*tmp_coord;
        coord = coords::r3(tmp_coord[0], tmp_coord[1], tmp_coord[2]);
    }
    return xyz;
}


coords::Representation_3D Scan2D::rotate_molecule_behind_a_ang(std::vector<std::size_t> const & abc, Scan2D::length_type const & deg) {
  using std::placeholders::_1;
  
  auto xyz = _coords.xyz();

  auto ba = xyz[abc[0] - 1] - xyz[abc[1] - 1];
  auto bc = xyz[abc[2] - 1] - xyz[abc[1] - 1];
  auto const & b = xyz[abc[1] - 1];

  RotationMatrix::Vector center{ b.x(),b.y(),b.z() };

  RotationMatrix::Vector axis = RotationMatrix::Vector{ ba.x(), ba.y(), ba.z() }.cross(RotationMatrix::Vector{bc.x(),bc.y(),bc.z()}).normalized();

  //coroutine_type::pull_type source{ std::bind(&Scan2D::go_along_backbone, this, _1, abc[0], abc[1]) };
  auto source = go_along_backbone(abc[0], abc[1]);

  for (auto const & aa : source) {
    std::size_t atom_number, bond_count;
    std::tie(atom_number, bond_count) = aa;

    length_type const rad = angle_type::from_deg(deg - change_from_atom_to_atom*static_cast<double>(bond_count)).radians();

    if (rad <= 0.) {
      break;
    }

    auto rot = RotationMatrix::rotate_around_axis_with_center(rad, axis, center);

    auto && coord = xyz[atom_number - 1];
    RotationMatrix::Vector tmp_coord{ coord.x(),coord.y(),coord.z() };
    tmp_coord = rot*tmp_coord;
    coord = coords::r3(tmp_coord[0], tmp_coord[1], tmp_coord[2]);
  }
  return xyz;

}

coords::Representation_3D Scan2D::transform_molecule_behind_a_bond(std::vector<std::size_t> const & ab, length_type const & length) {
  using std::placeholders::_1;
  auto xyz = _coords.xyz();

  auto ba = xyz[ab[0] - 1] - xyz[ab[1] - 1];

  auto axis = RotationMatrix::Vector{ ba.x(),ba.y(),ba.z() };
  auto distance = axis.norm();

  //coroutine_type::pull_type source{std::bind(&Scan2D::go_along_backbone, this, _1, ab[0], ab[1])};
  auto source = go_along_backbone(ab[0], ab[1]);

    for (auto const & bb : source) {
      std::size_t atom_number, bond_count;
      std::tie(atom_number, bond_count) = bb;

      auto change = (length - change_from_atom_to_atom*static_cast<double>(bond_count));

      if (change <= 0.) {
        break;
      }

      change -= distance;

      auto vec = axis.normalized() * change;

      RotationMatrix::Translation trans(vec);

      auto && coord = xyz[atom_number - 1];
      RotationMatrix::Vector tmp_coord{ coord.x(),coord.y(),coord.z() };
      tmp_coord = trans*tmp_coord;
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

Scan2D::bond_set Scan2D::go_along_backbone(std::size_t const & atom, std::size_t const & border) {
  auto const & atoms = _coords.atoms();
  bond_set ret;
  std::size_t recursion_count = 1;

  ret.insert(std::make_pair(atom, 0));

  std::function<void(std::vector<std::size_t>)> parse_neighbors = [&](std::vector<std::size_t> const & neigh) -> void {
    if(neigh.size()==1) return;
    for (auto const & n : neigh) {
      if (n == border) continue;
      if (ret.insert(std::make_pair(n, recursion_count)).second) {
        ++recursion_count;
        parse_neighbors(atoms.atom(n - 1).bonds());
        --recursion_count;
      }
    }
  };

  parse_neighbors(atoms.atom(atom - 1).bonds());

  return ret;
}

void Scan2D::make_scan() {

	prepare_scan();
	
	coords::output::formats::tinker output(_coords);
	parser->x_parser->set_coords(_coords.xyz());

	for (auto && x_step : axis->x_steps) {

		++x_circle;

		_coords.set_xyz(parser->x_parser->make_move(x_step,
          parser->x_parser->what->atoms
        ), true);
		parser->fix_atoms(_coords);

		write_energy_entry(_coords.o());
		parser->x_parser->set_coords(_coords.xyz());

		output.to_stream(logfile);

		go_along_y_axis(_coords);
		
	}

}

void Scan2D::prepare_scan() {
	auto xyz = _coords.xyz();

	auto & x_move = parser->x_parser->what->from_position;
	auto & y_move = parser->y_parser->what->from_position;

    _coords.set_xyz(parser->x_parser->make_move(x_move,
      parser->x_parser->what->atoms
    ),true);
    _coords.set_xyz(parser->y_parser->make_move(y_move,
      parser->x_parser->what->atoms
    ),true);

}

void Scan2D::go_along_y_axis(coords::Coordinates coords) {

	coords::output::formats::tinker output(coords);
	parser->y_parser->set_coords(coords.xyz());

	std::for_each(axis->y_steps.cbegin()+1, axis->y_steps.cend(), [&](auto && y_step){

		++y_circle;

		coords.set_xyz(parser->y_parser->make_move(y_step,
          parser->y_parser->what->atoms
        ), true);
		parser->fix_atoms(coords);

		this->write_energy_entry(coords.o());
		parser->y_parser->set_coords(coords.xyz());

		output.to_stream(logfile);

	});

	y_circle = 0;

}

void Scan2D::Normal_Bond_Input::set_coords(coords::Representation_3D const & xyz) {
	auto & atoms = what->atoms;
	bond = std::make_unique<Scan2D::cbond>(xyz[atoms[0] - 1u], xyz[atoms[1] - 1u]);
}

void Scan2D::Normal_Angle_Input::set_coords(coords::Representation_3D const & xyz) {
	auto & atoms = what->atoms;
	angle = std::make_unique<Scan2D::cangle>(xyz[atoms[0] - 1u], xyz[atoms[1] - 1u], xyz[atoms[2] - 1u]);
}

void Scan2D::Normal_Dihedral_Input::set_coords(coords::Representation_3D const & xyz) {
	auto & atoms = what->atoms;
	dihedral = std::make_unique<Scan2D::cdihedral>(xyz[atoms[0] - 1u], xyz[atoms[1] - 1u], xyz[atoms[2] - 1u], xyz[atoms[3] - 1u]);
}

std::vector<Scan2D::length_type> Scan2D::Normal_Bond_Input::make_axis() {

	auto current_position{ what->from_position };
	auto step_width{(what->to_position-current_position)/static_cast<length_type>(what->scans-1u)};
	auto made_vec = std::vector<length_type>(what->scans);

	for (auto && el : made_vec) {
		el = current_position;
		current_position += step_width;
	}

	return made_vec;
}

std::vector<Scan2D::length_type> Scan2D::Normal_Angle_Input::make_axis() {

	auto current_inclination{ what->from_position };
	auto step_width{ (what->to_position - current_inclination) / static_cast<length_type>(what->scans - 1u) };
	auto made_vec = std::vector<length_type>(what->scans);

	for (auto && el : made_vec) {
		el = current_inclination;
		current_inclination += step_width;
	}

	return made_vec;

}

std::vector<Scan2D::length_type> Scan2D::Normal_Dihedral_Input::make_axis() {

	auto current_azimuth{ what->from_position };
	auto step_width{ (what->to_position - current_azimuth) / static_cast<length_type>(what->scans - 1u) };
	auto made_vec = std::vector<length_type>(what->scans);

	for (auto && el : made_vec) {
		el = current_azimuth;
		current_azimuth += step_width;
	}

	return made_vec;

}

coords::Representation_3D Scan2D::Normal_Bond_Input::make_move(length_type const & new_pos, std::vector<std::size_t> const & atoms) {
    auto p = parent.lock();
    
    if (fabs(new_pos - say_val()) > p->max_change_rotation) {
      return p->transform_molecule_behind_a_bond(atoms, new_pos);
    }
    else {
      auto new_molecule = p->_coords.xyz();
      new_molecule[atoms.at(0)] = change_length_of_bond(*bond, new_pos);
      return new_molecule;
    }

	auto new_coord = Scan2D::change_length_of_bond(*bond, new_pos);
}

coords::Representation_3D Scan2D::Normal_Angle_Input::make_move(length_type const & new_pos, std::vector<std::size_t> const & atoms) {
    auto p = parent.lock();

    if (fabs(new_pos - say_val()) > p->max_change_rotation) {
        return p->rotate_molecule_behind_a_ang(atoms, new_pos);
    }
    else {
        auto new_molecule = p->_coords.xyz();
        new_molecule[atoms.at(0)] = rotate_a_to_new_angle(*angle, angle_type::from_deg(new_pos));
        return new_molecule;
    }
}

coords::Representation_3D Scan2D::Normal_Dihedral_Input::make_move(length_type const & new_pos, std::vector<std::size_t> const & atoms) {
    auto p = parent.lock();
    
    if (fabs(new_pos - say_val()) > p->max_change_rotation) {
       return p->rotate_molecule_behind_a_dih(atoms, new_pos);
    }
    else {
       auto new_molecule = p->_coords.xyz();
       new_molecule[atoms.at(0)] = rotate_a_to_new_dihedral(*dihedral, angle_type::from_deg(new_pos));
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

void Scan2D::XY_Parser::fix_atoms(coords::Coordinates & coords)const{
	
	for (auto && atom : x_parser->what->atoms) {
		coords.fix(atom-1);
	}
	for (auto && atom : y_parser->what->atoms) {
		coords.fix(atom-1);
	}
}
