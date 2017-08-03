#include "2DScan.h"

Scan2D::Scan2D(coords::Coordinates & coords) : _coords(coords) {
	auto & both_whats = Config::get().scan2d.AXES;

	if (both_whats.size() > 2) {
		throw std::runtime_error("You can't pass more than two axis!");
	}

	logfile.open(structures_file);
	energies.open(energie_file);

	auto x_input_parser = parse_input(both_whats.front());
	auto y_input_parser = parse_input(both_whats.back());

	x_input_parser->set_coords(_coords.xyz());
	y_input_parser->set_coords(_coords.xyz());

	auto x_changes = x_input_parser->make_axis();
	auto y_changes = y_input_parser->make_axis();

	XY_Parser parser(std::move(x_input_parser), std::move(y_input_parser));

	XY_steps axis(x_changes, y_changes);

	make_scan(parser, axis);

}

Scan2D::~Scan2D() {
	logfile.close();
	energies.close();
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

	splitted_vals[position_of_begin_val].compare("current") != 0 ?
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
			return std::make_unique<Normal_Bond_Input>();
		}
	}
	else if (kind.compare("angle") == 0) {
		if (size == 7) {
			return std::make_unique<Normal_Angle_Input>();
		}
	}
	else if (kind.compare("dihedral") == 0) {
		if (size == 8) {
			return std::make_unique<Normal_Dihedral_Input>();
		}
	}
	throw std::runtime_error("Your Input is wrong. Excpected 'bond', 'angle' or 'dihedral' by specifieing your 2D Scan!");

}

Scan2D::length_type Scan2D::get_length(Scan2D::bond const & ab) {

	auto const a_to_b = ab.a - ab.b;
	return scon::geometric_length(a_to_b);
}

Scan2D::angle_type Scan2D::get_angle(Scan2D::angle const & abc) {
	
	auto const b_to_a = abc.b - abc.a;
	auto const b_to_c = abc.b - abc.c;
	auto const AB_x_BC = scon::cross(b_to_c, b_to_a);
	auto const AB_O_BC = scon::dot(b_to_c, b_to_a);

	return Scan2D::angle_type::from_rad(atan2(scon::geometric_length(AB_x_BC), AB_O_BC));
}

Scan2D::angle_type Scan2D::get_dihedral(Scan2D::dihedral const & abcd) {
	
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

coords::Cartesian_Point Scan2D::change_length_of_bond(Scan2D::bond const & ab, length_type const & new_length) {

	auto direction = scon::normalized(ab.a - ab.b);

	return ab.b + direction * new_length;

}

coords::Cartesian_Point Scan2D::rotate_a_to_new_angle(Scan2D::angle const & abc, Scan2D::angle_type const & new_angle) {

	bond ab(abc.a, abc.b);
	
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

coords::Cartesian_Point Scan2D::rotate_a_to_new_dihedral(Scan2D::dihedral const & abcd, Scan2D::angle_type const & new_angle) {

	bond ab(abcd.a, abcd.b);
	angle abc(abcd.a, abcd.b, abcd.c);

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

void Scan2D::make_scan(Scan2D::XY_Parser const & parser, Scan2D::XY_steps const & steps) {

	prepare_scan(parser);
	
	coords::output::formats::tinker output(_coords);
	parser.x_parser->set_coords(_coords.xyz());

	for (auto && x_step : steps.x_steps) {

		++x_circle;

		auto new_xyz = _coords.xyz();

		auto atom_to_change = parser.x_parser->what->atoms[0];
		new_xyz[atom_to_change - 1u] = parser.x_parser->make_move(x_step);

		_coords.set_xyz(new_xyz, true);

		parser.fix_atoms(_coords);

		_coords.o();
		parser.x_parser->set_coords(_coords.xyz());

		go_along_y_axis(parser, steps.y_steps, _coords);
		
		output.to_stream(logfile);
	}

}

void Scan2D::prepare_scan(Scan2D::XY_Parser const & parser) {
	auto xyz = _coords.xyz();

	auto & x_atom = parser.x_parser->what->atoms[0];
	auto & y_atom = parser.y_parser->what->atoms[0];
	auto & x_move = parser.x_parser->what->from_position;
	auto & y_move = parser.y_parser->what->from_position;

	xyz[x_atom - 1] = parser.x_parser->make_move(x_move);
	xyz[y_atom - 1] = parser.y_parser->make_move(y_move);

	_coords.set_xyz(xyz, true);

}

void Scan2D::go_along_y_axis(Scan2D::XY_Parser const & parser, std::vector<length_type> const & y_steps, coords::Coordinates coords) {

	coords::output::formats::tinker output(coords);
	parser.y_parser->set_coords(coords.xyz());

	

	std::for_each(y_steps.cbegin()+1, y_steps.cend(), [&](auto && y_step){

		++y_circle;

		/*std::cout << "The " << x_circle << ". x step and the " <<
			y_circle << ". y step." << std::endl;*/

		energies <<
			parser.x_parser->say_val() << " " <<
			parser.y_parser->say_val() << " " <<
			coords.e() << "\n";

		auto new_xyz = coords.xyz();
		auto atom_to_change = parser.y_parser->what->atoms[0];

		new_xyz[atom_to_change - 1u] = parser.y_parser->make_move(y_step);
		coords.set_xyz(new_xyz, true);
		parser.fix_atoms(coords);

		coords.o(); 
		parser.y_parser->set_coords(coords.xyz());

		output.to_stream(logfile);

	});

	y_circle = 0;

	energies <<
		parser.x_parser->say_val() << " " <<
		parser.y_parser->say_val() << " " <<
		coords.e() << "\n" << std::endl;

}

void Scan2D::Normal_Bond_Input::set_coords(coords::Representation_3D const & xyz) {
	auto & atoms = what->atoms;
	bond = std::make_unique<Scan2D::bond>(xyz[atoms[0] - 1u], xyz[atoms[1] - 1u]);
}

void Scan2D::Normal_Angle_Input::set_coords(coords::Representation_3D const & xyz) {
	auto & atoms = what->atoms;
	angle = std::make_unique<Scan2D::angle>(xyz[atoms[0] - 1u], xyz[atoms[1] - 1u], xyz[atoms[2] - 1u]);
}

void Scan2D::Normal_Dihedral_Input::set_coords(coords::Representation_3D const & xyz) {
	auto & atoms = what->atoms;
	dihedral = std::make_unique<Scan2D::dihedral>(xyz[atoms[0] - 1u], xyz[atoms[1] - 1u], xyz[atoms[2] - 1u], xyz[atoms[3] - 1u]);
}

std::vector<Scan2D::length_type> Scan2D::Normal_Bond_Input::make_axis() {

	auto current_position{ what->from_position };
	auto step_width{(what->to_position-current_position)/static_cast<length_type>(what->scans-1u)};
	auto made_vec = std::vector<length_type>(what->scans);

	for (auto && el : made_vec) {
		el = current_position;
		current_position += step_width;
	}

	return std::move(made_vec);
}

std::vector<Scan2D::length_type> Scan2D::Normal_Angle_Input::make_axis() {

	auto current_inclination{ what->from_position };
	auto step_width{ (what->to_position - current_inclination) / static_cast<length_type>(what->scans - 1u) };
	auto made_vec = std::vector<length_type>(what->scans);

	for (auto && el : made_vec) {
		el = current_inclination;
		current_inclination += step_width;
	}

	return std::move(made_vec);

}

std::vector<Scan2D::length_type> Scan2D::Normal_Dihedral_Input::make_axis() {

	auto current_azimuth{ what->from_position };
	auto step_width{ (what->to_position - current_azimuth) / static_cast<length_type>(what->scans - 1u) };
	auto made_vec = std::vector<length_type>(what->scans);

	for (auto && el : made_vec) {
		el = current_azimuth;
		current_azimuth += step_width;
	}

	return std::move(made_vec);

}

coords::Cartesian_Point Scan2D::Normal_Bond_Input::make_move(length_type const & new_pos) {
	return Scan2D::change_length_of_bond(*bond, new_pos);
}

coords::Cartesian_Point Scan2D::Normal_Angle_Input::make_move(length_type const & new_pos) {
	return Scan2D::rotate_a_to_new_angle(*angle, angle_type::from_deg(new_pos));
}

coords::Cartesian_Point Scan2D::Normal_Dihedral_Input::make_move(length_type const & new_pos) {
	return Scan2D::rotate_a_to_new_dihedral(*dihedral, angle_type::from_deg(new_pos));
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