#include "2DScan.h"

Scan2D::Scan2D(coords::Coordinates & coords) : _coords(coords) {
	auto & both_whats = Config::get().scan2d.AXES;

	if (both_whats.size() > 2) {
		throw std::runtime_error("You can't pass more than two axis!");
	}

	parse_input(both_whats);



	/*std::cout << x_axis.what_kind << " with atoms: ";
	for (auto && el : x_axis.atoms) {
		std::cout << el << " ";
	}

	std::cout << "\n" << x_axis.to_position << " " << x_axis.scans << std::endl;

	std::cout << "\n\n";

	std::cout << y_axis.what_kind << "with atoms: ";
	for (auto && el : y_axis.atoms) {
		std::cout << el << " ";
	}

	std::cout << "\n" << y_axis.to_position << " " << y_axis.scans << std::endl;

	std::cout << _coords.xyz() << std::endl;*/

}

void Scan2D::Normal_Input::fill_what(std::vector<std::string> & splitted_vals, Scan2D::what & to_fill) {

	to_fill.what_kind = splitted_vals.front();

	std::transform(splitted_vals.begin() + 1, splitted_vals.end() - 2, std::back_inserter(to_fill.atoms),
		[](std::string const & atom) { return std::stoi(atom); }
	);

	auto double_back = splitted_vals.size() - 2;

	to_fill.to_position = std::stod(splitted_vals[double_back]);
	to_fill.scans = std::stoi(splitted_vals.back());

}

void Scan2D::parse_input(std::vector<std::string> const & input) {
	
	std::istringstream iss_x_axis(input.front());
	std::istringstream iss_y_axis(input.back());
	
	std::vector<std::string> splitted_x_vals{
		std::istream_iterator<std::string>{iss_x_axis},
		std::istream_iterator<std::string>{}
	};

	std::vector<std::string> splitted_y_vals{
		std::istream_iterator<std::string>{iss_y_axis},
		std::istream_iterator<std::string>{}
	};
	try {
		std::unique_ptr<Input_types> x_filler(Input_Factory(splitted_x_vals.size(), splitted_x_vals.front()));
		std::unique_ptr<Input_types> y_filler(Input_Factory(splitted_y_vals.size(), splitted_y_vals.front()));
		x_filler->fill_what(splitted_x_vals, x_axis);
		y_filler->fill_what(splitted_y_vals, y_axis);
	}
	catch (std::runtime_error e) {
		std::cout << e.what() << std::endl;
		exit(1);
	}

}

Scan2D::Input_types * Scan2D::Input_Factory(std::size_t size, std::string kind) {

	if(kind.compare("bond") == 0) {
		if (size == 5) {
			return new Normal_Input;
		}
	}
	else if (kind.compare("angle") == 0) {
		if (size == 6) {
			return new Normal_Input;
		}
	}
	else if (kind.compare("dihedral") == 0) {
		if (size == 7) {
			return new Normal_Input;
		}
	}
	throw std::runtime_error("Your Input is wrong. Excpected 'bond', 'angle' or 'dihedral' by specifieing your 2D Scan!");

}