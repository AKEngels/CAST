#ifndef H_2DSCAN
#define H_2DSCAN

#include<vector>
#include<sstream>
#include<iterator>
#include<memory>

#include"coords.h"
#include"coords_rep.h"
#include"coords_io.h"
#include"configuration.h"
#include"scon_vect.h"
#include"scon_angle.h"
#include"scon_spherical.h"

class Scan2D {
public:

	class dihedral {
	public:
		coords::Cartesian_Point const & a;
		coords::Cartesian_Point const & b;
		coords::Cartesian_Point const & c;
		coords::Cartesian_Point const & d;
		dihedral(coords::Cartesian_Point const & a, coords::Cartesian_Point const & b, coords::Cartesian_Point const & c, coords::Cartesian_Point const & d)
			: a(a), b(b), c(c), d(d)
		{}
	};

	class angle {
	public:
		coords::Cartesian_Point const & a;
		coords::Cartesian_Point const & b;
		coords::Cartesian_Point const & c;
		angle(coords::Cartesian_Point const & a, coords::Cartesian_Point const & b, coords::Cartesian_Point const & c)
			: a(a), b(b), c(c)
		{}
	};

	class bond {
	public:
		coords::Cartesian_Point const & a;
		coords::Cartesian_Point const & b;
		bond(coords::Cartesian_Point const & a, coords::Cartesian_Point const & b)
			: a(a), b(b)
		{}
	};

	Scan2D(coords::Coordinates & coords);

	static coords::float_type get_length(bond const & ab);
	static coords::float_type get_angle(angle const & abc);
	static coords::float_type get_dihedral(dihedral const & abcd);

	static coords::Cartesian_Point change_length_of_bond(Scan2D::bond const & ab, coords::float_type const & new_length);
	static coords::Cartesian_Point rotate_a_to_new_angle(angle const & abc, coords::float_type const & new_angle);
	static coords::Cartesian_Point rotate_a_to_new_dihedral(dihedral const & abcd, coords::float_type const & new_dihedral);


private:


	struct what {
		std::string what_kind;
		std::vector<std::size_t> atoms;
		double to_position;
		double from_position;
		std::size_t scans;
	};

	class Input_types {
	public:
		virtual void fill_what(std::vector<std::string> & splitted_vals, coords::Representation_3D const & xyz) = 0;
		virtual void set_coords(coords::Representation_3D const & xyz) = 0;
		virtual coords::float_type say_val() = 0;
		virtual std::vector<coords::float_type> make_axis() = 0;
		virtual coords::Cartesian_Point make_move(coords::float_type const & new_pos) = 0;
	public:
		std::unique_ptr<Scan2D::what> what;
	};

	class Normal_Input : public Input_types {
		virtual void fill_what(std::vector<std::string> & splitted_vals, coords::Representation_3D const & xyz) override;
	};

	class Normal_Bond_Input : public Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual coords::float_type say_val() override;
		virtual std::vector<coords::float_type> make_axis() override;
		virtual coords::Cartesian_Point make_move(coords::float_type const & new_pos) override;
	public:
		std::unique_ptr<Scan2D::bond> bond;
	};

	class Normal_Angle_Input : public Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual coords::float_type say_val() override;
		virtual std::vector<coords::float_type> make_axis() override;
		virtual coords::Cartesian_Point make_move(coords::float_type const & new_pos) override;
	public:
		std::unique_ptr<Scan2D::angle> angle;
	};

	class Normal_Dihedral_Input : public Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual coords::float_type say_val() override;
		virtual std::vector<coords::float_type> make_axis() override;
		virtual coords::Cartesian_Point make_move(coords::float_type const & new_pos) override;
	public:
		std::unique_ptr<Scan2D::dihedral> dihedral;
	};

	class XY_Parser {
	public:
		std::unique_ptr<Scan2D::Input_types> x_parser;
		std::unique_ptr<Scan2D::Input_types> y_parser;
		XY_Parser(std::unique_ptr<Scan2D::Input_types> x, std::unique_ptr<Scan2D::Input_types> y) 
			: x_parser(std::move(x)), y_parser(std::move(y)){}
		void fix_atoms(coords::Coordinates & coords)const;
	};

	class XY_steps {
	public:
		std::vector<coords::float_type> x_steps;
		std::vector<coords::float_type> y_steps;
		XY_steps(std::vector<coords::float_type> & x, std::vector<coords::float_type> & y)
			: x_steps(std::move(x)), y_steps(std::move(y)){}
	};

	std::unique_ptr<Scan2D::Input_types> parse_input(std::string const & input);
	std::unique_ptr<Scan2D::Input_types> Input_Factory(std::size_t size, std::string kind);

	void make_scan(XY_Parser const & parser, XY_steps const & steps);
	void prepare_scan(XY_Parser const & parser);
	void make_x_change(coords::float_type const & change, Scan2D::XY_Parser const & parser, coords::Representation_3D & y_steps);
	void go_along_y_axis(Scan2D::XY_Parser const & parser, std::vector<coords::float_type> const & y_steps, coords::Coordinates coords);


	coords::Coordinates & _coords;

	std::ofstream logfile;
	std::ofstream energies;

	int x_circle = 0;
	int y_circle = 0;
};

#endif