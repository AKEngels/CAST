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
private:
	coords::Coordinates & _coords;

public:

	template<typename T>
	using remove_cr = typename std::remove_const<typename std::remove_reference<T>::type>::type;

	template<typename T>
	struct get_spherical_types {
	};

	template<typename T, typename U>
	struct get_spherical_types<scon::sphericals<T, U>> {
		using length_type = T;
		using angle_type = U;
	};

	using length_type = typename get_spherical_types<remove_cr<decltype(_coords.intern(0))>>::length_type;
	using angle_type = typename get_spherical_types<remove_cr<decltype(_coords.intern(0))>>::angle_type;

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
	~Scan2D();

	static length_type get_length(bond const & ab);
	static angle_type get_angle(angle const & abc);
	static angle_type get_dihedral(dihedral const & abcd);

	static coords::Cartesian_Point change_length_of_bond(Scan2D::bond const & ab, length_type const & new_length);
	static coords::Cartesian_Point rotate_a_to_new_angle(angle const & abc, angle_type const & new_angle);
	static coords::Cartesian_Point rotate_a_to_new_dihedral(dihedral const & abcd, angle_type const & new_dihedral);


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
		virtual length_type say_val() = 0;
		virtual std::vector<length_type> make_axis() = 0;
		virtual coords::Cartesian_Point make_move(length_type const & new_pos) = 0;
	public:
		std::unique_ptr<Scan2D::what> what;
	};

	class Normal_Input : public Input_types {
		virtual void fill_what(std::vector<std::string> & splitted_vals, coords::Representation_3D const & xyz) override;
	};

	class Normal_Bond_Input : public Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual length_type say_val() override;
		virtual std::vector<length_type> make_axis() override;
		virtual coords::Cartesian_Point make_move(length_type const & new_pos) override;
	public:
		std::unique_ptr<Scan2D::bond> bond;
	};

	class Normal_Angle_Input : public Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual length_type say_val() override;
		virtual std::vector<length_type> make_axis() override;
		virtual coords::Cartesian_Point make_move(length_type const & new_pos) override;
	public:
		std::unique_ptr<Scan2D::angle> angle;
	};

	class Normal_Dihedral_Input : public Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual length_type say_val() override;
		virtual std::vector<length_type> make_axis() override;
		virtual coords::Cartesian_Point make_move(length_type const & new_pos) override;
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
		std::vector<length_type> x_steps;
		std::vector<length_type> y_steps;
		XY_steps(std::vector<length_type> & x, std::vector<length_type> & y)
			: x_steps(std::move(x)), y_steps(std::move(y)){}
	};

	std::unique_ptr<Scan2D::Input_types> parse_input(std::string const & input);
	std::unique_ptr<Scan2D::Input_types> Input_Factory(std::size_t size, std::string kind);

	void make_scan();
	void prepare_scan();
	//void make_x_change(length_type const & change, Scan2D::XY_Parser const & parser, coords::Representation_3D & y_steps);
	void go_along_y_axis(coords::Coordinates coords);

	std::unique_ptr<XY_Parser> parser;
	std::unique_ptr<XY_steps> axis;

	inline void write_energy_entry(double const & e) {

		energies << std::fixed << std::setprecision(5) <<
			std::setw(12) << parser->x_parser->say_val() <<
			std::setw(12) << parser->y_parser->say_val() <<
			std::setw(12) << e << "\n";

	}

	std::ofstream logfile;
	std::ofstream energies;

	int x_circle = 0;
	int y_circle = 0;

	std::string const energie_file = coords::output::filename("_ENERGIES", ".txt");
	std::string const structures_file = coords::output::filename("_STRUCTURES", ".arc");
};

#endif