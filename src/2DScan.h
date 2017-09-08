#ifndef H_2DSCAN
#define H_2DSCAN

#include<vector>
#include<sstream>
#include<iterator>
#include<memory>
#include<unordered_set>
//#include<boost/coroutine2/coroutine.hpp>

#include"coords.h"
#include"coords_rep.h"
#include"coords_io.h"
#include"configuration.h"
#include"scon_vect.h"
#include"scon_angle.h"
#include"scon_spherical.h"
#include"scon_mathmatrix.h"

class Scan2D : public std::enable_shared_from_this<Scan2D> {
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

    struct bond_hash;
    struct compare_bond_pair;
	void execute_scan();
	using length_type = typename get_spherical_types<remove_cr<decltype(_coords.intern(0))>>::length_type;
	using angle_type = typename get_spherical_types<remove_cr<decltype(_coords.intern(0))>>::angle_type;
    using bond_set = std::unordered_set<std::pair<std::size_t, std::size_t>, bond_hash, compare_bond_pair>;

    class dihedral {
    public:
      coords::Cartesian_Point & a;
      coords::Cartesian_Point & b;
      coords::Cartesian_Point & c;
      coords::Cartesian_Point & d;
      dihedral(coords::Cartesian_Point & a, coords::Cartesian_Point & b, coords::Cartesian_Point & c, coords::Cartesian_Point & d)
        : a(a), b(b), c(c), d(d)
      {}
    };

    class angle {
    public:
      coords::Cartesian_Point & a;
      coords::Cartesian_Point & b;
      coords::Cartesian_Point & c;
      angle(coords::Cartesian_Point & a, coords::Cartesian_Point & b, coords::Cartesian_Point & c)
        : a(a), b(b), c(c)
      {}
    };

    class bond {
    public:
      coords::Cartesian_Point & a;
      coords::Cartesian_Point & b;
      bond(coords::Cartesian_Point & a, coords::Cartesian_Point & b)
        : a(a), b(b)
      {}
    };

	class cdihedral {
	public:
		coords::Cartesian_Point const & a;
		coords::Cartesian_Point const & b;
		coords::Cartesian_Point const & c;
		coords::Cartesian_Point const & d;
        cdihedral(dihedral const & d)
            : a(d.a), b(d.b), c(d.c), d(d.d)
        {}
		cdihedral(coords::Cartesian_Point const & a, coords::Cartesian_Point const & b, coords::Cartesian_Point const & c, coords::Cartesian_Point const & d)
			: a(a), b(b), c(c), d(d)
		{}
	};

	class cangle {
	public:
		coords::Cartesian_Point const & a;
		coords::Cartesian_Point const & b;
		coords::Cartesian_Point const & c;
        cangle(angle const & a)
          : a(a.a), b(a.b), c(a.c)
        {}
		cangle(coords::Cartesian_Point const & a, coords::Cartesian_Point const & b, coords::Cartesian_Point const & c)
			: a(a), b(b), c(c)
		{}
	};

	class cbond {
	public:
		coords::Cartesian_Point const & a;
		coords::Cartesian_Point const & b;
        cbond(bond const & b)
          : a(b.a), b(b.a)
        {}
		cbond(coords::Cartesian_Point const & a, coords::Cartesian_Point const & b)
			: a(a), b(b)
		{}
	};

	Scan2D(coords::Coordinates & coords);
    Scan2D() = delete;
	~Scan2D();

	static length_type get_length(cbond const & ab);
	static angle_type get_angle(cangle const & abc);
	static angle_type get_dihedral(cdihedral const & abcd);

	static coords::Cartesian_Point change_length_of_bond(cbond const & ab, length_type const & new_length);
	static coords::Cartesian_Point rotate_a_to_new_angle(cangle const & abc, angle_type const & new_angle);
	static coords::Cartesian_Point rotate_a_to_new_dihedral(cdihedral const & abcd, angle_type const & new_dihedral);

    struct bond_hash {
      std::size_t operator()(std::pair<std::size_t, std::size_t> const & p)const {
        return std::hash<std::size_t>()(p.first);
      }
    };
    struct compare_bond_pair {
      bool operator()(std::pair<std::size_t, std::size_t> const & lhs, std::pair<std::size_t, std::size_t> const & rhs)const {
        return lhs.first == rhs.first;
      }
    };

private:
	struct what {
		std::string what_kind;
		std::vector<std::size_t> atoms;
		double to_position;
		double from_position;
		std::size_t scans;
        bool prepare_position = false;
	};

	class Input_types {
	public:
		virtual void fill_what(std::vector<std::string> & splitted_vals, coords::Representation_3D const & xyz) = 0;
		virtual void set_coords(coords::Representation_3D const & xyz) = 0;
		virtual length_type say_val() = 0;
		virtual std::vector<length_type> make_axis() = 0;
		virtual coords::Representation_3D make_move(length_type const & new_pos, std::vector<std::size_t> const & atoms) = 0;
	public:
        Input_types(std::weak_ptr<Scan2D> & p) : parent(p) {}
		std::unique_ptr<Scan2D::what> what;
        std::weak_ptr<Scan2D> parent;
	};

	class Normal_Input : public Input_types {
    public:
        Normal_Input(std::weak_ptr<Scan2D> & p) : Input_types(p) {}
		virtual void fill_what(std::vector<std::string> & splitted_vals, coords::Representation_3D const & xyz) override;
	};

	class Normal_Bond_Input : public Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual length_type say_val() override;
		virtual std::vector<length_type> make_axis() override;
		virtual coords::Representation_3D make_move(length_type const & new_pos, std::vector<std::size_t> const & atoms) override;
	public:
        Normal_Bond_Input(std::weak_ptr<Scan2D> p) : Normal_Input(p) {}
        std::unique_ptr<Scan2D::cbond> bond;
	};

	class Normal_Angle_Input : public Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual length_type say_val() override;
		virtual std::vector<length_type> make_axis() override;
		virtual coords::Representation_3D make_move(length_type const & new_pos, std::vector<std::size_t> const & atoms) override;
	public:
        Normal_Angle_Input(std::weak_ptr<Scan2D> p) : Normal_Input(p) {}
		std::unique_ptr<Scan2D::cangle> angle;
	};

	class Normal_Dihedral_Input : public Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual length_type say_val() override;
		virtual std::vector<length_type> make_axis() override;
		virtual coords::Representation_3D make_move(length_type const & new_pos, std::vector<std::size_t> const & atoms) override;
	public:
      Normal_Dihedral_Input(std::weak_ptr<Scan2D> p) : Normal_Input(p) {}
		std::unique_ptr<Scan2D::cdihedral> dihedral;
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

    double const change_from_atom_to_atom;
    double const max_change_rotation;

	std::ofstream logfile;
	std::ofstream energies;

	int x_circle = 0;
	int y_circle = 0;

	std::string const energie_file = coords::output::filename("_ENERGIES", ".txt");
	std::string const structures_file = coords::output::filename("_STRUCTURES", ".arc");

    coords::Representation_3D rotate_molecule_behind_a_dih(std::vector<std::size_t> const & abcd, Scan2D::length_type const & deg);
    coords::Representation_3D rotate_molecule_behind_a_ang(std::vector<std::size_t> const & abc, Scan2D::length_type const & deg);
    coords::Representation_3D transform_molecule_behind_a_bond(std::vector<std::size_t> const & ab, Scan2D::length_type const & length);


    //using coroutine_type = boost::coroutines2::coroutine<std::pair<std::size_t, std::size_t>>;

    bond_set go_along_backbone(std::vector<std::size_t> const & kind);
    //void go_along_backbone(coroutine_type::push_type & sink ,std::size_t const & atom, std::size_t const & border);

};

#endif
