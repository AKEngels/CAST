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

/**
  * @brief A class handling a two dimensional PES-scan.
  * @author Julian Erdmannsdörfer
  * @version 3.0
  *
  * The class inherits from enable_shared_from_this and therefore needs to be instantiated as a shared pointer.
  *
  * The types used in Scan2D are deduced from the types used in the coords class and thus do not need to be changed if the coords types are changed.
  *
  * The function execute_scan is the only one, besides the ctor, to be called from outside. The used variables got to be provided by the CAST input file.
*/

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
    /**
    * @brief takes preparation steps and executes the scan
    * @see make_scan()
    */
	void execute_scan();
	using length_type = typename get_spherical_types<remove_cr<decltype(_coords.intern(0))>>::length_type;
	using angle_type = typename get_spherical_types<remove_cr<decltype(_coords.intern(0))>>::angle_type;
    using bond_set = std::unordered_set<std::pair<std::size_t, std::size_t>, bond_hash, compare_bond_pair>;/**< alias for a unordered set, a hash-function, and a compare-function*/

    /**
    * @brief simple wrapper for four cartesian points as non const references
    * @see cdihedral
    */
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

    /**
    * @brief simple wrapper for three cartesian points as non const references
    * @see cangle
    */
    class angle {
    public:
      coords::Cartesian_Point & a;
      coords::Cartesian_Point & b;
      coords::Cartesian_Point & c;
      angle(coords::Cartesian_Point & a, coords::Cartesian_Point & b, coords::Cartesian_Point & c)
        : a(a), b(b), c(c)
      {}
    };

    /**
    * @brief simple wrapper for two cartesian points as non const references
    * @see cbond
    */
    class bond {
    public:
      coords::Cartesian_Point & a;
      coords::Cartesian_Point & b;
      bond(coords::Cartesian_Point & a, coords::Cartesian_Point & b)
        : a(a), b(b)
      {}
    };

    /**
    * @brief simple wrapper for four cartesian points as const references
    * @see dihedral
    */
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

    /**
    * @brief simple wrapper for three cartesian points as const references
    * @see angle
    */
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

    /**
    * @brief simple wrapper for two cartesian points as const references
    * @see bond
    */
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

    /**
    * @brief The only availible ctor to initialize some variables, set flags, and open files.
    * @param coords used mainly for the Cartesian representation
    */
	Scan2D(coords::Coordinates & coords);
    Scan2D() = delete;
    /**
    * @brief An empty dtor.
    */
	~Scan2D();

	static length_type get_length(cbond const & ab);
	static angle_type get_angle(cangle const & abc);
	static angle_type get_dihedral(cdihedral const & abcd);

	static coords::Cartesian_Point change_length_of_bond(cbond const & ab, length_type const & new_length);
	static coords::Cartesian_Point rotate_a_to_new_angle(cangle const & abc, angle_type const & new_angle);
	static coords::Cartesian_Point rotate_a_to_new_dihedral(cdihedral const & abcd, angle_type const & new_dihedral);
    /**
    * @brief The hash functor for the bond_set
    */
    struct bond_hash {
      /**
      * @param p the pait which is to hash
      * The only relevant value is the first one which is the atom index which defines if a atom is present or not
      * hashing only this index and keeping the second value unchanged should do the trick.
      */
      std::size_t operator()(std::pair<std::size_t, std::size_t> const & p)const {
        return std::hash<std::size_t>()(p.first);
      }
    };
    /**
    * @brief The compare functor for the bond_set
    */
    struct compare_bond_pair {
      /**
      * @param lhs The left hand side for the comparison
      * @param rhs The right hand side for the comparison
      * In both pairs the only relevant value is the one holding the index for the specific atom und thus only the first value is going to be compared to
      * identify atoms.
      */
      bool operator()(std::pair<std::size_t, std::size_t> const & lhs, std::pair<std::size_t, std::size_t> const & rhs)const {
        return lhs.first == rhs.first;
      }
    };

private:
    /**
    * @brief a simple wrapper struct to hold all information needed for an axis
    */
	struct what {
		std::string what_kind; /**<Variable to store which kind of constraint is used for this axis*/
		std::vector<std::size_t> atoms; /**<Variable to store the patricipating atom's indices*/
		double to_position; /**<Variable to strore the end position of the scan*/
		double from_position; /**<Variable to store the starting position of the scan*/
		std::size_t scans; /**<Variable to store the steps the scan takes from starting to end position including start and end structure*/
	};

    /**
    * @brief This struct stores all information to change the constraints value and thus move the atoms accordingly
    * The class' purpose is to handle the change in constraints and yields the possibility to rotate a greater part of the molecule behind the specific constraint in order to
    * not break the structure of a molecule.
    */
    struct Move_Handler {
      /**
      * @brief The ctor initializing all variables needed for a proper move except the new destination which is to defined withe the set_new_pos function
      * @see set_new_pos(length_type const & c)
      */
      Move_Handler(coords::Coordinates & coords, std::vector<std::size_t> const & atoms, std::shared_ptr<Scan2D> p) : _coords(coords), atoms(atoms), parent(std::move(p)) {}
      coords::Coordinates & _coords; /**<A reference to the used coords object*/
      std::shared_ptr<Scan2D> parent; /**<A shared pointer to the parent Scan2D class*/
      std::vector<std::size_t> const & atoms; /**<A constant reference to the std::vector holding the participating atom's indices*/
      length_type new_pos = 0.0; /**<Variable to hold the new position to which is to move*/

      /**
      * @brief Simple setter for the new_pos variable
      * @param c The value of the new position to which constraint is changed by moving the participating atoms
      */
      void set_new_pos(length_type const & c) { new_pos = c; }

      /**
      * @brief Rotate a dihedral constraint and a part of the molecule, too
      */
      coords::Representation_3D rotate_molecule_behind_a_dih(Scan2D::length_type const & deg) const;
      /**
      * @brief Rotate a angle constraint and a part of the molecule, too
      */
      coords::Representation_3D rotate_molecule_behind_a_ang(Scan2D::length_type const & deg) const;
      /**
      * @brief Rotate a bond constraint and a part of the molecule, too
      */
      coords::Representation_3D transform_molecule_behind_a_bond(Scan2D::length_type const & length) const;

      //using coroutine_type = boost::coroutines2::coroutine<std::pair<std::size_t, std::size_t>>;

      /**
      * @brief Scans the atoms behind the moved atom and returns the indices and the number of bonds between the atom which is to move and the respective one
      */
      bond_set go_along_backbone(std::vector<std::size_t> const & kind) const;
      //void go_along_backbone(coroutine_type::push_type & sink ,std::size_t const & atom, std::size_t const & border);
    };

    /**
    * @brief struct as a base class for all the possible constraints as a pure virtual class 
    */
	struct Input_types {
		virtual void fill_what(std::vector<std::string> & splitted_vals, coords::Representation_3D const & xyz) = 0;/**<Function to parse the needed information for an axis */
		virtual void set_coords(coords::Representation_3D const & xyz) = 0;/**<Function to set the xyz data for the, in the constraint, participating atoms*/
		virtual length_type say_val() = 0;/**<returns the current value of the angle or distace*/
		virtual std::vector<length_type> make_axis() = 0;/**<using the angle or distance stored in what.from_position and what.to_position*/
		virtual coords::Representation_3D make_move(Move_Handler const & mh) = 0;/**<uses the provided move handler to perform the corresponding move function for the respective constraint*/
        Input_types(std::weak_ptr<Scan2D> & p) : parent(p) {}
		std::unique_ptr<Scan2D::what> what;/**<Unique pointer to store the information for one axis */
        std::weak_ptr<Scan2D> parent;/**<A weak pointer to the instantiating parent */
	};

    /**
    * @brief The second level of inheritance for "normal" constraints
    * The normal constraints is the group of constraints related to internal coordinates e. i. bonds, angles, or diehdrals.
    * All three types have a pretty similar parsing in common and thus share the same method fill_what.
    */
	struct Normal_Input : Input_types {
        /**
        * @brief Ctor initializing the shared pointer to the inherited parent variable.
        * @param p weak pointer to be locked to be stored in the parent variable.
        */
        Normal_Input(std::weak_ptr<Scan2D> & p) : Input_types(p) {}
        /**
        * @brief Function to read the, in the input file, provided line to specify the one axis property.
        */
		virtual void fill_what(std::vector<std::string> & splitted_vals, coords::Representation_3D const & xyz) override;
	};

    /**
    * @brief Structure to store the information for a bond constraint
    */
	struct Normal_Bond_Input : Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual length_type say_val() override;
		virtual std::vector<length_type> make_axis() override;
		virtual coords::Representation_3D make_move(Move_Handler const & mh) override;
        Normal_Bond_Input(std::weak_ptr<Scan2D> p) : Normal_Input(p) {}
        std::unique_ptr<Scan2D::cbond> bond;
	};

    /**
    * @brief Structure to store the information for a angle constraint
    */
	struct Normal_Angle_Input : Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual length_type say_val() override;
		virtual std::vector<length_type> make_axis() override;
		virtual coords::Representation_3D make_move(Move_Handler const & mh) override;
        Normal_Angle_Input(std::weak_ptr<Scan2D> p) : Normal_Input(p) {}
		std::unique_ptr<Scan2D::cangle> angle;
	};

    /**
    * @brief Structure to store the information for a dihedral constraint
    */
	struct Normal_Dihedral_Input : Normal_Input {
		virtual void set_coords(coords::Representation_3D const & xyz) override;
		virtual length_type say_val() override;
		virtual std::vector<length_type> make_axis() override;
		virtual coords::Representation_3D make_move(Move_Handler const & mh) override;
       Normal_Dihedral_Input(std::weak_ptr<Scan2D> p) : Normal_Input(p) {}
		std::unique_ptr<Scan2D::cdihedral> dihedral;
	};

    /**
    * @brief Wrapper structure for both axes
    */
	class XY_Parser {
	public:
		std::unique_ptr<Scan2D::Input_types> x_parser;
		std::unique_ptr<Scan2D::Input_types> y_parser;
		XY_Parser(std::unique_ptr<Scan2D::Input_types> x, std::unique_ptr<Scan2D::Input_types> y) 
			: x_parser(std::move(x)), y_parser(std::move(y)){}
		void fix_atoms(coords::Coordinates & coords)const;
	};

    /**
    * @brief Wrapper structure for both axis-grids
    */
	class XY_Steps {
	public:
		std::vector<length_type> x_steps;
		std::vector<length_type> y_steps;
		XY_Steps(std::vector<length_type> & x, std::vector<length_type> & y)
			: x_steps(std::move(x)), y_steps(std::move(y)){}
	};

    /**
    * @brief function to parse the relevant input file line and build the corresponding axis.
    */
	std::unique_ptr<Scan2D::Input_types> parse_input(std::string const & input);

    /**
    * @brief factory method to return the respective axis for a specific input line.
    */
	std::unique_ptr<Scan2D::Input_types> Input_Factory(std::size_t size, std::string kind);

    /**
    * @brief Using the provided information one step after another is taken to scan all grid-points
    * For each x-step the structure is copied and the scan is done along the y-axis. The y-axis scan is done in go_along_y_axis
    * @see go_along_y_axis(coords::Coordinates coords)
    */
	void make_scan();

    /**
    * @brief This function takes all the neccessary preparation steps like moving the constraint in the right beginning structure.
    */
	void prepare_scan();
	//void make_x_change(length_type const & change, Scan2D::XY_Parser const & parser, coords::Representation_3D & y_steps);

    /**
    * @brief Performes the scan along the y-axis
    * @param coords Copy of the coords object
    */
	void go_along_y_axis(coords::Coordinates coords);

    /**
    * @brief The function calls the optimization routine of the coords object and sets the optimized Cartesian coordinates.
    */
    length_type optimize(coords::Coordinates & c);

	std::unique_ptr<XY_Parser> parser;/**<Instance of the wrapper tyoe XY_Parser as unique pointer to store the axes information*/
	std::unique_ptr<XY_Steps> axis;/**<Instance of the wrapper tyoe XY_Steps as unique pointer to store the axis-grids*/

    /**
    * @brief A function for writing the energy output in the outputfile.
    */
	inline void write_energy_entry(double const & e) {

      auto bla = parser->x_parser->say_val();

		energies << std::fixed << std::setprecision(5) <<
			std::setw(12) << parser->x_parser->say_val() << " " <<
			std::setw(12) << parser->y_parser->say_val() << " " <<
			std::setw(12) << e << "\n";

	}

    double const change_from_atom_to_atom;/**<The user provided value for the rotation in the preparation step*/
    double const max_change_rotation;/**<The threshhold which causes the rotation of the whole submolecule if it is overstepped*/

	std::ofstream logfile;/**<The file in which the evaluated structures are saved*/
	std::ofstream energies;/**<The file in which the evaluated energies are saved with the write_energy_entry function*/
	std::ofstream before;/**<A file for debug purpose saving the unoptimized structures*/

	std::string const energie_file = coords::output::filename("_ENERGIES", ".txt");/**<The name of the file for the energies*/
	std::string const structures_file = coords::output::filename("_STRUCTURES", ".arc");/**<The name of the file for the structures*/

};

#endif
