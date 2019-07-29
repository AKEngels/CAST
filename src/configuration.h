/**
CAST 3
configuration.h
Purpose: class for extraction of information from inputfile

@author Daniel Weber (modified by many)
@version 1.1
*/

#ifndef H_CONFIGURATION
#define H_CONFIGURATION


#include <cstddef>
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
#include<memory>
#include <map>
#include <utility>
#include <array>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cctype>

#if defined _OPENMP
#include <omp.h>
#endif

#include "Scon/scon.h"
#include "filemanipulation.h"
#include "Scon/scon_utility.h"
#include "Scon/scon_vect.h"
#include "coords_rep.h"
#include "configurationHelperfunctions.h"


/*! Namespace containing relevant configuration options
 */
namespace config
{
  /**get a vector of integers from a string, string is something like "5-7,19,42"*/
  std::vector<std::size_t> sorted_indices_from_cs_string(std::string str, bool minus_1 = false);

  /**function that reads a string that consists of numbers, seperated by comma, into a vector of doubles*/
  std::vector<double> doubles_from_string(std::string str);

  /**function that reads a string that consists of integers, seperated by comma, into a vector of integers*/
  std::vector<int> ints_from_string(std::string str);

  // Here we find some static members that only
  // exist once in CAST, like the version number or
  // some helper arrays containing the tasks etc.

  /** Name of the program*/
  static std::string const Programname("CAST");
  /** Version-Number of CAST*/
  static std::string const Version("3.2.0.2dev");


  /**Number of tasks*/
  static std::size_t const NUM_TASKS = 35;

  /** Names of all CAST tasks as strings*/
  static std::string const task_strings[NUM_TASKS] =
  {
    "SP", "GRAD", "TS", "LOCOPT", "REMOVE_EXPLICIT_WATER",
    "MC", "DIMER", "MD", "NEB", "GOSOL",
    "STARTOPT",  "INTERNAL", "ENTROPY", "PCAgen", "PCAproc",
    "DEVTEST", "UMBRELLA", "FEP", "PATHOPT",
    "GRID", "ALIGN", "PATHSAMPLING", "SCAN2D", "XB_EXCITON_BREAKUP",
    "XB_INTERFACE_CREATION", "XB_CENTER", "XB_COUPLINGS",
    "LAYER_DEPOSITION", "HESS", "WRITE_TINKER", "MODIFY_SK_FILES", "WRITE_GAUSSVIEW", 
    "MOVE_TO_ORIGIN", "WRITE_XYZ", "WRITE_PDB"
  };

  /*! contains enum with all tasks currently present in CAST
   *
   * Those taks are subsequently mapped using task_strings string[].
   */
  struct tasks
  {
    /*! contains all tasks currently present in CAST
    */
    enum T
    {
      ILLEGAL = -1,
      SP, GRAD, TS, LOCOPT, REMOVE_EXPLICIT_WATER,
      MC, DIMER, MD, NEB, GOSOL,
      STARTOPT, INTERNAL, ENTROPY, PCAgen, PCAproc,
      DEVTEST, UMBRELLA, FEP, PATHOPT,
      GRID, ALIGN, PATHSAMPLING, SCAN2D, XB_EXCITON_BREAKUP,
      XB_INTERFACE_CREATION, XB_CENTER, XB_COUPLINGS,
      LAYER_DEPOSITION, HESS, WRITE_TINKER, MODIFY_SK_FILES, WRITE_GAUSSVIEW,
      MOVE_TO_ORIGIN, WRITE_XYZ, WRITE_PDB
    };
  };

  /** number of Input Types */
  static std::size_t const NUM_INPUT = 4;
  /** Input Types */
  static std::string const input_strings[NUM_INPUT] =
  {
    "TINKER", "AMBER", "XYZ", "PDB"
  };

  /*! contains enum with all input_types currently supported in CAST
   *
   * Those input_types are subsequently mapped using input_strings string[].
   */
  struct input_types
  {
    /*! contains all input_types currently supported in CAST
    */
    enum T
    {
      ILLEGAL = -1,
      TINKER, AMBER, XYZ, PDB
    };
  };

  /**number of Output Types*/
  static std::size_t const NUM_OUTPUT = 4;
  /**Output Types*/
  static std::string const output_strings[NUM_OUTPUT] =
  {
    "TINKER", "XYZ", "MOLDEN", "ZMATRIX"
  };

  /*! contains enum with all output_types currently supported in CAST
   *
   * Those output_types are subsequently mapped using output_strings string[].
   */
  struct output_types
  {
    /*! contains all output_types currently supported in CAST
    */
    enum T
    {
      ILLEGAL = -1,
      TINKER, XYZ, MOLDEN, ZMATRIX
    };
  };

  /**number of Interface Types*/
  static std::size_t const NUM_INTERFACES = 15;

  /**Interface Types*/
  static std::string const
    interface_strings[NUM_INTERFACES] =
  {
    "AMBER", "AMOEBA", "CHARMM22", "OPLSAA", "TERACHEM", "MOPAC" , "DFTBABY", "GAUSSIAN", "QMMM", "DFTB", "CHEMSHELL", "PSI4", "ONIOM", "THREE_LAYER", "ORCA"
  };

  /*! contains enum with all energy interface_types currently supported in CAST
   *
   * Those interface_types are subsequently mapped using interface_strings string[].
   */
  struct interface_types
  {
    /*! contains all interface_types currently supported in CAST
    */
    enum T
    {
      ILLEGAL = -1,
      AMBER, AMOEBA, CHARMM22, OPLSAA, TERACHEM, MOPAC, DFTBABY, GAUSSIAN, QMMM, DFTB, CHEMSHELL, PSI4, ONIOM, THREE_LAYER, ORCA
    };
  };

  /**number of supported Mopac Versions*/
  static std::size_t const NUM_MOPAC_VERSION = 4;
  /**supported Mopac Versions*/
  static std::string const
    mopac_ver_string[NUM_MOPAC_VERSION] =
  {
    "2012", "2012MT", "7", "AVOID_HB"
  };

  /*! contains enum with all MOPAC versions currently supported as energy interfaces in CAST
   *
   * Those interface versions are subsequently mapped using mopac_ver_string string[].
   */
  struct mopac_ver_type
  {
    /*! contains all MOPAC versions currently supported in CAST
    */
    enum T
    {
      ILLEGAL = -1,
      MOPAC2012, MOPAC2012MT, MOPAC7, MOPAC7_HB
    };
  };

  /** number of Global optimization routines*/
  static std::size_t const NUM_GLOBOPT_ROUTINES = 2;
  /**Global optimization routines (TABUSEARCH, BASINHOPPING)*/
  static std::string const
    globopt_routines_str[NUM_GLOBOPT_ROUTINES] =
  {
    "TS", "BH"
  };

  /*! contains enum with all global optimization routines currently supported in CAST
   *
   * Those routines are subsequently mapped using NUM_GLOBOPT_ROUTINES string[].
   */
  struct globopt_routine_type
  {
    /*! contains all global optimization routines currently supported in CAST
    */
    enum T
    {
      ILLEGAL = -1,
      TABUSEARCH, BASINHOPPING
    };
  };

  // Global static stuff ends here...

  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////

  // ... now lets see about the members of the config namespace
  // They all have one instance as members in the global
  // config::Config object. This object contains all the
  // configoptions read from file for the current CAST run.

  /*! Struct containing all general information about the current CAST run
   */
  struct general
  {
    /** Name of the input file ("CAST.TXT")*/
    std::string inputFilename;
    /** Name of the force-field parameter file*/
    std::string paramFilename;
    /** Name of the output file*/
    std::string outputFilename;
    /**Type of the coordinate input (default: Tinker)*/
    input_types::T input;
    /** Type of the coordinate output (default: Tinker)*/
    output_types::T output;
    /** Current task*/
    config::tasks::T task;
    std::ofstream * trackstream;
    bool forcefield;

    /**Energy interface used for current run*/
    interface_types::T energy_interface;
    /**Energy interface used pre-optimization performed before the current run*/
    interface_types::T preopt_interface;
    /**Verbosity of the output of CAST (supposed to be between 0 and 5)*/
    std::size_t verbosity;
    /**are amber charges read from a seperate file?*/
    bool chargefile;

    /// Constructor with reasonable default parameters
    general(void) :
      paramFilename("oplsaa.prm"), outputFilename("%i.out"),
      input(input_types::TINKER), output(output_types::TINKER),
      task(config::tasks::SP), energy_interface(interface_types::OPLSAA),
      preopt_interface(interface_types::ILLEGAL),
      verbosity(1U), chargefile(false)
    { }
  };

  /**struct to collect all input information
  which doesn't fit anywhere else*/
  struct stuff
  {
    /**moving mode for task MOVE_TO_ORIGIN*/
    int moving_mode{ 0 };
    /**try to create energy type from amino acids if XYZ input is used*/
    bool xyz_atomtypes{ false };
  };

  struct periodics
  {
    // Periodic Box
    scon::c3<double> pb_box;
    // Are Periodic bounddries on?
    bool periodic;
    // Print periodic dummy atoms
    bool periodic_print;

    // Cut out atoms out of box when using periodics before calculation
    bool periodicCutout;
    // Tolerance for cut
    double cutout_distance_to_box;
    //
    unsigned int criterion;
    periodics(void) :
      pb_box(10.0, 10.0, 10.0), periodic(false), periodic_print(false),
      periodicCutout(false), cutout_distance_to_box(0.), criterion(0u)
    {
      if ((pb_box.x() <= cutout_distance_to_box
        || pb_box.y() <= cutout_distance_to_box
        || pb_box.z() <= cutout_distance_to_box) && periodicCutout)
      {
        throw std::runtime_error("Cutout distance cannot be bigger than box size for periodic boundries. Aborting.");
      }
    }
  };

  /*! Stream operator for config::periodics
   *
   * Prints configuration details for the current CAST run
   * Contains: Information about periodic box and periodic cutout functionality,
   */
  std::ostream & operator << (std::ostream &, periodics const &);

  /*
  ########  ####    ###     ######
  ##     ##  ##    ## ##   ##    ##
  ##     ##  ##   ##   ##  ##
  ########   ##  ##     ##  ######
  ##     ##  ##  #########       ##
  ##     ##  ##  ##     ## ##    ##
  ########  #### ##     ##  ######
  */

  /**namespace for biased potentials*/
  namespace biases
  {
    /**additional potential on distance of given atoms*/
    struct distance
    {
      /**force constant*/
      double force;
      /**ideal distance*/
      double ideal;
      /**current value (might change during run)*/
      double value;
      /**number of one atom*/
      std::size_t a;
      /**number of the other atom*/
      std::size_t b;
      /**constructor*/
      distance(void)
        : force(), ideal(), a(), b()
      { }
    };
    /**additional potential on angle between given atoms*/
    struct angle
    {
      /**force constant*/
      double force;
      /**ideal angle*/
      double ideal;
      /**current value (might change during run)*/
      double value;
      /**number of one atom*/
      std::size_t a;
      /**number of next atom (peak of angle)*/
      std::size_t b;
      /**number of the third atom*/
      std::size_t c;
      /**constructor*/
      angle(void)
        : force(), ideal(), a(), b(), c()
      { }
    };
    /**additional potential on a given dihedral*/
    struct dihedral
    {
      /**force constant*/
      double force;
      /**ideal dihedral angle*/
      ::coords::angle_type ideal;
      /**current value (might change during run)*/
      ::coords::angle_type value;
      /**atom 1*/
      std::size_t a;
      /**atom 2*/
      std::size_t b;
      /**atom 3*/
      std::size_t c;
      /**atom 4*/
      std::size_t d;
      /**constructor*/
      dihedral(void)
        : force(), ideal(), value(),
        a(), b(), c(), d()
      { }
    };
    /**sperical potential - prevents non-bonded systems from exploding*/
    struct spherical
    {
      /** distance to center where the additional potential starts*/
      double radius;
      /** force constant */
      double force;
      /** exponent of the potential function, 2 for harmonic potential, 4 is also possible */
      double exponent;
      /**constructor*/
      spherical()
        : radius(), force(), exponent()
      { }
    };
    /**cubic potential (similar to spherical but cubic)*/
    struct cubic
    {
      /**size of cubic box as cartesian point*/
      ::coords::Cartesian_Point dim;
      /**force constant*/
      double force;
      /**exponent of the potential function*/
      double exponent;
      /**constructor*/
      cubic()
        : dim(), force(), exponent()
      { }
    };
    /**threshold potential (very special, only used in task layerdeposiotion)*/
    struct thresholdstr
    {
      /**force constant*/
      double forceconstant;
      /**threshold distance*/
      double th_dist;
      /**constructor*/
      thresholdstr(void)
        : forceconstant(), th_dist()
      { }
    };
  }


  /*
     ######   #######   #######  ########  ########   ######
    ##    ## ##     ## ##     ## ##     ## ##     ## ##    ##
    ##       ##     ## ##     ## ##     ## ##     ## ##
    ##       ##     ## ##     ## ########  ##     ##  ######
    ##       ##     ## ##     ## ##   ##   ##     ##       ##
    ##    ## ##     ## ##     ## ##    ##  ##     ## ##    ##
     ######   #######   #######  ##     ## ########   ######
  */

  /**stuff for coords object that can be read in by inputfile CAST.txt*/
  struct coords
  {
    /**vector with amber charges in amber units, i.e. they must be divided by 18.2223 to get elementary charge
    (only filled if AMBER input is used or option chargefile is selected)*/
    std::vector<double> amber_charges;

    /**stuff for internal coordinates*/
    struct internals
    {
      /**???*/
      std::map<std::size_t, std::size_t> connect;
      /**dihedrals given here can't be main dihedrals*/
      std::vector<std::pair<std::size_t, std::size_t>> main_whitelist;
      /**dihedrals given here must be main dihedrals*/
      std::vector<std::pair<std::size_t, std::size_t>> main_blacklist;
    } internal;

    /**stuff for umbrella sampling*/
    struct umbrellas
    {
      /**use umbrella combination biases also for other tasks?*/
      bool use_comb{ false };

      /**struct for restrained torsional angle*/
      struct umbrella_tor
      {
        /**force constant*/
        double force;
        /**angle to which it is restrained*/
        double angle;
        /**array of atom indices*/
        std::size_t index[4U];
        /**???*/
        bool fix_all_torsions;
        /**constructor*/
        umbrella_tor(void) :
          force(0.0), index(), fix_all_torsions(false) { }
      };

      /**struct for restrained angle*/
      struct umbrella_angle
      {
        /**force constant*/
        double force;
        /**angle to which it is restrained*/
        double angle;
        /**array of atom indizes*/
        std::size_t index[3U];
        /**constructor*/
        umbrella_angle(void) :
          force(0.0), index() {}
      };

      /**struct for restrained distance*/
      struct umbrella_dist
      {
        /**force constant*/
        double force;
        /**distance to which it is restrained*/
        double dist;
        /**array of atom indices*/
        std::size_t index[2U];
        /**constructor*/
        umbrella_dist(void) :
          force(0.0), index() { }
      };

      /**struct for a restrained reaction coordinate that consists of several distances*/
      struct umbrella_comb
      {
        /**struct for one of these distances*/
        struct uscoord {
          /**atom index of first atom*/
          int index1;
          /**atom index of second atom*/
          int index2;
          /**factor for weighing potential on this distance, 
          if you want to define umbrella combination as difference of 2 distances
          the factor for one of them is -1*/
          int factor;
        };
        /**force constant*/
        double force_final;
        /**current force constant (raises during first half of equilibration)*/
        double force_current{ 0.0 };
        /**value (in Angstrom) to which it is restrained*/
        double value;
        /**vector of all dists that are included in reaction coordinate*/
        std::vector<uscoord> dists;
      };
    } umbrella;

    /**biased potentials*/
    struct coord_bias
    {
      /**biased potentials on distances*/
      std::vector<biases::distance>  distance;
      /**biased potentials on angles*/
      std::vector<biases::angle>     angle;
      /**biased potentials on dihedrals*/
      std::vector<biases::dihedral>  dihedral;
      /**spherical potential*/
      std::vector<biases::spherical> spherical;
      /**cubic potentials*/
      std::vector<biases::cubic>     cubic;
      /**threshold potentials*/
      std::vector<biases::thresholdstr>     threshold;
      /**biased pot on torsions for umbrella sampling*/
      std::vector<config::coords::umbrellas::umbrella_tor> utors;
      /**biased pot on angles for umbrella sampling*/
      std::vector<config::coords::umbrellas::umbrella_angle> uangles;
      /**biased pot on bonds for umbrella sampling*/
      std::vector<config::coords::umbrellas::umbrella_dist> udist;
      /**biased pot on combinations of bonds for umbrella sampling*/
      std::vector<config::coords::umbrellas::umbrella_comb> ucombs;
    } bias;


    struct conditionsForStructuresToBeConsideredEqual
    {
      double superposition; // every atom is 'superposed' by an atom with the same atomic number within this radius in angstroms
      ::coords::main_type main; // none of the main torsions differ more then this
      ::coords::internal_type intern; // no internal vector (bond, angle, dihedral) differs more than this
      ::coords::Cartesian_Point xyz; // no xyz position differs more than this
      conditionsForStructuresToBeConsideredEqual() :
        superposition(0.4), main(::coords::angle_type::from_deg(8.0)),
        intern(0.2, ::coords::angle_type::from_deg(1.0), ::coords::angle_type::from_deg(8.0)),
        xyz(0.1, 0.1, 0.1)
      {}
    } equals;

    /**vector with numbers of fixed atoms, indizes starting with 0 (i.e. these atoms are not allowed to move)*/
    std::vector<std::size_t> fixed;

    /**struct that contains radius and index of central atom to fix atoms around a sphere*/
		struct fix_sphere {
      /**true if a fix-sphere is active*/
			bool use{ false };
      /**radius in Angstrom*/
			double radius;
      /**index of atom that defines sphere center (starting with 0)*/
			int central_atom;
		} fix_sphere;

    /**vector with subsystems (used e.g. for IN and OUT in FEP)*/
    std::vector<std::vector<std::size_t>> subsystems;
    /**are rotations where only hydrogens move counting for main dihedrals?*/
    bool remove_hydrogen_rot;
    /**are internals starting new with every molecule?*/
    bool decouple_internals;

    /**constructor*/
    coords(void) :
      internal(), umbrella(), bias(), equals(), fixed(), subsystems(),
      remove_hydrogen_rot(true),
      decouple_internals(false)
    {}

  };


  /*
    ######## ##    ## ######## ########   ######   ##    ##
    ##       ###   ## ##       ##     ## ##    ##   ##  ##
    ##       ####  ## ##       ##     ## ##          ####
    ######   ## ## ## ######   ########  ##   ####    ##
    ##       ##  #### ##       ##   ##   ##    ##     ##
    ##       ##   ### ##       ##    ##  ##    ##     ##
    ######## ##    ## ######## ##     ##  ######      ##
  */

  /*!
   *
   * @todo: role of isotropic bool unclear... is berendesn barostat working?
   */
  struct energy
  {
    /**cutoff for non-bonded interactions in forcefield interfaces*/
    double cutoff;
    /**radius to start switching function to kick in; scales interactions smoothly to zero at cutoff radius*/
    double switchdist;

    /**???*/
    bool isotropic;
    /**???*/
    bool remove_fixed;

    /**struct for spackman correction*/
    struct spack
    {
      double cut;
      bool on, interp;
      spack(void) : cut(10.0), on(false), interp(true) { }
    } spackman;

    /**struct that contains information necessary for QM/MM calculation (also with ONIOM and THREE_LAYER)*/
    struct qmmm_conf
    {
      /**definition of QM systems
      every element of the vector is one QM system (in additive QMMM and THREE_LAYER only the first one is used)
      each QM system is defined by a vector of the atom indices which belong to the QM system*/
      std::vector<std::vector<size_t>> qm_systems;
      /**indices of SE atoms [only for three-layer]*/
      std::vector <size_t> seatoms;
      /**MM interface*/
      interface_types::T mminterface{ interface_types::T::OPLSAA };
      /**SE interface [only for three-layer]*/
      interface_types::T seinterface{ interface_types::T::DFTB };
      /**QM interface*/
      interface_types::T qminterface{ interface_types::T::MOPAC };
      /**is QM/MM interface active?*/
      bool use{ false };
	    /**should QM region be written into file?*/
	    bool qm_to_file{ false };
      /**vector of MM charges (external charges for inner calculation)*/
		  std::vector<PointCharge> mm_charges;
      /**energy types of link atoms (in the order of MM atom)
      every element of the vector corresponds to one QM system (in additive QMMM and THREE_LAYER only the first one is used)*/
      std::vector<std::vector<int>> linkatom_sets;
      /**cutoff for electrostatic interaction*/
			double cutoff{std::numeric_limits<double>::max()};
			/**central atom for cutoff (as atom index)
			one element for each QM system*/
			std::vector<std::size_t> centers;

			// stuff for three-layer:

			/**for atoms that are seperated from the inner region by a maximum of ... bonds the charges are set to zero for electronic embedding (1, 2 or 3)*/
			int zerocharge_bonds{ 1 };
			/**electronic embedding type for smallest system (0=EEx, 1=3-EE, 2=MM+SE) [only for three-layer]*/
			int emb_small{ 1 };
			/**central atom for cutoff in small system (as atom index)*/
			std::size_t small_center;
    } qmmm{};

    /**struct that contains information necessary for MOPAC calculation*/
    struct mopac_conf
    {
      /**command that is given to MOPAC*/
      std::string command;
      /**path to MOPAC*/
      std::string path;
      /**MOPAC version*/
      mopac_ver_type::T version;
      /**should MOPAC input be deleted after run?*/
      bool delete_input;
			/**charge of total system*/
			int charge;
      /**number of link atoms (needed for QM/MM)*/
      std::size_t link_atoms;
      /**constructor*/
			mopac_conf(void) : command("PM7 MOZYME"),
#if defined(MOPAC_EXEC_PATH)
				path(MOPAC_EXEC_PATH)
#elif defined(_MSC_VER)
				path("\"C:\\Program Files\\mopac\\MOPAC2012.exe\""),
#else
				path("/opt/mopac/MOPAC2012.exe"),
#endif
				version(mopac_ver_type::T::MOPAC2012MT),
				delete_input(true),
				charge(0), link_atoms(0)
      {}
    } mopac;

    /**struct that contains all information necessary for DFTBaby calculation*/
    struct dftbaby_conf
    {
      /**path to dftbaby*/
      std::string path;
      /**name of dftbaby gradient file (deleted again but necessary because otherwise
      dftbaby doesn't calculate gradients)*/
      std::string gradfile;
      /**total charge of the molecule*/
      int charge;
      /**state for which DFTB gradients are calculated (ground state = 0)*/
      int gradstate;
      /**verbosity for dftbaby*/
      int verbose;
      /**maximum number of SCF-iterations*/
      int maxiter;
      /**convergence threshold for relative change in SCF-calculation*/
      std::string conv_threshold;
      /**cutoff in bohr: orbitals that are further away don't interact*/
      float cutoff;
      /**long range correction on or off*/
      bool longrange;
      /**distance (in bohr) where long range correction is switched on*/
      float lr_dist;
      /**limit the TD-DFTB matrix (used for gradients) to the lowest ... eigenvalues*/
      int states;
      /**number of occupied orbitals taken into account for TD-DFTB*/
      int orb_occ;
      /**number of virtual orbitals taken into account for TD-DFTB*/
      int orb_virt;
      /**maximum number of iterations for TD-DFTB matrix diagonalisation*/
      int diag_maxiter;
      /**convergence threshold for TD-DFTB matrix diagonalisation*/
      std::string diag_conv;
      /**use own optimizer for optimization (otherwise steepest gradient)*/
      bool opt;

      /**constructor
      for most options if a value is set to 0, the default values from dftbaby are used
      exceptions: gradstate, verbose*/
      dftbaby_conf(void): path{""}, gradfile{"grad.xyz"}, charge{0}, gradstate{0}, verbose{0}, maxiter{0},
      conv_threshold{"0"}, cutoff{ 0.0f }, longrange{false}, lr_dist{0.0f},
      states{0}, orb_occ{0}, orb_virt{0}, diag_maxiter{0}, diag_conv{"0"}, opt{false} {}
    } dftbaby;

    /**struct that contains all information necessary for DFTB+ calculation*/
    struct dftb_conf
    {
      /**path to dftb+*/
      std::string path;
      /**path to slater-koster files*/
      std::string sk_files;
      /**verbosity for dftb+*/
      int verbosity;
      /**convergence threshold for SCC-DFTB calculation*/
      double scctol;
      /**maximum number of steps for SCC procedure*/
      int max_steps;
      /**total charge of the system*/
      int charge;
      /**use DFTB3 ?*/
      bool dftb3;
			/**number of K points in x-, y- and z-direction (only for periodic boundaries)*/
			std::vector<int> kpoints;
      /**optimizer (0 = CAST, 1 = Steepest Decent, 2 = Conjugate Gradient)*/
      int opt;
      /**maximal number of steps for optimization with DFTB+ optimizer*/
      int max_steps_opt;
			/**temperature for fermi filling (in K)*/
			double fermi_temp;
      /**constructor*/
      dftb_conf(void): verbosity(0), scctol(0.00001), max_steps(1000), charge(0),
				dftb3(false), kpoints({ 1,1,1 }), opt(2), max_steps_opt(5000), fermi_temp(0.0) {}
    } dftb;

    /**struct that contains all information necessary for ORCA calculation*/
		struct orca_conf
		{
			/**path to orca*/
			std::string path;
			/**number of processors used*/
      int nproc{ 1 };
      /**maximum amount of scratch memory per core (in MB)*/
      int maxcore{ 0 };

			/**method*/
			std::string method;
			/**basisset*/
      std::string basisset{ "" };
      /**further specifications for ORCA call*/
      std::string spec{ "" };

			/**total charge of the system*/
      int charge{ 0 };
			/**multiplicity*/
      int multiplicity{ 1 };

			/**optimizer (0 = CAST, 1 = ORCA)*/
      int opt{ 1 };
      
      /**verbosity (from 0 to 4)*/
      int verbose{ 1 };

      /**numbers of orbitals that should be plotted as cubefiles*/
      std::vector<size_t> cube_orbs;

			// stuff for casscf calculation

			/**add casscf section*/
      bool casscf{ false };
			/**number of electrons*/
			int nelec;
			/**number of orbitals*/
			int norb;
			/**number of roots*/
			int nroots;
			/**use Newton-Raphson algorithm?*/
      bool nr{ false };
			/**switch on NEVPT2?*/
      bool nevpt{ false };

			// stuff for implicit solvent (CPCM)

			/**switch on cpcm?*/
      bool cpcm{ false };
			/**dielectric constant*/
			double eps;
			/**refractive index*/
			double refrac;
			
		} orca;

    /**struct that contains all information necessary for gaussian calculation*/
    struct gaussian_conf
    {
      /**path to gaussian / command to open gaussian*/
      std::string path;
      /**link optinos for gaussian*/
      std::string link;
      /**charge of the molecule*/
      std::string charge;
      /**multiplicity of the molecule*/
      std::string multipl;
      /**method for energy calculation*/
      std::string method;
      /**basisset for energy calculation*/
      std::string basisset;
      /**further specifications for gaussian call*/
      std::string spec;
      /**name of checkpoint file*/
      std::string chk;
      /**should gaussian input be deleted after calculation?*/
      bool delete_input;
      /**should gaussian optimizer be used? (otherwise CAST optimizer)*/
      bool opt;
      /**if gaussian optimizer is used: optimize with steepest gradient? (better comparable with CAST optimization)*/
      bool steep;
      /**after this number of failed gaussian calls CAST breaks*/
      int maxfail;

			// stuff for implicit solvent (CPCM)

			/**switch on cpcm?*/
			bool cpcm;
			/**dielectric constant*/
			double eps;
			/**refractive index*/
			double epsinf;

      /**constructor*/
      gaussian_conf(void) : method{ "Hf/ " }, basisset{ "" }, spec{ "" }, chk{ "" }, delete_input { true }, opt{ true },
         steep{ true }, maxfail{1000u}, cpcm {false}
      {}
    } gaussian;

  /**struct that contains all information necessary for chemshell calculation*/
	struct chemshell_conf {
		std::string extra_pdb = "";
		std::string optional_inpcrd = "";
		std::string optional_prmtop = "";
		std::string path = "";
		std::string babel_path = "";

        std::string coords = "";
		std::string scheme = "";
		std::string qm_theory = "";
		std::string qm_ham = "";
		std::string qm_basis = "";
		std::string qm_charge = "";
		std::string qm_atoms = "";
		std::string com_residues = "";

		std::string maxcycle = "";
		std::string maxcyc = "";
		std::string tolerance = "";
		std::string mxlist = "";
		std::string cutoff = "";
        std::string scale14 = "";
        std::string active_radius = "";
        /*
        std::vector<std::string> tleap_sources;
        std::vector<std::string> tleap_loadamberparams;
        std::vector<std::string> tleap_loadoffs;
        */
		bool dispersion = false;
		bool delete_input = true;
	} chemshell;

  /**struct that contains all information necessary for PSI4 calculation*/
  struct psi4_conf{
    /**command to execute psi4*/
    std::string path = "";
    /**reserve memory (e.g. "4GB")*/
    std::string memory = "";
    /**basisset for calculation*/
    std::string basis = "";
    /**method for calculation*/
    std::string method = "";
    /**spin multiplicity for molecule*/
    std::string spin = "";
    /**charge of molecule*/
    std::string charge = "";
    std::string threads = "";
  }psi4;

  /**default constructor for struct energy*/
    energy() :
      cutoff(std::numeric_limits<double>::max()), switchdist(cutoff - 4.0),
      isotropic(true),
      remove_fixed(false),
      spackman(), mopac()
    { }
  };

   /*! Stream operator for config::energy
   *
   * Prints configuration details for the current CAST run
   * Contains: Information about main dihedrals,
   * black-/whitelists for main dihedrals,
   * Umbrella sampling information (if task == UMBRELLA),
   * Bias Potentials
   */
  std::ostream & operator << (std::ostream &, energy const &);

  /*
    ##     ##  #######  ##       ########  ##    ## ##    ##    ###
    ###   ### ##     ## ##       ##     ##  ##  ##  ###   ##   ## ##
    #### #### ##     ## ##       ##     ##   ####   ####  ##  ##   ##
    ## ### ## ##     ## ##       ##     ##    ##    ## ## ## ##     ##
    ##     ## ##     ## ##       ##     ##    ##    ##  #### #########
    ##     ## ##     ## ##       ##     ##    ##    ##   ### ##     ##
    ##     ##  #######  ######## ########     ##    ##    ## ##     ##
  */

  /**namespace for MD options that need an own struct*/
  namespace md_conf
  {
    /**integrator (velocity-verlet or beeman)*/
    struct integrators { enum T { VERLET, BEEMAN }; };

    /**options for spherical boundaries*/
    struct config_spherical
    {
      /**radius for starting the inner spherical potential*/
      double r_inner;
      /**radius for starting the outer spherical potential*/
      double r_outer;
      /**exponent for the inner spherical potential*/
      double e1;
      /**exponent for the outer spherical potential*/
      double e2;
      /**force constant for the inner spherical potential*/
      double f1;
      /**force constant for the outer spherical potential*/
      double f2;
      /**true if spherical potential is applied, false if not*/
      bool use;
      /**constructor*/
      config_spherical(void) :
        r_inner(20.0), r_outer(20.1), e1(2.0), e2(4.0),
        f1(10.0), f2(10.0), use(false)
      { }
    };

    /**contains information for one heatstep*/
    struct config_heat
    {
      /**temperature*/
      double raise;
      /**step number*/
      std::size_t offset;
      /**constructor*/
      config_heat(void) : raise(10.0), offset(100u) { }
      /**overwritten operator <:
      returns true if < is true for step number (offset)*/
      friend bool operator< (config_heat const &a, config_heat const &b) { return (a.offset < b.offset); }
      /**overwritten operator >:
      returns true if > is true for step number (offset)*/
      friend bool operator> (config_heat const &a, config_heat const &b) { return operator<(b, a); }
    };

    /**contains information for rattle algorithm*/
    struct config_rattle
    {
      /**contains information for one rattle bond*/
      struct rattle_constraint_bond
      {
        /**ideal bond length (from parameter file)*/
        double len;
        /**index of first atom (starting with 0)*/
        std::size_t a;
        /**index of second atom (starting with 0)*/
        std::size_t b;
      };
      /**maximum number of iterations in rattle algorithm*/
      std::size_t num_iter;
      /**tolerance critirion for rattle algorithm*/
      double tolerance;
      /**vectors of bonds that should be constrained*/
      std::vector<rattle_constraint_bond> specified_rattle;
      /**true if rattle algorithm is applied, false if not*/
      bool use;
      /**true all bonds with an H-atom should be constrained, false if only specified bonds*/
      bool all;
      /**use parameterfile to get constrained distances or rather define them by yourself?*/
      bool use_paramfile;
      /**distances for rattlepairs in the same order as specified rattle*/
      std::vector<double> dists;
      /**constructor*/
      config_rattle(void) : num_iter(100), tolerance(1.0e-6), use(false), all(true), use_paramfile(true)
      { }
    };
  }

  /**struct for MD options*/
  struct molecular_dynamics
  {
    /**temperature control active?*/
    bool temp_control;
    /**timestep in picoseconds*/
    double timeStep;
    /**initial temperature*/
    double T_init;
    /**final temperature*/
    double T_final;
    /**start MD again from beginning if molecule gets destroyed yes or no*/
    int broken_restart;

    //pressure things
    double pcompress, pdelay, ptarget;

    // Options for biased MD
    /**1 if a biased potential around an active site is applied, 0 if not*/
    std::size_t set_active_center;
    /**1 if the active site and the distances to the active site should be calculated new every step,
    0 if they should be calculated only once at the beginning of the simulation*/
    std::size_t adjustment_by_step;
    /**distance of inner cutoff for biased potential in angstrom*/
    double inner_cutoff;
    /**distance of outer cutoff for biased potential in angstrom*/
    double outer_cutoff;
    /**vector of atoms (tinker atom-numbers) that define the active site
    coordinates of active site are calculated as geometrical center*/
    std::vector<unsigned> active_center;

    /**number of MD steps*/
    std::size_t num_steps;
    /**number of snapshots*/
    std::size_t num_snapShots;
    /**number of snapshots in memory before written to file*/
    std::size_t max_snap_buffer;
    /**after this number of steps the list of non-bonded interactions is generated new*/
    std::size_t refine_offset;
    /**after this number of steps a restart file is generated*/
    std::size_t restart_offset;
    /**each trackoffset'th step is written to trace file*/
    std::size_t trackoffset;

    // Umbrella Sampling

    /**number of equilibration steps*/
    std::size_t usequil;
    /**offset for taking snapshots*/
    std::size_t usoffset;

    /**vector of heatsteps:
    each MDheat option is saved into one element of this vector*/
    std::vector<md_conf::config_heat> heat_steps;
    /**contains options for spherical boundaries if applied,
    otherwise the information that no spherical boundaries are applied*/
    md_conf::config_spherical spherical;
    /**contains information for rattle algorithm*/
    md_conf::config_rattle rattle;
    /**integrator that is used: VERLET (velocity-verlet) or BEEMAN (beeman) */
    md_conf::integrators::T integrator;
    /**Nos�-Hoover thermostat yes or no*/
    bool hooverHeatBath;
    /**remove translation and rotation after every step*/
    bool veloScale;
    /**free energy perturbation calculation yes or no*/
    bool fep;
    /**activate tracking yes or no*/
    bool track;
    /**perform local optimization with snapshots before they are written into file yes or no*/
    bool optimize_snapshots;
    /**pressure control yes or no?*/
    bool pressure;
    /**use a restart file for starting MD yes or no (does currently not work)*/
    bool resume;
    /**perform an umbrella sampling yes or no*/
    bool umbrella;
    /**perform local optimization before starting simulation yes or no*/
    bool pre_optimize;
    /**atom pairs to analyze*/
    std::vector<std::vector<size_t>> ana_pairs;
    /**analyze zones?*/
    bool analyze_zones;
    /**zone width (distance to active site where a new zone starts)*/
    double zone_width;
    //**scaling factor for nosehoover thermostat
    double nosehoover_Q;

    /**constructor*/
    molecular_dynamics(void) :
      temp_control{true}, timeStep{0.001}, T_init{0.0}, T_final{0.0},
      broken_restart{ 0 }, pcompress{0.000046}, pdelay{2.0}, ptarget{1.0},
      set_active_center{ 0 }, adjustment_by_step { 0 }, inner_cutoff{ 0.0 }, outer_cutoff{ 0.0 },
      active_center(), num_steps{10000}, num_snapShots{100}, max_snap_buffer{50},
      refine_offset{0}, restart_offset{0}, trackoffset{1}, usequil{0}, usoffset{ 0 }, 
       heat_steps(), spherical{}, rattle{},
      integrator(md_conf::integrators::VERLET),
      hooverHeatBath{false}, veloScale{true},  fep{false}, track{true},
      optimize_snapshots{false}, pressure{false},
      resume{false}, umbrella{false}, pre_optimize{false}, ana_pairs(), analyze_zones{false},
      zone_width{ 0.0 }, nosehoover_Q{ 0.1 }
    { }

  };

  /**contains information about FEP calculation if performed*/
  struct fep
  {
    /**final value for order parameter lambda,
    no need to set it to a value other than 1*/
    double lambda;
    /**size of the FEP windows*/
    double dlambda;
    /**controls lambda value for vdw-coupling*/
    double vdwcouple;
    /**controls lambda value for electrostatic coupling*/
    double eleccouple;
    /**value for vdw shifting parameter (softcore potential)*/
    double ljshift;
    /**value for coulomb shifting parameter (softcore potential)*/
    double cshift;
    /**number of MD steps in production run for every window*/
    std::size_t steps;
    /**number of MD steps in equilibration run for every window*/
    std::size_t equil;
    /**output frequency in alchemical.txt (does not affect calculation)*/
    std::size_t freq;
    /**perform graphical analysis?*/
    bool analyze;
    /**use Bennets acceptance ratio?*/
    bool bar;
    /**constructor*/
    fep(void) :
      lambda(1.0), dlambda(0.1), vdwcouple(1.0), eleccouple(1.0), ljshift(1.0), cshift(1.0),
      steps(10), equil(10), freq(1), analyze(true), bar(false)
    { }
  };


  /*
     #######  ########  ######## #### ##     ## #### ########    ###    ######## ####  #######  ##    ##
    ##     ## ##     ##    ##     ##  ###   ###  ##       ##    ## ##      ##     ##  ##     ## ###   ##
    ##     ## ##     ##    ##     ##  #### ####  ##      ##    ##   ##     ##     ##  ##     ## ####  ##
    ##     ## ########     ##     ##  ## ### ##  ##     ##    ##     ##    ##     ##  ##     ## ## ## ##
    ##     ## ##           ##     ##  ##     ##  ##    ##     #########    ##     ##  ##     ## ##  ####
    ##     ## ##           ##     ##  ##     ##  ##   ##      ##     ##    ##     ##  ##     ## ##   ###
     #######  ##           ##    #### ##     ## #### ######## ##     ##    ##    ####  #######  ##    ##
  */

  /**optimization options*/
  namespace optimization_conf
  {
    /**methods for local optimizations (currently only LBFGS)*/
    struct lo_types { enum T { LBFGS = 0 }; };
    /**methods for global optimizations (monte carlo with minimization, tabu-search)*/
    struct go_types { enum T { MCM, TABU }; };

    /**struct that contains configuration options for local optimisation via L-BFGS*/
    struct lo
    {
      /**convergence threshold for bfgs*/
      double grad;
      /**max number of steps for bfgs*/
      std::size_t maxstep;
      /**should trace written into file?*/
      bool trace;
      lo(void) : grad(0.001), maxstep(10000), trace(false) { }
    };

    /**struct that contains configuration options for monte-carlo*/
    struct mc
    {
      struct move_types { enum T { DIHEDRAL_OPT, DIHEDRAL, XYZ, WATER }; };
      // stepsize in cartesian space and temperature
      double cartesian_stepsize, dihedral_max_rot, move_frequency_probability;
      // move method (cartesian, dihedral or dihedral opt)
      move_types::T move;
      // use minimization after move (basin hopping / mcm)
      // tracking
      bool minimization;
      mc(void) :
        cartesian_stepsize(2.0), dihedral_max_rot(30.0),
        move_frequency_probability(0.75),
        move(move_types::DIHEDRAL), minimization(true)
      { }
    };

    /**struct that contains configuration options for tabu-search*/
    struct ts
    {
      std::size_t divers_iterations, divers_threshold, divers_limit;
      go_types::T divers_optimizer;
      bool mcm_first;
      ts(void) :
        divers_iterations(30), divers_threshold(25), divers_limit(50),
        divers_optimizer(go_types::MCM), mcm_first(false)
      { }
    };

    struct grd
    {
      double delta;
      scon::ang<double> main_delta;
      grd() : delta(3.0),
        main_delta(scon::ang<double>::from_deg(30.0)) {}
    };

    struct sel
    {
      struct fitness_types { enum T { INVALID = -1, LINEAR, EXPONENTIAL }; };
      double low_rank_fitness, high_rank_fitness;
      std::size_t included_minima;
      fitness_types::T fit_type;
      sel(void) :
        low_rank_fitness(0.5), high_rank_fitness(1.0),
        included_minima(10), fit_type(fitness_types::LINEAR)
      { }
    };

    struct evo
    {
      double chance_pointmutation,
        chance_crossingover;
      evo() :
        chance_pointmutation(0.3),
        chance_crossingover(0.4)
      { }
    };

    ///////////////////////////////////

    struct global
    {
      struct fallback_types { enum T { INVALID = -1, LAST_GLOBAL, FITNESS_ROULETTE }; };
      double temperature, temp_scale, delta_e;
      ts tabusearch;
      mc montecarlo;
      sel selection;
      evo evolution;
      grd grid;
      std::size_t iterations, fallback_limit, precision;
      fallback_types::T fallback;
      bool metropolis_local, pre_optimize, move_dehydrated;
      global(void) :
        temperature(298.15), temp_scale(1.0), delta_e(0.0),
        tabusearch(), montecarlo(), iterations(1000), fallback_limit(200), precision(4),
        fallback(fallback_types::FITNESS_ROULETTE),
        metropolis_local(true), pre_optimize(false), move_dehydrated(false)
      { }
    };

    struct local
    {
      std::ptrdiff_t method;
      lo bfgs;
      local(void) : method(lo_types::LBFGS) { }
    };

    ///////////////////////////////////

    /*! Stream operator for config::global
     *
     * Prints configuration details for global optimisation
     * Contains: Iterations, Temperature, Temperature Scaling,
     * TabuSearch-Iterations, MonteCarlo-Movetype and much more.
     */
    std::ostream& operator<< (std::ostream &, global const &);

    static std::size_t const NUM_FITNESS = 2;
    static std::string const fitness_strings[NUM_FITNESS] = { "LINEAR", "EXPONENTIAL" };

    static std::size_t const NUM_FALLBACKS = 2;
    static std::string const fallback_strings[NUM_FALLBACKS] = { "LAST_GLOBAL", "FITNESS_ROULETTE" };
  }

  struct optimization
  {
    optimization_conf::local local;
    optimization_conf::global global;
  };

  /*
    ######  ########    ###    ########  ########  #######  ########  ########
   ##    ##    ##      ## ##   ##     ##    ##    ##     ## ##     ##    ##
   ##          ##     ##   ##  ##     ##    ##    ##     ## ##     ##    ##
    ######     ##    ##     ## ########     ##    ##     ## ########     ##
         ##    ##    ######### ##   ##      ##    ##     ## ##           ##
   ##    ##    ##    ##     ## ##    ##     ##    ##     ## ##           ##
    ######     ##    ##     ## ##     ##    ##     #######  ##           ##
  */

  /**namespace for subtasks of STARTOPT (solvadd and ringsearch)*/
  namespace startopt_conf
  {
    /**config options for solveadd, i.e. solvation with water*/
    struct solvadd
    {
      // some enums for solveadd

      /**possible shapes of the water layer*/
      struct boundary_types { enum T { LAYER = 0, SPHERE = 1, BOX = 2 }; };
      /**possibilities for when an optimization will be performed
      (not at all, after every water shell, in the end, after every water shell and in the end)*/
      struct opt_types { enum T { NONE, SHELL, TOTAL, TOTAL_SHELL }; };

      // real options

      /**hydrogen bond length*/
      double defaultLenHB;
      /**size of the water layer (extended if more waters have to be added than fit into it)*/
      double maxDistance;
      /**length of bond between O and H*/
      double water_bond;
      /**angle in water*/
      ::coords::angle_type water_angle;
      /**number of water molecules to be added*/
      std::size_t maxNumWater;
      /**forcefield atom type of oxygen*/
      std::size_t ffTypeOxygen;
      /**forcefield atom type of hydrogen*/
      std::size_t ffTypeHydrogen;
      /**shape of the water layer*/
      boundary_types::T boundary;
      /**when and if will an optimization be performed?*/
      opt_types::T opt;
      /**should initial structure be fixed?*/
      bool fix_initial;
      /**???*/
      bool fix_intermediate;
      /**???*/
      globopt_routine_type::T go_type;
      /**constructor*/
      solvadd() :
        defaultLenHB(1.79), maxDistance(10.0),
        water_bond(0.95), water_angle(::coords::angle_type::from_deg(109.5)),
        maxNumWater(0), ffTypeOxygen(53), ffTypeHydrogen(54),
        boundary(boundary_types::LAYER), opt(opt_types::SHELL),
        fix_initial(true), fix_intermediate(true),
        go_type(globopt_routine_type::BASINHOPPING)
      { }
    };

    /**config options for ringsearch*/
    struct ringsearch
    {
      ::coords::float_type bias_force, chance_close;
      std::size_t population, generations;
      ringsearch(void)
        : bias_force(0.1), chance_close(0.75),
        population(12), generations(20)
      { }
    };

    /*! Stream operator for config::ringsearch
     *
     * Prints configuration details for RINGSEARCH subtask
     * of STARTOP task.
     */
    std::ostream& operator<< (std::ostream &, ringsearch const &);

    /*! Stream operator for config::solvadd
     *
     * Prints configuration details for SOLVADD subtask
     * of STARTOP task.
     */
    std::ostream& operator<< (std::ostream &, solvadd const &);
  }

  /**struct for config option for STARTOPT, contains struc solveadd and ringsearch*/
  struct startopt
  {
    /**possible subtasks*/
    struct types { enum T { RINGSEARCH, SOLVADD, RINGSEARCH_SOLVADD }; };

    /**options for solveadd*/
    startopt_conf::solvadd solvadd;
    /**options for ringsearch*/
    startopt_conf::ringsearch ringsearch;
    /**which subtask should be run?*/
    types::T type;
    /**number of structures???*/
    std::size_t number_of_structures;
    /**constructor*/
    startopt(void)
      : solvadd(), ringsearch(),
      type(types::SOLVADD),
      number_of_structures()
    { }
  };

  /*
    ########  #### ##     ## ######## ########
    ##     ##  ##  ###   ### ##       ##     ##
    ##     ##  ##  #### #### ##       ##     ##
    ##     ##  ##  ## ### ## ######   ########
    ##     ##  ##  ##     ## ##       ##   ##
    ##     ##  ##  ##     ## ##       ##    ##
    ########  #### ##     ## ######## ##     ##
  */

  /**struct for config options of task DIMER*/
  struct dimer
  {
    struct translation_types { enum T { CG, STRAIGHT }; };
    struct gradient_types { enum T { CALCULATE, EXTRAPOLATE }; };
    double distance, rotationConvergence, trans_F_rot_limit;
    std::size_t maxStep, maxRot;
    translation_types::T trans_type;
    gradient_types::T grad_type;
    dimer(void) :
      distance(0.01), rotationConvergence(5), trans_F_rot_limit(0.01),
      maxStep(100), maxRot(20),
      trans_type(translation_types::CG), grad_type(gradient_types::CALCULATE)
    { }
  };

  /*
     ####     ##    #########   ########          ##   #######
     ## ##    ##    ##          ##     ##        ##    ##    ##
     ##  ##   ##    ##          ##     ##       ##     ##    ##
     ##   ##  ##    ######      ########       ##      #######
     ##    ## ##    ##          ##     ##     ##       ##
     ##     ####    ##          ##     ##    ##        ##
     ##      ###    #########   ########    ##         ##
  */

  /**struct for config options of task NEB*/
  struct neb
  {
    std::string FINAL_STRUCTURE, OPTMODE;
    double SPRINGCONSTANT, TEMPERATURE, MCSTEPSIZE, BIASCONSTANT,
      VARIATION, PO_ENERGY_RANGE, BOND_PARAM, INT_IT;
    std::size_t IMAGES, MCITERATION, GLOBALITERATION,
      CONNECT_NEB_NUMBER, NUMBER_OF_DIHEDRALS, MCM_SAVEITER;
    bool NEB_CONN, CONSTRAINT_GLOBAL, TAU, CONN,
      MIXED_MOVE, INT_PATH, CLIMBING, IDPP, MAXFLUX, MAXFLUX_PATHOPT, COMPLETE_PATH, MULTIPLE_POINTS, INTERNAL_INTERPOLATION, MCM_OPT;
	neb() :
		FINAL_STRUCTURE{ "" }, OPTMODE("PROJECTED"),
		SPRINGCONSTANT(0.1), TEMPERATURE(298.15), MCSTEPSIZE(0.5),
		BIASCONSTANT(0.1), VARIATION(3.0), PO_ENERGY_RANGE(100.0),
		BOND_PARAM(2.2), INT_IT(0.5), IMAGES(12), MCITERATION(100),
		GLOBALITERATION(1), CONNECT_NEB_NUMBER(3), NUMBER_OF_DIHEDRALS(1),MCM_SAVEITER(1),
		NEB_CONN(false), CONSTRAINT_GLOBAL(false), TAU(true), CONN(true),MIXED_MOVE(false),
		INT_PATH(false), CLIMBING(true), IDPP(false), MAXFLUX(false), MAXFLUX_PATHOPT(false), COMPLETE_PATH(false), MULTIPLE_POINTS(false), INTERNAL_INTERPOLATION(false), MCM_OPT(true)
    {}
  };

  /**
   * ALIGN // KABSCH ALIGNMENT OF STRUCTURES
   * THIS TASK REMOVES TRANSLATION AND ROTATION
   */
  struct align
  {
    size_t dist_unit;
    size_t reference_frame_num;
    bool traj_align_translational;
    bool traj_align_rotational;
    bool traj_print_bool;
    double holm_sand_r0;
    std::string align_external_file;
    align(void) : dist_unit(0), reference_frame_num(0), traj_align_translational(true), traj_align_rotational(true), traj_print_bool(true), holm_sand_r0(20), align_external_file()
    {}
  };

  /**
  * PCA // Principal Component Analysis
  * THIS TASK PERFORMS PCA ON A TRAJECTORY
  */
  struct PCA
  {
    bool pca_alignment;
    size_t pca_ref_frame_num;
    size_t pca_start_frame_num;
    bool pca_read_vectors;
    bool pca_read_modes;
    bool pca_use_internal;
    bool pca_trunc_atoms_bool;
    bool pca_ignore_hydrogen;
    bool pca_print_probability_density;
    double pca_histogram_width;
    size_t pca_histogram_number_of_bins;
    size_t pca_offset;
    std::vector<size_t> pca_trunc_atoms_num;
    std::vector<size_t> pca_internal_dih;
    std::vector<size_t> pca_dimensions_for_histogramming;
    bool pca_histogram_all_marginal_degrees_of_freedom;
    std::vector<double> proc_desired_start;
    std::vector<double> proc_desired_stop;

    PCA(void) : pca_alignment(true), pca_ref_frame_num(0u), pca_start_frame_num(0u), pca_read_vectors(false), pca_read_modes(false),
      pca_use_internal(false), pca_trunc_atoms_bool(false), pca_ignore_hydrogen(false),
      pca_print_probability_density(true), pca_histogram_width(0.), pca_histogram_number_of_bins(32u), pca_offset(1u),
      pca_trunc_atoms_num(), pca_internal_dih(), pca_dimensions_for_histogramming(std::vector<size_t>{1u, 2u}),
      pca_histogram_all_marginal_degrees_of_freedom(false), proc_desired_start(), proc_desired_stop()

    {}
  };

  /**
  * ENTROPY // Entropy Calculations
  * THIS TASK PERFORMS CONFIGURATIONAL AND CONFORMATIONAL ENTROPY CACLULATIONS
  */
  struct entropy
  {
    bool entropy_alignment;
    double entropy_temp;
    size_t entropy_ref_frame_num;
    size_t entropy_start_frame_num;
    std::vector<size_t> entropy_method;
    size_t entropy_method_knn_k;
    bool entropy_remove_dof;
    bool entropy_use_internal;
    bool entropy_trunc_atoms_bool;
    size_t entropy_offset;
    int knnfunc;
    int knnnorm;
    std::vector<size_t> entropy_internal_dih;
    std::vector<size_t> entropy_trunc_atoms_num;
    entropy(void) : entropy_alignment(true), entropy_temp(300), entropy_ref_frame_num(0), entropy_start_frame_num(0), entropy_method(1, 6u),
      entropy_method_knn_k(4), entropy_remove_dof(true), entropy_use_internal(false), entropy_trunc_atoms_bool(false), entropy_offset(1),
      knnfunc(2), knnnorm(0), entropy_internal_dih(), entropy_trunc_atoms_num()
    {}
  };

  /**
   * IO // IO OPTIONS
   * THIS STRUCT KEEPS TRACK OF ADITIONAL IO-STUFF
   * AT THE MOMENT ONLY FOR AMBER INPUT
   */
  struct io
  {
    // mdcrd has greater priority than inpcr, inpcrd has greater priority than restrt
    std::string amber_mdcrd; //Filename
    std::string amber_mdvel; //Filename
    std::string amber_inpcrd; //Filename
    std::string amber_restrt; //Filename
    bool amber_trajectory_at_constant_pressure; // If true, box coordinates are ambigousily also written to coordinate file. We need to know this.
    io(void) : amber_mdcrd(), amber_mdvel(), amber_inpcrd(), amber_restrt(), amber_trajectory_at_constant_pressure(false) {}
  };

  /*
  2DScan Struct
  */
  struct scan2d {
	  std::vector<std::string> AXES;

      double change_from_atom_to_atom=0., max_change_to_rotate_whole_molecule=180.;
      bool constraints = false;
  };

  /*

      GBSA

  */
/*
  namespace gbsa_conf
  {
    struct method_types { enum T { VAC = -1, STILL = 0, HCT, OBC, GRYCUK, ACE, ONION, METHODNUM }; };
    struct surface_types { enum T { TINKER, SASASTILL, GAUSS, SURFACESNUM }; };
    struct radius_types { enum T { STD, VDW }; };

  }

  struct generalized_born
  {
    gbsa_conf::method_types::T method_type;
    gbsa_conf::surface_types::T surface_type;
    gbsa_conf::radius_types::T radius_type;
    generalized_born() :
      method_type(gbsa_conf::method_types::STILL),
      surface_type(gbsa_conf::surface_types::TINKER),
      radius_type(gbsa_conf::radius_types::STD)
    {}
  };
  */

  // EXCITON BREAKUP STUFF

  struct exbreak
  {
	  std::string masscenters; //Filename
	  std::string nscpairrates; //Filename
	  std::string pscpairexrates; //Filename
	  std::string pscpairchrates; //Filename
	  std::string pnscpairrates; //Filename
	  int nscnumber, pscnumber;
	  char interfaceorientation;
    double ReorgE_exc, ReorgE_ch, ReorgE_nSC, ReorgE_ct, ReorgE_rek,
       ct_triebkraft, rek_triebkraft,oscillatorstrength, wellenzahl;
  };

  struct interfcrea
  {
    std::string icfilename;
    input_types::T icfiletype;
    char        icaxis;
    double      icdist;
  };

  struct center
  {
    bool dimer;
    double distance;
  };

  struct couplings
  {
    double nbr_nSC, nbr_pSC, nbr_dimPairs;
    std::string ct_chara_all,
                pSCmultipl, pSCcharge, pSCmethod_el, pSCmethod_ex,
                nSCmultipl, nSCcharge, nSCmethod,
                hetmultipl, hetcharge, hetmethod;
  };

  struct layd
  {
    std::size_t amount, del_amount, sec_amount, sec_del_amount;
    char        laydaxis;
    double      layddist, sec_layddist;
    bool        hetero_option, replace;
    std::string layd_secname, reference1, reference2;
  };
  
  /* Constraints on internal coordinates
   */

  // Proposal: give all Constraints a shared pointer to the Internal Coorindates in order, to controll the constraint coordinates.
  //class ::InternalCoordinates::InternalCoordinate;

  enum class Constraint : int { NONE = 0, BONDS, ANGLE, DIHEDRAL, TRANSLATION_X, TRANSLATION_Y, TRANSLATION_Z, ROTATION_A, ROTATION_B, ROTATION_C };

  struct AbstractConstraint {
	  virtual Constraint getType() const = 0;
	  virtual bool isFreezed() const = 0;
	  virtual std::vector<std::size_t> const& getAtomIndices() const = 0;
  };

  struct NormalConstraint : AbstractConstraint {
	  NormalConstraint(std::vector<std::size_t> && indices, bool freeze) : atomIndices(std::move(indices)), isFreeze(freeze) {}
	  virtual std::vector<std::size_t> const& getAtomIndices() const override {
		  return atomIndices;
	  }
	  virtual bool isFreezed() const override { return isFreeze; }
	  std::vector<std::size_t> atomIndices;
	  bool isFreeze;
  };

  struct BondConstraint : NormalConstraint {
	  using NormalConstraint::NormalConstraint;
	  virtual Constraint getType() const override {
		  return Constraint::BONDS;
	  }
  };

  struct AngleConstraint : NormalConstraint {
	  using NormalConstraint::NormalConstraint;
	  virtual Constraint getType() const override {
		  return Constraint::ANGLE;
	  }
  };

  struct DihedralConstraint : NormalConstraint {
	  using NormalConstraint::NormalConstraint;
	  virtual Constraint getType() const override {
		  return Constraint::DIHEDRAL;
	  }
  };

  struct TranslationXConstraint : NormalConstraint {
	  using NormalConstraint::NormalConstraint;
	  virtual Constraint getType() const override {
		  return Constraint::TRANSLATION_X;
	  }
  };

  struct TranslationYConstraint : NormalConstraint {
	  using NormalConstraint::NormalConstraint;
	  virtual Constraint getType() const override {
		  return Constraint::TRANSLATION_Y;
	  }
  };

  struct TranslationZConstraint : NormalConstraint {
	  using NormalConstraint::NormalConstraint;
	  virtual Constraint getType() const override {
		  return Constraint::TRANSLATION_Z;
	  }
  };

  struct RotationAConstraint : NormalConstraint {
	  using NormalConstraint::NormalConstraint;
	  virtual Constraint getType() const override {
		  return Constraint::ROTATION_A;
	  }
  };

  struct RotationBConstraint : NormalConstraint {
	  using NormalConstraint::NormalConstraint;
	  virtual Constraint getType() const override {
		  return Constraint::ROTATION_B;
	  }
  };

  struct RotationCConstraint : NormalConstraint {
	  using NormalConstraint::NormalConstraint;
	  virtual Constraint getType() const override {
		  return Constraint::ROTATION_C;
	  }
  };


  struct constrained_internals
  {
    using constrain_vec = std::vector<std::shared_ptr<AbstractConstraint>>;
    
	constrained_internals() : constrain_bond_lengths{ false }, constrain_bond_angles{ false }, constrain_dihedrals{ false }, constrain_out_of_plane_bends{ false },
		constrain_translations{ false }, constrain_rotations{ false }, constrains{std::make_shared<constrain_vec>()}{};

    bool constrain_bond_lengths,
         constrain_bond_angles,
         constrain_dihedrals,
         constrain_out_of_plane_bends,
         constrain_translations,
         constrain_rotations;
    
    // Information on individual coordinates
    std::shared_ptr<constrain_vec> constrains;

	void handleConstraintInput(std::istringstream &);
  };

  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////

  //... now for some
  // important functions tinkering with the config
  // or reading it from file....


  /**
  * Helperfunction that matches a string
  * to an enum via a sorted array of strings
  *
  * @typename enum_type: Type of the enum which is to be returned
  * @param SIZE: Size if the sorted string array used for matching ("valarray")
  * @param valarray: Sorted array of strings in same order as target enum
  * @param value: Input string that is to be "converted" to enum value.
  */
  template<class enum_type, std::size_t SIZE,
    class CharT, class TraitT, class AllocT>
    inline enum_type enum_from_string(
      std::basic_string<CharT, TraitT, AllocT> const valarray[SIZE],
      std::basic_string<CharT, TraitT, AllocT> const & value)
  {
    for (std::size_t i(0U); i < SIZE; ++i)
    {
      if (value == valarray[i])
      {
        return static_cast<enum_type>(i);
      }
    }
    return static_cast<enum_type>(-1);
  }

  /*! Stream operator for config::general
  *
  * Prints contents of config::general in human
  * readable form. This contains: Where structure was read from,
  * which inputtype the structure originates from, where the
  * (forcefield) parameters originate from and which
  * energy interface is used.
  *
  * @todo: Remove line explaining parameterfile if MOPAC or TERACHEM energy interface is chosen
  */
  std::ostream & operator<< (std::ostream &, general const &);

  /*! Parses command line switches into cnofig object
  *
  * This function parses command lines switches.
  * They have priority over options from the inputfile
  * and therefore overwrite them. Pass over argc
  * and argv to this function.
  *
  * @param N: usually "argc"
  * @param V: usually "argv"
  */
  void parse_command_switches(std::ptrdiff_t const, char**);

  /*! Returns name of the config-file for the runtime of CAST
  *
  * This function determines the name
  * of the INPUTFILE which is to be read for
  * config-options. Default is either "CAST.txt"
  * or "INPUTFILE".
  * By starting the CAST executable with commandline option
  * -s or -setup the filename of a different
  * inputfile can be specified.
  *
  * Call like this:
  * CAST.exe -setup=filename.txt
  * or
  * CAST.exe -s filename.txt
  */
  std::string config_file_from_commandline(std::ptrdiff_t const, char**);

  /*! Parses one config-option and stores it in config-class
  *
  * This function parses one configoption
  * and puts the value into the corresponding
  * struct inside the Config class
  *
  * @param option: name of the configoption
  * @param value_string: corresponding value of the option
  */
  void parse_option(std::string const option, std::string const value);

  // Important function declarations end here...

  //... now some stream operators


  /*! Stream operator for config::eqval
   *
   * Prints reasoning for considering two structures
   * equal (important for TabuSearch etc.)
   */
  std::ostream & operator << (std::ostream &, coords::conditionsForStructuresToBeConsideredEqual const &);

  /*! Stream operator for config::coords
   *
   * Prints configuration details for the current CAST run
   * Contains: Information about main dihedrals,
   * black-/whitelists for main dihedrals,
   * Umbrella sampling information (if task == UMBRELLA),
   * Bias Potentials
   */
  std::ostream & operator << (std::ostream &, coords const &);

  /*! Stream operator for config::startopt
   *
   * Prints configuration details for STARTOP task,
   * mainly by using the ostream operators for
   * config::ringsearch and config::solvadd.
   */
  std::ostream& operator<< (std::ostream &, startopt const &);

}

//////////////////////////////////////
//////////////////////////////////////
//////////////////////////////////////
//////////////////////////////////////

/*! Class containing the global configuration options
 *
 * Constructor automatically reads options from file.
 * Members are public since access functions get() and set()
 * return pointer to singular static global obejct of this class
 *
 * @Note: There can only ever be one Config object! A pointer
 * to this object is contained in the private member m_instance
 */
class Config
{
public:
  config::general               general;
  config::coords                coords;
  config::energy                energy;
  config::startopt              startopt;
  config::optimization          optimization;
  config::fep                   fep;
  config::molecular_dynamics    md;
  config::dimer                 dimer;
  config::neb					          neb;
  config::align			            alignment;
  config::PCA					          PCA;
  config::entropy				        entropy;
  config::io                    io;
  config::scan2d				      	scan2d;
  config::exbreak				        exbreak;
  config::interfcrea            interfcrea;
  config::center                center;
  config::couplings             couplings;
  config::periodics             periodics;
  config::layd                  layd;
  config::constrained_internals constrained_internals;
  config::stuff                 stuff;

  /*! Constructor of Config object
   *
   * During construction the configuration-file
   * will be automatically parsed and the members
   * of the config-object will be filled in accordingly.
   *
   * @param filename: Filename of the configuration file to be read. Needs to be in the same folder as the CAST executable.
   */
  Config(std::string const &filename)
  {
    // There can only ever be one Config!
    if (m_instance) throw std::runtime_error("Configuration duplication.");

    m_instance = this;
    parse_file(filename);
  }

  /*! Obtain contents of Config
   *
   * This get() function is used to safely
   * obtain the contents of the global Config instance
   *
   * This function returns const& and can
   * therefore not be used to change values.
   */
  static Config const & get()
  {
    if (!m_instance) throw std::runtime_error("Configuration not loaded.");
    return *m_instance;
  }

  /*! Change contents of Config
  *
  * This set() function is used to alter
  * the contents of the global Config instance
  *
  * This function returns a non-const reference
  * and can therefore be used to change values.
  * To merely obtain read-access, use get() function.
  */
  static Config & set()
  {
    if (!m_instance) throw std::runtime_error("Configuration not loaded.");
    return *m_instance;
  }

  /**
   * Helper function that matches a task
   * as string to the corresponding enum via
   * the sorted "helper-array" config::task_strings
   * If you add a new task, add it to both the enum
   * config::tasks::T and config::task_strings
   *
   * @param S: task as string
   */
  static config::tasks::T            getTask(std::string const&);

  /**
   * Helper function that matches an energy interface
   * as string to the corresponding enum via
   * the sorted "helper-array" config::config::interface_strings
   * If you add a new energy interface, add it to both the enum
   * config::interface_types::T and config::interface_strings
   *
   * @param S: energy interface as string
   */
  static config::interface_types::T  getInterface(std::string const&);

  /*
   * Helper function that matches an outputformat
   * as string to the corresponding enum via
   * the sorted "helper-array" config::output_strings
   * If you add a new output type, add it to both the enum
   * config::output_types::T and config::output_strings
   *
   * @param S: output-type as string
   */
  static config::output_types::T     getOutFormat(std::string const&);

private:

  /*! Parse whole config-file for config-options
   *
   * This function parses a configuration file
   * and puts the options into the Config class
   *
   * @param filename: Full filename of the file
   */
  void parse_file(std::string const &filename);

  /*! Pointer to the single instance of the Config class
   *
   * There can only ever be one Config object.
   * A pointer to it is contained here.
   * If no object exists (yet), this will be a nullpointer.
   */
  static Config * m_instance;
};

#endif
