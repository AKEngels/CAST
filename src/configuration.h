/**
CAST 3
configuration.h
Purpose: class for extraction of information from inputfile

@author Daniel Weber (modified by many)
@version 1.1
*/

#pragma once 
#include <cstddef>
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
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

#include "scon.h"
#include "filemanipulation.h"
#include "scon_utility.h"
#include "scon_vect.h"
#include "coords_rep.h"
#include "configurationHelperfunctions.h"

/*! Namespace containing relevant configuration options
 */
namespace config
{

  // Here we find some static members that only
  // exist once in CAST, like the version number or
  // some helper arrays containing the tasks etc.

  /** Name of the program*/
  static std::string const Programname("CAST");
  /** Version-Number of CAST*/
  static std::string const Version("3.2.0.2dev");

  /**Number of tasks*/
  static std::size_t const NUM_TASKS = 23;
  /** Names of all CAST tasks as strings*/
  static std::string const task_strings[NUM_TASKS] =
  { 
    "SP", "GRAD", "TS", "LOCOPT", "REMOVE_EXPLICIT_WATER",
    "MC", "DIMER", "MD", "NEB", "GOSOL", 
    "STARTOPT",  "INTERNAL", "ENTROPY", "PCAgen", "PCAproc",
    "DEVTEST", "ADJUST", "UMBRELLA", "FEP", "PATHOPT",
    "GRID", "ALIGN", "PATHSAMPLING", 
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
      DEVTEST, ADJUST, UMBRELLA, FEP, PATHOPT,
      GRID, ALIGN, PATHSAMPLING, 
    };
  };

  /** number of Input Types */
  static std::size_t const NUM_INPUT = 2;
  /** Input Types */
  static std::string const input_strings[NUM_INPUT] =
  { 
    "TINKER", "AMBER" 
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
      TINKER, AMBER 
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
  static std::size_t const NUM_INTERFACES = 6;
  /**Interface Types*/
  static std::string const 
    interface_strings[NUM_INTERFACES] =
  { 
    "AMBER", "AMOEBA", "CHARMM22", "OPLSAA", "TERACHEM", "MOPAC" 
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
      AMBER, AMOEBA, CHARMM22, OPLSAA, TERACHEM, MOPAC 
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
    /**Energy interface used for current run*/
    interface_types::T energy_interface;
    /**Energy interface used pre-optimization performed before the current run*/
    interface_types::T preopt_interface;
    /**Verbosity of the output of CAST (supposed to be between 0 and 5)*/
    std::size_t verbosity;
    /// Constructor with reasonable default parameters
    general(void) :
      paramFilename("oplsaa.prm"), outputFilename("%i.out"),
      input(input_types::TINKER), output(output_types::TINKER),
      task(config::tasks::SP), energy_interface(interface_types::OPLSAA),
      preopt_interface(interface_types::ILLEGAL),
      verbosity(1U)
    { }
  };


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
		 /**???*/
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
		 /**???*/
		 double value;
		 /**number of one atom*/
		 std::size_t a;
		 /**number of next atom*/
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
	   /**???*/
	   ::coords::angle_type value;
	   /**atom 1*/
	   std::size_t a;
	   /**atom 2*/
	   std::size_t b;
	   /**atom 3*/
	   std::size_t c;
	   /**atom 4*/
	   std::size_t d;
	   /**???*/
       bool forward;
	   /**constructor*/
       dihedral(void)
         : force(), ideal(), value(),
         a(), b(), c(), d(), forward(false)
       { }
     };
       /**sperical potential - prevents non-bonded systems from exploding*/
     struct spherical
     {
		 /**distance to center where the additional potential starts*/
		 double radius;
		 /**force constant*/
		 double force;
		 /**exponent of the potential function, 2 for harmonic potential, 4 is also possible*/
		 double exponent;
		 /**constructor*/
       spherical()
         : radius(), force(), exponent()
       { }
     };
     /**cubic potential (similar to spherical but cubic)*/
     struct cubic
     {
		 /**???*/
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
    /**vector with amber charges (only filled if AMBER input is used)*/
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
      struct umbrella_tor
      {
        double force, angle;
        std::size_t index[4U];
        bool fix_all_torsions;
        umbrella_tor(void) :
          force(0.0), index(), fix_all_torsions(false) { }
      };
      struct umbrella_dist
      {
        double force, dist;
        std::size_t index[2U];
        umbrella_dist(void) :
          force(0.0), index() { }
      };
      std::vector<umbrella_tor> torsions;
      std::vector<umbrella_dist> distances;
      std::size_t steps, snap_offset;
      umbrellas(void) : steps(50), snap_offset(10) { }

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
	  /**biased pot on torsions for umbrella sampling*/
      std::vector<config::coords::umbrellas::umbrella_tor> utors;
	  /**biased pot on bonds for umbrella sampling*/
      std::vector<config::coords::umbrellas::umbrella_dist> udist;
    } bias;
	/**???*/
    struct eqval
    {
      double superposition;
      ::coords::main_type main;
      ::coords::internal_type intern;
      ::coords::Cartesian_Point xyz;
      eqval() :
        superposition(0.4), main(::coords::angle_type::from_deg(8.0)),
        intern(0.2, ::coords::angle_type::from_deg(1.0), ::coords::angle_type::from_deg(8.0)),
        xyz(0.1, 0.1, 0.1)
      {}
    } equals;
	/**vector with numbers of fixed atoms (i.e. these atoms are not allowed to move)*/
    std::vector<std::size_t> fixed;
	/**vector with subsystems*/
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


  //
  // TASK ADJUST
  //

  namespace adjust_conf
  {
    struct dihedral
    {
      std::size_t a, b, c, d;
      ::coords::angle_type value;
    };
  }

  struct adjust
  {
    std::vector<adjust_conf::dihedral> dihedrals;
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

    double cutoff, switchdist;
    scon::c3<double> pb_box;
    bool isotropic, pme, periodic, periodic_print, remove_fixed;


    struct spack
    {
      bool on, interp;
      double cut;
      spack(void) : on(false), interp(true), cut(10.0) { }
    } spackman;

    struct mopac_conf
    {
      std::string command, path;
      mopac_ver_type::T version;
      bool delete_input;
      mopac_conf(void) : command("PM7 MOZYME"),
#if defined(MOPAC_EXEC_PATH)
        path(MOPAC_EXEC_PATH)
#elif defined(_MSC_VER)
        path("\"C:\\Program Files\\mopac\\MOPAC2012.exe\""),
#else
        path("/opt/mopac/MOPAC2012.exe"),
#endif
        version(mopac_ver_type::T::MOPAC2012MT),
        delete_input(true)
      {}
    } mopac;

    energy() :
      cutoff(10000.0), switchdist(cutoff - 4.0),
      pb_box(10.0, 10.0, 10.0), isotropic(true),
      pme(false), periodic(false), periodic_print(false), remove_fixed(false),
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
    struct integrators { enum T { VERLET, BEEMAN}; };

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
		/**number of H-atom a (number from tinker file - 1)*/
		std::size_t a;
		/**number of atom b (number from tinker file - 1)*/
		std::size_t b;
      };
	  /**???*/
      std::size_t num_iter;
	  /**???*/
      double tolerance;
	  /**vectors of bonds that should be constrained*/
      std::vector<rattle_constraint_bond> specified_rattle;
	  /**true if rattle algorithm is applied, false if not*/
	  bool use;
	  /**true all bonds with an H-atom should be constrained, false if only specified bonds*/
	  bool all;
	  /**name of parameter file where bond lengths for constrained bonds are taken from*/
      std::string ratpar;
	  /**constructor*/
      config_rattle(void) : num_iter(100), tolerance(1.0e-6), use(false), all(true)
      { }
    };
  }

  struct molecular_dynamics
  {
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
    std::size_t usoffset, usequil;

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
	/**constructor*/
    molecular_dynamics(void) :
      timeStep(0.001), T_init(0.0), T_final(),
      pcompress(0.000046), pdelay(2.0), ptarget(1.0),
      num_steps(10000), num_snapShots(100), max_snap_buffer(50),
      refine_offset(0), restart_offset(0), usequil(), usoffset(),
      trackoffset(1), heat_steps(), spherical(), rattle(),
      integrator(md_conf::integrators::VERLET), 
      hooverHeatBath(false), veloScale(false), fep(false), track(true),
      optimize_snapshots(false), pressure(false),
      resume(false), umbrella(false), pre_optimize(false)
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
	  /**every freq'th MD step of production run is taken into account for energy calculation*/
	  std::size_t freq;
    /**constructor*/
    fep(void) :
      lambda(1.0), dlambda(0.1), vdwcouple(1.0), eleccouple(1.0), ljshift(1.0), cshift(1.0),
      steps(10), equil(10), freq(1000)
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

  namespace optimization_conf
  {
    struct lo_types { enum T { LBFGS = 0 }; };
    struct go_types { enum T { MCM, TABU }; };

    struct lo
    {
      double grad;
      std::size_t maxstep;
      lo(void) : grad(0.001), maxstep(10000) { }
    };

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

  namespace startopt_conf
  {

    struct solvadd
    {
      struct boundary_types { enum T { LAYER = 0, SPHERE = 1, BOX = 2 }; };
      struct opt_types { enum T { NONE, SHELL, TOTAL, TOTAL_SHELL }; };
      double defaultLenHB, maxDistance, water_bond;
      ::coords::angle_type water_angle;
      std::size_t maxNumWater, ffTypeOxygen, ffTypeHydrogen;
      boundary_types::T boundary;
      opt_types::T opt;
      bool fix_initial, fix_intermediate;
      globopt_routine_type::T go_type;
      solvadd() :
        defaultLenHB(1.79), maxDistance(10.0),
        water_bond(0.95), water_angle(::coords::angle_type::from_deg(109.5)),
        maxNumWater(0), ffTypeOxygen(53), ffTypeHydrogen(54),
        boundary(boundary_types::LAYER), opt(opt_types::SHELL),
        fix_initial(true), fix_intermediate(true),
        go_type(globopt_routine_type::BASINHOPPING)
      { }
    };

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

  struct startopt
  {
    struct types { enum T { RINGSEARCH, SOLVADD, RINGSEARCH_SOLVADD }; };
    startopt_conf::solvadd solvadd;
    startopt_conf::ringsearch ringsearch;
    types::T type;
    std::size_t number_of_structures;
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

  struct neb
  {
    std::string FINAL_STRUCTURE, OPTMODE;
    double SPRINGCONSTANT, TEMPERATURE, MCSTEPSIZE, BIASCONSTANT, 
           VARIATION, PO_ENERGY_RANGE, BOND_PARAM, INT_IT;
    std::size_t IMAGES, MCITERATION, GLOBALITERATION, 
                CONNECT_NEB_NUMBER, NUMBER_OF_DIHEDRALS;
    bool NEB_CONN, CONSTRAINT_GLOBAL, TAU, 
         MIXED_MOVE, INT_PATH, CLIMBING, IDPP,MAXFLUX;
    neb() :
      OPTMODE("PROJECTED"),
      SPRINGCONSTANT(0.1), TEMPERATURE(298.15), MCSTEPSIZE(0.5),
      BIASCONSTANT(0.1), VARIATION(3.0), PO_ENERGY_RANGE(100.0), 
      BOND_PARAM(2.2), INT_IT(0.5), IMAGES(12), MCITERATION(100),
      GLOBALITERATION(1), CONNECT_NEB_NUMBER(3), NUMBER_OF_DIHEDRALS(1),
      NEB_CONN(false), CONSTRAINT_GLOBAL(false), TAU(true), MIXED_MOVE(false), 
      INT_PATH(false), CLIMBING(true), IDPP(false), MAXFLUX(false)
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
    std::vector<double> proc_desired_start;
    std::vector<double> proc_desired_stop;

    PCA(void) : pca_alignment(true), pca_ref_frame_num(0u), pca_start_frame_num(0u), pca_read_vectors(false), pca_read_modes(false),
       pca_use_internal(false), pca_trunc_atoms_bool(false), pca_ignore_hydrogen(false),
      pca_print_probability_density(true), pca_histogram_width(0.), pca_histogram_number_of_bins(32u), pca_offset(1u), 
      pca_trunc_atoms_num(), pca_internal_dih(), pca_dimensions_for_histogramming(std::vector<size_t>{1u, 2u}),
      proc_desired_start(), proc_desired_stop()

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
    std::vector<size_t> entropy_internal_dih;
    std::vector<size_t> entropy_trunc_atoms_num;
    entropy(void) : entropy_alignment(true), entropy_temp(300), entropy_ref_frame_num(0), entropy_start_frame_num(0), entropy_method(1, 6u),
      entropy_method_knn_k(4), entropy_remove_dof(true), entropy_use_internal(false), entropy_trunc_atoms_bool(false), entropy_offset(1),
      entropy_internal_dih(), entropy_trunc_atoms_num()
    {}
  };

  /**
   * IO // IO OPTIONS
   * THIS STRUCT KEEPS TRACK OF ADITIONAL IO-STUFF
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

      GBSA

  */

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
   * equal (important for TaboSearch etc.)
   */
  std::ostream & operator << (std::ostream &, coords::eqval const &);

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
  config::generalized_born      gbsa;
  config::adjust                adjustment;
  config::align			            alignment;
  config::PCA					          PCA;
  config::entropy				        entropy;
  config::io                    io;

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
