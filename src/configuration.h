/**
CAST 3
configuration.h
Purpose: class for extraction of information from parameter file

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

namespace config
{
  // Program Name and Version
  static std::string const Programname("CAST");
  static std::string const Version("3.2.0.1.0.0dev");

  // Tasks
  static std::size_t const NUM_TASKS = 30;
  static std::string const 
    task_strings[NUM_TASKS] =
  { 
    "SP", "GRAD", "TS", "LOCOPT", "REMOVE_EXPLICIT_WATER"
    "MC", "DIMER", "MD", "NEB", "GOSOL", 
    "STARTOPT",  "INTERNAL", "ENTROPY", "PCAgen", "PCAproc",
    "DEVTEST", "ADJUST", "UMBRELLA", "FEP", "PATHOPT",
    "GRID", "ALIGN", "PATHSAMPLING", 
    
  };
  struct tasks
  {
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

  // Input Types
  static std::size_t const NUM_INPUT = 2;
  static std::string const 
    input_strings[NUM_INPUT] =
  { 
    "TINKER", "AMBER" 
  };
  struct input_types 
  { 
    enum T 
    { 
      ILLEGAL = -1, 
      TINKER, AMBER 
    }; 
  };

  // Output Types
  static std::size_t const NUM_OUTPUT = 4;
  static std::string const 
    output_strings[NUM_OUTPUT] =
  { 
    "TINKER", "XYZ", "MOLDEN", "ZMATRIX" 
  };
  struct output_types 
  { 
    enum T 
    {
      ILLEGAL = -1, 
      TINKER, XYZ, MOLDEN, ZMATRIX
    };
  };

  // Interface Types
  static std::size_t const NUM_INTERFACES = 6;
  static std::string const 
    interface_strings[NUM_INTERFACES] =
  { 
    "AMBER", "AMOEBA", "CHARMM22", "OPLSAA", "TERACHEM", "MOPAC" 
  };
  struct interface_types 
  { 
    enum T 
    { 
      ILLEGAL = -1, 
      AMBER, AMOEBA, CHARMM22, OPLSAA, TERACHEM, MOPAC 
    }; 
  };

  // Mopac Versions
  static std::size_t const NUM_MOPAC_VERSION = 4;
  static std::string const 
    mopac_ver_string[NUM_MOPAC_VERSION] = 
  { 
    "2012", "2012MT", "7", "AVOID_HB" 
  };
  struct mopac_ver_type 
  { 
    enum T 
    { 
      ILLEGAL = -1, 
      MOPAC2012, MOPAC2012MT, MOPAC7, MOPAC7_HB 
    }; 
  };

  // Global optimization routines
  static std::size_t const NUM_GLOBOPT_ROUTINES = 2;
  static std::string const 
    globopt_routines_str[NUM_GLOBOPT_ROUTINES] =
  { 
    "TS", "BH" 
  };
  struct globopt_routine_type 
  {
    enum T 
    { 
      ILLEGAL = -1, 
      TABUSEARCH, BASINHOPPING 
    }; 
  };

  // Implicit solvation method types
  static std::size_t const NUM_SOLV = 7;
  static std::string const 
    solv_strings[NUM_SOLV] =
  { 
    "VAC", "STILL", "HCT", "OBC", "GRYCUK", "ACE", "ONION" 
  };
  struct solvs
  {
    enum S 
    {
      ILLEGAL = -1, 
      VAC, STILL, HCT, OBC, GRYCUK, ACE, ONION 
    };
  };
  // implicit solvation surface types
  static std::size_t const NUM_SURF = 3;
  static std::string const 
    surf_strings[NUM_SURF] =
  { 
    "TINKER", "SASASTILL", "GAUSS" 
  };
  struct surfs
  {
    enum SA 
    { 
      ILLEGAL = -1, 
      TINKER, SASASTILL, GAUSS 
    };
  };

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


  struct requirements
  {
    bool req_parameter, got_input_structure, 
      got_energy_interface, got_parameters, config_file, got_task;
    requirements(void) :
      req_parameter(true), got_input_structure(false), got_energy_interface(false),
      got_parameters(false), config_file(false), got_task(false)
    { }
  };

  struct general
  {
    std::string inputFilename, paramFilename, outputFilename;
    input_types::T input;
    output_types::T output;
    config::tasks::T task;
    interface_types::T energy_interface, preopt_interface;
    std::size_t verbosity, profile_runs;
    std::ofstream * trackstream;
    config::solvs::S solvationmethod;
    config::surfs::SA surfacemethod;
    bool forcefield;
    general(void) :
      paramFilename("oplsaa.prm"), outputFilename("%i.out"),
      input(input_types::TINKER), output(output_types::TINKER),
      task(config::tasks::SP), energy_interface(interface_types::OPLSAA), 
      preopt_interface(interface_types::ILLEGAL),
      verbosity(1U), profile_runs(10U), trackstream(nullptr),
      solvationmethod(solvs::VAC), surfacemethod(surfs::TINKER), forcefield(true)
    { }
    void print(void);
  };

  std::ostream & operator << (std::ostream &, general const &);

  struct rdf
  {
    double width;
    rdf(void) : width(0.0) { }
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


  namespace biases
  {
    struct potential_types { enum T { QUADRATIC, BIQUADRATIC, PROGRESSIVE }; };
    struct distance
    {
      double force, ideal, value;
      std::size_t a, b;
      distance(void)
        : force(), ideal(), a(), b()
      { }
    };
    struct angle
    {
      double force, ideal, value;
      std::size_t a, b, c;
      angle(void)
        : force(), ideal(), a(), b()
      { }
    };
    struct dihedral
    {
      double force;
      ::coords::angle_type ideal, value;
      std::size_t a, b, c, d;
      bool forward;
      dihedral(void)
        : force(), ideal(), value(),
        a(), b(), forward(false)
      { }
    };
    inline std::ostream& operator<< (std::ostream & o, dihedral const &d)
    {
      o << d.a << ',' << d.b << ',' << d.c << ',' << d.d << ' ';
      o << d.force << ',' << d.ideal << ',' << d.value << ',' << d.forward;
      return o;
    }
    struct spherical
    {
      double radius, force, exponent;
      spherical()
        : radius(), force(), exponent()
      { }
    };
    struct cubic
    {
      ::coords::Cartesian_Point dim;
      double force, exponent;
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

  struct coords
  {

    struct internals
    {
      std::map<std::size_t, std::size_t> connect;
      std::vector<std::pair<std::size_t, std::size_t>> main_whitelist;
      std::vector<std::pair<std::size_t, std::size_t>> main_blacklist;
    } internal;


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
    struct coord_bias
    {
      std::vector<biases::distance>  distance;
      std::vector<biases::angle>     angle;
      std::vector<biases::dihedral>  dihedral;
      std::vector<biases::spherical> spherical;
      std::vector<biases::cubic>     cubic;
      std::vector<config::coords::umbrellas::umbrella_tor> utors;
      std::vector<config::coords::umbrellas::umbrella_dist> udist;
    } bias;
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
    std::vector<std::size_t> fixed;
    std::vector<std::vector<std::size_t>> subsystems;
    bool remove_hydrogen_rot, no_hydrot_mains, decouple_internals, nearest_internals;
    coords(void) :
      internal(), umbrella(), bias(), equals(), fixed(), subsystems(),
      remove_hydrogen_rot(true), no_hydrot_mains(false),
      decouple_internals(false), nearest_internals(false)
    {}

  };
  std::ostream & operator << (std::ostream &, coords::eqval const &);
  std::ostream & operator << (std::ostream &, coords const &);

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

  struct energy
  {

    double cutoff, switchdist, pb_cut, pmetresh;
    scon::c3<double> pb_box;
    int pmespline;
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
      cutoff(10000.0), switchdist(cutoff - 4.0), pb_cut(9.0), pmetresh(),
      pb_box(10.0, 10.0, 10.0), pmespline(), isotropic(true),
      pme(false), periodic(false), periodic_print(false), remove_fixed(false),
      spackman(), mopac()
    { }
  };

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

  namespace md_conf
  {
    struct integrators { enum T { VERLET, BEEMAN}; };

    struct config_spherical
    {
      double r_inner, r_outer, e1, e2, f1, f2;
      bool use;
      config_spherical(void) :
        r_inner(20.0), r_outer(20.1), e1(2.0), e2(4.0),
        f1(10.0), f2(10.0), use(false)
      { }
    };

    struct heat
    {
      struct point
      {
        std::size_t iteration;
        double temperature;
      };
      double raise;
      std::size_t offset;
      std::vector<point> points;
    };

    struct config_heat
    {
      double raise;
      std::size_t offset;
      config_heat(void) : raise(10.0), offset(100u) { }
      friend bool operator< (config_heat const &a, config_heat const &b) { return (a.offset < b.offset); }
      friend bool operator> (config_heat const &a, config_heat const &b) { return operator<(b, a); }
    };

    struct config_rattle
    {
      struct rattle_constraint_bond
      {
        double len;
        std::size_t a, b;
      };
      std::size_t num_iter;
      double tolerance;
      std::vector<rattle_constraint_bond> specified_rattle;
      bool use, all;
      std::string ratpar;
      config_rattle(void) : num_iter(100), tolerance(1.0e-6), use(false), all(true)
      { }
    };
  }

  struct fep
  {
    double lambda, dlambda, vdwcouple, eleccouple, ljshift, cshift;
    std::size_t steps, equil, freq, backward;
    bool couple;
    fep(void) :
      lambda(1.0), dlambda(0.0), vdwcouple(1.0), eleccouple(1.0), ljshift(1.0), cshift(1.0),
      steps(10), equil(10), freq(1000), backward(0), couple(false)
    { }
  };

  struct molecular_dynamics
  {
	  double timeStep, T_init, T_final, pcompress, pdelay, ptarget;
	unsigned set_active_center, adjustment_by_step;
	double inner_cutoff, outer_cutoff;
    std::size_t num_steps, num_snapShots, max_snap_buffer, refine_offset, restart_offset, usequil, usoffset, trackoffset;
    std::vector<md_conf::config_heat> heat_steps;
	std::vector<unsigned> active_center;
    md_conf::config_spherical spherical;
    md_conf::config_rattle rattle;
    md_conf::integrators::T integrator;
    bool fix, hooverHeatBath, veloScale, fep, track, silent, optimize_snapshots, pressure, resume, umbrella, pre_optimize;
    molecular_dynamics(void) :
      timeStep(0.001), T_init(293.15), T_final(293.15),
      pcompress(0.000046), pdelay(2.0), ptarget(1.0),
      num_steps(10000), num_snapShots(100), max_snap_buffer(50),
      refine_offset(200), restart_offset(5000), usequil(), usoffset(),
      trackoffset(1), heat_steps(), spherical(), rattle(),
      integrator(md_conf::integrators::VERLET), fix(false),
      hooverHeatBath(true), veloScale(false), fep(false), track(true),
      silent(false), optimize_snapshots(false), pressure(false),
      resume(false), umbrella(false), pre_optimize(false)
    { }

  };

  /*
    ########     ###    ######## ##     ##
    ##     ##   ## ##      ##    ##     ##
    ##     ##  ##   ##     ##    ##     ##
    ########  ##     ##    ##    #########
    ##        #########    ##    ##     ##
    ##        ##     ##    ##    ##     ##
    ##        ##     ##    ##    ##     ##
  */

  struct path
  {
    double maxDeltaE, maxDeltaX;
    std::string endpointFileName;
    path(void) :
      maxDeltaE(2.0), maxDeltaX(1.0), endpointFileName("PATH_END.xyz") { }
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
    //struct go_move  { enum mm { CARTESIAN, DIHEDRAL, DIHEDRAL_OPT }; };

    struct lo
    {
      double grad;
      std::size_t maxstep;
      lo(void) : grad(0.001), maxstep(10000) { }
    };

    struct local
    {
      std::ptrdiff_t method;
      lo bfgs;
      local(void) : method(lo_types::LBFGS) { }
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

    static std::size_t const NUM_FITNESS = 2;
    static std::string const fitness_strings[NUM_FITNESS] = { "LINEAR", "EXPONENTIAL" };

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

    std::ostream& operator<< (std::ostream &, sel const &);

    struct evo
    {
      double chance_pointmutation,
        chance_crossingover;
      evo() :
        chance_pointmutation(0.3),
        chance_crossingover(0.4)
      { }
    };

    static std::size_t const NUM_FALLBACKS = 2;
    static std::string const fallback_strings[NUM_FALLBACKS] = { "LAST_GLOBAL", "FITNESS_ROULETTE" };

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
    std::ostream& operator<< (std::ostream &, global const &);
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

    struct fold
    {
      struct helices { enum { alpha, threeten }; };
      struct sheets { enum { betaAParallel, betaParallel }; };
      struct turns
      {
        enum {
          turnBetaIa, turnBetaIb, turnBetaIIa,
          turnBetaIIb, turnBetaVIa, turnBetaVIb, turnBetaVIII
        };
      };
      std::string sequenceFile, outputFile;
      int helix, sheet, turn;
      fold() :
        helix(helices::alpha), sheet(sheets::betaAParallel),
        turn(turns::turnBetaIa)
      { }
    };

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

    std::ostream& operator<< (std::ostream &, solvadd const &);

    struct ringsearch
    {
      ::coords::float_type bias_force, chance_close;
      std::size_t population, generations;
      ringsearch(void)
        : bias_force(0.1), chance_close(0.75),
        population(12), generations(20)
      { }
    };

    std::ostream& operator<< (std::ostream &, ringsearch const &);

  }



  struct startopt
  {
    struct types { enum T { RINGSEARCH, SOLVADD, RINGSEARCH_SOLVADD }; };
    //startopt_conf::fold fold;
    startopt_conf::solvadd solvadd;
    startopt_conf::ringsearch ringsearch;
    types::T type;
    std::size_t number_of_structures;
    startopt(void)
      : /*fold(),*/ solvadd(), ringsearch(),
      type(types::SOLVADD),
      number_of_structures()
    { }
  };

  std::ostream& operator<< (std::ostream &, startopt const &);




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
    //double cdist_cutoff; <- CONTACT DISTANCE NOT YET IMPLEMENTED
    align(void) : dist_unit(0), reference_frame_num(0), traj_align_translational(true), traj_align_rotational(true), traj_print_bool(true), holm_sand_r0(20), align_external_file()//, cdist_cutoff(5) 
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
      pca_trunc_atoms_num(), pca_internal_dih(), pca_dimensions_for_histogramming(1u),
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

  void parse_command_switches(std::ptrdiff_t const, char**);
  std::string config_file_from_commandline(std::ptrdiff_t const, char**);
  void parse_option(std::string const option, std::string const value);
}

class Config
{
public:

  // Constructor, automatically parses file
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

  config::general               general;
  config::coords                coords;
  config::energy                energy;
  config::startopt              startopt;
  config::optimization          optimization;
  config::fep                   fep;
  config::molecular_dynamics    md;
  config::path                  path;
  //config::bias                  bias;
  config::dimer                 dimer;
  config::neb					          neb;
  config::generalized_born      gbsa;
  config::adjust                adjustment;
  config::align			            alignment;
  config::PCA					          PCA;
  config::entropy				        entropy;
  config::io                    io;


  static config::tasks::T            getTask(std::string const&);
  static config::interface_types::T  getInterface(std::string const&);
  static config::output_types::T     getOutFormat(std::string const&);
  static config::solvs::S            getSolv(std::string const&);
  static config::surfs::SA           getSurf(std::string const&);

private:

  void parse_file(std::string const &filename);
  static Config * m_instance;

};
