#pragma once 

#ifdef USE_PYTHON
#include <Python.h>
#endif
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <atomic>
#include <string>
#include <memory>



#include "configuration.h"
#include "coords.h"
#include "coords_io.h"
#include "scon_chrono.h"
#include "scon_vect.h"
#include "scon_serialization.h"
#include "scon_log.h"
#include "scon_utility.h"
#include "helperfunctions.h"

/**
*namespace for everything that has to do with molecular dynamics simulatinons
*/
namespace md
{
  /**boltzmann constant*/
  static const double kB = 0.83144725;
  /**conversion factor kcal to g*A^2/ps^2*/
  static const double convert = 418.4;
  /**negativ conversion factor kcal to g*A^2/ps^2*/
  static const double negconvert = -convert;
  /**pi*/
  static const double PI = 3.14159265358979323;
  /**gas constant*/
  static const double R = 1.9872066e-3;
  /**conversion factor: (kcal/mol) / (atm*A^3) */
  static const double presc = 6.85684112e4;   // 1.0/6.85684112e4 would make more sense in my opinion 
                                              // but the program wouldn't always give 0.00000 for pressure

  /**
  *collection of current simulation data
  */
  struct trace_data
  {
    /**energy of subsystems???*/
    std::vector<coords::float_type> Eia;
    /**temperature*/
    coords::float_type T;
    /**kinetic energy*/
    coords::float_type Ek;
    /**potential energy*/
    coords::float_type Ep;
    /**pressure*/
    coords::float_type P;
    /**step-number*/
    std::size_t i;
    /**snapshot-number*/
    std::size_t snapshot;
    /**
    * default constructor
    */
    trace_data() : Eia(), T(), Ek(), Ep(), P(), i(), snapshot() {}
    /**
    * another constructor that already takes all the public members as parameters
    */
    trace_data(std::vector<coords::float_type> const &E_ia,
      coords::float_type temp, coords::float_type E_kin,
      coords::float_type E_pot, coords::float_type press,
      std::size_t iteration, std::size_t snap_number) :
      Eia(E_ia), T(temp), Ek(E_kin), Ep(E_pot), P(press),
      i(iteration), snapshot(snap_number)
    { }
  };

  /**
  writes content of trace_data object d into a stream
  */
  std::ostream& operator<< (std::ostream &, trace_data const &);

  /**
  overload of << operator
  */
  template<class Strm>
  scon::binary_stream<Strm> & operator<< (scon::binary_stream<Strm> &str, trace_data const &t)
  {
    str << t.T << t.Ek << t.Ep << t.P << t.i << t.snapshot;
    str << t.Eia.size() << t.Eia;
    return str;
  }

  /**
  overload of >> operator
  */
  template<class Strm>
  scon::binary_stream<Strm> & operator >> (scon::binary_stream<Strm> &str, trace_data &t)
  {
    decltype(t.Eia.size()) x = 0;
    if (str >> t.T && str >> t.Ek && str >> t.Ep &&
      str >> t.P && str >> t.i && str >> t.snapshot && str >> x)
    {
      t.Eia.resize(x);
      str >> t.Eia;
    }
    return str;
  }

  /**
  class for writing trace_data into a file
  */
  class trace_writer
  {
    std::unique_ptr<std::ofstream> strm;
  public:
    /**default constructor*/
    trace_writer() : strm() {}
    /**another constructor
    * @param filename: name of the file where the information should be written
    */
    trace_writer(char const * const filename)
      : strm(new std::ofstream(filename, std::ios::out))
    {}
    /** function for writing the data
    * @param xyz: trace_data object where information should be taken
    */
    void operator() (trace_data  const & xyz);
  };

  /**
  class for collecting logging information (trace data and snapshots)
  */
  class Logger
  {

    coords::offset_buffered_cartesian_logfile snap_buffer;
    scon::vector_offset_buffered_callable<trace_data, trace_writer> data_buffer;
    std::size_t snapnum;

  public:

    /**writes snapshots
    @param coords: coords-object
    @param snap_offset: has something to do with MDsnapbuffer???
    */
    Logger(coords::Coordinates &coords, std::size_t snap_offset);

    /**looks every 5000 steps if temperature, pressure or energy is nan and throws an error if yes
    */
    bool operator() (std::size_t const iter,
      coords::float_type const T,
      coords::float_type const P,
      coords::float_type const Ek,
      coords::float_type const Ep,
      std::vector<coords::float_type> const Eia,
      coords::Representation_3D const & x);

    /**
    overload of << operator
    */
    template<class Strm>
    friend scon::binary_stream<Strm> & operator<< (scon::binary_stream<Strm> &str, Logger const &l)
    {
      str << l.snapnum;
      str << l.snap_buffer;
      str << l.data_buffer;
      return str;
    }

    /**
    overload of >> operator
    */
    template<class Strm>
    friend scon::binary_stream<Strm> & operator >> (scon::binary_stream<Strm> &str, Logger &l)
    {
      str >> l.snapnum;
      str >> l.snap_buffer;
      //str >> l.data_buffer;
      return str;
    }

  };


  /** Nose-Hover thermostat. Variable names and implementation are identical to the book of
  Frenkel and Smit, Understanding Molecular Simulation, Appendix E */
  struct nose_hoover
  {
    double v1, v2, x1, x2;
    double Q1, Q2, G1, G2;
    nose_hoover(void) :
      v1(0.0), v2(0.0), x1(0.0), x2(0.0),
      Q1(0.1), Q2(0.1), G1(0.0), G2(0.0)
    { }
  };

  /** collection of variables for FEP calculation
  */
  struct fepvar
  {
    /**lambda_el of former window for appearing atoms*/
    double mein;
    /**lambda_el of former window for disappearing atoms*/
    double meout;
    /**lambda_vdw of former window for appearing atoms*/
    double mvin;
    /**lambda_vdw of former window for disappearing atoms*/
    double mvout;
    /**lambda_el for appearing atoms*/
    double ein;
    /**lambda_el for disappearing atoms*/
    double eout;
    /**lambda_vdw for appearing atoms*/
    double vin;
    /**lambda_vdw for disappearing atoms*/
    double vout;
    /**lambda_el of next window for appearing atoms*/
    double dein;
    /**lambda_el of next window for disappearing atoms*/
    double deout;
    /**lambda_vdw of next window for appearing atoms*/
    double dvin;
    /**lambda_vdw of next window for disappearing atoms*/
    double dvout;
  };

  /**struct that contains information about an atom pair that is to be analyzed*/
  struct ana_pair
  {
    /**index of first atom (starting with 0)*/
    int a;
    /**index of second atom (starting with 0)*/
    int b;
    /**element symbol of first atom*/
    std::string symbol_a;
    /**element symbol of second atom*/
    std::string symbol_b;
    /**name of first atom (element and tinker atom index)*/
    std::string name_a;
    /**name of second atom (element and tinker atom index)*/
    std::string name_b;
    /**legend of this atom pair in the graph*/
    std::string legend;
    /**distances for every MD frame*/
    std::vector<double> dists;

    /**constructor
    @param p1: tinker atom index of first atom (i.e. starting with 1)
    @param p2: tinker atom index of second atom (i.e. starting with 1)*/
    ana_pair(int p1, int p2) { a = p1 - 1; b = p2 - 1; }

    /**returns a string with information about the atom pair*/
    std::string info()
    {
      std::string result = "Atoms: " + name_a + " , " + name_b + "\n";
      for (auto d : dists)
      {
        result += std::to_string(d) + " , ";
      }
      return result + "\n";
    }
  };

  /**information about a zone for which temperature is to by analyzed*/
  struct zone
  {
    /**legend for plotting*/
    std::string legend;
    /**atom indizes (starting with 0)*/
    std::vector<int> atoms;
    /**temperatures for every MD step*/
    std::vector<double> temperatures;
  };


  /** class for MD simulation
  */
  class simulation
  {

  private:

    /**pointer to coordinates*/
    coords::Coordinates & coordobj;

    /** Logger */
    Logger logging;

    // positions, forces, velocities
  /**positions*/
    coords::Representation_3D P;
    /**old positions*/
    coords::Representation_3D P_old;
    /**positions at the beginning of the simulation*/
    coords::Representation_3D P_start;
    /**forces*/
    coords::Representation_3D F;
    /**old forces*/
    coords::Representation_3D F_old;
    /**velocities*/
    coords::Representation_3D V;
    /** masses */
    std::vector<double> M;
    /** total Mass of the System */
    double M_total;
    /**tensor for kinetic energy*/
    coords::Tensor E_kin_tensor;
    /**tensor for kinetic energy*/
    double E_kin;
    /** desired temperature */
    double T;
    /** current temperature */
    double temp;
    /** Pressure stuff*/
    double press, presstemp;
    std::vector <double> press_box;
    /** timestep */
    double dt;
    /** degrees of freedom */
    std::size_t freedom;
    /**snapshot offset (gap between two snapshots) */
    std::size_t snapGap;
    /** geometric center */
    coords::Cartesian_Point C_geo;
    /**center of mass*/
    coords::Cartesian_Point C_mass;
    /** nose hoover thermostat values */
    md::nose_hoover nht;
    /** rattle constraints */
    std::vector<config::md_conf::config_rattle::rattle_constraint_bond> rattle_bonds;

    // stuff for biased potential
    /**distances to active site for every atom*/
    std::vector<double> distances;  // distances to active site for every atom
    /**atoms with a distance smaller than the inner cutoff*/
    std::vector<int> inner_atoms;   //
    /**atoms that move (distance smaller than outer cutoff)*/
    std::vector<int> movable_atoms; // 

      /** vector with lambda-values for every FEP window */
    std::vector<fepvar> window;
    /** Umbrella sampling vectors */
    std::vector<double> udatacontainer;

    /** save restarted status */
    bool restarted;

    /** initialization */
    void init(void);
    /** remove translational and rotational momentum of the whole system */
    void tune_momentum(void);

    /**function for calculation of current target temperature
    @param step: current MD step
    @param fep: true if in equilibration of production of FEP run, then temperature is kept constant
    */
    bool heat(std::size_t const step, bool fep);
    /** nose hoover thermostat (velocity scaling is done automatically in this function)*/
    void nose_hoover_thermostat(void);
    /** nose hoover thermostat only for some atoms when used together with biased potential or fixed atoms
    returns the temperature scaling factor for velocities (scaling has to be performed after this function)
    @param atoms: vector with atom indizes of those atoms that are to be used to calculate scaling factor*/
    double nose_hoover_thermostat_some_atoms(std::vector<int> atoms);

    /**sets coordinates to original values and assigns random velocities*/
    void restart_broken();

    /** function to control the temperature
    @param thermostat: determines if nose-hoover-thermostat or direct velocity scaling
    @param half: determines if half or fullstep (only relevant for console output)
    */
    double tempcontrol(bool thermostat, bool half);

    /** calculate distances to active center
    @param counter: current MD step
    */
    std::vector<double> init_active_center(int counter);
    /** scale down velocities according to distance to active site
    @param atom_number: atom number (starting with zero)
    @param inner_cutoff: radius of inner cutoff
    @param outer_cutoff: radius of outer cutoff
    */
    coords::Cartesian_Point adjust_velocities(int atom_number, double inner_cutoff, double outer_cutoff);

		/**function that checks if the two atoms of a rattlepair are bonded with each other,
		throws an error if not
		@param rctemp: rattlepair*/
		void check_rattlepair_for_bond(config::md_conf::config_rattle::rattle_constraint_bond &rctemp);
    /** rattle feature pre */
    void rattle_pre(void);
    /** rattle feature post */
    void rattle_post(void);

    /**select an integrator (velocity-verlet or beeman)
  @param fep: true if in equilibration of production of FEP run, then temperature is kept constant
  @param k_init: step where the MD starts (zero should be okay)
  */
    void integrate(bool fep = false, std::size_t const k_init = 0U);

    /**velocity-verlet or beeman integrator
    @param fep: true if in equilibration of production of FEP run, then temperature is kept constant
    @param k_init: step where the MD starts (zero should be okay)
    @param beeman: true if beeman integrator is used, false if velocity verlet integrator is used
    */
    void integrator(bool fep, std::size_t const k_init = 0U, bool beeman = false);

    /** tell user that he applies spherical boundary conditions
  */
    void boundary_adjustments(void);
    /** adjusting spherical boundary conditions (yet to be tested)
    */
    void spherical_adjust(void);

    /** Kinetic Energy update */
    void updateEkin(void);
    /** Kinetic Energy update for a selection of atoms
    @param atom_list: vector of atom numbers whose energy should be calculated
    */
    void updateEkin_some_atoms(std::vector<int> atom_list);

    /** Berendsen pressure coupling (doesn't work */
    void berendsen(double const);

    /**write a restartfile
    @param k: current MD step
    */
    void write_restartfile(std::size_t const k);

    /**vector with average temperature for every frame*/
    std::vector<double> temperatures;

    /**vector of atom pairs that are to be analyzed*/
    std::vector<ana_pair> ana_pairs;

    /**vector of zones to be plotted*/
    std::vector < zone> zones;

  public:

    /** constructor
  @param coords::Coordinates &: coords-object whose movement should be simulated
  */
    simulation(coords::Coordinates &);
    /** start simulation
  @param restart: set to true (default) if simulation should start from the beginning
  */
    void run(bool const restart = true);
    /**prints information about MD simulation before the run starts*/
    void print_init_info(void);

    /** set up constraints for H-X bonds if requested
    ideal bond lengths are taken from the foce field parameter file
    specified by RATTpar in the INPUTFILE */
    void rattlesetup(void);

    /** perform an umbrella sampling
    @param restart: set to true (default) if simulation should start from the beginning
    */
    void umbrella_run(bool const restart = true);

    /** If FEP calculation is requested: calculate lambda values for each window
    and print the scaling factors for van-der-Waals and electrostatics for each window */
    void fepinit(void);
    /** perform FEP calculation if requested */
    void feprun();
    /**Calculation of ensemble average and free energy change after every window if FEP calculation is performed
     calculation can be improved if at every step the current averages are stored
     currently calculation is performed at the end of each window */
    void freecalc();
    /**calculation of free energy from Bennets acceptance ratio
    @param window: current window*/
    void bar(int window);
    /** write the output FEP calculations into "alchemical.txt" and "FEP_Results.txt"*/
    void freewrite(int);
    /**function that returns a string 
    this string can be run as a pythonprogramme that adds all paths necessary for FEP analysis to pythonpath*/
    std::string get_pythonpath();
    /**function that performs an FEP analysis (histogram and overlap of probability distributions)
    @param dE_pots: vector of dE_pot values for a conformation. 
    this is used to transfer these values to the next window as the dE value of window i corresponds to the dE_back value of window i+1
    @param window: number of current window
    returns vector with dE_pot values (explanation see above)*/
    std::vector<double> fepanalyze(std::vector<double> dE_pots, int window);
    /**function to plot temperature
    @param temps: temperatures to be plotted*/
    void plot_temp(std::vector<double> temps);
    /**function to plot distances for atom pairs
    @param pairs: atom pairs to be plotted*/
    void plot_distances(std::vector<ana_pair> &pairs);
    /**function to plot temperatures for all zones*/
    void plot_zones();
    /**function that fills zones with atoms*/
    std::vector<zone> find_zones();
    /**bool that determines if the current run is a production run or an equilibration run*/
    bool prod;
    /**current free energy difference for forward transformation*/
    double FEPsum;
    /**current free energy difference for backwards transformation*/
    double FEPsum_back;
    /**current free energy difference for simple overlap sampling (SOS)*/
    double FEPsum_SOS;
    /**free energy change of current window from SOS (start value for BAR)*/
    double dG_SOS;
    /**free energy difference for bennets acceptance ratio (BAR)*/
    double FEPsum_BAR;
    /**<exp^(-1/kT)*dE/2> save for use after next window (for SOS)*/
    double de_ensemble_v_SOS;
    /**<w*exp^(-1/kT)*dE/2> save for use after next window (for BAR)*/
    double de_ensemble_v_BAR;

    //**overload for << operator*/
    template<class Strm>
    friend scon::binary_stream<Strm>& operator<< (scon::binary_stream<Strm> &strm, simulation const &sim)
    {
      std::array<std::size_t, 9u> const sizes = {
        sim.P.size(), sim.P_old.size(), sim.F.size(), 
        sim.F_old.size(), sim.V.size(), sim.M.size(),
        sim.rattle_bonds.size(), sim.window.size(), sim.udatacontainer.size()};
      // sizes
      for (auto const & s : sizes) 
        strm << s;

      // logger
      // Disabled 09.08.17 by Dustin Kaiser
      // as there is some bug in the I/O of the logger
      // data from binary streams
      //strm << sim.logging;


      // Non-Fundamental vectors
      for (auto const & x : sim.coordobj.xyz()) strm << x;
      for (auto const & p : sim.P) strm << p;
      for (auto const & p : sim.P_old) strm << p;
      for (auto const & f : sim.F) strm << f;
      for (auto const & f : sim.F_old) strm << f;
      for (auto const & v : sim.V) strm << v;
      for (auto const & m : sim.M) strm << m;
      // rest
      strm << sim.M_total << sim.E_kin_tensor <<
        sim.E_kin << sim.T << sim.temp << sim.press << sim.presstemp <<
        sim.dt << sim.freedom << sim.snapGap << sim.C_geo << sim.C_mass <<
        sim.nht << sim.rattle_bonds << sim.window << sim.udatacontainer;
      return strm;
    }

    //**another overload for >> operator*/
    template<class Strm>
    friend scon::binary_stream<Strm>& operator>> (scon::binary_stream<Strm> &strm, simulation &sim)
    {
      std::array<std::size_t, 9u> sizes;
      // sizes
      for (auto & s : sizes) 
        strm >> s;

      // logger
      // Disabled 09.08.17 by Dustin Kaiser
      // as there is some bug in the I/O of the logger
      // data from binary streams
      //strm >> sim.logging;
      //
      // non-fundamental vectors
      sim.P.resize(sizes[0]);
      sim.P_old.resize(sizes[1]);
      sim.F.resize(sizes[2]);
      sim.F_old.resize(sizes[3]);
      sim.V.resize(sizes[4]);
      sim.M.resize(sizes[5]);
      sim.rattle_bonds.resize(sizes[6]);
      sim.window.resize(sizes[7]);
      sim.udatacontainer.resize(sizes[8]);
      for (auto i : scon::index_range(sim.coordobj.xyz()))
      {
        coords::Cartesian_Point p_new;
        strm >> p_new;
        sim.coordobj.move_atom_to(i, p_new, true);
      }
      for (auto & p : sim.P) strm >> p;
      for (auto & p : sim.P_old) strm >> p;
      for (auto & f : sim.F) strm >> f;
      for (auto & f : sim.F_old) strm >> f;
      for (auto & v : sim.V) strm >> v;
      for (auto & m : sim.M) strm >> m;
      // rest
      strm >> sim.M_total >> sim.E_kin_tensor >>
        sim.E_kin >> sim.T >> sim.temp >> sim.press >> sim.presstemp >>
        sim.dt >> sim.freedom >> sim.snapGap >> sim.C_geo >> sim.C_mass >>
        sim.nht >> sim.rattle_bonds >> sim.window >> sim.udatacontainer;
      return strm;
    }


  };

}

