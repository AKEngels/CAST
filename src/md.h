/**
CAST 3
md.h
Purpose: header for molecular dynamics simulation

@version 1.0
*/

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
#include "Scon/scon_chrono.h"
#include "Scon/scon_vect.h"
#include "Scon/scon_serialization.h"
#include "Scon/scon_log.h"
#include "Scon/scon_utility.h"
#include "helperfunctions.h"
#include "md_analysis.h"
#include "md_logging.h"

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

    void setQ1(double Q1inp) { Q1 = Q1inp; }
    void setQ2(double Q2inp) { Q2 = Q2inp; }
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

  class CoordinatesUBIAS {
  public:
    CoordinatesUBIAS(coords::Coordinates* const coordinates) : coordinates(coordinates), broken_bonds() {}
    //umbrella
    void ubias(std::vector<double>& uout)
    {
      if (!coordinates->m_potentials.uempty())
        coordinates->m_potentials.umbrellaapply(coordinates->m_representation.structure.cartesian,
          coordinates->m_representation.gradient.cartesian,
          uout);
    }

    // very ugly approach to get the correct return type. No other way seen yet
    decltype(std::declval<coords::Coordinates>().atoms(std::declval<int>()))
      atoms(int const i) const {
      return coordinates->atoms(i);
    }
    decltype(std::declval<coords::Coordinates>().size())
      size() const {
      return coordinates->size();
    }
    decltype(std::declval<coords::Coordinates>().g())
      g() {
      return coordinates->g();
    }
    decltype(std::declval<coords::Coordinates>().xyz(std::declval<int>()))
      xyz(int const i) const {
      return coordinates->xyz(i);
    }
    decltype(std::declval<coords::Coordinates>().xyz())
      xyz() const {
      return coordinates->xyz();
    }
    bool validate_bonds();
    decltype(std::declval<coords::Coordinates>().getFep())
      getFep() const {
      return coordinates->getFep();
    }
    std::vector<std::vector<float>> const& getBrokenBonds() const { return broken_bonds; }
    auto center_of_mass() const {
      return coordinates->center_of_mass();
    }
    auto center_of_geometry() const {
      return coordinates->center_of_geometry();
    }
    template<typename ... Args>
    auto move_all_by(Args ... args) {
      return coordinates->move_all_by(args...);
    }
    template<typename ... Args>
    void move_atom_by(Args ... args) {
      coordinates->move_atom_by(args...);
    }
    template<typename ... Args>
    void move_atom_to(Args ... args) {
      coordinates->move_atom_to(args...);
    }
    /** add gradients to an atom (spherical boundaries)
    @param index: atom index
    @param g: gradients that should be added to the gradients of index*/
    void add_sp_gradients(std::size_t const index, coords::Cartesian_Point const& g)
    {
      coordinates->m_representation.gradient.cartesian[index] += g;
    }
    /**scale the coordinates of an atom (used for pressure control)
    @param index: index of atom that is to be moved
    @param p: factor by which the coordinates should be scaled
    @param force_move: if set to true also move fixed atoms*/
    void scale_atom_by(std::size_t const index, double& p, bool const force_move = false)
    {
      if (!atoms(index).fixed() || force_move)
      {
        coordinates->m_representation.structure.cartesian[index] *= p;
        coordinates->energy_valid = false;
      }
      coordinates->m_stereo.update(xyz());
    }
    /** set new atom coordinates if anisotropic pressure control is enabled*/
    void set_atom_aniso(std::size_t const index, double tvec[3][3], bool const force_move = false) {
      if (!atoms(index).fixed() || force_move)
      {
        double tcor;
        tcor = tvec[0][0] * coordinates->m_representation.structure.cartesian[index].x() +
          tvec[0][1] * coordinates->m_representation.structure.cartesian[index].y() + tvec[0][2] *
          coordinates->m_representation.structure.cartesian[index].z();
        coordinates->m_representation.structure.cartesian[index].x() = tcor;
        tcor = tvec[1][0] * coordinates->m_representation.structure.cartesian[index].x() +
          tvec[1][1] * coordinates->m_representation.structure.cartesian[index].y() + tvec[1][2] *
          coordinates->m_representation.structure.cartesian[index].z();
        coordinates->m_representation.structure.cartesian[index].y() = tcor;
        tcor = tvec[2][0] * coordinates->m_representation.structure.cartesian[index].x() +
          tvec[2][1] * coordinates->m_representation.structure.cartesian[index].y() + tvec[2][2] *
          coordinates->m_representation.structure.cartesian[index].z();
        coordinates->m_representation.structure.cartesian[index].z() = tcor;
        coordinates->energy_valid = false;
      }
      coordinates->m_stereo.update(xyz());
    }
    template<typename ... Args>
    void set_xyz(Args ... args) const {
      return coordinates->set_xyz(args...);
    }
    template<typename ... Args>
    auto g_xyz(Args ... args) const {
      return coordinates->g_xyz(args...);
    }

    /**updates the topology*/
    void energy_update(bool const skip_topology = false) { coordinates->energy_update(skip_topology); }
    /**returns matrix with interactions between subsystems*/
    coords::sub_ia_matrix_t& interactions()
    {
      return coordinates->interactions();
    }
    /**returns matrix with interactions between subsystems*/
    coords::sub_ia_matrix_t const& interactions() const
    {
      return coordinates->interactions();
    }
    /**adds V to virial coefficients
    @param V: values that should be added*/
    void add_to_virial(std::array<std::array<double, 3>, 3> & V)
    {
      coordinates->m_virial[0][0] += V[0][0];
      coordinates->m_virial[1][0] += V[1][0];
      coordinates->m_virial[2][0] += V[2][0];
      coordinates->m_virial[0][1] += V[0][1];
      coordinates->m_virial[1][1] += V[1][1];
      coordinates->m_virial[2][1] += V[2][1];
      coordinates->m_virial[0][2] += V[0][2];
      coordinates->m_virial[1][2] += V[1][2];
      coordinates->m_virial[2][2] += V[2][2];
    };
    /**returns the virial coefficients*/
    coords::virial_t const& virial() const { return coordinates->m_virial; }
    decltype(std::declval<coords::Coordinates>().pes())
      pes() const {
      return coordinates->pes();
    }
  protected:
    coords::Coordinates* const coordinates;

    /**vector of broken bonds (determined by validate_bonds())
    i.e. bondlength either too short or too long
    each element of the vector is a vector which contains the numbers of the two atoms that form the bond and the bond length*/
    std::vector<std::vector<float>> broken_bonds;
  };


  /** class for MD simulation
  */
  class simulation
  {

  private:

    /**coordinates object*/
    CoordinatesUBIAS coordobj;

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
    /**total kinetic energy*/
    double E_kin;
    /** desired temperature */
    double T;
    /** current temperature */
    double temp;
    /** Pressure stuff*/
    double press;
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
    std::vector<size_t> inner_atoms;   //
    /**atoms that move (distance smaller than outer cutoff)*/
    std::vector<size_t> movable_atoms; // 

      /** vector with lambda-values for every FEP window */
    std::vector<fepvar> window;
    /** Umbrella sampling vectors */
    std::vector<double> udatacontainer;

    /** save restarted status */
    bool restarted;

    /** initialization */
    void init(void);
    /** remove translational and rotational momentum of the whole system */
    void removeTranslationalAndRotationalMomentumOfWholeSystem(void);

    /**function for calculation of current target temperature
    @param step: current MD step
    @param fep: true if in equilibration of production of FEP run, then temperature is kept constant
    */
    bool heat(std::size_t const step, bool fep);
    /** nose hoover thermostat (velocity scaling is done automatically in this function)*/
    double nose_hoover_thermostat(void);
    /** nose hoover thermostat only for some atoms when used together with biased potential or fixed atoms
    returns the temperature scaling factor for velocities (scaling has to be performed after this function)
    @param atoms: vector with atom indizes of those atoms that are to be used to calculate scaling factor*/
    double nose_hoover_thermostat_some_atoms(std::vector<size_t> atoms);

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
    void check_rattlepair_for_bond(config::md_conf::config_rattle::rattle_constraint_bond& rctemp);
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

    /** Get new kinetic energy from current velocities of atoms
    @param atom_list: vector of atom numbers whose energy should be calculated*/
    void updateEkin(std::vector<size_t> atom_list);

    


    /** Berendsen pressure coupling (doesn't work */
    void berendsen(double const);

    /**write a restartfile
    @param k: current MD step*/
    void write_restartfile(std::size_t const k);



  public:

    /** constructor
    @param coords::Coordinates &: coords-object whose movement should be simulated
    */
    simulation(coords::Coordinates&);
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

    // ANALYZING STUFF 

    /**vector of atom pairs that are to be analyzed*/
    std::vector<md_analysis::ana_pair> ana_pairs;
    /**vector of zones to be plotted*/
    std::vector <md_analysis::zone> zones;

    /**calculate and get distances from active center*/
    std::vector<double> calc_distances_from_center() {
      return init_active_center(0);
    }

    coords::float_type md::simulation::getEkin(std::vector<std::size_t> atom_list) const;

    /**function to retrieve reference to coordinates object*/
    CoordinatesUBIAS& get_coords() { return coordobj; }

    // OPERATORS

    //**overload for << operator*/
    template<class Strm>
    friend scon::binary_stream<Strm>& operator<< (scon::binary_stream<Strm>& strm, simulation const& sim)
    {
      std::array<std::size_t, 9u> const sizes = {
        sim.P.size(), sim.P_old.size(), sim.F.size(),
        sim.F_old.size(), sim.V.size(), sim.M.size(),
        sim.rattle_bonds.size(), sim.window.size(), sim.udatacontainer.size() };
      // sizes
      for (auto const& s : sizes)
        strm << s;

      // logger
      // Disabled 09.08.17 by Dustin Kaiser
      // as there is some bug in the I/O of the logger
      // data from binary streams
      //strm << sim.logging;


      // Non-Fundamental vectors
      for (auto const& x : sim.coordobj.xyz()) strm << x;
      for (auto const& p : sim.P) strm << p;
      for (auto const& p : sim.P_old) strm << p;
      for (auto const& f : sim.F) strm << f;
      for (auto const& f : sim.F_old) strm << f;
      for (auto const& v : sim.V) strm << v;
      for (auto const& m : sim.M) strm << m;
      // rest
      strm << sim.M_total << sim.E_kin_tensor <<
        sim.E_kin << sim.T << sim.temp << sim.press <<
        sim.dt << sim.freedom << sim.snapGap << sim.C_geo << sim.C_mass <<
        sim.nht << sim.rattle_bonds << sim.window << sim.udatacontainer;
      return strm;
    }

    //**overload for >> operator*/
    template<class Strm>
    friend scon::binary_stream<Strm>& operator>> (scon::binary_stream<Strm>& strm, simulation& sim)
    {
      std::array<std::size_t, 9u> sizes;
      // sizes
      for (auto& s : sizes)
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
      for (auto& p : sim.P) strm >> p;
      for (auto& p : sim.P_old) strm >> p;
      for (auto& f : sim.F) strm >> f;
      for (auto& f : sim.F_old) strm >> f;
      for (auto& v : sim.V) strm >> v;
      for (auto& m : sim.M) strm >> m;
      // rest
      strm >> sim.M_total >> sim.E_kin_tensor >>
        sim.E_kin >> sim.T >> sim.temp >> sim.press >>
        sim.dt >> sim.freedom >> sim.snapGap >> sim.C_geo >> sim.C_mass >>
        sim.nht >> sim.rattle_bonds >> sim.window >> sim.udatacontainer;
      return strm;
    }


  };

}

