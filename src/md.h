#pragma once 

#include <vector>
#include <string>
#include <iostream>
#include <memory>

#include "configuration.h"
#include "coords.h"
#include "scon_vect.h"
#include "scon_serialization.h"
#include "scon_log.h"

/*
Q an J:
- refine nach e() ?
- Quelle Nosé-Hoover
- Inertia Tensor Fehler in Init?
*/

//class coordinates;

namespace md
{

  static const double kB = 0.83144725;
  static const double convert = 418.4; //kcal to g*A^2/ps^2
  static const double negconvert = -418.4; //kcal to g*A^2/ps^2
  static const double hnconvert = -209.2; //kcal to g*A^2/ps^2
  static const double PI = 3.14159265358979323;
  static const double R = 1.9872066e-3;
  static const double presc = 6.85684112e4;

  struct trace_data
  {
    std::vector<coords::float_type> Eia;
    coords::float_type T, Ek, Ep, P;
    std::size_t i, snapshot;
    trace_data() : Eia(), T(), Ek(), Ep(), P(), i(), snapshot() {}
    trace_data(std::vector<coords::float_type> const &E_ia,
      coords::float_type temp, coords::float_type E_kin, 
      coords::float_type E_pot, coords::float_type press,
      std::size_t iteration, std::size_t snap_number) :
      Eia(E_ia), T(temp), Ek(E_kin), Ep(E_pot), P(press),
      i(iteration), snapshot(snap_number)
    { }
  };

  std::ostream& operator<< (std::ostream &, trace_data const &);

  template<class Strm>
  scon::binary_stream<Strm> & operator<< (scon::binary_stream<Strm> &str, trace_data const &t)
  {
    str << t.T << t.Ek << t.Ep << t.P << t.i << t.snapshot;
    str << t.Eia.size() << t.Eia;
    return str;
  }

  template<class Strm>
  scon::binary_stream<Strm> & operator>> (scon::binary_stream<Strm> &str, trace_data &t)
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

  class trace_writer
  {
    std::unique_ptr<std::ofstream> strm;
  public:
    trace_writer() : strm() {}
    trace_writer(char const * const filename)
      : strm(new std::ofstream(filename, std::ios::out))
    {}
    void operator() (trace_data  const & xyz);
  };

  class Logger
  {

    coords::offset_buffered_cartesian_logfile snap_buffer;
    scon::vector_offset_buffered_callable<trace_data, trace_writer> data_buffer;
    std::size_t snapnum;

  public:

    Logger(coords::Coordinates &coords, std::size_t snap_offset);

    bool operator() (std::size_t const iter,
      coords::float_type const T, 
      coords::float_type const P,
      coords::float_type const Ek, 
      coords::float_type const Ep,
      std::vector<coords::float_type> const Eia,
      coords::Representation_3D const & x);

    template<class Strm>
    friend scon::binary_stream<Strm> & operator<< (scon::binary_stream<Strm> &str, Logger const &l)
    {
      str << l.snapnum;
      str << l.snap_buffer;
      str << l.data_buffer;
      return str;
    }

    template<class Strm>
    friend scon::binary_stream<Strm> & operator>> (scon::binary_stream<Strm> &str, Logger &l)
    {
      str >> l.snapnum;
      str >> l.snap_buffer;
      str >> l.data_buffer;
      return str;
    }

  };

  struct nose_hoover
  {
    double v1, v2, x1, x2;
    double Q1, Q2, G1, G2;
    nose_hoover(void) :
      v1(0.0), v2(0.0), x1(0.0), x2(0.0),
      Q1(0.1), Q2(0.1), G1(0.0), G2(0.0)
    { }
    //std::size_t bytesize (void) const { return sizeof(double)* 8; }
    //void copy_to_buffer(std::vector<char>&) const;
    //void deserialize_from_stream(std::istream &);
  };

  struct fepvar
  {
    double ein, eout, vin, vout;
    double dein, deout, dvin, dvout;
  };

  namespace thermostat
  {

    class velocity_rescaling
    {

      double temp_factor;

    public:

      velocity_rescaling(std::size_t const degrees_of_freedom)
        : temp_factor(2. / (R*degrees_of_freedom)) {}

      template<class V>
      double operator() (V & velocity_range, 
        double E_kin, double const target_temperature)
      {
        using std::sqrt;
        auto const is_temperature = E_kin * temp_factor;
        auto const F = sqrt(target_temperature / is_temperature);
        for (auto & v : velocity_range) v *= F;
        return E_kin * F * F;
      }

    };

    struct nose_hoover
    {
      double v1, v2, x1, x2, Q1, Q2, G1, G2;
      double d2, d4, d8;
      std::size_t freedom;

      nose_hoover(double const timestep,
        std::size_t const degrees_of_freedom) :
        v1(), v2(), x1(), x2(), Q1(0.1), Q2(.1), G1(), G2(),
        d2(timestep / double(2.)), d4(d2 / double(2.)),
        d8(d4 / double(2.)), freedom(degrees_of_freedom)
      { }

      template<class V>
      double operator() (V & velocity_range, double E_kin,
        double const target_temperature)
      {
        using std::exp;
        auto const TR = target_temperature*md::R;
        auto const fTR = TR*freedom;
        G1 = (2.0*E_kin - fTR) / Q1;
        G2 = (Q1*v1*v1 - TR) / Q2;
        v2 += G2*d4;
        v1 *= exp(-v2*d8);
        v1 += G1*d4;
        v1 *= exp(-v2*d8);
        x1 += v1*d2;
        x2 += v2*d2;
        auto tempscale = exp(-v1*d2);
        for (auto & v : velocity_range) v *= tempscale;
        E_kin *= tempscale*tempscale;
        v1 *= exp(-v2*d8);
        G1 = (2.0*E_kin - fTR) / Q1;
        v1 += G1*d4;
        v1 *= exp(-v2*d8);
        G2 = (Q1*v1*v1 - TR) / Q2;
        v2 += G2*d4;
        return E_kin;
      }
    };

  }

  namespace barostat
  {

    struct berendsen
    {

      double target, delay, compress;
      bool isotropic;

      double operator() (double const time, coords::Representation_3D & p,
        coords::Tensor const & Ek_T, coords::Tensor const & Vir_T, 
        coords::Cartesian_Point & box);

    };
  }
  
  namespace integrator
  {

    template<class Thermostat_T, class Barostat_T>
    class veloctiy_verlet
    {

      // Thermostat
      Thermostat_T thermo;

      // Barostat
      Barostat_T baro;

      // Logging
      Logger log;

      // Data
      coords::Tensor E_kin_tensor;
      coords::Tensor Virial_tensor;
      coords::Representation_3D x, g;

      std::size_t freedom;




    };

    using default_vv = veloctiy_verlet<thermostat::nose_hoover, barostat::berendsen>;


  }

  class simulation
  {

  private:

    simulation& operator= (simulation const &);

    // pointer to coordinates
    coords::Coordinates & coordobj;

    // Logger
    Logger logging;

    // positions, forces, velocities
    coords::Representation_3D P, P_old, F, F_old, V;
    // masses
    std::vector<double> M;
    // total Mass of the System
    double M_total;
    //! Kinetic Energy
    //std::array<std::array<double, 3>, 3> E_kin_tensor;
    coords::Tensor E_kin_tensor;
    double E_kin;
    //! Temperature
    double T, temp;
    //! Pressure
    double press, presstemp;
    //! timestep
    double dt;
    // degrees of freedom, snapshot offset
    std::size_t freedom, snapGap;
    // geometric center and center of mass
    coords::Cartesian_Point C_geo, C_mass;
    // trace temperature and energy evolvation
    //md::trace traces;
    //md::traces trace;
    // nose hoover thermostat values
    md::nose_hoover nht;
    // rattle constraints
    std::vector<config::md_conf::config_rattle::rattle_constraint_bond> rattle_bonds;

    //! Fep progress vector
    std::vector<fepvar> window;
    //! Umbrella sampling vectors
    std::vector<double> udatacontainer;

    // save restarted status
    bool restarted;

    // initialization
    void init(void);
    // remove translational and rotational momentum of the whole system
    void tune_momentum(void);

    //! heating
    bool heat(std::size_t const step);
    // nose hoover thermostat
    void nose_hoover_thermostat(void);

    //! rattle feature
    void rattle_pre(void);
    void rattle_post(void);

    // integrator selector
    void integrate(std::size_t const k_init = 0U);

    //! beeman integrator
    void beemanintegrator(std::size_t const k_init = 0U);
    //! velocity_verlet
    void velocity_verlet(std::size_t const k_init = 0u);

    //! adjusting boundary conditions
    void boundary_adjustments(void);
    void spherical_adjust(void);

    //! Kinetic Energy update
    void updateEkin(void);

    //! Berendsen pressure coupling
    void berendsen(double const &);

    //std::size_t data_bytesize(void) const;
    //void serialize_to_stream(std::ostream &, std::size_t const & current_step) const;
    //std::size_t deserialize_from_stream(std::istream &);

    //// read restart file and return iteration
    //std::size_t read_restartfile(void);

    void write_restartfile(std::size_t const);

  public:

    // construct
    simulation(coords::Coordinates &);
    // start
    void run(bool const restart = true);
    //umbrella sampling
    void umbrella_run(bool const restart = true);
    //! FEP calculation
    void rattlesetup(void);
    void fepinit(void);
    void feprun();
    void freecalc();
    void freewrite(std::size_t);
    bool prod;
    double FEPsum;
    void print_init_info(void);

    template<class Strm>
    friend scon::binary_stream<Strm> &
      operator<< (scon::binary_stream<Strm> &strm, simulation const &sim)
    {
      std::array<std::size_t, 9u> const sizes = { {sim.P.size(), sim.P_old.size(),
        sim.F.size(), sim.F_old.size(), sim.V.size(), sim.M.size(),
        sim.rattle_bonds.size(), sim.window.size(), sim.udatacontainer.size()} };
      // sizes
      for (auto const & s : sizes) strm << s;
      // logger
      strm << sim.logging;
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

    template<class Strm>
    friend scon::binary_stream<Strm> &
      operator>> (scon::binary_stream<Strm> &strm, simulation &sim)
    {
      std::array<std::size_t, 10u> sizes;
      // sizes
      for (auto & s : sizes) strm >> s;
      // logger
      strm >> sim.logging;
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

