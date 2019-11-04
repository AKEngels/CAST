/**
CAST 3
md_logging.h
Purpose: header for molecular dynamics simulation logging

@version 1.0
*/

#pragma once 

#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>


#include "coords.h"
#include "coords_io.h"
#include "Scon/scon_vect.h"
#include "Scon/scon_serialization.h"
#include "Scon/scon_log.h"
#include "Scon/scon_utility.h"

/**
*namespace for everything that has to do with molecular dynamics simulatinons
*/
namespace md
{

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
    trace_data(std::vector<coords::float_type> const& E_ia,
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
  std::ostream& operator<< (std::ostream&, trace_data const&);

  /**
  overload of << operator
  */
  template<class Strm>
  scon::binary_stream<Strm>& operator<< (scon::binary_stream<Strm>& str, trace_data const& t)
  {
    str << t.T << t.Ek << t.Ep << t.P << t.i << t.snapshot;
    str << t.Eia.size() << t.Eia;
    return str;
  }

  /**
  overload of >> operator
  */
  template<class Strm>
  scon::binary_stream<Strm>& operator >> (scon::binary_stream<Strm>& str, trace_data& t)
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
    trace_writer(char const* const filename)
      : strm(new std::ofstream(filename, std::ios::out))
    {}
    /** function for writing the data
    * @param xyz: trace_data object where information should be taken
    */
    void operator() (trace_data  const& xyz);
  };

  /**class for writing velocities into a file (like tinkerstructure but with velocities instead of coordinates)*/
  class velocity_logfile_drain
  {
    /**pointer to coordinates object*/
    coords::Coordinates* cp;
    /**unique pointer to std::ofstream object*/
    std::unique_ptr<std::ofstream> strm;

  public:
    /**default constructor*/
    velocity_logfile_drain() : cp(), strm() {}
    /**another constructor*/
    velocity_logfile_drain(coords::Coordinates& c, char const* const filename) :
      cp(&c), strm(std::make_unique<std::ofstream>(std::ofstream(filename, std::ios::out))) { }

    /**callback operator: writes structure (with given velocities velos) into ofstream*/
    void operator() (coords::Representation_3D&& velos)
    {
      if (!cp || !strm) return;
      // Save current state
      auto const tmp = (*cp).xyz();
      // plug velocities into coords (as xyz because then it can be better printed)
      (*cp).set_xyz(std::move(velos));
      // Print to stream
      *strm << *cp;
      // reset state
      (*cp).set_xyz(std::move(tmp));
    };
  };

  using offset_buffered_velocity_logfile =
    scon::vector_offset_buffered_callable<coords::Representation_3D, velocity_logfile_drain>;

  /**I guess this is a function to write the velostructures only after buffer is full*/
  inline offset_buffered_velocity_logfile make_buffered_velocity_log(coords::Coordinates& c, std::string file_suffix, std::size_t buffer_size,
    std::size_t log_offset)
  {
    return scon::offset_call_buffer<coords::Representation_3D, velocity_logfile_drain>(buffer_size, log_offset,
      velocity_logfile_drain{ c, coords::output::filename(file_suffix).c_str() });
  }

  /**
  class for collecting logging information (trace data, snapshots and velocities)
  */
  class Logger
  {
    /**object for writing snapshots*/
    coords::offset_buffered_cartesian_logfile snap_buffer;
    /**object for writing velocities*/
    offset_buffered_velocity_logfile velo_buffer;
    /**object for writing trace data*/
    scon::vector_offset_buffered_callable<trace_data, trace_writer> data_buffer;
    /**current snapshot number*/
    std::size_t snapnum;

  public:

    /**constructor
    @param coords: coords-object
    @param snap_offset: number of snapshots
    */
    Logger(coords::Coordinates& coords, std::size_t snap_offset);

    /**looks every 5000 steps if temperature, pressure or energy is nan and throws an error if yes
    if not writes tracefile, snapshots and velocities
    */
    bool operator() (std::size_t const iter,
      coords::float_type const T,
      coords::float_type const P,
      coords::float_type const Ek,
      coords::float_type const Ep,
      std::vector<coords::float_type> const Eia,
      coords::Representation_3D const& x,
      coords::Representation_3D const& v);

    /**
    overload of << operator
    */
    template<class Strm>
    friend scon::binary_stream<Strm>& operator<< (scon::binary_stream<Strm>& str, Logger const& l)
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
    friend scon::binary_stream<Strm>& operator >> (scon::binary_stream<Strm>& str, Logger& l)
    {
      str >> l.snapnum;
      str >> l.snap_buffer;
      //str >> l.data_buffer;
      return str;
    }

  };
}

