/**
CAST 3
md.h
Purpose: header for molecular dynamics simulation

@version 1.0
*/

#pragma once 

#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>


#include "coords.h"
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
		Logger(coords::Coordinates& coords, std::size_t snap_offset);

		/**looks every 5000 steps if temperature, pressure or energy is nan and throws an error if yes
		*/
		bool operator() (std::size_t const iter,
			coords::float_type const T,
			coords::float_type const P,
			coords::float_type const Ek,
			coords::float_type const Ep,
			std::vector<coords::float_type> const Eia,
			coords::Representation_3D const& x);

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

