#include "md_logging.h"
#include "md.h"
#include "helperfunctions.h"
#pragma once



std::ostream& md::operator<<(std::ostream& strm, trace_data const& d)
{
  strm << d.i << ",";
  strm << std::fixed << std::setprecision(3) << d.T << ",";
  strm << std::fixed << std::setprecision(5) << d.P << ",";
  strm << std::fixed << std::setprecision(5) << d.Ek << ",";
  strm << std::fixed << std::setprecision(5) << d.Ep << ",";
  strm << std::fixed << std::setprecision(5) << d.Ek + d.Ep << ",";
  if (d.snapshot > 0u) strm << d.snapshot;
  else strm << std::right << std::setw(20) << '-';
  for (auto iae : d.Eia)
  {
    strm << "," << std::fixed << std::setprecision(5) << iae;
  }
  strm << '\n';
  return strm;
}

void md::trace_writer::operator() (md::trace_data const& d)
{
  static std::atomic<bool> x(true);
  if (x && Config::get().general.verbosity > 1)
  {
    *strm << "It,";
    *strm << "T,";
    *strm << "P,";
    *strm << "kin.En.,";
    *strm << "pot.En.,";
    *strm << "tot.En.,";
    *strm << "Snapsh.";
    if (d.Eia.size() > 1u)
    {
      for (auto i : scon::index_range(d.Eia))
      {
        *strm << "," << "E_ia(# " << i << ')';
      }
    }
    *strm << '\n';
    x = false;
  }
  *strm << d;
}


md::Logger::Logger(coords::Coordinates& coords, std::size_t snap_offset) :
  snap_buffer(coords::make_buffered_cartesian_log(coords, "_MD_SNAP",
    Config::get().md.max_snap_buffer, snap_offset, Config::get().md.optimize_snapshots)),
  velo_buffer(make_buffered_velocity_log(coords, "_MD_VELO", Config::get().md.max_snap_buffer, snap_offset)),
  data_buffer(scon::offset_call_buffer<trace_data>(50u, Config::get().md.trackoffset,
    trace_writer{ coords::output::filename("_MD_TRACE", ".csv").c_str() })),
  snapnum()
{
}

bool md::Logger::operator()(std::size_t const iter, coords::float_type const T,
  coords::float_type const P, coords::float_type const Ek, coords::float_type const Ep,
  std::vector<coords::float_type> const Eia, coords::Representation_3D const& x, coords::Representation_3D const& v)
{
  if (iter % 5000u == 0)
  {
    if (std::isnan(T) || std::isnan(P) || std::isnan(Ek) || std::isnan(Ep))
    {
      std::cout << "NaN in simulation, please check your input options" << std::endl;
      throw std::logic_error("NaN in simulation, please check your input options");
    }
  }
  velo_buffer(v);   // writing velocities into file (snap_buffer(x) in next line writes snapshots into file)
  return data_buffer(trace_data(Eia, T, Ek, Ep, P, iter, snap_buffer(x) ? ++snapnum : 0u));
}


// Serialization Helper Function

static inline void append_to_buffer(std::vector<char>& buffer, void const* data, std::size_t const bytes)
{
  if (bytes > 0 && data)
  {
    buffer.resize(buffer.size() + bytes);
    std::vector<char>::iterator ins_pos_it(buffer.end() - bytes);
    memcpy(&(*ins_pos_it), data, bytes);
  }
}

static inline std::istream::pos_type istream_size_left(std::istream& S)
{
  std::istream::pos_type const T(S.tellg());
  S.seekg(0, std::istream::end);
  std::istream::pos_type const R(S.tellg());
  S.seekg(T, std::istream::beg);
  return R;
}


// Traces 


md::simulation::simulation(coords::Coordinates& coord_object) :
  coordobj(std::addressof(coord_object)),
  logging(coord_object, gap(Config::get().md.num_steps, Config::get().md.num_snapShots)),
  P(coord_object.xyz()), P_old(coord_object.xyz()),
  F(coord_object.g_xyz()), F_old(coord_object.g_xyz()),
  V(coord_object.xyz().size()), M(coord_object.xyz().size()),
  M_total(0.0), E_kin(0.0), desired_temp(Config::get().md.T_init), instantaneous_temp(0.0), press(0.0), dt(Config::get().md.timeStep),
  freedom(coord_object.size() * 3u), snapGap(0), C_geo(), C_mass(),
  thermostat(md::nose_hoover_arbitrary_length(std::vector<double>(Config::get().md.nosehoover_chainlength, Config::get().md.nosehoover_Q)),md::nose_hoover_2chained(Config::get().md.nosehoover_Q)),
  rattle_bonds(), window(), restarted(true)
{
  std::sort(Config::set().md.heat_steps.begin(), Config::set().md.heat_steps.end());

  if (Config::get().md.ana_pairs.size() > 0) md_analysis::create_ana_pairs(this);   // create atom pairs to analyze, fetch information and save
}

void md::simulation::write_restartfile(std::size_t const k)
{
  // Create binary buffer and copy data into vector char
  scon::binary_stream<std::vector<char>> buffer;
  buffer.v.reserve(scon::binary_size(k, *this));
  buffer << k << *this;
  // Create ofstream and insert buffer into stream
  std::ofstream restart_stream(
    (Config::get().general.outputFilename + "_MD_restart.cbf").c_str(),
    std::ofstream::out | std::ofstream::binary
  );
  restart_stream.write(buffer.v.data(), buffer.v.size());
}
