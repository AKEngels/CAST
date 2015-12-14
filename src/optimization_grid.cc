#include <limits>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <cmath>

#include "coords.h"
#include "optimization_global.h"
#include "scon_vect.h"
#include "scon_utility.h"
#include "configuration.h"
#include "scon_chrono.h"
#include "scon.h"

optimization::global::optimizers::main_grid::main_grid(coords::Coordinates &c, 
  coords::Ensemble_PES const &p, coords::angle_type const delta_grid)
  : optimizer(c,p), m_delta(delta_grid),
  m_num_offsets(static_cast<std::size_t>(coords::float_type(360) / delta_grid.degrees())+1u),
  m_offset(c.main().size(), 0u), 
  m_init_offset(c.main().size(), 0u),
  m_done(false)
{
  std::size_t const n = c.main().size();
  for (std::size_t k = 0; k < n; ++k)
  {
    m_init_offset[k] = static_cast<std::size_t>(std::round(coords::angle_type(c.main()[k]).radians() / m_delta.radians()));
  }
}

void optimization::global::optimizers::main_grid::increase_offset(std::size_t it)
{
  std::size_t const n = this->coordobj.main().size();
  if (it < n)
  {
    if (m_offset[it] < m_num_offsets) ++m_offset[i];
    if (m_offset[it] == m_num_offsets)
    {
      increase_offset(it + 1);
      m_offset[it] = 0;
    }
  }
  else
  {
    m_done = true;
  }
}


coords::Representation_Main optimization::global::optimizers::main_grid::next_main_from_offset()
{
  std::size_t const n = this->coordobj.main().size();
  coords::Representation_Main ret;
  ret.resize(n);
  increase_offset(0u);
  for (std::size_t k = 0; k < n; ++k)
  {
    ret[k] = coords::main_type::from_rad(coords::float_type(m_init_offset[k] + m_offset[k])*m_delta.radians());
  }
  return ret;
}


bool optimization::global::optimizers::main_grid::run(std::size_t const iterations, bool const)
{
  coordobj.set_pes(accepted_minima[min_index].pes);
  header_to_cout();
  std::size_t const iter_size(num_digits(Config::get().optimization.global.iterations) + 1);
  while (!m_done && i++<iterations)
  {
    scon::chrono::high_resolution_timer step_timer;
    coordobj.set_all_main(next_main_from_offset());
    //std::cout << "Mains: ";
    //for (auto a : coordobj.main()) std::cout << a << " ";
    //std::cout << "\n";
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << "GS  ";
      // iteration
      std::cout << std::right << std::setw(iter_size) << i + 1 << "/";
      std::cout << std::left << std::setw(iter_size) << Config::get().optimization.global.iterations;
      // Curennt minimum index
      std::cout << std::setw(10) << min_index;
      // Current minimum
      std::cout << std::setw(Config::get().optimization.global.precision + 10) << std::right;
      std::cout << std::setprecision(Config::get().optimization.global.precision);
      std::cout << std::showpoint << std::scientific << accepted_minima[min_index].pes.energy;
      // Transition
      std::cout << std::setw(Config::get().optimization.global.precision + 10) << std::right;
      std::cout << std::setprecision(Config::get().optimization.global.precision);
      std::cout << std::showpoint << std::scientific << coordobj.pes().energy;
    }
    if (coordobj.preoptimize()) coordobj.po();
    double ene = coordobj.o();
    coordobj.to_internal();
    coordobj.to_xyz();
    min_status::T const status(check_pes_of_coords());
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << std::setw(Config::get().optimization.global.precision + 10) << std::right;
      std::cout << std::setprecision(Config::get().optimization.global.precision);
      std::cout << std::showpoint << std::scientific << ene;
      std::cout << status;
      std::cout << std::setw(8) << std::left << accepted_minima.size();
      std::cout << std::setw(8) << std::left << range_minima.size();
    }
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << "   (" << std::setprecision(2) << std::showpoint;
      std::cout << std::fixed << T << " K, " << step_timer << ")" << lineend;
    }
    if (!restore(status))
    {
      if (Config::get().general.verbosity > 1U)
        std::cout << "Starting point selection limit reached (no non-tabu minimum accessible). Stop." << lineend;
      break;
    }
  }
  return found_new_minimum;
}
