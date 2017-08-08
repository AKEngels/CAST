#include "optimization.h"

template<> 
double const optimization::constants<double>::kB = 0.001987204118;
template<> 
float const optimization::constants<float>::kB = 0.001987204118f;


#include <string>
#include <fstream>
#include "coords.h"
#include "optimization_global.h"
#include "coords_io.h"
#include "configuration.h"
#include "startopt_solvadd.h"
#include "scon_utility.h"


bool optimization::global::Tabu_List::tabu (coords::PES_Point const &point, coords::Coordinates const &co) const
{
  if (empty() || point.energy < (front().pes.energy-1.0)) 
  {
    return false;
  }
  if ((back().pes.energy+1.0) < point.energy) 
  {
    return false;
  }
  double const low_bound(point.energy - 1.0), high_bound(point.energy + 1.0);
  base_type::size_type const n = this->size();
  for (base_type::size_type i(0u); i < n; ++i)
  {
    if ((*this)[i].pes.energy < high_bound && (*this)[i].pes.energy > low_bound)
    {
      if (co.equal_structure((*this)[i].pes, point))
      {
        return true;
      }
    }
  }
  //std::cout << "Did not find any structures in list to make it (" << point.energy << ") tabu\n";
  return false;
}


bool optimization::global::Tabu_List::has_superposition (coords::PES_Point const &point, coords::Coordinates const &coord) const
{
  for (auto const & tp : *this)
  {
    if (coord.check_superposition_xyz(tp.pes.structure.cartesian, point.structure.cartesian)) return true;
  }
  return false;
}


void optimization::global::Tabu_List::clear_above_e (double const minimum, double const kT)
{
  while (!empty() && back().pes.energy > minimum && std::exp(-(back().pes.energy-minimum)/kT) < 1.0E-5)
  {
    pop_back();
  }
}


std::unique_ptr<optimization::global::optimizer>
optimization::global::new_divers_optimizer(coords::Coordinates & c)
{
  return std::unique_ptr<optimization::global::optimizer>{ 
    new optimization::global::optimizers::monteCarlo(c) };
}

void optimization::global::optimizer::write_accepted(std::string const & suffix)
{
  if (accepted_minima.size() > 0)
  {
    std::ofstream final_accepted_log(coords::output::filename(suffix + "_accepted_final", ".log"), std::ios::out);
    std::ofstream final_accepted_arc(coords::output::filename(suffix + "_accepted_final"), std::ios::out);
    if (final_accepted_log && final_accepted_arc)
    {
      final_accepted_log << std::setw(10) << "Num";
      final_accepted_log << std::setw(12) << "Iteration";
      final_accepted_log << std::setw(20) << "Energy";
      final_accepted_log << std::setw(20) << "Interactions ...\n";
      std::size_t o = 0;
      for (auto const & am : accepted_minima)
      {
        final_accepted_log << std::setw(10) << ++o;
        final_accepted_log << std::setw(12) << am.iteration;
        final_accepted_log << std::setw(20) << am.pes.energy;
        for (auto && ia : am.pes.ia_matrix)
        {
          final_accepted_log << std::setw(20) << ia.energy;
        }
        final_accepted_log << '\n';
        coordobj.set_xyz(am.pes.structure.cartesian, true);
        final_accepted_arc << coordobj;
      }
    }
  }
}



optimization::global::result_drain::result_drain(coords::Coordinates &c, 
  std::string const & filename_suffix) : cp(&c),
  energy_stream(new std::ofstream(coords::output::filename(filename_suffix, ".log"), std::ios::out)),
  structure_stream(new std::ofstream(coords::output::filename(filename_suffix), std::ios::out)), n()
{
  if (energy_stream && *energy_stream)
  {
    *energy_stream << std::setw(10) << "Num";
    *energy_stream << std::setw(12) << "Iteration";
    *energy_stream << std::setw(20) << "Energy";
    *energy_stream << std::setw(20) << "Interactions ...\n";
  }
}


void optimization::global::result_drain::operator()(Result_Point && r)
{
  if (energy_stream && *energy_stream)
  {
    *energy_stream << std::setw(10) << ++n;
    *energy_stream << std::setw(12) << r.it;
    *energy_stream << std::setw(20) << r.pes.energy;
    for (auto & ia : r.pes.ia_matrix)
    {
      *energy_stream << std::setw(20) << ia.energy;
    }
    *energy_stream << '\n';
  }
  if (structure_stream && *structure_stream)
  {
    auto tmp = cp->pes();
    cp->set_pes(r.pes);
    *structure_stream << *cp;
    cp->set_pes(tmp);
  }
}

optimization::global::offset_buffered_resultfile 
optimization::global::make_buffered_result(coords::Coordinates &coord_obj,
  std::string const & file_suffix, std::size_t buffer_size)
{
  return scon::offset_call_buffer<Result_Point>(
    buffer_size, 1u, result_drain{ coord_obj, file_suffix });
}


optimization::global::optimizer::optimizer (coords::Coordinates &c, std::string const & output_name) :
  T(Config::get().optimization.global.temperature), kT(T*optimization::global::k),
  tabulist(), range_tabu(), accepted_minima(), range_minima(),
  init_stereo (c.stereos()), fallbacks(0), i(0), min_index (0), gmin_index (0),
  coordobj(c), accepted_log(make_buffered_result(c, output_name + "_accepted", 5)),
  opt_clock(), found_new_minimum(false)
{
}

optimization::global::optimizer::optimizer (
  coords::Coordinates &c, coords::Ensemble_PES const &initial_structures, 
  bool const minimized, std::string const & output_name) :
  T(Config::get().optimization.global.temperature), kT(T*optimization::global::k),
  tabulist(), range_tabu(), accepted_minima(), range_minima(), 
  init_stereo (c.stereos()), fallbacks(0), i(0), min_index (0), gmin_index (0),
  coordobj(c), accepted_log(make_buffered_result(c, output_name + "_accepted", 5)),
  opt_clock(), found_new_minimum(false)
{ 
  std::ofstream test("test.txt");
    test << "1" << '\n';
    test.close();

  if (Config::get().general.verbosity > 1)
  {
    std::cout << "Evaluating initial structures for global optimization.\n";
  }
  std::size_t const N(initial_structures.size());
  std::size_t const iter_size(scon::num_digits(N));
  coords::Ensemble_3d brokens;
  for (std::size_t kit = 0; kit < N; ++kit)
  {

    coordobj.set_xyz(initial_structures[kit].structure.cartesian);
    if (!minimized)
    {
      if (coordobj.preoptimize())
      {
        coordobj.po();
      }
      coordobj.o();
      coordobj.to_internal();
      coordobj.to_xyz();
    }
	
    else
    {
      std::ofstream test("test.txt");
      test << "1" << '\n';
      test.close();

      coordobj.set_pes(initial_structures[kit]);
    }


    std::ofstream test2("test2.txt");
    test2 << "2" << '\n';
    test2.close();

    min_status::T S(check_pes_of_coords());
    if (Config::get().general.verbosity > 1)
    {
      std::cout << "Status for initial structure #" << std::setw(iter_size) << std::left << kit + 1 << (minimized ? "" : " (minimized)");
      std::cout << ": " << S << "(Energy: " << coordobj.pes().energy << ", Minimum Index: " << min_index << ")\n";
      if (S == min_status::T::REJECT_BROKEN)
      {
        brokens.push_back(coordobj.xyz());
      }
    }
  }


  std::ofstream test3("test3.txt");
  test3 << "3" << '\n';
  test3.close();

  if (Config::get().general.verbosity > 1)
  {
    std::cout << '\n';
  }
  if (!brokens.empty())
  {
    auto filename = coords::output::formats::tinker::filename("_BROKEN_INIT");
    if (Config::get().general.verbosity > 1) std::cout << "Writing broken initial structures to '" << filename << "'\n";
    std::ofstream outputstream(filename.c_str(), std::ios_base::out);
    for (auto const & structure : brokens)
    {


      std::ofstream test4("test4.txt");
      test4 << "4" << '\n';
      test4.close();

      coordobj.set_xyz(structure);


      std::ofstream test5("test5.txt");
      test5 << "5" << '\n';
      test5.close();

      coordobj.e();


      std::ofstream test6("test6.txt");
      test6 << "6" << '\n';
      test6.close();

      outputstream << coords::output::formats::tinker(coordobj);
    }


    std::ofstream test7("test7.txt");
    test7 << "7" << '\n';
    test7.close();

  }
  if (!found_new_minimum)
  {
    throw std::runtime_error("No valid initial structure for global optimization.");
  }
}


optimization::global::optimizer::min_status::T optimization::global::optimizer::check_pes_of_coords()
{
  if (coordobj.pes().integrity)
  {
    if (Config::get().optimization.global.delta_e > 0.0 && (!found_new_minimum ||
      (coordobj.pes().energy - accepted_minima[gmin_index].pes.energy) < Config::get().optimization.global.delta_e))
    {
      updateRange(coordobj.pes());
    }
    if (!found_new_minimum || accept(coordobj.pes().energy))
    {
      if (init_stereo != coordobj.stereos())
      {
        return min_status::T::REJECT_STEREO;
      }
    }
    else return min_status::T::REJECT_ENERGY;
  }
  else return min_status::T::REJECT_BROKEN;
  if (!tabulist.tabu(coordobj.pes(), coordobj) && 
    !tabulist.has_superposition(coordobj.pes(), coordobj))
  {
    return new_minimum();
  }
  else return min_status::T::REJECT_TABU;
}

optimization::global::optimizer::min_status::T optimization::global::optimizer::new_minimum()
{
  bool global(false);
  accepted_minima.push_back(Tabu_Point(coordobj.pes()));
  //std::cout << "Mains of new minimum: " << scon::vector_delimeter(' ') << 
  //  "LM:" << coordobj.pes().structure.main << "\n";
  min_index = accepted_minima.size() - 1;
  if (!found_new_minimum || accepted_minima.back().pes.energy < accepted_minima[gmin_index].pes.energy)
  {
    global = true;
    gmin_index = min_index;
    tabulist.clear_above_e(accepted_minima[gmin_index].pes.energy, kT);
    write_range();
  }
  setTemp(T*Config::get().optimization.global.temp_scale);
  accepted_log(i, accepted_minima.back().pes);
  //accepted_iteration.push_back(i);
  scon::sorted::insert(tabulist, accepted_minima.back());
  found_new_minimum = true;
  return global ? min_status::ACCEPT_GLOBAL_MINIMUM : min_status::ACCEPT_MINIMUM;
}


void optimization::global::optimizer::header_to_cout()
{
  std::size_t const iter_size(scon::num_digits(
    Config::get().optimization.global.iterations) + 1);
  if (Config::get().general.verbosity > 1U)
  {
    std::cout << "M   ";
    std::cout << std::left << std::setw(2 * iter_size + 1) << "Iter" << ' ';
    std::cout << std::setw(10) << std::left << "C.Index" << ' ';
    std::cout << std::setw(Config::get().optimization.global.precision + 10) << "E.Current" << ' ';
    std::cout << std::setw(Config::get().optimization.global.precision + 10) << "E.Trans" << ' ';
    std::cout << std::setw(Config::get().optimization.global.precision + 10) << "E.Next" << ' ';
    std::cout << "Accepted?       ";
    std::cout << std::setw(10) << std::left << "N.Minima";
    std::cout << std::setw(10) << std::left << "N.Range";
    std::cout << "(Temperature, Timer)\n";
  }
}

std::ostream & optimization::global::operator<< (std::ostream &strm, optimization::global::optimizer::min_status::T const &S)
{
  typedef optimization::global::optimizer::min_status min_status;
  switch (S)
  {
    case min_status::T::ACCEPT_GLOBAL_MINIMUM:
    {
      strm << "ACCEPT (  GM  ) ";
      break;
    }
    case min_status::T::ACCEPT_MINIMUM:
    {
      strm << "ACCEPT (  OK  ) ";
      break;
    }
    case min_status::T::REJECT_BROKEN:
    {
      strm << "REJECT (broken) ";
      break;
    }
    case min_status::T::REJECT_ENERGY:
    {
      strm << "REJECT (energy) ";
      break;
    }
    case min_status::T::REJECT_STEREO:
    {
      strm << "REJECT (stereo) ";
      break;
    }
    case min_status::T::REJECT_TABU:
    {
      strm << "REJECT ( tabu ) ";
      break;
    }
  }
  return strm;
}

void optimization::global::optimizer::setTemp (const double temp)
{
  T = temp;
  kT = optimization::global::k*T;
}

namespace
{
  
}

bool optimization::global::optimizer::accept (double const E) const
{
  using std::abs;
  // Nan or anything? DO not accept!
  if (E != E) return false;
  double const E0(Config::get().optimization.global.metropolis_local ? 
    accepted_minima[min_index].pes.energy : accepted_minima[gmin_index].pes.energy);
  // More than two orders of magnitude difference?
  // Pobably not sane...
  if ((std::abs(E) / std::abs(E0)) > 100.0) return false;
  // R[0,1) < e^((E0-E)/kT) ?
  auto rand01 = scon::random::threaded_rand(std::uniform_real_distribution<double>{});
  return rand01 < std::exp(-(E-E0)/kT);
}


void optimization::global::optimizer::swap (optimizer &go)
{
  using std::swap;
  swap(T, go.T);
  swap(kT, go.kT);
  swap(found_new_minimum, go.found_new_minimum);
  swap(fallbacks, go.fallbacks);
  swap(i, go.i);
  swap(min_index, go.min_index);
  swap(gmin_index, go.gmin_index);
  tabulist.swap(go.tabulist);
  range_tabu.swap(go.range_tabu);
  accepted_minima.swap(go.accepted_minima);
  range_minima.swap(go.range_minima);
  //accepted_iteration.swap(go.accepted_iteration);
  range_iteration.swap(go.range_iteration);
  init_stereo.swap(go.init_stereo);
  swap(opt_clock, go.opt_clock);
  swap(accepted_log, go.accepted_log);
}

void optimization::global::optimizer::set_current_startpoint(std::size_t const new_index)
{
  if (new_index >= accepted_minima.size()) return;
  min_index = new_index;
  ++accepted_minima[min_index].visited;
  coordobj.set_pes(accepted_minima[min_index].pes);
}

bool optimization::global::optimizer::restore(min_status::T const S)
{
  if (S == min_status::T::ACCEPT_GLOBAL_MINIMUM 
    || S == min_status::ACCEPT_MINIMUM) return true;
  if (accepted_minima.empty()) return false;
  switch (Config::get().optimization.global.fallback)
  {
  case config::optimization_conf::global::fallback_types::LAST_GLOBAL:
  {
    return restore_last_global();
    break;
  }
  default: 
    return restore_roulette_selection();
    break;
  }
}

bool optimization::global::optimizer::restore_last_global()
{
  if (min_index >= accepted_minima.size()) return false;
  if (accepted_minima[min_index].visited
    >= Config::get().optimization.global.fallback_limit)
  {
    if (gmin_index >= accepted_minima.size() 
      || min_index == gmin_index
      || accepted_minima[gmin_index].visited 
        >= Config::get().optimization.global.fallback_limit)
    {
      return false;
    }
    set_current_startpoint(gmin_index);
  }
  else
  {
    set_current_startpoint(min_index);
  }
  return true;
}

bool optimization::global::optimizer::restore_roulette_selection()
{
  // If we do not have a new minimum we need to check back to a certain minimum
  bool new_index(false);
  std::size_t kit(0), m(0);
  std::size_t const accessible(
    std::min(Config::get().optimization.global.selection.included_minima, accepted_minima.size())
    );
  if (accessible > 1)
  {
    // Build a sorted list of all accepted minima
    std::vector<genetic_ranking::index_value_pair> iv_pairs;
    // check vector for every minimum whether it has been evaluated
    std::vector<bool> tried(accessible, false);
    std::size_t const N(accepted_minima.size());
    // Build sorted list of energy<->index pairs
    iv_pairs.reserve(N);
    for (std::size_t u(0); u < N; ++u)
    {
      auto skipit = accepted_minima[u].visited >= 
        Config::get().optimization.global.fallback_limit;
      iv_pairs.push_back({ u, accepted_minima[u].pes.energy, skipit });
      if (skipit) ++m;
    }
    std::sort(iv_pairs.begin(), iv_pairs.end());
    while (!new_index && kit < 100 && m < accessible)
    {
      // Linear fitness (1.0 for best minimum, 0.5 for minimum best+accessible) 
      genetic_ranking::Fitness_Linear lin_fit(accessible, 0.5, 1.0);
      // Roulette Selection using Linear fitness
      genetic_ranking::Roulette<genetic_ranking::Fitness_Linear> linear_roulette;
      // Selected minimum
      std::vector<std::size_t> selected = linear_roulette(lin_fit, accessible, 1);
      if (!selected.empty() && selected.front() < iv_pairs.size())
      {
        if (accepted_minima[iv_pairs[selected.front()].index].visited < Config::get().optimization.global.fallback_limit)
        {
          set_current_startpoint(iv_pairs[selected.front()].index);
          return true;
        }
        else
        {
          if (!tried[selected.front()])
          {
            ++m;
            tried[selected.front()] = true;
          }
        }
      }
      ++kit;
    }
    return false;
  }
  else
  {
    if (accepted_minima[0].visited
      >= Config::get().optimization.global.fallback_limit) return false;
    set_current_startpoint(0);
  }
  return true;
}

void optimization::global::optimizer::write_range(std::string const & suffix)
{
  std::size_t const M(range_minima.size());
  if (M > 0)
  {
    std::ofstream r_file(coords::output::filename(suffix + "_range").c_str(), std::ios_base::out),
      re_file(coords::output::filename(suffix + "_range", ".out").c_str(), std::ios_base::out);
    if (r_file && re_file)
    {
      re_file << std::setw(12) << "Num";
      re_file << std::setw(12) << "Iteration";
      re_file << std::setw(20) << "Energy";
      re_file << std::setw(20) << "Interactions ...\n";
      for (std::size_t it(0u); it < M; ++it)
      {
        re_file << std::setw(12) << it;
        re_file << std::setw(12) << range_iteration[it];
        re_file << std::setw(20) << range_minima[it].pes.energy;
        for (auto const & ia : range_minima[it].pes.ia_matrix)
        {
          re_file << std::setw(20) << ia.energy;
        }
        re_file << "\n";
        coordobj.set_pes(range_minima[it].pes);
        r_file << coordobj;
      }
    }
  }
}


void optimization::global::optimizer::updateRange (coords::PES_Point const &p)
{
  while (!range_minima.empty() && (range_minima.back().pes.energy - 
    accepted_minima[gmin_index].pes.energy) > 
    Config::get().optimization.global.delta_e)
  {
    range_minima.pop_back();
    range_iteration.pop_back();
  }
  if (!range_tabu.tabu(p, coordobj))
  {
    if (range_minima.empty()) 
    {
      range_minima.push_back(p);
      range_iteration.push_back(i);
    }
    else
    {
      std::ptrdiff_t const N(static_cast<std::ptrdiff_t>(range_minima.size()));
      std::ptrdiff_t insertion_point(0U);
      while (insertion_point < N && p.energy > range_minima[insertion_point].pes.energy) ++insertion_point;
      range_minima.insert(range_minima.begin()+insertion_point, p);
      range_iteration.insert(range_iteration.begin()+insertion_point, i);
    }
    scon::sorted::insert(range_tabu, Tabu_Point(p));
  }
}

coords::Coordinates optimization::global::optimizer::dehydrate(coords::Coordinates const & solv_coords)
{
  std::size_t const N(solv_coords.size() - 3 * Config::get().startopt.solvadd.maxNumWater);
  if (N < solv_coords.size())
  {
    coords::Atoms atms;
    coords::Representation_3D pos(N);
    for (std::size_t it(0); it < N; ++it)
    {
      coords::Atom a(solv_coords.atoms(it).number());
      for (auto bond : solv_coords.atoms(it).bonds())
      {
        a.bind_to(bond);
      }
      a.assign_to_system(0u);
      a.set_energy_type(solv_coords.atoms(it).energy_type());
      pos[it] = solv_coords.xyz()[it];
    }
    coords::PES_Point tmpos(pos);
    coords::Coordinates tmp;
    tmp.init_swap_in(atms, tmpos);
    return tmp;
  }
  return solv_coords;
}

coords::Coordinates optimization::global::optimizer::rehydrate(coords::Coordinates const & desolv_coords)
{
  startopt::preoptimizers::Solvadd sap(desolv_coords);
  sap.generate(coords::Ensemble_PES(1U, desolv_coords.pes()), 1U);
  return sap.final_coords();
}
