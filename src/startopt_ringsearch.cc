#include <vector>
#include <string>
#include <stdexcept>
#include <random>
#include <iterator>
#include <algorithm>
#include <utility>
#include "configuration.h"
#include "atomic.h"
#include "Scon/scon_vect.h"
#include "evolution.h"
#include "coords.h"
#include "startopt_ringsearch.h"
#include "configuration.h"
#include "Scon/scon_utility.h"

/*


  ring class


*/

startopt::ringsearch::ring::ring(std::size_t const _s)
  : size(_s), atoms(_s), dihedrals(_s - 3), dihedral_direction(_s - 3)
{ }

startopt::ringsearch::ring::ring(std::size_t const _s,
  std::vector<std::size_t> _v) :
  size(_s), atoms(_v.begin(), _v.begin() + size),
  dihedrals(size - 3), dihedral_direction(size - 3)
{ }

startopt::ringsearch::ring::ring(ring const& r) :
  size(r.size), atoms(r.atoms), dihedrals(r.dihedrals),
  dihedral_direction(r.dihedral_direction)
{ }

startopt::ringsearch::ring& startopt::ringsearch::ring::operator= (
  startopt::ringsearch::ring const& r)
{
  if (size == r.size)
  {
    atoms = r.atoms;
    dihedrals = r.dihedrals;
    dihedral_direction = r.dihedral_direction;
  }
  return *this;
}

/*


Search class


*/

startopt::ringsearch::Search::Search(coords::Coordinates& coords) :
  m_init_xyz(coords.xyz()), m_final_ensemble(),
  m_coord(coords), m_is_acceptor(coords.size(), false),
  m_rings(), m_overlap(), m_ringcontainer(7u)
{
  std::size_t const n = m_coord.size();
  for (std::size_t i(0U); i < n; ++i)
  {
    std::size_t const atomic_number_i = m_coord.atoms(i).number();
    if (atomic_number_i == 1 &&
      m_coord.atoms(i).bonds().size() == 1 &&
      atomic::number_is_heteroatom(
        m_coord.atoms(m_coord.atoms(i).bonds(0)).number()))
    {
      m_donors.push_back(i);
    }
    else if (atomic::number_is_heteroatom(atomic_number_i))
    {
      m_is_acceptor[i] = true;
    }
  }
  find_rings();
}

startopt::ringsearch::Search::~Search()
{
  m_coord.set_xyz(m_init_xyz, true);
  m_coord.to_internal();
}

std::size_t startopt::ringsearch::Search::find_rings()
{
  std::size_t const n = m_coord.size();
  std::vector<bool> tabu_vector(n, false);
  std::size_t const m = m_donors.size();
  for (std::size_t i = 0; i < m; ++i)
  {
    m_ringcontainer[0u] = m_donors[i];
    find_route(m_donors[i], 0u, std::vector<bool>(n, false));
  }
  if (Config::get().general.verbosity > 3)
  {
    std::cout << "Found " << m_rings.size()
      << " potential rings.\n";
    if (Config::get().general.verbosity > 4)
    {
      std::cout << "Ring-routes:\n";
      for (auto const& r : m_rings)
      {
        std::cout << "[" << r.size << "]: ";
        for (auto const& a : r.atoms)
        {
          std::cout << a << " ";
        }
        std::cout << "\n";
      }
    }
  }
  find_torsions();
  find_overlap();
  return m_rings.size();
}

void startopt::ringsearch::Search::find_route(std::size_t const i,
  std::size_t const size, std::vector<bool> tabulist)
{
  if (size > 5) return;
  tabulist[i] = true;
  const std::size_t increased = size + 1;
  for (auto const bound : m_coord.atoms(i).bonds())
  {
    m_ringcontainer[increased] = bound;
    if (tabulist[bound]) continue;
    else if (m_is_acceptor[bound] && increased > 3)
    {
      m_rings.push_back(ring(increased + 1u, m_ringcontainer));
    }
    find_route(bound, increased, tabulist);
  }
}

std::size_t startopt::ringsearch::Search::find_torsion(double& direction,
  std::size_t const index_b, std::size_t const index_c) const
{
  const std::size_t M = m_coord.atoms().mains().size();
  for (std::size_t i = 0u; i < M; ++i)
  {
    std::size_t const index = m_coord.atoms().intern_of_main_idihedral(i);
    std::size_t const ib = m_coord.atoms(m_coord.atoms(index).ibond()).i_to_a();
    std::size_t const ia = m_coord.atoms(m_coord.atoms(index).iangle()).i_to_a();
    if (index_b == ib && index_c == ia)
    {
      direction = 1.0;
      return i;
    }
    else if (index_b == ia && index_c == ib)
    {
      direction = -1.0;
      return i;
    }
  }
  std::cout << "Indices " << index_b << " / " << index_c << "\n";
  throw std::runtime_error("Error: Could not match a ring bond to a main torsion of the system.");
  //return 0;
}

void startopt::ringsearch::Search::find_torsions()
{
  const std::size_t M = m_rings.size();
  for (std::size_t i = 0; i < M; ++i)
  {
    m_rings[i].dihedrals[0] = find_torsion(m_rings[i].dihedral_direction[0], m_rings[i].atoms[1], m_rings[i].atoms[2]);
    m_rings[i].dihedrals[1] = find_torsion(m_rings[i].dihedral_direction[1], m_rings[i].atoms[2], m_rings[i].atoms[3]);
    if (m_rings[i].size > 5)
    {
      m_rings[i].dihedrals[2] = find_torsion(m_rings[i].dihedral_direction[2], m_rings[i].atoms[3], m_rings[i].atoms[4]);
    }
    if (m_rings[i].size > 6)
    {
      m_rings[i].dihedrals[3] = find_torsion(m_rings[i].dihedral_direction[3], m_rings[i].atoms[4], m_rings[i].atoms[5]);
    }
  }
}

void startopt::ringsearch::Search::find_overlap()
{
  std::size_t const N(m_rings.size()), M = (N * N - N) / 2;
  m_overlap.resize(N, false);
  for (std::size_t i(0u), row(1u), col(0u); i < M; ++i)
  {
    for (auto dih_a : m_rings[row].dihedrals)
    {
      for (auto dih_b : m_rings[col].dihedrals)
      {
        std::size_t const rel_a_of_a(m_coord.atoms(dih_a).iangle()),
          rel_d_of_a(m_coord.atoms(dih_a).idihedral());
        std::size_t const rel_a_of_b(m_coord.atoms(dih_b).iangle()),
          rel_d_of_b(m_coord.atoms(dih_b).idihedral());
        if ((rel_a_of_a == rel_a_of_b && rel_d_of_a == rel_d_of_b)
          || (rel_a_of_a == rel_d_of_b && rel_d_of_a == rel_a_of_b))
        {
          m_overlap(row, col) = true;
        }
      }
    }
    ++col;
    if (col == row)
    {
      col = 0;
      ++row;
    }
  }
}

void startopt::ringsearch::Search::bias_ring(std::size_t const index,
  coords::float_type const force)
{
  if (index < m_rings.size())
  {
    config::biases::dihedral bias;
    bias.force = force;
    auto const num_dihedrals = m_rings[index].size - 3u;
    auto const ring_dih_index = num_dihedrals - 2u;
    auto const n_atoms = m_coord.size();
    for (std::size_t i = 0; i < num_dihedrals; ++i)
    {
      bias.ideal = RING_DIH[ring_dih_index][i];
      bias.a = m_rings[index].atoms[0 + i];
      bias.b = m_rings[index].atoms[1 + i];
      bias.c = m_rings[index].atoms[2 + i];
      bias.d = m_rings[index].atoms[3 + i];
      // if a, b, c and d are valid atoms we add bias
      if (bias.a < n_atoms && bias.b < n_atoms &&
        bias.c < n_atoms && bias.d < n_atoms)
      {
        m_coord.potentials().add(bias);
      }
    }
  }
}

std::vector<std::size_t> startopt::ringsearch::Search::bias_rings(
  std::vector<bool> const& close_it,
  coords::float_type const force /*= 0.1*/)
{
  std::size_t const n = m_rings.size();
  auto actually_closed = std::vector<std::size_t>{};
  if (close_it.size() == n)
  {
    for (std::size_t i(0u); i < n; ++i)
    {
      if (close_it[i])
      {
        bool can_be_closed(true);
        for (std::size_t j(0u); j < i; ++j)
        {
          if (close_it[j] && m_overlap(i, j) == true)
          {
            can_be_closed = false;
            break;
          }
        }
        if (can_be_closed)
        {
          bias_ring(i, force);
          actually_closed.push_back(i);
        }
      }
    }
  }
  return actually_closed;
}

//************************************
// Method:    energy
// FullName:  startopt::ringsearch::Search::energy
// Access:    public 
// Returns:   coords::float_type 
//            (energy of relaxated structure with closed rings)
// Qualifier:
// Parameter: std::vector<bool> const & close_it 
//            (list of rings to be closed)
// Parameter: coords::float_type const force 
//            (force to be applied for bias potential)
//************************************

coords::float_type startopt::ringsearch::Search::energy(
  std::vector<bool> const& close_it,
  coords::float_type const force)
{
  // Set original coordinates
  m_coord.set_xyz(m_init_xyz, true);
  // Save original bias
  auto tmp_pot = m_coord.potentials();
  // Imply required closure bias potentials
  // returns actually closed rings for checkback
  /*auto closed_rings = */bias_rings(close_it, force);
  // Optimize with bias in place
  m_coord.o();
  // Remove ring closure bias
  m_coord.potentials().swap(tmp_pot);
  // Relaxate without artificial bias
  return m_coord.o();
}

/*


Genetic propagation function


*/

namespace
{

  struct rings_individual
  {
    std::vector<bool> close_state;
    coords::PES_Point pes;
    coords::float_type health;
    rings_individual(std::vector<bool> init_state) :
      close_state(std::move(init_state)), pes(), health() {}
    rings_individual(std::size_t n) : close_state(n), pes(), health() {}
  };

  template<class Rng>
  rings_individual randomized_rings_individual(
    std::size_t const n_rings, Rng& rng)
  {
    rings_individual r{ std::vector<bool>(n_rings, false) };
    std::generate_n(r.close_state.begin(), r.close_state.size(),
      [&rng]() -> bool
    {
      return std::bernoulli_distribution{
        Config::get().startopt.ringsearch.chance_close }(rng);
    });
    return r;
  }

  struct ring_ind_less
  {
    bool operator() (rings_individual const& a,
      rings_individual const& b)
    {
      return (a.pes.integrity && !b.pes.integrity) ||
        a.health < b.health;
    }
  };

  template<class Fitness>
  struct ring_population_fitness
  {
    startopt::ringsearch::Search& rso;
    Fitness fitness;
    // create object
    template<class ... Ts>
    ring_population_fitness(
      startopt::ringsearch::Search& search_object,
      Ts&& ... fits)
      : rso(search_object), fitness(std::forward<Ts>(fits)...) {}
    std::vector<coords::float_type> operator() (
      std::vector<rings_individual>& pop, std::size_t const)
    {
      auto const n = pop.size();
      auto values = std::vector<coords::float_type>(pop.size());
      auto mean = coords::float_type{ 0 };
      // get health and structures for population
      for (auto& individual : pop)
      {
        // get energy of individual
        individual.health = rso.energy(individual.close_state,
          Config::get().startopt.ringsearch.bias_force);
        mean += individual.health;
        individual.pes = rso.get_coords().pes();
      }
      // sort population by integrity and health
      std::sort(pop.begin(), pop.end(), ring_ind_less{});

      // Get fitness values
      auto fitness_values = std::vector<coords::float_type>(n);
      for (std::size_t i = 0; i < n; ++i)
      {
        fitness_values[i] = fitness(i);
      }
      // print individuals and energies
      if (Config::get().general.verbosity > 3)
      {
        mean /= n;
        std::cout << "Mean population energy value = " <<
          mean << "\n";
        if (Config::get().general.verbosity > 4)
        {
          for (auto& individual : pop)
          {
            for (bool s : individual.close_state)
            {
              std::cout << (s ? 1 : 0);
            }
            std::cout << " " << individual.health << "\n";
          }
        }
      }
      return fitness_values;
    }
  };

  template<class Fitness>
  ring_population_fitness<Fitness> make_fitness(
    startopt::ringsearch::Search& search_object,
    Fitness&& fit)
  {
    return{ search_object, std::forward<Fitness>(fit) };
  }

}

coords::Ensemble_PES startopt::ringsearch::genetic_propagation(
  coords::Coordinates& coords, std::size_t const population_count,
  std::size_t const generation_count)
{
  // typedefs
  using real_dist = std::uniform_real_distribution < coords::float_type >;
  using int_dist = std::uniform_int_distribution < std::size_t >;
  using bool_dist = std::bernoulli_distribution;
  using lin_rank_fitness = genetic::ranking::linear<coords::float_type>;
  using exp_rank_fitness = genetic::ranking::exponential<coords::float_type>;
  using ind_type = rings_individual;
  using pop_type = std::vector<ind_type>;
  // Mersenne Twister random engine
  auto rng = std::mt19937_64(std::random_device()());
  // linear ranking fitness
  auto lin_rank = lin_rank_fitness{ population_count,
    Config::get().optimization.global.selection.low_rank_fitness,
    Config::get().optimization.global.selection.high_rank_fitness };
  // exponential ranking fitness
  auto exp_rank = exp_rank_fitness{ population_count,
    Config::get().optimization.global.selection.low_rank_fitness,
    Config::get().optimization.global.selection.high_rank_fitness };
  // Ringsearch object
  Search rso(coords);
  // fitness callback
  auto pop_fitness = make_fitness(rso,
    [&](std::size_t const i) -> coords::float_type
  {
    return Config::get().optimization.global.selection.fit_type ==
      config::optimization_conf::sel::fitness_types::LINEAR ?
      lin_rank(i) : exp_rank(i);
  });
  // mating callback
  auto mating = [&](pop_type const& parents,
    std::vector<coords::float_type> const& fitness_values,
    std::size_t const) -> pop_type
  {
    // calculate sum of fitness values
    coords::float_type fitness_sum = 0;
    for (auto&& val : fitness_values) fitness_sum += val;
    // number of individuals
    auto const n = parents.size();
    // number of close states
    auto const m = n > 0 ? parents.front().close_state.size() : 0u;
    // child population
    pop_type children;
    children.reserve(n);
    children.push_back(parents.at(0u));
    for (std::size_t i = 1; i < n; ++i)
    {
      // create child
      children.emplace_back(std::vector<bool>(m, false));
      // select parents
      auto p1 = genetic::selection::roulette(fitness_values, fitness_sum);
      auto p2 = genetic::selection::roulette(fitness_values, fitness_sum);
      // combine child states from parents
      for (std::size_t j = 0; j < m; ++j)
      {
        // fitness weighted pick of parents for state i 
        auto const pick = real_dist{ 0,
          fitness_values[p1] + fitness_values[p2] }(rng);

        children.back().close_state[j] = parents.at(
          pick > fitness_values[p1] ? p1 : p2).close_state[j];
      }
    }
    return children;
  };
  // matation callback
  auto mutating = [&](pop_type& pop, std::size_t const)
  {
    auto&& go = Config::get().optimization.global;
    for (auto& i : pop)
    {
      // mutate as often as we hit the probability to do so
      while (bool_dist{ go.evolution.chance_pointmutation }(rng))
      {
        // select state j to change
        auto j = int_dist{ 0, i.close_state.size() - 1u }(rng);
        // flip state j with 50:50 chance
        i.close_state.at(j) = bool_dist{}(rng);
      }
    }
  };
  // initial population
  auto population = pop_type(population_count,
    ind_type{ std::vector<bool>(rso.rings().size(), false) });
  auto const n_rings = rso.rings().size();
  // generate randomized individuals
  std::generate_n(population.begin(), population_count, [&]() -> ind_type
  {
    return randomized_rings_individual(n_rings, rng);
  });
  // evolution
  population = genetic::evolve(population, generation_count,
    pop_fitness, mating, mutating);
  // return ensemble
  coords::Ensemble_PES pese;
  pese.reserve(population.size());
  for (auto& i : population)
  {
    pese.emplace_back(std::move(i.pes));
  }
  return pese;
}

/*


Preoptimizer class generate function


*/

void startopt::preoptimizers::R_evolution::generate(
  coords::Ensemble_PES const& init_ensemble, std::size_t const m)
{
  auto const n = init_ensemble.size();
  for (std::size_t i = 0u; i < n; ++i)
  {
    // set current structure to init_ensemble element i
    m_final_coords.set_pes(init_ensemble[i]);
    // propagate m individuals
    auto en = ringsearch::genetic_propagation(m_final_coords,
      m, Config::get().startopt.ringsearch.generations);
    // reserve space for final ensemble
    //m_ensemble.reserve(m_ensemble.size() + en.size());
    // move created elements into m_ensemble 
    m_ensemble.insert(m_ensemble.end(),
      std::make_move_iterator(en.begin()),
      std::make_move_iterator(en.end()));
  }
}
