#pragma once

#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <utility>
#include <limits>
#include <iostream>

namespace genetic
{

  namespace ranking
  {

    template<class Float>
    class linear
    {
      linear& operator= (linear const &);
      std::size_t const N;
      Float const inclination, offset;
    public:
      typedef Float float_type;
      linear(std::size_t const N_ind, Float const lower_limit, Float const upper_limit)
        : N(N_ind), inclination(-(upper_limit - lower_limit) / Float(N_ind - 1u)), offset(upper_limit)
      { }
      Float operator() (std::size_t const x) const
      {
        return (N == 0 || x <= N) ? inclination*x + offset : Float(0.);
      }
    };

    template<class Float>
    class exponential
    {
      exponential& operator= (exponential const &);
      std::size_t const N;
      Float const b, a, c;
    public:
      typedef Float float_type;
      exponential(std::size_t const N_ind, Float const L, Float const H)
        : N(N_ind),
        b((L / H - std::exp(-Float(N))) / (Float(1) - std::exp(-Float(N)))),
        a(1 - b), c(H)
      { }
      Float operator() (std::size_t const x) const
      {
        return (N == 0 || x < N) ? (a*std::exp(-Float(x)) + b)*c : Float(0.);
      }
    };

  }

  namespace selection
  {

    template<class float_type>
    std::size_t roulette(std::vector<float_type> const &values, float_type sum)
    {
      static std::mt19937_64 rng = std::mt19937_64{ std::random_device{}() };
      using dist_type = std::uniform_real_distribution<float_type>;
      auto const n = values.size();
      auto target_value = dist_type{ 0, sum }(rng);
      auto accu = float_type(0);
      for (std::size_t i = 0; i < n; ++i)
      {
        accu += values[i];
        if (accu >= target_value) return i;
      }
      throw std::logic_error("Roulette selection did not reach target value.");
    }



    template<class Fitness, class Float = typename Fitness::float_type>
    class Roulette
    {
    private:
      std::mt19937_64 rng;
    public:
      typedef Float float_type;
      Roulette() : rng(std::random_device()()) { }
      std::size_t operator()(Fitness fit, std::size_t const N_total)
      {
        std::vector<float_type> fitness_values(N_total);
        float_type total_fitness(0.0);
        for (std::size_t i(0); i < N_total; ++i)
        {
          fitness_values[i] = fit(i);
          total_fitness += fitness_values[i];
        }
        for (std::size_t i(0); i < N_total; ++i)
        {
          fitness_values[i] /= total_fitness;
        }
        std::uniform_real_distribution<float_type> distribution(0.0, 1.0);
        float_type const target_value(distribution(rng));
        float_type accumulated_value(0.0);
        for (std::size_t j(0); j < N_total; ++j)
        {
          accumulated_value += fitness_values[j];
          if (accumulated_value >= target_value)
          {
            return j;
          }
        }
        return 0;
      }
      std::vector<std::size_t> operator()(Fitness fit, std::size_t const N_total, std::size_t const N_select)
      {
        std::vector<float_type> fitness_values(N_total);
        float_type total_fitness(0.0);
        for (std::size_t i(0); i < N_total; ++i)
        {
          fitness_values[i] = fit(i);
          total_fitness += fitness_values[i];
        }
        for (std::size_t i(0); i < N_total; ++i)
        {
          fitness_values[i] /= total_fitness;
        }
        std::vector<bool> has_been_selected(N_total, false);
        std::uniform_real_distribution<float_type> distribution(0.0, 1.0);
        std::vector<std::size_t> selected;
        for (std::size_t i(0); i < N_select; ++i)
        {
          float_type const target_value(distribution(rng));
          float_type accumulated_value(0.0);
          for (std::size_t j(0); j < N_total; ++j)
          {
            accumulated_value += fitness_values[j];
            if (accumulated_value >= target_value)
            {
              selected.push_back(j);
              break;
            }
          }
        }
        return selected;
      }
    };

  }

  namespace mutators
  {
    template<class T>
    struct point_mutator
    { };
    template<>
    struct point_mutator<bool>
    {
      bool operator() (bool const &,
        bool const & = false,
        bool const & = true)
      {
        std::mt19937_64 mut_rng = std::mt19937_64(std::random_device()());
        std::uniform_int_distribution<std::size_t> m_dist(0, 1);
        return m_dist(mut_rng) == 0;
      }
    };
    template<>
    struct point_mutator<float>
    {
      float operator() (float const &,
        float const & lower = std::numeric_limits<float>::min(),
        float const & upper = std::numeric_limits<float>::max())
      {
        std::mt19937_64 mut_rng = std::mt19937_64(std::random_device()());
        std::uniform_real_distribution<float> m_dist(lower, upper);
        return m_dist(mut_rng);
      }
    };
  }

  template<class ValueType, class Mutator = mutators::point_mutator<ValueType>, class Float = double>
  class individuum
  {
    std::vector<ValueType> m_values;
    std::mt19937_64 m_rng;
    std::uniform_int_distribution<std::size_t> m_dist;
    Mutator m_mutator;
    Float m_health;
  public:
    typedef Float float_type;
    typedef ValueType value_type;
    individuum(std::vector<ValueType> const &values = std::vector<ValueType>(), Mutator const & mutator = Mutator())
      : m_values(values), m_rng(std::random_device()()),
      m_dist(0u, values.size()),
      m_mutator(mutator), m_health(Float())
    { }
    individuum(std::vector<ValueType> && values, Mutator const & mutator = Mutator())
      : m_values(std::forward<std::vector<ValueType>>(values)),
      m_rng(std::random_device()()),
      m_dist(0u, m_values.size() - 1u),
      m_mutator(mutator), m_health(Float())
    { }
    void mutate(ValueType const & lower_bound, ValueType const & upper_bound)
    {
      if (m_dist.max() < m_values.size())
      {
        std::size_t const i = m_dist(m_rng);
        m_values[i] = m_mutator(m_values[i], lower_bound, upper_bound);
      }
    }
    individuum<ValueType, Mutator, Float> mate(individuum<ValueType, Mutator, Float> const & rhs)
    {
      individuum<ValueType, Mutator, Float> ret(*this);
      std::uniform_int_distribution<unsigned int> lrs(0, 1);
      std::size_t const n = m_values.size();
      for (std::size_t i = 0; i < n; ++i)
      {
        if (lrs(m_rng) < 1) ret.m_values[i] = rhs.m_values[i];
      }
      return ret;
    }
    void cross()
    {
      std::size_t const n = m_values.size();
      std::size_t const cross_size = m_dist(m_rng), cross_a = m_dist(m_rng), cross_b = m_dist(m_rng);
      std::size_t end_a = std::min(cross_a + cross_size, n), end_b = std::min(cross_b + cross_size, n);
      std::size_t len = 0u;
      if (cross_a < cross_b)
      {
        end_a = std::min(end_a, cross_b);
        len = std::min(end_a - cross_a, end_b - cross_b);
      }
      else
      {
        end_b = std::min(end_b, cross_a);
        len = std::min(end_a - cross_a, end_b - cross_b);
      }
      for (std::size_t i = 0u; i < len; ++i)
      {
        typename std::vector<ValueType>::value_type tmp = std::move(m_values[cross_a + i]);
        m_values[cross_a + i] = std::move(m_values[cross_b + i]);
        m_values[cross_b + i] = std::move(tmp);
      }
    }
    std::vector<ValueType> const & values() const
    {
      return m_values;
    }
    float_type get_health() const
    {
      return m_health;
    }
    void set_health(float_type const value)
    {
      m_health = value;
    }
  };

  template<class V, class M, class F>
  inline std::ostream& operator<< (std::ostream & stream, individuum<V, M, F> const &i)
  {
    stream << i.get_health() << "; ";
    for (auto const & v : i.values()) stream << v << " ";
    return stream;
  }

  template<class V, class M, class F>
  inline bool operator< (individuum<V, M, F> const & a, individuum<V, M, F> const & b)
  {
    return a.get_health() < b.get_health();
  }

  template<
    class Individual,
    class Updater,
    class Fitness = ranking::linear<typename Individual::float_type>,
    class Selection = selection::Roulette<Fitness>
  >
  class evolution
  {
  public:
    typedef typename Individual::float_type float_type;
  private:
    std::vector<Individual> m_parents, m_children;
    std::mt19937_64 m_rng;
    std::uniform_real_distribution<float_type> m_dist;
    Selection m_selection;
    Fitness m_fitness;
    Updater m_update;
    evolution& operator= (evolution const &);
    bool draw(float_type const chance)
    {
      using std::max;
      using std::min;
      float_type const c = max(float_type(0.), min(chance, float_type(1.)));
      return m_dist(m_rng) <= c;
    }
  public:

    struct options
    {
      typename Individual::value_type mutation_lower_bound,
        mutation_upper_bound;
      float_type mutation_chance_point, crossing_chance,
        fitness_lower_bound, fitness_upper_bound;
      std::size_t fitness_included_individuals;
      options() :
        mutation_lower_bound(), mutation_upper_bound(),
        mutation_chance_point(0.25), crossing_chance(0.15),
        fitness_lower_bound(0.0), fitness_upper_bound(1.0),
        fitness_included_individuals(0)
      { }
    } const & m_options;

    typedef evolution<Individual, Updater, Fitness, Selection> this_type;

    evolution(std::vector<Individual> const & init_population,
      Updater const & updater = Updater(),
      options const & option_object = options()) :
      m_parents(init_population), m_children(init_population.size()),
      m_rng(std::random_device()()), m_dist(0.0, 1.0),
      m_selection(),
      m_fitness(option_object.fitness_included_individuals == 0
        ? m_parents.size() : option_object.fitness_included_individuals,
        option_object.fitness_lower_bound,
        option_object.fitness_upper_bound),
      m_update(updater),
      m_options(option_object)
    { }

    void propagate(std::size_t const steps)
    {
      using std::sort;
      std::size_t const n = m_parents.size();
      sort(m_parents.begin(), m_parents.end());
      for (std::size_t s(0u); s < steps; ++s)
      {

        //std::cout << "Parents:\n";
        //for (auto const & v : m_parents) std::cout << "(" << v << "), ";
        //std::cout << "\n";
        for (std::size_t i(0u); i < n; ++i)
        {
          std::vector<std::size_t> selected(m_selection(m_fitness, n, 2));
          //std::cout << "New child " << i + 1 << " from " << selected[0] << " and " << selected[1] << "\n";

          m_children[i] = m_parents[selected[0]].mate(m_parents[selected[1]]);
          //m_children[i] = m_parents[selected[0]];
          while (draw(m_options.mutation_chance_point))
          {
            m_children[i].mutate(m_options.mutation_lower_bound, m_options.mutation_upper_bound);
          }
          while (draw(m_options.crossing_chance))
          {
            m_children[i].cross();
          }
          // If updating current child fails -> select new parent
          //std::cout << "Updating\n";
          unsigned x = 0;
          while (!m_update(m_children[i], i))
          {
            std::vector<std::size_t> new_selected(m_selection(m_fitness, n, 2));
            //std::cout << "... Failed, Chaning to: " << new_selected[0] << ", " << new_selected[1] <<  " \n";
            m_children[i] = m_parents[new_selected[0]].mate(m_parents[new_selected[1]]);
            if (x > 1000)
            {
              throw std::runtime_error("Couldn't generate child in evolutionary propagation.");
            }
            ++x;
          }
        }
        // Sort children
        sort(m_children.begin(), m_children.end());
        // Swap best parent with worst child
        std::swap(m_parents.front(), m_children.back());
        // Resort
        sort(m_children.begin(), m_children.end());
        m_parents.swap(m_children);
        std::cout << "Mean population health: " << mean_health() << "\n";
      }
    }

    std::vector<Individual> const & population() const { return m_parents; }
    float_type mean_health() const
    {
      float_type ret = float_type();
      for (auto const & p : m_parents)
      {
        ret += p.get_health();
      }
      return ret / static_cast<float_type>(m_parents.size());
    }
  };


  //************************************
  // Method:    evolve
  // FullName:  genetic::evolve
  // Access:    public 
  // Returns:   Vector of individuals after
  //            generation_count generations of processing
  // Qualifier: -
  // Template:  Ind (Type of individual)
  // Template:  Fitness (Type of fitness functor)
  // Template:  Mutator (Type of mutation functor)
  // Template:  Mateor (Type of mating functor)
  // Parameter: std::vector<Ind> parents (initial population)
  // Parameter: std::size_t generation_count (number of generations)
  // Parameter: Fitness && fit (functor returning fitness values for mate functor)
  //            Paramter: std::vector<Ind>&
  // Parameter: Mateor && mate (functor mating parents
  //            Paramter: std::vector<Ind> const & (parent population)
  //            Paramter : ? (fitness values)
  // Parameter: Mutator && mutate (mutation functor)
  //            Paramter: Ind & (individuum to be mutated)
  //************************************

  template<class Ind, class Fitness, class Mutator, class Mateor>
  std::vector<Ind> evolve(std::vector<Ind> population, std::size_t n,
    Fitness && fit, Mateor && mate, Mutator && mutate)
  {
    // propagate for 'n' generations
    for (std::size_t gen = 0; gen < n; ++gen)
    {
      auto const g = gen;
      // build new population with population fitness
      population = mate(population, fit(population, g), g);
      // mutate population
      mutate(population, g);
    }
    // apply fitness function once again
    fit(population, n);
    return population;
  }

}
