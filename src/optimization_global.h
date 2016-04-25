#pragma once

#include <vector>
#include <random>
#include <string>
#include <memory>
#include "coords.h"
#include "scon_chrono.h"
#include "scon_log.h"

namespace optimization
{



  namespace global
  {

    struct Tabu_Point
    {
      coords::PES_Point pes;
      std::vector<coords::Gradients_Main> main_direction;
      std::vector<coords::Representation_Internal> intern_direction;
      std::vector<coords::Representation_3D> xyz_direction;
      std::size_t visited, iteration;
      Tabu_Point() : visited(0), iteration() {}
      Tabu_Point(coords::PES_Point const &pes_point, std::size_t const iter = 0) 
        : pes(pes_point), visited(0), iteration(iter)
      {}
      void swap(Tabu_Point &rhs)
      {
        pes.swap(rhs.pes);
        main_direction.swap(rhs.main_direction);
        intern_direction.swap(rhs.intern_direction);
        xyz_direction.swap(rhs.xyz_direction);
        std::swap(visited, rhs.visited);
      }
      operator coords::PES_Point () const { return pes; }
    };

    inline bool operator<  (const Tabu_Point& lhs, const Tabu_Point& rhs) { return lhs.pes.energy < rhs.pes.energy; }
    inline bool operator>  (const Tabu_Point& lhs, const Tabu_Point& rhs) { return  operator< (rhs, lhs); }
    inline bool operator<= (const Tabu_Point& lhs, const Tabu_Point& rhs) { return !operator> (lhs, rhs); }
    inline bool operator>= (const Tabu_Point& lhs, const Tabu_Point& rhs) { return !operator< (lhs, rhs); }

    struct Tabu_List
      : public std::vector<Tabu_Point>
    {
      typedef std::vector<Tabu_Point> base_type;
      bool tabu(coords::PES_Point const &, coords::Coordinates const &) const;
      bool has_superposition(coords::PES_Point const &, coords::Coordinates const &) const;
      void clear_above_e(coords::float_type const minimum, coords::float_type const kT);
    };

    static coords::float_type const k = static_cast<coords::float_type>(0.0019872966); // == kB (in kcal/K) * NA == R
    static coords::float_type const R_MAX = static_cast<coords::float_type>(RAND_MAX);

    struct Result_Point
    {
      coords::PES_Point pes;
      std::size_t it;
      Result_Point() : pes(), it() {}
      Result_Point (std::size_t iteration, coords::PES_Point const & init_pes) : 
        pes(init_pes), it(iteration)
      { }
    };

    class result_drain
    {
      coords::Coordinates * cp;
      std::unique_ptr<std::ofstream> energy_stream, structure_stream;
      std::size_t n;
    public:
      result_drain() : cp(), energy_stream(), structure_stream(), n() {}
      result_drain(coords::Coordinates &c, std::string const & filename);
      void operator() (Result_Point && xyz);
    };

    using offset_buffered_resultfile =
      scon::vector_offset_buffered_callable<Result_Point, result_drain>;

    offset_buffered_resultfile make_buffered_result(coords::Coordinates &coord_obj, 
      std::string const & file_suffix, std::size_t buffer_size);

    class optimizer
    {

      optimizer& operator= (optimizer const &);

      void write_accepted(std::string const &suffix = "");

    public:

      struct min_status 
      { enum T { 
        REJECT_ENERGY, 
        REJECT_BROKEN, 
        REJECT_STEREO, 
        REJECT_TABU, 
        ACCEPT_MINIMUM, 
        ACCEPT_GLOBAL_MINIMUM 
      }; };

      coords::float_type T, kT;
      Tabu_List tabulist, range_tabu;
      //Tabu_Point min_current, min_global;
      Tabu_List accepted_minima, range_minima;
      std::vector< std::size_t> /*accepted_iteration,*/ range_iteration;
      std::vector<coords::Stereo::pair> init_stereo;
      std::size_t fallbacks, i, min_index, gmin_index;
      coords::Coordinates &coordobj;
      offset_buffered_resultfile accepted_log;
      scon::chrono::high_resolution_timer opt_clock;
      bool found_new_minimum;

      // Constructors
      optimizer(coords::Coordinates & c, std::string const & output_name = "");
      optimizer(coords::Coordinates & c, 
        coords::Ensemble_PES const & initial_pes_points, 
        bool const minimized = false, std::string const & output_name = "");

      // Optimization function
      virtual bool run(std::size_t const iterations, bool const _reset = false) = 0;
      void write_range(std::string const &suffix = "");
      virtual ~optimizer() { write_accepted(); }

    protected:

      void header_to_cout();
      min_status::T check_pes_of_coords();
      void setTemp(const coords::float_type temp);
      bool accept(const coords::float_type new_energy) const;
      void swap(optimizer &);

      // mechanics to set next starting point
      // possibly restoring older / other minima
      void set_current_startpoint(std::size_t const);
      bool restore(min_status::T const );
      bool restore_last_global();
      bool restore_roulette_selection();

      min_status::T new_minimum(void);
      void updateRange(coords::PES_Point const &p);

      coords::Coordinates rehydrate(coords::Coordinates const & desolved);
      coords::Coordinates dehydrate(coords::Coordinates const & solved);

    };

    std::ostream & operator<< (std::ostream &, optimizer::min_status::T const &);
    std::unique_ptr<optimizer> new_divers_optimizer(coords::Coordinates &);

    namespace optimizers
    {

      class main_grid
        : public optimizer
      {
        coords::angle_type m_delta;
        std::size_t m_num_offsets;
        std::vector<std::size_t> m_offset, m_init_offset;
        bool m_done;
      public:
        main_grid(coords::Coordinates &c, coords::Ensemble_PES const &p, 
          coords::angle_type const delta_grid);
        bool run(std::size_t const iterations, bool const _reset = false);
        coords::Representation_Main next_main_from_offset();
        void increase_offset(std::size_t it);
      };

      class monteCarlo
        : public optimizer
      {
      public:
        // constructor
        monteCarlo(coords::Coordinates &c, std::string const & output_name = "", bool mc_in_neb = false);
		    monteCarlo(coords::Coordinates &c, coords::Ensemble_PES const &p, std::string const & output_name = "", bool mc_in_neb = false);
        // starter function
        bool run(std::size_t const iterations, bool const _reset = false);
        // move selector
        std::size_t move(coords::Coordinates &);
        // cartesian moves
        std::size_t move_xyz(coords::Coordinates &);
        // dihedral moves
        std::size_t move_main(coords::Coordinates &);
        // constrained dihedral move functions
        std::size_t move_main_strain(coords::Coordinates &);
        void bias_main_rot(coords::Coordinates &, std::size_t it, coords::angle_type rot);
        // move water molecules
        std::size_t move_water(coords::Coordinates &);
      private:
        std::vector<coords::Representation_1D> main_dih_rot_tabu;
		    bool neb_true;
      };

      class tabuSearch
        : public optimizer
      {
      public:
        std::unique_ptr<optimizer> divers_optimizer;
        //optimizer * divers_optimizer;
        // ..
        tabuSearch(coords::Coordinates & c, std::string const & output_name = "");
        tabuSearch(coords::Coordinates & c, coords::Ensemble_PES const &p, std::string const & output_name = "");
        // run tabu search
        bool run(std::size_t const iterations, bool const _reset = false);
        // perform ascent and descent
        void ascent(void);
        coords::float_type descent(void);
        //! Start a diversification search if tabu search is stuck
        bool diversification(void);

      };

    }

    namespace genetic_ranking
    {

      struct index_value_pair
      {
        std::size_t index;
        coords::float_type function_value;
        bool skip;
        bool operator< (index_value_pair const & r) const
        {
          return function_value < r.function_value;
        }
      };

      class Fitness_Linear
      {
        Fitness_Linear& operator= (Fitness_Linear const &);
        std::size_t const N;
        coords::float_type const inclination, offset;
      public:
        Fitness_Linear(std::size_t const N_ind, coords::float_type const lower_limit, coords::float_type const upper_limit)
          : N(N_ind), inclination(-(upper_limit - lower_limit) / (N_ind - 1)), offset(upper_limit)
        { }
        coords::float_type operator() (std::size_t const x) const
        {
          return x <= N ? inclination*x + offset : 0.0;
        }
      };

      class Fitness_Exponential
      {
        Fitness_Exponential& operator= (Fitness_Exponential const &);
        std::size_t const N;
        coords::float_type const b, a, c;
      public:
        Fitness_Exponential(std::size_t const N_ind, coords::float_type const L, coords::float_type const H)
          : N(N_ind),
          b((L / H - std::exp(-static_cast<coords::float_type>(N))) / (1 - std::exp(-static_cast<coords::float_type>(N)))),
          a(1-b), c(H)
        {
        }
        coords::float_type operator() (std::size_t const x) const
        {
          using std::exp;
          return x < N ? (a*exp(-static_cast<coords::float_type>(x)) + b)*c : 0.0;
        }
      };

      template<class Fitness>
      class Roulette
      {
      private:
        std::mt19937_64 rng;
      public:
        Roulette() : rng(std::random_device()()) { }
        std::vector<std::size_t> operator()(Fitness fit, std::size_t const N_total, std::size_t const N_select)
        {
          std::vector<coords::float_type> fitness_values(N_total);
          coords::float_type total_fitness(0.0);
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
          std::uniform_real_distribution<coords::float_type> distribution(0.0, 1.0);
          std::vector<std::size_t> selected;
          for (std::size_t i(0); i < N_select; ++i)
          {
            coords::float_type const target_value(distribution(rng));
            coords::float_type accumulated_value(0.0);
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

  }

  // Abstract interface base class, specialization for Forcefields, tinker, terachem etc.

}
