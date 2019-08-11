#ifndef cast_genetic_optimization_header

#define cast_genetic_optimization_header

#include <cstddef>
#include <cmath>
#include <vector>
#include <random>

namespace optimization
{

	namespace genetic
	{

		template<class F>
		struct index_value_pair
		{
			using float_type = F;
			std::size_t index;
			F f;
			bool skip;
			bool operator< (index_value_pair const& r) const
			{
				return f < r.f;
			}
		};

		template<class F>
		class Fitness_Linear
		{
			std::size_t N;
			F inclination, offset;
		public:
			using float_type = F;
			Fitness_Linear(std::size_t const N_ind, F const lower_limit = 0., F const upper_limit = 1.)
				: N(N_ind), inclination(-(upper_limit - lower_limit) / (N_ind - 1)), offset(upper_limit)
			{ }
			F operator() (std::size_t const x) const
			{
				return x <= N ? inclination * x + offset : F(0);
			}
		};

		template<class F>
		class Fitness_Exponential
		{
			std::size_t N;
			F b, a, c;
		public:
			using float_type = F;
			Fitness_Exponential(std::size_t const N_ind, F const lower_limit = 0., F const upper_limit = 1.)
				: N(N_ind), b((lower_limit / upper_limit - std::exp(-static_cast<F>(N)))
					/ (F(1) - std::exp(-static_cast<F>(N)))),
				a(F(1) - b), c(upper_limit)
			{ }
			F operator() (std::size_t const x) const
			{
				return x < N ? (a * std::exp(-static_cast<F>(x)) + b) * c : 0.0;
			}
		};

		template<class Fitness, class F = typename Fitness::float_type>
		class Rank_Roulette_Selection
		{
		private:
			std::mt19937_64 rng;
			Fitness fit;
		public:
			using float_type = F;
			Rank_Roulette_Selection(Fitness fitness)
				: rng(std::random_device()()), fit(fitness) { }
			std::vector<std::size_t> operator()(std::size_t const N_total, std::size_t const N_select)
			{
				std::vector<F> fitness_values(N_total);
				F total_fitness(F(0));
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
				std::uniform_real_distribution<F> distribution(0.0, 1.0);
				std::vector<std::size_t> selected;
				for (std::size_t i(0); i < N_select; ++i)
				{
					F const target_value(distribution(rng));
					F accumulated_value(0.0);
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

#endif
