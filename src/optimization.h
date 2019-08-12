#ifndef cast_general_optimization_header

#define cast_general_optimization_header

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

#include "Scon/scon_utility.h"
#include "Scon/scon_vect.h"
#include "coords_io.h"
#include "ls.h"
#include "optimization_genetic.h"
#include "representation.h"

namespace optimization {
	/**function to perform a local optimization
	@param coords: coordinates object
	@param ci: ???
	@param lo_structure_fn: filename for output structure (should end in .arc)
	@param lo_energies_fn: filename for energy and rmsd information (might end in
	.txt)*/
	void perform_locopt(coords::Coordinates& coords, coords::input::format& ci,
		std::string& lo_structure_fn, std::string& lo_energies_fn);

	enum status_types {
		REJECT_ENERGY = -10,
		REJECT_BROKEN,
		REJECT_STEREO,
		REJECT_TABU,
		REJECT_ERROR,
		SUCCESS = 0,
		ACCEPT_MINIMUM = 1,
		ACCEPT_GLOBAL_MINIMUM
	};

	template <class R, class G, class F, class D>
	struct Minimum : Point<R, G, F> {

		using type = R;
		using gradient_type = G;
		using float_type = F;
		using point_type = Point<R, G, F>;
		using state_type = State<R, G>;
		using directions_type = std::vector<D>;

		directions_type directions;
		std::size_t number_of_visits;

		// With Base

		Minimum(point_type const& init_p = point_type(),
			directions_type const& init_d = directions_type())
			: point_type(init_p), directions(init_d), number_of_visits() {}

		Minimum(point_type&& init_p,
			directions_type const& init_d = directions_type())
			: point_type(std::move(init_p)), directions(init_d), number_of_visits() {}

		Minimum(point_type&& init_p, directions_type&& init_d)
			: point_type(std::move(init_p)), directions(std::move(init_d)),
			number_of_visits() {}

		// With Base base

		Minimum(state_type const& init_s, F const init_f = F(),
			directions_type const& init_d = directions_type())
			: point_type(init_s, init_f), directions(init_d), number_of_visits() {}

		Minimum(state_type&& init_s, F const init_f = F(),
			directions_type const& init_d = directions_type())
			: point_type(std::move(init_s), init_f), directions(init_d),
			number_of_visits() {}

		Minimum(state_type&& init_s, F const init_f, directions_type&& init_d)
			: point_type(std::move(init_s), init_f), directions(std::move(init_d)),
			number_of_visits() {}

		// With base of base members

		Minimum(type const& init_x, gradient_type const& init_g = gradient_type(),
			F const init_f = F(),
			directions_type const& init_d = directions_type())
			: point_type(init_x, init_g, init_f), directions(init_d),
			number_of_visits() {}

		Minimum(type&& init_x, gradient_type const& init_g = gradient_type(),
			F const init_f = F(),
			directions_type const& init_d = directions_type())
			: point_type(std::move(init_x), init_g, init_f), directions(init_d),
			number_of_visits() {}

		Minimum(type&& init_x, gradient_type&& init_g, F const init_f = F(),
			directions_type const& init_d = directions_type())
			: point_type(std::move(init_x), std::move(init_g), init_f),
			directions(init_d), number_of_visits() {}

		Minimum(type&& init_x, gradient_type&& init_g, F const init_f,
			directions_type&& init_d)
			: point_type(std::move(init_x), std::move(init_g), init_f),
			directions(std::move(init_d)), number_of_visits() {}
	};

	template <class PointT, class DirectionT>
	using Point_Minimum =
		Minimum<typename PointT::type, typename PointT::gradient_type,
		typename PointT::float_type, DirectionT>;

	template <class CallbackT>
	using Callback_State =
		State<function_trait_detail::decayed_argument_type<CallbackT, 0U>,
		function_trait_detail::decayed_argument_type<CallbackT, 1U>>;

	template <class CallbackT>
	using Callback_Point =
		Point<function_trait_detail::decayed_argument_type<CallbackT, 0U>,
		function_trait_detail::decayed_argument_type<CallbackT, 1U>,
		function_trait_detail::return_type<CallbackT>>;

	template <class CallbackT>
	using Callback_Minimum =
		Minimum<function_trait_detail::decayed_argument_type<CallbackT, 0U>,
		function_trait_detail::decayed_argument_type<CallbackT, 1U>,
		function_trait_detail::return_type<CallbackT>,
		function_trait_detail::decayed_argument_type<CallbackT, 0U>>;

	template <typename T>
	inline bool success(T S) {
		if (S == T::SUCCESS || S == T::ACCEPT_MINIMUM ||
			S == T::ACCEPT_GLOBAL_MINIMUM) {
			return true;
		}
		else
			return false;
	}

	template <class F>
	struct constants {
		static F const kB; // kB in kcal/mol/K
	};

	struct binary_false {
		template <class T, class U = T>
		bool operator()(T const&, U const&) const {
			return false;
		}
	};

	template <class R, class G, class F, class SimilarityComparator>
	class UniqueStatesList {

	public:
		using type = R;
		using gradient_type = G;
		using float_type = F;
		using state_type = State<R, G>;
		using point_type = Point<R, G, F>;
		using compare_type = SimilarityComparator;

		compare_type similar;

	private:
		std::vector<point_type> m_list;
		float_type f_margin;

	public:
		UniqueStatesList(compare_type comparator_object,
			float_type const max_diff =
			std::numeric_limits<float_type>::max() / float_type(2),
			float_type const margin = 1.0)
			: similar(comparator_object), m_list(), f_margin(margin) {}

		std::vector<point_type> const& list() const { return m_list; }

		bool operator()(point_type const& s) {
			using std::begin;
			using std::end;
			using std::lower_bound;
			auto const list_begin = begin(m_list);
			auto const list_end = end(m_list);
			auto i = lower_bound(list_begin, list_end, s.f);
			if (i != list_end) {
				if (i == list_begin) {
					// if e is smaller than first f, v is not tabu
					// otherwise if v and x of i are similar, v is tabu
					if (!((s.f + f_margin) < (*i).f) && similar(s.x, (*i).x))
						return false;
				}
				else {
					// if e is smaller than f of i
					if ((s.f + f_margin) < (*i).f) {
						// if e is bigger than f of i-1, v is not tabu
						// otherwise if v and x of i-1 are similar, v is tabu
						if (!((s.f - f_margin) > (*(i - 1)).f) && similar(s.x, (*(i - 1)).x))
							return false;
					}
					else {
						// if e is in the range of f of i
						// and v and x of i are similiar, v is tabu
						if (similar(s.x, (*i).x))
							return false;
						// if v and x of i are not similar
						// and if e is not bigger than f of i-1
						// and if v and x of i-1 are similar, v is tabu
						else if (!((s.f - f_margin) > (*(i - 1)).f) &&
							similar(s.x, (*(i - 1)).x))
							return false;
					}
				}
			}
			else if (list_begin != list_end) {
				auto const f_diff = s.f - (*list_begin).f;
				if (!((s.f - f_margin) > (*(list_end - 1)).f) &&
					similar(s.x, (*(i - 1)).x))
					return false;
			}
			// Insert v into tabu list if not yet tabu
			m_list.insert(i, s);
			// update range
			return true;
		}
	};

	template <class CallbackT, class ComparatorT>
	using Callback_UniqueStatesList = UniqueStatesList<
		function_trait_detail::decayed_argument_type<CallbackT, 0U>,
		function_trait_detail::decayed_argument_type<CallbackT, 1U>,
		function_trait_detail::return_type<CallbackT>, ComparatorT>;

	template <class R, class G, class F, class SimilarityComparatorT>
	class Metropolis_Range_Acceptor {

	public:
		using type = R;
		using gradient_type = G;
		using float_type = F;
		using compare_type = SimilarityComparatorT;

		using point_type = Point<R, G, F>;

		std::vector<point_type> const& tabu() const { return m_tabu.list(); }
		std::vector<point_type> const& range() const { return m_range; }
		std::vector<point_type> const& accepted() const { return m_accepted; }

	private:
		UniqueStatesList<R, G, F, compare_type> m_tabu;
		std::vector<point_type> m_range;
		std::vector<point_type> m_accepted;

		float_type T, T_scaling, e_range, e_last;
		bool use_last, is_first;

		bool update_range() {
			while (!m_range.empty() &&
				(m_range.back().f - m_range.front().f) > e_range) {
				m_range.pop_back();
			}
			return true;
		}

		bool metropolis(float_type const e, float_type const e0) const {
			using scon::rand;
			using std::abs;
			using std::exp;
			// Nan or anything? DO not accept!
			if (e != e)
				return false;
			// More than two orders of magnitude difference
			// with an absolute difference of more than 50 ?
			// Pobably not sane...
			auto const ae = abs(e);
			auto const ae0 = abs(e0);
			if (abs(ae - ae0) > float_type(50) && (ae / ae0) > float_type(100))
				return false;
			// R[0,1) < e^((E0-E)/kT) ?
			return e < e0 || rand(float_type(0), float_type(1)) <
				exp((e0 - e) / (constants<float_type>::kB * T));
		}

	public:
		Metropolis_Range_Acceptor(compare_type compare,
			float_type const temperature = 298.15,
			float_type const range = 20.0,
			float_type const temperature_scaling = 0.95,
			bool logging = false)
			: m_tabu(compare), m_range(binary_false(), range),
			m_accepted(binary_false()), T(temperature),
			T_scaling(temperature_scaling), e_range(range), e_last(),
			is_first(true), use_last(false) {}

		status_types operator()(point_type const& s) {
			if (is_first)
				e_last = s.f + 1.0;
			auto const e0 = (use_last || range().list().empty())
				? e_last
				: range().list().front().f;
			if (is_first || metropolis(s.f, e0)) {
				is_first = false;
				bool const not_in_tabu = tabu(s);
				if (not_in_tabu) {
					e_last = s.f;
					m_accepted.push_back(s);
					scon::sorted::insert(m_range, s);
					update_range();
					T *= T_scaling;
					return status_types::ACCEPT_MINIMUM;
				}
				return status_types::REJECT_TABU;
			}
			return status_types::REJECT_ENERGY;
		}
	};

	template <class CallbackT, class ComparatorT>
	using Callback_Metropolis_Range_Acceptor = Metropolis_Range_Acceptor<
		function_trait_detail::decayed_argument_type<CallbackT, 0U>,
		function_trait_detail::decayed_argument_type<CallbackT, 1U>,
		function_trait_detail::return_type<CallbackT>, ComparatorT>;

	template <class F>
	struct MC_Point_Mover {
		F lim;
		MC_Point_Mover(F const limit) : lim(abs(limit)) {}
		template <class MinimumT, class CallbackT>
		status_types operator()(MinimumT& p, CallbackT const&) {
			using scon::randomize;
			auto m(p.x);
			randomize(m, -lim, lim);
			for (auto const& o : p.directions) {
				m -= o * dot(o, m);
			}
			p.directions.push_back(m);
			p.x += m;
			return status_types::SUCCESS;
		}
	};

	template <class R, class G, class F>
	struct Exponential_Roulette_Point_Selector {

		using type = R;
		using gradient_type = G;
		using float_type = F;
		using point_type = Point<R, G, F>;
		using fitness_type = genetic::Fitness_Exponential<float_type>;
		using selection_type =
			genetic::Rank_Roulette_Selection<fitness_type, float_type>;

		selection_type selection;
		std::size_t included;

		Exponential_Roulette_Point_Selector(std::size_t const included_points,
			float_type const lower_limit = 0.,
			float_type const upper_limit = 1.)
			: selection(fitness_type(included_points, lower_limit, upper_limit)) {}

		point_type const& operator()(scon::vector<point_type> const& v) {
			return v[selection(included, 1U)[0U]];
		}
	};

	template <class CallbackT, class MoverT, class OptimizerT, class AcceptorT,
		class SelectorT>
		class Global_Optimizer {

		public:
			using callback_type = CallbackT;
			using minimum_type = Callback_Minimum<callback_type>;
			using move_type = MoverT;
			using optimizer_type = OptimizerT;
			using acceptor_type = AcceptorT;
			using selector_type = SelectorT;

			using type = typename minimum_type::type;
			using gradient_type = typename minimum_type::gradient_type;
			using float_type = typename minimum_type::float_type;

			using state_type = typename minimum_type::state_type;
			using point_type = typename minimum_type::point_type;

			callback_type callback;
			move_type move;
			optimizer_type optimizer;
			acceptor_type accept;
			selector_type select;

		private:
			std::vector<minimum_type> accepted_points;
			std::size_t log_done, selection_limit, fallback_limit;

		public:
			Global_Optimizer(callback_type callback_object, move_type move_object,
				optimizer_type optimize_object,
				acceptor_type acceptor_object, selector_type selector_object)
				: callback(callback_object), move(move_object),
				optimizer(optimize_object), accept(acceptor_object),
				select(selector_object), accepted_points(), log_done(0u),
				selection_limit(20u), fallback_limit(20u) {}

			bool optimize(minimum_type p0, std::size_t const max_iterations = 1000u) {
				minimum_type p(p0);
				std::size_t pi = 0U;
				optimizer(p, callback);
				if (success(accept(p, 0u))) {
					accepted_points.push_back(p);
				}
				for (std::size_t i(1u); i <= max_iterations; ++i) {
					if (success(move(p, callback)) && success(optimizer(p, callback))) {
						if (success(accept(p, i))) {
							accepted_points.push_back(p);
							p0 = p;
							continue;
						}
					}
					if (!new_p(p))
						return false;
				}
				return true;
			}

			bool new_p(minimum_type& p) {
				minimum_type const* pp = &select(accepted_points);
				if (pp->number_of_visits < fallback_limit) {
					p = *pp;
					return true;
				}
				else {
					for (std::size_t i(0u); i < selection_limit; ++i) {
						pp = &select(accepted_points);
						if (pp->number_of_visits < fallback_limit) {
							p = *pp;
							++pp->number_of_visits;
							return true;
						}
					}
				}
				return false;
			}
	};

} // namespace optimization

#endif // cast_general_optimization_header
