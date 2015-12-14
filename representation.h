#ifndef representation_header

#define representation_header

#include <type_traits>
#include <functional>
#include <utility>

namespace optimization
{

  struct empty_log
  {
    void operator() (...) { }
  };

  struct bool_false_functor
  {
    bool operator() (...) { return false; }
  };

  template<class R, class G = R>
  struct State
  {

    using type = R;
    using gradient_type = G;

    type x;
    gradient_type g;

    State(type const & init_x = type(), gradient_type const &init_g = gradient_type())
      : x(init_x), g(init_g)
    { }

    State(type && init_x, gradient_type const &init_g = gradient_type())
      : x(std::move(init_x)), g(init_g)
    { }

    State(type && init_x, gradient_type && init_g)
      : x(std::move(init_x)),
      g(std::move(init_g))
    { }

  };

  template<class R, class G>
  inline bool operator== (State<R, G> const & lhs,
    State<R, G> const & rhs)
  {
    return lhs.x == rhs.x;
  }

  template<class R, class G>
  inline bool operator!= (State<R, G> const & lhs,
    State<R, G> const & rhs)
  {
    return !(lhs == rhs);
  }

  template<class R, class G, class F>
  struct Point : State<R, G>
  {

    using type = R;
    using gradient_type = G;
    using float_type = F;
    using state_type = State<R, G>;

    float_type f;

    // With Base

    Point(state_type const & init_s = state_type(), F const init_f = F())
      : state_type(init_s), f(init_f)
    { }

    Point(state_type && init_s, F const init_f = F())
      : state_type(std::move<state_type>(init_s)), f(init_f)
    { }

    // With Base members

    Point(type const & init_x, 
      gradient_type const &init_g = gradient_type(), 
      F const init_f = F())
      : state_type(init_x, init_g), f(init_f)
    { }

    Point(type && init_x,
      gradient_type const &init_g = gradient_type(),
      F const init_f = F())
      : state_type(std::move(init_x), init_g), f(init_f)
    { }

    Point(type && init_x,
      gradient_type &&init_g,
      F const init_f = F())
      : state_type(std::move(init_x),
      std::move<gradient_type>(init_g)),
      f(init_f)
    { }

  };

  template<class R, class G, class F>
  inline bool operator== (Point<R, G, F> const & lhs,
    Point<R, G, F> const & rhs)
  {
    return static_cast<State<R, G> const &>(lhs) 
      == static_cast<State<R, G> const &>(rhs);
  }
  template<class R, class G, class F>
  inline bool operator!= (Point<R, G, F> const & lhs,
    Point<R, G, F> const & rhs)
  {
    return !(lhs == rhs);
  }

  template<class R, class G, class F>
  inline bool operator< (Point<R, G, F> const & lhs,
    F const rhs)
  {
    return lhs.f < rhs;
  }

  template<class R, class G, class F>
  inline bool operator< (Point<R, G, F> const & lhs,
    Point<R, G, F> const & rhs)
  {
    return lhs.f < rhs.f;
  }
  template<class R, class G, class F>
  inline bool operator> (Point<R, G, F> const & lhs,
    Point<R, G, F> const & rhs)
  {
    return rhs < lhs;
  }
  template<class R, class G, class F>
  inline bool operator<= (Point<R, G, F> const & lhs,
    Point<R, G, F> const & rhs)
  {
    return !(rhs < lhs);
  }
  template<class R, class G, class F>
  inline bool operator>= (Point<R, G, F> const & lhs,
    Point<R, G, F> const & rhs)
  {
    return !(lhs < rhs);
  }

  template<class StateT, class F>
  using State_Point = Point <typename StateT::type, 
    typename StateT::gradient_type, F>;

}

#endif
