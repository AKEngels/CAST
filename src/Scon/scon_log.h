#ifndef scon_log_header_include_guard_
#define scon_log_header_include_guard_

#include "scon_traits.h"
#include "scon_utility.h"
#include "scon_serialization.h"

#include <cstddef>
#include <vector>
#include <iterator>
#include <utility>
#include <algorithm>
#include <type_traits>
#include <iostream>


namespace scon
{

  namespace _buffer
  {

    template<bool move> struct unary_flush;

    template<>
    struct unary_flush<true>
    {
      template<class F, class C>
      static inline void f(F&& f, C&& c)
      {
        for (auto&& d : c) f(std::move(d));
      }
    };

    template<>
    struct unary_flush<false>
    {
      template<class F, class C>
      static inline void f(F&& f, C&& c)
      {
        for (auto&& d : c) f(d);
      }
    };


    template<class F, class C>
    void unary_flush_container(F&& f, C&& c)
    {
      unary_flush<!std::is_lvalue_reference<
        scon::argument_type<F, 0u>>::value>::f(f, c);
      c.clear();
    }

    template<class F, class C>
    struct flushing_pair
    {

      using pair_type = std::pair<F, C>;

      std::pair<F, C> p;

      // Special members

      flushing_pair(flushing_pair&&) = default;
      flushing_pair(flushing_pair const&) = default;
      flushing_pair& operator=(flushing_pair const&) = default;
      flushing_pair& operator=(flushing_pair&&) = default;


      // Perfect forwarding construction
      SCON_PERFECT_FORWARDING_WRAPPER_CONSTRUCTORS(flushing_pair, pair_type, p)

        // Flush pair before destruction
        ~flushing_pair()
      {
        using ::scon::_buffer::unary_flush_container;
        unary_flush_container(p.first, p.second);
      }
    };
  } // namespace _buffer



  template<class Callable, class Container>
  struct offset_buffered_callable :
    _buffer::flushing_pair<Callable, Container>
  {

    using base_type = _buffer::flushing_pair<Callable, Container>;

    std::size_t mmax, moff, mi;

    template<class ... B>
    offset_buffered_callable(std::size_t const size,
      std::size_t const offset, B&& ... b)
      : base_type(std::forward<B>(b)...),
      mmax(size), moff(offset), mi()
    {
      this->p.second.reserve(size);
    }

    template<class ...U>
    bool operator() (U&& ... values)
    {
      if (moff > 0 && mi % moff == 0)
      {
        if (this->p.second.size() == mmax)
        { // Flush buffer if full
          using _buffer::unary_flush_container;
          unary_flush_container(this->p.first, this->p.second);
        }
        this->p.second.emplace_back(std::forward<U>(values)...);
        ++mi;
        return true;
      }
      ++mi;
      return false;
    }

  };

  template<class StreamBase, class Callable, class Container>
  typename std::enable_if<scon::has_stream_ops<bstream<StreamBase>,
    Callable>::value, bstream<StreamBase>&>::type
    operator<< (bstream<StreamBase>& bst,
      offset_buffered_callable<Callable, Container> const& obc)
  {
    bst << obc.mi << obc.mmax << obc.moff << obc.p.first << obc.p.second.size();
    for (auto const& e : obc.p.second)
    {
      bst << e;
    }
    return bst;
  }

  template<class StreamBase, class Callable, class Container>
  typename std::enable_if<!scon::has_stream_ops<bstream<StreamBase>,
    Callable>::value, bstream<StreamBase>&>::type
    operator<< (bstream<StreamBase>& bst,
      offset_buffered_callable<Callable, Container> const& obc)
  {
    bst << obc.mi << obc.mmax << obc.moff << obc.p.second.size();

    for (auto const& e : obc.p.second)
    {
      bst << e;
    }
    return bst;
  }

  template<class Callable, class Container, class StreamBase>
  typename std::enable_if<scon::has_stream_ops<bstream<StreamBase>,
    Callable>::value, bstream<StreamBase>&>::type operator>> (bstream<StreamBase>& bst,
      offset_buffered_callable<Callable, Container>& obc)
  {
    decltype(obc.p.second.size()) vsize;
    if (bst >> obc.mi && bst >> obc.mmax && bst >> obc.moff && bst >> obc.p.first && bst >> vsize)
    {
      obc.p.second.clear();
      obc.p.second.reserve(obc.mmax);
      obc.p.second.resize(vsize);
      for (auto& e : obc.p.second)
      {
        bst >> e;
      }
    }
    return bst;
  }

  template<class Callable, class Container, class StreamBase>
  typename std::enable_if<!scon::has_stream_ops<bstream<StreamBase>,
    Callable>::value, bstream<StreamBase>&>::type operator>> (bstream<StreamBase>& bst,
      offset_buffered_callable<Callable, Container>& obc)
  {
    decltype(obc.p.second.size()) vsize;
    if (bst >> obc.mi && bst >> obc.mmax && bst >> obc.moff && bst >> vsize)
    {
      obc.p.second.clear();
      obc.p.second.reserve(obc.mmax);
      obc.p.second.resize(vsize);
      for (auto& e : obc.p.second)
      {
        bst >> e;
      }
    }
    return bst;
  }



  template<class T, class Callable>
  using vector_offset_buffered_callable =
    offset_buffered_callable<Callable, std::vector<T>>;

  template<class T, class Callable>
  vector_offset_buffered_callable<T, Callable>
    offset_call_buffer(std::size_t const size, std::size_t const offset, Callable&& cb)
  {
    return vector_offset_buffered_callable<T, typename std::remove_reference<Callable>::type>(
      size, offset, std::forward<Callable>(cb), std::vector<T>());
  }

  template<class Callable, class Container>
  struct contraction_buffered_callable :
    _buffer::flushing_pair<Callable, Container>
  {

    using base_type = _buffer::flushing_pair<Callable, Container>;

    std::size_t mmax, msp, mi;

    template<class ... B>
    contraction_buffered_callable(std::size_t const size,
      std::size_t const offset, B&& ... b) :
      base_type(std::forward<B>(b)...),
      mmax(size % 2 == 0 ? size : size + 1), msp(offset), mi()
    {
      this->p.second.reserve(size);
    }

    template<class ...U>
    bool operator() (U&& ... values)
    {
      // Contact buffer if full
      if (mi % msp == 0 && this->p.second.size() == mmax)
      {
        contract();
      }
      // msp changed, new check required
      if (mi % msp == 0)
      {
        // Add new element
        this->p.second.emplace_back(std::forward<U>(values)...);
        ++mi;
        return true;
      }
      ++mi;
      return false;
    }

  private:

    void contract()
    {
      msp *= 2;
      std::size_t const N = this->p.second.size();
      std::size_t i(1u);
      for (std::size_t j(2u); j < N; j += 2, ++i)
      {
        this->p.second[i] = std::move(this->p.second[j]);
      }
      this->p.second.resize(i);
    }

  };

  template<class T, class Callable>
  using vector_contraction_buffered_callable =
    contraction_buffered_callable<Callable, std::vector<T>>;

  template<class T, class Callable>
  vector_contraction_buffered_callable<T, Callable>
    contracting_call_buffer(std::size_t const size, std::size_t const offset, Callable&& cb)
  {
    return vector_contraction_buffered_callable<T, Callable>(
      size, offset, std::forward<Callable>(cb), std::vector<T>());
  }

}


#endif // scon_log_header_include_guard_