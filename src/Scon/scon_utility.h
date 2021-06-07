#pragma once
#if !defined(SCON_UTILITY_HEADER)

#define SCON_UTILITY_HEADER

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <ctime>
#include <cstdint>
#include <cstddef>
#include <locale>
#include <algorithm>
#include <climits>
#include <iterator>
#include <limits>
#include <string>
#include <typeinfo>
#include <memory>
#include <type_traits>
#include <functional>
#include <thread>
#include <stdexcept>
#include <initializer_list>
#include "scon.h"
#include "scon_traits.h"
#include "scon_iterator.h"
#if !defined(_MSC_VER)
#include <cxxabi.h>
#endif

#if defined(_MSC_VER) && !defined(thread_local)  // Visual studio
#define thread_local __declspec( thread )
#elif defined(__GCC__) && (__GNUC__ < 4 || __GNUC_MINOR__ < 8) && !defined(thread_local)  // GCC
#define thread_local __thread
#elif  !defined(thread_local)
#define thread_local 
#endif

#if defined(_MSC_VER)
#include "../win_inc.h"
#endif

#ifdef _MSC_VER
#pragma warning (disable: 4996)
#endif

namespace scon
{

  namespace util
  {

#if defined(ULLONG_MAX)
    typedef unsigned long long maxi_uint_type;
#else
    typedef unsigned long maxi_uint_type;
#endif

    template<maxi_uint_type NUM, maxi_uint_type DENUM>
    struct ratio
    {
      typedef maxi_uint_type type;
      static const type num = NUM;
      static const type den = DENUM;
    };

    typedef ratio<1000000000000, 1> tera;
    typedef ratio<1000000000, 1> giga;
    typedef ratio<1000000, 1> mega;
    typedef ratio<1000, 1> kilo;
    typedef ratio<100, 1> hecto;
    typedef ratio<10, 1> deca;
    typedef ratio<1, 10> deci;
    typedef ratio<1, 100> centi;
    typedef ratio<1, 1000> milli;
    typedef ratio<1, 1000000> micro;
    typedef ratio<1, 1000000000> nano;
    typedef ratio<1, 1000000000000> pico;

  }



  template<class T>
  inline T clip_to_circular_range(T const& value, T const& lower, T const& upper)
  {
    using std::ceil;
    if (value < lower)
    {
      T const d = upper - lower;
      return value + ceil((lower - value) / d) * d;
    }
    else if (value > upper)
    {
      T const d = upper - lower;
      return value - ceil((value - upper) / d) * d;
    }
    return value;
  }

  // string stuff

  template <class T>
  std::string type_name()
  {
    typedef typename std::remove_reference<T>::type TR;
    std::unique_ptr<char, void(*)(void*)> own
    (
#ifndef _MSC_VER
      abi::__cxa_demangle(typeid(TR).name(), nullptr,
        nullptr, nullptr),
#else
      nullptr,
#endif
      std::free
    );
    std::string r = own != nullptr ? own.get() : typeid(TR).name();
    if (std::is_const<TR>::value)
      r += " const";
    if (std::is_volatile<TR>::value)
      r += " volatile";
    if (std::is_lvalue_reference<T>::value)
      r += " &";
    else if (std::is_rvalue_reference<T>::value)
      r += " &&";
    return r;
  }

  template<class charT>
  struct ci_str_equal
  {
    ci_str_equal(std::locale const& loc) : loc_(loc) { }
    bool operator()(charT ch1, charT ch2)
    {
      return std::toupper(ch1, loc_) == std::toupper(ch2, loc_);
    }
  private:
    ci_str_equal& operator= (ci_str_equal const&);
    const std::locale& loc_;
  };

  // find substring (case insensitive)
  template<class charT, class traitsT, class allocT,
    class U = std::basic_string<charT, traitsT, allocT>>
    typename std::basic_string<charT, traitsT, allocT>::size_type
    find_substr_ci(std::basic_string<charT, traitsT, allocT> const& str1,
      U const& needle, const std::locale& loc = std::locale())
  {
    std::basic_string<charT, traitsT, allocT> str2(needle);
    using std::begin;
    using std::end;
    auto it = std::search(begin(str1), end(str1),
      begin(str2), end(str2), ci_str_equal<charT>(loc));
    if (it != str1.end()) return static_cast<std::size_t>(it - str1.begin());
    else return std::basic_string<charT, traitsT, allocT>::npos; // not found
  }

  template<class T>
  T str_replace(T s, T const& toReplace, T const& replaceWith)
  {
    auto const f = s.find(toReplace);
    if (f != T::npos)
    {
      s.replace(f, toReplace.length(), replaceWith);
    }
    return s;
  }

  template<class T>
  T str_ireplace(T s, T const& toReplace, T const& replaceWith)
  {
    auto const f = find_substr_ci(s, toReplace);
    if (f != T::npos)
    {
      s.replace(f, toReplace.length(), replaceWith);
    }
    return s;
  }

  template<class T>
  struct sstream_of_char
  {
    typedef std::stringstream ss_type;
  };
  template<>
  struct sstream_of_char < wchar_t >
  {
    typedef std::wstringstream ss_type;
  };

  /*
   * Helperclass for handling filepaths
   *
   * Takes a string containg the full path as input
   * and can then return either only the filename,
   * only the base_path etc.
   *
   * Usually one will use the handy typedef
   * StringFilePath.
   */
  template<class T = std::string>
  struct FilePath
  {
    FilePath(T const& p) : path(p) {}
    T path_no_base(T const& delims = "/\\") const
    {
      typename T::size_type const p(path.find_last_of(delims));
      return p != T::npos ? path.substr(0, path.find_last_of(delims)) : T();
    }
    T base_name(T const& delims = "/\\") const
    {
      return path.substr(path.find_last_of(delims) + 1);
    }
    T extension(T const& delims = ".") const
    {
      T const base = base_name();
      return base.substr(base.find_last_of(delims) + 1);
    }
    T base_no_extension(T const& slash = "/\\", T const& dot = ".") const
    {
      typename T::size_type const
        s(path.find_last_of(slash)),
        p(path.find_last_of(dot));
      if (p != T::npos && (p > (s + 1U)) && ((p - s - 1U) > 0U))
      {
        return path.substr(s + 1, p - s - 1U);
      }
      return path;
    }
    T get_unique_path() const
    {
      T file(path);
      std::size_t next(0U);
      auto const path_base_no_extension = path_no_base() + base_no_extension();
      auto const ext = extension();
      while (std::ifstream(file.c_str(), std::ios_base::in).is_open() && next < 100)
      {
        typename sstream_of_char<typename T::value_type>::ss_type newfss;
        newfss << path_base_no_extension << "_" << ++next << "." << ext;
        file = newfss.str();
      }
      return file;
    }

    // The actual path
    T path;
  };

  typedef FilePath<std::string> StringFilePath;

  template<class It>
  inline std::size_t size_2d(It beg, It const end)
  {
    using std::begin;
    using std::end;
    std::ptrdiff_t n = std::ptrdiff_t();
    for (; beg != end; ++beg)
    {
      n += std::distance(begin(*beg), end(*beg));
    }
    return n >= 0 ? static_cast<std::size_t>(n) : 0u;
  }


  template<class C>
  std::size_t size_2d(C const& object)
  {
    using std::begin;
    using std::end;
    return size_2d(begin(object), end(object));
  }

  template <class T>
  inline std::size_t num_digits(T number)
  {
    std::size_t digits(0);
    while (number) {
      number /= 10;
      ++digits;
    }
    return digits;
  }

  template<class U, class T>
  void static_transform(T const& v, U& t)
  {
    using std::begin;
    using std::end;
    using std::transform;
#if defined(SCON_DEBUG)
    using std::distance;
    if (distance(begin(v), end(v)) != distance(begin(t), end(t)))
    {
      throw std::logic_error("static_transform requires euqally sized ranges.");
    }
#endif
    transform(begin(v), end(v), begin(t), [](scon::range_value<T> const& x)
      -> scon::range_value < U >
    {
      return static_cast<scon::range_value<U>>(x);
    });
  }

  template<class U, class T>
  void explicit_transform(T const& v, U& t)
  {
    using std::begin;
    using std::end;
    using std::transform;
#if defined(SCON_DEBUG)
    using std::distance;
    if (distance(begin(v), end(v)) != distance(begin(t), end(t)))
    {
      throw std::logic_error("static_transform requires euqally sized ranges.");
    }
#endif
    transform(begin(v), end(v), begin(t), [](scon::range_value<T> const& x)
      -> scon::range_value < U >
    {
      return scon::range_value<U>(x);
    });
  }

  template<class U, class T>
  U static_transform(T const& v)
  {
    using std::distance;
    using std::begin;
    using std::end;
    using std::transform;
    auto b = begin(v);
    auto const e = end(v);
    U t(std::distance(b, e));
    transform(b, e, begin(t), [](scon::range_value<T> const& x)
      -> scon::range_value<U>
    {
      return static_cast<scon::range_value<U>>(x);
    });
    return t;
  }

  template<class U, class T>
  U explicit_transform(T const& v)
  {
    using std::begin;
    using std::end;
    using std::transform;
    using std::distance;
    auto b = begin(v);
    auto const e = end(v);
    U t(std::distance(b, e));
    transform(b, e, begin(t), [](scon::range_value<T> const& x)
      -> scon::range_value < U >
    {
      return scon::range_value<U>(x);
    });
    return t;
  }

  namespace transformation
  {

    template<class Output, class Container, class Transformation>
    Output create(Container const& v, Transformation t)
    {
      using std::distance;
      using std::begin;
      using std::end;
      using std::transform;
      auto b = begin(v);
      auto e = end(v);
      Output x(distance(begin(v), end(v)));
      transform(b, e, begin(x), t);
      return x;
    }

  }

  template<class ConT>
  ConT concat(ConT const& a, ConT const& b)
  {
    using std::distance;
    using std::begin;
    using std::end;
    using std::copy;
    auto ba = begin(a);
    auto ea = end(a);
    auto bb = begin(b);
    auto eb = end(b);
    ConT r(distance(ba, ea) + distance(bb, eb));
    copy(bb, eb, copy(ba, ea, begin(r)));
    return r;
  }

  template<class ConT, class T>
  ConT pushed(ConT const& a, T const& b)
  {
    ConT r(a);
    r.push_back(b);
    return r;
  }

  namespace sorted
  {

    namespace sorted_detail
    {
      template<class It, class T>
      inline bool _exists_it(It first, It last, T const& value)
      {
        auto found = std::lower_bound(first, last, value);
        return found != last && !(value < *found);
      }

      template<class It, class T, class Compare>
      inline bool _exists_it(It first, It last, T const& value, Compare comp)
      {
        auto found = std::lower_bound(first, last, value, comp);
        return found != last && !comp(value, *found);
      }
    }

    template<class Container, class T, class Compare>
    std::size_t find(Container& object, T const& value, Compare comp)
    {
      using std::begin;
      using std::end;
      auto const first = begin(object);
      auto const last = end(object);
      auto const found = std::lower_bound(first, last, value, comp);
      return static_cast<std::size_t>(
        (found != last && !comp(value, *found)) ?
        std::distance(first, found) : std::distance(first, last)
        );
    }

    template<class Container, class T>
    std::size_t find(Container& object, T const& value)
    {
      using std::begin;
      using std::end;
      auto const first = begin(object);
      auto const last = end(object);
      auto const found = std::lower_bound(first, last, value);
      return static_cast<std::size_t>(
        (found != last && !(value < *found)) ?
        std::distance(first, found) : std::distance(first, last)
        );
    }

    template<class C, class T, class Compare>
    bool exists(C const& object, T const& value, Compare comp)
    {
      using std::begin;
      using std::end;
      return sorted_detail::_exists_it(begin(object), end(object), value, comp);
    }

    template<class C, class T>
    bool exists(C const& object, T const& value)
    {
      using std::begin;
      using std::end;
      return sorted_detail::_exists_it(begin(object), end(object), value);
    }

    template<class Container, class T>
    void insert(Container& object, T const& value)
    {
      using std::begin;
      using std::end;
      object.insert(std::lower_bound(begin(object),
        end(object), value), value);
    }

    template<class Container, class T, class Compare>
    void insert(Container& object, T const& value, Compare comp)
    {
      using std::begin;
      using std::end;
      object.insert(std::lower_bound(begin(object),
        end(object), value, comp), value);
    }

    template<class Container, class T>
    void insert_unique(Container& object, T const& value)
    {
      using std::begin;
      using std::end;
      auto first = begin(object);
      auto last = end(object);
      auto found = std::lower_bound(first, last, value);
      if (found == last || value < *found) object.insert(found, value);
    }

    template<class Container, class T, class Compare>
    void insert_unique(Container& object, T const& value, Compare comp)
    {
      using std::begin;
      using std::end;
      auto first = begin(object);
      auto last = end(object);
      auto found = std::lower_bound(first, last, value, comp);
      if (found == last || comp(value, *found)) object.insert(found, value);
    }

    template<class C, class T>
    void insert_noresize(C& object, T const& value)
    {
      using std::begin;
      using std::end;
      auto first = begin(object);
      auto last = end(object);
      auto const found = std::lower_bound(first, last, value);
      if (found != last)
      {
        for (auto i = last - 1; i != found; --i)
        {
          *i = *(i - 1);
        }
        *found = value;
      }
    }

    template<class Container, class T, class Compare>
    void insert_noresize(Container& object, T const& value, Compare&& comp)
    {
      using std::begin;
      using std::end;
      auto first = begin(object);
      auto last = end(object);
      auto const found = std::lower_bound(first, last, value, std::forward<Compare>(comp));
      if (found != last)
      {
        for (auto i = last - 1; i != found; --i)
        {
          *i = *(i - 1);
        }
        *found = value;
      }
    }


    template<class C, class T>
    void insert_unique_noresize(C& object, T const& value)
    {
      using std::begin;
      using std::end;
      auto first = begin(object);
      auto last = end(object);
      auto const found = std::lower_bound(first, last, value);
      if (found != last && value < *found)
      {
        for (auto i = last - 1; i != found; --i)
        {
          *i = *(i - 1);
        }
        *found = value;
      }
    }

    template<class C, class T, class Compare>
    void insert_unique_noresize(C& object, T const& value, Compare comp)
    {
      using std::begin;
      using std::end;
      auto first = begin(object);
      auto last = end(object);
      auto const found = std::lower_bound(first, last, value, comp);
      if (found != last && comp(value, *found))
      {
        for (auto i = last - 1; i != found; --i)
        {
          *i = *(i - 1);
        }
        *found = value;
      }
    }

    template<class C, class T>
    void remove_equal(C& object, T const& value)
    {
      using std::begin;
      using std::end;
      auto first = begin(object);
      auto last = end(object);
      auto er = std::equal_range(first, last, value);
      auto result = first;
      for (; result != er.first;)
      {
        ++result;
      }
      for (auto m_it = er.second; m_it != last; ++m_it)
      {
        *(result++) = *m_it;
      }
      object.resize(std::distance(first, result));
    }

  }


#if defined(SCON_CC11_RANDOM)

  struct random
  {
#if defined(SCON_CC11_THREAD)

    template<class T>
    inline std::mt19937_64 threaded_mt_engine(T x)
    {
      static std::mutex mtx;
      static std::random_device r;
      static auto mt = std::mt19937_64{
        std::seed_seq{ r(), r(), r(), r(), r(), r() } };
      std::unique_lock<std::mutex> lock(mtx);
      std::seed_seq seq{ mt(), mt(), mt(), mt(), x };
      return std::mt19937_64{ seq };
    }

    template<typename D>
    static typename D::result_type rand(D dist)
    {
      static auto lrng = std::mt19937_64{ std::random_device{}() };
      return dist(lrng);
    }
    template<class Distribution>
    static typename Distribution::result_type
      threaded_rand(Distribution dist)
    {
      static thread_local std::mt19937_64 tlrng =
        std::mt19937_64(std::random_device{}() +
          std::hash<std::thread::id>()(std::this_thread::get_id()));
      return dist(tlrng);
    }
#endif
  };

  template<class T, bool F = std::is_floating_point<T>::value>
  struct uniform_ditribution_trait
  {
    using _t = typename std::enable_if<
      std::is_arithmetic<T>::value, T>::type;
    using type = std::uniform_int_distribution < _t >;
  };

  template<class T>
  struct uniform_ditribution_trait < T, true >
  {
    using type = std::uniform_real_distribution < T >;
  };

  template<class T>
  struct uniform_distribution
  {
    using distribution_type =
      typename uniform_ditribution_trait<T>::type;
    using result_type = typename distribution_type::result_type;
    using param_type = typename distribution_type::param_type;
    distribution_type dist;
    uniform_distribution() : dist() { }
    uniform_distribution(T low, T high) : dist(low, high) { }
    template<class G>
    auto operator() (G& generator)
      -> decltype(dist(generator))
    {
      return dist(generator);
    }
    template<class G>
    auto operator() (G& generator, param_type const& p)
      -> decltype(dist(generator, p))
    {
      return dist(generator, p);
    }
  };

  template<class T>
  struct uniform_randomizer
  {
    uniform_distribution<T> d;
    uniform_randomizer() : d() { }
    uniform_randomizer(T const l)
      : d(-l, l)
    { }
    uniform_randomizer(T const l, T const h)
      : d(l, h) { }
    void operator ()
      (typename uniform_distribution<T>::result_type& v)
    {
      v = scon::random::rand(d);
    }
    typename uniform_distribution<T>::result_type operator()()
    {
      return scon::random::rand(d);
    }
  };

  template<typename T>
  inline T rand()
  {
    return random::rand(std::uniform_int_distribution<T>());
  }
  template<>
  inline float rand<float>()
  {
    return random::rand(std::uniform_real_distribution<float>());
  }
  template<>
  inline double rand<double>()
  {
    return random::rand(std::uniform_real_distribution<double>());
  }
  template<>
  inline long double rand<long double>()
  {
    return random::rand(std::uniform_real_distribution<long double>());
  }
  template<typename T>
  inline T rand(T const low, T const up)
  {
    return random::rand(std::uniform_int_distribution<T>(low, up));
  }
  template<>
  inline double rand<double>(double const low, double const up)
  {
    return random::rand(std::uniform_real_distribution<double>(low, up));
  }
  template<>
  inline float rand<float>(float const low, float const up)
  {
    return random::rand(std::uniform_real_distribution<float>(low, up));
  }
  template<>
  inline long double rand<long double>(long double const low, long double const up)
  {
    return random::rand(std::uniform_real_distribution<long double>(low, up));
  }
#if defined(SCON_CC11_THREAD)
  //template<typename T>
  //inline typename std::uniform_int_distribution<T>::result_type rand_thread(T const & low, T const & up)
  //{
  //  return random::threaded_random<std::uniform_int_distribution<T>>(std::uniform_int_distribution<T>(low, up));
  //}
  //template<>
  //inline std::uniform_real_distribution<double>::result_type rand_thread<double>(double const & low, double const & up)
  //{
  //  return random::threaded_random<std::uniform_real_distribution<double>>(std::uniform_real_distribution<double>(low, up));
  //}
  //template<>
  //inline std::uniform_real_distribution<float>::result_type rand_thread<float>(float const & low, float const & up)
  //{
  //  return random::threaded_random<std::uniform_real_distribution<float>>(std::uniform_real_distribution<float>(low, up));
  //}
#endif
#else

  template<typename T>
  inline T rand(T const& low, T const& up)
  {
    static const T rm(static_cast<T>(RAND_MAX));
    return low + static_cast<T>(rand()) / rm * (up - low);
  }

#endif

  template<class T, class U = T>
  inline typename std::enable_if<std::is_arithmetic<U>::value>::type
    rand(T & val, U const& low, U const& up)
  {
    val = rand<U>(low, up);
  }

  template<class T>
  inline typename std::enable_if<std::is_arithmetic<T>::value>::type
    randomize(T& val, T const& low, T const& up)
  {
    val = rand<T>(low, up);
  }

  template<class T>
  inline typename std::enable_if<std::is_arithmetic<T>::value>::type
    randomize(T& ref)
  {
    ref = scon::rand<T>();
  }

  template<class T>
  inline T randomized()
  {
    T tmp;
    randomize(tmp);
    return tmp;
  }

  template<class T, class ... CArgs>
  inline T randomized(CArgs&& ... arg)
  {
    T tmp(arg...);
    randomize(tmp);
    return tmp;
  }

  template<class T>
  void limit(T& value, T const& minimum, T const& maximum)
  {
    using std::max;
    using std::min;
    value = max(minimum, min(maximum, value));
  }



  template<class T>
  struct Randomizer
  {
    void operator () (T& v)
    {
      randomize(v);
    }
  };



  template<class T>
  void resize(std::size_t const N, T& v)
  {
    v.resize(N);
  }

  template<class T, class... U>
  void resize(std::size_t const N, T& a, U& ... others)
  {
    a.resize(N);
    resize(N, others...);
  }

  template<class T>
  void clear(T& v)
  {
    v.clear();
  }

  template<class T, class... U>
  void clear(T& a, U& ... others)
  {
    a.clear();
    clear(others...);
  }

  template<class T, class U>
  bool equal_size_ranges(T const& x, U const& y)
  {
    using std::distance;
    using std::begin;
    using std::end;
    return distance(begin(x), end(x)) == distance(begin(y), end(y));
  }

  template<class T, class U, class... V>
  bool equal_size_ranges(T const& x, U const& y, V const& ... others)
  {
    using std::distance;
    using std::begin;
    using std::end;
    return equal_size_ranges(x, y) && equal_size_ranges(y, others...);
  }

  template<class T>
  bool empty(T& v)
  {
    return v.empty();
  }

  template<class T, class... U>
  bool empty(T& a, U& ... others)
  {
    return a.empty() && empty(others...);
  }

  template<template<class...> class T, class ... TArgs>
  T<TArgs...> make(TArgs&& ... args)
  {
    return T<TArgs...>(std::forward<TArgs>(args)...);
  }

  template<class T>
  struct range_printer
  {
    T const* p;
    std::string d;
    range_printer(T const& t, std::string delimeter = " ") : p(&t), d(delimeter) {}
  };
  template<class T>
  inline std::ostream& operator<< (std::ostream& o, range_printer<T> const& rp)
  {
    bool first = true;
    for (auto const& elem : *(rp.p))
    {
      if (!first) o << rp.d;
      o << elem;
      first = false;
    }
    return o;
  }

  template<class T>
  range_printer<T> print_range(T const& t, std::string const delim = " ")
  {
    return range_printer<T>(t, delim);
  }



  namespace delimeted_tuple_detail
  {

    template<class Ch, class ChTr, class ... T>
    struct delimeted_tuple
    {
      std::tuple<T...> v;
      std::basic_string<Ch, ChTr> delim;
      template<class ... U>
      delimeted_tuple(std::basic_string<Ch, ChTr> const& delimeter, U&& ... us) :
        v(std::forward<U>(us)...), delim(delimeter) {}

    };

    template<std::size_t N, class ... Ts>
    struct TuplePrinter
    {
      template<class Ch, class ChTr>
      static void print(std::basic_ostream<Ch, ChTr>& strm,
        delimeted_tuple<Ch, ChTr, Ts...> const& t)
      {
        TuplePrinter<N - 1, Ts...>::print(strm, t);
        strm << std::get<N - 1>(t.v);
        if (N < sizeof...(Ts)) strm << t.delim;
      }
    };

    template<std::size_t N>
    struct TuplePrinter<N>
    {
      template<class Ch, class ChTr>
      static void print(std::basic_ostream<Ch, ChTr>& strm,
        delimeted_tuple<Ch, ChTr> const& t)
      { }
    };

    template<class ... Ts>
    struct TuplePrinter<1, Ts...>
    {
      template<class Ch, class ChTr>
      static void print(std::basic_ostream<Ch, ChTr>& strm,
        delimeted_tuple<Ch, ChTr, Ts...> const& t)
      {
        strm << std::get<0u>(t.v);
        if (sizeof...(Ts) > 1) strm << t.delim;
      }
    };

    template<class Ch, class ChTr, class ... T>
    std::basic_ostream<Ch, ChTr>& operator<< (std::basic_ostream<Ch, ChTr>& strm,
      delimeted_tuple<Ch, ChTr, T...> const& t)
    {
      TuplePrinter<sizeof...(T), T...>::template print<Ch, ChTr>(strm, t);
      return strm;
    }

  }

  template<class ... T>
  delimeted_tuple_detail::delimeted_tuple<std::string::value_type, std::string::traits_type, T...>
    delimeted(std::string const& delimeter, T&& ... ts)
  {
    return delimeted_tuple_detail::delimeted_tuple<std::string::value_type,
      std::string::traits_type, T...>(delimeter, std::forward<T>(ts)...);
  }

  namespace reverse_range_detail
  {
    using std::begin;
    using std::end;
    template<class T>
    struct reverted_iter_range
    {
      using begin_iterator = std::reverse_iterator<decltype(end(std::declval<T&>()))>;
      using end_iterator = std::reverse_iterator<decltype(begin(std::declval<T&>()))>;
      begin_iterator _begin;
      end_iterator _end;
      reverted_iter_range(T& r) : _begin(end(r)), _end(begin(r)) {}
    };
    template<class T>
    typename reverted_iter_range<T>::begin_iterator begin(reverted_iter_range<T> const& r)
    {
      return r._begin;
    }
    template<class T>
    typename reverted_iter_range<T>::end_iterator end(reverted_iter_range<T> const& r)
    {
      return r._end;
    }
  }

  template<class T>
  inline reverse_range_detail::reverted_iter_range<T> reverse_range(T& range)
  {
    return{ range };
  }

  //template<class T>
  //inline reverse_range_detail::reverted_iter_range<T> reverse_range(T && range)
  //{
  //  throw std::logic_error("rvalues not allow to revert range.");
  //}

  namespace range_index_detail
  {
    using std::begin;
    using std::end;
    struct range_indexer
    {
      std::size_t f, l;
      template<class R>
      range_indexer(R& r) : f(0u),
        l(static_cast<std::size_t>(std::distance(begin(r), end(r)))) {}
      range_indexer(std::size_t const start_index,
        std::size_t const end_index) : f(start_index), l(end_index)
      { }
    };
    inline scon::index_iterator begin(range_indexer const& ri) { return{ ri.f }; }
    inline scon::index_iterator end(range_indexer const& ri) { return{ ri.l }; }
  }

  template<class T>
  inline range_index_detail::range_indexer index_range(T& range)
  {
    return{ range };
  }

  inline range_index_detail::range_indexer index_range(
    std::size_t const first, std::size_t const last)
  {
    return{ first, last };
  }

  namespace detail
  {
    // termination function
    template<class T> void apply_unary(T&& value) { }
    // recursive function for applying all unarys
    template<class T, class F1, class... UnaryFs>
    void apply_unary(T&& value, F1&& f, UnaryFs&& ... fs)
    {
      // using function call syntax
      std::forward<F1>(f)(std::forward<T>(value));
      // using c++17 invoke
      //std::invoke(std::forward<F1>(f), std::forward<T>(value));
      // recurse and forward value and rest
      apply_unary(std::forward<T>(value), std::forward<UnaryFs>(fs)...);
    }
  }

  template<class Iter, class ... UnaryFs>
  std::tuple<UnaryFs...> for_each(Iter first, Iter last, UnaryFs&& ... fs)
  {
    for (; first != last; ++first)
    {
      // forward callable pack to apply
      detail::apply_unary(*first, std::forward<UnaryFs>(fs)...);
    }
    // return callable pack
    return std::forward_as_tuple(std::forward<UnaryFs>(fs)...);
  }

  /*! Function to call other programms
  *
  */

  int system_call(std::string const& command_line);


  /*! Function to seperate a string in letters, numbers and other
  *
  */


  enum charTypeT { other, alpha, digit };

  charTypeT charTypestring(char);

  std::string separateString(std::string);

  /**
   * Checks whether the ranges (begin1,  end1) and (begin2, end2) match, i.e. whether they are the same length and contain the same elements.
   */
  template<typename InputIt1, typename InputIt2>
  bool match(InputIt1 begin1, InputIt1 end1, InputIt2 begin2, InputIt2 end2) {
    if (std::distance(begin1, end1) != std::distance(begin2, end2)) return false;

    return std::all_of(begin1, end1, [begin2, end2](auto elem){
      return std::find(begin2, end2, elem) != end2;
    });
  }
}
#endif

