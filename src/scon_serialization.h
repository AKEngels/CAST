#if !defined(SCON_SERIALIZATION_HEADER)

#define SCON_SERIALIZATION_HEADER

//#include "scon_traits.h"
//#include "scon_utility.h"

#include <type_traits>
#include <iostream>
#include <cstddef>
#include <iterator>
#include <vector>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <cstring>
#include <utility>

#if __GNUG__ && __GNUC__ < 5
#define IS_TRIVIALLY_COPYABLE(T) std::has_trivial_copy_constructor<T>::value
#else
#define IS_TRIVIALLY_COPYABLE(T) std::is_trivially_copyable<T>::value
#endif

namespace scon
{

  template<class T>
  using is_trivially_serializable = std::integral_constant<bool,
    std::is_fundamental<T>::value || IS_TRIVIALLY_COPYABLE(T)/* && !scon::is_range<T>::value*/>;

  template<class Target>
  struct binary_stream { };

  template<>
  struct binary_stream<std::ostream>
  {
    std::ostream & strm;
    binary_stream(std::ostream &init_stream) : strm(init_stream) {}
    binary_stream& operator=(binary_stream const&) = delete;
    binary_stream<std::ostream>& write(char const * p, std::size_t const n)
    {
      strm.write(p, n);
      return *this;
    }
    operator bool() const { return static_cast<bool>(strm); }
  };

  template<>
  struct binary_stream<std::istream>
  {
    std::istream & strm;
    binary_stream(std::istream &init_stream) : strm(init_stream) {}
    binary_stream& operator=(binary_stream const&) = delete;
    binary_stream<std::istream>& read(char * const p, std::size_t const n)
    {
      strm.read(p, n);
      return *this;
    }
    operator bool() const { return static_cast<bool>(strm); }
  };

  template<>
  struct binary_stream<std::vector<char>>
  {
    std::vector<char> v;
    std::size_t x;
    bool ok;
    binary_stream() : v(), x(), ok(true) {}
    binary_stream(std::size_t const n) : v(), x(), ok(true) { v.reserve(n); }
    operator bool() const { return ok; }
    binary_stream<std::vector<char>>&
      write(char const * const p, std::size_t const n)
    {
      auto const s = v.size();
      v.resize(s + n);
      std::memcpy(v.data() + s, p, n);
      return *this;
    }
    binary_stream<std::vector<char>>&
      read(char * const p, std::size_t const n)
    {
      auto end = x + n;
      if (v.size() >= end)
      {
        std::memcpy(p, v.data() + x, n);
        x += n;
      }
      else
      {
        ok = false;
      }
      return *this;
    }
  };

  template<>
  struct binary_stream<std::size_t>
  {
    std::size_t N;
    binary_stream() : N() {}
    binary_stream& write(char const * const, std::size_t const n)
    {
      N += n;
      return *this;
    }
  };

  template<class T>
  using bstream = binary_stream<T>;

  using bostream = binary_stream<std::ostream>;
  using bistream = binary_stream<std::istream>;
  using bvstream = binary_stream< std::vector<char> >;

  using bstream_size = binary_stream<std::size_t>;

  template<class StreamBase, class Target>
  typename std::enable_if<is_trivially_serializable<Target>::value, binary_stream<StreamBase> &>::type
    operator<< (binary_stream<StreamBase> &bst, Target const &t)
  {
    auto p = reinterpret_cast<char const *>(std::addressof(t));
    return bst.write(p, sizeof(t));
  }

  template<class StreamBase, class Target>
  typename std::enable_if<is_trivially_serializable<Target>::value, binary_stream<StreamBase> &>::type
    operator<< (binary_stream<StreamBase> &bst, std::vector<Target> const &t)
  {
    if (!t.empty())
    {
      auto p = reinterpret_cast<char const *>(t.data());
      return bst.write(p, t.size()*sizeof(Target));
    }
    return bst;
  }


  template<class StreamBase, class Target>
  typename std::enable_if<is_trivially_serializable<Target>::value && 
    !std::is_const<Target>::value, binary_stream<StreamBase> &>::type
    operator>> (binary_stream<StreamBase> &bst, Target &t)
  {
    auto p = reinterpret_cast<char *>(std::addressof(t));
    return bst.read(p, sizeof(t));
  }

  template<class StreamBase, class Target>
  typename std::enable_if<is_trivially_serializable<Target>::value, binary_stream<StreamBase> &>::type
    operator>> (binary_stream<StreamBase> &bst, std::vector<Target> &t)
  {
    if (!t.empty())
    {
      auto p = reinterpret_cast<char *>(t.data());
      return bst.read(p, t.size()*sizeof(t));
    }
    return bst;
  }

  template<class T>
  std::size_t binary_size(T const &x)
  {
    binary_stream<std::size_t> size_stream;
    size_stream << x;
    return size_stream.N;
  }

  template<class T, class ... U>
  std::size_t binary_size(T const &x, U const & ... us)
  {
    return binary_size(x) + binary_size<U...>(us...);
  }

  //template<class T, bool ir = scon::is_range<T>::value>
  //struct _range_serializer { };

  //template<class T>
  //struct _range_serializer<T, true>
  //{

  //  template<class StreamBase, class Target>
  //  static typename std::enable_if<is_trivially_serializable<range_value<Target>>::value,
  //    binary_stream<StreamBase> &>::type
  //    to(binary_stream<StreamBase> &bst, Target const &t)
  //  {
  //    for (auto const & e : t)
  //    {
  //      bst << e;
  //    }
  //    return bst;
  //  }

  //  template<class StreamBase, class Target>
  //  static typename std::enable_if<is_range<range_value<Target>>::value,
  //    binary_stream<StreamBase> &>::type
  //    to(binary_stream<StreamBase> &bst, Target const &t)
  //  {
  //    for (auto const & e : t)
  //    {
  //      _range_serializer<range_value<Target>>::to(bst, e);
  //    }
  //    return bst;
  //  }

  //  template<class StreamBase, class Target>
  //  static typename std::enable_if<is_trivially_serializable<range_value<Target>>::value,
  //    binary_stream<StreamBase> &>::type
  //    from(binary_stream<StreamBase> &bst, Target &t)
  //  {
  //    for (auto & e : t)
  //    {
  //      bst >> e;
  //    }
  //    return bst;
  //  }

  //  template<class StreamBase, class Target>
  //  static typename std::enable_if<is_range<range_value<Target>>::value,
  //    binary_stream<StreamBase> &>::type
  //    to(binary_stream<StreamBase> &bst, Target &t)
  //  {
  //    for (auto & e : t)
  //    {
  //      _range_serializer<range_value<Target>>::from(bst, e);
  //    }
  //    return bst;
  //  }

  //};

  //template<class StreamBase, class Target>
  //typename std::enable_if<is_range<Target>::value, binary_stream<StreamBase> &>::type
  //  operator<< (binary_stream<StreamBase> &bst, Target const &t)
  //{
  //  return _range_serializer<Target>::to(bst, t);
  //}

  //template<class StreamBase, class Target>
  //typename std::enable_if<is_range<Target>::value, binary_stream<StreamBase> &>::type
  //  operator>> (binary_stream<StreamBase> &bst, Target &t)
  //{
  //  return _range_serializer<Target>::from(bst, t);
  //}


  template<class>
  struct sfinae_true : std::true_type {};

  namespace _detail {
    template<class Stream, class T>
    static auto test_out_op(int) -> 
      sfinae_true<decltype(std::declval<Stream>() << std::declval<T>())>;

    template<class , class T>
    static auto test_out_op(long) -> std::false_type;

    template<class Stream, class T>
    static auto test_in_op(int) ->
      sfinae_true<decltype(std::declval<Stream>() >> std::declval<typename std::remove_reference<T>::type &>())>;

    template<class, class T>
    static auto test_in_op(long)->std::false_type;
  } // detail::

  template<class Stream, class T>
  using has_stream_out = decltype(_detail::test_out_op<Stream, T>(0));

  template<class Stream, class T>
  using has_stream_in = decltype(_detail::test_out_op<Stream, T>(0));

  template<class S, class T>
  using has_stream_ops = std::integral_constant<bool, has_stream_out<S, T>::value && has_stream_in<S, T>::value>;

}


#endif