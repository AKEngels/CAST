#ifndef scon_c3_header_include_guard_

#define scon_c3_header_include_guard_

#include <type_traits>
#include <iostream>
#include <tuple>
#include <limits>

#include "scon_traits.h"

namespace scon
{

  namespace _c3
  {

    // template type _fc
    // T for fundamentals
    // T const & for all other types
    template<class Y>
    using _fc = typename std::conditional<
      std::is_fundamental<Y>::value,
      Y, Y const &>::type;

    template<class Y>
    using F = typename scon::float_helper<Y>::type;

  }

  template<class T = double>
  class c3
  {

  protected:

    T m_x, m_y, m_z;

  public:

    using this_type = c3 < T > ;
    using type = T;

    /*!
    *
    *  Constructors
    *
    */

    c3() : m_x(), m_y(), m_z() { }
    
    explicit constexpr c3(_c3::_fc<T> const v) :
      m_x(v), m_y(v), m_z(v) { }

    explicit constexpr c3(_c3::_fc<T> const x, _c3::_fc<T> const y, _c3::_fc<T> const z) :
      m_x(x), m_y(y), m_z(z)
    { }

    template <class U>
    explicit c3(c3<U> const &v) : 
      m_x(static_cast<T>(v.x())), 
      m_y(static_cast<T>(v.y())),
      m_z(static_cast<T>(v.z())) { }

    /*!
    *
    *  Elements by references
    *
    */

    _c3::_fc<T> x() const { return m_x; }
    _c3::_fc<T> y() const { return m_y; }
    _c3::_fc<T> z() const { return m_z; }

    T & x() { return m_x; }
    T & y() { return m_y; }
    T & z() { return m_z; }

    /*!
    *
    *  math ops
    *
    */
    
    this_type& operator+= (this_type const & v)
    {
      m_x += v.m_x;
      m_y += v.m_y;
      m_z += v.m_z;
      return *this;
    }
    this_type& operator-= (this_type const & v)
    {
      m_x -= v.m_x;
      m_y -= v.m_y;
      m_z -= v.m_z;
      return *this;
    }
    this_type& operator*= (this_type const & v)
    {
      m_x *= v.m_x;
      m_y *= v.m_y;
      m_z *= v.m_z;
      return *this;
    }
    this_type& operator/= (this_type const & v)
    {
      m_x /= v.m_x;
      m_y /= v.m_y;
      m_z /= v.m_z;
      return *this;
    }

    this_type& operator+= (_c3::_fc<T> v)
    {
      m_x += v;
      m_y += v;
      m_z += v;
      return *this;
    }
    this_type& operator-= (_c3::_fc<T> v)
    {
      m_x -= v;
      m_y -= v;
      m_z -= v;
      return *this;
    }
    this_type& operator*= (_c3::_fc<T> v)
    {
      m_x *= v;
      m_y *= v;
      m_z *= v;
      return *this;
    }
    this_type& operator/= (_c3::_fc<T> v)
    {
      m_x /= v;
      m_y /= v;
      m_z /= v;
      return *this;
    }

    std::ostream & write(std::ostream &stream) const
    {
      return stream
        .write(reinterpret_cast<char const *>(&m_x), sizeof(T))
        .write(reinterpret_cast<char const *>(&m_y), sizeof(T))
        .write(reinterpret_cast<char const *>(&m_z), sizeof(T));
    }

    std::istream & read(std::istream & stream)
    {
      stream.read(reinterpret_cast<char*>(&m_x), sizeof(T));
      stream.read(reinterpret_cast<char*>(&m_y), sizeof(T));
      stream.read(reinterpret_cast<char*>(&m_z), sizeof(T));
      return stream;
    }

  };

  /*
  
   ######   #######  ##     ## ########     ###    ########  ########
  ##    ## ##     ## ###   ### ##     ##   ## ##   ##     ## ##
  ##       ##     ## #### #### ##     ##  ##   ##  ##     ## ##
  ##       ##     ## ## ### ## ########  ##     ## ########  ######
  ##       ##     ## ##     ## ##        ######### ##   ##   ##
  ##    ## ##     ## ##     ## ##        ##     ## ##    ##  ##
   ######   #######  ##     ## ##        ##     ## ##     ## ########
  
  */

  template<class T>
  inline bool operator== (c3<T> const &l, c3<T> const &r)
  {
    return l.x() == r.x() && l.y() == r.y() && l.z() == r.z();
  }
  template<class T>
  inline bool operator!= (c3<T> const &l, c3<T> const &r)
  {
    return !(l == r);
  }

  template<class T>
  inline bool operator< (c3<T> const &l, c3<T> const &r)
  {
    using std::tie;
    return tie(l.x(), l.y(), l.z()) <
      tie(r.x(), r.y(), r.z());
  }

  template<class T>
  inline bool operator> (c3<T> const &l, c3<T> const &r)
  {
    return r < l;
  }
  template<class T>
  inline bool operator<= (c3<T> const &l, c3<T> const &r)
  {
    return !(r < l);
  }
  template<class T>
  inline bool operator>= (c3<T> const &l, c3<T> const &r)
  {
    return !(l < r);
  }

  /*
  
  ##     ##    ###    ######## ##     ##
  ###   ###   ## ##      ##    ##     ##
  #### ####  ##   ##     ##    ##     ##
  ## ### ## ##     ##    ##    #########
  ##     ## #########    ##    ##     ##
  ##     ## ##     ##    ##    ##     ##
  ##     ## ##     ##    ##    ##     ##

  */

  template<class T>
  c3<T> operator- (c3<T> const &a)
  {
    return c3<T>(-a.x(), -a.y(), -a.z());
  }

  template<class T>
  c3<T> operator+ (c3<T> a, c3<T> const &b)
  {
    a += b;
    return a;
  }
  template<class T>
  c3<T> operator- (c3<T> a, c3<T> const &b)
  {
    a -= b;
    return a;
  }
  template<class T>
  c3<T> operator* (c3<T> a, c3<T> const &b)
  {
    a *= b;
    return a;
  }
  template<class T>
  c3<T> operator/ (c3<T> a, c3<T> const &b)
  {
    a /= b;
    return a;
  }

  template<class T>
  c3<T> operator+ (c3<T> a, _c3::_fc<T> b)
  {
    a += b;
    return a;
  }
  template<class T>
  c3<T> operator- (c3<T> a, _c3::_fc<T> b)
  {
    a -= b;
    return a;
  }
  template<class T>
  c3<T> operator* (c3<T> a, _c3::_fc<T> b)
  {
    a *= b;
    return a;
  }
  template<class T>
  c3<T> operator/ (c3<T> a, _c3::_fc<T> b)
  {
    a /= b;
    return a;
  }

  /*

  ##     ## ######## #### ##       #### ######## ##    ## 
  ##     ##    ##     ##  ##        ##     ##     ##  ##  
  ##     ##    ##     ##  ##        ##     ##      ####   
  ##     ##    ##     ##  ##        ##     ##       ##    
  ##     ##    ##     ##  ##        ##     ##       ##    
  ##     ##    ##     ##  ##        ##     ##       ##    
   #######     ##    #### ######## ####    ##       ##    

  */

  template<class T>
  T min(c3<T> const &v)
  {
    using std::min;
    return min(v.x(), min(v.y(), v.z()));
  }

  template<class T>
  c3<T> max(c3<T> const &v)
  {
    using std::max;
    return max(v.x(), max(v.y(), v.z()));
  }

  template<class T>
  c3<T> min(c3<T> const &a, c3<T> const &b)
  {
    using std::min;
    return c3<T>(min(a.x(), b.x()), min(a.y(), b.y()), min(a.z(), b.z()));
  }

  template<class T>
  c3<T> max(c3<T> const &a, c3<T> const &b)
  {
    using std::max;
    return c3<T>(max(a.x(), b.x()), max(a.y(), b.y()), max(a.z(), b.z()));
  }

  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value, c3<T>>::type abs(c3<T> const &v)
  {
    using std::abs;
    return c3<T>(abs(v.x()), abs(v.y()), abs(v.z()));
  }

  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value, c3<T>>::type round(c3<T> const &v)
  {
    using std::round;
    return c3<T>(round(v.x()), round(v.y()), round(v.z()));
  }

  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value, c3<T>>::type floor(c3<T> const &v)
  {
    using std::floor;
    return c3<T>(floor(v.x()), floor(v.y()), floor(v.z()));
  }

  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value, c3<T>>::type ceil(c3<T> const &v)
  {
    using std::ceil;
    return c3<T>(ceil(v.x()), ceil(v.y()), ceil(v.z()));
  }


  /*

   ######  ######## ########  ########    ###    ##     ##  ######
  ##    ##    ##    ##     ## ##         ## ##   ###   ### ##    ##
  ##          ##    ##     ## ##        ##   ##  #### #### ##
   ######     ##    ########  ######   ##     ## ## ### ##  ######
        ##    ##    ##   ##   ##       ######### ##     ##       ##
  ##    ##    ##    ##    ##  ##       ##     ## ##     ## ##    ##
   ######     ##    ##     ## ######## ##     ## ##     ##  ######

  */

  namespace _c3
  {
    inline int _c3_delim_stream_index()
    {
      static int i = std::ios_base::xalloc();
      return i;
    }

    struct _c3_delimeter_holder
    {
      char d;
      _c3_delimeter_holder(char c) : d(c) { }
    };

    inline std::ostream& operator<< (std::ostream &s,
      _c3_delimeter_holder v)
    {
      s.iword(_c3_delim_stream_index()) = v.d;
      return s;
    }

    template<class S>
    char _c3_get_delimeter(S & strm)
    {
      auto delim = strm.iword(_c3_delim_stream_index());
      if (delim == 0 ||
        delim < std::numeric_limits<char>::min() ||
        delim > std::numeric_limits<char>::max()) return '\0';
      return static_cast<char>(delim);
    }
  }

  inline _c3::_c3_delimeter_holder c3_delimeter(char const c)
  {
    return _c3::_c3_delimeter_holder(c);
  }

  template<typename T>
  inline std::ostream& operator<< (std::ostream &s, c3<T> const &v)
  {
    std::streamsize const i = s.width();
    auto const del = _c3::_c3_get_delimeter(s);
    bool p = s && s << std::setw(i) && s << v.x();
    if (del != '\0') p = p && s << del;
    p = p && s << std::setw(i) && s << v.y();
    if (del != '\0') p = p && s << del;
    p = p && s << std::setw(i) && s << v.z();
    return s;
  }

  template<typename T>
  inline std::istream& operator>> (std::istream &s, c3<T> &v)
  {
    if ((s >> v.x()) && (s >> v.y())) s >> v.z();
    return s;
  }

  template<class T = double> using m3 = c3<c3<T>>;

}

#endif // scon_c3_header_include_guard_