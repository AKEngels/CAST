#if !defined(SCON_ANGLE_HEADER)
#define SCON_ANGLE_HEADER

#include <cmath>
#include <iostream>

#include "scon.h"
#include "scon_traits.h"
#include "scon_utility.h"

#ifndef SCON_PI_DEFINES
  #define SCON_PI_DEFINES
  #define SCON_PI     3.14159265358979323846264338327950288419716939937510
  #define SCON_PI_INV 0.31830988618379067153776752674502872406891929148091
  #define SCON_2PI    6.28318530717958647692528676655900576839433879875021
  #define SCON_PI180  0.01745329251994329576923690768488612713442871888541
  #define SCON_180PI 57.29577951308232087679815481410517033240547246656432
#endif

namespace scon
{

  template<class T>
  class ang
  {
    typename std::enable_if<
      std::is_floating_point<T>::value,
      T>::type m_value;

    static T to_range(T const value)
    {
      using std::fmod;
      T const m(fmod(value, T(SCON_2PI)));
      return to_range_s(m);
    }

    static T to_range_s(T const value)
    {
      return value <= -T(SCON_PI) ? value + T(SCON_2PI) :
        (value > T(SCON_PI) ? value - T(SCON_2PI) : value);
    }

  public:

    typedef T value_type;
    typedef ang<T> this_type;

    static this_type from_rad(T const value)
    {
      return this_type(value);
    }
    static this_type from_deg(T const value)
    {
      return this_type(value*T(SCON_PI180));
    }

    ang()
      : m_value()
    { }

    explicit ang(T const value)
      : m_value(to_range(value))
    { }

    explicit operator T () const { return m_value; }

    T degrees() const { return m_value*T(SCON_180PI); }
    T radians() const { return m_value; }

    T full_radians() const { return m_value < T() ? T(SCON_2PI) + m_value : m_value; }
    T full_degrees() const { return full_radians()*T(SCON_180PI); }

    this_type operator- () const
    {
      this_type tmp(*this);
      tmp.m_value = -tmp.m_value;
      return tmp;
    }

    this_type& operator+= (this_type const &rhs)
    {
      m_value = to_range_s(m_value + rhs.m_value);
      return *this;
    }
    this_type& operator-= (this_type const &rhs)
    {
      m_value = to_range_s(m_value - rhs.m_value);
      return *this;
    }
    this_type& operator*= (this_type const &rhs)
    {
      m_value = to_range(m_value * rhs.m_value);
      return *this;
    }
    this_type& operator/= (this_type const &rhs)
    {
      m_value = to_range(m_value / rhs.m_value);
      return *this;
    }
    this_type& operator+= (T rhs)
    {
      m_value = to_range(m_value + rhs);
      return *this;
    }
    this_type& operator-= (T rhs)
    {
      m_value = to_range(m_value + rhs);
      return *this;
    }
    this_type& operator*= (T rhs)
    {
      m_value = to_range(m_value + rhs);
      return *this;
    }
    this_type& operator/= (T rhs)
    {
      m_value = to_range(m_value + rhs);
      return *this;
    }

  };

  //// compare

  template<class T>
  bool operator== (ang<T> const &a, ang<T> const &b)
  {
    return a.radians() == b.radians();
  }
  template<class T>
  bool operator!= (ang<T> const &a, ang<T> const &b)
  {
    return !(a == b);
  }
  template<class T>
  bool operator< (ang<T> const &a, ang<T> const &b)
  {
    return a.radians() < b.radians();
  }
  template<class T>
  bool operator> (ang<T> const &a, ang<T> const &b)
  {
    return b < a;
  }
  template<class T>
  bool operator<= (ang<T> const &a, ang<T> const &b)
  {
    return !(b < a);
  }
  template<class T>
  bool operator>= (ang<T> const &a, ang<T> const &b)
  {
    return !(a < b);
  }

  // math

  template<class T>
  ang<T> operator+ (ang<T> a, ang<T> const &b)
  {
    a += b;
    return a;
  }
  template<class T>
  ang<T> operator- (ang<T> a, ang<T> const &b)
  {
    a -= b;
    return a;
  }
  template<class T>
  ang<T> operator* (ang<T> a, ang<T> const &b)
  {
    a *= b;
    return a;
  }
  template<class T>
  ang<T> operator/ (ang<T> a, ang<T> const &b)
  {
    a /= b;
    return a;
  }

  template<class T>
  ang<T> operator+ (ang<T> a, T b)
  {
    a += b;
    return a;
  }
  template<class T>
  ang<T> operator- (ang<T> a, T b)
  {
    a -= b;
    return a;
  }
  template<class T>
  ang<T> operator* (ang<T> a, T b)
  {
    a *= b;
    return a;
  }
  template<class T>
  ang<T> operator/ (ang<T> a, T b)
  {
    a /= b;
    return a;
  }


  // utility


  template<class T>
  inline T sin(ang<T> const & v)
  {
    using std::sin;
    return sin(v.radians());
  }

  template<class T>
  inline T cos(ang<T> const & v)
  {
    using std::cos;
    return cos(v.radians());
  }

  template<class T>
  inline T tan(ang<T> const &  v)
  {
    using std::tan;
    return tan(v.radians());
  }

  template<class T>
  inline ang<T> abs(ang<T> v)
  {
    if (v.radians() < T()) return -v;
    return v;
  }


  // helpers

  

  template<class T>
  inline void randomize(ang<T> & ref)
  {
    ref = ang<T>::from_rad(rand(-T(SCON_PI), T(SCON_PI)));
  }

  namespace _ang_detail
  {
    inline int _ang_stream_index()
    {
      static int i = std::ios_base::xalloc();
      return i;
    }
    template<class S>
    inline bool _stream_mode_degrees(S & s)
    {
      return s.iword(_ang_stream_index()) == 0;
    }
  }

  inline std::ostream& degrees(std::ostream& os)
  {
    os.iword(_ang_detail::_ang_stream_index()) = 0;
    return os;
  }

  inline std::ostream& radians(std::ostream& os)
  {
    os.iword(_ang_detail::_ang_stream_index()) = 1;
    return os;
  }

  inline std::istream& degrees(std::istream& is)
  {
    is.iword(_ang_detail::_ang_stream_index()) = 0;
    return is;
  }

  inline std::istream& radians(std::istream& is)
  {
    is.iword(_ang_detail::_ang_stream_index()) = 1;
    return is;
  }

  template<class T>
  inline std::ostream& operator<< (std::ostream &strm, ang<T> const &r)
  {
    strm << (_ang_detail::_stream_mode_degrees(strm) ? 
      r.degrees() : r.radians());
    return strm;
  }

  template<class T>
  inline std::istream& operator>> (std::istream &strm, ang<T> &r)
  {
    T v = T();
    if (strm >> v)
    {
      r = _ang_detail::_stream_mode_degrees(strm) ? 
        ang<T>::from_deg(v) : ang<T>::from_rad(v);
    }
    return strm;
  }

}

#endif
