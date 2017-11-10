#if !defined(SCON_SPHERICAL_COORDS_HEADER)
#define SCON_SPHERICAL_COORDS_HEADER

#include<iostream>
#include <type_traits>

#include "scon.h"
#include "scon_angle.h"
#include "scon_c3.h"

namespace scon
{

  template<class T, class U = ang<T>>
  struct sphericals
  {

    T m_radius;
    U m_inclination, m_azimuth;

  public:

    typedef sphericals<T,U>  this_type;
    typedef U              angle_type;

    // construction

    sphericals() : m_radius(), 
      m_inclination(), m_azimuth() { }

    explicit sphericals(T radius, 
      U const inclination, U const azimuth)
      : m_radius(radius), m_inclination(inclination), 
      m_azimuth(azimuth) { }

    //// Conversion from and to c3<T>
  
    //sphericals(c3<T> const &r)
    //  : m_radius(r.x()), m_inclination(r.y()),
    //  m_azimuth(r.z())
    //{ }
    //operator c3<T>()
    //{
    //  return c3<T>(m_radius, 
    //    m_inclination, 
    //    m_azimuth);
    //}

    // accessors

    T & radius() { return m_radius; }
    angle_type & inclination() { return m_inclination; }
    angle_type & azimuth() { return m_azimuth; }
    T radius() const { return m_radius; }
    angle_type const & inclination() const { return m_inclination; }
    angle_type const & azimuth() const { return m_azimuth; }

    // math operators

    this_type& operator+= (this_type const & v)
    {
      m_radius += v.m_radius;
      m_inclination += v.m_inclination;
      m_azimuth += v.m_azimuth;
      return *this;
    }
    this_type& operator-= (this_type const & v)
    {
      m_radius -= v.m_radius;
      m_inclination -= v.m_inclination;
      m_azimuth -= v.m_azimuth;
      return *this;
    }
    this_type& operator*= (this_type const & v)
    {
      m_radius *= v.m_radius;
      m_inclination *= v.m_inclination;
      m_azimuth *= v.m_azimuth;
      return *this;
    }
    this_type& operator/= (this_type const & v)
    {
      m_radius /= v.m_radius;
      m_inclination /= v.m_inclination;
      m_azimuth /= v.m_azimuth;
      return *this;
    }

    this_type& operator+= (T v)
    {
      m_radius += v;
      m_inclination += v;
      m_azimuth += v;
      return *this;
    }

    this_type& operator-= (T v)
    {
      m_radius -= v;
      m_inclination -= v;
      m_azimuth -= v;
      return *this;
    }

    this_type& operator*= (T v)
    {
      m_radius *= v;
      m_inclination *= v;
      m_azimuth *= v;
      return *this;
    }

    this_type& operator/= (T v)
    {
      m_radius /= v;
      m_inclination /= v;
      m_azimuth /= v;
      return *this;
    }

    this_type operator- () const
    {
      return this_type(-m_radius, -m_inclination, -m_azimuth);
    }

  };

  template<class T, class U>
  inline sphericals<T, U> operator+ (sphericals<T, U> a,
    sphericals<T, U> const & b)
  {
    a += b;
    return a;
  }
  template<class T, class U>
  inline sphericals<T, U> operator- (sphericals<T, U> a,
    sphericals<T, U> const & b)
  {
    a -= b;
    return a;
  }
  template<class T, class U>
  inline sphericals<T, U> operator* (sphericals<T, U> a,
    sphericals<T, U> const & b)
  {
    a *= b;
    return a;
  }
  template<class T, class U>
  inline sphericals<T, U> operator/ (sphericals<T, U> a,
    sphericals<T, U> const & b)
  {
    a /= b;
    return a;
  }

  template<class T, class U>
  inline sphericals<T, U> operator+ (sphericals<T, U> a, T b)
  {
    a += b;
    return a;
  }
  template<class T, class U>
  inline sphericals<T> operator- (sphericals<T, U> a, T b)
  {
    a -= b;
    return a;
  }
  template<class T, class U>
  inline sphericals<T, U> operator* (sphericals<T, U> a, T b)
  {
    a *= b;
    return a;
  }
  template<class T, class U>
  inline sphericals<T, U> operator/ (sphericals<T, U> a, T b)
  {
    a /= b;
    return a;
  }

  template<typename T, class U>
  std::ostream& operator<< (std::ostream &s, sphericals<T, U> const &v)
  {
    std::streamsize const i = s.width();
    bool p = s && s << std::setw(i) && s << v.radius();
    auto const del = _c3::_c3_get_delimeter(s);
    if (del != '\0') p = p && s << del;
    p = p && s << std::setw(i) && s << v.inclination();
    if (del != '\0') p = p && s << del;
    p = p && s << std::setw(i) && s << v.azimuth();
    return s;
  }

  template<typename T, class U>
  std::istream& operator>> (std::istream &s, sphericals<T,U> &v)
  {
    T t = T();
    U a = U(), d = U();
    if (s >> t && s >> a && s >> d)
    { v = sphericals<T>(t, a, d); }
    return s;
  }

}

#endif
