#if !defined(SCON_INTERNAL_HEADER)
#define SCON_INTERNAL_HEADER

#include<iostream>

#include "scon.h"
#include "scon_traits.h"
#include "scon_angle.h"
#include "scon_vect.h"

namespace scon
{

  template<class T>
  class internals
  {
    T m_x;
    angles<T> m_y, m_z;

  public:

    typedef internals<T>     this_type;
    typedef std::size_t          size_type;
    typedef std::ptrdiff_t       difference_type;
    typedef T                floatop_type;
    typedef ang<T> angle_type;
    typedef T                dot_type;
    typedef T                fundamental_type;

    internals()
      : m_x(static_cast<T>(0.)), m_y(T(0.)), m_z(T(0.))
    { }

    /*! copy constructors */

    internals(internals const& v)
      : m_x(v.m_x), m_y(v.m_y), m_z(v.m_z)
    { }

    // Implicit conversion from vect3d<T> as x, deg, deg

    internals(vect3d<T> const& v)
      : m_x(v.m_x),
      m_y(angle_type::from_deg(v.m_y)),
      m_z(angle_type::from_deg(v.m_z))
    { }

    operator vect3d<T>() const
    {
      return vect3d<T>(m_x, m_y.degrees(), m_z.degrees());
    }

    template <typename OT>
    internals(internals<OT> const& v) :
      m_x(static_cast<T>(v.x())),
      m_y(static_cast<T>(v.y().radians())),
      m_z(static_cast<T>(v.z().radians()))
    { }

    /*! constructor
    * assigning
    * @param1 to x
    * @param2 to y
    * @param3 to z */
    explicit internals(T const a, angle_type const b, angle_type const c)
      : m_x(a), m_y(b), m_z(c)
    { }


    void init(T const a, angle_type const b, angle_type const c)
    {
      m_x = a;
      m_y = b;
      m_z = c;
    }

    T const& x() const { return m_x; }
    angle_type const& y() const { return m_y; }
    angle_type const& z() const { return m_z; }

    T& x() { return m_x; }
    angle_type& y() { return m_y; }
    angle_type& z() { return m_z; }

    this_type& operator() () { return *this; }

    this_type& operator+= (this_type const& v)
    {
      m_x += v.m_x;
      m_y += v.m_y;
      m_z += v.m_z;
      return *this;
    }
    //! substraction by element
    this_type& operator-= (this_type const& v)
    {
      m_x -= v.m_x;
      m_y -= v.m_y;
      m_z -= v.m_z;
      return *this;
    }
    //! multiplication by element
    this_type& operator*= (this_type const& v)
    {
      m_x *= v.m_x;
      m_y *= v.m_y.radians();
      m_z *= v.m_z.radians();
      return *this;
    }
    //! division by element
    this_type& operator/= (this_type const& v)
    {
      m_x /= v.m_x;
      m_y /= v.m_y.radians();
      m_z /= v.m_z.radians();
      return *this;
    }

    this_type operator- () const
    {
      return this_type(-m_x, -m_y, -m_z);
    }

  };

  template<class T>
  inline internals<T> operator+ (internals<T> a, internals<T> const& b)
  {
    a += b;
    return a;
  }
  template<class T>
  inline internals<T> operator- (internals<T> a, internals<T> const& b)
  {
    a -= b;
    return a;
  }
  template<class T>
  inline internals<T> operator* (internals<T> a, internals<T> const& b)
  {
    a *= b;
    return a;
  }
  template<class T>
  inline internals<T> operator/ (internals<T> a, internals<T> const& b)
  {
    a /= b;
    return a;
  }

  template<class T>
#if defined (SCON_CC11_AUTO_DECLARATOR) && defined(SCON_CC11_F_LATE_DECLARATOR)
  inline auto dot(internals<T> const& a, internals<T> const& b)
    -> decltype(a.x()* b.x() + a.y() * b.y().rad + a.z() * b.z())
#else
  inline T dot(vect3d<T> const& a, vect3d<T> const& b)
#endif
  {
    return scon::dot(a.x(), b.x()) + scon::dot(a.y(), b.y()) + scon::dot(a.z(), b.z());
  }

  template<typename T>
  std::ostream& operator<< (std::ostream& s, const internals<T>& v)
  {
    s << "[" << v.x() << ",";
    s << v.y().degrees() << " d,";
    s << v.z().degrees() << " d]";
    return s;
  }

  template<typename T>
  std::istream& operator>> (std::istream& s, internals<T>& v)
  {
    T a(0.0), d(0.0);
    s >> v.x() >> a >> d;
    v.y() = scon::angles<T>::from_deg(a);
    v.z() = scon::angles<T>::from_deg(d);
    return s;
  }



}

#endif