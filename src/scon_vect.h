#ifndef SCON_VEC_HEADER

#define SCON_VEC_HEADER

#include <cstddef>
#include <iterator>
#include <utility>
#include <numeric>
#include <vector>
#include <iostream>
#include <iomanip>

#include "scon_traits.h"
#include "scon_utility.h"

// types

#include "scon_c3.h"
#include "scon_angle.h"
#include "scon_spherical.h"

#if defined(__GNUG__)

#if !defined(NDEBUG) && !defined(SCON_DEBUG)
#define SCON_DEBUG 1
#elif defined(DEBUG)
#define SCON_DEBUG 2
#endif

#elif defined(_MSC_VER)

#if defined(_DEBUG) && !defined(SCON_DEBUG)
#define SCON_DEBUG 2
#endif

#elif defined(__INTEL_COMPILER)

#if !defined(NDEBUG) && !defined(SCON_DEBUG)
#define SCON_DEBUG
#endif

#endif

#if defined(SCON_DEBUG)
#include <stdexcept>
#endif

namespace scon
{

  /*

  
  ##     ## ########  ######  ########  #######  ########
  ##     ## ##       ##    ##    ##    ##     ## ##     ##
  ##     ## ##       ##          ##    ##     ## ##     ##
  ##     ## ######   ##          ##    ##     ## ########
   ##   ##  ##       ##          ##    ##     ## ##   ##
    ## ##   ##       ##    ##    ##    ##     ## ##    ##
     ###    ########  ######     ##     #######  ##     ##

  
  */

  // Not to be used polymorphically !
  template<class T, class Allocator = std::allocator<T>>
  class vector : public std::vector<T, Allocator>
  { 

    typedef std::vector<T, Allocator> base_type;

  public:

    using allocator_type = typename base_type::allocator_type;
    using value_type = typename base_type::value_type;
    using size_type = typename base_type::size_type;
    using difference_type = typename base_type::difference_type;
    using iterator = typename base_type::iterator;
    using const_iterator = typename base_type::const_iterator;
    using reverse_iterator = typename base_type::reverse_iterator;
    using const_reverse_iterator = typename base_type::const_reverse_iterator;
    using pointer = typename base_type::pointer;
    using const_pointer = typename base_type::const_pointer;
    using reference = typename base_type::reference;
    using const_reference = typename base_type::const_reference;

    // constructor
    explicit vector(const Allocator& alloc = Allocator())
      : base_type(alloc)
    { }
    explicit vector(size_type count, T const & value = T(),
                    Allocator const & alloc = Allocator())
      : base_type(count, value, alloc)
    { }

    template< class InputIt >
    vector(InputIt first, InputIt last,
           Allocator const & alloc = Allocator())
      : base_type(first, last, alloc)
    { }

    vector(vector const & other)
      : base_type(other)
    { }

    vector(vector const & other, Allocator const & alloc)
      : base_type(other, alloc)
    { }

    vector(vector && other)
      : base_type(std::move(other))
    { }
    
    vector(vector && other, Allocator const & alloc)
      : base_type(std::move(other), alloc)
    { }

    vector(std::initializer_list<T> init,
           Allocator const & alloc = Allocator())
      : base_type(init, alloc)
    { }

    vector& operator=(vector const &other)
    {
      base_type::operator=(other);
      return *this;
    }

    vector& operator=(vector &&other)
    {
      base_type::operator=(std::move(other));
      return *this;
    }

    //template<class U, class Ulloc = std::allocator<U>>
    //vector(vector<U, Ulloc> const &u)
    //  : base_type(u.begin(), u.end())
    //{ }

    //using base_type::vector;
    // value assigment
    using base_type::assign;
    // allocator
    using base_type::get_allocator;
    // element access
    using base_type::operator[];
    using base_type::at;
    using base_type::front;
    using base_type::back;
    using base_type::data;
    // iterators
    using base_type::begin;
    using base_type::end;
    using base_type::rbegin;
    using base_type::rend;
    using base_type::cbegin;
    using base_type::cend;
    using base_type::crbegin;
    using base_type::crend;
    // capacity
    using base_type::capacity;
    using base_type::size;
    using base_type::max_size;
    using base_type::empty;
    using base_type::reserve;
    using base_type::shrink_to_fit;
    // modify
    using base_type::clear;
    using base_type::resize;
    using base_type::erase;
    using base_type::emplace;
    using base_type::emplace_back;
    using base_type::push_back;
    using base_type::pop_back;
    using base_type::insert;

    void swap(vector & other)
    {
      base_type::swap(other);
    }

  };

  template<class T, class Alloc>
  void swap(vector<T, Alloc> & a, vector<T, Alloc> & b)
  {
    a.swap(b);
  }

  namespace _vector_detail
  {
    inline int _vector_delim_stream_index()
    {
      static int i = std::ios_base::xalloc();
      return i;
    }
    inline int _vector_terminate_stream_index()
    {
      static int i = std::ios_base::xalloc();
      return i;
    }

    struct _vector_termination_holder
    {
      char d;
      _vector_termination_holder(char c) : d(c) { }
    };

    struct _vector_delimeter_holder
    {
      char d;
      _vector_delimeter_holder(char c) : d(c) { }
    };

    inline std::ostream& operator<< (std::ostream &s,
      _vector_delimeter_holder v)
    {
      s.iword(_vector_detail::_vector_delim_stream_index()) = v.d;
      return s;
    }

    inline std::ostream& operator<< (std::ostream &s,
      _vector_termination_holder v)
    {
      s.iword(_vector_terminate_stream_index()) = v.d;
      return s;
    }

    template<class S>
    char _vector_get_delimeter(S & strm)
    {
      auto delim = strm.iword(_vector_delim_stream_index());
      if (delim == 0 ||
        delim < std::numeric_limits<char>::min() ||
        delim > std::numeric_limits<char>::max()) return '\0';
      return static_cast<char>(delim);
    }

    template<class S>
    char _vector_get_terminator(S & strm)
    {
      auto delim = strm.iword(_vector_terminate_stream_index());
      if (delim == 0 ||
        delim < std::numeric_limits<char>::min() ||
        delim > std::numeric_limits<char>::max()) return '\0';
      return static_cast<char>(delim);
    }

  }

  inline _vector_detail::_vector_delimeter_holder 
    vector_delimeter(char const c)
  {
    return _vector_detail::_vector_delimeter_holder(c);
  }

  inline _vector_detail::_vector_termination_holder
    vector_terminator(char const c)
  {
    return _vector_detail::_vector_termination_holder(c);
  }

  template<class T, class Alloc>
  inline std::ostream& operator<< (std::ostream &s, vector<T, Alloc> const &v)
  {
    std::streamsize const i = s.width();
    auto delim = _vector_detail::_vector_get_delimeter(s);
    auto term = _vector_detail::_vector_get_terminator(s);
    bool delimate = delim != '\0', terminate = term != '\0';
    bool first(true), proceed(true);
    for (auto const & e : v)
    {
      if (delimate && !first && proceed) proceed = proceed && s << delim;
      proceed = proceed && s << std::setw(i) && s << e;
      if (!proceed) break;
      first = false;
    }
    if (proceed && terminate) s << term;
    return s;
  }



  /*


  ##     ## ######## ##       ########  ######## ########   ######
  ##     ## ##       ##       ##     ## ##       ##     ## ##    ##      ##
  ##     ## ##       ##       ##     ## ##       ##     ## ##            ##
  ######### ######   ##       ########  ######   ########   ######     ######
  ##     ## ##       ##       ##        ##       ##   ##         ##      ##
  ##     ## ##       ##       ##        ##       ##    ##  ##    ##      ##
  ##     ## ######## ######## ##        ######## ##     ##  ######


          ##     ## ######## #### ##       #### ######## ##    ##
          ##     ##    ##     ##  ##        ##     ##     ##  ##
          ##     ##    ##     ##  ##        ##     ##      ####
          ##     ##    ##     ##  ##        ##     ##       ##
          ##     ##    ##     ##  ##        ##     ##       ##
          ##     ##    ##     ##  ##        ##     ##       ##
           #######     ##    #### ######## ####    ##       ##

  */


  /*
   ######   ######## ##    ## ######## ########     ###    ##
  ##    ##  ##       ###   ## ##       ##     ##   ## ##   ##
  ##        ##       ####  ## ##       ##     ##  ##   ##  ##
  ##   #### ######   ## ## ## ######   ########  ##     ## ##
  ##    ##  ##       ##  #### ##       ##   ##   ######### ##
  ##    ##  ##       ##   ### ##       ##    ##  ##     ## ##
   ######   ######## ##    ## ######## ##     ## ##     ## ########
  */

  // cardinality

  template<class T>
  std::size_t cardinality(T const &)
  {
    return 1U;
  }

  template<class T>
  T arithmetic_mean(T const &v)
  {
    return v;
  }

  template<class T, class U = T, bool isa =
    std::is_arithmetic<T>::value && std::is_arithmetic<U>::value>
  struct dot_trait;

  // dot

  template<class T, class U>
  struct dot_trait < T, U, true >
  {
    typedef decltype(std::declval<T>()*std::declval<U>()) type;
  };

  template<class T, class U = T>
  using dot_type = typename dot_trait<T, U>::type;

  template<class T, class U = T>
  typename std::enable_if<std::is_arithmetic<T>::value
    && std::is_arithmetic<U>::value,
    dot_type<T, U >> ::type
    dot(T const &x, U const &y)
  {
    return x*y;
  }


  template<class T>
  using length_type = typename float_helper<typename dot_trait<T, T>::type>::type;

  // length

  template<class T>
  length_type<T> geometric_length(T const &a)
  {
    using std::sqrt;
    length_type<T> const d 
      = static_cast<length_type<T>>(dot(a, a));
    return sqrt(d);
  }

  template<class T>
  auto len(T const &a) -> decltype(geometric_length(a))
  {
    return geometric_length(a);
  }

  /*
       ######   #######  
      ##    ## ##     ## 
      ##              ## 
      ##        #######  
      ##              ## 
      ##    ## ##     ## 
       ######   #######   
  */

  // cardinality

  template<class T>
  typename std::enable_if<std::is_fundamental<T>::value,
    std::size_t>::type cardinality(c3<T> const &v)
  {
    return 3U;
  }

  template<class T>
  typename std::enable_if<std::is_fundamental<T>::value,
    T>::type arithmetic_mean(c3<T> const &v)
  {
    return (v.x()+v.y()+v.z()) / T(3);
  }

  template<class T>
  typename std::enable_if<!std::is_fundamental<T>::value, 
    std::size_t>::type cardinality(c3<T> const &v)
  {
    return cardinality(v.x()) + cardinality(v.y()) + cardinality(v.z());
  }

  template<class T>
  typename std::enable_if<std::is_fundamental<T>::value, 
    std::size_t>::type cardinality(scon::vector<c3<T>> const &v)
  {
    return std::size_t(3U * v.size());
  }

  // dot

  template<class T, class U>
  struct dot_trait < c3<T>, c3<U>, false >
  {
    typedef typename dot_trait<T, U>::type type;
  };

  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value, 
    T>::type dot(c3<T> const &a, c3<T> const &b)
  {
    return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
  }

  template<class T>
  typename std::enable_if<!std::is_arithmetic<T>::value, 
    typename dot_trait<T, T>::type>::type dot(c3<T> const &a, c3<T> const &b)
  {
    return dot(a.x(), b.x()) + dot(a.y(), b.y()) + dot(a.z(), b.z());
  }

  template<class T>
  c3<c3<T>> transpose(c3<c3<T>> const &matrix)
  {
    return c3<c3<T>>{
      { matrix.x().x(), matrix.y().x(), matrix.z().x() },
      { matrix.x().y(), matrix.y().y(), matrix.z().y() },
      { matrix.x().z(), matrix.y().z(), matrix.z().z() }
    };
  }


  /*
       ###    ##    ##  ######   ##       ######## 
      ## ##   ###   ## ##    ##  ##       ##       
     ##   ##  ####  ## ##        ##       ##       
    ##     ## ## ## ## ##   #### ##       ######   
    ######### ##  #### ##    ##  ##       ##       
    ##     ## ##   ### ##    ##  ##       ##       
    ##     ## ##    ##  ######   ######## ######## 
  */

  // cardinality

  template<class T>
  std::size_t cardinality(ang<T> const &v)
  {
    return 1U;
  }

  template<class T>
  std::size_t cardinality(scon::vector<ang<T>> const &v)
  {
    return static_cast<std::size_t>(v.size());
  }

  // dot

  template<class T>
  struct dot_trait < ang<T>, ang<T>, false >
  {
    using type = T;
  };
  template<class T>
  struct dot_trait < ang<T>, T, false >
  {
    using type = T;
  };
  template<class T>
  struct dot_trait < T, ang<T>, false >
  {
    using type = T;
  };

  template<class T>
  T dot(ang<T> const &a, ang<T> const &b)
  {
    return a.radians() * b.radians();
  }
  template<class T>
  T dot(ang<T> const &a, T b)
  {
    return a.radians() * b;
  }
  template<class T>
  T dot(T a, ang<T> const &b)
  {
    return a * b.radians();
  }

  /*
       ######  ########  ##     ## ######## ########  ####  ######
      ##    ## ##     ## ##     ## ##       ##     ##  ##  ##    ##
      ##       ##     ## ##     ## ##       ##     ##  ##  ##
       ######  ########  ######### ######   ########   ##  ##
            ## ##        ##     ## ##       ##   ##    ##  ##
      ##    ## ##        ##     ## ##       ##    ##   ##  ##    ##
       ######  ##        ##     ## ######## ##     ## ####  ######
  */

  template<class T>
  struct dot_trait < sphericals<T>, c3<T>, false >
  {
    typedef T type;
  };

  template<class T>
  struct dot_trait < c3<T>, sphericals<T>, false >
  {
    typedef T type;
  };

  template<class T>
  struct dot_trait < sphericals<T>, sphericals<T>, false >
  {
    typedef T type;
  };

  template<class T>
  T dot(sphericals<T> const & a, c3<T> const & b)
  {
    return dot(a.radius(), b.x())
      + dot(a.inclination(), b.y())
      + dot(a.azimuth(), b.z());
  }

  template<class T, bool B>
  T dot(c3<T> const & a, sphericals<T> const & b)
  {
    return dot(b, a);
  }

  template<class T>
  T dot(sphericals<T> const & a, sphericals<T> const & b)
  {
    return dot(a.radius(), b.radius())
      + dot(a.inclination(), b.inclination())
      + dot(a.azimuth(), b.azimuth());
  }

  /*
  ##     ## ########  ######  ########  #######  ########
  ##     ## ##       ##    ##    ##    ##     ## ##     ##
  ##     ## ##       ##          ##    ##     ## ##     ##
  ##     ## ######   ##          ##    ##     ## ########
   ##   ##  ##       ##          ##    ##     ## ##   ##
    ## ##   ##       ##    ##    ##    ##     ## ##    ##
     ###    ########  ######     ##     #######  ##     ##
  */

  // cardinality

  template<class T>
  typename std::enable_if<std::is_fundamental<T>::value,
    std::size_t>::type cardinality(scon::vector<T> const &v)
  {
    return static_cast<std::size_t>(v.size());
  }

  template<class T>
  typename std::enable_if<!std::is_fundamental<T>::value,
    std::size_t>::type cardinality(scon::vector<T> const &v)
  {
    std::size_t c(0u);
    auto const e = v.end();
    for (auto i = v.begin(); i != e; ++i) c += cardinality(*i);
    return c;
  }

  // dot

  template<class T, class U>
  struct dot_trait < scon::vector<T>, scon::vector<U>, false >
  {
    typedef typename dot_trait<T, U>::type type;
  };

  template<class T, class U = T>
  typename dot_trait<T, U>::type dot(scon::vector<T> const &x, scon::vector<U> const &y)
  {
    using scon::dot;
    using R = typename dot_trait<T, U>::type;
    R r = R();
    auto const ex = x.end();
    auto const ey = y.end();
    auto ix = x.begin();
    auto iy = y.begin();
#if defined(SCON_DEBUG)
    if (std::distance(ix, ex) != std::distance(iy, ey))
    {
      std::cout << "broken vectors: " << x << " , " << y << "\n";
      throw std::logic_error("Unqually sized vectors given for dot product.");
    }
#endif
    for (; ix != ex && iy != ey; ++ix, ++iy)
    {
      r += scon::dot(*ix, *iy);
    }
    return r;
  }

  template<class T, class A>
  typename vector<T, A>::value_type arithmetic_mean(vector<T, A> const &v)
  {
    using std::distance;
    auto bv = v.begin();
    auto const ev = v.end();
    if (distance(bv, ev) > 0)
    {
      typename vector<T, A>::value_type r = *(bv++);
      for (; bv != ev; ++bv)
      {
        r += *bv;
      }
      return r;
    }
    return typename vector<T, A>::value_type();
  }


  /*
       ######   #######      ##     ##    ###    ######## ##     ## 
      ##    ## ##     ##     ###   ###   ## ##      ##    ##     ## 
      ##              ##     #### ####  ##   ##     ##    ##     ## 
      ##        #######      ## ### ## ##     ##    ##    ######### 
      ##              ##     ##     ## #########    ##    ##     ## 
      ##    ## ##     ##     ##     ## ##     ##    ##    ##     ## 
       ######   #######      ##     ## ##     ##    ##    ##     ##   
  */

  template<class T>
  typename std::enable_if<std::is_floating_point<T>::value, T>::type
    distance(c3<T> a, c3<T> const &b)
  {
    using std::sqrt;
    a -= b;
    return sqrt(dot(a, a));
  }

  /**
   * Returns eucledean norm of vector
   */
  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value, 
    length_type<T>>::type geometric_length(c3<T> const &a)
  {
    using std::sqrt;
    return sqrt(static_cast<_c3::F<T>>(a.x()*a.x() + a.y()*a.y() + a.z()*a.z()));
  }


  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value, 
    c3<T>>::type cross(c3<T> const &a, c3<T> const &b)
  {
    return c3<T>(
      a.y()*b.z() - a.z()*b.y(),
      a.z()*b.x() - a.x()*b.z(),
      a.x()*b.y() - a.y()*b.x()
      );
  }

  //! return angle between two vectors
  template<class T>
  ang<T> angle(c3<T> const &a, c3<T> const &b)
  {
    using std::abs;
    using std::acos;
    T const norm(geometric_length(a)*geometric_length(b));
    T sc(dot(a, b));
    if (abs(norm) > T())
    {
      sc /= norm;
      return ang<T>::from_rad(acos(sc));
    }
    return ang<T>::from_rad(T());
  }

  template<class T>
  ang<T> angle_refined(c3<T> const &a, c3<T> const &b)
  {
    using std::atan2;
    return ang<T>::from_rad(atan2(geometric_length(cross(a, b)), dot(a, b)));
  }

  template<class T>
  ang<T> dihedral(c3<T> const &axis, c3<T> const &a, c3<T> const &b)
  {
    using std::abs;
    c3<T> const t(cross(axis, a)), u(cross(axis, b)), tu(cross(t, u));
    T norm = geometric_length(axis)*geometric_length(tu);
    if (abs(norm) > T()  && (dot(axis, tu) / norm) < T())
    {
      return -angle_refined(t, u);
    }
    return angle_refined(t, u);
  }

  template<class T>
  ang<T> angle_reference(c3<T> const &a, c3<T> const &b, 
    c3<T> const &orthonormal_reference)
  {
    using std::atan2;
    // this cross v
    auto this_x_v(cross(a, b));
    // Len of cross
    auto const x_len = geometric_length(this_x_v);
    // Get angle
    auto ret = ang<T>::from_rad(atan2(x_len, dot(a, b)));
    // Check reference
    if (dot(orthonormal_reference, this_x_v) < T()) ret = -ret;
    return ret;
  }

  template<class T>
  sphericals<T, ang<T>> spherical(c3<T> const & p, c3<T> const & origin,
    c3<T> const & zenith, c3<T> const & azimuth_reference)
  {
    using std::atan2;
    c3<T> const B(p - origin);
    c3<T> const ZxB(cross(zenith, B));
    c3<T> const ZxA(cross(zenith, azimuth_reference));
    auto const inclination = ang<T>::from_rad(atan2(geometric_length(ZxB), dot(zenith, B)));
    auto const azimuth = angle_reference(ZxA, ZxB, normalized(zenith));
    return sphericals<T>(geometric_length(B), inclination, azimuth);
  }

  //! compute vect3dor position based on the Natural Extension Reference Frame Algorithm
  template<class T, class U>
  c3<T> appendNERF(c3<T> const &P, c3<T> const &X, c3<T> const &Y, 
    T distance, U const & angle, U dihedral)
  {
    using std::sin;
    using std::cos;
    //dihedral.invert();
    dihedral = -dihedral;
    auto const anglesin = sin(angle);
    c3<T> const R(
      static_cast<T>(distance*cos(angle)),
      static_cast<T>(distance*cos(dihedral)*anglesin),
      static_cast<T>(distance*sin(dihedral)*anglesin));
    c3<T> m(normalized(Y));
    c3<T> n(normalized(cross(X, m)));
    c3<T> o(cross(n, m)); // normalized since n and m are normalized
    m *= R.x();
    o *= R.y();
    n *= R.z();
    return P + m + o + n;
  }

  template<class T>
  c3<T> append_spherical_NERF(c3<T> const &X, c3<T> zenith, c3<T> azimuth_ref,
    sphericals<T> const &spx)
  {
    T const inclination_sin = sin(spx.inclination());
    c3<T> const CC(spx.radius()*inclination_sin*cos(spx.azimuth()),
      spx.radius()*inclination_sin*sin(spx.azimuth()),
      spx.radius()*cos(spx.inclination()));
    normalize(zenith);
    c3<T> const ZxA(normalized(cross(zenith, azimuth_ref)));
    c3<T> const ZxAxZ(cross(ZxA, zenith));
    return X + ZxAxZ*CC.x() + ZxA*CC.y() + zenith*CC.z();
  }

  namespace vect_detail
  {
    template<class T>
    void rotate_x(c3<T> &v, T const theta)
    {
      using std::ceil;
      auto sin_t = std::sin(theta);
      auto cos_t = std::cos(theta);
      auto y = v.y();
      auto z = v.z();
      v.y() = y * cos_t - z * sin_t;
      v.z() = z * cos_t + y * sin_t;

      return c3<T>(ceil(v.x()), ceil(v.y()), ceil(v.z()));
    }
  }

  //template<class T>
  //typename std::enable_if<std::is_floating_point<T>::value, c3<T>>::type 
  //  rotated(c3<T> const &point, typename c3<T>::type const ang_x, 
  //    typename c3<T>::type const  ang_y, typename c3<T>::type const ang_z, 
  //    c3<T> const center = c3<T>{})
  //{
  //  using std::cos;
  //  using std::sin;
  //  auto p = point - center;
  //  auto ca = cos(scon::ang<T>::from_deg(ang_x));
  //  auto sa = sin(scon::ang<T>::from_deg(ang_x));
  //  auto cb = cos(scon::ang<T>::from_deg(ang_y));
  //  auto sb = sin(scon::ang<T>::from_deg(ang_y));
  //  auto cg = cos(scon::ang<T>::from_deg(ang_z));
  //  auto sg = sin(scon::ang<T>::from_deg(ang_z));
  //  auto X = p.x();
  //  auto Y = p.y();
  //  auto Z = p.z();
  //  c3<T> r;
  //  r.x() = ((((cb*cg)*X) + ((cb*-sg)*Y)) + (sb*Z));
  //  r.y() = ((((((-sa*-sb)*cg) + (ca*sg))*X) + ((((-sa*-sb)*-sg) + (ca*cg))*Y)) + ((-sa*cb)*Z));
  //  r.z() = ((((((ca*-sb)*cg) + (sa*sg))*X) + ((((ca*-sb)*-sg) + (sa*cg))*Y)) + ((ca*cb)*Z));
  //  return center + r;
  //}




  /*
       ######   #######     ##     ##     ######   #######
      ##    ## ##     ##     ##   ##     ##    ## ##     ##
      ##              ##      ## ##      ##              ##
      ##        #######        ###       ##        #######
      ##              ##      ## ##      ##              ##
      ##    ## ##     ##     ##   ##     ##    ## ##     ##
       ######   #######     ##     ##     ######   #######
  */

  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value,
    T>::type det(c3<c3<T>> const &m)
  {
    return  m.x().x()*m.y().y()*m.z().z() +
      m.x().z()*m.y().x()*m.z().y() +
      m.x().y()*m.y().z()*m.z().x() -
      m.x().z()*m.y().y()*m.z().x() -
      m.x().y()*m.y().x()*m.z().z() -
      m.x().x()*m.y().z()*m.z().y();
  }

  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value,
    void>::type transpose(c3<c3<T>> &m)
  {
    using std::swap;
    swap(m.x().y(), m.y().x());
    swap(m.x().z(), m.z().x());
    swap(m.y().z(), m.z().y());
  }

  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value,
    c3<c3<T>>>::type transposed(c3<c3<T>> const &m)
  {
    c3<c3<T>> r(m);
    transpose(r);
    return r;
  }

  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value,
    void>::type invert(c3<c3<T>> &m)
  {
    using std::abs;
    T const determinant = det(m);
    if (abs(determinant) > T())
    {
      c3<c3<T>> tmp(m); // copy
      // recalc matrix from copy
      m.x() = cross(tmp.y(), tmp.z());
      m.y() = cross(tmp.z(), tmp.x());
      m.z() = cross(tmp.x(), tmp.y());
      // resize
      T const dmv = T(1) / determinant;
      m.x() *= dmv;
      m.y() *= dmv;
      m.z() *= dmv;
      transpose(m);
    }
  }

  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value,
    c3<c3<T>>>::type inverted(c3<c3<T>> const &m)
  {
    c3<c3<T>> r(m);
    invert(r);
    return r;
  }

  //template<class T>
  //typename std::enable_if<std::is_floating_point<T>::value, c3<T>>::type
  //  operator*(c3<c3<T>> const & m1, c3<c3<T>> const & m2)
  //{
  //  return c3<c3<T>>{
  //    { m1.x().x()*m2.x().x() + m1.x().y()*m2.y().x() + m1.x().z()*m2.z().x(),
  //      m1.x().x()*m2.x().y(), m1.x().y()*m2.y().y(), m1.x().z()*m2.z().y(),
  //      m1.x().x()*m2.x().z(), m1.x().y()*m2.y().z(), m1.x().z()*m2.z().z() },

  //    { m1.y().x()*m2.x().x() + m1.y().y()*m2.y().x() + m1.y().z()*m2.z().x(),
  //      m1.y().x()*m2.x().y(), m1.y().y()*m2.y().y(), m1.y().z()*m2.z().y(),
  //      m1.y().x()*m2.x().z(), m1.y().y()*m2.y().z(), m1.y().z()*m2.z().z() },
  //    { TODO }
  //  };
  //}

  template<class T>
  typename std::enable_if<std::is_floating_point<T>::value, c3<T>>::type
    operator*(c3<c3<T>> const & matrix, c3<T> const & vec)
  {
    return c3<T>{
      matrix.x().x()*vec.x() + matrix.x().y()*vec.y() + matrix.x().z()*vec.z(),
      matrix.y().x()*vec.x() + matrix.y().y()*vec.y() + matrix.y().z()*vec.z(),
      matrix.z().x()*vec.x() + matrix.z().y()*vec.y() + matrix.z().z()*vec.z()
    };
  }

  namespace rotation_detail
  {
    template<class T>
    typename std::enable_if<std::is_floating_point<T>::value, c3<T>>::type
      tangle(c3<T> const & a, c3<T> const & b)
    {
      return c3<T>{
        a.y()*b.y() + a.z()*b.z(),
          a.x()*b.x() + a.z()*b.z(),
          a.x()*b.x() + a.y()*b.y()
      };
    }

    template<class T>
    typename std::enable_if<std::is_floating_point<T>::value, c3<T>>::type
      R_xz(c3<T> const & p, c3<T> const & axis)
    {
      using std::abs;
      auto const L = std::sqrt(axis.x()*axis.x() + axis.y()*axis.y());
      return c3<T>{
        axis.x()*p.x() / L + axis.y()*p.y() / L,
        -axis.y()*p.x() / L + axis.x()*p.y() / L,
        p.z()
      };
    }

  }


  template<class T>
  typename std::enable_if<std::is_floating_point<T>::value, c3<T>>::type
    rotated(c3<T> point, c3<T> axis, c3<T> const & center, ang<T> const &angle)
  {
    using scon::cos; using scon::sin; using std::abs;
    normalize(axis);
    if (abs(angle.radians()) < 1.e-16) return point;
    auto const S = sin(angle);
    auto const C = cos(angle);
    //std::cout << S << ", " << C << "\n";
    auto p2 = point * C;
    auto c2 = center * rotation_detail::tangle(axis, axis);
    auto m = dot(axis, point);
    auto d2 = axis * (rotation_detail::tangle(center, axis) - dot(axis, point)) * (1 - C);
    auto s2 = cross(center, axis) + cross(axis, point) * S;
    return c2 - d2 + p2 + s2;
  }
  

  template<class T>
  void randomize(c3<T> &v)
  {
    using scon::randomize;
    randomize(v.x());
    randomize(v.y());
    randomize(v.z());
  }

  template<class T>
  void randomize(c3<T> &v, T const &low, T const &high)
  {
    using scon::randomize;
    randomize(v.x(), low, high);
    randomize(v.y(), low, high);
    randomize(v.z(), low, high);
  }

  template<class T>
  struct uniform_c3_distribution
  {
    using type = typename uniform_ditribution_trait < T >::type;
    using result_type = c3 < T > ;
    using param_type = typename type::param_type;
    type d;
    uniform_c3_distribution() : d() { }
    uniform_c3_distribution(T const low, T const high) : d(low, high) { }
    uniform_c3_distribution(param_type const &p) : d(p) { }
    template<class G>
    c3<T> operator() (G & generator)
    {
      return c3<T>(d(generator), d(generator), d(generator));
    }
    template<class G>
    c3<T> operator() (G & generator, param_type const &p)
    {
      return c3<T>(d(generator, p),
        d(generator, p), d(generator, p));
    }
  };

  template<class T>
  struct uniform_ditribution_trait < c3<T>, true >
  {
    using type = uniform_c3_distribution < T >;
  };

  template<class T, class U>
  typename std::enable_if<is_range<T>::value, void>::type
    randomize(T &v, U const &low, U const & high)
  {
    using std::begin;
    using std::end;
    using std::for_each;
    for_each(begin(v), end(v), uniform_randomizer<range_value<T>>(low, high));
  }

  template<class T>
  typename std::enable_if<is_range<T>::value, void>::type
    randomize(T &v)
  {
    using std::begin;
    using std::end;
    using std::for_each;
    for_each(begin(v), end(v), Randomizer<range_value<T>>());
  }

  template<class T>
  void limit(c3<T> &v, T const mini, T const maxi)
  {
    using scon::limit;
    limit(v.x(), mini, maxi);
    limit(v.y(), mini, maxi);
    limit(v.z(), mini, maxi);
  }

  template<class V, class T>
  struct Limiter
  {
    T mini, maxi;
    Limiter(T const &minimum, T const &maximum)
      : mini(minimum), maxi(maximum)
    { }
    template<class U>
    void operator () (U & v)
    {
      limit(v, mini, maxi);
    }
  };

  template<class T, class U>
  typename std::enable_if<is_range<T>::value, void>::type
    limit(T &v, U const & minimum, U const & maximum)
  {
    using std::begin;
    using std::end;
    using std::for_each;
    //std::bind(, std::placeholders::_1, minimum, maximum);
    for_each(begin(v), end(v),
      Limiter<T, U>(minimum, maximum));
  }

  template<class T>
  c3<T> cog_align(c3<T> const &P, c3<T> const &cog_1, c3<T> const &cog_2)
  {
    return P - cog_1 + cog_2 - cog_1;
  }
  template<class T>
  c3<T> align_cog(c3<T> const &P, c3<T> const &cog_1, c3<T> const &cog_2)
  {
    return P - cog_1 + cog_2;
  }


  /*
      ########     ###    ##    ##  ######   ########    ##     ##    ###    ######## ##     ##
      ##     ##   ## ##   ###   ## ##    ##  ##          ###   ###   ## ##      ##    ##     ##
      ##     ##  ##   ##  ####  ## ##        ##          #### ####  ##   ##     ##    ##     ##
      ########  ##     ## ## ## ## ##   #### ######      ## ### ## ##     ##    ##    #########
      ##   ##   ######### ##  #### ##    ##  ##          ##     ## #########    ##    ##     ##
      ##    ##  ##     ## ##   ### ##    ##  ##          ##     ## ##     ##    ##    ##     ##
      ##     ## ##     ## ##    ##  ######   ########    ##     ## ##     ##    ##    ##     ##
  */

  // range / range

  template<class T>
  typename std::enable_if<is_range<T>::value, T>::type operator- (T v)
  {
    using std::begin;
    using std::end;
    auto const e = end(v);
    for (auto i(begin(v)); i != e; ++i) *i = -*i;
    return v;
  }

  template<class T, class U = T>
  typename std::enable_if<is_range<T>::value && is_range<U>::value, 
    T>::type& operator+= (T & v, U const &w)
  {
    using std::begin;
    using std::end;
    auto const ev = end(v);
    auto const ew = end(w);
    auto i(begin(v));
    auto j(begin(w));
    for (; i != ev && j != ew; ++i, ++j) *i += *j;
    return v;
  }

  template<class T, class U = T>
  typename std::enable_if<is_range<T>::value && is_range<U>::value, 
    T>::type operator+ (T v, U const &w)
  {
    v += w;
    return v;
  }

  template<class T, class U = T>
  typename std::enable_if<is_range<T>::value && is_range<U>::value, 
    T>::type& operator-= (T & v, U const &w)
  {
    using std::begin;
    using std::end;
    auto const ev = end(v);
    auto const ew = end(w);
    auto i(begin(v));
    auto j(begin(w));
    for (; i != ev && j != ew; ++i, ++j) *i -= *j;
    return v;
  }

  template<class T, class U = T>
  typename std::enable_if<is_range<T>::value && is_range<U>::value, 
    T>::type operator- (T v, U const &w)
  {
    v -= w;
    return v;
  }

  template<class T, class U = T>
  typename std::enable_if<is_range<T>::value && is_range<U>::value, 
    T>::type& operator*= (T & v, U const &w)
  {
    using std::begin;
    using std::end;
    auto const ev = end(v);
    auto const ew = end(w);
    auto i(begin(v));
    auto j(begin(w));
    for (; i != ev && j != ew; ++i, ++j) *i *= *j;
    return v;
  }

  template<class T, class U = T>
  typename std::enable_if<is_range<T>::value && is_range<U>::value, 
    T>::type operator* (T v, U const &w)
  {
    v *= w;
    return v;
  }

  template<class T, class U = T>
  typename std::enable_if<is_range<T>::value && is_range<U>::value, 
    T>::type& operator/= (T & v, U const &w)
  {
    using std::begin;
    using std::end;
    auto const ev = end(v);
    auto const ew = end(w);
    auto i(begin(v));
    auto j(begin(w));
    for (; i != ev && j != ew; ++i, ++j) *i /= *j;
    return v;
  }

  template<class T, class U = T>
  typename std::enable_if<is_range<T>::value && is_range<U>::value, 
    T>::type operator/ (T v, U const &w)
  {
    v /= w;
    return v;
  }
  
  // range / value

  template<class T, class U>
  typename std::enable_if<is_range<T>::value && !is_range<U>::value, 
    T>::type& operator+= (T & v, U const & w)
  {
    using std::begin;
    using std::end;
    auto const ev = end(v);
    for (auto i(begin(v)); i != ev; ++i) *i += w;
    return v;
  }

  template<class T, class U>
  typename std::enable_if<is_range<T>::value && !is_range<U>::value, 
    T>::type operator+ (T v, U const & w)
  {
    v += w;
    return v;
  }

  template<class T, class U>
  typename std::enable_if<is_range<T>::value && !is_range<U>::value, 
    T>::type& operator-= (T & v, U const & w)
  {
    using std::begin;
    using std::end;
    auto const ev = end(v);
    for (auto i(begin(v)); i != ev; ++i) *i -= w;
    return v;
  }

  template<class T, class U>
  typename std::enable_if<is_range<T>::value && !is_range<U>::value, 
    T>::type operator- (T v, U const & w)
  {
    v -= w;
    return v;
  }

  template<class T, class U>
  typename std::enable_if<is_range<T>::value && !is_range<U>::value, 
    T>::type& operator*= (T & v, U const & w)
  {
    using std::begin;
    using std::end;
    auto const ev = end(v);
    for (auto i(begin(v)); i != ev; ++i) *i *= w;
    return v;
  }

  template<class T, class U>
  typename std::enable_if<is_range<T>::value && !is_range<U>::value, 
    T>::type operator* (T v, U const & w)
  {
    v *= w;
    return v;
  }

  template<class T, class U>
  typename std::enable_if<is_range<T>::value && !is_range<U>::value, 
    T>::type& operator/= (T & v, U const & w)
  {
    using std::begin;
    using std::end;
    auto const ev = end(v);
    for (auto i(begin(v)); i != ev; ++i) *i /= w;
    return v;
  }

  template<class T, class U>
  typename std::enable_if<is_range<T>::value && !is_range<U>::value, 
    T>::type operator/ (T v, U const & w)
  {
    v /= w;
    return v;
  }

  template<class T>
  typename std::enable_if<is_range<T>::value, T>::type reciprocal(T v)
  {
    using std::begin;
    using std::end;
    auto const e = end(v);
    for (auto i(begin(v)); i != e; ++i) *i = T(1) / *i;
    return v;
  }







  /*
       ######   ######## ##    ## ######## ########  ####  ######
      ##    ##  ##       ###   ## ##       ##     ##  ##  ##    ##
      ##        ##       ####  ## ##       ##     ##  ##  ##
      ##   #### ######   ## ## ## ######   ########   ##  ##
      ##    ##  ##       ##  #### ##       ##   ##    ##  ##
      ##    ##  ##       ##   ### ##       ##    ##   ##  ##    ##
       ######   ######## ##    ## ######## ##     ## ####  ######
  */


  template<class T>
  typename float_helper<typename dot_trait<T, T>::type>::type root_mean_square(T const & c)
  { // root of the mean of the squares of c
    using F = typename float_helper<typename dot_trait<T, T>::type>::type;
    using std::sqrt;
    return sqrt(static_cast<F>(dot(c, c)) /
                static_cast<F>(cardinality(c)));
  }

  template<class T, class U = T>
  typename float_helper<T>::type root_mean_square_deviation(T const &x, U const &y)
  { // root of the mean of the squares of the deviation of two vecs
    return root_mean_square(x - y);
  }

  template<class T>
  T normalized(T x)
  { // normalization (result: geometric_length == 1)
    x /= geometric_length(x);
    return x;
  }

  template<class T>
  void normalize(T& x)
  { // normalization (result: geometric_length == 1)
    x /= geometric_length(x);
  }

  template<class T, class U = T>
  U projected_on_normal(T const & x, U const & direction)
  {
    return direction * dot(x, direction);
  }

  template<class T, class U = T>
  U projected_on(T const & x, U direction)
  {
    return projected_on_normal(x, normalized(direction));
  }

  template<class T, class U = T>
  void project_on_normal(T & x, U const & direction)
  {
    x = direction * dot(x, direction);
  }

  template<class T, class U = T>
  void project_on(T & x, U const & direction)
  {
    project_on_normal(x, normalized(direction));
  }

  // orthogonalization of a to b

  template<class T, class U = T>
  T orthogonalized_to_normal(T const & x, U const & direction)
  {
    return x - projected_on_normal(x, direction);
  }

  template<class T, class U = T>
  T orthogonalized_to(T const & x, U const & direction)
  {
    return x - projected_on(x, direction);
  }

  template<class T, class U = T>
  void orthogonalize_to_normal(T & x, U const & direction)
  {
    x -= projected_on_normal(x, direction);
  }

  template<class T, class U = T>
  void orthogonalize_to(T & x, U const & direction)
  {
    x -= projected_on_normal(x, normalized(direction));
  }

  template<class T>
#if defined(_MSC_VER)
  typename std::enable_if<
    std::is_trivially_copyable<T>::value, 
    std::ostream>::type & 
#else
  std::ostream &
#endif
    write(std::ostream &stream, scon::vector<T> const &v)
  {
    auto p = reinterpret_cast<char const *>(v.data());
    return stream.write(p, v.size()*sizeof(T));
  }

  template<class T>
#if defined(_MSC_VER)
  typename std::enable_if<
    std::is_trivially_copyable<T>::value,
    std::istream>::type &
#else
  std::istream &
#endif
    read(std::istream & stream, scon::vector<T> &v)
  {
    std::istream::pos_type const current_streampos(stream.tellg());
    stream.seekg(0, std::ios_base::end);
    std::istream::pos_type const pos_left(stream.tellg() - current_streampos);
    stream.seekg(current_streampos, std::ios_base::beg);
    if (pos_left < (sizeof(T)*v.size())) throw std::runtime_error("read vector from stream: Stream too short.");
    auto p = reinterpret_cast<char *>(v.data());
    stream.read(p, v.size()*sizeof(T));
    return stream;
  }


}

#endif
