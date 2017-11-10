#if !defined(SCON_ITERATOR_HEADER)
#define SCON_ITERATOR_HEADER

#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <algorithm>
#include <cstddef>
#include <type_traits>
#include "scon.h"
#include "scon_traits.h"

namespace scon
{

  namespace detail
  {
    template<class I, bool p = std::is_pointer<I>::value>
    struct itph
    {
      static inline typename std::iterator_traits < I >::pointer
        _ret(I const x) { return x; }
    };
    template<class I>
    struct itph < I, false >
    {
      static inline typename std::iterator_traits < I >::pointer
        _ret(I const x) { return x.operator->(); }
    };
  }

  /*

  Static / compile time stride iterator

  */

  template<class Iterator,
    typename Iterator::difference_type N>
  class static_stride_iterator
    : public std::iterator_traits<Iterator>
  {

    typedef std::iterator_traits<Iterator> It;

  public:

    typedef typename It::difference_type   difference_type;
    typedef typename It::iterator_category iterator_category;
    typedef typename It::pointer           pointer;
    typedef typename It::reference         reference;
    typedef typename It::value_type        value_type;

    Iterator _begin() const { return m_begin; }
    difference_type _offset() const { return offset; }

  private:

    using my_type = typename std::enable_if<
      std::is_same<typename It::iterator_category,
      std::random_access_iterator_tag>::value,
      static_stride_iterator<Iterator, N >> ::type;

    difference_type offset;
    Iterator m_begin;

  public:

    static_stride_iterator(Iterator it = Iterator())
      : offset(), m_begin(it)
    { }

    template<class Other>
    static_stride_iterator(static_stride_iterator<Other, N> const &rhs)
      : offset(rhs._offset()), m_begin(rhs._begin())
    { } // Allowed if Other is implicitly convertibe to Iterator

    reference operator* () const
    {
      return *(m_begin + offset);
    }
    pointer operator->() const
    {
      return detail::itph<Iterator>::_ret(m_begin + offset);
    }
    reference operator[] (difference_type const & i) const
    {
      return (m_begin + offset)[i*N];
    }
    // increment and decrement
    my_type & operator++ ()
    { // increment then return
      offset += N;
      return *this;
    }
    my_type & operator-- ()
    { // decrement then return
      offset -= N;
      return *this;
    }
    my_type operator++ (int)
    { // copy, increment this, return copy
      my_type ret(*this);
      ++*this;
      return ret;
    }
    my_type operator-- (int)
    { // copy, decrement this, return copy
      my_type ret(*this);
      --*this;
      return ret;
    }
    // Arithmetics using difference_type
    my_type & operator+= (difference_type d)
    {
      offset += d*N;
      return *this;
    }
    my_type & operator-= (difference_type d)
    {
      return (*this += -d);
    }
    my_type operator+ (difference_type d) const
    {
      my_type ret(*this);
      ret += d;
      return ret;
    }
    my_type operator- (difference_type d) const
    {
      return this->operator+(-d);
    }

    // Arithmetics with this type
    difference_type operator- (my_type const & rhs) const
    {
      return (offset - rhs.offset) / N;
    }

    // Comparison
    bool operator== (my_type const &rhs) const
    {
      return /*m_begin == rhs.m_begin &&*/ offset == rhs.offset;
    }
    bool operator!= (my_type const &rhs) const { return !operator==(rhs); }
    bool operator<  (my_type const &rhs) const { return offset < rhs.offset; }
    bool operator>= (my_type const &rhs) const { return !operator<(rhs); }
    bool operator>  (my_type const &rhs) const { return offset > rhs.offset; }
    bool operator<= (my_type const &rhs) const { return !operator>(rhs); }

  };

  template<class Iterator>
  class stride_iterator : 
    public std::iterator_traits<Iterator>
  {

    using my_type = typename std::enable_if <
      std::is_same<
      typename Iterator::iterator_category,
      std::random_access_iterator_tag
      >::value, stride_iterator < Iterator >
    > ::type;


  public:

    typedef std::iterator_traits<Iterator> iter_trait_type;
    typedef typename iter_trait_type::difference_type difference_type;
    typedef typename iter_trait_type::iterator_category iterator_category;
    typedef typename iter_trait_type::pointer pointer;
    typedef typename iter_trait_type::reference reference;
    typedef typename iter_trait_type::value_type value_type;

    stride_iterator() : N(), offset(), m_begin() { }
    stride_iterator(Iterator it, difference_type const stride)
      : N(stride), offset(0U), m_begin(it)
    { }

    template<class Other>
    stride_iterator(stride_iterator<Other> const &rhs)
      :N(rhs.N), offset(rhs.offset), m_begin(rhs.m_begin)
    { }

    /*
    Operators required for random_access_iterator category
    */
    reference operator* () const { return *(m_begin + offset); }
    pointer   operator->() const
    {
      return detail::itph<Iterator>::_ret(m_begin + offset);
    }
    reference operator[] (difference_type const & i) const
    {
      return (m_begin + offset)[i*N];
    }
    // increment and decrement
    my_type & operator++ ()
    { // increment then return
      offset += N;
      return *this;
    }
    my_type & operator-- ()
    { // decrement then return
      offset -= N;
      return *this;
    }
    my_type operator++ (int)
    { // copy, increment this, return copy
      my_type ret(*this);
      ++*this;
      return ret;
    }
    my_type operator-- (int)
    { // copy, decrement this, return copy
      my_type ret(*this);
      --*this;
      return ret;
    }
    // Arithmetics using difference_type
    my_type & operator+= (difference_type d)
    {
      offset += d*N;
      return *this;
    }
    my_type & operator-= (difference_type d)
    {
      return (*this += -d);
    }
    my_type operator+ (difference_type d) const
    {
      my_type ret(*this);
      ret += d;
      return ret;
    }
    my_type operator- (difference_type d) const
    {
      return this->operator+(-d);
    }

    // Arithmetics with this type
    difference_type operator- (my_type const & rhs) const
    {
      return (offset - rhs.offset) / N;
    }

    // Comparison
    bool operator== (my_type const &rhs) const { return /*m_begin == rhs.m_begin &&*/ offset == rhs.offset; }
    bool operator!= (my_type const &rhs) const { return !operator==(rhs); }
    bool operator<  (my_type const &rhs) const { return offset < rhs.offset; }
    bool operator>= (my_type const &rhs) const { return !operator<(rhs); }
    bool operator>  (my_type const &rhs) const { return offset > rhs.offset; }
    bool operator<= (my_type const &rhs) const { return !operator>(rhs); }

  public:

    difference_type N, offset;
    Iterator m_begin;

  };

  /*
  
  Iterator counting increases

  */

  template<class T>
  class count_iterator
  {

    T value;
    std::size_t offset;

    using my_type = count_iterator<T>;

  public:

    using iterator_category = std::random_access_iterator_tag;
    using reference = typename std::remove_reference<T>::type&;
    using pointer = typename std::remove_reference<T>::type*;
    using value_type = T;
    using difference_type = std::ptrdiff_t;

    count_iterator(T& val, std::size_t const i = 0)
      : value(val), offset(i) {}

    reference operator* () const
    {
      return value;
    }
    pointer operator->() const
    {
      return std::addressof(value);
    }
    reference operator[] (difference_type const) const
    {
      return value;
    }
    // increment and decrement
    my_type & operator++ ()
    { // increment then return
      ++offset;
      return *this;
    }
    my_type & operator-- ()
    { // decrement then return
      --offset;
      return *this;
    }
    my_type operator++ (int)
    { // copy, increment this, return copy
      my_type ret(*this);
      ++*this;
      return ret;
    }
    my_type operator-- (int)
    { // copy, decrement this, return copy
      my_type ret(*this);
      --*this;
      return ret;
    }
    // Arithmetics using difference_type
    my_type & operator+= (difference_type d)
    {
      offset += d;
      return *this;
    }
    my_type & operator-= (difference_type d)
    {
      return (*this += -d);
    }
    my_type operator+ (difference_type d) const
    {
      my_type ret(*this);
      ret += d;
      return ret;
    }
    my_type operator- (difference_type d) const
    {
      return this->operator+(-d);
    }

    // Arithmetics with this type
    difference_type operator- (my_type const & rhs) const
    {
      return offset - rhs.offset;
    }

    // Comparison
    bool operator== (my_type const &rhs) const
    {
      return std::addressof(value) == std::addressof(rhs.value)
        && offset == rhs.offset;
    }
    bool operator!= (my_type const &rhs) const { return !operator==(rhs); }
    bool operator<  (my_type const &rhs) const { return offset < rhs.offset; }
    bool operator>= (my_type const &rhs) const { return !operator<(rhs); }
    bool operator>  (my_type const &rhs) const { return offset > rhs.offset; }
    bool operator<= (my_type const &rhs) const { return !operator>(rhs); }

  };


  template<class _Matrix, class _Ref, class _Ptr, 
    bool ROWFIX = true, bool DIAG = false>
  class matrix_vector_iterator
  {

    typedef matrix_vector_iterator<_Matrix, _Ref, _Ptr, ROWFIX, DIAG> my_type;
    typedef matrix_vector_iterator < typename std::remove_const<_Matrix>::type,
      typename std::remove_const<_Ref>::type,
      typename std::remove_const<_Ptr>::type,
      ROWFIX, DIAG > non_const_type;

  public:

    typedef std::random_access_iterator_tag iterator_category;
    typedef typename _Matrix::value_type value_type;
    typedef typename _Matrix::difference_type difference_type;
    typedef difference_type distance_type;	// retained
    typedef _Ptr pointer;
    typedef _Ref reference;

    matrix_vector_iterator() : fixed(), free(), m_matrix() { }
    matrix_vector_iterator(_Matrix & m, difference_type const fix_id, difference_type const free_id = difference_type(0))
      : fixed(fix_id), free(DIAG?fix_id:free_id), m_matrix(&m)
    { }

    template<class OtherRef, class OtherPtr>
    matrix_vector_iterator(matrix_vector_iterator<typename std::remove_const<_Matrix>::type,
      OtherRef, OtherPtr, ROWFIX, DIAG> const & other)
      :fixed(other.fixed), free(other.free), m_matrix(other.m_matrix)
    { }


    /*
    Operators required for random_access_iterator category
    */
    reference operator* () const { return ROWFIX?(*m_matrix)(fixed, free):(*m_matrix)(free, fixed); }
    pointer   operator->() const { return &(this->operator*()); }
    reference operator[] (difference_type const & i) const
    {
      return ROWFIX?(*m_matrix)(fixed, free + i):(*m_matrix)(free + i, fixed);
    }

    // increment and decrement
    my_type & operator++ ()
    { // increment then return
      ++free;
      if (DIAG) ++fixed;
      return *this;
    }
    my_type & operator-- ()
    { // decrement then return
      --free;
      if (DIAG) --fixed;
      return *this;
    }
    my_type operator++ (int)
    { // copy, increment this, return copy
      my_type ret(*this);
      ++*this;
      return ret;
    }
    my_type operator-- (int)
    { // copy, decrement this, return copy
      my_type ret(*this);
      --*this;
      return ret;
    }
    // Arithmetics using difference_type
    my_type & operator+= (difference_type d)
    {
      free += d;
      if (DIAG) fixed += d;
      return *this;
    }
    my_type & operator-= (difference_type d)
    {
      return (*this += -d);
    }
    my_type operator+ (difference_type d) const
    {
      my_type ret(*this);
      ret += d;
      return ret;
    }
    my_type operator- (difference_type d) const
    {
      return this->operator+(-d);
    }

    // Arithmetics with this type
    difference_type operator- (my_type const & rhs) const
    {
      return free - rhs.free;
    }

    // Comparison
    bool operator== (my_type const &rhs) const
    {
      return m_matrix == rhs.m_matrix && fixed == rhs.fixed && free == rhs.free;
    }
    bool operator!= (my_type const &rhs) const { return !operator==(rhs); }
    bool operator<  (my_type const &rhs) const { return free < rhs.free; }
    bool operator>= (my_type const &rhs) const { return !operator<(rhs); }
    bool operator>  (my_type const &rhs) const { return free > rhs.free; }
    bool operator<= (my_type const &rhs) const { return !operator>(rhs); }

  public:

    difference_type fixed, free;
    _Matrix * m_matrix;

  };

  template<class ITER, class SIZE = std::size_t>
  class const_range_proxy
  {

    typedef const_range_proxy<ITER, SIZE> my_type;
    const_range_proxy& operator= (const_range_proxy const &);

  public:

    typedef ITER                               iterator;
    typedef SIZE                               size_type;
    typedef std::reverse_iterator<iterator>    reverse_iterator;
    typedef typename iterator::reference       reference;
    typedef typename iterator::pointer         pointer;
    typedef typename iterator::value_type      value_type;
    typedef typename iterator::difference_type difference_type;

  protected:

    iterator const m_begin;
    size_type const m_size;

  public:

    const_range_proxy(ITER const & begin, size_type const size)
      : m_begin(begin), m_size(size)
    { }

    iterator begin() const { return m_begin; }
    iterator end() const { return begin() + m_size; }

    iterator cbegin() const { return m_begin; }
    iterator cend() const { return end(); }

    reverse_iterator rbegin() const { return reverse_iterator(end()); }
    reverse_iterator rend() const { return reverse_iterator(begin()); }

    reverse_iterator crbegin() const { return reverse_iterator(cend()); }
    reverse_iterator crend() const { return reverse_iterator(cbegin()); }

    // content check
    size_type size() const { return m_size; }
    bool empty() const { return size() == 0; }

    reference front() const { return *begin(); }
    reference back() const { return *(end() - 1); }

    reference at(size_type const index) const
    {
      if (index >= size())
      {
        throw std::out_of_range("const_range_proxy::at");
      }
      return *(begin() + index);
    }

    reference operator[] (size_type const index) const { return *(begin() + index); }

    bool operator== (my_type const &other) const
    {
      using std::equal;
      if (size() != other.size()) return false;
      return equal(begin(), end(), other.begin());
    }
    bool operator!= (my_type const &other) const
    {
      return !(other == *this);
    }

  };


  template<class ITER, class CONST_ITER, class SIZE = std::size_t>
  class range_proxy
  {

    typedef range_proxy<ITER, CONST_ITER, SIZE> my_type;

  public:

    typedef ITER                                     iterator;
    typedef CONST_ITER                               const_iterator;
    typedef SIZE                                     size_type;
    typedef std::reverse_iterator<iterator>          reverse_iterator;
    typedef std::reverse_iterator<const_iterator>    const_reverse_iterator;
    typedef typename iterator::reference             reference;
    typedef typename const_iterator::reference       const_reference;
    typedef typename iterator::pointer               pointer;
    typedef typename const_iterator::pointer         const_pointer;
    typedef typename iterator::value_type            value_type;
    typedef typename const_iterator::value_type      const_value_type;
    typedef typename iterator::difference_type       difference_type;
    typedef typename const_iterator::difference_type const_difference_type;


  protected:

    iterator const m_begin;
    size_type const m_size;

    template<class X>
    void assign_from(X && rhs)
    {
      using std::begin;
      using std::end;
      auto first = begin(rhs);
      auto last = end(rhs);
      if (size() != static_cast<size_type>(std::distance(first, last)))
        throw std::logic_error("range_proxy assignment size mismatch");
      // copy if reference is given (lvalue, const value)
      // move if rvalue is passed
      if (std::is_reference<X>::value)
        std::copy(first, last, this->begin());
      else
        std::move(first, last, this->begin());
    }

  public:

    range_proxy(iterator const & iter,
      size_type const & size)
      : m_begin(iter), m_size(size)
    { }

    range_proxy(iterator const & first,
      iterator const & last)
      : m_begin(first), m_size(static_cast<size_type>(
        std::distance(first,last)))
    { }

    range_proxy(range_proxy const &) = default;
    range_proxy(range_proxy &&) = default;

    template<class X>
    range_proxy& operator= (X && rhs)
    {
      assign_from(std::forward<X>(rhs));
      return *this;
    }

    range_proxy& operator= (range_proxy const & rhs)
    { // call templated version explicitly
      assign_from(rhs);
      return *this;
    }

    range_proxy& operator= (range_proxy && rhs)
    { // call templated version explicitly
      assign_from(std::move(rhs));
      return *this;
    }

    // range
    iterator begin() { return m_begin; }
    iterator end() { return begin() + m_size; }

    const_iterator begin() const { return const_iterator(m_begin); }
    const_iterator end() const { return const_iterator(begin() + m_size); }

    const_iterator cbegin() const { return const_iterator(begin()); }
    const_iterator cend() const { return const_iterator(end()); }

    reverse_iterator rbegin() { return reverse_iterator(end()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }

    const_reverse_iterator rbegin() const { return const_reverse_iterator(cend()); }
    const_reverse_iterator rend() const { return const_reverse_iterator(cbegin()); }

    const_reverse_iterator crbegin() const { return const_reverse_iterator(cend()); }
    const_reverse_iterator crend() const { return const_reverse_iterator(cbegin()); }

    // front 
    reference front() { return *begin(); }
    reference back() { return *(end() - 1); }
    // back
    const_reference back() const { return *(cend() - 1); }
    const_reference front() const { return *cbegin(); }

    // size / content
    size_type size() const { return m_size; }
    bool empty() const { return size() == 0; }

    reference at(size_type const index)
    {
      if (index >= size())
      {
        throw std::out_of_range("range_proxy::at");
      }
      return *(begin() + index);
    }
    const_reference at(size_type const index) const
    {
      if (index >= size())
      {
        throw std::out_of_range("range_proxy::at");
      }
      return *(cbegin() + index);
    }

    reference operator[] (size_type const index)
    {
      return *(begin() + index);
    }
    const_reference operator[] (size_type const index) const
    {
      return *(cbegin() + index);
    }

    void swap(my_type &other)
    {
      if (size() != other.size())
        throw std::out_of_range("Swapping range proxys"
          " requires equal sized ranges");
      std::swap_ranges(this->begin(), this->end(), other.begin());
    }

    bool operator== (my_type const &other) const
    {
      if (size() != other.size()) return false;
      return std::equal(begin(), end(), other.begin());
    }
    bool operator!= (my_type const &other) const
    {
      return !(other == *this);
    }

  };

  template<class ITER, class CONST_ITER, class SIZE = std::size_t>
  struct range_proxy_traits
  {
    typedef ITER iterator;
    typedef CONST_ITER const_iterator;
    typedef range_proxy<iterator, const_iterator, SIZE> range;
    typedef const_range_proxy<const_iterator, SIZE> const_range;
  };

  struct index_iterator
  {

    using my_type = index_iterator;

    std::size_t _val;

  public:

    using difference_type = std::ptrdiff_t;
    using iterator_category = std::random_access_iterator_tag;
    using pointer = std::size_t const *;
    using reference = std::size_t const &;
    using value_type = std::size_t;

    index_iterator(std::size_t initial_value = 0) : _val(initial_value) {}

    std::size_t operator* () const
    {
      return _val;
    }
    pointer operator->() const
    {
      return &_val;
    }
    std::size_t operator[] (difference_type const & i) const
    {
      return _val + i;
    }
    // increment and decrement
    my_type & operator++ ()
    { // increment then return
      ++_val;
      return *this;
    }
    my_type & operator-- ()
    { // decrement then return
      --_val;
      return *this;
    }
    my_type operator++ (int)
    { // copy, increment this, return copy
      my_type ret(*this);
      ++*this;
      return ret;
    }
    my_type operator-- (int)
    { // copy, decrement this, return copy
      my_type ret(*this);
      --*this;
      return ret;
    }
    // Arithmetics using difference_type
    my_type & operator+= (difference_type d)
    {
      _val += d;
      return *this;
    }
    my_type & operator-= (difference_type d)
    {
      return (*this += -d);
    }
    my_type operator+ (difference_type d) const
    {
      my_type ret(*this);
      ret += d;
      return ret;
    }
    my_type operator- (difference_type d) const
    {
      return this->operator+(-d);
    }

    // Arithmetics with this type
    difference_type operator- (my_type const & rhs) const
    {
      return _val - rhs._val;
    }

    // Comparison
    bool operator== (my_type const &rhs) const
    {
      return _val == rhs._val;
    }
    bool operator!= (my_type const &rhs) const { return !operator==(rhs); }

    bool operator<  (my_type const &rhs) const { return _val < rhs._val; }
    bool operator>= (my_type const &rhs) const { return !operator<(rhs); }
    bool operator>  (my_type const &rhs) const { return rhs._val < _val; }
    bool operator<= (my_type const &rhs) const { return !operator>(rhs); }

  };


}

#endif
