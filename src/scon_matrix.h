#if !defined (SCON_MATRIX_HEADER)
#define SCON_MATRIX_HEADER

#include <vector>
#include <algorithm>
#include <cstddef>

#include "scon_traits.h"
#include "scon_iterator.h"
#include "scon_vect.h"


namespace scon
{

  // Tringular Index

  template <bool diag> inline std::size_t triangularIndex(std::size_t const a, std::size_t const b);

  template<> inline std::size_t triangularIndex<true>(std::size_t const a, std::size_t const b)
  {
    return a > b ? (a*a + a) / 2 + b : (b*b + b) / 2 + a;
  }
  template<> inline std::size_t triangularIndex<false>(std::size_t const a, std::size_t const b)
  {
    return a > b ? (a*a - a) / 2 + b : (b*b - b) / 2 + a;
  }

  /*

  DYNAMIC MATRIX

  */

  template<class T>
  using std_vector_wrapper = std::vector < T > ;

  template<class T, bool SYMMETRIC = false, template<class...> class ContainerT = std_vector_wrapper>
  class matrix
  {

  public: // Basic type interface

    typedef matrix<T, SYMMETRIC, ContainerT>                  my_type;
    typedef ContainerT<T>                                     container_type;
    typedef typename container_type::value_type               value_type;
    typedef typename container_type::size_type                size_type;
    typedef typename container_type::difference_type          difference_type;
    typedef typename container_type::reference                reference;
    typedef typename container_type::const_reference          const_reference;
    typedef typename container_type::pointer                  pointer;
    typedef typename container_type::const_pointer            const_pointer;
    typedef typename container_type::iterator                 iterator;
    typedef typename container_type::const_iterator           const_iterator;
    typedef typename container_type::reverse_iterator         reverse_iterator;
    typedef typename container_type::const_reverse_iterator   const_reverse_iterator;

  protected: // vector proxies (row, col, diagonals)

    // Proxy Types
    typedef range_proxy_traits < iterator,
      const_iterator, size_type > _row_proxy;
    typedef range_proxy_traits < stride_iterator<iterator>,
      stride_iterator<const_iterator>, size_type > _col_proxy;
    typedef range_proxy_traits < stride_iterator<iterator>,
      stride_iterator<const_iterator>, size_type > _diag_proxy;

  public: // vector types (row, col, diagonals)

    // storage is row major
    // row iterator is the container iterator
    typedef typename _row_proxy::range         row_type;
    typedef typename _row_proxy::const_range   const_row_type;
    // storage is row major !
    // column iteration is strided by number of columns
    typedef typename _col_proxy::range         col_type;
    typedef typename _col_proxy::const_range   const_col_type;
    // storage is row major !
    // diagonal iteration is strided by number of columns + 1
    typedef typename _diag_proxy::range        diag_type;
    typedef typename _diag_proxy::const_range  const_diag_type;

  protected: // data

    container_type m_data;
    std::size_t ROWS, COLS;

  public:  // function interface

    /* constructors */

    matrix()
      : m_data(), ROWS(0), COLS(0)
    { }

    matrix(size_type const rows,
      size_type const cols, T const &val = T())
      : m_data(rows*cols, val), ROWS(rows), COLS(cols)
    { }

    // Ownership transfer

    void swap(my_type & rhs)
    {
      if (this != &rhs)
      {
        using std::swap;
        swap(ROWS, rhs.ROWS);
        swap(COLS, rhs.COLS);
        swap(m_data, rhs.m_data);
      }
    }

    // Size interface

    size_type size() const { return m_data.size(); }
    bool empty() const { return m_data.empty(); }
    size_type rows() const { return ROWS; }
    size_type cols() const { return COLS; }

    void clear()
    {
      ROWS = COLS = 0u;
      m_data.clear();
    }

    void resize(size_type const rowcount,
      size_type const colcount, T const &val = T())
    {
      using std::min;
      if (rowcount == ROWS && colcount == COLS) return;
      // RAII -> get new matrix
      my_type tmp(rowcount, colcount, val);
      // Copy data
      size_type const copy_row_num(min(rowcount, ROWS));

      //Copy content
      if (colcount == COLS)
      {
        for (std::size_t i = 0; i < copy_row_num; ++i)
        {
          tmp.row(i) = row(i);
        }
      }
      else
      {
        for (std::size_t j = 0; j < min(colcount, COLS); ++j)
        {
          for (std::size_t i = 0; i < copy_row_num; ++i)
          {
            tmp(i,j) = (*this)(i,j);
          }
        }
      }

      // transfer ownership
      tmp.swap(*this);
    }

    // range for whole matrix

    iterator begin() { return m_data.begin(); }
    const_iterator begin() const { return m_data.begin(); }
    const_iterator cbegin() const { return m_data.cbegin(); }
    iterator end() { return m_data.end(); }
    const_iterator end() const { return m_data.end(); }
    const_iterator cend() const { return m_data.cend(); }

    reverse_iterator rbegin() { return m_data.rbegin(); }
    const_reverse_iterator rbegin() const { return m_data.rbegin(); }
    const_reverse_iterator crbegin() const { return m_data.crbegin(); }
    reverse_iterator rend() { return m_data.rend(); }
    const_reverse_iterator rend() const { return m_data.rend(); }
    const_reverse_iterator crend() const { return m_data.crend(); }

    // 2d access
    reference operator() (size_type const i) { return m_data[i]; }
    const_reference operator() (size_type const i) const { return m_data[i]; }
    reference operator() (size_type const r, size_type const c) { return m_data[COLS*r + c]; }
    const_reference operator() (size_type const r, size_type const c) const { return m_data[COLS*r + c]; }
    // 2d access with  
    reference at(size_type const i) { return m_data.at(i); }
    const_reference at(size_type const i) const { return m_data.at(i); }
    reference at (size_type const r, size_type const c) { return m_data.at(COLS*r + c); }
    const_reference at(size_type const r, size_type const c) const { return m_data.at(COLS*r + c); }

    // Matrix vectors

    // Rows
    row_type row(size_type const index)
    {
#if defined (SCON_DEBUG)
      if (index >= ROWS) throw std::out_of_range("matrix::row : index out of range");
#endif
      return row_type(typename row_type::iterator(m_data.begin() + COLS*index), COLS);
    }
    const_row_type row(size_type const index) const
    {
#if defined (SCON_DEBUG)
      if (index >= ROWS) throw std::out_of_range("matrix::row : index out of range");
#endif
      return const_row_type(typename const_row_type::iterator(m_data.cbegin() + COLS*index), COLS);
    }

    row_type operator[](size_type const index)
    {
#if defined (SCON_DEBUG)
      if (index >= ROWS) throw std::out_of_range("matrix::row : index out of range");
#endif
      return row_type(typename row_type::iterator(m_data.begin() + COLS*index), COLS);
    }
    const_row_type operator[](size_type const index) const
    {
#if defined (SCON_DEBUG)
      if (index >= ROWS) throw std::out_of_range("matrix::row : index out of range");
#endif
      return const_row_type(typename const_row_type::iterator(m_data.cbegin() + COLS*index), COLS);
    }

    // Columns
    col_type col(size_type const index)
    {
#if defined (SCON_DEBUG)
      if (index > COLS) throw std::out_of_range("matrix::col : index out of range");
#endif
      return col_type(typename col_type::iterator(m_data.begin() + index, COLS), ROWS);
    }
    const_col_type col(size_type const index) const
    {
#if defined (SCON_DEBUG)
      if (index > COLS) throw std::out_of_range("matrix::col : index out of range");
#endif
      return const_col_type(typename const_col_type::iterator(m_data.cbegin() + index, COLS), ROWS);
    }

    diag_type diag()
    {
#if defined (SCON_DEBUG)
      if (COLS != ROWS) throw std::logic_error("matrix::diag : matrix not square");
#endif
      return diag_type(typename diag_type::iterator(m_data.begin(), COLS + 1), ROWS);
    }
    const_diag_type diag() const
    {
#if defined (SCON_DEBUG)
      if (ROWS != COLS) throw std::logic_error("matrix::diag : matrix not square");
#endif
      return const_diag_type(typename const_diag_type::iterator(m_data.cbegin(), COLS + 1), ROWS);
    }

  };

  

  template<class T, template<class...> class ContainerT>
  class matrix<T, true, ContainerT>
  {

  public: // Type interface

    typedef matrix<T, true, ContainerT>                         my_type;
    typedef ContainerT<T>                                       container_type;
    typedef typename container_type::value_type                 value_type;
    typedef typename container_type::size_type                  size_type;
    typedef typename container_type::difference_type            difference_type;
    typedef typename container_type::reference                  reference;
    typedef typename container_type::const_reference            const_reference;
    typedef typename container_type::pointer                    pointer;
    typedef typename container_type::const_pointer              const_pointer;
    typedef typename container_type::iterator                   iterator;
    typedef typename container_type::const_iterator             const_iterator;
    typedef typename container_type::reverse_iterator           reverse_iterator;
    typedef typename container_type::const_reverse_iterator     const_reverse_iterator;

  protected: // vector proxies (row, col, diagonals)

    typedef matrix_vector_iterator
      <my_type, reference, pointer, true, false>                    _row_iter;
    typedef matrix_vector_iterator
      <const my_type, const_reference, const_pointer, true, false>  _const_row_iter;
    typedef matrix_vector_iterator
      <my_type, reference, pointer, false, false>                   _col_iter;
    typedef matrix_vector_iterator
      <const my_type, const_reference, const_pointer, false, false> _const_col_iter;
    typedef matrix_vector_iterator
      <my_type, reference, pointer, false, true>                    _diag_iter;
    typedef matrix_vector_iterator
      <const my_type, const_reference, const_pointer, false, true>  _const_diag_iter;

    // Proxy Types
    typedef range_proxy_traits<_row_iter, _const_row_iter, 
      size_type>   _row_proxy;
    typedef range_proxy_traits<_col_iter, _const_col_iter, 
      size_type>   _col_proxy;
    typedef range_proxy_traits<_diag_iter, _const_diag_iter,
      size_type>   _diag_proxy;

  public: // vector types (row, col, diagonals)

    // storage is row major
    // row iterator is the container iterator
    typedef typename _row_proxy::range        row_type;
    typedef typename _row_proxy::const_range  const_row_type;
    // storage is row major !
    // column iteration is strided by number of columns
    typedef typename _col_proxy::range        col_type;
    typedef typename _col_proxy::const_range  const_col_type;
    // storage is row major !
    // diagonal iteration is strided by number of columns + 1
    typedef typename _diag_proxy::range       diag_type;
    typedef typename _diag_proxy::const_range const_diag_type;

  protected: // data

    container_type m_data;
    std::size_t SIZE;

    // helpers

    static inline size_type sysize(size_type const & w)
    {
      return w*(w+1) / 2;
    }

    static inline size_type offset_sorted(size_type const &r, 
      size_type const &c)
    {
      return sysize(r) + c;
    }

    static inline size_type offset(size_type const &x, 
      size_type const &y)
    {
      return x >= y ? offset_sorted(x, y) : offset_sorted(y, x);
    }

  public:  // interface

    matrix()
      : m_data(), SIZE()
    { }

    matrix(size_type const width, value_type const &val = T())
      : m_data(sysize(width), val), SIZE(width)
    { }

    void swap(my_type & rhs)
    {
      if (this != &rhs)
      {
        using std::swap;
        swap(SIZE, rhs.SIZE);
        swap(m_data, rhs.m_data);
      }
    }

    // size interface

    size_type size() const { return m_data.size(); }
    bool empty() const { return m_data.empty(); }
    size_type rows() const { return SIZE; }
    size_type cols() const { return SIZE; }

    void clear()
    {
      SIZE = 0;
      m_data.clear();
    }

    void resize(size_type const new_size, T const &val = T())
    {
      if (new_size == SIZE) return;
      // RAII -> get new matrix
      my_type tmp(new_size, val);
      // Copy data
      size_type const min_size(std::min(size(), tmp.size()));
      std::copy(begin(), begin() + min_size, tmp.begin());
      // transfer ownership
      tmp.swap(*this);
    }

    // range for whole matrix

    iterator begin() { return m_data.begin(); }
    const_iterator begin() const { return m_data.begin(); }
    const_iterator cbegin() const { return m_data.cbegin(); }
    iterator end() { return m_data.end(); }
    const_iterator end() const { return m_data.end(); }
    const_iterator cend() const { return m_data.cend(); }

    reverse_iterator rbegin() { return m_data.rbegin(); }
    const_reverse_iterator rbegin() const { return m_data.rbegin(); }
    const_reverse_iterator crbegin() const { return m_data.crbegin(); }
    reverse_iterator rend() { return m_data.rend(); }
    const_reverse_iterator rend() const { return m_data.rend(); }
    const_reverse_iterator crend() const { return m_data.crend(); }

    reference operator() (size_type const i) { return m_data[i]; }
    const_reference operator() (size_type const i) const { return m_data[i]; }
    reference operator() (size_type const r, size_type const c) { return m_data[offset(r, c)]; }
    const_reference operator() (size_type const r, size_type const c) const { return m_data[offset(r, c)]; }

    reference at(size_type const i) { return m_data.at(i); }
    const_reference at(size_type const i) const { return m_data.at(i); }
    reference at(size_type const r, size_type const c) { return m_data.at(offset(r,c)); }
    const_reference at(size_type const r, size_type const c) const { return m_data.at(offset(r, c)); }

    // Matrix vectors

    // Rows
    row_type row(size_type const index)
    {
#if defined (SCON_DEBUG)
      if (index >= SIZE) throw std::out_of_range("matrix::row : index out of range");
#endif
      return row_type(typename row_type::iterator(*this, index), SIZE);
    }
    const_row_type row(size_type const index) const
    {
#if defined (SCON_DEBUG)
      if (index >= SIZE) throw std::out_of_range("matrix::row : index out of range");
#endif
      return const_row_type(typename const_row_type::iterator(*this, index), SIZE);
    }

    row_type operator[](size_type const index)
    {
#if defined (SCON_DEBUG)
      if (index >= SIZE) throw std::out_of_range("matrix::row : index out of range");
#endif
      return row_type(typename row_type::iterator(*this, index), SIZE);
    }
    const_row_type operator[](size_type const index) const
    {
#if defined (SCON_DEBUG)
      if (index >= SIZE) throw std::out_of_range("matrix::row : index out of range");
#endif
      return const_row_type(typename const_row_type::iterator(*this, index), SIZE);
    }

    // Columns
    col_type col(size_type const index)
    {
#if defined (SCON_DEBUG)
      if (index > SIZE) throw std::out_of_range("matrix::col : index out of range");
#endif
      return col_type(typename col_type::iterator(*this, index), SIZE);
    }
    const_col_type col(size_type const index) const
    {
#if defined (SCON_DEBUG)
      if (index > SIZE) throw std::out_of_range("matrix::col : index out of range");
#endif
      return const_col_type(typename const_col_type::iterator(*this, index), SIZE);
    }

    diag_type diag()
    {
      return diag_type(typename diag_type::iterator(*this, 0, 0), SIZE);
    }
    const_diag_type diag() const
    {
      return const_diag_type(typename const_diag_type::iterator(*this, 0, 0), SIZE);
    }

  };



  /*

    Equality checks

  */


  template<class T, bool S1, bool S2>
  bool operator== (matrix<T, S1> const & lhs, matrix<T, S2> const & rhs)
  {
    typedef typename matrix<T, true>::size_type size_type;
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) return false;
    for (size_type i = 0; i < lhs.rows(); ++i)
    {
      for (size_type j = 0; j <= lhs.cols(); ++j)
      {
        if (lhs(i, j) != rhs(i, j)) return false;
      }
    }
    return true;
  }

  template<class T, bool S1, bool S2>
  bool operator!= (matrix<T, S1> const & lhs, matrix<T, S2> const & rhs)
  {
    return !(lhs == rhs);
  }


  template<class T, bool S>
  bool operator== (matrix<T, S> const & lhs, matrix<T, S> const & rhs)
  {
    if (lhs.size() != rhs.size()) return false;
    return std::equal(lhs.begin(), lhs.end(), rhs.begin());
  }

  template<class T, bool S>
  bool operator!= (matrix<T, S> const & lhs, matrix<T, S> const & rhs)
  {
    return !(lhs == rhs);
  }



  /*
  
    STREAM OUTPUT
  
  */




  template<class Matrix>
  struct matrix_printer
  {
  private:
    matrix_printer& operator= (matrix_printer const &);
  public:
    Matrix const & matrix;
    std::streamsize margin_left, margin_right;
    std::size_t margin_top, margin_bottom;
    char fill_left, fill_right, spacer;
    bool print_indices;
    matrix_printer(Matrix const &m) :
      matrix(m),
      margin_left(0), margin_right(0), margin_top(0), margin_bottom(0),
      fill_left('\0'), fill_right('\0'), spacer('\0'),
      print_indices(false)
    {}
    matrix_printer operator() (Matrix const &m) const
    {
      matrix_printer tmp(m);
      tmp.margin_left = margin_left;
      tmp.margin_right = margin_right;
      tmp.margin_top = margin_top;
      tmp.margin_bottom = margin_bottom;
      tmp.fill_left = fill_left;
      tmp.fill_right = fill_right;
      tmp.spacer = spacer;
      tmp.print_indices = print_indices;
      return tmp;
    }
    void stream_fill_left(std::ostream & strm) const
    {
      if (fill_left != '\0')
        fill(strm, fill_left, margin_left);
    }
    void stream_fill_right(std::ostream & strm) const
    {
      if (fill_right != '\0')
        fill(strm, fill_right, margin_left);
    }
  private:
    void fill(std::ostream & strm, char const fillc, std::streamsize const & width) const
    {
      char const x = strm.fill();
      strm.fill(fillc);
      strm << std::setw(width) << fill_left;
      strm.fill(x);
    }
  };

  template<class M>
  inline std::ostream & operator<< (std::ostream & strm, matrix_printer<M> const &printer)
  {
    // save width for every element
    std::streamsize const W = strm.width();
    // Margin top
    for (std::size_t i = 0; i < printer.margin_top; ++i) strm << '\n';
    // Iterate matrix
    typedef typename M::size_type size_type;
    typedef typename M::row_type::const_iterator row_iterator;

    for (size_type i = 0; i < printer.matrix.rows(); ++i)
    {
      row_iterator const end = printer.matrix.row(i).cend();
      printer.stream_fill_left(strm);
      for (row_iterator it = printer.matrix.row(i).cbegin(); it != end; ++it)
      {
        strm << std::setw(W) << *it;
        if (printer.spacer != '\0') strm << printer.spacer;
      }
      printer.stream_fill_right(strm);
      strm << '\n';
    }
    // Margin bottom
    for (std::size_t i = 0; i < printer.margin_top; ++i) strm << '\n';
    // return stream
    return strm;
  }

  template<class T, bool B>
  inline std::ostream & operator<< (std::ostream & strm, scon::matrix<T,B> const &matrix)
  {
    typedef matrix_printer< scon::matrix<T,B> > p_type;
    strm << p_type(matrix);
    return strm;
  }

}


#endif
