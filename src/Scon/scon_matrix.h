#if !defined (SCON_MATRIX_HEADER)
#define SCON_MATRIX_HEADER

#include <vector>
#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <iterator>

#include "scon_traits.h"
#include "scon_iterator.h"
#include "scon_vect.h"


#if !defined(SCON_RESTRICT) && defined(_MSC_VER)
#define SCON_RESTRICT __restrict
#elif !defined(SCON_RESTRICT) && defined(__GNUG__)
#define SCON_RESTRICT __restrict__
#endif

namespace scon
{

	// Tringular Index

	template <bool diag> inline std::size_t triangularIndex(std::size_t const a, std::size_t const b);

	template<> inline std::size_t triangularIndex<true>(std::size_t const a, std::size_t const b)
	{
		return a > b ? (a * a + a) / 2 + b : (b * b + b) / 2 + a;
	}
	template<> inline std::size_t triangularIndex<false>(std::size_t const a, std::size_t const b)
	{
		return a > b ? (a * a - a) / 2 + b : (b * b - b) / 2 + a;
	}

	/*

	DYNAMIC MATRIX

	*/

	template<class T>
	using std_vector_wrapper = std::vector < T >;

	template<class T, bool SYMMETRIC = false,
		template<class...> class ContainerT = std_vector_wrapper>
	class matrix
	{

	public: // Basic type interface

		struct iterator_construct_type {};

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

		matrix() : m_data(), ROWS(0), COLS(0) { }

		matrix(size_type const rows,
			size_type const cols)
			: m_data(rows* cols), ROWS(rows), COLS(cols)
		{ }

		matrix(size_type const rows,
			size_type const cols, T const& val)
			: m_data(rows* cols, val), ROWS(rows), COLS(cols)
		{ }

		// create vector, row major
		matrix(std::size_t const in) :
			scon::matrix<T>(in, 1u) { };

		// Use iterator range to form diagonals of new matrix
		template<class InputIt>
		matrix(InputIt first, InputIt last, iterator_construct_type)
			: m_data(static_cast<std::size_t>(std::distance(first, last)),
				static_cast<std::size_t>(std::distance(first, last))),
			ROWS(static_cast<std::size_t>(std::distance(first, last))),
			COLS(static_cast<std::size_t>(std::distance(first, last)))
		{
			std::size_t index{ 0u };
			std::for_each(first, last,
				[this, &index](decltype(*first) const & val) -> void
			{
				(*this)(index, index) = val;
				++index;
			});
		}

		// construct from vector of vectors (with size checks)
		matrix(std::vector<std::vector<T>> const& in)
		{
			if (in.empty()) return;
			auto const n = in.size(), // rows
				m = in.front().size(); // cols
			bool same_size{ true };
			for (std::size_t i = 1u; i < n; ++i)
			{
				if (in[i].size() != m)
					throw std::runtime_error("Size mismatch "
						"in matrix construction.");
			}
			// resize matrix
			resize(n, m);
			// copy every submatrix into matrix rows
			for (std::size_t r = 0u; r < n; ++r)
			{
				std::copy(in[r].begin(), in[r].end(), row(r).begin());
			}
		}

		// put vector in diagonal: delegate to iterator constructor
		matrix(std::vector<T> const& input)
			: matrix(input.begin(), input.end(), iterator_construct_type{})
		{ }

		// create identity matrix
		static typename std::enable_if<std::is_arithmetic<T>::value, matrix>::type
			identity(std::size_t const num_rows, std::size_t const num_cols)
		{
			matrix r(num_rows, num_cols, T{});
			auto m = std::min(num_rows, num_cols);
			for (std::size_t i = 0; i < m; ++i)
			{
				r(i, i) = T{ 1 };
			}
			return r;
		}

		// Ownership transfer
		void swap(my_type& rhs)
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
			size_type const colcount, T const& val = T())
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
						tmp(i, j) = (*this)(i, j);
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
		reference operator() (size_type const r, size_type const c) { return m_data[COLS * r + c]; }
		const_reference operator() (size_type const r, size_type const c) const { return m_data[COLS * r + c]; }
		// 2d access with  
		reference at(size_type const i) { return m_data.at(i); }
		const_reference at(size_type const i) const { return m_data.at(i); }
		reference at(size_type const r, size_type const c) { return m_data.at(COLS * r + c); }
		const_reference at(size_type const r, size_type const c) const { return m_data.at(COLS * r + c); }

		// data access
		pointer data() { return m_data.data(); }
		const_pointer data() const { return m_data.data(); }

		// Matrix vectors

		// Rows
		row_type row(size_type const index)
		{
#if defined (SCON_DEBUG)
			if (index >= ROWS) throw std::out_of_range("matrix::row : index out of range");
#endif
			return row_type(typename row_type::iterator(m_data.begin() + COLS * index), COLS);
		}
		const_row_type row(size_type const index) const
		{
#if defined (SCON_DEBUG)
			if (index >= ROWS) throw std::out_of_range("matrix::row : index out of range");
#endif
			return const_row_type(typename const_row_type::iterator(m_data.cbegin() + COLS * index), COLS);
		}

		row_type operator[](size_type const index)
		{
#if defined (SCON_DEBUG)
			if (index >= ROWS) throw std::out_of_range("matrix::row : index out of range");
#endif
			return row_type(typename row_type::iterator(m_data.begin() + COLS * index), COLS);
		}
		const_row_type operator[](size_type const index) const
		{
#if defined (SCON_DEBUG)
			if (index >= ROWS) throw std::out_of_range("matrix::row : index out of range");
#endif
			return const_row_type(typename const_row_type::iterator(m_data.cbegin() + COLS * index), COLS);
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

		matrix& operator*= (T const& val)
		{
			// multiply ever element with val
			for (auto& e : *this)
			{
				e *= val;
			}
			return *this;
		}

	};

	/*

	math & utils

	*/

	// check whether matrix is a square matrix
	template<class T, template<class...> class C>
	bool is_square(matrix<T, false, C> const& m)
	{
		return m.rows() == m.cols();
	}


	// check whether matrix is a symmetric matrix
	template<class T, template<class...> class C>
	bool is_symmetric(matrix<T, false, C> const& m)
	{
		if (!is_square(m)) return false;
		auto const n = m.rows();
		for (std::size_t i = 0; i < n; ++i)
		{
			for (std::size_t j = 0; j < n; ++j)
			{
				if (m(i, j) != m(j, i)) return false;
			}
		}
		return true;
	}

	template<class T, template<class...> class C>
	scon::matrix<T> transposed(matrix<T, false, C> const& A)
	{
		auto const N = A.rows();
		auto const M = A.cols();
		auto B = scon::matrix<T>(M, N);
#if defined(_OPENMP)
		auto const oll = static_cast<std::ptrdiff_t>(N * M);
#pragma omp parallel for
		for (std::ptrdiff_t n = 0; n < oll; ++n)
#else
		for (std::size_t n = 0; n < N * M; ++n)
#endif
		{
			B(n) = A(M * (n % N) + (n / N));
		}
		return B;
	}

	// transpose matrix
	template<class T, template<class...> class C>
	void transpose(matrix<T, false, C>& m)
	{
		if (m.rows() == m.cols())
		{
			for (std::size_t r = 0; r < m.rows(); ++r)
			{
				for (std::size_t c = 0; c < r; ++c)
				{
					std::swap(m(r, c), m(c, r));
				}
			}
		}
		else
		{
			auto t = transposed(m);
			m.swap(t);
		}
	}

	// matrix multiplication
	template<class T, class U, template<class...> class C>
	typename std::enable_if < std::is_arithmetic<T>::value&& std::is_arithmetic<U>::value,
		scon::matrix<typename std::common_type<T, U>::type >> ::type
		operator*(scon::matrix<T, false, C> const& a, scon::matrix<U, false, C> const& b)
	{
		if (a.empty() || b.empty() || a.cols() != b.rows())
		{
			throw std::logic_error("Matrix multiplication dimension mismatch.");
		}
		// create result matrix r
		scon::matrix<T, false, C> r(a.rows(), b.cols(), 0.0);
		// make restricted pointer to const data of a
		T const* SCON_RESTRICT const ap = a.data();
		// make restricted pointer to const data of b
		//T const * SCON_RESTRICT const bp = b.data();;
		// make restricted pointer to data of r
		T* SCON_RESTRICT const rp = r.data();
#if defined(_OPENMP)
#pragma omp parallel
		{
#endif
			// create matrix copy inside parallel region: thread local
			scon::matrix<T> const my_b(b);
			// create restricted pointer to const data of new b
			T const* SCON_RESTRICT const bp = &my_b(0, 0);
#if defined(_OPENMP)
			// avoid signed unsigned mismatch in loop
			auto ar = static_cast<std::ptrdiff_t>(a.rows());
#pragma omp for 
			for (std::ptrdiff_t i = 0; i < ar; ++i)
#else
			for (std::size_t i = 0; i < a.rows(); ++i)
#endif
			{ // loop result rows
				// pointer to result_row
				T* SCON_RESTRICT const r_row_p = rp + i * b.cols();
				for (std::size_t k = 0; k < b.rows(); ++k)
				{ // loop b rows
					// current multiplicate value from a 
					auto const a_val = ap[i * a.cols() + k];
					// create pointer into the current row in b
					T const* SCON_RESTRICT const b_row_p = bp + k * b.cols();
					for (std::size_t j = 0; j < b.cols(); ++j)
					{ // loop result cols
						// add result to element in result row
						r_row_p[j] += a_val * b_row_p[j];
					}
				}
			}
#if defined(_OPENMP)
		}
#endif
		return r;
	}

	namespace detail
	{

		template<class T, bool B, template<class...> class C, class It>
		typename std::enable_if < std::is_arithmetic<T>::value,
			C<T>>::type
			matrix_mul_range(scon::matrix<T, B, C> const& A,
				It col_first, It col_last)
		{
			auto const n = A.rows();
			auto const m = A.cols();
			auto const rl = static_cast<std::size_t>(
				std::distance(col_first, col_last));
			if (rl != m)
				throw std::out_of_range("Improper sized column range for multiplication.");
			auto r = C<T>(n, T{});
			for (std::size_t i = 0; i < n; ++i)
			{
				for (std::size_t j = 0; j < m; ++j)
				{
					r[i] += A(i, j) * *(col_first + j);
				}
			}
			return r;
		}

		template<class T> struct is_matrix : std::false_type {};

		template<class T, bool b, template<class...> class C>
		struct is_matrix < matrix<T, b, C> > : std::true_type {};

	}

	// matrix vector multiplication
	template<class T, bool B, template<class...> class C, class U>
	typename std::enable_if <
		std::is_arithmetic<T>::value&& scon::is_range<U>::value
		&& !detail::is_matrix<typename std::remove_reference<U>::type>::value
		&& !std::is_base_of<matrix<T, B, C>, typename std::remove_reference<U>::type>::value &&
		!std::is_base_of<matrix<range_value<U>>, U>::value, C<T>>::type
		operator*(scon::matrix<T, B, C> const& a, U && b)
	{
		using std::begin;
		using std::end;
		return detail::matrix_mul_range(a, begin(b), end(b));
	}

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

		static inline size_type sysize(size_type const& w)
		{
			return w * (w + 1) / 2;
		}

		static inline size_type offset_sorted(size_type const& r,
			size_type const& c)
		{
			return sysize(r) + c;
		}

		static inline size_type offset(size_type const& x,
			size_type const& y)
		{
			return x >= y ? offset_sorted(x, y) : offset_sorted(y, x);
		}

	public:  // interface

		matrix()
			: m_data(), SIZE()
		{ }

		matrix(size_type const width, value_type const& val = T())
			: m_data(sysize(width), val), SIZE(width)
		{ }

		void swap(my_type& rhs)
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

		void resize(size_type const new_size, T const& val = T())
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

		// access operators
		reference operator() (size_type const i) { return m_data[i]; }
		const_reference operator() (size_type const i) const { return m_data[i]; }
		reference operator() (size_type const r, size_type const c) { return m_data[offset(r, c)]; }
		const_reference operator() (size_type const r, size_type const c) const { return m_data[offset(r, c)]; }
		reference at(size_type const i) { return m_data.at(i); }
		const_reference at(size_type const i) const { return m_data.at(i); }
		reference at(size_type const r, size_type const c) { return m_data.at(offset(r, c)); }
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

		matrix& operator*= (T const& val)
		{
			// multiply ever element with val
			for (auto& e : *this)
			{
				e *= val;
			}
			return *this;
		}


	};



	/*

		Equality checks

	*/


	// genral equality comparison
	// can compare symmetric and with non symmetric
	template<class T, bool S1, bool S2, template<class...> class C>
	bool operator== (matrix<T, S1, C> const& lhs, matrix<T, S2, C> const& rhs)
	{
		typedef typename matrix<T, true>::size_type size_type;
		if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) return false;
		// iterate all 
		for (size_type i = 0; i < lhs.rows(); ++i)
		{
			for (size_type j = 0; j < lhs.cols(); ++j)
			{
				if (lhs(i, j) != rhs(i, j)) return false;
			}
		}
		return true;
	}

	// specialization for symmetric only matrices.
	// can avoide half of the comparisons
	template<class T, template<class...> class C>
	bool operator== (matrix<T, true, C> const& lhs, matrix<T, true, C> const& rhs)
	{
		if (lhs.size() != rhs.size()) return false;
		return std::equal(lhs.begin(), lhs.end(), rhs.begin());
	}

	// non equality operator
	template<class T, bool S1, bool S2>
	bool operator!= (matrix<T, S1> const& lhs, matrix<T, S2> const& rhs)
	{
		return !(lhs == rhs);
	}

	// scalar multiplication
	template<class T, bool B, template<class...> class C>
	matrix<T, B, C> operator* (matrix<T, B, C> m, T const& v)
	{
		m *= v;
		return m;
	}

	// sym matrix is square by definition
	template<class T, template<class...> class C>
	bool is_square(matrix<T, true, C> const&)
	{
		return true;
	}
	// sym matrix is symmetric by definition
	template<class T, template<class...> class C>
	bool is_symmetric(matrix<T, true, C> const&)
	{
		return true;
	}

	// sym matrix is symmetric by definition
	// transpose does nothing
	template<class T, template<class...> class C>
	void transpose(matrix<T, true, C> const&)
	{ }



	/*

		STREAM OUTPUT

	*/




	template<class Matrix>
	struct matrix_printer
	{
	private:
		matrix_printer& operator= (matrix_printer const&);
	public:
		Matrix const& matrix;
		std::streamsize margin_left, margin_right;
		std::size_t margin_top, margin_bottom;
		char fill_left, fill_right, spacer;
		bool print_indices;
		matrix_printer(Matrix const& m) :
			matrix(m),
			margin_left(0), margin_right(0), margin_top(0), margin_bottom(0),
			fill_left(' '), fill_right(' '), spacer(' '),
			print_indices(false)
		{}
		matrix_printer operator() (Matrix const& m) const
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
		void stream_fill_left(std::ostream& strm) const
		{
			if (fill_left != '\0')
				fill(strm, fill_left, margin_left);
		}
		void stream_fill_right(std::ostream& strm) const
		{
			if (fill_right != '\0')
				fill(strm, fill_right, margin_left);
		}
	private:
		void fill(std::ostream& strm, char const fillc, std::streamsize const& width) const
		{
			char const x = strm.fill();
			strm.fill(fillc);
			strm << std::setw(width) << fill_left;
			strm.fill(x);
		}
	};

	template<class M>
	inline std::ostream& operator<< (std::ostream& strm, matrix_printer<M> const& printer)
	{
		// save width for every element
		std::streamsize const W = strm.width();
		// Margin top
		for (std::size_t i = 0; i < printer.margin_top; ++i) strm << '\n';
		// Iterate matrix
		typedef typename M::size_type size_type;
		//typedef typename M::row_type::const_iterator row_iterator;

		for (size_type i = 0; i < printer.matrix.rows(); ++i)
		{
			auto const end = printer.matrix.row(i).cend();
			printer.stream_fill_left(strm);
			for (auto it = printer.matrix.row(i).cbegin(); it != end; ++it)
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
	inline std::ostream& operator<< (std::ostream& strm, scon::matrix<T, B> const& matrix)
	{
		typedef matrix_printer< scon::matrix<T, B> > p_type;
		strm << p_type(matrix);
		return strm;
	}


	template<class T>
	struct LU_decomp
	{

		scon::matrix<T> lu;
		std::vector<std::size_t> P;
		typename std::enable_if<
			std::is_floating_point<T>::value, T>::type d;

		// construction

		LU_decomp(scon::matrix<T> const& A)
			: lu(A), P(A.rows()), d(1)
		{
			// check proper dimensionality
			if (A.rows() != A.cols())
				throw std::logic_error("LU decompostion requires square matrix.");
			// setup
			auto const tiny = T{ 1.e-40 };
			auto const n = A.rows();
			auto vv = std::vector<T>(n);
			for (std::size_t i = 0; i < n; ++i)
			{
				auto b = T{};
				for (std::size_t j = 0; j < n; ++j)
				{
					b = std::max(b, std::abs(lu(i, j)));
				}
				if (b == T{})
					throw std::runtime_error("Singular matrix in LU decompostion.");
				vv[i] = T{ 1 } / b;
			}

			for (std::size_t k = 0; k < n; ++k)
			{
				auto b = T{ 0 };
				auto imax = k;
				for (std::size_t i = k; i < n; ++i)
				{
					auto tmp = vv[i] * std::abs(lu(i, k));
					if (tmp > b)
					{
						b = tmp;
						imax = i;
					}
				}
				if (k != imax)
				{
					for (std::size_t j = 0; j < n; ++j)
					{
						std::swap(lu(imax, j), lu(k, j));
					}
					d = -d;
					vv[imax] = vv[k];
				}
				P[k] = imax;
				if (lu(k, k) == T{}) lu[k][k] = tiny;

				for (std::size_t i = k + 1; i < n; ++i)
				{
					auto tmp = lu(i, k) /= lu(k, k);
					for (std::size_t j = k + 1; j < n; ++j)
					{
						lu(i, j) -= tmp * lu(k, j);
					}
				}

			}

		}

		// return L and U

		scon::matrix<T> L() const
		{
			scon::matrix<T> l(lu);
			auto const n = l.rows();
			for (std::size_t i = 0; i < n; ++i)
			{
				for (std::size_t j = i; j < n; ++j)
				{
					l(i, j) = (i == j) ? T{ 1 } : T{ 0 };
				}
			}
			return l;
		}
		scon::matrix<T> U() const
		{
			scon::matrix<T> u(lu);
			auto const n = u.rows();
			for (std::size_t i = 1; i < n; ++i)
			{
				for (std::size_t j = 0; j < i; ++j)
				{
					u(i, j) = T{ 0 };
				}
			}
			return u;
		}

		// determinant of creation matrix A
		T determinant() const
		{
			auto r = d;
			auto const n = lu.rows();
			for (std::size_t i = 0; i < n; ++i) r *= lu(i, i);
			return r;
		}

		// solve linear equation set
		// A * x = b (param = b, return = x)
		std::vector<T> solve(std::vector<T> const& b) const
		{
			auto n = lu.rows();
			if (n != b.size())
				throw std::out_of_range("LU decompostion solve() got wrong vector size.");
			auto x = b;
			auto ii = std::size_t{};
			for (std::size_t i = 0; i < n; ++i)
			{
				auto ip = P[i];
				auto sum = x[ip];
				x[ip] = x[i];
				if (ii != std::size_t{ 0 })
				{
					for (std::size_t j = ii - 1; j < i; ++j)
					{
						sum -= lu(i, j) * x[j];
					}
				}
				else if (sum != T{}) ii = i + 1u;
				x[i] = sum;
			}
			auto const i_init = static_cast<std::ptrdiff_t>(n - 1u);
			for (auto i = i_init; i >= 0; --i)
			{
				auto sum = x[i];
				auto const j_init = static_cast<std::size_t>(i + 1);
				for (auto j = j_init; j < n; ++j)
				{
					sum -= lu(i, j) * x[j];
				}
				x[i] = sum / lu(i, i);
			}
			return x;
		}

		// solve multiple linear equation sets
		// A * X = B
		scon::matrix<T> solve(scon::matrix<T> const& b) const
		{
			auto const n = lu.rows();
			auto const m = b.cols();
			if (n != b.rows())
				throw std::out_of_range("LU decompostion solve() got wrong matrix size.");
			// result matrix
			auto x = scon::matrix<T>(n, m);
			auto xx = std::vector<T>(n, T{});
			for (std::size_t j = 0; j < m; ++j)
			{
				for (std::size_t i = 0; i < n; ++i) xx[i] = b(i, j);
				xx = solve(xx);
				for (std::size_t i = 0; i < n; ++i) x(i, j) = xx[i];
			}
			return x;
		}


		// get inverse of creation matrix
		scon::matrix<T> inverse() const
		{
			auto const n = lu.rows();
			auto inv_mat = scon::matrix<T>(n, n, T{ 0 });
			for (std::size_t i = 0; i < n; ++i)
			{
				inv_mat(i, i) = T{ 1 };
			}
			return solve(inv_mat);
		}

	};

	template<class T>
	LU_decomp<T> LU_decompose(scon::matrix<T> const& m)
	{
		return{ m };
	}

	template<class T>
	T determinant(scon::matrix<T> const& A)
	{
		return LU_decompose(A).determinant();
	}

	template<class T>
	T determinant(LU_decomp<T> const& LU)
	{
		return LU.determinant();
	}

	template<class T>
	scon::matrix<T> inverse(scon::matrix<T> const& A)
	{
		return LU_decompose(A).inverse();
	}

	template<class T>
	scon::matrix<T> inverse(LU_decomp<T> const& LU)
	{
		return LU.inverse();
	}

}


#endif
