#if !defined (SCON_MATRIXTSQC_HEADER)
#define SCON_MATRIXTSQC_HEADER

#include "scon_traits.h"
#include "scon_iterator.h"
#include "scon_vect.h"
#include <vector>
#include <array>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>

namespace scon
{

	template<typename T, std::size_t N>
	class tsqcmatrix
		: public std::array<T, ((N* N) + N) / 2>
	{
	public:
		typedef tsqcmatrix<T, N>                                 my_type;
		typedef std::array<T, ((N* N) + N) / 2>                   vector_type;
		typedef typename vector_type::value_type                 value_type;
		typedef typename vector_type::size_type                  size_type;
		typedef typename vector_type::difference_type            difference_type;
		typedef typename vector_type::reference                  reference;
		typedef typename vector_type::const_reference            const_reference;
		typedef typename vector_type::iterator                   iterator;
		typedef typename vector_type::const_iterator             const_iterator;
		typedef typename vector_type::reverse_iterator           reverse_iterator;
		typedef typename vector_type::const_reverse_iterator     const_reverse_iterator;
		/*	typedef constant_stride_iterator<value_type, N>          strided_iterator;
			typedef constant_stride_iterator<const value_type, N>    const_strided_iterator;*/

		size_type cols(void) const { return N; }
		size_type rows(void) const { return N; }

		inline reference operator() (size_type const row, size_type const col)
		{
			return *(this->begin() + ((row > col) ? offset(row, col) : offset(col, row)));
		}

		inline const_reference operator() (size_type const row, size_type const col) const
		{
			return *(this->begin() + ((row > col) ? offset(row, col) : offset(col, row)));
		}

		inline reference operator() (size_type const i)
		{
			return (*this)[i];
		}

		inline const_reference operator() (size_type const i) const
		{
			return (*this)[i];
		}

		reference at(size_type const row, size_type const col)
		{
			if (row > N || col > N) throw std::logic_error("scon::matrix::at() - out of range");
			return *(this->begin() + ((row > col) ? offset(row, col) : offset(col, row)));
		}
		const_reference at(size_type const row, size_type const col)  const
		{
			if (row > N || col > N) throw std::logic_error("scon::matrix::at() - out of range");
			return *(this->begin() + ((row > col) ? offset(row, col) : offset(col, row)));
		}

		my_type& operator += (my_type const& operand)
		{
			for (size_type i(0U); i < N; ++i) (*this)[i] += operand[i];
			return *this;
		}
		my_type& operator -= (my_type const& operand)
		{
			for (size_type i(0U); i < N; ++i) (*this)[i] -= operand[i];
			return *this;
		}
		my_type& operator *= (my_type const& operand)
		{
			for (size_type i(0U); i < N; ++i) (*this)[i] *= operand[i];
			return *this;
		}
		my_type& operator /= (my_type const& operand)
		{
			for (size_type i(0U); i < N; ++i) (*this)[i] /= operand[i];
			return *this;
		}

		my_type& operator += (T const& operand)
		{
			for (size_type i(0U); i < N; ++i) (*this)[i] += operand;
			return *this;
		}
		my_type& operator -= (T const& operand)
		{
			for (size_type i(0U); i < N; ++i) (*this)[i] -= operand;
			return *this;
		}
		my_type& operator *= (T const& operand)
		{
			for (size_type i(0U); i < N; ++i) (*this)[i] *= operand;
			return *this;
		}
		my_type& operator /= (T const& operand)
		{
			for (size_type i(0U); i < N; ++i) (*this)[i] /= operand;
			return *this;
		}

		my_type operator + (T const& operand)
		{
			my_type tmp(*this);
			tmp += operand;
			return tmp;
		}
		my_type operator - (T const& operand)
		{
			my_type tmp(*this);
			tmp -= operand;
			return tmp;
		}
		my_type operator * (T const& operand)
		{
			my_type tmp(*this);
			tmp *= operand;
			return tmp;
		}
		my_type operator / (T const& operand)
		{
			my_type tmp(*this);
			tmp /= operand;
			return tmp;
		}

		my_type operator + (my_type const& operand)
		{
			my_type tmp(*this);
			tmp += operand;
			return tmp;
		}
		my_type operator - (my_type const& operand)
		{
			my_type tmp(*this);
			tmp -= operand;
			return tmp;
		}
		my_type operator * (my_type const& operand)
		{
			my_type tmp(*this);
			tmp *= operand;
			return tmp;
		}
		my_type operator / (my_type const& operand)
		{
			my_type tmp(*this);
			tmp /= operand;
			return tmp;
		}

	private:

		static inline size_type offset(size_type const row, size_type const col)
		{
			return row * (row + 1) / 2 + col;
		}

	};

}

#endif
