#pragma once






/**
CAST 3
scon_mathmatrix.h
Purpose: Enabling matrix calculations. Uses Armadillo for enhanced speed when available. Otherwise uses slow internal routines.

Works on either internal scon::matrix types or arma::Mat types if
the flag USE_ARMADILLO is specified

@author Dustin Kaiser
@version 3.0
*/

/*
USAGE CONVENTIONS AS FOLLOWS:

mathmatrix(xyz, atom_nr) for mathmatrix of one frame in cartesian coordiantes
mathmatrix(dist/angle/dihedral, atom_nr) for mathmatrix of one frame in internal coordiantes
mathmatrix(coords, frames)

*/

/*
CODING CONVENTIONS AS FOLLOWS:

- Wraps underlying abstract matrix-obj of the "scon_matrix.h" - kind or arma::Mat kind
- There might be still some ambigous stuff going on regarding types.

*/



/////////////////////////////////
//                             //
//	D E F S                    //
//                             //
/////////////////////////////////

#include "coords.h"
typedef coords::float_type float_type;
typedef int int_type;
typedef size_t uint_type;
unsigned int constexpr printFunctionCallVerbosity = 5u;

///////////////////////////////
//                           //
//	I N C L U D E S          //
//                           //
///////////////////////////////

#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <cmath>
#include <limits>
#include <utility>

#ifdef CAST_USE_ARMADILLO
#include <armadillo>
#define CAST_ARMA_MATRIX_TYPE arma::Mat<T>
#else
#include "Eigen"
#define CAST_EIGEN_MATRIX_TYPE Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>
#endif


///////////////////////////////
//                           //
//	F L A G S                //
//                           //
///////////////////////////////

// Further flag important here:
// #define CAST_USE_ARMADILLO
// However, don't set this manually. This
// flag is set by Visual Studio or make according to your desired configuration
// It's all already automatized and integrated

///////////////////////////////


/////////////////////////////////
//                             //
//	D E F S                    //
//                             //
/////////////////////////////////

#include "coords.h"
typedef coords::float_type float_type;
typedef int int_type;
typedef size_t uint_type;


/**
 * @brief Class for handling matrix operations involving numerical entries.  Uses Armadillo for enhanced speed when available. Otherwise uses slow internal routines.
 * @author Dustin Kaiser
 * @version 3.0
 *
 * The class scon::mathmatrix is made to be mainly used for matrix operations
 * on floating point numbers such as multiplication, eigenvalue-decomposition and so on.
 * 
 * It is written as a generic template but it is entirely untested for types other than floating point numbers.
 *
 * In the file scon_mathmatrix_test.cc there are comprehensive unit tests for the methods of this class
 *
 * @note This class is two-faced. It is a wrapper around either 
 * - a Eigen matrix, a fast header only matrix math library for cpp
 * - or a armadillo::Mat. Armadillo is a C++ Framework which itself is a wrapper for LAPACK and BLAS, high speed fortran matrix routines.
 * If the preprocessor define "CAST_USE_ARMADILLO" is set, scon::mathmatrix will accompany an arma::Mat object. This of course
 * implies that CAST then hast to be compiled and linked with pre-existing LAPACK and BLAS libraries. Currently,
 * precompiled versions can be found in the CAST git repository in the folder optional_files. If you use the recomended
 * premake5 build automation to build CAST, everything should be fine and you do not need to worry about linking CAST with LAPACK.
 *
 * For environments where linking with BLAS and LAPACK is not possible, there is also the option of using scon::mathmatrix without armadillo.
 * For this end we have provided a stand-alone, internal version of each matrix-transformation. These are assured to yield similar results.
 * However, especially SVD-Decompositions on large matrices will be painsteakingly slow. It will be too slow for productive use by
 * computational chemists in case of the tasks PCA and ENTROPY.
 *
 * @warning Matrix access does NOT throw when out of bounds! Take caution here, only armadillo-enabled matrices on debug builds will
 * cause an exeption when a matrix is accessed out of bounds (for example, accesing a 2x2 matrix with matrix(4,4)). In ALL other scenarios,
 * CAST will continue without ANY error and you find yourself in undefined bahviour land.
 *
 * Transpose() and transposed() are not members of the mathmatrix class but available as free functions (found at the bottom of the scon_mathmatrix.h file).
 *
 * Usage conventions are as follows (these are guidelines and not enforced by design or assertions):
 * 
 * - mathmatrix(xyz, atom_nr) for mathmatrix of one frame in cartesian coordiantes.
 * - mathmatrix(dist/angle/dihedral, atom_nr) for mathmatrix of one frame in internal coordiantes.
 * - mathmatrix(coords, frames) for a matrix of a whole (MD) trajectory.
 *
 */
 
	template <typename T>
  class mathmatrix
#ifndef CAST_USE_ARMADILLO
    : public CAST_EIGEN_MATRIX_TYPE
#else
    : public CAST_ARMA_MATRIX_TYPE
#endif

	{
  private:
#ifdef CAST_USE_ARMADILLO
    using base_type = CAST_ARMA_MATRIX_TYPE;
#else
    using base_type = CAST_EIGEN_MATRIX_TYPE;
#endif


	public:

		/////////////////////////////////////
		/////                           /////
		/////  C O N S T R U C T O R S  /////
		/////                           /////
		/////////////////////////////////////

    /*! Construct empty mathmatrix
     *
     * Constructs an empty mathmatrix
    **/
#ifndef CAST_USE_ARMADILLO
    mathmatrix() : CAST_EIGEN_MATRIX_TYPE()
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Constructing empty matrix." << std::endl;
    };
#else
    mathmatrix() : CAST_ARMA_MATRIX_TYPE()
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Constructing empty matrix." << std::endl;
    };
#endif


    /*! Construct mathmatrix of certain size
     *
     * Elements are undefined
     * @param rows: Number of rows
     * @param cols: Number of columns
     */
#ifndef CAST_USE_ARMADILLO
    mathmatrix(uint_type rows, uint_type cols) : CAST_EIGEN_MATRIX_TYPE(rows, cols)
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Constructing empty " << rows << " x " << cols << " matrix." << std::endl;
    };
#else
    mathmatrix(uint_type rows, uint_type cols) : arma::Mat<T>(rows, cols) 
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Constructing empty " << rows << " x " << cols << " matrix." << std::endl;
    };
#endif

    /*! Construct mathmatrix from parent matrix
     *
     * @warning Do not call this manually. This is only used in internal operations.
     * @NOTE Used only during internal functions.
     */

#ifndef CAST_USE_ARMADILLO
    mathmatrix(CAST_EIGEN_MATRIX_TYPE const& in) : CAST_EIGEN_MATRIX_TYPE(in)
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Constructing matrix from eigen-matrix." << std::endl;
    };
    mathmatrix& operator=(CAST_EIGEN_MATRIX_TYPE const& in)
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Assignment from eigen-matrix." << std::endl;
      if (this != &in)
      {
        *this = mathmatrix(in);
      }
      return *this;
    };
#else
    mathmatrix(CAST_ARMA_MATRIX_TYPE in) : CAST_ARMA_MATRIX_TYPE(in)
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Constructing matrix from arma-matrix." << std::endl;
    };
#endif

    /*! Construct filled mathmatrix of certain size
     *
     * All values initialized to the same value
     * @param rows: Number of rows
     * @param cols: Number of columns
     * @param fill: Value to which all matrix elements will be initialized
     */
    mathmatrix(uint_type rows, uint_type cols, T fill)
      :mathmatrix(rows, cols)
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Constructing filled matrix." << std::endl;
      for (uint_type i = 0u; i < rows; i++)
        for (uint_type j = 0u; j < cols; j++)
          (*this)(i, j) = fill;
    };


    // base row and col proxy
    using base_type::row;
    using base_type::col;


    // element access from base class
    using base_type::operator();
    using base_type::operator[];
#ifdef CAST_USE_ARMADILLO
    using base_type::resize;
#else
    void resize(uint_type const rows, uint_type const cols)
    {
      this->conservativeResize(rows, cols);
    }
#endif


    using base_type::operator*=;
    using base_type::operator-=;
    using base_type::operator/=;
    using base_type::operator+=;
#ifndef CAST_USE_ARMADILLO
    using base_type::operator+;
    using base_type::operator-;
    using base_type::operator/;
    using base_type::operator*;
#else
    /*! mathmatrix += operator
     * 
     * @param in: Matrix to the right of the summation (this + in)
     * @return: Result of the addition of this + in
     */
    mathmatrix operator+(mathmatrix const& in) const
    {
        if (Config::get().general.verbosity >= printFunctionCallVerbosity)
          std::cout << "Function call: Operator+ for matrix-class" << std::endl;
    	if (!(this->rows() == in.rows() && this->cols() == in.cols() ))
    	{
    		throw("ERROR in mathmatrix Addition: Sizes of matrices do not match!");
    	}
        CAST_ARMA_MATRIX_TYPE const& base_this = *this;
        CAST_ARMA_MATRIX_TYPE const& base_in = in;
        return (mathmatrix(   base_this + base_in    ));
    };

    mathmatrix operator*(mathmatrix const& in) const
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Operator* for matrix-class" << std::endl;
      if (!(this->cols() == in.rows()))
      {
        throw("ERROR in mathmatrix multiplication: Sizes of matrices do not match!");
      }
      CAST_ARMA_MATRIX_TYPE const& base_this = *this;
      CAST_ARMA_MATRIX_TYPE const& base_in = in;
      CAST_ARMA_MATRIX_TYPE multiplied = base_this * base_in;
      return (mathmatrix(multiplied));
    }

    mathmatrix operator*(T const& in) const
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Operator* scalar for matrix-class" << std::endl;
      //Check if sizes match
      mathmatrix output = *this;
      for (uint_type i = 0; i < this->rows(); i++)
      {
        for (uint_type j = 0; j < this->cols(); j++)
        {
          output(i, j) *= in;
        }
      }
      return output;
    }
    mathmatrix operator/(mathmatrix const& in) const
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Operator/ for matrix-class" << std::endl;
      if (!(this->cols() == in.rows()))
      {
        throw("ERROR in mathmatrix divison: Sizes of matrices do not match!");
      }
      CAST_ARMA_MATRIX_TYPE const& base_this = *this;
      CAST_ARMA_MATRIX_TYPE const& base_in = in;
      return (mathmatrix(base_this / base_in));
    };

    /**
     * @brief Overload "-" Operator for mathmatrix
     */
    mathmatrix operator-(mathmatrix const& in) const
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Operator- for matrix-class" << std::endl;
    	//Check if sizes match
    	if ((in.rows() != this->rows()) || (in.cols() != this->cols()))
    	{
    		throw std::runtime_error("Error in Matrix mathmatrix subtraction, wrong input sizes");
    	}
    	mathmatrix output = *this;
    	for (uint_type i = 0; i < this->rows(); i++)
    	{
    		for (uint_type j = 0; j < this->cols(); j++)
    		{
    			output(i, j) -= in(i, j);
    		}
    	}
    	return output;
    }
#endif


    /*! Returns an "identity matrix" of certain size
     *
     * All diagonal elements will be initialized to one
     * All off-diagonal elements will be initialized to zero
     * @param num_rows: Number of rows
     * @param num_cols: Number of columns
     * @todo: Write as free function!!
     */
    static typename std::enable_if<std::is_arithmetic<T>::value, mathmatrix>::type
      Identity(std::size_t const num_rows, std::size_t const num_cols)
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Constructing identity matrix." << std::endl;
#ifndef CAST_USE_ARMADILLO
      return mathmatrix<T>(base_type::Identity(num_rows, num_cols));
#else
      return mathmatrix<T>(mathmatrix(num_rows, num_cols).eye());
#endif
    }


    // in case you are wondering:
    // transposed and some more stuff is available as free functions
    // eg transpose().

		/////////////////////////////////////
		/////                           /////
		/////  O P E R A T I O N S      /////
		/////                           /////
		/////////////////////////////////////


    /*! mathmatrix stream operator for writing output
     *
     * @param os: Stream to which data is passed on
     * @param object: mathmatrix whose contents are to be printed
     *
     * @return: Formatted contents of the matrix as strings (without a header)
     */
    friend std::ostream& operator<<(std::ostream& os, mathmatrix const & object)
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Operator>> for matrix-class" << std::endl;
      for (size_t i = 0u; i < object.rows(); ++i)
      {
        for (size_t j = 0u; j < object.cols(); ++j)
        {
          os << std::setw(18) << std::scientific << std::left << std::setprecision(8) << object(i, j) << " ";
        }
        os << "\n";
      }
      return os;
    };

		/**
		 * Checks for positive_definite matrix
		 */
		bool positive_definite_check()
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: positive_definite_check() for mathmatrix." << std::endl;
			bool positive_definite = true;
      if (!this->return_quadratic()) return false;
      for (unsigned int i = 1; i < (this->rows() + 1); i++)
      {
        if (mathmatrix( this->upper_left_submatrix(i, i) ).determ() <= 0)
        {
          return false;
        }
      }
			return true;
		};

		/**
		 * @brief Returns sign of the determinant (-1 / 1).
     * @return -1 if determinant is negative or zero, +1 if determinant is greater than zero.
		 */
    int_type det_sign() const
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: det_sign() for mathmatrix." << std::endl;
			return ((this->determ() <= 0) ? -1 : 1);
		};

		/**
		 * @brief Returns rank of the underlying matrix.
     * @return rank of the matrix.
     * @note A singular value decompostition is performed to determine the rank. This may be quite costly.
		 */
		size_t rank() const
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: rank() for mathmatrix." << std::endl;
#ifdef CAST_USE_ARMADILLO
      CAST_ARMA_MATRIX_TYPE const& base_this = *this;
      return static_cast<size_t>(arma::rank(base_this));
#else
      Eigen::ColPivHouseholderQR<CAST_EIGEN_MATRIX_TYPE> rank_colpivmat (*this);
      return static_cast<size_t>(rank_colpivmat.rank());
#endif
		};

		/**
		 * @brief Calculates determinant of the mathmatrix-obj
     *
     * Internal code uses a LU decompostion.
     * @return determinant
		 */
		float_type determ() const
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Starting determinant calculation of " << this->rows() << " x " << this->cols() << " matrix." << std::endl;

#ifdef CAST_USE_ARMADILLO
			return static_cast<float_type>(det(*this));
#else
      return static_cast<float_type>(static_cast<CAST_EIGEN_MATRIX_TYPE const*>(this)->determinant());
#endif
		}


		/**
		 * @brief Overload "/" Operator for mathmatrix and scalars
		 */
		mathmatrix operator/(T const& in) const
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: operator/ for mathmatrix." << std::endl;
			mathmatrix tempCopy(*this);
			for (uint_type i = 0; i < this->rows(); i++)
			{
				for (uint_type j = 0; j < this->cols(); j++)
				{
					tempCopy(i, j) = (*this)(i, j) / in;
				}
			}
			return tempCopy;
		}

		/**
		 * @brief Append one matrix to another, will check if sizes match, appends at the bottom end (rows are added)
		 */
		void append_bottom(const mathmatrix& I_will_be_the_bottom_part)
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: append_bottom for mathmatrix." << std::endl;
			//if this is transposed, append bottom means append right on underlying obj
			if (this->cols() != I_will_be_the_bottom_part.cols())
			{
				throw "Wrong Matrix size in mathmatrix:append()";
			}
			//Old size needs to be kept
			size_t holder = this->rows();

			this->resize(this->rows() + I_will_be_the_bottom_part.rows(), this->cols());

			//Add "in" to newly created space.
			for (uint_type i = 0; i < I_will_be_the_bottom_part.rows(); i++)
			{
				for (uint_type j = 0; j < this->cols(); j++)
				{
					(*this)(i + holder, j) = I_will_be_the_bottom_part(i, j);
				}
			}
		}

		/**
		 * @brief Append one matrix to another, will check if sizes match, appends on the top (rows are added)
		 */
		void append_top(const mathmatrix& I_will_be_the_top_part)
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: append_top for mathmatrix." << std::endl;
			if (this->cols() != I_will_be_the_top_part.cols())
			{
				throw "Wrong Matrix size in mathmatrix:append()";
			}
      const size_t thisOldRows = this->rows();
			this->resize(this->rows() + I_will_be_the_top_part.rows(), this->cols());

			//Move the entries in the parent matrix downward
			//We count downwards so that we dont overwrite
			for (uint_type i = thisOldRows - 1u; i > 0; i--)
			{
				for (uint_type j = 0; j < this->cols(); j++)
				{
					(*this)(i + I_will_be_the_top_part.rows(), j) = (*this)(i, j);
				}
			}

			//Add "I_will_be_the_top_part" to now absolete top space of the parent matrix (this).
			for (uint_type i = 0; i < I_will_be_the_top_part.rows(); i++)
			{
				for (uint_type j = 0; j < this->cols(); j++)
				{
					(*this)(i, j) = I_will_be_the_top_part(i, j);
				}
			}
		}

		/**
		 * @brief Append one matrix to another, will check if sizes match, appends left (columns are added)
		 */
		void append_left(const mathmatrix& I_will_be_the_left_part)
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: append_left for mathmatrix." << std::endl;
			if (this->rows() != I_will_be_the_left_part.rows())
			{
				throw "Wrong Matrix size in mathmatrix:append()";
			}
      const size_t thisOldCols = this->cols();

			this->resize(this->rows(), this->cols() + I_will_be_the_left_part.cols());

			//Move the entries in the parent matrix rightward
			//We count right so that we dont overwrite
			for (uint_type j = thisOldCols - 1u; j > 0; j--)
			{
				for (uint_type i = 0u; i < this->rows(); i++)
				{
					(*this)(i, j + I_will_be_the_left_part.cols()) = (*this)(i, j);
				}
			}

			//Add "I_will_be_the_left_part" to now absolete top space of the parent matrix (this).
			for (uint_type j = 0; j < I_will_be_the_left_part.cols(); j++)
			{
				for (uint_type i = 0; i < this->rows(); i++)
				{
					(*this)(i, j) = I_will_be_the_left_part(i, j);
				}
			}
		}

		/**
	   * @brief Append one matrix to another, will check if sizes match, appends left (columns are added)
  	 */
		void append_right(const mathmatrix& I_will_be_the_right_part)
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: append_right for mathmatrix." << std::endl;
			if (this->rows() != I_will_be_the_right_part.rows())
			{
				throw "Wrong Matrix size in mathmatrix:append()";
			}

			//Old size needs to be kept
      uint_type holder = this->cols();

			this->resize(this->rows(), this->cols() + I_will_be_the_right_part.cols());

			//Add "in" to newly created space.
			for (uint_type j = 0; j < I_will_be_the_right_part.cols(); j++)
			{
				for (uint_type i = 0; i < this->rows(); i++)
				{
					(*this)(i, j + holder - 1) = I_will_be_the_right_part(i, j);
				}
			}

		}

		/**
		 * @brief Sheds the specified rows from the matrix
		 */
		void shed_rows(size_t first_in, size_t last_in)
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: shed_rows for mathmatrix." << std::endl;
      if (first_in > last_in || last_in >= this->rows())
      {
        throw std::runtime_error("Index Out of Bounds in mathmatrix:shed_rows()");
      }
			mathmatrix newOne(this->rows() - (last_in - first_in + 1u), this->cols());
			for (size_t i = 0u; i < first_in; i++)
			{
				for (size_t j = 0u; j < this->cols(); j++)
				{
					newOne(i, j) = (*this)(i, j);
				}
			}
			for (size_t i = first_in; i < this->rows() - last_in - 1 + first_in; i++)
			{
				for (size_t j = 0u; j < this->cols(); j++)
				{
					newOne(i, j) = (*this)(i + (last_in - first_in) + 1, j);
				}
			}
			this->swap(newOne);
		}

		/**
		 * @brief Sheds the specified columns from the matrix
		 */
		void shed_cols(size_t first_in, size_t last_in)
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: shed_cols for mathmatrix." << std::endl;
      if (last_in >= this->cols() || first_in > last_in)
      {
        throw std::runtime_error("Index Out of Bounds in mathmatrix:shed_rows()");
      }
			mathmatrix newOne(this->rows(), this->cols() - (last_in - first_in + 1u));
			for (size_t j = 0u; j < this->rows(); j++)
			{
				for (size_t i = 0u; i < first_in; i++)
				{
					newOne(j, i) = (*this)(j, i);
				}
			}
			for (size_t j = 0u; j < this->rows(); j++)
			{
				for (size_t i = first_in; i < this->cols() - last_in - 1 + first_in; i++)
				{
					newOne(j, i) = (*this)(j, i + (last_in - first_in) + 1);
				}
			}
			this->swap(newOne);
		}

		/**
		 * @brief Returns number of rows
		 */
#ifndef CAST_USE_ARMADILLO
    using base_type::rows;
#else
    size_t rows() const
    {
      return this->n_rows;
    }
#endif

		/**
		 * @brief Returns number of columns
		 */
#ifndef CAST_USE_ARMADILLO
    using base_type::cols;
#else
    size_t cols() const
    {
      return this->n_cols;
    }
#endif

    /*! Performs Cholesky Decompostion on Matrix.
     * 
     * @NOTE: Code via https://rosettacode.org/wiki/Cholesky_decomposition#C 19.11.16
     * @param result: Upper triangular matrix as result of decompostion
     */
    void choleskyDecomposition(mathmatrix<T> & result) const
    {
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: choleskyDecomposition for mathmatrix." << std::endl;

      result = mathmatrix<T>(this->rows(), this->cols(), T(0));
      int n = static_cast<int>(this->rows());
      for (int i = 0; i < n; i++)
        for (int j = 0; j < (i + 1); j++) {
          T s = 0;
          for (int k = 0; k < j; k++)
            s += result(i, k) * result(j, k);
          result(i, j) = (i == j) ?
            ::sqrt((*this)(i, i) - s) :
            (1.0 / result(j, j) * ((*this)(i, j) - s));
        }
#if defined (_MSC_VER) && ! defined(CAST_USE_ARMADILLO)
      ::transpose(result);
#else
      transpose(result);
#endif
    }

		/**
		 * @brief Returns whether mathmatrix-obj is quadratic
		 */
		bool return_quadratic() const
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: return_quadratic for mathmatrix." << std::endl;
      return this->rows() == this->cols();
		}

		/**
		 * @brief Returns upper left submatrix. 
     * If no second argument for the function call,
		 * ie for columns_in is specified, then a quadratic submatrix with rows = columns = rows_in
		 * is yieled.
		 */
		mathmatrix upper_left_submatrix(uint_type rows_in, uint_type columns_in = 0) const
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: upper_left_submatrix for mathmatrix." << std::endl;
      if (columns_in == 0)
        columns_in = rows_in;
      if (rows_in <= this->rows() && columns_in <= this->cols())
      {
        mathmatrix copied(*this);
        copied.resize(rows_in, columns_in);
        return copied;
      }
      else
      {
        return *this;
      }
		}

    /*! Equality operator for armadillo mathmatrix
     *
     * Armadillo internally handels the operator== in a 
     * from my perspective very strange way. This is why we
     * use the approx_equal function with a very tight tolerance
     * of 0.1% internally to check for equality.
     *
     * @param in: Matrix that *this is compared to
     * @return: boolean that indicates if size and all elements are equal
     */
#ifndef CAST_TOLERANCE_FOR_MATRIX_COMPARISSON
#define CAST_TOLERANCE_FOR_MATRIX_COMPARISSON 0.001
#endif
    bool operator== (mathmatrix const& in) const
    {
#ifdef CAST_USE_ARMADILLO
      CAST_ARMA_MATRIX_TYPE const& a = *this;
      CAST_ARMA_MATRIX_TYPE const& b = in;
      return (approx_equal(a,b,"reldiff", CAST_TOLERANCE_FOR_MATRIX_COMPARISSON) );
#else
      if (this->rows() != in.rows() || this->cols() != in.cols()) return false;
      else
      {
        for (size_t i = 0u; i < this->rows(); i++)
        {
          for (size_t j = 0u; j < this->cols(); j++)
          {
            if (abs((*this)(i, j) - in(i, j)) > CAST_TOLERANCE_FOR_MATRIX_COMPARISSON * abs(in(i, j)))
            {
              return false;
            }
          }
        }
        return true;
      }
#endif
    }


		/**
		 * @brief Performs singular value decomposition on *this and writes results
		 * to the three resulting matrices U, s, V. 
     * 
     * Since the rank of a matrix
		 * is automatically calculated during our SVD-decomposition, a variable "rank"
		 * may be specified where the calculated rank may be stored for future use.
     *
     * justification: ptr is used for rank because it can represent a
     * bool intrinsically through nullptr and therefore stores information
     * wether rank should be passed through this function.
		 */
		void singular_value_decomposition(mathmatrix& U_in, mathmatrix& s_in, mathmatrix& V_in, int* rank = nullptr) const
		{
      U_in = *this;
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: Starting singular value decomposition of " << U_in.rows() << " x " << U_in.cols() << " matrix." << std::endl;

			//"Constructor"
			V_in.resize(this->cols(), this->cols());
			s_in.resize(this->cols(), 1u);

#ifndef CAST_USE_ARMADILLO
      CAST_EIGEN_MATRIX_TYPE const* this_eigenptr= static_cast<CAST_EIGEN_MATRIX_TYPE const*>(this);
      auto bdcsvd_obj = this_eigenptr->bdcSvd();
      bdcsvd_obj.compute(*this_eigenptr, Eigen::ComputeFullU | Eigen::ComputeFullV);
      U_in = mathmatrix<T>(bdcsvd_obj.matrixU());
      V_in = mathmatrix<T>(bdcsvd_obj.matrixV());
      for (size_t i = 0; i < U_in.rows(); i++)
      {
        s_in(i) = bdcsvd_obj.singularValues()(i);
      }
#else
			arma::Col<float_type> s;
			if(!svd_econ(U_in, s, V_in, *this)) throw std::runtime_error("Error in armadillo SVD: failed.");
			for (size_t i = 0; i < U_in.rows(); i++)
			{
			  s_in(i) = s(i);
			}
      if (rank != nullptr) *rank = arma::rank(*this);
#endif
		}

		/**
		 * @brief Performs eigenvalue decomposition on *this and returns eigenval and eigenvec
		 * as well as the rank of *this matrix.
		 *
		 * @param eigenval_in Eigenvalues will be written to this matrix
		 * @param eigenvec_in Eigenvectors will be written to this matrix
		 * @param rank_in Rank of the matrix is stored here
		 *
		 * justification: ptr is used for rank because in can represent a
		 * bool intrinsically through nullptr and therefore stores information
		 * weather rank should be passed through this function.
		 *
		 * @note Matrix is assumed to be symmetric.
		 */
		void eigensym(mathmatrix& eigenval_in, mathmatrix& eigenvec_in, int* rank_in = nullptr) const
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: eigensym for mathmatrix." << std::endl;
			//On basis of SVD, so take care that your matrix is symmetrical

			//And if you are unsure about your symmetry, try this beforehand:
			//this->symmetry_check();
			//(Although it might slow stuff down considerably)
#ifndef CAST_USE_ARMADILLO
			mathmatrix V;
			this->singular_value_decomposition(eigenvec_in, eigenval_in, V, rank_in);
#else
      eigenvec_in.resize(this->cols(), this->cols());
      eigenval_in.resize(this->cols(), 1u);
      arma::Col<float_type> s;
      arma::Mat<float_type> eigVecUnordered(eigenvec_in.rows(), eigenvec_in.cols());
      eig_sym(s, eigVecUnordered, *this);

      for (size_t i = 0; i < (*this).rows(); i++)
      {
        for (size_t j = 0; j < (*this).cols(); j++)
        {
          eigenvec_in(i, (*this).cols() - j - 1u) = eigVecUnordered(i, j);
        }
        eigenval_in((*this).rows() - i - 1u) = s(i);
      }
      if (rank_in != nullptr) *rank_in = arma::rank(*this);
#endif
		}

    mathmatrix inversed() const
    {
#ifndef CAST_USE_ARMADILLO
      CAST_EIGEN_MATRIX_TYPE const* this_eigenptr = static_cast<CAST_EIGEN_MATRIX_TYPE const*>(this);
      auto inversed = this_eigenptr->inverse();
      CAST_EIGEN_MATRIX_TYPE to_be_returned(inversed);
      return to_be_returned;
#else
      CAST_ARMA_MATRIX_TYPE const& base_this = *this;
      return mathmatrix(inv(base_this));
#endif
    }


		/**
		 * @breif Returns mathmatrix-obj as std vector of vector of float_type.
		 * USE THIS TO DEBUG, not in production code.
		 * This might be useful as the VS debugger cannot visualize content of arma arrays.
     *
     * @return Matrix in std::vector<vector<float-type> > form.
		 */
		std::vector <std::vector<T> > to_std_vector(void) const
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: to_std_vector for mathmatrix" << std::endl;
			std::vector <T> temp1(this->cols());
			std::vector < std::vector <T> > temp2(this->rows(), temp1);
			for (uint_type i = 0; i < this->rows(); i++)
			{
				for (uint_type j = 0; j < this->cols(); j++)
				{
					temp2[i][j] = (*this)(i, j);
				}
			}
			return temp2;
		}

		/**
		 * Updates the internal std::vector<std::vector<float_type> > array_debugview_internal array
		 * ONLY IF PREPROCESSOR FLAG "DEBUGVIEW" IS SET
		 */
		void update_debugview(void)
		{
      if (Config::get().general.verbosity >= printFunctionCallVerbosity)
        std::cout << "Function call: update_debugview." << std::endl;
#ifdef DEBUGVIEW
				std::vector <float_type> temp1(this->cols());
				std::vector < std::vector <float_type> > temp2(this->rows(), temp1);
				for (unsigned int i = 0; i < this->rows(); i++)
				{
					for (unsigned int j = 0; j < this->cols(); j++)
					{
						temp2[i][j] = (*this)(i, j);
					}
				}
				array_debugview_internal = temp2;
#endif
#ifndef DEBUGVIEW
				std::cout << "DEBUGVIEW Flag not enabled. update_debugview(void) may not be used.\nCheck your code.\n";
#endif

		}
	};


  //FREE FUNCTIONS
template<typename T>
void pow(mathmatrix<T> &matrix_in, T const& exp)
{
  for (unsigned int i = 0; i < matrix_in.rows(); i++)
  {
    for (unsigned int j = 0; j < matrix_in.cols(); j++)
    {
      matrix_in(i, j) = std::pow(matrix_in(i, j), exp);
    }
  }
}



  template<typename T>
  mathmatrix<T> transposed(mathmatrix<T> const& in)
  {
    mathmatrix<T> toBeReturned(in);
    transpose(toBeReturned);
    return toBeReturned;
  }

#ifdef CAST_USE_ARMADILLO
  template<typename T>
  void transpose(mathmatrix<T>& in)
  {
    in = mathmatrix<T>(static_cast<CAST_ARMA_MATRIX_TYPE&>(in).t());
  }
#else

  template<typename T>
  void transpose(mathmatrix<T>& in)
  {
    static_cast<CAST_EIGEN_MATRIX_TYPE&>(in).transposeInPlace();
  }
#endif

  /**
  * @brief Rotation Class. If armadillo is enabled it uses the LAPACK matrix routines otherwise it uses the Eigen matrices.
  * @author Julian Erdmannsdoerfer
  * @version 3.0
  *
  * 
  *
  */

  class RotationMatrix 
  {
  public:
#ifdef CAST_USE_ARMADILLO
    //Too much work for now.
#else
    using Translation = Eigen::Translation3d;
    using Rotation = Eigen::AngleAxisd;
    using Transformation = Eigen::Affine3d;
    using Vector = Eigen::Vector3d;

    static Transformation rotate_around_axis_with_center(double rad_deg, Vector axis, Vector center){
      Translation back(center);
      Translation to_center(-center);
      return Transformation(back * Rotation(rad_deg, axis.normalized()) * to_center);
    }
    static Transformation rotate_around_axis_in_center(double rad_deg, Vector axis) {
      Vector center{ 0., 0., 0. };
      Translation back(center);
      Translation to_center(-center);
      return Transformation(back * Rotation(rad_deg, axis.normalized()) * to_center);
    }
#endif
  };
//END HEADER
