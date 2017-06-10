///////////////////////////////
//                           //
//	I N C L U D E S          //
//                           //
///////////////////////////////

#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <cmath>
#include <limits>
#include <utility>

#ifdef CAST_USE_ARMADILLO
#include <armadillo>
#else
#include "scon_matrix.h"
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
 * - a scon::matrix, a fast but primtive C++ Standard compliant container used also elsewhere in CAST
 * - or a armadillo::Mat. Armadillo is a C++ Framework which itself is a wrapper for LAPACK and BLAS, high speed fortran matrix routines.
 * If the preprocessor define "CAST_USE_ARMADILLO" is set, scon::mathmatrix will accompany an arma::Mat object. This of course
 * implies that CAST then hast to be compiled and linked with pre-existing LAPACK and BLAS libraries. Currently,
 * precompiled versions can be found in the CAST git repository in the folder optional_files. If you use the recomended
 * premake5 build automation to build CAST, everything should be fine and you do not need to worry about linking CAST with LAPACK.
 *
 * For environments where linking with BLAS and LAPACK is not possible, there is also the option of using scon::mathmatrix without armadillo.
 * For this end we have provided a stand-alone, internal evrsion of each matrix-transformation. These are assured to yield similar results.
 * However, especially SVD-Decompositions on large matrices will be painsteakingly slow. It will be too slow for productive use by
 * computational chemists in case of the tasks PCA an ENTROPY.
 *
 * @warning Matrix access does NOT throw when out of bounds!!! Take caution here, only armadillo-enabled matrices on debug builds will
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
    : public scon::matrix<T>
#else
    : public arma::Mat<T>
#endif

	{
  private:
#ifndef CAST_USE_ARMADILLO
    using base_type = scon::matrix<T>;
#else
    using base_type = arma::Mat<T>;
#endif


	public:

		/////////////////////////////////////
		/////                           /////
		/////  C O N S T R U C T O R S  /////
		/////                           /////
		/////////////////////////////////////

#ifndef CAST_USE_ARMADILLO
    template<class ... Args>
    mathmatrix(Args && ... args) : base_type(std::forward<Args>(args)...) {}
#else
    /*! Construct empty mathmatrix
     *
     * Constructs an empty mathmatrix
     */
    mathmatrix() : arma::Mat<T>() 
    {
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: Constructing empty matrix." << std::endl;
    };

    /*! Construct mathmatrix of certain size
     *
     * Elements are undefined
     * @param rows: Number of rows
     * @param cols: Number of columns
     */
    mathmatrix(uint_type rows, uint_type cols) : arma::Mat<T>(rows, cols) 
    {
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: Constructing empty " << rows << " x " << cols << " matrix." << std::endl;
    };


    /*! Construct mathmatrix from armadillo matrix
     *
     * @warning Do not call this manually. This is only used in internal operations.
     * @NOTE Used only during internal functions when using armadillo.
     */
    mathmatrix(arma::Mat<T> in) : arma::Mat<T>(in) 
    {
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: Constructing matrix from arma-matrix." << std::endl;
    };

    /*! Construct filled mathmatrix of certain size
     *
     * All values initialized to the same value
     * @param rows: Number of rows
     * @param cols: Number of columns
     * @param fill: Value to which all matrix elements will be initialized
     */
    mathmatrix(uint_type rows, uint_type cols, T fill) : arma::Mat<T>(rows, cols)
    {
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: Constructing filled matrix." << std::endl;
      for (uint_type i = 0u; i < rows; i++)
        for (uint_type j = 0u; j < cols; j++)
          (*this)(i, j) = fill;
    };
#endif

    // pull in range functions
    // if we use scon::matrix
#ifndef CAST_USE_ARMADILLO
    using base_type::begin;
    using base_type::end;
    using base_type::cbegin;
    using base_type::cend;
    using base_type::rbegin;
    using base_type::rend;
    using base_type::crbegin;
    using base_type::crend;
#endif

    // base row and col proxy
    using base_type::row;
    using base_type::col;


    // element access from base class
    using base_type::operator();
    using base_type::operator[];

    using base_type::resize;

#ifndef CAST_USE_ARMADILLO
    using base_type::operator*=;


#endif

    // identity from base_classes
    // but wrapped in arma case to yield identical function name
#ifndef CAST_USE_ARMADILLO
    using base_type::identity;
#else

    /*! Returns an "identity matrix" of certain size
     *
     * All diagonal elements will be initialized to one
     * All off-diagonal elements will be initialized to zero
     * @param num_rows: Number of rows
     * @param num_cols: Number of columns
     * @todo: Write as free function!!
     */
    static typename std::enable_if<std::is_arithmetic<T>::value, mathmatrix>::type
      identity(std::size_t const num_rows, std::size_t const num_cols)
    {
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: Constructing identity matrix." << std::endl;
      return mathmatrix<T>(mathmatrix(num_rows, num_cols).eye());
    }
#endif

    // in case you are wondering:
    // transposed and some more stuff is available as free functions
    // eg transpose().

		/////////////////////////////////////
		/////                           /////
		/////  O P E R A T I O N S      /////
		/////                           /////
		/////////////////////////////////////

#ifdef CAST_USE_ARMADILLO
		/*! mathmatrix += operator
		 * 
     * @param in: Matrix to the right of the summation (this + in)
     * @return: Result of the addition of this + in
		 */
		mathmatrix operator+(mathmatrix const& in) const
		{
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: Operator+ for matrix-class" << std::endl;
			if (!(this->rows() == in.rows() && this->cols() == in.cols() ))
			{
				throw("ERROR in mathmatrix Addition: Sizes of matrices do not match!");
			}
      arma::Mat<T> const& base_this = *this;
      arma::Mat<T> const& base_in = in;
      return (mathmatrix(   base_this + base_in    ));
		};
#endif

    /*! mathmatrix stream operator for writing output
     *
     * @param os: Stream to which data is passed on
     * @param object: mathmatrix whose contents are to be printed
     *
     * @return: Formatted contents of the matrix as strings (without a header)
     */
    friend std::ostream& operator<<(std::ostream& os, mathmatrix const & object)
    {
      if (Config::get().general.verbosity > 4U)
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
      if (Config::get().general.verbosity > 4U)
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
      if (Config::get().general.verbosity > 4U)
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
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: rank() for mathmatrix." << std::endl;
#ifdef CAST_USE_ARMADILLO
      arma::Mat<T> const& base_this = *this;
      return static_cast<size_t>(arma::rank(base_this));
#else
			//Rewritten from NumRecipies, SVD

			//"Constructor"
			mathmatrix U_in = *this;
			mathmatrix V_in(this->cols(), this->cols());
			mathmatrix s_in(this->cols(), 1u);

			//Numerical Recipies Nomenclatur, I don't even.... who would name it like that? nnm strcpcps wtf
			int n = this->cols();
			int m = this->rows();

			float_type eps = std::numeric_limits<float_type>::epsilon();

			//Beginn SVD::decompose()
			{
				bool flag;
				int i, its, j, jj, k, l, nm;
				T anorm, c, f, g, h, s, scale, x, y, z;
				mathmatrix rv1(n, 1u);
				g = scale = anorm = 0.0;
				for (i = 0; i < n; i++) {
					l = i + 2;
					rv1(i) = scale*g;
					g = s = scale = 0.0;
					if (i < m) {
						for (k = i; k < m; k++) scale += std::abs(U_in(k, i));
						if (scale != 0.0) {
							for (k = i; k < m; k++) {
								U_in(k, i) /= scale;
								s += U_in(k, i) * U_in(k, i);
							}
							f = U_in(i, i);
							//g = -SIGN(sqrt(s), f); //((b) >= 0.0 ? fabs(a) : -fabs(a))
							g = -1 * (f >= 0.0 ? std::fabs(sqrt(s)) : -std::fabs(sqrt(s)));
							h = f*g - s;
							U_in(i, i) = f - g;
							for (j = l - 1; j < n; j++) {
								for (s = 0.0, k = i; k < m; k++) s += U_in(k, i) * U_in(k, j);
								f = s / h;
								for (k = i; k < m; k++)  U_in(k, j) += f* U_in(k, i);
							}
							for (k = i; k < m; k++) U_in(k, i) *= scale;
						}
					}
					s_in(i) = scale *g;
					g = s = scale = 0.0;
					if (i + 1 <= m && i + 1 != n) {
						for (k = l - 1; k < n; k++) scale += std::abs(U_in(i, k));
						if (scale != 0.0) {
							for (k = l - 1; k < n; k++) {
								U_in(i, k) /= scale;
								s += U_in(i, k) * U_in(i, k);
							}
							f = U_in(i, l - 1);
							g = -1 * (f >= 0.0 ? std::fabs(sqrt(s)) : -std::fabs(sqrt(s))); //g = -SIGN(sqrt(s), f)
							h = f*g - s;
							U_in(i, l - 1) = f - g;
							for (k = l - 1; k < n; k++) rv1(k) = U_in(i, k) / h;
							for (j = l - 1; j < m; j++) {
								for (s = 0.0, k = l - 1; k < n; k++) s += U_in(j, k) * U_in(i, k);
								for (k = l - 1; k < n; k++) U_in(j, k) += s*rv1(k);
							}
							for (k = l - 1; k < n; k++) U_in(i, k) *= scale;
						}
					}
					anorm = std::max(anorm, (std::abs(s_in(i)) + std::abs(rv1(i))));
				}
				for (i = n - 1; i >= 0; i--) {
					if (i < n - 1) {
						if (g != 0.0) {
							for (j = l; j < n; j++)
								V_in(j, i) = (U_in(i, j) / U_in(i, l)) / g;
							for (j = l; j < n; j++) {
								for (s = 0.0, k = l; k < n; k++) s += U_in(i, k) * V_in(k, j);
								for (k = l; k < n; k++) V_in(k, j) += s*V_in(k, i);
							}
						}
						for (j = l; j < n; j++) V_in(i, j) = V_in(j, i) = 0.0;
					}
					V_in(i, i) = 1.0;
					g = rv1(i);
					l = i;
				}
				for (i = std::min(m, n) - 1; i >= 0; i--) {
					l = i + 1;
					g = s_in(i);
					for (j = l; j < n; j++) U_in(i, j) = 0.0;
					if (g != 0.0) {
						g = 1.0 / g;
						for (j = l; j < n; j++) {
							for (s = 0.0, k = l; k < m; k++) s += U_in(k, i) * U_in(k, j);
							f = (s / U_in(i, i))*g;
							for (k = i; k < m; k++) U_in(k, j) += f*U_in(k, i);
						}
						for (j = i; j < m; j++) U_in(j, i) *= g;
					}
					else for (j = i; j < m; j++) U_in(j, i) = 0.0;
					++U_in(i, i);
				}
				for (k = n - 1; k >= 0; k--) {
					for (its = 0; its < 30; its++) {
						flag = true;
						for (l = k; l >= 0; l--) {
							nm = l - 1;
							if (l == 0 || std::abs(rv1(l)) <= eps*anorm) {
								flag = false;
								break;
							}
							if (std::abs(s_in(nm)) <= eps*anorm) break;
						}
						if (flag) {
							c = 0.0;
							s = 1.0;
							for (i = l; i<k + 1; i++) {
								f = s*rv1(i);
								rv1(i) = c*rv1(i);
								if (std::abs(f) <= eps*anorm) break;
								g = s_in(i);
								h = (std::abs(f) > std::abs(g) ? std::abs(f)*sqrt(1.0 + std::pow((std::abs(g) / std::abs(f)), 2)) :
									(std::abs(g) == 0.0 ? 0.0 : std::abs(g)*sqrt(1.0 + std::pow((std::abs(f) / std::abs(g)), 2))));

								/*
								pythag function , cahnged from nrutil.h from numrecipies
								(std::abs(a) > std::abs(b) ? std::abs(a)*sqrt(1.0 + std::pow((std::abs(b) / std::abs(a)),2)) :
								(std::abs(b) == 0.0 ? 0.0 : std::abs(b)*sqrt(1.0 + std::pow((std::abs(a) / std::abs(b)), 2))));
								*/

								s_in(i) = h;
								h = 1.0 / h;
								c = g*h;
								s = -f*h;
								for (j = 0; j < m; j++) {
									y = U_in(j, nm);
									z = U_in(j, i);
									U_in(j, nm) = y*c + z*s;
									U_in(j, i) = z*c - y*s;
								}
							}
						}
						z = s_in(k);
						if (l == k) {
							if (z < 0.0) {
								s_in(k) = -z;
								for (j = 0; j<n; j++) V_in(j, k) = -V_in(j, k);
							}
							break;
						}
						if (its == 29) throw("no convergence in 30 svdcmp iterations");
						x = s_in(l);
						nm = k - 1;
						y = s_in(nm);
						g = rv1(nm);
						h = rv1(k);
						f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0*h*y);
						g = (std::abs(f) > 1.0 ? std::abs(f)*sqrt(1.0 + std::pow((1.0 / std::abs(f)), 2)) : sqrt(1.0 + std::pow(f, 2)));
						f = ((x - z)*(x + z) + h*((y / (f + ((f) >= 0.0 ? std::fabs(g) : -std::fabs(g)))) - h)) / x;
						c = s = 1.0;
						for (j = l; j <= nm; j++) {
							i = j + 1;
							g = rv1(i);
							y = s_in(i);
							h = s*g;
							g = c*g;
							z = (std::abs(f) > std::abs(h) ? std::abs(f)*sqrt(1.0 + std::pow((std::abs(h) / std::abs(f)), 2)) :
								(std::abs(h) == 0.0 ? 0.0 : std::abs(h)*sqrt(1.0 + std::pow((std::abs(f) / std::abs(h)), 2))));
							rv1(j) = z;
							c = f / z;
							s = h / z;
							f = x*c + g*s;
							g = g*c - x*s;
							h = y*s;
							y *= c;
							for (jj = 0; jj<n; jj++) {
								x = V_in(jj, j);
								z = V_in(jj, i);
								V_in(jj, j) = x*c + z*s;
								V_in(jj, i) = z*c - x*s;
							}
							z = (std::abs(f) > std::abs(h) ? std::abs(f)*sqrt(1.0 + std::pow((std::abs(h) / std::abs(f)), 2)) :
								(std::abs(h) == 0.0 ? 0.0 : std::abs(h)*sqrt(1.0 + std::pow((std::abs(f) / std::abs(h)), 2))));
							s_in(j) = z;
							if (z) {
								z = 1.0 / z;
								c = f*z;
								s = h*z;
							}
							f = c*g + s*y;
							x = c*y - s*g;
							for (jj = 0; jj < m; jj++) {
								y = U_in(jj, j);
								z = U_in(jj, i);
								U_in(jj, j) = y*c + z*s;
								U_in(jj, i) = z*c - y*s;
							}
						}
						rv1(l) = 0.0;
						rv1(k) = f;
						s_in(k) = x;
					}
				}
			}
			//End SVD::decompose()

			//Beginn SVD::reorder()
			//I am not sure if this is necessary, if this
			//algorithm merely reorders it is surely not
			//necessary to determine the rank. However, I
			//did not fully read this procedure and am therefore
			//not 100% sure about wtf is even going on today.
			//End SVD::reorder()

			float_type tsh = 0.5*sqrt(m + n + 1.) * s_in(0) * eps;
			int j, nr = 0;
      for (j = 0; j < n; j++)
      {
        if (s_in(j) > sqrt(tsh) * 100 )
        {
          nr++;
        }
      }
			return nr;
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
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: Starting determinant calculation of " << this->rows() << " x " << this->cols() << " matrix." << std::endl;

#ifdef CAST_USE_ARMADILLO
			return static_cast<float_type>(det(*this));
#else
			//Via Numerical Recipies, LU Decomposition
			mathmatrix lu = *this;
			//const float_type TINY = 1.0e-40;
			int i, imax, j, k;
			int n = int(this->rows());
			std::vector<int> indx(n);
			float_type big, temp;
			mathmatrix vv(n, 1u);
			float_type d = 1.0;

			//Start LUdcmp
			{
				for (i = 0; i<n; i++) {
					big = 0.0;
					for (j = 0; j<n; j++)
						if ((temp = std::abs(lu(i, j))) > big) big = temp;
					if (big == 0.0)
					{
						//std::cout << "Singular matrix in LUdcmp (->calculation of determinant)";
						goto IF_SINGULAR;
					}
					vv(i) = 1.0 / big;
				}
				for (k = 0; k<n; k++)
				{
					big = 0.0;
					imax = k;
					for (i = k; i<n; i++) {
						temp = vv(i) * std::abs(lu(i, k));
						if (temp > big) {
							big = temp;
							imax = i;
						}
					}
					if (k != imax)
					{
						for (j = 0; j<n; j++)
						{
							temp = lu(imax, j);
							lu(imax, j) = lu(k, j);
							lu(k, j) = temp;
						}
						d = -d;
						vv(imax) = vv(k);
					}
					indx[k] = imax;

					//If pivot is zero, then the matrix is singular
					if (lu(k, k) == 0.0)
					{
						d = 0.0;
						break;
						//lu(k, k) = TINY;
					}

					for (i = k + 1; i<n; i++) {
						temp = lu(i, k) /= lu(k, k);
						for (j = k + 1; j<n; j++)
							lu(i, j) -= temp*lu(k, j);
					}
				}
			}
			//End LUdcmp

			//Get determinant
			if (d != 0.0)
			{
				for (int i2 = 0; i2 < n; i2++)
				{
					d *= lu(i2, i2);
				}
				return d;
			}
			else
			{
				IF_SINGULAR:
				return 0u;
			}
#endif
		}

		/**
		 * @brief Overload "-" Operator for mathmatrix
		 */
		mathmatrix operator-(mathmatrix const& in) const
		{
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

		/**
		 * @brief Overload "/" Operator for mathmatrix and scalars
		 */
		mathmatrix operator/(T const& in) const
		{
      if (Config::get().general.verbosity > 4U)
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
      if (Config::get().general.verbosity > 4U)
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
      if (Config::get().general.verbosity > 4U)
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
      if (Config::get().general.verbosity > 4U)
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
      if (Config::get().general.verbosity > 4U)
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
      if (Config::get().general.verbosity > 4U)
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
      if (Config::get().general.verbosity > 4U)
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
#ifndef CAST_USE_ARMADILLO
      result = mathmatrix<T>(this->rows(), this->cols(), T(0));
      int n = static_cast<int>(this->rows());
      for (int i = 0; i < n; i++)
        for (int j = 0; j < (i + 1); j++) {
          T s = 0;
          for (int k = 0; k < j; k++)
            s += result(i, k) * result(j, k);
          result(i, j) = (i == j) ?
            sqrt((*this)(i, i) - s) :
            (1.0 / result(j, j) * ((*this)(i, j) - s));
        }
      transpose(result);
#else
      arma::Mat<T> const& a = *this;
      result = mathmatrix<T>(arma::chol(a));
#endif
    }

		/**
		 * @brief Returns whether mathmatrix-obj is quadratic
		 */
		bool return_quadratic() const
		{
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: return_quadratic for mathmatrix." << std::endl;
#ifndef CAST_USE_ARMADILLO
			return is_square(*this);
#else
      return this->rows() == this->cols();
#endif
		}

		/**
		 * @brief Returns upper left submatrix. 
     * If no second argument for the function call,
		 * ie for columns_in is specified, then a quadratic submatrix with rows = columns = rows_in
		 * is yieled.
		 */
		mathmatrix upper_left_submatrix(uint_type rows_in, uint_type columns_in = 0) const
		{
      if (Config::get().general.verbosity > 4U)
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

#ifdef CAST_USE_ARMADILLO
    mathmatrix operator*(mathmatrix const& in) const
    {
      arma::Mat<T> const& a = *this;
      arma::Mat<T> const& b = in;
      return mathmatrix(a*b);
    };
#endif


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
    bool operator== (mathmatrix const& in) const
    {
#ifdef CAST_USE_ARMADILLO
      arma::Mat<T> const& a = *this;
      arma::Mat<T> const& b = in;
      return (approx_equal(a,b,"reldiff", 0.001) );
#else
      if (this->rows() != in.rows() || this->cols() != in.cols()) return false;
      else
      {
        for (size_t i = 0u; i < this->rows(); i++)
        {
          for (size_t j = 0u; j < this->cols(); j++)
          {
            if (abs((*this)(i, j) - in(i, j)) > 0.001 * abs(in(i, j)))
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
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: Starting singular value decomposition of " << U_in.rows() << " x " << U_in.cols() << " matrix." << std::endl;

			//"Constructor"
			V_in.resize(this->cols(), this->cols());
			s_in.resize(this->cols(), 1u);
#ifndef CAST_USE_ARMADILLO
			//Rewritten from NumRecipies
			//Numerical Recipies Nomenclatur, I don't even.... who would name it like that? nnm strcpcps wtf
			int n = int(this->cols());
			int m = int(this->rows());

			float_type eps = std::numeric_limits<float_type>::epsilon();

			//Beginn SVD::decompose()
			{
				bool flag;
				int i, its, j, jj, k, l = 0, nm;
				float_type anorm, c, f, g, h, s, scale, x, y, z;
				mathmatrix rv1(n);
				g = scale = anorm = 0.0;
				for (i = 0; i < n; i++)
				{
					l = i + 2;
					rv1(i) = scale * g;
					g = s = scale = 0.0;
					if (i < m) {
						for (k = i; k < m; k++) scale += std::abs(U_in(k, i));
						if (scale != 0.0) {
							for (k = i; k < m; k++) {
								U_in(k, i) /= scale;
								s += U_in(k, i) * U_in(k, i);
							}
							f = U_in(i, i);
							//g = -SIGN(sqrt(s), f);
							//SIGN(a, b) :::: b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
							g = -1 * (f >= 0.0 ? (sqrt(s) >= 0 ? sqrt(s) : -sqrt(s)) : (sqrt(s) >= 0 ? -sqrt(s) : sqrt(s)));

							h = f*g - s;
							U_in(i, i) = f - g;
							for (j = l - 1; j < n; j++) {
								for (s = 0.0, k = i; k < m; k++) s += U_in(k, i) * U_in(k, j);
								f = s / h;
								for (k = i; k < m; k++)  U_in(k, j) += f* U_in(k, i);
							}
							for (k = i; k < m; k++) U_in(k, i) *= scale;
						}
					}
					s_in(i) = scale *g;
					g = s = scale = 0.0;
					if (i + 1 <= m && i + 1 != n) {
						for (k = l - 1; k < n; k++) scale += std::abs(U_in(i, k));
						if (scale != 0.0) {
							for (k = l - 1; k < n; k++) {
								U_in(i, k) /= scale;
								s += U_in(i, k) * U_in(i, k);
							}
							f = U_in(i, l - 1);
							g = -1 * (f >= 0.0 ? (sqrt(s) >= 0 ? sqrt(s) : -sqrt(s)) : (sqrt(s) >= 0 ? -sqrt(s) : sqrt(s))); //g = -SIGN(sqrt(s), f)
							h = f*g - s;
							U_in(i, l - 1) = f - g;
							for (k = l - 1; k < n; k++) rv1(k) = U_in(i, k) / h;
							for (j = l - 1; j < m; j++) {
								for (s = 0.0, k = l - 1; k < n; k++) s += U_in(j, k) * U_in(i, k);
								for (k = l - 1; k < n; k++) U_in(j, k) += s*rv1(k);
							}
							for (k = l - 1; k < n; k++) U_in(i, k) *= scale;
						}
					}
					anorm = std::max(anorm, (std::abs(s_in(i)) + std::abs(rv1(i))));
				}
				for (i = n - 1; i >= 0; i--) {
					if (i < n - 1) {
						if (g != 0.0) {
							for (j = l; j < n; j++)
								V_in(j, i) = (U_in(i, j) / U_in(i, l)) / g;
							for (j = l; j < n; j++) {
								for (s = 0.0, k = l; k < n; k++) s += U_in(i, k) * V_in(k, j);
								for (k = l; k < n; k++) V_in(k, j) += s*V_in(k, i);
							}
						}
						for (j = l; j < n; j++) V_in(i, j) = V_in(j, i) = 0.0;
					}
					V_in(i, i) = 1.0;
					g = rv1(i);
					l = i;
				}
				for (i = std::min(m, n) - 1; i >= 0; i--) {
					l = i + 1;
					g = s_in(i);
					for (j = l; j < n; j++) U_in(i, j) = 0.0;
					if (g != 0.0) {
						g = 1.0 / g;
						for (j = l; j < n; j++) {
							for (s = 0.0, k = l; k < m; k++) s += U_in(k, i) * U_in(k, j);
							f = (s / U_in(i, i))*g;
							for (k = i; k < m; k++) U_in(k, j) += f*U_in(k, i);
						}
						for (j = i; j < m; j++) U_in(j, i) *= g;
					}
					else for (j = i; j < m; j++) U_in(j, i) = 0.0;
					++U_in(i, i);
				}
				for (k = n - 1; k >= 0; k--) {
					for (its = 0; its < 30; its++) {
						flag = true;
						for (l = k; l >= 0; l--) {
							nm = l - 1;
							if (l == 0 || std::abs(rv1(l)) <= eps*anorm) {
								flag = false;
								break;
							}
							if (std::abs(s_in(nm)) <= eps*anorm) break;
						}
						if (flag) {
							c = 0.0;
							s = 1.0;
							for (i = l; i<k + 1; i++) {
								f = s*rv1(i);
								rv1(i) = c*rv1(i);
								if (std::abs(f) <= eps*anorm) break;
								g = s_in(i);
								h = (std::abs(f) > std::abs(g) ? std::abs(f)*sqrt(1.0 + (std::abs(g) / std::abs(f)) * (std::abs(g) / std::abs(f))) :
									(std::abs(g) == 0.0 ? 0.0 : std::abs(g)*sqrt(1.0 + (std::abs(f) / std::abs(g)) * (std::abs(f) / std::abs(g)))));
								//h=pythag(f,g);
								/*
								pythag function , cahnged from nrutil.h from numrecipies
								(std::abs(a) > std::abs(b) ? std::abs(a)*sqrt(1.0 + std::pow((std::abs(b) / std::abs(a)),2)) :
								(std::abs(b) == 0.0 ? 0.0 : std::abs(b)*sqrt(1.0 + std::pow((std::abs(a) / std::abs(b)), 2))));
								*/

								s_in(i) = h;
								h = 1.0 / h;
								c = g*h;
								s = -f*h;
								for (j = 0; j < m; j++) {
									y = U_in(j, nm);
									z = U_in(j, i);
									U_in(j, nm) = y*c + z*s;
									U_in(j, i) = z*c - y*s;
								}
							}
						}
						z = s_in(k);
						if (l == k) {
							if (z < 0.0) {
								s_in(k) = -z;
								for (j = 0; j<n; j++) V_in(j, k) = -V_in(j, k);
							}
							break;
						}
						if (its == 29) throw("no convergence in 30 svdcmp iterations");
						x = s_in(l);
						nm = k - 1;
						y = s_in(nm);
						g = rv1(nm);
						h = rv1(k);
						f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0*h*y);
						g = (std::abs(f) > 1.0 ? std::abs(f)*sqrt(1.0 + std::pow((1.0 / std::abs(f)), 2)) : sqrt(1.0 + std::pow(f, 2)));
						f = ((x - z)*(x + z) + h*((y / (f +
							//SIGN(a, b) :::: b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
							(f >= 0.0 ? (g >= 0 ? g : -g) : (g >= 0 ? -g : g))
							//((f) >= 0.0 ? std::fabs(g) : -std::fabs(g))
							)) - h)) / x;
						c = s = 1.0;
						for (j = l; j <= nm; j++) {
							i = j + 1;
							g = rv1(i);
							y = s_in(i);
							h = s*g;
							g = c*g;
							z = (std::abs(f) > std::abs(h) ? std::abs(f)*sqrt(1.0 + std::pow((std::abs(h) / std::abs(f)), 2)) :
								(std::abs(h) == 0.0 ? 0.0 : std::abs(h)*sqrt(1.0 + std::pow((std::abs(f) / std::abs(h)), 2))));
							rv1(j) = z;
							c = f / z;
							s = h / z;
							f = x*c + g*s;
							g = g*c - x*s;
							h = y*s;
							y *= c;
							for (jj = 0; jj<n; jj++) {
								x = V_in(jj, j);
								z = V_in(jj, i);
								V_in(jj, j) = x*c + z*s;
								V_in(jj, i) = z*c - x*s;
							}
							z = (std::abs(f) > std::abs(h) ? std::abs(f)*sqrt(1.0 + std::pow((std::abs(h) / std::abs(f)), 2)) :
								(std::abs(h) == 0.0 ? 0.0 : std::abs(h)*sqrt(1.0 + std::pow((std::abs(f) / std::abs(h)), 2))));
							s_in(j) = z;
							if (z) {
								z = 1.0 / z;
								c = f*z;
								s = h*z;
							}
							f = c*g + s*y;
							x = c*y - s*g;
							for (jj = 0; jj < m; jj++) {
								y = U_in(jj, j);
								z = U_in(jj, i);
								U_in(jj, j) = y*c + z*s;
								U_in(jj, i) = z*c - y*s;
							}
						}
						rv1(l) = 0.0;
						rv1(k) = f;
						s_in(k) = x;
					}
				}
			}
			//End SVD::decompose()

			//Beginn SVD::reorder()
			{
				int i, j, k, s, inc = 1;
				float_type sw;
				mathmatrix su(m), sv(n);
				do { inc *= 3; inc++; } while (inc <= n);
				do
				{
					inc /= 3;
					for (i = inc; i < n; i++)
					{
						sw = s_in(i);
						for (k = 0; k < m; k++) su(k) = U_in(k, i);
						for (k = 0; k < n; k++) sv(k) = V_in(k, i);
						j = i;
						while (s_in(j - inc) < sw)
						{
							s_in(j) = s_in(j - inc);
							for (k = 0; k < m; k++) U_in(k, j) = U_in(k, j - inc);
							for (k = 0; k < n; k++) V_in(k, j) = V_in(k, j - inc);
							j -= inc;
							if (j < inc) break;
						}
						s_in(j) = sw;
						for (k = 0; k < m; k++) U_in(k, j) = su(k);
						for (k = 0; k < n; k++) V_in(k, j) = sv(k);
					}
				} while (inc > 1);
				for (k = 0; k < n; k++)
				{
					s = 0;
					for (i = 0; i < m; i++) if (U_in(i, k) < 0.) s++;
					for (j = 0; j < n; j++) if (V_in(j, k) < 0.) s++;
					if (s > (m + n) / 2) {
						for (i = 0; i < m; i++) U_in(i, k) = -U_in(i, k);
						for (j = 0; j < n; j++) V_in(j, k) = -V_in(j, k);
					}
				}
			}
			//End SVD::reorder()

			float_type tsh = 0.5*sqrt(m + n + 1.) * s_in(0) * eps;

			int j, nr = 0;
			for (j = 0; j<n; j++) if (s_in(j) > tsh) nr++;
			if (rank != nullptr)
			{
			  *rank = nr;
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
      if (Config::get().general.verbosity > 4U)
        std::cout << "Function call: eigensym for mathmatrix." << std::endl;
			//On basis of SVD, so take care that your matrix is symmetrical

			//And if you are unsure about your symmetry, try this beforehand:
			//this->symmetry_check();
			//(Although it might slow stuff down considerably.)
#ifndef CAST_USE_ARMADILLO
			mathmatrix V;
			this->singular_value_decomposition(eigenvec_in, eigenval_in, V, rank_in);
#else
      eigenvec_in.resize(this->cols(), this->cols());
      eigenval_in.resize(this->cols(), 1u);
      arma::Col<float_type> s;
      arma::Mat<float_type> eigVecUnordered(eigenvec_in.rows(), eigenvec_in.cols());
      eig_sym(s, eigVecUnordered, *this);

      // Swap the eigenvalues and vectors since they are ordered backwards (we need largest to
      // smallest).
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


		/**
		 * @breif Returns mathmatrix-obj as std vector of vector of float_type.
		 * USE THIS TO DEBUG, not in production code.
		 * This might be useful as the VS debugger cannot visualize content of arma arrays.
     *
     * @return Matrix in std::vector<vector<float-type> > form.
		 */
		std::vector <std::vector<T> > to_std_vector(void) const
		{
      if (Config::get().general.verbosity > 4U)
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

	};


  //FREE FUNCTIONS
#ifdef CAST_USE_ARMADILLO
  template<typename T>
  mathmatrix<T> transposed(mathmatrix<T> const& in)
  {
    return mathmatrix<T>(in.t());
  }

  template<typename T>
  void transpose(mathmatrix<T>& in)
  {
    in = transposed(in);
  }

#else
  using namespace scon;

  template<typename T>
  mathmatrix<T> transposed(mathmatrix<T> const &A)
  {
    matrix<T> const& in = A;
    return mathmatrix<T>(transposed(in));
  }

  // transpose matrix
  template<class T>
  void transpose(mathmatrix<T> &m)
  {
    matrix<T> & in = m;
   (transpose(in));
  }


  template<class T, class U = T>
  typename std::enable_if < std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,
    mathmatrix<typename std::common_type<T, U>::type >> ::type
    operator* (mathmatrix<T> const & a, mathmatrix<U> const &b)
  {
    scon::matrix<T> const &ma = a;
    matrix<U> const &mb = b;
    return ma*mb;
  }

  template<class T, class U = T>
  typename std::enable_if < std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,
    mathmatrix<typename std::common_type<T, U>::type >> ::type
    operator* (matrix<T> const & a, mathmatrix<U> const &b)
  {
    matrix<U> const &mb = b;
    return a *mb;
  }

  template<class T, class U = T>
  typename std::enable_if < std::is_arithmetic<T>::value && std::is_arithmetic<U>::value,
    mathmatrix<typename std::common_type<T, U>::type >> ::type
    operator* (mathmatrix<T> const & a, matrix<U> const &b)
  {
    matrix<T> const &ma = a;
    return ma * b;
  }
#endif

//END HEADER
