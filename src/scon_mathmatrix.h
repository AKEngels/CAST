/**
CAST 3
scon_mathmatrix.h
Purpose: Enabling matrix calculations

@author Dustin Kaiser
@version 2.0
*/

/*
USAGE CONVENTIONS AS FOLLOWS:

mathmatrix(xyz, atom_nr) for mathmatrix of one frame in cartesian coordiantes
mathmatrix(dist/angle/dihedral, atom_nr) for mathmatrix of one frame in internal coordiantes
mathmatrix(coords, frames)
mathmatrix(foo,1) == mathmatrix(foo) for single row vectors and such

*/

/* 
CODING CONVENTIONS AS FOLLOWS:

- Wraps underlying abstract matrix-obj of the "scon_matrix.h" - kind
- There might be still some ambigous stuff going on regarding types.
- If it breaks you get to keep all the peaces,

*/

#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <cmath>
#include <limits>
#include <utility>

#include "scon_matrix.h"


/////////////////////////////////
//                             //
//	D E F S                    //
//                             //
/////////////////////////////////

//typedef scon::matrix<coords::float_type> ArrayType;
#include "coords.h"
typedef coords::float_type float_type;
typedef int int_type;
//typedef scon::mathmatrix<float_type> Matrix_Class;
typedef unsigned int uint_type;

///////////////////////////////


///////////////////////////////
//                           //
//	F L A G S                //
//                           //
///////////////////////////////

// ENABLE DEBUGVIEW? 
// (became pretty obsolete after removing armadillo matrix library)
// -> Option to store array data as member of a private std::vector<std::vector<T>> for 
// debug pruproses.

//#define DEBUGVIEW

///////////////////////////////
namespace scon {
	template <typename T> 
  class mathmatrix : public scon::matrix<T>
	{
  private:
    using base_type = scon::matrix<T>;
#ifdef DEBUGVIEW
    std::vector<std::vector<T> > array_debugview_internal;
#else
#endif

	public:

		/////////////////////////////////////
		/////                           /////
		/////  C O N S T R U C T O R S  /////
		/////                           /////
		/////////////////////////////////////

    // pull in matrix constructors
    using base_type::base_type;

    // default constructor
    mathmatrix() = default;

    // construct from matrix
    mathmatrix(base_type const &m) : base_type(m) {}
    mathmatrix(base_type &&m) : base_type(std::move(m)) {}

    // pull in range functions
    using base_type::begin;
    using base_type::end;
    using base_type::cbegin;
    using base_type::cend;
    using base_type::rbegin;
    using base_type::rend;
    using base_type::crbegin;
    using base_type::crend;

    // base row and col proxy
    using base_type::row;
    using base_type::col;

    // element access
    using base_type::operator();
    using base_type::operator[];
    using base_type::operator*=;

    // identity
    using base_type::identity;

		/**
		 * Second "Copy Constructor", if boolean size_only is set to TRUE
		 * Only the size of the matrix will be copied (faster)
		 */
		mathmatrix(mathmatrix const& in, bool size_only) : 
      scon::matrix<T>(in.rows(), in.cols())
		{
			if (!size_only) std::copy(in.begin(), in.end(), this->begin());
		};

		/////////////////////////////////////
		/////                           /////
		/////  O P E R A T I O N S      /////
		/////                           /////
		/////////////////////////////////////


		/**
		 * Overload += operator
		 */
		mathmatrix operator+(mathmatrix const& in) const
		{
			if (!(this->rows() == in.rows() && this->cols() == in.cols() ))
			{
				throw("ERROR in mathmatrix Addition: Sizes of matrices do not match!");
			}
			mathmatrix out(*this);
			for (unsigned int i = 0; i < this->rows(); i++)
			{
				for (unsigned int j = 0; j < this->cols(); j++)
				{
					out(i, j) = out(i, j) + in(i, j);
				}
			}
			return out;
		};

    using base_type::resize;


		/**
		 * Overload * for matrix multiplication
		 *
		 * DEV-NOTE: MATRIX INTERNALS PRESENT
		 */
		mathmatrix operator* (mathmatrix const& in) const
		{
      // get reference to base class
      // (mathmatrix IS a matrix)
      base_type const & a = *this;
      base_type const & b = in;
      // return mathmatrix from multiplication of base types
      return mathmatrix{ a * b };
		};

		/**
		 * Overload * for scalar multiplication
		 */
		mathmatrix operator*(T const& in) const
		{
      base_type const & a = *this;
      return mathmatrix{ a*in };
		};

    

		/**
		 * Prints contents to std::cout
		 */
		void out() const
		{
      base_type const & a = *this;
      matrix_printer<base_type> mp(a);
      // 3 newlines on top
      mp.margin_top = 3u;
      // 3 newlines after
      mp.margin_bottom = 3u;
      // two newlines after each row
      mp.fill_right = '\n';
      mp.margin_right = 1u;
      // 
      std::cout << mp;
		};

		/**
		 * Checks for positive_definite matrix
		 */
		bool positive_definite_check()
		{
			bool positive_definite = true;
      if (!is_square(*this)) return false;
      for (unsigned int i = 1; i < (this->rows() + 1); i++)
      {
        if (mathmatrix{ resized(*this, i, i) }.det_sign() < 0)
        {
          return false;
        }
      }
			return true;
		};

		/**
		 * Returns dign of the determinant (-1 / 1)
		 */
		int det_sign() const
		{
			return ((this->determ() < 0) ? -1 : 1);
		};

		/**
		 * Returns rank of the underlying matrix
		 */
		unsigned int rank() const
		{
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
				float_type anorm, c, f, g, h, s, scale, x, y, z;
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
			for (j = 0; j<n; j++) if (s_in(j) > tsh) nr++;
			return nr;
		};

		/**
		 * Returns determinant of the mathmatrix-obj
		 */
		float_type determ() const
		{
			//Via Numerical Recipies, LU Decomposition
			mathmatrix lu = *this;
			const float_type TINY = 1.0e-40;
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
						//std::cerr << "Singular matrix in LUdcmp (->calculation of determinant)";
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
				for (int i = 0; i < n; i++)
				{
					d *= lu(i, i);
				}
				return d;
			}
			else
			{
				IF_SINGULAR:
				return 0u;
			}
		}

		/**
		 * Overload "-" Operator for mathmatrix
		 */
		mathmatrix operator-(mathmatrix const& in) const
		{
			//Check if sizes match
			if ((in.rows() != this->rows()) || (in.cols() != this->cols()))
			{
				throw ("Error in Matrix mathmatrix subtraction, wrong input sizes");
			}
			mathmatrix output = *this;
			for (unsigned int i = 0; i < this->rows(); i++)
			{
				for (unsigned int j = 0; j < this->cols(); j++)
				{
					output(i, j) -= in(i, j);
				}
			}
			return output;
		}

		/**
		 * Overload "/" Operator for mathmatrix and scalars
		 */
		mathmatrix operator/(T const& in) const
		{
			mathmatrix tempCopy(*this);
			for (unsigned int i = 0; i < this->rows(); i++)
			{
				for (unsigned int j = 0; j < this->cols(); j++)
				{
					tempCopy(i, j) = (*this)(i, j) / in;
				}
			}
			return tempCopy;
		}

		/**
		 * Append one matrix to another, will check if sizes match, appends at the bottom end (rows are added)
		 */
		void append_bottom(const mathmatrix& I_will_be_the_bottom_part)
		{
			//if this is transposed, append bottom means append right on underlying obj
			if (this->cols() != I_will_be_the_bottom_part.cols())
			{
				throw "Wrong Matrix size in mathmatrix:append()";
			}
			//Old size needs to be kept
			unsigned int holder = this->rows();

			this->resize(this->rows() + I_will_be_the_bottom_part.rows(), this->cols());

			//Add "in" to newly ceated space.
			for (unsigned int i = 0; i < I_will_be_the_bottom_part.rows(); i++)
			{
				for (unsigned int j = 0; j < this->cols(); j++)
				{
					(*this)(i + holder, j) = I_will_be_the_bottom_part(i, j);
				}
			}
		}

		/**
		 * Append one matrix to another, will check if sizes match, appends on the top (rows are added)
		 */
		void append_top(const mathmatrix& I_will_be_the_top_part)
		{
			if (this->columns() != I_will_be_the_top_part.columns())
			{
				throw "Wrong Matrix size in mathmatrix:append()";
			}
			this->resize(this->rows() + I_will_be_the_top_part.rows(), this->cols());

			//Move the entries in the parent matrix downward
			//We count downwards so that we dont overwrite
			for (unsigned int i = this->rows() - 1u; i > 0; i--)
			{
				for (unsigned int j = 0; j < this->cols(); j++)
				{
					(*this)(i + I_will_be_the_top_part.rows(), j) = (*this)(i, j);
				}
			}

			//Add "I_will_be_the_top_part" to now absolete top space of the parent matrix (this).
			for (unsigned int i = 0; i < I_will_be_the_top_part.rows(); i++)
			{
				for (unsigned int j = 0; j < this->cols(); j++)
				{
					(*this)(i, j) = I_will_be_the_top_part(i, j);
				}
			}
		}

		/**
		 * Append one matrix to another, will check if sizes match, appends left (columns are added)
		 */
		void append_left(const mathmatrix& I_will_be_the_left_part)
		{
			if (this->columns() != I_will_be_the_left_part.columns())
			{
				throw "Wrong Matrix size in mathmatrix:append()";
			}

			this->resize(this->rows(), this->cols() + I_will_be_the_left_part.cols());

			//Move the entries in the parent matrix downward
			//We count downwards so that we dont overwrite
			for (unsigned int j = this->cols() - 1u; j > 0; j--)
			{
				for (unsigned int i = 0u; i < this->rows(); i++)
				{
					(*this)(i, j + I_will_be_the_left_part.cols()) = (*this)(i, j);
				}
			}

			//Add "I_will_be_the_left_part" to now absolete top space of the parent matrix (this).
			for (unsigned int j = 0; j < I_will_be_the_left_part.cols(); j++)
			{
				for (unsigned int i = 0; i < this->rows(); i++)
				{
					(*this)(i, j) = (*this)(i, j);
				}
			}
		}

		/**
	   * Append one matrix to another, will check if sizes match, appends left (columns are added)
  	 */
		void append_right(const mathmatrix& I_will_be_the_right_part)
		{
			if (this->rows() != I_will_be_the_right_part.rows())
			{
				throw "Wrong Matrix size in mathmatrix:append()";
			}
			//Old size needs to be kept
			unsigned int holder = this->cols();

			this->resize(this->rows(), this->cols() + I_will_be_the_right_part.cols());

			//Add "in" to newly created space.
			for (unsigned int j = 0; j < I_will_be_the_right_part.cols(); j++)
			{
				for (unsigned int i = 0; i < this->rows(); i++)
				{
					(*this)(i, j + holder) = I_will_be_the_right_part(i, j);
				}
			}

		}

		/**
		 * Sheds the specified rows from the matrix
		 */
		void shed_rows(unsigned int first_in, unsigned int last_in)
		{
			mathmatrix newOne(this->rows() - (last_in - first_in + 1u), this->cols());
			for (unsigned int i = 0u; i < first_in; i++)
			{
				for (unsigned int j = 0u; j < this->cols(); j++)
				{
					newOne(i, j) = (*this)(i, j);
				}
			}
			for (unsigned int i = last_in + 1u; i < this->rows(); i++)
			{
				for (unsigned int j = 0u; j < this->cols(); j++)
				{
					newOne(i - last_in - 1u, j) = (*this)(i, j);
				}
			}
			this->swap(newOne);
		}

		/**
		 * Sheds the specified columns from the matrix
		 */
		void shed_cols(unsigned int first_in, unsigned int last_in)
		{
			mathmatrix newOne(this->rows(), this->cols() - (last_in - first_in + 1u));
			for (unsigned int j = 0u; j < this->rows(); j++)
			{
				for (unsigned int i = 0u; i < first_in; i++)
				{
					newOne(j, i) = (*this)(j, i);
				}
			}
			for (unsigned int j = 0u; j < this->rows(); j++)
			{
				for (unsigned int i = last_in + 1u; i < this->cols(); i++)
				{
					newOne(j, i - last_in - 1u) = (*this)(j, i);
				}
			}
			this->swap(newOne);
		}

		/**
		 * Returns number of rows
		 */
    using base_type::rows;

		/**
		 * Returns number of columns
		 */
    using base_type::cols;

		/**
		 * Returns whether mathmatrix-obj is quadratic
		 */
		bool return_quadratic() const
		{
			return this->quadratic();
		}

		/**
		 * Returns upper left submatrix. if no second argument for the function call,
		 * ie for columns_in is specified, then a quadratic submatrix with rows = columns = rows_in
		 * is yieled.
		 */
		mathmatrix upper_left_submatrix(unsigned int const& rows_in, unsigned int const& columns_in = 0) const
		{
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

		/**
		 * Performs singular value decomposition on *this and writes results
		 * to the three resulting matrices U, s, V. Since the rank of a matrix
		 * is automatically calculated during our SVD-decomposition, a variable "rank"
		 * may be specified where the calculated rank may be stored for future use.
     *
     * justification: ptr is used for rank because in can represent a
     * bool intrinsically through nullptr and therefore stores information
     * wether rank should be passed through this function.
		 */
		void singular_value_decomposition(mathmatrix& U_in, mathmatrix& s_in, mathmatrix& V_in, int* rank = nullptr) const
		{
			//Rewritten from NumRecipies

			//"Constructor"
			U_in = *this;
			V_in.resize(this->cols(), this->cols());
			s_in.resize(this->cols(), 1u);

			//Numerical Recipies Nomenclatur, I don't even.... who would name it like that? nnm strcpcps wtf
			int n = int(this->cols());
			int m = int(this->rows());

			float_type eps = std::numeric_limits<float_type>::epsilon();

			//Beginn SVD::decompose()
			{
				bool flag;
				int i, its, j, jj, k, l, nm;
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
		}

		/**
		 * Performs eigenvalue decomposition on *this and returns eigenval and eigenvec
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
		 * NOTE: Should we check if this is symmetric? I think we should -> "eig_sym",
		 * at least once we will do matrix-math internally in CAST.
		 */
		void eigensym(mathmatrix& eigenval_in, mathmatrix& eigenvec_in, int* rank_in = nullptr)
		{
			//On basis of SVD, so take care that your matrix is symmetrical
			//otherwise, well, you know, shit in -> shit out...

			//And if you are unsure about your symmetry, try this beforehand:
			//this->symmetry_check();
			//(Although it might slow stuff down considerably.)

			mathmatrix V;
			this->singular_value_decomposition(eigenvec_in, eigenval_in, V, rank_in);
		}

		/////////////////////////////////////
		/////                           /////
		/////   D E B U G  /  D E V     /////
		/////                           /////
		/////////////////////////////////////

		/**
		 * Returns mathmatrix-obj as vector of vector of float_type
		 * USE THIS TO DEBUG
		 * since the VS debugger cannot visualize content of arma arrays.
		 */
		std::vector <std::vector<T> > to_std_vector(void) const
		{
			std::vector <T> temp1(this->cols());
			std::vector < std::vector <T> > temp2(this->rows(), temp1);
			for (unsigned int i = 0; i < this->rows(); i++)
			{
				for (unsigned int j = 0; j < this->cols(); j++)
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

		/**
		 * Overload "==" Operator for mathmatrix
		 */
		bool operator==(mathmatrix const& in) const
		{
      base_type const& a = *this;
      base_type const& b = in;
      return a == b;
		}

	};

}





//END HEADER