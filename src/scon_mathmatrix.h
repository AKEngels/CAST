#pragma once

/**
CAST 3
scon_mathmatrix.h
Purpose: Enabling matrix calculations. Uses Armadillo for enhanced speed when
available. Otherwise uses slow internal routines.

Works on either internal scon::matrix types or arma::Mat types if
the flag USE_ARMADILLO is specified

@author Dustin Kaiser
@version 3.0
*/

/*
USAGE CONVENTIONS AS FOLLOWS:

mathmatrix(xyz, atom_nr) for mathmatrix of one frame in cartesian coordiantes
mathmatrix(dist/angle/dihedral, atom_nr) for mathmatrix of one frame in internal
coordiantes mathmatrix(coords, frames)

*/

/*
CODING CONVENTIONS AS FOLLOWS:

- Wraps underlying abstract matrix-obj of the "scon_matrix.h" - kind or
arma::Mat kind
- There might be still some ambigous stuff going on regarding types.

*/

///////////////////////////////
//                           //
//	I N C L U D E S          //
//                           //
///////////////////////////////

#include <cmath>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

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

#ifdef CAST_USE_ARMADILLO
#include <armadillo>
template <typename T>
using matrix_type = arma::Mat<T>;
#else
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
template <typename T>
using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
#endif

///////////////////////////////

/**
 * @brief Class for handling matrix operations involving numerical entries.
 * Uses Armadillo for enhanced speed when available. Otherwise uses slow
 * internal routines.
 * @author Dustin Kaiser
 * @version 3.0
 *
 * The class scon::mathmatrix is made to be mainly used for matrix operations
 * on floating point numbers such as multiplication, eigenvalue-decomposition
 * and so on.
 *
 * It is written as a generic template but it is entirely untested for types
 * other than floating point numbers.
 *
 * In the file scon_mathmatrix_test.cc there are comprehensive unit tests for
 * the methods of this class
 *
 * @note This class is two-faced. It is a wrapper around either
 * - a Eigen matrix, a fast header only matrix math library for cpp
 * - or a armadillo::Mat. Armadillo is a C++ Framework which itself is a wrapper
 * for LAPACK and BLAS, high speed fortran matrix routines. If the preprocessor
 * define "CAST_USE_ARMADILLO" is set, scon::mathmatrix will accompany an
 * arma::Mat object. This of course implies that CAST then hast to be compiled
 * and linked with pre-existing LAPACK and BLAS libraries. Currently,
 * precompiled versions can be found in the CAST git repository in the folder
 * optional_files. If you use the recomended premake5 build automation to build
 * CAST, everything should be fine and you do not need to worry about linking
 * CAST with LAPACK.
 *
 * For environments where linking with BLAS and LAPACK is not possible, there is
 * also the option of using scon::mathmatrix without armadillo. For this end we
 * have provided a stand-alone, internal version of each matrix-transformation.
 * These are assured to yield similar results. However, especially
 * SVD-Decompositions on large matrices will be painsteakingly slow. It will be
 * too slow for productive use by computational chemists in case of the tasks
 * PCA and ENTROPY.
 *
 * @warning Matrix access does NOT throw when out of bounds! Take caution here,
 * only armadillo-enabled matrices on debug builds will cause an exeption when a
 * matrix is accessed out of bounds (for example, accesing a 2x2 matrix with
 * matrix(4,4)). In ALL other scenarios, CAST will continue without ANY error
 * and you find yourself in undefined bahviour land.
 *
 * Transpose() and transposed() are not members of the mathmatrix class but
 * available as free functions (found at the bottom of the scon_mathmatrix.h
 * file).
 *
 * Usage conventions are as follows (these are guidelines and not enforced by
 * design or assertions):
 *
 * - mathmatrix(xyz, atom_nr) for mathmatrix of one frame in cartesian
 * coordiantes.
 * - mathmatrix(dist/angle/dihedral, atom_nr) for mathmatrix of one frame in
 * internal coordiantes.
 * - mathmatrix(coords, frames) for a matrix of a whole (MD) trajectory.
 *
 */
namespace scon {
template <typename T>
class mathmatrix
#ifndef CAST_USE_ARMADILLO
    : public matrix_type<T>
#else
    : public matrix_type<T>
#endif

{
private:
#ifdef CAST_USE_ARMADILLO
  using base_type = matrix_type<T>;
#else
  using base_type = matrix_type<T>;
#endif

public:
  using int_type = int;
  using uint_type = std::size_t;
  static auto constexpr printFunctionCallVerbosity = 5u;
  static auto constexpr mat_comp_tol = 0.001;

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
  using base_type::Matrix;
  mathmatrix(std::initializer_list<std::initializer_list<T>> ini);
#else
  using base_type::Mat;
#endif

  mathmatrix() = default;
  mathmatrix(mathmatrix const& other)
      : base_type(static_cast<base_type>(other)) {}
  mathmatrix(base_type const& other) : base_type(other) {}

  /*! Construct filled mathmatrix of certain size
   *
   * All values initialized to the same value
   * @param rows: Number of rows
   * @param cols: Number of columns
   * @param fill: Value to which all matrix elements will be initialized
   */
      mathmatrix(uint_type rows, uint_type cols, T fill);

      static mathmatrix zero(uint_type rows, uint_type cols);
      static mathmatrix eye(uint_type rows, uint_type cols);
      static mathmatrix fill_diag(uint_type rows, uint_type cols, T fill);

      // base row and col proxy
      using base_type::row;
      using base_type::col;


      // element access from base class
      using base_type::operator();
      using base_type::operator[];
#ifndef CAST_USE_ARMADILLO
      void resize(uint_type const rows, uint_type const cols);
#else
      using base_type::resize;
#endif


      using base_type::operator*=;
      using base_type::operator-=;
      using base_type::operator/=;
      using base_type::operator+=;
      /*! mathmatrix += operator
       *
       * @param in: Matrix to the right of the summation (this + in)
       * @return: Result of the addition of this + in
       */
      mathmatrix operator+(mathmatrix const& in)const;
      mathmatrix operator*(mathmatrix const& in) const;
      mathmatrix operator/(mathmatrix const& in) const;
      mathmatrix operator-(mathmatrix const& in) const;
      mathmatrix operator/(T const& in) const;
      mathmatrix operator*(T const& in) const;

      /*! Returns an "identity matrix" of certain size
       *
       * All diagonal elements will be initialized to one
       * All off-diagonal elements will be initialized to zero
       * @param num_rows: Number of rows
       * @param num_cols: Number of columns
       * @todo: Write as free function!!
       */
      static typename std::enable_if<std::is_arithmetic<T>::value, mathmatrix>::type
        Identity(std::size_t const num_rows, std::size_t const num_cols);


      // in case you are wondering:
      // transposed and some more stuff is available as free functions
      // eg transpose().

          /////////////////////////////////////
          /////                           /////
          /////  O P E R A T I O N S      /////
          /////                           /////
          /////////////////////////////////////

          /**
           * Checks for positive_definite matrix
           */
      bool positive_definite_check() const;

          /**
           * @brief Returns sign of the determinant (-1 / 1).
       * @return -1 if determinant is negative or zero, +1 if determinant is greater than zero.
           */
          int_type det_sign() const;

      /**
       * @brief Returns rank of the underlying matrix.
   * @return rank of the matrix.
   * @note A singular value decompostition is performed to determine the rank. This may be quite costly.
       */
      std::size_t rank() const;

      /**
       * @brief Calculates determinant of the mathmatrix-obj
   *
   * Internal code uses a LU decompostion.
   * @return determinant
       */
      T determ() const;




      /**
       * @brief Append one matrix to another, will check if sizes match, appends at the bottom end (rows are added)
       */
      void append_bottom(const mathmatrix& I_will_be_the_bottom_part);

      /**
       * @brief Append one matrix to another, will check if sizes match, appends on the top (rows are added)
       */
      void append_top(const mathmatrix& I_will_be_the_top_part);

      /**
       * @brief Append one matrix to another, will check if sizes match, appends left (columns are added)
       */
      void append_left(const mathmatrix& I_will_be_the_left_part);

      /**
     * @brief Append one matrix to another, will check if sizes match, appends left (columns are added)
   */
      void append_right(const mathmatrix& I_will_be_the_right_part);

      /**
       * @brief Sheds the specified rows from the matrix
       */
      void shed_rows(size_t first_in, size_t last_in);

      /**
       * @brief Sheds the specified columns from the matrix
       */
      void shed_cols(size_t first_in, size_t last_in);

      /**
       * @brief Returns number of rows
       */
#ifndef CAST_USE_ARMADILLO
    using base_type::rows;
#else
      std::size_t rows() const;
#endif

      /**
       * @brief Returns number of columns
       */
#ifndef CAST_USE_ARMADILLO
    using base_type::cols;
#else
      std::size_t cols() const;
#endif

    /*! Performs Cholesky Decompostion on Matrix.
     *
     * @NOTE: Code via https://rosettacode.org/wiki/Cholesky_decomposition#C 19.11.16
     * @param result: Upper triangular matrix as result of decompostion
     */
      void choleskyDecomposition(mathmatrix<T> & result) const;

    /**
     * @brief Returns whether mathmatrix-obj is quadratic
     */
    inline bool return_quadratic() const;

    /**
     * @brief Returns upper left submatrix.
 * If no second argument for the function call,
     * ie for columns_in is specified, then a quadratic submatrix with rows = columns = rows_in
     * is yieled.
     */
    mathmatrix upper_left_submatrix(uint_type rows_in, uint_type columns_in = 0) const;

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
    inline bool operator== (mathmatrix const& in) const;


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
    void singular_value_decomposition(mathmatrix& U_in, mathmatrix& s_in, mathmatrix& V_in) const;

    template<typename Comp>
    std::vector<T> sort_col_to_vec(Comp comp, std::size_t const & ind = 0) const;

    template<typename Comp>
    mathmatrix sort_col(Comp comp, std::size_t const & ind = 0) const;

    std::vector<std::size_t> sort_idx(std::size_t const & ind = 0) const;

    template<typename Comp>
    std::vector<std::size_t> find_idx(Comp comp, std::size_t const & ind = 0) const;

#ifndef CAST_USE_ARMADILLO
private:
  void removeRow(std::size_t const & rowToRemove);

  void removeColumn(std::size_t const & colToRemove);
public:
#endif

  mathmatrix submat(std::vector<std::size_t> const & columns, std::vector<std::size_t> const & rows) const;

    mathmatrix pinv() const;

    mathmatrix vectorise(std::size_t const & i = 0)const;

    void replace_idx_with(std::vector<std::size_t> const & idx, T const & val);

    std::vector<std::reference_wrapper<T>> elem(std::vector<std::size_t> const & idx)const;

    mathmatrix diagmat() const;

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
    std::pair<mathmatrix, mathmatrix> eigensym(/*bool const & sort = false*/) const;

    std::pair<mathmatrix, mathmatrix> diag();

        mathmatrix t() const;

        std::vector<T> col_to_std_vector(std::size_t const & iter = 0)const;

        /**
         * @brief Returns mathmatrix-obj as std vector of vector of float_type.
         * USE THIS TO DEBUG, not in production code.
         * This might be useful as the VS debugger cannot visualize content of arma arrays.
     *
     * @return Matrix in std::vector<vector<float-type> > form.
         */
        std::vector<std::vector<T>> to_std_vector() const;

        /**
         * Updates the internal std::vector<std::vector<float_type> > array_debugview_internal array
         * ONLY IF PREPROCESSOR FLAG "DEBUGVIEW" IS SET
         */
        void update_debugview(void)const;

        bool is_vec() const;

#ifndef CAST_USE_ARMADILLO
        using Quaternion = Eigen::Quaternion<T>;
#else
        class Quaternion {

        };
#endif
};

template <typename T>
mathmatrix<T> transpose(mathmatrix<T> const& in) {
  return in.t();
}

/**
 * @brief Rotation Class. If armadillo is enabled it uses the LAPACK matrix
 * routines otherwise it uses the Eigen matrices.
 * @author Julian Erdmannsdörfer
 * @version 3.0
 *
 *
 *
 */

class RotationMatrix {
public:
#ifdef CAST_USE_ARMADILLO
  // Too much work for now.
#else
  using Translation = Eigen::Translation3d;
  using Rotation = Eigen::AngleAxisd;
  using Transformation = Eigen::Affine3d;
  using Vector = Eigen::Vector3d;

  static Transformation
  rotate_around_axis_with_center(double rad_deg, Vector axis, Vector center) {
    Translation back(center);
    Translation to_center(-center);
    return Transformation(back * Rotation(rad_deg, axis.normalized()) *
                          to_center);
  }
  static Transformation rotate_around_axis_in_center(double rad_deg,
                                                     Vector axis) {
    Vector center{ 0., 0., 0. };
    Translation back(center);
    Translation to_center(-center);
    return Transformation(back * Rotation(rad_deg, axis.normalized()) *
                          to_center);
  }

#endif
};

#ifndef CAST_USE_ARMADILLO
template <typename T>
mathmatrix<T>::mathmatrix(std::initializer_list<std::initializer_list<T>> ini) {
  if (ini.size() == 0) {
    *this = mathmatrix();
    return;
  }

  *this = mathmatrix(ini.size(), (*ini.begin()).size());

  auto dat = data();

  for (auto const& i : ini) {
    for (auto const& ii : i) {
      *dat = ii;
      dat++;
    }
  }
}
template <typename T>
void mathmatrix<T>::resize(uint_type const rows, uint_type const cols) {
  this->conservativeResize(rows, cols);
}
#endif

template <typename T>
mathmatrix<T>::mathmatrix(uint_type rows, uint_type cols, T fill)
    : mathmatrix(rows, cols) {
  for (uint_type i = 0u; i < rows; i++)
    for (uint_type j = 0u; j < cols; j++)
      (*this)(i, j) = fill;
}

template <typename T>
mathmatrix<T> mathmatrix<T>::zero(uint_type rows, uint_type cols) {
  return mathmatrix(rows, cols, T());
}

template <typename T>
mathmatrix<T> mathmatrix<T>::eye(uint_type rows, uint_type cols) {
  mathmatrix ret(rows, cols, T());
  for (std::size_t i = 0; i < rows; ++i) {
    ret(i, i) = static_cast<T>(1.0);
  }
  return ret;
}

template <typename T>
mathmatrix<T> mathmatrix<T>::fill_diag(uint_type rows, uint_type cols, T fill) {
  mathmatrix ret(rows, cols, T());
  for (std::size_t i = 0; i < rows; ++i) {
    ret(i, i) = fill;
  }
  return ret;
}

template <typename T>
mathmatrix<T> mathmatrix<T>::operator+(mathmatrix const& in) const {
#ifndef CAST_USE_ARMADILLO
  return base_type::operator+(in);
#else
  return arma::operator+(static_cast<base_type>(*this),
                         static_cast<base_type>(in));
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::operator*(mathmatrix const& in) const {
#ifndef CAST_USE_ARMADILLO
  return base_type::operator*(in);
#else
  return arma::operator*(static_cast<base_type>(*this),
                         static_cast<base_type>(in));
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::operator/(mathmatrix const& in) const {
#ifndef CAST_USE_ARMADILLO
  return base_type::operator/(in);
#else
  return arma::operator/(static_cast<base_type>(*this),
                         static_cast<base_type>(in));
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::operator-(mathmatrix const& in) const {
#ifndef CAST_USE_ARMADILLO
  return base_type::operator-(in);
#else
  return arma::operator-(static_cast<arma::Mat<T>>(*this),
                         static_cast<arma::Mat<T>>(in));
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::operator/(T const& in) const {
#ifndef CAST_USE_ARMADILLO
  return base_type::operator/(in);
#else
  return arma::operator/(static_cast<base_type>(*this), in);
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::operator*(T const& in) const {
#ifndef CAST_USE_ARMADILLO
  return base_type::operator*(in);
#else
  return arma::operator*(static_cast<base_type>(*this), in);
#endif
}

template <typename T>
typename std::enable_if<std::is_arithmetic<T>::value, mathmatrix<T>>::type
mathmatrix<T>::Identity(std::size_t const num_rows,
                        std::size_t const num_cols) {
#ifndef CAST_USE_ARMADILLO
  return mathmatrix(base_type::Identity(num_rows, num_cols));
#else
  return mathmatrix(arma::eye<base_type>(num_rows, num_cols));
#endif
}

template <typename T>
bool mathmatrix<T>::positive_definite_check() const {
  if (!this->return_quadratic())
    return false;
  for (unsigned int i = 1; i < (this->rows() + 1); i++) {
    if (mathmatrix(this->upper_left_submatrix(i, i)).determ() <= 0) {
      return false;
    }
  }
  return true;
}

template <typename T>
int mathmatrix<T>::det_sign() const {
  return ((this->determ() <= 0) ? -1 : 1);
}

template <typename T>
std::size_t mathmatrix<T>::rank() const {
#ifndef CAST_USE_ARMADILLO
  Eigen::ColPivHouseholderQR<base_type> rank_colpivmat(*this);
  return static_cast<size_t>(rank_colpivmat.rank());
#else
  return static_cast<size_t>(arma::rank(static_cast<base_type>(*this)));
#endif
}

template <typename T>
T mathmatrix<T>::determ() const {
#ifdef CAST_USE_ARMADILLO
  return static_cast<T>(det(*this));
#else
  return static_cast<T>(static_cast<base_type const*>(this)->determinant());
#endif
}

template <typename T>
void mathmatrix<T>::append_bottom(const mathmatrix& I_will_be_the_bottom_part) {

  if (this->cols() != I_will_be_the_bottom_part.cols()) {
    throw std::runtime_error("Wrong Matrix size in mathmatrix:append()");
  }
  // Old size needs to be kept
  size_t holder = this->rows();

  this->resize(this->rows() + I_will_be_the_bottom_part.rows(), this->cols());

  // Add "in" to newly created space.
  for (uint_type i = 0; i < I_will_be_the_bottom_part.rows(); i++) {
    for (uint_type j = 0; j < this->cols(); j++) {
      (*this)(i + holder, j) = I_will_be_the_bottom_part(i, j);
    }
  }
}

template <typename T>
void mathmatrix<T>::append_top(const mathmatrix& I_will_be_the_top_part) {
  if (this->cols() != I_will_be_the_top_part.cols()) {
    throw std::runtime_error("Wrong Matrix size in mathmatrix:append()");
  }
  const size_t thisOldRows = this->rows();
  this->resize(this->rows() + I_will_be_the_top_part.rows(), this->cols());

  // Move the entries in the parent matrix downward
  // We count downwards so that we dont overwrite
  for (uint_type i = thisOldRows - 1u; i > 0; i--) {
    for (uint_type j = 0; j < this->cols(); j++) {
      (*this)(i + I_will_be_the_top_part.rows(), j) = (*this)(i, j);
    }
  }

  // Add "I_will_be_the_top_part" to now absolete top space of the parent matrix
  // (this).
  for (uint_type i = 0; i < I_will_be_the_top_part.rows(); i++) {
    for (uint_type j = 0; j < this->cols(); j++) {
      (*this)(i, j) = I_will_be_the_top_part(i, j);
    }
  }
}

template <typename T>
void mathmatrix<T>::append_left(const mathmatrix& I_will_be_the_left_part) {
  if (this->rows() != I_will_be_the_left_part.rows()) {
    throw std::runtime_error("Wrong Matrix size in mathmatrix:append()");
  }
  const size_t thisOldCols = this->cols();

  this->resize(this->rows(), this->cols() + I_will_be_the_left_part.cols());

  // Move the entries in the parent matrix rightward
  // We count right so that we dont overwrite
  for (uint_type j = thisOldCols - 1u; j > 0; j--) {
    for (uint_type i = 0u; i < this->rows(); i++) {
      (*this)(i, j + I_will_be_the_left_part.cols()) = (*this)(i, j);
    }
  }

  // Add "I_will_be_the_left_part" to now absolete top space of the parent
  // matrix (this).
  for (uint_type j = 0; j < I_will_be_the_left_part.cols(); j++) {
    for (uint_type i = 0; i < this->rows(); i++) {
      (*this)(i, j) = I_will_be_the_left_part(i, j);
    }
  }
}

template <typename T>
void mathmatrix<T>::append_right(const mathmatrix& I_will_be_the_right_part) {

  if (this->rows() != I_will_be_the_right_part.rows()) {
    throw std::runtime_error("Wrong Matrix size in mathmatrix:append()");
  }

  // Old size needs to be kept
  uint_type holder = this->cols();

  this->resize(this->rows(), this->cols() + I_will_be_the_right_part.cols());

  // Add "in" to newly created space.
  for (uint_type j = 0; j < I_will_be_the_right_part.cols(); j++) {
    for (uint_type i = 0; i < this->rows(); i++) {
      (*this)(i, j + holder - 1) = I_will_be_the_right_part(i, j);
    }
  }
}

template <typename T>
void mathmatrix<T>::shed_rows(size_t first_in, size_t last_in) {
  if (first_in > last_in || last_in >= this->rows()) {
    throw std::runtime_error("Index Out of Bounds in mathmatrix:shed_rows()");
  }
  mathmatrix newOne(this->rows() - (last_in - first_in + 1u), this->cols());
  for (size_t i = 0u; i < first_in; i++) {
    for (size_t j = 0u; j < this->cols(); j++) {
      newOne(i, j) = (*this)(i, j);
    }
  }
  for (size_t i = first_in; i < this->rows() - last_in - 1 + first_in; i++) {
    for (size_t j = 0u; j < this->cols(); j++) {
      newOne(i, j) = (*this)(i + (last_in - first_in) + 1, j);
    }
  }
  this->swap(newOne);
}

template <typename T>
void mathmatrix<T>::shed_cols(size_t first_in, size_t last_in) {
  if (last_in >= this->cols() || first_in > last_in) {
    throw std::runtime_error("Index Out of Bounds in mathmatrix:shed_rows()");
  }
  mathmatrix newOne(this->rows(), this->cols() - (last_in - first_in + 1u));
  for (size_t j = 0u; j < this->rows(); j++) {
    for (size_t i = 0u; i < first_in; i++) {
      newOne(j, i) = (*this)(j, i);
    }
  }
  for (size_t j = 0u; j < this->rows(); j++) {
    for (size_t i = first_in; i < this->cols() - last_in - 1 + first_in; i++) {
      newOne(j, i) = (*this)(j, i + (last_in - first_in) + 1);
    }
  }
  this->swap(newOne);
}

#ifdef CAST_USE_ARMADILLO
template <typename T>
std::size_t mathmatrix<T>::rows() const {
  return this->n_rows;
}
#endif

#ifdef CAST_USE_ARMADILLO
template <typename T>
std::size_t mathmatrix<T>::cols() const {
  return this->n_cols;
}
#endif

template <typename T>
void mathmatrix<T>::choleskyDecomposition(mathmatrix& result) const {
  result = mathmatrix(this->rows(), this->cols(), T(0));
  int n = static_cast<int>(this->rows());
  for (int i = 0; i < n; i++)
    for (int j = 0; j < (i + 1); j++) {
      T s = 0;
      for (int k = 0; k < j; k++)
        s += result(i, k) * result(j, k);
      result(i, j) = (i == j) ? ::sqrt((*this)(i, i) - s)
                              : (1.0 / result(j, j) * ((*this)(i, j) - s));
    }
#if defined(_MSC_VER) && !defined(CAST_USE_ARMADILLO)
  ::transpose(result);
#else
  transpose(result);
#endif
}
template <typename T>
bool mathmatrix<T>::return_quadratic() const {
  return this->rows() == this->cols();
}

template <typename T>
mathmatrix<T>
mathmatrix<T>::upper_left_submatrix(uint_type rows_in,
                                    uint_type columns_in = 0) const {
#ifndef CAST_USE_ARMADILLO
  return this->block(0, 0, rows_in, columns_in == 0 ? rows_in : columns_in);
#else
  return this->submat(0, 0, rows_in - 1,
                      columns_in == 0 ? rows_in - 1 : columns_in - 1);
#endif
}

template <typename T>
bool mathmatrix<T>::operator==(mathmatrix const& in) const {
#ifndef CAST_USE_ARMADILLO
  return this->isApprox(in, mat_comp_tol);
#else
  return (approx_equal(static_cast<base_type>(*this), static_cast<base_type>(b),
                       "reldiff", mat_comp_tol));
#endif
}

template <typename T>
void mathmatrix<T>::singular_value_decomposition(mathmatrix& U_in,
                                                 mathmatrix& s_in,
                                                 mathmatrix& V_in) const {

#ifndef CAST_USE_ARMADILLO
  auto svd = this->bdcSvd();
  svd.compute(static_cast<base_type>(*this),
              Eigen::ComputeFullU | Eigen::ComputeFullV);
  U_in = svd.matrixU();
  V_in = svd.matrixV();

  s_in = static_cast<base_type>(svd.singularValues());

#else
  arma::Col<T> s;
  if (!svd_econ(U_in, s, V_in, *this))
    throw std::runtime_error("Error in armadillo SVD: failed.");
  for (size_t i = 0; i < U_in.rows(); i++) {
    s_in(i) = s(i);
  }
#endif
}

template <typename T>
template <typename Comp>
std::vector<T> mathmatrix<T>::sort_col_to_vec(Comp comp,
                                              std::size_t const& ind) const {
  auto ret = col_to_std_vector(ind);

  std::sort(ret.begin(), ret.end(), comp);

  return ret;
}

template <typename T>
template <typename Comp>
mathmatrix<T> mathmatrix<T>::sort_col(Comp comp, std::size_t const& ind) const {

  auto sorted_vec = sort_col_to_vec(comp, ind);

  std::vector<std::vector<T>> ret;
  ret.reserve(val_vec.size());
  for (auto const& val : val_vec) {
    ret.emplace_back(std::vector<T>{
        val,
    });
  }

  return mathmatrix(ret);
}

template <typename T>
std::vector<std::size_t> mathmatrix<T>::sort_idx(std::size_t const& ind) const {
  auto val_vec = col_to_std_vector(ind);
  std::vector<std::size_t> ret(val_vec.size());
  std::iota(ret.begin(), ret.end(), 0);

  std::sort(ret.begin(), ret.end(),
            [&val_vec](std::size_t const& i, std::size_t const& j) {
              return val_vec.at(i) < val_vec.at(j);
            });

  return ret;
}

template <typename T>
template <typename Comp>
std::vector<std::size_t> mathmatrix<T>::find_idx(Comp comp,
                                                 std::size_t const& ind) const {
  auto val_vec = col_to_std_vector(ind);
  std::vector<std::size_t> ret(val_vec.size());
  std::iota(ret.begin(), ret.end(), 0);

  ret.erase(
      std::remove_if(ret.begin(), ret.end(),
                     [&](std::size_t const& i) { return !comp(val_vec[i]); }),
      ret.end());

  return ret;
}

#ifndef CAST_USE_ARMADILLO
template <typename T>
void mathmatrix<T>::removeRow(std::size_t const& rowToRemove) {
  auto numRows = rows() - 1;
  auto numCols = cols();

  if (rowToRemove < numRows) {
    base_type::block(rowToRemove, 0, numRows - rowToRemove, numCols) =
        base_type::block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
  }

  base_type::conservativeResize(numRows, numCols);
}

template <typename T>
void mathmatrix<T>::removeColumn(std::size_t const& colToRemove) {
  unsigned int numRows = rows();
  unsigned int numCols = cols() - 1;

  if (colToRemove < numCols) {
    base_type::block(0, colToRemove, numRows, numCols - colToRemove) =
        base_type::block(0, colToRemove + 1, numRows, numCols - colToRemove);
  }

  base_type::conservativeResize(numRows, numCols);
}
#endif

template <typename T>
mathmatrix<T>
mathmatrix<T>::submat(std::vector<std::size_t> const& columns,
                      std::vector<std::size_t> const& rows) const {

#ifndef CAST_USE_ARMADILLO
  mathmatrix ret(*this);
  for (auto const& i : columns) {
    ret.removeColumn(i + 1);
  }
  for (auto const& i : rows) {
    ret.removeRow(i + 1);
  }
  return ret;
#else
  return base_type::submat(columns, rows);
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::pinv() const {
#ifndef CAST_USE_ARMADILLO
  mathmatrix U, V, s;
  singular_value_decomposition(U, s, V);
  mathmatrix s_inv(cols(), rows(), T());

  for (std::size_t i = 0; i < s.rows(); ++i) {
    s_inv(i, i) = 1. / s(i);
  }

  return V * s_inv * U.t();
#else
  return arma::pinv(static_cast<base_type>(*this));
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::vectorise(std::size_t const& i) const {
#ifndef CAST_USE_ARMADILLO
  return mathmatrix(Eigen::Map<Eigen::RowVectorXd>(this->data(), this->size()));
#else
  return arma::vectorise(*this, i);
#endif
}

template <typename T>
void mathmatrix<T>::replace_idx_with(std::vector<std::size_t> const& idx,
                                     T const& val) {
#ifndef CAST_USE_ARMADILLO
  auto dat = data();
#else
  auto dat = memptr();
#endif
  for (auto const& i : idx) {
    *(dat + i) = val;
  }
}

template<typename T>
inline std::vector<std::reference_wrapper<T>> mathmatrix<T>::elem(std::vector<std::size_t> const & idx)const
{
  std::vector<std::reference_wrapper<T>> ret;
#ifndef CAST_USE_ARMADILLO
  auto dat = data();
#else
  auto dat = memptr();
#endif
  for (auto const& i : idx) {
    ret.emplace_back(*(dat + i));
  }
  return ret;
}

template <typename T>
mathmatrix<T> mathmatrix<T>::diagmat() const {
#ifndef CAST_USE_ARMADILLO
  if (cols() == 1) {
    return this->asDiagonal();
  }
  return this->diagonal().asDiagonal();
#else
  return arma::diagmat(*this);
#endif
}

template <typename T>
std::pair<mathmatrix<T>, mathmatrix<T>>
mathmatrix<T>::eigensym(/*bool const & sort = false*/) const {
#ifndef CAST_USE_ARMADILLO

  /*auto sort_eigenpairs = [](auto & EVals, auto & EVecs) {
  std::cout << "Hallo" << std::endl;
  std::vector<std::size_t> perms(EVals.cols());
  std::iota(perms.begin(), perms.end(), 0);

  std::sort(perms.begin(), perms.end(), [&](std::size_t const & i, std::size_t
  const & j) { std::cout << std::boolalpha << (EVals(i, 0) < EVals(j, 0)) <<
  std::endl; return EVals(i, 0) < EVals(j, 0);
  });

  auto EVals_temp = EVals, EVecs_temp = EVecs;
  for (auto i = 0; i < perms.size(); ++i) {

  EVals(i, 0) = EVals_temp(perms.at(i),0);
  EVecs.col(i) = EVecs_temp.col(i);
  }
  };*/

  Eigen::EigenSolver<base_type> es(static_cast<base_type>(*this));
  mathmatrix eigenval = es.eigenvalues().real();
  mathmatrix eigenvec = es.eigenvectors().real();

  // std::cout << eigenval << "\n\n";

  // sort_eigenpairs(eigenval, eigenvec);
#else
  arma::Col<T> eigVal;
  arma::Mat<T> eigVec;
  eig_sym(eigVal, eigVec, *this);

  mathmatrix eigenval = mathmatrix(eigVal);
  mathmatrix eigenvec = mathmatrix(eigVec);
#endif
  return std::make_pair(eigenval, eigenvec);
}

template <typename T>
std::pair<mathmatrix<T>, mathmatrix<T>> mathmatrix<T>::diag() {
  auto ret = eigensym();

  auto const& EVal = ret.first;
  auto const& EVec = ret.second;

  *this = EVec.t() * (*this) * EVec;
  return ret;
}

template <typename T>
mathmatrix<T> mathmatrix<T>::t() const {
#ifndef CAST_USE_ARMADILLO
  return this->transpose();
#else
  return arma::trans(static_cast<base_type>(*this));
#endif
}

template <typename T>
std::vector<T> mathmatrix<T>::col_to_std_vector(std::size_t const& iter) const {
  if (cols() != 1 && iter == 0) {
    throw std::runtime_error("This is not a vector it is a matrix!");
  }
  std::vector<T> ret;
  ret.reserve(rows());
#ifndef CAST_USE_ARMADILLO
  auto const& dat = data();
  auto const& r = rows();
  return std::vector<T>(dat + r * iter, dat + r * (iter + 1));
#else
  return arma::conv_to<std::vector<T>>::from(row(iter))
#endif
}

template <typename T>
std::vector<std::vector<T>> mathmatrix<T>::to_std_vector() const {
  std::vector<std::vector<T>> ret;
  ret.reserve(rows());
#ifndef CAST_USE_ARMADILLO
  auto dat = data();
  auto const& c = cols();
  for (auto i = 0; i < rows(); ++i) {
    ret.emplace_back(std::vector<double>(dat, dat + c));
    dat += c;
  }
#else
  for (auto i = 0; i < rows(); ++i) {
    ret.emplace_back(arma::conv_to<std::vector<T>>::from(row(i)));
  }
#endif
  return ret;
}

template <typename T>
void mathmatrix<T>::update_debugview(void) const {
  if (Config::get().general.verbosity >= printFunctionCallVerbosity)
    std::cout << "Function call: update_debugview." << std::endl;
#ifdef DEBUGVIEW
  std::vector<T> temp1(this->cols());
  std::vector<std::vector<T>> temp2(this->rows(), temp1);
  for (unsigned int i = 0; i < this->rows(); i++) {
    for (unsigned int j = 0; j < this->cols(); j++) {
      temp2[i][j] = (*this)(i, j);
    }
  }
  array_debugview_internal = temp2;
#endif
#ifndef DEBUGVIEW
  std::cout << "DEBUGVIEW Flag not enabled. update_debugview(void) may not be "
               "used.\nCheck your code.\n";
#endif
}

template <typename T>
bool mathmatrix<T>::is_vec() const {
  return cols() == 1;
}

} // namespace scon
// END HEADER
