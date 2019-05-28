#pragma once

/**
CAST 3
scon_mathmatrix.h
Purpose: Enabling matrix calculations. Uses Armadillo for enhanced speed when
available. Otherwise uses slow internal routines.

Works on either internal scon::matrix types or arma::Mat types if
the flag USE_ARMADILLO is specified

@author Julian Erdmannsdörfer, Dustin Kaiser
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
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "../configuration.h"

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
//WARNING!
//The Storage-Type is ColMajor by default. Default is used. Be careful changing this!
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
class mathmatrix : public matrix_type<T>{

private:
  using base_type = matrix_type<T>;

public:
  using int_type = int;
  using uint_type = std::size_t;
  static auto constexpr printFunctionCallVerbosity = 5u;
  static auto constexpr matCompTol(){
	  return 0.1;
  }
  static auto constexpr close_to_zero_tol = 1.e-7;

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
#else
  using base_type::Mat;
#endif

  mathmatrix(std::initializer_list<std::initializer_list<T>> const& ini);
  mathmatrix() = default;
  mathmatrix(mathmatrix const& other)
      : base_type(static_cast<base_type>(other)) {}
  mathmatrix(base_type const& other) : base_type(other) {}

  static mathmatrix col_from_vec(std::vector<T> const& col);
  static mathmatrix row_from_vec(std::vector<T> const& row);

  /*! Construct filled mathmatrix of certain size
   *
   * All values initialized to the same value
   * @param rows: Number of rows
   * @param cols: Number of columns
   * @param fill: Value to which all matrix elements will be initialized
   */
  mathmatrix(uint_type rows, uint_type cols, T fill);
  mathmatrix(T, T, T) = delete;

  static mathmatrix zero(uint_type rows, uint_type cols);
  /*! Returns an "identity matrix" of certain size
   *
   * All diagonal elements will be initialized to one
   * All off-diagonal elements will be initialized to zero
   * @param num_rows: Number of rows
   * @param num_cols: Number of columns
   * @todo: Write as free function!!
   */
  static typename std::enable_if<std::is_arithmetic<T>::value, mathmatrix>::type
  identity(std::size_t const num_rows, std::size_t const num_cols);
  static mathmatrix fill_diag(uint_type const& rows, uint_type const& cols,
                              T const& fill);

  // base row and col proxy
  /*using base_type::col;
  using base_type::row;*/

  // element access from base class
  template<typename LeftIntegral, typename RightIntegral = std::size_t>
  typename std::enable_if<std::is_integral<LeftIntegral>::value,
    typename std::enable_if<std::is_integral<RightIntegral>::value, T>::type
  >::type &
  operator()(LeftIntegral const row, RightIntegral const col = RightIntegral()) {
    if (checkIfIndexOutOfBounds(row, col)) {
      throw std::runtime_error("Error: out of bound exception in operator() of scon_mathmatrix");
    }
    return base_type::operator()(row,col);
  }
  template<typename LeftIntegral, typename RightIntegral = std::size_t>
  typename std::enable_if<std::is_integral<LeftIntegral>::value,
    typename std::enable_if<std::is_integral<RightIntegral>::value, T>::type
    >::type const 
  operator()(LeftIntegral const row, RightIntegral const col = RightIntegral()) const {
    if (checkIfIndexOutOfBounds(row, col)) {
      throw std::runtime_error("error: out of bound exception in operator() const of scon_mathmatrix");
    }
    return base_type::operator()(row, col);
  }
  /*T operator[](std::size_t const row, std::size_t const col) {
    if (checkIfIndexOutOfBounds()) {
      std::runtime_error("Error: out of bound exception in operator[] of scon_mathmatrix");
    }
    return base_type::operator[];
  }*/

#ifndef CAST_USE_ARMADILLO
  void resize(uint_type const rows, uint_type const cols);
#else
  using base_type::resize;
#endif
  
  //If these guys are enabled the code won't compile. Don't ask me why...
  /*using base_type::operator=;
  using base_type::operator*=;
  using base_type::operator-=;
  using base_type::operator/=;
  using base_type::operator+=;*/
  /*! mathmatrix += operator
   *
   * @param in: Matrix to the right of the summation (this + in)
   * @return: Result of the addition of this + in
   */
  mathmatrix operator+(mathmatrix const& in) const;
  mathmatrix operator*(mathmatrix const& in) const;
  mathmatrix operator/(mathmatrix const& in) const;
  mathmatrix operator-(mathmatrix const& in) const;
  mathmatrix operator/(T const& in) const;
  mathmatrix operator*(T const& in) const;
  mathmatrix operator+(T const& in) const;
  mathmatrix operator-(T const& in) const;

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
   * @return -1 if determinant is negative or zero, +1 if determinant is greater
   * than zero.
   */
  int_type det_sign() const;

  /**
   * @brief Returns rank of the underlying matrix.
   * @return rank of the matrix.
   * @note A singular value decompostition is performed to determine the rank.
   * This may be quite costly.
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
   * @brief Append one matrix to another, will check if sizes match, appends at
   * the bottom end (rows are added)
   */
  void append_bottom(const mathmatrix& I_will_be_the_bottom_part);

  /**
   * @brief Append one matrix to another, will check if sizes match, appends on
   * the top (rows are added)
   */
  void append_top(const mathmatrix& I_will_be_the_top_part);

  /**
   * @brief Append one matrix to another, will check if sizes match, appends
   * left (columns are added)
   */
  void append_left(const mathmatrix& I_will_be_the_left_part);

  /**
   * @brief Append one matrix to another, will check if sizes match, appends
   * left (columns are added)
   */
  void append_right(const mathmatrix& I_will_be_the_right_part);

  /**
   * @brief Sheds the specified rows from the matrix
   */
  void shed_rows(long const first_in, long const last_in = 0);

  /**
   * @brief Sheds the specified columns from the matrix
   */
  void shed_cols(long const first_in, long const last_in = 0);

  /**
   * @brief Returns number of rows
   */
#ifndef CAST_USE_ARMADILLO
  std::size_t rows() const { return base_type::rows();}
#else
  std::size_t rows() const;
#endif

  /**
   * @brief Returns number of columns
   */
#ifndef CAST_USE_ARMADILLO
  std::size_t cols() const { return base_type::cols();}
#else
  std::size_t cols() const;
#endif

  mathmatrix row(std::size_t const idx) const;
  mathmatrix col(std::size_t const idx) const;

  void set_row(std::size_t const, mathmatrix const&);
  void set_col(std::size_t const, mathmatrix const&);

  /*! Performs Cholesky Decompostion on Matrix.
   *
   * @NOTE: Code via
   * https://rosettacode.org/wiki/Cholesky_decomposition#C 19.11.16
   * @param result: Upper triangular matrix as result of decompostion
   */
  void choleskyDecomposition(mathmatrix<T>& result) const;

  /**
   * @brief Returns whether mathmatrix-obj is quadratic
   */
  inline bool return_quadratic() const;

  /**
   * @brief Returns upper left submatrix.
   * If no second argument for the function call,
   * ie for columns_in is specified, then a quadratic submatrix with rows =
   * columns = rows_in is yieled.
   */
  mathmatrix upper_left_submatrix(uint_type rows_in,
                                  uint_type columns_in = 0) const;

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
  inline bool operator==(mathmatrix const& in) const;

  /**
   * @brief Performs singular value decomposition on *this and writes results
   * to the three resulting matrices U, s, V.
   */
  void singular_value_decomposition(mathmatrix& U_in, mathmatrix& s_in,
                                    mathmatrix& V_in) const;

  /**
   * @brief Performs singular value decomposition on *this and returns results
   * to the three resulting matrices U, s, V in a tuple in this sepcific sequence.
   * 
   * @note When CAST is compiled against Eigen, JacobiSVD is used because it has been
   * observed that the faster BDCSVD yields inaccurate results (in addition to
   * sometimes returning matrices full of NaN). The effect becomes visible in
   * geometry optimization with constrained internal coordinates.
   *
   * @see singular_value_decomposition(mathmatrix&, mathmatrix&, mathmatrix&)
   */
  std::tuple<mathmatrix, mathmatrix, mathmatrix> svd() const;

  T norm()const{
    auto sum{T()};
    for(auto i = 0u; i<rows(); ++i){
      for(auto j = 0u; j<cols();++j){
        sum += ((*this)(i,j))*((*this)(i,j));
      }
    }
    return std::sqrt(sum);
  }

  T rmsd()const{
    return norm()/std::sqrt(static_cast<T>(rows()*cols()));
  }

  T max()const;

  void reshape(long new_rows, long new_columns);

  /**
   * @brief Sorts the column of a matrix and returns a standard vector with the sorted values.
   *
   * The sorting order of the vector is determined by the passed compare fcunction which got to take two values of
   * type T and returns a boolean vaule. If no specific column is passed the first one is used.
   * @see sort_col(Comp, std::size_t const&)
   * @see sort_idx(std::size_t const&)
   * @see sort_col_asc(std::size_t const&)
   * @see sort_col_disc(std::size_t const&)
   */
  template <typename Comp>
  std::vector<T> sort_col_to_vec(Comp comp, std::size_t const ind = 0) const;

  /**
   * @brief Sorts the column of a matrix and returns a column matrix with the sorted values.
   *
   * The sorting order of the vector is determined by the passed compare fcunction which got to take two values of
   * type T and returns a boolean vaule. If no specific column is passed the first one is used.
   * @see sort_col_to_vec(Comp, std::size_t const&)
   * @see sort_idx(std::size_t const&)
   * @see sort_col_asc(std::size_t const&)
   * @see sort_col_disc(std::size_t const&)
   */
  template <typename Comp>
  mathmatrix sort_col(Comp comp, std::size_t const ind = 0) const;

  /**
   * @brief Sorts the column of a matrix and returns a column matrix in ascending order.
   *
   * If no specific column is passed the first one is used.
   * @see sort_col_to_vec(Comp, std::size_t const&)
   * @see sort_col(Comp, std::size_t const&)
   * @see sort_idx(std::size_t const&)
   * @see sort_col_disc(std::size_t const&)
   */
  mathmatrix sort_col_asc(std::size_t const ind = 0) const;

  /**
   * @brief Sorts the column of a matrix and returns a column matrix in discending order.
   *
   * If no specific column is passed the first one is used.
   * @see sort_col_to_vec(Comp, std::size_t const&)
   * @see sort_col(Comp, std::size_t const&)
   * @see sort_idx(std::size_t const&)
   * @see sort_col_asc(std::size_t const&)
   */
  mathmatrix sort_col_disc(std::size_t const ind = 0) const;

  /**
   * @brief Sorts the column of a matrix and returns the indexes in which the entries got to be sorted.
   *
   * If no specific column is passed the first one is used. The order is always ascending.
   * @see sort_col_to_vec(Comp, std::size_t const&)
   * @see sort_col(Comp, std::size_t const&)
   * @see sort_col_asc(std::size_t const&)
   * @see sort_col_disc(std::size_t const&)
   */
  std::vector<std::size_t> sort_idx(std::size_t const ind = 0) const;

  /**
   * @brief Finds all elements in a column compared with a passed function.
   *
   * The passed function got to take one agrument of type T and returns a boolean
   * value. An example would be to pass this lambda [](T const& a){return a<1.0;}
   * which would generate a vector containing all indices of elements less than 1.0.
   * If no column is specified the first one is used.
   * @return vector of integer which holds the indices of found elements
   */
  template <typename Comp>
  std::vector<std::size_t> find_idx(Comp comp,
                                    std::size_t const ind = 0) const;

//#ifndef CAST_USE_ARMADILLO
//private:
//  void removeRow(std::size_t const& rowToRemove);
//
//  void removeColumn(std::size_t const& colToRemove);
//
//public:
//#endif

  /**
   * @brief Builds a new matrix out of the columns and rows passed by two standard vectors.
   *
   * The vectors determining the cloumns and rows got to be the indices starting at 0.
   * There is no need to keep the passed vectors sorted.
   * @return new matrix with the desired columns and rows.
   */
  mathmatrix submat(std::vector<std::size_t> const& columns,
                    std::vector<std::size_t> const& rows) const;

   /**
    * @brief Builds a new matrix spanning between passed rows and columns.
    *
    * The first and third argument determine which rows are taken for the new matirx as do
    * the second and fourth argument for the columns.
    * @return new matrix with the desired columns and rows.
    */
  mathmatrix submat(std::size_t const rb, std::size_t const cb,
                    std::size_t const re, std::size_t const ce) const;

//#ifdef CAST_USE_ARMADILLO
//  private:
//    mathmatrix _submat(std::vector<std::size_t> const& columns,
//      std::vector<std::size_t> const& rows) {
//      return base_type::submat(columns, rows);
//    }
//    public:
//#endif

  /**
   * @brief Compute the pseudoinverse of the matrix.
   *
   * If Eigen is used the computation of the pseudoinverse is done via a singular value decomposition.
   * Otherwise armadillo's intern routine is used.
   * @return The pseudoinverse of the matrix.
   */
  mathmatrix pinv() const;
  mathmatrix lppinv() const;

  /**
   * @brief transforms the matrix into a column matrix.
   *
   * The elements are taken row by row and stored in a column matrix.
   * @return resulting column matrix
   */
  mathmatrix vectorise_col() const;

  mathmatrix vectorise_row() const;

  /**
   * @brief Replaces the specific elements with the passed values.
   *
   * The passed vector got to consist of indices which are to change to the value passed.
   * The indices start again by 0. If a matrix and not a column matrix is changed the indices
   * are evaluated column major.
   */
  void replace_idx_with(std::vector<std::size_t> const& idx, T const& val);

  /**
   * @brief Picks specific elements of the matrix and returns them as references.
   *
   * The passed vector got to consist of indices of the elements desired to be returned as references.
   * The indices start again by 0. If a matrix and not a column matrix is changed the indices
   * are evaluated column major.
   * @return a vector of references_wrapper referencing the desired elements.
   */
  std::vector<std::reference_wrapper<T>>
  elem(std::vector<std::size_t> const& idx);

  /**
   * @brief Builds a diagonal matrix.
   *
   * If the passed matrix is a column matrix the elements are set as diagonal elements, otherwise
   * the diabonal elements are picked and the others are set to zero.
   * @return A diagonal matrix.
   */
  mathmatrix diagmat() const;

  /**
   * @brief Performs a eigenvalue decomposition.
   *
   * The resulting eigenvalues and eigenvectors are calculated and returned. The object itself
   * keeps unchanged.
   * @note Matrix is assumed to be symmetric.
   * @return the eigenvalues and eigenvectors as a std::pair.
   * @see diag()
   */
  std::pair<mathmatrix, mathmatrix>
      eigensym(bool const & sort = false);

  /**
   * @brief uses eigensym to diagonalize the matrix
   *
   * The object is changed to the eigenvalues as diagonal matrix.
   * @return the eigenvalues and eigenvectors as a std::pair.
   * @see eigensym()
   */
  std::pair<mathmatrix, mathmatrix> diag();

  /**
   * @brief returns the transpose of the matrix
   */
  mathmatrix t() const;

  /**
   * @brief tranforms the specific column into a vector.
   *
   * The passed index of the column is transformed into a standard vector.
   * If no value is passed the first column is used.
   * @return standard vector with the specific column.
   */
  std::vector<T> col_to_std_vector(std::size_t const iter = 0) const;

  /**
  * @brief tranforms the specific row into a vector.
  *
  * The passed index of the row is transformed into a standard vector.
  * If no value is passed the first row is used.
  * @return standard vector with the specific row.
  */
  std::vector<T> row_to_std_vector(std::size_t const iter = 0) const;

  /**
   * @brief Returns mathmatrix-obj as std vector of vector of T.
   * USE THIS TO DEBUG, not in production code.
   * This might be useful as the VS debugger cannot visualize content of arma
   * arrays.
   *
   * @return Matrix in std::vector<vector<float-type> > form.
   */
  std::vector<std::vector<T>> to_std_vector() const;

  /**
   * Updates the internal std::vector<std::vector<float_type> >
   * array_debugview_internal array ONLY IF PREPROCESSOR FLAG "DEBUGVIEW" IS SET
   */
  void update_debugview(void) const;

  /**
   * @brief returns true if a matrix is a column matrix
   */
  bool is_vec() const;

#ifndef CAST_USE_ARMADILLO
  using Quaternion = Eigen::Quaternion<T>;
#else
  class Quaternion {};
#endif
  private:
    bool checkIfIndexOutOfBounds(std::size_t const row, std::size_t const col) const {
      return row > rows() || col > cols();
    }
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

template <typename T>
mathmatrix<T>::mathmatrix(
    std::initializer_list<std::initializer_list<T>> const& ini) {
  if (ini.size() == 0) {
    *this = mathmatrix();
    return;
  }

  auto const& rows = ini.size();
  auto const& cols = (*ini.begin()).size();

  *this = mathmatrix(rows, cols);

  for (auto i = 0u; i < rows; ++i) {
    for (auto j = 0u; j < cols; ++j) {
      (*this)(i, j) = *((*(ini.begin() + i)).begin() + j);
    }
  }
}
#ifndef CAST_USE_ARMADILLO
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
mathmatrix<T> mathmatrix<T>::col_from_vec(std::vector<T> const& col) {
  auto const& size = col.size();

  mathmatrix ret(size, 1);

  for (auto i = 0u; i < size; ++i) {
    ret(i, 0) = col.at(i);
  }
  return ret;
}

template <typename T>
mathmatrix<T> mathmatrix<T>::row_from_vec(std::vector<T> const& row) {
  auto const& size = row.size();

  mathmatrix ret(1, size);

  for (auto i = 0u; i < size; ++i) {
    ret(0, i) = row.at(i);
  }
  return ret;
}

template <typename T>
mathmatrix<T> mathmatrix<T>::zero(uint_type rows, uint_type cols) {
  return mathmatrix(rows, cols, T());
}

template <typename T>
mathmatrix<T> mathmatrix<T>::fill_diag(uint_type const& rows,
                                       uint_type const& cols, T const& fill) {
  mathmatrix ret(rows, cols, T());
  for (uint_type i = 0; i < rows; ++i) {
    ret(i, i) = fill;
  }
  return ret;
}

template <typename T>
mathmatrix<T> mathmatrix<T>::operator+(mathmatrix const& in) const {
  if (cols() != in.cols() || rows() != in.rows()) {
    throw std::runtime_error(
        "Either the rows or the columns of the added matrices are not equal!");
  }
#ifndef CAST_USE_ARMADILLO
  return base_type::operator+(in);
#else
  return arma::operator+(static_cast<base_type>(*this),
                         static_cast<base_type>(in));
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::operator*(mathmatrix const& in) const {
  if (cols() != in.rows()) {
    throw std::runtime_error(
        "The number of colums of the left hand matrix got to be equal to the "
        "rows of the right hand matrix!");
  }

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
  if (cols() != in.cols() || rows() != in.rows()) {
    throw std::runtime_error("Either the rows or the columns of the "
                             "substracted matrices are not equal!");
  }
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
mathmatrix<T> mathmatrix<T>::operator+(T const& in) const {
#ifndef CAST_USE_ARMADILLO
  return this->array() + in;
#else
  return arma::operator+(static_cast<base_type>(*this), in);
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::operator-(T const& in) const {
#ifndef CAST_USE_ARMADILLO
  return this->array() - in;
#else
  return arma::operator-(static_cast<base_type>(*this), in);
#endif
}

template <typename T>
typename std::enable_if<std::is_arithmetic<T>::value, mathmatrix<T>>::type
mathmatrix<T>::identity(std::size_t const num_rows,
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

#ifndef CAST_USE_ARMADILLO
  auto old_mat = *this;
  resize(rows() + I_will_be_the_bottom_part.rows(), cols());
  *this << old_mat, I_will_be_the_bottom_part;
#else

  // Old size needs to be kept
  auto const holder = rows();

  resize(holder + I_will_be_the_bottom_part.rows(), cols());

  // Add "in" to newly created space.
  for (auto i = 0u; i < I_will_be_the_bottom_part.rows(); ++i) {
    for (auto j = 0u; j < cols(); ++j) {
      (*this)(i + holder, j) = I_will_be_the_bottom_part(i, j);
    }
  }
#endif
}

template <typename T>
void mathmatrix<T>::append_top(const mathmatrix& I_will_be_the_top_part) {
  if (this->cols() != I_will_be_the_top_part.cols()) {
    throw std::runtime_error("Wrong Matrix size in mathmatrix:append()");
  }

#ifndef CAST_USE_ARMADILLO
  auto old_mat = *this;
  resize(rows() + I_will_be_the_top_part.rows(), cols());
  *this << I_will_be_the_top_part, old_mat;
#else

  auto const otherRows = I_will_be_the_top_part.rows();
  auto const old_mat = *this;

  this->resize(rows() + otherRows, cols());

  for (auto i = 0u; i < cols(); ++i) {
    // Move the entries in the parent matrix rightward
    // We count right so that we dont overwrite
    for (auto j = 0u; j < old_mat.rows(); ++j) {
      (*this)(otherRows + j, i) = old_mat(j, i);
    }
    // Add "I_will_be_the_left_part" to now absolete top space of the parent
    // matrix (this).
    for (auto j = 0u; j < otherRows; ++j) {
      (*this)(j, i) = I_will_be_the_top_part(j, i);
    }
  }
#endif
}

template <typename T>
void mathmatrix<T>::append_left(const mathmatrix& I_will_be_the_left_part) {
  if (this->rows() != I_will_be_the_left_part.rows()) {
    throw std::runtime_error("Wrong Matrix size in mathmatrix:append()");
  }

#ifndef CAST_USE_ARMADILLO
  auto old_mat = *this;
  resize(rows(), cols() + I_will_be_the_left_part.cols());
  *this << I_will_be_the_left_part, old_mat;
#else

  auto const otherCols = I_will_be_the_left_part.cols();
  auto const old_mat = *this;

  this->resize(rows(), cols() + otherCols);

  for (auto i = 0u; i < rows(); ++i) {
    // Move the entries in the parent matrix rightward
    // We count right so that we dont overwrite
    for (auto j = 0u; j < old_mat.cols(); ++j) {
      (*this)(i, otherCols + j) = old_mat(i, j);
    }
    // Add "I_will_be_the_left_part" to now absolete top space of the parent
    // matrix (this).
    for (auto j = 0u; j < otherCols; ++j) {
      (*this)(i, j) = I_will_be_the_left_part(i, j);
    }
  }

#endif
}

template <typename T>
void mathmatrix<T>::append_right(const mathmatrix& I_will_be_the_right_part) {

  if (this->rows() != I_will_be_the_right_part.rows()) {
    throw std::runtime_error("Wrong Matrix size in mathmatrix:append()");
  }
#ifndef CAST_USE_ARMADILLO
  auto old_mat = *this;
  resize(rows(), cols() + I_will_be_the_right_part.cols());
  *this << old_mat, I_will_be_the_right_part;
#else
  // Old size needs to be kept
  auto const holder = cols();

  resize(rows(), holder + I_will_be_the_right_part.cols());

  // Add "in" to newly created space.
  for (auto i = 0u; i < I_will_be_the_right_part.cols(); ++i) {
    for (auto j = 0u; j < rows(); ++j) {
      (*this)(j, i + holder) = I_will_be_the_right_part(j, i);
    }
  }
#endif
}

template <typename T>
void mathmatrix<T>::shed_rows(long const first_in, long const last_in) {

  auto const last_in_ = last_in == 0 ? first_in : last_in;

  if (first_in < 0 || first_in > last_in_ || last_in >= static_cast<long>(this->rows())) {
    throw std::runtime_error("Index Out of Bounds in mathmatrix:shed_rows()");
  }


  mathmatrix newOne(this->rows() - (last_in_ - first_in + 1u), this->cols());
  for (auto i = 0u; i < first_in; ++i) {
    for (auto j = 0u; j < this->cols(); ++j) {
      newOne(i, j) = (*this)(i, j);
    }
  }
  for (auto i = first_in; i < static_cast<long>(this->rows()) - last_in_ - 1u + first_in; ++i) {
    for (auto j = 0u; j < this->cols(); ++j) {
      newOne(i, j) = (*this)(i + (last_in_ - first_in) + 1u, j);
    }
  }
  this->swap(newOne);
}

template <typename T>
void mathmatrix<T>::shed_cols(long const first_in, long const last_in) {

  auto const last_in_ = last_in == 0 ? first_in : last_in;

  if (first_in < 0 || first_in > last_in_ || last_in >= static_cast<long>(this->cols())) {
    throw std::runtime_error("Index Out of Bounds in mathmatrix:shed_cols()");
  }

  mathmatrix newOne(this->rows(), this->cols() - (last_in_ - first_in + 1u));
  for (auto j = 0u; j < this->rows(); ++j) {
    for (auto i = 0u; i < first_in; ++i) {
      newOne(j, i) = (*this)(j, i);
    }
  }
  for (auto j = 0u; j < this->rows(); ++j) {
    for (auto i = first_in; i < static_cast<long>(this->cols()) - last_in_ - 1u + first_in; ++i) {
      newOne(j, i) = (*this)(j, i + (last_in_ - first_in) + 1u);
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
inline mathmatrix<T> mathmatrix<T>::col(std::size_t const idx) const {
  if (idx > cols() - 1u) {
    throw std::runtime_error("The boundaries for the rows are exceeded. See function col(std::size_t const) in the mathmatrix class");
  }
#ifndef CAST_USE_ARMADILLO
  return base_type::col(idx);
#else
  auto const& nc = cols();
  mathmatrix ret(nc, 1);

  for (auto i = 0; i < nc; ++i) {
    ret(i, 0) = operator()(i, idx);
  }
  return ret;
#endif
}

template<typename T>
inline void mathmatrix<T>::set_row(std::size_t const nrow, mathmatrix const & other)
{
  if (other.cols() != cols() || other.rows() != 1) {
    throw std::runtime_error("By setting the row the sizes for both rows are different!");
  }
  for (auto i = 0u; i < cols(); ++i) {
    this->operator()(nrow, i) = other(0, i);
  }
}

template<typename T>
inline void mathmatrix<T>::set_col(std::size_t const ncol, mathmatrix const & other)
{
  if (other.rows() != rows() || other.cols() != 1) {
    throw std::runtime_error("By setting the col the sizes for both cols are different!");
  }
  for (auto i = 0u; i < rows(); ++i) {
    this->operator()(i, ncol) = other(i, 0);
  }
}

template <typename T>
inline mathmatrix<T> mathmatrix<T>::row(std::size_t const idx) const {
  if (idx > rows() - 1u) {
    throw std::runtime_error("The boundaries for the rows are exceeded. See function row(std::size_t const) in the mathmatrix class");
  }
#ifndef CAST_USE_ARMADILLO
  return base_type::row(idx);
#else
  auto const& nr = rows();
  mathmatrix ret(1, nr);

  for (auto i = 0; i < nr; ++i) {
    ret(0, i) = operator()(idx, i);
  }
  return ret;
#endif
}

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
mathmatrix<T>::submat(std::size_t const rb, std::size_t const cb,
                      std::size_t const re, std::size_t const ce) const {
#ifndef CAST_USE_ARMADILLO
  //Strange behaviour of Eigen if the submatrix spans from index 0. To circumvent:
  auto re_ = rb == 0 ? re + 1 : re;
  auto ce_ = cb == 0 ? ce + 1 : ce;

  return this->block(rb, cb, re_, ce_);
#else
  return base_type::submat(rb, cb, re, ce);
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::upper_left_submatrix(uint_type rows_in,
                                                  uint_type columns_in) const {
  return submat(0, 0, rows_in-1, (columns_in == 0 ? rows_in : columns_in)-1);
}

template <typename T>
bool mathmatrix<T>::operator==(mathmatrix<T> const& in) const {
#ifndef CAST_USE_ARMADILLO
  return this->isApprox(in, matCompTol());
#else
  return (arma::approx_equal(static_cast<base_type>(*this),
                             static_cast<base_type>(in), "absdiff",
                             matCompTol()));
#endif
}

//template<typename T>
//bool operator==(mathmatrix<T> const& lhs, mathmatrix<T> const& rhs) {
//  return lhs.operator==(rhs);
//}

template <typename T>
void mathmatrix<T>::singular_value_decomposition(mathmatrix& U_in,
                                                 mathmatrix& s_in,
                                                 mathmatrix& V_in) const {
  std::tie(U_in, s_in, V_in) = svd();
}

template <typename T>
inline std::tuple<mathmatrix<T>, mathmatrix<T>, mathmatrix<T>>
mathmatrix<T>::svd() const {

#ifndef CAST_USE_ARMADILLO
  mathmatrix s, U, V;
  Eigen::JacobiSVD<base_type> jacobiSvd(static_cast<base_type>(*this), Eigen::ComputeFullU | Eigen::ComputeFullV);
  s = static_cast<base_type>(jacobiSvd.singularValues());
  U = jacobiSvd.matrixU();
  V = jacobiSvd.matrixV();
  
#else
  arma::Col<T> s_arma;
  mathmatrix U, V;
  if (!svd_econ(U, s_arma, V, *this))
    throw std::runtime_error("Error in armadillo SVD: failed.");
  mathmatrix s(s_arma);
#endif
  return std::make_tuple(U, s, V);
}

template<typename T>
T mathmatrix<T>::max()const{
  #ifndef CAST_USE_ARMADILLO
  return this->maxCoeff();
  #else
  return base_type::max();
  #endif
}

template <typename T>
template <typename Comp>
std::vector<T> mathmatrix<T>::sort_col_to_vec(Comp comp,
                                              std::size_t const ind) const {
  auto ret = col_to_std_vector(ind);

  std::sort(ret.begin(), ret.end(), comp);

  return ret;
}

template <typename T>
template <typename Comp>
mathmatrix<T> mathmatrix<T>::sort_col(Comp comp, std::size_t const ind) const {

  auto sorted_vec = sort_col_to_vec(comp, ind);

  auto const& size = sorted_vec.size();

  mathmatrix ret(size, 1);

  for (auto i = 0u; i < size; ++i) {
    ret(i, 0) = sorted_vec.at(i);
  }

  return ret;
}

template <typename T>
mathmatrix<T> mathmatrix<T>::sort_col_asc(std::size_t const ind) const {
  return sort_col([](auto const& a, auto const& b) { return a < b; }, ind);
}

template <typename T>
inline mathmatrix<T>
mathmatrix<T>::sort_col_disc(std::size_t const ind) const {
  return sort_col([](auto const& a, auto const& b) { return a > b; }, ind);
}

template <typename T>
std::vector<std::size_t> mathmatrix<T>::sort_idx(std::size_t const ind) const {
  auto val_vec = col_to_std_vector(ind);
  std::vector<std::size_t> ret(val_vec.size());
  std::iota(ret.begin(), ret.end(), 0);

  std::sort(ret.begin(), ret.end(),
            [&val_vec](std::size_t const i, std::size_t const j) {
              return val_vec.at(i) < val_vec.at(j);
            });

  return ret;
}

template <typename T>
template <typename Comp>
std::vector<std::size_t> mathmatrix<T>::find_idx(Comp comp,
                                                 std::size_t const ind) const {
  auto val_vec = col_to_std_vector(ind);
  std::vector<std::size_t> ret(val_vec.size());
  std::iota(ret.begin(), ret.end(), 0);

  ret.erase(
      std::remove_if(ret.begin(), ret.end(),
                     [&](std::size_t const i) { return !comp(val_vec[i]); }),
      ret.end());

  return ret;
}

//#ifndef CAST_USE_ARMADILLO
//template <typename T>
//void mathmatrix<T>::removeRow(std::size_t const& rowToRemove) {
//  auto numRows = rows() - 1;
//  auto numCols = cols();
//
//  if (rowToRemove < numRows) {
//    base_type::block(rowToRemove, 0, numRows - rowToRemove, numCols) =
//        base_type::block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
//  }
//
//  base_type::conservativeResize(numRows, numCols);
//}
//
//template <typename T>
//void mathmatrix<T>::removeColumn(std::size_t const& colToRemove) {
//  unsigned int numRows = rows();
//  unsigned int numCols = cols() - 1;
//
//  if (colToRemove < numCols) {
//    base_type::block(0, colToRemove, numRows, numCols - colToRemove) =
//        base_type::block(0, colToRemove + 1, numRows, numCols - colToRemove);
//  }
//
//  base_type::conservativeResize(numRows, numCols);
//}
//#endif

template <typename T>
mathmatrix<T>
mathmatrix<T>::submat(std::vector<std::size_t> const& rows,
                      std::vector<std::size_t> const& columns) const {

//#ifndef CAST_USE_ARMADILLO

  auto reverse_vec = [&](std::vector<std::size_t> const& vec, std::size_t const size) {
    std::vector<std::size_t> ret(size);
    std::iota(ret.rbegin(), ret.rend(), 0);
    ret.erase(std::remove_if(ret.begin(), ret.end(), [&](std::size_t const x) {
      return std::find(vec.begin(), vec.end(),x) != vec.end();
    }), ret.end());
    return ret;
  };

  mathmatrix ret(*this);

  std::vector<std::size_t> del_rows = reverse_vec(rows, this->rows());
  std::vector<std::size_t> del_cols = reverse_vec(columns, cols());

  for (auto const& i : del_rows) {
    ret.shed_rows(i);
  }
  for (auto const& i : del_cols) {
    ret.shed_cols(i);
  }
  return ret;
//#else
//  auto ret = *this;
//  return ret._submat(rows, columns);
//#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::pinv() const {
#ifndef CAST_USE_ARMADILLO
  mathmatrix U, V, s;
  singular_value_decomposition(U, s, V);
  mathmatrix s_inv = zero(rows(), cols());

  for (auto i = 0u; i < s.rows(); ++i) {
    s_inv(i, i) = std::fabs(s(i)) > close_to_zero_tol ? 1. / s(i) : 0.0;
  }

  return V * s_inv * U.t();
#else
  return arma::pinv(static_cast<base_type>(*this));
#endif
}

template <typename T>
mathmatrix<T> mathmatrix<T>::lppinv() const {
#ifndef CAST_USE_ARMADILLO
  mathmatrix U, V, s;
  singular_value_decomposition(U, s, V);
  mathmatrix s_inv = zero(rows(), cols());

  for (auto i = 0u; i < s.rows(); ++i) {
    s_inv(i, i) = std::fabs(s(i)) > close_to_zero_tol ? 1. / s(i) : 0.0;
  }

  return V.t() * s_inv * U.t();
#else
  return arma::pinv(static_cast<base_type>(*this));
#endif
}

//template <typename T>
//mathmatrix<T> mathmatrix<T>::vectorise(std::size_t const& i) const {
//#ifndef CAST_USE_ARMADILLO
//  return base_type(Eigen::Map<Eigen::RowVectorXd>(this->data(), this->size()));
//#else
//  return arma::vectorise(*this, i);
//#endif
//}

template<typename T>
mathmatrix<T> mathmatrix<T>::vectorise_col() const {

  std::vector<T> ret;
  ret.reserve(rows()*cols());
  for (auto i = 0; i < rows(); ++i) {
    auto tmp_col = col_to_std_vector(i);
    ret.insert(ret.end(), tmp_col.begin(), tmp_col.end());
  }

  return col_from_vec(ret);
}

template<typename T>
mathmatrix<T> mathmatrix<T>::vectorise_row() const {

  std::vector<T> ret;
  ret.reserve(rows()*cols());
  for (auto i = 0u; i < rows(); ++i) {
    auto tmp_row = row_to_std_vector(i);
    ret.insert(ret.end(), tmp_row.begin(), tmp_row.end());
  }

  return col_from_vec(ret);
}

template <typename T>
void mathmatrix<T>::replace_idx_with(std::vector<std::size_t> const& idx,
                                     T const& val) {
#ifndef CAST_USE_ARMADILLO
  auto dat = base_type::data();
#else
  auto dat = base_type::memptr();
#endif
  for (auto const& i : idx) {
    *(dat + i) = val;
  }
}

template <typename T>
inline std::vector<std::reference_wrapper<T>>
mathmatrix<T>::elem(std::vector<std::size_t> const& idx) {
  std::vector<std::reference_wrapper<T>> ret;
#ifndef CAST_USE_ARMADILLO
  auto dat = base_type::data();
#else
  auto dat = base_type::memptr();
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
  return arma::diagmat(static_cast<base_type>(*this));
#endif
}

  template <typename T>
  void close_to_zero_to_zero(mathmatrix<T> & mat){
    for(auto r{0u}; r<mat.rows(); ++r){
      for(auto c{0u};c<mat.cols(); ++c){
        if(mat(r,c)<mathmatrix<T>::close_to_zero_tol){
          mat(r,c) = T();
        }
      }
    }
  }


template <typename T>
std::pair<mathmatrix<T>, mathmatrix<T>>
mathmatrix<T>::eigensym(bool const & sort) {
#ifndef CAST_USE_ARMADILLO

  //close_to_zero_to_zero(*this);

  Eigen::SelfAdjointEigenSolver<base_type> es(static_cast<base_type>(*this));
  mathmatrix eigenval = es.eigenvalues().real();
  mathmatrix eigenvec = es.eigenvectors().real();

  if (sort) {
    auto indices = eigenval.sort_idx();
    mathmatrix new_eigenvec(eigenvec.rows(), eigenvec.cols());
    mathmatrix new_eigenval(eigenval.rows(), eigenval.cols());
    for (auto i = 0u; i < indices.size(); ++i) {
      auto index = indices[i];
      new_eigenvec.set_col(i, eigenvec.col(index));
      new_eigenval(i,0) = eigenval(index,0);
    }
    return std::make_pair(new_eigenval, new_eigenvec);
  }

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
std::vector<T> mathmatrix<T>::col_to_std_vector(std::size_t const iter) const {
  if (cols() <= iter) {
    throw std::runtime_error("The required col is out of bounds!");
  }
#ifndef CAST_USE_ARMADILLO
  auto const& dat = base_type::data();
#else
  auto const& dat = base_type::memptr();
#endif
  auto const& r = rows();
  return std::vector<T>(dat + r * iter, dat + r * (iter + 1));
//#else
//  return arma::conv_to<std::vector<T>>::from(static_cast<arma::Col>(col(iter)));
//#endif
}

template <typename T>
std::vector<T> mathmatrix<T>::row_to_std_vector(std::size_t const iter) const {
  if (rows() <= iter) {
    throw std::runtime_error("The required row is out of bounds!");
  }
  std::vector<T> ret;
  ret.reserve(cols());
#ifndef CAST_USE_ARMADILLO
  auto const& dat = base_type::data();
#else
  auto const& dat = base_type::memptr();
#endif
  auto const& r = rows();
  for (auto i = 0u; i < cols(); ++i) {
    ret.emplace_back(*(dat + (r*i) + iter));
  }
  return ret;
  //#else
  //  return arma::conv_to<std::vector<T>>::from(static_cast<arma::Col>(col(iter)));
  //#endif
}

template <typename T>
std::vector<std::vector<T>> mathmatrix<T>::to_std_vector() const {
  std::vector<std::vector<T>> ret;
  ret.reserve(rows());
  for (auto i = 0u; i < rows(); ++i) {
    ret.emplace_back(row_to_std_vector(i));
  }
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
inline bool mathmatrix<T>::is_vec() const {
  return cols() == 1;
}

template<typename T>
void mathmatrix<T>::reshape(long new_rows, long new_cols){
  if(new_rows == -1 && new_cols == -1){
    throw std::runtime_error("Error in rehsape: You can't pass -1 for both columns and rows!");
  }
  else if(new_rows == -1){
    new_rows = rows()*cols()/new_cols;
  }
  else if(new_cols == -1){
    new_cols = rows()*cols()/new_rows;
  }
#ifndef CAST_USE_ARMADILLO
using RM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
RM rm = *this;
  *this = Eigen::Map<RM>(rm.data(), new_rows, new_cols);
#else
  this->reshape(new_rows, new_cols);
#endif
}

} // namespace scon
// END HEADER
