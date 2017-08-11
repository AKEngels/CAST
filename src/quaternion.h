#ifndef cast_quaternion_h_guard
#define cast_quaternion_h_guard

#pragma once

#include "coords.h"
#include "coords_rep.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>

/*!
 *  \addtogroup ic_util
 *  @{
 */
namespace ic_util {

/*!
\brief Quaternion class that only supports the minimal feature set required to
be operational.
\details Each quaternion is essentially an array with multiplication and an
output operator attached.
\tparam T Arithmetic type. Enforced by using type traits.
*/
template <typename T,
          typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
class Quaternion {
public:
  Quaternion(const std::vector<T>& vec) : q_{ create_quaternion(vec) } {}
  // Copy constructor.
  Quaternion(const Quaternion& quat) : q_{ quat.q_ } {}

  std::array<T, 4u> q_;

  /*!
  \brief Function for creating a quaternion from a 4-dimensional vector. Called
  during instantiation of the Quaternion class.
  \param v 4-dimensional vector used for construction.
  \return 4-dimensional array.
  */
  std::array<T, 4u> create_quaternion(const std::vector<T>& v) {
    if (v.size() == 4) {
      std::array<T, 4u> arr;
      // std::copy_n copies exactly 4 values from the range beginning at
      // v.begin() to the range beginning at arr.begin().
      std::copy_n(std::make_move_iterator(v.begin()), 4, arr.begin());
      return arr;
    } else {
      throw("Vector has not 4 elements.");
    }
  }

  /*!
  \brief Overloaded multiplication operator.
  \param Factor used for multiplication.
  \return Multiplied quaternion.
  */
  Quaternion operator*(T& num) {
    (*this).q_.at(0) *= num;
    (*this).q_.at(1) *= num;
    (*this).q_.at(2) *= num;
    (*this).q_.at(3) *= num;
    return *this;
  }

  /*!
  \brief Overloaded multiplication operator.
  \param Factor used for multiplication.
  \return Multiplied quaternion.
  */
  Quaternion operator*(int num) {
    (*this).q_.at(0) *= num;
    (*this).q_.at(1) *= num;
    (*this).q_.at(2) *= num;
    (*this).q_.at(3) *= num;
    return *this;
  }

  /*!
  \brief Templated and overloaded output operator.
  */
  template <typename T>
  friend std::ostream& operator<<(std::ostream& os, const Quaternion<T>& quat);
};

/*!
\brief Templated and overloaded output operator.
\tparam Arithmetic type. Same type as used for the instantiation of the
Quaternion object.
\param os std::ostream object for output.
\param quat Quaternion that shall be streamed to output.
\return Quaternion.
*/
template <typename T>
std::ostream& operator<<(std::ostream& os, const Quaternion<T>& quat) {
  return os << "Quaternion: " << quat.q_.at(0) << ", " << quat.q_.at(1) << ", "
            << quat.q_.at(2) << ", " << quat.q_.at(3);
}
}

/*! @} End of ic_util group*/

#endif // cast_quaternion_h_guard