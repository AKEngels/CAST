#ifndef cast_ic_util_h_guard
#define cast_ic_util_h_guard

#include<cstring>
#include<type_traits>
#include<vector>
#include<array>

#include"scon_mathmatrix.h"
#include "ic_atom.h"

/*!
 *  \addtogroup ic_util
 *  @{
 */
namespace ic_util{

  using coords::float_type;

  /*template <typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
  typename std::enable_if<std::is_arithmetic<T>::value, std::vector<std::vector<scon::mathmatrix<T>> >>::type
  F_matrix_derivs(ContainerType<CoordType<T>, ContainerArgs...> const& new_xyz);*/
  
  template<template <typename, typename, typename ...> class Map, typename Key, typename Value, typename DefValue, typename ... Args>
  inline Value getValueByKeyOrDefault(Map<Key, Value, Args...> const& map, Key const& key, DefValue const& defaultValue){
      typename Map<Key, Value, Args ...>::const_iterator it = map.find(key);
      if(it == map.end()){
          return defaultValue;
      }
      return it->second;
  }

  template<template <typename, typename...> class Vec, typename VecType, typename ... VecArgs>
  inline auto get_mean(Vec<VecType, VecArgs...> const & vec) {
    auto mean = std::accumulate(vec.begin(), vec.end(), VecType(), std::plus<VecType>());
    mean /= static_cast<float_type> (vec.size());
    return mean;
  }

  /*!
  \brief Converts a std::array to a std::vector.
  \tparam U Arithmetic type. Enforced by type traits.
  \tparam N Size of the array.
  \param arr Input array.
  \return std::vector.
  */
  template <template<typename, typename ...> class Conatiner, typename VectorType, typename ... VectorArguments> 
  inline typename std::enable_if<std::is_arithmetic<VectorType>::value, std::vector<VectorType>>::type
    arr_to_vec(const Conatiner<VectorType, VectorArguments ...>& container) {
    return std::vector<VectorType>(container.begin(), container.end());
  }
  template <typename U, std::size_t N>
  inline typename std::enable_if<std::is_arithmetic<U>::value, std::vector<U>>::type
  arr_to_vec(const std::array<U, N>& arr) {
    std::vector<U> res(arr.begin(), arr.end());
    return res;
  }

  template <typename T, template<typename> class CoordType, template <typename, typename ...> class ContainerType, typename ... ContainerArgs>
  inline typename std::enable_if<std::is_arithmetic<T>::value, scon::mathmatrix<T>>::type
  Rep3D_to_Mat(const ContainerType<CoordType<T>, ContainerArgs ...>& rep) {
    using Mat = scon::mathmatrix<T>;

    Mat A(rep.size(), 3);
    for (std::size_t i = 0; i != rep.size(); ++i) {
      A(i, 0) = rep.at(i).x();
      A(i, 1) = rep.at(i).y();
      A(i, 2) = rep.at(i).z();
    }
    return A;
  }


  // radius of gyration from Rep3D
  // coords have to be in Bohr
  template <typename T, template<typename> class CoordType, template<typename, typename ...> class ContainerType, typename ... ContainerArgs>
  inline typename std::enable_if<std::is_arithmetic<T>::value, T>::type
  rad_gyr(ContainerType<CoordType<T>, ContainerArgs...> const& rep) {
    auto mean = get_mean(rep);
    auto diff = rep - mean;
    auto sum{ 0.0 };
    for (auto& i : diff) {
      auto temp = std::pow(i.x(), 2) + std::pow(i.y(), 2) + std::pow(i.z(), 2);
      sum += temp;
    }
    sum /= rep.size();
    return std::sqrt(sum);
  }

  /*!
  \brief Normalizes a vector given as Cartesian_Point.
  \param a Vector, which shall be normalized.
  \return Normalized vector.
  */
  inline coords::Cartesian_Point normalize(const coords::Cartesian_Point& a) {
    return a / scon::geometric_length(a);
  }

  inline std::vector<std::pair<std::size_t, std::size_t>>
  bonds(std::vector<std::string> const& elem_vec,
        coords::Representation_3D const& cp_vec) {

    using scon::len;

    std::vector<std::pair<std::size_t, std::size_t>> connectedAtoms;

    for (auto i = 0u; i < elem_vec.size(); ++i) {
      for (auto j = i + 1u; j < elem_vec.size(); ++j) {
        auto maximumDistanceToBeBonded = 1.2 * (ic_atom::element_radius(elem_vec.at(i)) + ic_atom::element_radius(elem_vec.at(j)));
        if (scon::len(cp_vec.at(i) - cp_vec.at(j)) < maximumDistanceToBeBonded) {
          connectedAtoms.emplace_back(i+1u,j+1u);
        }
      }
    }

    return connectedAtoms;

  }

  /*!
  \brief Creates all possible 3-permutations from a vector containing 3 integers.
  \param vec Vector containing 3 integers.
  \return Vector of all created permutations, where each permutation is itself a
  3-dimensional vector.
  */
  inline std::vector<std::vector<std::size_t>>
  permutation_from_vec(std::vector<std::size_t>& vec) {
    std::vector<std::vector<std::size_t>> result;
    do {
      std::cout << vec.at(0) << " " << vec.at(1) << " " << vec.at(2) << "\n";
      result.emplace_back(std::vector<std::size_t>{ vec.at(0), vec.at(1), vec.at(2) });
    } while (std::next_permutation(vec.begin(), vec.end()));
    return result;
  }

  /*!
  \brief Performs standard container flattening.
  \tparam ContainerIt Type of container that is to be flattened.
  \tparam Result Type of resulting container; might be an iterator.
  \param st begin() iterator.
  \param fi end() iterator.
  \param res Usually an iterator to the flattened container.
  */
  template <typename ContainerIt, typename Result>
  void concatenate(ContainerIt st, ContainerIt fi, Result res) {
    while (st != fi) {
      res = std::move(st->begin(), st->end(), res);
      ++st;
    }
  }

  /*!
  \brief Flattens a std::vector of scon::c3<> vectors.
  \tparam T Arithmetic type. Enforced by type traits.
  \param vec std::vector that is to be flattened.
  \return std::vector.
  */
  template <typename T, template<typename> class CoordType, template<typename, typename...> class ContainerType, typename ... ContainerArgs>
  inline typename std::enable_if<std::is_arithmetic<T>::value, std::vector<T>>::type
  flatten_c3_vec(const ContainerType<CoordType<T>, ContainerArgs ...>& vec) {
    std::vector<T> result;
    result.reserve(3 * vec.size());
    for (auto& i : vec) {
      result.emplace_back(i.x());
      result.emplace_back(i.y());
      result.emplace_back(i.z());
    }
    return result;
  }

  template<typename Mat>
  coords::Representation_3D mat_to_rep3D(Mat&& mat){
    coords::Representation_3D result;
    for(auto i = 0u; i<mat.rows(); i+=3){
      result.emplace_back(coords::Cartesian_Point(
        mat(i, 0),mat(i + 1, 0),mat(i + 2, 0)
      ));
    }
    return result;
  }
}

/*! @} End of ic_util group*/

#endif // cast_ic_util_h_guard
