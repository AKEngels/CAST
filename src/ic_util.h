#ifndef cast_ic_util_h_guard
#define cast_ic_util_h_guard

#include<cstring>
#include<type_traits>
#include<vector>
#include<array>

#include"scon_mathmatrix.h"

/*!
 *  \addtogroup ic_util
 *  @{
 */
namespace ic_util{

  using float_type = coords::float_type;
  /*!
  \brief Converts a std::array to a std::vector.
  \tparam U Arithmetic type. Enforced by type traits.
  \tparam N Size of the array.
  \param arr Input array.
  \return std::vector.
  */
  template <typename U, std::size_t N>
  inline typename std::enable_if<std::is_arithmetic<U>::value, std::vector<U>>::type
  arr_to_vec(const std::array<U, N>& arr) {
    std::vector<U> res(arr.begin(), arr.end());
    return res;
  }

  template <typename T, template<typename> class CoordType, template <typename, typename ...> class ContainerType, typename ... ContainerArgs>
  inline typename std::enable_if<std::is_arithmetic<T>::value, scon::mathmatrix<T>>::type
  Rep3D_to_arma(const ContainerType<CoordType<T>, ContainerArgs ...>& rep) {
    using Mat = scon::mathmatrix<T>;

    Mat A(rep.size(), 3);
    for (std::size_t i = 0; i != rep.size(); ++i) {
      A(i, 0) = rep.at(i).x();
      A(i, 1) = rep.at(i).y();
      A(i, 2) = rep.at(i).z();
    }
    return A;
  }

  // mean function for coordinates
  inline coords::Cartesian_Point cp_mean(const coords::Representation_3D& rep) {
    coords::Cartesian_Point cp{};
    for (auto&& i : rep) {
      cp += i;
    }
    cp /= static_cast<float_type>(rep.size());
    return cp;
  }

  // radius of gyration from Rep3D
  // coords have to be in Bohr
  template <typename T>
  inline typename std::enable_if<std::is_arithmetic<T>::value, T>::type
  rad_gyr(const coords::Representation_3D& rep) {
    auto mean = cp_mean(rep);
    auto diff = rep - mean;
    auto sum{ 0.0 };
    for (auto& i : diff) {
      auto temp = std::pow(i.x(), 2) + std::pow(i.y(), 2) + std::pow(i.z(), 2);
      sum += temp;
    }
    sum /= rep.size();
    return std::sqrt(sum);
  }

  inline coords::Cartesian_Point
  normal_unit_vector(const coords::Cartesian_Point& a,
                     const coords::Cartesian_Point& b,
                     const coords::Cartesian_Point& c) {
    using scon::cross;
    using scon::len;

    auto a_vec = b - a;
    auto b_vec = b - c;
    auto t1 = cross(a_vec, b_vec);
    auto result = t1 / (len(t1));
    return result;
  }

  /*!
  \brief Creates all possible 2-permutations from a maximum number and sorts them
  according to size.
  \details Used as a helper function to create a vector of bonded atom pairs.
  \param num Maximum number used for permutations.
  \return Sorted std::vector of std::pair.
  */
  inline std::vector<std::pair<int, int>>
  permutation(const std::size_t& num) {
    std::vector<int> s(num);
    // std::iota fills the range from std::begin(s) to std::end(s) with
    // incremented integer values starting from 1.
    std::iota(std::begin(s), std::end(s), 1);
    std::vector<std::pair<int, int>> perm_vec;
    for (auto& i : s) {
      for (auto& f : s) {
        if (i < f) {
          auto par = std::make_pair(i, f);
          perm_vec.emplace_back(par);
        }
      }
    }
    return perm_vec;
  }

  /*!
  \brief Calculates the Euclidean distance between two Cartesian_Points.
  \tparam T Arithmetic type. Enforced by type traits.
  \param a First Cartesian_Point.
  \param b Second Cartesian_Point.
  \return Euclidean distance.
  */
  template <typename T>
  inline typename std::enable_if<std::is_arithmetic<T>::value, T>::type
  euclid_dist(const coords::Cartesian_Point& a, const coords::Cartesian_Point& b) {
    auto f{ (a - b) * (a - b) };
    auto dist = std::sqrt(f.x() + f.y() + f.z());
    return dist;
  }

  /*!
  \brief Creates a std::vector of Euclidean distances for each possible pair
  of atoms.
  \tparam T Arithmetic type. Enforced by type traits.
  \param cp_vec Structure of the system.
  \return std::vector of Euclidean distances.
  */
  template <typename T>
  inline typename std::enable_if<std::is_arithmetic<T>::value, std::vector<T>>::type
  dist_vec(const coords::Representation_3D& cp_vec) {
    std::vector<T> result;
    auto perm_vec = permutation(cp_vec.size());
    for (auto& i : perm_vec) {
      auto dist = euclid_dist<T>(cp_vec.at(i.first - 1), cp_vec.at(i.second - 1));
      result.emplace_back(dist);
    }
    return result;
  }

  /*!
  \brief Normalizes a vector given as Cartesian_Point.
  \param a Vector, which shall be normalized.
  \return Normalized vector.
  */
  inline coords::Cartesian_Point normalize(const coords::Cartesian_Point& a) {
    return a / scon::geometric_length(a);
  }

  template <typename Line, typename T>
  inline std::vector<std::pair<int, int>>
  bonds(const std::vector<Pdb::Atom<Line, T>>& vec,
        const coords::Representation_3D& cp_vec) {
    using ic_atom::radius_vec;

    auto radii = radius_vec(vec);
    auto vec_size = vec.size();
    auto perm = permutation(vec_size);
    std::vector<double> thres;
    for (const auto& i : perm) {
      auto t = 1.2 * (radii.at(i.first - 1) + radii.at(i.second - 1));
      thres.emplace_back(t);
    }
    auto dist = dist_vec<T>(cp_vec);
    auto thres_it{ thres.cbegin() };
    auto dist_it{ dist.cbegin() };
    std::vector<std::pair<int, int>> result;
    for (; dist_it != dist.cend() && thres_it != thres.cend();
         ++dist_it, ++thres_it) {
      if (*dist_it < *thres_it) {
        auto pos = std::distance(dist.cbegin(), dist_it);
        result.emplace_back(perm.at(pos));
      }
    }
    return result;
  }

  /*!
  \brief Creates all possible 3-permutations from a vector containing 3 integers.
  \param vec Vector containing 3 integers.
  \return Vector of all created permutations, where each permutation is itself a
  3-dimensional vector.
  */
  inline std::vector<std::vector<unsigned int>>
  permutation_from_vec(std::vector<unsigned int>& vec) {
    std::vector<std::vector<unsigned int>> result;
    if (vec.size() != 3) {
      throw("Function permutation_from_vec needs vector with size 3.");
    }
    do {
      std::vector<unsigned int> temp = { vec.at(0), vec.at(1), vec.at(2) };
      result.emplace_back(temp);
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
  template <typename T>
  inline typename std::enable_if<std::is_arithmetic<T>::value, std::vector<T>>::type
  flatten_c3_vec(const std::vector<scon::c3<T>>& vec) {
    std::vector<T> result;
    for (auto& i : vec) {
      result.emplace_back(i.x());
      result.emplace_back(i.y());
      result.emplace_back(i.z());
    }
    return result;
  }
}

/*! @} End of ic_util group*/

#endif // cast_ic_util_h_guard
