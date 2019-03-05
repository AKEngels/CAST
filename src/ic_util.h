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


//TODO make it independent of scon::mathmatrix!
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

  //not tested
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

  struct AtomConnector {
    using returnType = std::vector<std::pair<std::size_t, std::size_t>>;

    AtomConnector(std::vector<std::string> const& elem_vec,
      coords::Representation_3D const& cp_vec) : sequenceOfSymbols{ elem_vec }, cartesianRepresentation{ cp_vec }, firstAtomIndex{ 0u }, secondAtomIndex{ 0u } {}
    returnType operator()();

  private:
    double getThresholdForBeingNotConnected(std::string const& oneAtom, std::string const& otherAtom);
    void findTheFirstAtom(returnType & connectedAtoms);
    void findAtomWithIndexHigherThanFirst(returnType & connectedAtoms);
    void connectIfCloseEnough(returnType & connectedAtoms);
    bool areTheyCloseEnough();

    std::vector<std::string> const& sequenceOfSymbols;
    coords::Representation_3D const& cartesianRepresentation;

    std::size_t firstAtomIndex;
    std::size_t secondAtomIndex;
  };

  inline std::vector<std::pair<std::size_t, std::size_t>> AtomConnector::operator()() {
    std::vector<std::pair<std::size_t, std::size_t>> connectedAtoms;
    findTheFirstAtom(connectedAtoms);
    return connectedAtoms;
  }

  inline void AtomConnector::findTheFirstAtom(returnType & connectedAtoms){
    for (firstAtomIndex = 0u; firstAtomIndex < sequenceOfSymbols.size(); ++firstAtomIndex) {
      findAtomWithIndexHigherThanFirst(connectedAtoms);
    }
  }

  inline void AtomConnector::findAtomWithIndexHigherThanFirst(returnType & connectedAtoms){
    for (secondAtomIndex = firstAtomIndex+1u; secondAtomIndex < sequenceOfSymbols.size(); ++secondAtomIndex) {
      connectIfCloseEnough(connectedAtoms);
    }
  }

  inline void AtomConnector::connectIfCloseEnough(returnType & connectedAtoms) {
    if (areTheyCloseEnough()) {
      connectedAtoms.emplace_back(firstAtomIndex, secondAtomIndex);
    }
  }

  inline bool AtomConnector::areTheyCloseEnough() {
    double const threshold = getThresholdForBeingNotConnected(sequenceOfSymbols.at(firstAtomIndex), sequenceOfSymbols.at(secondAtomIndex));
    double const actualDistance = scon::len(cartesianRepresentation.at(firstAtomIndex) - cartesianRepresentation.at(secondAtomIndex));
    return actualDistance < threshold;
  }

  inline double AtomConnector::getThresholdForBeingNotConnected(std::string const& oneAtom, std::string const& otherAtom) {
    using ic_atom::element_radius;
    return 1.2 * (element_radius(oneAtom) + element_radius(otherAtom));
  }

  inline std::vector<std::pair<std::size_t, std::size_t>>
  bonds(std::vector<std::string> const& elem_vec,
        coords::Representation_3D const& cp_vec) {

    AtomConnector atomCreator(elem_vec, cp_vec);
    return atomCreator();

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

  namespace {
    struct MatrixConverterInterface {
      virtual coords::Representation_3D toRep3D() const = 0;
    };

    template<typename Mat>
    struct FlatMatrix : MatrixConverterInterface {
      FlatMatrix(Mat const& mat) : matrix(mat) {
        if (matrix.cols() != 1)
          throw std::runtime_error("Passed a Matrix as flattened Matrix which is not flat at all! Error in class FlatMatrix most likely thrown in MatToRep3D.");
      }
      Mat const& matrix;

      coords::Representation_3D toRep3D() const override {
        coords::Representation_3D result;
        for (auto i = 0u; i < matrix.rows(); i += 3) {
          result.emplace_back(matrix.operator()(i, 0u), matrix.operator()(i+1u, 0u), matrix.operator()(i+2u, 0u));
        }
        return result;
      }
    };
    template<typename Mat>
    struct NotFlatMatrix : MatrixConverterInterface {
      NotFlatMatrix(Mat const& mat) : matrix(mat) {
        if (matrix.cols() != 3)
          throw std::runtime_error("Passed a Matrix as to convert to a Rep3D object, but the matrix has not three rows! Error in class NotFlatMatrix  most likely thrown in MatToRep3D.");
      }
      Mat const& matrix;

      coords::Representation_3D toRep3D() const override {
        coords::Representation_3D result;
        for (auto i = 0u; i < matrix.rows(); ++i) {
          result.emplace_back(matrix(i, 0u), matrix(i, 1u), matrix(i, 2u));
        }
        return result;
      }
    };
  }

  template<typename Mat>
  void printMat(Mat const& mat) {
    for (auto i = 0u; i < mat.rows(); ++i) {
      for (auto j = 0u; j < mat.cols(); ++j) {
        std::cout << mat(i, j) << " ";
      }
      std::cout << "\n";
    }
  }

  template<typename Mat>
  coords::Representation_3D matToRep3D(Mat const& mat) {
    std::unique_ptr<MatrixConverterInterface> converter;
    if (mat.cols() == 1) {
      converter = std::make_unique<FlatMatrix<Mat>>(mat);
    }
    else {
      converter = std::make_unique<NotFlatMatrix<Mat>>(mat);
    }
    return converter->toRep3D();
  }
}

/*! @} End of ic_util group*/

#endif // cast_ic_util_h_guard
