#ifndef cast_ic_util_h_guard
#define cast_ic_util_h_guard

#include<cstring>
#include<type_traits>
#include<vector>
#include<array>

// TODO: find some possibility to include the values from "energy.h" without the unnecessary stuff.
namespace energy {
	internals::float_type constexpr bohr2Ang() {
		return 0.52917721067;
	}
	internals::float_type constexpr hartreePerBor2KcalPerMolAng(){
		return 627.5096080306 / bohr2Ang();
	}
}

/*!
 *  \addtogroup ic_util
 *  @{
 */
namespace ic_util{

	inline bool isSameSet(std::vector<std::size_t> lhs, std::vector<std::size_t> rhs) {
		if (lhs.size() != rhs.size()) return false;
		std::sort(lhs.begin(), lhs.end());
		std::sort(rhs.begin(), rhs.end());
		return std::includes(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
	}

  inline coords::Representation_3D grads_to_bohr(coords::Representation_3D const& grads) {
	  coords::Representation_3D bohr_grads;
	  bohr_grads.reserve(grads.size());
	  for (auto const& g : grads) {
		  bohr_grads.emplace_back(g / hartreePerBor2KcalPerMolAng());
	  }
	  return bohr_grads;
  }
  inline coords::Representation_3D rep3d_bohr_to_ang(coords::Representation_3D const& bohr) {
	  coords::Representation_3D ang;
	  ang.reserve(bohr.size());
	  for (auto const& b : bohr) {
		  ang.emplace_back(b * bohr2Ang());
	  }
	  return ang;
  }

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

  template <template<typename> class MatrixType, typename T, template<typename> class CoordType, template <typename, typename ...> class ContainerType, typename ... ContainerArgs>
  inline typename std::enable_if<std::is_arithmetic<T>::value, MatrixType<T>>::type
  Rep3D_to_Mat(const ContainerType<CoordType<T>, ContainerArgs ...>& rep) {
    using Mat = MatrixType<T>;

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

  

  /*!
\brief Number of elements currently supported by the code.
*/
  static constexpr auto N_elem{ 97u };

  /*!
  \brief Covalent radii of the first N_elem elements of the PSE.
  \details Corresponds to the symbol array.
  \see "Covalent radii revisited", Cordero et al., Dalton Trans., 21, 2008,
  pp. 2832-2838
  */
  static constexpr std::array<double, N_elem> radius{
	0.00, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41,
	1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.61,
	1.52, 1.50, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95,
	1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39,
	1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96,
	1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41,
	1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06,
	2.00, 1.96, 1.90, 1.87, 1.80, 1.69
  };

  /*!
  \brief List containing the first N_elem of the PSE.
  \details Corresponds to the radii of the radius array.
  */
  static constexpr std::array<char const*, N_elem> symbol{
	"DD", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
	"Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn",
	"Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
	"Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
	"Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
	"Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir",
	"Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
	"Pa", "U",  "Np", "Pu", "Am", "Cm"
  };

  

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

  template<typename Mat>
  Mat atomsNorm(Mat const& norm) {
	  Mat mat(norm.rows(), 1);
	  for (auto i = 0u; i < norm.rows(); ++i) {
		  mat(i, 0) = norm.row(i).norm();
	  }
	  return mat;
  }

}

/*! @} End of ic_util group*/

#endif // cast_ic_util_h_guard
