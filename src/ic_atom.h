#ifndef cast_ic_atom_h_guard
#define cast_ic_atom_h_guard

#pragma once

#include "pdb.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

/*!
 *  \addtogroup ic_atom
 *  @{
 */
namespace ic_atom {

/*!
\brief Number of elements currently supported by the code.
*/
static constexpr auto N_elem{ 97u };

/*!
\brief Factor for conversion from angstrom to bohr.
\details Simply divide angstrom by this factor.
*/
static constexpr auto bohr{ 0.52917721067 };

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
\brief Elements of the first period of the PSE.
*/
static constexpr std::array<char const*, 3u> period_one{ "DD", "H", "He" };

/*!
\brief Elements of the second period of the PSE.
*/
static constexpr std::array<char const*, 8u> period_two{
  "Li", "Be", "B", "C", "N", "O", "F", "Ne"
};

/*!
\brief Elements of the third period of the PSE.
*/
static constexpr std::array<char const*, 8u> period_three{ "Na", "Mg", "Al",
                                                           "Si", "P",  "S",
                                                           "Cl", "Ar" };

/*!
\brief Helper function for the radius_vec() function.
\param sym Element symbol.
\return Distance between symbol.begin() and the iterator that points to the
position of sym in the symbol array.
*/
inline auto atomic_num(const std::string& sym) {
  auto result = std::find(symbol.begin(), symbol.end(), sym);
  if (result == symbol.end()) {
    throw std::runtime_error("Element unknown.");
  }
  return std::distance(symbol.begin(), result);
}

/*!
\brief Creates a vector of radii from a vector of atoms.
\tparam Line Type used for PDB lines by the PDB parser.
\tparam T Arithmetic type required by the Pdb::Atom<> class for coordinates.
\param vec std::vector of atoms created by the PDB parser.
\return std::vector of radii.
*/
template <typename Line, typename T>
inline std::vector<double>
radius_vec(const std::vector<Pdb::Atom<Line, T>>& vec) {
  std::vector<double> result;
  for (auto& i : vec) {
    auto n = atomic_num(i.element);
    result.emplace_back(radius.at(n));
  }
  return result;
}
}

/*! @} End of ic_atom group*/

#endif // cast_ic_atom_h_guard