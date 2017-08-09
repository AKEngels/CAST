#pragma once
#include "coords.h"
#include <iostream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>

// Define Function to output molar mass of a coords object
double sys_mass(coords::Coordinates &sys)
{
  double m = 0;
  for (auto && a : sys.atoms())
  {
    m += a.mass();
  }
  return m;
}


// Energy print functions
void short_ene_stream(
  coords::Coordinates const &coords,
  std::ostream &strm, std::streamsize const w)
{
  strm << std::setw(w) << coords.pes().energy;
  for (auto && ia : coords.pes().ia_matrix)
  {
    strm << std::setw(w) << ia.energy;
  }
}

void short_ene_stream_h(
  coords::Coordinates const &coords,
  std::ostream &strm, std::streamsize const w)
{
  strm << std::setw(w) << "Energy";
  auto const n = coords.pes().ia_matrix.size();
  for (std::size_t i = 0; i < n; ++i)
  {
    strm << std::setw(w) << ("WW" + std::to_string(i));
  }
}

// VIA http://www.cplusplus.com/forum/general/209784/
// Computes the distance between two std::vectors
template <typename T>
T	eucledeanDistance(const std::vector<T>& a, const std::vector<T>& b)
{
  std::vector<T>	auxiliary;

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(auxiliary),//
    [](T element1, T element2) {return pow((element1 - element2), 2); });

  return std::sqrt(std::accumulate(auxiliary.begin(), auxiliary.end(), 0.0));

} // end template vectors_distance
