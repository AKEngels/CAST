#pragma once
#include "coords.h"

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