#pragma once
#include "coords.h"
#ifdef USE_PYTHON
#include <Python.h>
#endif

// Define Function to output molar mass of a coords object
inline double sys_mass(coords::Coordinates &sys)
{
  double m = 0;
  for (auto && a : sys.atoms())
  {
    m += a.mass();
  }
  return m;
}


// Energy print functions
inline void short_ene_stream(
  coords::Coordinates const &coords,
  std::ostream &strm, std::streamsize const w)
{
  strm << std::setw(w) << coords.pes().energy;
  for (auto && ia : coords.pes().ia_matrix)
  {
    strm << std::setw(w) << ia.energy;
  }
}

inline void short_ene_stream_h(
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

/**splits a string into a vector of strings
@ param text: string that is to be splitted
@ param sep: char where the string should be splitted*/
inline std::vector<std::string> split(const std::string &text, char sep) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }
  tokens.push_back(text.substr(start));
  return tokens;
}

/**calculates the distance between two points in Cartesian Space*/
inline double dist(coords::Cartesian_Point a, coords::Cartesian_Point b)
{
  return sqrt( (a.x()-b.x())*(a.x()-b.x()) + (a.y()-b.y())*(a.y()-b.y()) + (a.z()-b.z())*(a.z()-b.z()) );
}

#ifdef USE_PYTHON
/*
function that returns the path to a pythonmodule
(path has to be appended to pythonpath if you want to call this module)
@param modulename: name of the module
*/
inline std::string get_python_modulepath(std::string modulename)
{
  std::string find = "import " + modulename + "\nwith open('tmpfile.txt','w') as fn:\n    fn.write(" + modulename + ".__file__)";
  const char *c_find = find.c_str();
  PyRun_SimpleString(c_find);  //call a python programme to find the modulepath and write it to tmpfile
  std::ifstream file("tmpfile.txt");  //open tmpfile and read content
  std::string content;
  file >> content;
  file.close();
  remove("tmpfile.txt");
  return content.substr(0, content.size() - 14 - modulename.size());  //give back path without filename __init__.pyc and modulename
}
#endif

/**looks if vector v contains element x
returns true if yes and false if no */
inline bool is_in(std::vector<std::string> x, std::vector<std::vector<std::string>> v)
{
  if (std::find(v.begin(), v.end(), x) != v.end()) {
    return true;
  }
  else {
    return false;
  }
}

/**tests if a file exists
@param name: name of the file*/
inline bool file_exists(const std::string& name) {
  std::ifstream f(name.c_str());
  return f.good();
}
