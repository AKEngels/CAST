#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include<limits>
#include<algorithm>
#include "coords.h"
#include "scon_traits.h"

#ifdef USE_PYTHON
#include <Python.h>
#endif
#include "coords.h"

// Define Function to output molar mass of a coords object
inline double sys_mass(coords::Coordinates &sys)
{
  double m = 0;
  for (auto const& a : sys.atoms())
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
  strm << std::fixed << std::setw(w) << coords.pes().energy;
  for (auto const& ia : coords.pes().ia_matrix)
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
@ param sep: char where the string should be splitted
@ param remove: removes empty elements (i.e. if more than one seperator directly follow each other they are treated as one)*/
inline std::vector<std::string> split(const std::string &text, char sep, bool remove=false) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }
  tokens.push_back(text.substr(start));
  if (remove == true)  // remove empty elements
  {
	  std::vector<std::string> tokens2;
	  for (auto t : tokens)
	  {
		  if (t != "") tokens2.push_back(t);
	  }
	  tokens = tokens2;
  }
  return tokens;
}

/**calculates the distance between two points in Cartesian Space*/
inline double dist(coords::Cartesian_Point a, coords::Cartesian_Point b)
{
  return sqrt((a.x() - b.x())*(a.x() - b.x()) + (a.y() - b.y())*(a.y() - b.y()) + (a.z() - b.z())*(a.z() - b.z()));
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
returns true if yes and false if no  (overloaded function)*/
template<typename T, typename U, template<typename, typename ...> class Cont, typename ... ContArgs>
inline typename std::enable_if<scon::is_container<Cont<U, ContArgs...>>::value || std::is_same<Cont<U, ContArgs...>, std::string>::value, bool>::type
is_in(T const& x, Cont<U, ContArgs...> const& v) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

/**looks if vector v contains element x
returns true if yes and false if no (overloaded function)*/
template<typename T, typename U, std::size_t N>
inline bool is_in(T const& x, std::array<U, N> const& v) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

/**finds index of element x in vector v
if not inside it returns the maximum limit of an integer*/
template<typename T, template<typename, typename ...> class Cont, typename ... ContArgs>
inline typename std::enable_if<scon::is_container<Cont<T, ContArgs...>>::value || std::is_same<Cont<T, ContArgs...>, std::string>::value, int>::type
find_index(T const & x, Cont<T, ContArgs...> v) {
  auto found = std::find(v.begin(), v.end(), x);
  if (found != v.end()) return found - v.begin();
  else return std::numeric_limits<int>::max();
}

/**tests if a string is a number*/
inline bool check_if_number(std::string const & number) {
  return !number.empty() && std::find_if(number.cbegin(), number.cend(), [](char n) {
    return n != 'E' && n != 'e' && n != '-' && n != '+' && n != '.' && !std::isdigit(n); //check if the line contains digits, a minus or a dot to determine if its a floating point number
  }) == number.end();
}

/**tests if a (one-letter) string is a digit*/
inline bool isdigit(std::string s)
{
  return check_if_number(s);
}

/**tests if a file exists
@param name: name of the file*/
inline bool file_exists(const std::string& name) {
  std::ifstream f(name.c_str());
  return f.good();
}

/**returns last line of a file*/
inline std::string last_line(std::ifstream& in)
{
  std::string line;
  while (in >> std::ws && std::getline(in, line));
  return line;
}

/**tests if a file is empty
returns true if it is empty, false if not
@param filename: name of the file*/
inline bool file_is_empty(std::string &filename)
{
	std::ifstream file (filename);
	std::string some_string;
	file >> some_string;
	if (some_string.empty()) return true;
	else return false;
}

/**adds two vectors like in python, i.e. [A, B] + [C, D] = [A, B, C, D]
@param v1: first vector
@param v2: second vector
@param sort: if true sort resulting vector with the std::sort-function*/
template <typename T>
inline std::vector<T> add_vectors(std::vector<T> const &v1, std::vector<T> const &v2, bool sort = false)
{
  std::vector<T> v12;
  v12.reserve(v1.size() + v2.size());
  v12.insert(v12.end(), v1.begin(), v1.end());
  v12.insert(v12.end(), v2.begin(), v2.end());

	if (sort) std::sort(v12.begin(), v12.end());
  return v12;
}

/**tests if a vector contains at least one element twice (or more)
returns false if no element is in vector more than once, returns true otherwise
@param v: vector that is to be tested*/
template <typename T>
inline bool double_element(std::vector<T> const &v)
{
	for (auto i = 0u; i < v.size(); ++i)
	{
		for (auto j = 0u; j < i; ++j)
		{
			if (v[i] == v[j]) return true;
		}
	}
	return false;
}

/**function analogous to python range function (see https://stackoverflow.com/questions/13152252/is-there-a-compact-equivalent-to-python-range-in-c-stl) 
Attention! The order of start and stop is switched, so you can use:
range(10) = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
range(10, 5) = [5, 6, 7, 8, 9,]
range(11, 5, 2) = [5, 7, 9]*/
template <typename IntType>
std::vector<IntType> range(IntType stop, IntType start=0, IntType step=1)
{
	if (step == IntType(0))
	{
		throw std::invalid_argument("step for range must be non-zero");
	}

	std::vector<IntType> result;
	IntType i = start;
	while ((step > 0) ? (i < stop) : (i > stop))
	{
		result.push_back(i);
		i += step;
	}

	return result;
}

/**function to count an element in a vector
@param e: element that should be counted
@param vec: vector in which should be counted*/
template <typename T, typename U>
int count_element(T const &e, std::vector<U> const &vec)
{
	return std::count(vec.begin(), vec.end(), e);
}

#endif
