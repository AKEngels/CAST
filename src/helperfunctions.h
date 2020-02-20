/**
CAST 3
helperfunctions.h
Purpose: some functions that are helpful in general

@version 1.0
*/

#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include<limits>
#include<algorithm>
#include "coords.h"
#include "Scon/scon_traits.h"

#ifdef USE_PYTHON
#include <Python.h>
#endif


/**function to build up a vector with the element symbols of the bonding partners of an atom
@param a: atom
@param atoms: vector of atoms (needed to get the element symbol)*/
inline std::vector<std::string> get_bonding_symbols(coords::Atom& a, coords::Atoms& atoms)
{
  std::vector<std::string> result;
  for (auto b : a.bonds())
  {
    result.push_back(atoms.atom(b).symbol());
  }
  return result;
}

// Energy print functions
inline void short_ene_stream(
  coords::Coordinates const& coords,
  std::ostream& strm, std::streamsize const w)
{
  strm << std::fixed << std::setw(w) << coords.pes().energy;
  for (auto const& ia : coords.pes().ia_matrix)
  {
    strm << std::setw(w) << ia.energy;
  }
}

inline void short_ene_stream_h(
  coords::Coordinates const& coords,
  std::ostream& strm, std::streamsize const w)
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
inline std::vector<std::string> split(std::string const& text, char const sep, bool const remove = false) {
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

/**removes all spaces from string
@param str: inputstring*/
inline std::string remove_spaces(std::string const& str)
{
  std::string result{ "" };
  for (auto s : str)
  {
    if (s != ' ') result += s;
  }
  return result;
}

/**looks if a string consists only of whitespace characters
@param str: inputstring*/
inline bool only_whitespace(std::string const& str)
{
  for (char c : str) {
    if (c != ' ' && c != '\n' && c != '\r' && c != '\t') return false;
  }
  return true;
}

/**calculates the distance between two points in Cartesian Space*/
inline double dist(coords::Cartesian_Point const& a, coords::Cartesian_Point const& b)
{
  return sqrt((a.x() - b.x()) * (a.x() - b.x()) + (a.y() - b.y()) * (a.y() - b.y()) + (a.z() - b.z()) * (a.z() - b.z()));
}

#ifdef USE_PYTHON
/*
function that returns the path to a pythonmodule
(path has to be appended to pythonpath if you want to call this module)
@param modulename: name of the module
*/
inline std::string get_python_modulepath(std::string const& modulename)
{
  std::string find = "import " + modulename + "\nwith open('tmpfile.txt','w') as fn:\n    fn.write(" + modulename + ".__file__)";
  const char* c_find = find.c_str();
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

/**looks if element x is in any vector of vec (which is a vector of vectors)
returns true if yes and false if no*/
template <typename T>
bool is_in_any(T const& x, std::vector<std::vector<T>> const& vec)
{
  for (auto const& v : vec) {
    if (is_in(x, v)) return true;
  }
  return false;
}

/**finds index of element x in vector v
if not inside it returns the maximum limit of an integer*/
template<typename T, typename U>
inline std::size_t find_index(T const& x, std::vector<U> const& v) {
  auto found = std::find(v.begin(), v.end(), x);
  if (found != v.end()) return found - v.begin();
  else return std::numeric_limits<int>::max();
}

/**tests if a string is a number*/
inline bool check_if_number(std::string const& number) {
  size_t idx{ 0u };  // given to std::stod, gives afterward position in string which is behind double

  try { std::stod(number, &idx); }  // try to convert string to double
  catch (...) { return false; }  // if it doesn't work -> false

  if (idx == number.size()) return true;  // if it works and whole string has been converted -> true

  else // if not whole string has been converted
  {
    for (auto i = idx; i < number.size(); ++i)  // look if there is something else than whitespace behind
    {
      if (number[i] != ' ' && number[i] != '\n' && number[i] != '\t') return false; // if yes -> false
    }
    return true;                                                                    // if no -> true
  }
}

/**tests if a string is an integer*/
inline bool check_if_integer(std::string const& number) {
  if (check_if_number(number) && std::floor(std::stod(number)) == std::ceil(std::stod(number))) return true;
  else return false;
}

/**tests if a (one-letter) string is a digit*/
inline bool isdigit(std::string const& s)
{
  return check_if_number(s);
}

/**tests if a file exists
@param name: name of the file*/
inline bool file_exists(std::string const& name) {
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
inline bool file_is_empty(std::string const& filename)
{
  std::ifstream file(filename);
  std::string some_string;
  file >> some_string;
  if (some_string.empty()) return true;
  else return false;
}

/**adds two vectors like in python, i.e. [A, B] + [C, D] = [A, B, C, D]
@param v1: first vector
@param v2: second vector
@param sort: if true sort resulting vector with the std::sort-function*/
template <typename T, typename U>
inline std::vector<T> add_vectors(std::vector<T> const& v1, std::vector<U> const& v2, bool const sort = false)
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
inline bool double_element(std::vector<T> const& v)
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

/**sorts a vector and removes double elements from it
(taken from https://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector)
@param vec: vector, given as reference as it is changed*/
template <typename T>
void sort_and_remove_double_elements(std::vector<T>& vec)
{
  sort(vec.begin(), vec.end());
  vec.erase(unique(vec.begin(), vec.end()), vec.end());
}

/**tests if a two vectors contain the same element
(vectors not given as reference as they are changed in function and those changes should not be kept)
@param v1: first vector
@param v2: second vector*/
template <typename T>
inline bool double_element(std::vector<T> v1, std::vector<T> v2)
{
  sort_and_remove_double_elements(v1);
  sort_and_remove_double_elements(v2);
  auto combined_vector = add_vectors(v1, v2);
  return double_element(combined_vector);
}

/**This function takes a vector and a bunch of vectors and adds those vectors to the first which have common elements with it.
The resulting vector does contain each element only once.
As the function is called recursively also those vectors are added that have common elements with the newly added vectors.
Only those vectors are checked that are not marked as done in vector<bool> 'done' and vectors that are added are marked as done.
@param current: starting vector (given by reference as it is modified)
@param vecs: vector of vectors which are to compared with current vector
@param done: vectors of bools in which those elements that should not be added or that are already be added are marked as done.
(part of function combine_vectors())*/
template <typename T>
inline void add_vectors_with_common_elements(std::vector<T>& current, std::vector<std::vector<T>> const& vecs, std::vector<bool>& done)
{
  for (auto j{ 0u }; j < vecs.size(); ++j)
  {
    if (done[j] == false)
    {
      if (double_element(current, vecs[j]) == true)
      {
        current = add_vectors(current, vecs[j]);
        sort_and_remove_double_elements(current);
        done[j] = true;
        add_vectors_with_common_elements(current, vecs, done);
      }
    }
  }
}

/**
This function looks at a bunch of vectors and combines those which have at least one common element.
@param vec: vector of the vectors that should be investigated
returns a vector of the combined vectors (without any double elements)
*/
template <typename T>
std::vector<std::vector<T>> combine_vectors(std::vector<std::vector<T>> const& vec)
{
  // create a vector of bools that tracks which vectors are already included in the result vector
  std::vector<bool> done;
  done.resize(vec.size());
  for (auto&& d : done) d = false;

  // create result vector
  std::vector<std::vector<T>> result;
  for (auto i{ 0u }; i < vec.size(); ++i)   // for every element in vector...
  {
    if (done[i] == false)  //...that is not already included into result
    {
      auto current = vec[i];  // create a new vector
      done[i] = true;

      add_vectors_with_common_elements(current, vec, done); // add all other vectors that have common elements (also over several steps) to 'current' and track that they are included
      result.emplace_back(current);  // add current vector to result
    }
  }
  return result;
}

/**function analogous to python range function (see https://stackoverflow.com/questions/13152252/is-there-a-compact-equivalent-to-python-range-in-c-stl)
Attention! The order of start and stop is switched, so you can use:
range(10) = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
range(10, 5) = [5, 6, 7, 8, 9,]
range(11, 5, 2) = [5, 7, 9]*/
template <typename IntType>
std::vector<IntType> range(IntType const stop, IntType const start = 0, IntType const step = 1)
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
int count_element(T const& e, std::vector<U> const& vec)
{
  return std::count(vec.begin(), vec.end(), e);
}


/**function that does the same as (a < b)
but takes into account that doubles have a precision*/
inline bool is_smaller_than(double a, double b, double precision = 1e-10)
{
  if (fabs(a - b) < precision) return false;
  else if (a < b) return true;
  else return false;
}

/**function that finds the ideal number of bonding partners for an atom out of the parameterfile
@param atomtype: forcefield atomtype
@param paramfile: name of the forcefield parameterfile (default is what is given in inputfile)
returns number of bonds*/
inline unsigned int get_ideal_bond_number_from_parameterfile(int atomtype, std::string const& paramfile = Config::set().general.paramFilename)
{
  if (file_exists(paramfile) == false) {
    throw std::runtime_error("Parameterfile " + paramfile + " not found to get number of bonds for atomtype " + std::to_string(atomtype));
  }

  std::ifstream in_file(paramfile, std::ios_base::in);
  std::string line;
  std::vector<std::string> linevec;

  while (!in_file.eof())
  {
    std::getline(in_file, line);
    if (line.size() > 4 && line.substr(0, 4) == "atom")
    {
      linevec = split(line, ' ', true);
      if (std::stoi(linevec[1]) == atomtype) return std::stoi(line.substr(72, 5));
    }
  }
  throw std::runtime_error("Atomtype " + std::to_string(atomtype) + " not found");
}

// Usage: stream >> skipline;
// Skips a line
inline std::istream& skipline(std::istream& in)
{
  return in.ignore(std::numeric_limits < std::streamsize >::max(), '\n');
}

/**convert a vector to a string where the single elements are seperated by seperator 'sep'*/
template<typename T> 
inline std::string vec_to_string(std::vector<T> const &vec, std::string const& sep)
{
  std::string result;
  for (auto i{ 0u }; i < vec.size() - 1; ++i) result += std::to_string(vec[i]) + sep; // all elements except last one
  result += std::to_string(vec[vec.size() - 1]);  // last element of vector
  return result;
}

/**converts a number to a string that contains no dots (i.e. 1.7 gets "17")*/
inline std::string convert_number_to_string_without_dots(double number)
{
  std::string s = std::to_string(number);
  std::string n;
  for (char c : s) {
    if (c != '.') n += c;
  }
  return n;
}

/**function that checks if two values are equal within a given tolerance
returns true if they are about equal and false if not
@param one: first value
@param two: second value
@param precision: precision for comparison
This is the 'original' function which will be called by all overloaded functions earlier or later.*/
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
bool is_nearly_equal(T one, T two, double precision = 0.0000000001)
{
  if (std::abs(one - two) < precision) return true;
  else return false;
}

/**overloaded function that checks if two values are equal within a given tolerance
returns true if they are about equal and false if not
@param one: first value
@param two: second value
@param precision: precision for comparison*/
template <typename T>
bool is_nearly_equal(scon::c3<T> one, scon::c3<T> two, double precision = 0.0000000001)
{
  if (is_nearly_equal(one.x(), two.x(), precision) && is_nearly_equal(one.y(), two.y(), precision) && is_nearly_equal(one.z(), two.z(), precision)) {
    return true;
  }
  else return false;
}

/**overloaded function that checks if two values are equal within a given tolerance
returns true if they are about equal and false if not
@param one: first value
@param two: second value
@param precision: precision for comparison*/
template <typename T>
bool is_nearly_equal(std::vector<T> vec1, std::vector<T> vec2, double precision = 0.0000000001)
{
  if (vec1.size() != vec2.size()) return false;
  for (auto i{ 0u }; i < vec1.size(); ++i)
  {
    if (is_nearly_equal(vec1[i], vec2[i], precision) == false) return false;
  }
  return true;
}

inline std::size_t gap(std::size_t const largeNum, std::size_t const smallNum)
{
  if (smallNum == 0u) return 0u;
  return largeNum / std::min(largeNum, smallNum);
}

#endif
