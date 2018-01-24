#ifndef cast_pdb_h_guard
#define cast_pdb_h_guard

#include "coords_rep.h"

#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

/*!
 *  \addtogroup Pdb
 *  @{
 */
namespace Pdb {
/// header-only (i.e. no separate implementation file provided)

/*!
\brief Scoped enumeration type with the relevant Pdb entries as possible values.
*/
enum class Record : int {
  ATOM,
  HETATOM,
  NUMMDL,
  MODEL,
  ENDMDL,
  TER,
  END,
  UNKNOWN
};

/*!
\brief Helper function in order to retrieve information about the value of the
Record type.
\param v Record object.
\return String representation of the value of the Record object.
*/
inline std::string Record_string(Record v) {
  switch (v) {
  case Record::ATOM:
    return "ATOM";
  case Record::HETATOM:
    return "HETAOM";
  case Record::NUMMDL:
    return "NUMMDL";
  case Record::MODEL:
    return "MODEL";
  case Record::ENDMDL:
    return "ENDMDL";
  case Record::TER:
    return "TER";
  case Record::END:
    return "TER";
  case Record::UNKNOWN:
    return "Unknown Record_type.";
  default:
    return "Unknown Record_type.";
  }
}

/*!
\brief Overloaded output operator for the Record enumeration type.
\param os Output stream.
\param c Record object.
\return String representation for the Record object piped to standard output.
*/
inline std::ostream& operator<<(std::ostream& os, Record c) {
  return os << Record_string(c);
}

/*!
\brief Simple class representing the Atom concept.
\details The members are initialized using various helper functions.
\tparam Line Helper type representing an arbitrary line of the Pdb file.
\tparam T Type intended to hold the coordinate data.
*/
template <typename Line, typename T>
struct Atom {
public:
  Atom() = default;

  Atom(const Line& func, const std::string& file)
      : rec_name{ func.line_type(file) }, atom_serial{ func.atom_serial(file) },
        atom_name{ func.atom_name(file) }, alt_loc{ func.alt_loc(file) },
        res_name{ func.res_name(file) }, chain_id{ func.chain_id(file) },
        res_seq{ func.res_seq(file) }, insertion_code{ func.insertion_code(
                                           file) },
        cp{ std::move(func.cart_point(file)) }, element{ func.element(
                                                    file) } {};

  Record rec_name;
  unsigned int atom_serial;
  std::string atom_name;
  std::string alt_loc;
  std::string res_name;
  char chain_id;
  unsigned int res_seq;
  std::string insertion_code;
  std::array<T, 3> coord = {};
  coords::Cartesian_Point cp = {};
  std::string element;

private:
  template <typename _Line, typename _T>
  friend std::ostream& operator<<(std::ostream& os, const Atom<_Line, _T>& atom);
};

/*!
\brief Overloaded output operator for an individual Atom object.
\tparam Line Helper type representing an arbitrary line of the Pdb file.
\tparam T Type intended to hold the coordinate data.
\param os Output stream.
\param atom Atom object.
\return String representation of the Atom object piped to the standard output.
*/
template <typename Line, typename T>
std::ostream& operator<<(std::ostream& os, const Atom<Line, T>& atom) {
  return os << std::setw(7) << atom.rec_name << ", " << std::setw(5)
            << atom.atom_serial << ", " << std::setw(4) << atom.atom_name
            << ", " << std::setw(9) << atom.alt_loc << ", " << std::setw(3)
            << atom.res_name << ", " << std::setw(1) << atom.chain_id << ", "
            << std::setw(4) << atom.res_seq << ", " << std::setw(8)
            << atom.insertion_code << ", " << std::setw(8) << atom.cp << ", "
            << std::setw(2) << atom.element;
}

/*!
\brief Collection of helper functions for parsing the Pdb lines. For ease of
handling the functions are contained within one class.
*/
struct Line {
public:
  /*!
  \brief Function for recognizing the Record value of a Pdb line.
  \param line String representation of the Pdb line.
  \return Specific value of the Record type.
  */
  Record line_type(const std::string& line) const {
    if (line.compare(0, 4, "ATOM") == 0) {
      return Record::ATOM;
    } else if (line.compare(0, 7, "HETATOM") == 0) {
      return Record::HETATOM;
    } else if (line.compare(0, 6, "NUMMDL") == 0) {
      return Record::NUMMDL;
    } else if (line.compare(0, 5, "MODEL") == 0) {
      return Record::MODEL;
    } else if (line.compare(0, 6, "ENDMDL") == 0) {
      return Record::ENDMDL;
    } else if (line.compare(0, 3, "TER") == 0) {
      return Record::TER;
    } else if (line.compare(0, 3, "END") == 0) {
      return Record::END;
    } else {
      return Record::UNKNOWN;
    }
  }

  /*!
  \brief Function for checking whether the Pdb line is an ATOM or HETATOM entry.
  Only these are relevant.
  \param line String representation of the Pdb line.
  \return True, if the line is an ATOM or HETATOM entry; otherwise false.
  */
  bool line_check(const std::string& line) const {
    auto type = line_type(line);
    if (type == Record::ATOM || type == Record::HETATOM) {
      return true;
    } else {
      return false;
    }
  }

  /*!
  \brief Checks whether the string to_check contains the string name.
  \param to_check String that is to be checked.
  \param name String whose existence in the string to_check shall be proven.
  \return std::runtime_error if to_check does not contain name.
  */
  void field_check(const std::string& to_check, const std::string& name) const {
    if (std::all_of(to_check.begin(), to_check.end(), isspace)) {
      throw std::runtime_error("At least one field [ " + name +
                               " ] in the PDB file is empty.");
    }
  }

  /*!
  \brief Extracts the model information from the Pdb line.
  \param line The Pdb line.
  \return Serial number of the model.
  */
  unsigned int model_serial(const std::string& line) const {
    std::string d{ "model serial" };
    std::string f = line.substr(10, 4);
    field_check(f, d);
    return std::stoi(f);
  }

  /*!
  \brief Extracts the atom serial number from the Pdb line.
  \param line The Pdb line.
  \return Atom serial number.
  */
  unsigned int atom_serial(const std::string& line) const {
    std::string d{ "atom serial" };
    std::string f = line.substr(6, 5);
    field_check(f, d);
    return std::stoi(f);
  }

  /*!
  \brief Extracts the atom name from the Pdb line.
  \param line The Pdb line.
  \return Atom name.
  */
  std::string atom_name(const std::string& line) const {
    std::string d{ "atom name" };
    std::string f = line.substr(12, 4);
    field_check(f, d);
    return f;
  }

  /*!
  \brief Extracts the alternate location descriptor from the Pdb line.
  \param line The Pdb line.
  \return Alternate loaction.
  */
  std::string alt_loc(const std::string& line) const {
    std::string d{ "alternate location" };
    std::string f = line.substr(16, 1);
    if (std::all_of(f.begin(), f.end(), isspace)) {
      return "No altLoc";
    }
    return f;
  }

  /*!
  \brief Extracts the residue name from the Pdb line.
  \param line The Pdb line.
  \return Residue name.
  */
  std::string res_name(const std::string& line) const {
    std::string d{ "residue name" };
    std::string f = line.substr(17, 3);
    field_check(f, d);
    return f;
  }

  /*!
  \brief Extracts the chain ID from the Pdb line.
  \param line The Pdb line.
  \return Chain ID.
  */
  char chain_id(const std::string& line) const {
    std::string d{ "chain ID" };
    std::string f = line.substr(21, 1);
    field_check(f, d);
    return f.at(0);
  }

  /*!
  \brief Extracts the residue sequence number from the Pdb line.
  \param line The Pdb line.
  \return Residue sequence number.
  */
  unsigned int res_seq(const std::string& line) const {
    std::string d{ "residue sequence number" };
    std::string f = line.substr(22, 4);
    field_check(f, d);
    return std::stoi(f);
  }

  /*!
  \brief Extracts the insertion code from the Pdb line.
  \param line The Pdb line.
  \return Insertion code.
  */
  std::string insertion_code(const std::string& line) const {
    std::string d{ "insertion code" };
    std::string f = line.substr(26, 1);
    if (std::all_of(f.begin(), f.end(), isspace)) {
      return "No iCode";
    }
    return f;
  }

  /*!
  \brief Extracts the element information from the Pdb line.
  \param line The Pdb line.
  \return Element information.
  */
  std::string element(const std::string& line) const {
    std::string d{ "element" };
    std::string f = line.substr(76, 2);
    field_check(f, d);
    auto end = std::remove(f.begin(), f.end(), ' ');
    f.erase(end, f.end());
    return f;
  }

  /*!
  \brief Extracts the coordinate information from the Pdb line.
  \tparam T Numerical type intended for the storage of the coordinate
  information.
  \param line The Pdb line.
  \return 3-dimensional coordinate array.
  */
  template <typename T>
  std::array<T, 3> coord(const std::string& line) const {
    std::string d{ "coord" };
    std::string f1 = line.substr(30, 8);
    std::string f2 = line.substr(38, 8);
    std::string f3 = line.substr(46, 8);
    std::array<std::string, 3> a{ { f1, f2, f3 } };
    std::array<T, 3> c = {};
    for (const auto& i : a) {
      field_check(i, d);
    }
    std::transform(a.begin(), a.end(), c.begin(),
                   [](const std::string& s) { return std::stod(s); });
    return c;
  }

  /*!
  \brief Extracts the coordinate information from the Pdb line.
  \param line The Pdb line.
  \return coords::Cartesian_Point object.
  */
  coords::Cartesian_Point cart_point(const std::string& line) const {
    std::string d{ "coord" };
    std::string f1 = line.substr(30, 8);
    std::string f2 = line.substr(38, 8);
    std::string f3 = line.substr(46, 8);
    std::array<std::string, 3> a{ { f1, f2, f3 } };
    for (const auto& i : a) {
      field_check(i, d);
    }
    auto fx{ std::stod(f1) };
    auto fy{ std::stod(f2) };
    auto fz{ std::stod(f3) };
    coords::Cartesian_Point cp(fx, fy, fz);
    return cp;
  }
};

/*!
\brief Simple parser class intended to be used as a functor.
\tparam T Numerical type intended for the storage of the coordinate data.
*/
template <typename T>
class Parser {
public:
  using Atom_type = Atom<Line, T>;

  Parser<T>(const std::string& str) { this->operator()(str); }

  void operator()(const std::string& str) { read_file(str); }

  /*!
  \brief Function for parsing the content of a Pdb file in a
  std::vector<Atom_type>.
  \param pdb_file Name of the Pdb file to be parsed.
  \return void; the atom_vec member is filled by this function.
  */
  void read_file(const std::string& pdb_file) {
    std::string file_line;
    std::cout << "Reading PDB file: " << pdb_file << "\n";
    std::ifstream input_pdb(pdb_file.c_str());
    if (!input_pdb.is_open())
      throw std::ios::failure(pdb_file + " cannot be opened.");
    while (std::getline(input_pdb, file_line)) {
      Line line_func;
      if (!line_func.line_check(file_line)) {
        continue;
      }
      auto atom = Atom<Line, T>(line_func, file_line);
      atom_vec.emplace_back(atom);
      if (input_pdb.bad())
        throw std::ios::failure("Error reading file " + pdb_file);
    }
  }

  /*!
  \brief Fragments the std::vector of atoms into a std::vector of residues,
  where each residue is itself a std::vector.
  \param res_vec std::vector of atoms.
  \return std::vector of residue vectors.
  */
  std::vector<std::vector<Atom_type>>
  create_resids(const std::vector<Atom_type>& res_vec) {
    auto last = res_vec.back().atom_serial;
    std::vector<std::vector<Atom_type>> result;
    std::vector<Atom_type> temp;
    for (auto& i : res_vec) {
      auto index = std::distance(result.begin(), result.end()) + 1;
      if (last == i.atom_serial) {
        temp.emplace_back(i);
        result.emplace_back(temp);
      } else if (index == i.res_seq) {
        temp.emplace_back(i);
      } else {
        result.emplace_back(temp);
        temp.erase(temp.begin(), temp.end());
        temp.emplace_back(i);
      }
    }
    return result;
  }

  /*!
  \brief Uses the std::vector of atoms to create a std::vector of index
  std::vectors. Each index std::vector represents all the atom serial numbers
  that belong to one residue.
  \param vec std::vector of atoms.
  \return std::vector of index std::vectors.
  */
  std::vector<std::vector<std::size_t>>
  create_resids_indices(const std::vector<Atom_type>& vec) {
    std::vector<std::vector<std::size_t>> result;
    auto resids = create_resids(vec);
    for (auto& res : resids) {
      std::vector<std::size_t> temp;
      for (auto& i : res) {
        temp.emplace_back(i.atom_serial);
      }
      result.emplace_back(temp);
    }
    return result;
  }

  /*!
  \brief Creates a coords::Representation_3D object from a std::vector of atoms.
  \param vec std::vector of atoms.
  \return coords::Representation_3D object.
  */
  coords::Representation_3D
  create_rep_3D(const std::vector<Atom_type>& vec) const {
    coords::Representation_3D cp_vec;
    for (auto& i : vec) {
      cp_vec.emplace_back(i.cp);
    }
    return cp_vec;
  }

  /*!
  \brief Uses a std::vector of atoms to create a std::vector of residues, where
  each residue is represented as coords::Representation_3D object.
  \param vec std::vector of atoms.
  \return std::vector of coords::Representation_3D objects.
  */
  std::vector<coords::Representation_3D>
  create_resids_rep_3D(const std::vector<Atom_type>& vec) {
    auto resid_vec = create_resids(vec);
    std::vector<coords::Representation_3D> result;
    for (auto& i : resid_vec) {
      auto temp = create_rep_3D(i);
      result.emplace_back(temp);
    }
    return result;
  }

private:
  unsigned int model_number_;

public:
  std::vector<Atom_type> atom_vec;
};
}

/*! @} End of ic_util group*/

#endif // cast_pdb_h_guard
