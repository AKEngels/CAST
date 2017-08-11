#ifndef cast_pdb_h_guard
#define cast_pdb_h_guard

#pragma once

#include "coords_rep.h"

#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

namespace Pdb {
/// header-only (i.e. no separate implementation file provided)

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

inline std::ostream& operator<<(std::ostream& os, Record c) {
  return os << Record_string(c);
}

template <typename Line, typename T>
struct Atom {
public:
  Atom() : {};

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
  template <typename Line, typename T>
  friend std::ostream& operator<<(std::ostream& os, const Atom<Line, T>& atom);
};

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

struct Line {
public:
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

  bool line_check(const std::string& line) const {
    auto type = line_type(line);
    if (type == Record::ATOM || type == Record::HETATOM) {
      return true;
    } else {
      return false;
    }
  }

  void field_check(const std::string& to_check, const std::string& name) const {
    if (std::all_of(to_check.begin(), to_check.end(), isspace)) {
      throw std::runtime_error("At least one field [ " + name +
                               " ] in the PDB file is empty.");
    }
  }

  unsigned int model_serial(const std::string& line) const {
    std::string d{ "model serial" };
    std::string f = line.substr(10, 4);
    field_check(f, d);
    return std::stoi(f);
  }

  unsigned int atom_serial(const std::string& line) const {
    std::string d{ "atom serial" };
    std::string f = line.substr(6, 5);
    field_check(f, d);
    return std::stoi(f);
  }

  std::string atom_name(const std::string& line) const {
    std::string d{ "atom name" };
    std::string f = line.substr(12, 4);
    field_check(f, d);
    return f;
  }

  std::string alt_loc(const std::string& line) const {
    std::string d{ "alternate location" };
    std::string f = line.substr(16, 1);
    if (std::all_of(f.begin(), f.end(), isspace)) {
      return "No altLoc";
    }
    return f;
  }

  std::string res_name(const std::string& line) const {
    std::string d{ "residue name" };
    std::string f = line.substr(17, 3);
    field_check(f, d);
    return f;
  }

  char chain_id(const std::string& line) const {
    std::string d{ "chain ID" };
    std::string f = line.substr(21, 1);
    field_check(f, d);
    return f.at(0);
  }

  unsigned int res_seq(const std::string& line) const {
    std::string d{ "residue sequence number" };
    std::string f = line.substr(22, 4);
    field_check(f, d);
    return std::stoi(f);
  }

  std::string insertion_code(const std::string& line) const {
    std::string d{ "insertion code" };
    std::string f = line.substr(26, 1);
    if (std::all_of(f.begin(), f.end(), isspace)) {
      return "No iCode";
    }
    return f;
  }

  std::string element(const std::string& line) const {
    std::string d{ "element" };
    std::string f = line.substr(76, 2);
    field_check(f, d);
    auto end = std::remove(f.begin(), f.end(), ' ');
    f.erase(end, f.end());
    return f;
  }

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

template <typename T>
class Parser {
public:
  using Atom_type = Atom<Line, T>;

  Parser<T>(const std::string& str) { this->operator()(str); }

  void operator()(const std::string& str) { read_file(str); }

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

  coords::Representation_3D
  create_rep_3D(const std::vector<Atom_type>& vec) const {
    coords::Representation_3D cp_vec;
    for (auto& i : vec) {
      cp_vec.emplace_back(i.cp);
    }
    return cp_vec;
  }

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

#endif // cast_pdb_h_guard