#ifndef COORDS_IO_PDB_H
#define COORDS_IO_PDB_H

#include <set>
#include <vector>

#include "coords_io.h"
#include "energy.h"
#include "graph.h"
#include "ic_util.h"

namespace coords {
namespace input {
namespace formats {
using coords::float_type;
struct pdb : public format {
public:
  struct helper : public helper_base {

    /*!
    \brief Scoped enumeration type with the relevant Pdb entries as possible
    values.
    */
    enum class kindOfEntryInLine : int {
      ATOM,
      HETATOM,
      NUMMDL,
      MODEL,
      ENDMDL,
      TER,
      END,
      UNKNOWN
    };
    enum class terminals : int { none = 0, C, N };
    inline static std::set<std::string> const& residueNames() {
      static const std::set<std::string> ret{
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "CYX", "CYM", "CYP", "HID", "HIE", "HIP"
      };
      return ret;
    }
    /*!
    \brief Helper function in order to retrieve information about the value of
    the Record type. \param v Record object. \return String representation of
    the value of the Record object.
    */
    inline static std::map<kindOfEntryInLine, std::string> const&
    kindOfEntryInLineStringMap() {
      static const std::map<kindOfEntryInLine, std::string> ret{
        { kindOfEntryInLine::ATOM, "ATOM" },
        { kindOfEntryInLine::HETATOM, "HETATOM" },
        { kindOfEntryInLine::NUMMDL, "NUMMDL" },
        { kindOfEntryInLine::MODEL, "MODEL" },
        { kindOfEntryInLine::ENDMDL, "ENDMDL" },
        { kindOfEntryInLine::TER, "TER" },
        { kindOfEntryInLine::END, "TER" }
      };
      return ret;
    }
    inline static std::map<std::string, kindOfEntryInLine> const&
    stringKindOfEntryInLineMap() {
      static const std::map<std::string, kindOfEntryInLine> ret{
        { "ATOM", kindOfEntryInLine::ATOM },
        { "HETATOM", kindOfEntryInLine::HETATOM },
        { "HETATM", kindOfEntryInLine::HETATOM },
        { "NUMMDL", kindOfEntryInLine::NUMMDL },
        { "MODEL", kindOfEntryInLine::MODEL },
        { "ENDMDL", kindOfEntryInLine::ENDMDL },
        { "TER", kindOfEntryInLine::TER },
        { "TER", kindOfEntryInLine::END }
      };
      return ret;
    }

    static std::string recordString(kindOfEntryInLine v) {
      return ic_util::getValueByKeyOrDefault(kindOfEntryInLineStringMap(), v,
                                             "Unknown Record type.");
    }
    /**finds and returns atom type for atoms in protein sidechain (OPLSAA
    forcefield)
    @param atom_name: atom name from pdb file
    @param res_name: residue name from pdb file*/
    static inline int find_at_sidechain(std::string const& atom_name,
                                        std::string const& res_name);
    /**function that assigns atom types (oplsaa) to atoms of protein backbone
    (they are not suitable for force field calucations)
    @param atom_name: atom name from pdb file
    @param res_name: residue name from pdb file
    @param terminal: is residue N-terminal, C-terminal or not?*/
    static inline int find_energy_type(std::string atom_name,
                                       std::string res_name,
                                       terminals terminal);
    /**finds element symbol and energy type
    @param atom_name: atom name from pdb file
    @param res_name: residue name from pdb file
    returns element symbol*/
    static inline std::string find_element_symbol(std::string atom_name,
                                                  std::string res_name);
    /*!
    \brief Overloaded output operator for the Record enumeration type.
    \param os Output stream.
    \param c Record object.
    \return String representation for the Record object piped to standard
    output.
    */
    friend std::ostream& operator<<(std::ostream& os, kindOfEntryInLine c) {
      return os << recordString(c);
    }

    /*!
    \brief Simple class representing the Atom concept.
    \details The members are initialized using various helper functions.
    \tparam Line Helper type representing an arbitrary line of the Pdb file.
    \tparam T Type intended to hold the coordinate data.
    */
    template <typename T>
    struct Atom {
    public:
      Atom() = default;
      template <typename Line>
      Atom(const Line& func, const std::string& file)
          : rec_name{ func.lineType(file) },
            atom_serial{ func.atom_serial(file) }, atom_name{ func.atom_name(
                                                       file) },
            alt_loc{ func.alt_loc(file) }, res_name{ func.res_name(file) },
            chain_id{ func.chain_id(file) }, res_seq{ func.res_seq(file) },
            insertion_code{ func.insertion_code(file) },
            cp{ std::move(func.cart_point(file)) }, element{ func.element(
                                                        file) } {}

      kindOfEntryInLine rec_name;
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
      template <typename _T>
      friend std::ostream& operator<<(std::ostream& os, Atom<_T> const& atom) {
        return os << std::setw(7) << atom.rec_name << ", " << std::setw(5)
                  << atom.atom_serial << ", " << std::setw(4) << atom.atom_name
                  << ", " << std::setw(9) << atom.alt_loc << ", "
                  << std::setw(3) << atom.res_name << ", " << std::setw(1)
                  << atom.chain_id << ", " << std::setw(4) << atom.res_seq
                  << ", " << std::setw(8) << atom.insertion_code << ", "
                  << std::setw(8) << atom.cp << ", " << std::setw(2)
                  << atom.element;
      }
    };

    /*!
    \brief Collection of helper functions for parsing the Pdb lines. For ease of
    handling the functions are contained within one class.
    */
    struct Line {
    public:
      /*static constexpr std::array<std::string, 26> RESIDUE_NAMES{ { "ALA",
        "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
        "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "CYX",
        "CYM", "CYP", "HID", "HIE", "HIP" } };*/
      /*!
      \brief Function for recognizing the Record value of a Pdb line.
      \param line String representation of the Pdb line.
      \return Specific value of the Record type.
      */
    private:
      std::string getFirstWordOfLine(std::string const& line) const {
        std::stringstream lineStream{ line };
        std::string firstWordOfLine;
        lineStream >> firstWordOfLine;
        return firstWordOfLine;
      }

    public:
      kindOfEntryInLine lineType(std::string const& line) const {
        return ic_util::getValueByKeyOrDefault(
            stringKindOfEntryInLineMap(),
            removeBlanksFromString(getFirstWordOfLine(line)),
            kindOfEntryInLine::UNKNOWN);
      }

      /*!
      \brief Function for checking whether the Pdb line is an ATOM or HETATOM
      entry. Only these are relevant. \param line String representation of the
      Pdb line. \return True, if the line is an ATOM or HETATOM entry; otherwise
      false.
      */
      bool line_check(const std::string& line) const {
        auto type = lineType(line);
        return type == kindOfEntryInLine::ATOM ||
               type == kindOfEntryInLine::HETATOM;
      }

      /*!
      \brief Checks whether the string to_check contains the string name.
      \param to_check String that is to be checked.
      \param name String whose existence in the string to_check shall be proven.
      \return std::runtime_error if to_check does not contain name.
      */
      void field_check(const std::string& to_check,
                       const std::string& name) const {
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
        try {
          field_check(f, d);
        } catch (std::runtime_error const& e) {
          return 'E';
        }
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
        try {
          field_check(f, d);
        } catch (std::runtime_error const& e) {
          return find_element_symbol(atom_name(line), res_name(line));
        }
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
      std::array<coords::float_type, 3> coord(const std::string& line) const {
        std::string d{ "coord" };
        std::string f1 = line.substr(30, 8);
        std::string f2 = line.substr(38, 8);
        std::string f3 = line.substr(46, 8);
        std::array<std::string, 3> a{ { f1, f2, f3 } };
        std::array<coords::float_type, 3> c = {};
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
        auto c = coord(line);
        return coords::Cartesian_Point(c.at(0), c.at(1), c.at(2));
      }
    };
    /*!
    \brief Simple parser class intended to be used as a functor.
    \tparam T Numerical type intended for the storage of the coordinate data.
    */
    template <typename T>
    class Parser {
    public:
      using Atom_type = Atom<T>;

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
        std::ifstream input_pdb(pdb_file);
        if (!input_pdb.is_open())
          throw std::ios::failure(pdb_file + " cannot be opened.");
        while (std::getline(input_pdb, file_line)) {
          Line line_func;
          if (!line_func.line_check(file_line)) {
            continue;
          }
          auto atom = Atom<T>(line_func, file_line);
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
      static std::vector<std::vector<Atom_type>>
      create_resids(const std::vector<Atom_type>& res_vec) {
        std::vector<std::vector<Atom_type>> result;
        std::vector<Atom_type> temp;
        for (auto const& atom_i : res_vec) {
          if (result.size() + 1 != atom_i.res_seq) {
            result.emplace_back(std::move(temp));
          }
          temp.emplace_back(atom_i);
        }
        result.emplace_back(std::move(temp));
        return result;
      }

      std::vector<std::vector<Atom_type>> create_resids() const {
        return create_resids(atom_vec);
      }

      /*!
      \brief Uses the std::vector of atoms to create a std::vector of index
      std::vectors. Each index std::vector represents all the atom serial
      numbers that belong to one residue. \param vec std::vector of atoms.
      \return std::vector of index std::vectors.
      */
      static std::vector<std::vector<std::size_t>>
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

      std::vector<std::vector<std::size_t>> create_resids_indices() const {
        return create_resids_indices(atom_vec);
      }

      static coords::Representation_3D
        create_rep_3D(const std::vector<Atom_type>& vec) {
        return create_rep_3D_impl(vec);
      }

      coords::Representation_3D create_rep_3D() const {
        return create_rep_3D_impl(atom_vec);
      }

      static coords::Representation_3D
        create_rep_3D_bohr(const std::vector<Atom_type>& vec) {
        return create_rep_3D_bohr_impl(vec);
      }
      coords::Representation_3D create_rep_3D_bohr() const {
        return create_rep_3D_bohr_impl(atom_vec);
      }

      template <typename Vec>
      static std::vector<coords::Representation_3D>
        create_resids_rep_3D(Vec&& vec) {
        return create_resids_rep_3D_impl(std::forward<Vec>(vec),
          create_rep_3D_impl);
      }

      std::vector<coords::Representation_3D> create_resids_rep_3D() const {
        return create_resids_rep_3D(atom_vec);
      }

      template <typename Vec>
      static std::vector<coords::Representation_3D>
        create_resids_rep_3D_bohr(Vec&& vec) {
        return create_resids_rep_3D_impl(std::forward<Vec>(vec),
          create_rep_3D_bohr_impl);
      }

      std::vector<coords::Representation_3D> create_resids_rep_3D_bohr() const {
        return create_resids_rep_3D_bohr(atom_vec);
      }

      /*!
      \brief Uses a std::vector of pdb::Atom to form a std::vector of strings
      containing the element symbols \param vec std::vector of pdb::Atom \return
      std::vector of std::string
      */
      static std::vector<std::string>
        create_element_vec(std::vector<Atom_type> const& vec) {
        std::vector<std::string> result;
        for (auto const& atom : vec) {
          result.emplace_back(atom.element);
        }
        return result;
      }
      std::vector<std::string> create_element_vec() const {
        return create_element_vec(atom_vec);
      }

      /*!
      \brief Uses a std::vector of pdb::Atom to form a std::vector of
      coords::Atom \param vec std::vector of pdb::Atom \return std::vector of
      coords::Atom
      */
      static coords::Atoms
        create_coord_atoms(std::vector<Atom_type> const& vec) 
      {
        coords::Atoms atoms;
        for (auto const& atom : vec) {
          coords::Atom tmp_atom(atom.element);
          //tmp_atom.set_residue(atom.res_name);
          //tmp_atom.set_res_id(atom.chain_id);
          //tmp_atom.set_pdb_atom_name(atom.atom_name);
          atoms.add(std::move(tmp_atom));
        }
        return atoms;
      }

      coords::Atoms create_coord_atoms() const {
        return create_coord_atoms(atom_vec);
      }

      std::vector<Atom_type> atom_vec;
      

    private:

      
      /*!
      \brief Creates a coords::Representation_3D object from a std::vector of
      atoms. \param vec std::vector of atoms. \return coords::Representation_3D
      object.
      */
      static coords::Representation_3D
      create_rep_3D_impl(const std::vector<Atom_type>& vec) {
        coords::Representation_3D cp_vec;
        for (auto& i : vec) {
          cp_vec.emplace_back(i.cp);
        }
        return cp_vec;
      }

      static coords::Representation_3D
      create_rep_3D_bohr_impl(const std::vector<Atom_type>& vec) {
        auto rep3D = create_rep_3D(vec);
        for (auto&& coord : rep3D) {
          coord /= energy::bohr2ang;
        }
        return rep3D;
      }

      /*!
      \brief Uses a std::vector of atoms to create a std::vector of residues,
      where each residue is represented as coords::Representation_3D object.
      \param vec std::vector of atoms.
      \return std::vector of coords::Representation_3D objects.
      */
      template <typename Vec, typename Func>
      static std::vector<coords::Representation_3D>
      create_resids_rep_3D_impl(Vec&& vec, Func creator) {
        auto resid_vec = create_resids(std::forward<Vec>(vec));
        std::vector<coords::Representation_3D> result;
        for (auto& i : resid_vec) {
          auto temp = creator(i);
          result.emplace_back(temp);
        }
        return result;
      }

      unsigned int model_number_;
      
    };
    /**struct that contains information about a residue*/
    struct residue {
      /**name of the residue*/
      std::string res_name;
      /**atoms in the residue*/
      std::vector<coords::Atom> atoms;
      /**is an amino acid terminal or not
      possible values: no (not terminal), C (C terminal), N (N terminal)*/
      std::string terminal;
    };
    static inline void
    make_bonds(Atoms&, std::vector<std::pair<std::size_t, std::size_t>> const&);
    static inline void
    set_energy_type(Atoms&, std::vector<std::vector<std::size_t>> const&);
  };

public:
  Coordinates read(std::string) override;
  std::shared_ptr<coords::input::formats::pdb::helper::Parser<float_type>>
      parser;
};
} // namespace formats
} // namespace input
} // namespace coords



void coords::input::formats::pdb::helper::make_bonds(
    coords::Atoms& atoms,
    std::vector<std::pair<std::size_t, std::size_t>> const& bonds) {
  for (auto const& bond : bonds) {
    auto is_ion = [](std::string const& s) {
      return s.substr(s.size() - 1, 1) == "+" ||
             s.substr(s.size() - 1, 1) == "-";
    };
    // do not bond ions
    auto i = bond.first, j = bond.second;
    auto& atom1 = atoms.atom(i);
    auto& atom2 = atoms.atom(j);
    atom1.bind_to(j);
    atom2.bind_to(i);
  }
}

void coords::input::formats::pdb::helper::set_energy_type(
    coords::Atoms& atoms,
    std::vector<std::vector<std::size_t>> const& indices) 
{
  for (auto i = 0u; i < indices.size(); ++i) {
    auto const& res_ind = indices[i];
    for (auto const& ind : res_ind) {
      auto& atom = atoms.atom(ind - 1u);
      int energy_type = 0;
      if (energy_type == 0 && Config::get().general.energy_interface ==
                                  config::interface_types::T::OPLSAA) {
        if (Config::get().general.task == config::tasks::WRITE_TINKER) {
          throw std::runtime_error(
              "Yes, I know you just want to write a tinkerstructure and "
              "you don't need any energies. But it doesn't work like this. "
              "So just use GAUSSIAN or MOPAC as energy interface and all "
              "will be fine (even if you don't have access to any of these "
              "programmes).\n");
        }
        throw std::runtime_error("Assigment of atom types failed. Please use "
                                 "another energy interface.\n");
      }
      atom.set_energy_type(energy_type);
    }
  }
}

inline coords::Coordinates coords::input::formats::pdb::read(std::string file_name)
{
  if ((Config::get().general.energy_interface ==
       config::interface_types::T::AMBER) ||
      (Config::get().general.energy_interface ==
       config::interface_types::T::AMOEBA) ||
      (Config::get().general.energy_interface ==
       config::interface_types::T::CHARMM22)) {
    if (Config::get().general.task == config::tasks::WRITE_TINKER) {
      throw std::runtime_error("There is no way to use a PDB file in "
                               "order to write a Tinker File"
                               " because no force field parameters are "
                               "known. Please use another "
                               "input file.\n");
    }
    throw std::runtime_error(
        "ERROR: It is not possible to use PDB files with that "
        "interface "
        "because wrong atom types are assigned!\n");
  }

  parser =
      std::make_shared<coords::input::formats::pdb::helper::Parser<float_type>>(
          file_name);

  auto atoms = parser->create_coord_atoms();
  auto rep3D = parser->create_rep_3D();
  input_ensemble.emplace_back(rep3D);

  coords::input::formats::pdb::helper::make_bonds(
      atoms, ic_util::bonds(parser->create_element_vec(), rep3D));
  coords::input::formats::pdb::helper::set_energy_type(
      atoms, parser->create_resids_indices());

  Coordinates result;
  PES_Point pes(input_ensemble.at(0u));
  result.init_swap_in(atoms, pes);

  for (auto& p : input_ensemble) {
    p.gradient.cartesian.resize(p.structure.cartesian.size());
    result.set_xyz(p.structure.cartesian);
    result.to_internal_light();
    p = result.pes();
  }
  return result;
};

namespace {
bool elementExistsInResidueNames(std::string name) {
  auto const& residues = coords::input::formats::pdb::helper::residueNames();
  return residues.find(name) != residues.end();
}
} // namespace


std::string
coords::input::formats::pdb::helper::find_element_symbol(std::string atom_name,
                                                         std::string res_name) {
  std::string element;
  if (elementExistsInResidueNames(res_name)) // protein
  {
    element = atom_name.substr(0, 1);
  } else if (atom_name.substr(atom_name.size() - 1, 1) == "+" ||
             atom_name.substr(atom_name.size() - 1, 1) == "-") // ion
  {
    element = atom_name.substr(0, atom_name.size() - 1);
  } else if (res_name == "WAT" || res_name == "LIG") // water or ligand
  {
    element = "";
    for (auto s : atom_name) // every sign in atom name
    {
      if (!isdigit(s))
        element += s;
    }
  } else // unknown residue
  {
    std::cout << "WARNING: unknown residue: " << res_name
              << ". Guessing element symbol.\n";
    element = "";
    for (auto s : atom_name) // every sign in atom name
    {
      if (!isdigit(s))
        element += s;
    }
  }
  return element;
}
#endif
