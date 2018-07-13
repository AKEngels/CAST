#ifndef COORDS_IO_XYZ_H
#define COORDS_IO_XYZ_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <set>
#include <vector>

#include "coords_io.h"
#include "energy.h"
#include "graph.h"
#include "ic_util.h"

using Graph_type = boost::adjacency_list<boost::vecS, boost::vecS,
                                         boost::undirectedS, ic_util::Node>;
namespace coords {
namespace input {
namespace formats {
using coords::float_type;
struct xyz : public format {

public:
  struct helper {

    static coords::Representation_3D
    ang_from_bohr(coords::Representation_3D const& rep3D) {
      coords::Representation_3D result;
      for (auto const& a : rep3D) {
        result.emplace_back(a * energy::bohr2ang);
      }
      return result;
    }

    template <typename Line, typename T>
    struct Atom {
    public:
      Atom() = default;

      Atom(const Line& func, const std::string& file)
          : cp{ std::move(func.cart_point(file)) }, element{ func.element(
                                                        file) } {};
      coords::Cartesian_Point cp = {};
      std::string element;

    private:
      template <typename _Line, typename _T>
      friend std::ostream& operator<<(std::ostream& os,
                                      Atom<_Line, _T> const& atom) {
        return os << std::setw(8) << atom.cp << ", " << std::setw(2)
                  << atom.element;
      }
    };

    struct Line {
    public:
      /*!
           \brief Checks whether the string to_check contains the string name.
           \param to_check String that is to be checked.
           \param name String whose existence in the string to_check shall be
         proven. \return std::runtime_error if to_check does not contain name.
           */
      void field_check(const std::string& to_check,
                       const std::string& name) const {
        if (std::all_of(to_check.begin(), to_check.end(), isspace)) {
          throw std::runtime_error("At least one field [ " + name +
                                   " ] in the XYZ file is empty.");
        }
      }
      /*!
         \brief Extracts the element information from the Pdb line.
         \param line The Pdb line.
         \return Element information.
         */
      std::string element(const std::string& line) const {
        std::string d{ "element" };
        std::string f = line.substr(0, 2);

        auto end = std::remove(f.begin(), f.end(), ' ');
        f.erase(end, f.end());
        return f;
      }

      /*!
          \brief Extracts the coordinate information from the xyz line.
          \tparam T Numerical type intended for the storage of the coordinate
          information.
          \param line The xyz line.
          \return 3-dimensional coordinate array.
          */
      std::array<coords::float_type, 3> coord(const std::string& line) const {
        std::string d{ "coord" };
        std::string f1 = line.substr(5, 20);
        std::string f2 = line.substr(23, 34);
        std::string f3 = line.substr(34, 48);
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
          \brief Extracts the coordinate information from the xyz line.
          \param line The xyz line.
          \return coords::Cartesian_Point object.
          */
      coords::Cartesian_Point cart_point(const std::string& line) const {
        auto c = coord(line);
        return coords::Cartesian_Point(c.at(0), c.at(1), c.at(2));
      }
    }; // Line

    template <typename T>
    class Parser {
    public:
      using Atom_type = Atom<Line, T>;

      Parser<T>(const std::string& str) { this->operator()(str); }

      void operator()(const std::string& str) { read_file(str); }

      /*!
          \brief Function for parsing the content of a xyz file in a
          std::vector<Atom_type>.
          \param xyz_file Name of the xyz file to be parsed.
          \return void; the atom_vec member is filled by this function.
          */
      void read_file(const std::string& xyz_file) {
        std::string file_line;
        std::cout << "Reading XYZ file: " << xyz_file << "\n";
        std::ifstream input_xyz(xyz_file);
        if (!input_xyz.is_open())
          throw std::ios::failure(xyz_file + " cannot be opened.");

        std::getline(input_xyz, file_line); // first line: number of atoms
        double N = std::stoi(file_line);

        std::getline(input_xyz, file_line); // discard second line (comment)

        while (std::getline(input_xyz, file_line)) {
          Line line_func;

          auto atom = Atom<Line, T>(line_func, file_line);
          atom_vec.emplace_back(atom);
          if (input_xyz.bad())
            throw std::ios::failure("Error reading file " + xyz_file);
        }
      }

      /*!
      \brief Fragments the std::vector of atoms into a std::vector of
      residues, where each residue is itself a std::vector.
      \param res_vec std::vector of atoms.
      \return std::vector of residue vectors.
      */
      template <typename AtomVec, typename IndicesVec>
      static std::vector<std::vector<Atom_type>>
      create_resids(const AtomVec& atom_vec, const IndicesVec& res_index_vec) {

        std::vector<Atom_type> molecule;
        std::vector<std::vector<Atom_type>> molecules;

        for (auto i = 0; i < res_index_vec.size(); ++i) {
          for (int atom : res_index_vec[i]) {
            molecule.emplace_back(atom_vec[atom]);
          }
          molecules.emplace_back(molecule);
          molecule.clear();
        }

        return molecules;
      }

      template <typename Vec>
      static std::vector<std::vector<Atom_type>>
      create_resids(const Vec& res_index_vec) {

        return create_resids(res_index_vec, atom_vec);
      }

      /*!
          \brief Uses the std::vector of atoms to create a std::vector of index
          std::vectors. Each index std::vector represents all the atom serial
         numbers that belong to one residue.
         \param vec std::vector of atoms.
          \return std::vector of index std::vectors.
          */

      static std::vector<std::vector<std::size_t>>
      create_resids_indices(Graph_type& graph) {

        std::vector<std::size_t> molecule_vec;

        std::vector<std::vector<std::size_t>> result, temp_result;
        std::vector<int> component(boost::num_vertices(graph));

        int num = connected_components(graph, &component[0]);

        if (num > 1) {
          result = find_connected_components(graph);
        }

        temp_result = find_amino_acids(graph);
        for (auto a = 0; a < temp_result.size(); ++a) {
          for (auto b : temp_result[a]) {
            molecule_vec.emplace_back(b);
          }
          result.emplace_back(molecule_vec);
          molecule_vec.clear();
        }

        return result;
      }

      // function to find connected subsystems in graph

      static std::vector<std::vector<size_t>>
      find_connected_components(Graph_type& graph) {
        std::vector<std::size_t> molecule_vec;

        std::vector<std::vector<std::size_t>> result;

        std::vector<int> component(boost::num_vertices(graph));
        auto num_components = connected_components(graph, &component[0]);

        for (auto j = 0; j < num_components; ++j) {

          for (auto i = 0; i < component.size(); ++i) {
            if (component[i] == j) {

              molecule_vec.emplace_back(i);
            }
          }
          result.emplace_back(molecule_vec);
          molecule_vec.clear();
        }
        return result;
      }
      // function to find amino acids in peptides
      static std::vector<std::vector<size_t>>
      find_amino_acids(Graph_type& graph) {

        auto number_components_graph = number_of_components(graph);
        
        std::vector<std::size_t> molecule_vec, chiral_C_vec, carboxy_C_vec,
            amine_vec;

        std::vector<std::vector<std::size_t>> result, temp_result;

        boost::graph_traits<Graph_type>::edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
          int u = boost::source(*ei, graph);
          int v = boost::target(*ei, graph);
          
          //the "terminal condition" is important for serin
          if (graph[u].element == "O" && graph[v].element == "C" &&
              terminal(graph, u)) {

            carboxy_C_vec = adjacents(graph, v);

            for (auto w : carboxy_C_vec) {

              if (graph[w].element == "C") {

                chiral_C_vec = adjacents(graph, w);

                for (auto x : chiral_C_vec) {
                  if (graph[x].element == "N") {

                    // backbone found
                    molecule_vec.emplace_back(u);
                    molecule_vec.emplace_back(v);
                    molecule_vec.emplace_back(w);
                    molecule_vec.emplace_back(x);

                    // terminal -OH group
                    for (auto z : carboxy_C_vec) {
                      if (graph[z].element == "O" && !terminal(graph, z)) {
                        molecule_vec.emplace_back(z);
                        std::vector<size_t> O_vec;
                        O_vec = adjacents(graph, z);
                        for (auto zo : O_vec) {
                          if (graph[zo].element == "H") {
                            molecule_vec.emplace_back(zo);
                          }
                        }
                      }
                    }

                    //finds hydrogen(s) from (terminal) amin group
                    amine_vec = adjacents(graph, x);

                    for (auto y : amine_vec) {

                      if (graph[y].element == "H") {
                        molecule_vec.emplace_back(y);
                      }
                    }

                    // finds sidechain
                    for (auto a : chiral_C_vec) {

                      //finds hydrogens bonded to the chiral C and therefore glycin
                      if (graph[a].element == "H") {
                        molecule_vec.emplace_back(a);
                      }

                      //to get the complete sidechain, the edge between the first carbon in the sidechain and the chiral carbon from backbone is removed.
                      //with "find_connected_components" we get the connected subsystems and search for the vector which contains the first carbon from sidechain
                      else if (a != v && a != x && graph[a].element != "H") {

                        boost::remove_edge(a, w, graph);

                        auto num_comp = number_of_components(graph);

                        // special case prolin
                        //to use the "trick" mentioned above, we have to remove a second edge (between nitrogen and the carbon from ring)
                        if (num_comp == number_components_graph) { 
                          for (auto y : amine_vec) {
                            if (graph[y].element == "C") {
                              auto Pro_C_vec = adjacents(graph, y);
                              if (Pro_C_vec.size() == 4) {

                                boost::remove_edge(y, x, graph);
                                temp_result = find_connected_components(graph);
                                for (auto b = 0; b < temp_result.size(); ++b) {
                                  for (auto c : temp_result[b]) {
                                    if (c == a) {
                                      for (auto d : temp_result[b]) {
                                        molecule_vec.emplace_back(d);
                                      }
                                    }
                                  }
                                }
                                temp_result.clear();

                                boost::add_edge(y, x, graph);
                              }
                            }
                          }
                        } // prolin case close

                        else {

                          temp_result = find_connected_components(graph);
                          for (auto b = 0; b < temp_result.size(); ++b) {
                            for (auto c : temp_result[b]) {
                              if (c == a) {
                                for (auto d : temp_result[b]) {
                                  molecule_vec.emplace_back(d);
                                }
                              }
                            }
                          }
                          temp_result.clear();

                          
                        }
                        boost::add_edge(a, w, graph);
                        auto comp_graph = number_of_components(graph);

                        if (comp_graph != number_components_graph) {
                          std::cout << "ERROR: graph is ruined! Amino acids "
                                       "may be wrong.\n";
                        }
                      }
                    }
                    result.emplace_back(molecule_vec);
                    molecule_vec.clear();
                  }
                }
                chiral_C_vec.clear();
                carboxy_C_vec.clear();
                amine_vec.clear();
              }
            }
          }
          // the "terminal condition" is important for serin
          if (graph[v].element == "O" && graph[u].element == "C" &&
              terminal(graph, v)) {

            carboxy_C_vec = adjacents(graph, u);

            for (auto w : carboxy_C_vec) {

              if (graph[w].element == "C") {

                chiral_C_vec = adjacents(graph, w);

                for (auto x : chiral_C_vec) {
                  if (graph[x].element == "N") {

                    // backbone found
                    molecule_vec.emplace_back(u);
                    molecule_vec.emplace_back(v);
                    molecule_vec.emplace_back(w);
                    molecule_vec.emplace_back(x);

                    // terminal -OH group
                    for (auto z : carboxy_C_vec) {
                      if (graph[z].element == "O" && !terminal(graph, z)) {
                        molecule_vec.emplace_back(z);
                        std::vector<size_t> O_vec;
                        O_vec = adjacents(graph, z);
                        for (auto zo : O_vec) {
                          if (graph[zo].element == "H") {
                            molecule_vec.emplace_back(zo);
                          }
                        }
                      }
                    }

                    //finds hydrogen(s) from (terminal) amin group
                    amine_vec = adjacents(graph, x);
                                        
                    for (auto y : amine_vec) {

                      if (graph[y].element == "H") {
                        molecule_vec.emplace_back(y);
                      }
                    }

                    //finds sidechain
                    for (auto a : chiral_C_vec) {

                      // finds hydrogens bonded to the chiral C and therefore glycin
                      if (graph[a].element == "H") {
                        molecule_vec.emplace_back(a);
                      }

                      // to get the complete sidechain, the edge between the
                      // first carbon in the sidechain and the chiral carbon from
                      // backbone is removed. With "find_connected_components" we
                      // get the connected subsystems and search for the vector
                      // which contains the first carbon from sidechain
                      else if (a != u && a != x && graph[a].element != "H") {

                        boost::remove_edge(a, w, graph);

                        auto num_comp = number_of_components(graph);

                        // special case prolin
                        // to use the "trick" mentioned above, we have to remove a second edge (between nitrogen and the carbon from ring)
                        if (num_comp == number_components_graph) { 
                          for (auto y : amine_vec) {
                            if (graph[y].element == "C") {
                              auto Pro_C_vec = adjacents(graph, y);
                              if (Pro_C_vec.size() == 4) {

                                boost::remove_edge(y, x, graph);
                                temp_result = find_connected_components(graph);
                                for (auto b = 0; b < temp_result.size(); ++b) {
                                  for (auto c : temp_result[b]) {
                                    if (c == a) {
                                      for (auto d : temp_result[b]) {
                                        molecule_vec.emplace_back(d);
                                      }
                                    }
                                  }
                                }
                                temp_result.clear();

                                boost::add_edge(y, x, graph);
                              }
                            }
                          }
                        } // prolin case closed

                        else {

                          temp_result = find_connected_components(graph);
                          for (auto b = 0; b < temp_result.size(); ++b) {
                            for (auto c : temp_result[b]) {
                              if (c == a) {
                                for (auto d : temp_result[b]) {
                                  molecule_vec.emplace_back(d);
                                }
                              }
                            }
                          }
                          temp_result.clear();

                          
                        }
                        boost::add_edge(a, w, graph);
                        auto comp_graph = number_of_components(graph);

                        if (comp_graph != number_components_graph) {
                          std::cout << "ERROR: graph is ruined! Amino acids "
                                       "may be wrong\n";
                        }
                      }
                    }

                    result.emplace_back(molecule_vec);
                    molecule_vec.clear();
                  }
                }
                chiral_C_vec.clear();
                carboxy_C_vec.clear();
                amine_vec.clear();
              }
            }
          }
        }
        return result;
      }

      static bool terminal(const Graph_type& graph, size_t const& x) {
        std::vector<size_t> number_adjacents;
        Graph_type::adjacency_iterator ai, ai_end;
        for (tie(ai, ai_end) = boost::adjacent_vertices(x, graph); ai != ai_end;
             ++ai) {

          number_adjacents.emplace_back(*ai);
        }
        return (number_adjacents.size() == 1);
      }

      static std::vector<size_t> adjacents(const Graph_type& graph,
                                           const size_t& x) {

        std::vector<size_t> adjacents;
        Graph_type::adjacency_iterator ai, ai_end;
        for (tie(ai, ai_end) = boost::adjacent_vertices(x, graph); ai != ai_end;
             ++ai) {

          adjacents.emplace_back(*ai);
        }
        return adjacents;
      }

      static int number_of_components(const Graph_type& graph) {

        std::vector<int> component(boost::num_vertices(graph));

        int num = connected_components(graph, &component[0]);
        return num;
      }
      /*!
      \brief Creates a coords::Representation_3D object from a std::vector of
      atoms. \param vec std::vector of atoms. \return coords::Representation_3D
      object.  */
    private:
      static coords::Representation_3D
      create_rep_3D_impl(const std::vector<Atom_type>& vec) {
        coords::Representation_3D cp_vec;
        for (auto& i : vec) {
          cp_vec.emplace_back(i.cp);
        }
        return cp_vec;
      }

    public:
      static coords::Representation_3D
      create_rep_3D(const std::vector<Atom_type>& vec) {
        return create_rep_3D_impl(vec);
      }

      coords::Representation_3D create_rep_3D() const {
        return create_rep_3D_impl(atom_vec);
      }

    private:
      static coords::Representation_3D
      create_rep_3D_bohr_impl(const std::vector<Atom_type>& vec) {
        auto rep3D = create_rep_3D(vec);
        for (auto&& coord : rep3D) {
          coord /= energy::bohr2ang;
        }
        return rep3D;
      }

    public:
      static coords::Representation_3D
      create_rep_3D_bohr(const std::vector<Atom_type>& vec) {
        return create_rep_3D_bohr_impl(vec);
      }
      coords::Representation_3D create_rep_3D_bohr() const {
        return create_rep_3D_bohr_impl(atom_vec);
      }
      /*!
         \brief Uses a std::vector of atoms to create a std::vector of residues,
         where each residue is represented as coords::Representation_3D object.
         \param vec std::vector of atoms.
         \return std::vector of coords::Representation_3D objects.
         */
    private:
      template <typename Vec, typename Func, typename Vec_res>
      static std::vector<coords::Representation_3D>
      create_resids_rep_3D_impl(Vec&& vec, Func creator,
                                Vec_res const& res_index_vec) {
        auto resid_vec = create_resids(std::forward<Vec>(vec), res_index_vec);
        std::vector<coords::Representation_3D> result;
        for (auto& i : resid_vec) {
          auto temp = creator(i);
          result.emplace_back(temp);
        }
        return result;
      }

    public:
      template <typename Vec>
      static std::vector<coords::Representation_3D>
      create_resids_rep_3D(Vec&& vec) {
        return create_resids_rep_3D_impl(std::forward<Vec>(vec),
                                         create_rep_3D_impl);
      }

      std::vector<coords::Representation_3D> create_resids_rep_3D() const {
        return create_resids_rep_3D(atom_vec);
      }

      template <typename Vec, typename Vec_res>
      static std::vector<coords::Representation_3D>
      create_resids_rep_3D_bohr(Vec&& vec, const Vec_res& res_index_vec) {
        return create_resids_rep_3D_impl(
            std::forward<Vec>(vec), create_rep_3D_bohr_impl, res_index_vec);
      }

      template <typename Vec>
      std::vector<coords::Representation_3D>
      create_resids_rep_3D_bohr(const Vec& res_index_vec) const {
        return create_resids_rep_3D_bohr(atom_vec, res_index_vec);
      }

      /*!
          \brief Uses a std::vector of xyz::Atom to form a std::vector of
         strings containing the element symbols
         \param vec std::vector of xyz::Atom
         \return std::vector of std::string
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
      create_cooord_atoms(std::vector<Atom_type> const& vec) {
        coords::Atoms atoms;
        for (auto const& atom : vec) {
          coords::Atom tmp_atom(atom.element);
          atoms.add(std::move(tmp_atom));
        }
        return atoms;
      }

      coords::Atoms create_cooord_atoms() const {
        return create_cooord_atoms(atom_vec);
      }

    private:
      unsigned int model_number_;

    public:
      std::vector<Atom_type> atom_vec;

    }; // Parser

    static inline void make_bonds(Atoms&,
                                  std::vector<std::pair<int, int>> const&);

  }; // helper

public:
  inline Coordinates read(std::string const&) override;
  std::shared_ptr<coords::input::formats::xyz::helper::Parser<float_type>>
      parser;
}; // namespace formats
} // namespace formats
} // namespace input
} // namespace coords

void coords::input::formats::xyz::helper::make_bonds(
    coords::Atoms& atoms, std::vector<std::pair<int, int>> const& bonds) {
  for (auto const& bond : bonds) {
    auto i = bond.first - 1u, j = bond.second - 1u;
    auto& atom1 = atoms.atom(i);
    auto& atom2 = atoms.atom(j);
    atom1.bind_to(j);
    atom2.bind_to(i);
  }
}
coords::Coordinates coords::input::formats::xyz::read(std::string const& file) {
  if ((Config::get().general.energy_interface ==
       config::interface_types::T::AMBER) ||
      (Config::get().general.energy_interface ==
       config::interface_types::T::AMOEBA) ||
      (Config::get().general.energy_interface ==
       config::interface_types::T::CHARMM22) ||
      (Config::get().general.energy_interface ==
       config::interface_types::T::OPLSAA)) {
    std::cout << "ERROR: It is not possible to use XYZ files with a forcefield "
                 "interface because no atom types are assigned!\n";
    if (Config::get().general.task == config::tasks::WRITE_TINKER) {
      std::cout << "Yes, I know you just want to write a tinkerstructure and "
                   "you don't need any energies. But it doesn't work like "
                   "this. So just use GAUSSIAN or MOPAC as energy interface "
                   "and all will be fine (even if you don't have access to any "
                   "of these programmes).\n";
    }
    std::exit(0);
  }

  parser =
      std::make_shared<coords::input::formats::xyz::helper::Parser<float_type>>(
          file);

  auto atoms = parser->create_cooord_atoms();
  auto rep3D = parser->create_rep_3D();
  input_ensemble.emplace_back(rep3D);

  coords::input::formats::xyz::helper::make_bonds(
      atoms, ic_util::bonds(parser->create_element_vec(), rep3D));

  Coordinates coord_object;

  PES_Point pes(input_ensemble[0u]);

  coord_object.init_swap_in(atoms,
                            pes); // fill atoms and positions into coord_object

  for (auto& p :
       input_ensemble) // do some important stuff (see coords_io_AMBER.cc)
  {
    p.gradient.cartesian.resize(p.structure.cartesian.size());
    coord_object.set_xyz(p.structure.cartesian);
    coord_object.to_internal_light();
    p = coord_object.pes();
  }

  return coord_object;
}

#endif