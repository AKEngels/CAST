#ifndef COORDS_IO_XYZ_H
#define COORDS_IO_XYZ_H

#include <boost/cstdlib.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <set>
#include <vector>

#include "coords_io.h"
#include "energy.h"
#include "graph.h"
#include "ic_util.h"

using Graph_type = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::undirectedS, ic_util::Node,
    boost::property<boost::edge_color_t, boost::default_color_type>>;
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

      /* this is not working for molecules where rings share edges
      only conjugated molecules (sp2)
         source is the last atom visited, target the first one*/
      struct detect_loops : public boost::dfs_visitor<> {

        template <typename Function>
        detect_loops(Function callback)
            : boost::dfs_visitor<>(), callback(callback) {}

        std::function<void(std::vector<int> const&)> callback;

        template <class Edge, class Graph>
        void back_edge(Edge e, const Graph& g) {

          // std::cout << source(e, g) << " -- " << target(e, g) << "\n";

          std::vector<int> back_edge_indices;
          back_edge_indices.emplace_back(source(e, g));
          back_edge_indices.emplace_back(target(e, g));

          callback(back_edge_indices);
        }
      };

      // static std::vector<std::vector<int>> get_back_edges(Graph_type const&
      // graph) {

      //  using Vertex_type =
      //  boost::graph_traits<Graph_type>::vertex_descriptor;

      //  std::vector<std::vector<int>> back_edges;
      //  auto back_edges_callback = [&back_edges](std::vector<int> const& edge)
      //  {
      //    back_edges.emplace_back(edge);
      //  };

      //  detect_loops vis(back_edges_callback);
      //  undirected_dfs(graph,
      //                 boost::root_vertex(Vertex_type(0))
      //                     .visitor(vis)
      //                     .edge_color_map(get(boost::edge_color, graph)));
      //  return back_edges;

      //}

      static std::vector<std::vector<size_t>>
      find_aromatic_subunit(Graph_type graph) {

        std::vector<std::vector<size_t>> aromatic_subunits;

        using vertex_type = boost::graph_traits<Graph_type>::vertex_descriptor;

        std::vector<std::vector<int>> back_edges;
        auto back_edges_callback = [&back_edges](std::vector<int> const& edge) {
          back_edges.emplace_back(edge);
        };

        detect_loops vis(back_edges_callback);
        undirected_dfs(graph,
                       boost::root_vertex(vertex_type(0))
                           .visitor(vis)
                           .edge_color_map(get(boost::edge_color, graph)));
        /*       auto back_edges = get_back_edges(graph);*/

        if (back_edges.size() == 0)
          return aromatic_subunits;

        else {
          std::vector<std::vector<size_t>> cycles_in_molecule;
          for (auto backedge : back_edges) {
            auto one_cycle =
                find_atoms_in_ring(graph, backedge[0], backedge[1]);
            cycles_in_molecule.emplace_back(one_cycle);
          }

          std::vector<std::pair<int, int>> edges_to_ruin;

          for (auto cycle = 1; cycle < cycles_in_molecule.size(); ++cycle) {

            auto adj_first_atom_ring =
                adjacents(graph, cycles_in_molecule[cycle][0]);

            for (auto bridge_atom : adj_first_atom_ring) {
              if (bridge_atom != cycles_in_molecule[cycle][1] &&
                  bridge_atom != cycles_in_molecule[cycle][5] &&
                  graph[bridge_atom].element != "H") {
                std::pair<int, int> ruin_edge(bridge_atom,
                                              cycles_in_molecule[cycle][0]);
                edges_to_ruin.emplace_back(ruin_edge);
              } else
                continue;
            }
          }
          for (auto one_edge_to_ruin : edges_to_ruin) {
            boost::remove_edge(one_edge_to_ruin.first, one_edge_to_ruin.second,
                               graph);
          }
          auto all_subunits = find_connected_components(graph);

          for (auto one_edge_to_ruin : edges_to_ruin) {
            boost::add_edge(one_edge_to_ruin.first, one_edge_to_ruin.second,
                            graph);
          }
          for (auto unit : all_subunits) {
            for (auto i = 0; i < back_edges.size(); ++i) {
              for (auto atom_unit : unit) {
                if (atom_unit == back_edges[i][0]) {

                  aromatic_subunits.emplace_back(unit);
                }
              }
            }
          }
          return aromatic_subunits;
        }
      }

      /*gets the atoms of a (six-membered) ring
      first atom in vector is the first atom visited, which is the C binding to
      three Cs connecting bridge and ring. next one is adj, after that adj of
      adj, fourth is the opposit atom binding to three Cs. Last entry is the
      atom binding to the first one.*/
      static std::vector<size_t> find_atoms_in_ring(Graph_type graph,
                                                    int source, int target) {

        int first_atom_ring, second_atom_ring;
        std::vector<size_t> atoms_of_ring;

        auto adj_source = adjacents(graph, source);
        auto adj_target = adjacents(graph, target);
        if (graph[adj_source[0]].element == "C" &&
            graph[adj_source[1]].element == "C" &&
            graph[adj_source[2]].element == "C") {
          // source is atom connecting "bridge" and ring

          first_atom_ring = source;
          second_atom_ring = target;

        } else if (graph[adj_target[0]].element == "C" &&
                   graph[adj_target[1]].element == "C" &&
                   graph[adj_target[2]].element == "C") {
          // target is atom connecting "bridge" and ring

          first_atom_ring = target;
          second_atom_ring = source;

        } else {

          atoms_of_ring = first_cycle(graph, source, target);
          return atoms_of_ring;
        }

        /* while loop could be more efficent, especially for NOT six-membered
         rings. should already work for five-membered ones. aditionally while
         loop searches for atoms between (and including) both Cs connecting ring
         and bridge*/

        double atomC_binds_three_Cs = 0.5;
        int atom_current, atom_before;
        atom_current = second_atom_ring;
        atom_before = first_atom_ring;
        atoms_of_ring.emplace_back(atom_before);
        atoms_of_ring.emplace_back(atom_current);

        while (atom_current != atomC_binds_three_Cs) {

          auto adj_atom_current = adjacents(graph, atom_current);

          if (graph[adj_atom_current[0]].element == "C" &&
              graph[adj_atom_current[1]].element == "C" &&
              graph[adj_atom_current[2]].element == "C") {

            atomC_binds_three_Cs = atom_current;
          } else {

            for (auto atom_next : adj_atom_current) {
              if (atom_next != atom_before && graph[atom_next].element != "H") {

                atoms_of_ring.emplace_back(atom_next);
                atom_before = atom_current;
                atom_current = atom_next;
                break;
              }
            }
          }
        }
        if (first_atom_ring == atomC_binds_three_Cs) {
          atoms_of_ring.pop_back();
          return atoms_of_ring;

        }

        else {

          auto adj_first_atom = adjacents(graph, first_atom_ring);
          auto adj_C_binds_three_Cs = adjacents(graph, atomC_binds_three_Cs);

          for (auto atom_next_to_first : adj_first_atom) {
            for (auto atom_next_to_C_CCC : adj_C_binds_three_Cs) {
              if (edge_exist(graph, atom_next_to_first, atom_next_to_C_CCC)) {

                atoms_of_ring.emplace_back(atom_next_to_C_CCC);
                atoms_of_ring.emplace_back(atom_next_to_first);
                return atoms_of_ring;
              } else if (atom_next_to_first == atom_next_to_C_CCC) {
                atoms_of_ring.emplace_back(atom_next_to_first);
                return atoms_of_ring;
              }
            }
          }
        }
        std::cout << "Ring search failed probably. This ring seems to contain "
                  << atoms_of_ring.size() << " atoms. Is this right?"
                  << std::endl;
        return atoms_of_ring;
      }

      static bool edge_exist(Graph_type graph, int one_atom, int another_atom) {
        auto adj_one_atom = adjacents(graph, one_atom);
        for (auto adj_atom : adj_one_atom) {
          if (adj_atom == another_atom) {
            return true;
          } else
            continue;
        }
        return false;
      }

      // funktion is working not only for six-membered rings, but for every one
      // you can imagine IF no edges are shared with another ring last entry in
      // vector cycle is the atom binding to three Cs connecting ring and bridge
      static std::vector<size_t> first_cycle(Graph_type graph, int atom_a,
                                             int atom_f) {
        std::vector<size_t> cycle;

        double C_binds_three_Cs = 0.5;
        int atom_before = atom_f;
        int atom_current = atom_a;

        cycle.emplace_back(atom_a);
        cycle.emplace_back(atom_f);

        while (atom_current != C_binds_three_Cs) {

          auto adj_atom_current = adjacents(graph, atom_current);

          if (graph[adj_atom_current[0]].element == "C" &&
              graph[adj_atom_current[1]].element == "C" &&
              graph[adj_atom_current[2]].element == "C") {

            C_binds_three_Cs = atom_current;
            cycle.pop_back();
          } else
            for (auto atom_next : adj_atom_current) {

              if (atom_next != atom_before && graph[atom_next].element != "H") {

                cycle.emplace_back(atom_next);
                atom_before = atom_current;
                atom_current = atom_next;
                break;
              }
            }
        }
        atom_before = atom_a;
        atom_current = atom_f;

        while (atom_current != atom_a) {

          auto adj_atom_current = adjacents(graph, atom_current);

          if (graph[adj_atom_current[0]].element == "C" &&
              graph[adj_atom_current[1]].element == "C" &&
              graph[adj_atom_current[2]].element == "C") {

            return cycle;
          }
          for (auto atom_next : adj_atom_current) {

            if (atom_next != atom_before && graph[atom_next].element != "H") {

              cycle.emplace_back(atom_next);
              atom_before = atom_current;
              atom_current = atom_next;
              break;
            }
          }
        }

        return cycle;
      }

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

        return create_resids(atom_vec, res_index_vec);
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

        std::vector<std::size_t> resids_vec;

        std::vector<std::vector<std::size_t>> result, temp_result_aa,
            temp_result_cycle;
        std::vector<int> component(boost::num_vertices(graph));

        int num = connected_components(graph, &component[0]);

        if (num > 1) {
          result = find_connected_components(graph);
        }

        temp_result_aa = find_amino_acids(graph);
        for (auto aa : temp_result_aa) {

          result.emplace_back(aa);
        }

        temp_result_cycle = find_aromatic_subunit(graph);
        for (auto cycle_subunit : temp_result_cycle) {
          result.emplace_back(cycle_subunit);
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

        std::vector<std::vector<size_t>> backbones, amino_acids, subgraphs;
        backbones = find_backbone_in_aminoacid(graph);

        if (backbones.size() == 0)
          return amino_acids;

        else
          for (auto one_backbone = 0; one_backbone < backbones.size() - 1;
               ++one_backbone) {

            auto carboxyC = backbones[one_backbone][1];
            auto n_nextaa = backbones[one_backbone + 1][3];

            boost::remove_edge(carboxyC, n_nextaa, graph);
          }

        subgraphs = find_connected_components(graph);

        // subgraphs contains the amino acids and additionally seperated
        // molecules if excisting else no problem
        // needs to be checked!
        if (subgraphs.size() != backbones.size()) {

          for (auto sub = 0; sub < subgraphs.size(); ++sub) {
            for (auto a = 0; a < subgraphs[sub].size(); ++a) {
              for (auto b = 0; b < backbones.size(); ++b) {
                if (subgraphs[sub][a] == backbones[b][0]) {

                  amino_acids.push_back(subgraphs[sub]);
                }
              }
            }
          }
        }

        for (auto one_backbone = 0; one_backbone < backbones.size() - 1;
             ++one_backbone) {

          auto carboxyC = backbones[one_backbone][1];
          auto n_nextaa = backbones[one_backbone + 1][3];

          boost::add_edge(carboxyC, n_nextaa, graph);
        }

        return amino_acids;
      }

      static std::vector<std::vector<size_t>>
      find_backbone_in_aminoacid(const Graph_type& graph) {

        std::vector<std::vector<size_t>> backbones_indices;
        std::vector<size_t> one_backbone;
        int o, c1, c2, n;

        boost::graph_traits<Graph_type>::edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
          int u = boost::source(*ei, graph);
          int v = boost::target(*ei, graph);

          if (graph[u].element == "O" && graph[v].element == "C" &&
              terminal(graph, u)) {

            o = u;
            c1 = v;
          } else if (graph[v].element == "O" && graph[u].element == "C" &&
                     terminal(graph, v)) {

            o = v;
            c1 = u;
          }

          else
            continue;

          auto carboxy_C_vec = adjacents(graph, c1);
          for (auto c2 : carboxy_C_vec) {

            if (graph[c2].element == "C") {

              auto chiral_C_vec = adjacents(graph, c2);

              for (auto n : chiral_C_vec) {
                if (graph[n].element == "N") {
                  one_backbone.emplace_back(o);
                  one_backbone.emplace_back(c1);
                  one_backbone.emplace_back(c2);
                  one_backbone.emplace_back(n);
                  backbones_indices.emplace_back(one_backbone);
                  one_backbone.clear();
                }
              }
            }
          }
        }
        return backbones_indices;
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

    }; // class Parser

    static inline void make_bonds(Atoms&,
                                  std::vector<std::pair<int, int>> const&);

  }; // struct helper

public:
  inline Coordinates read(std::string const&) override;
  std::shared_ptr<coords::input::formats::xyz::helper::Parser<float_type>>
      parser;
}; // struct xyz
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