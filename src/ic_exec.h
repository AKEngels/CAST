#pragma once

#include "coords.h"
#include "coords_rep.h"
#include "InternalCoordinates/BondGraph.h"
#include "InternalCoordinates/InternalCoordinates.h"
#include "InternalCoordinates/ConstrainedInternalCoordinates.h"
#include "InternalCoordinates/InternalCoordinateDecorator.h"
#include "InternalCoordinates/TranslationRotationInternalCoordinates.h"
#include "InternalCoordinates/Optimizer.h"
#include "InternalCoordinates/InternalCoordinateUtilities.h"

#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <type_traits>
#include "Scon/scon_angle.h"

#include <iostream>
#include <iomanip>

template <typename Inserter>
class mst_visitor : public boost::default_bfs_visitor{
public:
  explicit mst_visitor(Inserter inserter):
    m_inserter{inserter}
  {}

  template <typename EdgeDescriptor, typename Graph>
  void tree_edge(EdgeDescriptor edge, Graph /*g*/){
    //std::cout << edge << '\n';
    m_inserter = std::make_pair(edge.m_source, edge.m_target);
  }

private:
  Inserter m_inserter;
};

struct z_matrix_node : public ic_util::Node{
  template <typename Graph>
  void print_z_matrix_entry(std::ostream & out, Graph const& g, coords::Representation_3D const& c) const{
    auto print_coord = [&out, &g, &c](auto coord_pair){
      auto val = coord_pair->first.val(c);
      auto is_distance = std::is_same<decltype(coord_pair->first), InternalCoordinates::BondDistance>::value;

      out << ' ' << g[coord_pair->second].m_position << ' '
          << (is_distance ? coord_pair->first.val(c) : scon::ang<decltype(val)>::from_rad(val).degrees());
    };

    out << element;
    if (m_distance) {
      print_coord(m_distance);
      if (m_angle){
        print_coord(m_angle);
        if (m_dihedral){
          print_coord(m_dihedral);
        }
      }
    }
    out << '\n';
  }

  std::optional<std::pair<InternalCoordinates::BondDistance, std::size_t>> m_distance;
  std::optional<std::pair<InternalCoordinates::BondAngle, std::size_t>> m_angle;
  std::optional<std::pair<InternalCoordinates::DihedralAngle, std::size_t>> m_dihedral;

  // Position in the Z matrix
  std::size_t m_position;
};

template <typename Vertex, typename Graph>
std::optional<Vertex> get_parent_vertex(Vertex u, Graph g){
  auto pair = boost::inv_adjacent_vertices(u, g);

  // Look whether we have a parent vertex
  if (pair.first == pair.second)
    return std::nullopt;
  else
    return *pair.first;
}

template <typename Graph>
class z_matrix_visitor : public boost::default_dfs_visitor {
  using VertexDescriptor = typename boost::graph_traits<Graph>::vertex_descriptor;

public:
  explicit z_matrix_visitor(Graph &g, std::vector<VertexDescriptor> &z_matrix_order) :
          m_graph{g},
          m_z_matrix_order{z_matrix_order}
          {}

  void discover_vertex(VertexDescriptor u, Graph /*g*/) {
    m_graph[u].m_position = m_curr_vertex;
    m_z_matrix_order.emplace_back(u);

    auto anchor = u;
    std::vector<VertexDescriptor> curr_anchors;
    for (std::size_t i = 0; i < 3; ++i) {
      auto parent = get_parent_vertex(anchor, m_graph);
      if (parent) {
        curr_anchors.emplace_back(*parent);
        // Construct the appropriate internal coordinates
        switch (i) {
          case 0:
            m_graph[u].m_distance = std::make_pair(
                    InternalCoordinates::BondDistance(m_graph[curr_anchors[0]], m_graph[u]), curr_anchors[0]);
            break;
          case 1:
            m_graph[u].m_angle = std::make_pair(
                    InternalCoordinates::BondAngle(m_graph[curr_anchors[1]], m_graph[curr_anchors[0]], m_graph[u]),
                    curr_anchors[1]);
            break;
          case 2:
            m_graph[u].m_dihedral = std::make_pair(
                    InternalCoordinates::DihedralAngle(m_graph[curr_anchors[2]], m_graph[curr_anchors[1]],
                                                       m_graph[curr_anchors[0]], m_graph[u]), curr_anchors[2]);
            break;
          default:
            break;
        }
        anchor = *parent;
      } else break;
    }
    ++m_curr_vertex;
  }

private:
  std::size_t m_curr_vertex{1};
  Graph &m_graph;
  std::vector<VertexDescriptor> &m_z_matrix_order;
};

template <typename Graph>
boost::optional<typename boost::graph_traits<Graph>::vertex_descriptor>
find_root_vertex(Graph const& g){
  using VertexDescriptor = typename boost::graph_traits<Graph>::vertex_descriptor;

  auto get_vertex_count = [](auto vertex_range){
    return vertex_range.second - vertex_range.first;
  };

  auto vertices = boost::vertices(g);
  for (auto it = vertices.first; it != vertices.second; ++it){
    auto adj_vertices = boost::adjacent_vertices(*it, g);

    if (get_vertex_count(adj_vertices) == 2){
      for(auto adj_vertex_it = adj_vertices.first; adj_vertex_it != adj_vertices.second; ++adj_vertex_it){
        if (get_vertex_count(boost::adjacent_vertices(*adj_vertex_it, g)) == 1)
          return static_cast<VertexDescriptor>(*adj_vertex_it);
      }
    }
  }
  return boost::none;
}

class ic_testing
{
public:
  coords::Coordinates* cPtr;
  coords::Representation_3D inp_struc_cartesian;
  coords::Representation_Internal inp_struc_internal;
  coords::Representation_Main inp_struc_main;

  void ic_execution(coords::Coordinates& coords)
  {
    //auto const& p = *coords.parser.get();


    //OLD
    //auto cp_vec = p.create_rep_3D_bohr();
    //NEW
    auto cp_vec2 = coords::Representation_3D(coords.xyz());
    auto cp_vec2_bohr = cp_vec2;
    for (auto&& coord_ : cp_vec2_bohr)
    {
      coord_ /= energy::bohr2ang;
    }

    // create residue index vector from Parser atom vector
    //OLD
    //auto index_vec = p.create_resids_indices();
    //NEW
    auto index_vec2 = coords.molecules(); //--> Probably +1 has to be added, indexation wrong
    std::vector<std::vector<std::size_t>> index_vec3;
    for (auto&& i : index_vec2)
    {
      for (auto&& j : i)
      {
        j += 1;
      }
      index_vec3.push_back(i);
    }


    //OLD
    //auto el_vec = p.create_element_vec();
    //NEW
    std::vector<std::string> el_vec2;
    for (auto&& i : coords.atoms())
    {
      el_vec2.emplace_back(i.symbol());
    }

    struct graphinfo
    {
      std::size_t atom_serial;
      std::string atom_name, element;
      coords::Cartesian_Point cp;
    };
    std::vector<graphinfo> curGraphinfo;
    for (std::size_t i = 0u; i < coords.xyz().size(); i++)
    {
      graphinfo tempinfo;
      tempinfo.cp = coords.xyz(i);
      tempinfo.atom_serial = i + 1u;
      tempinfo.atom_name = coords.atoms(i).symbol();
      tempinfo.element = tempinfo.atom_name;
      curGraphinfo.emplace_back(tempinfo);
    }


    // create vector of bonds
    auto bonds = ic_util::bonds(el_vec2, cp_vec2);

    // create graph from bonds vector and atom vector
    ic_util::Graph<ic_util::Node> graph = ic_util::make_graph(bonds, curGraphinfo);

    auto const root_index = find_root_vertex(graph);
    if (root_index) {
      std::cout << "Index of root vertex: " << *root_index << '\n';

      // Build a minimum spanning tree from breadth-first search
      std::vector<std::pair<std::size_t, std::size_t>> tree_edges;
      auto inserter = std::back_inserter(tree_edges);
      mst_visitor v(inserter);
      boost::breadth_first_search(graph, *root_index, boost::visitor(v));

      boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, z_matrix_node> spanning_tree{
          tree_edges.begin(),
          tree_edges.end(),
          curGraphinfo.size()
      };

      for (std::size_t i = 0; i < curGraphinfo.size(); ++i) {
        spanning_tree[i].atom_serial = curGraphinfo.at(i).atom_serial;
        spanning_tree[i].atom_name = curGraphinfo.at(i).atom_name;
        spanning_tree[i].element = curGraphinfo.at(i).element;
        spanning_tree[i].cp = curGraphinfo.at(i).cp;
      }

      std::vector<boost::graph_traits<decltype(spanning_tree)>::vertex_descriptor> z_matrix_order;
      z_matrix_visitor vis{spanning_tree, z_matrix_order};
      boost::depth_first_search(spanning_tree, boost::visitor(vis).root_vertex(*root_index));

      std::ofstream s("spanning-tree.txt");
      boost::write_graphviz(s, spanning_tree,
                            boost::make_label_writer(boost::get(&ic_util::Node::atom_name, spanning_tree)));

      std::ofstream o("z-matrix.txt");
      for (auto curr_vertex : z_matrix_order)
        spanning_tree[curr_vertex].print_z_matrix_entry(o, spanning_tree, cp_vec2);
      o.close();
    }
    else
      std::cout << "No suitable root vertex found\n";

    // output graphviz file from graph
    graph.visualize_graph("Graphviz");

    /*InternalCoordinates::CartesiansForInternalCoordinates cartesians(cp_vec2_bohr);

    auto manager = std::make_shared<internals::ConstraintManager>(Config::get().constrained_internals.constraints);

    manager->constrainAllDistances(Config::get().constrained_internals.constrain_bond_lengths)
      .constrainAllAngles(Config::get().constrained_internals.constrain_bond_angles)
      .constrainAllDihedrals(Config::get().constrained_internals.constrain_dihedrals)
      .constrainAllOOPs(Config::get().constrained_internals.constrain_out_of_plane_bends)
      .constrainAllTranslations(Config::get().constrained_internals.constrain_translations)
      .constrainAllRotations(Config::get().constrained_internals.constrain_rotations);

    // create initial internal coordinates system
    std::unique_ptr<internals::ICDecoratorBase> decorator{ nullptr };
    decorator = std::make_unique<internals::ICRotationDecorator>(std::move(decorator));
    decorator = std::make_unique<internals::ICTranslationDecorator>(std::move(decorator));
    decorator = std::make_unique<internals::ICOutOfPlaneDecorator>(std::move(decorator));
    decorator = std::make_unique<internals::ICDihedralDecorator>(std::move(decorator));
    decorator = std::make_unique<internals::ICAngleDecorator>(std::move(decorator));
    decorator = std::make_unique<internals::ICBondDecorator>(std::move(decorator));

    decorator->buildCoordinates(cartesians, graph, index_vec3, *manager);
    internals::ConstrainedInternalCoordinates icSystem{ *decorator };

    if (Config::get().general.verbosity > 0)
    {
      std::cout << "CAST delocalized internals read in the following info:\n";
      for (auto const& pic : icSystem.primitive_internals)
      {
        std::cout << pic->info(cp_vec2_bohr) << "\n";
      }
      std::cout << "Starting...\n" << std::endl;
    }
    Optimizer optimizer(icSystem, cartesians);
    optimizer.optimize(coords);*/

  }
};
