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

#include <boost/graph/breadth_first_search.hpp>

#include <iostream>
#include <iomanip>

template <typename Inserter>
class mst_visitor : public boost::default_bfs_visitor{
public:
  mst_visitor(Inserter inserter):
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

    // Build a minimum spanning tree from breadth-first search
    std::vector<std::pair<std::size_t, std::size_t>> tree_edges;
    auto inserter = std::back_inserter(tree_edges);
    mst_visitor<decltype(inserter)> v(inserter);
    boost::breadth_first_search(graph, 0, boost::visitor(v));

    boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, ic_util::Node> spanning_tree{
        tree_edges.begin(),
        tree_edges.end(),
        curGraphinfo.size()
    };

    for (auto i = 0u; i < curGraphinfo.size(); ++i) {
      spanning_tree[i].atom_serial = curGraphinfo.at(i).atom_serial;
      spanning_tree[i].atom_name = curGraphinfo.at(i).atom_name;
      spanning_tree[i].element = curGraphinfo.at(i).element;
      spanning_tree[i].cp = curGraphinfo.at(i).cp;
    }

    std::ofstream s("spanning-tree.txt");
    boost::write_graphviz(s, spanning_tree, boost::make_label_writer(boost::get(&ic_util::Node::atom_name, spanning_tree)));

    // output graphviz file from graph
    /*graph.visualize_graph("Graphviz");

    InternalCoordinates::CartesiansForInternalCoordinates cartesians(cp_vec2_bohr);

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
