#pragma once

#include "coords.h"
#include "coords_rep.h"
#include "graph.h"
#include "InternalCoordinates.h"
#include "ConstrainedInternalCoordinates.h"
#include "InternalCoordinateDecorator.h"
#include "Optimizer.h"

#include <iostream>
#include <iomanip>

class ic_testing
{
public: 
	coords::Coordinates* cPtr;
	coords::Representation_3D inp_struc_cartesian;
	coords::Representation_Internal inp_struc_internal;
	coords::Representation_Main inp_struc_main;

  void ic_execution(coords::Coordinates & coords)
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
  
    // output graphviz file from graph
    graph.visualize_graph("Graphviz");
  
    InternalCoordinates::CartesiansForInternalCoordinates cartesians(cp_vec2_bohr);
    
	auto manager = std::make_shared<internals::ConstraintManager>(Config::get().constrained_internals.constraints);

	manager->constrainAllDistances(Config::get().constrained_internals.constrain_bond_lengths)
		.constrainAllAngles(Config::get().constrained_internals.constrain_bond_angles)
		.constrainAllDihedrals(Config::get().constrained_internals.constrain_dihedrals)
		.constrainAllOOPs(Config::get().constrained_internals.constrain_out_of_plane_bends)
		.constrainAllTranslations(Config::get().constrained_internals.constrain_translations)
		.constrainAllRotations(Config::get().constrained_internals.constrain_rotations);

    // create initial internal coordinates system
    auto icSystem = std::make_shared<internals::ConstrainedInternalCoordinates>();
    std::shared_ptr<internals::InternalCoordinatesBase> decorator = icSystem;
    decorator = std::make_shared<internals::ICRotationDecorator>(decorator);
    decorator = std::make_shared<internals::ICTranslationDecorator>(decorator);
    //decorator = std::make_shared<internals::ICOutOfPlaneDecorator>(decorator);
    decorator = std::make_shared<internals::ICDihedralDecorator>(decorator);
    decorator = std::make_shared<internals::ICAngleDecorator>(decorator);
    decorator = std::make_shared<internals::ICBondDecorator>(decorator);
    decorator->buildCoordinates(cartesians, graph, index_vec3, *manager);
    
    std::cout << "CAST delocalized internals read in the following info:\n";
    for (auto const & pic : icSystem->primitive_internals) 
    {
      std::cout << pic->info(cp_vec2_bohr) << "\n";
    }

    std::cout << "Starting...\n" << std::endl;
	  Optimizer optimizer(*icSystem, cartesians);
    optimizer.optimize(coords);
	    
  }
};
