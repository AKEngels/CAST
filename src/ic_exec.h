#pragma once

#include "coords.h"
#include "coords_rep.h"
#include "graph.h"
#include "InternalCoordinates.h"
#include "TranslationRotationInternalCoordinates.h"
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
        template<typename T>
	void ic_execution(coords::DL_Coordinates<T> & coords){

            // test of my code
            // create a Parser object from the input file
            /*Pdb::Parser<coords::float_type> pAlaGly("ethanol.pdb");
            Pdb::Parser<coords::float_type> p0("test.pdb");*/
            auto const& p = *coords.parser.get();

            auto cp_vec = p.create_rep_3D_bohr();

            // create residue vector from Parser atom vector
            //auto residue_vec = p.create_resids_rep_3D_bohr();

            // create residue index vector from Parser atom vector
            auto index_vec = p.create_resids_indices();

            auto el_vec = p.create_element_vec();

            // create vector of bonds
            auto bonds = ic_util::bonds(p.create_element_vec(), coords::input::formats::pdb::helper::ang_from_bohr(cp_vec));

            // create graph from bonds vector and atom vector
            ic_util::Graph<ic_util::Node> graph = ic_util::make_graph(bonds, p.atom_vec);

            // output graphviz file from graph
            graph.visualize_graph("Graphviz");

            InternalCoordinates::CartesiansForInternalCoordinates cartesians(cp_vec);
            
            // create initial internal coordinates system
            //internals::TRIC icSystem(residue_vec, index_vec, cartesians, graph);
            //internals::PrimitiveInternalsTransRot icSystem(residue_vec, index_vec, cartesians, graph);
            
            auto icSystem = std::make_shared<internals::TRIC>();
            std::shared_ptr<internals::InternalCoordinatesBase> decorator = icSystem;
            decorator = std::make_shared<internals::ICRotationDecorator>(decorator);
            decorator = std::make_shared<internals::ICTranslationDecorator>(decorator);
            decorator = std::make_shared<internals::ICOutOfPlaneDecorator>(decorator);
            decorator = std::make_shared<internals::ICDihedralDecorator>(decorator);
            decorator = std::make_shared<internals::ICAngleDecorator>(decorator);
            decorator = std::make_shared<internals::ICBondDecorator>(decorator);
            decorator->buildCoordinates(cartesians, graph, index_vec);
            icSystem->delocalize_ic_system(cartesians);
        
            
            //auto write_with_zero = [](auto&& ofs, auto&& mat) {
            //  for (auto r = 0; r < mat.rows(); ++r) {
            //    for (auto c = 0; c < mat.cols(); ++c) {
            //      ofs << std::setw(10) << std::setprecision(5) << std::fixed << (std::fabs(mat(r, c)) > 1.e-6 ? mat(r, c) : 0.0);
            //    }
            //    ofs << "\n";
            //  }
            //};

            for (auto const & pic : icSystem->primitive_internals) {
              std::cout << pic->info(cp_vec) << "\n";
            }

	    Optimizer optimizer(*icSystem, cartesians);
            optimizer.optimize(coords);
	    
        }
};