#pragma once

#include "ic_exec.h"
#include "coords.h"
#include "coords_rep.h"
#include "graph.h"
#include "ic_core.h"
#include "ic_rotation.h"
#include "ic_util.h"
#include "quaternion.h"
#include "scon_angle.h"
#include "scon_spherical.h"
#include "scon_vect.h"

#include "scon_mathmatrix.h"
#include "coords_io_pdb.h"
#include <array>
#include <iostream>
#include <iomanip>
#include <vector>

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
            auto residue_vec = p.create_resids_rep_3D_bohr();

            // create residue index vector from Parser atom vector
            auto index_vec = p.create_resids_indices();

            auto el_vec = p.create_element_vec();

            // create vector of bonds
            auto bonds = ic_util::bonds(p.create_element_vec(), coords::input::formats::pdb::helper::ang_from_bohr(cp_vec));

            // create graph from bonds vector and atom vector
            auto  graph = ic_util::make_graph(bonds, p.atom_vec);

            // output graphviz file from graph
            graph.visualize_graph("Graphviz");

            // create initial internal coordinates system
            ic_core::system icSystem(residue_vec, index_vec, cp_vec);

            icSystem.create_ic_system(graph.g);

            //auto write_with_zero = [](auto&& ofs, auto&& mat) {
            //  for (auto r = 0; r < mat.rows(); ++r) {
            //    for (auto c = 0; c < mat.cols(); ++c) {
            //      ofs << std::setw(10) << std::setprecision(5) << std::fixed << (std::fabs(mat(r, c)) > 1.e-6 ? mat(r, c) : 0.0);
            //    }
            //    ofs << "\n";
            //  }
            //};

            icSystem.delocalize_ic_system();

            std::cout << "Rotation:\n";
            for (auto& i : icSystem.rotation_vec_)
            {
                    auto j = i.rot_val(cp_vec);
                    std::cout << j.at(0) << "||" << j.at(1) << "||" << j.at(2) << "||" << std::endl;
            }

            for (auto const & pic : icSystem.primitive_internals) {
              std::cout << pic->info(cp_vec);
            }
            icSystem.optimize(coords);
        }
};