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

void ic_testing::ic_execution(coords::DL_Coordinates & coords) {

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
  auto bonds = ic_util::bonds(p.create_element_vec(), coords::input::formats::pdb_helper::ang_from_bohr(cp_vec));

  // create graph from bonds vector and atom vector
  ic_util::Graph<decltype(p.atom_vec)::value_type> graph(bonds, p.atom_vec);

  // output graphviz file from graph
  graph.visualize_graph("Graphviz");

  // create initial internal coordinates system
  ic_core::system icSystem(residue_vec, index_vec, cp_vec);

  icSystem.create_ic_system(graph.g);

  auto write_with_zero = [](auto&& ofs, auto&& mat) {
    for (auto r = 0; r < mat.rows(); ++r) {
      for (auto c = 0; c < mat.cols(); ++c) {
        ofs << std::setw(10) << std::setprecision(5) << std::fixed << (std::fabs(mat(r, c)) > 1.e-6 ? mat(r, c) : 0.0);
      }
      ofs << "\n";
    }
  };

  std::ofstream offpb("prim_Bmat.dat");
  write_with_zero(offpb, icSystem.Bmat());
  std::ofstream offpg("prim_Gmat.dat");
  write_with_zero(offpg, icSystem.Gmat());

  auto G_matrix = icSystem.delocalize_ic_system();

  std::ofstream offib("ic_Bmat.dat");
  write_with_zero(offib,icSystem.ic_Bmat());
  std::ofstream offig("ic_Gmat.dat");
  write_with_zero(offig,icSystem.ic_Gmat());

  std::cout << "DLC matrix: \n";
  write_with_zero(std::cout, G_matrix);
  std::cout << "\n\n";
  auto&& del_hessian = icSystem.initial_delocalized_hessian();

  std::cout << "DelHessian:\n";
  write_with_zero(std::cout, del_hessian);
  std::cout << "\n\n";

  /*std::cout << "Ginversed:\n" << G_matrix_inv << "\n\n";
  std::cout << "Gmatrix:\n" << G_matrix << "\n\n";*/

  //std::cout << icSystem.angle_vec_.size() << "||" << icSystem.distance_vec_.size() << "||" << icSystem.dihed_vec_.size() << "||" << icSystem.oop_vec_.size() << "||" << icSystem.rotation_vec_.size() << "||" << icSystem.trans_x_vec_.size() << std::endl;
  /*for (auto& i : icSystem.distance_vec_)
  {
	  std::cout << i.dist() << std::endl;
  }*/
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












  // test matrix stuff
  /*auto matrix_trial = ic_util::Rep3D_to_arma<coords::float_type>(cp_vec);
  matrix_trial.print();
  std::cout << "\n";*/
  // create trial and target structure
  // coords::Representation_3D trial = cp_vec;
  // coords::Representation_3D target = ca_vec;
  /*auto matrix_target = ic_util::Rep3D_to_arma<coords::float_type>(target);
  matrix_target.print();
  // test correlation matrix
  auto correl_mat =
      ic_rotation::correlation_matrix<coords::float_type>(trial, target);
  std::cout << "\n";
  std::cout << "Correlation matrix: \n";
  correl_mat.print();
  // test F matrix
  std::cout << "F matrix: \n";
  auto F = ic_rotation::F_matrix<coords::float_type>(trial, target);
  F.print();
  // test quaternion creation
  auto quater = ic_rotation::quaternion<coords::float_type>(trial, target);
  std::cout << "Quaternion: \n";
  std::cout << quater.first << ", " << quater.second << std::endl;
  // test expmap creation
  auto expmap = ic_rotation::exponential_map<coords::float_type>(trial, target);
  std::cout << "Expmap: \n";
  std::cout << expmap.at(0) << ", " << expmap.at(1) << ", " << expmap.at(2)
            << "\n";
  // test of correlation matrix derivatives
  auto correl_mat_ders =
      ic_rotation::correlation_matrix_derivs<coords::float_type>(correl_mat,
                                                                 target);
  std::cout << correl_mat_ders.size() << std::endl;
  // test of F matrix derivatives
  auto f_mat_ders =
      ic_rotation::F_matrix_derivs<coords::float_type>(trial, target);
  // test of quaternion derivatives
  auto q_mat_der =
      ic_rotation::quaternion_derivs<coords::float_type>(trial, target);
  std::cout << "Quaternion derivatives: \n";
  for (auto& j : q_mat_der)
    std::cout << j << std::endl;
  // test of expmap derivatives
  std::cout << "Expmap derivs: " << std::endl;
  auto expmap_derivs =
      ic_rotation::exponential_derivs<coords::float_type>(trial, target);
  for (auto& k : expmap_derivs)
    std::cout << k << std::endl;*/
  // test rotation class
  /*ic_core::rotation rot_vec(target);
  auto expmap_from_rotvec = rot_vec.rot_val(trial);
  std::cout << expmap_from_rotvec.at(0) << " || " << expmap_from_rotvec.at(1)
            << " || " << expmap_from_rotvec.at(2) << std::endl;
  auto expmapdd =
      ic_rotation::exponential_map<coords::float_type>(trial, target);
  std::cout << expmapdd.at(0) << " || " << expmapdd.at(1) << " || "
            << expmapdd.at(2) << std::endl;
  auto expder = rot_vec.rot_der(trial);
  for (auto& q : expder)
    std::cout << q << std::endl;
  auto expder1 =
      ic_rotation::exponential_derivs<coords::float_type>(trial, target);
  for (auto& w : expder1)
    std::cout << w << std::endl;*/
}
