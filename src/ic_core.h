////////////////////////////////
// new functions to be tested //
////////////////////////////////
#ifndef cast_ic_core_h_guard
#define cast_ic_core_h_guard

#include "coords.h"
#include "coords_rep.h"
#include "ic_atom.h"
#include "ic_rotation.h"
#include "pdb.h"
#include "scon_angle.h"
#include "scon_spherical.h"
#include "scon_vect.h"

#include <algorithm>
#include "scon_mathmatrix.h"
#include <array>
#include <boost/graph/adjacency_list.hpp>
#include <cmath>
#include <iterator>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include<memory>

#include "coords_io_pdb.h"

namespace ic_core {

using coords::float_type;

struct internal_coord {
  virtual float_type val(coords::Representation_3D const& xyz) const = 0;
  virtual std::vector<float_type> der_vec(coords::Representation_3D const& xyz) const = 0;
  virtual float_type hessian_guess(coords::Representation_3D const& xyz) const = 0;
  virtual std::string info(coords::Representation_3D const & xyz) const = 0;
};

struct distance : public internal_coord {
  distance(const unsigned int& index_a, const unsigned int& index_b,
           const std::string& elem_a, const std::string& elem_b)
      : index_a_{ index_a }, index_b_{ index_b },
        elem_a_{ elem_a }, elem_b_{ elem_b } {}

  std::size_t index_a_;
  std::size_t index_b_;
  std::string elem_a_;
  std::string elem_b_;

  float_type val(coords::Representation_3D const& xyz) const override;
  std::pair<scon::c3<float_type>, scon::c3<float_type>> der(coords::Representation_3D const& xyz) const;
  std::vector<float_type> der_vec(coords::Representation_3D const& xyz) const override;
  float_type hessian_guess(coords::Representation_3D const& xyz) const override;
  std::string info(coords::Representation_3D const& xyz) const override;
};

struct angle : internal_coord {
  angle(const unsigned int& index_a,
        const unsigned int& index_b, const unsigned int& index_c,
        const std::string& elem_a, const std::string& elem_b,
        const std::string& elem_c)
      : index_a_{ index_a }, index_b_{ index_b },
        index_c_{ index_c }, elem_a_{ elem_a }, elem_b_{ elem_b }, elem_c_{
          elem_c
        } {}

  std::size_t index_a_;
  std::size_t index_b_;
  std::size_t index_c_;
  std::string elem_a_;
  std::string elem_b_;
  std::string elem_c_;

  coords::float_type val(coords::Representation_3D const& xyz) const override;
  std::tuple<scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>> der(coords::Representation_3D const& xyz) const;
  std::vector<float_type> der_vec(coords::Representation_3D const& xyz) const override;
  float_type hessian_guess(coords::Representation_3D const& xyz) const override;
  std::string info(coords::Representation_3D const& xyz) const override;
};
//
//template<bool isLValue>
//class dihedral_points;
//
//class dihedral_points<true> {
//  using rwcp = std::reference_wrapper<coords::Cartesian_Point const>;
//  rwcp a_;
//  rwcp b_;
//  rwcp c_;
//  rwcp d_;
//};
//
//class dihedral_points<false> {
//  coords::Cartesian_Point a_;
//  coords::Cartesian_Point b_;
//  coords::Cartesian_Point c_;
//  coords::Cartesian_Point d_;
//};
//
//template<bool isLValue>
//class dihedral_indices;
//
//class dihedral_indices<true> {
//  using rwst = std::reference_wrapper<std::size_t const>;
//  rwst index_a_;
//  rwst index_b_;
//  rwst index_c_;
//  rwst index_d_;
//};
//
//class dihedral_indices<false> {
//  std::size_t index_a_;
//  std::size_t index_b_;
//  std::size_t index_c_;
//  std::size_t index_d_;
//};

//template<bool CP_isLVal, bool Ind_isLVal>
struct dihedral : public internal_coord /*: public dihedral_points<CP_isLVal>, public dihedral_indices<Ind_isLVal> */{
  dihedral(const unsigned int& index_a, const unsigned int& index_b,
           const unsigned int& index_c, const unsigned int& index_d)
      : index_a_{ index_a },
        index_b_{ index_b }, index_c_{ index_c }, index_d_{ index_d } {}

  coords::float_type val(coords::Representation_3D const& xyz) const override;
  std::tuple<scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>, scon::c3<float_type>> 
    der(coords::Representation_3D const& xyz) const;
  std::vector<float_type> der_vec(coords::Representation_3D const& xyz) const override;
  float_type hessian_guess(coords::Representation_3D const& xyz) const override;
  std::string info(coords::Representation_3D const& xyz) const override;

  std::size_t get_ind_a() { return index_a_; }
  std::size_t get_ind_b() { return index_b_; }
  std::size_t get_ind_c() { return index_c_; }
  std::size_t get_ind_d() { return index_d_; }

  std::size_t index_a_;
  std::size_t index_b_;
  std::size_t index_c_;
  std::size_t index_d_;
};

struct out_of_plane : public internal_coord {
  out_of_plane(const unsigned int& index_a,
               const unsigned int& index_b, unsigned int& index_c,
               const unsigned int& index_d)
      : index_a_{ index_a },
        index_b_{ index_b }, index_c_{ index_c }, index_d_{ index_d } {}

  std::size_t index_a_;
  std::size_t index_b_;
  std::size_t index_c_;
  std::size_t index_d_;

  coords::float_type val(coords::Representation_3D const& xyz) const override;
  std::vector<scon::c3<float_type>> der(coords::Representation_3D const& xyz) const;
  std::vector<float_type> der_vec(coords::Representation_3D const& xyz) const override;
  float_type hessian_guess(coords::Representation_3D const& xyz) const override;
  std::string info(coords::Representation_3D const& xyz) const override;
};

struct trans : public internal_coord {
  trans(std::vector<std::size_t> const& index_vec)
    : indices_(index_vec) {}

  virtual float_type val(coords::Representation_3D const&) const = 0;
  std::vector<std::size_t> indices_;

  template<typename Func>
  coords::Representation_3D der(std::size_t const & sys_size, Func const & coord) const {
    using rep3D = coords::Representation_3D;
    using cp = coords::Cartesian_Point;

    rep3D result(sys_size, cp(0.,0.,0.));

    auto const & s{ indices_.size() };
    for (auto const & i : indices_) {
      result.at(i-1) = coord(s);
    }
    return result;
  }
  float_type hessian_guess(coords::Representation_3D const& xyz) const override {
    return 0.05;
  }
};

struct trans_x : trans{
  trans_x(const std::vector<std::size_t>& index_vec)
      : trans(index_vec) {}

  float_type val(const coords::Representation_3D& xyz) const override {
    auto coord_sum{ 0.0 };
    for (auto& i : indices_) {
      coord_sum += xyz.at(i-1u).x();
    }
    return coord_sum / indices_.size();
  }

  std::vector<float_type> der_vec(coords::Representation_3D const& rep)const override;
  std::string info(coords::Representation_3D const& xyz) const override;
};

struct trans_y : trans {
  trans_y(const std::vector<std::size_t>& index_vec)
      : trans(index_vec) {}

  float_type val(const coords::Representation_3D& xyz) const override {
    auto coord_sum{ 0.0 };
    for (auto& i : indices_) {
      coord_sum += xyz.at(i - 1u).y();
    }
    return coord_sum / indices_.size();
  }

  std::vector<float_type> der_vec(coords::Representation_3D const& rep)const override;
  std::string info(coords::Representation_3D const& xyz) const override;
};

struct trans_z : trans {
  trans_z(const std::vector<std::size_t>& index_vec)
      : trans(index_vec) {}

  float_type val(const coords::Representation_3D& xyz) const override {
    auto coord_sum{ 0.0 };
    for (auto& i : indices_) {
      coord_sum += xyz.at(i - 1u).z();
    }
    return coord_sum / indices_.size();
  }

  std::vector<float_type> der_vec(coords::Representation_3D const& rep)const override;
  std::string info(coords::Representation_3D const& xyz) const override;
};

struct rotation {
  rotation(const coords::Representation_3D& target,
    const std::vector<std::size_t>& index_vec);

  std::vector<std::size_t> indices_;

  coords::Representation_3D reference_;
  float_type rad_gyr_;

  std::array<float_type, 3u> rot_val(const coords::Representation_3D&) const;
  std::vector<scon::mathmatrix<float_type>> rot_der(const coords::Representation_3D&) const;
  scon::mathmatrix<float_type> rot_der_mat(std::size_t const&, const coords::Representation_3D&) const;
  float_type radius_gyration(const coords::Representation_3D&);

  static coords::Representation_3D xyz0;
};

class system {
public:
  system(const std::vector<coords::Representation_3D>& res_init,
         const coords::Representation_3D& rep_init,
         const std::vector<std::vector<std::size_t>>& res_index)
      : res_vec_{ res_init }, rep_{ rep_init }, res_index_vec_{ res_index } {}

  std::vector<std::unique_ptr<internal_coord>> primitive_internals;
  std::vector<rotation> rotation_vec_;

private:
  const std::vector<coords::Representation_3D> res_vec_;
  const std::vector<std::vector<std::size_t>> res_index_vec_;
  const coords::Representation_3D rep_;

public:
  void append_primitives(std::vector<std::unique_ptr<internal_coord>> && pic) {
    primitive_internals.insert(primitive_internals.end(), 
      std::make_move_iterator(pic.begin()), 
      std::make_move_iterator(pic.end()));
  }

  std::vector<std::unique_ptr<internal_coord>>
  create_trans_x(const std::vector<std::vector<std::size_t>>&) const;

  std::vector<std::unique_ptr<internal_coord>>
  create_trans_y(const std::vector<std::vector<std::size_t>>&) const;

  std::vector<std::unique_ptr<internal_coord>>
  create_trans_z(const std::vector<std::vector<std::size_t>>&) const;

  std::tuple<std::vector<std::unique_ptr<internal_coord>>, std::vector<std::unique_ptr<internal_coord>>, std::vector<std::unique_ptr<internal_coord>>>
    create_translations(const std::vector<std::vector<std::size_t>>&) const;

  std::vector<rotation>
  create_rotations(const coords::Representation_3D&,
                   const std::vector<std::vector<std::size_t>>&);

  /*template<typename Graph>
  IC_System create_system(Graph const &);*/

  template <typename Graph>
  std::vector<std::unique_ptr<internal_coord>> create_distances(const Graph&) const;

  template <typename Graph>
  std::vector<std::unique_ptr<internal_coord>> create_angles(const Graph&) const;

  template <typename Graph>
  std::vector<std::unique_ptr<internal_coord>> create_oops(const coords::Representation_3D&,
                                        const Graph&) const;

  template <typename Graph>
  std::vector<std::unique_ptr<internal_coord>> create_dihedrals(const Graph&) const;

  /*template <typename Graph>
  void create_ic_system(const std::vector<coords::Representation_3D>&,
                        const coords::Representation_3D&, const Graph&);*/

  template <typename Graph>
  void create_ic_system(const Graph&);

  scon::mathmatrix<float_type> del_mat;

  scon::mathmatrix<float_type> delocalize_ic_system(coords::Representation_3D const &);
  scon::mathmatrix<float_type> initial_hessian(coords::Representation_3D const &);
  scon::mathmatrix<float_type> delocalize_hessian(scon::mathmatrix<float_type> const &);
  scon::mathmatrix<float_type> G_mat_inversion(scon::mathmatrix<float_type> const &);
  scon::mathmatrix<float_type> Bmat(const coords::Representation_3D& trial) const;

  scon::mathmatrix<float_type> calc(coords::Representation_3D const&) const;
  scon::mathmatrix<float_type> calcGrad(coords::Representation_3D const&, coords::Representation_3D const&) const;

};

template <typename Graph>
inline std::vector<std::unique_ptr<internal_coord>>
system::create_distances(const Graph& g) const {
  using boost::edges;
  using boost::source;
  using boost::target;

  std::vector<std::unique_ptr<internal_coord>> result;
  auto ed = edges(g);
  for (auto it = ed.first; it != ed.second; ++it) {
    auto u = source(*it, g);
    auto v = target(*it, g);
    auto u_index = g[u].atom_serial;
    auto v_index = g[v].atom_serial;
    auto u_elem = g[u].element;
    auto v_elem = g[v].element;
    result.emplace_back(std::make_unique<distance>(u_index, v_index, u_elem, v_elem));
  }
  return result;
}

template <typename Graph>
inline std::vector<std::unique_ptr<internal_coord>>
system::create_angles(const Graph& g) const {
  using boost::adjacent_vertices;
  using boost::vertices;

  std::vector<std::unique_ptr<internal_coord>> result;
  auto vert = vertices(g);
  for (auto it = vert.first; it != vert.second; ++it) {
    auto a_vert = adjacent_vertices(*it, g);
    for (auto it2 = a_vert.first; it2 != a_vert.second; ++it2) {
      for (auto it3 = a_vert.first; it3 != a_vert.second; ++it3) {
        if (g[*it2].atom_serial < g[*it3].atom_serial) {
          auto a_index = g[*it2].atom_serial;
          auto b_index = g[*it].atom_serial;
          auto c_index = g[*it3].atom_serial;
          auto a_elem = g[*it2].element;
          auto b_elem = g[*it].element;
          auto c_elem = g[*it3].element;
          result.emplace_back(std::make_unique<angle>(a_index, b_index, c_index, a_elem, b_elem, c_elem));
        }
      }
    }
  }
  return result;
}

template <typename Graph>
inline std::vector<std::unique_ptr<internal_coord>>
system::create_oops(const coords::Representation_3D& coords, const Graph& g) const {
  using boost::adjacent_vertices;
  using boost::vertices;
  using scon::dot;

  std::vector<std::unique_ptr<internal_coord>> result;
  auto vert = vertices(g);
  for (auto it = vert.first; it != vert.second; ++it) {
    auto a_vert = adjacent_vertices(*it, g);
    for (auto it2 = a_vert.first; it2 != a_vert.second; ++it2) {
      for (auto it3 = a_vert.first; it3 != a_vert.second; ++it3) {
        for (auto it4 = a_vert.first; it4 != a_vert.second; ++it4) {
          if (g[*it2].atom_serial < g[*it3].atom_serial &&
              g[*it3].atom_serial < g[*it4].atom_serial) {
            auto core = g[*it].atom_serial;
            auto core_cp = coords.at(core - 1);
            auto a = g[*it2].atom_serial;
            auto b = g[*it3].atom_serial;
            auto c = g[*it4].atom_serial;
            std::vector<unsigned int> vec = { a, b, c };
            auto permutation_vec = ic_util::permutation_from_vec(vec);
            for (auto& i : permutation_vec) {
              auto u_cp = coords.at(i.at(0) - 1);
              auto v_cp = coords.at(i.at(1) - 1);
              auto w_cp = coords.at(i.at(2) - 1);
              auto n_vec1 = ic_util::normal_unit_vector(u_cp, v_cp, w_cp);
              auto n_vec2 = ic_util::normal_unit_vector(core_cp, u_cp, v_cp);
              auto dot_n_vecs = dot(n_vec1, n_vec2);
              if (0.95 < dot_n_vecs && dot_n_vecs < 1.05) {
                result.emplace_back(std::make_unique<out_of_plane>(core, a, b, c));
              }
            }
          }
        }
      }
    }
  }
  return result;
}

template <typename Graph>
inline std::vector<std::unique_ptr<internal_coord>>
system::create_dihedrals(const Graph& g) const {
  using boost::adjacent_vertices;
  using boost::edges;
  using boost::source;
  using boost::target;

  std::vector<std::unique_ptr<internal_coord>> result;
  auto ed = edges(g);
  for (auto it = ed.first; it != ed.second; ++it) {
    auto u = source(*it, g);
    auto v = target(*it, g);
    auto u_vert = adjacent_vertices(u, g);
    for (auto u_vert_it = u_vert.first; u_vert_it != u_vert.second;
         ++u_vert_it) {
      auto v_vert = adjacent_vertices(v, g);
      for (auto v_vert_it = v_vert.first; v_vert_it != v_vert.second;
           ++v_vert_it) {
        if (g[*u_vert_it].atom_serial != g[*v_vert_it].atom_serial &&
            g[*u_vert_it].atom_serial != g[v].atom_serial &&
            g[u].atom_serial != g[*v_vert_it].atom_serial) {
          auto a_index = g[*u_vert_it].atom_serial;
          auto b_index = g[u].atom_serial;
          auto c_index = g[v].atom_serial;
          auto d_index = g[*v_vert_it].atom_serial;
          ;
          result.emplace_back(std::make_unique<dihedral>(a_index, b_index, c_index, d_index));
        }
      }
    }
  }
  return result;
}

// template <typename Graph>
// inline void
// system::create_ic_system(const std::vector<coords::Representation_3D>& res_vec,
//                          const coords::Representation_3D& coords,
//                          const Graph& g) {
//   trans_x_vec_.emplace_back(trans_x::create_trans(res_vec));
//   trans_y_vec_.emplace_back(trans_y::create_trans(res_vec));
//   trans_z_vec_.emplace_back(trans_z::create_trans(res_vec));
//   rotation_vec_ = create_rotations(res_vec);
//   distance_vec_ = create_distances(coords, g);
//   angle_vec_ = create_angles(coords, g);
//   oop_vec_ = create_oops(coords, g);
//   dihed_vec_ = create_dihedrals(coords, g);
// }

//template <typename Graph>
//inline void system::create_ic_system(const Graph& g) {
//  trans_x_vec_ = create_trans_x(res_vec_, res_index_vec_);
//  trans_y_vec_ = create_trans_y(res_vec_, res_index_vec_);
//  trans_z_vec_ = create_trans_z(res_vec_, res_index_vec_);
//  rotation_vec_ = create_rotations(res_vec_, res_index_vec_);
//  distance_vec_ = create_distances(rep_, g);
//  angle_vec_ = create_angles(rep_, g);
//  oop_vec_ = create_oops(rep_, g);
//  dihed_vec_ = create_dihedrals(rep_, g);
//}

inline std::tuple<std::vector<std::unique_ptr<internal_coord>>, std::vector<std::unique_ptr<internal_coord>>, std::vector<std::unique_ptr<internal_coord>>>
system::create_translations(const std::vector<std::vector<std::size_t>>& res_index_vec) const {
  return std::make_tuple(
    create_trans_x(res_index_vec),
    create_trans_y(res_index_vec),
    create_trans_z(res_index_vec)
  );
}

template <typename Graph>
inline void system::create_ic_system(const Graph& g) {
  append_primitives(create_distances(g));
  append_primitives(create_angles(g));
  append_primitives(create_oops(rep_, g));
  append_primitives(create_dihedrals(g));

  std::vector<std::unique_ptr<internal_coord>> trans_x, trans_y, trans_z;
  std::tie(trans_x, trans_y, trans_z) = create_translations(res_index_vec_);

  append_primitives(std::move(trans_x));
  append_primitives(std::move(trans_y));
  append_primitives(std::move(trans_z));

  rotation_vec_ = create_rotations(rep_, res_index_vec_);
}
//template <typename Graph>
//inline system::IC_System system::create_system(Graph const& g) {
//  IC_System new_ic_system;
//
//  new_ic_system.append_primitives(create_distances(g));
//  new_ic_system.append_primitives(create_angles(g));
//  new_ic_system.append_primitives(create_oops(rep_, g));
//  new_ic_system.append_primitives(create_dihedrals(g));
//
//  std::vector<intrtnal_coord> trans_x, trans_y, trans_z;
//  std::tie(trans_x, trans_y, trans_z) = create_translations(res_index_vec_);
//
//  new_ic_system.append_primitives(trans_x);
//  new_ic_system.append_primitives(trans_y);
//  new_ic_system.append_primitives(trans_z);
//
//  new_ic_system.configure_rotations(create_rotations(res_vec_, res_index_vec_));
//
//  return new_ic_system;
//}
}
namespace coords {
  class DL_Coordinates : public Coordinates {
  public:
    std::shared_ptr<coords::input::formats::pdb_helper::Parser<float_type>> parser;
    DL_Coordinates(Coordinates const& coords, std::unique_ptr<coords::input::formats::pdb> format) : Coordinates(coords) {
      if (!format) {
        throw std::runtime_error("You need to pass a format with fragments i. e. a pdb format.\n");
      }
      parser = format->parser;
    }
  };
}
#endif // cast_ic_core_h_guard
