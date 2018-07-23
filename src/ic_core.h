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
#include <memory>

#include"InternalCoordinates.h"

#include "coords_io_pdb.h"

namespace coords {
  template<typename CoordKind>
  class DL_Coordinates : public Coordinates {
  public:
    std::shared_ptr<typename CoordKind::helper::template Parser<double>> parser;
    DL_Coordinates(Coordinates const& coords, std::unique_ptr<CoordKind> format) : Coordinates(coords) {
      if (!format) {
        throw std::runtime_error("You need to pass a format with fragments i. e. a pdb format.\n");
      }
      parser = format->parser;
    }
  };
}
namespace ic_core {

using coords::float_type;

coords::Representation_3D grads_to_bohr(coords::Representation_3D const& grads);
coords::Representation_3D rep3d_bohr_to_ang(coords::Representation_3D const& bohr);

//TODO Replace target as a instance of CartesiansForInternalCoordinates
std::shared_ptr<InternalCoordinates::Rotator> build_rotation(const coords::Representation_3D& target,
  const std::vector<std::size_t>& index_vec);

class system {
public:
  system(const std::vector<coords::Representation_3D>& res_init,
         const std::vector<std::vector<std::size_t>>& res_index,
         const coords::Representation_3D& xyz_init)
      : res_vec_{ res_init }, res_index_vec_{ res_index }, xyz_{ xyz_init } {}

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> primitive_internals;
  std::vector<std::shared_ptr<InternalCoordinates::Rotator>> rotation_vec_;

private:
  
  const std::vector<coords::Representation_3D> res_vec_;
  const std::vector<std::vector<std::size_t>> res_index_vec_;
  coords::Representation_3D xyz_;
  scon::mathmatrix<float_type> B_matrix;
  scon::mathmatrix<float_type> G_matrix;
  scon::mathmatrix<float_type> del_mat;
  scon::mathmatrix<float_type> hessian;

  template<typename VertIter>
  static std::vector<std::vector<std::size_t>> possible_sets_of_3(VertIter const& vbegin, VertIter const& vend);

  bool new_B_matrix = true;
  bool new_G_matrix = true;


public:
  void append_primitives(std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> && pic) {
    primitive_internals.insert(primitive_internals.end(),
      std::make_move_iterator(pic.begin()),
      std::make_move_iterator(pic.end()));
  }

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
  create_trans_x(const std::vector<std::vector<std::size_t>>&) const;

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
  create_trans_y(const std::vector<std::vector<std::size_t>>&) const;

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
  create_trans_z(const std::vector<std::vector<std::size_t>>&) const;

  std::tuple<std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>, std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>, std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>>
    create_translations(const std::vector<std::vector<std::size_t>>&) const;

  std::vector<std::shared_ptr<InternalCoordinates::Rotator>>
  create_rotations(const coords::Representation_3D&,
                   const std::vector<std::vector<std::size_t>>&);

  /*template<typename Graph>
  IC_System create_system(Graph const &);*/

  template <typename Graph>
  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> create_distances(const Graph&) const;

  template <typename Graph>
  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> create_angles(const Graph&) const;

  template <typename Graph>
  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> create_oops(const coords::Representation_3D&,
                                        const Graph&) const;

  template <typename Graph>
  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> create_dihedrals(const Graph&) const;

  /*template <typename Graph>
  void create_ic_system(const std::vector<coords::Representation_3D>&,
                        const coords::Representation_3D&, const Graph&);*/

  template <typename Graph>
  void create_ic_system(const Graph&);

  scon::mathmatrix<float_type>& delocalize_ic_system();
  scon::mathmatrix<float_type> guess_hessian();
  scon::mathmatrix<float_type>& initial_delocalized_hessian();
  std::vector<std::vector<float_type>> deriv_vec();
  scon::mathmatrix<float_type>& Bmat();
  scon::mathmatrix<float_type>& Gmat();
  scon::mathmatrix<float_type>& ic_Bmat();
  scon::mathmatrix<float_type>& ic_Gmat();


  scon::mathmatrix<float_type> calc() const;
  scon::mathmatrix<float_type> calculate_internal_grads(scon::mathmatrix<float_type> const&);
  template<typename Gint>
  scon::mathmatrix<float_type> get_internal_step(Gint&&);
  template<typename Dint>
  void apply_internal_change(Dint&&);
  template<typename Dcart>
  coords::Representation_3D& take_Cartesian_step(Dcart&& d_cart);
  template<typename XYZ>
  coords::Representation_3D& set_xyz(XYZ&& new_xyz);
  
  scon::mathmatrix<float_type> calc_prims(coords::Representation_3D const& xyz) const;

  template<typename XYZ>
  scon::mathmatrix<float_type> calc(XYZ&& xyz) const;
  template<typename XYZ>
  scon::mathmatrix<float_type> calc_diff(XYZ&& lhs, XYZ&& rhs) const;

  void optimize(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
  
private:
    void initializeOptimization(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
    void setCartesianCoordinatesForGradientCalculation(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
    void prepareOldVariablesPtr(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
    void evaluateNewCartesianStructure(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
    void applyHessianChange();
    void setNewToOldVariables();
    scon::mathmatrix<float_type> getInternalGradientsButReturnCartesianOnes(coords::DL_Coordinates<coords::input::formats::pdb> & coords);
    
    
    class ConvergenceCheck{
    public:
        ConvergenceCheck(int step, scon::mathmatrix<float_type> & gxyz, system const& sys) :
            step{step},
            cartesianGradients{gxyz},
            internalCoordinateSystem{sys},
            energyDiff{0.0},
            gradientRms{0.0},
            displacementRms{0.0},
            gradientMax{0.0},
            displacementMax{0.0}{}
        
        void writeAndCalcEnergyDiffs();
        void writeAndCalcGradientRmsd();
        void writeAndCalcDisplacementRmsd();
        bool checkConvergence()const;
        bool operator()();
    private:
        int step;
        scon::mathmatrix<float_type> & cartesianGradients;
        system const& internalCoordinateSystem;
        
        static auto constexpr threshEnergy = 1.e-6;
        static auto constexpr threshGradientRms = 0.0003;
        static auto constexpr threshDisplacementRms = 0.00045;
        static auto constexpr threshGradientMax = 0.0012;
        static auto constexpr threshDisplacementMax = 0.0018;
        
        float_type energyDiff;
        float_type gradientRms;
        float_type displacementRms;
        float_type gradientMax;
        float_type displacementMax;
    };
    
    static std::pair<float_type,float_type> gradientRmsValAndMax(scon::mathmatrix<float_type> const& grads);
    std::pair<float_type,float_type> displacementRmsValAndMax()const;
    
    static std::pair<float_type,float_type> displacementRmsValAndMaxTwoStructures(coords::Representation_3D const& oldXyz, coords::Representation_3D const& newXyz);
    
    struct SystemVariables{
        float_type systemEnergy;
        scon::mathmatrix<coords::float_type> systemGradients;
        coords::Representation_3D systemCartesianRepresentation;
    };
    
    SystemVariables currentVariables;
    std::unique_ptr<SystemVariables> oldVariables;

};

template <typename Graph>
inline std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
system::create_distances(const Graph& g) const {
  using boost::edges;
  using boost::source;
  using boost::target;

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
  auto ed = edges(g);
  for (auto it = ed.first; it != ed.second; ++it) {
    auto u = source(*it, g);
    auto v = target(*it, g);
    auto u_index = g[u].atom_serial;
    auto v_index = g[v].atom_serial;
    auto u_elem = g[u].element;
    auto v_elem = g[v].element;
    result.emplace_back(std::make_unique<InternalCoordinates::BondDistance>(u_index, v_index, u_elem, v_elem));
  }
  return result;
}

template <typename Graph>
inline std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
system::create_angles(const Graph& g) const {
  using boost::adjacent_vertices;
  using boost::vertices;

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
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
          result.emplace_back(std::make_unique<InternalCoordinates::BondAngle>(a_index, b_index, c_index, a_elem, b_elem, c_elem));
        }
      }
    }
  }
  return result;
}

template<typename VertIter>
inline std::vector<std::vector<std::size_t>> system::possible_sets_of_3(VertIter const& vbegin, VertIter const& vend){
  std::vector<std::vector<std::size_t>> result;
  for(auto first = vbegin; first < vend-2; ++first){
    for(auto second = first+1; second< vend-1; ++second){
      for(auto third = first+2; third < vend; ++third){
        result.emplace_back(std::vector<std::size_t>{*first, *second, *third});
      }
    }
  }
  return result;
}

template <typename Graph>
inline std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
system::create_oops(const coords::Representation_3D& coords, const Graph& g) const {
  using boost::adjacent_vertices;
  using boost::vertices;
  using scon::dot;

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
  auto vert = vertices(g);
  for (auto it = vert.first; it != vert.second; ++it) {
    auto core = g[*it].atom_serial;
    auto core_cp = coords.at(core - 1);
    auto vert_i = adjacent_vertices(*it, g);
    if((vert_i.second - vert_i.first) >= 3){
      auto permutations = possible_sets_of_3(vert_i.first, vert_i.second);
      for(auto & combination : permutations){
        std::sort(combination.begin(), combination.end());


      //auto permutation_vec = ic_util::permutation_from_vec(neighbours);
      //for (auto& permutation : permutation_vec) {
        auto u_cp = coords.at(combination.at(0));
        auto v_cp = coords.at(combination.at(1));
        auto w_cp = coords.at(combination.at(2));
        auto n_vec1 = ic_util::normal_unit_vector(u_cp, v_cp, w_cp);
        auto n_vec2 = ic_util::normal_unit_vector(core_cp, u_cp, v_cp);
        auto dot_n_vecs = dot(n_vec1, n_vec2);
        if (0.95 < std::fabs(dot_n_vecs)) {
          result.emplace_back(std::make_unique<InternalCoordinates::OutOfPlane>(core, combination.at(0)+1, combination.at(1)+1, combination.at(2)+1));
        }
      }
      //}
    }
  }
  return result;
}

template <typename Graph>
inline std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
system::create_dihedrals(const Graph& g) const {
  using boost::adjacent_vertices;
  using boost::edges;
  using boost::source;
  using boost::target;

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
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
          result.emplace_back(std::make_unique<InternalCoordinates::DihedralAngle>(a_index, b_index, c_index, d_index));
        }
      }
    }
  }
  return result;
}

inline std::tuple<std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>, std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>, std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>>
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
  append_primitives(create_oops(xyz_, g));
  append_primitives(create_dihedrals(g));

  std::vector<std::unique_ptr<InternalCoordinateImpl>> trans_x, trans_y, trans_z;
  std::tie(trans_x, trans_y, trans_z) = create_translations(res_index_vec_);

  append_primitives(std::move(trans_x));
  append_primitives(std::move(trans_y));
  append_primitives(std::move(trans_z));

  rotation_vec_ = create_rotations(xyz_, res_index_vec_);
}

template<typename Gint>
inline scon::mathmatrix<float_type> system::get_internal_step(Gint&& g_int){
  return -1.*hessian.pinv()*g_int;
}

template<typename Dint>
inline void system::apply_internal_change(Dint&& d_int){
  using ic_util::flatten_c3_vec;

  auto old_xyz = xyz_;
  coords::Representation_3D first_struct, last_good_xyz;
  auto d_int_left = std::forward<Dint>(d_int);
  auto micro_iter{ 0 }, fail_count{ 0 };
  auto damp{1.};
  auto old_rmsd{ 0.0 }, old_inorm{ 0.0 };
  for(; micro_iter < 50; ++micro_iter){

    //std::cout << "Bmat:\n" << ic_Bmat().t() << "\n\n";
    //std::cout << "Ginv:\n" << ic_Gmat().pinv() << "\n\n";

    take_Cartesian_step(damp*ic_Bmat().t()*ic_Gmat().pinv()*d_int_left); //should it not be G^-1*B^T?
    //std::cout << "Cartesian:\n" << xyz_ << std::endl;

    auto d_now = calc_diff(xyz_, old_xyz);
    //std::cout << "Diff internal coordinates:\n" << d_now << std::endl;

    auto d_int_remain = d_int_left - d_now;
    auto cartesian_rmsd = ic_util::Rep3D_to_Mat(old_xyz - xyz_).rmsd();
    auto internal_norm = d_int_remain.norm();
    //std::cout << "Left change internal coordinates:\n" << d_int_remain << "\n\n";
    //std::cout << "internal norm: " << internal_norm << "\n\n";
    if (micro_iter == 0){
      first_struct = xyz_;
      last_good_xyz = xyz_;
      old_rmsd = cartesian_rmsd;
      old_inorm = internal_norm;
    }
    else {
      if(internal_norm>old_inorm){
        damp /=2.;
        ++fail_count;
      }
      else{
        fail_count = 0;
        damp = std::min(1.2*damp,1.);
        old_rmsd = cartesian_rmsd;
        old_inorm = internal_norm;
        last_good_xyz = xyz_;
      }
    }
    if(cartesian_rmsd < 1.e-6 || internal_norm < 1.e-6){
      std::cout << "Took " << micro_iter << " steps to converge.\n";
      return;
    }
    else if(fail_count>=10){
        std::cout << "Failed ten times to converge.\n";
        xyz_ = first_struct;
        return;
    }

    old_xyz = xyz_;
    d_int_left = std::move(d_int_remain);
  }
  std::cout << "Took all " << micro_iter+1 << " steps, still not converged.\n";
}

template<typename Dcart>
coords::Representation_3D& system::take_Cartesian_step(Dcart&& d_cart){
  auto d_cart_rep3D = ic_util::mat_to_rep3D(std::forward<Dcart>(d_cart));
  return set_xyz(xyz_ + d_cart_rep3D);
}

template<typename XYZ>
coords::Representation_3D& system::set_xyz(XYZ&& new_xyz){
  new_B_matrix = true;
  new_G_matrix = true;
  return xyz_ = std::forward<XYZ>(new_xyz);
}

template<typename XYZ>
scon::mathmatrix<float_type> system::calc(XYZ&& xyz) const{
  auto prims = calc_prims(std::forward<XYZ>(xyz));
  return (prims * del_mat).t();
}

template<typename XYZ>
scon::mathmatrix<float_type> system::calc_diff(XYZ&& lhs, XYZ&& rhs) const{
  auto lprims = calc_prims(std::forward<XYZ>(lhs));
  auto rprims = calc_prims(std::forward<XYZ>(rhs));
  auto diff = lprims - rprims;

  for(auto i = 0u;i<primitive_internals.size(); ++i){
    if(dynamic_cast<InternalCoordinates::DihedralAngle*>(primitive_internals.at(i).get())){
      if(std::fabs(diff(0, i)) > SCON_PI){
        if(diff(0, i)<0.0){
          diff(0, i) += 2.*SCON_PI;
        }
        else{
          diff(0, i) -= 2.*SCON_PI;
        }
      }
    }
  }
  //std::cout << "Diff:\n" << diff.t() << "\n";
  return (diff * del_mat).t();
}
}
#endif // cast_ic_core_h_guard
