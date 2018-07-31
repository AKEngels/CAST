#ifndef PRIMITIVE_INTERNAL_COORDINATES_H
#define PRIMITIVE_INTERNAL_COORDINATES_H

#include"InternalCoordinates.h"
#include"coords.h"
#include "ic_core.h"
#include "graph.h"


namespace internals {
  class system {
  protected:
    using CartesianType = InternalCoordinates::CartesiansForInternalCoordinates;
    using BondGraph = ic_util::Graph<ic_util::Node>;
    using InternalVec = std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>;
  public:
    system(const std::vector<coords::Representation_3D>& res_init,
      const std::vector<std::vector<std::size_t>>& res_index,
      CartesianType const& xyz_init)
      : res_vec_{ res_init }, subSystemIndices{ res_index }, xyz_{ xyz_init } {}

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> primitive_internals;
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> rotation_vec_;

  protected:

    const std::vector<coords::Representation_3D> res_vec_;
    const std::vector<std::vector<std::size_t>> subSystemIndices;
    CartesianType xyz_;
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> registeredRotators;

    scon::mathmatrix<coords::float_type> B_matrix;
    scon::mathmatrix<coords::float_type> G_matrix;
    scon::mathmatrix<coords::float_type> hessian;

    std::vector<std::vector<coords::float_type>> deriv_vec();
    template<typename XYZ>
    coords::Representation_3D& set_xyz(XYZ&& new_xyz);

    static std::vector<std::vector<std::size_t>> possible_sets_of_3(BondGraph::adjacency_iterator const vbegin, BondGraph::adjacency_iterator const vend);
    std::shared_ptr<InternalCoordinates::Rotator> build_rotation(InternalCoordinates::CartesiansForInternalCoordinates & target,
      std::vector<std::size_t> const& index_vec);

    bool new_B_matrix = true;
    bool new_G_matrix = true;

    class InternalCoordinatesCreator {
    protected:
      using BondGraph = ic_util::Graph<ic_util::Node>;
    public:
      InternalCoordinatesCreator(BondGraph const& graph) : bondGraph{ graph } {}
      virtual std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> getInternals() = 0;
    protected:
      ic_util::Graph<ic_util::Node> const& bondGraph;
    };

    class DistanceCreator : public InternalCoordinatesCreator {
    public:
      DistanceCreator(BondGraph const& graph) : InternalCoordinatesCreator{ graph }, source{ 0u }, target{ 0u }, edgeIterators{ boost::edges(bondGraph) } {}

      virtual std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> getInternals() override;

    protected:
      bool nextEdgeDistances();
      std::size_t source, target;
      std::pair<ic_util::Graph<ic_util::Node>::edge_iterator, ic_util::Graph<ic_util::Node>::edge_iterator> edgeIterators;
    };

    class AngleCreator : public InternalCoordinatesCreator {
    public:
      AngleCreator(BondGraph const& graph) : InternalCoordinatesCreator{ graph }, leftAtom{ 0u }, middleAtom{ 0u }, rightAtom{ 0u }, vertexIterators{ boost::vertices(graph) } {}
      virtual std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> getInternals() override;
    protected:
      bool nextVertex();
      void addAngleForAllNeighbors();
      void spanLeftAndRightNeighborsForAngle(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & neighbors);
      bool findLeftAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & neighbors);

      bool findRightAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & neighborsLeft);

      std::size_t leftAtom, middleAtom, rightAtom;
      std::pair<ic_util::Graph<ic_util::Node>::vertex_iterator, ic_util::Graph<ic_util::Node>::vertex_iterator> vertexIterators;
      std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> * pointerToResult;
    };

    class DihedralCreator : public DistanceCreator {
    public:
      DihedralCreator(BondGraph const& graph) : DistanceCreator{ graph }, outerLeft{ 0u }, outerRight{ 0u } {}
      virtual std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> getInternals() override;
    protected:
      void findLeftAndRightAtoms();
      bool findLeftAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & sourceNeighbors);
      bool findRightAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator> & targetNeighbors);

      std::size_t outerLeft, outerRight;

      std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> * pointerToResult;
    };

  public:
    void append_primitives(InternalVec && pic);

    CartesianType const& getXyz() const { return xyz_; }

    InternalVec create_distances(BondGraph const&) const;
    InternalVec create_angles(BondGraph const&) const;
    InternalVec create_oops(coords::Representation_3D const&, BondGraph const&) const;
    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> create_dihedrals(BondGraph const&) const;

    InternalVec create_trans_x() const;
    InternalVec create_trans_y() const;
    InternalVec create_trans_z() const;
    std::tuple<InternalVec, InternalVec, InternalVec>
      create_translations()const;

    std::tuple<InternalVec, InternalVec, InternalVec>
      createRotationABC(std::vector<InternalCoordinates::Rotations> & rotations);
    std::tuple<InternalVec, InternalVec, InternalVec>
      create_rotations();

    void create_ic_system(BondGraph const&);

    virtual scon::mathmatrix<coords::float_type> calc(coords::Representation_3D const& xyz) const;//F
    virtual scon::mathmatrix<coords::float_type> calc_diff(coords::Representation_3D const& lhs, coords::Representation_3D const& rhs) const;//F

    virtual scon::mathmatrix<coords::float_type>& guess_hessian();//F
    virtual scon::mathmatrix<coords::float_type>& Bmat();//F
    virtual scon::mathmatrix<coords::float_type>& Gmat();//F

    scon::mathmatrix<coords::float_type>& getHessian() { return hessian; }
    scon::mathmatrix<coords::float_type> const& getHessian() const { return hessian; }; 
    template<typename Hessian>
    /*typename std::enable_if<std::is_same<Hessian, scon::Mathmatrix<coords::float_type>::value>::type*/ void
      setHessian(Hessian && newHessian) { hessian std::forward<Hessian>(newHessian); }

    scon::mathmatrix<coords::float_type> calculate_internal_grads(scon::mathmatrix<coords::float_type> const&);//F
    template<typename Gint>
    scon::mathmatrix<coords::float_type> get_internal_step(Gint&&);//F
    template<typename Dint>
    void apply_internal_change(Dint&&);//F
    template<typename Dcart>
    coords::Representation_3D& take_Cartesian_step(Dcart&& d_cart);
  };

  template<typename XYZ>
  coords::Representation_3D& system::set_xyz(XYZ&& new_xyz) {
    new_B_matrix = true;
    new_G_matrix = true;
    return xyz_.setCartesianCoordnates(std::forward<XYZ>(new_xyz));
  }

  template<typename Gint>
  inline scon::mathmatrix<coords::float_type> system::get_internal_step(Gint&& g_int) {
    return -1.*hessian.pinv()*g_int;
  }

  template<typename Dint>
  inline void system::apply_internal_change(Dint&& d_int) {
    using ic_util::flatten_c3_vec;

    auto old_xyz = xyz_;
    coords::Representation_3D first_struct, last_good_xyz;
    auto d_int_left = std::forward<Dint>(d_int);
    auto micro_iter{ 0 }, fail_count{ 0 };
    auto damp{ 1. };
    auto old_rmsd{ 0.0 }, old_inorm{ 0.0 };
    for (; micro_iter < 50; ++micro_iter) {

      take_Cartesian_step(damp*Bmat().t()*Gmat().pinv()*d_int_left); //should it not be G^-1*B^T?
                                                                     //std::cout << "Cartesian:\n" << xyz_ << std::endl;

      auto d_now = calc_diff(xyz_, old_xyz);
      //std::cout << "Diff internal coordinates:\n" << d_now << std::endl;

      auto d_int_remain = d_int_left - d_now;
      auto cartesian_rmsd = ic_util::Rep3D_to_Mat(old_xyz - xyz_).rmsd();
      auto internal_norm = d_int_remain.norm();
      //std::cout << "Left change internal coordinates:\n" << d_int_remain << "\n\n";
      //std::cout << "internal norm: " << internal_norm << "\n\n";
      if (micro_iter == 0) {
        first_struct = xyz_;
        last_good_xyz = xyz_;
        old_rmsd = cartesian_rmsd;
        old_inorm = internal_norm;
      }
      else {
        if (internal_norm > old_inorm) {
          damp /= 2.;
          ++fail_count;
        }
        else {
          fail_count = 0;
          damp = std::min(1.2*damp, 1.);
          old_rmsd = cartesian_rmsd;
          old_inorm = internal_norm;
          last_good_xyz = xyz_;
        }
      }
      if (cartesian_rmsd < 1.e-6 || internal_norm < 1.e-6) {
        std::cout << "Took " << micro_iter << " steps to converge.\n";
        return;
      }
      else if (fail_count >= 10) {
        std::cout << "Failed ten times to converge.\n";
        xyz_ = first_struct;
        return;
      }

      old_xyz = xyz_;
      d_int_left = std::move(d_int_remain);
    }
    std::cout << "Took all " << micro_iter + 1 << " steps, still not converged.\n";
  }

  template<typename Dcart>
  coords::Representation_3D& system::take_Cartesian_step(Dcart&& d_cart) {
    auto d_cart_rep3D = ic_util::mat_to_rep3D(std::forward<Dcart>(d_cart));
    return set_xyz(xyz_ + d_cart_rep3D);
  }
}

#endif