/**
CAST 3
PrimitiveInternalCoordinates.h
Purpose: Definition of primitive Internal Coordinate Systems


@author Julian Erdmannsdörfer, Michael Prem
@version 3.0
*/

#ifndef PRIMITIVE_INTERNAL_COORDINATES_H
#define PRIMITIVE_INTERNAL_COORDINATES_H

#include"InternalCoordinates.h"
#include"coords.h"
#include "ic_core.h"
#include "graph.h"


namespace internals {
  class PrimitiveInternalCoordinates {
  protected:
    using CartesianType = InternalCoordinates::CartesiansForInternalCoordinates;
    using BondGraph = ic_util::Graph<ic_util::Node>;
    using InternalVec = std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>;
  public:
    PrimitiveInternalCoordinates(const std::vector<coords::Representation_3D>& res_init,
      const std::vector<std::vector<std::size_t>>& res_index,
      CartesianType & xyz_init, BondGraph const& graph)
      : res_vec_{ res_init }, subSystemIndices{ res_index } {
      create_ic_system(graph, xyz_init);
    }
    PrimitiveInternalCoordinates() = default;
    virtual ~PrimitiveInternalCoordinates() = default;

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> primitive_internals;
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> rotation_vec_;

    void requestNewBAndG(){
     new_B_matrix = true;
     new_G_matrix = true;
    }

  protected:

    const std::vector<coords::Representation_3D> res_vec_;
    const std::vector<std::vector<std::size_t>> subSystemIndices;
    //CartesianType xyz_;
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> registeredRotators;

    scon::mathmatrix<coords::float_type> B_matrix;
    scon::mathmatrix<coords::float_type> G_matrix;
    scon::mathmatrix<coords::float_type> hessian;

    std::vector<std::vector<coords::float_type>> deriv_vec(CartesianType const& cartesians);
    

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
      virtual ~InternalCoordinatesCreator() = default;
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
      create_rotations(CartesianType & cartesians);

    void create_ic_system(BondGraph const&, CartesianType &);

    virtual scon::mathmatrix<coords::float_type> calc(coords::Representation_3D const& xyz) const;//F
    virtual scon::mathmatrix<coords::float_type> calc_diff(coords::Representation_3D const& lhs, coords::Representation_3D const& rhs) const;//F

    virtual scon::mathmatrix<coords::float_type> guess_hessian(CartesianType const&) const;//F
    virtual scon::mathmatrix<coords::float_type>& Bmat(CartesianType const& cartesians);//F
    virtual scon::mathmatrix<coords::float_type>& Gmat(CartesianType const& cartesians);//F
    virtual scon::mathmatrix<coords::float_type> transposeOfBmat(CartesianType const& cartesian);
    virtual scon::mathmatrix<coords::float_type> pseudoInverseOfGmat(CartesianType const& cartesian);

  };


  class InternalToCartesianConverter {
  public:
    InternalToCartesianConverter(PrimitiveInternalCoordinates & internals,
      InternalCoordinates::CartesiansForInternalCoordinates & cartesians) : internalCoordinates{ internals }, cartesianCoordinates{ cartesians }, inverseHessian(0u,0u) {}

    scon::mathmatrix<coords::float_type> calculateInternalGradients(scon::mathmatrix<coords::float_type> const&);
    virtual scon::mathmatrix<coords::float_type> getInternalStep(scon::mathmatrix<coords::float_type> const&, scon::mathmatrix<coords::float_type> const&);

    virtual std::pair<coords::float_type, coords::float_type> getDeltaYPrimeAndSol(scon::mathmatrix<coords::float_type> const& internalStep, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian);
    
    virtual void applyInternalChange(scon::mathmatrix<coords::float_type>);//F
    template<typename XYZ>
    coords::Representation_3D& set_xyz(XYZ&& new_xyz);
    virtual InternalCoordinates::CartesiansForInternalCoordinates const& getCartesianCoordinates() const { return cartesianCoordinates; }
    virtual InternalCoordinates::CartesiansForInternalCoordinates & getCartesianCoordinates() { return cartesianCoordinates; }
    
    void invertNormalHessian(scon::mathmatrix<double> const& hessian);
  protected:
    PrimitiveInternalCoordinates & internalCoordinates;
    InternalCoordinates::CartesiansForInternalCoordinates & cartesianCoordinates;
    //write in another class so that it gets deleted after used by all functions
    scon::mathmatrix<coords::float_type> inverseHessian;
    coords::float_type getDeltaYPrime(scon::mathmatrix<coords::float_type> const& internalStep);
    coords::float_type InternalToCartesianConverter::getSol(scon::mathmatrix<coords::float_type> const& internalStep, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian);
  private:
    template<typename Dcart>
    coords::Representation_3D& takeCartesianStep(Dcart&& d_cart);
  };

  class StepRestrictor {
  public:
    StepRestrictor(InternalToCartesianConverter & converter, coords::float_type const target) : converter{ converter }, target{ target }, restrictedStep{}, restrictedSol{ 0.0 }, v0{0.0} {}
    coords::float_type operator()(scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const & hessian);
    scon::mathmatrix<coords::float_type> const& getRestrictedStep() const { return restrictedStep; }
    scon::mathmatrix<coords::float_type> & getRestrictedStep() { return restrictedStep; }
    void setInitialV0(coords::float_type const initialV0) { v0 = initialV0; }
  protected:
    scon::mathmatrix<coords::float_type> alterHessian(scon::mathmatrix<coords::float_type> const & hessian, coords::float_type const alteration) const;
    coords::float_type getStepNorm() const { return restrictedStep.norm(); }

    InternalToCartesianConverter & converter;
    coords::float_type target;
    scon::mathmatrix<coords::float_type> restrictedStep;
    coords::float_type restrictedSol, v0;
  };

  class InternalToCartesianStep {
  public:
    InternalToCartesianStep(InternalToCartesianConverter & converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian, coords::float_type const trustRadius) 
      : converter{ converter }, gradients{ gradients }, hessian{ hessian }, trustRadius{ trustRadius } {}
    coords::float_type operator()(coords::float_type const trial);
  protected:
    InternalToCartesianConverter & converter;
    scon::mathmatrix<coords::float_type> const& gradients;
    scon::mathmatrix<coords::float_type> const& hessian;
    coords::float_type trustRadius;
  };

  class AppropriateStepFinder {
  public:
    AppropriateStepFinder(){}
  };
  
  inline void InternalToCartesianConverter::applyInternalChange(scon::mathmatrix<coords::float_type> d_int_left) {
    using ic_util::flatten_c3_vec;

    auto old_xyz = cartesianCoordinates;
    coords::Representation_3D first_struct, last_good_xyz;
    auto micro_iter{ 0 }, fail_count{ 0 };
    auto damp{ 1. };
    auto old_inorm{ 0.0 };
    for (; micro_iter < 50; ++micro_iter) {
      takeCartesianStep(damp*internalCoordinates.transposeOfBmat(cartesianCoordinates)*internalCoordinates.pseudoInverseOfGmat(cartesianCoordinates)*d_int_left);

      auto d_now = internalCoordinates.calc_diff(cartesianCoordinates, old_xyz);

      auto d_int_remain = d_int_left - d_now;
      auto cartesian_rmsd = ic_util::Rep3D_to_Mat(old_xyz - cartesianCoordinates).rmsd();
      auto internal_norm = d_int_remain.norm();
      //std::cout << "Left change internal coordinates:\n" << d_int_remain << "\n\n";
      //std::cout << "internal norm: " << internal_norm << "\n\n";
      if (micro_iter == 0) {
        first_struct = cartesianCoordinates;
        last_good_xyz = cartesianCoordinates;
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
          old_inorm = internal_norm;
          last_good_xyz = cartesianCoordinates;
        }
      }
      if (cartesian_rmsd < 1.e-6 || internal_norm < 1.e-6) {
        std::cout << "Took " << micro_iter << " steps to converge.\n";
        return;
      }
      else if (fail_count >= 10) {
        std::cout << "Failed ten times to converge.\n";
        cartesianCoordinates = first_struct;
        return;
      }

      old_xyz = cartesianCoordinates;
      d_int_left = std::move(d_int_remain);
    }
    std::cout << "Took all " << micro_iter + 1 << " steps, still not converged.\n";
  }


  template<typename Dcart>
  coords::Representation_3D& InternalToCartesianConverter::takeCartesianStep(Dcart&& d_cart) {
    auto d_cart_rep3D = ic_util::mat_to_rep3D(std::forward<Dcart>(d_cart));
    return set_xyz(cartesianCoordinates + d_cart_rep3D);
  }

  template<typename XYZ>
  coords::Representation_3D& InternalToCartesianConverter::set_xyz(XYZ&& new_xyz) {
    internalCoordinates.requestNewBAndG();
    return cartesianCoordinates.setCartesianCoordnates(std::forward<XYZ>(new_xyz));
  }
}

#endif
