/*#include "PrimitiveInternalsTransRot.h"

namespace internals{
  
  PrimitiveInternalsTransRot::PrimitiveInternalsTransRot(
    const std::vector<coords::Representation_3D>& res_init,
    const std::vector<std::vector<std::size_t>>& res_index,
    CartesianType & xyz_init,
    BondGraph const& graph
  ) :
    PrimitiveInternalCoordinates{ res_init, res_index, xyz_init, graph }
  {
    create_translations_rotations(xyz_init);
  }


  void PrimitiveInternalsTransRot::create_translations_rotations(CartesianType & cartesians){
    append_primitives(create_trans_x());
    append_primitives(create_trans_y());
    append_primitives(create_trans_z());
    
    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> rotationA, rotationB, rotationC;
    std::tie(rotationA, rotationB, rotationC) = create_rotations(cartesians);

    append_primitives(std::move(rotationA));
    append_primitives(std::move(rotationB));
    append_primitives(std::move(rotationC));
  }
  
  std::shared_ptr<InternalCoordinates::Rotator> PrimitiveInternalsTransRot::build_rotation(
    InternalCoordinates::CartesiansForInternalCoordinates & target,
    std::vector<std::size_t> const & index_vec
  ) {
    coords::Representation_3D reference;
    for (auto const & ind : index_vec) {
      reference.emplace_back(target.at(ind-1));
    }
    return InternalCoordinates::Rotator::buildRotator(target, index_vec);
  }
  
  PrimitiveInternalCoordinates::InternalVec PrimitiveInternalsTransRot::create_trans_x() const {

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    for (auto const & indices : subSystemIndices) {
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationX>(indices));
    }
    return result;
  }

  PrimitiveInternalCoordinates::InternalVec PrimitiveInternalsTransRot::create_trans_y() const {

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    for (auto const & indices : subSystemIndices) {
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationY>(indices));
    }
    return result;
  }

  PrimitiveInternalCoordinates::InternalVec PrimitiveInternalsTransRot::create_trans_z() const {

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    for (auto const & indices : subSystemIndices) {
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationZ>(indices));
    }
    return result;
  }

  /*std::tuple<PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec>
    PrimitiveInternalsTransRot::create_translations() const {
    return std::make_tuple(
      create_trans_x(),
      create_trans_y(),
      create_trans_z()
    );
  }*/
  
  /*std::tuple<PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec> PrimitiveInternalsTransRot::create_rotations(CartesianType & cartesians) {
    std::vector<InternalCoordinates::Rotations> result;
    for (auto const& indices : subSystemIndices) {
      result.emplace_back(build_rotation(cartesians, indices)->makeRotations());
    }
    //return createRotationABC(result);
    
    InternalVec resultA, resultB, resultC;
    for (auto & rotation : result) {
      resultA.emplace_back(std::move(rotation.rotationA));
      resultB.emplace_back(std::move(rotation.rotationB));
      resultC.emplace_back(std::move(rotation.rotationC));
      registeredRotators.emplace_back(rotation.rotator);
    }
    return { std::move(resultA), std::move(resultB), std::move(resultC) };
  }
  
  /*scon::mathmatrix<coords::float_type> PrimitiveInternalsTransRot::calc_diff(coords::Representation_3D const & lhs, coords::Representation_3D const & rhs) const {
    //TODO remove these from here
    
    auto lprims = PrimitiveInternalCoordinates::calc(lhs);
    //TODO remove these from here
    for (auto & r : registeredRotators) {
      r->requestNewValueEvaluation();
    }
    auto rprims = PrimitiveInternalCoordinates::calc(rhs);
    auto diff = lprims - rprims;

    for (auto i = 0u; i < primitive_internals.size(); ++i) {
      if (dynamic_cast<InternalCoordinates::DihedralAngle*>(primitive_internals.at(i).get())) {
        if (std::fabs(diff(0, i)) > SCON_PI) {
          if (diff(0, i) < 0.0) {
            diff(0, i) += 2.*SCON_PI;
          }
          else {
            diff(0, i) -= 2.*SCON_PI;
          }
        }
      }
    }
    //std::cout << "Diff:\n" << diff.t() << "\n";
    return diff;
  }*/
  
  /*void PrimitiveInternalsTransRot::prepare_rotations() const{
    for (auto & r : registeredRotators) {
      r->requestNewValueEvaluation();
    }
  }
}*/
