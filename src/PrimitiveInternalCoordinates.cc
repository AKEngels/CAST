#include "PrimitiveInternalCoordinates.h"
namespace internals {
  std::vector<std::vector<std::size_t>> PrimitiveInternalCoordinates::possible_sets_of_3(BondGraph::adjacency_iterator const vbegin, BondGraph::adjacency_iterator const vend) {
    std::vector<std::vector<std::size_t>> result;
    for (auto first = vbegin; first < vend - 2; ++first) {
      for (auto second = first + 1; second < vend - 1; ++second) {
        for (auto third = first + 2; third < vend; ++third) {
          result.emplace_back(std::vector<std::size_t>{*first, *second, *third});
        }
      }
    }
    return result;
  }

  std::shared_ptr<InternalCoordinates::Rotator> PrimitiveInternalCoordinates::build_rotation(InternalCoordinates::CartesiansForInternalCoordinates & target, std::vector<std::size_t> const & index_vec) {
    coords::Representation_3D reference;
    for (auto const & ind : index_vec) {
      reference.emplace_back(target.at(ind));
    }
    return InternalCoordinates::Rotator::buildRotator(target, index_vec);
  }

  void PrimitiveInternalCoordinates::append_primitives(PrimitiveInternalCoordinates::InternalVec&& pic) {
    primitive_internals.insert(primitive_internals.end(),
      std::make_move_iterator(pic.begin()),
      std::make_move_iterator(pic.end()));
  }

  PrimitiveInternalCoordinates::InternalVec PrimitiveInternalCoordinates::create_trans_x() const {

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    for (auto const & indices : subSystemIndices) {
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationX>(indices));
    }
    return result;
  }

  PrimitiveInternalCoordinates::InternalVec PrimitiveInternalCoordinates::create_trans_y() const {

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    for (auto const & indices : subSystemIndices) {
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationY>(indices));
    }
    return result;
  }

  PrimitiveInternalCoordinates::InternalVec PrimitiveInternalCoordinates::create_trans_z() const {

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    for (auto const & indices : subSystemIndices) {
      result.emplace_back(std::make_unique<InternalCoordinates::TranslationZ>(indices));
    }
    return result;
  }

  std::tuple<PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec>
    PrimitiveInternalCoordinates::create_translations() const {
    return std::make_tuple(
      create_trans_x(),
      create_trans_y(),
      create_trans_z()
    );
  }

  std::tuple<PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec>
    PrimitiveInternalCoordinates::createRotationABC(std::vector<InternalCoordinates::Rotations> & rotations) {
    InternalVec resultA, resultB, resultC;
    for (auto & rotation : rotations) {
      resultA.emplace_back(std::move(rotation.rotationA));
      resultB.emplace_back(std::move(rotation.rotationB));
      resultC.emplace_back(std::move(rotation.rotationC));
      registeredRotators.emplace_back(rotation.rotator);
    }
    return { std::move(resultA), std::move(resultB), std::move(resultC) };
  }

  std::tuple<PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec> PrimitiveInternalCoordinates::create_rotations() {
    std::vector<InternalCoordinates::Rotations> result;
    for (auto const& indices : subSystemIndices) {
      result.emplace_back(build_rotation(xyz_, indices)->makeRotations());
    }
    return createRotationABC(result);
  }

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> PrimitiveInternalCoordinates::DistanceCreator::getInternals() {
    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    while (nextEdgeDistances()) {
      result.emplace_back(std::make_unique<InternalCoordinates::BondDistance>(bondGraph[source], bondGraph[target]));
    }
    return result;
  }

  bool PrimitiveInternalCoordinates::DistanceCreator::nextEdgeDistances() {
    if (edgeIterators.first == edgeIterators.second) return false;
    source = boost::source(*edgeIterators.first, bondGraph);
    target = boost::target(*edgeIterators.first, bondGraph);
    ++edgeIterators.first;
    return true;
  }

  inline std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
    PrimitiveInternalCoordinates::create_distances(const BondGraph& g) const {

    DistanceCreator dc(g);
    return dc.getInternals();
  }

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> PrimitiveInternalCoordinates::AngleCreator::getInternals() {
    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    pointerToResult = &result;
    while (nextVertex()) {
      addAngleForAllNeighbors();
    }
    return result;
  }

  bool PrimitiveInternalCoordinates::AngleCreator::nextVertex() {
    if (vertexIterators.first == vertexIterators.second) return false;
    middleAtom = *vertexIterators.first;
    vertexIterators.first++;
    return true;
  }

  void PrimitiveInternalCoordinates::AngleCreator::addAngleForAllNeighbors() {
    auto allNeighbors = boost::adjacent_vertices(middleAtom, bondGraph);
    spanLeftAndRightNeighborsForAngle(allNeighbors);
  }

  void PrimitiveInternalCoordinates::AngleCreator::spanLeftAndRightNeighborsForAngle(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& neighbors) {
    while (findLeftAtom(neighbors)) {
      auto copyOfLeftNeighbors = neighbors;
      while (findRightAtom(copyOfLeftNeighbors)) {
        pointerToResult->emplace_back(std::make_unique<InternalCoordinates::BondAngle>(
          bondGraph[leftAtom], bondGraph[middleAtom], bondGraph[rightAtom]));
      }
    }
  }

  bool PrimitiveInternalCoordinates::AngleCreator::findLeftAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& neighbors) {
    if (neighbors.first == neighbors.second) return false;
    leftAtom = *neighbors.first;
    ++neighbors.first;
    return true;
  }

  bool PrimitiveInternalCoordinates::AngleCreator::findRightAtom(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& neighborsLeft) {
    if (neighborsLeft.first == neighborsLeft.second) return false;
    rightAtom = *neighborsLeft.first;
    ++neighborsLeft.first;
    return true;
  }

  inline std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
    PrimitiveInternalCoordinates::create_angles(const BondGraph& g) const {
    AngleCreator ac(g);
    return ac.getInternals();
  }

  //This function surely does not work.
  inline std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
    PrimitiveInternalCoordinates::create_oops(const coords::Representation_3D& coords, const BondGraph& g) const {
    using boost::adjacent_vertices;
    using boost::vertices;
    using scon::dot;

    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    auto vert = vertices(g);
    for (auto it = vert.first; it != vert.second; ++it) {
      auto core = g[*it].atom_serial;
      auto core_cp = coords.at(core - 1);
      auto vert_i = adjacent_vertices(*it, g);
      if ((vert_i.second - vert_i.first) >= 3) {
        auto permutations = possible_sets_of_3(vert_i.first, vert_i.second);
        for (auto & combination : permutations) {
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
            result.emplace_back(std::make_unique<InternalCoordinates::OutOfPlane>(g[*it], g[combination.at(0)], g[combination.at(1)], g[combination.at(2)]));
          }
        }
        //}
      }
    }
    return result;
  }

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> PrimitiveInternalCoordinates::DihedralCreator::getInternals() {
    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> result;
    pointerToResult = &result;
    while (nextEdgeDistances()) {
      findLeftAndRightAtoms();
    }
    return result;
  }

  void PrimitiveInternalCoordinates::DihedralCreator::findLeftAndRightAtoms() {
    auto leftVertices = boost::adjacent_vertices(source, bondGraph);
    while (findLeftAtoms(leftVertices)) {
      auto rightVertices = boost::adjacent_vertices(target, bondGraph);
      while (findRightAtoms(rightVertices)) {
        pointerToResult->emplace_back(std::make_unique<InternalCoordinates::DihedralAngle>(
          bondGraph[outerLeft], bondGraph[source], bondGraph[target], bondGraph[outerRight]));
      }
    }
  }

  bool PrimitiveInternalCoordinates::DihedralCreator::findLeftAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& sourceNeighbors) {
    if (sourceNeighbors.first == sourceNeighbors.second) return false;
    outerLeft = *sourceNeighbors.first;
    ++sourceNeighbors.first;
    if (outerLeft == target) return findLeftAtoms(sourceNeighbors);
    return true;
  }

  bool PrimitiveInternalCoordinates::DihedralCreator::findRightAtoms(std::pair<BondGraph::adjacency_iterator, BondGraph::adjacency_iterator>& targetNeighbors) {
    if (targetNeighbors.first == targetNeighbors.second) return false;
    outerRight = *targetNeighbors.first;
    ++targetNeighbors.first;
    if (outerRight == source) return findRightAtoms(targetNeighbors);
    return true;
  }

  inline std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
    PrimitiveInternalCoordinates::create_dihedrals(const BondGraph& g) const {
    DihedralCreator dc(g);
    return dc.getInternals();
  }

  void PrimitiveInternalCoordinates::create_ic_system(BondGraph const& g) {
    append_primitives(create_distances(g));
    append_primitives(create_angles(g));
    append_primitives(create_oops(xyz_, g));
    append_primitives(create_dihedrals(g));

    //TODO own function for translations
    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> trans_x, trans_y, trans_z;
    std::tie(trans_x, trans_y, trans_z) = create_translations();

    append_primitives(std::move(trans_x));
    append_primitives(std::move(trans_y));
    append_primitives(std::move(trans_z));

    //TODO own function for rotations
    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> rotationA, rotationB, rotationC;
    std::tie(rotationA, rotationB, rotationC) = create_rotations();

    append_primitives(std::move(rotationA));
    append_primitives(std::move(rotationB));
    append_primitives(std::move(rotationC));
  }

  scon::mathmatrix<coords::float_type>& PrimitiveInternalCoordinates::guess_hessian() {
    using Mat = scon::mathmatrix<coords::float_type>;
    using scon::cross;
    using scon::dot;
    using scon::len;

    std::vector<coords::float_type> values;
    for (auto const & pic : primitive_internals) {
      values.emplace_back(pic->hessian_guess(xyz_));
    }

    hessian = Mat::col_from_vec(values).diagmat();

    return hessian;
  }

  scon::mathmatrix<coords::float_type>& PrimitiveInternalCoordinates::Bmat() {
    if (!new_B_matrix) {
      return B_matrix;
    }
    using Mat = scon::mathmatrix<coords::float_type>;

    auto ders = deriv_vec();

    std::size_t n_rows = ders.size(), n_cols = ders.at(0).size();
    B_matrix = Mat(n_rows, n_cols);

    for (std::size_t i{ 0 }; i < n_rows; ++i) {
      B_matrix.set_row(i, Mat::row_from_vec(ders.at(i)));
    }

    new_B_matrix = false;
    return B_matrix;
  }

  std::vector<std::vector<coords::float_type>> PrimitiveInternalCoordinates::deriv_vec() {
    std::vector<std::vector<coords::float_type>> result;

    for (auto const& pic : primitive_internals) {
      result.emplace_back(pic->der_vec(xyz_));
    }

    return result;
  }

  scon::mathmatrix<coords::float_type>& PrimitiveInternalCoordinates::Gmat() {
    if (!new_G_matrix) {
      return G_matrix;
    }
    PrimitiveInternalCoordinates::Bmat();
    G_matrix = B_matrix * B_matrix.t();
    new_G_matrix = false;
    return G_matrix;
  }

  scon::mathmatrix<coords::float_type> PrimitiveInternalCoordinates::calc(coords::Representation_3D const& xyz) const {
    std::vector<coords::float_type> primitives;
    primitives.reserve(primitive_internals.size());

    for (auto const & pic : primitive_internals) {
      primitives.emplace_back(pic->val(xyz));
    }

    return scon::mathmatrix<coords::float_type>::row_from_vec(primitives);
  }

  scon::mathmatrix<coords::float_type> PrimitiveInternalCoordinates::calc_diff(coords::Representation_3D const & lhs, coords::Representation_3D const & rhs) const {
    //TODO remove these from here
    for (auto & r : registeredRotators) {
      r->requestNewValueEvaluation();
    }
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
  }

  scon::mathmatrix<coords::float_type> PrimitiveInternalCoordinates::calculate_internal_grads(scon::mathmatrix<coords::float_type> const& g) {
    return Gmat().pinv() * Bmat() * g;
  }
}