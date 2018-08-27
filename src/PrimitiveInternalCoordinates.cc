#include "PrimitiveInternalCoordinates.h"
#include "Optimizer.h"

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
      reference.emplace_back(target.at(ind-1));
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

  std::tuple<PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec, PrimitiveInternalCoordinates::InternalVec> PrimitiveInternalCoordinates::create_rotations(CartesianType & cartesians) {
    std::vector<InternalCoordinates::Rotations> result;
    for (auto const& indices : subSystemIndices) {
      result.emplace_back(build_rotation(cartesians, indices)->makeRotations());
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
      auto core_cp = coords.at(core - 1u);
      auto vert_i = adjacent_vertices(*it, g);
      auto sym = g[*it].element;
      auto dist = (vert_i.second - vert_i.first);
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

  void PrimitiveInternalCoordinates::create_ic_system(BondGraph const& g, CartesianType & cartesians) {
    append_primitives(create_distances(g));
    append_primitives(create_angles(g));
    append_primitives(create_oops(cartesians, g));
    append_primitives(create_dihedrals(g));

    //TODO own function for translations
    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> trans_x, trans_y, trans_z;
    std::tie(trans_x, trans_y, trans_z) = create_translations();

    append_primitives(std::move(trans_x));
    append_primitives(std::move(trans_y));
    append_primitives(std::move(trans_z));

    //TODO own function for rotations
    std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> rotationA, rotationB, rotationC;
    std::tie(rotationA, rotationB, rotationC) = create_rotations(cartesians);

    append_primitives(std::move(rotationA));
    append_primitives(std::move(rotationB));
    append_primitives(std::move(rotationC));
  }

  scon::mathmatrix<coords::float_type> PrimitiveInternalCoordinates::guess_hessian(InternalCoordinates::CartesiansForInternalCoordinates const& cartesians) const {
    using Mat = scon::mathmatrix<coords::float_type>;

    std::vector<coords::float_type> values;
    for (auto const & pic : primitive_internals) {
      values.emplace_back(pic->hessian_guess(cartesians));
    }

    return Mat::col_from_vec(values).diagmat();
  }

  scon::mathmatrix<coords::float_type>& PrimitiveInternalCoordinates::Bmat(CartesianType const& cartesians) {
    if (!new_B_matrix) {
      return B_matrix;
    }
    using Mat = scon::mathmatrix<coords::float_type>;

    auto ders = deriv_vec(cartesians);

    std::size_t n_rows = ders.size(), n_cols = ders.at(0).size();
    B_matrix = Mat(n_rows, n_cols);

    for (std::size_t i{ 0 }; i < n_rows; ++i) {
      B_matrix.set_row(i, Mat::row_from_vec(ders.at(i)));
    }

    new_B_matrix = false;
    return B_matrix;
  }

  scon::mathmatrix<coords::float_type> PrimitiveInternalCoordinates::transposeOfBmat(CartesianType const& cartesian) {
    return Bmat(cartesian).t();
  }

  scon::mathmatrix<coords::float_type> PrimitiveInternalCoordinates::pseudoInverseOfGmat(CartesianType const& cartesian) {
    return Gmat(cartesian).pinv();
  }

  std::vector<std::vector<coords::float_type>> PrimitiveInternalCoordinates::deriv_vec(CartesianType const& cartesians) {
    std::vector<std::vector<coords::float_type>> result;

    for (auto const& pic : primitive_internals) {
      result.emplace_back(pic->der_vec(cartesians));
    }

    return result;
  }

  scon::mathmatrix<coords::float_type>& PrimitiveInternalCoordinates::Gmat(CartesianType const& cartesians) {
    if (!new_G_matrix) {
      return G_matrix;
    }
    PrimitiveInternalCoordinates::Bmat(cartesians);
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

  scon::mathmatrix<coords::float_type> InternalToCartesianConverter::calculateInternalGradients(scon::mathmatrix<coords::float_type> const& gradients) {
    return internalCoordinates.pseudoInverseOfGmat(cartesianCoordinates) * internalCoordinates.Bmat(cartesianCoordinates) * gradients;
  }

  std::pair<coords::float_type, coords::float_type> InternalToCartesianConverter::cartesianNormOfOtherStructureAndCurrent(coords::Representation_3D const& otherCartesians) const {
    return Optimizer::displacementRmsValAndMaxTwoStructures(cartesianCoordinates, otherCartesians);
  }

  scon::mathmatrix<coords::float_type> AppropriateStepFinder::alterHessian(coords::float_type const alteration) const {
    return hessian + alteration * scon::mathmatrix<coords::float_type>::identity(hessian.rows(), hessian.cols());
  }

  coords::float_type StepRestrictor::randomizeAlteration(std::size_t const step){
    static RandomNumberForHessianAlteration randomNumberForHessianAlteration;
    v0 += randomNumberForHessianAlteration.getRandomNumberBetweenZeroAndOne() * static_cast<coords::float_type>(step) / 100.;
  }

  coords::float_type StepRestrictor::operator()(AppropriateStepFinder & finder){
    restrictedStep = finder.getInternalStep();
    auto deltaYPrime = finder.getDeltaYPrime(restrictedStep);
    auto lastInternalNorm = 0.0;
    auto internalStepNorm = getStepNorm();
    auto i = 0u;
    for (; i < 1000; ++i) {
      v0 += (1. - internalStepNorm / target)*(internalStepNorm / deltaYPrime);
      restrictedStep = finder.getInternalStep(finder.alterHessian(v0));
      internalStepNorm = getStepNorm();
      if (std::fabs(internalStepNorm - target) / target < 0.001) {
        return restrictedSol = finder.getSol(restrictedStep);
      }
      else if(i>10 && ((std::fabs(lastInternalNorm - internalStepNorm) / internalStepNorm)<0.001)){
        return restrictedSol = finder.getSol(restrictedStep);
      }
      else if((i+1u) % 100u == 0){
	std::cout << "Trust Step did not converge after " << i+1 << " steps. Starting to randomize.\n";
	randomizeAlteration(i);
      }
      deltaYPrime = finder.getDeltaYPrime(restrictedStep);
      lastInternalNorm = internalStepNorm;
    }
    std::cout << "Took over 1000 steps to retrict the trust step. Breaking up optimization and return current value.\n";
    return restrictedSol = finder.getSol(restrictedStep);
    //throw std::runtime_error("Took over 1000 steps to retrict the trust step in InternalToCartesianConverter::restrictStep. Breaking up optimization.");
  }

  coords::float_type InternalToCartesianStep::operator()(StepRestrictor & restrictor){
    if (restrictor.targetIsZero()) {
      return -trustRadius;
    }
    auto expect = restrictor(finder);
    return finder.applyInternalChangeAndGetNorm(restrictor) - trustRadius;
  }

  bool BrentsMethod::useBisection() const {
    auto firstCondition = !((result > (3.*leftLimit + rightLimit) / 4.) && (result < rightLimit));
    auto secondCondition = (bisectionWasUsed && std::fabs(result - rightLimit) >= std::fabs(rightLimit - middle) / 2.);
    auto thirdCondition = (!bisectionWasUsed && std::fabs(result - rightLimit) >= std::fabs(result - oldMiddle)/2.);
    auto fourthCondition = (bisectionWasUsed && std::fabs(rightLimit - middle) < delta);
    auto fifthCondition = (!bisectionWasUsed && std::fabs(middle-oldMiddle) < delta);
    return firstCondition || secondCondition || thirdCondition || fourthCondition || fifthCondition;
  }

  coords::float_type BrentsMethod::operator()(InternalToCartesianStep & internalToCartesianStep){
    auto leftRestrictor = finder.generateStepRestrictor(leftLimit);
    auto rightRestrictor = finder.generateStepRestrictor(rightLimit);

    if (valueLeft * valueRight > 0.0) {
      throw std::runtime_error("In Brents method both limits are either positive or negative. Thus no zero can be found.");
    }
    if (std::fabs(valueLeft) < std::fabs(valueRight)) {
      std::swap(leftRestrictor, rightRestrictor);
      std::swap(valueLeft, valueRight);
      std::swap(leftLimit, rightLimit);
    }
    middle = leftLimit;
    auto valueMiddle = valueLeft;

    auto epsilon = 0.01 > 1e-2*std::fabs(valueLeft - valueRight) ? 1e-2*std::fabs(valueLeft - valueRight) : 0.01;

    for(;;) {
      if (valueLeft != valueMiddle && valueRight != valueMiddle) {
        result = leftLimit * valueRight * valueMiddle / ((valueLeft - valueRight) * (valueLeft - valueMiddle));
        result += rightLimit * valueLeft * valueMiddle / ((valueMiddle - valueLeft) * (valueRight - valueRight));
        result += middle * valueLeft * valueRight / ((valueMiddle - valueLeft) * (valueMiddle - valueRight));
      }
      else {
        result = rightLimit - valueRight * (rightLimit - leftLimit) / (valueRight - valueLeft);
      }
      if (useBisection()) {
        result = (leftLimit + rightLimit) / 2.;
        bisectionWasUsed = true;
      }
      else {
        bisectionWasUsed = false;
      }
      auto resultRestrictor = finder.generateStepRestrictor(result);
      auto valueResult = internalToCartesianStep(resultRestrictor);

      if (valueResult / trustStep <= threshold) {
        resultRestrictor.registerBestGuess();
        return result;
      }

      oldMiddle = middle;
      middle = rightLimit;

      if (valueLeft * valueResult < 0.0) {
        std::swap(rightRestrictor, resultRestrictor);
        std::swap(valueRight, valueResult);
        std::swap(rightLimit, result);
      }
      else {
        std::swap(leftRestrictor, resultRestrictor);
        std::swap(valueLeft, valueResult);
        std::swap(leftLimit, result);
      }

      if (std::fabs(leftLimit - rightLimit) < epsilon) {
        std::cout << "Warning: In Brents method the left and right interval are getting too close. Now returning.\n";
        resultRestrictor.registerBestGuess();
        return result;
      }

      if (std::fabs(valueLeft) < std::fabs(valueRight)) {
        std::swap(leftRestrictor, rightRestrictor);
        std::swap(valueLeft, valueRight);
        std::swap(leftLimit, rightLimit);
      }
    }
    return 0.0;
  }

  StepRestrictorFactory::StepRestrictorFactory(AppropriateStepFinder & finder) : finalStep{ finder.getAddressOfInternalStep() }, finalCartesians{ finder.getAddressOfCartesians() } {}

  void AppropriateStepFinder::appropriateStep(coords::float_type const trustRadius){
    auto cartesianNorm = applyInternalChangeAndGetNorm(getInternalStep());
    if (cartesianNorm > 1.1*trustRadius) {
      InternalToCartesianStep internalToCartesianStep(*this, trustRadius);
      BrentsMethod brent(*this, 0.0, bestStepSoFar.norm(), trustRadius, cartesianNorm);//very ugly
      brent(internalToCartesianStep);
    }
  }

  coords::float_type AppropriateStepFinder::applyInternalChangeAndGetNorm(scon::mathmatrix<coords::float_type> const& internalStep) {
    bestCartesiansSoFar = converter.applyInternalChange(internalStep);
    bestStepSoFar = internalStep;//TODO move the internal step
    return converter.cartesianNormOfOtherStructureAndCurrent(bestCartesiansSoFar).first;
  }
  
  coords::float_type AppropriateStepFinder::applyInternalChangeAndGetNorm(StepRestrictor & restrictor) {
    restrictor.setCartesians(converter.applyInternalChange(restrictor.getRestrictedStep()));
    return converter.cartesianNormOfOtherStructureAndCurrent(restrictor.getCartesians()).first;
  }

  coords::float_type AppropriateStepFinder::getDeltaYPrime(scon::mathmatrix<coords::float_type> const & internalStep) const {
    scon::mathmatrix<coords::float_type> deltaPrime = -1.*inverseHessian*internalStep;
    return (internalStep.t()*deltaPrime)(0, 0) / internalStep.norm();
  }

  coords::float_type AppropriateStepFinder::getSol(scon::mathmatrix<coords::float_type> const & internalStep) const {
    auto transposedInternalStep = internalStep.t();
    return (0.5 * transposedInternalStep * hessian * internalStep)(0, 0) + (transposedInternalStep*gradients)(0, 0);
  }

  coords::float_type AppropriateStepFinder::getSolBestStep() const{
    return getSol(bestStepSoFar);
  }

  scon::mathmatrix<coords::float_type> AppropriateStepFinder::getInternalStep() const {
    return -1.*inverseHessian*gradients;
  }

  scon::mathmatrix<coords::float_type> AppropriateStepFinder::getInternalStep(scon::mathmatrix<coords::float_type> const& hessian) const {
    return -1.*hessian.pinv()*gradients;
  }

  void InternalToCartesianConverter::takeCartesianStep(scon::mathmatrix <coords::float_type> && cartesianChange, InternalCoordinates::temporaryCartesian & cartesians) const {
    cartesians.coordinates += ic_util::mat_to_rep3D(std::move(cartesianChange));
    cartesians.stolenNotify();
    internalCoordinates.requestNewBAndG();
  }

  coords::Representation_3D InternalToCartesianConverter::applyInternalChange(scon::mathmatrix<coords::float_type> d_int_left) const {
    using ic_util::flatten_c3_vec;

    InternalCoordinates::temporaryCartesian actual_xyz = cartesianCoordinates;
    InternalCoordinates::temporaryCartesian old_xyz = cartesianCoordinates;;
    coords::Representation_3D first_struct, last_good_xyz;
    auto micro_iter{ 0 }, fail_count{ 0 };
    auto damp{ 1. };
    auto old_inorm{ 0.0 };
    for (; micro_iter < 50; ++micro_iter) {
      takeCartesianStep(damp*internalCoordinates.transposeOfBmat(actual_xyz.coordinates)*internalCoordinates.pseudoInverseOfGmat(actual_xyz.coordinates)*d_int_left, actual_xyz);

      auto d_now = internalCoordinates.calc_diff(actual_xyz.coordinates, old_xyz.coordinates);

      auto d_int_remain = d_int_left - d_now;
      auto cartesian_rmsd = ic_util::Rep3D_to_Mat(old_xyz.coordinates - actual_xyz.coordinates).rmsd();
      auto internal_norm = d_int_remain.norm();
      //std::cout << "Left change internal coordinates:\n" << d_int_remain << "\n\n";
      //std::cout << "internal norm: " << internal_norm << "\n\n";
      if (micro_iter == 0) {
        first_struct = actual_xyz.coordinates;
        last_good_xyz = actual_xyz.coordinates;
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
          last_good_xyz = actual_xyz.coordinates;
        }
      }
      if (cartesian_rmsd < 1.e-6 || internal_norm < 1.e-6) {
        std::cout << "Took " << micro_iter << " steps to converge.\n";
        return actual_xyz.coordinates;
      }
      else if (fail_count >= 10) {
        std::cout << "Failed ten times to converge.\n";
        return first_struct;
      }

      old_xyz = actual_xyz;
      d_int_left = std::move(d_int_remain);
    }
    std::cout << "Took all " << micro_iter + 1 << " steps, still not converged.\n";
    return actual_xyz.coordinates;
  }

  StepRestrictor StepRestrictorFactory::makeStepRestrictor(coords::float_type const target) {
    return StepRestrictor{ finalStep, finalCartesians, target };
  }

  void StepRestrictor::registerBestGuess() {
    *stepCallbackReference = std::move(restrictedStep);
    *cartesianCallbackReference = std::move(correspondingCartesians);
  }

  StepRestrictor AppropriateStepFinder::generateStepRestrictor(coords::float_type const target) {
    return stepRestrictorFactory.makeStepRestrictor(target);
  }
}
