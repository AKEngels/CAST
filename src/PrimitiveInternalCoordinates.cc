#include "PrimitiveInternalCoordinates.h"

#include "InternalCoordinateDecorator.h"
#include "Optimizer.h"

#include "Scon/scon_mathmatrix.h"
#include "InternalCoordinateUtilities.h"
#include "InternalCoordinates.h"

namespace internals {

	struct GradientsAndHessians {
		GradientsAndHessians(scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian) :
			gradients{ gradients }, hessian{ hessian }, inverseHessian{ hessian.pinv() }{}

		GradientsAndHessians(scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian, scon::mathmatrix<coords::float_type> && inverseHessian) :
			gradients{ gradients }, hessian{ hessian }, inverseHessian{ std::move(inverseHessian) }{}

		scon::mathmatrix<coords::float_type> gradients;
		scon::mathmatrix<coords::float_type> hessian;
		scon::mathmatrix<coords::float_type> inverseHessian;
	};

PrimitiveInternalCoordinates::PrimitiveInternalCoordinates() = default;

PrimitiveInternalCoordinates::PrimitiveInternalCoordinates(ICDecoratorBase &decorator) : B_matrix{ std::make_unique<scon::mathmatrix<coords::float_type>>() }, G_matrix{ std::make_unique<scon::mathmatrix<coords::float_type>>() }, hessian{std::make_unique<scon::mathmatrix<coords::float_type>>()} {
	decorator.appendCoordinates(*this);
}

PrimitiveInternalCoordinates::~PrimitiveInternalCoordinates() = default;

void PrimitiveInternalCoordinates::appendPrimitives(InternalVec&& pic) {
  primitive_internals.insert(primitive_internals.end(),
                             std::make_move_iterator(pic.begin()),
                             std::make_move_iterator(pic.end()));
}

void PrimitiveInternalCoordinates::appendRotators(
    std::vector<std::shared_ptr<InternalCoordinates::Rotator>> const&
        rotators) {
  for (auto const& curr_rot : rotators) {
    registeredRotators.emplace_back(curr_rot);
  }
}

std::unique_ptr<AppropriateStepFinder>
PrimitiveInternalCoordinates::constructStepFinder(
    InternalToCartesianConverter const& converter,
    scon::mathmatrix<coords::float_type> const& gradients,
    scon::mathmatrix<coords::float_type> const& hessianMatrix,
    CartesianType const& /*cartesians*/
) {
  return std::make_unique<AppropriateStepFinder>(converter, gradients, hessianMatrix);
}

// This function surely does not work.
/*inline std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>>
  PrimitiveInternalCoordinates::create_oops(const coords::Representation_3D&
coords, const BondGraph& g) const { using boost::adjacent_vertices; using
boost::vertices; using scon::dot;

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
          result.emplace_back(std::make_unique<InternalCoordinates::OutOfPlane>(g[*it],
g[combination.at(0)], g[combination.at(1)], g[combination.at(2)]));
        }
      }
      //}
    }
  }
  return result;
}//*/

scon::mathmatrix<coords::float_type>
PrimitiveInternalCoordinates::guess_hessian(
    InternalCoordinates::CartesiansForInternalCoordinates const& cartesians)
    const {
  using Mat = scon::mathmatrix<coords::float_type>;

  std::vector<coords::float_type> values;
  for (auto const& pic : primitive_internals) {
    values.emplace_back(cartesians.getInternalHessianGuess(*pic));
  }

  return Mat::col_from_vec(values).diagmat();
}

scon::mathmatrix<coords::float_type>&
PrimitiveInternalCoordinates::Bmat(CartesianType const& cartesians) {
  // TODO activate it again!!!!!
  /*if (!new_B_matrix) {
    return B_matrix;
  }*/
  using Mat = scon::mathmatrix<coords::float_type>;

  auto ders = deriv_vec(cartesians);

  std::size_t n_rows = ders.size(), n_cols = ders.at(0).size();
  *B_matrix = Mat(n_rows, n_cols);

  for (std::size_t i{ 0 }; i < n_rows; ++i) {
    B_matrix->set_row(i, Mat::row_from_vec(ders.at(i)));
  }

  new_B_matrix = false;
  // std::cout << "Bmat:\n" << std::fixed << std::setprecision(15) << B_matrix
  // << "\n\n";
  return *B_matrix;
}

scon::mathmatrix<coords::float_type>
PrimitiveInternalCoordinates::transposeOfBmat(CartesianType const& cartesian) {
  return Bmat(cartesian).t();
}

scon::mathmatrix<coords::float_type>
PrimitiveInternalCoordinates::pseudoInverseOfGmat(
    CartesianType const& cartesian) {
  return Gmat(cartesian).pinv();
}

scon::mathmatrix<coords::float_type>
PrimitiveInternalCoordinates::projectorMatrix(CartesianType const& cartesian) {
  return Gmat(cartesian) * pseudoInverseOfGmat(cartesian);
}

std::vector<std::vector<coords::float_type>>
PrimitiveInternalCoordinates::deriv_vec(CartesianType const& cartesians) {
  std::vector<std::vector<coords::float_type>> result;

  for (auto const& pic : primitive_internals) {
    result.emplace_back(cartesians.getInternalDerivativeVector(*pic));
  }

  return result;
}

scon::mathmatrix<coords::float_type>&
PrimitiveInternalCoordinates::Gmat(CartesianType const& cartesians) {
  // TODO activate it again!!!!!
  /*if (!new_G_matrix) {
    return G_matrix;
  }*/
  PrimitiveInternalCoordinates::Bmat(cartesians);
  *G_matrix = *B_matrix * B_matrix->t();
  new_G_matrix = false;
  return *G_matrix;
}

scon::mathmatrix<coords::float_type>
PrimitiveInternalCoordinates::calc(CartesianType const& xyz) const {
  std::vector<coords::float_type> primitives;
  primitives.reserve(primitive_internals.size());

  for (auto const& pic : primitive_internals) {
    primitives.emplace_back(xyz.getInternalValue(*pic));
  }

  return scon::mathmatrix<coords::float_type>::row_from_vec(primitives);
}

scon::mathmatrix<coords::float_type> PrimitiveInternalCoordinates::calc_diff(
	CartesianType const& lhs, CartesianType const& rhs) const {
  //// TODO remove these from here
  //for (auto& r : registeredRotators) {
  //  r->requestNewValueEvaluation();
  //}
  //auto lprims = PrimitiveInternalCoordinates::calc(lhs);
  //// TODO remove these from here
  //for (auto& r : registeredRotators) {
  //  r->requestNewValueEvaluation();
  //}
  //auto rprims = PrimitiveInternalCoordinates::calc(rhs);
  //auto diff = lprims - rprims;

	std::vector<coords::float_type> diff;
	diff.reserve(primitive_internals.size());

	for (auto const& pic : primitive_internals) {
		diff.emplace_back(lhs.getInternalDifference(rhs, *pic));
	}

	return scon::mathmatrix<coords::float_type>::row_from_vec(diff);
  // std::cout << "Diff:\n" << diff.t() << "\n";
  //return diff;
}

InternalToCartesianConverter::~InternalToCartesianConverter() = default;

scon::mathmatrix<coords::float_type>
InternalToCartesianConverter::calculateInternalGradients(
    scon::mathmatrix<coords::float_type> const& gradients) {
  return internalCoordinates.pseudoInverseOfGmat(cartesianCoordinates) *
         internalCoordinates.Bmat(cartesianCoordinates) * gradients;
}

std::pair<coords::float_type, coords::float_type>
InternalToCartesianConverter::cartesianNormOfOtherStructureAndCurrent(
    coords::Representation_3D const& otherCartesians) const {
  return cartesianCoordinates.displacementRmsValAndMaxTwoStructures(otherCartesians);
}

void InternalToCartesianConverter::takeCartesianStep(scon::mathmatrix<coords::float_type>&& d_cart){
	auto d_cart_rep3D = ic_util::matToRep3D(std::move(d_cart));
	set_xyz(cartesianCoordinates + d_cart_rep3D);
}

scon::mathmatrix<coords::float_type>
AppropriateStepFinder::alterHessian(coords::float_type const alteration) const {
  return matrices->hessian + scon::mathmatrix<coords::float_type>::identity(
	  matrices->hessian.rows(), matrices->hessian.cols()) *
                       alteration;
}

void StepRestrictor::randomizeAlteration(std::size_t const step) {
  static RandomNumberForHessianAlteration randomNumberForHessianAlteration;
  v0 += randomNumberForHessianAlteration.getRandomNumberBetweenZeroAndOne() *
        static_cast<coords::float_type>(step) / 100.;
}

//void StepRestrictor::swap(StepRestrictor & other) {
//	stepCallbackReference;
//	std::swap(stepCallbackReference, other.ste);
//	std::swap(cartesianCallbackReference, other.cartesianCallbackReference);
//	std::swap(target, other.target);
//	restrictedStep.swap(other.restrictedStep);
//	std::swap(correspondingCartesians, other.correspondingCartesians);
//	std::swap(restrictedSol, other.restrictedSol);
//	std::swap(v0, other.v0);
//}

StepRestrictor::StepRestrictor(StepRestrictor && other) { swap(std::move(other)); }

coords::float_type StepRestrictor::operator()(AppropriateStepFinder& finder) {
  *restrictedStep = finder.getInternalStep();
  auto deltaYPrime = finder.getDeltaYPrime(*restrictedStep);
  auto lastInternalNorm = 0.0;
  auto internalStepNorm = getStepNorm();
  auto i = 0u;
  for (; i < 1000; ++i) {
    v0 += (1. - internalStepNorm / target) * (internalStepNorm / deltaYPrime);
    *restrictedStep = finder.getInternalStep(finder.alterHessian(v0));
    internalStepNorm = getStepNorm();
    if (std::fabs(internalStepNorm - target) / target < 0.001) {
      return restrictedSol = finder.getSol(*restrictedStep);
    } else if (i > 10 && ((std::fabs(lastInternalNorm - internalStepNorm) /
                           internalStepNorm) < 0.001)) {
      return restrictedSol = finder.getSol(*restrictedStep);
    } else if ((i + 1u) % 100u == 0) {
      std::cout << "Trust Step did not converge after " << i + 1
                << " steps. Starting to randomize.\n";
      randomizeAlteration(i);
    }
    deltaYPrime = finder.getDeltaYPrime(*restrictedStep);
    lastInternalNorm = internalStepNorm;
  }
  std::cout << "Took over 1000 steps to retrict the trust step. Breaking up "
               "optimization and return current value.\n";
  return restrictedSol = finder.getSol(*restrictedStep);
  // throw std::runtime_error("Took over 1000 steps to retrict the trust step in
  // InternalToCartesianConverter::restrictStep. Breaking up optimization.");
}

coords::float_type InternalToCartesianStep::
operator()(StepRestrictor& restrictor) {
  if (restrictor.targetIsZero()) {
    return -trustRadius;
  }
  restrictor(finder);
  return finder.applyInternalChangeAndGetNorm(restrictor) - trustRadius;
}

bool BrentsMethod::useBisection() const {
  auto firstCondition =
      !((result > (3. * leftLimit + rightLimit) / 4.) && (result < rightLimit));
  auto secondCondition =
      (bisectionWasUsed &&
       std::fabs(result - rightLimit) >= std::fabs(rightLimit - middle) / 2.);
  auto thirdCondition =
      (!bisectionWasUsed &&
       std::fabs(result - rightLimit) >= std::fabs(result - oldMiddle) / 2.);
  auto fourthCondition =
      (bisectionWasUsed && std::fabs(rightLimit - middle) < delta);
  auto fifthCondition =
      (!bisectionWasUsed && std::fabs(middle - oldMiddle) < delta);
	if (Config::get().general.verbosity > 0)
	{
		std::cout << "Condition 1: " << std::boolalpha << firstCondition << "\n";
		std::cout << "Condition 2: " << std::boolalpha << secondCondition << "\n";
		std::cout << "Condition 3: " << std::boolalpha << thirdCondition << "\n";
		std::cout << "Condition 4: " << std::boolalpha << fourthCondition << "\n";
		std::cout << "Condition 5: " << std::boolalpha << fifthCondition << "\n";
	}
  return firstCondition || secondCondition || thirdCondition ||
         fourthCondition || fifthCondition;
}

coords::float_type BrentsMethod::
operator()(InternalToCartesianStep& internalToCartesianStep) {
  auto leftRestrictor = finder.generateStepRestrictor(leftLimit);
  auto rightRestrictor = finder.generateStepRestrictor(rightLimit);

  if (valueLeft * valueRight > 0.0) {
    throw std::runtime_error("In Brents method both limits are either positive "
                             "or negative. Thus no zero can be found.");
  }
  if (std::fabs(valueLeft) < std::fabs(valueRight)) {
    leftRestrictor.swap(std::move(rightRestrictor));
    std::swap(valueLeft, valueRight);
    std::swap(leftLimit, rightLimit);
  }
  middle = leftLimit;
  auto valueMiddle = valueLeft;

  auto epsilon = 0.01 > 1e-2 * std::fabs(valueLeft - valueRight)
                     ? 1e-2 * std::fabs(valueLeft - valueRight)
                     : 0.01;

  for (;;) {
		if (Config::get().general.verbosity > 0) {
			std::cout << std::fixed << std::setprecision(15) << valueLeft << "\n";
			std::cout << std::fixed << std::setprecision(15) << valueRight << "\n";
			std::cout << std::fixed << std::setprecision(15) << valueMiddle << "\n";
		}
    if (valueLeft != valueMiddle && valueRight != valueMiddle) {
      result = leftLimit * valueRight * valueMiddle /
               ((valueLeft - valueRight) * (valueLeft - valueMiddle));
      result += rightLimit * valueLeft * valueMiddle /
                ((valueMiddle - valueLeft) * (valueRight - valueRight));
      result += middle * valueLeft * valueRight /
                ((valueMiddle - valueLeft) * (valueMiddle - valueRight));
      std::cout << "Inverse Quadratic interpolation\n";
    } else {
      result = rightLimit -
               valueRight * (rightLimit - leftLimit) / (valueRight - valueLeft);
			if (Config::get().general.verbosity > 0) std::cout << "Secant method\n";
    }
		if (Config::get().general.verbosity > 0) std::cout << "Before Bisection" << std::fixed << std::setprecision(15)
              << result << "\n";
    if (useBisection()) {
      result = (leftLimit + rightLimit) / 2.;
      bisectionWasUsed = true;
      std::cout << "Bisection\n";
    } else {
      bisectionWasUsed = false;
    }
		if (Config::get().general.verbosity > 0) std::cout << "After Bisection" << std::fixed << std::setprecision(15)
              << result << "\n";
    auto resultRestrictor = finder.generateStepRestrictor(result);
    auto valueResult = internalToCartesianStep(resultRestrictor);

    if (valueResult / trustStep <= threshold) {
      resultRestrictor.registerBestGuess();
      return result;
    }

    oldMiddle = middle;
    middle = rightLimit;

    if (valueLeft * valueResult < 0.0) {
	  rightRestrictor.swap(std::move(resultRestrictor));
      std::swap(valueRight, valueResult);
      std::swap(rightLimit, result);
    } else {
      leftRestrictor.swap(std::move(resultRestrictor));
      std::swap(valueLeft, valueResult);
      std::swap(leftLimit, result);
    }

    if (std::fabs(leftLimit - rightLimit) < epsilon) {
      std::cout << "Warning: In Brents method the left and right interval are "
                   "getting too close. Now returning.\n";
      resultRestrictor.registerBestGuess();
      return result;
    }

    if (std::fabs(valueLeft) < std::fabs(valueRight)) {
      leftRestrictor.swap(std::move(rightRestrictor));
      std::swap(valueLeft, valueRight);
      std::swap(leftLimit, rightLimit);
    }
  }
  return 0.0;
}

StepRestrictorFactory::StepRestrictorFactory(AppropriateStepFinder& finder)
    : finalStep{ finder.getAddressOfInternalStep() }, finalCartesians{
        finder.getAddressOfCartesians()
      } {}


StepRestrictorFactory::~StepRestrictorFactory() = default;

AppropriateStepFinder::AppropriateStepFinder(InternalToCartesianConverter const& converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian) :
	matrices{ std::make_unique<GradientsAndHessians>(gradients, hessian) }, converter{ converter }, bestStepSoFar{ std::make_shared<scon::mathmatrix<coords::float_type>>() }, bestCartesiansSoFar{ std::make_shared<coords::Representation_3D>() }, stepRestrictorFactory{ *this }{}

AppropriateStepFinder::AppropriateStepFinder(InternalToCartesianConverter const& converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian, scon::mathmatrix<coords::float_type> && invertedHessian) :
	matrices{ std::make_unique<GradientsAndHessians>(gradients, hessian, std::move(invertedHessian)) }, converter{ converter }, bestStepSoFar{ std::make_shared<scon::mathmatrix<coords::float_type>>() }, bestCartesiansSoFar{ std::make_shared<coords::Representation_3D>() }, stepRestrictorFactory{ *this } {}


void AppropriateStepFinder::appropriateStep(
    coords::float_type const trustRadius) {
  auto cartesianNorm = applyInternalChangeAndGetNorm(getInternalStep());
  if (cartesianNorm > 1.1 * trustRadius) {
    if (Config::get().general.verbosity > 0) std::cout << "\033[31mTrust radius exceeded.\033[0m\n";
    InternalToCartesianStep internalToCartesianStep(*this, trustRadius);
    BrentsMethod brent(*this, 0.0, bestStepSoFar->norm(), trustRadius,
                       cartesianNorm); // very ugly
    brent(internalToCartesianStep);
  }
}

coords::float_type AppropriateStepFinder::applyInternalChangeAndGetNorm(
    scon::mathmatrix<coords::float_type> const& internalStep) {
  *bestCartesiansSoFar = converter.applyInternalChange(internalStep);
  *bestStepSoFar = internalStep; // TODO move the internal step
  return converter.cartesianNormOfOtherStructureAndCurrent(*bestCartesiansSoFar)
      .first;
}

coords::float_type AppropriateStepFinder::applyInternalChangeAndGetNorm(
    StepRestrictor& restrictor) {
  restrictor.setCartesians(
      converter.applyInternalChange(restrictor.getRestrictedStep()));
  return converter
      .cartesianNormOfOtherStructureAndCurrent(restrictor.getCartesians())
      .first;
}

coords::float_type AppropriateStepFinder::getDeltaYPrime(
    scon::mathmatrix<coords::float_type> const& internalStep) const {
  scon::mathmatrix<coords::float_type> deltaPrime =
      matrices->inverseHessian * internalStep * -1.0;
  return (internalStep.t() * deltaPrime)(0, 0) / internalStep.norm();
}

coords::float_type AppropriateStepFinder::getSol(
    scon::mathmatrix<coords::float_type> const& internalStep) const {
  auto transposedInternalStep = internalStep.t();
  return (transposedInternalStep * matrices->hessian * internalStep * 0.5)(0, 0) +
         (transposedInternalStep * matrices->gradients)(0, 0);
}

coords::float_type AppropriateStepFinder::getSolBestStep() const {
  return getSol(*bestStepSoFar);
}

scon::mathmatrix<coords::float_type>
AppropriateStepFinder::getInternalStep() const {
  return matrices->inverseHessian * matrices->gradients * -1.f;
}

scon::mathmatrix<coords::float_type> AppropriateStepFinder::getInternalStep(
    scon::mathmatrix<coords::float_type> const& hessian) const {
  return hessian.pinv() * matrices->gradients * -1.f;
}

InternalToCartesianStep::~InternalToCartesianStep() = default;

scon::mathmatrix<coords::float_type> InternalToCartesianConverter::calculateInternalValues()const {
	return internalCoordinates.calc(cartesianCoordinates);
}

void InternalToCartesianConverter::reset() {
  internalCoordinates.requestNewBAndG();
  cartesianCoordinates.reset();
}

void InternalToCartesianConverter::takeCartesianStep(
    scon::mathmatrix<coords::float_type>&& cartesianChange,
    InternalCoordinates::temporaryCartesian& cartesians) const {
  cartesians.coordinates += ic_util::matToRep3D(std::move(cartesianChange));
  cartesians.stolenNotify();
  internalCoordinates.requestNewBAndG();
}

coords::Representation_3D InternalToCartesianConverter::applyInternalChange(
    scon::mathmatrix<coords::float_type> d_int_left) const {
  using ic_util::flatten_c3_vec;
  // std::cout << "Internal Step:\n" << std::fixed << std::setprecision(15) <<
  // d_int_left << "\n\n";
  InternalCoordinates::temporaryCartesian actual_xyz = cartesianCoordinates;
  actual_xyz.stolenNotify();
  InternalCoordinates::temporaryCartesian old_xyz = cartesianCoordinates;

  coords::Representation_3D first_struct, last_good_xyz;
  auto micro_iter{ 0 }, fail_count{ 0 };
  auto damp{ 1. };
  auto old_inorm{ 0.0 };

  for (; micro_iter < 50; ++micro_iter) {
    /*std::cout << "Bmat MicroIteration " << micro_iter << "\n\n"
      << std::fixed << std::setprecision(15) <<
    internalCoordinates.transposeOfBmat(actual_xyz.coordinates) << "\n\n"
      << "Gmat inverser MicroIteration " << micro_iter << "\n\n"
      << std::fixed << std::setprecision(15) <<
    internalCoordinates.pseudoInverseOfGmat(actual_xyz.coordinates) << "\n\n";
    std::cout << "Cartesians in Microiteration " << micro_iter << ":\n\n"
            << actual_xyz.coordinates << "\n\n";*/

    takeCartesianStep(
        internalCoordinates.transposeOfBmat(actual_xyz.coordinates) *
            internalCoordinates.pseudoInverseOfGmat(actual_xyz.coordinates) *
            damp * d_int_left,
        actual_xyz);

    auto d_now = internalCoordinates
                     .calc_diff(actual_xyz.coordinates, old_xyz.coordinates)
                     .t();

    // std::cout << d_int_left.cols() << " " << d_int_left.rows() << " " <<
    // d_now.cols() << " " << d_now.rows() << std::endl;

    auto d_int_remain = d_int_left - d_now;
    auto cartesian_rmsd =
        ic_util::Rep3D_to_Mat<scon::mathmatrix>(old_xyz.coordinates - actual_xyz.coordinates)
            .rmsd(); // Optimizer::displacementRmsValAndMaxTwoStructures(old_xyz.coordinates,
                     // actual_xyz.coordinates).first;
    auto internal_norm = d_int_remain.norm();

    if (internal_norm != internal_norm) { /// Test for NaN
      throw std::runtime_error("Internal norm is NaN!");
    }
    // std::cout << "Left change internal coordinates:\n" << d_int_remain <<
    // "\n\n"; std::cout << "internal norm: " << internal_norm << "\n\n";
    if (micro_iter == 0) {
			if (Config::get().general.verbosity > 0)
			{
				std::cout << "Iter " << micro_iter + 1u
					<< " Internal Norm: " << std::scientific << std::setprecision(5)
					<< internal_norm << " Cartesian RMSD: " << std::scientific
					<< std::setprecision(5) << cartesian_rmsd
					<< " Damp: " << std::scientific << std::setprecision(5) << damp
					<< std::endl;
			}
      first_struct = actual_xyz.coordinates;
      last_good_xyz = actual_xyz.coordinates;
      old_inorm = internal_norm;
    } else {
      if (internal_norm > old_inorm) {
				if (Config::get().general.verbosity > 0)
				{
					std::cout << "Iter " << micro_iter + 1u
						<< " Internal Norm: " << std::scientific
						<< std::setprecision(5) << internal_norm
						<< " Last Internal Norm: " << std::scientific
						<< std::setprecision(5) << old_inorm
						<< " Cartesian RMSD: " << std::scientific
						<< std::setprecision(5) << cartesian_rmsd
						<< " Damp: " << std::scientific << std::setprecision(5)
						<< damp << std::endl;
				}
        damp /= 2.;
        ++fail_count;
      } else {
				if (Config::get().general.verbosity > 0)
				{
					std::cout << "Iter " << micro_iter + 1u
						<< " Internal Norm: " << std::scientific
						<< std::setprecision(5) << internal_norm
						<< " Last Internal Norm: " << std::scientific
						<< std::setprecision(5) << old_inorm
						<< " Cartesian RMSD: " << std::scientific
						<< std::setprecision(5) << cartesian_rmsd
						<< " Damp: " << std::scientific << std::setprecision(5)
						<< damp << std::endl;
				}
        fail_count = 0;
        damp = std::min(1.2 * damp, 1.);
        old_inorm = internal_norm;
        last_good_xyz = actual_xyz.coordinates;
      }
    }
    if (cartesian_rmsd < 1.e-6 || internal_norm < 1.e-6) {
			if (Config::get().general.verbosity > 0) {
				std::cout << "Applying internal changes took " << micro_iter + 1u << " steps to converge." << std::endl;
			}
      return actual_xyz.coordinates;
    } else if (fail_count >= 10) {
			if (Config::get().general.verbosity > 0)
			{
				std::cout << "Applying internal changes failed ten times to converge. "
					"Aborting after " << micro_iter + 1u << "steps." << std::endl;
			}
      return first_struct;
    }

    old_xyz = actual_xyz;
    d_int_left = std::move(d_int_remain);
  }
	if (Config::get().general.verbosity > 0)
	{
		std::cout << "Applying internal changes took all " << micro_iter + 1u
			<< " steps, still not converged. Returning best step."
			<< std::endl;
	}
  return actual_xyz.coordinates;
}

StepRestrictor::StepRestrictor(std::shared_ptr<scon::mathmatrix<coords::float_type>> step, std::shared_ptr<coords::Representation_3D> cartesians, coords::float_type const target) :
	stepCallbackReference{ step }, cartesianCallbackReference{ cartesians }, target{ target }, restrictedStep{std::make_unique<scon::mathmatrix<coords::float_type>>()}, correspondingCartesians{std::make_unique<coords::Representation_3D>()}, restrictedSol{ 0.0 }, v0{ 0.0 } {}

StepRestrictor::~StepRestrictor() = default;

coords::float_type StepRestrictor::getStepNorm() const { return restrictedStep->norm(); }

StepRestrictor
StepRestrictorFactory::makeStepRestrictor(coords::float_type const target) {
  return StepRestrictor(finalStep, finalCartesians, target );
}

void StepRestrictor::registerBestGuess() {
  *stepCallbackReference = std::move(*restrictedStep);
  *cartesianCallbackReference = std::move(*correspondingCartesians);
}

StepRestrictor
AppropriateStepFinder::generateStepRestrictor(coords::float_type const target) {
  return stepRestrictorFactory.makeStepRestrictor(target);
}

AppropriateStepFinder::~AppropriateStepFinder() = default;
BrentsMethod::~BrentsMethod() = default;

} // namespace internals
