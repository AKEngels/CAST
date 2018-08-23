#include "optimization_test.h"
#include "TestFiles\ExpectedValuesForTrustRadius.h"
#include "../Optimizer.h"

namespace {
  double constexpr doubleNearThreshold = 1.e-10;
}

AppropriateStepFinderMock::AppropriateStepFinderMock(internals::InternalToCartesianConverter const& converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian) :
  AppropriateStepFinder(converter, gradients, hessian) {}


FinderImplementation::FinderImplementation() : internals{}, cartesians{}, converter{ internals, cartesians }, emptyGradients{}, emptyHessian{ { 1., 0. },{ 0.,1. }, },
finder{ converter, emptyGradients, emptyHessian } {}

StepRestrictorTest::StepRestrictorTest() : FinderImplementation{}, expectedStep {}, restrictor{ &expectedStep, ExpectedValuesForTrustRadius::initialTarget() } {}

TEST_F(StepRestrictorTest, TestIfTargetIsZero) {
  scon::mathmatrix<coords::float_type> someDummy{};
  EXPECT_TRUE(internals::StepRestrictor(&someDummy, 0.0).targetIsZero());
}

TEST_F(StepRestrictorTest, restrictStepTest) {
  EXPECT_CALL(finder, getInternalStep())
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::internalStepInitial()));
  EXPECT_CALL(finder, getInternalStep(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::expectedTrustStep()));
  EXPECT_CALL(finder, getDeltaYPrime(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::initialDeltaYPrime()))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::finalDeltaYPrime()));
  EXPECT_CALL(finder, getSol(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::finalSol()));
  EXPECT_CALL(finder, alterHessian(testing::_)).Times(testing::AtLeast(1));

  restrictor.setInitialV0(ExpectedValuesForTrustRadius::initialAlterationOfDiagonals());
  auto sol = restrictor(finder);
  EXPECT_EQ(restrictor.getRestrictedStep(), ExpectedValuesForTrustRadius::expectedTrustStep());
  EXPECT_NEAR(sol, ExpectedValuesForTrustRadius::expectedSol(), doubleNearThreshold);
}

InternalToCartesianStepTest::InternalToCartesianStepTest() : FinderImplementation{}, fakeStep{}, toCartesianNorm { finder, ExpectedValuesForTrustRadius::initialTrustRadius() } {}

TEST_F(InternalToCartesianStepTest, convertFromInternalToCartesianTest_NonZeroTrial) {
  StepRestrictorMock NonZeroRestrictor{ &fakeStep, ExpectedValuesForTrustRadius::initialTarget() };
  
  EXPECT_CALL(NonZeroRestrictor, execute(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::expectedSol()));
  EXPECT_CALL(NonZeroRestrictor, getRestrictedStep())
    .WillOnce(testing::ReturnRefOfCopy(ExpectedValuesForTrustRadius::expectedTrustStep()));

  
  EXPECT_CALL(finder, applyInternalChangeAndGetNorm(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::restrictedStepFromInternalToCartesian()));

  EXPECT_NEAR(toCartesianNorm(NonZeroRestrictor), ExpectedValuesForTrustRadius::expectedRestrictedCartesianNorm(), doubleNearThreshold);
}

TEST_F(InternalToCartesianStepTest, convertFromInternalToCartesianTest_ZeroTrial) {
  internals::InternalToCartesianStep toCartesianNorm{ finder, ExpectedValuesForTrustRadius::initialTrustRadius() };

  StepRestrictorMock TargetZeroRestrictor{ &fakeStep, 0.0 };

  EXPECT_NEAR(toCartesianNorm(TargetZeroRestrictor), -ExpectedValuesForTrustRadius::initialTrustRadius(), doubleNearThreshold);
}

BrentsMethodTest::BrentsMethodTest() : FinderImplementation{}, brent{ finder, 0.0, ExpectedValuesForTrustRadius::initialInternalNorm(), ExpectedValuesForTrustRadius::initialTrustRadius() } {}

TEST_F(BrentsMethodTest, testFirstRestrictedStepOfTwoMethanol) {
  InternalToCartesianStepMock internalToCartesian{ finder, ExpectedValuesForTrustRadius::initialTrustRadius() };

  EXPECT_CALL(internalToCartesian, execute(testing::_))
    .WillOnce(testing::Return(-0.282842712474619))
    .WillOnce(testing::Return(0.031306767697298))
    .WillOnce(testing::Return(0.000914776842788));

  EXPECT_NEAR(brent(internalToCartesian), ExpectedValuesForTrustRadius::expectedBrent(), doubleNearThreshold);
}