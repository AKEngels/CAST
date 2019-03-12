#ifdef GOOGLE_MOCK

#include "optimization_test.h"
#include "TestFiles/ExpectedValuesForTrustRadius.h"
#include "../Optimizer.h"

namespace {
  double constexpr doubleNearThreshold = 1.e-10;

  inline void isCartesianPointNear(coords::r3 const& lhs, coords::r3 const& rhs) {
    EXPECT_NEAR(lhs.x(), rhs.x(), doubleNearThreshold);
    EXPECT_NEAR(lhs.y(), rhs.y(), doubleNearThreshold);
    EXPECT_NEAR(lhs.z(), rhs.z(), doubleNearThreshold);
  }
}

AppropriateStepFinderMock::AppropriateStepFinderMock(internals::InternalToCartesianConverter const& converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian) :
  AppropriateStepFinder(converter, gradients, hessian) {}

AppropriateStepFinderAvoidingInverseOfHessian::AppropriateStepFinderAvoidingInverseOfHessian(internals::InternalToCartesianConverter const& converter, 
  scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian) 
  : AppropriateStepFinder(converter, gradients, hessian, ExpectedValuesForTrustRadius::initialInverseHessianForTrust()) {}

FinderImplementation::FinderImplementation() : internals{}, cartesians{}, converter{ internals, cartesians }, emptyGradients{}, emptyHessian{ { 1., 0. },{ 0.,1. }, },
finder{ converter, emptyGradients, emptyHessian } {}

StepRestrictorTest::StepRestrictorTest() : FinderImplementation{}, expectedStep{}, expectedCartesians{}, restrictor{ &expectedStep, &expectedCartesians, ExpectedValuesForTrustRadius::initialTarget() } {}

TEST_F(StepRestrictorTest, TestIfTargetIsZero) {
  EXPECT_TRUE(internals::StepRestrictor(&expectedStep, &expectedCartesians, 0.0).targetIsZero());
}

TEST_F(StepRestrictorTest, restrictStepTest) {
  EXPECT_CALL(finder, getInternalStep())
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::internalStepInitial()));
  EXPECT_CALL(finder, getInternalStep(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::expectedTrustStep()));
  EXPECT_CALL(finder, getDeltaYPrime(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::initialDeltaYPrime()))/*
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::finalDeltaYPrime()))*/;
  EXPECT_CALL(finder, getSol(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::finalSol()));
  EXPECT_CALL(finder, alterHessian(testing::_)).Times(testing::AtLeast(1));

  restrictor.setInitialV0(ExpectedValuesForTrustRadius::initialAlterationOfDiagonals());
  auto sol = restrictor(finder);
  EXPECT_EQ(restrictor.getRestrictedStep(), ExpectedValuesForTrustRadius::expectedTrustStep());
  EXPECT_NEAR(sol, ExpectedValuesForTrustRadius::expectedSol(), doubleNearThreshold);
}

TEST_F(StepRestrictorTest, registerBestGuessTest) {
  auto bla = restrictor.getRestrictedStep();
}

InternalToCartesianStepTest::InternalToCartesianStepTest() : FinderImplementation{}, fakeStep{}, fakeCartesians{}, toCartesianNorm{ finder, ExpectedValuesForTrustRadius::initialTrustRadius() } {}

TEST_F(InternalToCartesianStepTest, convertFromInternalToCartesianTest_NonZeroTrial) {
  StepRestrictorMock NonZeroRestrictor{ &fakeStep, &fakeCartesians, ExpectedValuesForTrustRadius::initialTarget() };
  
  EXPECT_CALL(NonZeroRestrictor, execute(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::expectedSol()));
  
  EXPECT_CALL(finder, applyInternalChangeAndGetNorm(testing::Matcher<internals::StepRestrictor &>(testing::_)))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::restrictedStepFromInternalToCartesian()));

  EXPECT_NEAR(toCartesianNorm(NonZeroRestrictor), ExpectedValuesForTrustRadius::expectedRestrictedCartesianNorm(), doubleNearThreshold);
}

TEST_F(InternalToCartesianStepTest, convertFromInternalToCartesianTest_ZeroTrial) {
  internals::InternalToCartesianStep toCartesianNorm{ finder, ExpectedValuesForTrustRadius::initialTrustRadius() };

  StepRestrictorMock TargetZeroRestrictor{ &fakeStep, &fakeCartesians, 0.0 };

  EXPECT_NEAR(toCartesianNorm(TargetZeroRestrictor), -ExpectedValuesForTrustRadius::initialTrustRadius(), doubleNearThreshold);
}

BrentsMethodTest::BrentsMethodTest() : FinderImplementation{}, brent{ finder, 0.0, ExpectedValuesForTrustRadius::initialInternalNorm(), ExpectedValuesForTrustRadius::initialTrustRadius(), ExpectedValuesForTrustRadius::initialCartesianNorm() } {}

TEST_F(BrentsMethodTest, testFirstRestrictedStepOfTwoMethanol) {
  InternalToCartesianStepMock internalToCartesian{ finder, ExpectedValuesForTrustRadius::initialTrustRadius() };

  EXPECT_CALL(internalToCartesian, execute(testing::_))
    .WillOnce(testing::Return(0.000914776842788));

  EXPECT_NEAR(brent(internalToCartesian), ExpectedValuesForTrustRadius::expectedBrent(), doubleNearThreshold);
}

AppropriateStepFinderTest::AppropriateStepFinderTest() : internals{}, cartesians{ ExpectedValuesForTrustRadius::initialCartesians() / energy::bohr2ang },
converter{ internals, cartesians }, gradients{ ExpectedValuesForTrustRadius::initialGradients() }, hessian{ ExpectedValuesForTrustRadius::initialHessianForTrust() }, 
finder{ converter, gradients, hessian } {}

TEST_F(AppropriateStepFinderTest, getDeltaYPrimeTest) {
  EXPECT_NEAR(finder.getDeltaYPrime(ExpectedValuesForTrustRadius::internalStepInitial()), ExpectedValuesForTrustRadius::initialDeltaYPrime(), doubleNearThreshold);
}

TEST_F(AppropriateStepFinderTest, getSol) {
  EXPECT_NEAR(finder.getSol(ExpectedValuesForTrustRadius::internalStepInitial()), ExpectedValuesForTrustRadius::initialSol(), doubleNearThreshold);
}

TEST_F(AppropriateStepFinderTest, getInternalStepWithMemberHessian) {
  EXPECT_EQ(finder.getInternalStep(), ExpectedValuesForTrustRadius::internalStepInitial());
}

TEST_F(AppropriateStepFinderTest, getInternalStepWithPassedHessian) {
  EXPECT_EQ(finder.getInternalStep(hessian), ExpectedValuesForTrustRadius::internalStepInitial());
}

TEST_F(AppropriateStepFinderTest, applyChangeAndGetNormWithRestrictorTest) {
  scon::mathmatrix<coords::float_type> step;
  coords::Representation_3D cartesians;
  StepRestrictorMock restrictor{ &step, &cartesians, ExpectedValuesForTrustRadius::initialTarget() };

  EXPECT_CALL(restrictor, getRestrictedStep())
    .WillOnce(testing::ReturnRefOfCopy(ExpectedValuesForTrustRadius::expectedTrustStep()));

  EXPECT_CALL(converter, applyInternalChange(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::cartesiansAfterwards() / energy::bohr2ang));

  EXPECT_NEAR(finder.applyInternalChangeAndGetNorm(restrictor), ExpectedValuesForTrustRadius::restrictedStepFromInternalToCartesian(), doubleNearThreshold);
}

TEST_F(AppropriateStepFinderTest, applyChangeAndGetNormWithStepVectorTest) {
  scon::mathmatrix<coords::float_type> step;
  coords::Representation_3D cartesians;

  EXPECT_CALL(converter, applyInternalChange(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::cartesiansAfterwards() / energy::bohr2ang));

  EXPECT_NEAR(finder.applyInternalChangeAndGetNorm(ExpectedValuesForTrustRadius::internalStepInitial()), ExpectedValuesForTrustRadius::restrictedStepFromInternalToCartesian(), doubleNearThreshold);
}

TEST_F(AppropriateStepFinderTest, alterHessianTest) {
  double constexpr alteration = 0.25;
  EXPECT_EQ(finder.alterHessian(alteration), hessian + scon::mathmatrix<coords::float_type>::identity(hessian.cols(), hessian.rows()) * alteration);
}

TEST_F(AppropriateStepFinderTest, appropriateStepTest) {
  EXPECT_CALL(converter, applyInternalChange(testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::cartesiansTooFar() / energy::bohr2ang))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::cartesiansAfterwards() / energy::bohr2ang));

  finder.appropriateStep(ExpectedValuesForTrustRadius::initialTrustRadius());

  auto cartesians = finder.extractCartesians();
  auto finalCartesians = ExpectedValuesForTrustRadius::cartesiansAfterwards() / energy::bohr2ang;
  for (auto i = 0u; i < cartesians.size(); ++i) {
    isCartesianPointNear(cartesians.at(i), finalCartesians.at(i));
  }

  EXPECT_EQ(finder.extractBestStep(), ExpectedValuesForTrustRadius::expectedTrustStep());
}
#endif