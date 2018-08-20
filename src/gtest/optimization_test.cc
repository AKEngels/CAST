#include "optimization_test.h"
#include "TestFiles\ExpectedValuesForTrustRadius.h"

namespace {
  double constexpr doubleNearThreshold = 1.e-10;
}

OptimizerTest::OptimizerTest() : gradients{ ExpectedValuesForTrustRadius::initialGradients() }, hessian{ ExpectedValuesForTrustRadius::initialHessianForTrust() }, cartesians { ExpectedValuesForTrustRadius::initialCartesians() }, testSystem{},
 converter{ testSystem, cartesians }, conv{ testSystem, cartesians }, restrictor(converter, ExpectedValuesForTrustRadius::initialTarget()), toCartesianNorm(conv, gradients, hessian, ExpectedValuesForTrustRadius::initialTrustRadius())
{}

void OptimizerTest::restrictStepTest(){
  restrictor.setInitialV0(ExpectedValuesForTrustRadius::initialAlterationOfDiagonals());
  auto sol = restrictor(gradients, hessian);
  EXPECT_EQ(restrictor.getRestrictedStep(), ExpectedValuesForTrustRadius::expectedTrustStep());
  EXPECT_NEAR(sol, ExpectedValuesForTrustRadius::expectedSol(), doubleNearThreshold);
}

TEST_F(OptimizerTest, restrictStepTest) {
  restrictStepTest();
}

void OptimizerTest::restrictCartesianStepTest() {
  EXPECT_CALL(conv, getInternalStep(testing::_, testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::internalStepInitial()))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::expectedTrustStep()));
  EXPECT_CALL(conv, getDeltaYPrimeAndSol(testing::_, testing::_, testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::initialSolAndPrime()))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::finalSolAndPrime()));
  EXPECT_CALL(conv, getCartesianCoordinates())
    .WillOnce(testing::ReturnRefOfCopy(InternalCoordinates::CartesiansForInternalCoordinates(ExpectedValuesForTrustRadius::initialCartesians() / energy::bohr2ang)))
    .WillOnce(testing::ReturnRefOfCopy(InternalCoordinates::CartesiansForInternalCoordinates(ExpectedValuesForTrustRadius::cartesiansAfterwards() / energy::bohr2ang)));
  EXPECT_CALL(conv, applyInternalChange(testing::_));

  EXPECT_NEAR(toCartesianNorm(ExpectedValuesForTrustRadius::initialTarget()), ExpectedValuesForTrustRadius::expectedRestrictedCartesianNorm(), doubleNearThreshold);
}

TEST_F(OptimizerTest, restrictCartesianStepTest) {
  restrictCartesianStepTest();
}