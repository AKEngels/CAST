#include "optimization_test.h"
#include "TestFiles\ExpectedValuesForTrustRadius.h"

namespace {
  double constexpr doubleNearThreshold = 1.e-10;
}

OptimizerTest::OptimizerTest() : gradients{ ExpectedValuesForTrustRadius::initialGradients() }, hessian{ ExpectedValuesForTrustRadius::initialHessianForTrust() }, cartesians { ExpectedValuesForTrustRadius::initialCartesians() }, testSystem{},
converter{ testSystem, cartesians }, restrictor(converter, ExpectedValuesForTrustRadius::initialTarget()), toCartesianNorm(converter, gradients, hessian, ExpectedValuesForTrustRadius::initialTrustRadius()),
brent{ 0.0, ExpectedValuesForTrustRadius::initialInternalNorm(), ExpectedValuesForTrustRadius::initialTrustRadius() }
{}

void OptimizerTest::restrictStepTest(){
  EXPECT_CALL(converter, getInternalStep(testing::_, testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::internalStepInitial()))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::expectedTrustStep()));
  EXPECT_CALL(converter, getDeltaYPrimeAndSol(testing::_, testing::_, testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::initialSolAndPrime()))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::finalSolAndPrime()));
  restrictor.setInitialV0(ExpectedValuesForTrustRadius::initialAlterationOfDiagonals());
  auto sol = restrictor(gradients, hessian);
  EXPECT_EQ(restrictor.getRestrictedStep(), ExpectedValuesForTrustRadius::expectedTrustStep());
  EXPECT_NEAR(sol, ExpectedValuesForTrustRadius::expectedSol(), doubleNearThreshold);
}

TEST_F(OptimizerTest, restrictStepTest) {
  restrictStepTest();
}

void OptimizerTest::restrictCartesianStepTest() {

  StepRestrictorMock restrictorMock{ converter, ExpectedValuesForTrustRadius::initialTarget() };
  EXPECT_CALL(restrictorMock, execute(testing::_, testing::_))
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::expectedSol()));
  EXPECT_CALL(restrictorMock, getRestrictedStep())
    .WillOnce(testing::ReturnRefOfCopy(ExpectedValuesForTrustRadius::expectedTrustStep()));
  EXPECT_CALL(restrictorMock, targetIsZero())
    .WillOnce(testing::Return(false));

  EXPECT_CALL(converter, getCartesianCoordinates())
    .WillOnce(testing::ReturnRefOfCopy(InternalCoordinates::CartesiansForInternalCoordinates(ExpectedValuesForTrustRadius::initialCartesians() / energy::bohr2ang)))
    .WillOnce(testing::ReturnRefOfCopy(InternalCoordinates::CartesiansForInternalCoordinates(ExpectedValuesForTrustRadius::cartesiansAfterwards() / energy::bohr2ang)));
  EXPECT_CALL(converter, applyInternalChange(testing::_));

  EXPECT_NEAR(toCartesianNorm(restrictorMock), ExpectedValuesForTrustRadius::expectedRestrictedCartesianNorm(), doubleNearThreshold);
}

void OptimizerTest::restrictCartesianStepWithZeroTargetTest() {

  StepRestrictorMock restrictorMock{ converter, 0.0 };
  EXPECT_CALL(restrictorMock, getTarget())
    .WillOnce(testing::Return(ExpectedValuesForTrustRadius::initialTarget()));
  EXPECT_CALL(restrictorMock, targetIsZero())
    .WillOnce(testing::Return(true));

  EXPECT_NEAR(toCartesianNorm(restrictorMock), ExpectedValuesForTrustRadius::initialTarget(), doubleNearThreshold);
}

TEST_F(OptimizerTest, restrictCartesianStepWithZeroTargetTest) {
  restrictCartesianStepWithZeroTargetTest();
}

void OptimizerTest::BrentsTrustStepTest() {
  InternalToCartesianStepMock internalCartesianStep{ converter, gradients, hessian, ExpectedValuesForTrustRadius::initialTrustRadius() };

  EXPECT_CALL(internalCartesianStep, execute(testing::_))
    .WillOnce(testing::Return(-0.282842712474619))
    .WillOnce(testing::Return(0.031306767697298))
    .WillOnce(testing::Return(0.000914776842788));

  EXPECT_NEAR(brent(internalCartesianStep).first, ExpectedValuesForTrustRadius::expectedBrent(), doubleNearThreshold);
}

TEST_F(OptimizerTest, BrentsTrustStepTest) {
  BrentsTrustStepTest();
}