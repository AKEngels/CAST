#include "optimization_test.h"
#include "TestFiles\ExpectedValuesForTrustRadius.h"

namespace {
  double constexpr doubleNearThreshold = 1.e-10;
}

OptimizerTest::OptimizerTest() : cartesians{ ExpectedValuesForTrustRadius::initialCartesians() }, testSystem{}, converter{ testSystem, cartesians }, 
  restrictor(converter, ExpectedValuesForTrustRadius::initialTarget()), toCartesianNorm(converter, ExpectedValuesForTrustRadius::initialTarget()) {}

void OptimizerTest::restrictStepTest(){
  restrictor.setInitialV0(ExpectedValuesForTrustRadius::initialAlterationOfDiagonals());
  auto sol = restrictor(
    ExpectedValuesForTrustRadius::initialGradients(),
    ExpectedValuesForTrustRadius::initialHessianForTrust()
  );
  EXPECT_EQ(restrictor.getRestrictedStep(), ExpectedValuesForTrustRadius::expectedTrustStep());
  EXPECT_NEAR(sol, ExpectedValuesForTrustRadius::expectedSol(), doubleNearThreshold);
}

TEST_F(OptimizerTest, restrictStepTest) {
  restrictStepTest();
}

void OptimizerTest::restrictCartesianStepTest() {
  EXPECT_NEAR(toCartesianNorm(ExpectedValuesForTrustRadius::initialTarget()), ExpectedValuesForTrustRadius::expectedRestrictedCartesianNorm(), doubleNearThreshold);
}

TEST_F(OptimizerTest, restrictCartesianStepTest) {
  restrictCartesianStepTest();
}