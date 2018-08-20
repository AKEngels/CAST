#include "optimization_test.h"
#include "TestFiles\ExpectedValuesForTrustRadius.h"

namespace {
  double constexpr doubleNearThreshold = 1.e-10;
}

OptimizerTest::OptimizerTest() : cartesians{ ExpectedValuesForTrustRadius::initialCartesians() }, testSystem{}, converter{ testSystem, cartesians } {}

void OptimizerTest::restrictStepTest(){
  auto expectedValues = converter.restrictStep(
    ExpectedValuesForTrustRadius::initialTarget(),
    ExpectedValuesForTrustRadius::initialAlterationOfDiagonals(),
    ExpectedValuesForTrustRadius::initialCartesians(),
    ExpectedValuesForTrustRadius::initialGradients(),
    ExpectedValuesForTrustRadius::initialHessianForTrust()
  );
  EXPECT_EQ(expectedValues.first, ExpectedValuesForTrustRadius::expectedTrustStep());
  EXPECT_NEAR(expectedValues.second, ExpectedValuesForTrustRadius::expectedSol(), doubleNearThreshold);
}

TEST_F(OptimizerTest, restrictStepTest) {
  restrictStepTest();
}
