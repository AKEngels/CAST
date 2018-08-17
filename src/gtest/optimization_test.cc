#include "optimization_test.h"
#include "TestFiles\ExpectedValuesForTrustRadius.h"
#include "primitive_internals_test.h"

namespace {
  double constexpr doubleNearThreshold = 1.e-10;
}

OptimizerTest::OptimizerTest() : converter{ ExpectedValuesForTrustRadius::initialCartesians(), MockPrimitiveInternals{} } {}

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
