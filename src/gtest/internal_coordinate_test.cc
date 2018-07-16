#include <gtest/gtest.h>

#include "../ic_core.h"
#include "../coords.h"

class InternalCoordinatesTest : public testing::Test {
public:
  InternalCoordinatesTest()
      : moleculeCartesianRepresentation{ coords::r3{ 0.000, 0.000, 0.000 },
                                         coords::r3{ -0.684, 1.184, -0.483 },
                                         coords::r3{ -0.024, 2.327, -0.017 },
                                         coords::r3{ -0.000, -0.000, 1.089 },
                                         coords::r3{ 1.027, -0.000, -0.363 },
                                         coords::r3{ -0.513, -0.889, -0.363 },
                                         coords::r3{ -0.684, 1.184, -1.572 },
                                         coords::r3{ -1.710, 1.184, -0.120 },
                                         coords::r3{ -0.470, 3.100, -0.332 } },
        bond(1, 2, "C", "C") {}

  double testBondLength(){
    return bond.val(moleculeCartesianRepresentation);
  }

private:
  coords::Representation_3D moleculeCartesianRepresentation;
  ic_core::distance bond;
};

TEST_F(InternalCoordinatesTest, testBondLength) {
  auto bondLength = testBondLength();
  EXPECT_EQ(testBondLength(), 1.4501727483303497);
}