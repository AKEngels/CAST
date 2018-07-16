#include <gtest/gtest.h>

#include <vector>
#include <tuple>

#include "../coords.h"
#include "../ic_core.h"
//#include "../scon_mathmatrix.h"

class InternalCoordinatesTest : public testing::Test {
public:
  InternalCoordinatesTest()
      : moleculeCartesianRepresentation{ coords::r3{ -6.053, -0.324, -0.108 },
                                         coords::r3{ -4.677, -0.093, -0.024 },
                                         coords::r3{ -6.262, -1.158, -0.813 },
                                         coords::r3{ -6.582, 0.600, -0.424 },
                                         coords::r3{ -6.431, -0.613, 0.894 },
                                         coords::r3{ -4.387, 0.166, -0.937 },
                                         coords::r3{ -6.146, 3.587, -0.024 },
                                         coords::r3{ -4.755, 3.671, -0.133 },
                                         coords::r3{ -6.427, 2.922, 0.821 },
                                         coords::r3{ -6.587, 3.223, -0.978 },
                                         coords::r3{ -6.552, 4.599, 0.179 },
                                         coords::r3{ -4.441, 2.753, -0.339 } },
        elementSymbols{ "C", "O", "H", "H", "H", "H",
                        "C", "O", "H", "H", "H", "H" },
        bond(1, 2, elementSymbols.at(0), elementSymbols.at(1)),
        angle(1, 2, 3, elementSymbols.at(0), elementSymbols.at(1),
              elementSymbols.at(2)),
        dihedralAngle(1, 2, 3, 4) {}

  double testBondLength() { return bond.val(moleculeCartesianRepresentation); }
  double testAngleValue() { return angle.val(moleculeCartesianRepresentation); }
  double testDihedralValue() {
    return dihedralAngle.val(moleculeCartesianRepresentation);
  }

  std::pair<coords::r3, coords::r3> testBondDerivatives() {

  }

private:
  coords::Representation_3D moleculeCartesianRepresentation;
  std::vector<std::string> elementSymbols;
  ic_core::distance bond;
  ic_core::angle angle;
  ic_core::dihedral dihedralAngle;
};

TEST_F(InternalCoordinatesTest, testBondLength) {
  auto bondLength = testBondLength();
  EXPECT_EQ(testBondLength(), 1.397781456451616);
}

TEST_F(InternalCoordinatesTest, testAngleValue) {
  auto angleVal = testAngleValue();
  EXPECT_EQ(testAngleValue(), 0.52901078997179596);
}

TEST_F(InternalCoordinatesTest, testDihedralValue) {
  auto angleVal = testDihedralValue();
  EXPECT_EQ(testDihedralValue(), 0.56342327253755931);
}