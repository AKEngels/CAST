#include <gtest/gtest.h>

#include <vector>
#include <tuple>

#include "../coords.h"
#include "../ic_core.h"
//#include "../scon_mathmatrix.h"

namespace {
  double constexpr doubleNearThreshold = 1.e-10;
  coords::r3 constexpr r3NearThreshold{ doubleNearThreshold, doubleNearThreshold, doubleNearThreshold };
  inline void isCartesianPointNear(coords::r3 const& lhs, coords::r3 const& rhs) {
    EXPECT_NEAR(lhs.x(), rhs.x(), doubleNearThreshold);
    EXPECT_NEAR(lhs.y(), rhs.y(), doubleNearThreshold);
    EXPECT_NEAR(lhs.z(), rhs.z(), doubleNearThreshold);
  }
  
}

struct TwoMethanolMolecules {
  TwoMethanolMolecules() : moleculeCartesianRepresentation{ 
    coords::r3{ -6.053, -0.324, -0.108 },
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
    "C", "O", "H", "H", "H", "H" } {
    moleculeCartesianRepresentation /= energy::bohr2ang;
  }

  std::vector<std::string> elementSymbols;
  coords::Representation_3D moleculeCartesianRepresentation;
};

class InternalCoordinatesTest : public testing::Test {
public:
  
  InternalCoordinatesTest()
      : twoMethanolMolecules{std::make_unique<TwoMethanolMolecules>()}, 
        bond(1, 2, twoMethanolMolecules->elementSymbols.at(0), twoMethanolMolecules->elementSymbols.at(1)),
        angle(1, 2, 3, twoMethanolMolecules->elementSymbols.at(0), twoMethanolMolecules->elementSymbols.at(1),
          twoMethanolMolecules->elementSymbols.at(2)),
        dihedralAngle(1, 2, 3, 4) {}

  double testBondLength() { return bond.val(twoMethanolMolecules->moleculeCartesianRepresentation); }
  double testAngleValue() { return angle.val(twoMethanolMolecules->moleculeCartesianRepresentation); }
  double testDihedralValue() {
    return dihedralAngle.val(twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  std::pair<coords::r3, coords::r3> testBondDerivatives() {
    return bond.der(twoMethanolMolecules->moleculeCartesianRepresentation);
  }
  std::tuple<coords::r3, coords::r3, coords::r3> testAngleDerivatives() {
    return angle.der(twoMethanolMolecules->moleculeCartesianRepresentation);
  }
  std::tuple<coords::r3, coords::r3, coords::r3, coords::r3> testDihedralDerivatives() {
    return dihedralAngle.der(twoMethanolMolecules->moleculeCartesianRepresentation);
  }

private:
  std::unique_ptr<TwoMethanolMolecules> twoMethanolMolecules;
  ic_core::distance bond;
  ic_core::angle angle;
  ic_core::dihedral dihedralAngle;
};

//class InternalCoordinateDistanceTest : 

TEST_F(InternalCoordinatesTest, testBondLength) {
  auto bla = testBondLength();
  EXPECT_NEAR(testBondLength(), 2.6414241359371124, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTest, testAngleValue) {
  auto bla = testAngleValue();
  EXPECT_NEAR(testAngleValue(), 0.52901078997179563, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTest, testDihedralValue) {
  auto bla = testDihedralValue();
  EXPECT_NEAR(testDihedralValue(), 0.56342327253755953, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTest, testBondDerivatives) {
  std::pair<coords::r3, coords::r3> testValuesForDerivatives;

  auto bondDerivatives = testBondDerivatives();

  isCartesianPointNear(bondDerivatives.first, coords::r3{ -0.98441712304088669, -0.16526188620817209, -0.060095231348426204 });
  isCartesianPointNear(bondDerivatives.second, coords::r3{ 0.98441712304088669, 0.16526188620817209, 0.060095231348426204 });
}

TEST_F(InternalCoordinatesTest, testAngleDerivatives) {
  std::tuple<coords::r3, coords::r3, coords::r3> testValuesForDerivatives;

  auto angleDerivatives = testAngleDerivatives();

  isCartesianPointNear(std::get<0>(angleDerivatives), coords::r3{ -0.062056791850036874, 0.27963965489457271, 0.24754030125005314 });
  isCartesianPointNear(std::get<1>(angleDerivatives), coords::r3{ -0.10143003133725584, -0.13768014867512546, -0.11073454633478358 });
  isCartesianPointNear(std::get<2>(angleDerivatives), coords::r3{ 0.16348682318729271, -0.14195950621944725, -0.13680575491526956 });
}

TEST_F(InternalCoordinatesTest, testDihedralDerivatives) {
  std::tuple<coords::r3, coords::r3, coords::r3> testValuesForDerivatives;

  auto dihedralDerivatives = testDihedralDerivatives();

  isCartesianPointNear(std::get<0>(dihedralDerivatives), coords::r3{ 0.047760930904702938, -0.49023624122103959, 0.56578012853796311 });
  isCartesianPointNear(std::get<1>(dihedralDerivatives), coords::r3{ -0.056149623980491933, 0.17150456631353736, -0.11870115223680294 });
  isCartesianPointNear(std::get<2>(dihedralDerivatives), coords::r3{ -0.084250009050912372, 0.23597149360212683, -0.14926917153430774 });
  isCartesianPointNear(std::get<3>(dihedralDerivatives), coords::r3{ 0.092638702126701375, 0.082760181305375408, -0.29780980476685243 });
}