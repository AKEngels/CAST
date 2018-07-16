#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "../coords.h"
#include "../ic_core.h"
//#include "../scon_mathmatrix.h"

namespace {
double constexpr doubleNearThreshold = 1.e-10;
coords::r3 constexpr r3NearThreshold{ doubleNearThreshold, doubleNearThreshold,
                                      doubleNearThreshold };
inline void isCartesianPointNear(coords::r3 const& lhs, coords::r3 const& rhs) {
  EXPECT_NEAR(lhs.x(), rhs.x(), doubleNearThreshold);
  EXPECT_NEAR(lhs.y(), rhs.y(), doubleNearThreshold);
  EXPECT_NEAR(lhs.z(), rhs.z(), doubleNearThreshold);
}

} // namespace

struct TwoMethanolMolecules {
  TwoMethanolMolecules()
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
                        "C", "O", "H", "H", "H", "H" } {
    moleculeCartesianRepresentation /= energy::bohr2ang;
  }

  std::vector<std::string> elementSymbols;
  coords::Representation_3D moleculeCartesianRepresentation;
};

class InternalCoordinatesTest : public testing::Test {
public:
  InternalCoordinatesTest()
      : twoMethanolMolecules{ std::make_unique<TwoMethanolMolecules>() } {}

protected:
  std::unique_ptr<TwoMethanolMolecules> twoMethanolMolecules;
};

class InternalCoordinatesDistancesTest : public InternalCoordinatesTest {
public:
  InternalCoordinatesDistancesTest()
      : InternalCoordinatesTest(),
        firstAtomDerivatives{ -0.98441712304088669, -0.16526188620817209,
                              -0.060095231348426204 },
        secondAtomDerivatives{ 0.98441712304088669, 0.16526188620817209,
                               0.060095231348426204 },
        bond(1, 2, twoMethanolMolecules->elementSymbols.at(0),
             twoMethanolMolecules->elementSymbols.at(1)),
        derivativeVector(3 * 12, 0.) {
    derivativeVector.at(0) = firstAtomDerivatives.x();
    derivativeVector.at(1) = firstAtomDerivatives.y();
    derivativeVector.at(2) = firstAtomDerivatives.z();
    derivativeVector.at(3) = secondAtomDerivatives.x();
    derivativeVector.at(4) = secondAtomDerivatives.y();
    derivativeVector.at(5) = secondAtomDerivatives.z();
  }

  double testBondLength() {
    return bond.val(twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  std::pair<coords::r3, coords::r3> testBondDerivatives() {
    return bond.der(twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  void derivativeVectorTest() {
    auto derivatives = bond.der_vec(twoMethanolMolecules->moleculeCartesianRepresentation);
    for (auto i = 1u; i < derivatives.size(); ++i) {
      EXPECT_NEAR(derivatives.at(i), derivativeVector.at(i), doubleNearThreshold);
    }
  }

  double hessianGuessTest() {
    return bond.hessian_guess(
        twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  std::string returnInfoTest() {
    return bond.info(twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  coords::r3 const firstAtomDerivatives;
  coords::r3 const secondAtomDerivatives;

private:
  ic_core::distance bond;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesAnglesTest : public InternalCoordinatesTest {
public:
  InternalCoordinatesAnglesTest()
      : InternalCoordinatesTest(), leftAtomsDerivative{ -0.062056791850036874,
                                                        0.27963965489457271,
                                                        0.24754030125005314 },
        middleAtomsDerivative{ -0.10143003133725584, -0.13768014867512546,
                               -0.11073454633478358 },
        rightAtomsDerivative{ 0.16348682318729271, -0.14195950621944725,
                              -0.13680575491526956 },
        angle(1, 2, 3, twoMethanolMolecules->elementSymbols.at(0),
              twoMethanolMolecules->elementSymbols.at(1),
              twoMethanolMolecules->elementSymbols.at(2)),
        derivativeVector(3 * 12, 0.) {
    derivativeVector.at(0) = leftAtomsDerivative.x();
    derivativeVector.at(1) = leftAtomsDerivative.y();
    derivativeVector.at(2) = leftAtomsDerivative.z();
    derivativeVector.at(3) = middleAtomsDerivative.x();
    derivativeVector.at(4) = middleAtomsDerivative.y();
    derivativeVector.at(5) = middleAtomsDerivative.z();
    derivativeVector.at(6) = rightAtomsDerivative.x();
    derivativeVector.at(7) = rightAtomsDerivative.y();
    derivativeVector.at(8) = rightAtomsDerivative.z();
  }

  double testAngleValue() {
    return angle.val(twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  std::tuple<coords::r3, coords::r3, coords::r3> testAngleDerivatives() {
    return angle.der(twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  void derivativeVectorTest() {
    auto derivatives = angle.der_vec(twoMethanolMolecules->moleculeCartesianRepresentation);
    for (auto i = 1u; i < derivatives.size(); ++i) {
      EXPECT_NEAR(derivatives.at(i), derivativeVector.at(i), doubleNearThreshold);
    }
  }

  double hessianGuessTest() {
    return angle.hessian_guess(
      twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  std::string returnInfoTest() {
    return angle.info(twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  coords::r3 const leftAtomsDerivative;
  coords::r3 const middleAtomsDerivative;
  coords::r3 const rightAtomsDerivative;

private:
  ic_core::angle angle;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesDihedralsTest : public InternalCoordinatesTest {
public:
  InternalCoordinatesDihedralsTest()
      : InternalCoordinatesTest(), leftLeftDerivative{ 0.047760930904702938,
                                                       -0.49023624122103959,
                                                       0.56578012853796311 },
        leftMiddleDerivative{ -0.056149623980491933, 0.17150456631353736,
                              -0.11870115223680294 },
        rightMiddleDerivative{ -0.084250009050912372, 0.23597149360212683,
                               -0.14926917153430774 },
        rightRightDerivative{ 0.092638702126701375, 0.082760181305375408,
                              -0.29780980476685243 },
        dihedralAngle(1, 2, 3, 4), derivativeVector(3 * 12, 0.) {
    derivativeVector.at(0) = leftLeftDerivative.x();
    derivativeVector.at(1) = leftLeftDerivative.y();
    derivativeVector.at(2) = leftLeftDerivative.z();
    derivativeVector.at(3) = leftMiddleDerivative.x();
    derivativeVector.at(4) = leftMiddleDerivative.y();
    derivativeVector.at(5) = leftMiddleDerivative.z();
    derivativeVector.at(6) = rightMiddleDerivative.x();
    derivativeVector.at(7) = rightMiddleDerivative.y();
    derivativeVector.at(8) = rightMiddleDerivative.z();
    derivativeVector.at(9) = rightRightDerivative.x();
    derivativeVector.at(10) = rightRightDerivative.y();
    derivativeVector.at(11) = rightRightDerivative.z();
  }

  coords::r3 leftLeftDerivative;
  coords::r3 leftMiddleDerivative;
  coords::r3 rightMiddleDerivative;
  coords::r3 rightRightDerivative;

  double testDihedralValue() {
    return dihedralAngle.val(
        twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
  testDihedralDerivatives() {
    return dihedralAngle.der(
        twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  void derivativeVectorTest() {
    auto derivatives = dihedralAngle.der_vec(twoMethanolMolecules->moleculeCartesianRepresentation);
    for (auto i = 1u; i < derivatives.size(); ++i) {
      EXPECT_NEAR(derivatives.at(i), derivativeVector.at(i), doubleNearThreshold);
    }
  }

  double hessianGuessTest() {
    return dihedralAngle.hessian_guess(
      twoMethanolMolecules->moleculeCartesianRepresentation);
  }

  std::string returnInfoTest() {
    return dihedralAngle.info(twoMethanolMolecules->moleculeCartesianRepresentation);
  }

private:
  ic_core::dihedral dihedralAngle;
  std::vector<double> derivativeVector;
};

TEST_F(InternalCoordinatesDistancesTest, testBondLength) {
  auto bla = testBondLength();
  EXPECT_NEAR(testBondLength(), 2.6414241359371124, doubleNearThreshold);
}

TEST_F(InternalCoordinatesDistancesTest, testBondDerivatives) {
  std::pair<coords::r3, coords::r3> testValuesForDerivatives;

  auto bondDerivatives = testBondDerivatives();

  isCartesianPointNear(bondDerivatives.first, firstAtomDerivatives);
  isCartesianPointNear(bondDerivatives.second, secondAtomDerivatives);
}

TEST_F(InternalCoordinatesDistancesTest, derivativeVectorTest) {
  derivativeVectorTest();
}

TEST_F(InternalCoordinatesDistancesTest, hessianGuessTest) {
  EXPECT_NEAR(hessianGuessTest(), 0.45990191969419725, doubleNearThreshold);
}

TEST_F(InternalCoordinatesDistancesTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Bond: 2.64142 || 1 || 2 ||");
}

TEST_F(InternalCoordinatesAnglesTest, testAngleValue) {
  auto bla = testAngleValue();
  EXPECT_NEAR(testAngleValue(), 0.52901078997179563, doubleNearThreshold);
}

TEST_F(InternalCoordinatesAnglesTest, testAngleDerivatives) {
  std::tuple<coords::r3, coords::r3, coords::r3> testValuesForDerivatives;

  auto angleDerivatives = testAngleDerivatives();

  isCartesianPointNear(std::get<0>(angleDerivatives), leftAtomsDerivative);
  isCartesianPointNear(std::get<1>(angleDerivatives), middleAtomsDerivative);
  isCartesianPointNear(std::get<2>(angleDerivatives), rightAtomsDerivative);
}

TEST_F(InternalCoordinatesAnglesTest, derivativeVectorTest) {
  derivativeVectorTest();
}

TEST_F(InternalCoordinatesAnglesTest, hessianGuessTest) {
  EXPECT_NEAR(hessianGuessTest(), 0.16, doubleNearThreshold);
}

TEST_F(InternalCoordinatesAnglesTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Angle: 30.3101 || 1 || 2 || 3 ||");
}

TEST_F(InternalCoordinatesDihedralsTest, testDihedralValue) {
  auto bla = testDihedralValue();
  EXPECT_NEAR(testDihedralValue(), 0.56342327253755953, doubleNearThreshold);
}

TEST_F(InternalCoordinatesDihedralsTest, testDihedralDerivatives) {
  std::tuple<coords::r3, coords::r3, coords::r3> testValuesForDerivatives;

  auto dihedralDerivatives = testDihedralDerivatives();

  isCartesianPointNear(std::get<0>(dihedralDerivatives), leftLeftDerivative);
  isCartesianPointNear(std::get<1>(dihedralDerivatives), leftMiddleDerivative);
  isCartesianPointNear(std::get<2>(dihedralDerivatives), rightMiddleDerivative);
  isCartesianPointNear(std::get<3>(dihedralDerivatives), rightRightDerivative);
}

TEST_F(InternalCoordinatesDihedralsTest, derivativeVectorTest) {
  derivativeVectorTest();
}

TEST_F(InternalCoordinatesDihedralsTest, hessianGuessTest) {
  EXPECT_NEAR(hessianGuessTest(), 0.023, doubleNearThreshold);
}

TEST_F(InternalCoordinatesDihedralsTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Dihedral: 32.2818 || 1 || 2 || 3 || 4 ||");
}