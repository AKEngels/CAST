#ifdef GOOGLE_MOCK

#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "internal_coordinate_test.h"
#include "../ic_util.h"

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

SubsystemOfTwoMethanolMolecules::SubsystemOfTwoMethanolMolecules()
  : subSystem{ {coords::r3{ -6.053, -0.324, -0.108 },
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
                                     coords::r3{ -4.441, 2.753, -0.339 } }, {
      "C", "O", "H", "H", "H", "H", "C", "O", "H", "H", "H", "H"
    } } {
  subSystem.cartesianRepresentation /= energy::bohr2ang;
}

RotatetdMethanolMolecules::RotatetdMethanolMolecules()
    : initialMethanolSystem{ { coords::r3{ -0.321, -0.087, 0.12733 },
                               coords::r3{ 1.055, 0.144, 0.21133 },
                               coords::r3{ -0.53, -0.921, -0.57767 },
                               coords::r3{ -0.85, 0.837, -0.18867 },
                               coords::r3{ -0.699, -0.376, 1.12933 },
                               coords::r3{ 1.345, 0.403, -0.70167 } },
                             { "C", "O", "H", "H", "H", "H" } },
      rotatedMethanolSystem{ {
                                 coords::r3{ 0.12733, -0.087, 0.321 },
                                 coords::r3{ 0.21133, 0.144, -1.055 },
                                 coords::r3{ -0.57767, -0.921, 0.53 },
                                 coords::r3{ -0.18867, 0.837, 0.85 },
                                 coords::r3{ 1.12933, -0.376, 0.699 },
                                 coords::r3{ -0.70167, 0.403, -1.345 },
                             },
                             { "C", "O", "H", "H", "H", "H" } } {
  initialMethanolSystem.cartesianRepresentation /= energy::bohr2ang;
  rotatedMethanolSystem.cartesianRepresentation /= energy::bohr2ang;
}

InternalCoordinatesDistancesTest::InternalCoordinatesDistancesTest()
    : InternalCoordinatesTest(), firstAtomDerivatives{ -0.98441712304088669,
                                                       -0.16526188620817209,
                                                       -0.060095231348426204 },
      secondAtomDerivatives{ 0.98441712304088669, 0.16526188620817209,
                             0.060095231348426204 },
      bond(1, 2, twoMethanolMolecules->getOneRepresentation().elementSymbols.at(0),
           twoMethanolMolecules->getOneRepresentation().elementSymbols.at(1)),
      derivativeVector(3u * 12u, 0.) {
  derivativeVector.at(0) = firstAtomDerivatives.x();
  derivativeVector.at(1) = firstAtomDerivatives.y();
  derivativeVector.at(2) = firstAtomDerivatives.z();
  derivativeVector.at(3) = secondAtomDerivatives.x();
  derivativeVector.at(4) = secondAtomDerivatives.y();
  derivativeVector.at(5) = secondAtomDerivatives.z();
}

double InternalCoordinatesDistancesTest::testBondLength() {
  return bond.val(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

std::pair<coords::r3, coords::r3> InternalCoordinatesDistancesTest::testBondDerivatives() {
  return bond.der(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

void InternalCoordinatesDistancesTest::derivativeVectorTest() {
  auto derivatives =
      bond.der_vec(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
  for (auto i = 1u; i < derivatives.size(); ++i) {
    EXPECT_NEAR(derivatives.at(i), derivativeVector.at(i), doubleNearThreshold);
  }
}

double InternalCoordinatesDistancesTest::hessianGuessTest() {
  return bond.hessian_guess(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

std::string InternalCoordinatesDistancesTest::returnInfoTest() {
  return bond.info(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

InternalCoordinatesAnglesTest::InternalCoordinatesAnglesTest()
    : InternalCoordinatesTest(), leftAtomsDerivative{ -0.062056791850036874,
                                                      0.27963965489457271,
                                                      0.24754030125005314 },
      middleAtomsDerivative{ -0.10143003133725584, -0.13768014867512546,
                             -0.11073454633478358 },
      rightAtomsDerivative{ 0.16348682318729271, -0.14195950621944725,
                            -0.13680575491526956 },
      angle(1, 2, 3, twoMethanolMolecules->getOneRepresentation().elementSymbols.at(0),
            twoMethanolMolecules->getOneRepresentation().elementSymbols.at(1),
            twoMethanolMolecules->getOneRepresentation().elementSymbols.at(2)),
      derivativeVector(3u * 12u, 0.) {
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

double InternalCoordinatesAnglesTest::testAngleValue() {
  return angle.val(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

std::tuple<coords::r3, coords::r3, coords::r3> InternalCoordinatesAnglesTest::testAngleDerivatives() {
  return angle.der(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

void InternalCoordinatesAnglesTest::derivativeVectorTest() {
  auto derivatives =
      angle.der_vec(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
  for (auto i = 1u; i < derivatives.size(); ++i) {
    EXPECT_NEAR(derivatives.at(i), derivativeVector.at(i), doubleNearThreshold);
  }
}

double InternalCoordinatesAnglesTest::hessianGuessTest() {
  return angle.hessian_guess(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

std::string InternalCoordinatesAnglesTest::returnInfoTest() {
  return angle.info(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

InternalCoordinatesDihedralsTest::InternalCoordinatesDihedralsTest()
    : InternalCoordinatesTest(), leftLeftDerivative{ 0.047760930904702938,
                                                     -0.49023624122103959,
                                                     0.56578012853796311 },
      leftMiddleDerivative{ -0.056149623980491933, 0.17150456631353736,
                            -0.11870115223680294 },
      rightMiddleDerivative{ -0.084250009050912372, 0.23597149360212683,
                             -0.14926917153430774 },
      rightRightDerivative{ 0.092638702126701375, 0.082760181305375408,
                            -0.29780980476685243 },
      dihedralAngle(1, 2, 3, 4), derivativeVector(3u * 12u, 0.) {
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

double InternalCoordinatesDihedralsTest::testDihedralValue() {
  return dihedralAngle.val(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
InternalCoordinatesDihedralsTest::testDihedralDerivatives() {
  return dihedralAngle.der(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

void InternalCoordinatesDihedralsTest::derivativeVectorTest() {
  auto derivatives = dihedralAngle.der_vec(
      twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
  for (auto i = 1u; i < derivatives.size(); ++i) {
    EXPECT_NEAR(derivatives.at(i), derivativeVector.at(i), doubleNearThreshold);
  }
}

double InternalCoordinatesDihedralsTest::hessianGuessTest() {
  return dihedralAngle.hessian_guess(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

std::string InternalCoordinatesDihedralsTest::returnInfoTest() {
  return dihedralAngle.info(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

InternalCoordinatesTranslationTest::InternalCoordinatesTranslationTest() : InternalCoordinatesTest(), translation{ { 1u,2u,3u,4u,5u,6u } } {}

void InternalCoordinatesTranslationTest::testTranslationDerivativeTest() {

  auto flattenedVector = ic_util::flatten_c3_vec(translation.der(6, [](auto const & s) {
    auto doubleS = static_cast<double>(s);
    return coords::r3(1. / doubleS, 2. / doubleS, 3. / doubleS);
  }));

  std::vector<double> someTestValues;

  for (auto i = 0; i < 6; ++i) {
    someTestValues.emplace_back(1. / static_cast<double>(translation.indices_.size()));
    someTestValues.emplace_back(2. / static_cast<double>(translation.indices_.size()));
    someTestValues.emplace_back(3. / static_cast<double>(translation.indices_.size()));
  }

  for (auto i = 1u; i < flattenedVector.size(); ++i) {
    EXPECT_NEAR(flattenedVector.at(i), someTestValues.at(i), doubleNearThreshold);
  }

}

double InternalCoordinatesTranslationTest::hessianGuessTest() {
  return translation.hessian_guess(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

InternalCoordinatesTranslationXTest::InternalCoordinatesTranslationXTest() : InternalCoordinatesTest(), translation{ {1u,2u,3u,4u,5u,6u} },
derivativeVector(3u*12u,0.) {
  auto constexpr sizeOfMethanolWhichIsDescribed = 3u * 12u / 2u;
  for (auto i = 0u; i < sizeOfMethanolWhichIsDescribed; i += 3) {
    derivativeVector.at(i) = 0.16666666666666666;
  }
}

double InternalCoordinatesTranslationXTest::testTranslationValue() {
  return translation.val(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

void InternalCoordinatesTranslationXTest::derivativeVectorTest() {
  auto derivatives = translation.der_vec(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
  for (auto i = 0u; i < derivatives.size(); ++i) {
    EXPECT_NEAR(derivatives.at(i), derivativeVector.at(i), doubleNearThreshold);
  }
}

std::string InternalCoordinatesTranslationXTest::returnInfoTest() {
  return translation.info(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

InternalCoordinatesTranslationYTest::InternalCoordinatesTranslationYTest() : InternalCoordinatesTest(), translation{ { 1u,2u,3u,4u,5u,6u } },
derivativeVector(3u * 12u, 0.) {
  auto constexpr sizeOfMethanolWhichIsDescribed = 3u * 12u / 2u;
  for (auto i = 1u; i < sizeOfMethanolWhichIsDescribed; i += 3) {
    derivativeVector.at(i) = 0.16666666666666666;
  }
}

double InternalCoordinatesTranslationYTest::testTranslationValue() {
  return translation.val(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

void InternalCoordinatesTranslationYTest::derivativeVectorTest() {
  auto derivatives = translation.der_vec(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
  for (auto i = 0u; i < derivatives.size(); ++i) {
    EXPECT_NEAR(derivatives.at(i), derivativeVector.at(i), doubleNearThreshold);
  }
}

std::string InternalCoordinatesTranslationYTest::returnInfoTest() {
  return translation.info(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

InternalCoordinatesTranslationZTest::InternalCoordinatesTranslationZTest() : InternalCoordinatesTest(), translation{ { 1u,2u,3u,4u,5u,6u } },
derivativeVector(3u * 12u, 0.) {
  auto constexpr sizeOfMethanolWhichIsDescribed = 3u * 12u / 2u;
  for (auto i = 2u; i < sizeOfMethanolWhichIsDescribed; i += 3) {
    derivativeVector.at(i) = 0.16666666666666666;
  }
}

double InternalCoordinatesTranslationZTest::testTranslationValue() {
  return translation.val(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

void InternalCoordinatesTranslationZTest::derivativeVectorTest() {
  auto derivatives = translation.der_vec(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
  for (auto i = 0u; i < derivatives.size(); ++i) {
    EXPECT_NEAR(derivatives.at(i), derivativeVector.at(i), doubleNearThreshold);
  }
}

std::string InternalCoordinatesTranslationZTest::returnInfoTest() {
  return translation.info(
    twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

TEST_F(InternalCoordinatesDistancesTest, testBondLength) {
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

TEST_F(InternalCoordinatesTranslationTest, testTranslationDerivativeTest) {
  testTranslationDerivativeTest();
}

TEST_F(InternalCoordinatesTranslationTest, hessianGuessTest) {
  EXPECT_NEAR(hessianGuessTest(), 0.05, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTranslationXTest, testTranslationValue) {
  EXPECT_NEAR(testTranslationValue(), -10.831910151124269, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTranslationXTest, derivativeVectorTest) {
  derivativeVectorTest();
}

TEST_F(InternalCoordinatesTranslationXTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Trans X: -10.8319");
}

TEST_F(InternalCoordinatesTranslationYTest, testTranslationValue) {
  EXPECT_NEAR(testTranslationValue(), -0.44786509173350525, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTranslationYTest, derivativeVectorTest) {
  derivativeVectorTest();
}

TEST_F(InternalCoordinatesTranslationYTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Trans Y: -0.447865");
}

TEST_F(InternalCoordinatesTranslationZTest, testTranslationValue) {
  auto bla = testTranslationValue();
  EXPECT_NEAR(testTranslationValue(), -0.44471554819107562, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTranslationZTest, derivativeVectorTest) {
  derivativeVectorTest();
}

TEST_F(InternalCoordinatesTranslationZTest, returnInfoTest) {
  auto bla = returnInfoTest();
  EXPECT_EQ(returnInfoTest(), "Trans Z: -0.444716");
}
#endif