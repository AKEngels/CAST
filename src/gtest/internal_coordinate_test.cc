#ifdef GOOGLE_MOCK

#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "internal_coordinate_test.h"
#include "../ic_util.h"
#include "../ic_rotation.h"
#include "../scon_mathmatrix.h"
#include "TestFiles/WriteOutputForFMatrixDerivatives.h"

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
                                 coords::r3{ -0.11699, -0.10048, 0.32100 },
                                 coords::r3{ 0.08440, -0.24140, -1.05500 },
                                 coords::r3{ -0.74011, 0.79636, 0.53000 },
                                 coords::r3{ 0.85731, -0.03439, 0.85000 },
                                 coords::r3{ -0.65548, -0.99354, 0.69900 },
                                 coords::r3{ 0.57087, 0.57345, -1.34500 },
                             },
                             { "C", "O", "H", "H", "H", "H" } } {
  initialMethanolSystem.cartesianRepresentation /= energy::bohr2ang;
  rotatedMethanolSystem.cartesianRepresentation /= energy::bohr2ang;
}

InternalCoordinatesDistancesTest::InternalCoordinatesDistancesTest()
    : InternalCoordinatesTestSubsystem(), firstAtomDerivatives{ -0.98441712304088669,
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
    : InternalCoordinatesTestSubsystem(), leftAtomsDerivative{ -0.062056791850036874,
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
    : InternalCoordinatesTestSubsystem(), leftLeftDerivative{ 0.047760930904702938,
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

InternalCoordinatesTranslationTest::InternalCoordinatesTranslationTest() : InternalCoordinatesTestSubsystem(), translation{ { 1u,2u,3u,4u,5u,6u } } {}

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

InternalCoordinatesTranslationXTest::InternalCoordinatesTranslationXTest() : InternalCoordinatesTestSubsystem(), translation{ {1u,2u,3u,4u,5u,6u} },
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

InternalCoordinatesTranslationYTest::InternalCoordinatesTranslationYTest() : InternalCoordinatesTestSubsystem(), translation{ { 1u,2u,3u,4u,5u,6u } },
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

InternalCoordinatesTranslationZTest::InternalCoordinatesTranslationZTest() : InternalCoordinatesTestSubsystem(), translation{ { 1u,2u,3u,4u,5u,6u } },
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

InternalCoordinatesRotationTest::InternalCoordinatesRotationTest()
    : InternalCoordinatesTestRotatedMolecules(),
      rotation(twoMethanolMolecules->getTwoRepresentations()
                   .first.cartesianRepresentation,
        std::vector<std::size_t>{ 1, 2, 3, 4, 5, 6 }) {}

double InternalCoordinatesRotationTest::testRadiusOfGyration()
{
  return rotation.radiusOfGyration(twoMethanolMolecules->getTwoRepresentations()
    .first.cartesianRepresentation);
}

void InternalCoordinatesRotationTest::testRotationValue(){
  auto rotationsForXyz = rotation.valueOfInternalCoordinate(twoMethanolMolecules->getTwoRepresentations().second.cartesianRepresentation);
  std::array<double, 3u> expectedValues = { -3.1652558307984515, -2.4287895503201611, 3.1652568986800538 };
  for (auto i = 0u; i < rotationsForXyz.size(); ++i) {
    EXPECT_NEAR(rotationsForXyz.at(i), expectedValues.at(i), doubleNearThreshold);
  }
}

void CorrelationTests::testCorrelationMatrix() {
  auto correlationMatrix = ic_rotation::correlation_matrix(
      twoMethanolMolecules->getTwoRepresentations()
          .first.cartesianRepresentation,
      twoMethanolMolecules->getTwoRepresentations()
          .second.cartesianRepresentation);

  
  scon::mathmatrix<double> expectedValues{
    { 3.62870631512344, 6.77810459455176, -3.11427207160589 },
    { 3.03723523986651, -0.65557765262001, -7.29126193241831 },
    { -16.1308283716317, -2.71895876607409, 3.87293765244782 }
  };
  EXPECT_EQ(correlationMatrix, expectedValues);
}

void CorrelationTests::testFMatrix() {

  // Maybe Mock the F_matrix method to not calculate the correlation Matrix
  // explicitly?
  auto F = ic_rotation::F_matrix(twoMethanolMolecules->getTwoRepresentations()
                                     .first.cartesianRepresentation,
                                 twoMethanolMolecules->getTwoRepresentations()
                                     .second.cartesianRepresentation);

  scon::mathmatrix<double> expectedValues{
    { 6.84606631495125, - 4.57230316634422, - 13.0165563000258,  3.74086935468525 },
    { -4.57230316634422,  0.41134631529563,  9.81533983441827, - 19.2451004432376 },
    { -13.0165563000258,  9.81533983441827, - 8.15722162019127, - 10.0102206984924 },
    { 3.74086935468525, - 19.2451004432376, - 10.0102206984924, 0.899808989944389 }
  };
  EXPECT_EQ(F, expectedValues);

  auto eigenVectorsAndValuesOfF = F.eigensym();

  scon::mathmatrix<double> expectedEigenValues{ { -18.7538392250112,},
                                                { -18.0459463050214,},
                                                { 6.10328376463374,},
                                                { 30.6965017653989,} };

  scon::mathmatrix<double> expectedEigenVectors{
    { -0.199380447172903, 0.342091987934248, -0.811125403595648,
      0.430460321885915 },
    { -0.549673535066269, -0.495167463863083, -0.371435099023426,
      -0.560984986245274 },
    { -0.401611672479942, 0.783069660011025, 0.200535739330358,
      -0.430459509535228 },
    { -0.70485069813454, -0.156793373880449, 0.404841900136933,
      0.56098517550819 }
  };

  EXPECT_EQ(eigenVectorsAndValuesOfF.first, expectedEigenValues);
  EXPECT_EQ(eigenVectorsAndValuesOfF.second, expectedEigenVectors);

}

void CorrelationTests::testExopentialMap() {
  auto exponentialMap = ic_rotation::exponential_map(twoMethanolMolecules->getTwoRepresentations().first.cartesianRepresentation, twoMethanolMolecules->getTwoRepresentations().second.cartesianRepresentation);
  std::array<double, 3u> expectedValues = { -1.3993943532121675, -1.0737945251652472, 1.3993948253343480 };
  for (auto i = 0u; i < exponentialMap.size(); ++i) {
    EXPECT_NEAR(exponentialMap.at(i), expectedValues.at(i), doubleNearThreshold);
  }
}

void CorrelationTests::testQuaternionForTwoMolecules() {
  auto quaternion = ic_rotation::quaternion(twoMethanolMolecules->getTwoRepresentations().first.cartesianRepresentation, twoMethanolMolecules->getTwoRepresentations().second.cartesianRepresentation);

  std::array<double, 4u> expectedValuesForTheQuaternion{ 0.43046032188591526, -0.56098498624527438, -0.43045950953522794, 0.56098517550818972 };
  for (auto i = 0u; i < expectedValuesForTheQuaternion.size(); ++i) {
    EXPECT_NEAR(quaternion.second.at(i), expectedValuesForTheQuaternion.at(i), doubleNearThreshold);
  }

  double expectedValueForTheHighestEigenvalueOfF = 30.696501765398882;
  EXPECT_NEAR(quaternion.first, expectedValueForTheHighestEigenvalueOfF, doubleNearThreshold);
}

void CorrelationTests::testCorrelationMatrixDerivatives() {
  auto const& cartesians = twoMethanolMolecules->getTwoRepresentations().first.cartesianRepresentation;
  auto derives = ic_rotation::correlation_matrix_derivs(cartesians);

  std::vector<scon::mathmatrix<double> > containerForRowsOfCartesianRepresentation;
  
  for (auto const& atom : cartesians) {
    containerForRowsOfCartesianRepresentation.emplace_back(scon::mathmatrix<double>{ {atom.x(), atom.y(), atom.z()} });
  }

  auto constexpr numberOfCartesiansPerRow = 3u;

  for (auto i = 0u; i < containerForRowsOfCartesianRepresentation.size(); ++i) {
    auto & expectedRow = containerForRowsOfCartesianRepresentation.at(i);
    for (auto j = 0u; j < numberOfCartesiansPerRow; ++j) {
      EXPECT_EQ(derives.at(i).at(j).row(j), expectedRow);
    }
  }

}

void CorrelationTests::testFMatrixDerivatives(){
  provideExpectedValuesForFMatrixDerivativesInFile("DerivativesOfF.dat");
  std::ifstream matrixFileWithDerivatives("DerivativesOfF.dat");
  auto Fderivatives = ic_rotation::F_matrix_derivs(twoMethanolMolecules->getTwoRepresentations().first.cartesianRepresentation);

  for (auto i = 0u; i < Fderivatives.size(); ++i) {
    for (auto j = 0u; j < Fderivatives.at(i).size(); ++j) {
      auto expectedDerivative = readNextFderivative(matrixFileWithDerivatives);
      std::cout << expectedDerivative << std::endl;
      EXPECT_EQ(Fderivatives.at(i).at(j), expectedDerivative);
    }
  }
}

scon::mathmatrix<double> CorrelationTests::readNextFderivative(std::istream & inputFileStream)
{
  auto constexpr sizeOfDerivatives = 4u;
  scon::mathmatrix<double> nextDerivative(sizeOfDerivatives, sizeOfDerivatives);

  for (auto i = 0u; i < sizeOfDerivatives; ++i) {
    std::string line;
    std::getline(inputFileStream, line);
    std::stringstream lineStream(line);
    double row1, row2, row3, row4;

    lineStream >> nextDerivative(i,0) >> nextDerivative(i, 1) >> nextDerivative(i, 2) >> nextDerivative(i, 3);
  }

  return nextDerivative;
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
  EXPECT_NEAR(testTranslationValue(), -0.44471554819107562, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTranslationZTest, derivativeVectorTest) {
  derivativeVectorTest();
}

TEST_F(InternalCoordinatesTranslationZTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Trans Z: -0.444716");
}

TEST_F(InternalCoordinatesRotationTest, testRadiusOfGyration) {
  EXPECT_NEAR(testRadiusOfGyration(), 2.2618755203155767, doubleNearThreshold);
}

TEST_F(InternalCoordinatesRotationTest, testRotationValue) {
  testRotationValue();
}

TEST_F(CorrelationTests, testCorrelationMatrix) {
  testCorrelationMatrix();
}

TEST_F(CorrelationTests, testExopentialMap) {
  testExopentialMap();
}

TEST_F(CorrelationTests, testFMatrix) {
  testFMatrix();
}

TEST_F(CorrelationTests, testQuaternionForTwoMolecules) {
  testQuaternionForTwoMolecules();
}

TEST_F(CorrelationTests, testCorrelationMatrixDerivatives) {
  testCorrelationMatrixDerivatives();
}

TEST_F(CorrelationTests, testFMatrixDerivatives) {
  testFMatrixDerivatives();
}

#endif


