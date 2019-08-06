#ifdef GOOGLE_MOCK

#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "internal_coordinate_test.h"
#include "../../ic_util.h"
#include "../../ic_rotation.h"
#include "../../graph.h"
#include "../../Scon/scon_mathmatrix.h"
#include "ExpectedValuesForInternalCoordinatesTest.h"

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

using namespace ExpectedValuesForInternalCoordinates;

SubsystemOfTwoMethanolMolecules::SubsystemOfTwoMethanolMolecules()
  : subSystem{ {createSystemOfTwoMethanolMolecules() }, 
                createSequenceOfSymbolsForTwoMethanolMolecules()} {
  subSystem.cartesianRepresentation /= energy::bohr2ang;
}

RotatetdMethanolMolecules::RotatetdMethanolMolecules()
    : initialMethanolSystem{ createInitialMethanolForRotationSystem(),
  createSequenceOfSymbolsForInitialMethanolForRotationSystem() },
      rotatedMethanolSystem{ createRotatedMethanolForRotationSystem(),
  createSequenceOfSymbolsForRotatedMethanolForRotationSystem() } {
  initialMethanolSystem.cartesianRepresentation /= energy::bohr2ang;
  rotatedMethanolSystem.cartesianRepresentation /= energy::bohr2ang;
}

InternalCoordinatesDistancesTest::InternalCoordinatesDistancesTest()
    : InternalCoordinatesTestSubsystem(), firstAtomDerivatives{ firstBondAtomDerivative() },
      secondAtomDerivatives{ secondBondAtomDerivative() },
  bond(ic_util::Node{ 1, "C", "C", coords::Cartesian_Point() }, ic_util::Node{ 2,"H","H", coords::Cartesian_Point() }),
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

std::string InternalCoordinatesDistancesTest::returnInfoTest() {
  return bond.info(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

InternalCoordinatesAnglesTest::InternalCoordinatesAnglesTest()
    : InternalCoordinatesTestSubsystem(), leftAtomsDerivative{ leftAngleAtomDerivative() },
      middleAtomsDerivative{ middleAngleAtomDerivative() },
      rightAtomsDerivative{ rightAngleAtomDerivative() },
      angle(ic_util::Node{ 1, "H", "H", coords::Cartesian_Point() }, ic_util::Node{ 2, "C", "C", coords::Cartesian_Point() }, 
				ic_util::Node{ 3, "H", "H", coords::Cartesian_Point() }),
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

std::string InternalCoordinatesAnglesTest::returnInfoTest() {
  return angle.info(twoMethanolMolecules->getOneRepresentation().cartesianRepresentation);
}

InternalCoordinatesDihedralsTest::InternalCoordinatesDihedralsTest()
    : InternalCoordinatesTestSubsystem(), leftLeftDerivative{ leftLeftDihedralDerivative() },
      leftMiddleDerivative{ leftMiddleDihedralDerivative() },
      rightMiddleDerivative{ rightMiddleDihedralDerivative() },
      rightRightDerivative{ rightRightDihedralDerivative() },
  dihedralAngle(ic_util::Node{ 1, std::string(), std::string(), coords::Cartesian_Point() }, 
		ic_util::Node{ 2, std::string(), std::string(), coords::Cartesian_Point() },
		ic_util::Node{ 3, std::string(), std::string(), coords::Cartesian_Point() },
		ic_util::Node{ 4, std::string(), std::string(), coords::Cartesian_Point() }),
	derivativeVector(3u * 12u, 0.) {
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

std::string InternalCoordinatesDihedralsTest::returnInfoTest() {
  return dihedralAngle.info(
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

InternalCoordinatesRotatorTest::InternalCoordinatesRotatorTest()
    : InternalCoordinatesTestRotatedMolecules(), cartesianCoordinates(twoMethanolMolecules->getTwoRepresentations()
      .first.cartesianRepresentation),
      rotation(InternalCoordinates::Rotator::buildRotator(cartesianCoordinates,
        std::vector<std::size_t>{ 1u,2u,3u,4u,5u,6u })) {}

void InternalCoordinatesRotatorTest::testRotationValue(){
  auto rotationsForXyz = rotation->valueOfInternalCoordinate(twoMethanolMolecules->getTwoRepresentations().second.cartesianRepresentation);
  std::array<double, 3u> expectedValues = { -3.1652558307984515, -2.4287895503201611, 3.1652568986800538 };
  for (auto i = 0u; i < rotationsForXyz.size(); ++i) {
    EXPECT_NEAR(rotationsForXyz.at(i), expectedValues.at(i), doubleNearThreshold);
  }
}

void InternalCoordinatesRotatorTest::testRotationDerivatives(){
  auto rotationDerivatives = rotation->rot_der_mat(twoMethanolMolecules->getTwoRepresentations().second.cartesianRepresentation);
  EXPECT_EQ(rotationDerivatives, expectedRotationDerivatives());
}


void InternalCoordinatesRotatorTest::testRadiusOfGyration() {
  EXPECT_NEAR(ic_util::rad_gyr(ExpectedValuesForInternalCoordinates::createInitialMethanolForRotationSystem() / energy::bohr2ang), 2.2618755203155767, doubleNearThreshold);
}

void InternalCoordinatesRotationsTest::checkIfVectorsAreSame(std::vector<double> const & lhs, std::vector<double> const & rhs){
  ASSERT_EQ(lhs.size(), rhs.size());
  for (auto i = 0u; i < lhs.size(); ++i) {
    EXPECT_NEAR(lhs.at(i), rhs.at(i), doubleNearThreshold);
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

  scon::mathmatrix<double> expectedValues = expectedValuesForF();

  EXPECT_EQ(F, expectedValues);

  auto eigenVectorsAndValuesOfF = F.eigensym();

  scon::mathmatrix<double> expectedEigenValues = expectedEigenValuesForF();

  scon::mathmatrix<double> expectedEigenVectors = expectedEigenVectorsForF();

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
  auto Fderivatives = ic_rotation::F_matrix_derivs(twoMethanolMolecules->getTwoRepresentations().first.cartesianRepresentation);
  auto expectedFderivatives = provideExpectedValuesForFMatrixDerivativesInFile();
  for (auto i = 0u; i < Fderivatives.size(); ++i) {
    for (auto j = 0u; j < Fderivatives.at(i).size(); ++j) {
      EXPECT_EQ(Fderivatives.at(i).at(j), expectedFderivatives.at(i * Fderivatives.at(i).size() + j));
    }
  }
}

void CorrelationTests::testQuaternionDerivatives(){
  auto quaternionDerivatives = provideExpectedValuesForQuaternionDerivatives();
  for (auto i = 0u; i < quaternionDerivatives.size(); ++i) {
    EXPECT_EQ(quaternionDerivatives.at(i), quaternionDerivatives.at(i));
  }
}

void CorrelationTests::testExponentialMapDerivatives(){
  auto exponentialMapDerivatives = provideExpectedValuesForExponentialMapDerivatives();
  for (auto i = 0u; i < exponentialMapDerivatives.size(); ++i) {
    EXPECT_EQ(exponentialMapDerivatives.at(i), exponentialMapDerivatives.at(i));
  }

}

scon::mathmatrix<double> CorrelationTests::readNextFderivative(std::istream & inputFileStream){
  auto constexpr sizeOfDerivatives = 4u;
  return ReadMatrixFiles(inputFileStream).readNLinesOfFileWithMNumbers(sizeOfDerivatives, sizeOfDerivatives);
}

scon::mathmatrix<double> CorrelationTests::readNextQuaternionDerivative(std::istream & inputFileStream) {
  auto constexpr sizeOfDerivativeRows = 3u;
  auto constexpr sizeOfDerivativeCols = 4u;
  return ReadMatrixFiles(inputFileStream).readNLinesOfFileWithMNumbers(sizeOfDerivativeRows, sizeOfDerivativeCols);
}

scon::mathmatrix<double> CorrelationTests::readNextExponentialMapderivative(std::istream & inputFileStream) {
  auto constexpr sizeOfDerivatives = 3u;
  return ReadMatrixFiles(inputFileStream).readNLinesOfFileWithMNumbers(sizeOfDerivatives, sizeOfDerivatives);
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

TEST_F(InternalCoordinatesDistancesTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Bond: 2.64142 (a. u.) 1.39778 (Angstrom) || 1 || 2 || Constrained: false");
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

TEST_F(InternalCoordinatesAnglesTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Angle: 30.3101 || 1 || 2 || 3 || Constrained: false");
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

TEST_F(InternalCoordinatesDihedralsTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Dihedral: 32.2818 || 1 || 2 || 3 || 4 || Constrained: false");
}

/*TEST_F(InternalCoordinatesTranslationTest, testTranslationDerivativeTest) {
  testTranslationDerivativeTest();
}*/

TEST_F(InternalCoordinatesTranslationXTest, testTranslationValue) {
  EXPECT_NEAR(testTranslationValue(), -10.831910151124269, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTranslationXTest, derivativeVectorTest) {
  derivativeVectorTest();
}

TEST_F(InternalCoordinatesTranslationXTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Trans X: -10.8319 | Constrained: false");
}

TEST_F(InternalCoordinatesTranslationYTest, testTranslationValue) {
  EXPECT_NEAR(testTranslationValue(), -0.44786509173350525, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTranslationYTest, derivativeVectorTest) {
  derivativeVectorTest();
}

TEST_F(InternalCoordinatesTranslationYTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Trans Y: -0.447865 | Constrained: false");
}

TEST_F(InternalCoordinatesTranslationZTest, testTranslationValue) {
  EXPECT_NEAR(testTranslationValue(), -0.44471554819107562, doubleNearThreshold);
}

TEST_F(InternalCoordinatesTranslationZTest, derivativeVectorTest) {
  derivativeVectorTest();
}

TEST_F(InternalCoordinatesTranslationZTest, returnInfoTest) {
  EXPECT_EQ(returnInfoTest(), "Trans Z: -0.444716 | Constrained: false");
}

TEST_F(InternalCoordinatesRotatorTest, testRadiusOfGyration) {
  testRadiusOfGyration();
}

TEST_F(InternalCoordinatesRotatorTest, testRotationValue) {
  testRotationValue();
}

TEST_F(InternalCoordinatesRotatorTest, testRotationDerivatives) {
  testRotationDerivatives();
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

TEST_F(CorrelationTests, testQuaternionDerivatives) {
  testQuaternionDerivatives();
}

TEST_F(CorrelationTests, testExponentialMapDerivatives) {
  testExponentialMapDerivatives();
}

TEST_P(InternalCoordinatesHessianTests, testHessianGuessesForAllInternalCoordinates) {
  auto const& expectedValues = GetParam();
  internalCoordinate->hessian_guess(twoMethanolMolecules->getTwoRepresentations().first.cartesianRepresentation);
  EXPECT_NEAR(expectedValues.expectedValue, internalCoordinate->hessian_guess(twoMethanolMolecules->getTwoRepresentations().first.cartesianRepresentation), doubleNearThreshold);
}

INSTANTIATE_TEST_CASE_P(BondDistances, InternalCoordinatesHessianTests, testing::Values(
  
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondDistance>(ic_util::Node{ 1,"H","H", coords::Cartesian_Point() }, ic_util::Node{ 2,"H","H", coords::Cartesian_Point() }), 0.072180537605377640 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondDistance>(ic_util::Node{ 1,"C","C", coords::Cartesian_Point() }, ic_util::Node{ 2,"H","H", coords::Cartesian_Point() }), 0.14450082351214560 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondDistance>(ic_util::Node{ 1,"H","H", coords::Cartesian_Point() }, ic_util::Node{ 2,"C","C", coords::Cartesian_Point() }), 0.14450082351214560 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondDistance>(ic_util::Node{ 1,"Si","Si", coords::Cartesian_Point() }, ic_util::Node{ 2,"H","H", coords::Cartesian_Point() }), 0.22290342750129621 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondDistance>(ic_util::Node{ 1,"H","H", coords::Cartesian_Point() }, ic_util::Node{ 2,"Si","Si", coords::Cartesian_Point() }), 0.22290342750129621 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondDistance>(ic_util::Node{ 1,"Si","Si", coords::Cartesian_Point() }, ic_util::Node{ 2,"C","C", coords::Cartesian_Point() }), 1.2361326955675498 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondDistance>(ic_util::Node{ 1,"C","C", coords::Cartesian_Point() }, ic_util::Node{ 2,"Si","Si", coords::Cartesian_Point() }), 1.2361326955675498 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondDistance>(ic_util::Node{ 1,"C","C", coords::Cartesian_Point() }, ic_util::Node{ 2,"C","C", coords::Cartesian_Point() }), 0.45990191969419802 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondDistance>(ic_util::Node{ 1,"Si","Si", coords::Cartesian_Point() }, ic_util::Node{ 2,"Si","Si", coords::Cartesian_Point() }), 9.1964705962169191 }
));

INSTANTIATE_TEST_CASE_P(BondAndDihedralAngles, InternalCoordinatesHessianTests, testing::Values(
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondAngle>(ic_util::Node{ 1, "H", "H", coords::Cartesian_Point() }, ic_util::Node{ 2, "C", "C", coords::Cartesian_Point() }, ic_util::Node{ 3, "H", "H", coords::Cartesian_Point() }), 0.16 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondAngle>(ic_util::Node{ 1, "H", "H", coords::Cartesian_Point() }, ic_util::Node{ 2, "C", "C", coords::Cartesian_Point() }, ic_util::Node{ 3, "C", "C", coords::Cartesian_Point() }), 0.16 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondAngle>(ic_util::Node{ 1, "C", "C", coords::Cartesian_Point() }, ic_util::Node{ 2, "C", "C", coords::Cartesian_Point() }, ic_util::Node{ 3, "H", "H", coords::Cartesian_Point() }), 0.16 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::BondAngle>(ic_util::Node{ 1, "C", "C", coords::Cartesian_Point() }, ic_util::Node{ 2, "C", "C", coords::Cartesian_Point() }, ic_util::Node{ 3, "C", "C", coords::Cartesian_Point() }), 0.25 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::DihedralAngle>(ic_util::Node{1, std::string(), std::string(), coords::Cartesian_Point()}, 
		ic_util::Node{2, std::string(), std::string(), coords::Cartesian_Point() }, ic_util::Node{3, std::string(), std::string(), coords::Cartesian_Point()}, 
		ic_util::Node{4, std::string(), std::string(), coords::Cartesian_Point()}), 0.023 }
));

INSTANTIATE_TEST_CASE_P(Translations, InternalCoordinatesHessianTests, testing::Values(
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::TranslationX>(std::vector<std::size_t>{1,2,3}), 0.05 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::TranslationY>(std::vector<std::size_t>{1,2,3}), 0.05 },
  DifferentInternalCoordinates{ std::make_shared<InternalCoordinates::TranslationZ>(std::vector<std::size_t>{1,2,3}), 0.05 }
));

INSTANTIATE_TEST_CASE_P(Rotations, InternalCoordinatesHessianTests, testing::Values(
  DifferentInternalCoordinates{ rotations(Rotation::A), 0.05 },
  DifferentInternalCoordinates{ rotations(Rotation::B), 0.05 }
));

std::unique_ptr<InternalCoordinates::InternalCoordinate> & InternalCoordinatesRotationsTest::getRotation(ExpectedValuesForInternalCoordinates::Rotation const kindOfRotation) {
  if (kindOfRotation == Rotation::A) {
    return rotations.rotationA;
  }
  else if (kindOfRotation == Rotation::B) {
    return rotations.rotationB;
  }
  else if (kindOfRotation == Rotation::C) {
    return rotations.rotationC;
  }
	else throw std::runtime_error("Error: strange kind of rotation.");
}

TEST_P(InternalCoordinatesRotationsTest, testValuesForAllRotations) {

  auto const& rotation = getRotation(GetParam().kindOfRotation);

  auto const& expectedValue = GetParam().expectedValue;
  if (GetParam().rotateMolecule) {
    cartesianCoordinates.setCartesianCoordnates(twoMethanolMolecules->getTwoRepresentations()
      .second.cartesianRepresentation);
  }
  if (GetParam().evaluateValues) {
    EXPECT_NEAR(rotation->val(cartesianCoordinates), expectedValue, doubleNearThreshold);
  }
  if (GetParam().evaluateDerivatives) {
    checkIfVectorsAreSame(rotation->der_vec(cartesianCoordinates),
      GetParam().expectedDerivatives);
  }
  EXPECT_EQ(GetParam().evaluateValues, rotations.rotator->areValuesUpToDate());
  EXPECT_EQ(GetParam().evaluateDerivatives, rotations.rotator->areDerivativesUpToDate());
}

INSTANTIATE_TEST_CASE_P(RotationA, InternalCoordinatesRotationsTest, testing::Values(
  ExpectedValuesForRotations{ Rotation::A, 0.0, false, false, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::A, 0.0, true, false, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::A, 0.0, false, true, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::A, -3.1652558307984515, true, true, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::A, 0.0, false, true, true, expectedDerivativesRotationANoRotation() },
  ExpectedValuesForRotations{ Rotation::A, -3.1652558307984515, true, true, true, expectedDerivativesRotatoionA() }
));

INSTANTIATE_TEST_CASE_P(RotationB, InternalCoordinatesRotationsTest, testing::Values(
  ExpectedValuesForRotations{ Rotation::B, 0.0, false, false, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::B, 0.0, true, false, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::B, 0.0, false, true, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::B, -2.4287895503201611, true, true, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::B, 0.0, false, true, true, expectedDerivativesRotationBNoRotation() },
  ExpectedValuesForRotations{ Rotation::B, -2.4287895503201611, true, true, true, expectedDerivativesRotatoionB() }
));

INSTANTIATE_TEST_CASE_P(RotationC, InternalCoordinatesRotationsTest, testing::Values(
  ExpectedValuesForRotations{ Rotation::C, 0.0, false, false, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::C, 0.0, true, false, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::C, 0.0, false, true, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::C, 3.1652568986800538, true, true, false, notCalculatingDerivatives() },
  ExpectedValuesForRotations{ Rotation::C, 0.0, false, true, true, expectedDerivativesRotationCNoRotation() },
  ExpectedValuesForRotations{ Rotation::C, 3.1652568986800538, true, true, true, expectedDerivativesRotatoionC() }
));

std::string InternalCoordinatesRotationInfoTest::infoOfRotationA(){
  return rotations.rotationA->info(twoMethanolMolecules->getTwoRepresentations().second.cartesianRepresentation);
}

std::string InternalCoordinatesRotationInfoTest::infoOfRotationB() {
  return rotations.rotationB->info(twoMethanolMolecules->getTwoRepresentations().second.cartesianRepresentation);
}

std::string InternalCoordinatesRotationInfoTest::infoOfRotationC() {
  return rotations.rotationC->info(twoMethanolMolecules->getTwoRepresentations().second.cartesianRepresentation);
}

TEST_F(InternalCoordinatesRotationInfoTest, infoOfRotationA) {
  EXPECT_EQ(infoOfRotationA(), "Rotation A: -3.16526 | Constrained: false");
}

TEST_F(InternalCoordinatesRotationInfoTest, infoOfRotationB) {
  EXPECT_EQ(infoOfRotationB(), "Rotation B: -2.42879 | Constrained: false");
}

TEST_F(InternalCoordinatesRotationInfoTest, infoOfRotationC) {
  EXPECT_EQ(infoOfRotationC(), "Rotation C: 3.16526 | Constrained: false");
}

RotatorObserverTest::RotatorObserverTest() : cartesianCoordinates{ ExpectedValuesForInternalCoordinates::createInitialMethanolForRotationSystem() }, rotator{ InterestedRotator::buildInterestedRotator(cartesianCoordinates) } {}

void RotatorObserverTest::testInitiallyFlagIsSetToFalse(){
  EXPECT_EQ(rotator->isFlagSet(), false);
}

void RotatorObserverTest::testWhenGeometryIsUpdatedThenFlagIsTrue(){
  auto changedStructure = ExpectedValuesForInternalCoordinates::createRotatedMethanolForRotationSystem();
  cartesianCoordinates.setCartesianCoordnates(std::move(changedStructure));
  EXPECT_EQ(rotator->isFlagSet(), true);
}

std::shared_ptr<InterestedRotator> InterestedRotator::buildInterestedRotator(InternalCoordinates::CartesiansForInternalCoordinates & cartesians) {
  auto newInstance = std::make_shared<InterestedRotator>(InterestedRotator());
  newInstance->registerCartesians(cartesians);
  return newInstance;
}

void InterestedRotator::registerCartesians(InternalCoordinates::CartesiansForInternalCoordinates & cartesians) {
  auto observer = std::make_shared<InternalCoordinates::RotatorObserver>();
  observer->setNewRotator(shared_from_this());
  cartesians.registerObserver(observer);
}

TEST_F(RotatorObserverTest, testInitiallyFlagIsSetToFalse) {
  testInitiallyFlagIsSetToFalse();
}

TEST_F(RotatorObserverTest, testWhenGeometryIsUpdatedThenFlagIsTrue) {
  testWhenGeometryIsUpdatedThenFlagIsTrue();
}
#endif
