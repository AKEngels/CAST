#ifdef GOOGLE_MOCK

#ifndef H_INTERNAL_COORDINATE_TEST
#define H_INTERNAL_COORDINATE_TEST

#include "../InternalCoordinates.h"
#include "../coords.h"
#include "../ic_core.h"

struct Molecule {
  coords::Representation_3D cartesianRepresentation;
  std::vector<std::string> elementSymbols;
};

struct MethanolMoleculesImpl {
  virtual Molecule & getOneRepresentation() = 0;
  virtual std::pair<Molecule&, Molecule&> getTwoRepresentations() = 0;
};

struct SubsystemOfTwoMethanolMolecules : MethanolMoleculesImpl {
  SubsystemOfTwoMethanolMolecules();
  Molecule & getOneRepresentation() override { return subSystem; }
  std::pair<Molecule&, Molecule&> getTwoRepresentations() override { throw std::runtime_error("Wrong Implementation. Use one of the other Methanol implementations."); }

  Molecule subSystem;
};

struct RotatetdMethanolMolecules : MethanolMoleculesImpl {
  RotatetdMethanolMolecules();
  Molecule & getOneRepresentation() override { throw std::runtime_error("Wrong Implementation. Use one of the other Methanol implementations."); }
  std::pair<Molecule&, Molecule&> getTwoRepresentations() override {
    return { initialMethanolSystem, rotatedMethanolSystem };
  }

  std::vector<std::string> elementSymbols;
  Molecule initialMethanolSystem;
  Molecule rotatedMethanolSystem;
};

class InternalCoordinatesTestSubsystem : public testing::Test {
public:
  InternalCoordinatesTestSubsystem()
      : twoMethanolMolecules{ std::make_unique<SubsystemOfTwoMethanolMolecules>() } {}

protected:
  std::unique_ptr<MethanolMoleculesImpl> twoMethanolMolecules;
};

class InternalCoordinatesTestRotatedMolecules : public testing::Test {
public:
  InternalCoordinatesTestRotatedMolecules()
    : twoMethanolMolecules{ std::make_unique<RotatetdMethanolMolecules>() } {}

protected:
  std::unique_ptr<MethanolMoleculesImpl> twoMethanolMolecules;
};

class InternalCoordinatesDistancesTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesDistancesTest();

  double testBondLength();

  std::pair<coords::r3, coords::r3> testBondDerivatives();

  void derivativeVectorTest();

  std::string returnInfoTest();

  coords::r3 const firstAtomDerivatives;
  coords::r3 const secondAtomDerivatives;

private:
  InternalCoordinates::BondDistance bond;
  std::vector<double> derivativeVector;
};

struct DifferentInternalCoordinates {
  std::shared_ptr<InternalCoordinates::InternalCoordinate> internalCoordinate;
  double expectedValue;
  friend std::ostream& operator<<(std::ostream & os, DifferentInternalCoordinates const& internalCoordinate) {
    return os << "Hessian Guess schould be: " << internalCoordinate.expectedValue;
  }
};

class InternalCoordinatesHessianTests : public testing::WithParamInterface<DifferentInternalCoordinates>, public testing::Test {
public:
  InternalCoordinatesHessianTests() : twoMethanolMolecules{ std::make_unique<RotatetdMethanolMolecules>() }, 
    internalCoordinate(GetParam().internalCoordinate) {}
  std::unique_ptr<MethanolMoleculesImpl> twoMethanolMolecules;

  std::shared_ptr<InternalCoordinates::InternalCoordinate> internalCoordinate;

};

class InternalCoordinatesAnglesTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesAnglesTest();

  double testAngleValue();

  std::tuple<coords::r3, coords::r3, coords::r3> testAngleDerivatives();

  void derivativeVectorTest();

  std::string returnInfoTest();

  coords::r3 const leftAtomsDerivative;
  coords::r3 const middleAtomsDerivative;
  coords::r3 const rightAtomsDerivative;

private:
  InternalCoordinates::BondAngle angle;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesDihedralsTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesDihedralsTest();

  double testDihedralValue();

  std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
    testDihedralDerivatives();

  void derivativeVectorTest();

  std::string returnInfoTest();

  coords::r3 leftLeftDerivative;
  coords::r3 leftMiddleDerivative;
  coords::r3 rightMiddleDerivative;
  coords::r3 rightRightDerivative;

private:
  InternalCoordinates::DihedralAngle dihedralAngle;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesTranslationTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesTranslationTest();

  //TODO Test it with the TEST_P thingy
  void testTranslationDerivativeTest();

private:
  InternalCoordinates::TranslationX translation;
};

class InternalCoordinatesTranslationXTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesTranslationXTest();

  double testTranslationValue();

  void derivativeVectorTest();

  std::string returnInfoTest();

private:
  InternalCoordinates::TranslationX translation;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesTranslationYTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesTranslationYTest();

  double testTranslationValue();

  void derivativeVectorTest();

  std::string returnInfoTest();

private:
  InternalCoordinates::TranslationY translation;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesTranslationZTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesTranslationZTest();

  double testTranslationValue();

  void derivativeVectorTest();

  std::string returnInfoTest();

private:
  InternalCoordinates::TranslationZ translation;
  std::vector<double> derivativeVector;
};

struct ExpectedValuesForRotations{
  double expectedValue;
  bool rotateMolecule;
  bool evaluateValues;
  bool evaluateDerivatives;
  std::vector<double> expectedDerivatives;
};

class InternalCoordinatesRotationsTest : public InternalCoordinatesTestRotatedMolecules, public testing::WithParamInterface<ExpectedValuesForRotations> {
public:
  InternalCoordinatesRotationsTest() : InternalCoordinatesTestRotatedMolecules(), cartesianCoordinates(twoMethanolMolecules->getTwoRepresentations()
    .first.cartesianRepresentation), rotations{ InternalCoordinates::Rotator::buildRotator(cartesianCoordinates, std::vector<std::size_t>{1,2,3,4,5,6})->makeRotations() }{}

  InternalCoordinates::CartesiansForInternalCoordinates cartesianCoordinates;
  InternalCoordinates::Rotations rotations;

  void checkIfVectorsAreSame(std::vector<double> const& lhs, std::vector<double> const& rhs);
};

class InternalCoordinatesRotationInfoTest : public InternalCoordinatesTestRotatedMolecules {
public:
  InternalCoordinatesRotationInfoTest() : InternalCoordinatesTestRotatedMolecules(), cartesianCoordinates(twoMethanolMolecules->getTwoRepresentations()
    .first.cartesianRepresentation), rotations{ InternalCoordinates::Rotator::buildRotator(cartesianCoordinates, std::vector<std::size_t>{1,2,3,4,5,6})->makeRotations() } {}

  std::string infoOfRotationA();
  std::string infoOfRotationB();
  std::string infoOfRotationC();

  InternalCoordinates::CartesiansForInternalCoordinates cartesianCoordinates;
  InternalCoordinates::Rotations rotations;
};

class InternalCoordinatesRotatorTest : public InternalCoordinatesTestRotatedMolecules {
public:
  InternalCoordinatesRotatorTest();
  void testRotationValue();
  void testRotationDerivatives();
  void testRadiusOfGyration();

private:
  InternalCoordinates::CartesiansForInternalCoordinates cartesianCoordinates;
  std::shared_ptr<InternalCoordinates::Rotator> rotation;
};

class CorrelationTests : public testing::Test {
public:
  //Make the input molecules not being in the origin
  CorrelationTests()
    : twoMethanolMolecules{ std::make_unique<RotatetdMethanolMolecules>() } {}

  void testExopentialMap();
  void testCorrelationMatrix();
  void testFMatrix();
  void testQuaternionForTwoMolecules();
  void testCorrelationMatrixDerivatives();
  void testFMatrixDerivatives();
  void testQuaternionDerivatives();
  void testExponentialMapDerivatives();

private:
  class ReadMatrixFiles {
  public:
    ReadMatrixFiles(std::istream & inputStream) : inputStream{ inputStream }{}
    scon::mathmatrix<double> readNLinesOfFileWithMNumbers(std::size_t const n, std::size_t const m);
  private:
    std::istream & inputStream;
  };
  scon::mathmatrix<double> readNextFderivative(std::istream & inputFileStream);
  scon::mathmatrix<double> readNextQuaternionDerivative(std::istream & inputFileStream);
  scon::mathmatrix<double> readNextExponentialMapderivative(std::istream & inputFileStream);
  std::unique_ptr<MethanolMoleculesImpl> twoMethanolMolecules;
};

class TranslationRotationCoordinatesTest : testing::Test {
private:
  std::unique_ptr<MethanolMoleculesImpl> twoMethanolMolecules;
};

class InterestedRotator : public InternalCoordinates::AbstractRotatorListener, public std::enable_shared_from_this<InterestedRotator>{
public:
  static std::shared_ptr<InterestedRotator> buildInterestedRotator(InternalCoordinates::CartesiansForInternalCoordinates & cartesians) ;
  void setAllFlag()override { updateFlag = true; }
  bool isFlagSet() { return updateFlag; }
private:
  InterestedRotator() : updateFlag{ false } {}
  void registerCartesians(InternalCoordinates::CartesiansForInternalCoordinates & cartesians);
  bool updateFlag;
};

class RotatorObserverTest : public testing::Test {
public:
  RotatorObserverTest();
  void testInitiallyFlagIsSetToFalse();
  void testWhenGeometryIsUpdatedThenFlagIsTrue();
private:
  InternalCoordinates::CartesiansForInternalCoordinates cartesianCoordinates;
  std::shared_ptr<InterestedRotator> rotator;
};

#endif

#endif