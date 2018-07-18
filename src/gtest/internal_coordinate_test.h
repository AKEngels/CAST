#ifdef GOOGLE_MOCK

#ifndef H_INTERNAL_COORDINATE_TEST
#define H_INTERNAL_COORDINATE_TEST

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

  double hessianGuessTest();

  std::string returnInfoTest();

  coords::r3 const firstAtomDerivatives;
  coords::r3 const secondAtomDerivatives;

private:
  ic_core::BondDistance bond;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesAnglesTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesAnglesTest();

  double testAngleValue();

  std::tuple<coords::r3, coords::r3, coords::r3> testAngleDerivatives();

  void derivativeVectorTest();

  double hessianGuessTest();

  std::string returnInfoTest();

  coords::r3 const leftAtomsDerivative;
  coords::r3 const middleAtomsDerivative;
  coords::r3 const rightAtomsDerivative;

private:
  ic_core::BondAngle angle;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesDihedralsTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesDihedralsTest();

  double testDihedralValue();

  std::tuple<coords::r3, coords::r3, coords::r3, coords::r3>
    testDihedralDerivatives();

  void derivativeVectorTest();

  double hessianGuessTest();

  std::string returnInfoTest();

  coords::r3 leftLeftDerivative;
  coords::r3 leftMiddleDerivative;
  coords::r3 rightMiddleDerivative;
  coords::r3 rightRightDerivative;

private:
  ic_core::DihedralAngle dihedralAngle;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesTranslationTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesTranslationTest();

  void testTranslationDerivativeTest();

  double hessianGuessTest();

private:
  ic_core::TranslationX translation;
};

class InternalCoordinatesTranslationXTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesTranslationXTest();

  double testTranslationValue();

  void derivativeVectorTest();

  std::string returnInfoTest();

private:
  ic_core::TranslationX translation;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesTranslationYTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesTranslationYTest();

  double testTranslationValue();

  void derivativeVectorTest();

  std::string returnInfoTest();

private:
  ic_core::TranslationY translation;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesTranslationZTest : public InternalCoordinatesTestSubsystem {
public:
  InternalCoordinatesTranslationZTest();

  double testTranslationValue();

  void derivativeVectorTest();

  std::string returnInfoTest();

private:
  ic_core::TranslationZ translation;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesRotationTest : public InternalCoordinatesTestRotatedMolecules {
public:
  InternalCoordinatesRotationTest();
  double testRadiusOfGyration();
  void testRotationValue();
  void testRotationDerivatives();

private:
  ic_core::Rotation rotation;
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

#endif

#endif