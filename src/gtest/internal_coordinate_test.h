#ifdef GOOGLE_MOCK

#ifndef H_INTERNAL_COORDINATE_TEST
#define H_INTERNAL_COORDINATE_TEST

#include "../coords.h"
#include "../ic_core.h"

struct TwoMethanolMolecules {
  TwoMethanolMolecules();

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

class InternalCoordinatesAnglesTest : public InternalCoordinatesTest {
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

class InternalCoordinatesDihedralsTest : public InternalCoordinatesTest {
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

class InternalCoordinatesTranslationTest : public InternalCoordinatesTest {
public:
  InternalCoordinatesTranslationTest();

  void testTranslationDerivativeTest();

  double hessianGuessTest();

private:
  ic_core::TranslationX translation;
};

class InternalCoordinatesTranslationXTest : public InternalCoordinatesTest {
public:
  InternalCoordinatesTranslationXTest();

  double testTranslationValue();

  void derivativeVectorTest();

  std::string returnInfoTest();

private:
  ic_core::TranslationX translation;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesTranslationYTest : public InternalCoordinatesTest {
public:
  InternalCoordinatesTranslationYTest();

  double testTranslationValue();

  void derivativeVectorTest();

  std::string returnInfoTest();

private:
  ic_core::TranslationY translation;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesTranslationZTest : public InternalCoordinatesTest {
public:
  InternalCoordinatesTranslationZTest();

  double testTranslationValue();

  void derivativeVectorTest();

  std::string returnInfoTest();

private:
  ic_core::TranslationZ translation;
  std::vector<double> derivativeVector;
};

class InternalCoordinatesRotationTest : public InternalCoordinatesTest {
public:
private:
  ic_core::Rotation rotation;
};

#endif

#endif