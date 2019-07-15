#ifdef GOOGLE_MOCK

#include "primitive_internals_test.h"
#include"../../InternalCoordinates.h"
#include"../../graph.h"
#include "../../InternalCoordinateDecorator.h"

using namespace ExpectedValuesForInternalCoordinates;

namespace {

  ic_util::Graph<ic_util::Node> createTestGraph(){
    std::vector<ic_util::Node> atomVector{ 
      ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 2, "O", "O" }, ic_util::Node{ 3, "H", "H" }, ic_util::Node{ 4, "H", "H" }, ic_util::Node{ 5, "H", "H" }, ic_util::Node{ 6, "H", "H" },
      ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 8, "O", "O" }, ic_util::Node{ 9, "H", "H" }, ic_util::Node{ 10, "H", "H" }, ic_util::Node{ 11, "H", "H" }, ic_util::Node{ 12, "H", "H" }
    };
    std::vector<std::pair<std::size_t, std::size_t>> connectivity{ 
      { 0u, 1u },{ 0u, 2u },{ 0u, 3u },{ 0u, 4u },{ 1u, 5u },
      { 6u, 7u },{ 6u, 8u },{ 6u, 9u },{ 6u, 10u },{ 7u, 11u },
    };
    return ic_util::make_graph(connectivity, atomVector);
  }

  coords::Representation_3D createFirstResidue() {
    return coords::Representation_3D{ coords::r3{ -6.053, -0.324, -0.108 },
      coords::r3{ -4.677, -0.093, -0.024 },
      coords::r3{ -6.262, -1.158, -0.813 },
      coords::r3{ -6.582, 0.600, -0.424 },
      coords::r3{ -6.431, -0.613, 0.894 },
      coords::r3{ -4.387, 0.166, -0.937 }
      } / energy::bohr2ang;
  }

  coords::Representation_3D createSecondResidue() {
    return coords::Representation_3D{ coords::r3{ -6.146, 3.587, -0.024 },
      coords::r3{ -4.755, 3.671, -0.133 },
      coords::r3{ -6.427, 2.922, 0.821 },
      coords::r3{ -6.587, 3.223, -0.978 },
      coords::r3{ -6.552, 4.599, 0.179 },
      coords::r3{ -4.441, 2.753, -0.339 }
    } / energy::bohr2ang;
  }

  std::vector<std::size_t> createFirstResidueIndices() {
    return { 1u, 2u, 3u, 4u, 5u, 6u };
  }

  std::vector<std::size_t> createSecondResidueIndices() {
    return { 7u, 8u, 9u, 10u, 11u, 12u };
  }

  std::vector<InternalCoordinates::BondDistance> expectedBondsForTwoMethanol() {
    return { 
      InternalCoordinates::BondDistance{ ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 2,"O","O" } },
      InternalCoordinates::BondDistance{ ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 3,"H","H" } },
      InternalCoordinates::BondDistance{ ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 4,"H","H" } },
      InternalCoordinates::BondDistance{ ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 5,"H","H" } },
      InternalCoordinates::BondDistance{ ic_util::Node{ 2, "O", "O" }, ic_util::Node{ 6,"H","H" } },

      InternalCoordinates::BondDistance{ ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 8,"O","O" } },
      InternalCoordinates::BondDistance{ ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 9,"H","H" } },
      InternalCoordinates::BondDistance{ ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 10,"H","H" } },
      InternalCoordinates::BondDistance{ ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 11,"H","H" } },
      InternalCoordinates::BondDistance{ ic_util::Node{ 8, "O", "O" }, ic_util::Node{ 12,"H","H" } },
    };
  }

  std::vector<InternalCoordinates::BondAngle> expectedAnglesForTwoMethanol() {
    return {
      InternalCoordinates::BondAngle{ ic_util::Node{ 2, "O", "O" }, ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 3, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 2, "O", "O" }, ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 4, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 2, "O", "O" }, ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 5, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 3, "H", "H" }, ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 4, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 3, "H", "H" }, ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 5, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 4, "H", "H" }, ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 5, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 2, "O", "O" }, ic_util::Node{ 6, "H", "H" } },

      InternalCoordinates::BondAngle{ ic_util::Node{ 8, "O", "O" }, ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 9, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 8, "O", "O" }, ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 10, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 8, "O", "O" }, ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 11, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 9, "H", "H" }, ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 10, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 9, "H", "H" }, ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 11, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 10, "H", "H" }, ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 11, "H", "H" } },
      InternalCoordinates::BondAngle{ ic_util::Node{ 7, "C", "C" }, ic_util::Node{ 8, "O", "O" }, ic_util::Node{ 12, "H", "H" } },
    };
  }

  std::vector<InternalCoordinates::DihedralAngle> expectedDihedralsForTwoMethanol() {
    return {
      InternalCoordinates::DihedralAngle{ ic_util::Node{3}, ic_util::Node{1}, ic_util::Node{2}, ic_util::Node{6} },
      InternalCoordinates::DihedralAngle{ ic_util::Node{4}, ic_util::Node{1}, ic_util::Node{2}, ic_util::Node{6} },
      InternalCoordinates::DihedralAngle{ ic_util::Node{5}, ic_util::Node{1}, ic_util::Node{2}, ic_util::Node{6} },
      InternalCoordinates::DihedralAngle{ ic_util::Node{9}, ic_util::Node{7}, ic_util::Node{8}, ic_util::Node{12} },
      InternalCoordinates::DihedralAngle{ ic_util::Node{10}, ic_util::Node{7}, ic_util::Node{8}, ic_util::Node{12} },
      InternalCoordinates::DihedralAngle{ ic_util::Node{11}, ic_util::Node{7}, ic_util::Node{8}, ic_util::Node{12} }
    };
  }

  std::vector<std::unique_ptr<InternalCoordinates::Translations>> expectedTranslationsForTwoMethanol() {
	  std::vector<std::unique_ptr<InternalCoordinates::Translations>> translations;
	  translations.emplace_back(std::make_unique<InternalCoordinates::TranslationX>(createFirstResidueIndices()));
	  translations.emplace_back(std::make_unique<InternalCoordinates::TranslationY>(createFirstResidueIndices()));
	  translations.emplace_back(std::make_unique<InternalCoordinates::TranslationZ>(createFirstResidueIndices()));
	  translations.emplace_back(std::make_unique<InternalCoordinates::TranslationX>(createSecondResidueIndices()));
	  translations.emplace_back(std::make_unique<InternalCoordinates::TranslationY>(createSecondResidueIndices()));
	  translations.emplace_back(std::make_unique<InternalCoordinates::TranslationZ>(createSecondResidueIndices()));
    return translations;
  }

  std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> expectedRotationsForTwoMethanol() {
    InternalCoordinates::CartesiansForInternalCoordinates cartesians{createSystemOfTwoMethanolMolecules()};

	auto rotatorFirst = InternalCoordinates::Rotator::buildRotator(cartesians, createFirstResidueIndices());
	auto rotatorSecond = InternalCoordinates::Rotator::buildRotator(cartesians, createSecondResidueIndices());

	std::vector<std::unique_ptr<InternalCoordinates::InternalCoordinate>> rotations;

	auto rotationsFirst = rotatorFirst->makeRotations();

	rotations.emplace_back(std::move(rotationsFirst.rotationA));
	rotations.emplace_back(std::move(rotationsFirst.rotationB));
	rotations.emplace_back(std::move(rotationsFirst.rotationC));

	auto rotationsSecond = rotatorSecond->makeRotations();

	rotations.emplace_back(std::move(rotationsSecond.rotationA));
	rotations.emplace_back(std::move(rotationsSecond.rotationB));
	rotations.emplace_back(std::move(rotationsSecond.rotationC));

    return rotations;
  }

  double constexpr doubleNearThreshold = 1.e-10;

  inline void isCartesianPointNear(coords::r3 const& lhs, coords::r3 const& rhs) {
    EXPECT_NEAR(lhs.x(), rhs.x(), doubleNearThreshold);
    EXPECT_NEAR(lhs.y(), rhs.y(), doubleNearThreshold);
    EXPECT_NEAR(lhs.z(), rhs.z(), doubleNearThreshold);
  }

}

PrimitiveInternalSetTest::PrimitiveInternalSetTest() : cartesians{ createSystemOfTwoMethanolMolecules() }, testSystem(), molecules({ createFirstResidueIndices(), createSecondResidueIndices() }),
  systemGraph{ createTestGraph() }{}

void PrimitiveInternalSetTest::createDistances() {
  auto decorator = std::make_unique<internals::ICBondDecorator>(nullptr);
	internals::NoConstraintManager manager;
	decorator->buildCoordinates(cartesians, systemGraph, molecules, manager);
  testSystem = std::make_shared<internals::PrimitiveInternalCoordinates>(*decorator);
}

void PrimitiveInternalSetTest::createAngles() {
	auto decorator = std::make_unique<internals::ICAngleDecorator>(nullptr);
	internals::NoConstraintManager manager;
	decorator->buildCoordinates(cartesians, systemGraph, molecules, manager);
  testSystem = std::make_shared<internals::PrimitiveInternalCoordinates>(*decorator);
}

void PrimitiveInternalSetTest::createDihedrals() {
	auto decorator = std::make_unique<internals::ICDihedralDecorator>(nullptr);
	internals::NoConstraintManager manager;
  decorator->buildCoordinates(cartesians, systemGraph, molecules, manager);
  testSystem = std::make_shared<internals::PrimitiveInternalCoordinates>(*decorator);
}

void PrimitiveInternalSetTest::createTranslations() {
	auto decorator = std::make_unique<internals::ICTranslationDecorator>(nullptr);
	internals::NoConstraintManager manager;
  decorator->buildCoordinates(cartesians, systemGraph, molecules, manager);
  testSystem = std::make_shared<internals::PrimitiveInternalCoordinates>(*decorator);
}

void PrimitiveInternalSetTest::createRotations() {
	auto decorator = std::make_unique<internals::ICRotationDecorator>(nullptr);
	internals::NoConstraintManager manager;
  decorator->buildCoordinates(cartesians, systemGraph, molecules, manager);
  testSystem = std::make_shared<internals::PrimitiveInternalCoordinates>(*decorator);
}

void PrimitiveInternalSetTest::distanceCreationTest() {
	createDistances();
  auto & allBonds = testSystem->primitive_internals;
  auto expectedBonds = expectedBondsForTwoMethanol();
  for (auto i = 0u; i < allBonds.size(); ++i) {
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::BondDistance*>(allBonds.at(i).get()), expectedBonds.at(i));
  }
}

TEST_F(PrimitiveInternalSetTest, distanceCreationTest) {
  distanceCreationTest();
}

void PrimitiveInternalSetTest::bondAngleCreationTest() {
	createAngles();
  auto & allAngles = testSystem->primitive_internals;
  auto expectedAngles = expectedAnglesForTwoMethanol();
  for (auto i = 0u; i < allAngles.size(); ++i) {
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::BondAngle*>(allAngles.at(i).get()), expectedAngles.at(i));
  }
}

TEST_F(PrimitiveInternalSetTest, bondAngleCreationTest) {
  bondAngleCreationTest();
}

void PrimitiveInternalSetTest::dihedralCreationTest() {
	createDihedrals();
  auto & allDihedrals = testSystem->primitive_internals;
  auto expectedDihedrals = expectedDihedralsForTwoMethanol();
  for (auto i = 0u; i < allDihedrals.size(); ++i) {
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::DihedralAngle*>(allDihedrals.at(i).get()), expectedDihedrals.at(i));
  }
}

TEST_F(PrimitiveInternalSetTest, dihedralCreationTest) {
  dihedralCreationTest();
}

void PrimitiveInternalSetTest::tarnslationXCreationTest() {
	createTranslations();
  auto & trans = testSystem->primitive_internals;
  auto expectedTranslations = expectedTranslationsForTwoMethanol();
  for (auto i = 0u; i < expectedTranslations.size(); ++i) {
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::Translations*>(trans.at(i).get()), *expectedTranslations.at(i));
  }
}

TEST_F(PrimitiveInternalSetTest, tarnslationsCreationTest) {
  tarnslationXCreationTest();
}

void PrimitiveInternalSetTest::rotationsCreationTest() {
	createRotations();
  auto & rotations = testSystem->primitive_internals;
  auto expectedRotations = expectedRotationsForTwoMethanol();
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::RotationA*>(expectedRotations.at(0).get()), *dynamic_cast<InternalCoordinates::RotationA*>(rotations.at(0).get()));
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::RotationB*>(expectedRotations.at(1).get()), *dynamic_cast<InternalCoordinates::RotationB*>(rotations.at(1).get()));
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::RotationC*>(expectedRotations.at(2).get()), *dynamic_cast<InternalCoordinates::RotationC*>(rotations.at(2).get()));
	EXPECT_EQ(*dynamic_cast<InternalCoordinates::RotationA*>(expectedRotations.at(3).get()), *dynamic_cast<InternalCoordinates::RotationA*>(rotations.at(3).get()));
	EXPECT_EQ(*dynamic_cast<InternalCoordinates::RotationB*>(expectedRotations.at(4).get()), *dynamic_cast<InternalCoordinates::RotationB*>(rotations.at(4).get()));
	EXPECT_EQ(*dynamic_cast<InternalCoordinates::RotationC*>(expectedRotations.at(5).get()), *dynamic_cast<InternalCoordinates::RotationC*>(rotations.at(5).get()));
}

TEST_F(PrimitiveInternalSetTest, rotationsCreationTest) {
  rotationsCreationTest();
}

MatricesTest::MatricesTest() : cartesians{ createSystemOfTwoMethanolMolecules() / energy::bohr2ang }, testSystem(), molecules({ createFirstResidueIndices(), createSecondResidueIndices() }),
systemGraph{ createTestGraph() }{
	std::unique_ptr<internals::ICAbstractDecorator> decorator = std::make_unique<internals::ICRotationDecorator>(nullptr);
	decorator = std::make_unique<internals::ICTranslationDecorator>(std::move(decorator));
	decorator = std::make_unique<internals::ICDihedralDecorator>(std::move(decorator));
	decorator = std::make_unique<internals::ICAngleDecorator>(std::move(decorator));
	decorator = std::make_unique<internals::ICBondDecorator>(std::move(decorator));
	internals::NoConstraintManager manager;
	decorator->buildCoordinates(cartesians, systemGraph, molecules, manager);
	testSystem = std::make_shared<internals::PrimitiveInternalCoordinates>(*decorator);
}

void MatricesTest::bMatrixTest(){
  EXPECT_EQ(testSystem->Bmat(cartesians), exampleBmatrixForTwoMethanols());
}

TEST_F(MatricesTest, bMatrixTest) {
  bMatrixTest();
}

void MatricesTest::gMatrixTest() {
  EXPECT_EQ(testSystem->Gmat(cartesians), exampleGmatrixForTwoMethanols());
}

TEST_F(MatricesTest, gMatrixTest) {
  gMatrixTest();
}

void MatricesTest::hessianGuessTest() {
  EXPECT_EQ(testSystem->guess_hessian(cartesians), exampleGuessHessianForTwoMethanols());
}

TEST_F(MatricesTest, hessianGuessTest) {
  hessianGuessTest();
}

DelocalizedMatricesTest::DelocalizedMatricesTest() : cartesians{ createSystemOfTwoMethanolMolecules() / energy::bohr2ang }, testSystem(), molecules({ createFirstResidueIndices(), createSecondResidueIndices() }), systemGraph(createTestGraph()) {
	std::unique_ptr<internals::ICAbstractDecorator> decorator = std::make_unique<internals::ICRotationDecorator>(nullptr);
	decorator = std::make_unique<internals::ICTranslationDecorator>(std::move(decorator));
	//decorator = std::make_unique<internals::ICOutOfPlaneDecorator>(std::move(decorator));
	decorator = std::make_unique<internals::ICDihedralDecorator>(std::move(decorator));
	decorator = std::make_unique<internals::ICAngleDecorator>(std::move(decorator));
	decorator = std::make_unique<internals::ICBondDecorator>(std::move(decorator));
	internals::NoConstraintManager manager;
	decorator->buildCoordinates(cartesians, systemGraph, molecules, manager);
	testSystem = std::make_shared<internals::TRIC>(*decorator, cartesians);
}

void DelocalizedMatricesTest::delocalizedMatrixTest() {
  EXPECT_EQ(testSystem->getDelMat(), exampleDelocalizedMatrixForTwoMethanols());
}

TEST_F(DelocalizedMatricesTest, delocalizedMatrixTest) {
  delocalizedMatrixTest();
}

void DelocalizedMatricesTest::delocalizedBMatrixTest() {
  EXPECT_EQ(testSystem->Bmat(cartesians), exampleDelocalizedBMatrixForTwoMethanols());
}

TEST_F(DelocalizedMatricesTest, delocalizedBMatrixTest) {
  delocalizedBMatrixTest();
}

void DelocalizedMatricesTest::delocalizedGMatrixTest() {
  EXPECT_EQ(testSystem->Gmat(cartesians), exampleDelocalizedGMatrixForTwoMethanols());
}

TEST_F(DelocalizedMatricesTest, delocalizedGMatrixTest) {
  delocalizedGMatrixTest();
}

void DelocalizedMatricesTest::delocalizedInitialHessianTest() {
  EXPECT_EQ(testSystem->guess_hessian(cartesians), exampleDelocalizedInitialHessianForTwoMethanols());
}

TEST_F(DelocalizedMatricesTest, delocalizedInitialHessianTest) {
  delocalizedInitialHessianTest();
}

ConverterMatricesTest::ConverterMatricesTest() : cartesians{ createSystemOfTwoMethanolMolecules() / energy::bohr2ang }, testSystem{}, converter{ testSystem, cartesians } {}

void ConverterMatricesTest::calculateInternalGradsTest() {
  EXPECT_CALL(testSystem, Bmat(testing::_)).WillRepeatedly(testing::ReturnRefOfCopy(exampleDelocalizedBMatrixForTwoMethanols()));
  EXPECT_CALL(testSystem, pseudoInverseOfGmat(testing::_)).WillRepeatedly(testing::Return(inverseOfGForTheFirstStep()));

  EXPECT_EQ(internalGradientsOfTwoMethanolMolecules(), converter.calculateInternalGradients(gradientsOfTwoMethanolMolecules()));
}

TEST_F(ConverterMatricesTest, calculateInternalGradsTest) {
  calculateInternalGradsTest();
}

void ConverterMatricesTest::applyInternalChangeTest() {
  EXPECT_CALL(testSystem, transposeOfBmat(testing::_))
    .WillOnce(testing::Return(transposeOfBForTheFirstStep()))
    .WillOnce(testing::Return(transposeOfBForTheSecondStep()))
    .WillOnce(testing::Return(transposeOfBForTheThirdStep()));
  EXPECT_CALL(testSystem, pseudoInverseOfGmat(testing::_))
    .WillOnce(testing::Return(inverseOfGForTheFirstStep()))
    .WillOnce(testing::Return(inverseOfGForTheSecondStep()))
    .WillOnce(testing::Return(inverseOfGForTheThirdStep()));
  EXPECT_CALL(testSystem, calc_diff(testing::_, testing::_))
    .WillOnce(testing::Return(differncesInInternalCoordinatesAfterFirstStep()))
    .WillOnce(testing::Return(differncesInInternalCoordinatesAfterSecondStep()))
    .WillOnce(testing::Return(differncesInInternalCoordinatesAfterThirdStep()));

  auto newCartesians = converter.applyInternalChange(internalInitialStepOfTwoMethanolMolecules());

  auto expectedChangeAfterFirstStep = cartesianChangeOfTwoMethanolMoleculesAfterFirstStep();

  for (auto i = 0u; i < expectedChangeAfterFirstStep.size(); ++i) {
    isCartesianPointNear(newCartesians.at(i), expectedChangeAfterFirstStep.at(i));
  }
}

TEST_F(ConverterMatricesTest, applyInternalChangeTest) {
  applyInternalChangeTest();
}

void MatricesTest::calculatePrimitiveInternalValuesTest() {
  EXPECT_EQ(testSystem->calc(cartesianChangeOfTwoMethanolMoleculesAfterFirstStep()), expectedPrimitiveValuesForTwoMethanol());
}

TEST_F(MatricesTest, calculatePrimitiveInternalValuesTest) {
  calculatePrimitiveInternalValuesTest();
}

void DelocalizedMatricesTest::internalDifferencesTest() {
  EXPECT_EQ(testSystem->calc_diff(cartesianChangeOfTwoMethanolMoleculesAfterFirstStep(),createSystemOfTwoMethanolMolecules() / energy::bohr2ang), expectedInternalChangeBetweenInitialAndFirstStepStructure());
}

TEST_F(DelocalizedMatricesTest, internalDifferencesTest) {
  internalDifferencesTest();
}

void DelocalizedMatricesTest::internalValuesForTricTest() {
  EXPECT_EQ(tricValuesAfterTheFirstStep(), testSystem->calc(cartesianChangeOfTwoMethanolMoleculesAfterFirstStep()));
}

TEST_F(DelocalizedMatricesTest, internalValuesForTricTest) {
  internalValuesForTricTest();
}
#endif
