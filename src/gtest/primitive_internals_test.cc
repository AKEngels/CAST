#ifdef GOOGLE_MOCK

#include "primitive_internals_test.h"
#include"../InternalCoordinates.h"

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
    return { 0u, 1u, 2u, 3u, 4u, 5u };
  }

  std::vector<std::size_t> createSecondResidueIndices() {
    return { 6u, 7u, 8u, 9u, 10u, 11u };
  }

  std::vector<InternalCoordinates::BondDistance> expectedBondsForTwoMethanol() {
    return { 
      InternalCoordinates::BondDistance{ 1, 2, "C", "O" },
      InternalCoordinates::BondDistance{ 1, 3, "C", "H" },
      InternalCoordinates::BondDistance{ 1, 4, "C", "H" },
      InternalCoordinates::BondDistance{ 1, 5, "C", "H" },
      InternalCoordinates::BondDistance{ 2, 6, "O", "H" },

      InternalCoordinates::BondDistance{ 7, 8, "C", "O" },
      InternalCoordinates::BondDistance{ 7, 9, "C", "H" },
      InternalCoordinates::BondDistance{ 7, 10, "C", "H" },
      InternalCoordinates::BondDistance{ 7, 11, "C", "H" },
      InternalCoordinates::BondDistance{ 8, 12, "O", "H" }
    };
  }

  std::vector<InternalCoordinates::BondAngle> expectedAnglesForTwoMethanol() {
    return {
      InternalCoordinates::BondAngle{ 2, 1, 3, "O", "C", "H" },
      InternalCoordinates::BondAngle{ 2, 1, 4, "O", "C", "H" },
      InternalCoordinates::BondAngle{ 2, 1, 5, "O", "C", "H" },
      InternalCoordinates::BondAngle{ 3, 1, 4, "H", "C", "H" },
      InternalCoordinates::BondAngle{ 3, 1, 5, "H", "C", "H" },
      InternalCoordinates::BondAngle{ 4, 1, 5, "H", "C", "H" },
      InternalCoordinates::BondAngle{ 1, 2, 6, "C", "O", "H" },

      InternalCoordinates::BondAngle{ 8, 7, 9, "O", "C", "H" },
      InternalCoordinates::BondAngle{ 8, 7, 10, "O", "C", "H" },
      InternalCoordinates::BondAngle{ 8, 7, 11, "O", "C", "H" },
      InternalCoordinates::BondAngle{ 9, 7, 10, "H", "C", "H" },
      InternalCoordinates::BondAngle{ 9, 7, 11, "H", "C", "H" },
      InternalCoordinates::BondAngle{ 10, 7, 11, "H", "C", "H" },
      InternalCoordinates::BondAngle{ 7, 8, 12, "C", "O", "H" }
    };
  }

  std::vector<InternalCoordinates::DihedralAngle> expectedDihedralsForTwoMethanol() {
    return {
      InternalCoordinates::DihedralAngle{ 3, 1, 2, 6 },
      InternalCoordinates::DihedralAngle{ 4, 1, 2, 6 },
      InternalCoordinates::DihedralAngle{ 5, 1, 2, 6 },
      InternalCoordinates::DihedralAngle{ 9, 7, 8, 12 },
      InternalCoordinates::DihedralAngle{ 10, 7, 8, 12 },
      InternalCoordinates::DihedralAngle{ 11, 7, 8, 12 }
    };
  }

  std::vector<InternalCoordinates::TranslationX> expectedTranslationsForTwoMethanol() {
    return {
      InternalCoordinates::TranslationX{ createFirstResidueIndices() },
      InternalCoordinates::TranslationX{ createSecondResidueIndices() }
    };
  }

  std::vector<std::shared_ptr<InternalCoordinates::Rotator>> expectedRotationsForTwoMethanol() {
    InternalCoordinates::CartesiansForInternalCoordinates cartesians{createSystemOfTwoMethanolMolecules()};
    return { InternalCoordinates::Rotator::buildRotator(cartesians, createFirstResidueIndices()),
      InternalCoordinates::Rotator::buildRotator(cartesians, createSecondResidueIndices()) };
  }

  double constexpr doubleNearThreshold = 1.e-10;

  inline void isCartesianPointNear(coords::r3 const& lhs, coords::r3 const& rhs) {
    EXPECT_NEAR(lhs.x(), rhs.x(), doubleNearThreshold);
    EXPECT_NEAR(lhs.y(), rhs.y(), doubleNearThreshold);
    EXPECT_NEAR(lhs.z(), rhs.z(), doubleNearThreshold);
  }

}

PrimitiveInternalSetTest::PrimitiveInternalSetTest() : testSystem({ createFirstResidue(), createSecondResidue() }, { createFirstResidueIndices(), createSecondResidueIndices() }, createSystemOfTwoMethanolMolecules()),
  systemGraph{ createTestGraph() }{}

void PrimitiveInternalSetTest::distanceCreationTest() {
  auto allBonds = testSystem.create_distances(systemGraph);
  auto expectedBonds = expectedBondsForTwoMethanol();
  for (auto i = 0u; i < allBonds.size(); ++i) {
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::BondDistance*>(allBonds.at(i).get()), expectedBonds.at(i));
  }
}

TEST_F(PrimitiveInternalSetTest, distanceCreationTest) {
  distanceCreationTest();
}

void PrimitiveInternalSetTest::bondAngleCreationTest() {
  auto allAngles = testSystem.create_angles(systemGraph);
  auto expectedAngles = expectedAnglesForTwoMethanol();
  for (auto i = 0u; i < allAngles.size(); ++i) {
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::BondAngle*>(allAngles.at(i).get()), expectedAngles.at(i));
  }
}

TEST_F(PrimitiveInternalSetTest, bondAngleCreationTest) {
  bondAngleCreationTest();
}

void PrimitiveInternalSetTest::dihedralCreationTest() {
  auto allDihedrals = testSystem.create_dihedrals(systemGraph);
  auto expectedDihedrals = expectedDihedralsForTwoMethanol();
  for (auto i = 0u; i < allDihedrals.size(); ++i) {
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::DihedralAngle*>(allDihedrals.at(i).get()), expectedDihedrals.at(i));
  }
}

TEST_F(PrimitiveInternalSetTest, dihedralCreationTest) {
  dihedralCreationTest();
}

void PrimitiveInternalSetTest::tarnslationXCreationTest() {
  auto transX = testSystem.create_trans_x({ createFirstResidueIndices(), createSecondResidueIndices() });
  auto expectedTorsions = expectedTranslationsForTwoMethanol();
  for (auto i = 0u; i < expectedTorsions.size(); ++i) {
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::Translations*>(transX.at(i).get()), expectedTorsions.at(i));
  }
}

TEST_F(PrimitiveInternalSetTest, tarnslationXCreationTest) {
  tarnslationXCreationTest();
}

void PrimitiveInternalSetTest::tarnslationYCreationTest() {
  auto transY = testSystem.create_trans_y({ createFirstResidueIndices(), createSecondResidueIndices() });
  auto expectedTorsions = expectedTranslationsForTwoMethanol();
  for (auto i = 0u; i < expectedTorsions.size(); ++i) {
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::Translations*>(transY.at(i).get()), expectedTorsions.at(i));
  }
}

TEST_F(PrimitiveInternalSetTest, tarnslationYCreationTest) {
  tarnslationYCreationTest();
}

void PrimitiveInternalSetTest::tarnslationZCreationTest() {
  auto transZ = testSystem.create_trans_z({ createFirstResidueIndices(), createSecondResidueIndices() });
  auto expectedTorsions = expectedTranslationsForTwoMethanol();
  for (auto i = 0u; i < expectedTorsions.size(); ++i) {
    EXPECT_EQ(*dynamic_cast<InternalCoordinates::Translations*>(transZ.at(i).get()), expectedTorsions.at(i));
  }
}

TEST_F(PrimitiveInternalSetTest, tarnslationZCreationTest) {
  tarnslationYCreationTest();
}

void PrimitiveInternalSetTest::rotationsCreationTest() {
  auto translations = testSystem.create_rotations(createSystemOfTwoMethanolMolecules(), { createFirstResidueIndices(), createSecondResidueIndices() });
  auto expectedRotations = expectedRotationsForTwoMethanol();
  for (auto i = 0u; i < expectedRotations.size(); ++i) {
    EXPECT_EQ(*translations.at(i).get(), *expectedRotations.at(i).get());
  }
}

TEST_F(PrimitiveInternalSetTest, rotationsCreationTest) {
  rotationsCreationTest();
}

MatricesTest::MatricesTest() : testSystem({ createFirstResidue(), createSecondResidue() }, { createFirstResidueIndices(), createSecondResidueIndices() }, createSystemOfTwoMethanolMolecules()/energy::bohr2ang) {
  ic_util::Graph<ic_util::Node> graph{ createTestGraph() };
  testSystem.create_ic_system(graph);
}

void MatricesTest::bMatrixTest(){
  std::istringstream iss(exampleBmatrixForTwoMethanols());

  auto constexpr rowsOfBmatrix = 42u;
  auto constexpr colsOfBmatrix = 36u;

  auto expectedValuesForTheBmatrix = ReadMatrixFiles(iss).readNLinesOfFileWithMNumbers(rowsOfBmatrix, colsOfBmatrix);

  EXPECT_EQ(testSystem.Bmat(), expectedValuesForTheBmatrix);
}

TEST_F(MatricesTest, bMatrixTest) {
  bMatrixTest();
}

void MatricesTest::gMatrixTest() {
  std::istringstream iss(exampleGmatrixForTwoMethanols());

  auto constexpr rowsOfGmatrix = 42u;
  auto constexpr colsOfGmatrix = 42u;

  auto expectedValuesForTheGmatrix = ReadMatrixFiles(iss).readNLinesOfFileWithMNumbers(rowsOfGmatrix, colsOfGmatrix);

  EXPECT_EQ(testSystem.Gmat(), expectedValuesForTheGmatrix);
}

TEST_F(MatricesTest, gMatrixTest) {
  gMatrixTest();
}

void MatricesTest::hessianGuessTest() {
  std::istringstream iss(exampleGuessHessianForTwoMethanols());

  auto constexpr rowsOfHessian = 42u;
  auto constexpr colsOfHessian = 42u;

  auto expectedValuesForHessian = ReadMatrixFiles(iss).readNLinesOfFileWithMNumbers(rowsOfHessian, colsOfHessian);

  EXPECT_EQ(testSystem.guess_hessian(), expectedValuesForHessian);
}

TEST_F(MatricesTest, hessianGuessTest) {
  hessianGuessTest();
}

DelocalizedMatricesTest::DelocalizedMatricesTest() : MatricesTest{} {
  testSystem.delocalize_ic_system();
  testSystem.initial_delocalized_hessian();
}

void DelocalizedMatricesTest::delocalizedMatrixTest() {
  std::istringstream iss(exampleDelocalizedMatrixForTwoMethanols());

  auto constexpr rowsOfDelocalizedMatrix = 42u;
  auto constexpr colsOfDelocalizedMatrix = 36u;

  auto expectedValuesForDelocalizedMatrix = ReadMatrixFiles(iss).readNLinesOfFileWithMNumbers(rowsOfDelocalizedMatrix, colsOfDelocalizedMatrix);

  EXPECT_EQ(testSystem.getDelMat(), expectedValuesForDelocalizedMatrix);
}

TEST_F(DelocalizedMatricesTest, delocalizedMatrixTest) {
  delocalizedMatrixTest();
}

void DelocalizedMatricesTest::delocalizedBMatrixTest() {
  std::istringstream iss(exampleDelocalizedBMatrixForTwoMethanols());

  auto constexpr rowsOfDelocalizedBMatrix = 36u;
  auto constexpr colsOfDelocalizedBMatrix = 36u;

  auto expectedValuesForDelocalizedBMatrix = ReadMatrixFiles(iss).readNLinesOfFileWithMNumbers(rowsOfDelocalizedBMatrix, colsOfDelocalizedBMatrix);

  EXPECT_EQ(testSystem.ic_Bmat(), expectedValuesForDelocalizedBMatrix);
}

TEST_F(DelocalizedMatricesTest, delocalizedBMatrixTest) {
  delocalizedBMatrixTest();
}

void DelocalizedMatricesTest::delocalizedGMatrixTest() {
  std::istringstream iss(exampleDelocalizedGMatrixForTwoMethanols());

  auto constexpr rowsOfDelocalizedGMatrix = 36u;
  auto constexpr colsOfDelocalizedGMatrix = 36u;

  auto expectedValuesForDelocalizedGMatrix = ReadMatrixFiles(iss).readNLinesOfFileWithMNumbers(rowsOfDelocalizedGMatrix, colsOfDelocalizedGMatrix);

  EXPECT_EQ(testSystem.ic_Gmat(), expectedValuesForDelocalizedGMatrix);
}

TEST_F(DelocalizedMatricesTest, delocalizedGMatrixTest) {
  delocalizedGMatrixTest();
}

void DelocalizedMatricesTest::delocalizedInitialHessianTest() {
  std::istringstream iss(exampleDelocalizedInitialHessianForTwoMethanols());

  auto constexpr rowsOfDelocalizedInitialHessian = 36u;
  auto constexpr colsOfDelocalizedInitialHessian = 36u;

  auto expectedValuesForDelocalizedInitialHessian = ReadMatrixFiles(iss).readNLinesOfFileWithMNumbers(rowsOfDelocalizedInitialHessian, colsOfDelocalizedInitialHessian);

  EXPECT_EQ(testSystem.initial_delocalized_hessian(), expectedValuesForDelocalizedInitialHessian);
}

TEST_F(DelocalizedMatricesTest, delocalizedInitialHessianTest) {
  delocalizedInitialHessianTest();
}

void DelocalizedMatricesTest::calculateInternalGradsTest() {
  EXPECT_EQ(internalGradientsOfTwoMethanolMolecules(), testSystem.calculate_internal_grads(gradientsOfTwoMethanolMolecules()));
}

TEST_F(DelocalizedMatricesTest, calculateInternalGradsTest) {
  calculateInternalGradsTest();
}

void DelocalizedMatricesTest::getInternalStepTest() {
  std::ofstream ofs("InternalStep");
  ofs << std::setprecision(15) << testSystem.get_internal_step(internalGradientsOfTwoMethanolMolecules());

  EXPECT_EQ(internalInitialStepOfTwoMethanolMolecules(), testSystem.get_internal_step(internalGradientsOfTwoMethanolMolecules()));
}

TEST_F(DelocalizedMatricesTest, getInternalStepTest) {
  getInternalStepTest();
}

void DelocalizedMatricesTest::applyInternalChangeTest() {
  testSystem.apply_internal_change(internalInitialStepOfTwoMethanolMolecules());

  auto expectedChangeAfterFirstStep = cartesianChangeOfTwoMethanolMoleculesAfterFirstStep();

  for (auto i = 0u; i < expectedChangeAfterFirstStep.size(); ++i) {
    isCartesianPointNear(testSystem.getXyz().at(i), expectedChangeAfterFirstStep.at(i));
  }
}

TEST_F(DelocalizedMatricesTest, applyInternalChangeTest) {
  applyInternalChangeTest();
}

void DelocalizedMatricesTest::calculatePrimitiveInternalValuesTest() {
  EXPECT_EQ(testSystem.calc_prims(cartesianChangeOfTwoMethanolMoleculesAfterFirstStep()), expectedPrimitiveValuesForTwoMethanol());
}

TEST_F(DelocalizedMatricesTest, calculatePrimitiveInternalValuesTest) {
  calculatePrimitiveInternalValuesTest();
}

void DelocalizedMatricesTest::internalDifferencesTest() {
  EXPECT_EQ(testSystem.calc_diff(cartesianChangeOfTwoMethanolMoleculesAfterFirstStep(),createSystemOfTwoMethanolMolecules() / energy::bohr2ang), expectedInternalChangeBetweenInitialAndFirstStepStructure());
}

TEST_F(DelocalizedMatricesTest, internalDifferencesTest) {
  internalDifferencesTest();
}

#endif
