#ifdef GOOGLE_MOCK

#ifndef PRIMITIVE_INTERNALS_TEST_H
#define PRIMITIVE_INTERNALS_TEST_H

#include<gtest/gtest.h>
#include<vector>
#include<utility>

#include"../graph.h"
#include"../coords.h"
#include"../PrimitiveInternalCoordinates.h"
#include"../TranslationRotationInternalCoordinates.h"
#include"../InternalCoordinates.h"
#include"TestFiles/ExpectedValuesForInternalCoordinatesTest.h"

class PrimitiveInternalSetTest : public testing::Test {
public:
  PrimitiveInternalSetTest();
  InternalCoordinates::CartesiansForInternalCoordinates cartesians;
  internals::PrimitiveInternalCoordinates testSystem;
  ic_util::Graph<ic_util::Node> systemGraph;
  void distanceCreationTest();
  void bondAngleCreationTest();
  void dihedralCreationTest();
  void tarnslationXCreationTest();
  void tarnslationYCreationTest();
  void tarnslationZCreationTest();
  void rotationsCreationTest();
};

class MatricesTest : public testing::Test {
public:
  MatricesTest();
  InternalCoordinates::CartesiansForInternalCoordinates cartesians;
  internals::PrimitiveInternalCoordinates testSystem;

  void bMatrixTest();
  void gMatrixTest();
  void hessianGuessTest();
  void calculatePrimitiveInternalValuesTest();
};

class DelocalizedMatricesTest : public testing::Test {
public:
  DelocalizedMatricesTest();
  void delocalizedMatrixTest();
  void delocalizedBMatrixTest();
  void delocalizedGMatrixTest();
  void delocalizedInitialHessianTest();
  void internalDifferencesTest();
  void internalValuesForTricTest();


  InternalCoordinates::CartesiansForInternalCoordinates cartesians;
  internals::TRIC testSystem;
};

class ConverterMatricesTest : public testing::Test {
public:
  ConverterMatricesTest();
  void calculateInternalGradsTest();
  void getInternalStepTest();
  void applyInternalChangeTest();
  InternalCoordinates::CartesiansForInternalCoordinates cartesians;
  internals::TRIC testSystem;
  internals::InternalToCartesianConverter converter;
};

#endif

#endif