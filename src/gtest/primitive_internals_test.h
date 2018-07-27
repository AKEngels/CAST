#ifdef GOOGLE_MOCK

#ifndef PRIMITIVE_INTERNALS_TEST_H
#define PRIMITIVE_INTERNALS_TEST_H

#include<gtest/gtest.h>
#include<vector>
#include<utility>

#include"../graph.h"
#include"../coords.h"
#include"../ic_core.h"
#include"TestFiles/ExpectedValuesForInternalCoordinatesTest.h"

class PrimitiveInternalSetTest : public testing::Test {
public:
  PrimitiveInternalSetTest();
  ic_core::system testSystem;
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
  ic_core::system testSystem;

  void bMatrixTest();
  void gMatrixTest();
  void hessianGuessTest();
};

class DelocalizedMatricesTest : public MatricesTest {
public:
  DelocalizedMatricesTest();
  void delocalizedMatrixTest();
  void delocalizedBMatrixTest();
  void delocalizedGMatrixTest();
  void delocalizedInitialHessianTest();
  void calculateInternalGradsTest();
  void getInternalStepTest(); 
  void applyInternalChangeTest();
  void calculatePrimitiveInternalValuesTest();
  void internalDifferencesTest();
  void internalValuesForTricTest();
};

#endif

#endif