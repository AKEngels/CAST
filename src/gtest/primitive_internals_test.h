#ifdef GOOGLE_MOCK

#ifndef PRIMITIVE_INTERNALS_TEST_H
#define PRIMITIVE_INTERNALS_TEST_H

#include<gtest/gtest.h>
#include<gmock/gmock.h>
#include<vector>
#include<utility>

#include"../graph.h"
#include"../coords.h"
#include"../PrimitiveInternalCoordinates.h"
#include"../TranslationRotationInternalCoordinates.h"
#include"../InternalCoordinates.h"
#include"TestFiles/ExpectedValuesForInternalCoordinatesTest.h"

class MockPrimitiveInternals : public internals::PrimitiveInternalCoordinates {
public:
  MockPrimitiveInternals(){}
  MOCK_METHOD1(calc, scon::mathmatrix<coords::float_type>(coords::Representation_3D const&  xyz));
  MOCK_METHOD2(calc_diff, scon::mathmatrix<coords::float_type>(coords::Representation_3D const&  lhs, coords::Representation_3D const&  rhs));
  MOCK_CONST_METHOD1(guess_hessian, scon::mathmatrix<coords::float_type>(CartesianType const&));
  MOCK_METHOD1(Bmat, scon::mathmatrix<coords::float_type>&(CartesianType const& cartesians));
  MOCK_METHOD1(Gmat, scon::mathmatrix<coords::float_type>&(CartesianType const& cartesians));
};

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

class ConverterMatricesTestTest : public testing::Test {
public:
  ConverterMatricesTestTest();
  void calculateInternalGradsTest(); 
  void getInternalStepTest();
  InternalCoordinates::CartesiansForInternalCoordinates cartesians;
  MockPrimitiveInternals testSystem;
  internals::InternalToCartesianConverter converter;
};

#endif

#endif