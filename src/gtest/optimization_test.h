#ifdef GOOGLE_MOCK

#ifndef OPTIMIZATION_TEST
#define OPTIMIZATION_TEST

#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include "../PrimitiveInternalCoordinates.h"
#include "primitive_internals_test.h"

class InternalToCartesianConverterMock : public internals::InternalToCartesianConverter {
public:
  InternalToCartesianConverterMock(internals::PrimitiveInternalCoordinates & internals,
    InternalCoordinates::CartesiansForInternalCoordinates & cartesians) : InternalToCartesianConverter{ internals, cartesians } {}
  MOCK_METHOD2(getInternalStep, scon::mathmatrix<coords::float_type>(scon::mathmatrix<coords::float_type> const&, scon::mathmatrix<coords::float_type> const&));
  MOCK_METHOD3(getDeltaYPrimeAndSol, std::pair<coords::float_type, coords::float_type>(scon::mathmatrix<coords::float_type> const& internalStep, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian));
  MOCK_METHOD1(applyInternalChange, void(scon::mathmatrix<coords::float_type>));
  MOCK_CONST_METHOD0(getCartesianCoordinates, InternalCoordinates::CartesiansForInternalCoordinates const&());
  MOCK_METHOD0(getCartesianCoordinates, InternalCoordinates::CartesiansForInternalCoordinates &());
};

class OptimizerTest : public testing::Test {
public:
  OptimizerTest();
  void restrictStepTest();
  void restrictCartesianStepTest();
  scon::mathmatrix<double> gradients;
  scon::mathmatrix<double> hessian;
  InternalCoordinates::CartesiansForInternalCoordinates cartesians;
  MockPrimitiveInternals testSystem;
  internals::InternalToCartesianConverter converter;
  internals::StepRestrictor restrictor;
  InternalToCartesianConverterMock conv;
  internals::InternalToCartesianStep toCartesianNorm;
};

#endif
#endif