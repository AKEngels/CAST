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

class StepRestrictorMock : public internals::StepRestrictor {
public:
  StepRestrictorMock(scon::mathmatrix<coords::float_type> * const step, coords::float_type const target) : StepRestrictor{ step, target } {}

  MOCK_METHOD1(execute, coords::float_type(internals::AppropriateStepFinder &));
  coords::float_type operator()(internals::AppropriateStepFinder & finder) override {
    return execute(finder);
  }

  MOCK_CONST_METHOD0(getRestrictedStep, scon::mathmatrix<coords::float_type> const&());
  MOCK_METHOD0(getRestrictedStep, scon::mathmatrix<coords::float_type> &());
};

class InternalToCartesianStepMock : public internals::InternalToCartesianStep {
public:
  InternalToCartesianStepMock(internals::AppropriateStepFinder & finder, coords::float_type const trustRadius)
    : InternalToCartesianStep(finder, trustRadius) {}

  MOCK_METHOD1(execute, coords::float_type(internals::StepRestrictor &));
  coords::float_type operator()(internals::StepRestrictor & restrictor) override {
    return execute(restrictor);
  }

  /*MOCK_METHOD1(execute, std::pair<coords::float_type, internals::StepRestrictor>(coords::float_type const));
  std::pair<coords::float_type, internals::StepRestrictor> operator()(coords::float_type const target) override {
    return execute(target);
  }*/
};

class AppropriateStepFinderMock : public internals::AppropriateStepFinder {
public:
  AppropriateStepFinderMock(internals::InternalToCartesianConverter const& converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian);

  MOCK_CONST_METHOD1(getDeltaYPrime, coords::float_type(scon::mathmatrix<coords::float_type> const&));
  MOCK_CONST_METHOD1(getSol, coords::float_type(scon::mathmatrix<coords::float_type> const&));
  MOCK_CONST_METHOD0(getInternalStep, scon::mathmatrix<coords::float_type>());
  MOCK_CONST_METHOD1(getInternalStep, scon::mathmatrix<coords::float_type>(scon::mathmatrix<coords::float_type> const&));
  MOCK_METHOD1(applyInternalChangeAndGetNorm, coords::float_type(scon::mathmatrix<coords::float_type> const&));
  MOCK_CONST_METHOD1(alterHessian, scon::mathmatrix<coords::float_type>(coords::float_type const));
};

class FinderImplementation {
private:
  MockPrimitiveInternals internals;
  InternalCoordinates::CartesiansForInternalCoordinates cartesians;
  InternalToCartesianConverterMock converter{ internals, cartesians };
  scon::mathmatrix<coords::float_type> emptyGradients;
  scon::mathmatrix<coords::float_type> emptyHessian;
public:
  FinderImplementation();
  AppropriateStepFinderMock finder;
};

class StepRestrictorTest : public testing::Test, public FinderImplementation {
public:
  StepRestrictorTest();
  scon::mathmatrix<coords::float_type> expectedStep;
  internals::StepRestrictor  restrictor;
};

class InternalToCartesianStepTest : public testing::Test, public FinderImplementation {
public:
  InternalToCartesianStepTest();
  scon::mathmatrix<coords::float_type> fakeStep;
  internals::InternalToCartesianStep toCartesianNorm;
};

class BrentsMethodTest : public testing::Test, public FinderImplementation {
public:
  BrentsMethodTest();
  internals::BrentsMethod brent;
};

#endif
#endif