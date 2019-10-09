#ifdef GOOGLE_MOCK

#ifndef OPTIMIZATION_TEST
#define OPTIMIZATION_TEST

#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include "../../InternalCoordinates/PrimitiveInternalCoordinates.h"
#include "../InternalCoordinates/primitive_internals_test.h"

class InternalToCartesianConverterMock : public internals::InternalToCartesianConverter {
public:
  InternalToCartesianConverterMock(internals::PrimitiveInternalCoordinates& internals,
    InternalCoordinates::CartesiansForInternalCoordinates& cartesians) : InternalToCartesianConverter{ internals, cartesians } {}

  MOCK_METHOD2(getInternalStep, scon::mathmatrix<coords::float_type>(scon::mathmatrix<coords::float_type> const&, scon::mathmatrix<coords::float_type> const&));
  MOCK_METHOD3(getDeltaYPrimeAndSol, std::pair<coords::float_type, coords::float_type>(scon::mathmatrix<coords::float_type> const& internalStep, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian));
  MOCK_CONST_METHOD1(applyInternalChange, coords::Representation_3D(scon::mathmatrix<coords::float_type>));
  MOCK_CONST_METHOD0(getCartesianCoordinates, InternalCoordinates::CartesiansForInternalCoordinates const& ());
  MOCK_METHOD0(getCartesianCoordinates, InternalCoordinates::CartesiansForInternalCoordinates& ());
};

class StepRestrictorMock : public internals::StepRestrictor {
public:
  StepRestrictorMock(std::shared_ptr<scon::mathmatrix<coords::float_type>> const step, std::shared_ptr<coords::Representation_3D> cartesians, coords::float_type const target) : StepRestrictor{ step, cartesians, target } {}

  MOCK_METHOD1(execute, coords::float_type(internals::AppropriateStepFinder&));
  coords::float_type operator()(internals::AppropriateStepFinder& finder) override {
    return execute(finder);
  }

  MOCK_CONST_METHOD0(getRestrictedStep, scon::mathmatrix<coords::float_type> const& ());
  MOCK_METHOD0(getRestrictedStep, scon::mathmatrix<coords::float_type>& ());
};

class InternalToCartesianStepMock : public internals::InternalToCartesianStep {
public:
  InternalToCartesianStepMock(internals::AppropriateStepFinder& finder, coords::float_type const trustRadius)
    : InternalToCartesianStep(finder, trustRadius) {}

  MOCK_METHOD1(execute, coords::float_type(internals::StepRestrictor&));
  coords::float_type operator()(internals::StepRestrictor& restrictor) override {
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
  MOCK_METHOD1(applyInternalChangeAndGetNorm, coords::float_type(internals::StepRestrictor&));
  MOCK_CONST_METHOD1(alterHessian, scon::mathmatrix<coords::float_type>(coords::float_type const));
};

class FinderImplementation {
private:
  MockPrimitiveInternals internals;
  InternalCoordinates::CartesiansForInternalCoordinates cartesians;
  InternalToCartesianConverterMock converter;
  scon::mathmatrix<coords::float_type> emptyGradients;
  scon::mathmatrix<coords::float_type> emptyHessian;
public:
  FinderImplementation();
  AppropriateStepFinderMock finder;
};

class AppropriateStepFinderAvoidingInverseOfHessian : public internals::AppropriateStepFinder {
public:
  AppropriateStepFinderAvoidingInverseOfHessian(internals::InternalToCartesianConverter const& converter, scon::mathmatrix<coords::float_type> const& gradients, scon::mathmatrix<coords::float_type> const& hessian);
};

class StepRestrictorTest : public testing::Test, public FinderImplementation {
public:
  StepRestrictorTest();
  std::shared_ptr<scon::mathmatrix<coords::float_type>> expectedStep;
  std::shared_ptr<coords::Representation_3D> expectedCartesians;
  internals::StepRestrictor  restrictor;
};

class InternalToCartesianStepTest : public testing::Test, public FinderImplementation {
public:
  InternalToCartesianStepTest();
  std::shared_ptr<scon::mathmatrix<coords::float_type>> fakeStep;
  std::shared_ptr<coords::Representation_3D> fakeCartesians;
  internals::InternalToCartesianStep toCartesianNorm;
};

class BrentsMethodTest : public testing::Test, public FinderImplementation {
public:
  BrentsMethodTest();
  internals::BrentsMethod brent;
};

class AppropriateStepFinderTest : public testing::Test {
public:
  AppropriateStepFinderTest();
  MockPrimitiveInternals internals;
  InternalCoordinates::CartesiansForInternalCoordinates cartesians;
  InternalToCartesianConverterMock converter;
  scon::mathmatrix<coords::float_type> gradients;
  scon::mathmatrix<coords::float_type> hessian;
  AppropriateStepFinderAvoidingInverseOfHessian finder;
};

#endif
#endif