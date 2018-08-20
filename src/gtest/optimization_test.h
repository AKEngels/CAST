#ifdef GOOGLE_MOCK

#ifndef OPTIMIZATION_TEST
#define OPTIMIZATION_TEST

#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include "../PrimitiveInternalCoordinates.h"
#include "primitive_internals_test.h"

class OptimizerTest : public testing::Test {
public:
  OptimizerTest();
  void restrictStepTest();
  void restrictCartesianStepTest();
  InternalCoordinates::CartesiansForInternalCoordinates cartesians;
  MockPrimitiveInternals testSystem;
  internals::InternalToCartesianConverter converter;
  internals::StepRestrictor restrictor;
  internals::InternalToCartesianStep toCartesianNorm;
};

#endif
#endif