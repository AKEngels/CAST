#ifdef GOOGLE_MOCK

#ifndef OPTIMIZATION_TEST
#define OPTIMIZATION_TEST

#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include "../PrimitiveInternalCoordinates.h"

class OptimizerTest : public testing::Test {
public:
  OptimizerTest();
  void restrictStepTest();
  internals::InternalToCartesianConverter converter;
};

#endif
#endif