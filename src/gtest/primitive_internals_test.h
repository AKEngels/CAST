#ifdef GOOGLE_MOCK

#ifndef PRIMITIVE_INTERNALS_TEST_H
#define PRIMITIVE_INTERNALS_TEST_H

#include<gtest/gtest.h>
#include<vector>
#include<utility>

#include"../graph.h"
#include"../coords.h"
#include"../ic_core.h"

class PrimitiveInternalSetTest : public testing::Test {
public:
  PrimitiveInternalSetTest();
  ic_core::system testSystem;
  ic_util::Graph<ic_util::Node> systemGraph;
  void distanceCreationTest();
  void bondAngleCreationTest();
};



#endif

#endif