/**
CAST 3
configurationTest
Purpose: Tests functions in configuration.cc

@author Dustin Kaiser
@version 1.0
*/


#ifdef GOOGLE_MOCK
#pragma once
#include "../src/configuration.h"
#include "../src/configurationHelperfunctions.h"
#include "gtest/gtest.h"

TEST(sorted_indices_from_cs_string, withProperInput)
{
  std::string input = "1, 2, 3, 5";
  std::vector<std::size_t> sorted = config::sorted_indices_from_cs_string(input);
  std::vector<std::size_t> compare = { 1u, 2u, 3u, 5u };
  ASSERT_EQ(sorted, compare);
}

TEST(sorted_indices_from_cs_string, FailsWithFloatingPointInput)
{
  std::string input = "1, 2.5, 3";
  ASSERT_THROW(config::sorted_indices_from_cs_string(input), std::runtime_error);
}

TEST(sorted_indices_from_cs_string, onlyCapturesIdenticalNumbersOnce)
{
  std::string input = "1, 3, 5, 5, 1,";
  std::vector<std::size_t> sorted = config::sorted_indices_from_cs_string(input);
  std::vector<std::size_t> compare = { 1u, 3u, 5u };
  ASSERT_EQ(sorted, compare);
}

TEST(sorted_indices_from_cs_string, throwsWhenNegativeNumbersArePresent)
{
  std::string input = "1, 3, -5,";
  ASSERT_THROW(config::sorted_indices_from_cs_string(input), std::runtime_error);
  std::string input2 = "1, -10.5";
  ASSERT_THROW(config::sorted_indices_from_cs_string(input2), std::runtime_error);
}
TEST(sorted_indices_from_cs_string, takesRangeAndIsUnique)
{
  std::string input = "0 - 1, 3-7, 4,";
  std::vector<std::size_t> sorted = config::sorted_indices_from_cs_string(input);
  std::vector<std::size_t> compare = { 0u, 1u, 3u, 4u, 5u, 6u, 7u };
  ASSERT_EQ(sorted, compare);
}

#endif
