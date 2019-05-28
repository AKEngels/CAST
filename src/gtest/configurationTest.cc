/**
CAST 3
configurationTest
Purpose: Tests functions in configuration.cc

@author Dustin Kaiser
@version 1.0
*/


#ifdef GOOGLE_MOCK

#include "../configuration.h"
#include "../configurationHelperfunctions.h"
#include <gtest/gtest.h>

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


TEST(doubles_from_string, withProperInput)
{
	std::string input = "1.0, 2.2, 3.5, 5.2";
	std::vector<double> sorted = config::doubles_from_string(input);
	std::vector<double> compare = { 1.0, 2.2, 3.5, 5.2 };
	ASSERT_EQ(sorted, compare);
}

TEST(doubles_from_string, IntegersAreAlsoOk)
{
	std::string input = "1, 2.5, 3";
	std::vector<double> sorted = config::doubles_from_string(input);
	std::vector<double> compare = { 1, 2.5, 3 };
	ASSERT_EQ(sorted, compare);
}

TEST(doubles_from_string, identicalNumbersDontMatter)
{
	std::string input = "1.1, 3.5, 5.5, 5.5, 1.1";
	std::vector<double> sorted = config::doubles_from_string(input);
	std::vector<double> compare = { 1.1, 3.5, 5.5, 5.5, 1.1 };
	ASSERT_EQ(sorted, compare);
}

TEST(doubles_from_string, NegativeNumbersOk)
{
	std::string input = "1.1, 3.3, -5.5";
	std::vector<double> sorted = config::doubles_from_string(input);
	std::vector<double> compare = { 1.1, 3.3, -5.5 };
	ASSERT_EQ(sorted, compare);

	std::string input_2 = "1, -10.5";
	std::vector<double> sorted_2 = config::doubles_from_string(input_2);
	std::vector<double> compare_2 = { 1, -10.5 };
	ASSERT_EQ(sorted_2, compare_2);
}

TEST(doubles_from_string, rangeThrowsError)
{
	std::string input = "0 - 1, 3-7, 4";
	ASSERT_THROW(config::doubles_from_string(input), std::runtime_error);
}

TEST(doubles_from_string, emtpyFieldThrowsError)
{
	std::string input = "0.1, 1, 3.7,, 4";
	ASSERT_THROW(config::doubles_from_string(input), std::runtime_error);

	input = "0.1, 1, 3.7, , 4";
	ASSERT_THROW(config::doubles_from_string(input), std::runtime_error);

	input = "0.1, 1, 3.7, 4,";
	ASSERT_THROW(config::doubles_from_string(input), std::runtime_error);
}

#endif
