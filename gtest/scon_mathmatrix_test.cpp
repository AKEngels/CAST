/**
CAST 3
scon_mathmatrix_tests
Purpose: Tests matrix procedures

@author Dustin Kaiser
@version 1.0
*/


#ifdef GOOGLE_MOCK
#include "gtest/gtest.h"
#pragma once

#include "../src/scon_mathmatrix.h"
#include <iostream>
#include <string>
#include <sstream>

TEST(mathmatrix, constructedEmpty)
{
  mathmatrix<float> one;
  ASSERT_EQ(one.rows(), 0);
  ASSERT_EQ(one.cols(), 0);
}

TEST(mathmatrix, constructedAndFilled)
{
  mathmatrix<float> one(3u, 3u, 2.f);
  float count = 0.f;
  for (size_t i = 0u; i < 3u; i++)
    for (size_t j = 0u; j < 3u; j++)
      count += one(i, j);
  ASSERT_FLOAT_EQ(count, 18.f);
}

TEST(mathmatrix, constructedIdentityMatrixQuadratic)
{
  mathmatrix<float> one(3u, 3u);
  mathmatrix<float> two = one.identity(4u, 4u);
  ASSERT_EQ(two.rows(), 4u);
  ASSERT_EQ(two.cols(), 4u);
  ASSERT_FLOAT_EQ(two(2, 2), 1.f);
  ASSERT_FLOAT_EQ(two(1, 2), 0.f);
  ASSERT_FLOAT_EQ(two(3, 3), 1.f);
  ASSERT_FLOAT_EQ(two(3, 2), 0.f);
  ASSERT_FLOAT_EQ(two(1, 1), 1.f);
  ASSERT_FLOAT_EQ(two(0, 1), 0.f);
}

TEST(mathmatrix, constructedIdentityMatrixRectangular)
{

  mathmatrix<float> one(3u, 3u);
  mathmatrix<float> two = one.identity(2u, 4u);

  ASSERT_EQ(two.rows(), 2u);
  ASSERT_EQ(two.cols(), 4u);
  ASSERT_FLOAT_EQ(two(0, 0), 1.f);
  ASSERT_FLOAT_EQ(two(1, 1), 1.f);
  ASSERT_FLOAT_EQ(two(0, 2), 0.f);
  ASSERT_FLOAT_EQ(two(0, 3), 0.f);
  ASSERT_FLOAT_EQ(two(1, 2), 0.f);
  ASSERT_FLOAT_EQ(two(1, 3), 0.f);

  //

  mathmatrix<float> three = one.identity(4u, 2u);
  ASSERT_EQ(three.rows(), 4u);
  ASSERT_EQ(three.cols(), 2u);
  ASSERT_FLOAT_EQ(three(0, 0), 1.f);
  ASSERT_FLOAT_EQ(three(1, 1), 1.f);
  ASSERT_FLOAT_EQ(three(2, 0), 0.f);
  ASSERT_FLOAT_EQ(three(3, 0), 0.f);
  ASSERT_FLOAT_EQ(three(3, 1), 0.f);
  ASSERT_FLOAT_EQ(three(2, 1), 0.f);
}

TEST(mathmatrix, outputStreamIsAsItShouldBe)
{
  mathmatrix<float> one(3u, 3u, 3.f);
  one(2, 2) = 5.f;
  one(1, 0) = -4.;
  one(0, 0) = -3.125;

  std::stringstream stream;
  stream << one;
  std::string temp = stream.str();
  std::string thisIsCorrect = "-3.12500000e+00    3.00000000e+00     3.00000000e+00     \n-4.00000000e+00    3.00000000e+00     3.00000000e+00     \n3.00000000e+00     3.00000000e+00     5.00000000e+00     \n";
  ASSERT_EQ(temp, thisIsCorrect);
}

TEST(mathmatrix, upperLeftSubmatrixCorrectlyTruncated)
{
  mathmatrix<float> one(3u, 3u, 1.f);
  mathmatrix<float> two = one.upper_left_submatrix(2, 2);
  mathmatrix<float> three = one.upper_left_submatrix(1, 2);
  ASSERT_EQ(two.rows(), 2);
  ASSERT_EQ(two.cols(), 2);

  ASSERT_EQ(three.rows(), 1);
  ASSERT_EQ(three.cols(), 2);
}

TEST(mathmatrix, positiveDefiniteCheck)
{
  // via https://en.wikipedia.org/wiki/Positive-definite_matrix
  mathmatrix<float> one(3u, 3u, 0.f);
  one(0, 0) = 2.;
  one(1, 0) = -1.;
  one(1, 1) = 2;
  one(2, 2) = 2;
  one(0, 1) = -1;
  one(2, 1) = -1.;
  one(1, 2) = -1;

  ASSERT_EQ(one.positive_definite_check(), true);
  mathmatrix<float> two(3u, 3u, 1.f);
  ASSERT_EQ(two.positive_definite_check(), false);

}

TEST(mathmatrix, detSign)
{
  // via https://en.wikipedia.org/wiki/Positive-definite_matrix
  mathmatrix<float> one(3u, 3u, 0.f);
  one(0, 0) = 2.;
  one(1, 0) = -1.;
  one(1, 1) = 2;
  one(2, 2) = 2;
  one(0, 1) = -1;
  one(2, 1) = -1.;
  one(1, 2) = -1;
  ASSERT_EQ(one.det_sign(), 1);
  mathmatrix<float> two(3u, 3u, 1.f);
  ASSERT_EQ(two.det_sign(), -1);
}

TEST(mathmatrix, rank)
{
  // via https://en.wikipedia.org/wiki/Positive-definite_matrix
  mathmatrix<float> one(3u, 3u, 0.f);
  one(0, 0) = 2.;
  one(1, 0) = -1.;
  one(1, 1) = 2;
  one(2, 2) = 2;
  one(0, 1) = -1;
  one(2, 1) = -1.;
  one(1, 2) = -1;

  mathmatrix<float> three(2u, 4u, 1.f);
  three(0, 2) = 0;
  three(1, 2) = 0;
  three(1, 0) = -1;
  three(1, 1) = -1;
  three(0, 3) = 2;
  three(1, 3) = -2;
  ASSERT_EQ(three.rank(), 1);

  ASSERT_EQ(one.rank(), 3);
  mathmatrix<float> two(3u, 3u, 1.f);
  two(0, 1) = 2;
  two(1, 0) = -2;
  two(2, 0) = 3;
  two(1, 1) = -3;
  two(2, 1) = 5;
  two(2, 2) = 0;
  ASSERT_EQ(two.rank(), 2);


}


#endif