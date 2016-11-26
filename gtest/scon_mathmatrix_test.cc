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

TEST(mathmatrix, determ)
{
  mathmatrix<float> one(3u, 3u, 0.f);
  one(0, 0) = 2.;
  one(1, 0) = -1.;
  one(1, 1) = 2;
  one(2, 2) = 2;
  one(0, 1) = -1;
  one(2, 1) = -1.;
  one(1, 2) = -1;
  ASSERT_FLOAT_EQ(one.determ(), 4.f);

  mathmatrix<float> two(3u, 3u, 1.f);
  ASSERT_EQ(two.determ(), 0);

  mathmatrix<float> three(4u, 4u, 5.f);
  three(0, 0) = 3;
  three(0, 2) = 8;
  three(0, 3) = 15;
  three(1, 0) = -4;
  three(1, 1) = 6;
  three(1, 2) = -3.5;
  three(1, 3) = 9;
  three(2, 1) = -5;
  three(2, 3) = 4;
  three(3, 0) = 1;
  three(3, 1) = 2;
  three(3, 2) = 3;
  three(3, 3) = 4;
  ASSERT_NEAR(three.determ(), 0.5, 0.001);
}

TEST(mathmatrix, minusOperatorWorksReasonably)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  one(0, 0) = 3;
  one(0, 2) = 8;
  one(0, 3) = 15;
  one(1, 0) = -4;
  one(1, 1) = 6;
  one(1, 2) = -3.5;
  one(1, 3) = 9;
  one(2, 1) = -5;
  one(2, 3) = 4;
  one(3, 0) = 1;
  one(3, 1) = 2;
  one(3, 2) = 3;
  one(3, 3) = 4;

  mathmatrix<float> two(4u, 4u, 3.f);

  mathmatrix<float> three(4u, 4u, 0.f);
  three(0, 1) = 2;
  three(0, 2) = 5;
  three(0, 3) = 12;
  three(1, 0) = -7;
  three(1, 1) = 3;
  three(1, 2) = -6.5;
  three(1, 3) = 6;
  three(2, 0) = 2;
  three(2, 1) = -8;
  three(2, 2) = 2;
  three(2, 3) = 1;
  three(3, 0) = -2;
  three(3, 1) = -1;
  three(3, 2) = 0;
  three(3, 3) = 1;
  ASSERT_EQ(one - two, three);

}

TEST(mathmatrix, minusOperatorThrowsAtSizeMismatch)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  mathmatrix<float> two(3u, 2u, 5.f);

  ASSERT_ANY_THROW(one - two);

}

TEST(mathmatrix, choleskyDecomposition)
{
  mathmatrix<float> one(3u, 3u, 0.f);
  one(0, 0) = 2.;
  one(1, 0) = -1.;
  one(1, 1) = 2;
  one(2, 2) = 2;
  one(0, 1) = -1;
  one(2, 1) = -1.;
  one(1, 2) = -1;

  mathmatrix<float> reference(3u, 3u, 0.f);
  reference(0, 0) = 1.41421354e+00;
  reference(0, 1) = -7.07106769e-01;
  reference(1, 1) = 1.22474492e+00;
  reference(1, 2) = -8.16496551e-01;
  reference(2, 2) = 1.15470052e+00;

  mathmatrix<float> result;
  one.choleskyDecomposition(result);

  for (int i = 0; i < result.rows(); i++)
    for (int j = i; j < result.rows(); j++)
      ASSERT_FLOAT_EQ(result(i, j), reference(i, j));

}

TEST(mathmatrix, divisonByFloatOperatorWorksReasonably)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  one(0, 0) = 3;
  one(0, 2) = 8;
  one(0, 3) = 15;
  one(1, 0) = -4;
  one(1, 1) = 6;
  one(1, 2) = -3.5;
  one(1, 3) = 9;
  one(2, 1) = -5;
  one(2, 3) = 4;
  one(3, 0) = 1;
  one(3, 1) = 2;
  one(3, 2) = 3;
  one(3, 3) = 4;

  mathmatrix<float> two(one);
  for (size_t i = 0u; i < two.rows(); i++)
    for (size_t j = 0u; j < two.cols(); j++)
      two(i, j) = two(i, j) / 3.5;

  ASSERT_EQ(one / 3.5 , two);
}

TEST(mathmatrix, appendFunctionsThrowUponSizeMissmatch)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  mathmatrix<float> two(1u, 2u, 5.f);

  ASSERT_ANY_THROW(one.append_bottom(two));
  ASSERT_ANY_THROW(one.append_left(two));
  ASSERT_ANY_THROW(one.append_right(two));
  ASSERT_ANY_THROW(one.append_top(two));
}

TEST(mathmatrix, appendTopWorksCorrectly)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  mathmatrix<float> two(1u, 4u, 0.f);

  one.append_top(two);
  ASSERT_EQ(one.rows(), 5u);
  ASSERT_EQ(one.cols(), 4u);
  ASSERT_FLOAT_EQ(one(0, 0), 0.f);
  ASSERT_FLOAT_EQ(one(1, 1), 5.f);
}

TEST(mathmatrix, appendBottomWorksCorrectly)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  mathmatrix<float> two(1u, 4u, 0.f);

  one.append_bottom(two);
  ASSERT_EQ(one.rows(), 5u);
  ASSERT_EQ(one.cols(), 4u);
  ASSERT_FLOAT_EQ(one(0, 0), 5.f);
  ASSERT_FLOAT_EQ(one(4, 3), 0.f);
}

TEST(mathmatrix, throwsWhenArgumentsToRoundBracketOperatorIsOutOfRange)
{
  mathmatrix<float> one(2u, 2u, 1.f);
  ASSERT_ANY_THROW(one(2, 2));
  ASSERT_ANY_THROW(one(-1, 1));
}

TEST(mathmatrix, appendLeftWorksCorrectly)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  mathmatrix<float> two(4u, 1u, 0.f);

  one.append_left(two);
  ASSERT_EQ(one.rows(), 4u);
  ASSERT_EQ(one.cols(), 5u);
  ASSERT_FLOAT_EQ(one(0, 0), 0.f);
  ASSERT_FLOAT_EQ(one(1, 1), 5.f);
}

TEST(mathmatrix, appendRightWorksCorrectly)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  mathmatrix<float> two(4u, 1u, 0.f);

  one.append_right(two);
  ASSERT_EQ(one.rows(), 4u);
  ASSERT_EQ(one.cols(), 5u);
  ASSERT_FLOAT_EQ(one(0, 0), 5.f);
  ASSERT_FLOAT_EQ(one(3, 4), 0.f);
}

TEST(mathmatrix, shedFunctionsThrowWhenOutOfBounds)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  ASSERT_ANY_THROW(one.shed_rows(7, 7));
  ASSERT_ANY_THROW(one.shed_rows(-1, -1));
  ASSERT_ANY_THROW(one.shed_rows(3, 5));
  ASSERT_ANY_THROW(one.shed_cols(7, 7));
  ASSERT_ANY_THROW(one.shed_cols(-1, 1));
  ASSERT_ANY_THROW(one.shed_cols(2, 5));
}

TEST(mathmatrix, shedRowsWorksCorrectly)
{
  mathmatrix<float> one(4u, 4u, 0.f);
  for (size_t i = 0u; i < one.rows(); i++)
  {
    for (size_t j = 0u; j < one.cols(); j++)
      one(i, j) = static_cast<float>(i);
  }

  one.shed_rows(2, 2);
  ASSERT_EQ(one.rows(), 3u);
  ASSERT_FLOAT_EQ(one(0, 0), 0.);
  ASSERT_FLOAT_EQ(one(1, 0), 1.);
  ASSERT_FLOAT_EQ(one(2, 0), 3.);

  mathmatrix<float> two(7u, 7u, 0.f);
  for (size_t i = 0u; i < two.rows(); i++)
  {
    for (size_t j = 0u; j < two.cols(); j++)
      two(i, j) = static_cast<float>(i);
  }

  two.shed_rows(2, 4);
  ASSERT_EQ(two.rows(), 4u);
  ASSERT_FLOAT_EQ(two(0, 0), 0.);
  ASSERT_FLOAT_EQ(two(1, 0), 1.);
  ASSERT_FLOAT_EQ(two(2, 0), 5.);
  ASSERT_FLOAT_EQ(two(3, 1), 6.);
}

TEST(mathmatrix, shedColsWorksCorrectly)
{
  mathmatrix<float> one(4u, 4u, 0.f);
  for (size_t i = 0u; i < one.rows(); i++)
  {
    for (size_t j = 0u; j < one.cols(); j++)
      one(i, j) = static_cast<float>(j);
  }

  one.shed_cols(2, 2);
  ASSERT_EQ(one.cols(), 3u);
  ASSERT_FLOAT_EQ(one(0, 0), 0.);
  ASSERT_FLOAT_EQ(one(0, 1), 1.);
  ASSERT_FLOAT_EQ(one(0, 2), 3.);

  mathmatrix<float> two(7u, 7u, 0.f);
  for (size_t i = 0u; i < two.rows(); i++)
  {
    for (size_t j = 0u; j < two.cols(); j++)
      two(i, j) = static_cast<float>(j);
  }

  two.shed_cols(2, 4);
  ASSERT_EQ(two.cols(), 4u);
  ASSERT_FLOAT_EQ(two(0, 0), 0.);
  ASSERT_FLOAT_EQ(two(0, 1), 1.);
  ASSERT_FLOAT_EQ(two(0, 2), 5.);
  ASSERT_FLOAT_EQ(two(1, 3), 6.);
}

TEST(mathmatrix, returnQuadraticWorksCorrectly)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  mathmatrix<float> two(4u, 1u, 0.f);

  ASSERT_EQ(one.return_quadratic(), true);
  ASSERT_EQ(two.return_quadratic(), false);

}

TEST(mathmatrix, upperLeftSubmatrixWorksCorrectly)
{
  mathmatrix<float> one(4u, 4u, 0.f);
  for (size_t i = 0u; i < one.rows(); i++)
  {
    for (size_t j = 0u; j < one.cols(); j++)
      one(i, j) = static_cast<float>(i);
  }

  mathmatrix<float> sub = one.upper_left_submatrix(2);
  ASSERT_EQ(sub.return_quadratic(), true);
  ASSERT_FLOAT_EQ(sub(0, 0), 0.);
  ASSERT_FLOAT_EQ(sub(1, 1), 1.);
  ASSERT_EQ(sub.rows(), 2u);

  mathmatrix<float> one1(6u, 6u, 0.f);
  for (size_t i = 0u; i < one1.rows(); i++)
  {
    for (size_t j = 0u; j < one1.cols(); j++)
      one1(i, j) = static_cast<float>(i);
  }

  mathmatrix<float> sub1 = one1.upper_left_submatrix(2, 3);
  ASSERT_FLOAT_EQ(sub1(0, 0), 0.);
  ASSERT_FLOAT_EQ(sub1(1, 1), 1.);
  ASSERT_EQ(sub1.rows(), 2u);
  ASSERT_EQ(sub1.cols(), 3u);
}

TEST(mathmatrix, matrixMultiplicationThrowsAtSizeMissmatch)
{
  mathmatrix<float> one(4u, 3u, 5.f);

  mathmatrix<float> two(4u, 5u, 3.f);


  ASSERT_ANY_THROW(one*two);
}

TEST(mathmatrix, matrixMultiplicationAndEqualityOperatorWorkingReasonably)
{
  mathmatrix<float> one(4u, 4u, 5.f);
  one(0, 0) = 3;
  one(0, 2) = 8;
  one(0, 3) = 15;
  one(1, 0) = -4;
  one(1, 1) = 6;
  one(1, 2) = -3.5;
  one(1, 3) = 9;
  one(2, 1) = -5;
  one(2, 3) = 4;
  one(3, 0) = 1;
  one(3, 1) = 2;
  one(3, 2) = 3;
  one(3, 3) = 4;

  mathmatrix<float> two(4u, 3u, 3.f);
  two(1, 1) = -0.2;

  mathmatrix<float> three(4u, 3u, 0.f);
  three(0, 0) = 93;
  three(0, 1) = 77;
  three(0, 2) = 93;
  three(1, 0) = 22.5;
  three(1, 1) = 3.3;
  three(1, 2) = 22.5;
  three(2, 0) = 27;
  three(2, 1) = 43;
  three(2, 2) = 27;
  three(3, 0) = 30;
  three(3, 1) = 23.6;
  three(3, 2) = 30;

  mathmatrix<float> result = one * two;
  for (size_t i = 0u; i < three.rows(); i++)
  {
    for (size_t j = 0u; j < three.cols(); j++)
    {
      ASSERT_FLOAT_EQ(result(i,j), three(i,j));
    }
  }



  // OPERATOR==
  ASSERT_TRUE(result == three);
}
#endif