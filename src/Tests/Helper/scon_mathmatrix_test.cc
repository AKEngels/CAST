#ifdef GOOGLE_MOCK

#include <gtest/gtest.h>

#include <algorithm>
#include <random>
//#define CAST_USE_ARMADILLO
#undef eigen_assert
#define eigen_assert(x) \
  if (!(x)) { throw (std::runtime_error("Something went wrong with Eigen!")); }

#include "../../Scon/scon_mathmatrix.h"
#include <iostream>

#include <stdexcept>

auto is_eq = [](auto const& a, auto const& b) {
  return a == b;
};

TEST(SconMathmatrix, IsEqual) {
  scon::mathmatrix<double> A{
    { 1.,2.,3. },
  { 1.,2.,3. },
  { 1.,2.,3. }
  };
  scon::mathmatrix<double> B{
    { 1.,2.,3. },
  { 1.,2.,3. },
  { 1.,2.,3. }
  };

  ASSERT_EQ(A, B);
}

TEST(SconMathmatrix, Resize) {
  scon::mathmatrix<double> A(3, 3);

  ASSERT_EQ(A.cols(), 3);
  ASSERT_EQ(A.rows(), 3);

  A.resize(5, 5);

  ASSERT_EQ(A.cols(), 5);
  ASSERT_EQ(A.rows(), 5);
}

TEST(SconMathmatrix, Constrctors) {

  scon::mathmatrix<double> empty;
  ASSERT_EQ(empty.cols(), 0);
  ASSERT_EQ(empty.rows(), 0);

  scon::mathmatrix<double> A(3, 3);
  A(0, 0) = 2.0; A(0, 1) = 2.0; A(0, 2) = 2.0;
  A(1, 0) = 2.0; A(1, 1) = 2.0; A(1, 2) = 2.0;
  A(2, 0) = 2.0; A(2, 1) = 2.0; A(2, 2) = 2.0;

  scon::mathmatrix<double> B{
    { 2.,2.,2. },
  { 2.,2.,2. },
  { 2.,2.,2. }
  };

  ASSERT_EQ(A.cols(), 3);
  ASSERT_EQ(A.rows(), 3);
  ASSERT_EQ(B.cols(), 3);
  ASSERT_EQ(B.rows(), 3);

  ASSERT_EQ(A, B);

  scon::mathmatrix<double> C(std::size_t(3), std::size_t(3), 2.);

  ASSERT_EQ(A, C);

  ASSERT_EQ(C, B);

  C = scon::mathmatrix<double>{
    { 1.,1.,1. },
  { 2.,2.,2. },
  { 3.,3.,3. }
  };

  auto D(C);
  ASSERT_EQ(D, C);

  D = A;
  ASSERT_EQ(D, A);
}

TEST(SconMathmatrix, Identity) {

  scon::mathmatrix<double> A{
    { 1.,0.,0. },
  { 0.,1.,0. },
  { 0.,0.,1. }
  };

  auto B = scon::mathmatrix<double>::identity(3, 3);

  ASSERT_EQ(A, B);
}

TEST(SconMathmatrix, Zeroes) {

  scon::mathmatrix<double> A{
    { 0.,0.,0. },
  { 0.,0.,0. },
  { 0.,0.,0. }
  };

  auto B = scon::mathmatrix<double>::zero(3, 3);

  ASSERT_EQ(A, B);
}

TEST(SconMathmatrix, FillDiag) {

  scon::mathmatrix<double> A{
    { 3.,0.,0. },
  { 0.,3.,0. },
  { 0.,0.,3. }
  };

  auto B = scon::mathmatrix<double>::fill_diag(3, 3, 3.);

  ASSERT_EQ(A, B);
}

TEST(SconMathmatrix, IsVec) {
  scon::mathmatrix<double> A{
    std::initializer_list<double>{1.,},
    std::initializer_list<double>{2.,},
    std::initializer_list<double>{3.,},
  };

  EXPECT_EQ(true, A.is_vec());
}

TEST(SconMathmatrix, PositiveDefinite) {

  scon::mathmatrix<double> A{
    { 2.,-1.,0. },
  { -1.,2.,-1. },
  { 0.,-1.,2. }
  };

  ASSERT_EQ(A.positive_definite_check(), true);

  A = scon::mathmatrix<double>(std::size_t(3), std::size_t(3), 1.0);

  ASSERT_EQ(A.positive_definite_check(), false);
}

TEST(SconMathmatrix, Rank) {

  scon::mathmatrix<double> A{
    { 2.,-1.,0. },
  { -1.,2.,-1. },
  { 0.,-1.,2. }
  };

  ASSERT_EQ(A.rank(), 3);

  scon::mathmatrix<double> B{
    { 1.,2.,1. },
  { -2.,-3.,1. },
  { 3.,5.,0. }
  };

  ASSERT_EQ(B.rank(), 2);

  scon::mathmatrix<double> C{
    { 1.,1.,0.,2. },
  { -1.,-1.,0.,-2. }
  };

  ASSERT_EQ(C.rank(), 1);

}

TEST(SconMathmatrix, Det) {

  scon::mathmatrix<double> A{
    { 2.,-1.,0. },
  { -1.,2.,-1. },
  { 0.,-1.,2. }
  };

  ASSERT_EQ(A.determ(), 4.);

  scon::mathmatrix<double> B{
    { 3., 1., 1. },
  { 1.,3.,1. },
  { 1.,1.,3. }
  }, C{
    { 1.,1.,3. },
  { 1.,3.,1. },
  { 1.,1.,3. }
  };

  ASSERT_EQ(B.det_sign(), 1);
  ASSERT_EQ(C.det_sign(), -1);

}

TEST(SconMathmatrix, Plus) {

  scon::mathmatrix<double> A{
    { 2.,-1.,0. },
  { -1.,2.,-1. },
  { 0.,-1.,2. }
  };

  scon::mathmatrix<double> Aplus1{
    { 3.,0.,1. },
  { 0.,3.,0. },
  { 1.,0.,3. }
  };

  ASSERT_EQ(Aplus1, A + 1.);

  scon::mathmatrix<double> AplusA{
    { 4.,-2.,0. },
  { -2.,4.,-2. },
  { 0.,-2.,4. }
  };

  ASSERT_EQ(AplusA, A + A);

  scon::mathmatrix<double> B{
    { 1.,2.,3.,4. },
  { 1.,2.,3.,4. },
  { 1.,2.,3.,4. },
  { 1.,2.,3.,4. },
  };
  ASSERT_ANY_THROW(A + B);

}

TEST(SconMathmatrix, Minus) {

  scon::mathmatrix<double> A{
    { 2.,-1.,0. },
  { -1.,2.,-1. },
  { 0.,-1.,2. }
  };

  scon::mathmatrix<double> Aminus1{
    { 1.,-2.,-1. },
  { -2.,1.,-2. },
  { -1.,-2.,1. }
  };

  ASSERT_EQ(Aminus1, A - 1.);

  scon::mathmatrix<double> minusA{
    { -2.,1.,0. },
  { 1.,-2.,1. },
  { 0.,1.,-2. }
  };

  scon::mathmatrix<double> AminusminusA{
    { 4.,-2.,0. },
  { -2.,4.,-2. },
  { 0.,-2.,4. }
  };

  ASSERT_EQ(AminusminusA, A - minusA);

  scon::mathmatrix<double> B{
    { 1.,2.,3.,4. },
  { 1.,2.,3.,4. },
  { 1.,2.,3.,4. },
  { 1.,2.,3.,4. },
  };
  ASSERT_ANY_THROW(A - B);

}

TEST(SconMathmatrix, Mul) {

  scon::mathmatrix<double> A{
    { 2.,-1.,0. },
  { -1.,2.,-1. },
  { 0.,-1.,2. }
  };

  scon::mathmatrix<double> AtimesA{
    { 5., -4., 1. },
  { -4.,6.,-4. },
  { 1.,-4.,5. }
  };

  ASSERT_EQ(AtimesA, A * A);

  scon::mathmatrix<double> aCol{
    std::initializer_list<double>{1.,},
    std::initializer_list<double>{2.,},
    std::initializer_list<double>{3.,}
  };

  scon::mathmatrix<double> Aa{
    std::initializer_list<double>{0.,},
    std::initializer_list<double>{0.,},
    std::initializer_list<double>{4.,}
  };

  ASSERT_EQ(Aa, A * aCol);

  scon::mathmatrix<double> aRow{ { 1.,2.,3. } };

  scon::mathmatrix<double> aA{ { 0.,0.,4. } };

  ASSERT_EQ(aA, aRow * A);

  scon::mathmatrix<double> B{
    { 1.,2.,3.,4. },
  { 1.,2.,3.,4. },
  { 1.,2.,3.,4. },
  { 1.,2.,3.,4. },
  };
  ASSERT_ANY_THROW(A * B);
}

TEST(SconMathmatrix, Div) {

  scon::mathmatrix<double> A{
    { 2.,-1.,0. },
  { -1.,2.,-1. },
  { 0.,-1.,2. }
  };

  scon::mathmatrix<double> Adiv2{
    { 1., -.5, 0. },
  { -.5,1.,-.5 },
  { 0.,-.5,1. }
  };

  ASSERT_EQ(Adiv2, A / 2.);

}

TEST(SconMathmatrix, Append) {

  auto resetA = []() {
    return scon::mathmatrix<double>{
      { 7., 8., 9. },
      { 12.,13.,14. },
      { 17.,18.,19. }
    };
  };

  scon::mathmatrix<double> A = resetA(),
    appendRight{
    std::initializer_list<double>{10.,},
    std::initializer_list<double>{15.,},
    std::initializer_list<double>{20.,}
  },
    appendLeft{
    std::initializer_list<double>{6.,},
    std::initializer_list<double>{11.,},
    std::initializer_list<double>{16.,}
  },
    appendTop{
      { 1.,2.,3.,4.,5. }
  }, appendBottom{
    { 21.,22.,23.,24.,25. }
  },
    rightAppended{
      { 7.,8.,9., 10. },
  { 12.,13.,14.,15. },
  { 17.,18.,19.,20. }
  },
    leftAppended{
      { 6.,7.,8.,9., 10. },
  { 11.,12.,13.,14.,15. },
  { 16.,17.,18.,19.,20. }
  },
    topAppended{
      { 1.,2.,3.,4.,5. },
  { 6.,7.,8.,9., 10. },
  { 11.,12.,13.,14.,15. },
  { 16.,17.,18.,19.,20. }
  },
    bottomAppended{
      { 1.,2.,3.,4.,5. },
  { 6.,7.,8.,9., 10. },
  { 11.,12.,13.,14.,15. },
  { 16.,17.,18.,19.,20. },
  { 21.,22.,23.,24.,25. }
  };

  A.append_right(appendRight);

  ASSERT_EQ(A, rightAppended);

  A.append_left(appendLeft);

  ASSERT_EQ(A, leftAppended);

  A.append_top(appendTop);

  ASSERT_EQ(A, topAppended);

  A.append_bottom(appendBottom);

  ASSERT_EQ(A, bottomAppended);

  A = resetA();

  topAppended = scon::mathmatrix<double>{
    { 7.,8.,9. },
  { 12.,13.,14. },
  { 17.,18.,19. },
  { 7.,8.,9. },
  { 12.,13.,14. },
  { 17.,18.,19. }
  };
  bottomAppended = topAppended;

  auto B = A;

  A.append_top(B);

  ASSERT_EQ(A, topAppended);

  A = resetA();

  A.append_bottom(B);

  ASSERT_EQ(A, bottomAppended);

  A = resetA();

  rightAppended = scon::mathmatrix<double>{
    { 7.,8.,9.,7.,8.,9. },
  { 12.,13.,14.,12.,13.,14. },
  { 17.,18.,19.,17.,18.,19. }
  };
  leftAppended = rightAppended;

  A.append_right(B);

  ASSERT_EQ(A, rightAppended);

  A = resetA();

  A.append_left(B);

  ASSERT_EQ(A, leftAppended);
}

TEST(SconMathmatrix, ErrorIfOutOfBound) {
  scon::mathmatrix<double> A(std::size_t(3), std::size_t(3), 1.);
  ASSERT_ANY_THROW(A(5, 5));
}

TEST(SconMathmatrix, ErrorIfAppendFails) {
  scon::mathmatrix<double> A(std::size_t(3), std::size_t(3), 1.),
    B(std::size_t(2), std::size_t(2), 1.);
  ASSERT_ANY_THROW(A.append_bottom(B));
  ASSERT_ANY_THROW(A.append_top(B));
  ASSERT_ANY_THROW(A.append_left(B));
  ASSERT_ANY_THROW(A.append_right(B));
}

TEST(SconMathmatrix, ShedErrorIfOutOfBound)
{
  scon::mathmatrix<double> A(std::size_t(4u), std::size_t(4u), 5.);
  ASSERT_ANY_THROW(A.shed_rows(7, 7));
  ASSERT_ANY_THROW(A.shed_rows(-1, -1));
  ASSERT_ANY_THROW(A.shed_rows(3, 5));
  ASSERT_ANY_THROW(A.shed_cols(7, 7));
  ASSERT_ANY_THROW(A.shed_cols(-1, 1));
  ASSERT_ANY_THROW(A.shed_cols(2, 5));
}

TEST(SconMathmatrix, Shed)
{
  auto resetA = []() {
    return scon::mathmatrix<double>{
      {1., 2., 3., 4.},
      { 5.,6.,7.,8. },
      { 9.,10.,11.,12 },
      { 13.,14.,15.,16. }
    };
  };

  scon::mathmatrix<double> A = resetA(),
    B{
      { 1.,2.,3.,4. },
  { 13.,14.,15.,16. }
  };

  A.shed_rows(1, 2);

  ASSERT_EQ(A.rows(), 2u);
  ASSERT_EQ(A, B);

  A = resetA();

  B = scon::mathmatrix<double>{
    { 1.,4. },
  { 5.,8. },
  { 9.,12. },
  { 13.,16. },
  };

  A.shed_cols(1, 2);

  ASSERT_EQ(A.cols(), 2u);
  ASSERT_EQ(A, B);

  A = resetA();

  A.shed_cols(2);

  B = scon::mathmatrix<double>{
    { 1., 2., 4. },
  { 5.,6.,8. },
  { 9.,10.,12 },
  { 13.,14.,16. }
  };

  ASSERT_EQ(A.cols(), 3u);
  ASSERT_EQ(A, B);

  A = resetA();

  A.shed_rows(2);

  B = scon::mathmatrix<double>{
    { 1., 2., 3., 4. },
  { 5.,6.,7.,8. },
  { 13.,14.,15.,16. }
  };

  ASSERT_EQ(A.rows(), 3u);
  ASSERT_EQ(A, B);

}

TEST(SconMathmatrix, ReturnQuadratic)
{
  scon::mathmatrix<double> A(std::size_t(4), std::size_t(4), 5.);
  scon::mathmatrix<double> B(std::size_t(4), std::size_t(1), 0.);

  ASSERT_EQ(A.return_quadratic(), true);
  ASSERT_EQ(B.return_quadratic(), false);

}

TEST(SconMathmatrix, Transpose) {
  scon::mathmatrix<double> A{
    { 1.,1.,1. },
  { 2.,2.,2. },
  { 3.,3.,3. }
  },
    B{
      { 1.,2.,3. },
  { 1.,2.,3. },
  { 1.,2.,3. }
  }, a{
    std::initializer_list<double>{1.,},
    std::initializer_list<double>{2.,},
    std::initializer_list<double>{3.,},
  }, b{
    { 1.,2.,3. }
  };

  auto C = A.t();
  ASSERT_EQ(C, B);

  auto c = a.t();
  ASSERT_EQ(c, b);
}

TEST(SconMathmatrix, UpperLeftSubmatrix)
{
  scon::mathmatrix<double> A{
    { 1.,1.,1.,1. },
  { 2.,2.,2.,2. },
  { 3.,3.,3.,3. },
  { 4.,4.,4.,4. }
  }, testA{
    { 1.,1. },
  { 2.,2. }
  };

  scon::mathmatrix<double> subA = A.upper_left_submatrix(2);

  ASSERT_EQ(subA.return_quadratic(), true);
  ASSERT_EQ(subA, testA);
  ASSERT_EQ(subA.rows(), 2u);

  scon::mathmatrix<double> B{
    { 1.,1.,1.,1.,1.,1. },
  { 2.,2.,2.,2.,2.,2. },
  { 3.,3.,3.,3.,3.,3. },
  { 4.,4.,4.,4.,4.,4. },
  { 5.,5.,5.,5.,5.,5. },
  { 6.,6.,6.,6.,6.,6. },
  }, testB{
    { 1.,1.,1.,1. },
  { 2.,2.,2.,2. },
  { 3.,3.,3.,3. }
  };

  scon::mathmatrix<double> subB = B.upper_left_submatrix(3, 4);

  ASSERT_EQ(subB, testB);
  ASSERT_EQ(subB.rows(), 3u);
  ASSERT_EQ(subB.cols(), 4u);
}

TEST(SconMathmatrix, SVD) {
  scon::mathmatrix<double> A{
    { 2.,-1.,0. },
  { -1.,2.,-1. },
  { 0.,-1.,2. }
  }, Stest{
    std::initializer_list<double>{3.41421,},
    std::initializer_list<double>{2.00000,},
    std::initializer_list<double>{0.58579,}
  },
    I = scon::mathmatrix<double>::identity(3, 3),
    U, s, V;



  std::tie(U, s, V) = A.svd();

  ASSERT_EQ(s, Stest);
  ASSERT_EQ(A, U * s.diagmat() * V.t());
  ASSERT_EQ(I, V.t() * V);
  ASSERT_EQ(I, U.t() * U);

  A.singular_value_decomposition(U, s, V);

  ASSERT_EQ(s, Stest);
  ASSERT_EQ(A, U * s.diagmat() * V.t());
  ASSERT_EQ(I, V.t() * V);
  ASSERT_EQ(I, U.t() * U);

}

TEST(SconMathmatrix, Eigensym) {
  scon::mathmatrix<double> A{
    { 2.,-1.,0. },
  { -1.,2.,-1. },
  { 0.,-1.,2. }
  }, lambdaTest{
    std::initializer_list<double>{0.58579,},
    std::initializer_list<double>{2.00000,},
    std::initializer_list<double>{3.41421,}
  },
    I = scon::mathmatrix<double>::identity(3, 3),
    EVal, EVec;



  std::tie(EVal, EVec) = A.eigensym();

  auto B = EVec.t() * A * EVec;

  ASSERT_EQ(lambdaTest, EVal);
  ASSERT_EQ(EVal.diagmat(), B);
  ASSERT_EQ(I, EVec.t() * EVec);

  B = lambdaTest.diagmat();

  A.diag();

  ASSERT_EQ(A, B);

}

class TestRowsAndCols {
private:
  static const scon::mathmatrix<double> A;
  FRIEND_TEST(SconMathmatrix, RowCol);
  FRIEND_TEST(SconMathmatrix, DoesColWork);
  FRIEND_TEST(SconMathmatrix, DoesRowWork);
  FRIEND_TEST(SconMathmatrix, DoesColThrow);
  FRIEND_TEST(SconMathmatrix, DoesRowThrow);
};

const scon::mathmatrix<double> TestRowsAndCols::A{
  { 1.,2.,3. },
{ 4.,5.,6. },
{ 7.,8.,9. }
};


TEST(SconMathmatrix, DoesColWork) {
  scon::mathmatrix<double> ac{
    std::initializer_list<double>{1.,},
    std::initializer_list<double>{4.,},
    std::initializer_list<double>{7.,}
  };
  ASSERT_EQ(TestRowsAndCols::A.col(0), ac);
}

TEST(SconMathmatrix, DoesRowWork) {
  scon::mathmatrix<double> ar{
    { 1.,2.,3. }
  };
  ASSERT_EQ(TestRowsAndCols::A.row(0), ar);
}

TEST(SconMathmatrix, DoesColThrow) {
  ASSERT_ANY_THROW(TestRowsAndCols::A.col(5));
}

TEST(SconMathmatrix, DoesRowThrow) {
  ASSERT_ANY_THROW(TestRowsAndCols::A.row(5));
}

TEST(SconMathmatrix, ColFromVec) {
  std::vector<double> Avec{ 1.,2.,3.,4.,5.,6.,7.,8.,9. };

  scon::mathmatrix<double> A = scon::mathmatrix<double>::col_from_vec(Avec);

  scon::mathmatrix<double> sorted{
    std::initializer_list<double>{1.,},
    std::initializer_list<double>{2.,},
    std::initializer_list<double>{3.,},
    std::initializer_list<double>{4.,},
    std::initializer_list<double>{5.,},
    std::initializer_list<double>{6.,},
    std::initializer_list<double>{7.,},
    std::initializer_list<double>{8.,},
    std::initializer_list<double>{9.,}
  };


  ASSERT_EQ(A, sorted);
}

TEST(SconMathmatrix, ToStdVec) {
  std::vector<double> Avec{ 1.,2.,3.,4.,5.,6.,7.,8.,9. };
  scon::mathmatrix<double> A{
    std::initializer_list<double>{1.,},
    std::initializer_list<double>{2.,},
    std::initializer_list<double>{3.,},
    std::initializer_list<double>{4.,},
    std::initializer_list<double>{5.,},
    std::initializer_list<double>{6.,},
    std::initializer_list<double>{7.,},
    std::initializer_list<double>{8.,},
    std::initializer_list<double>{9.,}
  };

  auto vec = A.col_to_std_vector();

  EXPECT_EQ(vec, Avec);

  A = scon::mathmatrix<double>{
    { 1., 2., 3., 4., 5., 6., 7., 8., 9. }
  };

  auto row = A.row_to_std_vector();

  EXPECT_EQ(row, Avec);

  A = scon::mathmatrix<double>{
    { 1.,2.,3., },
  { 4.,5.,6., },
  { 7.,8.,9., },
  };

  vec = A.col_to_std_vector(2);
  row = A.row_to_std_vector(2);

  Avec = std::vector<double>{ 3., 6., 9. };
  EXPECT_EQ(Avec, vec);

  Avec = std::vector<double>{ 7., 8., 9. };
  EXPECT_EQ(Avec, row);

  std::vector<std::vector<double>> Avecvec{
    { 1.,2.,3., },
  { 4.,5.,6., },
  { 7.,8.,9., },
  };

  auto vecvec = A.to_std_vector();

  EXPECT_EQ(vecvec, Avecvec);
}

TEST(SconMathmatrix, SortCol) {
  std::vector<double> Avec{ 1.,2.,3.,4.,5.,6.,7.,8.,9. };

  std::shuffle(Avec.begin(), Avec.end(), std::mt19937(std::random_device()()));

  scon::mathmatrix<double> A = scon::mathmatrix<double>::col_from_vec(Avec), sorted{
    std::initializer_list<double>{1.,},
    std::initializer_list<double>{2.,},
    std::initializer_list<double>{3.,},
    std::initializer_list<double>{4.,},
    std::initializer_list<double>{5.,},
    std::initializer_list<double>{6.,},
    std::initializer_list<double>{7.,},
    std::initializer_list<double>{8.,},
    std::initializer_list<double>{9.,}
  };

  std::vector<double> asc_vec{ 1.,2.,3.,4.,5.,6.,7.,8.,9. };
  auto sorted_vec = A.sort_col_to_vec([](auto const& a, auto const& b) {
    return a < b;
  });
  EXPECT_EQ(sorted_vec, asc_vec);
  EXPECT_PRED2(is_eq, A.sort_col([](auto const& a, auto const& b) {
    return a < b;
  }), sorted);
  EXPECT_PRED2(is_eq, A.sort_col_asc(), sorted);

  sorted = scon::mathmatrix<double>{
    std::initializer_list<double>{9.,},
    std::initializer_list<double>{8.,},
    std::initializer_list<double>{7.,},
    std::initializer_list<double>{6.,},
    std::initializer_list<double>{5.,},
    std::initializer_list<double>{4.,},
    std::initializer_list<double>{3.,},
    std::initializer_list<double>{2.,},
    std::initializer_list<double>{1.,}
  };

  EXPECT_PRED2(is_eq, A.sort_col_disc(), sorted);
}

TEST(SconMathmatrix, SortIndex) {
  std::vector<std::size_t> Avec_sorted(9);
  std::iota(Avec_sorted.begin(), Avec_sorted.end(), 0);

  std::vector<std::size_t> Avec = Avec_sorted;

  std::shuffle(Avec.begin(), Avec.end(), std::mt19937(std::random_device()()));

  std::vector<std::size_t> Avec_indices = Avec_sorted;

  std::sort(Avec_indices.begin(), Avec_indices.end(), [&](auto const& a, auto const& b) {
    return Avec.at(a) < Avec.at(b);
  });

  scon::mathmatrix<double> A{
    std::initializer_list<double>{static_cast<double>(Avec.at(0)),},
    std::initializer_list<double>{static_cast<double>(Avec.at(1)),},
    std::initializer_list<double>{static_cast<double>(Avec.at(2)),},
    std::initializer_list<double>{static_cast<double>(Avec.at(3)),},
    std::initializer_list<double>{static_cast<double>(Avec.at(4)),},
    std::initializer_list<double>{static_cast<double>(Avec.at(5)),},
    std::initializer_list<double>{static_cast<double>(Avec.at(6)),},
    std::initializer_list<double>{static_cast<double>(Avec.at(7)),},
    std::initializer_list<double>{static_cast<double>(Avec.at(8)),}
  };

  auto idx = A.sort_idx();
  EXPECT_EQ(idx.size(), Avec_indices.size());
  EXPECT_EQ(idx, Avec_indices);
}

TEST(SconMathmatrix, FindIdx) {
  scon::mathmatrix<double> A{
    std::initializer_list<double>{.5,},
    std::initializer_list<double>{3.,},
    std::initializer_list<double>{.2,},
    std::initializer_list<double>{1.,},
    std::initializer_list<double>{11.,},
    std::initializer_list<double>{.25,},
    std::initializer_list<double>{27.,},
    std::initializer_list<double>{-9.,},
    std::initializer_list<double>{0.,}
  };

  std::vector<std::size_t> foundTest{ 0,2,5,7,8, };

  auto found = A.find_idx([](auto const& a) {
    return a < 1.;
  });

  EXPECT_EQ(found, foundTest);
}

TEST(SconMathmatrix, SubmatTest) {
  scon::mathmatrix<double> A{
    { 1.,  2.,  3.,  4.,  5., },
  { 6.,  7.,  8.,  9., 10., },
  { 11., 12., 13., 14., 15., },
  { 16., 17., 18., 19., 20., },
  { 21., 22., 23., 24., 25., },
  }, subVec{
    { 1.,  2.,  4., },
  { 11., 12., 14., },
  { 16., 17., 19., },
  }, subInd{
    { 7.,  8.,  9., },
  { 12., 13., 14., },
  { 17., 18., 19., },
  };

  std::vector<std::size_t> cols{ 0,1,3 }, rows{ 0,2,3 };

  auto TestVec = A.submat(rows, cols);

  EXPECT_EQ(TestVec, subVec);

  TestVec = A.submat(1, 1, 3, 3);

  EXPECT_EQ(TestVec, subInd);

  TestVec = A.submat(0, 0, 2, 2);

  subInd = scon::mathmatrix<double>{
    { 1.,  2.,  3., },
  { 6.,  7.,  8., },
  { 11., 12., 13., },
  };

  EXPECT_EQ(TestVec, subInd);

  TestVec = A.submat(0, 1, 2, 3);

  subInd = scon::mathmatrix<double>{
    { 2.,  3.,  4., },
  { 7.,  8.,  9., },
  { 12., 13., 14., },
  };

  EXPECT_EQ(TestVec, subInd);

  TestVec = A.submat(1, 0, 3, 2);

  subInd = scon::mathmatrix<double>{
    { 6.,  7.,  8., },
  { 11., 12., 13., },
  { 16., 17., 18., },
  };

  EXPECT_EQ(TestVec, subInd);

}

TEST(SconMathmatrix, PseudoInverse) {

  scon::mathmatrix<double> A{
    { 1.,  2.,  3., },
  { 4.,  5.,  6., },
  { 7.,  8.,  9., },
  },
  //result from octave
  pinv{
    { -6.3889e-001, -1.6667e-001,  3.0556e-001 },
  { -5.5556e-002,  3.6675e-017,  5.5556e-002 },
  { 5.2778e-001,  1.6667e-001, -1.9444e-001 }
  };

  auto B = A.pinv();

  EXPECT_EQ(A, A * B * A);
  EXPECT_EQ(B, B * A * B);
  EXPECT_EQ(B, pinv);

}

//TEST(SconMathmatrix, Vectorise) {
//
//	scon::mathmatrix<double> A{
//		{ 1.,  2.,  3., },
//	{ 4.,  5.,  6., },
//	{ 7.,  8.,  9., },
//	}, vec_test{
//		std::initializer_list<double>{1.,},
//		std::initializer_list<double>{2.,},
//		std::initializer_list<double>{3.,},
//		std::initializer_list<double>{4.,},
//		std::initializer_list<double>{5.,},
//		std::initializer_list<double>{6.,},
//		std::initializer_list<double>{7.,},
//		std::initializer_list<double>{8.,},
//		std::initializer_list<double>{9.,}
//	};
//
//	auto vec = A.vectorise();
//
//	EXPECT_EQ(vec, vec_test);
//}

TEST(SconMathmatrix, ReplaceItemsWithIndex) {
  scon::mathmatrix<double> A{
    { 1.,  2.,  3., },
  { 4.,  5.,  6., },
  { 7.,  8.,  9., },
  };
  scon::mathmatrix<double> replaced{
    { 1.,  2.,  3., },
  { 4.,  0.,  0., },
  { 0.,  8.,  9., },
  };

  std::vector<std::size_t> indices{ 2, 4, 7 };

  A.replace_idx_with(indices, 0.0);

  EXPECT_EQ(A, replaced);

  A = scon::mathmatrix<double>{
    { 1., 2., 3., 4., 5., 6., 7., 8., 9., },
  };
  replaced = scon::mathmatrix<double>{
    { 1., 2., 0., 4., 0., 6., 7., 0., 9., },
  };

  A.replace_idx_with(indices, 0.0);

  EXPECT_EQ(A, replaced);

}

TEST(SconMathmatrix, Elem) {
  scon::mathmatrix<double> A{
    { 1.,  2.,  3., },
  { 4.,  5.,  6., },
  { 7.,  8.,  9., },
  };

  std::vector<std::size_t> indices{ 2, 4, 7 };

  std::vector<double> elements_ind{
    A(indices.at(0) % A.rows(), indices.at(0) / A.rows()),
    A(indices.at(1) % A.rows(), indices.at(1) / A.rows()),
    A(indices.at(2) % A.rows(), indices.at(2) / A.rows()),
  }, elements_container;

  elements_container.reserve(indices.size());
  auto elements = A.elem(indices);

  for (auto const& el : elements) {
    elements_container.emplace_back(el.get());
  }

  EXPECT_EQ(elements_container, elements_ind);

  scon::mathmatrix<double> B{
    { 1.,  2.,  3., },
  { 4.,  6.,  7., },
  { 8.,  8.,  9., },
  };

  for (auto& el : elements) {
    el.get() += 1.0;
  }

  EXPECT_EQ(A, B);
}

TEST(SconMathmatrix, DiagMat) {
  scon::mathmatrix<double> diag_elms{
    std::initializer_list<double>{1.,},
    std::initializer_list<double>{2.,},
    std::initializer_list<double>{3.,},
  }, diagmat{
    { 1., 0., 0., },
  { 0., 2., 0., },
  { 0., 0., 3., },
  }, A{
    { 1., 1., 1., },
  { 1., 2., 1., },
  { 1., 1., 3., },
  };

  EXPECT_EQ(diagmat, diag_elms.diagmat());
  EXPECT_EQ(diagmat, A.diagmat());

}

TEST(SconMathmatrix, Solve) {
  scon::mathmatrix<double> A{
          {1., 0., -1.,  1.},
          {0., 1.,  2.,  0.},
          {0., 0., -3.,  2.},
          {0., 1.,  4., -1.}
  };
  auto y = scon::mathmatrix<double>::col_from_vec({1., 0.,  2., -2.});
  auto x = scon::mathmatrix<double>::col_from_vec({1., 4., -2., -2.});

  EXPECT_EQ(x, A.solve(y));
}

#endif