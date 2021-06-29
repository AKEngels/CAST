#ifdef GOOGLE_MOCK

#include <gtest/gtest.h>

#include <algorithm>
#include <random>
//#define CAST_USE_ARMADILLO
#undef eigen_assert
#define eigen_assert(x) \
  if (!(x)) { throw (std::runtime_error("Something went wrong with Eigen!")); }

#include "../.././Scon/scon_mathmatrix.h"
#include "../../coords.h"
#include <iostream>

#include <stdexcept>

#include "../../entropy.h"

// tests use the test system butanol.arc

TEST(entropy, andersenDarlingTestWorking)
{
  //
  std::vector<coords::float_type> stdVecIn = { 20.000 ,21.000 ,21.500 ,19.000 ,19.000 ,20.400 ,18.300 ,19.900 ,\
  18.700 ,18.000 ,17.700 ,19.100 ,19.700 ,18.100 ,18.400 ,17.500 ,18.900 ,19.000 ,20.500 ,17.300 ,18.300 ,18.400 ,\
  18.600 ,19.800 ,20.200 ,18.500 ,18.500 ,18.000 ,20.900 ,18.100 ,19.400 ,20.500 ,20.400 ,16.100 ,18.700 ,18.800 ,\
  17.300 ,18.100 ,19.900 ,19.600 };
  auto drawMatrix = scon::mathmatrix<coords::float_type>::col_from_vec(stdVecIn);
  auto andersonDarlingValues = entropy::normalityCheck(drawMatrix, drawMatrix.cols());
  ASSERT_NEAR(andersonDarlingValues.at(0u),0.327934,0.001);
}
#endif