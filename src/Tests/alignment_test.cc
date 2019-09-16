#ifdef GOOGLE_MOCK

#include <gtest/gtest.h>

#include <algorithm>
#include <random>
//#define CAST_USE_ARMADILLO
#undef eigen_assert
#define eigen_assert(x) \
  if (!(x)) { throw (std::runtime_error("Something went wrong with Eigen!")); }

#include "../Scon/scon_mathmatrix.h"
#include "../alignment.h"
#include "../coords.h"
#include <iostream>

#include <stdexcept>

#include "../coords_io.h"
#include "../tinker_parameters.h"
#include "../energy_int_aco.h"

// tests use the test system butanol.arc

TEST(alignment, kabschAlignmentLeavesStructureUnchanged)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  const coords::Coordinates coords(ci->read("test_files/butanol.arc"));
  const coords::Coordinates com_aligned = align::centerOfGeometryAligned(coords);

  const coords::Coordinates after(align::kabschAligned(com_aligned, com_aligned));
  constexpr double maxDiffAngstrom = 10e-5;
  for (std::size_t i = 0u; i < after.xyz().size(); i++)
  {
    ASSERT_NEAR(after.xyz().at(i).x(), com_aligned.xyz().at(i).x(), maxDiffAngstrom);
    ASSERT_NEAR(after.xyz().at(i).y(), com_aligned.xyz().at(i).y(), maxDiffAngstrom);
    ASSERT_NEAR(after.xyz().at(i).z(), com_aligned.xyz().at(i).z(), maxDiffAngstrom);
  }
}
#endif