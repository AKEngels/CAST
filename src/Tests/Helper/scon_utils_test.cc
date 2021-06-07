#ifdef GOOGLE_MOCK
#include <gtest/gtest.h>

#include "../../Scon/scon_utility.h"

TEST(SconUtilities, Matches) {
  std::vector a{1, 2, 3}, b{2, 1, 3}, c{1, 2, 4}, s{1, 2};

  EXPECT_TRUE(scon::match(a.begin(), a.end(), a.begin(), a.end()));
  EXPECT_TRUE(scon::match(a.begin(), a.end(), b.begin(), b.end()));
  EXPECT_FALSE(scon::match(a.begin(), a.end(), c.begin(), c.end()));

  EXPECT_FALSE(scon::match(a.begin(), a.end(), s.begin(), s.end()));
  EXPECT_FALSE(scon::match(s.begin(), s.end(), a.begin(), a.end()));
}
#endif