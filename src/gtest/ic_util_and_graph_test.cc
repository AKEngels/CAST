#ifdef GOOGLE_MOCK

#include "ic_util_and_graph_test.h"
#include "../ic_util.h"

namespace {
  double constexpr doubleNearThreshold = 1.e-10;
  inline void isCartesianPointNear(coords::r3 const& lhs, coords::r3 const& rhs) {
    EXPECT_NEAR(lhs.x(), rhs.x(), doubleNearThreshold);
    EXPECT_NEAR(lhs.y(), rhs.y(), doubleNearThreshold);
    EXPECT_NEAR(lhs.z(), rhs.z(), doubleNearThreshold);
  }
}
TestCreateGraph::TestCreateGraph() : atomVector{
  {1, "C" ,"C", coords::r3{ -0.321, -0.087, 0.12733 } },
{2, "O", "O", coords::r3{ 1.055, 0.144, 0.21133 } },
{3, "H" , "H" , coords::r3{ -0.53, -0.921, -0.57767 } },
{4, "H" , "H" , coords::r3{ -0.85, 0.837, -0.18867 } },
{5,"H","H", coords::r3{ -0.699, -0.376, 1.12933 } },
{6, "H" , "H" , coords::r3{ 1.345, 0.403, -0.70167 } }
}, connectivity{
  {1,2},
{ 1,3 },
{ 1,4 },
{ 1,5 },
{ 2,6 }
} {}

#endif

KeyOrDefaultTest::KeyOrDefaultTest() : testMap{ {"one", 1u }, {"three", 3u}, {"ten", 10u} }{}

void KeyOrDefaultTest::testKeyOrDefault(){
  auto value = ic_util::getValueByKeyOrDefault(testMap, GetParam().key, GetParam().defaultValue);
  EXPECT_EQ(GetParam().expectedValue == value, GetParam().isExpectedValue);
}

TEST_P(KeyOrDefaultTest, testKeyOrDefault) {
  testKeyOrDefault();
}

INSTANTIATE_TEST_CASE_P(getValueByKeyOrDefault, KeyOrDefaultTest, testing::Values(
  expectedValuesForMap{ "one", 1u, 0u, true},
  expectedValuesForMap{ "two", 2u, 0u, false }
));

GetMeanTest::GetMeanTest() : randomCoordinates{
  coords::r3{ 0.243315, 0.96908646, 0.93456417 },
  coords::r3{ 0.66872486, 0.06794672, 0.0607708 },
  coords::r3{ 0.84176226, 0.31045712, 0.17150052 },
  coords::r3{ 0.46912322, 0.40808598, 0.39124196 },
  coords::r3{ 0.90971707, 0.7951251, 0.37516915 },
  coords::r3{ 0.67829251, 0.22864326, 0.1908322 }
}, methanolMolecule{
  coords::r3{ -0.321, -0.087, 0.12733 },
  coords::r3{ 1.055, 0.144, 0.21133 },
  coords::r3{ -0.53, -0.921, -0.57767 },
  coords::r3{ -0.85, 0.837, -0.18867 },
  coords::r3{ -0.699, -0.376, 1.12933 },
  coords::r3{ 1.345, 0.403, -0.70167 }
}, vectorOfDoubles{
  0.0636934 ,  0.12765765,  0.75803137,  0.70638564,  0.48240959,
  0.47581971,  0.52241619,  0.08248583,  0.40740382,  0.18346276
}{}

void GetMeanTest::testRandomMolecule() {
  auto mean = ic_util::get_mean(randomCoordinates);
  isCartesianPointNear(mean, coords::r3{ 0.63515582, 0.46322410666666669, 0.35401313333333334 });
}

TEST_F(GetMeanTest, testRandomMolecule) {
  testRandomMolecule();
}

void GetMeanTest::testMethanolMolecule() {
  auto mean = ic_util::get_mean(methanolMolecule);
  isCartesianPointNear(mean, coords::r3{ 0., 0., -3.3333333333551707e-06 });
}

TEST_F(GetMeanTest, testMethanolMolecule) {
  testMethanolMolecule();
}

void GetMeanTest::testVectorOfDoubles() {
  auto mean = ic_util::get_mean(vectorOfDoubles);
  EXPECT_NEAR(mean, 0.380976596,doubleNearThreshold);
}

TEST_F(GetMeanTest, testVectorOfDoubles) {
  testVectorOfDoubles();
}

TEST(ConvertContainerToVector, ArrayToVector) {
  std::array<double, 3u> container{ 1.,2.,3. };
  auto sameAsVector = ic_util::arr_to_vec(container);
  EXPECT_EQ(is_vector<decltype(sameAsVector)>::value, true);
}

