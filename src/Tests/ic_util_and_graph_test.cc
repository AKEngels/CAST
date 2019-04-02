#ifdef GOOGLE_MOCK

#include<array>
#include<list>
#include<set>
#include<array>

#include "ic_util_and_graph_test.h"
#include "../ic_util.h"

#include "../Scon/scon_traits.h"

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
  {0,1},
{ 0,2 },
{ 0,3 },
{ 0,4 },
{ 1,5 }
} {}

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

void ConvertContainerToVectorTest::arrayToVector(){
  std::array<double, 3u> container{ 1.,2.,3. };
  auto sameAsVector = ic_util::arr_to_vec(container);
  EXPECT_EQ(scon::is_vector<decltype(sameAsVector)>::value, true);
  EXPECT_EQ(sameAsVector, referenceVector);
}

TEST_F(ConvertContainerToVectorTest, arrayToVector) {
  arrayToVector();
}

void ConvertContainerToVectorTest::vectorToVector() {
  std::vector<double> container{ 1.,2.,3. };
  auto sameAsVector = ic_util::arr_to_vec(container);
  EXPECT_EQ(scon::is_vector<decltype(sameAsVector)>::value, true);
  EXPECT_EQ(sameAsVector, referenceVector);
}

TEST_F(ConvertContainerToVectorTest, vectorToVector) {
  vectorToVector();
}

void ConvertContainerToVectorTest::listToVector() {
  std::list<double> container{ 1.,2.,3. };
  auto sameAsVector = ic_util::arr_to_vec(container);
  EXPECT_EQ(scon::is_vector<decltype(sameAsVector)>::value, true);
  EXPECT_EQ(sameAsVector, referenceVector);
}

TEST_F(ConvertContainerToVectorTest, listToVector) {
  listToVector();
}

void ConvertContainerToVectorTest::dequeToVector() {
  std::deque<double> container{ 1.,2.,3. };
  auto sameAsVector = ic_util::arr_to_vec(container);
  EXPECT_EQ(scon::is_vector<decltype(sameAsVector)>::value, true);
  EXPECT_EQ(sameAsVector, referenceVector);
}

TEST_F(ConvertContainerToVectorTest, dequeToVector) {
  dequeToVector();
}

void ConvertContainerToVectorTest::setToVector(){
  std::set<double> container{ 1.,2.,3. };
  auto sameAsVector = ic_util::arr_to_vec(container);
  EXPECT_EQ(scon::is_vector<decltype(sameAsVector)>::value, true);
  EXPECT_EQ(sameAsVector, referenceVector);
}

TEST_F(ConvertContainerToVectorTest, setToVector) {
  setToVector();
}

void ConvertContainerToVectorTest::multisetToVector() {
  std::multiset<double> container{ 1.,2.,3. };
  auto sameAsVector = ic_util::arr_to_vec(container);
  EXPECT_EQ(scon::is_container<decltype(sameAsVector)>::value, true);
  EXPECT_EQ(sameAsVector, referenceVector);
}

TEST_F(ConvertContainerToVectorTest, multisetToVector) {
  multisetToVector();
}

TEST(IcUtilityFreeFunctions, Rep3D_to_Mat) {
  coords::Representation_3D cartesianCoordinates{
    coords::r3{ -0.321, -0.087, 0.12733 },
    coords::r3{ 1.055, 0.144, 0.21133 },
    coords::r3{ -0.53, -0.921, -0.57767 },
    coords::r3{ -0.85, 0.837, -0.18867 },
    coords::r3{ -0.699, -0.376, 1.12933 },
    coords::r3{ 1.345, 0.403, -0.70167 }
  };
  auto cartesianMatrix = ic_util::Rep3D_to_Mat(cartesianCoordinates);
  scon::mathmatrix<double> expectedValues { { -0.321, -0.087, 0.12733 },
  { 1.055, 0.144, 0.21133 },
  { -0.53, -0.921, -0.57767 },
  { -0.85, 0.837, -0.18867 },
  { -0.699, -0.376, 1.12933 },
  { 1.345, 0.403, -0.70167 }};
  EXPECT_EQ(cartesianMatrix, expectedValues);
}

TEST(IcUtilityFreeFunctions, matToRep3D) {
  scon::mathmatrix<double> matrixToTransform{ {-0.321,}, {-0.087,}, {0.12733,} , {1.055,}, {0.144,}, {0.21133,}, {-0.53,}, {-0.921,}, {-0.57767,},
  {-0.85, }, {0.837,}, {-0.18867,}, {-0.699,}, {-0.376,}, {1.12933,}, {1.345,}, {0.403,}, {-0.70167,}};
  auto evaluatedValues = ic_util::matToRep3D(matrixToTransform);

  coords::Representation_3D expectedValues{
    coords::r3{ -0.321, -0.087, 0.12733 },
    coords::r3{ 1.055, 0.144, 0.21133 },
    coords::r3{ -0.53, -0.921, -0.57767 },
    coords::r3{ -0.85, 0.837, -0.18867 },
    coords::r3{ -0.699, -0.376, 1.12933 },
    coords::r3{ 1.345, 0.403, -0.70167 }
  };
  
  
  EXPECT_EQ(evaluatedValues, expectedValues);
}

TEST(IcUtilityFreeFunctions, bonds) {
  std::vector<std::string> symbols{ "C", "O", "H", "H", "H", "H" };
  coords::Representation_3D cartesianCoordinates{
    coords::r3{ -0.321, -0.087, 0.12733 },
    coords::r3{ 1.055, 0.144, 0.21133 },
    coords::r3{ -0.53, -0.921, -0.57767 },
    coords::r3{ -0.85, 0.837, -0.18867 },
    coords::r3{ -0.699, -0.376, 1.12933 },
    coords::r3{ 1.345, 0.403, -0.70167 }
  };
  std::vector<std::pair<std::size_t, std::size_t>> expectedBonds{
    { 0u, 1u }, { 0u, 2u }, { 0u, 3u }, { 0u, 4u }, { 1u, 5u }
  };

  auto bonds = ic_util::bonds(symbols, cartesianCoordinates);
  EXPECT_EQ(bonds, expectedBonds);

  coords::Representation_3D twoMethanolMolecules{ coords::r3{ -6.053, -0.324, -0.108 },
    coords::r3{ -4.677, -0.093, -0.024 },
    coords::r3{ -6.262, -1.158, -0.813 },
    coords::r3{ -6.582, 0.600, -0.424 },
    coords::r3{ -6.431, -0.613, 0.894 },
    coords::r3{ -4.387, 0.166, -0.937 },
    coords::r3{ -6.146, 3.587, -0.024 },
    coords::r3{ -4.755, 3.671, -0.133 },
    coords::r3{ -6.427, 2.922, 0.821 },
    coords::r3{ -6.587, 3.223, -0.978 },
    coords::r3{ -6.552, 4.599, 0.179 },
    coords::r3{ -4.441, 2.753, -0.339 } };

  std::vector<std::string> symbolsForTwo{ "C", "O", "H", "H", "H", "H", "C", "O", "H", "H", "H", "H" };

  auto bondsForTwo = ic_util::bonds(symbolsForTwo, twoMethanolMolecules);

  std::vector<std::pair<std::size_t, std::size_t>> expectedBondsForTwo{
    { 0u, 1u },{ 0u, 2u },{ 0u, 3u },{ 0u, 4u },{ 1u, 5u }, { 6u, 7u },{ 6u, 8u }, {6u, 9u}, { 6u, 10u}, { 7u, 11u}
  };

  EXPECT_EQ(bondsForTwo, expectedBondsForTwo);

}

TEST(IcUtilityFreeFunctions, flatten_c3_vec) {
  coords::Representation_3D cartesianCoordinates{
    coords::r3{ -0.321, -0.087, 0.12733 },
    coords::r3{ 1.055, 0.144, 0.21133 },
    coords::r3{ -0.53, -0.921, -0.57767 },
    coords::r3{ -0.85, 0.837, -0.18867 },
    coords::r3{ -0.699, -0.376, 1.12933 },
    coords::r3{ 1.345, 0.403, -0.70167 }
  };
  std::vector<double> expectedValue{ -0.321, -0.087, 0.12733 , 1.055, 0.144, 0.21133, -0.53, -0.921, -0.57767,
     -0.85, 0.837, -0.18867, -0.699, -0.376, 1.12933, 1.345, 0.403, -0.70167 };
  auto evaluatedValues = ic_util::flatten_c3_vec(cartesianCoordinates);

  EXPECT_EQ(evaluatedValues, expectedValue);

}

BuildUpGraphTest::BuildUpGraphTest() : expectedGraph{ makeExpectedGraph() }, actualGraph{ic_util::make_graph(createTestEdgesVector(), createTestAtomVector())}{}

BuildUpGraphTest::Graph BuildUpGraphTest::makeExpectedGraph() {

  BuildUpGraphTest::Graph graph;

  auto firstAtom = boost::add_vertex(graph);
  auto secondAtom = boost::add_vertex(graph);
  auto thirdAtom = boost::add_vertex(graph);
  auto fourthAtom = boost::add_vertex(graph);
  auto fifthAtom = boost::add_vertex(graph);
  auto sixthAtom = boost::add_vertex(graph);

  graph[firstAtom] = ic_util::Node{ 1, "C", "C" };
  graph[secondAtom] = ic_util::Node{ 2, "O", "O" };
  graph[thirdAtom] = ic_util::Node{ 3, "H", "H" };
  graph[fourthAtom] = ic_util::Node{ 4, "H", "H" };
  graph[fifthAtom] = ic_util::Node{ 5, "H", "H" };
  graph[sixthAtom] = ic_util::Node{ 6, "H", "H" };

  boost::add_edge(firstAtom, secondAtom, graph);
  boost::add_edge(firstAtom, thirdAtom, graph);
  boost::add_edge(firstAtom, fourthAtom, graph);
  boost::add_edge(firstAtom, fifthAtom, graph);
  boost::add_edge(secondAtom, sixthAtom, graph);

  return graph;
}

std::vector<ic_util::Node> BuildUpGraphTest::createTestAtomVector() {
  return { ic_util::Node{ 1, "C", "C" }, ic_util::Node{ 2, "O", "O" }, ic_util::Node{ 3, "H", "H" }, ic_util::Node{ 4, "H", "H" }, ic_util::Node{ 5, "H", "H" }, ic_util::Node{ 6, "H", "H" } };
}

std::vector<std::pair<std::size_t, std::size_t>> BuildUpGraphTest::createTestEdgesVector() {
  return { { 0u, 1u },{ 0u, 2u },{ 0u, 3u },{ 0u, 4u },{ 1u, 5u } };
}

void BuildUpGraphTest::testIfAtomsAreSetRight(){
  auto expectedIterators = boost::vertices(expectedGraph);
  auto actualIterators = boost::vertices(actualGraph);

  for (; expectedIterators.first != expectedIterators.second || actualIterators.first != actualIterators.second; ++expectedIterators.first, ++actualIterators.first) {
    EXPECT_EQ(expectedGraph[*expectedIterators.first].atom_serial, actualGraph[*actualIterators.first].atom_serial);
    EXPECT_EQ(expectedGraph[*expectedIterators.first].atom_name, actualGraph[*actualIterators.first].atom_name);
    EXPECT_EQ(expectedGraph[*expectedIterators.first].element, actualGraph[*actualIterators.first].element);
  }

}

TEST_F(BuildUpGraphTest, testIfAtomsAreSetRight) {
  testIfAtomsAreSetRight();
}


void BuildUpGraphTest::testIfEdgesAreSetRight() {
  auto expectedIterators = boost::edges(expectedGraph);
  auto actualIterators = boost::edges(actualGraph);

  for (; expectedIterators.first != expectedIterators.second || actualIterators.first != actualIterators.second; ++expectedIterators.first, ++actualIterators.first) {
    auto expectedSource = boost::source(*expectedIterators.first,expectedGraph);
    auto expectedTarget = boost::target(*expectedIterators.first, expectedGraph);
    auto actualSource = boost::source(*expectedIterators.first, actualGraph);
    auto actualTarget = boost::target(*expectedIterators.first, actualGraph);

    EXPECT_EQ(expectedSource, actualSource);
    EXPECT_EQ(expectedTarget, actualTarget);
  }
}

TEST_F(BuildUpGraphTest, testIfEdgesAreSetRight) {
  testIfEdgesAreSetRight();
}

#endif
