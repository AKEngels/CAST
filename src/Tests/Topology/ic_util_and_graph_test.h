#ifdef GOOGLE_MOCK

#ifndef IC_UTIL_AND_GRAPH_TEST_H
#define IC_UTIL_AND_GRAPH_TEST_H

#include<map>
#include<type_traits>
#include<vector>

#include<gtest/gtest.h>
#include <boost/graph/adjacency_list.hpp>
#include"../../graph.h"
#include"../../coords.h"

struct expectedValuesForMap {
  std::string key;
  std::size_t expectedValue;
  std::size_t defaultValue;
  bool isExpectedValue;
};

class KeyOrDefaultTest : public testing::Test, public testing::WithParamInterface<expectedValuesForMap> {
public:
  KeyOrDefaultTest();
  void testKeyOrDefault();
private:
  std::map<std::string, std::size_t> testMap;
};

class TestCreateGraph : public testing::Test {
public:
  TestCreateGraph();

private:
  std::vector<ic_util::Node> atomVector; 
  std::vector<std::pair<int, int>> connectivity;
};

class GetMeanTest : public testing::Test{
public:
  GetMeanTest();
  void testRandomMolecule();
  void testMethanolMolecule();
  void testVectorOfDoubles();
private:
  coords::Representation_3D randomCoordinates;
  coords::Representation_3D methanolMolecule;
  std::vector<double> vectorOfDoubles;
};

class ConvertContainerToVectorTest : public testing::Test {
public:
  ConvertContainerToVectorTest() : referenceVector{ 1.,2.,3. }{}

  void arrayToVector();
  void vectorToVector();
  void listToVector();
  void dequeToVector();
  void setToVector();
  void multisetToVector();

private:
  std::vector<double> referenceVector;
};



class BuildUpGraphTest : public testing::Test {
public:
  using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, ic_util::Node>;
  BuildUpGraphTest();
  Graph const expectedGraph;
  ic_util::Graph<ic_util::Node> const actualGraph;
  void testIfAtomsAreSetRight();
  void testIfEdgesAreSetRight();
private:
    std::vector<ic_util::Node> createTestAtomVector();
    std::vector<std::pair<std::size_t, std::size_t>> createTestEdgesVector();
    Graph makeExpectedGraph();

};

#endif

#endif
