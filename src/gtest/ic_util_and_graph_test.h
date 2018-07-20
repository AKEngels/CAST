#ifdef GOOGLE_MOCK

#ifndef IC_UTIL_AND_GRAPH_TEST_H
#define IC_UTIL_AND_GRAPH_TEST_H

#include<map>

#include<gtest/gtest.h>
#include"../graph.h"
#include"../coords.h"
#include<type_traits>
#include<vector>

namespace is_vector_imp {
  template<typename T> struct is_vector : std::false_type {};
  template<typename T> struct is_vector<std::vector<T>> : std::true_type {};
}

template<typename T>
struct is_vector {
  static constexpr bool value = is_vector_imp::is_vector<std::decay_t<T>>::value;
};

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

class ConvertArrayToVectorTest : public testing::Test {
public:

};

#endif

#endif