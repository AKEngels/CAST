/**
CAST 3
Purpose: Tests stuff for helperfunctions

@author Susanne Sauer
@version 1.0
*/

#ifdef GOOGLE_MOCK
#include "../coords_io.h"
#include "gtest/gtest.h"

TEST(helperfuncs, test_system_mass)
{
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
	coords::Coordinates coords(ci->read("test_files/butanol.arc"));

	auto mass = sys_mass(coords);

	ASSERT_EQ(mass, 74.1216);
}

// energy printing functions are not tested yet

TEST(helperfuncs, test_split)
{
  std::string teststring = "this is a string";

  auto stringvec = split(teststring, ' ');
  std::vector<std::string> result = { "this", "is", "a", "string" };

  ASSERT_EQ(stringvec, result);
}

TEST(helperfuncs, test_split_remove)
{
  std::string teststring = "this is  a string";

  auto stringvec = split(teststring, ' ', true);
  std::vector<std::string> result = { "this", "is", "a", "string" };

  ASSERT_EQ(stringvec, result);
}

TEST(helperfuncs, test_split_no_remove)
{
  std::string teststring = "this is  a string";

  auto stringvec = split(teststring, ' ', false);
  std::vector<std::string> result = { "this", "is", "", "a", "string" };

  ASSERT_EQ(stringvec, result);
}

TEST(helperfuncs, test_dist)
{
  coords::Cartesian_Point a;
  a.x() = 1.4;
  a.y() = 3.2;
  a.z() = 5.6;

  coords::Cartesian_Point b;
  b.x() = 2.4;
  b.y() = 4.2;
  b.z() = 6.6;

  auto distance = dist(a, b);
  ASSERT_EQ(distance, std::sqrt(3));
}

TEST(helperfuncs, test_isin_numbers)
{
  std::vector<double> testvec = { 1.2, 3.4, 81.5, 101.7 };
  bool result = is_in(1.2, testvec);
  ASSERT_TRUE(result);

  result = is_in(5.2, testvec);
  ASSERT_FALSE(result);
}

TEST(helperfuncs, test_isin_strings)
{
  std::vector<std::string> testvec = {"these", "are", "some", "stupid", "strings"};
  bool result = is_in("stupid", testvec);
  ASSERT_TRUE(result);

  result = is_in("bullshit", testvec);
  ASSERT_FALSE(result);
}

TEST(helperfuncs, test_find_index)
{
  std::vector<double> testvec = { 1.2, 3.4, 81.5, 101.7 };
  auto result = find_index(1.2, testvec);
  ASSERT_EQ(result, 0);

  std::vector<std::string> testvec2 = { "these", "are", "some", "stupid", "strings" };
  std::string s = "stupid";
  result = find_index(s, testvec2);
  ASSERT_EQ(result, 3);
}

TEST(helperfuncs, test_find_index_not_in)
{
  std::vector<double> testvec = { 1.2, 3.4, 81.5, 101.7 };
  auto result = find_index(7.2, testvec);
  ASSERT_EQ(result, std::numeric_limits<int>::max());
}

TEST(helperfuncs, test_is_number)
{
  bool result = check_if_number("10");
  ASSERT_TRUE(result);

  result = check_if_number("5.2");
  ASSERT_TRUE(result);
}

TEST(helperfuncs, test_is_number_with_space)
{
	bool result = check_if_number(" 10");
	ASSERT_TRUE(result);

	result = check_if_number("5.2 ");
	ASSERT_TRUE(result);
}

TEST(helperfuncs, test_is_no_number)
{
  bool result = check_if_number("blabla");
  ASSERT_FALSE(result);
}

TEST(helperfuncs, test_is_range)
{
	bool result = check_if_number("1 - 3");
	ASSERT_FALSE(result);

	result = check_if_number("1-3");
	ASSERT_FALSE(result);
}

TEST(helperfuncs, test_is_number_with_sign)
{
  bool result = check_if_number("+10");
  ASSERT_TRUE(result);

  result = check_if_number("-5.2");
  ASSERT_TRUE(result);
}

TEST(helperfuncs, test_is_number_scientific)
{
  bool result = check_if_number("7.3e-5");
  ASSERT_TRUE(result);

  result = check_if_number("-5.2E+3");
	ASSERT_TRUE(result);
}

TEST(helperfuncs, test_file_exists)
{
  bool result = file_exists("test_files/butanol.arc");
	ASSERT_TRUE(result);
}

TEST(helperfuncs, test_file_exists_not)
{
  bool result = file_exists("no_file.dat");
  ASSERT_FALSE(result);
}

TEST(helperfuncs, test_last_line)
{
  std::ifstream in_file("test_files/butanol.arc");
  auto line = last_line(in_file);
  ASSERT_EQ(line, "15  H      1.559120    1.699330   -0.626680    97    12");
}

TEST(helperfuncs, test_file_is_empty)
{
  std::string s = "test_files/empty.txt";
  bool result = file_is_empty(s);
	ASSERT_TRUE(result);
}

TEST(helperfuncs, test_file_is_not_empty)
{
  std::string s = "test_files/butanol.arc";
  bool result = file_is_empty(s);
  ASSERT_FALSE(result);
}

TEST(helperfuncs, test_not_existent_file_is_empty_returns_true)
{
  std::string s = "empty.txt";
  bool result = file_is_empty(s);
	ASSERT_TRUE(result);
}

TEST(helperfuncs, test_add_vectors)
{
  std::vector<int> vec1 = { 5,8,9,10,11,12,13,14 };
  std::vector<int> vec2 = { 0, 1, 2, 3, 4, 6, 7 };

  auto sum = add_vectors(vec1, vec2);
  std::vector<int> result = { 5,8,9,10,11,12,13,14 ,0, 1, 2, 3, 4, 6, 7 };

  ASSERT_EQ(sum, result);
}

TEST(helperfuncs, test_add_and_sort_vectors)
{
  std::vector<int> vec1 = { 5,8,9,10,11,12,13,14 };
  std::vector<int> vec2 = { 0, 1, 2, 3, 4, 6, 7 };

  auto sum = add_vectors(vec1, vec2, true);
  std::vector<int> result = { 0, 1, 2, 3, 4,5, 6, 7 ,8,9,10,11,12,13,14 };

  ASSERT_EQ(sum, result);
}

TEST(helperfuncs, test_double_element_numbers)
{
  std::vector<int> vec = { 5,8,9,10,11,5,12,13,14 };
  bool result = double_element(vec);
	ASSERT_TRUE(result);

  vec = { 8,9,10,11,5,12,13,14 };
  result = double_element(vec);
  ASSERT_FALSE(result);
}

TEST(helperfuncs, test_double_element_strings)
{
  std::vector<std::string> vec = { "these", "are", "some", "stupid", "stupid", "strings" };
  bool result = double_element(vec);
	ASSERT_TRUE(result);

  vec = { "these", "are", "some", "stupid", "strings" };
  result = double_element(vec);
  ASSERT_FALSE(result);
}

TEST(helperfuncs, test_range)
{
  auto range1 = range(10);
  std::vector<int> result = { 0,1,2,3,4,5,6,7,8,9 };
  ASSERT_EQ(range1, result);

  auto range2 = range(10, 5);
  result = {5,6,7,8,9 };
  ASSERT_EQ(range2, result);

  auto range3 = range(11, 5, 2);
  result = { 5,7,9 };
  ASSERT_EQ(range3, result);

  auto range4 = range(5, 11, -2);
  result = { 11,9,7 };
  ASSERT_EQ(range4, result);
}

TEST(helperfuncs, test_count_element)
{
  std::vector<int> testvec = { 0,1,2,7,3,4,5,6,7,8,9,8,7 };

  int count2 = count_element(2, testvec);
  ASSERT_EQ(count2, 1);

  int count7 = count_element(7, testvec);
  ASSERT_EQ(count7, 3);

  int count8 = count_element(8, testvec);
  ASSERT_EQ(count8, 2);

  int count42 = count_element(42, testvec);
  ASSERT_EQ(count42, 0);
}

TEST(helperfuncs, test_count_element_string)
{
  std::vector<std::string> testvec = { "these", "are", "some", "stupid", "stupid", "strings" };

  int count_stupid = count_element("stupid", testvec);
  ASSERT_EQ(count_stupid, 2);

  int count_are = count_element("are", testvec);
  ASSERT_EQ(count_are, 1);

  int count_bullshit = count_element("bullshit", testvec);
  ASSERT_EQ(count_bullshit, 0);
}

TEST(helperfuncs, test_is_smaller_than_really_smaller)
{
	bool result = is_smaller_than(5.3, 10.22);
	ASSERT_TRUE(result);
}

TEST(helperfuncs, test_is_smaller_than_bigger)
{
	bool result = is_smaller_than(10.22, 5.3);
	ASSERT_FALSE(result);
}

TEST(helperfuncs, test_is_smaller_than_just_smaller)
{
	bool result = is_smaller_than(0.3, 0.30000000011);
	ASSERT_TRUE(result);
}

TEST(helperfuncs, test_is_smaller_than_smaller_but_inside_precision)
{
	bool result = is_smaller_than(0.3, 0.30000000009);
	ASSERT_FALSE(result);
}

TEST(helperfuncs, test_is_smaller_userdefined_precision)
{
	bool result = is_smaller_than(0.3, 0.302, 0.001);
	ASSERT_TRUE(result);

	result = is_smaller_than(0.3, 0.3009, 0.001);
	ASSERT_FALSE(result);
}

TEST(helperfuncs, test_existing_atomtype)
{
	int number_of_bonds = get_ideal_bond_number_from_parameterfile(80);
	ASSERT_EQ(number_of_bonds, 4);
}

TEST(helperfuncs, test_non_existing_atomtype)
{
	ASSERT_THROW(get_ideal_bond_number_from_parameterfile(4242), std::runtime_error);
}

TEST(helperfuncs, test_non_existing_parameterfile)
{
	ASSERT_THROW(get_ideal_bond_number_from_parameterfile(80, "test_files/oplsaa_mod2.prm"), std::runtime_error);
}


#endif
