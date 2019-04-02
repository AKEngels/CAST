#ifdef GOOGLE_MOCK
#include "gtest/gtest.h"
#include <type_traits>
#include <array>
#include <vector>
#include <functional>
#include "../Scon/scon_traits.h"

TEST(SconTraitHelpers, BooleanHelper)
{
  EXPECT_TRUE(scon::Boolean<true>::value == std::true_type::value);
  EXPECT_TRUE(scon::Boolean<true>::value == true);
  EXPECT_TRUE(scon::Boolean<false>::value == std::false_type::value);
  EXPECT_TRUE(scon::Boolean<false>::value == false);
}

TEST(SconTraitHelpers, ConditionalHelper)
{
  EXPECT_TRUE((std::is_same<scon::Conditional<true, char, int>, char>::value));
  EXPECT_TRUE((std::is_same<scon::Conditional<false, char, int>, int>::value));
}

TEST(SconTraitHelpers, AndHelper)
{
  // 1x false
  EXPECT_FALSE((scon::And<false>::value));
  // 1x true
  EXPECT_TRUE((scon::And<true>::value));
  // 2x false
  EXPECT_FALSE((scon::And<false, false>::value));
  // 1x true, 1x false
  EXPECT_FALSE((scon::And<true, false>::value));
  EXPECT_FALSE((scon::And<false, true>::value));
  // 2x true
  EXPECT_TRUE((scon::And<true, true>::value));
  // 3x false
  EXPECT_FALSE((scon::And<false, false, false>::value));
  // 1x true, 2x false
  EXPECT_FALSE((scon::And<true, false, false>::value));
  EXPECT_FALSE((scon::And<false, true, false>::value));
  EXPECT_FALSE((scon::And<false, false, true>::value));
  // 2x true, 1x false
  EXPECT_FALSE((scon::And<true, true, false>::value));
  EXPECT_FALSE((scon::And<true, false, true>::value));
  EXPECT_FALSE((scon::And<false, true, true>::value));
  // 3x true
  EXPECT_TRUE((scon::And<true, true, true>::value));

  // NOTE: Cannot accomplish "full" test
  // since one cannot test any arbitrary combination

}

TEST(SconTraitHelpers, NotHelper)
{
  EXPECT_FALSE((scon::Not<std::true_type>::value));
  EXPECT_TRUE((scon::Not<std::false_type>::value));
}

TEST(SconTraitHelpers, AnyHelper)
{
  // None true -> false
  EXPECT_FALSE((scon::Any<std::false_type>::value));
  EXPECT_FALSE((scon::Any<std::false_type, std::false_type>::value));
  EXPECT_FALSE((scon::Any<std::false_type, std::false_type, std::false_type>::value));

  // Either 1, 2 or 3 true -> true
  EXPECT_TRUE((scon::Any<std::true_type>::value));

  EXPECT_TRUE((scon::Any<std::true_type, std::false_type>::value));
  EXPECT_TRUE((scon::Any<std::false_type, std::true_type>::value));
  EXPECT_TRUE((scon::Any<std::true_type, std::true_type>::value));

  EXPECT_TRUE((scon::Any<std::true_type, std::false_type, std::false_type>::value));
  EXPECT_TRUE((scon::Any<std::false_type, std::true_type, std::false_type>::value));
  EXPECT_TRUE((scon::Any<std::false_type, std::false_type, std::true_type>::value));

  EXPECT_TRUE((scon::Any<std::true_type, std::true_type, std::false_type>::value));
  EXPECT_TRUE((scon::Any<std::false_type, std::true_type, std::true_type>::value));
  EXPECT_TRUE((scon::Any<std::true_type, std::false_type, std::true_type>::value));

  EXPECT_TRUE((scon::Any<std::true_type, std::true_type, std::true_type>::value));
}

TEST(SconTraitHelpers, AllHelper)
{
  // All True, -> true
  EXPECT_TRUE((scon::All<std::true_type>::value));
  EXPECT_TRUE((scon::All<std::true_type, std::true_type>::value));
  EXPECT_TRUE((scon::All<std::true_type, std::true_type, std::true_type>::value));

  // One or more false -> false
  EXPECT_FALSE((scon::All<std::false_type>::value));

  EXPECT_FALSE((scon::All<std::true_type, std::false_type>::value));
  EXPECT_FALSE((scon::All<std::false_type, std::true_type>::value));
  EXPECT_FALSE((scon::All<std::false_type, std::false_type>::value));

  EXPECT_FALSE((scon::All<std::true_type, std::false_type, std::false_type>::value));
  EXPECT_FALSE((scon::All<std::false_type, std::true_type, std::false_type>::value));
  EXPECT_FALSE((scon::All<std::false_type, std::false_type, std::true_type>::value));

  EXPECT_FALSE((scon::All<std::true_type, std::true_type, std::false_type>::value));
  EXPECT_FALSE((scon::All<std::false_type, std::true_type, std::true_type>::value));
  EXPECT_FALSE((scon::All<std::true_type, std::false_type, std::true_type>::value));

  EXPECT_FALSE((scon::All<std::false_type, std::false_type, std::false_type>::value));
}

TEST(SconTraitHelpers, UnqualifiedTypeHelper)
{
  // All True, -> true
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const>, int>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile>, int>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile>, int>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int&>, int>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const &>, int>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile &>, int>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile &>, int>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int&&>, int>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const &&>, int>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile &&>, int>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile &&>, int>::value));

  EXPECT_TRUE((std::is_same<scon::unqualified_type<int * const>, int *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int * volatile>, int *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int * const volatile>, int *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int * const &>, int *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int * volatile &>, int *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int * const volatile &>, int *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int * const &&>, int *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int * volatile &&>, int *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int * const volatile && >, int *>::value));

  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const * const>, int const *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const * volatile>, int const *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const * const volatile>, int const *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const * const &>, int const *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const * volatile &>, int const *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const * const volatile &>, int const *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const * const &&>, int const *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const * volatile &&>, int const *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const * const volatile &&>, int const *>::value));

  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile * const>, int volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile * volatile>, int volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile * const volatile>, int volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile * const &>, int volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile * volatile &>, int volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile * const volatile &>, int volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile * const &&>, int volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile * volatile &&>, int volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int volatile * const volatile &&>, int volatile *>::value));

  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile * const>, int const volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile * volatile>, int const volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile * const volatile>, int const volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile * const &>, int const volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile * volatile &>, int const volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile * const volatile &>, int const volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile * const &&>, int const volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile * volatile &&>, int const volatile *>::value));
  EXPECT_TRUE((std::is_same<scon::unqualified_type<int const volatile * const volatile &&>, int const volatile *>::value));

}

TEST(SconTraitHelpers, FloatHelper)
{
  EXPECT_TRUE((std::is_same<scon::float_helper<float>::type, float>::value));
  EXPECT_TRUE((std::is_same<scon::float_helper<double>::type, double>::value));
  EXPECT_TRUE((std::is_same<scon::float_helper<long double>::type, long double>::value));
  // every other type ...
  EXPECT_TRUE((std::is_same<scon::float_helper<int>::type, double>::value));
  EXPECT_TRUE((std::is_same<scon::float_helper<unsigned>::type, double>::value));
  EXPECT_TRUE((std::is_same<scon::float_helper<char>::type, double>::value));
  EXPECT_TRUE((std::is_same<scon::float_helper<long>::type, double>::value));
  EXPECT_TRUE((std::is_same<scon::float_helper<long long>::type, double>::value));
  EXPECT_TRUE((std::is_same<scon::float_helper<std::vector<int>>::type, double>::value));
  // extendible to the end of days...
}

// No test for Noreference and Decayed as they're just wrappers for std:: stuff

namespace some_namespace
{
  struct some_test_value {};
  struct some_other_value {};
  struct some_test_iterator { some_test_value & operator*(); };
  struct some_test_iterator_2 { some_other_value & operator*(); };
  struct st_no_range_1 {};
  struct st_no_range_2 { some_test_iterator begin(); };
  struct st_range_1 { some_test_iterator begin() const; some_test_iterator end() const; };
  struct st_range_2 { };
  some_test_iterator_2 begin(st_range_2 const &);
  some_test_iterator_2 end(st_range_2 const &);
}

TEST(SconTraits, IsRangeTrait)
{
  // begin and end members -> indirect via std::begin / std::end ?
  EXPECT_TRUE((scon::is_range<some_namespace::st_range_1>::value));
  // begin and end free functions -> direct
  EXPECT_TRUE((scon::is_range<some_namespace::st_range_2>::value));
  // ranges anyway
  EXPECT_TRUE((scon::is_range<std::vector<int>>::value));
  EXPECT_TRUE((scon::is_range<std::array<int, 3u>>::value));
  // non-ranges
  EXPECT_FALSE((scon::is_range<int>::value));
  EXPECT_FALSE((scon::is_range<float>::value));
  EXPECT_FALSE((scon::is_range<some_namespace::st_no_range_1>::value));
  EXPECT_FALSE((scon::is_range<some_namespace::st_no_range_2>::value));
}

TEST(SconTraits, RangeValueTrait)
{
  // value type of vector<int> should be int, shouldn't it?
  EXPECT_TRUE((std::is_same<scon::range_value<std::vector<int>>, int>::value));
  EXPECT_TRUE((std::is_same<scon::range_value<std::vector<int> const>, int const>::value));
  EXPECT_TRUE((std::is_same<scon::range_value<std::vector<int> const &>, int const>::value));
  EXPECT_TRUE((std::is_same<scon::range_value<some_namespace::st_range_1>, some_namespace::some_test_value>::value));
  EXPECT_TRUE((std::is_same<scon::range_value<some_namespace::st_range_1&>, some_namespace::some_test_value>::value));
  EXPECT_TRUE((std::is_same<scon::range_value<some_namespace::st_range_1 const>, some_namespace::some_test_value>::value));
  EXPECT_TRUE((std::is_same<scon::range_value<some_namespace::st_range_1 const &>, some_namespace::some_test_value>::value));
  EXPECT_TRUE((std::is_same<scon::range_value<some_namespace::st_range_2>, some_namespace::some_other_value>::value));
  EXPECT_TRUE((std::is_same<scon::range_value<some_namespace::st_range_2 &>, some_namespace::some_other_value>::value));
  EXPECT_TRUE((std::is_same<scon::range_value<some_namespace::st_range_2 const>, some_namespace::some_other_value>::value));
  EXPECT_TRUE((std::is_same<scon::range_value<some_namespace::st_range_2 const &>, some_namespace::some_other_value>::value));

}

namespace
{
	bool some_bool_ret_test() { 
		bool somebool{ false };   // doesn't matter if true or false
		return somebool; 
	}
	float somefloat;            // just a float to reference to
	float & some_float_ref_ret_test() { 
		return somefloat; 
	}
  void some_void_ret_int_int(int, int) { }
  struct some_functor {    char operator()() const; };
  int some_int_returning_function() { return 1; }
}

TEST(SconTraits, ReturnTypeTrait)
{
  // general syntax

  EXPECT_TRUE((std::is_same<scon::return_type<int(void)>, int>::value));
  EXPECT_TRUE((std::is_same<scon::return_type<double(float, int&)>, double>::value));

  // lambda expressions
  auto char_returning_lambda = []() -> char { return ' '; };
  auto int_and_float_ref_taking_bool_returning_lambda = [](int, float&) -> bool { return true; };
  EXPECT_TRUE((std::is_same<scon::return_type<decltype(char_returning_lambda)>, char>::value));
  EXPECT_TRUE((std::is_same<scon::return_type<decltype(int_and_float_ref_taking_bool_returning_lambda)>, bool>::value));

  // std::function objects
  auto int_returning_std_function = std::function<int()>{ &some_int_returning_function };
  EXPECT_TRUE((std::is_same<scon::return_type<decltype(int_returning_std_function)>, int>::value));
  auto void_returning_int_int_std_function = std::function<void(int,int)>{ &some_void_ret_int_int };
  EXPECT_TRUE((std::is_same<scon::return_type<decltype(int_returning_std_function)>, int>::value));

  // function object
  auto char_returning_functor = some_functor{};
  EXPECT_TRUE((std::is_same<scon::return_type<decltype(char_returning_functor)>, char>::value));

  // Note; bind expressions not possible using scon::return_type:
  // empty operator() required which is not necessarily there

  // function pointers
  EXPECT_TRUE((std::is_same<scon::return_type<decltype(&some_bool_ret_test)>, bool>::value));
  EXPECT_TRUE((std::is_same<scon::return_type<decltype(&some_float_ref_ret_test)>, float&>::value));

  // member function pointer
  int& (std::vector<int>::* fpn) (std::vector<int>::size_type) = &std::vector<int>::operator[];
  EXPECT_TRUE((std::is_same<scon::return_type<decltype(fpn)>, int&>::value));

}

TEST(SconTraits, ArgumentTypeTrait)
{
  // general syntax

  EXPECT_TRUE((std::is_same<scon::argument_type<double(float, int&), 0u>, float>::value));
  EXPECT_TRUE((std::is_same<scon::argument_type<double(float, int&), 1u>, int&>::value));

  auto int_and_floatref_taking_lambda = [](int , float & ) -> void {};
  EXPECT_TRUE((std::is_same<scon::argument_type<decltype(int_and_floatref_taking_lambda), 0u>, int>::value));
  EXPECT_TRUE((std::is_same<scon::argument_type<decltype(int_and_floatref_taking_lambda), 1u>, float&>::value));

  int& (std::vector<int>::* fpn) (std::vector<int>::size_type) = &std::vector<int>::operator[];
  EXPECT_TRUE((std::is_same<scon::argument_type<decltype(fpn), 1u>, std::vector<int>::size_type>::value));
}

namespace
{

  struct some_A
  {
    some_A(std::vector<int>);
    explicit some_A(unsigned);
    some_A(unsigned, float);
  };
  struct convertible_to_some_A
  { operator some_A(); };
  struct explicitly_convertible_to_some_A
  { explicit operator some_A(); };

}

TEST(SconTraits, ConvertibleFromtTrait)
{
  // multi argument construction -> "conversion" possible
  EXPECT_TRUE((scon::is_convertible_from<some_A, unsigned, float>::value));
  // single argument construction with non explicit constructor -> possible
  EXPECT_TRUE((scon::is_convertible_from<some_A, std::vector<int>>::value));
  // constructor is explicit, no implicit conversion unsigned -> some_A
  EXPECT_FALSE((scon::is_convertible_from<some_A, unsigned>::value));
  // conversion operator -> implicitly convertible
  EXPECT_TRUE((scon::is_convertible_from<some_A, convertible_to_some_A>::value));
  // explicit conversion operator -> not convertible
  EXPECT_FALSE((scon::is_convertible_from<some_A, explicitly_convertible_to_some_A>::value));
}

TEST(SconTraits, UnsignedShouldBeExplicitlyConstructible) {
  EXPECT_TRUE((scon::is_explicitly_constructible<some_A, unsigned>::value));
}

//Not sure if this test gets the intention right.
TEST(SconTraits, ExplicitlyConvertibleShouldBeExplicitlyConstructible) {
  EXPECT_TRUE((scon::is_explicitly_constructible<some_A, explicitly_convertible_to_some_A>::value));
}

TEST(SconTraits, StdVectorValueShouldNotBeExplicitlyConstructible) {
  EXPECT_FALSE((scon::is_explicitly_constructible<some_A, std::vector<int>>::value));
}

TEST(SconTraits, ConvertibleShouldNotBeExplicitlyConstructible) {
  EXPECT_FALSE((scon::is_explicitly_constructible<some_A, convertible_to_some_A>::value));
}

#endif
