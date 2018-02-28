#if !defined(SCON_TRAITS_HEADER)
#define SCON_TRAITS_HEADER

#include <iterator>
#include <type_traits>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <array>

#include "scon.h"
#include "function_trait.h"

namespace scon
{
  namespace is_container_impl {
    template <typename T>                struct is_container :std::false_type {};
    template <typename T>                struct is_container<T[]>                          :std::true_type{};
    template <typename T, std::size_t N> struct is_container<T[N]>                         :std::true_type{};
    template <typename T, std::size_t N> struct is_container<std::array         <T, N>>    :std::true_type{};
    template <typename... Args>          struct is_container<std::vector        <Args...>> :std::true_type{};
    template <typename... Args>          struct is_container<std::set           <Args...>> :std::true_type{};
    template <typename... Args>          struct is_container<std::map           <Args...>> :std::true_type{};
    template <typename... Args>          struct is_container<std::unordered_set <Args...>> :std::true_type{};
    template <typename... Args>          struct is_container<std::unordered_map <Args...>> :std::true_type{};
  }

  template<typename T> struct is_container {
    static constexpr bool value = is_container_impl::is_container<std::decay_t<T>>::value;
  };

  namespace trait_detail
  {
    //enum class enabled_ty {};
    using enabled_ty = int;
  }

  template<bool B> struct Boolean : std::true_type { };
  template<> struct Boolean<false> : std::false_type {};

  // Conditional

  template <bool If, typename Then, typename Else>
  using Conditional = typename std::conditional<If, Then, Else>::type;

  // And

  template<bool... Args> struct And;

  template<bool A, bool... Args>
  struct And<A, Args...> :
    Boolean<A && And<Args...>::value>{};

  template<bool A> struct And<A> :
    Boolean<A> {};

  // Not

  template<class T> using Not = Boolean<!T::value>;

  // Any

  template <typename... T>
  struct Any : Boolean<false> {};

  template <typename Head, typename... Tail>
  struct Any<Head, Tail...> : 
    Conditional<Head::value, Boolean<true>, Any<Tail...>>{};

  // All

  template <typename... T>
  struct All : Boolean<true> {};

  template <typename Head, typename... Tail>
  struct All<Head, Tail...> : 
    Conditional<Head::value, All<Tail...>, Boolean<false>>{};

  // Enable and Disable

  template<class ... Conditions>
  using EnableIf = typename std::enable_if<All<Conditions...>::value, 
    trait_detail::enabled_ty>::type;

  template<class ... Condition>
  using DisableIf = typename std::enable_if<Not<Any<Condition...>>::value, 
    trait_detail::enabled_ty>::type;

  // Qualification conversion

  template<typename T>
  using unqualified_type = typename std::remove_cv<
    typename std::remove_reference<T>::type
  >::type;

  template<typename T>
  using Noreference = typename std::remove_reference<T>::type;

  template<class T>
  using Decayed = typename std::decay<T>::type;

  // Range Trait Details

  namespace trait_detail
  {

    using std::begin;
    using std::end;
    template <int I> using _ol = std::integral_constant<int, I>*;

    // is_range: is begin(v) and end(v) valid for v of type T

    template<class T>
    struct is_range
    {
    private:
      // helper function declarations usign expression sfinae
      template <class U, _ol<0> = nullptr>
      static std::false_type b(...);
      template <class U, _ol<1> = nullptr>
      static auto b(U &v) -> decltype(begin(v), std::true_type());
      template <class U, _ol<0> = nullptr>
      static std::false_type e(...);
      template <class U, _ol<1> = nullptr>
      static auto e(U &v) -> decltype(end(v), std::true_type());
      // return types
      using b_return = decltype(b<T>(std::declval<T&>()));
      using e_return = decltype(e<T>(std::declval<T&>()));
    public:
      static const bool value = b_return::value && e_return::value;
    };

    template<class T, bool b = is_range<T>::value>
    struct range_begin_iterator_type { };

    template<class T>
    struct range_begin_iterator_type < T, true >
    {
      using type = decltype(begin(std::declval<T&>()));
    };

    template<class T, bool b = is_range<T>::value>
    struct range_end_iterator_type { };

    template<class T>
    struct range_end_iterator_type < T, true >
    {
      using type = decltype(end(std::declval<T&>()));
    };

    // range_value_type: dereferenced begin(v) for v of type T

    template<class T, bool b = is_range<T>::value>
    struct range_value_type { };

    template<class T>
    struct range_value_type < T, true >
    {
      using _drb = decltype(*begin(std::declval<T&>()));
      using type = typename std::remove_reference<_drb>::type;
    };

  }

  // alias for is_range and range_value_type

  template <class T> using is_range = 
    trait_detail::is_range < T >;

  template <class T> using range_value = 
    typename trait_detail::range_value_type<T>::type;

  template <class T> using range_begin_iterator =
    typename trait_detail::range_begin_iterator_type<T>::type;

  template <class T> using range_end_iterator =
    typename trait_detail::range_end_iterator_type<T>::type;

  // alias for function traits

  template<class T>
  using return_type = function_trait_detail::return_type<T>;

  template<class T, std::size_t N>
  using argument_type = function_trait_detail::argument_type<T, N>;

  // Floating point operation helpers

  template <class T> struct float_helper       { using type = double; };
  template <> struct float_helper<float>       { using type = float; };
  template <> struct float_helper<long double> { using type = long double; };

  // Construction and conversion traits

  template<typename T, typename... Args>
  struct is_convertible_from : 
    std::is_constructible<T, Args...> {};

  template<typename T, typename U>
  struct is_convertible_from<T, U> :
    std::is_convertible<U, T> {};

  template<typename T, typename U>
  struct is_explicitly_constructible :
    Boolean<std::is_constructible<T, U>::value && 
    !std::is_convertible<U, T>::value> {};


}

#define SCON_PERFECT_FORWARDING_WRAPPER_CONSTRUCTORS(classname, wrapped_type, wrapped_name) \
template<typename... Args, \
  scon::EnableIf< scon::is_convertible_from< wrapped_type , Args... > >  = 0 > \
classname (Args&&... args) : wrapped_name (std::forward<Args>(args)...) { } \
\
template<typename T,\
  scon::EnableIf< scon::is_explicitly_constructible< wrapped_type , T > >  = 0>\
explicit classname (T&& arg) : wrapped_name (std::forward<T>(arg)) { }\
\
template<class T,\
  scon::EnableIf< std::is_convertible< std::initializer_list<T>, wrapped_type > >  = 0 >\
classname(std::initializer_list<T> list) : wrapped_name(list) { }\
\
template<class T, class ... Args,\
  scon::EnableIf< std::is_constructible< wrapped_type, std::initializer_list<T>, Args...> >  = 0 >\
classname(std::initializer_list<T> list, Args && ... args) : wrapped_name(list, std::forward<Args>(args)...) { }


#endif
