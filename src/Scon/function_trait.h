#ifndef function_trait_header

#define function_trait_header


#include <tuple>
#include <type_traits>
#include <cstddef>

namespace function_trait_detail
{

	template<class F>
	struct function_traits
	{
	private:
		using call_type = function_traits < decltype(&F::operator()) >;
	public:
		using return_type = typename call_type::return_type;

		static const std::size_t arity = call_type::arity - 1;

		template <std::size_t N>
		struct argument
		{
			static_assert(N < arity, "error: invalid parameter index.");
			using type = typename call_type::template argument<N + 1>::type;
		};
	};

	template<class F>
	struct function_traits<F&> : public function_traits < F >
	{ };

	template<class F>
	struct function_traits<F&&> : public function_traits < F >
	{ };

	// function pointer

	template<class R, class... Args>
	struct function_traits<R(*)(Args...)> : public function_traits < R(Args...) >
	{ };

	template<class R, class... Args>
	struct function_traits < R(Args...) >
	{
		using return_type = R;

		static const std::size_t arity = sizeof...(Args);

		template <std::size_t N>
		struct argument
		{
			static_assert(N < arity, "error: invalid parameter index.");
			using type = typename std::tuple_element<N, std::tuple<Args...>>::type;
		};
	};

	template<class C, class R, class... Args>
	struct function_traits<R(C::*)(Args...)>
		: public function_traits < R(C&, Args...) >
	{ };

	template<class C, class R, class... Args>
	struct function_traits<R(C::*)(Args...) const>
		: public function_traits < R(C&, Args...) >
	{ };

	template<class C, class R>
	struct function_traits<R(C::*)> : public function_traits < R(C&) >
	{ };

	// ...

	template<class T>
	using return_type = typename function_traits<T>::return_type;

	template<class T, std::size_t N>
	using argument_type = typename function_traits<T>::template argument<N>::type;

	template<class T, std::size_t N>
	using decayed_argument_type = typename std::decay<argument_type<T, N>>::type;

	template<class T, class U>
	struct has_dot
	{
	private:
		template <class V, class W> static char d(...);
		template <class V, class W> static auto d(V& x, W& y) -> decltype(dot(x, y));
	public:
		static const bool value = sizeof(d<T, U>(std::declval<T&>(), std::declval<U&>())) > sizeof(char);
	};

}


#endif
