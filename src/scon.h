/* SCON LIB CENTRAL HEADER */

#if !defined(SCON_HEADER)
#define SCON_HEADER

#if defined(__GNUG__)

  #if !defined(NDEBUG) && !defined(SCON_DEBUG)
    #define SCON_DEBUG 1
  #elif defined(DEBUG)
    #define SCON_DEBUG 2
  #endif

  #if defined(__i386__) || defined(__i486__) || defined(__i586__) || defined(__i686__)
    #define SCON_x86
  #else
    #define SCON_x64
  #endif

  #ifndef GCC_VERSION
    #define GCC_VERSION (__GNUC__ * 10000 \
      + __GNUC_MINOR__ * 100 \
      + __GNUC_PATCHLEVEL__)
  #endif

  #define SCON_CC11_RVALUE_INIT
  #define SCON_CC11_EXTERN_TEMPLATES

  #if GCC_VERSION >= 40300
    #define SCON_CC11_RVALUE_REF
    #define SCON_CC11_DECLTYPE
    #define SCON_CC11_STATIC_ASSERT
    #define SCON_CC11_F_TEMPLATE_DEFAULT
    #define SCON_CC11_VARIADIC_TEMPLATES
    #define SCON_CC11_USING_DECLARATIONS
    #define SCON_CC11_RIGHT_ANGLE_BRACKETS
  #endif

  #if GCC_VERSION >= 40400
    #define SCON_CC11_INITIALIZER_LISTS
    #define SCON_CC11_AUTO_VARIABLES
    #define SCON_CC11_AUTO_DECLARATOR
    #define SCON_CC11_F_LATE_DECLARATOR
    #define SCON_CC11_SFINAE_EXPRESSIONS
    #define SCON_CC11_STRONGLY_TYPED_ENUMS
    #define SCON_CC11_UNICODE_CHAR
    #define SCON_CC11_DEF_DEL_FUNCTIONS
    #define SCON_CC11_SIZEOF_STATIC
    #define SCON_CC11_INLINE_NAMESPACES
    #define SCON_CC11_CHRONO
    #define SCON_CC11_SFINAE_TRAITS
    #define SCON_CC11_ARRAY
  #endif

  #if GCC_VERSION >= 40500
    #define SCON_CC11_LAMBDA
    #define SCON_CC11_EXPLICIT_CONVERSION_OPERATORS
    #define SCON_CC11_UNICODE_LITERALS
    #define SCON_CC11_RAW_STRING_LITERALS
    #define SCON_CC11_UNIVERSAL_CHAR_LITERALS
    #define SCON_CC11_LOCAL_UNNAMED_TEMPLATE_ARGUMENTS
    #define SCON_CC11_STANDARD_LAYOUT_TYPES
  #endif

  #if GCC_VERSION >= 40600
    #define SCON_CC11_NULLPTR
    #define SCON_CC11_ENUM_FORWARD_DECL
    #define SCON_CC11_GENRALIZED_CONSTEXPR
    #define SCON_CC11_RANGE_BASED_FOR
    #define SCON_CC11_MOVE
    #define SCON_CC11_NOEXCEPT
    #define SCON_CC11_UNRESTRICTED_UNIONS
  #endif

  #if GCC_VERSION >= 40700
    #define SCON_CC11_NONSTATIC_MEMBER_INIT
    #define SCON_CC11_TEMPLATE_ALIAS
    #define SCON_CC11_DELEGATING_CONSTRUCTORS
    #define SCON_CC11_USER_LITERALS
    #define SCON_CC11_EXTENDED_FRIEND
    #define SCON_CC11_EXPLICIT_VIRTUAL_OVERRRIDE
    #define SCON_CC11_OVERRIDE_FINAL
    #define SCON_CC11_RANDOM
    #define SCON_CC11_THREAD
  #endif

  #if GCC_VERSION >= 40800
    #define SCON_CC11_GENERAL_ATTRIBUTES
    #define SCON_CC11_ALIGNMENT
    #define SCON_CC11_INHERIT_CONSTRUCTORS
    #define SCON_CC11_RVALUEREF_THIS
    #define SCON_CC11_DECLTYPE_CALL
  #endif

#elif defined(_MSC_VER)

  #if defined(_M_X64) || defined(_M_AMD64) || defined(_WIN64)
    #define SCON_x64
  #else
    #define SCON_x86
  #endif

  #if defined(_DEBUG) && !defined(SCON_DEBUG)
    #define SCON_DEBUG 2
  #endif


  #include <yvals.h>

  #if _MSC_VER >= 1600
    #define SCON_CC11_F_LATE_DECLARATOR
    #define SCON_CC11_LOCAL_UNNAMED_TEMPLATE_ARGUMENTS
    #define SCON_CC11_EXTENDED_FRIEND
    #define SCON_CC11_NULLPTR
    #define SCON_CC11_EXTERN_TEMPLATES
    #define SCON_CC11_STATIC_ASSERT
    #define SCON_CC11_RVALUE_REF
    #define SCON_CC11_DECLTYPE
    #define SCON_CC11_AUTO_VARIABLES
    #define SCON_CC11_AUTO_DECLARATOR
    #define SCON_CC11_RIGHT_ANGLE_BRACKETS
  #endif

  #if _MSC_VER >= 1700
    #define SCON_CC11_STRONGLY_TYPED_ENUMS
    #define SCON_CC11_LAMBDA
    #define SCON_CC11_RANGE_BASED_FOR
    #define SCON_CC11_STANDARD_LAYOUT_TYPES
    #define SCON_CC11_ENUM_FORWARD_DECL
    #define SCON_CC11_RAW_STRING_LITERALS
    #define SCON_CC11_MOVE
    #define SCON_CC11_OVERRIDE_FINAL
    #define SCON_CC11_CHRONO
    #define SCON_CC11_SFINAE_TRAITS
    #define SCON_CC11_ARRAY
    #define SCON_CC11_RANDOM
  #endif

  #if _MSC_VER >= 1800
    #define SCON_CC11_DECLTYPE_CALL
    #define SCON_CC11_DEF_DEL_FUNCTIONS
    #define SCON_CC11_EXPLICIT_CONVERSION_OPERATORS
    #define SCON_CC11_DELEGATING_CONSTRUCTORS
    #define SCON_CC11_TEMPLATE_ALIAS
    #define SCON_CC11_F_TEMPLATE_DEFAULT
    #define SCON_CC11_INITIALIZER_LISTS
    #define SCON_CC11_VARIADIC_TEMPLATES
    #define SCON_CC11_NONSTATIC_MEMBER_INIT
    #define SCON_CC11_THREAD
  #endif

#elif defined(__INTEL_COMPILER)

  #if !defined(NDEBUG) && !defined(SCON_DEBUG)
    #define SCON_DEBUG
  #endif

  #if __INTEL_COMPILER >= 1100
    #define SCON_CC11_RVALUE_INIT
    #define SCON_CC11_STATIC_ASSERT
    #define SCON_CC11_AUTO_DECLARATOR
    #define SCON_CC11_RIGHT_ANGLE_BRACKETS
    #define SCON_CC11_EXTERN_TEMPLATES
  #endif

  #if __INTEL_COMPILER >= 1200
    #define SCON_CC11_EXTENDED_FRIEND
    #define SCON_CC11_AUTO_VARIABLES
    #define SCON_CC11_RVALUE_REF
    #define SCON_CC11_DECLTYPE
    #define SCON_CC11_DECLTYPE_CALL
    #define SCON_CC11_DEF_DEL_FUNCTIONS
    #define SCON_CC11_LAMBDA
    #define SCON_CC11_LOCAL_UNNAMED_TEMPLATE_ARGUMENTS
  #endif

  #if __INTEL_COMPILER >= 1206
    #define SCON_CC11_UNIVERSAL_CHAR_LITERALS
    #define SCON_CC11_VARIADIC_TEMPLATES
    #define SCON_CC11_TEMPLATE_ALIAS
    #define SCON_CC11_NULLPTR
    #define SCON_CC11_F_LATE_DECLARATOR
    #define SCON_CC11_SFINAE_EXPRESSIONS
    #define SCON_CC11_GENERAL_ATTRIBUTES
    #define SCON_CC11_F_TEMPLATE_DEFAULT
    #define SCON_CC11_CHRONO
    #define SCON_CC11_SFINAE_TRAITS
    #define SCON_CC11_ARRAY
    #define SCON_CC11_RANDOM
  #endif

  #if __INTEL_COMPILER >= 1300
    #define SCON_CC11_EXPLICIT_ENUM_BASE
    #define SCON_CC11_SCOPED_ENUM
    #define SCON_CC11_EXPLICIT_CONVERSION_OPERATORS
    #define SCON_CC11_RANGE_BASED_FOR
    #define SCON_CC11_ADDITIONAL_TYPE_TRAITS
  #endif

  #if __INTEL_COMPILER >= 1400
    #define SCON_CC11_DELEGATING_CONSTRUCTORS
    #define SCON_CC11_GENRALIZED_CONSTEXPR
    #define SCON_CC11_STRONGLY_TYPED_ENUMS
    #define SCON_CC11_RVALUEREF_THIS
    #define SCON_CC11_RAW_STRING_LITERALS
    #define SCON_CC11_INLINE_NAMESPACES
    #define SCON_CC11_UNRESTRICTED_UNIONS
    #define SCON_CC11_INITIALIZER_LISTS
    #define SCON_CC11_NONSTATIC_MEMBER_INIT
    #define SCON_CC11_EXPLICIT_VIRTUAL_OVERRRIDE
    #define SCON_CC11_MOVE
    #define SCON_CC11_OVERRIDE_FINAL
  #endif

#elif defined(__clang__)

  #error "CLANG SUPPORT NOT YET INCLUDED FOR SCON."

  #ifndef CLANG_VERSION
  #define CLANG_VERSION (__clang_major__ * 10000 \
    + __clang_minor__ * 100 \
    + __clang_patchlevel__)
  #endif

  #if CLANG_VERSION >= 20900
  #endif
  #if CLANG_VERSION >= 30000
  #endif
  #if CLANG_VERSION >= 30100
  #endif
  #if CLANG_VERSION >= 30200
  #endif
  #if CLANG_VERSION >= 30300
  #endif

#endif // compiler defines

#if defined(SCON_CC11_RANDOM)
#include <random>
#include <limits>
#else
#include <cstdlib>
#endif

#if defined(SCON_CC11_THREAD)
#include <thread>
#include <mutex>
#endif

namespace scon
{
  static const char endl = '\n';
}

#endif // ifdef SCON_HEADER
