project(CAST LANGUAGES CXX) 
cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)

# Enable LTO for optimized builds
if (${CMAKE_BUILD_TYPE} STREQUAL "Release" OR ${CMAKE_BUILD_TYPE} STREQUAL "MinSizeRel")
    message(STATUS "Building with link-time optimization")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION On)
endif()

# Compile ALGLIB
file(GLOB ALGLIB_SRC "submodules/ALGLIB/src/*.cpp")
add_library(ALGLIB STATIC ${ALGLIB_SRC})

# Compile GTest and GMock
add_library(GoogleTest STATIC "submodules/googletest/googletest/src/gtest-all.cc" "submodules/googletest/googlemock/src/gmock-all.cc")
target_include_directories(GoogleTest SYSTEM PRIVATE "submodules/googletest/googletest/include" "submodules/googletest/googlemock/include")
target_include_directories(GoogleTest PRIVATE "submodules/googletest/googletest" "submodules/googletest/googlemock")

# Set up CAST options
option(CAST_USE_ARMADILLO "Use Armadillo as backend for matrix computations" Off)
option(CAST_USE_OPTPP "Use Opt++ library for nonlinear optimization" Off)
option(CAST_USE_PYTHON "Enable Python integration" Off)

# Add CAST targets
set(CAST_EXEC ${PROJECT_NAME})
set(CAST_TEST_EXEC ${CAST_EXEC}_test)

file(GLOB_RECURSE CAST_SRC CONFIGURE_DEPENDS src/*.cc src/*.h)
file(GLOB_RECURSE CAST_TEST_SRC CONFIGURE_DEPENDS src/Tests/*.cc src/Tests/*.h)

# Remove all source files that are needed in either CAST or unit tests.
list(REMOVE_ITEM CAST_SRC ${CMAKE_SOURCE_DIR}/src/main.cc)
list(FILTER CAST_SRC EXCLUDE REGEX ${CMAKE_SOURCE_DIR}/src/Tests)

add_library(CAST_COMMON OBJECT ${CAST_SRC})

add_executable(${CAST_EXEC} src/main.cc)
target_link_libraries(${CAST_EXEC} PUBLIC CAST_COMMON)

add_executable(${CAST_TEST_EXEC} ${CAST_TEST_SRC})
target_link_libraries(${CAST_TEST_EXEC} PUBLIC CAST_COMMON GoogleTest)
target_compile_definitions(${CAST_TEST_EXEC} PUBLIC "GOOGLE_MOCK")
target_include_directories(${CAST_TEST_EXEC} SYSTEM PUBLIC "submodules/googletest/googletest/include" "submodules/googletest/googlemock/include")

# Drop exceptions when building in debug mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(CAST_COMMON PUBLIC "CAST_DEBUG_DROP_EXCEPTIONS")
endif()

# Enable OpenMP when available
find_package(OpenMP)
if (${OpenMP_CXX_FOUND})
    target_link_libraries(CAST_COMMON PUBLIC OpenMP::OpenMP_CXX)
endif()

# Check if we compile for x64
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    target_compile_definitions(CAST_COMMON PUBLIC "COMPILEX64")
endif()

# Set up compiler-dependent warning options and libraries
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(CAST_COMMON PUBLIC "-Wall" "-Wextra" "-pedantic")
    target_link_libraries(CAST_COMMON PUBLIC "stdc++fs")
endif()

target_include_directories(CAST_COMMON SYSTEM PUBLIC "submodules/eigen" "submodules/boost" "submodules/ALGLIB")

# Use Armadillo if requested
if (${CAST_USE_ARMADILLO})
    message(STATUS "Using Armadillo")
    target_compile_definitions(CAST_COMMON PUBLIC "ARMA_DONT_USE_WRAPPER" "CAST_USE_ARMADILLO")
    target_include_directories(CAST_COMMON SYSTEM PUBLIC "optional_files/includes/armadillo" "optional_files/includes")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        link_directories("optional_files/windows_precompiled_libs")
        if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
            target_link_libraries(CAST_COMMON PUBLIC "lapack_win64_MT" "blas_win64_MT")
        else()
            target_link_libraries(CAST_COMMON PUBLIC "lapack_x86rel" "blas_x86rel")
        endif()
    else()
        target_link_libraries(CAST_COMMON PUBLIC "openblas" "lapack")
    endif()
endif()

# Use OPT++ if requested
if (${CAST_USE_OPTPP})
    message(STATUS "Using OPT++")
    find_library(OPTPP_LIB "opt" REQUIRED)
    find_library(NEWMAT_LIB "newmat" REQUIRED)
    target_compile_definitions(CAST_COMMON PUBLIC "HAVE_NAMESPACES" "USE_OPTPP")
    target_include_directories(CAST_COMMON SYSTEM PUBLIC "submodules/optpp/include" "submodules/optpp/newmat11")
    target_link_libraries(CAST_COMMON PUBLIC ${OPTPP_LIB} ${NEWMAT_LIB} lapack openblas)
endif()

# Use Python if requested
if (${CAST_USE_PYTHON})
    message(STATUS "Using Python")
    target_compile_definitions(CAST_COMMON PUBLIC "USE_PYTHON")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        set(PYTHON27_DIR "C:/Python27")
        target_include_directories(CAST_COMMON SYSTEM PUBLIC "${PYTHON27_DIR}/include")
        target_link_directories(CAST_COMMON PUBLIC "${PYTHON27_DIR}/libs")
        target_link_libraries(CAST_COMMON PUBLIC python27)
    else()
        target_include_directories(CAST_COMMON SYSTEM PUBLIC "/usr/include/python2.7")
        target_link_libraries(CAST_COMMON PUBLIC "python2.7")
    endif()
endif()

target_link_libraries(CAST_COMMON PUBLIC "ALGLIB")
