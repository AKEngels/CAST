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

# Set up CAST
set(CAST_EXEC ${PROJECT_NAME})

option(CAST_USE_ARMADILLO "Use Armadillo as backend for matrix computations" Off)
option(CAST_USE_OPTPP "Use Opt++ library for nonlinear optimization" Off)

file(GLOB_RECURSE CAST_SRC CONFIGURE_DEPENDS src/*.cc src/*.h)

add_executable(${CAST_EXEC} ${CAST_SRC})
add_executable(${CAST_EXEC}_test ${CAST_SRC})
target_compile_definitions(${CAST_EXEC}_test PUBLIC "GOOGLE_MOCK")
target_include_directories(${CAST_EXEC}_test PUBLIC "submodules/googletest/googletest/include" "submodules/googletest/googlemock/include")
target_link_libraries(${CAST_EXEC}_test PUBLIC GoogleTest)

set(CAST_COMPILE_OPTIONS "")
set(CAST_DEFINITIONS "")
set(CAST_INCLUDE_DIRS "")
set(CAST_LIBRARIES "")

# Enable OpenMP when available
find_package(OpenMP)
if (${OpenMP_CXX_FOUND})
    list(APPEND CAST_LIBRARIES OpenMP::OpenMP_CXX)
endif()

# Check if we compile for x64
if (${CMAKE_SIZEOF_VOID_P} STREQUAL "8")
    list(APPEND CAST_DEFINITIONS "COMPILEX64")
endif()

# Set up warning options
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
    list(APPEND CAST_COMPILE_OPTIONS "/W3")
else()
    list(APPEND CAST_COMPILE_OPTIONS "-Wall" "-Wextra" "-pedantic")
endif()

list(APPEND CAST_INCLUDE_DIRS "submodules/eigen" "submodules/boost" "submodules/ALGLIB")

# Use Armadillo if requested
if (${CAST_USE_ARMADILLO})
    message(STATUS "Using Armadillo")
    list(APPEND CAST_DEFINITIONS "ARMA_DONT_USE_WRAPPER" "CAST_USE_ARMADILLO")
    list(APPEND CAST_INCLUDE_DIRS "optional_files/includes/armadillo" "optional_files/includes")
    list(APPEND CAST_LIBRARIES "openblas" "lapack")
endif()

# Use OPT++ if requested
if (${CAST_USE_OPTPP})
    message(STATUS "Using OPT++")
    find_library(OPTPP_LIB "opt" REQUIRED)
    find_library(NEWMAT_LIB "newmat" REQUIRED)
    list(APPEND CAST_DEFINITIONS "HAVE_NAMESPACES" "USE_OPTPP")
    list(APPEND CAST_INCLUDE_DIRS "submodules/optpp/include" "submodules/optpp/newmat11")
    list(APPEND CAST_LIBRARIES ${OPTPP_LIB} ${NEWMAT_LIB} lapack openblas)
endif()

list(APPEND CAST_LIBRARIES "ALGLIB" "stdc++fs")

target_compile_options(${CAST_EXEC}      PUBLIC ${CAST_COMPILE_OPTIONS})
target_compile_options(${CAST_EXEC}_test PUBLIC ${CAST_COMPILE_OPTIONS})
target_compile_definitions(${CAST_EXEC}      PUBLIC ${CAST_DEFINITIONS})
target_compile_definitions(${CAST_EXEC}_test PUBLIC ${CAST_DEFINITIONS})
target_include_directories(${CAST_EXEC}      SYSTEM PUBLIC ${CAST_INCLUDE_DIRS})
target_include_directories(${CAST_EXEC}_test SYSTEM PUBLIC ${CAST_INCLUDE_DIRS})
target_link_libraries(${CAST_EXEC}      PUBLIC ${CAST_LIBRARIES})
target_link_libraries(${CAST_EXEC}_test PUBLIC ${CAST_LIBRARIES})
