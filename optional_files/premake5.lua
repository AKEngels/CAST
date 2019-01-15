-- premake5.lua
-- Get Premake5 via https://premake.github.io/
--
-- CAST automatic build configuration
-- Targeted at VS2017 and our own Linux-Cluster "ECPC"
--
-- Build for windows: "premake5 vs2017"
-- Open CAST.sln in project/
--
-- Build for Linux on ECPC "premake5 gmake"
-- Run "make config=armadillo_release_x64 CXX=g++-5" in project/
--
-- Build for Linux on Smurf: "premake5 gmake --mpi"
-- Run "make config=armadillo_release_x64 CXX='/apps/mpich/2.1.4.1p1/bin/mpicxx -cxx=/apps/gcc-6.1/bin/g++-6.1 -static-libstdc++ -static-libgcc'" in project/
--

newoption {
   trigger     = "mpi",
   description = "Target SMURF cluster with MPI"
}

function os.winSdkVersion()
    local reg_arch = iif( os.is64bit(), "\\Wow6432Node\\", "\\" )
          if os.is("windows") then
                    local sdk_version = os.getWindowsRegistry( "HKLM:SOFTWARE" .. reg_arch .."Microsoft\\Microsoft SDKs\\Windows\\v10.0\\ProductVersion" )
                if sdk_version ~= nil then return sdk_version end
    else return "nothing" end
end


workspace "CAST"
  configurations { "Debug", "Release", "Testing", "Armadillo_Testing", "Armadillo_Release", "Armadillo_Debug", "Python_Release", "Python_Debug" }
    location "../optional_files/project"
    platforms { "x86", "x64"}
  filter { "platforms:x86" }
    architecture "x32"
  filter { "platforms:x64" }
    defines { "COMPILEX64"}
    architecture "x64"
  filter{}

project "CAST"
  kind "ConsoleApp"
  language "C++"
  targetdir "../optional_files/build/"
  files { "../src/**.h", "../src/**.cc", "../src/gtest/**.cc", "../submodules/cubature/hcubature.c", "../submodules/cubature/pcubature.c","../submodules/cubature/cubature.h" }

  vpaths { 
  ["Headers"] = "../src/**.h" , 
  ["Sources"] = "../src/*.cc", 
  ["Testing"] = "../src/gtest/**.cc", 
  ["cubature"] = { "../submodules/cubature/hcubature.c", "../submodules/cubature/pcubature.c","../submodules/cubature/cubature.h"} }


  configuration "gmake"
    includedirs { "../submodules/eigen/Eigen/", "../submodules/cubature/" }
    linkoptions { "-fopenmp" }
                targetname "CAST_undefined.exe"
    filter { "options:mpi" }
      defines { "USE_MPI" }
    filter { "action:gmake" }
      buildoptions { "-Wextra", "-Wall", "-std=c++0x", "-pedantic", "-fopenmp", "-static", }
    filter { "configurations:Release", "action:gmake" }
      optimize "Full"
    filter { "configurations:Release",  "platforms:x86", "action:gmake"}
      targetname "CAST_linux_x86_release"
    filter { "configurations:Release",  "platforms:x64", "action:gmake"}
      targetname "CAST_linux_x64_release"
    filter { "configurations:Python_Release",  "platforms:x64", "action:gmake"}
      optimize "Full"
      targetname "CAST_linux_x64_python_release"
      includedirs { "/apps/python27/include/python2.7" }
      defines {"USE_PYTHON"}
      links {"python2.7", "util", "lapack"}
      linkoptions {"-Xlinker", "-export-dynamic", "-Wl,-rpath,/apps/lapack-3.4.2/lib"}
    filter { "configurations:Python_Release",  "platforms:x86", "action:gmake"}
      optimize "Full"
      targetname "CAST_linux_x86_python_release"
      includedirs { "/apps/python27/include/python2.7" }
      defines {"USE_PYTHON"}
      links {"python2.7", "util", "lapack"}
      linkoptions {"-Xlinker", "-export-dynamic", "-Wl,-rpath,/apps/lapack-3.4.2/lib"}

    filter { "configurations:Armadillo_Testing", "action:gmake" }
      optimize "Debug"
      defines { "GOOGLE_MOCK", "USE_ARMADILLO", "ARMA_DONT_USE_WRAPPER" }
      includedirs { "./includes/gtest/", "../optional_files/includes/armadillo/", "../submodules/cubature/" }
      buildoptions { "-I ../optional_files/includes -I ../includes -lgfortran" }
      linkoptions { "../linux_precompiled_libs/libgmock.a ../linux_precompiled_libs/libopenblas.a ../linux_precompiled_libs/liblapack.a -lgfortran" }
      flags { "LinkTimeOptimization" }
    filter { "configurations:Armadillo_Testing",  "platforms:x86", "action:gmake"}
      targetname "CAST_linux_x86_armadillo_testing"
    filter { "configurations:Armadillo_Testing",  "platforms:x64", "action:gmake"}
      targetname "CAST_linux_x64_armadillo_testing"

    filter { "configurations:Testing", "action:gmake" }
      optimize "Debug"
      defines { "GOOGLE_MOCK" }
      includedirs { "./includes/gtest/", "../submodules/eigen/Eigen/", "../submodules/cubature/"}
      buildoptions { "-I ../optional_files/includes" }
      linkoptions { "../linux_precompiled_libs/libgmock.a -I ../optional_files/includes/" }
      flags { "LinkTimeOptimization" }
    filter { "configurations:Testing",  "platforms:x86", "action:gmake"}
      targetname "CAST_linux_x86_testing"
    filter { "configurations:Testing",  "platforms:x64", "action:gmake"}
      targetname "CAST_linux_x64_testing"

    filter { "configurations:Debug", "action:gmake" }
      defines { "CAST_DEBUG_DROP_EXCEPTIONS" }
      optimize "Debug"
      includedirs { "../submodules/eigen/Eigen/", "../submodules/cubature/"}
      flags { "Symbols" }
    filter { "configurations:Debug",  "platforms:x86", "action:gmake"}
      targetname "CAST_linux_x86_debug"
    filter { "configurations:Debug",  "platforms:x64", "action:gmake"}
      targetname "CAST_linux_x64_debug"
    filter { "configurations:Python_Debug",  "platforms:x86", "action:gmake"}
      optimize "Debug"
      targetname "CAST_linux_x86_python_debug"
      includedirs { "/apps/python27/include/python2.7" }
      defines {"USE_PYTHON"}
      links {"python2.7", "util", "lapack"}
      linkoptions {"-Xlinker", "-export-dynamic", "-Wl,-rpath,/apps/lapack-3.4.2/lib"}
    filter { "configurations:Python_Debug",  "platforms:x64", "action:gmake"}
      optimize "Debug"
      targetname "CAST_linux_x64_python_debug"
      includedirs { "/apps/python27/include/python2.7" }
      defines {"USE_PYTHON"}
      links {"python2.7", "util", "lapack"}
      linkoptions {"-Xlinker", "-export-dynamic", "-Wl,-rpath,/apps/lapack-3.4.2/lib"}

    filter { "configurations:Armadillo_Debug", "action:gmake" }
      includedirs { "../optional_files/includes/armadillo/", "../submodules/cubature/"}
      buildoptions { "-I ../includes -DARMA_DONT_USE_WRAPPER -lgfortran" }
      linkoptions { "-lblas  -I ../optional_files/includes/ -DARMA_DONT_USE_WRAPPER -llapack  -lgfortran" }
      defines { "CAST_DEBUG_DROP_EXCEPTIONS" }
      optimize "Debug"
      flags { "Symbols" }
      defines { "CAST_USE_ARMADILLO" }
    filter { "configurations:Armadillo_Debug",  "platforms:x86", "action:gmake"}
      targetname "CAST_linux_x86_armadillo_debug"
    filter { "configurations:Armadillo_Debug",  "platforms:x64", "action:gmake"}
      targetname "CAST_linux_x64_armadillo_debug"

    filter { "configurations:Armadillo_Release", "action:gmake" }
      includedirs { "./includes/armadillo/", "../submodules/cubature/"}
      buildoptions { "-I ../optional_files/includes -DARMA_DONT_USE_WRAPPER -lgfortran" }
      linkoptions { "-lblas -I ../optional_files/includes/ -DARMA_DONT_USE_WRAPPER -llapack -lgfortran" }
      optimize "Full"
      flags { "LinkTimeOptimization" }
      defines { "CAST_USE_ARMADILLO"}
    filter { "configurations:Armadillo_Release",  "platforms:x86", "action:gmake"}
      targetname "CAST_linux_x64_armadillo_release"
    filter { "configurations:Armadillo_Release",  "platforms:x64", "action:gmake"}
      targetname "CAST_linux_x64_armadillo_release"



  configuration "vs2017"
    systemversion(os.winSdkVersion() .. ".0")
    targetname "CAST_undefined.exe"
    debugdir "../optional_files/build"
    flags { "MultiProcessorCompile" }
    filter { "configurations:Release", "action:vs2017" }
      optimize "Full"
      flags { "LinkTimeOptimization" }
    filter { "configurations:Release",  "platforms:x86", "action:vs2017"}
      targetname "CAST_win_x86_release"
      defines {"EIGEN_NO_DEBUG"}
      includedirs { "../submodules/eigen/Eigen/", "../submodules/cubature/"}
    filter { "configurations:Release",  "platforms:x64", "action:vs2017"}
      targetname "CAST_win_x64_release"
      defines {"EIGEN_NO_DEBUG"}
      includedirs { "../submodules/eigen/Eigen/", "../submodules/cubature/"}

    filter { "configurations:Python_Release",  "platforms:x64", "action:vs2017"}
      optimize "Full"
      targetname "CAST_win_x64_python_release"
      defines {"EIGEN_NO_DEBUG", "USE_PYTHON"}
      includedirs { "../submodules/eigen/Eigen/", "C:/Python27/include"}
      libdirs {"C:/Python27/libs"}
      links {"python27"}
    filter { "configurations:Python_Release",  "platforms:x86", "action:vs2017"}
      optimize "Full"
      targetname "CAST_win_x86_python_release"
      defines {"EIGEN_NO_DEBUG", "USE_PYTHON"}
      includedirs { "../submodules/eigen/Eigen/", "C:/Python27/include"}
      libdirs {"C:/Python27/libs"}
      links {"python27"}

    filter { "configurations:Armadillo_Release", "action:vs2017"}
      includedirs { "../optional_files/includes/armadillo/", "../submodules/cubature/"}
      libdirs { "../optional_files/windows_precompiled_libs/" }
      optimize "Full"
      defines { "CAST_USE_ARMADILLO" }
      flags { "LinkTimeOptimization" }
    filter { "configurations:Armadillo_Release",  "platforms:x86", "action:vs2017"}
      targetname "CAST_win_x86_release_lapack"
      links { "blas_x86rel", "lapack_x86rel" }
    filter { "configurations:Armadillo_Release",  "platforms:x64", "action:vs2017"}
      targetname "CAST_win_x64_release_lapack"
      links { "blas_win64_MT", "lapack_win64_MT" }

    filter { "configurations:Armadillo_Debug", "action:vs2017"}
      includedirs { "../optional_files/includes/armadillo/", "../submodules/cubature/"}
      libdirs { "../optional_files/windows_precompiled_libs/" }
      defines { "CAST_DEBUG_DROP_EXCEPTIONS" }
      optimize "Debug"
      flags { "Symbols" }
      defines { "CAST_USE_ARMADILLO" }
    filter { "configurations:Armadillo_Debug",  "platforms:x86", "action:vs2017"}
      targetname "CAST_win_x86_debug_lapack"
      links { "blas_x86rel", "lapack_x86rel" }
    filter { "configurations:Armadillo_Debug",  "platforms:x64", "action:vs2017"}
      targetname "CAST_win_x64_debug_lapack"
      links { "blas_win64_MT", "lapack_win64_MT" }


    filter { "configurations:Testing", "action:vs2017" }
      optimize "Debug"
      includedirs { "../optional_files/includes/gtest/", "../submodules/eigen/Eigen/", "../submodules/cubature/"}
      defines { "GOOGLE_MOCK" }
      libdirs { "../optional_files/windows_precompiled_libs/" }
      links { "gmock" }
    linkoptions {"/DEBUG"}
    filter { "configurations:Testing",  "platforms:x86", "action:vs2017"}
    targetname "CAST_win_x86_testing"
    filter { "configurations:Testing",  "platforms:x64", "action:vs2017"}
    targetname "CAST_win_x64_testing"

    filter { "configurations:Armadillo_Testing", "action:vs2017" }
      optimize "Debug"
      defines { "GOOGLE_MOCK", "CAST_USE_ARMADILLO" }
      includedirs { "../optional_files/includes/armadillo/", "../optional_files/includes/gtest/", "../submodules/cubature/"}
      libdirs { "../optional_files/windows_precompiled_libs/" }
    linkoptions {"/DEBUG"}
      filter { "configurations:Armadillo_Testing",  "platforms:x86", "action:vs2017"}
        targetname "CAST_win_x86_testing_lapack"
        links { "blas_x86rel", "lapack_x86rel", "gmock" }
      filter { "configurations:Armadillo_Testing",  "platforms:x64", "action:vs2017"}
        targetname "CAST_win_x64_testing_lapack"
        links { "blas_win64_MT", "lapack_win64_MT", "gmock" }

    filter { "configurations:Debug", "action:vs2017" }
      removefiles { "./src/tests/**.cc", "./src/test/**.h"}
      defines { "CAST_DEBUG_DROP_EXCEPTIONS" }
      includedirs { "../submodules/eigen/Eigen/", "../submodules/cubature/"}
      optimize "Debug"
      flags { "Symbols" }
    filter { "configurations:Debug",  "platforms:x86", "action:vs2017"}
      targetname "CAST_win_x86_debug"
    filter { "configurations:Debug",  "platforms:x64", "action:vs2017"}
      targetname "CAST_win_x64_debug"
    filter { "configurations:Python_Debug",  "platforms:x86", "action:vs2017"}
      optimize "Debug"
      targetname "CAST_win_x86_python_debug"
      includedirs {"../submodules/eigen/Eigen/", "C:/Python27/include"}
      defines {"USE_PYTHON"}
      libdirs {"C:/Python27/libs"}
      links {"python27"}
    filter { "configurations:Python_Debug",  "platforms:x64", "action:vs2017"}
      optimize "Debug"
      targetname "CAST_win_x64_python_debug"
      includedirs {"../submodules/eigen/Eigen/", "C:/Python27/include"}
      defines {"USE_PYTHON"}
      libdirs {"C:/Python27/libs"}
      links {"python27"}
