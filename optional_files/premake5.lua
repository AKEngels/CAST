-- premake5.lua
-- Get Premake5 via https://premake.github.io/
--
-- CAST automatic build configuration
-- Targeted at VS2015 and our own Linux-Cluster "ECPC"
--
-- Build for windows: "premake5 vs2015"
-- Open CAST.sln in project/
--
-- Build for Linux on ECPC without MPI: "premake5 gmake"
-- Run "make config=release_x64 CXX=g++-5" in project/
--
-- Build for Linux on ECPC wit MPI: "premake5 gmake --mpi"
-- Run "make config=release_x64 CXX=mpic++" in project/
--

newoption {
   trigger     = "mpi",
   description = "Add MPI Flag to gmake makefiles"
}

workspace "CAST"
	configurations { "Debug", "Release", "Testing", "Armadillo_Release", "Armadillo_Debug" }
    location "../project"
   	platforms { "x86", "x64"}
	filter { "platforms:x86" }
		architecture "x32"
	filter { "platforms:x64" }
		architecture "x64"
	filter{}

project "CAST"
	kind "ConsoleApp"
	language "C++"
	targetdir "../build/"
	files { "../src/**.h", "../src/**.cc" }

	vpaths { ["Headers"] = "../src/**.h" , ["Sources"] = "../src/**.cc"}


	configuration "gmake"
		linkoptions { "-fopenmp" }
		filter { "options:mpi" }
			buildoptions { "-Wextra", "-Wall", "-std=c++0x", "-pedantic", "-fopenmp", "-static", "-DTERACHEM_MPI" }
		filter { "action:gmake" }
		  buildoptions { "-Wextra", "-Wall", "-std=c++0x", "-pedantic", "-fopenmp", "-static", }
		filter { "configurations:Release", "action:gmake" }
			optimize "Full"
		filter { "configurations:Release",  "platforms:x86", "action:gmake"}
		targetname "CAST_linux_x86_release"
		filter { "configurations:Release",  "platforms:x64", "action:gmake"}
		targetname "CAST_linux_x64_release"

		filter { "configurations:Testing", "action:gmake" }
			optimize "Debug"
			defines { "GOOGLE_MOCK" }
		filter { "configurations:Testing",  "platforms:x86", "action:gmake"}
			targetname "CAST_linux_x86_testing"
		filter { "configurations:Testing",  "platforms:x64", "action:gmake"}
		  targetname "CAST_linux_x64_testing"

		filter { "configurations:Debug", "action:gmake" }
			defines { "CAST_DEBUG_DROP_EXCEPTIONS" }
			optimize "Debug"
			flags { "Symbols" }
		filter { "configurations:Debug",  "platforms:x86", "action:gmake"}
		targetname "CAST_linux_x86_debug"
		filter { "configurations:Debug",  "platforms:x64", "action:gmake"}
		targetname "CAST_linux_x64_debug"

    filter { "configurations:Armadillo_Debug", "action:gmake" }
			includedirs { "../optional_files/includes/armadillo/"}
      buildoptions { "-I ../optional_files/includes -DARMA_DONT_USE_WRAPPER -lgfortran" }
			linkoptions { "../optional_files/linux_precompiled_libs/libopenblas.a -I ../optional_files/includes/ -DARMA_DONT_USE_WRAPPER ../optional_files/linux_precompiled_libs/liblapack.a -lgfortran" }
			defines { "CAST_DEBUG_DROP_EXCEPTIONS" }
			optimize "Debug"
			flags { "Symbols" }
			defines { "USE_ARMADILLO" }
		filter { "configurations:Armadillo_Debug",  "platforms:x86", "action:gmake"}
			targetname "CAST_linux_x64_armadillo_debug"
		filter { "configurations:Armadillo_Debug",  "platforms:x64", "action:gmake"}
			targetname "CAST_linux_x64_armadillo_debug"

    filter { "configurations:Armadillo_Release", "action:gmake" }
    includedirs { "../optional_files/includes/armadillo/"}
      buildoptions { "-I ../optional_files/includes -DARMA_DONT_USE_WRAPPER -lgfortran" }
      linkoptions { "../optional_files/linux_precompiled_libs/libopenblas.a -I ../optional_files/includes/ -DARMA_DONT_USE_WRAPPER ../optional_files/linux_precompiled_libs/liblapack.a -lgfortran" }
  		optimize "Full"
			flags { "LinkTimeOptimization" }
  		defines { "USE_ARMADILLO" }
  	filter { "configurations:Armadillo_Release",  "platforms:x86", "action:gmake"}
  		targetname "CAST_linux_x64_armadillo_release"
  	filter { "configurations:Armadillo_Release",  "platforms:x64", "action:gmake"}
  		targetname "CAST_linux_x64_armadillo_release"

	configuration "vs2015"
		targetname "CAST.exe"
		buildoptions { "/openmp" }
		flags { "MultiProcessorCompile" }

		filter { "configurations:Release", "action:vs2015" }
			optimize "Full"
			flags { "LinkTimeOptimization" }
		filter { "configurations:Release",  "platforms:x86", "action:vs2015"}
			targetname "CAST_win_x86_release"
		filter { "configurations:Release",  "platforms:x64", "action:vs2015"}
			targetname "CAST_win_x64_release"

		filter { "configurations:Armadillo_Release", "action:vs2015"}
			includedirs { "../optional_files/includes/armadillo/"}
			libdirs { "../optional_files/windows_precompiled_libs/" }
			links { "blas_win64_MT", "lapack_win64_MT" }
			optimize "Full"
			defines { "USE_ARMADILLO" }
			flags { "LinkTimeOptimization" }
		filter { "configurations:Armadillo_Release",  "platforms:x86", "action:vs2015"}
			targetname "CAST_win_x86_release_lapack"
		filter { "configurations:Armadillo_Release",  "platforms:x64", "action:vs2015"}
			targetname "CAST_win_x64_release_lapack"

		filter { "configurations:Armadillo_Debug", "action:vs2015"}
      includedirs { "../optional_files/includes/armadillo/"}
      libdirs { "../optional_files/windows_precompiled_libs/" }			libdirs { "../libs/" }
			links { "blas_win64_MT", "lapack_win64_MT" }
			defines { "CAST_DEBUG_DROP_EXCEPTIONS" }
			optimize "Debug"
			flags { "Symbols" }
			defines { "USE_ARMADILLO" }
		filter { "configurations:Armadillo_Release",  "platforms:x86", "action:vs2015"}
			targetname "CAST_win_x86_release_lapack"
		filter { "configurations:Armadillo_Release",  "platforms:x64", "action:vs2015"}
			targetname "CAST_win_x64_release_lapack"


		filter { "configurations:Testing", "action:vs2015" }
			optimize "Debug"
			defines { "GOOGLE_MOCK" }
		filter { "configurations:Testing",  "platforms:x86", "action:vs2015"}
		targetname "CAST_win_x86_testing"
		filter { "configurations:Testing",  "platforms:x64", "action:vs2015"}
		targetname "CAST_win_x64_testing"

		filter { "configurations:Debug", "action:vs2015" }
			removefiles { "./src/tests/**.cc", "./src/test/**.h"}
			defines { "CAST_DEBUG_DROP_EXCEPTIONS" }
			optimize "Debug"
			flags { "Symbols" }
		filter { "configurations:Debug",  "platforms:x86", "action:vs2015"}
		targetname "CAST_win_x86_debug"
		filter { "configurations:Debug",  "platforms:x64", "action:vs2015"}
		targetname "CAST_win_x64_debug"
