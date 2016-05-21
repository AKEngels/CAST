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
	configurations { "Debug", "Release", "Testing", "Armadillo_Release" }
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
		filter { "configurations:Release" }
			removefiles { "../src/tests/**.cc", "./src/test/**.h"}
			optimize "Full"
		filter { "configurations:Release",  "platforms:x86", "action:gmake"}
		targetname "CAST_linux_x86_release"	
		filter { "configurations:Release",  "platforms:x64", "action:gmake"}
		targetname "CAST_linux_x64_release"	
		
		filter { "configurations:Testing" }
			optimize "Debug"
			defines { "GOOGLE_MOCK" }
		filter { "configurations:Testing",  "platforms:x86", "action:gmake"}
		targetname "CAST_linux_x86_testing"	
		filter { "configurations:Testing",  "platforms:x64", "action:gmake"}
		targetname "CAST_linux_x64_testing"			
		
		filter { "configurations:Debug" }
			removefiles { "./src/tests/**.cc", "./src/test/**.h"}
			defines { "CAST_DEBUG_DROP_EXCEPTIONS" }
			optimize "Debug"
			flags { "Symbols" }
		filter { "configurations:Debug",  "platforms:x86", "action:gmake"}
		targetname "CAST_linux_x86_debug"
		filter { "configurations:Debug",  "platforms:x64", "action:gmake"}
		targetname "CAST_linux_x64_debug"	
		
	configuration "vs2015"
		targetname "CAST.exe" 
		buildoptions { "/openmp" }
		flags { "MultiProcessorCompile" }

		filter { "configurations:Release" }
			removefiles { "./src/tests/**.cc", "./src/test/**.h"}
			optimize "Full"
			flags { "LinkTimeOptimization" }
		filter { "configurations:Release",  "platforms:x86", "action:vs2015"}
			targetname "CAST_win_x86_release"
		filter { "configurations:Release",  "platforms:x64", "action:vs2015"}
			targetname "CAST_win_x64_release"	

		filter { "configurations:Armadillo_Release"}
			includedirs { "../includes/armadillo/"}
			libdirs { "../libs/" }
			links { "blas_win64_MT", "lapack_win64_MT" }
			removefiles { "./src/tests/**.cc", "./src/test/**.h"}
			optimize "Full"
			defines { "USE_ARMADILLO" }
			flags { "LinkTimeOptimization" }
		filter { "configurations:Armadillo_Release",  "platforms:x86", "action:vs2015"}
			targetname "CAST_win_x86_release_lapack"
		filter { "configurations:Armadillo_Release",  "platforms:x64", "action:vs2015"}
			targetname "CAST_win_x64_release_lapack"
		
		
		filter { "configurations:Testing" }
			optimize "Debug"
			defines { "GOOGLE_MOCK" }
		filter { "configurations:Testing",  "platforms:x86", "action:vs2015"}
		targetname "CAST_win_x86_testing"
		filter { "configurations:Testing",  "platforms:x64", "action:vs2015"}
		targetname "CAST_win_x64_testing"	
			
		filter { "configurations:Debug" }
			removefiles { "./src/tests/**.cc", "./src/test/**.h"}
			defines { "CAST_DEBUG_DROP_EXCEPTIONS" }
			optimize "Debug"
			flags { "Symbols" }
		filter { "configurations:Debug",  "platforms:x86", "action:vs2015"}
		targetname "CAST_win_x86_debug"
		filter { "configurations:Debug",  "platforms:x64", "action:vs2015"}
		targetname "CAST_win_x64_debug"	
