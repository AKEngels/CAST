-- some local variables
local python27_dir = "C:/Python27" --path to python on windows

newoption {
   trigger     = "mpi",
   description = "Target SMURF cluster with MPI"
}

workspace "CAST"
	configurations { "Debug", "Release", "Armadillo_Debug", "Armadillo_Release", "Testing", "Armadillo_Testing", "Python_Release", "Python_Debug"}
		location "project"
		platforms { "x86", "x64" }
		filter "platforms:x86"
			architecture "x32"
		filter "platforms:x64"
			defines "COMPILEX64"
			architecture "x64"
		filter{}

	-- Do not change the GoogleTest settings if you have no idea what you are doing!
    project "GoogleTest"
        kind"StaticLib"
        language"C++"
        cppdialect"C++14"
        targetdir"libs/%{cfg.buildcfg}"
        location"project/libs/GoogleTest"

		-- setup for google test AND google Mock
		files{ "../submodules/googletest/googletest/src/gtest-all.cc","../submodules/googletest/googlemock/src/gmock-all.cc" }
		-- google test include directories
		sysincludedirs { "../submodules/googletest/googletest/include" }
		includedirs{ "../submodules/googletest/googletest" }
		-- google mock include directories
		sysincludedirs { "../submodules/googletest/googlemock/include" }
		includedirs{ "../submodules/googletest/googlemock" }

        filter "action:vs*"
            system("windows")
		symbols"On"

	project "ALGLIB"
        kind "StaticLib"
        language"C++"
        targetdir"libs/%{cfg.buildcfg}"
        location"project/libs/ALGLIB"
        warnings"Off"
        files { "../submodules/ALGLIB/src/**.h", "../submodules/ALGLIB/src/**.cpp"}

        filter "*Debug"
			symbols "On"
		filter "*Testing"
			symbols "On"

	project "CAST"
		kind "ConsoleApp"
		language "C++"
		targetdir "build"
		files { "../src/**.h", "../src/**.cc"}
		vpaths {
			["Headers/*"] = "../src/**.h",
			["Sources/*"] = "../src/**.cc",
		}
		cppdialect "C++17"
		warnings "Extra"
		sysincludedirs "../submodules/boost"
		sysincludedirs "../submodules/ALGLIB"
		libdirs "../submodules/boost/stage/lib"
		links "ALGLIB"

		--enable if Armadillo Transformations are implemented
		--filter "not Armadillo_*"
		sysincludedirs "../submodules/eigen"
		filter { "not Armadillo_*", "not *Debug" }
			defines "EIGEN_NO_DEBUG"

		filter "*Release"
			optimize "Full"
			flags "LinkTimeOptimization"

		filter "*Debug"
			symbols "On"
			defines "CAST_DEBUG_DROP_EXCEPTIONS"
			debugdir "./build"

		filter "Armadillo_*"
			sysincludedirs { "includes/armadillo", "includes" }
			defines { "ARMA_DONT_USE_WRAPPER", "CAST_USE_ARMADILLO" }

		filter "Python_*"
			defines "USE_PYTHON"

		filter "*Testing"
			symbols "On"
			defines "GOOGLE_MOCK"
			sysincludedirs {"../submodules/googletest/googletest/include", "../submodules/googletest/googlemock/include"}
			links"GoogleTest"

		filter "action:gmake"
			buildoptions { "-Wextra", "-Wall", "-pedantic", "-static", "-fopenmp" }
			linkoptions {"-fopenmp","-lstdc++fs"}

		filter { "options:mpi", "action:gmake" }
			linkoptions "-fopenmp"
			defines "USE_MPI"

		filter {"Release", "platforms:x86", "action:gmake"}
			targetname "CAST_linux_x86_release"
		filter {"Release", "platforms:x64"}
			targetname "CAST_linux_x64_release"

		filter { "Debug", "platforms:x86", "action:gmake" }
			targetname "CAST_linux_x86_debug"
		filter { "Debug", "platforms:x64" }
			targetname "CAST_linux_x64_debug"

		filter {"Armadillo_*", "action:gmake"}
			links { "openblas", "lapack", "pthread" , "gfortran"}
			libdirs { "linux_precompiled_libs" }

		filter {"Armadillo_Release", "platforms:x86", "action:gmake" }
			targetname "CAST_linux_x86_armadillo_release"
		filter {"Armadillo_Release", "platforms:x64", "action:gmake" }
			targetname "CAST_linux_x64_armadillo_release"

		filter {"Armadillo_Debug", "platforms:x86", "action:gmake" }
			targetname "CAST_linux_x86_armadillo_debug"
		filter {"Armadillo_Debug", "platforms:x64", "action:gmake" }
			targetname "CAST_linux_x64_armadillo_debug"

		filter {"Testing", "platforms:x86", "action:gmake" }
			targetname "CAST_linux_x86_testing"
            buildoptions { "-fprofile-arcs", "-ftest-coverage"}
            linkoptions { "-fprofile-arcs" }
		filter {"Testing", "platforms:x64", "action:gmake" }
			targetname "CAST_linux_x64_testing"
            buildoptions { "-fprofile-arcs", "-ftest-coverage"}
            linkoptions { "-fprofile-arcs" }
		filter {"Armadillo_Testing", "platforms:x86", "action:gmake" }
			targetname "CAST_linux_x86_armadillo_testing"
		filter {"Armadillo_Testing", "platforms:x64", "action:gmake" }
			targetname "CAST_linux_x64_armadillo_testing"

		filter {"Python_*", "action:gmake"}
			links { "python2.7", "util", "lapack" }
			linkoptions {  "-export-dynamic", "-pthread", "-ldl", --[["-Wl"--]] }
			libdirs { "linux_precompiled_libs" }
            sysincludedirs "/usr/include/python2.7"

		filter {"Python_Release", "platforms:x86", "action:gmake" }
			targetname "CAST_linux_x86_python_release"
		filter {"Python_Release", "platforms:x64", "action:gmake" }
			targetname "CAST_linux_x64_python_release"

		filter {"Python_Debug", "platforms:x86", "action:gmake" }
			targetname "CAST_linux_x86_python_debug"
		filter {"Python_Debug", "platforms:x64", "action:gmake" }
			targetname "CAST_linux_x64_python_debug"

		filter "action:vs*"
            system("windows")

			buildoptions "/openmp"
			flags "MultiProcessorCompile"

		--enable if Armadillo Transformations are implemented
		--filter { "not Armadillo_*", "action:vs*" }
		    buildoptions "/bigobj"

		filter { "Release", "platforms:x86", "action:vs*" }
			targetname "CAST_win_x86_release"
		filter { "Release", "platforms:x64", "action:vs*" }
			targetname "CAST_win_x64_release"

		filter { "Debug", "platforms:x86", "action:vs*" }
			targetname "CAST_win_x86_debug"
		filter { "Debug", "platforms:x64", "action:vs*" }
			targetname "CAST_win_x64_debug"

		filter {"Armadillo_*", "action:vs*"}
			libdirs "windows_precompiled_libs"
		filter {"Armadillo_*", "platforms:x86", "action:vs*"}
			links {"lapack_x86rel", "blas_x86rel"}
		filter {"Armadillo_*", "platforms:x64", "action:vs*"}
			links {"lapack_win64_MT", "blas_win64_MT"}

		filter { "Armadillo_Release", "platforms:x86", "action:vs*" }
			targetname "CAST_win_x86_armadillo_release"
		filter { "Armadillo_Release", "platforms:x64", "action:vs*" }
			targetname "CAST_win_x64_armadillo_release"

		filter { "Armadillo_Debug", "platforms:x86", "action:vs*"}
			targetname "CAST_win_x86_armadillo_debug"
		filter { "Armadillo_Debug", "platforms:x64", "action:vs*" }
			targetname "CAST_win_x64_armadillo_debug"

		filter {"Python_*", "action:vs*"}
			includedirs(python27_dir .. "/include")
			libdirs(python27_dir .. "/libs")
			links "python27"

		filter {"Python_Release", "platforms:x86", "action:vs*"}
			targetname "CAST_win_x86_python_release"
		filter {"Python_Release", "platforms:x64", "action:vs*"}
			targetname "CAST_win_x64_python_release"

		filter {"Python_Debug", "platforms:x86", "action:vs*"}
			targetname "CAST_win_x86_python_debug"
		filter {"Python_Debug", "platforms:x64", "action:vs*"}
			targetname "CAST_win_x64_python_debug"

		filter { "Testing", "platforms:x86", "action:vs*" }
			targetname "CAST_win_x86_testing"
		filter { "Testing", "platforms:x64", "action:vs*" }
			targetname "CAST_win_x64_testing"

		filter { "Armadillo_Testing", "platforms:x86", "action:vs*"  }
			targetname "CAST_win_x86_armadillo_testing"
		filter { "Armadillo_Testing", "platforms:x64", "action:vs*" }
			targetname "CAST_win_x64_armadillo_testing"
