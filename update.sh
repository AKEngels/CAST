git pull
cd optional_files
/apps/premake5/premake5 gmake
cd project
make config=release_x64 CXX=g++-5
make config=python_release_x64 CXX=g++-5
make config=armadillo_release_x64 CXX=g++-5
make config=testing_x64 CXX=g++-5
cd ..
cd build
./CAST_linux_x64_testing