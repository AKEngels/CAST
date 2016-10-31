#ifdef GOOGLE_MOCK
#include <gtest/gtest.h>

//TESTING MAIN
//////////////////////////
//                      //
//  Test Driven Devel   //
//                      //
//////////////////////////

int main(int argc, char** argv) {

	testing::InitGoogleTest(&argc, argv);
	int result = RUN_ALL_TESTS();
  
  // If you want to run TESTS live from Visual Studio,
  // this will help
  std::cin >> result;

  return result;
}
#endif