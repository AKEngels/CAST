#ifdef GOOGLE_MOCK
#include <gtest/gtest.h>
#include "../configuration.h"
#pragma once

//TESTING MAIN
//////////////////////////
//                      //
//  Test Driven Devel   //
//                      //
//////////////////////////

int main(int argc, char** argv) {

  // Parse config file and command line 
  auto config_filename = config::config_file_from_commandline(argc, argv);
  Config main_configuration(config_filename);
  config::parse_command_switches(argc, argv);
  Config::set().general.verbosity = 3u;

	testing::InitGoogleTest(&argc, argv);
	int result = RUN_ALL_TESTS();
  
  // If you want to run TESTS live from Visual Studio,
  // this will help
  std::cin >> result;

  return result;
}
#endif