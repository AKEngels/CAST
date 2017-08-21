#ifdef GOOGLE_MOCK
#include <gtest/gtest.h>
#include "../configuration.h"
#ifdef _MSC_VER
#include "../win_inc.h"
#endif
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
#ifdef _MSC_VER 
  // make window stay open in debug session on windows
  if (IsDebuggerPresent()) std::system("pause");
#endif

  return result;
}
#endif