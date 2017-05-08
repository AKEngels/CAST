#include "scon_utility.h"

namespace scon
{
/*! Function to call other programms
*  
*/

  int system_call(std::string const & command_line)
  {
#if defined (_MSC_VER)
    // get a modifiable character sequence of the command: 
#if defined _UNICODE
    using cur_char = wchar_t;
    using cur_string = std::wstring;
#else
    using cur_char = char;
    using cur_string = std::string;
#endif
    cur_string wide_cmd(command_line.begin(), command_line.end());
    std::vector<cur_char> w_cmdl(wide_cmd.c_str(), wide_cmd.c_str() + wide_cmd.length() + 1u);
    STARTUPINFO si;
    PROCESS_INFORMATION pi;
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));

    if (CreateProcess(NULL, w_cmdl.data(), NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, NULL, &si, &pi))
    {
      WaitForSingleObject(pi.hProcess, INFINITE);
      CloseHandle(pi.hProcess);
      CloseHandle(pi.hThread);
      return 0;
    }
    else /*if call failed*/ return 666;
#else
    return system(command_line.c_str());
#endif
  }


}

std::string scon::stringseparation::separateString(std::string inString)
{
  std::string oString = ""; charTypeT st = other;
  for (auto c : inString) {
    if ((st == alpha && charTypestring(c) == digit) || (st == digit && charTypestring(c) == alpha)
      || (st == other && charTypestring(c) == digit) || (st == digit && charTypestring(c) == other)
      || (st == other && charTypestring(c) == alpha) || (st == alpha && charTypestring(c) == other))
      oString.push_back(' ');
    oString.push_back(c); st = charTypestring(c);
  }
  return oString;
}
