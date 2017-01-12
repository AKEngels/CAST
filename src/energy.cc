#include <stdexcept>
#include "configuration.h"
#include "energy.h"
#include "energy_int_aco.h"
#include "energy_int_mopac.h"
#include "energy_int_terachem.h"
#include "energy_int_amoeba.h"
#include "energy_int_gaussian.h"
#include "coords.h"
#include "scon_utility.h"


/*! Creates a specific energy interface
*
* This function returns a pointer to
* a new energy interface. It should not be called
* manually but is used by the energy::new_interface
* function. Which interface is created depends on the 
* specifications in the global Config instance
*
* @param coordinates: Pointer to coordinates object for which energy interface will perform
* @return: Base-Class Pointer to the energy interface. Nullpointer is returned if something went wrong.
*/
static inline energy::interface_base * get_interface (coords::Coordinates * coordinates, config::interface_types::T const &inf)
{
  switch (inf)
  {
  case config::interface_types::T::ILLEGAL:
    {
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "Illegal Energy interface.\n";
      }
      return nullptr;
    }  
  case config::interface_types::T::AMOEBA:
    {
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "AMOEBA force field energy interface chosen for energy calculations\n";
      }
      return new energy::interfaces::amoeba::amoeba_ff(coordinates);
    }
  case config::interface_types::T::MOPAC:
    {
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "Mopac choosen for energy calculations.\n";
      }
      return new energy::interfaces::mopac::sysCallInterface(coordinates);
    }
  case config::interface_types::T::GAUSSIAN:
    {
     if (Config::get().general.verbosity >= 3)
     {
       std::cout << "Gaussian choosen for energy calculations.\n";
     }
     return new energy::interfaces::gaussian::sysCallInterfaceGauss(coordinates);
  }
#if defined(USE_MPI)
  case config::interface_types::T::TERACHEM:
    {
      if (Config::get().general.verbosity >= 3) std::cout << "Terachem choosen for energy calculations.\n";
      return new energy::interfaces::terachem::mpiInterface(coordinates);
    }
#endif
  default:
    {
    if (Config::get().general.verbosity >= 3) std::cout << "Default (force field) interface choosen for energy calculations.\n";
      return new energy::interfaces::aco::aco_ff(coordinates);
    }
  }
}

/*! Creates a new energy interface
 *
 * This function returns a pointer to
 * a fresh, new energy interface. Which interface
 * is created depends on the specifications in the global
 * Config instance
 *
 * @param coordinates: Pointer to coordinates object for which energy interface will perform
 * @return: Base-Class Pointer to the energy interface. Nullpointer is returned if something went wrong.
 */
energy::interface_base* energy::new_interface (coords::Coordinates * coordinates)
{
  return get_interface(coordinates, Config::get().general.energy_interface);
}

energy::interface_base* energy::pre_interface(coords::Coordinates * coordinates)
{
  energy::interface_base * const r = get_interface(coordinates, Config::get().general.preopt_interface);
  //std::cout << "New Preinterface " << (r ? r : 0) << "." << std::endl;
  return r;
}

/*! Override of virtual void swap
*
* Virtual void swap is decleared in header, but is overrided in energy.cc, so swap is
* useable with all general members of the class. For additional members 
* of derived classes it must be decleared there additionaly to
* "void swap(sysCallInterface&)"
*/
void energy::interface_base::swap (interface_base &other)    
{
  std::swap(energy,    other.energy);
  std::swap(pb_max,    other.pb_max);
  std::swap(pb_min,    other.pb_min);
  std::swap(pb_dim,    other.pb_dim);
  std::swap(periodic,  other.periodic);
  std::swap(integrity, other.integrity);
  std::swap(optimizer, other.optimizer);
  std::swap(interactions, other.interactions);
}


/*! Override of virtual void to_stream
*
* Virtual void to_stream is decleared in header, but is overrided in energy.cc, so to_stream 
* is useable with all general members of the class fot output. For additional members
* of derived classes it must be decleared there additionaly.
*/
void energy::interface_base::to_stream (std::ostream &stream) const
{
  stream << "Energy: " << energy << ", Periodic: " << periodic << ", Integrity: " << integrity << ", Optimizer: " << optimizer << '\n';
  stream << "Periodics:  Max: " << pb_max << ", Min: " << pb_min << ", Dim: " << pb_dim << '\n';
}

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
  else /*if MOPAC call failed*/ return 666;
#else
  return system(command_line.c_str());
#endif
}

int scon::system_call(std::string const & command_line)
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