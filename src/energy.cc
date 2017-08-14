#include <stdexcept>
#include "configuration.h"
#include "energy.h"
#include "energy_int_aco.h"
#include "energy_int_mopac.h"
#include "energy_int_terachem.h"
#include "energy_int_amoeba.h"
#include "energy_int_dftb.h"
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
  case config::interface_types::T::DFTB:
   {
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "DFTB choosen for energy calculations.\n";
      }
      return new energy::interfaces::dftb::sysCallInterface(coordinates);
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

//energy::interface_base::interface_base(void)  //Constructor
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