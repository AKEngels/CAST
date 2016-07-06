#include <stdexcept>
#include "configuration.h"
#include "energy.h"
#include "energy_int_aco.h"
#include "energy_int_mopac.h"
#include "energy_int_terachem.h"
#include "energy_int_amoeba.h"
#include "energy_int_qmmm.h"
#include "coords.h"

static inline energy::interface_base * get_interface (coords::Coordinates * cp, config::interface_types::T const &inf)
{
  switch (inf)
  {
  case config::interface_types::T::ILLEGAL:
    {
      if (Config::get().general.verbosity > 29)
      {
        std::cout << "Illegal Energy interface.\n";
      }
      return nullptr;
    }  
  case config::interface_types::T::AMOEBA:
    {
      if (Config::get().general.verbosity > 29)
      {
        std::cout << "Not (yet) existent energy interface choosen.\n";
      }
      return new energy::interfaces::amoeba::amoeba_ff(cp);
    }
  case config::interface_types::T::MOPAC:
    {
      if (Config::get().general.verbosity > 29) 
      {
        std::cout << "Mopac choosen for energy calculations.\n";
      }
      std::cout << "GET MOPAC!\n";
      return new energy::interfaces::mopac::sysCallInterface(cp);
    }
  case config::interface_types::T::QMMM:
  {
    if (Config::get().general.verbosity > 29)
    {
      std::cout << "QMMM-Interface choosen for energy calculations.\n";
    }
    return new energy::interfaces::qmmm::QMMM(cp);
  }
#if defined(USE_MPI)
  case config::interface_types::T::TERACHEM:
    {
      if (Config::get().general.verbosity > 29) std::cout << "Terachem choosen for energy calculations.\n";
      return new energy::interfaces::terachem::mpiInterface(cp);
    }
#endif
  default:
    {
    if (Config::get().general.verbosity > 29) std::cout << "Default (force field) energy interface choosen.\n";
      return new energy::interfaces::aco::aco_ff(cp);
    }
  }
}

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

void energy::interface_base::to_stream (std::ostream &stream) const
{
  stream << "Energy: " << energy << ", Periodic: " << periodic << ", Integrity: " << integrity << ", Optimizer: " << optimizer << '\n';
  stream << "Periodics:  Max: " << pb_max << ", Min: " << pb_min << ", Dim: " << pb_dim << '\n';
}

