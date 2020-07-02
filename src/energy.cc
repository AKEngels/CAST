#include <stdexcept>
#include "configuration.h"
#include "energy.h"
#include "energy_int_aco.h"
#include "energy_int_mopac.h"
#include "energy_int_terachem.h"
#include "energy_int_amoeba.h"
#include "energy_int_qmmm_a.h"
#include "energy_int_qmmm_s.h"
#include "energy_int_3layer.h"
#ifdef USE_PYTHON
#include "energy_int_dftbaby.h"
#endif
#include "energy_int_dftb.h"
#include "energy_int_gaussian.h"
#include "energy_int_chemshell.h"
#include "energy_int_psi4.h"
#include "energy_int_orca.h"
#include "coords.h"
#include "Scon/scon_utility.h"


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

static inline energy::interface_base* get_interface(coords::Coordinates* coordinates, config::interface_types::T const& inf)
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
        std::cout << "Mopac chosen for energy calculations.\n";
      }
      return new energy::interfaces::mopac::sysCallInterface(coordinates);
    }
    case config::interface_types::T::QMMM_A:
    {
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "Additive QMMM-Interface choosen for energy calculations.\n";
      }
      return new energy::interfaces::qmmm::QMMM_A(coordinates);
    }
    case config::interface_types::T::QMMM_S:
    {
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "Subtractive QMMM-Interface choosen for energy calculations.\n";
      }
      return new energy::interfaces::qmmm::QMMM_S(coordinates);
    }
    case config::interface_types::T::THREE_LAYER:
    {
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "Three-Layer-Interface choosen for energy calculations.\n";
      }
      return new energy::interfaces::qmmm::THREE_LAYER(coordinates);
    }

    case config::interface_types::T::TERACHEM:
    {
#if defined(USE_MPI)
      if (Config::get().general.verbosity >= 3) std::cout << "Terachem choosen for energy calculations.\n";
      return new energy::interfaces::terachem::mpiInterface(coordinates);
#else
      std::cout << "You need to include MPI for this.\n";
      std::exit(0);
#endif
#ifdef USE_PYTHON
    case config::interface_types::T::DFTBABY:
    {
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "DFTBaby choosen for energy calculations.\n";
      }
      return new energy::interfaces::dftbaby::sysCallInterface(coordinates);
    }
#endif
    case config::interface_types::T::DFTB:
    {
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "DFTBplus choosen for energy calculations.\n";
      }
      return new energy::interfaces::dftb::sysCallInterface(coordinates);
    }
    case config::interface_types::T::GAUSSIAN:
    {
      if (Config::get().general.verbosity >= 3)
      {
        std::cout << "Gaussian chosen for energy calculations.\n";
      }
      return new energy::interfaces::gaussian::sysCallInterfaceGauss(coordinates);
    }
    case config::interface_types::T::CHEMSHELL:
    {
      if (Config::get().general.verbosity >= 3) {
        std::cout << "Chemshell chosen for energy calculations.\n";
      }
      return new energy::interfaces::chemshell::sysCallInterface(coordinates);
    }
    case config::interface_types::T::PSI4:
    {
      if (Config::get().general.verbosity >= 3) {
        std::cout << "Psi4 chosen for energy calculations.\n";
      }
      return new energy::interfaces::psi4::sysCallInterface(coordinates);
    }
    case config::interface_types::T::ORCA:
    {
      if (Config::get().general.verbosity >= 3) {
        std::cout << "ORCA chosen for energy calculations.\n";
      }
      return new energy::interfaces::orca::sysCallInterface(coordinates);
    }
#if defined(USE_MPI)
    case config::interface_types::T::TERACHEM:
    {
      if (Config::get().general.verbosity >= 3) std::cout << "Terachem choosen for energy calculations.\n";
      return new energy::interfaces::terachem::mpiInterface(coordinates);
    }
#endif
  }
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
energy::interface_base* energy::new_interface(coords::Coordinates* coordinates)
{
  return get_interface(coordinates, Config::get().general.energy_interface);
}

energy::interface_base* energy::pre_interface(coords::Coordinates* coordinates)
{
  if (Config::get().general.preopt_interface == config::interface_types::T::ILLEGAL) { return nullptr; }
  energy::interface_base* const r = get_interface(coordinates, Config::get().general.preopt_interface);
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
void energy::interface_base::swap(interface_base& other)
{
  std::swap(energy, other.energy);
  std::swap(periodic, other.periodic);
  std::swap(integrity, other.integrity);
  std::swap(optimizer, other.optimizer);
  std::swap(interactions, other.interactions);
}

void energy::interface_base::print_G_tinkerlike(std::ostream& S, bool const endline) const {
  S << " Cartesian Gradient Breakdown over Individual Atoms :" << std::endl << std::endl;
  S << "  Type      Atom              dE/dX       dE/dY       dE/dZ          Norm" << std::endl << std::endl;
  for (std::size_t k = 0; k < coords->size(); ++k)
  {
    S << " Anlyt";
    S << std::right << std::setw(10) << k + 1U;
    S << "       ";
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).x();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).y();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).z();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4);
    S << std::sqrt(
      coords->g_xyz(k).x() * coords->g_xyz(k).x()
      + coords->g_xyz(k).y() * coords->g_xyz(k).y()
      + coords->g_xyz(k).z() * coords->g_xyz(k).z()) << (endline ? "\n" : "");
  }
}

/*! Override of virtual void to_stream
*
* Virtual void to_stream is decleared in header, but is overrided in energy.cc, so to_stream
* is useable with all general members of the class fot output. For additional members
* of derived classes it must be decleared there additionaly.
*/
void energy::interface_base::to_stream(std::ostream& stream) const
{
  stream << "Energy: " << energy << ", Periodic: " << periodic << ", Integrity: " << integrity << ", Optimizer: " << optimizer << '\n';
}

void energy::interface_base::boundary(coords::Cartesian_Point& r)
{
  coords::Cartesian_Point const halfbox(Config::get().periodics.pb_box / 2.0);
  if (r.x() > halfbox.x())
  {
    r.x() -= Config::get().periodics.pb_box.x();
  }
  else if (r.x() < -halfbox.x())
  {
    r.x() += Config::get().periodics.pb_box.x();
  }
  if (r.y() > halfbox.y())
  {
    r.y() -= Config::get().periodics.pb_box.y();
  }
  else if (r.y() < -halfbox.y())
  {
    r.y() += Config::get().periodics.pb_box.y();
  }
  if (r.z() > halfbox.z())
  {
    r.z() -= Config::get().periodics.pb_box.z();
  }
  else if (r.z() < -halfbox.z())
  {
    r.z() += Config::get().periodics.pb_box.z();
  }
}

// definition of static variable, is necessary 
// (see https://stackoverflow.com/questions/195207/unresolved-external-symbol-on-static-class-members)
std::vector<energy::PointCharge> energy::interface_base::external_charges;