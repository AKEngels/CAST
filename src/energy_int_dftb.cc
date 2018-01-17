#ifdef USE_PYTHON
#include "energy_int_dftb.h"


/*
dftb sysCall functions
*/





energy::interfaces::dftb::sysCallInterface::sysCallInterface(coords::Coordinates * cp) :
  energy::interface_base(cp),
  e_bs(0.0), e_coul(0.0), e_rep(0.0), e_tot(0.0)
{
    
}

energy::interfaces::dftb::sysCallInterface::sysCallInterface(sysCallInterface const & rhs, coords::Coordinates *cobj) :
  interface_base(cobj),
  e_bs(rhs.e_bs), e_coul(rhs.e_coul), e_rep(rhs.e_rep), e_tot(rhs.e_tot)
{
  interface_base::operator=(rhs);
}

energy::interface_base * energy::interfaces::dftb::sysCallInterface::clone(coords::Coordinates * coord_object) const
{
  sysCallInterface * tmp = new sysCallInterface(*this, coord_object);
  return tmp;
}

energy::interface_base * energy::interfaces::dftb::sysCallInterface::move(coords::Coordinates * coord_object)
{
  sysCallInterface * tmp = new sysCallInterface(*this, coord_object);
  return tmp;
}

void energy::interfaces::dftb::sysCallInterface::swap(interface_base &rhs)
{
  swap(dynamic_cast<sysCallInterface&>(rhs));
}

void energy::interfaces::dftb::sysCallInterface::swap(sysCallInterface &rhs)
{
  interface_base::swap(rhs);
}

energy::interfaces::dftb::sysCallInterface::~sysCallInterface(void)
{

}

/*
Energy class functions that need to be overloaded
*/

// Energy function
double energy::interfaces::dftb::sysCallInterface::e(void)
{
  
  return e_tot;
}

// Energy+Gradient function
double energy::interfaces::dftb::sysCallInterface::g(void)
{
  
  return e_tot;
}

// Hessian function
double energy::interfaces::dftb::sysCallInterface::h(void)
{
  

  return e_tot;
}

// Optimization
double energy::interfaces::dftb::sysCallInterface::o(void)
{
    
  return e_tot;
}

// Output functions
void energy::interfaces::dftb::sysCallInterface::print_E(std::ostream &S) const
{
  S << "Total Energy:      ";
  S << std::right << std::setw(16) << std::fixed << std::setprecision(8) << e_tot;
}

void energy::interfaces::dftb::sysCallInterface::print_E_head(std::ostream &S, bool const endline) const
{
  S << "Energies\n";
  S << std::right << std::setw(24) << "E_bs";
  S << std::right << std::setw(24) << "E_coul";
  S << std::right << std::setw(24) << "E_lr";
  S << std::right << std::setw(24) << "E_rep";
  S << std::right << std::setw(24) << "SUM\n\n";
}

void energy::interfaces::dftb::sysCallInterface::print_E_short(std::ostream &S, bool const endline) const
{
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << e_bs;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << e_coul;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << e_lr;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << e_rep;
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << e_tot << '\n';
  S << "\n";
}

void energy::interfaces::dftb::sysCallInterface::print_G_tinkerlike(std::ostream &S, bool const) const 
{ 
  S << " Cartesian Gradient Breakdown over Individual Atoms :" << std::endl << std::endl;
  S << "  Type      Atom              dE/dX       dE/dY       dE/dZ          Norm" << std::endl << std::endl;
  for(std::size_t k=0; k < coords->size(); ++k)
  {
    S << " Anlyt";
    S << std::right << std::setw(10) << k+1U;
    S << "       ";
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).x();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).y();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).z();
    S << std::right << std::fixed << std::setw(12) << std::setprecision(4);
    S << std::sqrt(
      coords->g_xyz(k).x() * coords->g_xyz(k).x()
    + coords->g_xyz(k).y() * coords->g_xyz(k).y()
    + coords->g_xyz(k).z() * coords->g_xyz(k).z()) << std::endl;
  }
}

void energy::interfaces::dftb::sysCallInterface::to_stream(std::ostream&) const { }

bool energy::interfaces::dftb::sysCallInterface::check_bond_preservation(void) const
{
  std::size_t const N(coords->size());
  for (std::size_t i(0U); i < N; ++i)
  { // cycle over all atoms i
    if (!coords->atoms(i).bonds().empty())
    {
      std::size_t const M(coords->atoms(i).bonds().size());
      for (std::size_t j(0U); j < M && coords->atoms(i).bonds(j) < i; ++j)
      { // cycle over all atoms bound to i
        double const L(geometric_length(coords->xyz(i) - coords->xyz(coords->atoms(i).bonds(j))));
        if (L > 2.2) return false;
      }
    }
  }
  return true;
}

std::vector<coords::float_type>
energy::interfaces::dftb::sysCallInterface::charges() const
{
  std::vector<coords::float_type> charges;
  return charges;
}
#endif