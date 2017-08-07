#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <cstdlib>
#include <utility>
#include "atomic.h"
#include "energy_int_dftb.h"
#include "configuration.h"
#include "coords.h"
#include "coords_io.h"
#if defined (_MSC_VER)
#include "win_inc.h"
#endif

#ifdef _MSC_VER
#pragma warning (disable: 4996)
#endif

/*
dftb sysCall functions
*/

energy::interfaces::dftb::sysCallInterface::sysCallInterface(coords::Coordinates * cp) :
  energy::interface_base(cp),
  hof_kcal_mol(0.0), hof_kj_mol(0.0), e_total(0.0),
  e_electron(0.0), e_core(0.0), id(Config::get().general.outputFilename), failcounter(0u)
{
  std::stringstream ss;
  ss << (std::size_t(std::rand()) | (std::size_t(std::rand()) << 15));
  id.append("_tmp_").append(ss.str());
  optimizer = true;
}

energy::interfaces::dftb::sysCallInterface::sysCallInterface(sysCallInterface const & rhs, coords::Coordinates *cobj) :
  interface_base(cobj),
  hof_kcal_mol(rhs.hof_kcal_mol), hof_kj_mol(rhs.hof_kj_mol), e_total(rhs.e_total),
  e_electron(rhs.e_electron), e_core(rhs.e_core), id(rhs.id), failcounter(rhs.failcounter)
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
  std::swap(failcounter, rhs.failcounter);
}

energy::interfaces::dftb::sysCallInterface::~sysCallInterface(void)
{
  
}





namespace
{
  int dftb_system_call(std::string const & command_line)
  {

  }
}



/*
Energy class functions that need to be overloaded
*/

// Energy function
double energy::interfaces::dftb::sysCallInterface::e(void)
{
  integrity = true;
  grad_var = false;
  std::cout<<"No energy yet\n";
  return energy;
}

// Energy+Gradient function
double energy::interfaces::dftb::sysCallInterface::g(void)
{
  integrity = true;
  grad_var = true;
  std::cout<<"no gradients yet\n";
  return energy;
}

// Energy+Gradient+Hessian function
double energy::interfaces::dftb::sysCallInterface::h(void)
{
  integrity = true;
  grad_var = false;
  std::cout<<"no hessian yet\n";
  return energy;
}

// Optimization
double energy::interfaces::dftb::sysCallInterface::o(void)
{
  integrity = true;
  grad_var = false;
  std::cout<<"no optimizer yet\n";
  return energy;
}

// Output functions
void energy::interfaces::dftb::sysCallInterface::print_E(std::ostream &S) const
{
  std::cout<<"DFTB energy\n";
}

void energy::interfaces::dftb::sysCallInterface::print_E_head(std::ostream &S, bool const endline) const
{
  std::cout<<"DFTB head\n";
}

void energy::interfaces::dftb::sysCallInterface::print_E_short(std::ostream &S, bool const endline) const
{
  std::cout<<"DFTB short\n";
}

void energy::interfaces::dftb::sysCallInterface::print_G_tinkerlike(std::ostream &, bool const) const { }

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