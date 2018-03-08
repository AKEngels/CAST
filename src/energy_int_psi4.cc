#include "energy_int_psi4.h"


void energy::interfaces::psi4::sysCallInterface::swap(interface_base & other){}
energy::interface_base * energy::interfaces::psi4::sysCallInterface::clone(coords::Coordinates* coord_object) const{ return new sysCallInterface(*this, coord_object); }
energy::interface_base * energy::interfaces::psi4::sysCallInterface::move(coords::Coordinates* coord_object) { return new sysCallInterface(*this, coord_object); }

void energy::interfaces::psi4::sysCallInterface::update(bool const skip_topology){}

coords::float_type energy::interfaces::psi4::sysCallInterface::e(void){}
coords::float_type energy::interfaces::psi4::sysCallInterface::g(void){}
coords::float_type energy::interfaces::psi4::sysCallInterface::h(void){}
coords::float_type energy::interfaces::psi4::sysCallInterface::o(void){}

void energy::interfaces::psi4::sysCallInterface::print_E(std::ostream&) const{}
void energy::interfaces::psi4::sysCallInterface::print_E_head(std::ostream&, bool const endline) const{}
void energy::interfaces::psi4::sysCallInterface::print_E_short(std::ostream&, bool const endline) const{}
void energy::interfaces::psi4::sysCallInterface::print_G_tinkerlike(std::ostream&, bool const aggregate) const{}
void energy::interfaces::psi4::sysCallInterface::to_stream(std::ostream&) const{}
