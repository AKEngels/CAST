#include "energy_int_terachem.h"

#if defined(USE_MPI)

#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <climits>

#include "atomic.h"

#include "filemanipulation.h"
#include "configuration.h"
#include "coords.h"

/*
  Mopac sysCall functions
*/
cmpi::CMPI energy::interfaces::terachem::mpiInterface::mpo;

bool energy::interfaces::terachem::mpiInterface::option_init_done(false);

void energy::interfaces::terachem::mpiInterface::init_terachem (void)
{
  if(Config::get().general.verbosity >= 3) std::cout << "Reading TeraChem options from 'CAST_TERACHEM_OPTIONS.txt'.\n";
  LBL_FileReader options("CAST_TERACHEM_OPTIONS.txt");
  std::string option_string;
  for (std::size_t i(0U); i < options.data.size() && i < 100; i+=2)
  {
    if (i < options.data.size()-1)
    {
      options.data[i].resize(128, ' ');
      option_string += options.data[i];
      options.data[i+1].resize(128, ' ');
      option_string += options.data[i+1];
    }
    //option_string += '\n';
  }
  option_string += "end";
  option_string.resize(option_string.length()+253, ' ');
  //option_string += '\n';
  std::vector<char> option_buffer(option_string.begin(), option_string.end());
  mpiInterface::mpo.send<char>(option_buffer, 0, 2);
  if (Config::get().general.verbosity >= 4)
  {
    std::cout << "OptBuffer x " << option_buffer.size() << ":";
    for (auto p : option_buffer) std::cout << p;
    std::cout << '\n';
    std::cout << (options.data.size()/2) << " options transfered to TeraChem.\n";
  }
  mpiInterface::option_init_done = true;
}


energy::interfaces::terachem::mpiInterface::mpiInterface (coords::Coordinates * cp)
: energy::interface_base(cp), qm_population_charges(cp->size())
{
  if (!mpiInterface::mpo.connected())
  {
    mpiInterface::mpo.connect("terachem_port");
  }
  if (!mpiInterface::option_init_done)
  {
    mpiInterface::init_terachem();
    std::size_t const N(coords->size());
    // Check and convert number of atoms
    if (N > 10000U || N > INT_MAX)
    {
      std::runtime_error("ERR_ENERGY_TERACHEM_MAXIMUM_ATOMS_EXCEEDED");
    }
    int num_atoms(static_cast<int>(N));
     //! atom symbol buffer
    std::vector<char> char_buffer(2U*N+1, 0U);
    // build atom-symbol C-style-string
    for (std::size_t i=0; i<N; ++i)
    {
      std::string const & smb = atomic::symbolMap[coords->atoms(i).number()];
      std::size_t const Q(i*2U);
      char_buffer[Q] = smb[0];
      if (smb.length() > 1) char_buffer[Q+1] = smb[1];
    }
    // send number of atoms and atom symbols for initial orbital checks
    if (Config::get().general.verbosity >= 4) std::cout << "Sending initial number of atoms to TeraChem.\n";
    mpiInterface::mpo.send(num_atoms,  0, 2);
    if (Config::get().general.verbosity >= 4) std::cout << "Sending initial atom symbols to TeraChem.\n";
    mpiInterface::mpo.send(char_buffer, 0, 2);
  }
  optimizer = true;
}
energy::interfaces::terachem::mpiInterface::~mpiInterface (void)
{
  mpi_send_data(0);
}


energy::interfaces::terachem::mpiInterface::mpiInterface (mpiInterface const & rhs, coords::Coordinates *cobj)
  : interface_base(cobj), qm_population_charges(rhs.qm_population_charges)
{
  interface_base::operator=(rhs);
  qm_population_charges.resize(cobj->size());
}

energy::interfaces::terachem::mpiInterface::mpiInterface (mpiInterface && rhs, coords::Coordinates *cobj)
  : interface_base(cobj), qm_population_charges(std::move(rhs.qm_population_charges))
{
  interface_base::swap(rhs);
  qm_population_charges.resize(coords->size());
}

energy::interface_base * energy::interfaces::terachem::mpiInterface::clone (coords::Coordinates * coord_object) const
{
  mpiInterface * tmp = new mpiInterface(*this, coord_object);
  return tmp;
}

energy::interface_base * energy::interfaces::terachem::mpiInterface::move (coords::Coordinates * coord_object)
{
  mpiInterface * tmp = new mpiInterface(std::move(*this), coord_object);
  return tmp;
}

void energy::interfaces::terachem::mpiInterface::swap (interface_base &rhs)
{
  swap(dynamic_cast<mpiInterface&>(rhs));
}

void energy::interfaces::terachem::mpiInterface::swap (mpiInterface &rhs)
{
  interface_base::swap(rhs);
  qm_population_charges.swap(rhs.qm_population_charges);
}

// MPI send calculation request
void energy::interfaces::terachem::mpiInterface::mpi_send_data (const int tag)
{
  std::size_t const N(coords->size());
  // Check and convert number of atoms
  if (N > 10000U || N > INT_MAX) std::runtime_error("ERR_ENERGY_TERACHEM_MAXIMUM_ATOMS_EXCEEDED");
  int num_atoms(static_cast<int>(N));
   //! atom symbol buffer
  std::vector<char> char_buffer(2U*N+1, ' ');
  char_buffer.back() = '\0';
  //! atom coordinates buffer
  std::vector<double> coords_buffer(N*3U, 0.0);
  // build atom-symbol C-style-string and coord buffer
  for (std::size_t i=0; i<N; ++i)
  {
    std::string const & smb = atomic::symbolMap[coords->atoms(i).number()];
    std::size_t const P(i*3U), Q(i*2U);
    char_buffer[Q] = smb[0];
    if (smb.length() > 1) char_buffer[Q+1] = smb[1];
    coords_buffer[P] = coords->xyz(i).x();
    coords_buffer[P+1U] = coords->xyz(i).y();
    coords_buffer[P+2U] = coords->xyz(i).z();
  }
  // send number of atoms
  if (Config::get().general.verbosity >= 4) std::cout << "Sending Number of atoms to TeraChem.\n";
  mpiInterface::mpo.send(num_atoms,  0, tag);
  // send atom symbols
  if (Config::get().general.verbosity >= 4) std::cout << "Sending atom symbols to TeraChem.\n";
  mpiInterface::mpo.send(char_buffer, 0, tag);
  // send coordinates
  if (Config::get().general.verbosity >= 4) std::cout << "Sending coords to TeraChem.\n";
  //for (auto coord : coords_buffer) std::cout << coord << '\n';
  mpiInterface::mpo.send(coords_buffer, 0, tag);
  qm_population_charges.resize(N);
}

// MPI recieve calculated data
void energy::interfaces::terachem::mpiInterface::mpi_recv_energy (void)
{
  if (Config::get().general.verbosity >= 4) std::cout << "Recieving energy from TeraChem.\n";
  mpiInterface::mpo.recv(energy, MPI_ANY_SOURCE, MPI_ANY_TAG);
  if (Config::get().general.verbosity >= 4) std::cout << "Recieved energy with tag " << mpiInterface::mpo.status().MPI_TAG << " from TeraChem.\n";
  if (mpiInterface::mpo.status().MPI_TAG == 1)
  {
    integrity = false;
    return;
  }
  if (Config::get().general.verbosity >= 4) std::cout << "Recieving QM population charges from TeraChem.\n";
  mpiInterface::mpo.recv(qm_population_charges, MPI_ANY_SOURCE, MPI_ANY_TAG);
  if (Config::get().general.verbosity >= 4) std::cout << "Recieved QM population charges from TeraChem.\n";
  std::vector<double> buffer(4);
  if (Config::get().general.verbosity >= 4) std::cout << "Recieving QM dipole from TeraChem.\n";
  mpiInterface::mpo.recv(buffer, MPI_ANY_SOURCE, MPI_ANY_TAG);
  if (Config::get().general.verbosity >= 4) std::cout << "Recieving MM dipole from TeraChem.\n";
  mpiInterface::mpo.recv(buffer, MPI_ANY_SOURCE, MPI_ANY_TAG);
  if (Config::get().general.verbosity >= 4) std::cout << "Recieving TOT dipole from TeraChem.\n";
  mpiInterface::mpo.recv(buffer, MPI_ANY_SOURCE, MPI_ANY_TAG);
}

void energy::interfaces::terachem::mpiInterface::mpi_recv_gradients (void)
{
  std::size_t const N = coords->size();
  std::vector<double> grad_buffer(N*3U, 0.0);
  if (Config::get().general.verbosity >= 4) std::cout << "Recieving gradients.\n";
  mpiInterface::mpo.recv(grad_buffer, MPI_ANY_SOURCE, MPI_ANY_TAG);
  if (mpiInterface::mpo.status().MPI_TAG == 1) throw std::runtime_error("ERR_ENERGY_TERACHEM : SCF did not converge");
  for (std::size_t i=0; i<N; ++i)
  {
    std::size_t const P(i*3U);
    coords->update_g_xyz(i, coords::Cartesian_Point(grad_buffer[P]*(energy::Hartree_Bohr2Kcal_MolAng), grad_buffer[P+1U]*(energy::Hartree_Bohr2Kcal_MolAng), grad_buffer[P+2U]*(energy::Hartree_Bohr2Kcal_MolAng)));
  }
}

void energy::interfaces::terachem::mpiInterface::mpi_recv_positions (void)
{
  std::size_t const N = coords->size();
  std::vector<double> pos_buffer(N*3U, 0.0);
  if (Config::get().general.verbosity >= 4) std::cout << "Recieving positions.\n";
  mpiInterface::mpo.recv(pos_buffer, MPI_ANY_SOURCE, MPI_ANY_TAG);
  for (std::size_t i(0U); i<N; ++i)
  {
    std::size_t const P(i*3U);
    coords->move_atom_to(i, coords::Cartesian_Point(pos_buffer[P], pos_buffer[P+1U], pos_buffer[P+2U]));
  }
}


/*
  Energy class functions that need to be overloaded

*/

// Energy function
coords::float_type energy::interfaces::terachem::mpiInterface::e (void)
{
  integrity = true;
  if (Config::get().general.verbosity >= 4) std::cout << "Data transfer for engery call to -> TeraChem.\n";
  mpi_send_data(2);
  if (Config::get().general.verbosity >= 4) std::cout << "Data transfer 1 for engery call from <- TeraChem.\n";
  mpi_recv_energy();
  if (Config::get().general.verbosity >= 4) std::cout << "Data transfer 2 for engery call from <- TeraChem.\n";
  mpi_recv_gradients();
  return energy*=energy::au2kcal_mol;
}

// Energy+Gradient function
coords::float_type energy::interfaces::terachem::mpiInterface::g (void)
{
  integrity = true;
  if (Config::get().general.verbosity >= 4) std::cout << "Data transfer for gradient call to -> TeraChem.\n";
  mpi_send_data(2);
  if (Config::get().general.verbosity >= 4) std::cout << "Data transfer 1 for gradient call from <- TeraChem.\n";
  mpi_recv_energy();
  if (Config::get().general.verbosity >= 4) std::cout << "Data transfer 2 for gradient call from <- TeraChem.\n";
  mpi_recv_gradients();
  return energy*=energy::au2kcal_mol;
}

// Energy+Gradient+Hessian function
coords::float_type energy::interfaces::terachem::mpiInterface::h (void)
{
  throw std::runtime_error("NO HESSIAN YET!");
  return 0.0;
}

// Optimization
coords::float_type energy::interfaces::terachem::mpiInterface::o (void)
{
  integrity = true;
  if (Config::get().general.verbosity >= 4) std::cout << "Data transfer for opt call to -> TeraChem.\n";
  mpi_send_data(3);
  if (Config::get().general.verbosity >= 4) std::cout << "Data transfer 1 (recv energy) for opt call from <- TeraChem.\n";
  mpi_recv_energy();
  //if (Config::get().general.verbosity >= 4) std::cout << "Data transfer 2 (recv gradients) for opt call from <- TeraChem.\n";
  //mpi_recv_gradients();
  if (Config::get().general.verbosity >= 4) std::cout << "Data transfer 3 (recv positions) for opt call from <- TeraChem.\n";
  mpi_recv_positions();
  return energy*=energy::au2kcal_mol;
}

// Output functions
void energy::interfaces::terachem::mpiInterface::print_E (std::ostream &S) const
{
  S << "Total Energy: ";
  S << std::right << std::setw(26) << std::fixed << std::setprecision(8) << energy << std::endl;
}

void energy::interfaces::terachem::mpiInterface::print_E_head (std::ostream &S, bool const endline) const
{
  if(!coords->potentials().empty()) S << std::right << std::setw(26) << "BIAS";
  S << std::right << std::setw(26) << "E";
  if (endline) S << std::endl;
}

void energy::interfaces::terachem::mpiInterface::print_E_short (std::ostream &S, bool const endline) const
{
  if(!coords->potentials().empty()) S << std::right << std::setw(26) << std::scientific << std::setprecision(5) << coords->potentials().energy();
  S << std:: right << std::setw(26) << std::scientific << std::setprecision(5) << energy;
  if (endline) S << std::endl;
  S << std::fixed;
}

void energy::interfaces::terachem::mpiInterface::update (bool const) { }

void energy::interfaces::terachem::mpiInterface::to_stream (std::ostream &S) const
{
  interface_base::to_stream(S);
}

#endif
