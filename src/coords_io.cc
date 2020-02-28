#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstddef>
#include <cstdio>
#include "atomic.h"
#include "configuration.h"
#include "coords_io.h"
#include "Scon/scon_utility.h"
#include "helperfunctions.h"

#if defined(_MSC_VER) && !defined(CAST_SSCANF_COORDS_IO)
#define CAST_SSCANF_COORDS_IO sscanf_s
#elif !defined(CAST_SSCANF_COORDS_IO)
#define CAST_SSCANF_COORDS_IO sscanf
#endif

struct TinkerCoordFileLine
{
  std::string data;
public:
  friend std::istream& operator>> (std::istream& is, TinkerCoordFileLine& l)
  {
    std::getline(is, l.data);
    l.data.erase(std::remove(l.data.begin(), l.data.end(), '\r'), l.data.end());
    return is;
  }
};


/*! Creates new coords::input::format
 *
 * Creates new coords::input::format
 * according to the input type specified
 * in the configfile
 */
coords::input::format* coords::input::new_format(void)
{
  switch (Config::get().general.input)
  {
  case config::input_types::TINKER:
    //TINKER
    return new formats::tinker;
    break;
  case config::input_types::AMBER:
    //AMBER
    return new formats::amber;
    break;
  case config::input_types::XYZ:
    //XYZ
    return new formats::xyz;
    break;
  case config::input_types::PDB:
    //PDB
    return new formats::pdb;
    break;
  default:
  {
    return new formats::tinker;
  }
  }
}

/*! Creates new coords::input::format
*
* Creates new coords::input::format
* according to the input type specified
* in the configfile for interface creation
*/
coords::input::format* coords::input::additional_format(void)
{
  switch (Config::get().interfcrea.icfiletype)
  {
  case config::input_types::TINKER:
    //TINKER
    return new formats::tinker;
    break;
  case config::input_types::AMBER:
    //AMBER
    return new formats::amber;
    break;
  case config::input_types::XYZ:
    //XYZ
    return new formats::xyz;
    break;
  case config::input_types::PDB:
    //PDB
    return new formats::pdb;
    break;
  default:
  {
    return new formats::tinker;
  }
  }
}

coords::input::format* coords::input::new_interf_format(void)
{
  return new formats::tinker;
}


std::ostream& coords::operator<< (std::ostream& stream, coords::Coordinates const& coord)
{
  if (Config::get().general.output == config::output_types::TINKER)
    stream << coords::output::formats::tinker(coord);
  else if (Config::get().general.output == config::output_types::XYZ)
    stream << coords::output::formats::xyz(coord);
  else if (Config::get().general.output == config::output_types::XYZ)
    stream << coords::output::formats::xyz_dftb(coord);
  else if (Config::get().general.output == config::output_types::MOLDEN)
    stream << coords::output::formats::moldenxyz(coord);
  else if (Config::get().general.output == config::output_types::ZMATRIX)
    stream << coords::output::formats::zmatrix(coord);
  stream << std::flush;
  return stream;
}

/*! Read coordinates from a tinker .arc file
 *
 * This function is used to read coordinates from
 * a tinker xyz file (.arc or sometimes .xyz)
 *
 * @return: Coordinates object containing the coordinates
 * @param file: Filename of the .arc tinker file
 */
coords::Coordinates coords::input::formats::tinker::read(std::string file) {
  // Create empty coordinates object!
  Coordinates coord_object;
  std::ifstream coord_file_stream(file.c_str(), std::ios_base::in);
  std::size_t number_of_structures_before{ input_ensemble.size() };  // if this is not zero the new structure are just added

  if (coord_file_stream) {
    std::size_t N(0U);
    std::string line;
    std::getline(coord_file_stream, line);
    std::istringstream first_line_stream(line);
    first_line_stream >> N;
    // coord_object.m_topology.resize(N);
    Atoms atoms;
    if (N == 0U)
      throw std::logic_error("Error in reading tinker (arc) coordinate frame: Zero atoms present. It seems like the file is not valid (Reading file at: '" + file +
        "').");
    Representation_3D positions;
    std::vector<std::size_t> index_of_atom(N);
    bool indexation_not_contiguous(false), has_in_out_subsystems(false);


    // loop fetching atoms and positions
    for (std::size_t i(1U); std::getline(coord_file_stream, line); ++i) {
      if (i <= N) {
        std::istringstream linestream(line);
        std::size_t const nbmax(7u);
        tinker::line tfl;
        linestream >> tfl;
        Atom atom(tfl.symbol);
        index_of_atom[tfl.index - 1] = positions.size();
        if (positions.size() != (tfl.index - 1)) {
          indexation_not_contiguous = true;
        }
        for (std::size_t j(0u); j < nbmax && tfl.bonds[j] > 0u; ++j) {
          atom.bind_to(tfl.bonds[j] - 1u);
          // coord_object.topo((tfl.bonds[j] - 1u), i - 1);
        }
        positions.push_back(tfl.position);
        if (scon::find_substr_ci(line, "in") != std::string::npos) {
          atom.set_sub_type(Atom::ST_IN);
          atom.assign_to_system(1u);
          has_in_out_subsystems = true;
        }
        else if (scon::find_substr_ci(line, "out") != std::string::npos) {
          atom.set_sub_type(Atom::ST_OUT);
          atom.assign_to_system(2u);
          has_in_out_subsystems = true;
        }
        atom.set_energy_type(tfl.tinker_type);
        atoms.add(atom);
        if (i == N) {
          input_ensemble.push_back(positions);
          positions.clear();
        }
      }
      else {
        if (i % (N + 1u) != 0) {
          double x(0), y(0), z(0);
          CAST_SSCANF_COORDS_IO(line.c_str(), "%*u %*s %lf %lf %lf", &x, &y,
            &z);
          positions.emplace_back(x, y, z);
          if ((i - (input_ensemble.size()-number_of_structures_before) * (N + 1u)) ==
            N) { // if we are at the end of a structure
            if (positions.size() != atoms.size())
              throw std::logic_error(
                "The size of an additionally provided"
                " structure does not match the number of atoms.");
            input_ensemble.push_back(positions);
            positions.clear();
          }
        }
      }
    }

    // dividing subsystems
    if (!has_in_out_subsystems) {
      auto n_susy = Config::get().coords.subsystems.size();
      for (std::size_t i = 0; i < n_susy; ++i) {
        for (auto a : Config::get().coords.subsystems[i]) {
          if (0 < a && a < (N - 1u)) {
            atoms.atom(a - 1u).assign_to_system(i + 1u);
          }
        }
      }
    }

    if (indexation_not_contiguous) {
      std::cout << "Indexation not contiguous. Rebinding atoms.\n";
      for (std::size_t i(0U); i < atoms.size(); ++i) {
        for (auto bonding_partner : atoms.atom(i).bonds()) {
          std::cout << "Partner of atom which is now " << i + 1
            << " detached from " << bonding_partner
            << " and rebound to " << index_of_atom[bonding_partner]
            << '\n';
          atoms.atom(i).detach_from(bonding_partner);
          atoms.atom(i).bind_to(index_of_atom[bonding_partner]);
        }
      }
    }

    if (input_ensemble.empty())
      throw std::logic_error("No structures found.");
    coords::PES_Point x(input_ensemble[0u]);
    if (!Config::get().coords.fixed.empty()) {
      for (auto fix : Config::get().coords.fixed) {
        if (fix < atoms.size())
          atoms.atom(fix).fix(true);
      }
    }
    coord_object.init_swap_in(atoms, x);
    for (auto& p : input_ensemble) {
      p.gradient.cartesian.resize(p.structure.cartesian.size());
      coord_object.set_xyz(p.structure.cartesian, true);
      coord_object.to_internal_light();
      p = coord_object.pes();
    }

    if (Config::get().general.chargefile) // read charges from chargefile
    {
      std::vector<double> charges;
      std::ifstream charge_stream("charges.txt", std::ios_base::in);
      std::string read;
      while (charge_stream) {
        if (charge_stream >> read) {
          // ignore atom number
        }
        if (charge_stream >> read) {
          // ignore atom type
        }
        if (charge_stream >> read) {
          charges.push_back(std::stod(read));
        }
      }
      if (charges.size() == coord_object.size()) {
        Config::set().coords.atom_charges = charges;
        if (Config::get().general.verbosity > 3) {
          std::cout << "Reading charges from chargefile successful.\n";
        }
      }
      else // not the correct number of charges in chargefile
      {
        std::cout << "Reading charges from chargefile failed.\n";
        throw std::logic_error("Reading the structure input file failed.");
      }
    }

  }
  else
    throw std::logic_error("Reading the structure from file '" + file + "' failed.");
  return coord_object;
}


static void tinker_dummy_to_stream(std::ostream& stream, std::size_t index, coords::float_type x, coords::float_type y, coords::float_type z)
{
  stream << std::right << std::setw(6) << index << "  ";
  stream << std::left << "XX ";
  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << x;
  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << y;
  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << z;
  stream << std::right << std::setw(6) << 0 << "\n";
}

void coords::output::formats::tinker::to_stream(std::ostream& stream) const
{
  bool const pp = Config::get().periodics.periodic && Config::get().periodics.periodic_print;
  std::size_t const N(ref.size());
  stream << (pp ? N + 8 : N) << '\n';
  //std::size_t index_width(1), tens(N);
  //while(tens >= 10U)
  //{
  //  tens /= 10U;
  //  ++index_width;
  //}
  for (std::size_t i(0U); i < N; ++i)
  {
    stream << std::right << std::setw(6) << i + 1U << "  ";
    stream << std::left << std::setw(3) << ref.atoms(i).symbol().substr(0U, 2U);
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z();
    stream << std::right << std::setw(6) << ref.atoms(i).energy_type();
    std::size_t const bSize(ref.atoms(i).bonds().size());
    for (std::size_t j(0U); j < bSize; ++j)
    {
      stream << std::right << std::setw(6) << ref.atoms(i).bonds()[j] + 1U;
    }
    if (ref.atoms(i).sub_type() == coords::Atom::sub_types::ST_IN) stream << " IN";
    else if (ref.atoms(i).sub_type() == coords::Atom::sub_types::ST_OUT) stream << " OUT";
    stream << '\n';
  }
  if (pp)
  {
    coords::Cartesian_Point const halfbox = Config::get().periodics.pb_box / 2.;
    tinker_dummy_to_stream(stream, N + 1, halfbox.x(), halfbox.y(), -halfbox.z());
    tinker_dummy_to_stream(stream, N + 2, -halfbox.x(), -halfbox.y(), -halfbox.z());
    tinker_dummy_to_stream(stream, N + 3, halfbox.x(), -halfbox.y(), -halfbox.z());
    tinker_dummy_to_stream(stream, N + 4, -halfbox.x(), halfbox.y(), -halfbox.z());
    tinker_dummy_to_stream(stream, N + 5, halfbox.x(), halfbox.y(), halfbox.z());
    tinker_dummy_to_stream(stream, N + 6, -halfbox.x(), -halfbox.y(), halfbox.z());
    tinker_dummy_to_stream(stream, N + 7, halfbox.x(), -halfbox.y(), halfbox.z());
    tinker_dummy_to_stream(stream, N + 8, -halfbox.x(), halfbox.y(), halfbox.z());
  }
}


void coords::output::formats::moldenxyz::to_stream(std::ostream& stream) const
{
  std::size_t const N(ref.size());
  stream << N << '\n';
  stream << "Energy = " << ref.energyinterface()->energy;
  for (std::size_t i(0U); i < N; ++i)
  {
    stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z();
    stream << '\n';
  }
}

void coords::output::formats::xyz::to_stream(std::ostream& stream) const
{
  std::size_t const N(ref.size());
  for (std::size_t i(0U); i < N; ++i)
  {
    stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];//Made the precision from 6 to 7 according to Lee-Pings Code
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << ref.xyz(i).z();
    stream << '\n';
  }
}

void coords::output::formats::xyz_cast::to_stream(std::ostream& stream) const
{
  std::size_t const N(ref.size());
  stream << std::to_string(N) << "\nCreated_Using_CAST\n";
  for (std::size_t i(0U); i < N; ++i)
  {
    stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];//Made the precision from 6 to 7 according to Lee-Pings Code
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << ref.xyz(i).z();
    stream << '\n';
  }
}

/**writes the geometry of the coord object in the gen format (description of format see manual of dftb+, appendix C)*/
void coords::output::formats::xyz_gen::to_stream(std::ostream& stream) const
{
  // create a vector with all element symbols of the input structure 
  std::vector<std::string> existing_symbols;
  for (auto a : ref.atoms())
  {
    if (is_in(a.symbol(), existing_symbols) == false)
    {
      existing_symbols.push_back(a.symbol());
    }
  }

  // write structure
  std::size_t const N(ref.size());
  if (Config::get().periodics.periodic) stream << N << "  S\n";  // supercell in cartesian coordinates
  else stream << N << "  C\n";                                   // cluster (non-periodic)
  for (auto s : existing_symbols)
  {
    stream << s << " ";
  }
  stream << "\n";
  for (std::size_t i(0U); i < N; ++i)      // atomic coordinates
  {
    stream << std::left << std::setw(5) << i + 1 << std::left << std::setw(5) << find_index(ref.atoms(i).symbol(), existing_symbols) + 1;
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z();
    stream << '\n';
  }

  if (Config::get().periodics.periodic)   // stuff for periodic cell
  {
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0;   // coordinate origin (is ignored by parser)
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0;
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0;
    stream << '\n';

    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << Config::get().periodics.pb_box.x(); // lattice vector in x-direction
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0;
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0;
    stream << '\n';

    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0; // lattice vector in y-direction
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << Config::get().periodics.pb_box.y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0;
    stream << '\n';

    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0; // lattice vector in z-direction
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0;
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << Config::get().periodics.pb_box.z();
    stream << '\n';
  }
}

void coords::output::formats::xyz_dftb::to_stream(std::ostream& stream) const
{
  std::size_t const N(ref.size());
  stream << N << '\n';
  if (Config::get().energy.dftbaby.charge == 0)
  {
    stream << Config::get().general.inputFilename << "\n";
  }
  else
  {
    stream << "charge=" + std::to_string(Config::get().energy.dftbaby.charge) << "\n";
  }

  for (std::size_t i(0U); i < N; ++i)
  {
    stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z();
    stream << '\n';
  }
}



void coords::output::formats::zmatrix::to_stream(std::ostream& stream) const
{
  std::size_t const N(ref.size());
  if (stream.good() && N > 0)
  {
    stream << "zmat angstroms\n";
    stream << std::left << std::setw(10) << 1U;
    stream << std::left << std::setw(10) << ref.atoms(0u).i_to_a() + 1;
    stream << std::left << std::setw(4) << ref.atoms(ref.atoms(0u).i_to_a()).symbol() << '\n';
    if (N > 1)
    {
      stream << std::left << std::setw(10) << 2U;
      std::size_t const A = ref.atoms(1u).i_to_a();
      stream << std::left << std::setw(10) << A + 1;
      stream << std::left << std::setw(4) << ref.atoms(A).symbol();
      stream << std::right << std::setw(10) << ref.atoms(1u).ibond() + 1;
      stream << std::right << std::setw(10) << "bnd" << 1u;
      stream << '\n';
    }
    if (N > 2)
    {
      stream << std::left << std::setw(10) << 3U;
      std::size_t const A = ref.atoms(2u).i_to_a();
      stream << std::left << std::setw(10) << A + 1;
      stream << std::left << std::setw(4) << ref.atoms(A).symbol();
      stream << std::right << std::setw(10) << ref.atoms(2u).ibond() + 1;
      stream << std::right << std::setw(10) << "bnd" << 2u;
      stream << std::right << std::setw(10) << ref.atoms(2u).iangle() + 1;
      stream << std::right << std::setw(10) << "ang" << 2u;
      stream << '\n';
    }

    for (std::size_t i = 3; i < N; ++i)
    {
      stream << std::left << std::setw(10) << i + 1U;
      std::size_t const A = ref.atoms(i).i_to_a();
      stream << std::left << std::setw(10) << A + 1;
      stream << std::left << std::setw(4) << ref.atoms(A).symbol();
      stream << std::right << std::setw(10) << ref.atoms(i).ibond() + 1;
      stream << std::right << std::setw(10) << "bnd" << i;
      stream << std::right << std::setw(10) << ref.atoms(i).iangle() + 1;
      stream << std::right << std::setw(10) << "ang" << i;
      stream << std::right << std::setw(10) << ref.atoms(i).idihedral() + 1;
      stream << std::right << std::setw(10) << "dih" << i;
      stream << '\n';
    }
    stream << "variables\n";
    if (N > 1)
    {
      stream << "bnd" << std::left << std::setw(10) << 1u;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(1u).radius();
      stream << '\n';
    }
    if (N > 2)
    {
      stream << "bnd" << std::left << std::setw(10) << 2u;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(2u).radius();
      stream << '\n';
      stream << "ang" << std::left << std::setw(10) << 2u;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(2u).inclination();
      stream << '\n';
    }
    for (std::size_t i = 3; i < N; ++i)
    {
      stream << "bnd" << std::left << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(i).radius();
      stream << '\n';
      stream << "ang" << std::left << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(i).inclination();
      stream << '\n';
      stream << "dih" << std::left << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(i).azimuth();
      stream << '\n';
    }
    stream << "constants\n";
    for (std::size_t i = 1; i < N; ++i)
    {
      stream << "g_bnd" << std::left << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.g_intern(i).x();
      stream << '\n';
      if (i > 1)
      {
        stream << "g_ang" << std::left << std::setw(10) << i;
        stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.g_intern(i).y();
        stream << '\n';
      }
      if (i > 2)
      {
        stream << "g_dih" << std::left << std::setw(10) << i;
        stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.g_intern(i).z();
        stream << '\n';
      }
    }
    stream << "end\n";
  }
  else throw std::runtime_error("ERR_FILE_WRITE: stream bad");
}

void coords::output::formats::xyz_mopac::to_stream(std::ostream& stream) const
{
  std::size_t const N(ref.size());
  for (std::size_t i(0U); i < N; ++i)
  {
    if (ref.atoms(i).fixed()) {
      stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x() << " 0";
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y() << " 0";
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z() << " 0";
      stream << '\n';
    }
    else
    {
      stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x() << " 1";
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y() << " 1";
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z() << " 1";
      stream << '\n';
    }
  }
}

void coords::output::formats::gaussview::to_stream(std::ostream &gstream) const
{
  if (Config::get().general.energy_interface == config::interface_types::ONIOM || Config::get().general.energy_interface == config::interface_types::QMMM)
  {
    gstream << "# ONIOM(HF/6-31G:UFF)\n\n";                      // method
    gstream << "some title\n\n";
    gstream << "0 1 0 1 0 1\n";                                  // charge and multiplicity
    gstream << std::fixed << std::setprecision(3);
    for (auto i = 0u; i < ref.size(); ++i)                    // atom information
    {
      auto symbol = ref.atoms(i).symbol();
      auto x = ref.xyz(i).x();
      auto y = ref.xyz(i).y();
      auto z = ref.xyz(i).z();
      auto system = "L";
      if (is_in(i, Config::get().energy.qmmm.qm_systems[0])) system = "H";
      gstream << symbol << "\t" << x << "\t" << y << "\t" << z << "\t" << system << "\n";
    }
    gstream << "\n";
  }
  else if (Config::get().general.energy_interface == config::interface_types::THREE_LAYER)
  {
    gstream << "# ONIOM(B3LYP/6-31G:HF/STO-3G:UFF)\n\n";                      // method
    gstream << "some title\n\n";
    gstream << "0 1 0 1 0 1 0 1 0 1\n";                                       // charge and multiplicity
    gstream << std::fixed << std::setprecision(3);
    for (auto i = 0u; i < ref.size(); ++i)                    // atom information
    {
      auto symbol = ref.atoms(i).symbol();
      auto x = ref.xyz(i).x();
      auto y = ref.xyz(i).y();
      auto z = ref.xyz(i).z();
      auto system = "L";
      if (is_in(i, Config::get().energy.qmmm.qm_systems[0])) system = "H";
      else if (is_in(i, Config::get().energy.qmmm.seatoms)) system = "M";
      gstream << symbol << "\t" << x << "\t" << y << "\t" << z << "\t" << system << "\n";
    }
    gstream << "\n";
  }
  else  // input for "normal molecule"
  {
    gstream << "# HF/6-31G\n\n";                      // method
    gstream << "some title\n\n";
    gstream << "0 1\n";                               // charge and multiplicity
    gstream << coords::output::formats::xyz(ref);  // coordinates
    gstream << "\n";
  }
}
