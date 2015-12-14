#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstddef>
#include "atomic.h"
#include "global.h"
#include "configuration.h"
#include "coords_io.h"
#include "scon_utility.h"

struct TinkerCoordFileLine
{
  std::string data;
public:
  friend std::istream & operator>> (std::istream &is, TinkerCoordFileLine &l)
  {
    std::getline(is, l.data);
    l.data.erase(std::remove(l.data.begin(), l.data.end(), '\r'), l.data.end());
    return is;
  }
};

coords::input::format* coords::input::new_format (void)
{
  switch (Config::get().general.input)
  {
  default:
    {
      return new formats::tinker;
    }
  }
  //return new formats::tinker;
}

coords::Coordinates coords::input::formats::tinker::read (std::string file)
{
  Coordinates coord_object;
  std::ifstream coord_file_stream(file.c_str(), std::ios_base::in);
  if (coord_file_stream)
  {
    std::size_t N(0U);
    std::string line;
    std::getline(coord_file_stream, line);
    std::istringstream first_line_stream(line);
    first_line_stream >> N;
	  //coord_object.m_topology.resize(N);
    Atoms atoms;
    if (N == 0U) throw std::logic_error(std::string("ERR_COORD: Expecting no atoms from ").append(file));
    Representation_3D positions;
    std::vector<std::size_t> index_of_atom(N);
    bool indexation_not_contiguous(false);

    for (std::size_t i(1U); std::getline(coord_file_stream, line); ++i)
    {
      //std::cout << "Line " << i << " mod: " << i%(N+1u) << lineend;
      std::istringstream linestream(line);
      if (i <= N)
      {
        std::size_t const nbmax(7u);
        tinker::line tfl;
        linestream >> tfl;
        Atom atom(tfl.symbol);
        index_of_atom[tfl.index-1] = positions.size();
        if (positions.size() != (tfl.index-1)) 
        {
          indexation_not_contiguous = true;
        }
    		for (std::size_t j(0u); j<nbmax && tfl.bonds[j] > 0u; ++j) 
        { 
          atom.bind_to(tfl.bonds[j] - 1u); 
          //coord_object.topo((tfl.bonds[j] - 1u), i - 1); 
        }
		    positions.push_back(tfl.position);
        if (scon::find_substr_ci(line, "in") != std::string::npos) 
        {  
          atom.set_sub_type(Atom::ST_IN);
          atom.assign_to_system(1u);
        }
        else if (scon::find_substr_ci(line, "out") != std::string::npos)
        {
          atom.set_sub_type(Atom::ST_OUT);
          atom.assign_to_system(2u);
        }
        atom.set_energy_type(tfl.tinker_type);
        atoms.add(atom);
        if (i == N) 
        {
          input_ensemble.push_back(positions);
          positions.clear();
        }
      }
      else
      {
        if (i%(N+1u) != 0)
        {
          tinker::line tfl;
          linestream >> tfl.index >> tfl.symbol >> tfl.position.x() >> tfl.position.y() >> tfl.position.z();
          positions.push_back(tfl.position);
          if ((i-input_ensemble.size()*(N+1u)) == N)
          { // if we are at the end of a structure 
            if(positions.size() != atoms.size()) 
              throw std::logic_error("The size of an additionally provided structure does not match the number of atoms.");
            input_ensemble.push_back(positions);
            positions.clear();
          }
        }
      }
    } // for

    if (indexation_not_contiguous)
    {
      std::cout << "Indexation not contiguous. Rebinding atoms." << lineend;
      for (std::size_t i(0U); i<atoms.size(); ++i)
      {
        for (auto bonding_partner : atoms.atom(i).bonds())
        {
          std::cout << "Partner of atom which is now " << i+1 << " detached from " << bonding_partner << " and rebound to " << index_of_atom[bonding_partner] << lineend;
          atoms.atom(i).detach_from(bonding_partner);
          atoms.atom(i).bind_to(index_of_atom[bonding_partner]);
        }
      }
    }

    if (input_ensemble.empty()) throw std::logic_error("No structures found.");
    coords::PES_Point x(input_ensemble[0u]);
    if (!Config::get().coords.fixed.empty())
    {
      for (auto fix : Config::get().coords.fixed)
      {
        if (fix < atoms.size()) atoms.atom(fix).fix(true);
      }
    }
    coord_object.init_swap_in(atoms, x);
    for (auto & p : input_ensemble)
    {
      p.gradient.cartesian.resize(p.structure.cartesian.size());
      coord_object.set_xyz(p.structure.cartesian);
      coord_object.to_internal();
      p = coord_object.pes();
    }
    
  }
  else throw std::logic_error("Reading the structure input file failed.");
  return coord_object;
}


std::ostream& coords::operator<< (std::ostream &stream, coords::Coordinates const & coord)
{
  if (Config::get().general.output == config::output_types::TINKER) 
    stream << coords::output::formats::tinker(coord);
  else if (Config::get().general.output == config::output_types::XYZ) 
    stream << coords::output::formats::xyz(coord);
  else if (Config::get().general.output == config::output_types::MOLDEN) 
    stream << coords::output::formats::moldenxyz(coord);
  else if (Config::get().general.output == config::output_types::ZMATRIX) 
    stream << coords::output::formats::zmatrix(coord);
  return stream;
}



static void tinker_dummy_to_stream(std::ostream & stream, std::size_t index, coords::float_type x, coords::float_type y, coords::float_type z)
{
  stream << std::right << std::setw(6) << index << "  ";
  stream << std::left << "XX ";
  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << x;
  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << y;
  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << z;
  stream << std::right << std::setw(6) << 0 << "\n";
}

void coords::output::formats::tinker::to_stream (std::ostream & stream) const
{
  bool const pp = Config::get().energy.periodic && Config::get().energy.periodic_print;
  std::size_t const N(ref.size());
  stream << (pp ? N+8 : N) << lineend;
  //std::size_t index_width(1), tens(N);
  //while(tens >= 10U)
  //{
  //  tens /= 10U;
  //  ++index_width;
  //}
  for (std::size_t i(0U); i<N; ++i)
  {
    stream << std::right << std::setw(6) << i+1U << "  ";
    stream << std::left  << std::setw(3) << ref.atoms(i).symbol().substr(0U, 2U);
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z();
    stream << std::right << std::setw(6) << ref.atoms(i).energy_type();
    std::size_t const bSize(ref.atoms(i).bonds().size());
    for (std::size_t j(0U); j<bSize; ++j)
    {
      stream << std::right << std::setw(6) << ref.atoms(i).bonds()[j]+1U;
    }
    if (ref.atoms(i).sub_type() == coords::Atom::sub_types::ST_IN) stream << " IN";
    else if (ref.atoms(i).sub_type() == coords::Atom::sub_types::ST_OUT) stream << " OUT";
    stream << lineend;
  }
  if (pp)
  {
    coords::Cartesian_Point const halfbox = Config::get().energy.pb_box / 2.;
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


void coords::output::formats::moldenxyz::to_stream(std::ostream & stream) const
{
  std::size_t const N(ref.size());
  stream << N << lineend;
  stream << "Energy = " << ref.energyinterface()->energy;
  for (std::size_t i(0U); i<N; ++i)
  {
    stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z();
    stream << lineend;
  }
}

void coords::output::formats::xyz_mopac7::to_stream(std::ostream & stream) const
{
	std::size_t const N(ref.size());
	for (std::size_t i(0U); i < N; ++i)

	{
		if (ref.atoms(i).fixed()) {
			stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
			stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x() << " 0";
			stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y() << " 0";
			stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z() << " 0";
			stream << lineend;
		}
		else
		{
			stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
			stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x() << " 1";
			stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y() << " 1";
			stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z() << " 1";
			stream << lineend;
		}
	}
}

void coords::output::formats::xyz::to_stream (std::ostream & stream) const
{
  std::size_t const N(ref.size());
  //stream << N << lineend;
  for (std::size_t i(0U); i<N; ++i)
  {
   /* stream << std::left  << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z();
    stream << lineend;*/
	  if (ref.atoms(i).fixed()) {
		  stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
		  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x() << " 0";
		  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y() << " 0";
		  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z() << " 0";
		  stream << lineend;
	  }
	  else
	  {
		  stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
		  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x() << " 1";
		  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y() << " 1";
		  stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z() << " 1";
		  stream << lineend;
	  }
  }
}


void coords::output::formats::zmatrix::to_stream (std::ostream & stream) const
{
  std::size_t const N(ref.size());
  if (stream.good() && N > 0)
  {
    stream << "zmat angstroms" << lineend;
    stream << std::left << std::setw(10) << 1U;
    stream << std::left << std::setw(10) << ref.atoms(0u).i_to_a()+1;
    stream << std::left  << std::setw(4)  << ref.atoms(ref.atoms(0u).i_to_a()).symbol() << lineend;
    if (N > 1)
    {
      stream << std::left << std::setw(10) << 2U;
      std::size_t const A = ref.atoms(1u).i_to_a();
      stream << std::left << std::setw(10) << A+1;
      stream << std::left  << std::setw(4)  << ref.atoms(A).symbol();
      stream << std::right  << std::setw(10)  << ref.atoms(1u).ibond()+1;
      stream << std::right  << std::setw(10)  << "bnd" << 1u;
      stream << lineend;
    }
    if (N > 2)
    {
      stream << std::left << std::setw(10) << 3U;
      std::size_t const A = ref.atoms(2u).i_to_a();
      stream << std::left << std::setw(10) << A+1;
      stream << std::left  << std::setw(4)  << ref.atoms(A).symbol();
      stream << std::right  << std::setw(10)  <<  ref.atoms(2u).ibond()+1;
      stream << std::right  << std::setw(10)  << "bnd" << 2u;
      stream << std::right  << std::setw(10)  <<  ref.atoms(2u).iangle()+1;
      stream << std::right  << std::setw(10)  << "ang" << 2u;
      stream << lineend;
    }

    for (std::size_t i=3; i<N; ++i)
    {
      stream << std::left << std::setw(10) << i+1U;
      std::size_t const A = ref.atoms(i).i_to_a();
      stream << std::left << std::setw(10) << A+1;
      stream << std::left  << std::setw(4)  << ref.atoms(A).symbol();
      stream << std::right  << std::setw(10)  <<  ref.atoms(i).ibond()+1;
      stream << std::right  << std::setw(10)  << "bnd" << i;
      stream << std::right  << std::setw(10)  <<  ref.atoms(i).iangle()+1;
      stream << std::right  << std::setw(10)  << "ang" << i;
      stream << std::right  << std::setw(10)  <<  ref.atoms(i).idihedral()+1;
      stream << std::right  << std::setw(10)  << "dih" << i;
      stream << lineend;
    }
    stream << "variables" << lineend;
    if (N > 1)
    {
      stream  << "bnd" << std::left  << std::setw(10) << 1u;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6)  << ref.intern(1u).radius();
      stream << lineend;
    }
    if (N > 2)
    {
      stream  << "bnd" << std::left  << std::setw(10) << 2u;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(2u).radius();
      stream << lineend;
      stream  << "ang" << std::left  << std::setw(10) << 2u;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6)  << ref.intern(2u).inclination();
      stream << lineend;
    }
    for (std::size_t i=3; i<N; ++i)
    {
      stream  << "bnd" << std::left  << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(i).radius();
      stream << lineend;
      stream  << "ang" << std::left  << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(i).inclination();
      stream << lineend;
      stream  << "dih" << std::left  << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6)  << ref.intern(i).azimuth();
      stream << lineend;
    }
    stream << "constants" << lineend;
    for (std::size_t i=1; i<N; ++i)
    {
      stream  << "g_bnd" << std::left  << std::setw(10) << i;
      stream << std::right << std::fixed  << std::setw(20) << std::setprecision(6)  << ref.g_intern(i).x();
      stream << lineend;
      if (i > 1)
      {
        stream  << "g_ang" << std::left  << std::setw(10) << i;
        stream << std::right << std::fixed << std::setw(20) << std::setprecision(6)  << ref.g_intern(i).y();
        stream << lineend;
      }
      if (i > 2)
      {
        stream  << "g_dih" << std::left  << std::setw(10) << i;
        stream << std::right << std::fixed << std::setw(20) << std::setprecision(6)  << ref.g_intern(i).z();
        stream << lineend;
      }
    }
    stream << "end" << lineend;
  } 
  else throw std::runtime_error("ERR_FILE_WRITE: stream bad");
}

void coords::output::formats::xyz_mopac::to_stream(std::ostream &stream) const
{
  std::size_t const N(ref.size());
  //stream << N << lineend;
  for (std::size_t i(0U); i<N; ++i)
  {
    /* stream << std::left  << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y();
    stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z();
    stream << lineend;*/
    if (ref.atoms(i).fixed()) {
      stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x() << " +0";
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y() << " +0";
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z() << " +0";
      stream << lineend;
    }
    else
    {
      stream << std::left << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x() << " +1";
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y() << " +1";
      stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z() << " +1";
      stream << lineend;
    }
  }
}
