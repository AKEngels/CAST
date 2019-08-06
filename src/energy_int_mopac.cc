#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <cstdlib>
#include <utility>
#include <cstdio>

#include "atomic.h"
#include "energy_int_mopac.h"
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
Mopac sysCall functions
*/

energy::interfaces::mopac::sysCallInterface::sysCallInterface(coords::Coordinates * cp) :
  energy::interface_base(cp),
  hof_kcal_mol(0.0), hof_kj_mol(0.0), e_total(0.0),
  e_electron(0.0), e_core(0.0), failcounter(0u)
{
	id = Config::get().general.outputFilename;
  std::stringstream ss;
  std::srand(std::time(0));
  ss << (std::size_t(std::rand()) | (std::size_t(std::rand()) << 15));
  id.append("_tmp_").append(ss.str());
  optimizer = true;
	charge = Config::get().energy.mopac.charge;
}

energy::interfaces::mopac::sysCallInterface::sysCallInterface(sysCallInterface const & rhs, coords::Coordinates *cobj) :
  interface_base(cobj),
  hof_kcal_mol(rhs.hof_kcal_mol), hof_kj_mol(rhs.hof_kj_mol), e_total(rhs.e_total),
  e_electron(rhs.e_electron), e_core(rhs.e_core), failcounter(rhs.failcounter)
{
	id = rhs.id;
	charge = rhs.charge;
  interface_base::operator=(rhs);
}

energy::interface_base * energy::interfaces::mopac::sysCallInterface::clone(coords::Coordinates * coord_object) const
{
  sysCallInterface * tmp = new sysCallInterface(*this, coord_object);
  return tmp;
}

energy::interface_base * energy::interfaces::mopac::sysCallInterface::move(coords::Coordinates * coord_object)
{
  sysCallInterface * tmp = new sysCallInterface(*this, coord_object);
  return tmp;
}

void energy::interfaces::mopac::sysCallInterface::swap(interface_base &rhs)
{
  swap(dynamic_cast<sysCallInterface&>(rhs));
}

void energy::interfaces::mopac::sysCallInterface::swap(sysCallInterface &rhs)
{
  interface_base::swap(rhs);
  std::swap(failcounter, rhs.failcounter);
}

energy::interfaces::mopac::sysCallInterface::~sysCallInterface(void)
{
  std::string rem_file(id);
  if (Config::get().energy.mopac.delete_input)
  {
    std::remove(std::string(id).append(".xyz").c_str());
    std::remove(std::string(id).append(".out").c_str());
    std::remove(std::string(id).append(".arc").c_str());
    std::remove(std::string(id).append("_sys.out").c_str());
    std::remove(std::string(id).append(".xyz.out").c_str());
    if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB)
    {
      std::remove("FOR005");
      std::remove("FOR006");
      std::remove("FOR012");
    }

  }
}

void energy::interfaces::mopac::sysCallInterface::read_charges() 
{
  auto file = id + ".xyz.aux";
  std::ifstream auxstream{ file };
  if (!auxstream)
  {
    throw std::runtime_error("AUX file not found.");
  }
	atomic_charges.clear();
  std::string auxline;
  while (std::getline(auxstream, auxline))
  {
    if (auxline.find("ATOM_CHARGES") != std::string::npos)
    {
      break;
    }
  }
  double charge{};
  while (auxstream >> charge)
  {
		atomic_charges.push_back(charge);
  }
  if (atomic_charges.size() != coords->size())
  {
    throw std::logic_error("Found " + std::to_string(atomic_charges.size()) +
      " charges instead of " + std::to_string(coords->size()) + " charges.");
  }
}

// write mol.in (see http://openmopac.net/manual/QMMM.html)
void energy::interfaces::mopac::sysCallInterface::write_mol_in()
{
	auto elec_factor = 332.0;
	std::cout << std::setprecision(6);

	std::ofstream molstream{ "mol.in" };
	if (molstream)
	{
		auto const n_qm = coords->size() - Config::get().energy.mopac.link_atoms; // number of QM atoms
		molstream << '\n';
		molstream << n_qm << " "<<Config::get().energy.mopac.link_atoms <<"\n";  // number of link atoms
		for (std::size_t i = 0; i < coords->size(); ++i)  // for every atom (QM + link atom)
		{
			double qi{};
			for (std::size_t j = 0; j < Config::get().energy.qmmm.mm_charges.size(); ++j) // for every MM atom
			{
				auto current_charge = Config::get().energy.qmmm.mm_charges[j];
				double dist_x = coords->xyz(i).x() - current_charge.x;
				double dist_y = coords->xyz(i).y() - current_charge.y;
				double dist_z = coords->xyz(i).z() - current_charge.z;
				double d = std::sqrt(dist_x*dist_x + dist_y * dist_y + dist_z * dist_z); // distance between QM and MM atom

				qi += current_charge.scaled_charge / d;
			}
			qi *= elec_factor;
			molstream << "0 0 0 0 " << qi << "\n";
		}
		molstream.close();
	}
	else
	{
		throw std::runtime_error("Cannot write mol.in file.");
	}
}

void energy::interfaces::mopac::sysCallInterface::print_mopacInput(bool const grad, bool const hess, bool const opt)
{
	if (Config::get().energy.qmmm.mm_charges.size() != 0) write_mol_in();

  std::string outstring(id);
  outstring.append(".xyz");

  /*! Special Input for modified MOPAC7
  *
  * MOPAC7_HB is a modified version of MOPAC7 which
  * neglects hydrogenbonding and needs a special input.
  */
  if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB)
  {
    outstring = "FOR005";
  }

  std::ofstream out_file(outstring.c_str(), std::ios_base::out);
  if (out_file)
  {
    if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7)
    {
      out_file << Config::get().energy.mopac.command << " AUX ";
			if (Config::get().energy.qmmm.mm_charges.size() != 0) out_file << " QMMM ";
      out_file << (opt ?  " LINMIN" : " 1SCF");
      out_file << (grad ? " GRADIENTS" : "") << (hess ? " HESSIAN" : "");
			if (charge != 0) out_file << " CHARGE=" << charge;
    }
    else if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB)
    {
      std::size_t found(0U);
      found = Config::get().energy.mopac.command.find("+");
      if (found > 0 && found < 100)
      {
        std::string str_t1, str_t2;
        str_t1 = Config::get().energy.mopac.command.substr(0, found + 1);
        str_t2 = Config::get().energy.mopac.command.substr(found + 1);
				if (Config::get().energy.qmmm.mm_charges.size() != 0) out_file << " QMMM ";
        out_file << (opt ? (coords->size() > 250 ? "EF " : "EF ") : "1SCF ") << (grad ? "GRADIENTS " : " ");
        out_file << (hess ? "HESSIAN " : " ") << str_t1 << '\n';
        out_file << str_t2;
				if (charge != 0) out_file << " CHARGE=" << charge;
      }
      else
      {
        out_file << Config::get().energy.mopac.command << " AUX ";
				if (Config::get().energy.qmmm.mm_charges.size() != 0) out_file << " QMMM ";
        out_file << (opt ? (coords->size() > 250 ? " EF" : " EF") : " 1SCF");
        out_file << (grad ? " GRADIENTS" : "") << (hess ? " HESSIAN" : "");
				if (charge != 0) out_file << " CHARGE=" << charge;
      }
    }
    else                                                                       // "normal" MOPAC, e.g. version 2016
    {
      out_file << Config::get().energy.mopac.command << " AUX ";
			if (Config::get().energy.qmmm.mm_charges.size() != 0) out_file << " QMMM ";
      out_file << (opt ? (coords->size() > 250 ? " LBFGS" : " EF") : " 1SCF");
      out_file << (grad ? " GRADIENTS" : "") << (hess ? " HESSIAN" : "");
			if (charge != 0) out_file << " CHARGE=" << charge;
    }
    if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC2012MT)
    {
#if defined _OPENMP
      out_file << " THREADS=" << omp_get_num_threads();
#else
      out_file << " THREADS=1";
#endif
    }

    out_file << '\n';
    out_file << '\n' << '\n';

    out_file << coords::output::formats::xyz_mopac(*coords);

		if (Config::get().periodics.periodic)  // translation vectors for unit cell
		{
			out_file << std::left << std::setw(3) << "Tv";
			out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << Config::get().periodics.pb_box.x() << " 1";
			out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0 << " 1";
			out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0 << " 1";
			out_file << '\n';

			out_file << std::left << std::setw(3) << "Tv";
			out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0 << " 1";
			out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << Config::get().periodics.pb_box.y() << " 1";
			out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0 << " 1";
			out_file << '\n';

			out_file << std::left << std::setw(3) << "Tv";
			out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0 << " 1";
			out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << 0.0 << " 1";
			out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << Config::get().periodics.pb_box.z() << " 1";
			out_file << '\n';
		}
  }
  else throw std::runtime_error("Writing MOPAC Inputfile failed.");
}

void energy::interfaces::mopac::sysCallInterface::read_mopacOutput(bool const grad, bool const, bool const opt)
{
  hof_kcal_mol = hof_kj_mol = energy = e_total = e_electron = e_core = 0.0;
  std::string in_string = id + ".out";
  if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB) { in_string = "FOR006"; remove("FOR012"); }
  std::ifstream in_file(in_string.c_str(), std::ios_base::in);
  //std::size_t fixcounter(0);
  // ...

  if (!in_file)
  {
    std::string const alt_infile(std::string(id).append(".xyz.out"));
    remove(std::string(id).append(".xyz.arc").c_str());
    if (Config::get().general.verbosity > 4)
    {
      std::cout << "Input file '" << in_string << "' not found, trying '" << alt_infile << "'.\n";
    }
    in_file.open(alt_infile.c_str(), std::ios_base::in);

	if (in_file)
		in_string = alt_infile;
	else
		throw std::runtime_error(std::string("MOPAC OUTPUT NOT PRESENT; ID: ").append(id));


  }
  // ...
  double constexpr ev_to_kcal(23.060547411538901355);
  bool const mozyme(Config::get().energy.mopac.command.find("MOZYME") != std::string::npos);
  bool done(false);
  coords::Representation_3D g_tmp(coords->size()), xyz_tmp(coords->size());
  if (in_file)
  {
		if (Config::get().energy.qmmm.use) read_charges();    // read atomic charges

    std::string buffer;
    while (!in_file.eof())
    {
      std::getline(in_file, buffer);

      if (buffer.find("FINAL HEAT OF FORMATION") != std::string::npos)
      {
        done = true;

        std::sscanf(buffer.c_str(), "%*s %*s %*s %*s %*s %lf %*s %*s %lf %*s", &hof_kcal_mol, &hof_kj_mol);

      }

      else if (buffer.find("TOTAL ENERGY") != std::string::npos)
        std::sscanf(buffer.c_str(), "%*s %*s %*s %lf %*s", &e_total);

      else if (buffer.find("ELECTRONIC ENERGY") != std::string::npos)
        std::sscanf(buffer.c_str(), "%*s %*s %*s %lf %*s", &e_electron);

      else if (buffer.find("CORE-CORE REPULSION") != std::string::npos)
        std::sscanf(buffer.c_str(), "%*s %*s %*s %lf %*s", &e_core);

      else if (grad && buffer.find("FINAL  POINT  AND  DERIVATIVES") != std::string::npos)
      {
        done = true;
        std::getline(in_file, buffer); // empty line
        std::getline(in_file, buffer); // Header
        if (mozyme)
        {
          std::getline(in_file, buffer); // Headers 2
          std::getline(in_file, buffer); // empty line
        }
        std::size_t const atoms = coords->size();
        std::size_t i(0), j;
        coords::Cartesian_Point p, g;
        for (; i < atoms && !in_file.eof(); ++i)
        {
          // 				  if(!integrity){i=j;}
          // 				        integrity=true;
          if (coords->atoms(i).fixed()) { p.x() = coords->xyz(i).x(); p.y() = coords->xyz(i).y(); p.z() = coords->xyz(i).z(); continue; }
          if (mozyme)
          {
            std::getline(in_file, buffer);
            std::istringstream bss(buffer);
            std::size_t tx;
            std::string sx;
            bss >> tx >> sx >> g.x() >> g.y() >> g.z() >> p.x() >> p.y() >> p.z();
          }
          else if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7)
          {
            std::getline(in_file, buffer);
            std::istringstream bss_x(buffer);
            std::size_t tx;
            std::string sx[3];
            bss_x >> tx >> j >> sx[0] >> sx[1] >> p.x() >> g.x();

            if (i != (j - 1))
            {
              if (Config::get().general.verbosity >= 2) std::cout << "Index " << i << " does not match line index " << j - 1 << '\n';
              integrity = false;
              return;
            }
            std::getline(in_file, buffer);
            std::istringstream bss_y(buffer);
            bss_y >> tx >> j >> sx[0] >> sx[1] >> p.y() >> g.y();
            if (i != (j - 1))
            {
              if (Config::get().general.verbosity >= 2) std::cout << "MOPAC: Index " << i << " does not match line index " << j - 1 << '\n';
              integrity = false;
              return;
            }
            std::getline(in_file, buffer);
            std::istringstream bss_z(buffer);
            bss_z >> tx >> j >> sx[0] >> sx[1] >> p.z() >> g.z();
            if (i != (j - 1))
            {
              if (Config::get().general.verbosity >= 2) std::cout << "MOPAC: Index " << i << " does not match line index " << j - 1 << '\n';
              integrity = false;
              return;
            }


          }
          else
          {
            if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB && coords->xyz(i).x() == 0.0) {}
            else {
              std::getline(in_file, buffer);
              std::istringstream bss_x(buffer);
              std::size_t tx;
              std::string sx[3];
              bss_x >> tx >> j >> sx[0] >> sx[1] >> sx[2] >> p.x() >> g.x();
              // 						std::cout << g.x() << '\n';
              /*if (i != (j - 1) && Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB)
              {

              integrity = false;
              if(sx[2] == "X") {g_tmp[j-1].x()=g.x();}
              }else*/ if (i != (j - 1) && Config::get().energy.mopac.version != config::mopac_ver_type::MOPAC7_HB)
              {
                if (Config::get().general.verbosity >= 2) std::cout << "MOPAC: Index " << i << " does not match line index " << j - 1 << '\n';
                integrity = false;
                return;
              }
            }
            if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB && coords->xyz(i).y() == 0.0) {}
            else {
              std::getline(in_file, buffer);
              std::istringstream bss_y(buffer);
              std::size_t tx;
              std::string sx[3];
              bss_y >> tx >> j >> sx[0] >> sx[1] >> sx[2] >> p.y() >> g.y();
              /*if (i != (j - 1) && Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB)
              {

              integrity = false;
              if(sx[2] == "Y") {g_tmp[j-1].y()=g.y();}
              }else*/ if (i != (j - 1) && Config::get().energy.mopac.version != config::mopac_ver_type::MOPAC7_HB)
              {
                if (Config::get().general.verbosity >= 2) std::cout << "MOPAC: Index " << i << " does not match line index " << j - 1 << '\n';
                integrity = false;
                return;
              }
            }if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB && coords->xyz(i).z() == 0.0) {}
            else {
              std::getline(in_file, buffer);
              std::istringstream bss_z(buffer);
              std::size_t tx;
              std::string sx[3];
              bss_z >> tx >> j >> sx[0] >> sx[1] >> sx[2] >> p.z() >> g.z();
              /*if (i != (j - 1) && Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB)
              {

              integrity = false;
              if(sx[2] == "Z") {g_tmp[j-1].z()=g.z();}
              }else*/ if (i != (j - 1) && Config::get().energy.mopac.version != config::mopac_ver_type::MOPAC7_HB)
              {
                if (Config::get().general.verbosity >= 2) std::cout << "MOPAC: Index " << i << " does not match line index " << j - 1 << '\n';
                integrity = false;
                return;
              }
            }
          }

          // update gradients
          g_tmp[i] = g;

          //coords->update_g_xyz(i, g);
          // update optimized structure if needed
          if (opt)
          {
            xyz_tmp[i] = p;
            //coords->move_atom_to(i, p);
          }
        } // for atoms

				if (Config::get().energy.qmmm.mm_charges.size() != 0)  // if QM/MM: add coulomb gradient due to external charges
				{
					double constexpr elec_factor = 332.06;
					grad_ext_charges.clear();
					grad_ext_charges.resize(Config::get().energy.qmmm.mm_charges.size());

					for (auto i = 0u; i < coords->size(); ++i) // for every QM atom
					{
						double charge_i = charges()[i];

						for (auto j = 0u; j < Config::get().energy.qmmm.mm_charges.size(); ++j) // for every external charge
						{
							auto current_charge = Config::get().energy.qmmm.mm_charges[j];

							coords::r3 r_ij;
							r_ij.x() = current_charge.x - coords->xyz(i).x();
							r_ij.y() = current_charge.y - coords->xyz(i).y();
							r_ij.z() = current_charge.z - coords->xyz(i).z();
							coords::float_type d = len(r_ij);

							coords::float_type b = (charge_i*current_charge.scaled_charge) / d * elec_factor;
							coords::float_type db = b / d;
							auto c_gradient_ij = r_ij * db / d;
							g_tmp[i] += c_gradient_ij;            // gradient on QM atom
							grad_ext_charges[j] -= c_gradient_ij; // gradient on external charge (= MM atom)
						}

					}

				}
      }
      else if (buffer.find("CARTESIAN COORDINATES") != std::string::npos && Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7_HB)
      {
        for (std::size_t n(0); n < 3; n++)
        {
          std::getline(in_file, buffer);
        }
        for (std::size_t i = 0; i < coords->size(); ++i)
        {
          std::size_t tx;
          coords::Cartesian_Point l;
          std::string sx[4];
          std::getline(in_file, buffer);
          std::istringstream bss_z(buffer);

          bss_z >> tx >> sx[0] >> l.x() >> l.y() >> l.z();
          if (opt)
          {
            xyz_tmp[i] = l;
            //coords->move_atom_to(i, p);
          }
        }
      }



      // if grad && final point

      //Cartesian Coordinates for MOPAC7

      if (Config::get().energy.mopac.version == config::mopac_ver_type::MOPAC7)
      {

        for (std::size_t n(0); n < 9; n++)
        {
          std::getline(in_file, buffer);
        }

        /*  std::cout << "BUFFER " << buffer << '\n';*/
        for (std::size_t i = 0; i < coords->size(); ++i)
        {
          std::size_t tx;
          std::string sx[4];
          std::getline(in_file, buffer);
          coords::Cartesian_Point p;
          std::istringstream bss_z(buffer);
          if (!coords->atoms(i).fixed()) bss_z >> tx >> sx[0] >> p.x() >> sx[1] >> p.y() >> sx[2] >> p.z() >> sx[3];
          if (coords->atoms(i).fixed()) bss_z >> tx >> sx[0] >> p.x() >> p.y() >> p.z();
          if (opt)
          {
            xyz_tmp[i] = p;
            //coords->move_atom_to(i, p);
          }
        }

      }


    } // while !eof
    energy = hof_kcal_mol;
    e_total *= ev_to_kcal;
    e_electron *= ev_to_kcal;
    e_core *= ev_to_kcal;


    if (!done)
    {
#ifndef _CAST_ALL_WARNINGS_ARE_ERRORS
      if (Config::get().general.verbosity >= 1)
		  std::cout << "ATTENTION: MOPAC calculation was not finished correctly.\n";
	  // Keep (=Copy) failed output file if verbosity is >=3
	  if (Config::get().general.verbosity >= 3)
	  {
		  std::ifstream in_file_failed(in_string.c_str(), std::ios_base::in);
		  std::ofstream keep_file_failed(std::string(in_string + "_FAILED").c_str(), std::ios_base::out);
		  if (in_file_failed && keep_file_failed)
		  {
			  while (!in_file_failed.eof())
			  {
				  std::string str_buffer;
				  std::getline(in_file_failed, str_buffer);
				  keep_file_failed << str_buffer;
			  }
		  }
		  else
		  {
			  std::cout << "Could not copy failed MOPAC output file. Run with verbosity < 3 (this only masks the error, but it still exists) or fix the issue.\n";
			  throw std::runtime_error("Could not copy failed MOPAC output file");
		  }
	  }
#else
	  throw std::runtime_error("ATTENTION: MOPAC calculation was not finished correctly.\n");
#endif
      energy = e_total = e_electron = e_core = 0.0;
      integrity = false;
      return;
    }
    else
    {
      if (opt)
      {
        coords->set_xyz(std::move(xyz_tmp));
      }
      coords->swap_g_xyz(g_tmp);
    }


  } // if open
  else
  {
    throw std::runtime_error(std::string("MOPAC OUTPUT NOT PRESENT; ID: ").append(id));
  }
  auto bpv = true;
	if (coords->check_bond_preservation() == false) bpv = false;
	else if (coords->check_for_crashes() == false) bpv = false;
  if (Config::get().general.verbosity > 3 && !bpv)
  {
    std::cout << "Broken bond detected, structure integrity not warranted anymore. "
      " Current conformation will be treated as broken structure.\n";
  }
  integrity = bpv;
}

int energy::interfaces::mopac::sysCallInterface::callMopac()
{
  auto mopac_call = Config::get().energy.mopac.path + " " + id + ".xyz";
  mopac_call.append(" > output_mopac.txt 2>&1");
  auto ret = scon::system_call(mopac_call);
  if (ret != 0)
  {
    ++failcounter;
    std::cout << "Warning: MOPAC call to '" << mopac_call <<
      "' did not return 0 and probably failed.\n(A total of "
      << failcounter << " MOPAC call" << (failcounter != 1 ? "s have" : " has")
      << " failed so far.)\n";
    if (failcounter > 100u)
    {
      throw std::runtime_error("More than 100 MOPAC calls have failed.");
    }
  }
  return ret;
}

/*
Energy class functions that need to be overloaded
*/

// Energy function
double energy::interfaces::mopac::sysCallInterface::e(void)
{
  integrity = coords->check_structure();
  grad_var = false;
	if (integrity == true)
	{
		print_mopacInput(false, false, false);
		if (callMopac() == 0) read_mopacOutput(false, false, false);
		else
		{
			if (Config::get().general.verbosity >= 2)
			{
				std::cout << "MOPAC call return value was not 0. Treating structure as broken.\n";
			}
			integrity = false;
		}
		return energy;
	}
	else return 0;  // energy = 0 if structure contains NaN
}

// Energy+Gradient function
double energy::interfaces::mopac::sysCallInterface::g(void)
{
	integrity = coords->check_structure();
	grad_var = true;
	if (integrity == true)
	{
	  print_mopacInput(true, false, false);
	  if (callMopac() == 0) read_mopacOutput(true, false, false);
	  else
	  {
		  ++failcounter;
		  std::cout << "MOPAC call failed. " << failcounter << " MOPAC calls have failed so far.\n";
		  if (failcounter > 1000u)
		  {
			  std::cout << "More than 1000 MOPAC calls have failed. Aborting." << std::endl;
			  throw std::runtime_error("MOPAC call failed.");
		  }
	  }
	  return energy;
  }
  else return 0;  // energy = 0 if structure contains NaN
}

// Energy+Gradient+Hessian function
double energy::interfaces::mopac::sysCallInterface::h(void)
{
  integrity = coords->check_structure();
  grad_var = false;
	if (integrity == true)
	{
		print_mopacInput(true, true, false);
		if (callMopac() == 0) read_mopacOutput(true, true, false);
		else
		{
			++failcounter;
			std::cout << "MOPAC call failed. A total of " << failcounter << " MOPAC calls have failed so far.\n";
			if (failcounter > 1000u)
			{
				std::cout << "More than 1000 MOPAC calls have failed. Aborting." << std::endl;
				throw std::runtime_error("MOPAC call failed.");
			}
		}
		return energy;
	}
	else return 0;  // energy = 0 if structure contains NaN
}

// Optimization
double energy::interfaces::mopac::sysCallInterface::o(void)
{
  integrity = coords->check_structure();
  grad_var = false;
	if (integrity == true)
	{
		print_mopacInput(true, false, true);
		if (callMopac() == 0) read_mopacOutput(true, false, true);
		else
		{
			++failcounter;
			std::cout << "MOPAC call failed. A total of " << failcounter << " MOPAC calls have failed so far.\n";
			if (failcounter > 1000u)
			{
				std::cout << "More than 1000 MOPAC calls have failed. Aborting." << std::endl;
				throw std::runtime_error("MOPAC call failed.");
			}
		}
		return energy;
	}
	else return 0;  // energy = 0 if structure contains NaN
}

// Output functions
void energy::interfaces::mopac::sysCallInterface::print_E(std::ostream &S) const
{
  S << "Heat of Formation: ";
  S << std::right << std::setw(16) << std::fixed << std::setprecision(8) << hof_kcal_mol;
  S << "Total Energy:      ";
  S << std::right << std::setw(16) << std::fixed << std::setprecision(8) << e_total;
  S << "Electronic Energy: ";
  S << std::right << std::setw(16) << std::fixed << std::setprecision(8) << e_electron;
  S << "Core-Core Energy:  ";
  S << std::right << std::setw(16) << std::fixed << std::setprecision(8) << e_core << '\n';
}

void energy::interfaces::mopac::sysCallInterface::print_E_head(std::ostream &S, bool const endline) const
{
  S << std::right << std::setw(16) << "HOF";
  S << std::right << std::setw(16) << "TOT";
  S << std::right << std::setw(16) << "EL";
  S << std::right << std::setw(16) << "CORE";
  if (endline) S << '\n';
}

void energy::interfaces::mopac::sysCallInterface::print_E_short(std::ostream &S, bool const endline) const
{
  S << std::right << std::setw(16) << std::scientific << std::setprecision(5) << hof_kcal_mol;
  S << std::right << std::setw(16) << std::scientific << std::setprecision(5) << e_total;
  S << std::right << std::setw(16) << std::scientific << std::setprecision(5) << e_electron;
  S << std::right << std::setw(16) << std::scientific << std::setprecision(5) << e_core;
  if (endline) S << '\n';
}

void energy::interfaces::mopac::sysCallInterface::to_stream(std::ostream&) const { }
