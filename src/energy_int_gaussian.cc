#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <cstdlib>
#include <utility>
#include "atomic.h"
#include "energy_int_gaussian.h"
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
Gaussian sysCall functions
*/

energy::interfaces::gaussian::sysCallInterfaceGauss::sysCallInterfaceGauss(coords::Coordinates * cp) :
  energy::interface_base(cp),
  hof_kcal_mol(0.0), hof_kj_mol(0.0), e_total(0.0),
  e_electron(0.0), e_core(0.0), id(Config::get().general.outputFilename), failcounter(0u)
{
  std::stringstream ss;
  ss << (std::size_t(std::rand()) | (std::size_t(std::rand()) << 15));
  id.append("_tmp_").append(ss.str());
  optimizer = true;
}

energy::interfaces::gaussian::sysCallInterfaceGauss::sysCallInterfaceGauss(sysCallInterfaceGauss const & rhs, coords::Coordinates *cobj) :
  interface_base(cobj),
  hof_kcal_mol(rhs.hof_kcal_mol), hof_kj_mol(rhs.hof_kj_mol), e_total(rhs.e_total),
  e_electron(rhs.e_electron), e_core(rhs.e_core), id(rhs.id), failcounter(rhs.failcounter)
{
  interface_base::operator=(rhs);
}

energy::interface_base * energy::interfaces::gaussian::sysCallInterfaceGauss::clone(coords::Coordinates * coord_object) const
{
  sysCallInterfaceGauss * tmp = new sysCallInterfaceGauss(*this, coord_object);
  return tmp;
}

energy::interface_base * energy::interfaces::gaussian::sysCallInterfaceGauss::move(coords::Coordinates * coord_object)
{
  sysCallInterfaceGauss * tmp = new sysCallInterfaceGauss(*this, coord_object);
  return tmp;
}

void energy::interfaces::gaussian::sysCallInterfaceGauss::swap(interface_base &rhs)
{
  swap(dynamic_cast<sysCallInterfaceGauss&>(rhs));
}

void energy::interfaces::gaussian::sysCallInterfaceGauss::swap(sysCallInterfaceGauss &rhs)
{
  interface_base::swap(rhs);
  std::swap(failcounter, rhs.failcounter);
}

energy::interfaces::gaussian::sysCallInterfaceGauss::~sysCallInterfaceGauss(void)
{
  std::string rem_file(id);
  if (Config::get().energy.gaussian.delete_input)
  {
    remove(std::string(id).append(".gjf").c_str());
    remove(std::string(id).append(".log").c_str());
  }
}

void energy::interfaces::gaussian::sysCallInterfaceGauss::print_gaussianInput()
{
  std::string outstring(id);
  outstring.append(".gjf");

  std::ofstream out_file(outstring.c_str(), std::ios_base::out);
  if (out_file)
  {
    if (Config::get().energy.gaussian.link.length() != 0) { //if no link commands are issued the line wil be skipped
      out_file << "%" << Config::get().energy.gaussian.link;
      out_file << '\n';
    }
    out_file << "# " << Config::get().energy.gaussian.command;
    out_file << '\n';
    out_file << '\n';
    out_file << Config::get().general.outputFilename;
    out_file << '\n';
    out_file << '\n';
    out_file << Config::get().energy.gaussian.charge << " ";
    out_file << Config::get().energy.gaussian.multipl;
    out_file << '\n';
    out_file << coords::output::formats::xyz(*coords);
    out_file << '\n';

  }
  else std::runtime_error("Writing Gaussian Inputfile failed.");
}

void::energy::interfaces::gaussian::sysCallInterfaceGauss::read_gaussianOutput(bool const grad, bool const opt)
{
  double const au2kcal_mol(627.5095), eV2kcal_mol(23.061078);  //1 au = 627.5095 kcal/mol
  hof_kcal_mol = hof_kj_mol = energy = e_total = e_electron = e_core = 0.0;
  float hof_au(0.0), e_total_au(0.0);

  auto in_string = id + ".log";

  std::ifstream in_file(in_string.c_str(), std::ios_base::in);
  std::vector <float> occMO, virtMO, excitE;
  std::ofstream mos("MOs.txt", std::ios_base::out), test("test.txt", std::ios_base::out); //ofstream for mo testoutput

  bool done(false),test_lastMOs(false);//to controll if reading was successfull
  coords::Representation_3D g_tmp(coords->size()), xyz_tmp(coords->size());
  std::size_t const atoms = coords->size();

  if (in_file)
  {
    


    std::string buffer;
    while (!in_file.eof())
    {
      std::getline(in_file, buffer);

      if (buffer.find("Alpha  occ. eigenvalues --") != std::string::npos) //ascertain if before mo energieds were read and deleting older data
      {
        if (test_lastMOs == true)
        {
          occMO.erase(occMO.begin(), occMO.end());
          virtMO.erase(virtMO.begin(), virtMO.end());
          test_lastMOs = false;
        }
      }

      if (buffer.find("Alpha  occ. eigenvalues --") != std::string::npos) //reading Mo energies
      {
        for (int i = 0; buffer.length() > (29 + i * 10); i++) //in gaussian output orbital energies are presented in rows of 5
        {
          occMO.push_back(std::stof(buffer.substr(29 + i * 10)));
        }
      }



      if (buffer.find("Alpha virt. eigenvalues --") != std::string::npos)
      {
        for (int i = 0; buffer.length() > (29 + i * 10); i++)
        {
          virtMO.push_back(std::stof(buffer.substr(29 + i * 10)));
        }
        test_lastMOs = true;
      }



      if (buffer.find(" Excited State   ") != std::string::npos)//fetches excitation energies from gaussian output
      {
        excitE.push_back(std::stof(buffer.substr(38)));
      }

      if (buffer.find(" SCF Done:") != std::string::npos)
      {
        e_total_au = std::stof(buffer.substr(buffer.find_first_of("=") + 1));
      }

      if (grad && buffer.find("Old X    -DE/DX   Delta X") != std::string::npos) //fetches last calculated gradients from output
      {
        coords::Cartesian_Point g;


        g_tmp.erase(g_tmp.begin(), g_tmp.end());

        std::getline(in_file, buffer);
        std::getline(in_file, buffer);

        for (std::size_t i(0); i < atoms && !in_file.eof(); ++i)
        {
          std::sscanf(buffer.c_str(), "%*s %*s %*s %*s %*s %lf %*s", &g.x());
          std::getline(in_file, buffer);
          std::sscanf(buffer.c_str(), "%*s %*s %*s %*s %*s %lf %*s", &g.y());
          std::getline(in_file, buffer);
          std::sscanf(buffer.c_str(), "%*s %*s %*s %*s %*s %lf %*s", &g.z());

          std::getline(in_file, buffer);

          g_tmp[i] = g;

        }
      }//end gradient reading

      if (grad && buffer.find("Center     Atomic      Atomic             Coordinates (Angstroms)") != std::string::npos)//reads last coordinates from file
      {
        coords::Cartesian_Point p;

        xyz_tmp.erase(xyz_tmp.begin(), xyz_tmp.end());

        std::getline(in_file, buffer);
        std::getline(in_file, buffer);


        for (std::size_t i(0); i < atoms && !in_file.eof(); ++i)
        {
          std::getline(in_file, buffer);

          test << buffer << '\n';

          std::sscanf(buffer.c_str(), "%*s %*s %*s %lf %lf %lf", &p.x(), &p.y(), &p.z());


          if (opt)
          {
            xyz_tmp[i] = p;
          }
        }
      }//end coordinater reading
    }//end while(!in_file.eof())

    } //end if(in_file)

     std::sort(occMO.begin(), occMO.end(), std::greater <float>()); //sort occupied mos highest to lowest
     std::sort(virtMO.begin(), virtMO.end()); //sort unoccupied mos lowest to highest
     std::sort(excitE.begin(), excitE.end()); //sort excitation energies lowest to highest

     //converting energie units to kcal/mol

     for (int i = 0; i < occMO.size(); i++)
     {
       occMO[i] *= au2kcal_mol;
     }

     for (int i = 0; i < virtMO.size(); i++)
     {
       virtMO[i] *= au2kcal_mol;
     }

     for (int i = 0; i < excitE.size(); i++)
     {
       excitE[i] *= eV2kcal_mol;
     }

    e_total = e_total_au * au2kcal_mol;

    for (int i=0; i < occMO.size(); i++)
    {
      mos << occMO[i] << "    " << virtMO[i] << '\n';
    }

    mos << '\n';

    for (int i=0; i < xyz_tmp.size(); i++)
    {
      mos << xyz_tmp[i].x() << "  " << xyz_tmp[i].y() << "  " << xyz_tmp[i].z() << "  " << " | " << g_tmp[i].x() << " " << g_tmp[i].y() << " " << g_tmp[i].z() << " " << '\n';
    }

    mos << '\n';

    for (int i = 0; i < excitE.size(); i++)
    {
      mos << excitE[i] << '\n';
    }

    mos << '\n' << e_total << " | " << e_total_au << '\n';

  }


int energy::interfaces::gaussian::sysCallInterfaceGauss::callGaussian()
{
  auto gaussian_call = Config::get().energy.gaussian.path + " " + id + ".gjf";
//#ifdef _MSC_VER
//  gaussian_call.append(" > nul 2>&1");
//#else
//  gaussian_call.append(" /dev/null 2>&1");
//#endif

  auto ret = scon::system_call(gaussian_call);
  if (ret != 0)
  {
    ++failcounter;
    std::cout << "Warning Gaussian call to '" << gaussian_call
      << " ' did not return 0 and probably failes.\n( A total of "
      << failcounter << " Gaussian call" << (failcounter != 1 ? "s have" : " has")
      << " failed so far.)\n";
    if (failcounter > 100u)
    {
      throw std::runtime_error("More than 100 Gaussian calls have failed.");
    }
  }
  return ret;
}

//Energy functions
double energy::interfaces::gaussian::sysCallInterfaceGauss::e(void)
{
  integrity = true;
  print_gaussianInput();
  if (callGaussian() == 0) read_gaussianOutput();
  else
  {
    if (Config::get().general.verbosity >=2)
    {
      std::cout << "Gaussian call return value was not 0. Treating structure as broken.\n";
    }
    integrity = false;
  }
  return energy;
}

double energy::interfaces::gaussian::sysCallInterfaceGauss::g(void)
{
  if (Config::get().general.verbosity >= 2)
  {
    std::cout << "Gradient not implemented in CAST as yet.";
  }
  return energy;
}

double energy::interfaces::gaussian::sysCallInterfaceGauss::h(void)
{
  if (Config::get().general.verbosity >= 2)
  {
    std::cout << "Hessian not implemented in CAST as yet.";
  }
  integrity = true;
  print_gaussianInput();
  if (callGaussian() == 0) read_gaussianOutput();
  else
  {
    if (Config::get().general.verbosity >= 2)
    {
      std::cout << "Gaussian call return value was not 0. Treating structure as broken.\n";
    }
    integrity = false;
  }
  return energy;
}

double energy::interfaces::gaussian::sysCallInterfaceGauss::o(void)
{
  if (Config::get().general.verbosity >= 2)
  {
    std::cout << "Optimization not implemented in CAST as yet.";
  }
  integrity = true;
  print_gaussianInput();
  if (callGaussian() == 0) read_gaussianOutput();
  else
  {
    if (Config::get().general.verbosity >= 2)
    {
      std::cout << "Gaussian call return value was not 0. Treating structure as broken.\n";
    }
    integrity = false;
  }
  return energy;
}

//Outputfunctions

void energy::interfaces::gaussian::sysCallInterfaceGauss::print_E(std::ostream &S) const
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

void energy::interfaces::gaussian::sysCallInterfaceGauss::print_E_head(std::ostream &S, bool const endline) const
{
  S << std::right << std::setw(16) << "HOF";
  S << std::right << std::setw(16) << "TOT";
  S << std::right << std::setw(16) << "EL";
  S << std::right << std::setw(16) << "CORE";
  if (endline) S << '\n';
}

void energy::interfaces::gaussian::sysCallInterfaceGauss::print_E_short(std::ostream &S, bool const endline) const
{
  S << std::right << std::setw(16) << std::scientific << std::setprecision(5) << hof_kcal_mol;
  S << std::right << std::setw(16) << std::scientific << std::setprecision(5) << e_total;
  S << std::right << std::setw(16) << std::scientific << std::setprecision(5) << e_electron;
  S << std::right << std::setw(16) << std::scientific << std::setprecision(5) << e_core;
  if (endline) S << '\n';
}

void energy::interfaces::gaussian::sysCallInterfaceGauss::print_G_tinkerlike(std::ostream &, bool const) const { }

void energy::interfaces::gaussian::sysCallInterfaceGauss::to_stream(std::ostream&) const { }

bool energy::interfaces::gaussian::sysCallInterfaceGauss::check_bond_preservation(void) const
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