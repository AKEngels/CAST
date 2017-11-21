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
#include <iterator>
#if defined (_MSC_VER)
#include "win_inc.h"
#endif

#ifdef _MSC_VER
#pragma warning (disable: 4996)
#endif

/*
Gaussian sysCall functions
*/

/**
* Global static instance of the sysCallInterfaceGauss-object.
*/

energy::interfaces::gaussian::sysCallInterfaceGauss::sysCallInterfaceGauss(coords::Coordinates * cp) :
  energy::interface_base(cp),
  hof_kcal_mol(0.0), hof_kj_mol(0.0), e_total(0.0),
  e_electron(0.0), e_core(0.0), id(Config::get().general.outputFilename), failcounter(0u)
{
  std::stringstream ss;
  std::srand(std::time(0));
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

void energy::interfaces::gaussian::sysCallInterfaceGauss::print_gaussianInput(char calc_type)
{
  std::string outstring(id);
  outstring.append(".gjf");

  std::ofstream out_file(outstring.c_str(), std::ios_base::out);

  if (out_file)
  {
    if (Config::get().energy.gaussian.link.length() != 0) { //if no link commands are issued the line wil be skipped

      std::istringstream iss(Config::get().energy.gaussian.link);

      std::vector<std::string> splitted_str(
        std::istream_iterator<std::string>{iss},
        std::istream_iterator<std::string>{}
      );

      for (std::size_t i = 0; i < splitted_str.size(); i++)
      {
        out_file << '%' << splitted_str[i] << '\n';
      }
      
    }
    out_file << "# " << Config::get().energy.gaussian.method << " " << Config::get().energy.gaussian.basisset << " " << Config::get().energy.gaussian.spec << " ";

    switch (calc_type) {// to ensure the needed gaussian keywords are used in gausian inputfile for the specified calculation
      case 'o' :
        if (Config::get().energy.gaussian.steep)
        {
          out_file << " Opt=(Steep) ";//Cartesian,
        }
        else
        {
          out_file << " Opt"; //=Cartesian 
        }
        break;
      case 'g' :
        out_file << " Force";
        break;
    }

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
    out_file.close();
  }
  else std::runtime_error("Writing Gaussian Inputfile failed.");
}

void energy::interfaces::gaussian::sysCallInterfaceGauss::read_gaussianOutput(bool const grad, bool const opt)
{
  //std::ofstream mos("MOs.txt", std::ios_base::out); //ofstream for mo testoutput keep commented if not needed

  double const au2kcal_mol(627.5095), eV2kcal_mol(23.061078);  //1 au = 627.5095 kcal/mol
  double const HartreePerBohr2KcalperMolperAngstr = 627.5095 * (1 / 0.52918);
  hof_kcal_mol = hof_kj_mol = energy = e_total = e_electron = e_core = 0.0;

  auto in_string = id + ".log";

  std::ifstream in_file(in_string.c_str(), std::ios_base::in);
  

  bool test_lastMOs(false);//to controll if reading was successfull
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

      if (buffer.find("Excited to excited state transition electric dipole moments (Au):") != std::string::npos)
      {

        std::getline(in_file, buffer);

        bool el_dipm = true;
        int tmp_i, tmp_j;
        coords::Cartesian_Point tmp_ex_ex_trans;

        std::getline(in_file, buffer);
        while (el_dipm)
        {

          std::sscanf(buffer.c_str(), "%i %i %lf %lf %lf %*s %*s", &tmp_i, &tmp_j, &tmp_ex_ex_trans.x(), &tmp_ex_ex_trans.y(), &tmp_ex_ex_trans.z());
          state_i.push_back(tmp_i);
          state_j.push_back(tmp_j);
          ex_ex_trans.push_back(tmp_ex_ex_trans);

          std::getline(in_file, buffer);

          if (buffer.find("Excited to excited state transition velocity") != std::string::npos) {el_dipm = false;}
        }

      }

      if (buffer.find("Ground to excited state transition electric") != std::string::npos)
      {
        std::getline(in_file, buffer);
        bool gz_az_dipm = true;
        int tmp_gz_i;
        coords::Cartesian_Point tmp_gz_ex_trans;

        std::getline(in_file, buffer);
        while (gz_az_dipm)
        { 
          std::sscanf(buffer.c_str(), "%i %lf %lf %lf %*s %*s", &tmp_gz_i, &tmp_gz_ex_trans.x(), &tmp_gz_ex_trans.y(), &tmp_gz_ex_trans.z());
          gz_i_state.push_back(tmp_gz_i);
          gz_ex_trans.push_back(tmp_gz_ex_trans);

          std::getline(in_file, buffer);

          if (buffer.find("Ground to excited state transition velocity dipole moments (Au):") != std::string::npos) { gz_az_dipm = false; }
        }

      }



      if (buffer.find(" Excited State   ") != std::string::npos)//fetches excitation energies from gaussian output
      {
        excitE.push_back(std::stof(buffer.substr(38)));
      }

      if (buffer.find(" SCF Done:") != std::string::npos)
      {
        e_total = std::stof(buffer.substr(buffer.find_first_of("=") + 1));
      }



      if (grad && buffer.find("Old X    -DE/DX   Delta X") != std::string::npos) //fetches last calculated gradients from output
      {
        coords::Cartesian_Point g;

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

          g_tmp[i] = g * HartreePerBohr2KcalperMolperAngstr;

        }
      }//end gradient reading

      if (grad && buffer.find("Center     Atomic      Atomic             Coordinates (Angstroms)") != std::string::npos)//reads last coordinates from file
      {
        coords::Cartesian_Point p;

        std::getline(in_file, buffer);
        std::getline(in_file, buffer);

        for (std::size_t i(0); i < atoms && !in_file.eof(); ++i)
        {
          std::getline(in_file, buffer);

          std::sscanf(buffer.c_str(), "%*s %*s %*s %lf %lf %lf", &p.x(), &p.y(), &p.z());


          if (grad || opt)
          {
            xyz_tmp[i] = p;
          }
        }
      }//end coordinater reading


    }//end while(!in_file.eof())

    if (grad || opt)
    {
      coords->set_xyz(std::move(xyz_tmp));
    }

    if (grad || opt)
    {
      coords->swap_g_xyz(g_tmp);
    }

    } //end if(in_file)

     std::sort(occMO.begin(), occMO.end(), std::greater <double>()); //sort occupied mos highest to lowest
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

    e_total *= au2kcal_mol;
    energy = e_total;

    //test output for interface, shound be outcommented

   /* for (unsigned int i = 0; i < xyz_tmp.size(); i++)
    { mos << xyz_tmp[i] << '\n'; }*/

   /* mos << "\n occ" << "       " << "virt \n";*/

   /* for (unsigned int i = 0; i < occMO.size(); i++)
    { mos << occMO[i] << "    " << virtMO[i] << '\n'; }*/

    /*mos << "\n total energy: " << e_total << '\n';*/

   /* mos << "\n Excitation energies \n";*/

    //for (unsigned int i = 0; i < excitE.size(); i++)
    //{ mos << excitE[i] << '\n'; }

    //mos << "\n Gradients \n";

   /* for (unsigned int i = 0; i < g_tmp.size(); i++)
    { mos << g_tmp[i] << '\n'; }*/

    /*for (unsigned int i = 0; i < state_i.size(); i++)
    { mos << state_i[i] << "   " << state_j[i] << "   " << ex_ex_trans[i] << '\n'; }*/

   /* for (unsigned int i = 0; i < gz_i_state.size(); i++)
    { mos << gz_i_state[i] << "   " << gz_ex_trans[i] << '\n'; }*/

   /* mos.close();*/
  }

int energy::interfaces::gaussian::sysCallInterfaceGauss::callGaussian()
{
  std::string gaussian_call = Config::get().energy.gaussian.path + " " + id + ".gjf";

  int ret = scon::system_call(gaussian_call);
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
  print_gaussianInput('e');

  if (callGaussian() == 0)
  {
    read_gaussianOutput(false, false);
  }
  else
  {
    if (Config::get().general.verbosity >=2)
    {
      std::cout << "Gaussian call return value was not 0. Treating structure as broken.e\n";
    }
    integrity = false;
  }
  return energy;
}

double energy::interfaces::gaussian::sysCallInterfaceGauss::g(void)
{
  std::string tmp_id = id;
  id = id + "_G_";

  integrity = true;
  print_gaussianInput('g');
  if (callGaussian() == 0) read_gaussianOutput(true, false);
  else
  {
    if (Config::get().general.verbosity >= 2)
    {
      std::cout << "Gaussian call return value was not 0. Treating structure as broken.g\n";
    }
    integrity = false;
  }

  id = tmp_id;

  return energy;
}

double energy::interfaces::gaussian::sysCallInterfaceGauss::h(void)
{
  if (Config::get().general.verbosity >= 2)
  {
    std::cout << "Hessian not implemented in CAST as yet.";
  }
  integrity = true;
  print_gaussianInput('h');
  if (callGaussian() == 0) read_gaussianOutput();
  else
  {
    if (Config::get().general.verbosity >= 2)
    {
      std::cout << "Gaussian call return value was not 0. Treating structure as broken.h\n";
    }
    integrity = false;
  }
  return energy;
}

double energy::interfaces::gaussian::sysCallInterfaceGauss::o(void)
{
  std::string tmp_id = id;
  id = id +"_O_";

  integrity = true;
  print_gaussianInput('o');
  if (callGaussian() == 0) read_gaussianOutput(true, true);
  else
  {
    if (Config::get().general.verbosity >= 2)
    {
      std::cout << "Gaussian call return value was not 0. Treating structure as broken.o\n";
    }
    integrity = false;
  }

  id=tmp_id;

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
  S << std::right << std::setw(20) << std::fixed << std::setprecision(16) << hof_kcal_mol;
  S << std::right << std::setw(20) << std::fixed << std::setprecision(16) << e_total;
  S << std::right << std::setw(20) << std::fixed << std::setprecision(16) << e_electron;
  S << std::right << std::setw(20) << std::fixed << std::setprecision(16) << e_core;
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
