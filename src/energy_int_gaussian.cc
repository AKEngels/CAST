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
#pragma warning (disable: 4996)
#endif
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

/*
Gaussian sysCall functions
*/

/**
* Global static instance of the sysCallInterfaceGauss-object.
*/

energy::interfaces::gaussian::sysCallInterfaceGauss::sysCallInterfaceGauss(coords::Coordinates* cp) :
  energy::interface_base(cp),
  hof_kcal_mol(0.0), hof_kj_mol(0.0), e_total(0.0),
  e_electron(0.0), e_core(0.0), failcounter(0u)
{
  id = Config::get().general.outputFilename;
  std::stringstream ss;
  std::srand(static_cast<unsigned>(std::time(0)));
  ss << (std::size_t(std::rand()) | (std::size_t(std::rand()) << 15));
  id.append("_tmp_").append(ss.str());
  optimizer = Config::get().energy.gaussian.opt;
  charge = std::stoi(Config::get().energy.gaussian.charge);
}

energy::interfaces::gaussian::sysCallInterfaceGauss::sysCallInterfaceGauss(sysCallInterfaceGauss const& rhs, coords::Coordinates* cobj) :
  interface_base(cobj),
  hof_kcal_mol(rhs.hof_kcal_mol), hof_kj_mol(rhs.hof_kj_mol), e_total(rhs.e_total),
  e_electron(rhs.e_electron), e_core(rhs.e_core), failcounter(rhs.failcounter)
{
  id = rhs.id;
  charge = rhs.charge;
  interface_base::operator=(rhs);
}

energy::interface_base* energy::interfaces::gaussian::sysCallInterfaceGauss::clone(coords::Coordinates* coord_object) const
{
  sysCallInterfaceGauss* tmp = new sysCallInterfaceGauss(*this, coord_object);
  return tmp;
}

energy::interface_base* energy::interfaces::gaussian::sysCallInterfaceGauss::move(coords::Coordinates* coord_object)
{
  sysCallInterfaceGauss* tmp = new sysCallInterfaceGauss(*this, coord_object);
  return tmp;
}

void energy::interfaces::gaussian::sysCallInterfaceGauss::swap(interface_base& rhs)
{
  swap(dynamic_cast<sysCallInterfaceGauss&>(rhs));
}

void energy::interfaces::gaussian::sysCallInterfaceGauss::swap(sysCallInterfaceGauss& rhs)
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
    if (Config::get().energy.gaussian.chk.length() != 0) {    // if checkpoint file specified
      out_file << "%chk=" << Config::get().energy.gaussian.chk << "\n";
    }
    out_file << "# " << Config::get().energy.gaussian.method << " " << Config::get().energy.gaussian.basisset;
    out_file << " " << Config::get().energy.gaussian.spec << " ";
    if (Config::get().energy.gaussian.cpcm == true) out_file << "scrf(cpcm,solvent=generic,read) ";
    if (Config::get().energy.qmmm.mm_charges.size() != 0) out_file << "Charge ";
    if (Config::get().energy.qmmm.use == true) out_file << "NoSymm ";

    switch (calc_type) {// to ensure the needed gaussian keywords are used in gausian inputfile for the specified calculation
    case 'o':
      if (Config::get().energy.gaussian.steep)
      {
        out_file << " Opt=(Steep) ";//Cartesian,
      }
      else
      {
        out_file << " Opt"; //=Cartesian
      }
      break;
    case 'g':
      out_file << " Force";
      if (Config::get().energy.qmmm.mm_charges.size() != 0) out_file << " Prop=(Field,Read) Density";
      break;
    }

    out_file << '\n';
    out_file << '\n';
    out_file << Config::get().general.outputFilename;
    out_file << '\n';
    out_file << '\n';
    out_file << charge << " ";
    out_file << Config::get().energy.gaussian.multipl;
    out_file << '\n';
    out_file << coords::output::formats::xyz(*coords);
    if (Config::get().periodics.periodic)
    {
      out_file << std::left << std::setw(3) << "TV";   // Translation Vector for unit cell (x-direction)
      out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << Config::get().periodics.pb_box.x();
      out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << 0.0;
      out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << 0.0;
      out_file << '\n';

      out_file << std::left << std::setw(3) << "TV";   // Translation Vector for unit cell (y-direction)
      out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << 0.0;
      out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << Config::get().periodics.pb_box.y();
      out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << 0.0;
      out_file << '\n';

      out_file << std::left << std::setw(3) << "TV";   // Translation Vector for unit cell (z-direction)
      out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << 0.0;
      out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << 0.0;
      out_file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(7) << Config::get().periodics.pb_box.z();
      out_file << '\n';
    }
    out_file << '\n';
    if (Config::get().energy.qmmm.mm_charges.size() != 0)  // if desired: writing additional point charges (from MM atoms)
    {
      for (auto& c : Config::get().energy.qmmm.mm_charges)
      {
        out_file << c.x << " " << c.y << " " << c.z << " " << c.scaled_charge << "\n";
      }
      out_file << '\n';
    }
    if (Config::get().energy.gaussian.cpcm == true)    // parameters for CPCM
    {
      out_file << "eps=" << Config::get().energy.gaussian.eps << "\n";
      out_file << "epsinf=" << Config::get().energy.gaussian.epsinf << "\n";
      out_file << "\n";
    }
    if (Config::get().energy.gaussian.method == "DFTB=read")  // slater koster files for DFTB
    {
      std::vector<std::vector<std::string>> pairs = find_pairs(*coords);
      for (auto p : pairs)
      {
        std::string filename = p[0] + "-" + p[1] + ".skf";
        if (file_exists(filename) == false)
        {
          throw std::runtime_error("ERROR! Slater Koster file " + filename + " does not exist. Please download it from dftb.org and convert it with the task MODIFY_SK_FILES!");
        }
        out_file << "@./" << filename << " /N\n";
      }
      out_file << '\n';
    }
    else if (Config::get().energy.gaussian.method == "DFTBA")  // stuff for DFTBA
    {
      out_file << "@GAUSS_EXEDIR:dftba.prm\n\n";
    }
    if (calc_type == 'g' && Config::get().energy.qmmm.mm_charges.size() != 0)   // writing points for electric field (positions of MM atoms)
    {
      for (auto& c : Config::get().energy.qmmm.mm_charges)
      {
        out_file << c.x << " " << c.y << " " << c.z << "\n";
      }
      out_file << '\n';
    }
    out_file.close();
  }
  else throw std::runtime_error("Writing Gaussian Inputfile failed.");
}

bool energy::interfaces::gaussian::sysCallInterfaceGauss::read_gaussianOutput(bool const grad, bool const opt, bool const qmmm)
{
  //std::ofstream mos("MOs.txt", std::ios_base::out); //ofstream for mo testoutput keep commented if not needed

  hof_kcal_mol = hof_kj_mol = energy = e_total = e_electron = e_core = 0.0;
  double mm_el_energy(0.0);
  int atoms(coords->size());

  auto in_string = id + ".log";
  std::ifstream in_file(in_string.c_str(), std::ios_base::in);

  bool test_lastMOs(false);//to controll if reading was successfull
  coords::Representation_3D g_tmp(coords->size()), xyz_tmp(coords->size());

  if (in_file)
  {
    std::string buffer;
    bool normalGaussianTermination = false;
    while (!in_file.eof())
    {

      std::getline(in_file, buffer);

      if (buffer.find("NAtoms=") != std::string::npos)   // get total number of atoms (in case of bonded QM/MM including link atoms)
      {                                                  // without QM/MM it's the same as coords->size()
        std::vector<std::string> stringvec = split(buffer, ' ', true);
        atoms = std::stoi(stringvec[1]);
      }

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
        for (auto i = 0u; buffer.length() > (29 + i * 10); i++) //in gaussian output orbital energies are presented in rows of 5
        {
          occMO.push_back(std::stof(buffer.substr(29 + i * 10)));
        }
      }

      if (buffer.find("Alpha virt. eigenvalues --") != std::string::npos)
      {
        for (auto i = 0u; buffer.length() > (29 + i * 10); i++)
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

          if (buffer.find("Excited to excited state transition velocity") != std::string::npos) { el_dipm = false; }
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
        e_total = std::stod(buffer.substr(buffer.find_first_of("=") + 1));
      }

      if (qmmm)
      { // for QM/MM calculation: get energy of the interaction between the external charges (i.e. the MM atoms)
        if (buffer.find("Self energy of the charges") != std::string::npos)
        {
          mm_el_energy = std::stod(buffer.substr(buffer.find_first_of("=") + 1, 21));
        }

        // get the electric field
        // the electric field at MM atoms due to QM atoms is used to calculate gradients of electrostatic interaction on MM atoms
        if (buffer.find("-------- Electric Field --------") != std::string::npos)
        {
          coords::Cartesian_Point p;
          std::vector<coords::Cartesian_Point> el_field_tmp;
          std::getline(in_file, buffer);
          std::getline(in_file, buffer);

          std::getline(in_file, buffer);
          while (buffer.substr(0, 15) != " --------------")
          {
            p.x() = std::stod(buffer.substr(24, 14));
            p.y() = std::stod(buffer.substr(38, 14));
            p.z() = std::stod(buffer.substr(52, 14));

            el_field_tmp.push_back(p * energy::Hartree_Bohr2Kcal_MolAng);
            std::getline(in_file, buffer);
          }
          calc_grads_from_field(el_field_tmp);   // calculate the gradients on external charges from electric field
        }
      }

      if (grad && buffer.find("Old X    -DE/DX   Delta X") != std::string::npos) //fetches last calculated gradients from output
      {
        coords::Cartesian_Point g;
        double temp;

        std::getline(in_file, buffer);
        std::getline(in_file, buffer);

        for (std::size_t i(0); i < (*this->coords).size() && !in_file.eof(); ++i)
        {
          std::sscanf(buffer.c_str(), "%*s %*s %lf %*s %*s %*s %*s", &temp);
          g.x() = -temp;
          std::getline(in_file, buffer);
          std::sscanf(buffer.c_str(), "%*s %*s %lf %*s %*s %*s %*s", &temp);
          g.y() = -temp;
          std::getline(in_file, buffer);
          std::sscanf(buffer.c_str(), "%*s %*s %lf %*s %*s %*s %*s", &temp);
          g.z() = -temp;

          std::getline(in_file, buffer);
          g_tmp[i] = g * energy::Hartree_Bohr2Kcal_MolAng;
        }
      }//end gradient reading

      if (grad && buffer.find("Center     Atomic      Atomic             Coordinates (Angstroms)") != std::string::npos)//reads last coordinates from file
      {
        coords::Cartesian_Point p;

        std::getline(in_file, buffer);
        std::getline(in_file, buffer);

        for (auto i(0); i < atoms && !in_file.eof(); ++i)
        {
          std::getline(in_file, buffer);

          std::sscanf(buffer.c_str(), "%*s %*s %*s %lf %lf %lf", &p.x(), &p.y(), &p.z());


          if (grad || opt)
          {
            xyz_tmp[i] = p;
          }
        }
      }//end coordinater reading

      if (buffer.find("Mulliken charges:") != std::string::npos)  // read charges (restricted calculation)
      {
        double charge;
        atom_charges.clear();

        std::getline(in_file, buffer); // discard next line
        for (std::size_t i(0); i < coords->size(); ++i)
        {
          std::getline(in_file, buffer);
          std::sscanf(buffer.c_str(), "%*s %*s %lf", &charge);
          atom_charges.push_back(charge);
        }
      }

      if (buffer.find("Mulliken charges and spin densities:") != std::string::npos)  // read charges (unrestricted calculation)
      {
        double charge;
        atom_charges.clear();

        std::getline(in_file, buffer); // discard next line
        for (std::size_t i(0); i < coords->size(); ++i)
        {
          std::getline(in_file, buffer);
          std::sscanf(buffer.c_str(), "%*s %*s %lf %*s", &charge);
          atom_charges.push_back(charge);
        }
      }

      if (buffer.find("Normal termination of Gaussian") != std::string::npos)
      {
        normalGaussianTermination = true;
      }

    }//end while(!in_file.eof())

    if (normalGaussianTermination == false)
    {
      if (Config::set().energy.gaussian.delete_input == false)
      {                                   // save logfile for failed gaussian calls
        failcounter++;
        std::string oldname = id + ".log";
        std::string newname = "failed_gaussian_call_" + std::to_string(failcounter) + ".log";
        rename(oldname.c_str(), newname.c_str());
      }
      return false; // GAUSSIAN DID NOT TERMINATE PROPERLY
    }

    if (qmmm)
    {
      e_total = e_total - mm_el_energy;
    }

    if (grad && opt)
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

  for (auto i = 0u; i < occMO.size(); i++)
  {
    occMO[i] *= energy::au2kcal_mol;
  }

  for (std::size_t i = 0; i < virtMO.size(); i++)
  {
    virtMO[i] *= energy::au2kcal_mol;
  }

  for (std::size_t i = 0; i < excitE.size(); i++)
  {
    excitE[i] *= energy::eV2kcal_mol;
  }

  e_total *= energy::au2kcal_mol;
  energy = e_total;

  //test output for interface, should be commented, only for debug

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
  return true;
}

int energy::interfaces::gaussian::sysCallInterfaceGauss::callGaussian()
{
  std::string gaussian_call = "export GAUSS_SCRDIR=" + fs::current_path().string() + " && " + Config::get().energy.gaussian.path + " " + id + ".gjf";

  const int ret = scon::system_call(gaussian_call);
  if (ret != 0)
  {
    ++failcounter;
    std::cout << "Warning Gaussian call to '" << gaussian_call
      << " ' did not return 0 and probably failed.\n( A total of "
      << failcounter << " Gaussian call" << (failcounter != 1 ? "s have" : " has")
      << " failed so far.)\n";

    if (Config::set().energy.gaussian.delete_input == false)
    {                                   // save logfile for failed gaussian calls
      std::string oldname = id + ".log";
      std::string newname = "failed_gaussian_call_" + std::to_string(failcounter) + ".log";
      rename(oldname.c_str(), newname.c_str());
    }

    if (failcounter > Config::get().energy.gaussian.maxfail && Config::get().energy.qmmm.use == false)
    {
      throw std::runtime_error("More than " + std::to_string(Config::get().energy.gaussian.maxfail) + " Gaussian calls have failed.");
    }
  }
  return ret;
}

//Energy functions
double energy::interfaces::gaussian::sysCallInterfaceGauss::e(void)
{
  integrity = coords->check_structure();
  if (integrity == true)
  {
    print_gaussianInput('e');

    if (callGaussian() == 0)
    {
      if (!read_gaussianOutput(false, false, Config::get().energy.qmmm.use))
      {
        throw std::runtime_error("Gaussian calculation did not terminate normally.");
      }
    }
    else
    {
      if (Config::get().general.verbosity >= 2)
      {
        std::cout << "Gaussian call (e) return value was not 0. Treating structure as broken.\n";
      }
      integrity = false;
      return 0.0;
    }
    return energy;
  }
  else return 0;
}

double energy::interfaces::gaussian::sysCallInterfaceGauss::g(void)
{
  integrity = coords->check_structure();

  if (integrity == true)
  {
    std::string tmp_id = id;
    id = id + "_G_";


    print_gaussianInput('g');
    if (callGaussian() == 0)
    {
      if (!read_gaussianOutput(false, false, Config::get().energy.qmmm.use))
      {
        throw std::runtime_error("Gaussian calculation did not terminate normally.");
      }
    }
    else
    {
      if (Config::get().general.verbosity >= 2)
      {
        std::cout << "Gaussian call (g) return value was not 0. Treating structure as broken.\n";
      }
      integrity = false;
    }

    id = tmp_id;

    return energy;
  }
  else return 0;
}

double energy::interfaces::gaussian::sysCallInterfaceGauss::h(void)
{

  throw std::runtime_error("Hessian not implemented in CAST as yet.");

  /*integrity = true;
  print_gaussianInput('h');
  if (callGaussian() == 0) read_gaussianOutput();
  else
  {
    if (Config::get().general.verbosity >= 2)
    {
      std::cout << "Gaussian call (h) return value was not 0. Treating structure as broken.\n";
    }
    integrity = false;
  }
  return energy;*/
}

double energy::interfaces::gaussian::sysCallInterfaceGauss::o(void)
{
  integrity = coords->check_structure();

  if (integrity == true)
  {
    std::string tmp_id = id;
    id = id + "_O_";


    print_gaussianInput('o');
    if (callGaussian() == 0) read_gaussianOutput(true, true);
    else
    {
      if (Config::get().general.verbosity >= 2)
      {
        std::cout << "Gaussian call (o) return value was not 0. Treating structure as broken.\n";
      }
      integrity = false;
    }

    id = tmp_id;

    return energy;
  }
  else return 0;
}

//Outputfunctions

void energy::interfaces::gaussian::sysCallInterfaceGauss::print_E(std::ostream& S) const
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

void energy::interfaces::gaussian::sysCallInterfaceGauss::print_E_head(std::ostream& S, bool const endline) const
{
  S << std::right << std::setw(24) << "HOF";
  S << std::right << std::setw(24) << "EL";
  S << std::right << std::setw(24) << "CORE";
  S << std::right << std::setw(24) << "TOT";
  if (endline) S << '\n';
}

void energy::interfaces::gaussian::sysCallInterfaceGauss::print_E_short(std::ostream& S, bool const endline) const
{
  S << '\n';
  S << std::right << std::setw(24) << hof_kcal_mol;
  S << std::right << std::setw(24) << e_electron;
  S << std::right << std::setw(24) << e_core;
  S << std::right << std::setw(24) << e_total;
  if (endline) S << '\n';
}

void energy::interfaces::gaussian::sysCallInterfaceGauss::to_stream(std::ostream&) const { }

void energy::interfaces::gaussian::sysCallInterfaceGauss::calc_grads_from_field(std::vector<coords::Cartesian_Point> const& electric_field)
{
  if (electric_field.size() - coords->size() != Config::get().energy.qmmm.mm_charges.size())
  {
    throw std::logic_error("Electric field has not the same size as the external charges. Can't calculate gradients.");
  }

  grad_ext_charges.clear();
  for (auto i = coords->size(); i < electric_field.size(); ++i)  // calculate gradients on external charges from electric field
  {
    coords::Cartesian_Point E = electric_field[i];
    double q = Config::get().energy.qmmm.mm_charges[i - coords->size()].scaled_charge;

    double x = q * E.x();
    double y = q * E.y();
    double z = q * E.z();
    coords::Cartesian_Point new_grad;
    new_grad.x() = -x;
    new_grad.y() = -y;
    new_grad.z() = -z;

    grad_ext_charges.push_back(new_grad);
  }
}

std::vector<double>
energy::interfaces::gaussian::sysCallInterfaceGauss::charges() const
{
  return atom_charges;
}

std::vector<coords::Cartesian_Point>
energy::interfaces::gaussian::sysCallInterfaceGauss::get_g_ext_chg() const
{
  return grad_ext_charges;
}

