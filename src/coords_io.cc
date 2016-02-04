#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstddef>
#include <cstdio>
#include "atomic.h"
#include "global.h"
#include "configuration.h"
#include "coords_io.h"
#include "scon_utility.h"

#if defined(_MSC_VER) && !defined(CAST_SSCANF_COORDS_IO)
#define CAST_SSCANF_COORDS_IO sscanf_s
#elif !defined(CAST_SSCANF_COORDS_IO)
#define CAST_SSCANF_COORDS_IO sscanf
#endif

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
    default:
    {
      return new formats::tinker;
    }
  }
}

// Amber is aimed at FORTRAN parsers so this code is just, horrible... just horrible...
// See http://ambermd.org/prmtop.pdf
// Since AMBER Stores actual coordinates in a seperate file, this procedure will
// read from the config-namespace. We search for the name of a possibly existing second file
// containing the raw coordinates. Pay attantion to this!
coords::Coordinates coords::input::formats::amber::read(std::string file)
{
  Coordinates coord_object;
  std::ifstream config_file_stream(file.c_str(), std::ios_base::in);
  std::string line;
  std::cout << "NOTE: You are specifiying AMBER-Files as input.\n This feature is highly, let me repeated that, HIGHLY experimental in CAST\n";
  std::cout << "All that fancy AMBER-stuff like residues and such is OMITTED. We collect only the atoms and their positions.\n";
  std::cout << "Keep in mind: If anything breaks, you get to keep all the pieces!\n";

  //First, let's see which sections are present at all. Then we can judge from where we'll take our information.
  if (config_file_stream)
  {
    //%VERSION
    std::getline(config_file_stream, line);
    if (line.substr(1, 7) != std::string("VERSION"))
    {
      std::cout << "In reading AMBER prmtop: First Line needs to be %VERSION.\n";
      throw("In reading AMBER prmtop: First Line needs to be %VERSION.\n");
    }

    while (std::getline(config_file_stream, line))
    {
      if (line.substr(0, 5) == "%FLAG")
      {
        for (unsigned int i(0U); i < sections_size; ++i)
        {
          // remove "%FORMAT%" and search
          if (line.substr(6).find(amber_sections[i]) != std::string::npos)
          {
            sectionsPresent.push_back(i);
            break;
          }
        }
      }
    }


    if (sectionsPresent.size() == 0)
    {
      throw("found no sections when reading AMBER input sorry \n");
    }
  }

  // Reopen file for second run, sicne we now know which sections are present
  config_file_stream.clear();
  config_file_stream.close();
  config_file_stream.open(file.c_str(), std::ios_base::in);

  // First, lets see if we can infer atom types from section
  // "AMBER_Atom_Types". If it isn't present, we will have
  // to take "ATOM_NAME". If it is present, we can skip past "ATOM_NAME"

  // FUTURE: Use Atomic_Number, however, in the workflow of this
  // current working group this is not present. (ie: this scetion is just not present
  // in the files for which processing with CAST is desired).
  bool skipPastAtomName = std::find(sectionsPresent.begin(), sectionsPresent.end(), 34u) != sectionsPresent.end();

  //NOW LETS GO!
  if (config_file_stream)
  {
    //get %VERSION -> discard
    std::getline(config_file_stream, line);

    unsigned int currentSectionID = 515845u; //self explaining, random high number as starting value
    std::string currentSectionFormat; //FORTRAN style format specifiers

    // LOOP, reading file
    while (std::getline(config_file_stream, line))
    {
      //Check if line starts with "%", else continue
      if (!(line.substr(0, 1) == std::string("%")))
      {
        continue;
      }

      //Check if %COMMENT, if yes then get nextline
      if (line.substr(0, 8) == std::string("%COMMENT"))
      {
        continue;
      }

      //Search for a valid section
      for (unsigned int i(0U); i < sections_size; ++i)
      {
        // remove "%FLAG%" and search
        if (line.substr(6).find(amber_sections[i]) != std::string::npos)
        {
          //Get the format specifier for debug purposes.
          std::getline(config_file_stream, currentSectionFormat);
          currentSectionID = i;
          break;
        }
        if (i == sections_size - 1u)
        {
          goto FAILED;
        }
      }

      //If section == TITLE
      if (currentSectionID == 0u)
      {
        std::getline(config_file_stream, this->name);
      }

      //If section == POINTERS <- very important stuff in there.
      if (currentSectionID == 1u)
      {
        std::getline(config_file_stream, line);
        numberOfAtoms = std::stoi(line.substr(0, 8));
        for (unsigned int i = 0u; i < 10u; i++)
        {
          pointers_raw[i + 0u] = std::stoi(line.substr(8u * i, 8));
        }
        std::getline(config_file_stream, line);
        for (unsigned int i = 0u; i < 10u; i++)
        {
          pointers_raw[i + 10u] = std::stoi(line.substr(8u * i, 8));
        }
        std::getline(config_file_stream, line);
        for (unsigned int i = 0u; i < 10u; i++)
        {
          pointers_raw[i + 20u] = std::stoi(line.substr(8u * i, 8));
        }
        std::getline(config_file_stream, line);
        pointers_raw[30u] = std::stoi(line.substr(0, 8));
      }

      // Now we know some stuff, like how many atoms there are. 
      // So we can start allocating.

      //If section == BONDS_INC_HYDROGEN
      if (currentSectionID == 24u)
      {
        unsigned long state = 1u; //We dont care bout the third value, its AMBER stuff...
        for (unsigned int countlines = 0u; countlines < std::ceil(float(pointers_raw[2]) * 3.f / 10.f) - 1u; countlines++)
        {
          std::getline(config_file_stream, line);
          // A line has 10 members in the format
          for (unsigned int i = 0; i < 10u; state++, i++)
          {
            //We dont care bout the third value, its AMBER stuff...
            if (state % 3u != 0)
            {
              // Citing the amber file format manual:
              // "For run-time efficiency, the atom indexes are actually
              // indexes into a coordinate array, so the actual atom index A is calculated
              // from the coordinate array index N by A = N / 3 + 1. (N is the value in the
              // topology file)"
              bondsWithHydrogen.push_back(std::stoi(line.substr(i * 8u, 8)) / 3u + 1u);
            }
          }
        }
      }

      //If section == BONDS_WITHOUT_HYDROGEN
      if (currentSectionID == 25u)
      {
        unsigned long state = 1u; //We dont care bout the third value, its AMBER stuff...
        for (unsigned int countlines = 0u; countlines < std::ceil(float(pointers_raw[3]) * 3.f / 10.f) - 1u; countlines++)
        {
          std::getline(config_file_stream, line);
          // A line has 10 members in the format
          for (unsigned int i = 0; i < 10u; state++, i++)
          {
            //We dont care bout the third value, its AMBER stuff...
            if (state % 3u != 0)
            {
              // Citing the amber file format manual:
              // "For run-time efficiency, the atom indexes are actually
              // indexes into a coordinate array, so the actual atom index A is calculated
              // from the coordinate array index N by A = N / 3 + 1. (N is the value in the
              // topology file)"
              bondsWithoutHydrogen.push_back(std::stoi(line.substr(i * 8u, 8)) / 3u + 1u);
            }
          }
        }
      }

      //If section == ATOM_NAME
      if (!skipPastAtomName && currentSectionID == 2u)
      {
        std::cout << "Did not find section AMBER_Atom_Types in .prmtop file!\n";
        std::cout << "Continuing, now we will try to get the names from the section ATOM_NAME.\n";
        for (unsigned int countlines = 0u; countlines < std::ceil(float(numberOfAtoms) / 20.f); countlines++)
        {
          std::getline(config_file_stream, line); // Discard the format specifier
          // A Line has 20 members in this format
          for (unsigned int i = 0u; i < 20u; i++)
          {
            std::string symbol = line.substr(i * 4u, 2u);

            // First, let's process all those AMBER Atom Types with
            // just incredibly retarded names. damn.
            // via: http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html
            if (symbol == "CU")
            {
              symbol = "Cu";
            }
            else if (symbol == "C�") symbol = "Co";
            else if (symbol == "IM") symbol = "Cl";
            else if (symbol == "QC") symbol = "Cs";
            else if (symbol == "QK") symbol = "K";
            else if (symbol == "QL") symbol = "Li";
            else if (symbol == "QN") symbol = "Na";
            else if (symbol == "QR") symbol = "Na";

            //wtf is lone pair? certainly not an atom!
            else if (symbol == "LP") continue;
            else if (symbol == "nu") continue;

            //Si? Starts with an s like sulfur, so we have to check..gnaah.....
            else if (symbol == "si" || symbol == "Si" || symbol == "SI") { symbol = "Si"; }

            //Now, lets see
            else if (symbol.substr(0, 1) == "h" || symbol.substr(0, 1) == "H") symbol = "H";
            else if (symbol.substr(0, 1) == "c" || symbol.substr(0, 1) == "C") symbol = "C";
            else if (symbol.substr(0, 1) == "o" || symbol.substr(0, 1) == "O") symbol = "O";
            else if (symbol.substr(0, 1) == "n" || symbol.substr(0, 1) == "N") symbol = "N";
            else if (symbol.substr(0, 1) == "s" || symbol.substr(0, 1) == "S") symbol = "S";
            else if (symbol.substr(0, 1) == "p" || symbol.substr(0, 1) == "P") symbol = "P";
            else
            {
              std::cout << "Could not identify AMER_ATOM_TYPE: " << line.substr(i * 4u, 4u) << "\nstopping!\n";
              goto FAILED;
            }

            Atom current(symbol);
            atoms.add(current);
          }
        }
      }

      //If section == AMBER_Atom_Types
      if (currentSectionID == 34u)
      {
        // Discard the format specifier
        // Read all lines except last one
        for (unsigned int countlines = 0u; countlines < std::ceil(float(numberOfAtoms) / 20.f) - 1; countlines++)
        {
          std::getline(config_file_stream, line);
          // A Line has 20 members in this format
          for (unsigned int i = 0u; i < 20u; i++)
          {
            std::string symbol = line.substr(i * 4u, 2u);

            // First, let's process all those AMBER Atom Types with
            // just incredibly retarded names. damn.
            // via: http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html
            if (symbol == "CU")
            {
              symbol = "Cu";
            }
            else if (symbol == "C�") symbol = "Co";
            else if (symbol == "IM") symbol = "Cl";
            else if (symbol == "QC") symbol = "Cs";
            else if (symbol == "QK") symbol = "K";
            else if (symbol == "QL") symbol = "Li";
            else if (symbol == "QN") symbol = "Na";
            else if (symbol == "QR") symbol = "Na";
            else if (symbol == "Na") symbol = "Na";
            else if (symbol == "IP")
            {
              symbol = "Na";
              std::cout << "Atom Type \"IP\" used. We are guessing it stands for sodium.\n This is not documented in the official "
                << "AMBER Atom Types webpage (\"http://ambermd.org/antechamber/gaff.html\" or "
                << "\"http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html\").\n"
                << "So you need to check if it really IS meant to be sodium.\n"
                << "Please talk to a CAST programmer about this.\n";
            }

            //wtf is lone pair? certainly not an atom!
            else if (symbol == "LP") continue;
            else if (symbol == "nu") continue;

            //Si? Starts with an s like sulfur, so we have to check..gnaah.....
            else if (symbol == "si" || symbol == "Si" || symbol == "SI") { symbol = "Si"; }

            //Now, lets see
            else if (symbol.substr(0, 1) == "h" || symbol.substr(0, 1) == "H") symbol = "H";
            else if (symbol.substr(0, 1) == "c" || symbol.substr(0, 1) == "C") symbol = "C";
            else if (symbol.substr(0, 1) == "o" || symbol.substr(0, 1) == "O") symbol = "O";
            else if (symbol.substr(0, 1) == "n" || symbol.substr(0, 1) == "N") symbol = "N";
            else if (symbol.substr(0, 1) == "s" || symbol.substr(0, 1) == "S") symbol = "S";
            else if (symbol.substr(0, 1) == "p" || symbol.substr(0, 1) == "P") symbol = "P";

            //Maybe we are already done?
            //FAILED
            else if (!symbol.empty())
            {
              std::cout << "Could not identify AMER_ATOM_TYPE: " << line.substr(i * 4u, 4u) << "\nstopping!\n";
              std::cout << "You should probably go talk to a CAST programmer about this.\n";
              goto FAILED;
            }

            Atom current(symbol);
            atoms.add(current);
          }
        }
        //Process last line
        std::getline(config_file_stream, line);
        for (unsigned int i = 0u; i < numberOfAtoms % 20; i++)
        {
          std::string symbol = line.substr(i * 4u, 2u);

          // First, let's process all those AMBER Atom Types with
          // just incredibly retarded names. damn.
          // via: http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html
          if (symbol == "CU")
          {
            symbol = "Cu";
          }
          else if (symbol == "C�") symbol = "Co";
          else if (symbol == "IM") symbol = "Cl";
          else if (symbol == "QC") symbol = "Cs";
          else if (symbol == "QK") symbol = "K";
          else if (symbol == "QL") symbol = "Li";
          else if (symbol == "QN") symbol = "Na";
          else if (symbol == "QR") symbol = "Na";
          else if (symbol == "Na") symbol = "Na";
          else if (symbol == "IP")
          {
            symbol = "Na";
            std::cout << "Atom Type \"IP\" used. We are guessing it stands for sodium.\n This is not documented in the official "
              << "AMBER Atom Types webpage (\"http://ambermd.org/antechamber/gaff.html\" or "
              << "\"http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html\").\n"
              << "So you need to check if it really IS meant to be sodium.\n"
              << "Please talk to a CAST programmer about this.\n";
          }

          //wtf is lone pair? certainly not an atom!
          else if (symbol == "LP") continue;
          else if (symbol == "nu") continue;

          //Si? Starts with an s like sulfur, so we have to check..gnaah.....
          else if (symbol == "si" || symbol == "Si" || symbol == "SI") { symbol = "Si"; }

          //Now, lets see
          else if (symbol.substr(0, 1) == "h" || symbol.substr(0, 1) == "H") symbol = "H";
          else if (symbol.substr(0, 1) == "c" || symbol.substr(0, 1) == "C") symbol = "C";
          else if (symbol.substr(0, 1) == "o" || symbol.substr(0, 1) == "O") symbol = "O";
          else if (symbol.substr(0, 1) == "n" || symbol.substr(0, 1) == "N") symbol = "N";
          else if (symbol.substr(0, 1) == "s" || symbol.substr(0, 1) == "S") symbol = "S";
          else if (symbol.substr(0, 1) == "p" || symbol.substr(0, 1) == "P") symbol = "P";

          //Maybe we are already done?
          //FAILED
          else if (!symbol.empty())
          {
            std::cout << "Could not identify AMER_ATOM_TYPE: " << line.substr(i * 4u, 4u) << "\nstopping!\n";
            std::cout << "You should probably go talk to a CAST programmer about this.\n";
            goto FAILED;
          }

          Atom current(symbol);
          atoms.add(current);

        }
      }

    }

    // Now lets implement who's bonded to who. 
    // Good thing we kept track of this with our two std::vector members
    for (unsigned int i = 0u; i < bondsWithHydrogen.size(); i = i + 2u)
    {
      atoms.atom(i).bind_to(i + 1u);
    }
    for (unsigned int i = 0u; i < bondsWithoutHydrogen.size(); i = i + 2u)
    {
      atoms.atom(i).bind_to(i + 1u);
    }

    // OK, let's now fetch the actual coordinates
    Representation_3D positions;
    positions.reserve(numberOfAtoms); //Reserve space
    if (!Config::get().io.amber_mdcrd.empty())
    {
      std::ifstream coord_file_stream(Config::get().io.amber_mdcrd.c_str(), std::ios_base::in);
      // Now, discard the title
      std::getline(coord_file_stream, line);

      //In mdcrd we need to check if a linebreak occurs and box coordiantes are written.
      if (!Config::get().io.amber_trajectory_at_constant_pressure)
      {
        unsigned long state = 0u; //Counts each processed floating point number.
        while (std::getline(coord_file_stream, line))
        {
          // A line has 10 members in the format
          for (unsigned int i = 0; i < 10u; state++, i++)
          {
            // Needed: Check if substr is empty
            if (line.substr(i * 8u, 8u).empty())
            {
              if (!(state % 3u == 0))
              {
                std::cout << "Encounterd unexpected linebreak or EOF in reading AMBER mdcrd file.\n";
                goto FAILED;
              }
              else
              {
                ///////////////////////////////
                // EXIT WHILE AND DO LOOP!!! //    <- THIS HERE IS IMPORTANT!!!!!!!
                ///////////////////////////////
                goto DONE;
              }
            }

            if (state % 3u == 0)
            {
              position.x() = std::stod(line.substr(i * 8u, 8u));
            }
            else if (state % 3u == 1)
            {
              position.y() = std::stod(line.substr(i * 8u, 8u));
            }
            else if (state % 3u == 2)
            {
              position.z() = std::stod(line.substr(i * 8u, 8u));
              positions.push_back(position);
              //Check if we reached end of structure
              if (((state + 1) / 3) % numberOfAtoms == 0)
              {
                input_ensemble.push_back(positions);
                if (positions.size() != atoms.size())
                {
                  throw std::logic_error("The size of an provided structure does not match the number of atoms.");
                }
                positions.clear();
                state++;
                break;
              }
            }
          }
        }
      }
      else
      {
        std::cout << "AMBER trajectories with constant pressure not yet supported, talk to your admin.";
        goto FAILED;
      }
    }
    else if (!Config::get().io.amber_inpcrd.empty())
    {
      std::cout << "Reading from inpcrd not yet supported.";
      goto FAILED;
    }
    else if (!Config::get().io.amber_restrt.empty())
    {
      std::cout << "Reading from restrt not yet supported.";
      goto FAILED;
    }
    else
    {
      std::cout << "No file containing AMBER coordinates was specified.";
      goto FAILED;
    }

    // THIS HAPPENS WHEN DONE
  DONE:
    if (input_ensemble.empty()) throw std::logic_error("No structures found.");
    coords::PES_Point x(input_ensemble[0u]);
    coord_object.init_swap_in(atoms, x);
    //i do not know what this is for.
    //however, it slows io down considerably. Thats why we are not doing it (whatever it is) atm.
    /*for (auto & p : input_ensemble)
    {
      p.gradient.cartesian.resize(p.structure.cartesian.size());
      coord_object.set_xyz(p.structure.cartesian);
      coord_object.to_internal();
      p = coord_object.pes();
    }*/
  }
  else
  {
    // IF FAILED!!
  FAILED:
    throw std::logic_error("Reading the structure input file failed.");
    //return Coordinates();
  }
  return coord_object;
}



coords::Coordinates coords::input::formats::tinker::read(std::string file)
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
    if (N == 0U) throw std::logic_error("ERR_COORD: Expecting no atoms from '" + file + "'.");
    Representation_3D positions;
    std::vector<std::size_t> index_of_atom(N);
    bool indexation_not_contiguous(false),
      has_in_out_subsystems(false);

    // loop fetching atoms and positions
    for (std::size_t i(1U); std::getline(coord_file_stream, line); ++i)
    {
      //std::cout << "Line " << i << " mod: " << i%(N+1u) << lineend;
      //std::istringstream linestream(line);
      if (i <= N)
      {
        std::istringstream linestream(line);
        std::size_t const nbmax(7u);
        tinker::line tfl;
        linestream >> tfl;
        Atom atom(tfl.symbol);
        index_of_atom[tfl.index - 1] = positions.size();
        if (positions.size() != (tfl.index - 1))
        {
          indexation_not_contiguous = true;
        }
        for (std::size_t j(0u); j < nbmax && tfl.bonds[j] > 0u; ++j)
        {
          atom.bind_to(tfl.bonds[j] - 1u);
          //coord_object.topo((tfl.bonds[j] - 1u), i - 1); 
        }
        positions.push_back(tfl.position);
        if (scon::find_substr_ci(line, "in") != std::string::npos)
        {
          atom.set_sub_type(Atom::ST_IN);
          atom.assign_to_system(1u);
          has_in_out_subsystems = true;
        }
        else if (scon::find_substr_ci(line, "out") != std::string::npos)
        {
          atom.set_sub_type(Atom::ST_OUT);
          atom.assign_to_system(2u);
          has_in_out_subsystems = true;
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
        if (i % (N + 1u) != 0)
        {
          double x(0), y(0), z(0);
          CAST_SSCANF_COORDS_IO(line.c_str(), "%*lu %*s %lf %lf %lf", &x, &y, &z);
          positions.emplace_back(x, y, z);
          /*std::size_t curr_ind{ 0 };
          std::string curr_sym;
          if (linestream >> curr_ind >> curr_sym >> x >> y >> z)
          {
            positions.emplace_back( x, y, z );
          }
          else
          {
            throw std::logic_error("Cannot obtain x,y,z for " + std::to_string(i) + ".");
          }*/
          if ((i - input_ensemble.size()*(N + 1u)) == N)
          { // if we are at the end of a structure 
            if (positions.size() != atoms.size())
              throw std::logic_error("The size of an additionally provided structure does not match the number of atoms.");
            input_ensemble.push_back(positions);
            positions.clear();
          }
        }
      }
    } // for

    // dividing subsystems
    if (!has_in_out_subsystems)
    {
      auto n_susy = Config::get().coords.subsystems.size();
      for (std::size_t i = 0; i < n_susy; ++i)
      {
        for (auto a : Config::get().coords.subsystems[i])
        {
          if (0 < a && a < (N - 1u))
          {
            atoms.atom(a - 1u).assign_to_system(i + 1u);
          }
        }
      }


    }

    if (indexation_not_contiguous)
    {
      std::cout << "Indexation not contiguous. Rebinding atoms." << lineend;
      for (std::size_t i(0U); i < atoms.size(); ++i)
      {
        for (auto bonding_partner : atoms.atom(i).bonds())
        {
          std::cout << "Partner of atom which is now " << i + 1 << " detached from " << bonding_partner << " and rebound to " << index_of_atom[bonding_partner] << lineend;
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
      coord_object.to_internal_light();
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

void coords::output::formats::tinker::to_stream(std::ostream & stream) const
{
  bool const pp = Config::get().energy.periodic && Config::get().energy.periodic_print;
  std::size_t const N(ref.size());
  stream << (pp ? N + 8 : N) << lineend;
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
  for (std::size_t i(0U); i < N; ++i)
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

void coords::output::formats::xyz::to_stream(std::ostream & stream) const
{
  std::size_t const N(ref.size());
  //stream << N << lineend;
  for (std::size_t i(0U); i < N; ++i)
  {
    /* stream << std::left  << std::setw(3) << atomic::symbolMap[ref.atoms(i).number()];
     stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).x();
     stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).y();
     stream << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << ref.xyz(i).z();
     stream << lineend;*/
    if (ref.atoms(i).fixed())
    {
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


void coords::output::formats::zmatrix::to_stream(std::ostream & stream) const
{
  std::size_t const N(ref.size());
  if (stream.good() && N > 0)
  {
    stream << "zmat angstroms" << lineend;
    stream << std::left << std::setw(10) << 1U;
    stream << std::left << std::setw(10) << ref.atoms(0u).i_to_a() + 1;
    stream << std::left << std::setw(4) << ref.atoms(ref.atoms(0u).i_to_a()).symbol() << lineend;
    if (N > 1)
    {
      stream << std::left << std::setw(10) << 2U;
      std::size_t const A = ref.atoms(1u).i_to_a();
      stream << std::left << std::setw(10) << A + 1;
      stream << std::left << std::setw(4) << ref.atoms(A).symbol();
      stream << std::right << std::setw(10) << ref.atoms(1u).ibond() + 1;
      stream << std::right << std::setw(10) << "bnd" << 1u;
      stream << lineend;
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
      stream << lineend;
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
      stream << lineend;
    }
    stream << "variables" << lineend;
    if (N > 1)
    {
      stream << "bnd" << std::left << std::setw(10) << 1u;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(1u).radius();
      stream << lineend;
    }
    if (N > 2)
    {
      stream << "bnd" << std::left << std::setw(10) << 2u;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(2u).radius();
      stream << lineend;
      stream << "ang" << std::left << std::setw(10) << 2u;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(2u).inclination();
      stream << lineend;
    }
    for (std::size_t i = 3; i < N; ++i)
    {
      stream << "bnd" << std::left << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(i).radius();
      stream << lineend;
      stream << "ang" << std::left << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(i).inclination();
      stream << lineend;
      stream << "dih" << std::left << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.intern(i).azimuth();
      stream << lineend;
    }
    stream << "constants" << lineend;
    for (std::size_t i = 1; i < N; ++i)
    {
      stream << "g_bnd" << std::left << std::setw(10) << i;
      stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.g_intern(i).x();
      stream << lineend;
      if (i > 1)
      {
        stream << "g_ang" << std::left << std::setw(10) << i;
        stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.g_intern(i).y();
        stream << lineend;
      }
      if (i > 2)
      {
        stream << "g_dih" << std::left << std::setw(10) << i;
        stream << std::right << std::fixed << std::setw(20) << std::setprecision(6) << ref.g_intern(i).z();
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
  for (std::size_t i(0U); i < N; ++i)
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
