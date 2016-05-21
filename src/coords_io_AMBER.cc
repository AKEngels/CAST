/**
CAST 3
coords_io_AMEBR.cpp
Purpose: Reading from AMBER .prmtop and .crd files

@author Dustin Kaiser
@version 1.0
*/

#include "coords_io.h"

namespace amberUtil
{
  // because of the currently very entangled
  // processes in CAST, we simply assign (wrong)
  // TINKER atom types to the atoms according to
  // their AMBER atom types. These fake TINKER types
  // correspond however to the correct TINKER classes
  // in the amber.prm file,
  // so everything is well
  size_t toTinkerType(std::string const& in)
  {
    if (in == "CT") return 2u;
    else if (in == "C") return 3u;
    else if (in == "CA") return 115u;
    else if (in == "CM") return 1082u;
    else if (in == "CC") return 169u;
    else if (in == "CV") return 189u;
    else if (in == "CW") return 204u;
    else if (in == "CR") return 206u;
    else if (in == "CB") return 1018u;
    else if (in == "C*") return 146u;
    else if (in == "CN") return 152u;
    else if (in == "CK") return 1021u;
    else if (in == "CQ") return 1023u;
    else if (in == "N")  return 7u;
    else if (in == "NA") return 150u;
    else if (in == "NB") return 193u;
    else if (in == "NC") return 1022u;
    else if (in == "N*") return 1047u;
    else if (in == "N2") return 1057u;
    else if (in == "N3") return 285u;
    else if (in == "OW") return 2001u;
    else if (in == "OH") return 145u;
    else if (in == "OS") return 1001u;
    else if (in == "O")  return 1060u;
    else if (in == "S")  return 95u;
    else if (in == "SH") return 85u;
    else if (in == "P")  return 1230u;
    else if (in == "H ")  return 90u;
    else if (in == "HW") return 2002u;
    else if (in == "HO") return 774u;
    else if (in == "HS") return 86u;
    else if (in == "HA") return 117u;
    else if (in == "HC") return 129u;
    else if (in == "H1") return 143u;
    else if (in == "H2") return 1009u;
    else if (in == "H3") return 1113u;
    else if (in == "HP") return 284u;
    else if (in == "H4") return 148u;
    else if (in == "H5") return 86u;
    else return 0u;
  }
}


// Amber is aimed at FORTRAN parsers so this code is just, horrible... just horrible...
// See http://ambermd.org/prmtop.pdf
// Since AMBER Stores actual coordinates in a seperate file, this procedure will
// read from the config-namespace. We search for the name of a possibly existing second file
// containing the raw coordinates. Pay attantion to this!
//////////////////////////////
/////
///// This whole section needs to be rewritten cleanly ASAP
///// I dont know if I was like drunk when I wrote this, but, omg, convoluted as sh*t.
///// In case anyone has to dig through this, I am terribly sorry.
/////
//////////////////////////////
coords::Coordinates coords::input::formats::amber::read(std::string file)
{
  Coordinates coord_object;
  std::ifstream config_file_stream(file.c_str(), std::ios_base::in);
  std::string line;
  std::cout << "NOTE: You are specifiying AMBER-Files as input.\nThis feature is highly, let me repeat that, HIGHLY experimental in CAST.\n";
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

    bool ignoreFFtypes = false;

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
          //Special case, we should get rid of this really fcking soon omg -.-'
          //It is by design because of the find function
          if (line.substr(6, 9) == "SOLVENT_P") {
            currentSectionID = 666; break; // this... is... just... horrible
          }

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
            const std::string amberAtomType(symbol);

            // First, let's process all those AMBER Atom Types with
            // just incredibly retarded names. damn.
            // via: http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html
            if (symbol == "CU")
            {
              symbol = "Cu";
            }
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

            //Weird stuff

            //3C 2C 1C
            else if (symbol.substr(1, 1) == "c" || symbol.substr(1, 1) == "C") symbol = "C";

            //Maybe we are already done?
            //FAILED
            else if (!symbol.empty())
            {
              std::cout << "Could not identify AMBER_ATOM_TYPE: " << line.substr(i * 4u, 4u) << "\nstopping!\n";
              std::cout << "You should probably go talk to a CAST programmer about this.\n";
              goto FAILED;
            }

            Atom current(symbol);
            if (!ignoreFFtypes)
            {
              size_t type = amberUtil::toTinkerType(amberAtomType);
              if (type != 0u) current.set_energy_type(type);
              else
              {
                ignoreFFtypes = true;
                std::cout << "AMBER atom type " << amberAtomType << " could not be matched to TINKER atom class. Atom types omitted, do not use fore field energy interfaces." << std::endl;
              }
            }
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
            std::cout << "Could not identify AMBER_ATOM_TYPE: " << line.substr(i * 4u, 4u) << "\nstopping!\n";
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
      atoms.atom(bondsWithHydrogen[i] - 1u).bind_to(bondsWithHydrogen[i + 1u] - 1u);
      atoms.atom(bondsWithHydrogen[i + 1u] - 1u).bind_to(bondsWithHydrogen[i] - 1u);

    }
    for (unsigned int i = 0u; i < bondsWithoutHydrogen.size(); i = i + 2u)
    {
      atoms.atom(bondsWithoutHydrogen[i] - 1u).bind_to(bondsWithoutHydrogen[i + 1u] - 1u);
      atoms.atom(bondsWithoutHydrogen[i + 1u] - 1u).bind_to(bondsWithoutHydrogen[i] - 1u);
    }

    // OK, let's now fetch the actual coordinates
    Representation_3D positions;
    positions.reserve(numberOfAtoms); //Reserve space
    if (!Config::get().io.amber_mdcrd.empty())
    {
      std::ifstream coord_file_stream(Config::get().io.amber_mdcrd.c_str(), std::ios_base::in);
      // Now, discard the title
      std::getline(coord_file_stream, line);

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
            //Check if we reached end of structure, conditional is true if yes
            if (((state + 1) / 3) % numberOfAtoms == 0)
            {
              input_ensemble.push_back(positions);
              if (positions.size() != atoms.size())
              {
                throw std::logic_error("The size of an provided structure does not match the number of atoms.");
              }
              positions.clear();
              state++;

              //In mdcrd we need to check if a linebreak occurs and box coordiantes are written.
              if (Config::get().io.amber_trajectory_at_constant_pressure)
              {
                // Discard line containing box size
                std::getline(coord_file_stream, line);
              }
              break;
            }
          }
        }
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

    for (auto & p : input_ensemble)
    {
      p.gradient.cartesian.resize(p.structure.cartesian.size());
      coord_object.set_xyz(p.structure.cartesian);
      coord_object.to_internal_light();
      p = coord_object.pes();
    }
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