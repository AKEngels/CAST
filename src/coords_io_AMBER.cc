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
  // in the amber99.prm file,
  // so everything is well
	//NOT COMPLETE
	size_t toTinkerType(std::string const & in)
	{
		if (in == "CT ") return 2u;          //tetrahedral C
		else if (in == "CX ") return 2u;     //tetrahedral C
		else if (in == "C3 ") return 2u;     //tetrahedral C
		else if (in == "2C ") return 2u;     //tetrahedral C
		else if (in == "3C ") return 2u;     //tetrahedral C
		else if (in == "C8 ") return 2u;     //tetrahedral C
		else if (in == "CL-") return 2011u; //Cloride ion
		else if (in == "IM ") return 2011u;  //Cloride ion
		else if (in == "C  ") return 3u;      //sp2 C (carbonyl)
		else if (in == "CO ") return 3u;     //sp2 C (carbonyl)
		else if (in == "CA ") return 115u;   //sp2 C (aromatic)
		else if (in == "CM ") return 1082u;  //sp2 C (pyridine)
		else if (in == "CC ") return 169u;   //sp2 C (aromatic)
		else if (in == "CV ") return 189u;   //sp2 C (aromatic)
		else if (in == "CW ") return 204u;   //sp2 C (aromatic)
		else if (in == "CR ") return 206u;   //sp2 C (aromatic)
		else if (in == "CB ") return 1018u;  //sp2 C (aromatic)
		else if (in == "C* ") return 146u;   //sp2 C (aromatic)
		else if (in == "CN ") return 152u;   //sp2 C (aromatic)
		else if (in == "CK ") return 1021u;  //sp2 C (aromatic)
		else if (in == "CQ ") return 1023u;  //sp2 C (aromatic)
		else if (in == "N  ")  return 7u;     //amide N
		else if (in == "NA ") return 150u;   //N in ring, bound to H
		else if (in == "Cs+") return 2007u; //Cs ion
		else if (in == "QC ") return 2007u;  //Cs ion
		else if (in == "K+ ") return 2005u;  //K ion
		else if (in == "QK ") return 2005u;  //K ion
		else if (in == "Li+") return 2003u; //Li ion
		else if (in == "QL ") return 2003u;  //Li ion
		else if (in == "Na+") return 2004u; //Na ion
		else if (in == "QN ") return 2004u;  //Na ion
		else if (in == "QR ") return 2006u;  //Rb ion
		else if (in == "IP ") return 2004u;  //Na ion
		else if (in == "NB ") return 193u;   //N in ring, with LP
		else if (in == "NC ") return 1022u;  //N in ring, with LP
		else if (in == "N* ") return 1047u;  //N in ring, bound to alkyle group
		else if (in == "N2 ") return 1057u;  //basic NH2 group
		else if (in == "N3 ") return 285u;   //sp3 N
		else if (in == "OW ") return 2001u;  //water O
		else if (in == "OH ") return 63u;    //alcohole O
		else if (in == "OS ") return 1001u;  //O in ether or ester
		else if (in == "O  ")  return 1060u;  //carbonyle O
		else if (in == "O2 ") return 219u;   //carboyle or phosphate (non-bounded) O
		else if (in == "S  ")  return 95u;    //S without H
		else if (in == "SH ") return 85u;    //S with H
		else if (in == "P  ")  return 1230u;  //P in phosphate
		else if (in == "H  ")  return 90u;    //amide or imino H
		else if (in == "HW ") return 2002u;  //water H
		else if (in == "HO ") return 64u;    //alcohole H
		else if (in == "HS ") return 86u;    //H in SH
		else if (in == "H5 ") return 1026u;  //???
		else if (in == "HA ") return 117u;   //H bound to aromatic C
		else if (in == "HC ") return 129u;   //H bound to aliphatic C
		else if (in == "H1 ") return 143u;   //H bound to aliphatic C
		else if (in == "H2 ") return 1009u;  //H in NH2
		else if (in == "H3 ") return 1113u;  //H in NH3+
		else if (in == "HP ") return 284u;   //H bound to P???
		else if (in == "H4 ") return 148u;   //??? 
		else if (in == "c3 ") return 3003u;  // GAFF  -> charge from atom type no. 2 (CT)
		else if (in == "ca ") return 3004u;  // GAFF  -> charge from atom type no. 115 (C) 
		else if (in == "n  ") return 3035u;  // GAFF  -> charge = example from prmtop
		else if (in == "c  ") return 3000u;  // GAFF  -> charge = example from prmtop
		else if (in == "oh ") return 3049u;  // GAFF  -> charge = example from prmtop
		else if (in == "ho ") return 3026u;  // GAFF  -> charge = example from prmtop
		else if (in == "o  ") return 3048u;  // GAFF  -> charge = example from prmtop
		else if (in == "hc ") return 3024u;  // GAFF  -> charge from atom type no. 129 (HC) 
		else if (in == "ha ") return 3023u;  // GAFF  -> charge = example from prmtop
		else if (in == "h1 ") return 3018u;  // GAFF  -> charge = example from prmtop
		else if (in == "cl ") return 3032u;  // GAFF  -> charge = example from prmtop
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
  std::cout << "NOTE: You are specifiying AMBER-Files as input.\n";
  std::cout << "All fancy AMBER-stuff like residues and such is OMITTED. CAST only collects the atoms and their positions.\n";
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
  
  // -> Currently not used, uncomment if you want to use it...
  // bool skipPastAtomName = std::find(sectionsPresent.begin(), sectionsPresent.end(), 34u) != sectionsPresent.end();

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
		  // read last line of bondsWithHydrogen if existent
		  std::getline(config_file_stream, line);
		  if (line.substr(0, 1) != "%")
		  {
			  for (unsigned int i = 0; i < line.size(); i = i + 8)
			  {
				  if (state % 3u != 0)
				  {
					  bondsWithHydrogen.push_back(std::stoi(line.substr(i * 1u, 8)) / 3u + 1u);
				  }
				  state = state + 1;
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
		// read last line of bondsWithoutHydrogen if existent
		std::getline(config_file_stream, line);
		if (line.substr(0, 1) != "%")
		{
			for (unsigned int i = 0; i < line.size(); i = i + 8)
			{
				if (state % 3u != 0)
				{
					bondsWithoutHydrogen.push_back(std::stoi(line.substr(i * 1u, 8)) / 3u + 1u);
				}
				state = state + 1;
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
            std::string symbol = line.substr(i * 4u, 3u);
			const std::string amberAtomType(symbol);

            // First, let's process all those AMBER Atom Types with
            // just incredibly retarded names. damn.
            // via: http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html

			if (symbol == "CU") symbol = "Cu";
            else if (symbol == "IM") symbol = "CL";
            else if (symbol == "QC") symbol = "Cs";
            else if (symbol == "QK") symbol = "K";
            else if (symbol == "QL") symbol = "Li";
            else if (symbol == "QN") symbol = "NA";
            else if (symbol == "QR") symbol = "NA";
            else if (symbol == "Na") symbol = "NA";
            else if (symbol == "IP")
            {
              symbol = "NA";
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
                std::cout << "AMBER atom type " << amberAtomType << " could not be matched to TINKER atom class. All atom types are therefore omitted, do not use force field energy interfaces." << std::endl;
              }
            }
            atoms.add(current);
          }
        }
        //Process last line
        std::getline(config_file_stream, line);
        for (unsigned int i = 0u; i < numberOfAtoms % 20; i++)
        {
          std::string symbol = line.substr(i * 4u, 3u);
		  const std::string amberAtomType(symbol);

          // First, let's process all those AMBER Atom Types with
          // just incredibly retarded names. damn.
          // via: http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html
          
		  if (symbol == "CU") symbol = "Cu";
          else if (symbol == "IM") symbol = "CL";
          else if (symbol == "QC") symbol = "Cs";
          else if (symbol == "QK") symbol = "K";
          else if (symbol == "QL") symbol = "Li";
          else if (symbol == "QN") symbol = "NA";
          else if (symbol == "QR") symbol = "NA";
          else if (symbol == "Na") symbol = "NA";
          else if (symbol == "IP")
          {
            symbol = "NA";
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
		  if (!ignoreFFtypes)
		  {
			  size_t type = amberUtil::toTinkerType(amberAtomType);
			  if (type != 0u) current.set_energy_type(type);
			  else
			  {
				  ignoreFFtypes = true;
				  std::cout << "AMBER atom type " << amberAtomType << " could not be matched to TINKER atom class. All atom types are therefore omitted, do not use force field energy interfaces." << std::endl;
			  }
		  }
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

      // Now, discard the title (if existent) and the line with the atom number
      std::getline(coord_file_stream, line);
	  int atom_number;
	  bool title_exists = false;
      try
	  {
		  atom_number = std::stoi(line);
	  }
	  catch(...)
	  {
		  title_exists = true;
	  }
	  if (title_exists == true)
	  {
		  std::getline(coord_file_stream, line);
	  }

      unsigned long state = 0u; //Counts each processed floating point number.
      while (std::getline(coord_file_stream, line))
      {
		  if (line.size() > 36)
		  {
			  double tempx(0.0), tempy(0.0), tempz(0.0), tempx_2(0.0), tempy_2(0.0), tempz_2(0.0);
			  sscanf_s(line.c_str(), "%lf %lf %lf %lf %lf %lf", &tempx, &tempy, &tempz, &tempx_2, &tempy_2, &tempz_2);
			  position.x() = tempx;
			  position.y() = tempy;
			  position.z() = tempz;
			  positions.push_back(position);
			  position.x() = tempx_2;
			  position.y() = tempy_2;
			  position.z() = tempz_2;
			  positions.push_back(position);
		  }
		  else   // coordinate line only half
		  {
			  double tempx(0.0), tempy(0.0), tempz(0.0);
			  sscanf_s(line.c_str(), "%lf %lf %lf", &tempx, &tempy, &tempz);
			  position.x() = tempx;
			  position.y() = tempy;
			  position.z() = tempz;
			  positions.push_back(position);
		  }
	  }
	  // delete the box parameters in the last line if there are some
	  if (positions.size() == atoms.size() + 2)
	  {
		  positions.pop_back();
		  positions.pop_back();
	  }
	  
	  input_ensemble.push_back(positions);
	  goto DONE;
	}

		  
        // A line has 10 members in the format
   //     for (unsigned int i = 0; i < 10u; state++, i++)
   //     {
   //       // Needed: Check if substr is empty
   //       if (line.substr(i * 8u, 8u).empty())
   //       {
   //         if (!(state % 3u == 0))
   //         {
   //           std::cout << "Encounterd unexpected linebreak or EOF in reading AMBER mdcrd file.\n";
   //           goto FAILED;
   //         }
   //         else
   //         {
   //           ///////////////////////////////
   //           // EXIT WHILE AND DO LOOP!!! //    <- THIS HERE IS IMPORTANT!!!!!!!
   //           ///////////////////////////////
   //           goto DONE;
   //         }
   //       }
   //       if (state % 3u == 0)
   //       {
   //         position.x() = std::stod(line.substr(i * 12u, 12u));
   //       }
   //       else if (state % 3u == 1)
   //       {
   //         position.y() = std::stod(line.substr(i * 12u, 12u));
   //       }
   //       else if (state % 3u == 2)
   //       {
   //         position.z() = std::stod(line.substr(i * 12u, 12u));
   //         positions.push_back(position);
			//std::cout << positions << "\n";
   //         //Check if we reached end of structure, conditional is true if yes
   //         if (((state + 1) / 3) % numberOfAtoms == 0)
   //         {
   //           input_ensemble.push_back(positions);
   //           if (positions.size() != atoms.size())
   //           {
   //             throw std::logic_error("The size of an provided structure does not match the number of atoms.");
   //           }
   //           positions.clear();
   //           state++;

   //           //In mdcrd we need to check if a linebreak occurs and box coordiantes are written.
   //           if (Config::get().io.amber_trajectory_at_constant_pressure)
   //           {
   //             // Discard line containing box size
   //             std::getline(coord_file_stream, line);
   //           }
   //           break;
   //         }
   //       }
   //     }

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