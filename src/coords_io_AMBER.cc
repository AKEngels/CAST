/**
CAST 3
coords_io_AMEBR.cpp
Purpose: Reading from AMBER .prmtop and .rst (.crd) files

@author Dustin Kaiser, Susanne Sauer
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
  // in the amber99_gaff.prm file,
  // so everything is well
  // all atom types included that are in the amber_gaff.prm file 
  // however there might be alternative symbols for the atom types (like 2C instead of CT) that are not included yet
	size_t toTinkerType(std::string const & in)
	{
    if (in == "N  ") return    1u;
    else if (in == "CT " || in == "2C " || in == "3C " || in == "C8 ") return 2u;
    else if (in == "C  " || in == "CO ") return 3u;
    else if (in == "H  ") return    4u;
    else if (in == "O  ") return    5u;
    else if (in == "H1 ") return    6u;
    else if (in == "HC ") return   14u;
    else if (in == "OH ") return   63u;
    else if (in == "HO ") return   64u;
    else if (in == "SH ") return   85u;
    else if (in == "HS ") return   86u;
    else if (in == "S  ") return   95u;
    else if (in == "CA ") return  115u;
    else if (in == "HA ") return  117u;
    else if (in == "C* ") return  146u;
    else if (in == "CW ") return  147u;
    else if (in == "H4 ") return  148u;
    else if (in == "CB ") return  149u;
    else if (in == "NA ") return  150u;
    else if (in == "CN ") return  152u;
    else if (in == "CC ") return  169u;
    else if (in == "CR ") return  174u;
    else if (in == "H5 ") return  175u;
    else if (in == "CV ") return  189u;
    else if (in == "NB ") return  193u;
    else if (in == "O2 ") return  219u;
    else if (in == "HP ") return  284u;
    else if (in == "N3 ") return  285u;
    else if (in == "N2 ") return  299u;
    else if (in == "OS ") return 1001u;
    else if (in == "H2 ") return 1009u;
    else if (in == "N* ") return 1017u;
    else if (in == "CK ") return 1021u;
    else if (in == "NC ") return 1022u;
    else if (in == "CQ ") return 1023u;
    else if (in == "CM ") return 1082u;
    else if (in == "P  ") return 1230u;
    else if (in == "OW ") return 2001u;
    else if (in == "HW ") return 2002u;
    else if (in == "Li+") return 2003u;
    else if (in == "Na+") return 2004u;
    else if (in == "K+ ") return 2005u;
    else if (in == "Rb+") return 2006u;
    else if (in == "Cs+") return 2007u;
    else if (in == "Mg+") return 2008u;
    else if (in == "Ca+") return 2009u;
    else if (in == "Zn+") return 2010u;
    else if (in == "Ba+") return 2011u;
    else if (in == "Cl-") return 2012u;
    else if (in == "c  ") return 3000u;
    else if (in == "c1 ") return 3001u;
    else if (in == "c2 ") return 3002u;
    else if (in == "c3 ") return 3003u;
    else if (in == "ca ") return 3004u;
    else if (in == "cp ") return 3005u;
    else if (in == "cq ") return 3006u;
    else if (in == "cc ") return 3007u;
    else if (in == "cd ") return 3008u;
    else if (in == "ce ") return 3009u;
    else if (in == "cf ") return 3010u;
    else if (in == "cg ") return 3011u;
    else if (in == "ch ") return 3012u;
    else if (in == "cx ") return 3013u;
    else if (in == "cy ") return 3014u;
    else if (in == "cu ") return 3015u;
    else if (in == "cv ") return 3016u;
    else if (in == "cz ") return 3017u;
    else if (in == "h1 ") return 3018u;
    else if (in == "h2 ") return 3019u;
    else if (in == "h3 ") return 3020u;
    else if (in == "h4 ") return 3021u;
    else if (in == "h5 ") return 3022u;
    else if (in == "ha ") return 3023u;
    else if (in == "hc ") return 3024u;
    else if (in == "hn ") return 3025u;
    else if (in == "ho ") return 3026u;
    else if (in == "hp ") return 3027u;
    else if (in == "hs ") return 3028u;
    else if (in == "hw ") return 3029u;
    else if (in == "hx ") return 3030u;
    else if (in == "f  ") return 3031u;
    else if (in == "cl ") return 3032u;
    else if (in == "br ") return 3033u;
    else if (in == "i  ") return 3034u;
    else if (in == "n  ") return 3035u;
    else if (in == "n1 ") return 3036u;
    else if (in == "n2 ") return 3037u;
    else if (in == "n3 ") return 3038u;
    else if (in == "n4 ") return 3039u;
    else if (in == "na ") return 3040u;
    else if (in == "nb ") return 3041u;
    else if (in == "nc ") return 3042u;
    else if (in == "nd ") return 3043u;
    else if (in == "ne ") return 3044u;
    else if (in == "nf ") return 3045u;
    else if (in == "nh ") return 3046u;
    else if (in == "no ") return 3047u;
    else if (in == "o  ") return 3048u;
    else if (in == "oh ") return 3049u;
    else if (in == "os ") return 3050u;
    else if (in == "ow ") return 3051u;
    else if (in == "p2 ") return 3052u;
    else if (in == "p3 ") return 3053u;
    else if (in == "p4 ") return 3054u;
    else if (in == "p5 ") return 3055u;
    else if (in == "pb ") return 3056u;
    else if (in == "pc ") return 3057u;
    else if (in == "pd ") return 3058u;
    else if (in == "pe ") return 3059u;
    else if (in == "pf ") return 3060u;
    else if (in == "px ") return 3061u;
    else if (in == "py ") return 3062u;
    else if (in == "s  ") return 3063u;
    else if (in == "s2 ") return 3064u;
    else if (in == "s4 ") return 3065u;
    else if (in == "s6 ") return 3066u;
    else if (in == "sh ") return 3067u;
    else if (in == "ss ") return 3068u;
    else if (in == "sx ") return 3069u;
    else if (in == "sy ") return 3070u;
    else if (Config::get().general.energy_interface <= 3)
      std::cout << "Warning: Could not match AMBER atom type "
        << in << "to any TINKER atom type.\n\nThis will cause havok when force-field energy interfaces are specified (this is currently the case). "
        << "In case you did not expect this error, please expect everything to break very soon. Do not consider results after this message valid. "
        << "Please note that nevertheless you may use this AMBER input structure to perform calculations using CAST with energy-interfaces other than force-fields. "
        << "To do this, adjust your configuration file accordingly.\n"
        << std::endl;
    return 0u;
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
  std::vector<double> charges;
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
          if (line.substr(6, 9) == "SOLVENT_P") 
          {
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

      //Let's save the atom charges
      if (currentSectionID == 3u)
      {
        for (unsigned int countlines = 0u; countlines < std::ceil(float(numberOfAtoms) / 5.f) - 1; countlines++)
        {
          std::getline(config_file_stream, line);
          // A Line has 5 members in this format
          for (unsigned int i = 0u; i < 5u; i++)
          {
            std::string charge_str = line.substr(i * 16u, 16u);
            int exponent = std::stoi(charge_str.substr(13u, 3u));
            double koeff = std::stod(charge_str.substr(0u, 12u));
            charges.push_back((koeff * pow(10, exponent)));
          }
        }
        // read last line if existent
        std::getline(config_file_stream, line);
        if (line.substr(0, 1) != "%")
        {
          for (unsigned int i = 0; i < line.size(); i = i + 16)
          {
            std::string charge_str = line.substr(i * 1u, 16);
            int exponent = std::stoi(charge_str.substr(13u, 3u));
            double koeff = std::stod(charge_str.substr(0u, 12u));
            charges.push_back((koeff * pow(10, exponent)));
          }
        }
        Config::set().coords.amber_charges = charges;
      }

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

      //If section == AMBER_ATOM_TYPE
      if (currentSectionID == 34u)
      {
        // Discard the format specifier
        // Read all lines except last one
        for (unsigned int countlines = 0u; countlines < std::ceil((float(numberOfAtoms) + 0.5) / 20.f) - 1; countlines++)
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

		      	if (symbol == "CU ") symbol = "Cu";
            else if (symbol == "IM ") symbol = "CL";
            else if (symbol == "QC ") symbol = "Cs";
            else if (symbol == "QK ") symbol = "K";
            else if (symbol == "QL ") symbol = "Li";
            else if (symbol == "QN ") symbol = "NA";
            else if (symbol == "QR ") symbol = "NA";
            else if (symbol == "Na+") symbol = "NA";
            else if (symbol == "IP ")
            {
              symbol = "NA";
              std::cout << "Atom Type \"IP\" used. We are guessing it stands for sodium.\n This is not documented in the official "
                << "AMBER Atom Types webpage (\"http://ambermd.org/antechamber/gaff.html\" or "
                << "\"http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html\").\n"
                << "So you need to check if it really IS meant to be sodium.\n"
                << "Please talk to a CAST programmer about this.\n";
            }

            //wtf is lone pair? certainly not an atom!
            else if (symbol == "LP ") continue;
            else if (symbol == "nu ") continue;

            //Si? Starts with an s like sulfur, so we have to check..gnaah.....
            else if (symbol == "si " || symbol == "Si " || symbol == "SI ") { symbol = "Si"; }

		      	//Cl? Starts with a c like carbon, so we have to check..gnaah.....
		      	else if (symbol == "cl " || symbol == "Cl " || symbol == "CL ") { symbol = "CL"; }

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
          
		      if (symbol == "CU ") symbol = "Cu";
          else if (symbol == "IM ") symbol = "CL";
          else if (symbol == "QC ") symbol = "Cs";
          else if (symbol == "QK ") symbol = "K";
          else if (symbol == "QL ") symbol = "Li";
          else if (symbol == "QN ") symbol = "NA";
          else if (symbol == "QR ") symbol = "NA";
          else if (symbol == "Na+") symbol = "NA";
          else if (symbol == "IP ")
          {
            symbol = "NA";
            std::cout << "Atom Type \"IP\" used. We are guessing it stands for sodium.\n This is not documented in the official "
              << "AMBER Atom Types webpage (\"http://ambermd.org/antechamber/gaff.html\" or "
              << "\"http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html\").\n"
              << "So you need to check if it really IS meant to be sodium.\n"
              << "Please talk to a CAST programmer about this.\n";
          }

          //wtf is lone pair? certainly not an atom!
          else if (symbol == "LP ") continue; 
          else if (symbol == "nu ") continue;

          //Si? Starts with an s like sulfur, so we have to check..gnaah.....
          else if (symbol == "si " || symbol == "Si" || symbol == "SI") { symbol = "Si"; }

		      //Cl? Starts with a c like carbon, so we have to check..gnaah.....
		      else if (symbol == "cl " || symbol == "Cl " || symbol == "CL ") { symbol = "CL"; }

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
      if (atoms.atom(bondsWithHydrogen[i] - 1u).energy_type() == 2002 && atoms.atom(bondsWithHydrogen[i + 1u] - 1u).energy_type() == 2002)
      {
        // do not connect the two hydrogens of a water molecule, even it there is a bond in the prmtop file
      }
      else
      {
        atoms.atom(bondsWithHydrogen[i] - 1u).bind_to(bondsWithHydrogen[i + 1u] - 1u);
        atoms.atom(bondsWithHydrogen[i + 1u] - 1u).bind_to(bondsWithHydrogen[i] - 1u);
      } 
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
      std::ifstream coord_file_stream(Config::get().io.amber_restrt.c_str(), std::ios_base::in);
      if (! coord_file_stream.good()) 
        throw std::runtime_error("Can't read AMBER restart file " 
          + Config::get().io.amber_restrt 
          + " which was specified in the inputfile.\n");
      // Now, discard the title and the line with the atom number (and simulation time)
      std::getline(coord_file_stream, line);
      std::getline(coord_file_stream, line);

      while (std::getline(coord_file_stream, line))
      {
        std::string x = line.substr(0, 12);
        std::string y = line.substr(12, 12);
        std::string z = line.substr(24, 12);
        position.x() = std::stod(x);
        position.y() = std::stod(y);
        position.z() = std::stod(z);
        positions.push_back(position);

        if (line.size() > 36)   // normally two coordinates in one line
        {
          x = line.substr(36, 12);
          y = line.substr(48, 12);
          z = line.substr(60, 12);
          position.x() = std::stod(x);
          position.y() = std::stod(y);
          position.z() = std::stod(z);
          positions.push_back(position);
        }

      }

      // delete the box parameters in the last line if there are some (size and angles)
      if (positions.size() == atoms.size() + 2)
      {
        positions.pop_back();
        positions.pop_back();
      }
      // delete the box parameters in the last line if there are some (only box size)
      if (positions.size() == atoms.size() + 1)
      {
        positions.pop_back();
      }

      input_ensemble.push_back(positions);
      goto DONE;
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

	if (Config::get().general.chargefile)   // read charges from chargefile
	{
		charges.clear();  // delete former charges
		std::ifstream charge_stream("charges.txt", std::ios_base::in);
		std::string read;
		while (charge_stream)
		{
			if (charge_stream >> read)
			{
				// ignore atom number
			}
			if (charge_stream >> read)
			{
				// ignore atom type
			}
			if (charge_stream >> read)
			{
				charges.push_back(std::stod(read));
			}
		}
		if (charges.size() == coord_object.size())
		{
			Config::set().coords.amber_charges = charges;
			if (Config::get().general.verbosity > 3)
			{
				std::cout << "Reading charges from chargefile successful.\n";
			}
		}
		else   // not the correct number of charges in chargefile
		{
			std::cout << "Reading charges from chargefile failed. Using charges from prmtopfile instead.\n";
		}
	}

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