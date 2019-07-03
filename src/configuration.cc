#include "configuration.h"
#include "helperfunctions.h"

#include "helperfunctions.h"

/**
 * Global static instance of the config-object.
 * There can only ever be one config-object.
 */
Config * Config::m_instance = nullptr;

/**
* Helper function that sorts numerical
* integer type numbers into a vector. Every number
* is only inserted once ("uniquely").
*
* @param str: String containing number range (ex. "0-2, 7, 2, 3, 3, 13-15")
* @return: vector containing sorted unique numbers (ex. "0, 1, 2, 3, 7, 13, 14, 15")
*/
std::vector<std::size_t> config::sorted_indices_from_cs_string(std::string str, bool minus_1)
{
  // remove all spaces
  str.erase(std::remove_if(str.begin(), str.end(),
    [](char c) -> bool {return std::isspace(c) > 0; }),
    str.end());
  // replace all commas with spaces
  std::replace(str.begin(), str.end(), ',', ' ');
  std::string d;
  std::stringstream iss(str);
  std::vector<std::size_t> re;
  // get each seperated value as single string d
  while (iss >> d)
  {
    auto dash_pos = d.find('-');
    // if element is a range (contains '-')
    if (dash_pos != std::string::npos)
    {
      d[dash_pos] = ' ';
      std::stringstream pss{ d };
      std::size_t first(0), last(0);
      if (pss >> first && pss >> last)
      {
        if (first <= last && (!minus_1 || first > 0))
        {
          for (auto i = first; i <= last; ++i)
          {
            re.push_back(minus_1 ? i - 1u : i);
          }
        }
        else
        {
          throw std::runtime_error("Invalid range for indices: '" +
            std::to_string(first) + " - " + std::to_string(last) + "'.");
        }
      }

      else
      {
        throw std::runtime_error("Cannot read index range from '" + d + "'.");
      }
    }
    // throw if non-numeric character is found
    else if (d.find_first_not_of("0123456789") != std::string::npos)
    {
      throw std::runtime_error("Non numeric character found in '" + d + "'.");
    }
    // read number from stringstream of d
    else
    {
      std::stringstream pss{ d };
      std::size_t value;
      if (pss >> value && (!minus_1 || value > 0))
      {
        re.push_back(minus_1 ? value- 1u : value);
      }
      else
      {
        throw std::runtime_error("Cannot read index from '" + d + "'.");
      }
    }
  }
  // sort resulting numbers
  std::sort(re.begin(), re.end());
  // remove duplicates and return rest
  return std::vector<std::size_t>{re.begin(), std::unique(re.begin(), re.end())};
}

std::vector<double> config::doubles_from_string(std::string str)
{
  std::vector<double> result;
  std::vector<std::string> stringvec = split(str, ',');
  for (auto i : stringvec)
  {
		if (check_if_number(i) == true) result.emplace_back(std::stod(i));
		else throw std::runtime_error(i + " can't be converted to double.");
  }
  return result;
}


template<typename T>
static T clip(T value, T const LOW, T const HIGH)
{
  using std::min;
  using std::max;
  return min(HIGH, max(LOW, value));
}


template<typename ENUM_T, std::size_t ARR_SIZE>
static ENUM_T enum_type_from_string_arr(std::string const & S, std::string const (&arr)[ARR_SIZE])
{
  for (std::size_t i = 0; i < ARR_SIZE; ++i)
  {
    if (S.find(arr[i]) != S.npos) return static_cast<ENUM_T>(i);
  }
  return static_cast<ENUM_T>(-1);
}


config::tasks::T Config::getTask(std::string const & S)
{
  for (std::size_t i = 0; i < config::NUM_TASKS; ++i)
  {
    if (S.find(config::task_strings[i]) != S.npos)
      return static_cast<config::tasks::T>(i);
  }
  return config::tasks::ILLEGAL;
}

config::interface_types::T Config::getInterface(std::string const & S)
{
  for (std::size_t i = 0; i < config::NUM_INTERFACES; ++i)
  {
    if (S.find(config::interface_strings[i]) != S.npos)
      return static_cast<config::interface_types::T>(i);
  }
  return config::interface_types::ILLEGAL;
}

config::output_types::T Config::getOutFormat(std::string const & S)
{
  for (std::size_t i = 0; i < config::NUM_INTERFACES; ++i)
  {
    if (S.find(config::output_strings[i]) != S.npos)
      return static_cast<config::output_types::T>(i);
  }
  return config::output_types::ILLEGAL;
}


/*


  ########     ###    ########   ######  ########
  ##     ##   ## ##   ##     ## ##    ## ##
  ##     ##  ##   ##  ##     ## ##       ##
  ########  ##     ## ########   ######  ######
  ##        ######### ##   ##         ## ##
  ##        ##     ## ##    ##  ##    ## ##
  ##        ##     ## ##     ##  ######  ########


   ######   #######  ##     ## ##       #### ##    ## ########
  ##    ## ##     ## ###   ### ##        ##  ###   ## ##
  ##       ##     ## #### #### ##        ##  ####  ## ##
  ##       ##     ## ## ### ## ##        ##  ## ## ## ######
  ##       ##     ## ##     ## ##        ##  ##  #### ##
  ##    ## ##     ## ##     ## ##        ##  ##   ### ##
   ######   #######  ##     ## ######## #### ##    ## ########


*/

std::string config::config_file_from_commandline(std::ptrdiff_t const N, char **V)
{
  if (N > 1)
  {
    for (std::ptrdiff_t i(1); i<N; ++i)
    {
      std::string argument(V[i]);
      if (argument.substr(0U, 7U) == "-setup=")
      {
        std::string ffile(argument.substr(7U, argument.length()));
        if (file_exists_readable(ffile)) return ffile;
      }
      else if (argument.substr(0U, 2U) == "-s" && argument.length() == 2U && N > i + 1)
      {
        std::string nextArgument(V[i + 1]);
        if (file_exists_readable(nextArgument)) return nextArgument;
      }
    }
  }
  if (file_exists_readable("CAST.txt")) return std::string("CAST.txt");
  return std::string("INPUTFILE");
}

void config::parse_command_switches(std::ptrdiff_t const N, char **V)
{
  if (N > 1)
  {
    for (std::ptrdiff_t i(1U); i<N; ++i)
    {
      // argument to string
      std::string argument(V[i]);
      std::string::size_type const equality_pos(argument.find("="));

      // This section transfers the commandline switches
      // into a form which can be parsed by the function
      // "parse option", see also parse_option(std::string option, std::string value)
      //
      // This only works for switches of the form -option=value
      // Example, -cutoff=5.0
      if (argument.substr(0U, 1U) == "-" && equality_pos != argument.npos &&
        argument.length() >(equality_pos + 1U))
      {
        std::string option(argument.substr(1U, (equality_pos - 1U)));
        std::string value(argument.substr((equality_pos + 1U), (argument.length() - equality_pos - 1U)));

        // if the value is passed to this function containing quotation marks,
        // they will be removed for your conveniance
        if (value.substr(0U, 1U) == "\"")
        {
          std::string::size_type const last_quotationmark(argument.find_last_of("\""));
          if (last_quotationmark == 0U)
          {
            bool found_end(false);
            while (!found_end && (i<(N - 1)))
            {
              value.append(" ");
              std::string this_arg(V[++i]);
              std::string::size_type quot_pos(this_arg.find("\""));
              if (quot_pos != this_arg.npos)
              {
                value.append(this_arg.substr(0U, quot_pos));
                found_end = true;
              }
              else value.append(this_arg);
            }
          }
          else value = value.substr(1U, last_quotationmark - 1U);

        }
        config::parse_option(option, value);
      }

      // The following switches can also take the form:
      // "-n yourChosenName" or "-p paramfile.prm"
      //
      // This only works for name, paramfile, outname, task and inputtype
      if (argument.substr(0, 2) == "-n" && argument.length() == 2U && N > i + 1)
      {
        config::parse_option("name", V[++i]);
      }
      else if (argument.substr(0, 2) == "-p" && argument.length() == 2U && N > i + 1)
      {
        config::parse_option("paramfile", V[++i]);
      }
      else if (argument.substr(0, 2) == "-o" && argument.length() == 2U && N > i + 1)
      {
        config::parse_option("outname", V[++i]);
      }
      else if (argument.substr(0, 2) == "-t" && argument.length() == 2U && N > i + 1)
      {
        config::parse_option("task", V[++i]);
      }
      else if (argument.substr(0, 5) == "-type" && argument.length() == 5U && N > i + 1)
      {
        config::parse_option("inputtype", V[++i]);
      }

      // Enable spackman by entering -spack
      else if (argument.substr(0, 6) == "-spack")
      {
        Config::set().energy.spackman.on = true;
      }
    }
  }

  std::string::size_type f = Config::set().general.outputFilename.find("%i");
  if (f != std::string::npos)
  {
    Config::set().general.outputFilename =
      Config::set().general.outputFilename.substr(0, f) +
      scon::StringFilePath(Config::get().general.inputFilename).base_no_extension() +
      Config::set().general.outputFilename.substr(f + 2);
  }
}


/*


########     ###    ########   ######  ########
##     ##   ## ##   ##     ## ##    ## ##
##     ##  ##   ##  ##     ## ##       ##
########  ##     ## ########   ######  ######
##        ######### ##   ##         ## ##
##        ##     ## ##    ##  ##    ## ##
##        ##     ## ##     ##  ######  ########


*/

/*! Helperfunction to keep track of the output-filename
 *
 * Small helper function that keeps
 * track of wether an "outname" has
 * been specified and read up until now.
 *
 * Reasoning: If no "outname" is specified in the
 * configfile, "inputname" + "_out" will be taken
 * as the outname. This function is used only in
 * config::parse_option
 */
static bool outname_check(int x = 0)
{
  static int chk = 0;
  chk += x;
  return chk < 1;
}

/**
* function to find all appearing or disappearing atoms in an FEP input file
* returns a vector with tinker atom numbers
*/
std::vector<unsigned> FEP_get_inout()
{
	std::vector<unsigned> temp;
	std::ifstream config_file_stream(Config::set().general.inputFilename.c_str(), std::ios_base::in);
	std::string line;
	while (std::getline(config_file_stream, line))
	{
		std::vector<std::string> fields;
		std::istringstream iss(line);
		std::string sub;
		while (iss >> sub)
		{
			fields.push_back(sub);
		}
		if (fields[fields.size()-1] == "IN" || fields[fields.size() - 1] == "OUT")
		{
			std::string atom_number_str = fields[0];
			unsigned atom_number = std::stoi(atom_number_str);
			temp.push_back(atom_number);
		}
	}
	return temp;
}

void config::parse_option(std::string const option, std::string const value_string)
{
  std::istringstream cv(value_string);

  /////////////////////
  //// Config::general
  ////////////////////
  if (option == "name")
  {
	  Config::set().general.inputFilename = value_string;

	  // If no outname is specified,
	  // the output-file will have the same name as the inputfile
	  // but with "_out" added to the name.
	  //
	  // outname_check is a small function used to test if
	  // an outname has already been specified by keeping a
	  // static counter.
	  if (outname_check(0))
	  {
		  std::string path = value_string;
		  // Remove quote signs
		  if (path.front() == '"' && path.back() == '"')
		  {
			  path = path.substr(1u, path.length() - 2U);
		  }
		  scon::FilePath<std::string> in_file_path(path);
		  Config::set().general.outputFilename = in_file_path.base_no_extension() + "_out";
	  }
  }
  // Name of the outputfile
  else if (option == "outname")
  {
	  // outname_check is a small function used to test if
	  // an outname has already been specified by keeping a
	  // static counter.
	  outname_check(1);

	  // If the "outname" starts with an +,
	  // we append the (input)name by this string
	  if (value_string[0] == '+')
	  {
		  Config::set().general.outputFilename += value_string.substr(1);
	  }
	  else
		  Config::set().general.outputFilename = value_string;
  }

  // Filename of the ForceField parameter-file
  // Has to be in same folder as executable
  // Default: oplsaa.prm
  else if (option == "paramfile")
    Config::set().general.paramFilename = value_string;

  // Option to read charges from seperate file "charges.txt"
  // has to be in the same folder as executable
  // makes sense in combination with amber force field
  else if (option == "chargefile")
  {
	  if (value_string == "1")
	  {
		  Config::set().general.chargefile = true;
	  }
	  else
	  {
		  Config::set().general.chargefile = false;
	  }
  }

  // Input format.
  // Default: TINKER
  else if (option == "inputtype")
    Config::set().general.input = enum_from_string<input_types::T, NUM_INPUT>(input_strings, value_string);

  /////////////////////
  //// Config::energy
  ////////////////////

  // This is an option only valid for force-fields
  // Cutoff radius for non bonded interactions
  // {supported for any internal FF interface}
  // If smaller than 9, it will be set to 10
  // Default: 10000
  else if (option == "cutoff")
  {
    cv >> Config::set().energy.cutoff;
    Config::set().energy.cutoff = Config::get().energy.cutoff < 9.0 ? 10.0 : Config::get().energy.cutoff;
    if (Config::get().periodics.periodic)
    {
      double const min_cut(min(abs(Config::get().periodics.pb_box)) / 2.0);
      if (min_cut > 9.0)
        Config::set().energy.cutoff = min_cut;
    }
    Config::set().energy.switchdist = Config::get().energy.cutoff - 4.0;
  }

  // Radius to start switching function to kick in; scales interactions smoothly to zero at cutoff radius
  // Default: Cutoff - 4.0
  else if (option == "switchdist")
  {
    cv >> Config::set().energy.switchdist;
  }

	else if (option == "xyz_atomtypes")
	{
		Config::set().stuff.xyz_atomtypes = bool_from_iss(cv);
	}

  // Set Verbosity
  // Default: 1
  else if (option == "verbosity")
  {
    cv >> Config::set().general.verbosity;
  }

  //! Task to be performed by CAST
  else if (option == "task")
  {
    config::tasks::T task(Config::getTask(value_string));
    if (task != config::tasks::ILLEGAL)
    {
      Config::set().general.task = task;
    }
    else
    {
      // In case the task was specified numerically (why the hell would
      // you even do that dude?), we try to infer the task by matching
      // the number to the enum
      std::ptrdiff_t identifier(from_string<std::ptrdiff_t>(value_string));
      if (identifier >= 0 && identifier < static_cast<std::ptrdiff_t>(config::NUM_TASKS))
      {
        Config::set().general.task = static_cast<config::tasks::T>(identifier);
      }
    }
  }

  // Output type for structures
  else if (option == "outputtype")
  {
    config::output_types::T type(Config::getOutFormat(value_string));
    if (type != config::output_types::ILLEGAL)
    {
      Config::set().general.output = type;
    }
    else
    {
      // In case the output was specified numerically (why the hell would
      // you even do that dude?), we try to infer the outputtype by matching
      // the number to the enum
      std::ptrdiff_t identifier(from_string<std::ptrdiff_t>(value_string));
      if (identifier >= 0 && identifier < static_cast<std::ptrdiff_t>(config::NUM_OUTPUT))
      {
        Config::set().general.output = static_cast<config::output_types::T>(identifier);
      }
    }

  }

  // Energy calculation interface
  else if (option == "interface")
  {
    interface_types::T inter = Config::getInterface(value_string);
    if (inter != interface_types::ILLEGAL)
    {
      Config::set().general.energy_interface = inter;
    }
    else
    {
      throw std::runtime_error("Configuration contained illegal interface.");
    }
		if (inter == interface_types::QMMM || inter == interface_types::ONIOM || inter == interface_types::THREE_LAYER)
		{
			Config::set().energy.qmmm.use = true;
		}
  }

  // Preoptimizazion energy calculation interface
  else if (option == "preinterface")
  {
    interface_types::T inter = Config::getInterface(value_string);
    Config::set().general.preopt_interface = inter;
  }

  else if (option == "cores")
  {
#if defined _OPENMP
    long const T(from_iss<long>(cv));
    if (T > 0) omp_set_num_threads(T);
#else
    if (Config::get().general.verbosity > 2u)
    std::cout << "CAST was compiled without multithreading. Ignoring the config-option \"cores\"." << std::endl;
#endif
  }

  else if (option == "MOVEmode")
  {
    Config::set().stuff.moving_mode = std::stoi(value_string);
  }

  //! Qmmm-Option
  else if (option.substr(0, 4u) == "QMMM")
  {
    if (option.substr(4u) == "qmatoms")
    {
      Config::set().energy.qmmm.qmatoms =
        sorted_indices_from_cs_string(value_string, true);
    }
    else if (option.substr(4u) == "seatoms")
    {
      Config::set().energy.qmmm.seatoms =
        sorted_indices_from_cs_string(value_string, true);
    }
    else if (option.substr(4u) == "mminterface")
    {
      interface_types::T inter = Config::getInterface(value_string);
      if (inter != interface_types::ILLEGAL)
      {
        Config::set().energy.qmmm.mminterface = inter;
      }
      else
      {
        std::cout << "Configuration contained illegal QMMM MM-interface." << std::endl;
      }
    }
    else if (option.substr(4u) == "seinterface")
    {
      interface_types::T inter = Config::getInterface(value_string);
      if (inter != interface_types::ILLEGAL)
      {
        Config::set().energy.qmmm.seinterface = inter;
      }
      else
      {
        std::cout << "Configuration contained illegal interface for middle layer." << std::endl;
      }
    }
    else if (option.substr(4u) == "qminterface")
    {
      interface_types::T inter = Config::getInterface(value_string);
      if (inter != interface_types::ILLEGAL)
      {
        Config::set().energy.qmmm.qminterface = inter;
      }
      else
      {
        std::cout << "Configuration contained illegal QMMM QM-interface." << std::endl;
      }
    }
	  else if (option.substr(4u) == "writeqmintofile")
	  {
	    if (value_string == "1") Config::set().energy.qmmm.qm_to_file = true;
	  }
	  else if (option.substr(4u) == "linkatomtype")
	  {
	    Config::set().energy.qmmm.linkatom_types.push_back(std::stoi(value_string));
	  }
		else if (option.substr(4u) == "cutoff")
		{
			Config::set().energy.qmmm.cutoff = std::stod(value_string);
		}
		else if (option.substr(4u) == "small_charges")
		{
			Config::set().energy.qmmm.emb_small = std::stoi(value_string);
		}
		else if (option.substr(4u) == "zerocharge_bonds")
		{
			Config::set().energy.qmmm.zerocharge_bonds = std::stoi(value_string);
		}
  }


  //!SPACKMAN
  else if (option == "Spackman")
  {
    std::size_t a(0U), b(0U);
    if (cv >> a && cv >> b && cv >> Config::set().energy.spackman.cut)
    {
      if (a > 0) Config::set().energy.spackman.on = true;
      if (b > 0) Config::set().energy.spackman.interp = true;
    }

  }
  // Options for NEB and pathopt
  else if (option.substr(0, 11) == "NEB-PATHOPT")
  {
	  if (option.substr(11, 6) == "-FINAL")
		  Config::set().neb.FINAL_STRUCTURE = value_string;
	  else if (option.substr(11, 7) == "-IMAGES")
		  cv >> Config::set().neb.IMAGES;
	  else if (option.substr(11, 7) == "-SPRING")
		  cv >> Config::set().neb.SPRINGCONSTANT;
	  else if (option.substr(11, 5) == "-TEMP")
		  cv >> Config::set().neb.TEMPERATURE;
	  else if (option.substr(11, 5) == "-ITER")
		  cv >> Config::set().neb.MCITERATION;
	  else if (option.substr(11, 9) == "-GLOBITER")
		  cv >> Config::set().neb.GLOBALITERATION;
	  else if (option.substr(11, 5) == "-MODE")
		  Config::set().neb.OPTMODE = value_string;
	  else if (option.substr(11, 13) == "-BIASCONSTANT")
		  cv >> Config::set().neb.BIASCONSTANT;
	  else if (option.substr(11, 7) == "-MAXVAR")
		  cv >> Config::set().neb.VARIATION;
	  else if (option.substr(11, 13) == "-ENERGY_RANGE")
		  cv >> Config::set().neb.PO_ENERGY_RANGE;
	  else if (option.substr(11, 9) == "-STEPSIZE")
		  cv >> Config::set().neb.MCSTEPSIZE;
	  else if (option.substr(11, 8) == "-NEBCONN")
		  cv >> Config::set().neb.NEB_CONN;
	  else if (option.substr(11, 15) == "-NUMBER_NEBCONN")
		  cv >> Config::set().neb.CONNECT_NEB_NUMBER;
	  else if (option.substr(11, 8) == "-MIXMOVE")
		  cv >> Config::set().neb.MIXED_MOVE;
	  else if (option.substr(11, 18) == "-CONSTRAINT_GLOBAL")
		  cv >> Config::set().neb.CONSTRAINT_GLOBAL;
	  else if (option.substr(11, 28) == "-CONSTRAINT_NUMBER_DIHEDRALS")
		  cv >> Config::set().neb.NUMBER_OF_DIHEDRALS;
	  else if (option.substr(11, 11) == "-BOND_PARAM")
		  cv >> Config::set().neb.BOND_PARAM;
	  else if (option.substr(11, 4) == "-TAU")
		  cv >> Config::set().neb.TAU;
	  else if (option.substr(11, 9) == "-INT_PATH")
		  cv >> Config::set().neb.INT_PATH;
	  else if (option.substr(11, 9) == "-CLIMBING")
		  cv >> Config::set().neb.CLIMBING;
	  else if (option.substr(11, 7) == "-INT_IT")
		  cv >> Config::set().neb.INT_IT;
	  else if (option.substr(11, 9) == "-NEB-IDPP")
		  cv >> Config::set().neb.IDPP;
	  else if (option.substr(11, 8) == "-MAXFLUX")
		  cv >> Config::set().neb.MAXFLUX;
	  else if (option.substr(11, 11) == "-MF_PATHOPT")
		  cv >> Config::set().neb.MAXFLUX_PATHOPT;
	  else if (option.substr(11, 13) == "-NEB-COMPLETE")
		  cv >> Config::set().neb.COMPLETE_PATH;
	  else if (option.substr(11, 20) == "-NEB-MULTIPLE_POINTS")
		  cv >> Config::set().neb.MULTIPLE_POINTS;
     else if (option.substr(11, 27) == "-NEB-INTERNAL_INTERPOLATION")
      cv >> Config::set().neb.INTERNAL_INTERPOLATION;
	 else if (option.substr(11, 12) == "-NEB-MCM_OPT")
		 cv >> Config::set().neb.MCM_OPT;
	 else if (option.substr(11, 12) == "-NEB-MC_SAVE")
		 cv >> Config::set().neb.MCM_SAVEITER;
	 else if (option.substr(11, 5) == "-CONN")
		 cv >> Config::set().neb.CONN;
  }

  // MOPAC options
  else if (option.substr(0, 5) == "MOPAC")
  {
		if (option.substr(5, 3) == "key")
			Config::set().energy.mopac.command = value_string;
		else if (option.substr(5, 4) == "link")
			Config::set().energy.gaussian.link = value_string;
		else if (option.substr(5, 4) == "path")
			Config::set().energy.mopac.path = value_string;
		else if (option.substr(5, 6) == "delete")
			Config::set().energy.mopac.delete_input = bool_from_iss(cv);
		else if (option.substr(5, 7) == "version")
		{
			// Matching the "value string"
			// to the enum config::mopac_ver_type
			Config::set().energy.mopac.version =
				enum_type_from_string_arr<
				config::mopac_ver_type::T,
				config::NUM_MOPAC_VERSION
				>(value_string, config::mopac_ver_string);
		}
		else if (option.substr(5, 6) == "charge")
			Config::set().energy.mopac.charge = std::stoi(value_string);
  }

  //DFTBaby options
  else if (option.substr(0,7) == "DFTBABY")
  {
    if (option.substr(7,3) == "ath")
        Config::set().energy.dftbaby.path = value_string;
    else if (option.substr(7,8) == "gradfile")
        Config::set().energy.dftbaby.gradfile = value_string;
    else if (option.substr(7,9) == "gradstate")
        Config::set().energy.dftbaby.gradstate = std::stoi(value_string);
    else if (option.substr(7,7) == "verbose")
        Config::set().energy.dftbaby.verbose = std::stoi(value_string);
    else if (option.substr(7,6) == "cutoff")
        Config::set().energy.dftbaby.cutoff = std::stof(value_string);
    else if (option.substr(7,7) == "lr_dist")
        Config::set().energy.dftbaby.lr_dist = std::stof(value_string);
    else if (option.substr(7,7) == "maxiter")
        Config::set().energy.dftbaby.maxiter = std::stoi(value_string);
    else if (option.substr(7,4) == "conv")
        Config::set().energy.dftbaby.conv_threshold = value_string;
    else if (option.substr(7,6) == "states")
        Config::set().energy.dftbaby.states = std::stoi(value_string);
    else if (option.substr(7,7) == "occ_orb")
        Config::set().energy.dftbaby.orb_occ = std::stoi(value_string);
    else if (option.substr(7,8) == "virt_orb")
        Config::set().energy.dftbaby.orb_virt = std::stoi(value_string);
    else if (option.substr(7,9) == "diag_conv")
         Config::set().energy.dftbaby.diag_conv = value_string;
    else if (option.substr(7,12) == "diag_maxiter")
        Config::set().energy.dftbaby.diag_maxiter = std::stoi(value_string);
    else if (option.substr(7,6) == "charge")
        Config::set().energy.dftbaby.charge = std::stoi(value_string);
    else if (option.substr(7,7) == "lr_corr")
    {
       if (value_string == "1")
          Config::set().energy.dftbaby.longrange = true;
    }
    else if (option.substr(4,3) == "opt")
    {
      if (value_string == "1") Config::set().energy.dftbaby.opt = true;
    }
  }

  //DFTB+ Options
  else if (option.substr(0, 5) == "DFTB+")
  {
    if (option.substr(5, 4) == "path"){
      Config::set().energy.dftb.path = value_string;
    }
    else if (option.substr(5, 7) == "skfiles"){
      Config::set().energy.dftb.sk_files = value_string;
    }
    else if (option.substr(5, 7) == "verbose"){
      Config::set().energy.dftb.verbosity = std::stoi(value_string);
    }
    else if (option.substr(5, 6) == "scctol"){
      Config::set().energy.dftb.scctol = std::stod(value_string);
    }
    else if (option.substr(5, 13) == "max_steps_scc"){
      Config::set().energy.dftb.max_steps = std::stoi(value_string);
    }
    else if (option.substr(5, 6) == "charge"){
      Config::set().energy.dftb.charge = std::stoi(value_string);
    }
    else if (option.substr(5, 1) == "3"){
      if (value_string == "1") Config::set().energy.dftb.dftb3 = true;
    }
    else if (option.substr(5, 9) == "optimizer"){
      Config::set().energy.dftb.opt = std::stoi(value_string);
    }
    else if (option.substr(5, 13) == "max_steps_opt"){
      Config::set().energy.dftb.max_steps_opt = std::stoi(value_string);
    }
		else if (option.substr(5, 10) == "fermi_temp") {
			Config::set().energy.dftb.fermi_temp = std::stod(value_string);
		}
  }

	// ORCA options
	else if (option.substr(0, 4) == "ORCA")
	{
		if (option.substr(4, 4) == "path") {
			Config::set().energy.orca.path = value_string;
		}
		else if (option.substr(4, 5) == "nproc") {
			Config::set().energy.orca.nproc = std::stoi(value_string);
		}
    else if (option.substr(4, 7) == "maxcore") {
      Config::set().energy.orca.maxcore = std::stoi(value_string);
    }
		else if (option.substr(4, 6) == "method") {
			Config::set().energy.orca.method = value_string;
		}
		else if (option.substr(4, 8) == "basisset") {
			Config::set().energy.orca.basisset = value_string;
		}
    else if (option.substr(4, 4) == "spec") {
      Config::set().energy.orca.spec = value_string;
    }
		else if (option.substr(4, 6) == "charge") {
			Config::set().energy.orca.charge = std::stoi(value_string);
		}
		else if (option.substr(4, 12) == "multiplicity") {
			Config::set().energy.orca.multiplicity = std::stoi(value_string);
		}
		else if (option.substr(4, 3) == "opt") {
			Config::set().energy.orca.opt = std::stoi(value_string);
		}
    else if (option.substr(4, 7) == "verbose") {
      Config::set().energy.orca.verbose = std::stoi(value_string);
    }
    else if (option.substr(4, 4) == "cube") {
      Config::set().energy.orca.cube_orbs = sorted_indices_from_cs_string(value_string);
    }
		else if (option.substr(4, 6) == "casscf") {
			if (value_string == "1") Config::set().energy.orca.casscf = true;
		}
		else if (option.substr(4, 5) == "nelec") {
			Config::set().energy.orca.nelec = std::stoi(value_string);
		}
		else if (option.substr(4, 4) == "norb") {
			Config::set().energy.orca.norb = std::stoi(value_string);
		}
		else if (option.substr(4, 6) == "nroots") {
			Config::set().energy.orca.nroots = std::stoi(value_string);
		}
		else if (option.substr(4, 2) == "nr") {
			if (value_string == "1") Config::set().energy.orca.nr = true;
		}
		else if (option.substr(4, 5) == "nevpt") {
			if (value_string == "1") Config::set().energy.orca.nevpt = true;
		}
		else if (option.substr(4, 4) == "cpcm") {
			if (value_string == "1") Config::set().energy.orca.cpcm = true;
		}
		else if (option.substr(4, 3) == "eps") {
			Config::set().energy.orca.eps = std::stod(value_string);
		}
		else if (option.substr(4, 6) == "refrac") {
			Config::set().energy.orca.refrac = std::stod(value_string);
		}
	}

  //Gaussian options
  else if (option.substr(0, 8) == "GAUSSIAN")
  {
    if (option.substr(8, 6) == "method")
      Config::set().energy.gaussian.method = value_string;
    else if (option.substr(8, 8) == "basisset")
      Config::set().energy.gaussian.basisset = value_string;
    else if (option.substr(8, 14) == "specifications")
      Config::set().energy.gaussian.spec = value_string;
    else if (option.substr(8, 3) == "chk")
      Config::set().energy.gaussian.chk = value_string;
    else if (option.substr(8, 4) == "link")
      Config::set().energy.gaussian.link = value_string;
    else if (option.substr(8, 6) == "charge")
      Config::set().energy.gaussian.charge = value_string;
    else if (option.substr(8, 12) == "multiplicity")
      Config::set().energy.gaussian.multipl = value_string;
    else if (option.substr(8, 4) == "path")
      Config::set().energy.gaussian.path = value_string;
    else if (option.substr(8, 3) == "opt")
      Config::set().energy.gaussian.opt = bool_from_iss(cv);
    else if (option.substr(8, 5) == "steep")
      Config::set().energy.gaussian.steep = bool_from_iss(cv);
    else if (option.substr(8, 6) == "delete")
      Config::set().energy.gaussian.delete_input = bool_from_iss(cv);
    else if (option.substr(8, 7) == "maxfail")
      Config::set().energy.gaussian.maxfail = std::stoi(value_string);
		else if (option.substr(8, 4) == "cpcm") {
			if (value_string == "1") Config::set().energy.gaussian.cpcm = true;
		}
		else if (option.substr(8, 6) == "epsinf") {
			Config::set().energy.gaussian.epsinf = std::stod(value_string);
		}
		else if (option.substr(8, 3) == "eps") {
			Config::set().energy.gaussian.eps = std::stod(value_string);
		}
  }
  else if (option.substr(0, 9) == "CHEMSHELL") {
	  auto sub_option = option.substr(10);
	  if (sub_option == "path") {
		  Config::set().energy.chemshell.path = value_string;
	  }
      else if (sub_option == "coords") {
          Config::set().energy.chemshell.coords = value_string;
      }
	  else if (sub_option == "pdb") {
		  Config::set().energy.chemshell.extra_pdb = value_string;
	  }
	  else if (sub_option == "inpcrd") {
		  Config::set().energy.chemshell.optional_inpcrd = value_string;
	  }
	  else if (sub_option == "prmtop") {
		  Config::set().energy.chemshell.optional_prmtop = value_string;
	  }
	  else if (sub_option == "babel_path") {
		  Config::set().energy.chemshell.babel_path = value_string;
	  }
	  else if (sub_option == "scheme") {
		  Config::set().energy.chemshell.scheme = value_string;
	  }
	  else if (sub_option == "qm_theory") {
		  Config::set().energy.chemshell.qm_theory = value_string;
	  }
	  else if (sub_option == "qm_hamiltonian") {
		  Config::set().energy.chemshell.qm_ham = value_string;
	  }
	  else if (sub_option == "qm_basisset") {
		  Config::set().energy.chemshell.qm_basis = value_string;
	  }
	  else if (sub_option == "qm_charge") {
		  Config::set().energy.chemshell.qm_charge = value_string;
	  }
	  else if (sub_option == "combine_residues") {
		  Config::set().energy.chemshell.com_residues = value_string;
	  }
	  else if (sub_option == "qm_atoms") {
		  Config::set().energy.chemshell.qm_atoms = value_string;
	  }
	  else if (sub_option == "delete") {
		  Config::set().energy.chemshell.delete_input = bool_from_iss(cv);
	  }
	  else if (sub_option == "maxcyc") {
		  Config::set().energy.chemshell.maxcyc = value_string;
	  }
	  else if (sub_option == "maxcycle") {
		  Config::set().energy.chemshell.maxcycle = value_string;
	  }
	  else if (sub_option == "tolerance") {
		  Config::set().energy.chemshell.tolerance = value_string;
	  }
	  else if (sub_option == "mxlist") {
		  Config::set().energy.chemshell.mxlist = value_string;
	  }
	  else if (sub_option == "cutoff") {
		  Config::set().energy.chemshell.cutoff = value_string;
	  }
	  else if (sub_option == "dispersion_correction") {
		  Config::set().energy.chemshell.dispersion = bool_from_iss(cv);
	  }
      else if (sub_option == "scale14") {
        Config::set().energy.chemshell.scale14 = value_string;
      }
      else if (sub_option == "active_radius") {
        Config::set().energy.chemshell.active_radius = value_string;
      }
  }
  else if (option.substr(0, 5) == "PSI4-") {
    auto sub_option = option.substr(5);
    if(sub_option == "path"){
      Config::set().energy.psi4.path = value_string;
    }
    else if(sub_option == "memory"){
      Config::set().energy.psi4.memory = value_string;
    }
    else if(sub_option == "basis"){
      Config::set().energy.psi4.basis = value_string;
    }
    else if(sub_option == "method"){
      Config::set().energy.psi4.method = value_string;
    }
    else if(sub_option == "spin"){
      Config::set().energy.psi4.spin = value_string;
    }
    else if(sub_option == "charge"){
      Config::set().energy.psi4.charge = value_string;
    }
    else if(sub_option == "threads"){
      Config::set().energy.psi4.threads = value_string;
    }
  }

  // convergence threshold for bfgs
  // Default 0.001
  else if (option == "BFGSgrad")
    cv >> Config::set().optimization.local.bfgs.grad;

  // max number of steps for bfgs
  // Default: 10000
  else if (option == "BFGSmaxstep")
    cv >> Config::set().optimization.local.bfgs.maxstep;
  
  //should trace written into file?
  else if (option == "BFGStrace")
    Config::set().optimization.local.bfgs.trace = bool_from_iss(cv);

  //! STARTOPT
  else if (option == "SOtype")
  {
    Config::set().startopt.type = static_cast<config::startopt::types::T>(from_iss<int>(cv));
  }
  //! STARTOPT
  else if (option == "SOstructures")
  {
    cv >> Config::set().startopt.number_of_structures;
  }


  //! SIMULATION OPTIONS

  else if (option == "Temperature")
  {
    double T(from_iss<double>(cv));
    if (T > 0.0)
    {
      Config::set().optimization.global.temperature = T;
    }
  }

  else if (option == "Tempscale")
  {
    Config::set().optimization.global.temp_scale
      = clip<double>(from_iss<double>(cv), 0.0, 1.0);
  }

  else if (option == "Iterations")
  {
    std::size_t const I(from_iss<std::size_t>(cv));
    if (I > 0)
    {
      Config::set().optimization.global.iterations = I;
    }
    cv >> Config::set().optimization.global.iterations;
  }

  //! SOLVADD
  else if (option.substr(0, 2) == "SA")
  {
    //! default hydrogen bond length
    if (option.substr(2, 2) == "hb")
    {
      cv >> Config::set().startopt.solvadd.defaultLenHB;
    }
    //! maximum number of water atoms
    else if (option.substr(2, 5) == "limit")
    {
      cv >> Config::set().startopt.solvadd.maxNumWater;
    }
    //! maximum number of water atoms
    else if (option.substr(2, 8) == "boundary")
    {
      Config::set().startopt.solvadd.boundary = enum_from_iss<config::startopt_conf::solvadd::boundary_types::T>(cv);
    }
    //! maximum distance of water and initial structure
    else if (option.substr(2, 6) == "radius")
    {
      cv >> Config::set().startopt.solvadd.maxDistance;
    }
    //! forcefield types
    else if (option.substr(2, 5) == "types")
    {
      cv >> Config::set().startopt.solvadd.ffTypeOxygen >> Config::set().startopt.solvadd.ffTypeHydrogen;
    }
    else if (option.substr(2, 3) == "opt")
    {
      Config::set().startopt.solvadd.opt = enum_from_iss<config::startopt_conf::solvadd::opt_types::T>(cv);
    }
    else if (option.substr(2, 7) == "go_type")
    {
      Config::set().startopt.solvadd.go_type =
        enum_type_from_string_arr<
        config::globopt_routine_type::T,
        config::NUM_GLOBOPT_ROUTINES
        >(value_string, config::globopt_routines_str);
    }
    else if (option.substr(2, 7) == "fixinit")
    {
      Config::set().startopt.solvadd.fix_initial = bool_from_iss(cv);
    }
  }
  //! md
  else if (option.substr(0, 2) == "MD")
  {
	  if (option.substr(2, 5) == "steps")
	  {
		  cv >> Config::set().md.num_steps;
	  }
	  else if (option.substr(2, 10) == "integrator")
	  {
		  Config::set().md.integrator = enum_from_iss<config::md_conf::integrators::T>(cv);
	  }
	  else if (option.substr(2, 11) == "trackoffset")
	  {
		  cv >> Config::set().md.trackoffset;
	  }
	  else if (option.substr(2, 5) == "track")
	  {
		  Config::set().md.track = bool_from_iss(cv);
	  }
	  else if (option.substr(2, 8) == "snap_opt")
	  {
		  Config::set().md.optimize_snapshots = bool_from_iss(cv);
	  }
	  else if (option.substr(2, 11) == "snap_buffer")
	  {
		  cv >> Config::set().md.max_snap_buffer;
	  }
	  else if (option.substr(2, 5) == "snap")
	  {
		  cv >> Config::set().md.num_snapShots;
	  }
	  else if (option.substr(2, 9) == "veloscale")
	  {
		  cv >> Config::set().md.veloScale;
	  }
      else if (option.substr(2, 12) == "temp_control")
      {
        if (cv.str() == "0") Config::set().md.temp_control = false;
      }
	  else if (option.substr(2, 10) == "thermostat")
	  {
		  Config::set().md.hooverHeatBath = bool_from_iss(cv);
	  }
	  else if (option.substr(2, 14) == "restart_offset")
	  {
		  cv >> Config::set().md.restart_offset;
	  }
	  else if (option.substr(2, 13) == "refine_offset")
	  {
		  cv >> Config::set().md.refine_offset;
	  }
	  else if (option.substr(2, 6) == "resume")
	  {
		  Config::set().md.resume = bool_from_iss(cv);
	  }
	  else if (option.substr(2, 8) == "timestep")
	  {
		  cv >> Config::set().md.timeStep;
	  }
	  else if (option.substr(2, 17) == "restart_if_broken")
	  {
		  cv >> Config::set().md.broken_restart;
	  }
	  else if (option.substr(2, 5) == "press")
	  {
		  cv >> Config::set().md.pressure;
	  }
	  else if (option.substr(2, 9) == "pcompress")
	  {
		  cv >> Config::set().md.pcompress;
	  }
	  else if (option.substr(2, 6) == "pdelay")
	  {
		  cv >> Config::set().md.pdelay;
	  }
	  else if (option.substr(2, 7) == "ptarget")
	  {
		  cv >> Config::set().md.ptarget;
	  }
	  else if (option.substr(2, 12) == "pre_optimize")
	  {
		  Config::set().md.pre_optimize = bool_from_iss(cv);
	  }
	  else if (option.substr(2, 4) == "heat")
	  {
		  config::md_conf::config_heat heatconf;
		  if (cv >> heatconf.offset && cv >> heatconf.raise)
		  {
			  Config::set().md.heat_steps.push_back(heatconf);
			  Config::set().md.T_init = Config::get().md.heat_steps.front().raise;
			  Config::set().md.T_final = Config::get().md.heat_steps.back().raise;
		  }
	  }
	  else if (option.substr(2, 9) == "spherical")
	  {
		  if (bool_from_iss(cv) && cv >> Config::set().md.spherical.r_inner
			  && cv >> Config::set().md.spherical.r_outer
			  && cv >> Config::set().md.spherical.f1
			  && cv >> Config::set().md.spherical.f2
			  && cv >> Config::set().md.spherical.e1
			  && cv >> Config::set().md.spherical.e2)
		  {
			  Config::set().md.spherical.use = true;
		  }
	  }
	  else if (option.substr(2, 10) == "rattlebond")
	  {
		  std::size_t a, b;
		  cv >> a >> b;

		  config::md_conf::config_rattle::rattle_constraint_bond rattlebondBuffer;
		  rattlebondBuffer.a = a-1;
		  rattlebondBuffer.b = b-1;
		  Config::set().md.rattle.specified_rattle.push_back(rattlebondBuffer);
	  }
		else if (option.substr(2, 10) == "rattledist")
		{
			double value;
			cv >> value;
			Config::set().md.rattle.dists.push_back(value);
		}
		else if (option.substr(2, 20) == "rattle_use_paramfile")
		{
			Config::set().md.rattle.use_paramfile = bool_from_iss(cv);
		}
	  else if (option.substr(2, 6) == "rattle")
	  {
		  int val(0);
		  if (cv >> val)
		  {
			  Config::set().md.rattle.use = (val != 0);
			  if (val == 2) Config::set().md.rattle.all = false;
		  }
	  }
	  else if (option.substr(2, 16) == "biased_potential")
	  {
		  cv >> Config::set().md.set_active_center;
	  }
	  else if (option.substr(2, 11) == "active_site")
	  {
		  unsigned act_cent_atom;
		  if (cv >> act_cent_atom)
		  {
			  if (act_cent_atom == 0 && Config::get().general.task == tasks::FEP)
			  {     // for FEP calculation: if active_site is set to zero:
				    // calculate active site out of all appearing or disappearing atoms
				  Config::set().md.active_center = FEP_get_inout();
			  }
			  else
			  {
				  Config::set().md.active_center.push_back(act_cent_atom);
			  }
		  }
	  }
	  else if (option.substr(2, 6) == "cutoff")
	  {
		  double inner, outer;
		  cv >> inner >> outer;

		  Config::set().md.inner_cutoff = inner;
		  Config::set().md.outer_cutoff = outer;
	  }
	  else if (option.substr(2, 18) == "adjust_by_step")
	  {
		  cv >> Config::set().md.adjustment_by_step;
	  }
    else if (option.substr(2, 8) == "ana_pair")
    {
      std::vector<size_t> vec = sorted_indices_from_cs_string(value_string);
      Config::set().md.ana_pairs.push_back(vec);
    }
    else if (option.substr(2, 13) == "analyze_zones")
    {
      if (value_string == "1") Config::set().md.analyze_zones = true;
    }
    else if (option.substr(2, 9) == "zonewidth")
    {
      Config::set().md.zone_width = std::stod(value_string);
    }
    else if (option.substr(2, 12) == "nosehoover_Q")
    {
      Config::set().md.nosehoover_Q = std::stod(value_string);
    }
  }

  //! dimer

  else if (option.substr(0, 5) == "DIMER")
  {
    if (option.substr(5, 14) == "rotconvergence")
    {
      cv >> Config::set().dimer.rotationConvergence;
    }
    else if (option.substr(5, 8) == "distance")
    {
      cv >> Config::set().dimer.distance;
    }
    else if (option.substr(5, 5) == "maxit")
    {
      cv >> Config::set().dimer.maxRot >> Config::set().dimer.maxStep;
    }
    else if (option.substr(5, 8) == "tflimit")
    {
      cv >> Config::set().dimer.trans_F_rot_limit;
    }
    //else if (option.substr(5,10) == "trans_type")
    //{
    //  Config::set().dimer.trans_type = enum_from_iss<config::dimer::translation_types::T>(cv);
    //}
  }

  else if (option.substr(0, 9) == "Periodics")
  {
    Config::set().periodics.periodic = bool_from_iss(cv, option.substr(0, 9));
    cv >> Config::set().periodics.pb_box.x() >> Config::set().periodics.pb_box.y() >> Config::set().periodics.pb_box.z();
  }
  else if (option.substr(0, 9) == "Periodicp")
  {
    Config::set().periodics.periodic_print = bool_from_iss(cv, option.substr(0, 9));
  }
  else if (option.substr(0, 15) == "PeriodicCutout")
  {
    Config::set().periodics.periodicCutout = bool_from_iss(cv, option.substr(0, 15));
  }
  else if (option.substr(0, 23) == "PeriodicCutoutCriterion")
  {
    cv >> Config::set().periodics.criterion;
  }
  else if (option.substr(0, 23) == "PeriodicCutoutDistance")
  {
    cv >> Config::set().periodics.cutout_distance_to_box;
  }

  else if (option.substr(0, 3) == "FEP")
  {
    if (option.substr(3, 6) == "lambda")
    {
      cv >> Config::set().fep.lambda;
    }
    else if (option.substr(3, 7) == "dlambda")
    {
      cv >> Config::set().fep.dlambda;
    }
    else if (option.substr(3, 9) == "vdwcouple")
    {
      cv >> Config::set().fep.vdwcouple;
    }
    else if (option.substr(3, 10) == "eleccouple")
    {
      cv >> Config::set().fep.eleccouple;
    }
    else if (option.substr(3, 6) == "vshift")
    {
      cv >> Config::set().fep.ljshift;
    }
    else if (option.substr(3, 6) == "cshift")
    {
      cv >> Config::set().fep.cshift;
    }
    else if (option.substr(3, 5) == "equil")
    {
      cv >> Config::set().fep.equil;
    }
    else if (option.substr(3, 5) == "steps")
    {
      cv >> Config::set().fep.steps;
    }
    else if (option.substr(3, 4) == "freq")
    {
      cv >> Config::set().fep.freq;
    }
    else if (option.substr(3, 7) == "analyze")
    {
      std::string a;
      cv >> a;
      if (a == "0") Config::set().fep.analyze = false;
    }
    else if (option.substr(3, 3) == "bar")
    {
      std::string a;
      cv >> a;
      if (a == "1") Config::set().fep.bar = true;
    }
  }

  else if (option.substr(0, 5) == "COORD")
  {
    if (option.substr(5, 7) == "eq_main")
    {
      cv >> Config::set().coords.equals.main;
    }
    else if (option.substr(5, 6) == "eq_int")
    {
      cv >> Config::set().coords.equals.intern;
    }
    else if (option.substr(5, 6) == "eq_xyz")
    {
      cv >> Config::set().coords.equals.xyz;
    }
  }

  else if (option.substr(0, 2) == "GO")
  {
    if (option.substr(2, 10) == "metrolocal")
    {
      Config::set().optimization.global.metropolis_local = bool_from_iss(cv);
    }
    //! Energy range for global optimization tracking
    else if (option.substr(2, 6) == "erange")
    {
      Config::set().optimization.global.delta_e = std::abs(from_iss<double>(cv));
    }
    else if (option.substr(2, 9) == "main_grid")
    {
      cv >> Config::set().optimization.global.grid.main_delta;
    }
    else if (option.substr(2, 9) == "precision")
    {
      Config::set().optimization.global.precision
        = clip<std::size_t>(from_iss<std::size_t>(cv), 4, 30);
    }
    else if (option.substr(2, 8) == "startopt")
    {
      Config::set().optimization.global.pre_optimize = bool_from_iss(cv);
    }
    else if (option.substr(2, 10) == "move_dehyd")
    {
      Config::set().optimization.global.move_dehydrated = bool_from_iss(cv);
    }
    else if (option.substr(2, 14) == "fallback_limit")
    {
      cv >> Config::set().optimization.global.fallback_limit;
    }
    else if (option.substr(2, 18) == "fallback_fr_minima")
    {
      Config::set().optimization.global.selection.included_minima
        = clip<std::size_t>(from_iss<std::size_t>(cv), 1, 50);
    }
    else if (option.substr(2, 18) == "fallback_fr_bounds")
    {
      if (cv >> Config::set().optimization.global.selection.low_rank_fitness)
        cv >> Config::set().optimization.global.selection.high_rank_fitness;
    }
    else if (option.substr(2, 15) == "fallback_fr_fit")
    {
      Config::set().optimization.global.selection.fit_type =
        enum_type_from_string_arr<
        config::optimization_conf::sel::fitness_types::T,
        config::optimization_conf::NUM_FITNESS
        >(value_string, config::optimization_conf::fitness_strings);
    }
    else if (option.substr(2, 8) == "fallback")
    {
      Config::set().optimization.global.fallback =
        enum_type_from_string_arr<
        config::optimization_conf::global::fallback_types::T,
        config::optimization_conf::NUM_FALLBACKS
        >(value_string, config::optimization_conf::fallback_strings);
    }

  }

  else if (option.substr(0, 2) == "RS")
  {
    if (option.substr(2, 10) == "bias_force")
    {
      cv >> Config::set().startopt.ringsearch.bias_force;
    }
    else if (option.substr(2, 12) == "chance_close")
    {
      cv >> Config::set().startopt.ringsearch.chance_close;
    }
    else if (option.substr(2, 10) == "population")
    {
      cv >> Config::set().startopt.ringsearch.population;
    }
    else if (option.substr(2, 11) == "generations")
    {
      cv >> Config::set().startopt.ringsearch.generations;
    }
  }

  // TABU SEARCH

  else if (option.substr(0, 2) == "TS")
  {
    if (option.substr(2, 8) == "mc_first")
    {
      Config::set().optimization.global.tabusearch.mcm_first = bool_from_iss(cv);
    }
    else if (option.substr(2, 16) == "divers_threshold")
    {
      cv >> Config::set().optimization.global.tabusearch.divers_threshold;
    }
    else if (option.substr(2, 12) == "divers_limit")
    {
      cv >> Config::set().optimization.global.tabusearch.divers_limit;
    }
    else if (option.substr(2, 11) == "divers_iter")
    {
      cv >> Config::set().optimization.global.tabusearch.divers_iterations;
    }
  }

  //MONTE CARLO

  else if (option.substr(0, 2) == "MC")
  {
    if (option.substr(2, 9) == "step_size")
    {
      cv >> Config::set().optimization.global.montecarlo.cartesian_stepsize;
    }
    else if (option.substr(2, 12) == "max_dihedral")
    {
      cv >> Config::set().optimization.global.montecarlo.dihedral_max_rot;
    }
    else if (option.substr(2, 12) == "minimization")
    {
      Config::set().optimization.global.montecarlo.minimization = bool_from_iss(cv);
    }
    else if (option.substr(2, 8) == "movetype")
    {
      Config::set().optimization.global.montecarlo.move
        = enum_from_iss<config::optimization_conf::mc::move_types::T>(cv);
    }
    else if (option.substr(2, 5) == "scale")
    {
      cv >> Config::set().optimization.global.temp_scale;
    }
  }

  else if (option.substr(0, 2) == "US")
  {
    if (option.substr(2, 3) == "use")
    {
      Config::set().coords.umbrella.use_comb = bool_from_iss(cv);
    }
    if (option.substr(2, 5) == "equil")
    {
      cv >> Config::set().md.usequil;
    }
    if (option.substr(2, 4) == "snap")
    {
      cv >> Config::set().md.usoffset;
    }
    if (option.substr(2, 7) == "torsion")
    {
      config::coords::umbrellas::umbrella_tor ustorBuffer;
      if (cv >> ustorBuffer.index[0] && cv >> ustorBuffer.index[1]
        && cv >> ustorBuffer.index[2] && cv >> ustorBuffer.index[3]
        && cv >> ustorBuffer.force && cv >> ustorBuffer.angle)
      {
        --ustorBuffer.index[0];
        --ustorBuffer.index[1];
        --ustorBuffer.index[2];
        --ustorBuffer.index[3];
        Config::set().coords.bias.utors.push_back(ustorBuffer);
      }
    }
    if (option.substr(2, 5) == "angle")
    {
      config::coords::umbrellas::umbrella_angle usangleBuffer;
      if (cv >> usangleBuffer.index[0] && cv >> usangleBuffer.index[1]
        && cv >> usangleBuffer.index[2]
        && cv >> usangleBuffer.force && cv >> usangleBuffer.angle)
      {
        --usangleBuffer.index[0];
        --usangleBuffer.index[1];
        --usangleBuffer.index[2];
        Config::set().coords.bias.uangles.push_back(usangleBuffer);
      }
    }
    if (option.substr(2, 4) == "dist")
    {
      config::coords::umbrellas::umbrella_dist usdistBuffer;
      if (cv >> usdistBuffer.index[0] && cv >> usdistBuffer.index[1]
        && cv >> usdistBuffer.force && cv >> usdistBuffer.dist)
      {
        --usdistBuffer.index[0];
        --usdistBuffer.index[1];
        Config::set().coords.bias.udist.push_back(usdistBuffer);
      }
    }
    if (option.substr(2, 4) == "comb")
    {
      config::coords::umbrellas::umbrella_comb uscombBuffer;
      int number_of_dists;
      std::string buffer;

      cv >> number_of_dists >> uscombBuffer.force_final >> uscombBuffer.value;
      for (int i=0; i < number_of_dists; ++i)
      {
        config::coords::umbrellas::umbrella_comb::uscoord dist;
        cv >> buffer >> dist.index1 >> dist.index2 >> dist.factor >> buffer;  // buffer is ( and )
        dist.index1 -= 1;   // because user input starts with 1 and we need numbers starting form 0
        dist.index2 -= 1;
        uscombBuffer.dists.emplace_back(dist);
      }
      Config::set().coords.bias.ucombs.push_back(uscombBuffer);
    }
  }

  //! Fixation excluding
  else if (option.substr(0, 10) == "FIXexclude")
  {
    Config::set().energy.remove_fixed = bool_from_iss(cv, option.substr(0, 10));
  }
  //! Fixation
  else if (option.substr(0, 8) == "FIXrange")
  {
    std::vector<size_t> indicesFromString = sorted_indices_from_cs_string(value_string);
    for (auto &i : indicesFromString) i = i - 1;  // convert atom indizes from tinker numbering (starting with 1) to numbering starting with 0
    Config::set().coords.fixed = indicesFromString;
  }
	else if (option.substr(0, 9) == "FIXsphere")
	{
		int atom_number;
		double radius;
		cv >> atom_number >> radius;
		Config::set().coords.fix_sphere.radius = radius;
		Config::set().coords.fix_sphere.central_atom = atom_number - 1; // convert atom indizes from tinker numbering (starting with 1) to numbering starting with 0
		Config::set().coords.fix_sphere.use = true;
	}

  //! Connect two atoms internally
  else if (option.substr(0, 10) == "INTERNconnect")
  {
    std::size_t a, b;
    if (cv >> a && cv >> b)
    {
      Config::set().coords.internal.connect[a] = b;
      Config::set().coords.internal.connect[b] = a;
    }
  }
  else if (option.substr(0u, 4u) == "MAIN")
  {
    if (option.substr(4u, 9u) == "blacklist")
    {
      std::pair<std::size_t, std::size_t> p;
      if (cv >> p.first && cv >> p.second)
      {
        --p.first;
        --p.second;
        scon::sorted::insert_unique(Config::set().coords.internal.main_blacklist, p);
      }

    }
    else if (option.substr(4u, 9u) == "whitelist")
    {
      std::pair<std::size_t, std::size_t> p;
      if (cv >> p.first && cv >> p.second)
      {
        --p.first;
        --p.second;
        scon::sorted::insert_unique(Config::set().coords.internal.main_whitelist, p);
      }
    }
  }

  else if (option.substr(0, 10) == "REMOVEHROT")
  {
    Config::set().coords.remove_hydrogen_rot = bool_from_iss(cv, option.substr(0, 10));
  }

  else if (option.substr(0, 4) == "BIAS")
  {
    if (option.substr(4, 9) == "spherical")
    {

      config::biases::spherical buffer;
      if (cv >> buffer.radius
        && cv >> buffer.force
        && cv >> buffer.exponent)
      {

        Config::set().coords.bias.spherical.push_back(buffer);
      }
    }

    else if (option.substr(4, 5) == "cubic")
    {

      config::biases::cubic buffer;
      if (cv >> buffer.dim
        && cv >> buffer.force
        && cv >> buffer.exponent)
      {

        Config::set().coords.bias.cubic.push_back(buffer);
      }
    }

    else if (option.substr(4, 3) == "dih")
    {

      config::biases::dihedral biasBuffer;
      if (cv >> biasBuffer.a && cv >> biasBuffer.b
        && cv >> biasBuffer.c && cv >> biasBuffer.d
        && cv >> biasBuffer.ideal && cv >> biasBuffer.force)
      {

        --biasBuffer.a;
        --biasBuffer.b;
        --biasBuffer.c;
        --biasBuffer.d;
        Config::set().coords.bias.dihedral.push_back(biasBuffer);
      }
    }

		else if (option.substr(4, 5) == "angle")
		{
			config::biases::angle biasBuffer;
			if (cv >> biasBuffer.a && cv >> biasBuffer.b && cv >> biasBuffer.c 
				&& cv >> biasBuffer.ideal && cv >> biasBuffer.force)
			{
				--biasBuffer.a;
				--biasBuffer.b;
				--biasBuffer.c;
				Config::set().coords.bias.angle.push_back(biasBuffer);
			}
		}

    else if (option.substr(4, 4) == "dist")
    {

      config::biases::distance biasBuffer;
      if (cv >> biasBuffer.a && cv >> biasBuffer.b
        && cv >> biasBuffer.ideal && cv >> biasBuffer.force)
      {
        --biasBuffer.a;
        --biasBuffer.b;
        Config::set().coords.bias.distance.push_back(biasBuffer);
      }
    }
  } // BIAS

  else if (option.substr(0, 18) == "thresholdpotential")
  {
    config::biases::thresholdstr thrBuffer;
    if (cv >> thrBuffer.th_dist && cv >> thrBuffer.forceconstant)
    {
      Config::set().coords.bias.threshold.push_back(thrBuffer);
    }
  }

  else if (option.substr(0, 9) == "Subsystem")
  {
    // indices from current option value
    std::vector<size_t> ssi = sorted_indices_from_cs_string(value_string);
    // check whether one of the atoms in that subsystem
    // is contained in any already given one
    for (auto & susy : Config::get().coords.subsystems)
    { // cycle present subsystems
      for (auto ssa : susy)
      { // cycle atoms of that subsystem
        if (scon::sorted::exists(ssi, ssa))
        { // check if that atom is in ssi
          throw std::runtime_error(
            "Atom '" + std::to_string(ssa) +
            "' is already part of another subsystem.");
        }
      }
    }
    // push if valid
    Config::set().coords.subsystems.push_back(std::move(ssi));
  }

  else if (option.substr(0, 7) == "2DSCAN-")
  {
    auto const command = option.substr(7);
    if (command == "bond" || command == "angle" || command == "dihedral") {
      Config::set().scan2d.AXES.emplace_back(option.substr(7) + " " + value_string);
    }
    else if (command.substr(0,7) == "PREDEF-") {
      auto const predef_command = command.substr(7);
      if (predef_command == "change_from_atom_to_atom") {
        Config::set().scan2d.change_from_atom_to_atom = std::stod(value_string);
      }
      else if (predef_command == "max_change_to_rotate_whole_molecule") {
        Config::set().scan2d.max_change_to_rotate_whole_molecule = std::stod(value_string);
      }
    }
  }

  //Trajectory Alignment and Analasys options
  else if (option == "dist_unit")
  {
    cv >> Config::set().alignment.dist_unit;
  }
  else if (option == "holm_sand_r0")
  {
    cv >> Config::set().alignment.holm_sand_r0;
  }
  else if (option == "ref_frame_num")
  {
    cv >> Config::set().alignment.reference_frame_num;
  }
  else if (option == "align_external_file")
  {
    cv >> Config::set().alignment.align_external_file;
  }
  else if (option == "traj_align_translational")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE" || holder == "t" || holder == "T" || holder == "1")
    {
      Config::set().alignment.traj_align_translational = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE" || holder == "f" || holder == "F" || holder == "0")
    {
      Config::set().alignment.traj_align_translational = false;
    }
  }
  else if (option == "traj_align_rotational")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE" || holder == "t" || holder == "T" || holder == "1")
    {
      Config::set().alignment.traj_align_rotational = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE" || holder == "f" || holder == "F" || holder == "0")
    {
      Config::set().alignment.traj_align_rotational = false;
    }
  }
  else if (option == "traj_print_bool")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE" || holder == "t" || holder == "T" || holder == "1")
    {
      Config::set().alignment.traj_print_bool = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE" || holder == "f" || holder == "F" || holder == "0")
    {
      Config::set().alignment.traj_print_bool = false;
    }
  }
  // PCA Options
  else if (option == "pca_alignment")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().PCA.pca_alignment = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().PCA.pca_alignment = false;
    }
  }
  else if (option == "pca_read_modes")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().PCA.pca_read_modes = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().PCA.pca_read_modes = false;
    }
  }
  else if (option == "pca_read_vectors")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().PCA.pca_read_vectors = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().PCA.pca_read_vectors = false;
    }
  }
  else if (option == "pca_use_internal")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().PCA.pca_use_internal = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().PCA.pca_use_internal = false;
    }
  }

  else if (option == "pca_trunc_atoms_bool")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().PCA.pca_trunc_atoms_bool = true;

    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().PCA.pca_trunc_atoms_bool = false;
    }
  }
  else if (option == "pca_trunc_atoms_num")
  {
    std::vector<std::string> holder;
    while (cv)
    {
      std::string temp2;
      cv >> temp2;
      holder.push_back(temp2);
    }
    holder.pop_back();
    Config::set().PCA.pca_trunc_atoms_num = configuration_range_int<size_t>(holder);
  }
  else if (option == "pca_start_frame_num")
  {
    cv >> Config::set().PCA.pca_start_frame_num;
  }
  else if (option == "pca_offset")
  {
    cv >> Config::set().PCA.pca_offset;
  }
  else if (option == "pca_ref_frame_num")
  {
    cv >> Config::set().PCA.pca_ref_frame_num;
  }
  else if (option == "pca_internal_dih" && Config::get().PCA.pca_use_internal)
  {
    std::vector<std::string> holder;
    while (cv)
    {
      std::string temp2;
      cv >> temp2;
      holder.push_back(temp2);
    }
    holder.pop_back();
    Config::set().PCA.pca_internal_dih = configuration_range_int<size_t>(holder);
  }
  else if (option == "pca_ignore_hydrogen")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().PCA.pca_ignore_hydrogen = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().PCA.pca_ignore_hydrogen = false;
    }
  }
  else if (option == "pca_print_probability_density")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().PCA.pca_print_probability_density = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().PCA.pca_print_probability_density = false;
    }
  }
  else if (option == "pca_histogram_width")
  {
    cv >> Config::set().PCA.pca_histogram_width;
  }
  else if (option == "pca_histogram_number_of_bins")
  {
    cv >> Config::set().PCA.pca_histogram_number_of_bins;
  }
  else if (option == "pca_dimensions_for_histogramming")
  {
    std::vector<std::string> holder;
    while (cv)
    {
      std::string temp2;
      cv >> temp2;
      holder.push_back(temp2);
    }
    holder.pop_back();
    if (holder.at(0u) != "all")
      Config::set().PCA.pca_dimensions_for_histogramming = configuration_range_int<size_t>(holder);
    else
      Config::set().PCA.pca_histogram_all_marginal_degrees_of_freedom = true;
  }
  else if (option == "proc_desired_start")
  {
    Config::set().PCA.proc_desired_start = configuration_range_float<double>(cv);
  }
  else if (option == "proc_desired_stop")
  {
    Config::set().PCA.proc_desired_stop = configuration_range_float<double>(cv);
  }

  // entropy Options

  else if (option == "entropy_alignment")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().entropy.entropy_alignment = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().entropy.entropy_alignment = false;
    }
  }
  else if (option == "entropy_use_internal")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().entropy.entropy_use_internal = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().entropy.entropy_use_internal = false;
    }
  }
  else if (option == "entropy_trunc_atoms_bool")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().entropy.entropy_trunc_atoms_bool = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().entropy.entropy_trunc_atoms_bool = false;
    }
  }
  else if (option == "entropy_trunc_atoms_num" && Config::get().entropy.entropy_trunc_atoms_bool)
  {
    std::vector<std::string> holder;
    while (cv)
    {
      std::string temp2;
      cv >> temp2;
      holder.push_back(temp2);
    }
    holder.pop_back();
    Config::set().entropy.entropy_trunc_atoms_num = configuration_range_int<size_t>(holder);
  }
  else if (option == "entropy_internal_dih" && Config::get().entropy.entropy_use_internal)
  {
    std::vector<std::string> holder;
    while (cv)
    {
      std::string temp2;
      cv >> temp2;
      holder.push_back(temp2);
    }
    holder.pop_back();
    Config::set().entropy.entropy_internal_dih = configuration_range_int<size_t>(holder);
  }
  else if (option == "entropy_method")
  {
    std::vector<std::string> holder;
    while (cv)
    {
      std::string temp2;
      cv >> temp2;
      holder.push_back(temp2);
    }
    holder.pop_back();
    Config::set().entropy.entropy_method = configuration_range_int<size_t>(holder);
  }
  else if (option == "entropy_method_knn_k")
  {
    cv >> Config::set().entropy.entropy_method_knn_k;
  }
  else if (option == "entropy_offset")
  {
    cv >> Config::set().entropy.entropy_offset;
  }
  else if (option == "entropy_start_frame_num")
  {
    cv >> Config::set().entropy.entropy_start_frame_num;
  }
  else if (option == "entropy_ref_frame_num")
  {
    cv >> Config::set().entropy.entropy_ref_frame_num;
  }
  else if (option == "entropy_temp")
  {
    cv >> Config::set().entropy.entropy_temp;
  }
  else if (option == "entropy_knn_func")
  {
    std::string holder;
    cv >> holder;
    if (holder == "goria" || holder == "Goria" || holder == "GORIA")
    {
      Config::set().entropy.knnfunc = 0;
    }
    else if (holder == "lombardi" || holder == "Lombardi" || holder == "LOMBARDI")
    {
      Config::set().entropy.knnfunc = 1;
    }
  }
  else if (option == "entropy_knn_norm")
  {
    std::string holder;
    cv >> holder;
    if (holder == "maximum" || holder == "Maximum" || holder == "MAXIMUM")
    {
      Config::set().entropy.knnnorm = 1;
    }
  }
  else if (option == "entropy_remove_dof")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().entropy.entropy_remove_dof = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().entropy.entropy_remove_dof = false;
    }
  }
  // NOT IMPLEMENTED AS OF NOW!
  // I/O Atoms index options
  //else if (option == "atomexclude")
  //{
    //Config::set().general.bool_atomsexclude = true;
    //Config::set().general.atomexclude = configuration_makearray<unsigned int>(cv);
  //}
  //

  //IO Options
  else if (option == "amber_mdcrd")
  {
    cv >> Config::set().io.amber_mdcrd;
  }
  else if (option == "amber_mdvel")
  {
    cv >> Config::set().io.amber_mdvel;
  }
  else if (option == "amber_inpcrd")
  {
    cv >> Config::set().io.amber_inpcrd;
  }
  else if (option == "amber_restrt")
  {
    cv >> Config::set().io.amber_restrt;
  }
  else if (option == "amber_trajectory_at_constant_pressure")
  {
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
    {
      Config::set().io.amber_trajectory_at_constant_pressure = true;
    }
    else if (holder == "false" || holder == "False" || holder == "FALSE")
    {
      Config::set().io.amber_trajectory_at_constant_pressure = false;
    }
  }

  /* Inputoptions for exciton_breakup
  */
  else if (option.substr(0u,2u) == "EX")
  {
	  if (option.substr(2u,11u) == "masscenters")
	  {
		  Config::set().exbreak.masscenters = value_string;
	  }
	  else if (option.substr(2u,7u) == "numbern")
	  {
		  cv >> Config::set().exbreak.nscnumber;
	  }
	  else if (option.substr(2u,7u) == "numberp")
	  {
		  cv >> Config::set().exbreak.pscnumber;
	  }
	  else if (option.substr(2u, 11u) == "planeinterf")
	  {
		  cv >> Config::set().exbreak.interfaceorientation;
	  }
	  else if (option.substr(2u, 12u) == "nscpairrates")
	  {
		  Config::set().exbreak.nscpairrates = value_string;
	  }
	  else if (option.substr(2u, 14u) == "pscpairexrates")
	  {
		  Config::set().exbreak.pscpairexrates = value_string;
	  }
	  else if (option.substr(2u, 14u) == "pscpairchrates")
	  {
		  Config::set().exbreak.pscpairchrates = value_string;
	  }
	  else if (option.substr(2u, 13u) == "pnscpairrates")
	  {
		  Config::set().exbreak.pnscpairrates = value_string;
	  }
    else if (option.substr(2u, 10u) == "ReorgE_exc")
    {
      cv >> Config::set().exbreak.ReorgE_exc;
    }
    else if (option.substr(2u, 9u) == "ReorgE_ch")
    {
      cv >> Config::set().exbreak.ReorgE_ch;
    }
    else if (option.substr(2u, 10u) == "ReorgE_nSC")
    {
      cv >> Config::set().exbreak.ReorgE_nSC;
    }
    else if (option.substr(2u, 9u) == "ReorgE_ct")
    {
      cv >> Config::set().exbreak.ReorgE_ct;
    }
   else if (option.substr(2u, 10u) == "ReorgE_rek")
    {
     cv >> Config::set().exbreak.ReorgE_rek;
    }
    else if (option.substr(2u, 13u) == "ct_triebkraft")
    {
      cv >> Config::set().exbreak.ct_triebkraft;
    }
    else if (option.substr(2u, 14u) == "rek_triebkraft")
    {
      cv >> Config::set().exbreak.rek_triebkraft;
    }
    else if (option.substr(2u, 18u) == "oscillatorstrength")
    {
      cv >> Config::set().exbreak.oscillatorstrength;
    }
    else if (option.substr(2u, 10u) == "wellenzahl")
    {
      cv >> Config::set().exbreak.wellenzahl;
    }
  }

  /* Inputoptions for interfacecreation
  */
  else if (option.substr(0u, 2u) == "IC")
  {
    if (option.substr(2u, 4u) == "name")
    {
      Config::set().interfcrea.icfilename = value_string;
    }
    else if (option.substr(2u, 9u) == "inputtype")
    {
      Config::set().interfcrea.icfiletype = enum_from_string<input_types::T, NUM_INPUT>(input_strings, value_string);
    }
    else if (option.substr(2u, 4u) == "axis")
    {
      cv >> Config::set().interfcrea.icaxis;
    }
    else if (option.substr(2u, 8u) == "distance")
    {
      cv >> Config::set().interfcrea.icdist;
    }
  }

  /* Inputoptions for Center
  */
  else if (option.substr(0u, 6u) == "CENTER")
  {
    if (option.substr(6u, 5u) == "dimer")
    {
      Config::set().center.dimer = bool_from_iss(cv);
    }
    else if (option.substr(6u, 8u) == "distance")
    {
      cv >> Config::set().center.distance;
    }
  }

  /* Inputoptions for Couplings
  */
  else if (option.substr(0u, 9u) == "Couplings")
  {
    if (option.substr(9u, 11u) == "dimernumber")
    {
      cv >> Config::set().couplings.nbr_dimPairs;
    }
    else if (option.substr(9u, 9u) == "nSCnumber")
    {
      cv >> Config::set().couplings.nbr_nSC;
    }
    else if (option.substr(9u, 9u) == "pSCnumber")
    {
      cv >> Config::set().couplings.nbr_pSC;
    }
    else if (option.substr(9u, 13u) == "CTcharastates")
    {
      cv >> Config::set().couplings.ct_chara_all;
    }
    else if (option.substr(9u, 6u) == "pSCdim")
    {
      if (option.substr(15u, 12u) == "Multiplicity")
      {
        cv >> Config::set().couplings.pSCmultipl;
      }
      else if (option.substr(15u, 6u) == "Charge")
      {
        cv >> Config::set().couplings.pSCcharge;
      }
      else if (option.substr(15u, 12u) == "ElCalcmethod")
      {
        while (!cv.eof())
        {
          std::string tmp;
          cv >> tmp;
          Config::set().couplings.pSCmethod_el.append(tmp);
          Config::set().couplings.pSCmethod_el.append(" ");
        }
      }
      else if (option.substr(15u, 14u) == "ExciCalcmethod")
      {
        while (!cv.eof())
        {
          std::string tmp;
          cv >> tmp;
          Config::set().couplings.pSCmethod_ex.append(tmp);
          Config::set().couplings.pSCmethod_ex.append(" ");
        }
      }
    }
    else if (option.substr(9u, 6u) == "nSCdim")
    {
      if (option.substr(15u, 12u) == "Multiplicity")
      {
        cv >> Config::set().couplings.nSCmultipl;
      }
      else if (option.substr(15u, 6u) == "Charge")
      {
        cv >> Config::set().couplings.nSCcharge;
      }
      else if (option.substr(15u, 13u) == "holCalcmethod")
      {
        while (!cv.eof())
        {
          std::string tmp;
          cv >> tmp;
          Config::set().couplings.nSCmethod.append(tmp);
          Config::set().couplings.nSCmethod.append(" ");
        }
      }
    }
    else if (option.substr(9u, 9u) == "heterodim")
    {
      if (option.substr(18u, 12u) == "Multiplicity")
      {
        cv >> Config::set().couplings.hetmultipl;
      }
      else if (option.substr(18u, 6u) == "Charge")
      {
        cv >> Config::set().couplings.hetcharge;
      }
      else if (option.substr(18u, 10u) == "Calcmethod")
      {
        while (!cv.eof())
        {
          std::string tmp;
          cv >> tmp;
          Config::set().couplings.hetmethod.append(tmp);
          Config::set().couplings.hetmethod.append(" ");
        }
      }
    }
  }

  /* Inputoptions for Layer_Deposition
  */
  else if (option.substr(0u, 4u) == "LayD")
  {
    if (option.substr(4u, 6u) == "layers")
    {
      cv >> Config::set().layd.amount;
    }
    else if (option.substr(4u, 10u) == "del_number")
    {
      cv >> Config::set().layd.del_amount;
    }
    else if (option.substr(4u, 4u) == "axis")
    {
      cv >> Config::set().layd.laydaxis;
    }
    else if (option.substr(4u, 8u) == "distance")
    {
      cv >> Config::set().layd.layddist;
    }
    else if (option.substr(4u, 9u) == "het_struc")
    {
      Config::set().layd.hetero_option = bool_from_iss(cv);
    }
    else if (option.substr(4u, 8u) == "het_dist")
    {
      cv >> Config::set().layd.sec_layddist;
    }
    else if (option.substr(4u, 10u) == "het_layers")
    {
      cv >> Config::set().layd.sec_amount;
    }
    else if (option.substr(4u, 14u) == "het_del_number")
    {
      cv >> Config::set().layd.sec_del_amount;
    }
    else if (option.substr(4u, 8u) == "het_name")
    {
      Config::set().layd.layd_secname = value_string;
    }
    else if (option.substr(4u, 7u) == "replace")
    {
      Config::set().layd.replace = bool_from_iss(cv);
    }
    else if (option.substr(4u, 10u) == "reference1")
    {
      Config::set().layd.reference1 = value_string;
    }
    else if (option.substr(4u, 10u) == "reference2")
    {
      Config::set().layd.reference2 = value_string;
    }
  }
  
  /* Options for constraint internal coordinates
   */
  else if (option.substr(0u, 10u) == "constraint")
  {
    bool constraint;
    std::string holder;
    cv >> holder;
    if (holder == "true" || holder == "True" || holder == "TRUE")
      constraint = true;
    else if (holder == "false" || holder == "False" || holder == "FALSE")
      constraint = false;
    
    if (option.substr(11u, 12u) == "bond_lengths")
    {
      Config::set().constrained_internals.constrain_bond_lengths = constraint;
    }
    else if (option.substr(11u, 11u) == "bond_angles")
    {
      Config::set().constrained_internals.constrain_bond_angles = constraint;
    }
    else if (option.substr(11u, 9u) == "dihedrals")
    {
      Config::set().constrained_internals.constrain_dihedrals = constraint;
    }
    else if (option.substr(11u, 18u) == "out_of_plane_bends")
    {
      Config::set().constrained_internals.constrain_out_of_plane_bends = constraint;
    }
    else if (option.substr(11u, 12u) == "translations")
    {
      Config::set().constrained_internals.constrain_translations = constraint;
    }
    else if (option.substr(11u, 9u) == "rotations")
    {
      Config::set().constrained_internals.constrain_rotations = constraint;
    }
    else if (option.substr(11u, 10u) == "coordinate")
    {
      std::vector<std::size_t> atom_indices;
      std::istream_iterator<std::size_t> eos, it(cv);
      for(; it != eos; ++it)
      {
        // Only append if atom_indices does not already contain *it
        if (!is_in(*it, atom_indices))
          atom_indices.push_back(*it);
      }
      
      if (atom_indices.size() == 2)
      {
        Config::set().constrained_internals.constrained_bond_lengths.push_back(std::make_pair(std::move(atom_indices), constraint));
      }
      else if (atom_indices.size() == 3)
      {
        Config::set().constrained_internals.constrained_bond_angles.push_back(std::make_pair(std::move(atom_indices), constraint));
      }
      else if (atom_indices.size() == 4)
      {
       Config::set().constrained_internals.constrained_dihedrals.push_back(std::make_pair(std::move(atom_indices), constraint)); 
      }
    }
  }
}



/*


  ########     ###    ########   ######  ########
  ##     ##   ## ##   ##     ## ##    ## ##
  ##     ##  ##   ##  ##     ## ##       ##
  ########  ##     ## ########   ######  ######
  ##        ######### ##   ##         ## ##
  ##        ##     ## ##    ##  ##    ## ##
  ##        ##     ## ##     ##  ######  ########


  #### ##    ## ########  ##     ## ######## ######## #### ##       ########
   ##  ###   ## ##     ## ##     ##    ##    ##        ##  ##       ##
   ##  ####  ## ##     ## ##     ##    ##    ##        ##  ##       ##
   ##  ## ## ## ########  ##     ##    ##    ######    ##  ##       ######
   ##  ##  #### ##        ##     ##    ##    ##        ##  ##       ##
   ##  ##   ### ##        ##     ##    ##    ##        ##  ##       ##
  #### ##    ## ##         #######     ##    ##       #### ######## ########


*/

void Config::parse_file(std::string const & filename)
{
  //TODO not even save commented lines by LBL_FileReader
  auto data = LBL_FileReader(filename).data;
  std::size_t const N(data.size());
  for (std::size_t i = 0; i < N; i++)
  {
    // In the configfile, first there will be an option, then a varying number of whitespaces
    // and then the value of the option
    std::string option_string(data[i].substr(0U, data[i].find_first_of(" ")));
    std::string value_string(data[i]);

    // erase whitespaces
    value_string.erase(0, value_string.find_first_of(" "));
    value_string.erase(0, value_string.find_first_not_of(" "));

    // Values beginning with an "#" are ignored.
    if (option_string.size() > 0 && option_string[0] != '#')
    {
      config::parse_option(option_string, value_string);
    }
  }
}

std::ostream & config::operator << (std::ostream &strm, general const &g)
{
  strm << "Reading structure(s) from '" << g.inputFilename;
  strm << "' (type: " << config::input_strings[g.input] << ")";
  strm << " and parameters from '" << g.paramFilename << "'.\n";
  strm << "Energy calculations will be performed using '" << interface_strings[g.energy_interface] << "'.\n";
  return strm;
}

std::ostream & config::operator<< (std::ostream &strm, coords::conditionsForStructuresToBeConsideredEqual const &equals)
{
  strm << "Two structures will be considered to be equal if either\n";
  strm << " - none of the main torsions differ more then " << equals.main << ", or\n";
  strm << " - no internal vector (bond, angle, dihedral) differs more than " << equals.intern << ", or\n";
  strm << " - no xyz position differs more than " << equals.xyz << ", or\n";
  strm << " - every atom is 'superposed' by an atom with the same atomic number within ";
  strm << equals.superposition << " angstroms.\n";
  return strm;
}

std::ostream & config::operator<< (std::ostream &strm, coords const & coords_in)
{
  if (coords_in.remove_hydrogen_rot)
  {
    strm << "The torsional rotation of a hydrogen atom "
      << "will not be considered as main torsion.\n";
  }
  if (!coords_in.internal.main_blacklist.empty() && coords_in.internal.main_whitelist.empty())
  {
    strm << "Torsional rotation axes ";
    std::size_t k = 1;
    for (auto const & po : coords_in.internal.main_blacklist)
    {
      strm << (k == 1 ? (k % 10 == 0 ? ",\n" : "") : ", ") << po.first << "-" << po.second;
      ++k;
    }
    strm << (coords_in.internal.main_blacklist.size() > 1 ? " are" : " is")
      << " are not considered as main torsions.\n";
  }
  else if (!coords_in.internal.main_whitelist.empty())
  {
    strm << "Torsional rotation around ax" << (coords_in.internal.main_whitelist.size() > 1 ? "es " : "is ");
    std::size_t k = 1;
    for (auto const & po : coords_in.internal.main_whitelist)
    {
      strm << (k == 1 ? (k % 10 == 0 ? ",\n" : "") : ", ") << po.first << "-" << po.second;
      ++k;
    }
    strm << (coords_in.internal.main_whitelist.size() > 1 ? " are" : " is")
      << " exclusively considered for main torsions.\n";
  }
  if (! coords_in.fixed.empty())
  {
    auto const fsi = coords_in.fixed.size();
    if (fsi == 1)
    {
      strm << "1 atom will be fixed: " << coords_in.fixed.front() << '\n';
    }
    else
    {
      strm << fsi << " atoms will be fixed:";
      std::size_t first(0u), last(0u);
      auto const fsim1 = fsi - 1u;
      while (last < fsi)
      {
        while (last < fsim1 && coords_in.fixed[last + 1u] == (coords_in.fixed[last] + 1u)) { ++last; }
        if (last > first)
        {
          strm << " [" << coords_in.fixed[first] + 1 << " to " << coords_in.fixed[last] + 1 << "]";
        }
        else
        {
          strm << " " << coords_in.fixed[first] + 1;
        }
        ++last;
        first = last;
      }
      strm << '\n';
    }
  }

  for (auto const & torsion : coords_in.bias.dihedral)
  {
    strm << "Dihedral " << torsion.a + 1 << "->" << torsion.b + 1 << "->";
    strm << torsion.c + 1 << "->" << torsion.d + 1 << " will be forced to be ";
    strm << torsion.ideal << " deg with force =  " << torsion.force << ".\n";
  }

  for (auto const & dist : coords_in.bias.distance)
  {
    strm << "Distance " << dist.a << "<->" << dist.b;
    strm << " will be forced to be ";
    strm << dist.ideal << " A. Force =  " << dist.force << "\n";
  }

  for (auto const & angle : coords_in.bias.angle)
  {
    strm << "Angle " << angle.a << "->" << angle.b << "<-" << angle.c;
    strm << " will be forced to be ";
    strm << angle.ideal << " A. Force =  " << angle.force << "\n";
  }

  for (auto const & sphere : coords_in.bias.spherical)
  {
    strm << "Spherical boundary with radius " << sphere.radius;
    strm << " will be applied; Force =  " << sphere.force;
    strm << ", Exponent = " << sphere.exponent << "\n";
  }

  for (auto const & cube : coords_in.bias.cubic)
  {
    strm << "Cubic boundary with box size " << cube.dim;
    strm << " will be applied; Force =  " << cube.force;
    strm << ", Exponent = " << cube.exponent << "\n";
  }

  return strm;
}

std::ostream & config::operator<< (std::ostream &strm, energy const &p)
{
  if (p.cutoff < 1000.0) strm << "Cutoff radius of " << p.cutoff << " Angstroms (switching to zero starting at " << p.switchdist << " Angstroms) applied.\n";
  else strm << "No cutoff radius applied.\n";
  if (p.remove_fixed)
  {
    strm << "Nonbonded terms between fixed atoms will be excluded in internal forcefield calculations.\n";
  }
  if (p.spackman.on)
  {
    strm << "Spackman correction applied.\n";
  }
  if (Config::get().general.energy_interface == interface_types::MOPAC)
  {
    strm << "Mopac path is '" << p.mopac.path << "' and command is '" << p.mopac.command << "'.\n";
  }
  else if (Config::get().general.energy_interface == interface_types::GAUSSIAN)
  {
    strm << "Gaussian path is '" << p.gaussian.path << "' and command is '# "
         << p.gaussian.method << " " << p.gaussian.basisset << " " << p.gaussian.spec << "'.\n";
  }
  else if (Config::get().general.energy_interface == interface_types::CHEMSHELL) {
	  strm << "Chemshell path is '" << p.chemshell.path << "'.\n"; //<- Not done here!!!!!!
  }
  else if(Config::get().general.energy_interface == interface_types::PSI4){
    strm << "Psi4 path is '" << p.psi4.path << "'\n"; //<- Not done here!!!!!!
  }
	else if (Config::get().general.energy_interface == interface_types::DFTB) {
		strm << "Path to DFTB+ is '" << p.dftb.path << "'\n"; 
	}
	else if (Config::get().general.energy_interface == interface_types::ORCA) {
		strm << "Path to ORCA is '" << p.orca.path << "'\n";
	}
  return strm;
}

std::ostream & config::operator<< (std::ostream &strm, periodics const &p)
{
  if (p.periodic)
  {
    strm << "Periodics box [ x, y, z ] " << p.pb_box << " applied.\n";
    if (p.periodicCutout)
    {
      strm << "Molecules ";
      if (p.criterion == 1u)
        strm << "of which the center of mass is ";
      else if (p.criterion == 0u)
        strm << "which have atoms that are ";
      if (p.cutout_distance_to_box == 0.)
        strm << "outside the periodic box ";
      else
        strm << "closer than " << p.cutout_distance_to_box << "Angstrom to periodic box edges ";
      strm << "are removed.\n";
    }
  }
  return strm;
}

std::ostream& config::optimization_conf::operator<< (std::ostream &strm, global const &opt)
{
  strm << "At most " << opt.iterations << " global optimization iterations at " << opt.temperature;
  strm << "K (multiplied by " << opt.temp_scale << " for each accepted minimum) will be performed.\n";
  strm << "All structures within " << opt.delta_e << " kcal/mol above the lowest minimum will be saved.\n";
  strm << "TabuSearch will use " << opt.tabusearch.divers_iterations << " iterations of";
  switch (opt.tabusearch.divers_optimizer)
  {
    case go_types::MCM:
    {
      strm << " MCM";
      break;
    }
    default:
    {
      break;
    }
  }
  strm << " (at max. " << opt.tabusearch.divers_limit << " times)";
  strm << " for diversification after " << opt.tabusearch.divers_threshold << " TS steps failed to accept a minimum.\n";
  strm << "Monte Carlo will be performed by random";
  switch (opt.montecarlo.move)
  {
    case mc::move_types::XYZ:
    {
      strm << " cartesian";
      strm << " (max. distortion " << opt.montecarlo.cartesian_stepsize << " angstroms)";
      break;
    }
    case mc::move_types::DIHEDRAL:
    {
      strm << " dihedral";
      strm << " (max. distortion " << opt.montecarlo.dihedral_max_rot << " degrees)";
      break;
    }
    case mc::move_types::WATER:
    {
      strm << " water";
      break;
    }
    default:
    {
      strm << " dihedral (strained optimization, ";
      strm << "max. distortion " << opt.montecarlo.dihedral_max_rot << " degrees)";
      break;
    }
  }
  strm << " movement" << (opt.montecarlo.minimization ? " with" : " without") << " subsequent local optimization.\n";
  if (opt.pre_optimize) strm << "The system will be pre-optimized (startopt) prior to optimization.\n";
  strm << "If the current iteration fails to find a new (accepted) minimum, \n";
  if (opt.fallback == global::fallback_types::LAST_GLOBAL)
  {
    strm << "the routine will use the last minimum as its starting point at most " << opt.fallback_limit << " times.\n";
    strm << "If the limit is reached, the lowest minimum available will be selected.\n";
    strm << "No minimum will be used more than " << opt.fallback_limit << " times.\n";
  }
  else if (opt.fallback == global::fallback_types::FITNESS_ROULETTE)
  {
    strm << "a new starting point will be selected using roulette selection,\n";
    strm << "utilizing a rank-based fitness function among the ";
    strm << opt.selection.included_minima << " lowest, accepted minima.\n";
    strm << "The minimum which is lowest in energy (rank 1) will have a fitness value of  " << opt.selection.high_rank_fitness;
    strm << ", while rank " << opt.selection.included_minima << " will have fitness " << opt.selection.low_rank_fitness << '\n';
    strm << "Fitness interpolation type is '" << fitness_strings[opt.selection.fit_type >= 0 ? opt.selection.fit_type : 0] << "'.\n";
    strm << "No minimum will be selected more than " << opt.fallback_limit << " times, ";
    strm << " while not more than 100 approaches to find a valid minimum will be performed.\n";
  }
  else strm << "invalid action is performed.\n";
  return strm;
}

std::ostream& config::startopt_conf::operator<< (std::ostream &strm, solvadd const &sa)
{
  strm << "SolvAdd will fill ";
  if (sa.boundary == startopt_conf::solvadd::boundary_types::BOX)
  {
    strm << "a water box with a required side length of " << sa.maxDistance << "\n";
  }
  else if (sa.boundary == startopt_conf::solvadd::boundary_types::SPHERE)
  {
    strm << "a water sphere with a required radius of " << sa.maxDistance << "\n";
  }
  else if (sa.boundary == startopt_conf::solvadd::boundary_types::LAYER)
  {
    strm << "a water layer with a required thickness of " << sa.maxDistance << "\n";
  }
  if (sa.maxNumWater == 0)
  {
    strm << "The boundary distance is binding.\n";
  }
  else
  {
    strm << "The solvation shell will contain " << sa.maxNumWater << " water molecules in the end.\n";
    strm << "if the boundary does not allow placing all water molecules, the boundary is expanded in steps of 2A.\n";
  }
  strm << "The placed water molecules have a bond length of " << sa.water_bond << " A and a bond angle of " << sa.water_angle << " degrees\n";
  strm << "Oxygen will be FF type " << sa.ffTypeOxygen << " whereas hydrogen will be " << sa.ffTypeHydrogen << "\n";
  if (sa.opt != startopt_conf::solvadd::opt_types::NONE)
  {
    if (sa.opt == startopt_conf::solvadd::opt_types::SHELL || sa.opt == startopt_conf::solvadd::opt_types::TOTAL_SHELL)
    {
      strm << "A local optimization will be performed after each iteration ";
      if (sa.fix_intermediate) strm << "(where the water of the previous iterations will be fixed)";
      if (sa.opt == startopt_conf::solvadd::opt_types::TOTAL_SHELL)
      {
        strm << " and once the hydration is complete";
      }
      strm << ".\n";
    }
    if (sa.opt == startopt_conf::solvadd::opt_types::TOTAL)
    {
      strm << "The system will be relaxed after placing all water molecules.\n";
    }

    if (sa.fix_initial) strm << "The initial geometry will be fixed.\n";
  }
  return strm;
}

std::ostream& config::startopt_conf::operator<< (std::ostream &strm, ringsearch const &rso)
{
  strm << "Ringsearch will propagate ";
  strm << rso.population << " individual structures for ";
  strm << rso.generations << " generations (evolutionary).\n";
  strm << "The chance to initially close a ring will be " << rso.chance_close;
  strm << " and a force of " << rso.bias_force << " will be applied to close a ring.\n";
  return strm;
}

std::ostream& config::operator<< (std::ostream &strm, startopt const &sopt)
{
  if (sopt.type == config::startopt::types::T::RINGSEARCH ||
    sopt.type == config::startopt::types::T::RINGSEARCH_SOLVADD)
  {
    strm << "Startopt with Ringsearch.\n";
    strm << sopt.ringsearch;
  }
  if (sopt.type == config::startopt::types::T::SOLVADD ||
    sopt.type == config::startopt::types::T::RINGSEARCH_SOLVADD)
  {
    strm << "Startopt with Solvadd.\n";
    strm << sopt.solvadd;
  }
  return strm;
}
