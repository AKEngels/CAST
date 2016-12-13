#include "configuration.h"

/**
 * Global static instance of the config-object. 
 * There can only ever be one config-object.
 */
Config * Config::m_instance = nullptr;

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

config::solvs::S Config::getSolv(std::string const & S)
{
  for (std::size_t i = 0; i < config::NUM_SOLV; ++i)
  {
    if (S.find(config::solv_strings[i]) != S.npos)
      return static_cast<config::solvs::S>(i);
  }
  return config::solvs::VAC;
}

config::surfs::SA Config::getSurf(std::string const & S)
{
  for (std::size_t i = 0; i < config::NUM_SURF; ++i)
  {
    if (S.find(config::surf_strings[i]) != S.npos)
      return static_cast<config::surfs::SA>(i);
  }
  return config::surfs::TINKER;
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
  // Default: oplsaa.prm
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
    if (Config::get().energy.periodic)
    {
      double const min_cut(min(abs(Config::get().energy.pb_box)) / 2.0);
      if (min_cut > 9.0) Config::set().energy.cutoff = min_cut;
    }
    Config::set().energy.switchdist = Config::get().energy.cutoff - 4.0;
  }
  // Turn particle mesh ewald on
  // CURRENTLY OUT OF ORDER
  // Default: Off
  else if (option == "PME")
  {
    cv >> Config::set().energy.pme;
  }

  // Radius to start switching function to kick in; scales interactions smoothly to zero at cutoff radius
  // Default: Cutoff - 4.0
  else if (option == "switchdist")
  {
    cv >> Config::set().energy.switchdist;
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

  //! Methods for implicit solvation

  // Method for solvation
  // Default: VAC (i guess vacuum?)
  else if (option == "solvmethod")
  {
    config::solvs::S solvmethod(Config::getSolv(value_string));
    Config::set().general.solvationmethod = solvmethod;
  }
  else if (option == "surface")
  {
    config::surfs::SA surfmethod(Config::getSurf(value_string));
    Config::set().general.surfacemethod = surfmethod;
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
      std::cout << "Configuration contained illegal interface." << std::endl;
      std::cout << "Using default energy interface: OPLSAA." << std::endl;
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
  }

  // MOPAC options
  else if (option.substr(0, 5) == "MOPAC")
  {
    if (option.substr(5, 3) == "key")
      Config::set().energy.mopac.command = value_string;
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
  }

  // convergence threshold for bfgs
  // Default 0.001
  else if (option == "BFGSgrad")
    cv >> Config::set().optimization.local.bfgs.grad;

  // max number of steps for bfgs
  // Default: 10000
  else if (option == "BFGSmaxstep")
    cv >> Config::set().optimization.local.bfgs.maxstep;


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
		  rattlebondBuffer.a = a;
		  rattlebondBuffer.b = b;
		  Config::set().md.rattle.specified_rattle.push_back(rattlebondBuffer);
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
	  else if (option.substr(2, 7) == "rattpar")
	  {
		  cv >> Config::set().md.rattle.ratpar;
		  std::cout << "Rattle pars:  " << Config::get().md.rattle.ratpar << std::endl;
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
    Config::set().energy.periodic = bool_from_iss(cv);
    if (cv >> Config::set().energy.pb_box.x()
      && cv >> Config::set().energy.pb_box.y()
      && cv >> Config::set().energy.pb_box.z() )
    {
      double const min_cut(min(abs(Config::get().energy.pb_box)) / 2.0);
      if (Config::set().energy.periodic
        && Config::get().energy.cutoff > min_cut)
      {
        Config::set().energy.cutoff = min_cut;
      }
    }
  }
  else if (option.substr(0, 9) == "Periodicp")
  {
    Config::set().energy.periodic_print = bool_from_iss(cv);
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
  }

  //! Fixation excluding
  else if (option.substr(0, 10) == "FIXexclude")
  {
    Config::set().energy.remove_fixed = bool_from_iss(cv);
  }
  //! Fixation 
  else if (option.substr(0, 8) == "FIXrange")
  {
    std::size_t start(0), end(0);
    if (cv >> start && cv >> end && start > 0 && end > start)
    {
      for (std::size_t a(start - 1u); a < end; ++a)
      {
        //std::cout << "RangeFIXING: " << a << "\n";
        scon::sorted::insert_unique(Config::set().coords.fixed, a);
      }
    }
  }

  else if (option.substr(0, 7) == "ATOMFIX")
  {
    auto fixed = from_iss<std::size_t>(cv) - 1u;
    //std::cout << "ATOMIXing: " << fixed << "\n";
    scon::sorted::insert_unique(Config::set().coords.fixed, fixed);
  }

  else if (option.substr(0, 5) == "ADJUST")
  {
    if (option.substr(5, 3) == "dih")
    {
      config::adjust_conf::dihedral ald;
      if (cv >> ald.a && cv >> ald.b
        && cv >> ald.c && cv >> ald.d && cv >> ald.value)
      {
        --ald.a; --ald.b; --ald.c; --ald.d;
        Config::set().adjustment.dihedrals.push_back(ald);
      }
    }
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
    Config::set().energy.remove_fixed = bool_from_iss(cv);
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
    Config::set().coords.remove_hydrogen_rot = bool_from_iss(cv);
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
        biasBuffer.forward = bool_from_iss(cv);
        Config::set().coords.bias.dihedral.push_back(biasBuffer);
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

  // D U S T I N S T U F F 

  //Trajectory Alignment and Analasys options
  else if (option == "dist_unit")
  {
    cv >> Config::set().alignment.dist_unit;
  }
  else if (option == "holm_sand_r0")
  {
    cv >> Config::set().alignment.holm_sand_r0;
  }
  //else if (option == "cdist_cutoff")
  //{
    //cv >> Config::set().alignment.cdist_cutoff;
  //}
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
    Config::set().PCA.pca_trunc_atoms_num = configuration_range_int<size_t>(cv);
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
    Config::set().PCA.pca_internal_dih = configuration_range_int<size_t>(cv);
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
    Config::set().PCA.pca_dimensions_for_histogramming = configuration_range_int<size_t>(cv);
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
    Config::set().entropy.entropy_trunc_atoms_num = configuration_range_int<size_t>(cv);
  }
  else if (option == "entropy_internal_dih" && Config::get().entropy.entropy_use_internal)
  {
    Config::set().entropy.entropy_internal_dih = configuration_range_int<size_t>(cv);
  }
  else if (option == "entropy_method")
  {
    Config::set().entropy.entropy_method = configuration_range_int<size_t>(cv);
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

std::ostream & config::operator<< (std::ostream &strm, coords::eqval const &equals)
{
  strm << "Two structures will be considered to be equal if either\n";
  strm << " - none of the main torsions differ more then " << equals.main << ", or\n";
  strm << " - no internal vector (bond, angle, dihedral) differs more than " << equals.intern << ", or\n";
  strm << " - no xyz position differs more than " << equals.xyz << ", or\n";
  strm << " - every atom is 'superposed' by an atom with the same atomic number within ";
  strm << equals.superposition << " angstroms.\n";
  return strm;
}

std::ostream & config::operator<< (std::ostream &strm, coords const &p)
{
  if (p.remove_hydrogen_rot)
  {
    strm << "The torsional rotation of a hydrogen atom "
      << "will not be considered as main torsion.\n";
  }
  if (!p.internal.main_blacklist.empty() && p.internal.main_whitelist.empty())
  {
    strm << "Torsional rotation axes ";
    std::size_t k = 1;
    for (auto const & po : p.internal.main_blacklist)
    {
      strm << (k == 1 ? (k % 10 == 0 ? ",\n" : "") : ", ") << po.first << "-" << po.second;
      ++k;
    }
    strm << (p.internal.main_blacklist.size() > 1 ? " are" : " is")
      << " are not considered as main torsions.\n";
  }
  else if (!p.internal.main_whitelist.empty())
  {
    strm << "Torsional rotation around ax" << (p.internal.main_whitelist.size() > 1 ? "es " : "is ");
    std::size_t k = 1;
    for (auto const & po : p.internal.main_whitelist)
    {
      strm << (k == 1 ? (k % 10 == 0 ? ",\n" : "") : ", ") << po.first << "-" << po.second;
      ++k;
    }
    strm << (p.internal.main_whitelist.size() > 1 ? " are" : " is")
      << " exclusively considered for main torsions.\n";
  }
  if (!p.fixed.empty())
  {
    auto const fsi = p.fixed.size();
    if (fsi == 1)
    {
      strm << "1 atom will be fixed: " << p.fixed.front() << '\n';
    }
    else
    {
      strm << fsi << " atoms will be fixed:";
      std::size_t first(0u), last(0u);
      auto const fsim1 = fsi - 1u;
      while (last < fsi)
      {
        while (last < fsim1 && p.fixed[last + 1u] == (p.fixed[last] + 1u)) { ++last; }
        if (last > first)
        {
          strm << " [" << p.fixed[first] + 1 << " to " << p.fixed[last] + 1 << "]";
        }
        else
        {
          strm << " " << p.fixed[first] + 1;
        }
        ++last;
        first = last;
      }
      strm << '\n';
    }
  }

  if (Config::get().general.task == tasks::UMBRELLA && (!p.umbrella.torsions.empty() || !p.umbrella.distances.empty()))
  {
    strm << "Umbrella Sampling with " << " steps and snapshots every " << p.umbrella.snap_offset << " steps.\n";
    if (!p.umbrella.torsions.empty())
    {
      strm << "Umbrella torsions:\n";
      for (auto const & torsion : p.umbrella.torsions)
      {
        strm << "[UT] Indices: " << torsion.index[0] << ", " << torsion.index[1] << ", " << torsion.index[2] << ", " << torsion.index[3];
        strm << ". Start: " << " - End: " << ". Step: " << ". \n";
      }
    }
    if (!p.umbrella.distances.empty())
    {
      strm << "Umbrella distances:\n";
      for (auto const & dist : p.umbrella.distances)
      {
        strm << "[UD] Indices: " << dist.index[0] << ", " << dist.index[1];
      }
    }
  }

  for (auto const & torsion : p.bias.dihedral)
  {
    strm << "Dihedral " << torsion.a + 1 << "->" << torsion.b + 1 << "->";
    strm << torsion.c + 1 << "->" << torsion.d + 1 << " will be forced to be ";
    strm << torsion.ideal << " deg with force =  " << torsion.force << ".\n";
  }

  for (auto const & dist : p.bias.distance)
  {
    strm << "Distance " << dist.a << "<->" << dist.b;
    strm << " will be forced to be ";
    strm << dist.ideal << " A. Force =  " << dist.force << "\n";
  }

  for (auto const & angle : p.bias.angle)
  {
    strm << "Angle " << angle.a << "->" << angle.b << "<-" << angle.c;
    strm << " will be forced to be ";
    strm << angle.ideal << " A. Force =  " << angle.force << "\n";
  }

  for (auto const & sphere : p.bias.spherical)
  {
    strm << "Spherical boundary with radius " << sphere.radius;
    strm << " will be applied; Force =  " << sphere.force;
    strm << ", Exponent = " << sphere.exponent << "\n";
  }

  for (auto const & cube : p.bias.cubic)
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
  if (p.periodic)
  {
    strm << "Periodics box [ x, y, z] " << p.pb_box << " applied.\n";
  }
  if (p.spackman.on)
  {
    strm << "Spackman correction applied.\n";
  }
  if (Config::get().general.energy_interface == interface_types::MOPAC)
  {
    strm << "Mopac path is '" << p.mopac.path << "' and command is '" << p.mopac.command << "'.\n";
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