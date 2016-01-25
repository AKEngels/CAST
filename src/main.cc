#include <cstdlib>
#include <fstream>
#include <memory>
#include "global.h"
#include "configuration.h"
#include "coords_io.h"
#include "scon_chrono.h"
// Task items
#include "startopt_solvadd.h"
#include "startopt_ringsearch.h"
#include "md.h"
#include "optimization_global.h"
#include "pathopt.h"
#include "Path_perp.h"
#include "reaccoord.h"
#include "scon_log.h"
#include "matop.h" //For ALIGN, PCAgen, ENTROPY, PCAproc
//#include "gbsa.h"
#include <omp.h>


#if defined(_MSC_VER)
#include <process.h>
#define pid_func _getpid
#else 
#include <unistd.h>
#define pid_func getpid
#endif

//#define CAST_DEBUG_DROP_EXCEPTIONS


int main(int argc, char **argv)
{

  std::ios::sync_with_stdio(false);
  std::cout.sync_with_stdio(false);

#ifndef CAST_DEBUG_DROP_EXCEPTIONS
  try
  {
#endif

    /*
    Preparation
    */

    std::cout << scon::c3_delimeter(',');

    // start execution and initialization timer
    scon::chrono::high_resolution_timer exec_timer, init_timer;

    // initialize (old) RNG
    srand((unsigned int)time(NULL)+pid_func());

    // Parse config file and command line 
    auto config_filename = config::config_file_from_commandline(argc, argv);
    Config main_configuration(config_filename);
    config::parse_command_switches(argc, argv);

    // Print configuration
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << "-------------------------------------------------" << lineend;
      std::cout << "Configuration ('" << config_filename << "')" << lineend;
      std::cout << "-------------------------------------------------" << lineend;
      std::cout << Config::get().general;
      std::cout << Config::get().coords;
      std::cout << Config::get().energy;
    }

    /*

    Initialization

    */

    // read coordinate input file

    std::unique_ptr<coords::input::format> ci(coords::input::new_format());
    
    //coords::input::format * ci(coords::input::new_format());
    coords::Coordinates coords(ci->read(Config::get().general.inputFilename));

	  // setting the methods for implicit solvation
	  //GB::born::set(coords);
	  //GB::born::SET_METHOD();
	  //GB::born::SET_SURFACE();

    if (Config::get().general.verbosity > 1U)
    {
      std::cout << "-------------------------------------------------" << lineend;
      std::cout << "Initialization" << lineend;
      std::cout << "-------------------------------------------------" << lineend;
      std::cout << "Loaded " << ci->size() << " structure" << (ci->size() == 1 ? "" : "s");
      std::cout << ". (" << ci->atoms() << " atom" << (ci->atoms() == 1 ? "" : "s");
      std::cout << (ci->size() > 1 ? " each" : "") << ")" << lineend;
      std::size_t const susysize(coords.subsystems().size());
      if (susysize > 1U)
      {
        std::cout << susysize << " subsystems: ";
        for (std::size_t i(0U); i < susysize; ++i)
        {
          std::size_t const atms(coords.subsystems(i).size());
          std::cout << "[#" << i + 1 << " with " << atms << " atom" << (atms == 1 ? ".]" : "s.]");
        }
        std::cout << lineend;
      }
    }
    // ... 
    if (Config::get().energy.periodic)
    {
      for (auto & pes : *ci)
      {
        coords.set_xyz(pes.structure.cartesian);
        coords.move_all_by(-coords.center_of_mass());
        pes = coords.pes();
      }
    }


	// Initialize PME STUFF
	/* if (Config::get().energy.pme == true)
	{
	if (Config::get().energy.periodic == false)
	{
	std::cout << "PME can only be used with Periodic Boundary Conditions! Check your INPUTFILE!" << std::endl;
	exit(0);
	}
	std::cout << "Initializing PME parameters" << std::endl;
	int size = coords.size();
	coords.pme_stuff(size);
	}*/



    // stop and print initialization time
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << "-------------------------------------------------" << lineend;
      std::cout << "Initialization done after " << init_timer << lineend;
    }
    /*

    Preoptimization

    */
    if (coords.preoptimize())
    {
      if (Config::get().general.verbosity > 1U)
      {
        std::cout << "-------------------------------------------------" << lineend;
        std::cout << "Preoptimization:" << lineend;
        std::cout << "-------------------------------------------------" << lineend;
      }
      std::size_t i(0);
      for (auto const & pes : *ci)
      {
        coords.set_xyz(pes.structure.cartesian);
        coords.pe();
        if (Config::get().general.verbosity > 1U)
        {
          std::cout << "Preoptimization initial: " << ++i << lineend;
          coords.e_tostream_short(std::cout, coords.preinterface());
        }
        coords.po();
        if (Config::get().general.verbosity > 1U)
        {
          std::cout << "Preoptimization post-opt: " << i << lineend;
          coords.e_tostream_short(std::cout, coords.preinterface());
        }
      }
    }

    /*

    Tasks

    */

    std::cout << "-------------------------------------------------" << lineend;
    std::cout << "Task '" << config::task_strings[Config::get().general.task];
    std::cout << "' (" << Config::get().general.task << ") computation:" << lineend;
    std::cout << "-------------------------------------------------" << lineend;
    // start task timer
    scon::chrono::high_resolution_timer task_timer;
    // select task
    switch (Config::get().general.task)
    {
      // DEVTEST: Room for Development testing
    case config::tasks::DEVTEST:
    {
      break;
    }
    case config::tasks::SP:
      { // singlepoint
        coords.e_head_tostream_short(std::cout);
        std::size_t i(0u);
        auto sp_energies_fn = coords::output::filename("SP", ".txt");
        std::ofstream sp_estr(sp_energies_fn, std::ios_base::out);
        if (!sp_estr) throw std::runtime_error("Cannot open '" + sp_energies_fn + "' to write SP energies.");
        for (auto const & pes : *ci)
        {
          using namespace std::chrono;
          coords.set_xyz(pes.structure.cartesian);
          auto start = high_resolution_clock::now();
          auto en = coords.e();
          auto tim = duration_cast<duration<double>>
            (high_resolution_clock::now() - start);
          std::cout << tim.count() << " ns\n";
          std::cout << "Structure " << ++i << lineend;
          sp_estr << i << ' ' << en << ' ' << tim.count() << '\n';
          coords.e_tostream_short(std::cout);
        }
        break;
      }
    case config::tasks::ADJUST:
      { // alignment / change strucutre
        coords.e_head_tostream_short(std::cout);
        std::size_t i(0u);
        std::ofstream outputstream(coords::output::filename("_ADJUSTED").c_str(), std::ios_base::out);
        for (auto const & pes : *ci)
        {
          coords.set_xyz(pes.structure.cartesian);
          coords.to_internal();
          for (auto const &d : Config::get().adjustment.dihedrals)
          {
            std::size_t const di = coords.atoms().intern_of_dihedral(d.a, d.b, d.c, d.d);
            if (di < coords.size() && coords.atoms(di).i_to_a() == d.d)
            {
              std::cout << "Setting dihedral " << di << " to " << d.value << "\n";
              coords.set_dih(di, d.value, true, true);
            }
          }
          coords.to_xyz();
          coords.e();
          std::cout << "Structure " << ++i << lineend;
          coords.e_tostream_short(std::cout);
          outputstream << coords;
        }
        break;
      }
    case config::tasks::GRAD:
      { // gradients 
        coords.e_head_tostream_short(std::cout);
        std::size_t i(0u);
        std::ofstream gstream(coords::output::filename("_GRAD", ".txt").c_str());
        for (auto const & pes : *ci)
        {
          coords.set_xyz(pes.structure.cartesian);
          coords.g();
          std::cout << "Structure " << ++i << lineend;
          coords.e_tostream_short(std::cout);
          coords.energyinterface()->print_G_tinkerlike(gstream);
        }
        break;
      }
    case config::tasks::PROFILE:
      { // gradient profiling 
        coords.e_head_tostream_short(std::cout);
        coords.set_xyz((*ci).PES()[0].structure.cartesian);
        for (std::size_t i(0U); i < Config::get().general.profile_runs; ++i)
        {
          coords.g();
        }
        std::cout << "Energy " << lineend;
        coords.e_tostream_short(std::cout);
        break;
      }
    case config::tasks::LOCOPT:
      { // local optimization
        coords.e_head_tostream_short(std::cout);
        auto lo_structure_fn = coords::output::filename("_LOCOPT");
        std::ofstream locoptstream(lo_structure_fn, std::ios_base::out);
        if (!locoptstream) throw std::runtime_error("Cannot open '" + lo_structure_fn + "' for LOCOPT structures.");
        auto lo_energies_fn = coords::output::filename("_LOCOPT", ".txt");
        std::ofstream loclogstream(lo_energies_fn, std::ios_base::out );
        if (!loclogstream) throw std::runtime_error("Cannot open '" + lo_structure_fn + "' for LOCOPT energies.");
        std::size_t i(0U);
        for (auto const & pes : *ci)
        {
          using namespace std::chrono;
          auto start = high_resolution_clock::now();
          coords.set_xyz(pes.structure.cartesian);
          coords.e();
          std::cout << "Initial: " << ++i << lineend;
          coords.e_tostream_short(std::cout);
          loclogstream << i << ' ' << coords.pes().energy << ' ';
          coords.o();
          auto tim = duration_cast<duration<double>>
            (high_resolution_clock::now() - start);
          loclogstream << coords.pes().energy << ' ' << tim.count() << '\n';
          std::cout << "Post-Opt: " << i << "(" << tim.count() << " s)" << '\n';
          coords.e_tostream_short(std::cout);
          locoptstream << coords;
        }
        break;
      }
    case config::tasks::CENTER:
      { // Center
      std::ofstream centerstream(coords::output::filename("_CENTER").c_str(), std::ios_base::out);
        for (auto const & pes : *ci)
        {
          coords.set_xyz(pes.structure.cartesian);
          coords.move_all_by(-coords.center_of_mass(), true);
          coords::float_type maxdist = 0.0;
          for (auto const & p : coords.xyz())
          {
            coords::float_type const l = geometric_length(p);
            maxdist = l > maxdist ? l : maxdist;
          }
          std::cout << "Maximum distance to center is " << maxdist << lineend;
          centerstream << coords;
        }
        break;
      }
    case config::tasks::TS:
      { // Gradient only tabu search
        std::cout << Config::get().coords.equals;
        std::cout << "-------------------------------------------------" << lineend;
        std::cout << Config::get().optimization.global;
        std::cout << "-------------------------------------------------" << lineend;
        if (Config::get().optimization.global.pre_optimize)
        {
          startopt::apply(coords, ci->PES());
        }
        optimization::global::optimizers::tabuSearch gots(coords, ci->PES());
        gots.run(Config::get().optimization.global.iterations, true);
        gots.write_range("_TS");
        break;
      }
    case config::tasks::MC:
      { // MonteCarlo
        std::cout << Config::get().coords.equals;
        std::cout << "-------------------------------------------------" << lineend;
        std::cout << Config::get().optimization.global;
        std::cout << "-------------------------------------------------" << lineend;
        if (Config::get().optimization.global.pre_optimize)
        {
          startopt::apply(coords, ci->PES());
        }
        optimization::global::optimizers::monteCarlo mc(coords, ci->PES());
        mc.run(Config::get().optimization.global.iterations, true);
        mc.write_range("_MCM");
        break;
      }
    case config::tasks::GRID:
      { // Grid Search
        std::cout << Config::get().coords.equals;
        std::cout << "-------------------------------------------------" << lineend;
        std::cout << Config::get().optimization.global;
        std::cout << "-------------------------------------------------" << lineend;
        optimization::global::optimizers::main_grid mc(coords, ci->PES(), 
          Config::get().optimization.global.grid.main_delta);
        mc.run(Config::get().optimization.global.iterations, true);
        mc.write_range("_GRID");
        break;
      }
    case config::tasks::INTERNAL:
      { // Internal
        for (auto const & pes : *ci)
        {
          coords.set_xyz(pes.structure.cartesian);
          coords.to_internal();
          std::cout << coords::output::formats::zmatrix(coords);
          std::size_t const TX(coords.atoms().mains().size());
          for (std::size_t i(0U); i < TX; ++i)
          {
            std::size_t const j(coords.atoms().intern_of_main_idihedral(i));
            std::size_t const bound_intern(coords.atoms(j).ibond());
            std::size_t const angle_intern(coords.atoms(j).iangle());
            std::cout << "Main " << i << " along " << coords.atoms(bound_intern).i_to_a();
            std::cout << " and " << coords.atoms(angle_intern).i_to_a();
            std::cout << " : " << coords.main(i) << lineend;
          }
          for (auto const & e : coords.main()) std::cout << e << lineend;
        }
        break;
      }
    case config::tasks::DIMER:
      { // Dimer
        coords.e_head_tostream_short(std::cout);
        std::size_t i(0U);
        std::ofstream dimerstream(coords::output::filename("_DIMERTRANS").c_str(), std::ios_base::out);
        for (auto const & pes : *ci)
        {
          coords.set_xyz(pes.structure.cartesian);
          coords.o();
          dimerstream << coords;
          std::cout << "Pre-Dimer Minimum" << ++i << lineend;
          coords.e_tostream_short(std::cout);
          coords.dimermethod_dihedral();
          dimerstream << coords;
          std::cout << "Dimer Transition" << i << lineend;
          coords.e_tostream_short(std::cout);
          coords.o();
          dimerstream << coords;
          std::cout << "Post-Dimer Minimum" << i << lineend;
          coords.e_tostream_short(std::cout);
        }
        break;
      }
    case config::tasks::MD:
      { // Molecular Dynamics
        if (Config::get().md.pre_optimize) coords.o();
        md::simulation mdObject(coords);
        mdObject.run();
        break;
      }
    case config::tasks::FEP:
      { // Free energy perturbation
        md::simulation mdObject(coords);
        mdObject.fepinit();
        mdObject.run();
        mdObject.feprun();
        break;
      }
    case config::tasks::UMBRELLA:
      {
        Config::set().md.umbrella = true;
        md::simulation mdObject(coords);
        mdObject.umbrella_run();
        break;
      }
    case config::tasks::STARTOPT:
      { // Preoptimization
        //std::cout << "PreApply.\n";
        startopt::apply(coords, ci->PES());
        //std::cout << "PostApply.\n";
        std::ofstream gstream(coords::output::filename("_SO").c_str());
        for (auto const & pes : ci->PES())
        {
          //std::cout << "PreSet.\n";
          coords.set_pes(pes,true);
          //std::cout << "PostSet.\n";
          gstream << coords;
        }
        break;
      }
    case config::tasks::GOSOL:
      { // Combined Solvation + Global Optimization
        std::cout << Config::get().startopt.solvadd;
        std::cout << "-------------------------------------------------" << lineend;
        std::cout << Config::get().coords.equals;
        std::cout << "-------------------------------------------------" << lineend;
        std::cout << Config::get().optimization.global;
        std::cout << "-------------------------------------------------" << lineend;
        startopt::preoptimizers::GOSol sopt(coords, ci->PES());
        sopt.run(Config::get().startopt.solvadd.maxNumWater);
        break;
      }
    case config::tasks::NEB:
      {
        std::ptrdiff_t counter = 0;
        coords::Coordinates const coord_obj(coords);
        for (auto const & pes : *ci)
        {
          coords.set_xyz(pes.structure.cartesian);
          coords.mult_struc_counter++;
          neb nobj(&coords);
          nobj.preprocess(counter, Config::get().neb.IMAGES);
        }
        break;
      }
    case config::tasks::PATHOPT:
      {
        std::ptrdiff_t counter = 0;
        coords::Coordinates const coord_obj(coords);
        for (auto const & pes : *ci)
        {
          coords.set_xyz(pes.structure.cartesian);
          coords.mult_struc_counter++;
          neb nobj(&coords);
          nobj.preprocess(counter, Config::get().neb.IMAGES);
          pathx Pobj(&nobj, &coords);
          Pobj.pathx_ini();
        }
        break;
      }
    case config::tasks::PATHSAMPLING:
      {
        coords::Coordinates const coord_obj(coords);
        path_perp path_perpobj(&coords);
        path_perpobj.pathx_ini();
			  break;
		  }
	  case config::tasks::REACTIONCOORDINATE:
		  {
			  coords::Coordinates coords2(coords),coords3(coords);
			  reaccoords reac_obj(&coords,&coords2);
			  std::fstream RMSD("RMSD_REAC.dat",std::ios::app);

			  coords3.set_xyz((*ci).PES()[0].structure.cartesian);

			   for (auto const & pes : *ci)
			   {
				  coords.mult_struc_counter++;
				  coords.set_xyz(pes.structure.cartesian);
				  reac_obj.rmsd_align(coords3);
				  coords.set_xyz(reac_obj.rmsd_align(coords3));
				  reac_obj.bonds();
				  reac_obj.angles();
			   }
			 
			   for (auto const & pes : *ci)
			   {
				  coords2.mult_struc_counter++;
				  coords2.set_xyz(pes.structure.cartesian);
          RMSD << "RMSD : " << scon::root_mean_square_deviation(coords2.xyz(), coords.xyz()) << lineend;

			   }
			   
			  reac_obj.bonds_alteration();
			  reac_obj.angles_alteration();

			  break;
		  }
    case config::tasks::ALIGN:
      {
        /*
         * THIS TASK ALIGNES A SIMULATION TRAJECTORY
         *
         * This task will perform a translational- and rotational fit of the conformations
         * obtained from a molecular simulation according to Kabsch's method. 
         * Furthermore, molecular distance meassures may be computed afterwards
         *
         * Now even fancier through OpenMP magic
         */
        alignment(ci, coords);
        std::cout << "Everything is done. Have a nice day." << std::endl;
        break;
      }
    case config::tasks::PCAgen:
      {
        /**
         * THIS TASK PERFORMS PRINCIPAL COMPONENT ANALYSIS ON A SIMULATION TRAJECTORY
         *
         * This task will perform a principal component analysis (PCA) on a molecular simulation
         * trajectory. Prior translational- and rotational fit of the conformations
         * obtained is possible. Options can be specified in the INPUTFILE.
         *
         * Further processing can be done via PCAproc - Task
         */
        pca_gen(ci, coords);
        std::cout << "Everything is done. Have a nice day." << std::endl;
        break;
      }
    case config::tasks::PCAproc:
      {
        /**
         * THIS TASK PERFORMS Processing of previously obtained PRINCIPAL COMPONENTs
         * To be precise, it will write out the structures coresponding to user specified PC-Ranges.
         * see also: Task PCAgen
         */
        pca_proc(ci, coords);
        std::cout << "Everything is done. Have a nice day." << std::endl;
        break;
      }
    case config::tasks::ENTROPY:
      {
        /**
         * THIS TASK PERFORMS CONFIGURATIONAL ENTROPY CALCULATIONS ON A SIMULATION TRAJECTORY
         *
         * This task will perform verious configurational or conformational entropy calculations
         * on a molecular simualtion trajectory. Prior translational- and rotational fit of the
         * conformations obtained is possible. Options can be specified in the INPUTFILE
         *
         */
        entropy(ci, coords);
        std::cout << "Everything is done. Have a nice day." << std::endl;
        break;
      }
    default:
      {
      
      }

    }
    // stop and print task and execution time
    std::cout << lineend << "Task " << config::task_strings[Config::get().general.task];
    std::cout << " took " << task_timer << " to complete." << lineend;
    std::cout << "Execution of " << config::Programname << " (" << config::Version << ")";
    std::cout << " ended after " << exec_timer << lineend;

#ifndef CAST_DEBUG_DROP_EXCEPTIONS
  }
#if !defined(COMPILEX64)
  catch (std::bad_alloc &)
  {
    std::cout << "Memory allocation failure. CAST probably ran out of memory. Try using 64bit compiled " << config::Programname << " instead." << lineend;
  }
#else
  catch (std::bad_alloc &)
  {
    std::cout << "Memory allocation failure. Input structure probably too large." << lineend;
  }
#endif
  catch (std::exception & e)
  {
    std::cout << "An exception occured. The execution of " << config::Programname << " failed. " << lineend;
    std::cout << "Error: " << e.what() << lineend;
  }
#endif
  return 0;
}
