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
#include <deque> //For PCAproc
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
      using namespace matop;
      using namespace matop::align;

      coords::Coordinates coordsReferenceStructure(coords), coordsTemporaryStructure(coords);
      auto temporaryPESpoint = ci->PES()[Config::get().alignment.reference_frame_num].structure.cartesian;

      //Alignment to external reference frame (different file)
      if (!Config::get().alignment.align_external_file.empty())
      {
        std::unique_ptr<coords::input::format> externalReferenceStructurePtr(coords::input::new_format());
        coords::Coordinates externalReferenceStructure(externalReferenceStructurePtr->read(Config::get().alignment.align_external_file));
        if (Config::get().alignment.reference_frame_num >= externalReferenceStructurePtr->PES().size())
        {
          throw std::out_of_range("Requested reference frame number not within reference structure ensemble.");
        }
        temporaryPESpoint = externalReferenceStructurePtr->PES()[Config::get().alignment.reference_frame_num].structure.cartesian;
      }
      //Constructs two coordinate objects and sets reference frame according to INPUTFILE
      coordsReferenceStructure.set_xyz(temporaryPESpoint);
      Matrix_Class matrixReferenceStructure = matop::transfer_to_matr(coordsReferenceStructure);

      //Construct and Allocate arrays for stringoutput (necessary for OpenMP)
      double mean_value = 0;
      std::string *hold_str, *hold_coords_str;
      hold_str = new std::string[ci->size()];
      hold_coords_str = new std::string[ci->size()];

      //Perform translational alignment for reference frame
      if (Config::get().alignment.traj_align_bool)
      {
        centerOfMassAlignment(coordsReferenceStructure);
      }

#ifdef _OPENMP
      #pragma omp parallel for firstprivate(coords, coordsReferenceStructure, coordsTemporaryStructure, matrixReferenceStructure) reduction(+:mean_value) shared(hold_coords_str, hold_str)
#endif
      for (int i = 0; i < (int) ci->size(); i++)
      {
        if (i != Config::get().alignment.reference_frame_num)
        {
          auto temporaryPESpoint2 = ci->PES()[i].structure.cartesian;
          coordsTemporaryStructure.set_xyz(temporaryPESpoint2);
          //Create temporary objects for current frame

          if (Config::get().alignment.traj_align_bool)
          {
            kabschAlignment(coordsTemporaryStructure, coordsReferenceStructure, true);
          }

          if (Config::get().alignment.traj_print_bool)
          {
            if (Config::get().alignment.dist_unit == 0)
            //RMSD 
            {
              std::stringstream temporaryStringstream;
              double currentRootMeanSquareDevaition = root_mean_square_deviation(coordsTemporaryStructure.xyz(), coordsReferenceStructure.xyz());
              temporaryStringstream << std::setw(13) << i << " ";
              temporaryStringstream << std::setw(13) << currentRootMeanSquareDevaition << "\n";
              mean_value += currentRootMeanSquareDevaition;
              hold_str[i] = temporaryStringstream.str();
            }
            else if (Config::get().alignment.dist_unit == 1)
            //dRMSD
            {
              std::stringstream temporaryStringstream;
              temporaryStringstream << i << " ";
              double value = (double) drmsd_calc(coordsTemporaryStructure, coordsReferenceStructure);
              temporaryStringstream << std::setw(13) << value << "\n";
              mean_value += value;
              hold_str[i] = temporaryStringstream.str();
            }
            else if (Config::get().alignment.dist_unit == 2)
            //Holm&Sander Distance
            {
              std::stringstream temporaryStringstream;
              double value = (double) holmsander_calc(coordsTemporaryStructure, coordsReferenceStructure, Config::get().alignment.holm_sand_r0);
              temporaryStringstream << std::setw(13) << i << " " << value << "\n";
              mean_value += value;
              hold_str[i] = temporaryStringstream.str();
            }
          }
          //Molecular distance measure calculation

          std::stringstream hold_coords;
          hold_coords << coordsTemporaryStructure;
          hold_coords_str[i] = hold_coords.str();
          //Formatted string-output
        }

        else if ( i == Config::get().alignment.reference_frame_num)
        {
          std::stringstream hold_coords;
          hold_coords << coordsReferenceStructure;
          hold_coords_str[i] = hold_coords.str();
          //Formatted string-output (first to array because of OpenMP parallelization)
        }
      }

      std::ofstream distance(coords::output::filename("_distances").c_str(), std::ios::app);
      std::ofstream outputstream(coords::output::filename("_aligned").c_str(), std::ios::app);
      for (size_t i = 0; i < ci->size(); i++)
      {
        if (Config::get().alignment.traj_print_bool)
        {
          distance << hold_str[i];
        }
        if (Config::get().alignment.traj_align_bool)
        {
          outputstream << hold_coords_str[i];
        }
      }
      distance << "\n";
      distance << "Mean value: " << (mean_value / (ci->size() - 1)) << "\n";
      //Formatted string-output

      delete[] hold_str;
      delete[] hold_coords_str;
      //Cleaning Up

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
        using namespace matop;
        using namespace matop::pca;
        
        Matrix_Class pca_modes, eigenvalues, eigenvectors, matrix_aligned;

        if (!Config::get().PCA.pca_read_modes)
        {
          coords::Coordinates coords_ref(coords);
          auto holder = ci->PES()[Config::get().PCA.pca_ref_frame_num].structure.cartesian;
          coords_ref.set_xyz(holder);
          //Constructs two coordinate objects and sets reference frame according to INPUTFILE

          //Perform translational alignment for reference frame
          if (Config::get().PCA.pca_alignment)
          {
            align::centerOfMassAlignment(coords_ref);
          }

          //bool has_it_started = false;
          const size_t FRAME_SIZE = (size_t) ci->size();
          //Initializing some stuff

          for (size_t i = 0; i < FRAME_SIZE; i++)
          {
            if ((Config::get().PCA.pca_start_frame_num <= i) && ((Config::get().PCA.pca_start_frame_num % Config::get().PCA.pca_offset) == (i % Config::get().PCA.pca_offset)))
            {
              auto holder2 = ci->PES()[i].structure.cartesian;
              coords.set_xyz(holder2);
              //Initializing current frame

              if (Config::get().PCA.pca_alignment && !Config::get().PCA.pca_use_internal)
              {
                align::centerOfMassAlignment(coords); //Alignes center of mass
                align::kabschAlignment(coords, coords_ref); //Rotates
              }
              //Translational and rotational alignment

              if (Config::get().PCA.pca_use_internal)
              {
                coords.to_internal();
              }
              //Conversion to internal coordinates if desired

              if (Config::get().PCA.pca_start_frame_num < i)
              {
                if (Config::get().PCA.pca_use_internal)
                {
                  matrix_aligned.append_bottom(transformToOneline(coords, Config::get().PCA.pca_internal_dih, true));
                }
                else
                {
                  //Works for full and truncated PCA
                  matrix_aligned.append_bottom(transformToOneline(coords, Config::get().PCA.pca_trunc_atoms_num, false));
                }
              }
              else if (Config::get().PCA.pca_start_frame_num == i)
              {
                if (Config::get().PCA.pca_use_internal)
                {
                  //First a little check if the user-specified values are reasonable
                  matrix_aligned = transformToOneline(coords, Config::get().PCA.pca_internal_dih, true);
                }
                else
                {
                  matrix_aligned = transformToOneline(coords, Config::get().PCA.pca_trunc_atoms_num, false);
                }
              }
              //Building one huge [frames] x [coordinates] matrix by appending for every frame
            }
          }

          transpose(matrix_aligned);
          // NECESSARY because of implementation details, don't worry about it for now
          // Matrix is now [coordinates] x [frames] !!! !!!
          // THIS IS THE CONVENTION FOR USAGE, STICK TO IT!

          if (!Config::get().PCA.pca_use_internal)
          {
            coords_ref.set_xyz(holder);
            if (!Config::get().PCA.pca_trunc_atoms_bool)
            {
              massweight(matrix_aligned, coords_ref, false);
            }
            else
            {
              massweight(matrix_aligned, coords_ref, false, Config::get().PCA.pca_trunc_atoms_num);
            }
          }
          //Mass-weightening coordinates if cartesians are used
          prepare_pca(matrix_aligned, eigenvalues, eigenvectors);
        }

        if(Config::get().PCA.pca_read_vectors)
        {
          std::string iAmNotImportant_YouMayDiscardMe;
          readEigenvectorsAndModes(eigenvectors, pca_modes, iAmNotImportant_YouMayDiscardMe);
        }

        if (!Config::get().PCA.pca_read_modes)
        {
          pca_modes = transposed(eigenvectors) * matrix_aligned;
        }

        ///////////////////////////////////////

        std::string additionalInformation;
        if (Config::get().PCA.pca_use_internal)
        {
          additionalInformation += "int ";
          for (size_t i = 0u; i < Config::get().PCA.pca_internal_dih.size(); i++)
          {
            additionalInformation += "d" + std::to_string(Config::get().PCA.pca_internal_dih[i]) + " ";
          }
        }
        else
        {
          additionalInformation += "car ";
          for (size_t i = 0u; Config::get().PCA.pca_trunc_atoms_bool && i < Config::get().PCA.pca_trunc_atoms_num.size(); i++)
          {
            additionalInformation += std::to_string(Config::get().PCA.pca_trunc_atoms_num[i]) + " ";
          }
        }

        ///////////////////////////////////////

        if (Config::get().PCA.pca_read_vectors)
        {
          output_pca_modes(eigenvalues, eigenvectors, pca_modes, "pca_modes_new.dat", additionalInformation);
        }
        else
        {
          output_pca_modes(eigenvalues, eigenvectors, pca_modes, "pca_modes.dat", additionalInformation);
        }

        if (Config::get().PCA.pca_print_probability_density)
        {
          output_probability_density(pca_modes);
        }

        ///////////////////////////////////////

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

        //TODO: The first structure on output from this task "seems" to be off. Revie plz

        using namespace matop;
        using namespace matop::pca;
        Matrix_Class eigenvectors, trajectory;
        std::vector<size_t> structuresToBeWrittenToFile;
        std::string additionalInformation;
        readEigenvectorsAndModes(eigenvectors, trajectory, additionalInformation);
        if (Config::get().PCA.proc_desired_start.size() > trajectory.rows() || Config::get().PCA.proc_desired_stop.size() > trajectory.rows())
        {
          std::cerr << "Desired PCA-Ranges have higher dimensionality then modes. Omitting the last values.\n";
        }

        for (size_t j = 0u; j < trajectory.cols(); j++)
        {
          bool isStructureInRange = true;
          for (size_t i = 0u; i < trajectory.rows() && i < std::max(Config::get().PCA.proc_desired_stop.size(), Config::get().PCA.proc_desired_start.size()); i++)

          {
            if (i < Config::get().PCA.proc_desired_start.size())
            {
              if (trajectory(i, j) < Config::get().PCA.proc_desired_start[i])
              {
                isStructureInRange = false;
              }
            }
            if (i < Config::get().PCA.proc_desired_stop.size())
            {
              if (trajectory(i, j) > Config::get().PCA.proc_desired_stop[i])
              {
                isStructureInRange = false;
              }
            }
          }
          if (isStructureInRange)
          {
            structuresToBeWrittenToFile.push_back(j);
          }
        }

        if (Config::get().general.verbosity >= 3u) std::cout << "Found " << structuresToBeWrittenToFile.size() << " structures in desired range.\n";

        //Undoing PCA
        trajectory = eigenvectors * trajectory;

        std::ofstream outstream(coords::output::filename("_pca_selection").c_str(), std::ios::app);

        // Case: Cartesian Coordinates.
        if (additionalInformation.substr(0, 3) == "car")
        {
          //Additional Information Processing -> Read "DOFS that were used" from file. Put their identifying numbers in vector.
          stringstream ss(additionalInformation.substr(4));
          std::string buffer;
          std::vector<size_t> tokens;
          std::deque<bool> alreadyFoundStructures(ci->size(), false);

          while (ss >> buffer) tokens.push_back((size_t)std::stoi(buffer));

          if (tokens.size() != 0u)
          {
            undoMassweight(trajectory, coords, false, tokens);
            for (size_t i = 0u; i < structuresToBeWrittenToFile.size(); i++)
            {
              Matrix_Class out_mat(3, trajectory.rows() / 3u);
              for (size_t j = 0u; j < trajectory.rows(); j = j + 3)
              {
                out_mat(0, j / 3u) = trajectory(j, structuresToBeWrittenToFile[i]);
                out_mat(1, j / 3u) = trajectory(j + 1u, structuresToBeWrittenToFile[i]);
                out_mat(2, j / 3u) = trajectory(j + 2u, structuresToBeWrittenToFile[i]);
              }
              coords::Coordinates current(coords);

              // For every partial (truncated) structure that is inside the user-defined 
              // range regarding its PCA-Modes,
              // we now search the matching full structure in the input trajectory.
              bool structureFound = false;
              auto structureCartesian = ci->PES()[0].structure.cartesian;
              int structureNumber = -1;

              for (size_t k = 0u; k < ci->size() && !structureFound; k++)
              {
                if (alreadyFoundStructures[k]) continue;

                structureCartesian = ci->PES()[k].structure.cartesian; //Current structure
                structureFound = true;
                structureNumber = (int) k;
                // Remeber, we are now iterating over certain atoms (those to which
                // the PCA was truncated.
                for (size_t l = 0u; l < tokens.size(); l++)
                {
                  // If abs() of diff of every coordinate is smaller than 0.5% of coordinate (or, if this value
                  // is very small, the arbitrary cutoff 2e-4), consider it a match.
                  // However, we look for "not-matching" and break the loop. If everything matches, we continue. 
                  // Thats why we negate the criterion in the if clause (!)
                  float_type xCompare = 0.005 * std::max(std::abs(structureCartesian[tokens[l]].x()), 2e-4);
                  float_type yCompare = 0.005 * std::max(std::abs(structureCartesian[tokens[l]].y()), 2e-4);
                  float_type zCompare = 0.005 * std::max(std::abs(structureCartesian[tokens[l]].z()), 2e-4);

                  if (!(std::abs(out_mat(0, l) - structureCartesian[tokens[l]].x()) <= xCompare &&
                    std::abs(out_mat(1, l) - structureCartesian[tokens[l]].y()) <= yCompare &&
                    std::abs(out_mat(2, l) - structureCartesian[tokens[l]].z()) <= zCompare))
                  {
                    structureFound = false;
                    break;
                  }
                }
              }
              if (structureFound)
              {
                // Horray, we found it, now write it out!
                alreadyFoundStructures[structureNumber] = true;
                current.set_xyz(structureCartesian);
                current.to_internal();
                outstream << current;
              }
              else
              {
                std::cerr << "Could not find structure restored from PCA-Modes in ensemble of structures from original coordinates.\n";
                std::cerr << "This means that there was no provided structure with a deviance of less than 0.5% to the current restored structure.\n\n";
              }
            }
          }
          else
          {
            // Here, we merely restore the coordinates from the PCA-modes
            // since no trucnation took plage, and write them out.
            undoMassweight(trajectory, coords, false);
            for (size_t i = 0u; i < structuresToBeWrittenToFile.size(); i++)
            {
              Matrix_Class out_mat(3, trajectory.rows() / 3u);
              for (size_t j = 0u; j < trajectory.rows(); j = j + 3)
              {
                out_mat(0, j / 3u) = trajectory(j, structuresToBeWrittenToFile[i]);
                out_mat(1, j / 3u) = trajectory(j + 1u, structuresToBeWrittenToFile[i]);
                out_mat(2, j / 3u) = trajectory(j + 2u, structuresToBeWrittenToFile[i]);
              }
              coords::Coordinates out(coords);
              out.set_xyz(transfer_to_3DRepressentation(out_mat));
              out.to_internal();
              outstream << out;
            }
          }
        }

        // Case: Internal Coordiantes
        else if (additionalInformation.substr(0, 3) == "int")
        {
          // [0]: bond distance tokens, [1]: bond angle tokens [2]: dihedrals tokens.
          // Here we keep track of which DOFs are to be considered
          std::deque<bool> tokens(coords.atoms().size(), false);
          std::deque<bool> alreadyFoundStructures(ci->size(), false);

          //Additional Information Processing -> Read "DOFS that were used" from file. Store in "tokens".
          stringstream ss(additionalInformation.substr(4));
          std::string buffer;
          // Read the additional information
          while (ss >> buffer)
          {
            if (buffer.substr(0, 1) == "d") tokens[(size_t) std::stoi(buffer.substr(1))] = true;
          }

          // This is gonna be complicated. Im sorry.
          // Iterate over structures chosen from PCA-Ensemble
          for (size_t i = 0u; i < structuresToBeWrittenToFile.size(); i++)
          {
            bool structureFound = false;
            //auto structureCartesian = ci->PES()[0].structure.cartesian;
            int structureNumber = -1;

            //Iterate over input strucutures -> find matching structure
            for (size_t k = 0u; k < ci->size() && !structureFound; k++)
            {
              if (alreadyFoundStructures[k]) continue;
              coords.set_internal(ci->PES()[k].structure.intern);
              size_t quicksearch = 0u;
              structureFound = true;
              structureNumber = (int)k;
              //Iterating over atoms, see if they all match
              for (size_t j = 0u; j < tokens.size(); j++)
              {
                // If abs() of diff of every coordinate is smaller than 0.1% of coordinate, consider it a match.
                // However, we look for "not-matching" and break the loop. If everything matches, we continue. 
                // Thats why we negate the criterion in the if clause (!)
                if (tokens[j] == true)
                {
                  auto compareFromPCA1 = trajectory(quicksearch, structuresToBeWrittenToFile[i]);
                  auto compareFromTrajectory1 = std::cos(coords.intern(j).azimuth().radians());
                  auto compareFromPCA2 = trajectory(quicksearch + 1u, structuresToBeWrittenToFile[i]);
                  auto compareFromTrajectory2 = std::sin(coords.intern(j).azimuth().radians());
                  bool found1 = std::abs(compareFromTrajectory1 - compareFromPCA1) <= 0.001 * std::abs(compareFromPCA1) \
                    || std::abs(compareFromTrajectory1 - compareFromPCA1) < 0.0000001;
                  bool found2 = std::abs(compareFromTrajectory2 - compareFromPCA2) <= 0.001 * std::abs(compareFromPCA2) \
                    || std::abs(compareFromTrajectory2 - compareFromPCA2) < 0.0000001;
                  //If structures did not match
                  if (!(found1 && found2))
                  {
                    structureFound = false;
                    break;
                  }
                  quicksearch += 2u;
                }
              }
              //If match was found, write out.
              if (structureFound)
              {
                alreadyFoundStructures[structureNumber] = true;
                coords.set_pes( ci->PES()[k]);
                outstream << coords;
              }
            }
            if (!structureFound)
            {
              std::cerr << "Could not find structure restored from PCA-Modes in ensemble of structures from original coordinates.\n";
              std::cerr << "You probably made a mistake somewhere in your INPUTFILE.\nDo not consider structures written out after this message as valid.\n";
            } 
          }
        }
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
        using namespace matop;
        using namespace matop::entropy;

        //Initialize the reference frame (for alignment etc)
        coords::Coordinates coords_ref(coords);
        auto holder = (*ci).PES()[Config::get().entropy.entropy_ref_frame_num].structure.cartesian;
        coords_ref.set_xyz(holder);

        //Translational alignment of the reference frame
        if (Config::get().entropy.entropy_alignment)
        {
          align::centerOfMassAlignment(coords_ref);
        }

        Matrix_Class matrix_aligned;
        const size_t FRAME_SIZE = int(ci->size());
        if (Config::get().entropy.entropy_alignment && Config::get().entropy.entropy_use_internal)
        {
          std::cerr << "Alignment is (in this case) redundant since internal coordinates are used. Alignment is skipped. Check your INPUTFILE please.\n";
          std::cerr << "Continuing anyway...";
        }
        //Initializing and checking...

        for (size_t i = 0; i < FRAME_SIZE; i++)
        {
          // Meet-your-Maker-Note: Keep it like this with the seeminlgy stupid "is it started", please.
          if (Config::get().entropy.entropy_start_frame_num <= i && (Config::get().entropy.entropy_start_frame_num % Config::get().entropy.entropy_offset == i % Config::get().entropy.entropy_offset))
          {
            auto holder2 = ci->PES()[i].structure.cartesian;
            coords.set_xyz(holder2);
            //Initializing current frame

            if (Config::get().entropy.entropy_alignment && !Config::get().entropy.entropy_use_internal)
            {
              align::centerOfMassAlignment(coords); //Alignes center of mass
              align::kabschAlignment(coords, coords_ref); //Rotates
            }
            //Translational and rotational alignment

            if (Config::get().entropy.entropy_use_internal)
            {
              coords.to_internal();
            }
            //Conversion to internal coordinates if desired

            if (Config::get().entropy.entropy_start_frame_num < i)
            {
              if (Config::get().entropy.entropy_use_internal)
              {
                matrix_aligned.append_bottom(transformToOneline(coords, Config::get().entropy.entropy_internal_dih, true));
              }
              else
              {
                matrix_aligned.append_bottom(transformToOneline(coords, Config::get().entropy.entropy_trunc_atoms_num, false));
              }
            }
            else if (Config::get().entropy.entropy_start_frame_num == i)
            {
              if (Config::get().entropy.entropy_use_internal)
              {
                matrix_aligned = transformToOneline(coords, Config::get().entropy.entropy_internal_dih, true);
              }
              else
              {
                matrix_aligned = transformToOneline(coords, Config::get().entropy.entropy_trunc_atoms_num, false);
              }
            }
            //Building one huge [coordinates] x [frames] matrix by appending for every frame

          }
        }
        transpose(matrix_aligned);
        //NECESSARY because of implementation details, don't worry about it for now; rows are DOFs, columns are frames FROM HERE ON!

        if (!Config::get().entropy.entropy_use_internal)
        {
          massweight(matrix_aligned, coords_ref, true);
        }
        //Mass-weightening cartesian coordinates

        for (size_t u = 0u; u < Config::get().entropy.entropy_method.size(); u++)
        {
          Matrix_Class& workobj = matrix_aligned;
          int m = (int) Config::get().entropy.entropy_method[u];
          if (m == 1 || m == 0)
          {
            /*double entropy_value = */karplus_wrapper(workobj);
          }
          if (m == 2)
          {
            /*double entropy_value = */knapp_m_wrapper(workobj);
          }
          if (m == 3 || m == 0)
          {
            /*double entropy_value = */knapp_wrapper(workobj);
          }
          if (m == 4 || m == 0)
          {
            /*double entropy_value = */hnizdo_wrapper(workobj);
          }
          if (m == 5 || m == 0)
          {
            /*double entropy_value = */hnizdo_m_wrapper(workobj);
          }
          if (m == 6 || m == 0)
          {
            /*double entropy_value = */schlitter_wrapper(workobj);
          }
        }
        //Perform desired calculations

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
