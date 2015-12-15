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
#include "matop.h"
//#include "gbsa.h"
//#include <omp.h>


#if defined(_MSC_VER)
#include <process.h>
#define pid_func _getpid
#else 
#include <unistd.h>
#define pid_func getpid
#endif

//#define CAST_DEBUG_DROP_EXCEPTIONS


void apply_startopt(coords::Coordinates & c, coords::Ensemble_PES & e)
{
  startopt::Preoptimizer * optimizer(nullptr);
  std::size_t multi(Config::get().startopt.number_of_structures / e.size());
  std::cout << Config::get().startopt;
  std::cout << "-------------------------------------------------" << lineend;
  if (Config::get().startopt.type == config::startopt::types::T::SOLVADD)
  {
    optimizer = new startopt::preoptimizers::Solvadd(c);
  }
  else if (Config::get().startopt.type == config::startopt::types::T::RINGSEARCH ||
           Config::get().startopt.type == config::startopt::types::T::RINGSEARCH_SOLVADD)
  {
    optimizer = new startopt::preoptimizers::R_evolution(c);
    multi = Config::get().startopt.ringsearch.population;
  }
  //std::cout << c;
  if (optimizer != nullptr)
  {
    // Generate structures using initial preoptimizer (ringsearch or solvadd)
    //std::cout << "PreGenerate.\n";
    optimizer->generate(e, (multi > 0u?multi:1u));
    //std::cout << "PostGenerate.\n";
    if (Config::get().startopt.type == config::startopt::types::T::RINGSEARCH_SOLVADD)
    {
      std::size_t sa_multi(Config::get().startopt.number_of_structures / optimizer->PES().size());
      startopt::preoptimizers::Solvadd sa(optimizer->final_coords());
      sa.generate(optimizer->PES(), sa_multi);
      delete optimizer;
      optimizer = &sa;
    }
    c.swap(optimizer->final_coords());
    e = optimizer->PES();
    if (Config::get().startopt.type != config::startopt::types::T::RINGSEARCH_SOLVADD)
      delete optimizer;
  }
}


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

    //std::cout << "-------------------------------------------------" << lineend;
    //std::cout << "Structure:" << lineend;
    //std::cout << "-------------------------------------------------" << lineend;

    //std::cout << coords;

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

      scon::matrix<coords::float_type> test(5, 4);
      test.resize(4u, 4u);
      test.resize(7u, 8u);
      
      break;
    }
    case config::tasks::SP:
      { // singlepoint
        coords.e_head_tostream_short(std::cout);
        std::size_t i(0u);
        //std::ofstream outputstream(coords::output::filename("_SP").c_str(), std::ios_base::out);
        for (auto const & pes : *ci)
        {
          coords.set_xyz(pes.structure.cartesian);
          auto start = std::chrono::high_resolution_clock::now();
          coords.e();
          std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>
            (std::chrono::high_resolution_clock::now() - start).count() << " ns\n";
          std::cout << "Structure " << ++i << lineend;
          coords.e_tostream_short(std::cout);
          //outputstream << coords;
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
        std::ofstream locoptstream(coords::output::filename("_LOCOPT").c_str(), std::ios_base::out);
        std::size_t i(0U);
        for (auto const & pes : *ci)
        {
          coords.set_xyz(pes.structure.cartesian);
          coords.e();
          std::cout << "Initial: " << ++i << lineend;
          coords.e_tostream_short(std::cout);
          coords.o();
          std::cout << "Post-Opt: " << i << lineend;
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
          apply_startopt(coords, ci->PES());
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
          apply_startopt(coords, ci->PES());
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
        coords.e_head_tostream_short(std::cout);
        coords.o();
        Config::set().general.trackstream = &dimerstream;
        coords.dimermethod_dihedral();
        std::cout << "Optimized energy after transition: " << coords.o() << lineend;
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
        apply_startopt(coords, ci->PES());
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
          Pobj.pathx_ini(Config::get().neb.IMAGES);
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
       * obtained from a molecular simulation. Furthermore, molecular distance meassures
       * may be computed afterwards
       *
       * Now even fancier through OpenMP magic
       */
      using namespace matop;
      using namespace matop::align;


      coords::Coordinates coords_ref(coords), coords_temp(coords);
      auto holder = ci->PES()[Config::get().alignment.reference_frame_num].structure.cartesian;
      coords_ref.set_xyz(holder);
      Matrix_Class ref = matop::transfer_to_matr(coords_ref);
      //Constructs two coordinate objects and sets reference frame according to INPUTFILE

      double mean_value = 0;
      std::string *hold_str, *hold_coords_str;
      hold_str = new std::string[ci->size()];
      hold_coords_str = new std::string[ci->size()];
      //Construct and Allocate arrays for stringoutput (necessary for OpenMP)

      if (Config::get().alignment.traj_align_bool)
      {
        matop::align::align_center_of_mass(ref, coords_ref);
        coords_ref.set_xyz(matop::transfer_to_3DRepressentation(ref));
      }
      //Perform translational alignment for reference frame
#ifdef _OPENMP
      #pragma omp parallel for firstprivate(coords, ref) reduction(+:mean_value) shared(hold_coords_str, hold_str)
#endif
      for (int i = 0; i < (int) ci->size(); i++)
      {
        if ((unsigned int) i != Config::get().alignment.reference_frame_num)
        {
          auto holder2 = ci->PES()[i].structure.cartesian;
          coords_temp.set_xyz(holder2);
          Matrix_Class matr_structure = transfer_to_matr(coords_temp);
          //Create temporary objects for current frame

          if (Config::get().alignment.traj_align_bool)
          {
            align_center_of_mass(matr_structure, coords_temp);
            rotate(matr_structure, ref);
            //matr_structure = (matr_structure.rotated(ref));
          }
          coords_temp.set_xyz(transfer_to_3DRepressentation(matr_structure));
          //Alignment taking place

          if (Config::get().alignment.traj_print_bool)
          {
            if (Config::get().alignment.dist_unit == 0)
              //RMSD 
            {
              std::stringstream hold_distance;
              hold_distance << "Frame " << i << "/" << ci->size() << " ";
              hold_distance << "Distance: " << root_mean_square_deviation(coords_temp.xyz(), coords_ref.xyz()) << "\n";
              mean_value += root_mean_square_deviation(coords_temp.xyz(), coords_ref.xyz());
              hold_str[i] = hold_distance.str();
            }
            else if (Config::get().alignment.dist_unit == 1)
              //dRMSD
            {
              std::stringstream hold_distance;
              hold_distance << "Frame " << i << "/" << ci->size() << " ";
              long double value = drmsd_calc(matr_structure, ref);
              value = sqrt(double(ci->size()) * (double(ci->size()) + 1) * value);
              hold_distance << "Distance: " << value << "\n";
              mean_value += value;
              hold_str[i] = hold_distance.str();
            }
            else if (Config::get().alignment.dist_unit == 2)
              //Holm&Sander Distance
            {
              std::stringstream hold_distance;
              long double value = holmsander_calc(matr_structure, ref);
              hold_distance << "Frame " << i << "/" << ci->size() << " " << "Distance: " << value << "\n";
              mean_value += value;
              hold_str[i] = hold_distance.str();
            }
          }
          //Molecular distance measure calculation

          std::stringstream hold_coords;
          coords.set_xyz(transfer_to_3DRepressentation(matr_structure));
          hold_coords << coords;
          hold_coords_str[i] = hold_coords.str();
          //Formatted string-output
        }

        else if ((unsigned int) i == Config::get().alignment.reference_frame_num)
        {
          std::stringstream hold_coords;
          hold_coords << coords_ref;
          hold_coords_str[i] = hold_coords.str();
          //Formatted string-output (first to array because of OpenMP parallelization)
        }
      }

      std::ofstream distance(coords::output::filename("_distances").c_str(), std::ios::app);
      std::ofstream outputstream(coords::output::filename("_aligned").c_str(), std::ios::app);
      for (unsigned int i = 0; i < ci->size(); i++)
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
          Matrix_Class ref = transfer_to_matr(coords_ref);
          //Constructs two coordinate objects and sets reference frame according to INPUTFILE

          if (Config::get().PCA.pca_alignment)
          {
            align::align_center_of_mass(ref, coords);
            coords_ref.set_xyz(transfer_to_3DRepressentation(ref));;
          }
          //Perform translational alignment for reference frame

          bool has_it_started = false;
          const unsigned int FRAME_SIZE = int(ci->size());
          //Initializing some stuff

          for (unsigned int i = 0; i < FRAME_SIZE; i++)
          {
            if ((Config::get().PCA.pca_start_frame_num <= i) && ((Config::get().PCA.pca_start_frame_num % Config::get().PCA.pca_offset) == (i % Config::get().PCA.pca_offset)))
            {
              auto holder2 = ci->PES()[i].structure.cartesian;
              coords.set_xyz(holder2);
              Matrix_Class matr_structure = transfer_to_matr(coords);
              //Initializing current frame

              if (Config::get().PCA.pca_alignment && !Config::get().PCA.pca_use_internal)
              {
                align::align_center_of_mass(matr_structure, coords); //Alignes center of mass
                align::rotate(matr_structure, ref); //Rotates
              }
              //Translational and rotational alignment

              if (Config::get().PCA.pca_use_internal)
              {
                coords.set_xyz(transfer_to_3DRepressentation(matr_structure));
                coords.to_internal();
                matr_structure = transfer_to_matr_internal(coords);
              }
              //Conversion to internal coordinates if desired

              if (Config::get().PCA.pca_start_frame_num < i)
              {
                if (Config::get().PCA.pca_use_internal)
                {
                  matrix_aligned.append_bottom(transform_3n_nf_internal_pca(matr_structure));
                }
                else
                {
                  if (Config::get().PCA.pca_trunc_atoms_bool)
                  {
                    matrix_aligned.append_bottom(transform_3n_nf_trunc_pca(matr_structure));
                  }
                  else
                  {
                    matrix_aligned.append_bottom(transform_3n_nf(matr_structure));
                  }
                }
              }
              else if (Config::get().PCA.pca_start_frame_num == i)
              {
                if (Config::get().PCA.pca_use_internal)
                {
                  matrix_aligned = transform_3n_nf_internal_pca(matr_structure);
                }
                else
                {
                  matrix_aligned = transform_3n_nf(matr_structure);
                  if (Config::get().PCA.pca_trunc_atoms_bool)
                  {
                    matrix_aligned = transform_3n_nf_trunc_pca(matr_structure);
                  }
                }
              }
              //Building one huge [frames] x [coordinates] matrix by appending for every frame
            }
          }

          matrix_aligned.transpose();
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

          prepare_pca(matrix_aligned, eigenvalues, eigenvectors, pca_modes);
        }

        if(Config::get().PCA.pca_read_vectors)
        {
          readEigenvectorsAndModes(eigenvectors, pca_modes);
        }

        if (!Config::get().PCA.pca_read_modes)
        {
          pca_modes = eigenvectors.transposed() * matrix_aligned;
        }

        ///////////////////////////////////////

        if (Config::get().PCA.pca_read_vectors)
        {
          output_pca_modes(eigenvalues, eigenvectors, pca_modes, "pca_modes_new.dat");
        }
        else
        {
          output_pca_modes(eigenvalues, eigenvectors, pca_modes);
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

        using namespace matop;
        using namespace matop::pca;
        Matrix_Class eigenvectors, trajectory;
        std::vector<unsigned int> structuresToBeWrittenToFile;
        readEigenvectorsAndModes(eigenvectors, trajectory);
        if ( Config::get().PCA.proc_desired_start.size() > trajectory.return_rows() || Config::get().PCA.proc_desired_stop.size() > trajectory.return_rows() )
        {
          std::cerr << "Desired PCA-Ranges have higher dimensionality then modes. Omitting the last values.\n";
        }


        for (unsigned int j = 0u; j < trajectory.return_cols(); j++)
        {
          bool isStructureInRange = true;
          for (unsigned int i = 0u; i < trajectory.return_rows() && i < std::max(Config::get().PCA.proc_desired_stop.size(), 
               Config::get().PCA.proc_desired_start.size()); i++)
          {
            if (i < Config::get().PCA.proc_desired_start.size())
            {
              if (trajectory(i,j) < Config::get().PCA.proc_desired_start[i])
              {
                isStructureInRange = false;
              }
            }
            if (i < Config::get().PCA.proc_desired_stop.size())
            {
              if (trajectory(i,j) > Config::get().PCA.proc_desired_stop[i])
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

        trajectory = eigenvectors * trajectory;
        std::ofstream outstream(coords::output::filename("_pca_selection").c_str(), std::ios::app);
        for (unsigned int i = 0u; i < structuresToBeWrittenToFile.size(); i++)
        {
          Matrix_Class out_mat(3, trajectory.return_rows() / 3u);
          for (unsigned int j = 0u; j < trajectory.return_rows(); j = j + 3)
          {
            out_mat(0, j / 3u) = trajectory(j, i);
            out_mat(1, j / 3u) = trajectory(j + 1u, i);
            out_mat(2, j / 3u) = trajectory(j + 2u, i);
          }
          coords::Coordinates out(coords);
          out.set_xyz(transfer_to_3DRepressentation(out_mat));
          outstream << out;
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

        coords::Coordinates coords_ref(coords);
        auto holder = (*ci).PES()[Config::get().entropy.entropy_ref_frame_num].structure.cartesian;
        coords_ref.set_xyz(holder);
        Matrix_Class ref = transfer_to_matr(coords_ref);
        //Initialize the reference frame (for alignment etc)

        if (Config::get().entropy.entropy_alignment)
        {
          align::align_center_of_mass(ref,coords_ref);
        }
        //Translational alignment of the reference frame

        Matrix_Class matrix_aligned;
        const unsigned int FRAME_SIZE = int(ci->size());
        if (Config::get().entropy.entropy_alignment && Config::get().entropy.entropy_use_internal)
        {
          std::cerr << "Alignment is (in this case) redundant when internal coordinates are used. Alignment is skipped. Check your INPUTFILE please.\n";
          std::cerr << "Continuing anyway...";
        }
        //Initializing and checking...

        for (unsigned int i = 0; i < FRAME_SIZE; i++)
        {
          // Meet-your-Maker-Note: Keep it like this with the seeminlgy stupid "is it started", please.
          if (Config::get().entropy.entropy_start_frame_num <= i && (Config::get().entropy.entropy_start_frame_num % Config::get().entropy.entropy_offset == i % Config::get().entropy.entropy_offset))
          {
            auto holder2 = ci->PES()[i].structure.cartesian;
            coords.set_xyz(holder2);
            Matrix_Class matr_structure = transfer_to_matr(coords);
            //Initializing current frame

            if (Config::get().entropy.entropy_alignment && !Config::get().entropy.entropy_use_internal)
            {
              align::align_center_of_mass(matr_structure, coords); //Alignes center of mass
              align::rotate(matr_structure, ref); //Rotates
            }
            //Translational and rotational alignment

            if (Config::get().entropy.entropy_use_internal)
            {
              coords.set_xyz(transfer_to_3DRepressentation(matr_structure));
              coords.to_internal();
              matr_structure = transfer_to_matr_internal(coords);
            }
            //Conversion to internal coordinates if desired

            if (Config::get().entropy.entropy_start_frame_num < i)
            {
              if (Config::get().entropy.entropy_use_internal)
              {
                matrix_aligned.append_bottom(transform_3n_nf_internal_entropy(matr_structure));
              }
              else
              {
                if (Config::get().entropy.entropy_trunc_atoms_bool)
                {
                  matrix_aligned.append_bottom(transform_3n_nf_trunc_entropy(matr_structure));
                }
                else
                {
                  matrix_aligned.append_bottom(transform_3n_nf(matr_structure));
                }
              }
            }
            else if (Config::get().entropy.entropy_start_frame_num == i)
            {
              if (Config::get().entropy.entropy_use_internal)
              {
                matrix_aligned = transform_3n_nf_internal_entropy(matr_structure);
              }
              else
              {
                matrix_aligned = transform_3n_nf(matr_structure);
                if (Config::get().entropy.entropy_trunc_atoms_bool)
                {
                  matrix_aligned = transform_3n_nf_trunc_entropy(matr_structure);
                }
              }
            }
            //Building one huge [coordinates] x [frames] matrix by appending for every frame

          }
        }
        matrix_aligned.transpose();
        //NECESSARY because of implementation details, don't worry about it for now; rows are DOFs, columns are frames FROM HERE ON!

        if (!Config::get().entropy.entropy_use_internal)
        {
          massweight(matrix_aligned, coords_ref, true);
        }
        //Mass-weightening cartesian coordinates

        //std::cout << matrix_aligned;
        for (unsigned int u = 0u; u < Config::get().entropy.entropy_method.size(); u++)
        {
          Matrix_Class& workobj = matrix_aligned;
          int m = Config::get().entropy.entropy_method[u];
          if (m == 1 || m == 0)
          {
            double entropy_value = karplus_wrapper(workobj);
          }
          if (m == 2)
          {
            double entropy_value = knapp_m_wrapper(workobj);
          }
          if (m == 3 || m == 0)
          {
            double entropy_value = knapp_wrapper(workobj);
          }
          if (m == 4 || m == 0)
          {
            double entropy_value = hnizdo_wrapper(workobj);
          }
          if (m == 5 || m == 0)
          {
            double entropy_value = hnizdo_m_wrapper(workobj);
          }
          if (m == 6 || m == 0)
          {
            double entropy_value = schlitter_wrapper(workobj);
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
    throw;
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
