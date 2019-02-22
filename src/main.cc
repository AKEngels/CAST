
//////////   //////////   //////////   ////////// 
//           //      //   //               //
//           //      //   //               //
//           //      //   //////////       //
//           //////////           //       //
//           //      //           //       //
//           //      //           //       //
//////////   //      //   //////////       //
/////conformational analysis and search tool/////

// If we run our tests, we will
// take the testing main from gtest/testing_main.cc
#ifndef GOOGLE_MOCK

//////////////////////////
//                      //
//       L I B S        //
//                      //
//////////////////////////
#ifdef USE_PYTHON
#include <Python.h>
#endif
#include <cstdlib>
#include <fstream>
#include <memory>
#include <omp.h>

//////////////////////////
//                      //
//    H E A D E R S     //
//                      //
//////////////////////////
#include "configuration.h"
#include "coords_io.h"
#include "scon_chrono.h"
#include "helperfunctions.h"
#include "scon_log.h"
#ifdef _MSC_VER
#include "win_inc.h"
#endif
// Task items
#include "startopt_solvadd.h"
#include "startopt_ringsearch.h"
#include "md.h"
#include "optimization_global.h"
#include "pathopt.h"
#include "Path_perp.h"
#include "matop.h" //For ALIGN, PCAgen, ENTROPY, PCAproc
#include "PCA.h"
#include "2DScan.h"
#include "exciton_breakup.h"
#include "Center.h"
#include "Couplings.h"
#include "periodicCutout.h"
#include "replaceMonomers.h"
#include "modify_sk.h"


//////////////////////////
//                      //
//  R A N D O M         //
//  N U M B E R         //
//  G E N E R A T O R   //
//                      //
//////////////////////////
#if defined(_MSC_VER)
#include <process.h>
#define pid_func _getpid
#else 
#include <unistd.h>
#define pid_func getpid
#endif


// Enable this define to drop exceptions
//
// If CAST_DEBUG_DROP_EXCEPTIONS is set,
// exceptions will be dropped (this is good for debugging
// and default in Debug configuration.
// Otherwise CAST will crash and print the string
// attached to the exception. This is default behaviour 
// in release mode.
//
//#define CAST_DEBUG_DROP_EXCEPTIONS

int main(int argc, char **argv)
{
#ifdef USE_PYTHON   
  Py_Initialize();  // if python enabled: initialise here
#endif

#ifndef CAST_DEBUG_DROP_EXCEPTIONS
  try
  {
#endif

    //////////////////////////
    //                      //
    //      Preparation     //
    //                      //
    //////////////////////////

    std::cout << scon::c3_delimeter(',');

    // start execution and initialization timer
    scon::chrono::high_resolution_timer exec_timer, init_timer;

    // initialize (old) Random Number Generator
    srand((unsigned int)time(NULL) + pid_func());

    // Parse config file and command line 
    auto config_filename = config::config_file_from_commandline(argc, argv);
    Config main_configuration(config_filename);
    config::parse_command_switches(argc, argv);

    // Print configuration
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << "\n";
      std::cout << "  |-----------------------------------------------------|\n";
      std::cout << "  |                                                     |\n";
      std::cout << "  |  //////////   //////////   //////////   //////////  |\n";
      std::cout << "  |  //           //      //   //               //      |\n";
      std::cout << "  |  //           //      //   //               //      |\n";
      std::cout << "  |  //           //      //   //////////       //      |\n";
      std::cout << "  |  //           //////////           //       //      |\n";
      std::cout << "  |  //           //      //           //       //      |\n";
      std::cout << "  |  //           //      //           //       //      |\n";
      std::cout << "  |  //////////   //      //   //////////       //      |\n";
      std::cout << "  |                                                     |\n";
      std::cout << "  |       conformational analysis and search tool       |\n";
      std::cout << "  |                                                     |\n";
      std::cout << "  |-----------------------------------------------------|\n\n\n";

      std::cout << "-------------------------------------------------------\n";
      std::cout << "Configuration ('" << config_filename << "')\n";
      std::cout << "-------------------------------------------------------\n";
      std::cout << Config::get().general;
      std::cout << Config::get().coords;
      std::cout << Config::get().energy;
      std::cout << Config::get().periodics;
    }

    //////////////////////////
    //                      //
    //    Initialization    //
    //                      //
    //////////////////////////

    if (Config::get().general.energy_interface == config::interface_types::T::DFTBABY)
    {  
#ifdef USE_PYTHON
#else
      printf("It is not possible to use DFTBaby without python!\n");
      std::exit(0);
#endif
      std::remove("output_dftb.txt"); // delete dftbaby output files from former run
      std::remove("tmp_struc_trace.xyz");
    }

    // read coordinate input file
    // "ci" contains all the input structures
    std::unique_ptr<coords::input::format> ci(coords::input::new_format());
    coords::Coordinates coords(ci->read(Config::get().general.inputFilename));

    // Print "Header"
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << "-------------------------------------------------\n";
      std::cout << "Initialization\n";
      std::cout << "-------------------------------------------------\n";
      std::cout << "Loaded " << ci->size() << " structure" << (ci->size() == 1 ? "" : "s");
      std::cout << ". (" << ci->atoms() << " atom" << (ci->atoms() == 1 ? "" : "s");
      std::cout << " and " << sys_mass(coords) << " g/mol";
      std::cout << (ci->size() > 1 ? " each" : "") << ")\n";
      std::size_t const susysize(coords.subsystems().size());
      if (susysize > 1U)
      {
        std::cout << susysize << " subsystems: ";
        for (std::size_t i(0U); i < susysize; ++i)
        {
          std::size_t const atms(coords.subsystems(i).size());
          std::cout << "[#" << i + 1 << " with " << atms << " atom" << (atms == 1 ? ".]" : "s.]");
        }
        std::cout << '\n';
      }
    }

    // If Periodic Boundry Conditions are used, translate all structures
    // so that their center of mass is on the origin of the coordinate system
    if (Config::get().periodics.periodic)
    {
      for (auto & pes : *ci)
      {
        coords.set_xyz(pes.structure.cartesian);
        coords.move_all_by(-coords.center_of_mass());
        pes = coords.pes();
      }
      // If Cutout option is on, cut off all atoms outside of box + radius
      if (Config::get().periodics.periodicCutout)
      {
        coords::Coordinates newCoords(coords);
        for (auto & pes : *ci)
        {
          newCoords.set_xyz(pes.structure.cartesian);
          newCoords = periodicsHelperfunctions::periodicCutout(coords);
          pes = newCoords.pes();
        }
        newCoords.set_xyz(ci->structure(0u).structure.cartesian);
        coords = newCoords;
      }
    }



    // stop and print initialization time
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << "-------------------------------------------------\n";
      std::cout << "Initialization done after " << init_timer << '\n';
    }

    //////////////////////////
    //                      //
    //    Preoptimization   //
    //                      //
    //////////////////////////

    if (coords.preoptimize())
    {
      if (Config::get().general.verbosity > 1U)
      {
        std::cout << "-------------------------------------------------\n";
        std::cout << "Preoptimization:\n";
        std::cout << "-------------------------------------------------\n";
      }
      std::size_t i(0);
      for (auto & pes : *ci)
      {
        coords.set_xyz(pes.structure.cartesian);
        coords.pe();
        if (Config::get().general.verbosity > 1U)
        {
          std::cout << "Preoptimization initial: " << ++i << '\n';
          coords.e_head_tostream_short(std::cout, coords.preinterface());
          coords.e_tostream_short(std::cout, coords.preinterface());
        }
        coords.po();
        pes.structure.cartesian = coords.xyz();
        if (Config::get().general.verbosity > 1U)
        {
          std::cout << "Preoptimization post-opt: " << i << '\n';
          coords.e_tostream_short(std::cout, coords.preinterface());
        }
      }
    }


    //////////////////////////
    //                      //
    //       T A S K S      //
    //                      //
    //////////////////////////

    // print which task
    std::cout << "-------------------------------------------------\n";
    std::cout << "Task '" << config::task_strings[Config::get().general.task];
    std::cout << "' (" << Config::get().general.task << ") computation:\n";
    std::cout << "-------------------------------------------------\n";

    // start task timer
    scon::chrono::high_resolution_timer task_timer;

    // select task
    switch (Config::get().general.task)
    {
    case config::tasks::DEVTEST:
    {
      // DEVTEST: Room for Development testing
      break;
    }
    case config::tasks::SP:
      { // singlepoint
        std::size_t i(0u);
        auto sp_energies_fn = coords::output::filename("_SP", ".txt");
        std::ofstream sp_estr(sp_energies_fn, std::ios_base::out);
        if (!sp_estr) throw std::runtime_error("Cannot open '" + 
          sp_energies_fn + "' to write SP energies.");
        sp_estr << std::setw(16) << "#";
        short_ene_stream_h(coords, sp_estr, 16);
        sp_estr << std::setw(16) << 't';
        sp_estr << '\n';
        i = 0;
        for (auto const& pes : *ci)
        {
          using namespace std::chrono;
          coords.set_xyz(pes.structure.cartesian, true);
          auto start = high_resolution_clock::now();
          coords.e();
          auto tim = duration_cast<duration<double>>
            (high_resolution_clock::now() - start);
          std::cout << "Structure " << ++i << " (" << tim.count() << " s)" << '\n';
          short_ene_stream(coords, sp_estr, 16);
          sp_estr << std::setw(16) << tim.count() << '\n';
          coords.e_head_tostream_short(std::cout);
          coords.e_tostream_short(std::cout);
        }
        break;
      }
    case config::tasks::GRAD:
    {
      // calculate gradient
      std::size_t i(0u);
      std::ofstream gstream(coords::output::filename("_GRAD", ".txt").c_str());
      for (auto const & pes : *ci)
      {
        coords.set_xyz(pes.structure.cartesian);
        coords.g();
        std::cout << "Structure " << ++i << '\n';
        coords.e_head_tostream_short(std::cout);
        coords.e_tostream_short(std::cout);
        coords.energyinterface()->print_G_tinkerlike(gstream);
      }
      break;
    }
    case config::tasks::HESS:
    {
      // calculate hessian matrix
      coords.e_head_tostream_short(std::cout);
      std::size_t i(0u);
      std::ofstream gstream(coords::output::filename("_HESS", ".txt").c_str());
      for (auto const & pes : *ci)
      {
        coords.set_xyz(pes.structure.cartesian);
        coords.h();
        std::cout << "Structure " << ++i << '\n';
        coords.e_tostream_short(std::cout);
        coords.h_tostream(gstream);
      }
      break;
    }
    case config::tasks::LOCOPT:
    {
      // local optimization
      std::remove("trace.arc"); // delete trace.arc file from former run
      coords.e_head_tostream_short(std::cout);
      auto lo_structure_fn = coords::output::filename("_LOCOPT");
      std::ofstream locoptstream(lo_structure_fn, std::ios_base::out);
      if (!locoptstream) throw std::runtime_error("Cannot open '" + lo_structure_fn + "' for LOCOPT structures.");
      auto lo_energies_fn = coords::output::filename("_LOCOPT", ".txt");
      std::ofstream loclogstream(lo_energies_fn, std::ios_base::out);
      if (!loclogstream) throw std::runtime_error("Cannot open '" + lo_structure_fn + "' for LOCOPT energies.");
      loclogstream << std::setw(16) << "#";
      short_ene_stream_h(coords, loclogstream, 16);
      short_ene_stream_h(coords, loclogstream, 16);
      loclogstream << std::setw(16) << "t";
      loclogstream << '\n';
      std::size_t i(0U);
      for (auto const & pes : *ci)
      {
        using namespace std::chrono;
        auto start = high_resolution_clock::now();
        coords.set_xyz(pes.structure.cartesian);
        coords.e();
        std::cout << "Initial: " << ++i << '\n';
        coords.e_tostream_short(std::cout);
        loclogstream << std::setw(16) << i;
        short_ene_stream(coords, loclogstream, 16);
        coords::Representation_3D oldC = coords.xyz();
        coords.o();
        coords::Representation_3D newC = coords.xyz();
        auto tim = duration_cast<duration<double>>
          (high_resolution_clock::now() - start);
        short_ene_stream(coords, loclogstream, 16);
        loclogstream << std::setw(16) << tim.count() << '\n';
        std::cout << "Post-Opt: " << i << "(" << tim.count() << " s)\n";
        coords.e_tostream_short(std::cout);
        locoptstream << coords;

        // calculate RMSD
        double sum_d_square=0, sum_d_square_not_fixed = 0;
        for (auto i = 0u; i < coords.size(); i++)
        {
          sum_d_square += dist(oldC[i], newC[i]) * dist(oldC[i], newC[i]);
          if (is_in(i, Config::get().coords.fixed) == false)
          {
            sum_d_square_not_fixed += dist(oldC[i], newC[i]) * dist(oldC[i], newC[i]);
          }
        }
        double rmsd = std::sqrt(sum_d_square / coords.size());
        double rmsd_not_fixed = std::sqrt(sum_d_square_not_fixed / (coords.size()- Config::get().coords.fixed.size()));
        std::cout << "RMSD between starting and optimized structure is " << rmsd << " angstrom.\n";
        if (Config::get().coords.fixed.size() != 0) std::cout << "If taking into account only non-fixed atoms it is " << rmsd_not_fixed << " angstrom.\n";
        loclogstream << "\nRMSD: " << rmsd << "\nRMSD(only_non_fixed): " << rmsd_not_fixed << "\n";
      }
      break;
    }
    case config::tasks::TS:
    {
      // Gradient only tabu search
      std::cout << Config::get().coords.equals;
      std::cout << "-------------------------------------------------\n";
      std::cout << Config::get().optimization.global;
      std::cout << "-------------------------------------------------\n";
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
    {
      // MonteCarlo Simulation
      std::cout << Config::get().coords.equals;
      std::cout << "-------------------------------------------------\n";
      std::cout << Config::get().optimization.global;
      std::cout << "-------------------------------------------------\n";
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
    {
      // Grid Search
      std::cout << Config::get().coords.equals;
      std::cout << "-------------------------------------------------\n";
      std::cout << Config::get().optimization.global;
      std::cout << "-------------------------------------------------\n";
      optimization::global::optimizers::main_grid mc(coords, ci->PES(),
        Config::get().optimization.global.grid.main_delta);
      mc.run(Config::get().optimization.global.iterations, true);
      mc.write_range("_GRID");
      break;
    }
    case config::tasks::INTERNAL:
    {
      // Explicitly shows CAST conversion to internal coordiantes
      // Beware when chaning this, PCA-task depend on this output and need to be adjusted accordingly.
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
          std::cout << " : " << coords.main(i) << '\n';
        }

        for (auto const & e : coords.main()) std::cout << e << '\n';

      }
      break;
    }
    case config::tasks::DIMER:
    {
      // Dimer method
      coords.e_head_tostream_short(std::cout);
      std::size_t i(0U);
      std::ofstream dimerstream(coords::output::filename("_DIMERTRANS").c_str(), std::ios_base::out);
      for (auto const & pes : *ci)
      {
        coords.set_xyz(pes.structure.cartesian);
        coords.o();
        dimerstream << coords;
        std::cout << "Pre-Dimer Minimum" << ++i << '\n';
        coords.e_tostream_short(std::cout);
        coords.dimermethod_dihedral();
        dimerstream << coords;
        std::cout << "Dimer Transition" << i << '\n';
        coords.e_tostream_short(std::cout);
        coords.o();
        dimerstream << coords;
        std::cout << "Post-Dimer Minimum" << i << '\n';
        coords.e_tostream_short(std::cout);
      }
      break;
    }
    case config::tasks::MD:
    {
      // Molecular Dynamics Simulation
      if (Config::get().md.pre_optimize) coords.o();
      md::simulation mdObject(coords);
      mdObject.run();
      break;
    }
    case config::tasks::FEP:
    {
      // Free energy perturbation
      md::simulation mdObject(coords);
      mdObject.fepinit();
      mdObject.run();
      mdObject.feprun();
      break;
    }
    case config::tasks::UMBRELLA:
    {
      // Umbrella Sampling
      Config::set().md.umbrella = true;
      if (Config::get().md.pre_optimize) coords.o();
      md::simulation mdObject(coords);
      mdObject.umbrella_run();
      break;
    }
    case config::tasks::STARTOPT:
    {
      // Preoptimization
      //std::cout << "PreApply.\n";
      startopt::apply(coords, ci->PES());
      //std::cout << "PostApply.\n";
      std::ofstream gstream(coords::output::filename("_SO").c_str());
      for (auto const & pes : ci->PES())
      {
        //std::cout << "PreSet.\n";
        coords.set_pes(pes, true);
        //std::cout << "PostSet.\n";
        gstream << coords;
      }
      break;
    }
    case config::tasks::GOSOL:
    { // Combined Solvation + Global Optimization
      std::cout << Config::get().startopt.solvadd;
      std::cout << "-------------------------------------------------\n";
      std::cout << Config::get().coords.equals;
      std::cout << "-------------------------------------------------\n";
      std::cout << Config::get().optimization.global;
      std::cout << "-------------------------------------------------\n";
      startopt::preoptimizers::GOSol sopt(coords, ci->PES());
      sopt.run(Config::get().startopt.solvadd.maxNumWater);
      break;
    }
    case config::tasks::NEB:
    {
      std::ptrdiff_t counter = 0;
      std::vector<coords::Representation_3D> input_pathway;
      coords::Representation_3D start_struc, final_struc;
      ptrdiff_t image_connect = ptrdiff_t(Config::get().neb.CONNECT_NEB_NUMBER);

      for (auto const & pes : *ci)
      {
        coords.set_xyz(pes.structure.cartesian);
        coords.mult_struc_counter++;
        if (Config::get().neb.COMPLETE_PATH)
        {
          input_pathway.push_back(pes.structure.cartesian);
        }
        else if (!Config::get().neb.MULTIPLE_POINTS)
        {

          neb nobj(&coords);
          nobj.preprocess(counter);
        }
      }
      if (Config::get().neb.COMPLETE_PATH && !(Config::get().neb.MULTIPLE_POINTS))
      {
        neb nobj(&coords);
        nobj.preprocess(input_pathway, counter);
      }
      else if ((Config::get().neb.MULTIPLE_POINTS))
      {
        for (size_t i = 0; i < (input_pathway.size() - 1); ++i)
        {
          start_struc = input_pathway[i];
          final_struc = input_pathway[i + 1];
          neb nobj(&coords);
          nobj.preprocess(counter, image_connect, counter, start_struc, final_struc, true);
        }
      }
      break;
    }
    case config::tasks::PATHOPT:
    {
      std::ptrdiff_t counter = 0;
      for (auto const & pes : *ci)
      {
        coords.set_xyz(pes.structure.cartesian);
        coords.mult_struc_counter++;
        neb nobj(&coords);
        nobj.preprocess(counter);
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
    case config::tasks::WRITE_TINKER:
    {
      std::ofstream gstream(coords::output::filename("", ".arc").c_str());
      gstream << coords::output::formats::tinker(coords);
      break;
    }
		case config::tasks::WRITE_GAUSSVIEW:
		{
			std::ofstream gstream(coords::output::filename("", ".gjf").c_str());
			if (Config::get().general.energy_interface == config::interface_types::ONIOM || Config::get().general.energy_interface == config::interface_types::QMMM)
			{
				gstream << "# ONIOM(HF/6-31G:UFF)\n\n";                      // method
				gstream << "some stupid title\n\n";
				gstream << "0 1 0 1 0 1\n";                                  // charge and multiplicity
        gstream << std::fixed << std::setprecision(3);
				for (auto i = 0u; i < coords.size(); ++i)                    // atom information
				{
					auto symbol = coords.atoms(i).symbol();
					auto x = coords.xyz(i).x();
					auto y = coords.xyz(i).y();
					auto z = coords.xyz(i).z();
					auto system = "L";
					if (is_in(i, Config::get().energy.qmmm.qmatoms)) system = "H";
					gstream << symbol << "\t" << x << "\t" << y << "\t" << z << "\t" << system << "\n";
				}
				gstream << "\n";
			}
			else if (Config::get().general.energy_interface == config::interface_types::THREE_LAYER)
			{
				gstream << "# ONIOM(B3LYP/6-31G:HF/STO-3G:UFF)\n\n";                      // method
				gstream << "some stupid title\n\n";
				gstream << "0 1 0 1 0 1 0 1 0 1\n";                                       // charge and multiplicity
        gstream << std::fixed << std::setprecision(3);
				for (auto i = 0u; i < coords.size(); ++i)                    // atom information
				{
					auto symbol = coords.atoms(i).symbol();
					auto x = coords.xyz(i).x();
					auto y = coords.xyz(i).y();
					auto z = coords.xyz(i).z();
					auto system = "L";
					if (is_in(i, Config::get().energy.qmmm.qmatoms)) system = "H";
					else if (is_in(i, Config::get().energy.qmmm.seatoms)) system = "M";
					gstream << symbol << "\t" << x << "\t" << y << "\t" << z << "\t" << system << "\n";
				}
				gstream << "\n";
			}
			else  // input for "normal molecule"
			{
				gstream << "# HF/6-31G\n\n";                      // method
				gstream << "some stupid title\n\n";
				gstream << "0 1\n";                               // charge and multiplicity
				gstream << coords::output::formats::xyz(coords);  // coordinates
				gstream << "\n";
			}
			break;
		}
    case config::tasks::MODIFY_SK_FILES:
    {
      std::vector<std::vector<std::string>> pairs = find_pairs(coords);
      for (auto p : pairs)
      {
        modify_file(p);
      }
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

       // Create empty pointer since we do not know yet if PCA eigenvectors etc.
       // will be generated from coordinates or read from file
      pca::PrincipalComponentRepresentation* pcaptr = nullptr;

      // Create new PCA eigenvectors and modes
      if (!Config::get().PCA.pca_read_modes && !Config::get().PCA.pca_read_vectors)
      {
        pcaptr = new pca::PrincipalComponentRepresentation(ci, coords);
        pcaptr->writePCAModesFile("pca_modes.dat");
      }
      // Read modes and eigenvectors from (properly formated) file "pca_modes.dat"
      else if (Config::get().PCA.pca_read_modes && Config::get().PCA.pca_read_vectors) pcaptr = new pca::PrincipalComponentRepresentation("pca_modes.dat");
      else
      {
        pcaptr = new pca::PrincipalComponentRepresentation(ci, coords);
        // Read PCA-Modes from file but generate new eigenvectors from input coordinates
        if (Config::get().PCA.pca_read_modes) pcaptr->readModes("pca_modes.dat");
        // Read PCA-Eigenvectors from file but generate new modes using the eigenvectors
        // and the input coordinates
        else if (Config::get().PCA.pca_read_vectors)
        {
          pcaptr->readEigenvectors("pca_modes.dat");
          pcaptr->generatePCAModesFromPCAEigenvectorsAndCoordinates();
        }
      }

      // If modes or vectors have changed, write them to new file
      if (Config::get().PCA.pca_read_modes != Config::get().PCA.pca_read_vectors) pcaptr->writePCAModesFile("pca_modes_new.dat");

      // Create Histograms
      // ATTENTION: This function read from Config::PCA
      pcaptr->writeHistogrammedProbabilityDensity("pca_histogrammed.dat");

      // Write Stock's Delta, see DOI 10.1063/1.2746330
      // ATTENTION: This function read from Config::PCA
      pcaptr->writeStocksDelta("pca_stocksdelta.dat");

      // Cleanup
      delete pcaptr;
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
      pca::ProcessedPrincipalComponentRepresentation pcaproc("pca_modes.dat");
      pcaproc.determineStructures(ci, coords);
      pcaproc.writeDeterminedStructures(coords);
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
    case config::tasks::REMOVE_EXPLICIT_WATER:
    {
      /**
      * THIS TASK REMOVES EXPLICIT WATER FROM STRUCTURES AND WRITES THE TRUNCATED STRUCTURES TO FILE
      *
      */

      std::ofstream out(coords::output::filename("_noexplwater").c_str(), std::ios::app);
      std::string* hold_str = new std::string[ci->size()];
#ifdef _OPENMP
      auto const n_omp = static_cast<std::ptrdiff_t>(ci->size());
#pragma omp parallel for firstprivate(coords) shared(hold_str)
      for (std::ptrdiff_t iter = 0; iter < n_omp; ++iter)
#else
      for (std::size_t iter = 0; iter < ci->size(); ++iter)
#endif
      {
        auto holder = ci->PES()[iter].structure.cartesian;
        coords.set_xyz(holder);

        std::vector<size_t> atomsToBePurged;
        coords::Atoms truncAtoms;
        coords::Representation_3D positions;
        for (size_t i = 0u; i < coords.atoms().size(); i++)
        {
          coords::Atom atom(coords.atoms().atom(i));
          if (atom.number() != 8u && atom.number() != 1u)
          {
            truncAtoms.add(atom);
            positions.push_back(coords.xyz(i));
          }
          else if (atom.number() == 1u)
          {
            // Check if hydrogen is bound to something else than Oxygen
            bool checker = true;
            for (size_t j = 0u; j < atom.bonds().size(); j++)
            {
              if (coords.atoms().atom(atom.bonds()[j]).number() == 8u) checker = false;
            }
            if (checker)
            {
              truncAtoms.add(atom);
              positions.push_back(coords.xyz(i));
            }
          }
          else if (atom.number() == 8u)
          {
            //checker checks if only hydrogens are bound to this current oxygen
            bool checker = true;
            for (size_t j = 0u; j < atom.bonds().size(); j++)
            {
              if (coords.atoms().atom(atom.bonds()[j]).number() != 1u) checker = false;
            }
            if (!checker)
            {
              truncAtoms.add(atom);
              positions.push_back(coords.xyz(i));
              for (auto const& bond : atom.bonds())
              {
                if (coords.atoms().atom(bond).number() == 1u)
                {
                  truncAtoms.add(coords.atoms().atom(bond));
                  positions.push_back(coords.xyz(bond));
                }
              }
            }
          }
        }
        coords::Coordinates newCoords;
        coords::PES_Point x(positions);
        newCoords.init_in(truncAtoms, x);
        std::stringstream temporaryStringstream;
        temporaryStringstream << newCoords;
        hold_str[iter] = temporaryStringstream.str();
      }
      for (size_t i = 0; i < ci->size(); i++)

      {
        out << hold_str[i];
      }
		    break;
      }
	  case config::tasks::SCAN2D:
	  {
		  auto scan = std::make_shared<Scan2D>(coords);
		scan->execute_scan();
		  break;
	  }
	    case config::tasks::XB_EXCITON_BREAKUP:
	    {
		  /**
		  * THIS TASK SIMULATES THE EXCITON_BREAKUP ON AN 
		  * INTERFACE OF TWO ORGANIC SEMICONDUCTORS: 
		  * (AT THE MOMENT ONLY ORGANIC SEMICONDUCTOR/FULLERENE INTERFACE)
		  * NEEDS SPECIALLY PREPEARED INPUT
		  */  
		  exciton_breakup(Config::get().exbreak.pscnumber, Config::get().exbreak.nscnumber, Config::get().exbreak.interfaceorientation, Config::get().exbreak.masscenters, 
						 Config::get().exbreak.nscpairrates, Config::get().exbreak.pscpairexrates, Config::get().exbreak.pscpairchrates, Config::get().exbreak.pnscpairrates);
      break;
	  }
      case config::tasks::XB_INTERFACE_CREATION:
      {
      /**
      * THIS TASK CREATES A NEW COORDINATE SET FROM TWO PRECURSORS
      */
        //creating second coords object
        std::unique_ptr<coords::input::format> add_strukt_uptr(coords::input::additional_format());
        coords::Coordinates add_coords(add_strukt_uptr->read(Config::get().interfcrea.icfilename));
        coords::Coordinates newCoords(coords);

 
        newCoords = periodicsHelperfunctions::interface_creation(Config::get().interfcrea.icaxis, Config::get().interfcrea.icdist, coords, add_coords);

        coords = newCoords;

        std::ofstream new_structure(Config::get().general.outputFilename, std::ios_base::out);
        new_structure << coords;

        break;
      }
      case config::tasks::XB_CENTER:
      {
        /**
        * THIS  TASK CALCULATES THE CENTERS OF MASSES FOR ALL MONOMERS IN THE STRUCTURE AND IF WANTED GIVES STRUCTURE FILES FOR DIMERS
        * WITHIN A DEFINED DISTANCE BETWEEN THE MONOMERS
        */

        center(coords);
        break;
      }
      case config::tasks::XB_COUPLINGS:
      {
        couplings::coupling coup;

        coup.kopplung();

        break;
      }
      case config::tasks::LAYER_DEPOSITION:
      {
        //Generating layer with random defects
        coords::Coordinates newCoords(coords);
        coords::Coordinates inp_add_coords(coords);
        coords::Coordinates add_coords;

        for (auto & pes : *ci)
        {
          newCoords.set_xyz(pes.structure.cartesian);
          newCoords = periodicsHelperfunctions::delete_random_molecules(coords, Config::get().layd.del_amount);
          pes = newCoords.pes();
        }
        newCoords.set_xyz(ci->structure(0u).structure.cartesian);
        coords = newCoords;

        

        for (std::size_t i = 0u; i < 1; i++) //this loops purpose is to ensure mdObject1 is destroyed before further changes to coords happen and the destructor goes bonkers. Not elegant but does the job.
        {
          //Molecular Dynamics Simulation
          if (Config::get().md.pre_optimize) coords.o();
          md::simulation mdObject1(coords);
          mdObject1.run();
        }

        for (std::size_t i = 1; i < Config::get().layd.amount; i++)
        {
          add_coords = inp_add_coords;
          add_coords = periodicsHelperfunctions::delete_random_molecules(add_coords, Config::get().layd.del_amount);
  
          for (auto & pes : *ci)
          {
            newCoords.set_xyz(pes.structure.cartesian);
            newCoords = periodicsHelperfunctions::interface_creation(Config::get().layd.laydaxis, Config::get().layd.layddist, coords, add_coords);
            pes = newCoords.pes();
          }
          newCoords.set_xyz(ci->structure(0u).structure.cartesian);
          coords = newCoords;
                                                                                                                            

          for (std::size_t j = 0; j < (coords.size() - add_coords.size()); j++)//fix all atoms already moved by md
          {
            coords.set_fix(j, true);
          }

          // Molecular Dynamics Simulation
          if (Config::get().md.pre_optimize) coords.o();
          md::simulation mdObject2(coords);
          mdObject2.run();
        }

        std::size_t mon_amount_type1 = coords.molecules().size();//save number of molecules of first kind for later use in replacement

        //option if a heterogenous structure shall be created
        if (Config::get().layd.hetero_option == true)
        {
          std::unique_ptr<coords::input::format> sec_strukt_uptr(coords::input::additional_format());
          coords::Coordinates sec_coords(sec_strukt_uptr->read(Config::get().layd.layd_secname));
          coords::Coordinates add_sec_coords;

          for (std::size_t i = 0; i < Config::get().layd.sec_amount; i++)
          {
            add_sec_coords = sec_coords;
            add_sec_coords = periodicsHelperfunctions::delete_random_molecules(add_sec_coords, Config::get().layd.sec_del_amount);

            for (auto & pes : *ci)
            {
              newCoords.set_xyz(pes.structure.cartesian);
              newCoords = periodicsHelperfunctions::interface_creation(Config::get().layd.laydaxis, Config::get().layd.sec_layddist, coords, add_sec_coords);
              pes = newCoords.pes();
            }
            newCoords.set_xyz(ci->structure(0u).structure.cartesian);
            coords = newCoords;

            for (std::size_t j = 0; j < (coords.size() - add_sec_coords.size()); j++)//fix all atoms already moved by md
            {
              coords.set_fix(j, true);
            }

            // Molecular Dynamics Simulation
            if (Config::get().md.pre_optimize) coords.o();
            md::simulation mdObject3(coords);
            mdObject3.run();
          }
        }

        //option if monomers in structure shall be replaced
        if (Config::get().layd.replace == true)
        {
          std::unique_ptr<coords::input::format> add_strukt_uptr(coords::input::additional_format());
          coords::Coordinates add_coords1(add_strukt_uptr->read(Config::get().layd.reference1));
          std::unique_ptr<coords::input::format> add_strukt_uptr2(coords::input::additional_format());
          coords::Coordinates add_coords2(add_strukt_uptr2->read(Config::get().layd.reference2));

          coords = monomerManipulation::replaceMonomers(coords, add_coords1, add_coords2, mon_amount_type1);
        }

        std::ofstream output(Config::get().general.outputFilename, std::ios_base::out);       
        output << coords;
        break;
      }

      default:
      {
      
      }
    }
#ifdef USE_PYTHON
      Py_Finalize(); //  close python
#endif 

    // stop and print task and execution time
    std::cout << '\n' << "Task " << config::task_strings[Config::get().general.task];
    std::cout << " took " << task_timer << " to complete.\n";
    std::cout << "Execution of " << config::Programname << " (" << config::Version << ")";
    std::cout << " ended after " << exec_timer << '\n';

    //////////////////////////
    //                      //
    //       EXCEPTION      //
    //       HANDLING       //
    //                      //
    //////////////////////////
#ifndef CAST_DEBUG_DROP_EXCEPTIONS
  }
#if defined COMPILEX64 || defined __LP64__ || defined _WIN64 
  catch (std::bad_alloc &)
  {
    std::cout << "Memory allocation failure. Input structure probably too large.\n";
  }
#else
  catch (std::bad_alloc &)
  {
    std::cout << "Memory allocation failure. CAST probably ran out of memory. Try using 64bit compiled " << config::Programname << " instead.\n";
  }
#endif
  catch (std::exception & e)
  {
    std::cout << "An exception occured. The execution of " << config::Programname << " failed. \n";
    std::cout << "Error: " << e.what() << '\n';
  }
#endif
#ifdef _MSC_VER 
  // make window stay open in debug session on windows
  if (IsDebuggerPresent()) std::system("pause");
#endif
  return 0;
  }
#endif
