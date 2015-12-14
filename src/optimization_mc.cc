#include <limits>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include "coords.h"
#include "optimization_global.h"
#include "scon_vect.h"
#include "scon_utility.h"
#include "configuration.h"
#include "scon_chrono.h"
#include "scon.h"

optimization::global::optimizers::monteCarlo::monteCarlo(
  coords::Coordinates &c, std::string const & output_name, bool mc_in_neb)
  : optimizer(c, output_name), neb_true(mc_in_neb)
{
  neb_true = mc_in_neb;
}

optimization::global::optimizers::monteCarlo::monteCarlo(
  coords::Coordinates &c, coords::Ensemble_PES const &p, 
  std::string const & output_name, bool mc_in_neb)
: optimizer(c, p, false, output_name), neb_true(false)
{
	neb_true = false;
}




bool optimization::global::optimizers::monteCarlo::run(std::size_t const iterations, bool const reset)
{
  coords::float_type ene(0.0);
  if (reset)
  {
    i = 0;
    opt_clock.restart();
    header_to_cout();
  }
  coordobj.set_pes(accepted_minima[min_index].pes);
  init_stereo = coordobj.stereos();
  // Print header
  std::size_t const iter_size(num_digits(Config::get().optimization.global.iterations) + 1);
  
  std::string const method(std::string("MC").append(Config::get().optimization.global.montecarlo.minimization ? "M " : "  "));
  for (; i<iterations; ++i)
  {
    scon::chrono::high_resolution_timer step_timer;
    // move 
    if (Config::get().optimization.global.move_dehydrated)
    {
      std::size_t const current_size(coordobj.size());
      // create temporary dehydrated coords 
      // if dehydrated movement is required
      coords::Coordinates tmp(dehydrate(coordobj));
      // move in dehydrated coord space
      move(tmp);
      // set coords to rehydrated version of dehydrated coords
      coordobj = rehydrate(tmp);
      if (coordobj.size() != current_size) throw std::logic_error("Rehydrated coords after dehydrated movement exhibit wrong number of atoms.");
    } else move(coordobj);
    coordobj.to_internal();
    coordobj.to_xyz();
   // if (coordobj.preoptimize()) coordobj.po();
    // Print Method, Iteration, Start Minimum and Transition energy

    if (Config::get().general.verbosity > 1U)
    {
      ene = coordobj.e();
      std::cout << method;
      // iteration
      std::cout << std::right << std::setw(iter_size) << i + 1 << "/";
      std::cout << std::left << std::setw(iter_size) << Config::get().optimization.global.iterations << ' ';
      // Curennt minimum index
      std::cout << std::setw(10) << min_index << ' ';
      // Current minimum
      std::cout << std::setw(Config::get().optimization.global.precision + 10) << std::left;
      std::cout << std::setprecision(Config::get().optimization.global.precision);
      std::cout << std::showpoint << std::scientific << accepted_minima[min_index].pes.energy << ' ';
      // Transition
      std::cout << std::setw(Config::get().optimization.global.precision + 10) << std::left;
      std::cout << std::setprecision(Config::get().optimization.global.precision);
      std::cout << std::showpoint << std::scientific << ene << ' ';
    }
    // Optimization if MC+Minimization
    if (Config::get().optimization.global.montecarlo.minimization)
    {
      ene = coordobj.o();
    }
    // if verbosity < 2 we did not calc the energy yet
    else if (Config::get().general.verbosity < 2U)
    {
      ene = coordobj.e();
    }
    coordobj.to_internal();
    coordobj.to_xyz();
    // Check whether pes point is to be accepted
    min_status::T status(check_pes_of_coords());
    // Print final energy
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << std::setw(Config::get().optimization.global.precision + 10) << std::left;
      std::cout << std::setprecision(Config::get().optimization.global.precision);
      std::cout << std::showpoint << std::scientific << ene << ' ';
      std::cout << status;
      std::cout << std::setw(10) << std::left << accepted_minima.size();
      std::cout << std::setw(10) << std::left << range_minima.size();
    }
    if (Config::get().general.verbosity > 1U)
    {
      std::cout << "(" << std::setprecision(2) << std::showpoint;
      std::cout << std::fixed << T << " K, " << step_timer << ")" << lineend;
    }
    if (!restore(status))
    {
      if (Config::get().general.verbosity > 1U)
      {
          std::cout << "Starting point selection limit reached (no non-banned minimum accessible). Stop." << lineend;
      }
      break;
    }
  }
  --i;
  return found_new_minimum;
}


void optimization::global::optimizers::monteCarlo::move(coords::Coordinates &c)
{
  switch (Config::get().optimization.global.montecarlo.move)
  {
    case config::optimization_conf::mc::move_types::XYZ:
    {
      move_xyz(c);
      break;
    }
    case config::optimization_conf::mc::move_types::DIHEDRAL:
    {
      move_main(c);
      break;
    }
    default:
    {
      move_main_strain(c);
      break;
    }
  }
}


void optimization::global::optimizers::monteCarlo::move_xyz(coords::Coordinates &movecoords)
{
  const std::size_t N = movecoords.size();
  for (std::size_t ak = 0; ak < N; ++ak)
  {
    if (!movecoords.atoms(ak).fixed())
    {
      coords::Cartesian_Point move(normalized(scon::randomized<coords::Cartesian_Point>()));
      move *= Config::get().optimization.global.montecarlo.cartesian_stepsize * scon::rand<coords::float_type>(0, 1);
      movecoords.move_atom_by(ak, move);

    }
  }
}


void optimization::global::optimizers::monteCarlo::move_main(coords::Coordinates &movecoords)
{
  std::size_t const N(movecoords.main().size());
  coords::Representation_Main main_tors(N);
  // choose number of torsions to modify
  //coords::float_type const NUM_MAINS(static_cast<coords::float_type>(N));
  std::size_t const NUM_MOD(std::min((static_cast<std::size_t>(-std::log(scon::rand<coords::float_type>(0.0,1.0))) + 1U), N));
  // apply those torsions
  if (Config::get().general.verbosity > 14U) std::cout << "Changing " << NUM_MOD << " of " << N << " mains." << lineend;
  for (std::size_t ii(0U); ii<NUM_MOD; ++ii)
  {
    std::size_t K = scon::rand<std::size_t>(0u, N - 1u);
    coords::float_type const F = scon::rand<coords::float_type>(-Config::get().optimization.global.montecarlo.dihedral_max_rot, 
                                        Config::get().optimization.global.montecarlo.dihedral_max_rot);
    if (Config::get().general.verbosity > 14U) std::cout << "Changing main " << K << " by " << F << lineend;
    main_tors[K] = coords::main_type::from_deg(F);
  }
  // orthogonalize movement to main directions
  for (auto tabu_direction : accepted_minima[min_index].main_direction)
  {
    tabu_direction.resize(main_tors.size());
    scon::orthogonalize_to_normal(main_tors, tabu_direction);
    //main_tors.orthogonalize_toNormed(tabu_direction);
  }
  for (std::size_t ii = 0; ii < N; ++ii)
  {
    movecoords.rotate_main(ii, main_tors[ii]);
  }
  //main_tors.print(cout);
  movecoords.to_xyz();
}


void optimization::global::optimizers::monteCarlo::move_main_strain(coords::Coordinates &movecoords)
{
  std::size_t const N(movecoords.main().size());
  coords::Representation_Main main_tors(N);
  // choose number of torsions to modify
  coords::float_type const NUM_MAINS(static_cast<coords::float_type>(N));
  std::size_t const NUM_MOD(std::min((static_cast<std::size_t>(-std::log(scon::rand<coords::float_type>(0.0, 1.0))) + 1U), N));
  // apply those torsions
  if (Config::get().general.verbosity > 9U) std::cout << "Changing " << NUM_MOD << " of " << N << " mains." << lineend;
  for (std::size_t ii(0U); ii<NUM_MOD; ++ii)
  {
    std::size_t const K(static_cast<std::size_t>(NUM_MAINS*scon::rand<coords::float_type>(0.0, 1.0)));
    coords::float_type const F(Config::get().optimization.global.montecarlo.dihedral_max_rot*((d_rand()*2.0) - 1.0));
    if (Config::get().general.verbosity > 9U) std::cout << "Changing main " << K << " by " << F << lineend;
    main_tors[K] = coords::main_type::from_deg(F);
  }
  // orthogonalize movement to main directions
  for (auto tabu_direction : accepted_minima[min_index].main_direction)
  {
    tabu_direction.resize(main_tors.size());
    scon::orthogonalize_to_normal(main_tors, tabu_direction);
    //main_tors.orthogonalize_toNormed(tabu_direction);
  }
  // add bias potentials to main torsions and dependant torsions to be current value + main_tors[i]
  for (std::size_t ii(0U); ii<N; ++ii)
  {
    if (abs(main_tors[ii]) > coords::main_type::from_deg(0)) bias_main_rot(movecoords, ii, main_tors[ii]);
  }
  // optimize using biased potentials
  movecoords.o();
  movecoords.to_internal();
  // obtain direction of movement and make it tabu
  //accepted_minima[min_index].main_direction.push_back(scon::normalized(main_tors - movecoords.main()));
  movecoords.potentials().clear();
  movecoords.potentials().append_config();
}


void optimization::global::optimizers::monteCarlo::bias_main_rot(coords::Coordinates &movecoords, std::size_t index, coords::angle_type const rot)
{
  if (index >= movecoords.main().size()) throw std::out_of_range("index out of range for main");
  std::size_t const intern(movecoords.atoms().intern_of_main_idihedral(index));
  config::biases::dihedral bias;
  bias.force = 10.0;
  for (auto const idih_index:movecoords.atoms(intern).dependant_idihedrals())
  {
    bias.a = movecoords.atoms(idih_index).i_to_a();
    bias.b = movecoords.atoms(movecoords.atoms(idih_index).ibond()).i_to_a();
    bias.c = movecoords.atoms(movecoords.atoms(idih_index).iangle()).i_to_a();
    bias.d = movecoords.atoms(movecoords.atoms(idih_index).idihedral()).i_to_a();
    coords::angle_type phi = scon::dihedral(movecoords.xyz()[bias.c] - movecoords.xyz()[bias.b], 
                                            movecoords.xyz()[bias.a], movecoords.xyz()[bias.d]);
    bias.ideal = phi + rot;
    movecoords.potentials().add(bias);
  }
}