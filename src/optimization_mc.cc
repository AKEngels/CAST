#include <limits>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <algorithm>
#include "coords_io.h"
#include "coords.h"
#include "optimization_global.h"
#include "scon_vect.h"
#include "scon_utility.h"
#include "configuration.h"
#include "scon_chrono.h"
#include "scon.h"
#include "startopt_solvadd.h"

optimization::global::optimizers::monteCarlo::monteCarlo(
  coords::Coordinates &c, std::string const & output_name, bool mc_in_neb)
  : optimizer(c, output_name), neb_true(mc_in_neb)
{ }

optimization::global::optimizers::monteCarlo::monteCarlo(
  coords::Coordinates &c, coords::Ensemble_PES const &p, 
  std::string const & output_name, bool mc_in_neb)
: optimizer(c, p, false, output_name), neb_true(mc_in_neb)
{ }

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
  std::size_t const iter_size(scon::num_digits(Config::get().optimization.global.iterations) + 1);
  std::string const method(std::string("MC").append(Config::get().optimization.global.montecarlo.minimization ? "M " : "  "));
  std::unique_ptr<std::ofstream> pref, transf, postf;
  //if (Config::get().optimization.global.track_all)
  //{
  //  pref.reset(new std::ofstream{ coords::output::filename("MC_pre").c_str() });
  //  transf.reset(new std::ofstream{ coords::output::filename("MC_trans").c_str() });
  //  postf.reset(new std::ofstream{ coords::output::filename("MC_post").c_str() });
  //}
  for (; i<iterations; ++i)
  {
    scon::chrono::high_resolution_timer step_timer;
    if (pref && *pref) { *pref << coordobj; }
    // move 
    std::size_t movecount(0);
    if (Config::get().optimization.global.move_dehydrated)
    {
      std::size_t const current_size(coordobj.size());
      // create temporary dehydrated coords 
      // if dehydrated movement is required
      coords::Coordinates tmp(dehydrate(coordobj));
      // move in dehydrated coord space
      movecount = move(tmp);
      // set coords to rehydrated version of dehydrated coords
      coordobj = rehydrate(tmp);
      if (coordobj.size() != current_size) throw std::logic_error("Rehydrated coords after dehydrated movement exhibit wrong number of atoms.");
    } else movecount = move(coordobj);
    coordobj.to_internal();
    coordobj.to_xyz();
    if (transf && *transf) { *transf << coordobj; }
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
      coordobj.to_internal();
      coordobj.to_xyz();
    }
    // if verbosity < 2 we did not calc the energy yet
    else if (Config::get().general.verbosity < 2U)
    {
      ene = coordobj.e();
    }
    if (postf && *postf) { *postf << coordobj; }
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
      std::cout << std::fixed << T << " K, " << step_timer << ", Changed " << movecount << ")\n";
    }
    if (!restore(status))
    {
      if (Config::get().general.verbosity > 1U)
      {
          std::cout << "Starting point selection limit reached (no non-banned minimum accessible). Stop.\n";
      }
      break;
    }
  }
  --i;
  return found_new_minimum;
}


std::size_t optimization::global::optimizers::monteCarlo::move(coords::Coordinates &c)
{
  using config::optimization_conf::mc;
  switch (Config::get().optimization.global.montecarlo.move)
  {
    case mc::move_types::XYZ:
    {
      return move_xyz(c);
      break;
    }
    case mc::move_types::DIHEDRAL:
    {
      return move_main(c);
      break;
    }
    case mc::move_types::WATER:
    {
      return move_water(c);
      break;
    }
    default:
    {
      return move_main_strain(c);
      break;
    }
  }
}


std::size_t optimization::global::optimizers::monteCarlo::move_xyz(coords::Coordinates &movecoords)
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
  return N;
}


std::size_t optimization::global::optimizers::monteCarlo::move_main(coords::Coordinates &movecoords)
{
  auto const n = movecoords.main().size();
  coords::Representation_Main main_tors(n);
  // choose number of torsions to modify
  auto dist = std::geometric_distribution<std::size_t>{ 0.5 };
  auto m = scon::random::threaded_rand(dist);
  m = std::min<std::size_t>(m, n);
  // apply those torsions
  if (Config::get().general.verbosity > 4U)
  {
    std::cout << "Changing " << m << " of " << n << " mains.\n";
  }
  auto maxrot = Config::get().optimization.global.montecarlo.dihedral_max_rot;
  for (std::size_t ii(0U); ii<m; ++ii)
  {
    std::size_t kk = scon::rand<std::size_t>(0u, n - 1u);
    coords::float_type const F = 
      scon::rand<coords::float_type>(-maxrot, maxrot);
    if (Config::get().general.verbosity > 4U)
    {
      std::cout << "Changing main " << k << " by " << F << '\n';
    }
    main_tors[kk] = coords::main_type::from_deg(F);
  }
  // orthogonalize movement to main directions
  for (auto tabu_direction : accepted_minima[min_index].main_direction)
  {
    tabu_direction.resize(main_tors.size());
    scon::orthogonalize_to_normal(main_tors, tabu_direction);
    //main_tors.orthogonalize_toNormed(tabu_direction);
  }
  //std::cout << scon::vector_delimeter('\n');
  {
    std::ofstream pref("pre.arc", std::ios::app);
    pref << movecoords;
  }
  auto tmpi = movecoords.intern();
  for (std::size_t ii = 0; ii < n; ++ii)
  {
    movecoords.rotate_main(ii, main_tors[ii], false);
  }
  //std::cout << (tmpi - movecoords.intern()) << '\n';
  //main_tors.print(cout);
  //std::cout << movecoords.xyz() << "\n";
  movecoords.to_xyz();
  {
    std::ofstream postf("post.arc", std::ios::app);
    postf << movecoords;
  }
  return m;
}


std::size_t optimization::global::optimizers::monteCarlo::move_main_strain(coords::Coordinates &movecoords)
{
  auto max_rot = std::abs(
    Config::get().optimization.global.montecarlo.dihedral_max_rot);
  auto const n = movecoords.main().size();
  coords::Representation_Main main_tors(n);
  // distributions
  auto select_distribution = std::uniform_int_distribution<std::size_t>{ 0, n-1u };
  auto rot_distribution = std::uniform_real_distribution<coords::float_type>{ -max_rot, max_rot };
  // choose number of torsions to modify
  auto m = std::min<std::size_t>(
    scon::random::threaded_rand(std::geometric_distribution<std::size_t>{ 0.5 }), n);
  // apply those torsions
  if (Config::get().general.verbosity > 4U)
  {
    std::cout << "Changing " << m << " of " << n << " mains.\n";
  }
  for (std::size_t ii(0U); ii<m; ++ii)
  {
    // obtain random rotation
    auto rotation = scon::random::threaded_rand(rot_distribution);
    // get angle type from rotation
    auto rot_main = coords::main_type::from_deg(rotation);
    // select id of main to be modified by 'rotation' degrees
    auto o = scon::random::threaded_rand(select_distribution);
    if (Config::get().general.verbosity > 4U)
    {
      std::cout << "Changing main " << o << " by " << rot_main << '\n';
    }
    // set main k to be rotated by 'rot_main'
    main_tors[o] = rot_main;
  }
  //// orthogonalize movement to main directions
  //for (auto tabu_direction : accepted_minima[min_index].main_direction)
  //{
  //  tabu_direction.resize(main_tors.size());
  //  scon::orthogonalize_to_normal(main_tors, tabu_direction);
  //  //main_tors.orthogonalize_toNormed(tabu_direction);
  //}
  // add bias potentials to main torsions 
  // and dependant torsions to be current value + main_tors[i]
  for (std::size_t ii(0U); ii<n; ++ii)
  {
    if (abs(main_tors[ii]) > coords::main_type::from_deg(0))
    {
      bias_main_rot(movecoords, ii, main_tors[ii]);
    }
  }
  // optimize using biased potentials
  movecoords.o();
  movecoords.to_internal();
  // obtain direction of movement and make it tabu
  //accepted_minima[min_index].main_direction.push_back(scon::normalized(main_tors - movecoords.main()));
  movecoords.potentials().clear();
  movecoords.potentials().append_config();
  return m;
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

namespace
{
  template<class IV>
  bool molecule_is_moveable_water(IV const &indices, coords::Atoms const &atm)
  {
    // more than 3 atoms in molecule? -> not water
    if (indices.size() != 3u) return false;
    // any fixed atoms? -> not moveable
    if (atm.atom(indices[0]).fixed() || 
      atm.atom(indices[1]).fixed() || 
      atm.atom(indices[2]).fixed()) return false;
    std::size_t nh{ 0 }, no{ 0 };
    for (auto i : indices)
    {
      if (atm.atom(i).number() == 1 && nh < 2) ++nh;
      else if (atm.atom(i).number() == 8 && no < 1) ++no;
      else return false;
    }
    return true;
  }

  template<class IV, class MV, class D3V>
  IV remove_burried_waters(
    IV const &waters, MV const & molecules, D3V const &xyz)
  {
    scon::linked::Cells < coords::float_type, coords::Cartesian_Point, 
      coords::Representation_3D > cells{ xyz, 3.0, false,{}, 3.0 };
    IV ret;
    // cycle water molecules
    for (auto w : waters)
    { 
      // obtain geometric center of water molecule
      auto p = (xyz.at(molecules.at(w).at(0)) + xyz.at(molecules.at(w).at(1)) + xyz.at(molecules.at(w).at(2))) / 3.;
      bool exposed{ false };
      // cycle all directions around geomatric center
      for (double i = 0.0; i < 360.0; i += 4.0)
      {
        for (double j = 0.0; j < 360.0; j += 4.0)
        {
          auto inc = scon::ang<double>::from_deg(i);
          auto azi = scon::ang<double>::from_deg(j);
          auto sini = sin(inc);
          // get position 2.0 A from center in current direction
          coords::Cartesian_Point x{ 2.0*sini*cos(azi), 2.0*sini*sin(azi), 2.0*cos(inc) };
          // translate to make vector originate at p
          x += p;
          // get linked cell box of the position
          auto b = cells.box_of_point(x);
          exposed = true;
          // cycle all atoms in adjacent boxes to check for steric interference
          for (auto n : b.adjacencies())
          {
            if (n >= 0)
            { 
              auto a = static_cast<std::size_t>(n);
              // if n is valid atom index
              if (a >= xyz.size()) continue; 
              // if atom is not part of current water and nearer than 2.0 A
              // this direction is "burried"
              if (a != molecules.at(w).at(0) && a != molecules.at(w).at(1) && 
                a != molecules.at(w).at(2) && len(xyz[a] - x) < 2.0)
              {
                exposed = false;
                break;
              }
            }
          }
          // if current water is exposed in current direction we stop 
          // and add to return vector
          if (exposed)
          {
            ret.push_back(w);
            exposed = true;
            break;
          }
        }
        if (exposed)
        {
          break;
        }
      }
    }
    return ret;
  }

  struct windices
  {
    std::size_t h[2], o;
    windices() : h(), o() {}
  };

  template<class MI>
  windices get_windices(std::size_t const water_mol_index, 
    MI const & molecules, coords::Atoms const &a)
  {
    windices r;
    auto const &wm = molecules[water_mol_index];

    if (a.atom(wm[0]).number() == 1u)
    {
      r.h[0] = wm[0];
      if (a.atom(wm[1]).number() == 1u)
      {
        r.h[1] = wm[1];
        r.o = wm[2];
      }
      else
      {
        r.h[1] = wm[2];
        r.o = wm[1];
      }
    }
    else 
    {
      r.h[0] = wm[1];
      r.h[1] = wm[2];
      r.o = wm[0];
    }
    return r;
  }

}

std::size_t optimization::global::optimizers::monteCarlo::move_water(coords::Coordinates &c)
{
  auto const n = c.size();
  std::vector<bool> tabu(n, false);
  // linked cell object for adjacent atoms
  scon::linked::Cells < coords::float_type, coords::Cartesian_Point,
    coords::Representation_3D > cells{ c.xyz(), 3.0, false,{}, 3.0 };
  // create and screen sites
  auto sites = startopt::solvadd::accessible_sites(
    startopt::solvadd::build_hb_sites(c.xyz(), c.atoms(), tabu), c.xyz());
  if (sites.empty())
  {
    throw std::logic_error("Random water movement requires accessible hydrogen bonding sites.");
  }
  // shuffle sites
  std::mt19937_64 mte{ std::random_device{}() };
  std::shuffle(std::begin(sites), std::end(sites), mte);
  // get water all molecules
  std::vector<std::size_t> wmol;
  for (auto j : scon::index_range(c.molecules()))
  {
    if (molecule_is_moveable_water(c.molecules()[j], c.atoms()))
    { // if molecule is water we add its index to wmol
      // and its geometric center to wpos
      wmol.push_back(j);
    }
  }
  if (wmol.empty())
  {
    throw std::logic_error("Random water movement requires moveable "
      "(non fixed) water molecules to be present.");
  }
  // remove burried waters to obtain candidates for movement
  wmol = remove_burried_waters(wmol, c.molecules(), c.xyz());
  if (wmol.empty())
  {
    throw std::logic_error("Random water movement requires unburried "
      "(free 2 A sphere in 2 A from) water molecules to be present.");
  }
  // shuffle water molecules
  std::shuffle(std::begin(wmol), std::end(wmol), mte);
  // select number of water molecules to be moved
  auto m = std::geometric_distribution<std::size_t>{ 
    Config::get().optimization.global.montecarlo.move_frequency_probability 
  }(mte);
  m += 1;
  m = std::min<std::size_t>(m, wmol.size());
  // get water atoms to be potentially moved
  std::vector<std::size_t> potentially_moved_atoms;
  for (std::size_t j{ 0 }; j < m; ++j)
  {
    auto wat = wmol.at(j);
    scon::sorted::insert_unique(potentially_moved_atoms, c.molecules().at(wat).at(0));
    scon::sorted::insert_unique(potentially_moved_atoms, c.molecules().at(wat).at(1));
    scon::sorted::insert_unique(potentially_moved_atoms, c.molecules().at(wat).at(2));
  }
  // Print number of moved water molecules
  std::size_t actually_moved_molecules(0);
  for (std::size_t j{ 0 }; j < m; ++j)
  {
    for (auto & s : sites)
    {
      // if site is tabu or site is "part of" moved water 
      // go on without moving to this site
      if (s.tabu || scon::sorted::exists(potentially_moved_atoms, s.atom))
      {
        continue;
      }
      auto b = cells.box_of_point(s.p);
      auto f = startopt::solvadd::fit_water_into_site(s, 
      { b.adjacencies().begin(), b.adjacencies().end() }, c.xyz(), c.atoms());
      if (f.first)
      {
        auto wi = get_windices(wmol[j], c.molecules(), c.atoms());
        //std::cout << "((" << wi.o << "," << wi.h[0] << ", " << wi.h[1] << "))";
        //std::cout << "Moving " << wi.o << " at " << c.xyz()[wi.o] << " to " << f.second.o << "\n";
        //std::cout << "Moving " << wi.h[0] << " at " << c.xyz()[wi.h[0]] << " to " << f.second.h[0] << "\n";
        //std::cout << "Moving " << wi.h[1] << " at " << c.xyz()[wi.h[1]] << " to " << f.second.h[1] << "\n";
        c.move_atom_to(wi.o, f.second.o);
        c.move_atom_to(wi.h[0], f.second.h[0]);
        c.move_atom_to(wi.h[1], f.second.h[1]);
        s.tabu = true;
        ++actually_moved_molecules;
        break;
      }
    }
  }
  return actually_moved_molecules;
}
