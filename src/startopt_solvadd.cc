#include <cmath>
#include <ctime>
#include <string>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>

#include "atomic.h"
#include "startopt_solvadd.h"
#include "scon_chrono.h"
#include "scon_utility.h"
#include "histogram.h"
#include "optimization_global.h"
#include "coords_io.h"

bool startopt::solvadd::water::check_geometry (void) const
{
  using std::abs;
  coords::Cartesian_Point const b1(h[0] - o), b2(h[1] - o);
  // Config::get().startopt.solvadd.water_bond
  double const d1(std::abs(Config::get().startopt.solvadd.water_bond - len(b1)));
  double const d2(std::abs(Config::get().startopt.solvadd.water_bond - len(b2)));
  double const da(abs((Config::get().startopt.solvadd.water_angle - angle(b1, b2))).degrees());
  double const d3{ std::abs(len(h[0]-h[1])) };
  return (d1 > 0.05 || d2 > 0.05 || da > 0.05 || d3 < 1.0) ? false : true;
}

double ratio(std::size_t const num_atoms, std::size_t const num_water)
{
  auto const n_pp = (num_atoms*num_atoms - num_atoms) / 2;
  auto const n_ww = (9 * num_water*num_water - 3 * num_water) / 2;
  auto n_wp = num_atoms*num_water * 3;
  return (static_cast<double>(n_wp - n_ww) / static_cast<double>(n_pp + n_ww + n_wp));
}

std::size_t startopt::solvadd::num_w_max_w_solute_ia(std::size_t const a)
{
  std::size_t w(0u);
  auto rat = ratio(a, w);
  auto rat_next = ratio(a, w + 1u);
  while (rat_next > rat)
  {
    rat = rat_next;
    ++w;
    rat_next = ratio(a, w + 1u);
  }
  return w;
}

void startopt::preoptimizers::Solvadd::generate (
  coords::Ensemble_PES const & init_ensemble, std::size_t const multiplier)
{
  std::size_t const init_ensemble_size(init_ensemble.size());
  std::size_t const pre_per_structure(multiplier > 0u ? multiplier : 1u);
  auto const b_tmp = m_boundary;
  //std::cout << "Optimum waters: " << solvadd::wp_optimum(coords.size()) << std::endl;
  //Multihistogram<double> sahist(2U, 2.0, 0.0);
  std::size_t const twopps(100U*pre_per_structure);
  bool const SHELLOPT((Config::get().startopt.solvadd.opt == config::startopt_conf::solvadd::opt_types::SHELL || 
                      Config::get().startopt.solvadd.opt == config::startopt_conf::solvadd::opt_types::TOTAL_SHELL));
  for (std::size_t e =0u; e < init_ensemble_size; ++e)
  {
    std::size_t p(0U), pi(0U);
    while(p < pre_per_structure && pi < twopps)
    {
      m_boundary = b_tmp;
      std::size_t w(0u);
      m_sites.clear();
      tabu_atoms.clear();
      m_solvated_atoms = coords.atoms();
      m_interlinks.clear();
      m_solvated_positions = init_ensemble[e].structure.cartesian;
      m_solvated_cells.update();
      tabu_atoms.assign(coords.size(), false);
      if (Config::get().startopt.solvadd.fix_initial) m_solvated_atoms.fix_all();
      bool check(true);
      while (check)
      {
        std::size_t added_total(0U);
        build_sites();
        for (auto const & site : m_sites) 
        {
          tabu_atoms[site.atom] = true;
          if ((w+added_total) >= Config::get().startopt.solvadd.maxNumWater && 
            Config::get().startopt.solvadd.maxNumWater != 0U) break;
          if (populate_site(site)) ++added_total; 
        }
        std::size_t purged(0U);
        if (added_total > 0U && SHELLOPT)
        {
          // populate solvated_coords object with solvation data
          populate_coords(added_total);
          //// Set current ensemble positions
          //for (std::size_t a = 0; a < init_ensemble[e].size(); ++a)
          //{
          //  solvated_coords.move_atom_to(a, init_ensemble[e].structure.cartesian[a]);
          //}
          // optimize solvated_coords
          solvated_coords.o();
          // get optimized positions
          m_solvated_positions = solvated_coords.xyz();
        }
        // remove out of boundary waters
        purged = purge_coords();
        w += added_total;
        w -= purged;
        m_solvated_cells.update();
        if (SHELLOPT) 
        {
          tabu_atoms.resize(m_solvated_atoms.size());
          tabu_atoms.assign(m_solvated_atoms.size(), false);
        }
        if (purged >= added_total && w < Config::get().startopt.solvadd.maxNumWater) 
        {
          push_boundary(); 
          tabu_atoms.assign(m_solvated_atoms.size(), false);
          continue;
        } 
        check = !m_sites.empty() && added_total != 0 && purged < added_total && 
          (Config::get().startopt.solvadd.maxNumWater == 0U || 
            w < Config::get().startopt.solvadd.maxNumWater);
        if (Config::get().general.verbosity > 2U)
        {
          std::cout << "Waters: " << w << " (Added: " 
            << added_total << " , Purged: " << purged 
            << " , Diff: " << added_total - purged << ")" << std::endl;
        }
      }

      // we need to maintain a stable number of water molecules throughout the ensemble
      std::size_t const iter(e*pre_per_structure+p);
      if (iter > 0u && m_solvated_positions.size() != m_final_coords.size() && Config::get().general.verbosity > 2U)
      {
        auto currentNumberOfWaters = (m_solvated_positions.size() - coords.size()) / 3u;
        auto desiredNumberOfWaters = (m_final_coords.size() - coords.size()) / 3u;
        //auto d = ds > df ? ds-df : df-ds;
        std::cout << "Number of waters " << currentNumberOfWaters << " does not match required number " << desiredNumberOfWaters << " - retrying.\n";
      }

      if (iter < 1U || m_solvated_positions.size() == m_final_coords.size())
      {
        populate_coords(w); 
        if (Config::get().startopt.solvadd.opt == 
            config::startopt_conf::solvadd::opt_types::TOTAL_SHELL || 
            Config::get().startopt.solvadd.opt == 
            config::startopt_conf::solvadd::opt_types::TOTAL)
        {
          solvated_coords.o();
        }
        else
        {
          //solvated_coords.e();
        }
        solvated_coords.to_internal();
        solvated_coords.to_xyz();
        
        m_solvated_positions = solvated_coords.xyz();
        m_ensemble.push_back(solvated_coords.pes());
        if (Config::get().general.verbosity > 1U)
        {
          std::cout << "Solvation " << iter << " / " 
            << (init_ensemble_size*pre_per_structure) << " has "<< w 
            << " waters. Energy: " << solvated_coords.pes().energy << std::endl;
        }
        if (iter < 1U) 
        {
          m_final_coords = solvated_coords;
          if (Config::get().general.verbosity > 2U)
          {
            auto d = (m_final_coords.size() - this->coords.size()) / 3;
            std::cout << "First structure has " << d << " water molecules.\n";
          }
        }
        ++p;
      }
      ++pi;
    }
  }
  m_final_coords.swap(solvated_coords);
  //scon::chrono::high_resolution_progresstimer::duration pgtime = pgt.stopover();
  //std::cout << std::endl << "SA took " << scon::chrono::to_seconds(pgtime) << "s (" << pgtime.count() << " clock ticks) to complete." << std::endl;
  //sahist.distribute();
  //std::ofstream histogram_stream(std::string(Config::get().general.outputFilename).append("_SA_histogram.txt").c_str());
  //histogram_stream << sahist;
}

void startopt::preoptimizers::Solvadd::build_sites (void)
{
  std::size_t const N(m_solvated_atoms.size());
  m_sites.clear();
  for (std::size_t i(0u); i<N; ++i)
  {
    if (tabu_atoms[i]) continue;

    switch (m_solvated_atoms.atom(i).number())
    {
    case 1: 
      {
        build_site_group_1(i);
        break;
      }
    case 7:
    case 15:
      {
        build_site_group_15(i);
        break;
      }
    case 8:
    case 16:
      {
        build_site_group_16(i);
        break;
      }
    case 9:
    case 17:
    case 35:
    case 53:
      {
        build_site_group_17(i);
        break;
      }
    default: break;
    }
  }
  std::random_device rd;
  auto engine = std::default_random_engine{ rd() };
  std::shuffle(m_sites.begin(), m_sites.end(), engine);
}

void startopt::preoptimizers::Solvadd::build_site_group_1  (std::size_t atom)
{
  if(!m_solvated_atoms.atom(atom).bonds().empty()) 
  {
    std::size_t bound = m_solvated_atoms.atom(atom).bonds(0u), an = m_solvated_atoms.atom(bound).number();
    if (an == 7u || an == 8u || an == 15u || an == 16u)
    {
      solvadd::site s;
      s.donor = true;
      s.v = m_solvated_positions[atom] - m_solvated_positions[bound];
      s.v *= Config::get().startopt.solvadd.defaultLenHB / geometric_length(s.v);
      //s.v.resize(Config::get().startopt.solvadd.defaultLenHB);
      s.p = m_solvated_positions[atom] + s.v;
      s.atom = atom;
      m_sites.push_back(s);
    }
  }
}

void startopt::preoptimizers::Solvadd::build_site_group_15 (std::size_t atom)
{
  config::startopt_conf::solvadd const & soc = Config::get().startopt.solvadd;
  coords::size_1d const & b(m_solvated_atoms.atom(atom).bonds());
  std::size_t nb(b.size());
  coords::Representation_3D const & p = m_solvated_positions;
  if(nb == 3) 
  {
    solvadd::site s;
    //! cross two Vectors in the plane of the three ligands for perpendicular vector
    coords::Cartesian_Point perpendicular = normalized(cross(p[b[1]] - p[b[0]], p[b[2]] - p[b[0]]));
    //! set length to default hydrogen-bond length + water bond (reaching Oxygen)
    perpendicular *= soc.defaultLenHB + soc.water_bond;
    //! compute relative positions "Atom+p" and "Atom-p"
    coords::Cartesian_Point relPos[2] = { p[atom] + perpendicular, p[atom] - perpendicular };
    //! get distances between these relative positions and the bonded atoms A, B, C
    double dist1[3] = {distance(relPos[0], p[b[0]]),
                       distance(relPos[0], p[b[1]]),
                       distance(relPos[0], p[b[2]])},
           dist2[3] = {distance(relPos[1], p[b[0]]),
                       distance(relPos[1], p[b[1]]),
                       distance(relPos[1], p[b[2]])},
           meandist[2] = {(dist1[0]+dist1[1]+dist1[2])/3.0,
                          (dist2[0]+dist2[1]+dist2[2])/3.0},
           meandiff = meandist[0] - meandist[1],
           abmd = abs(meandiff);
    /*! if both rel posis have nearly the same mean distance to bonded atoms -> sp2
     * -> no electron-pair donation
     * -> no extension */
    if(abmd < 0.02) return;
    //! if meandiff is smaller than 0 we need the relPos[1] combination
    else if(meandiff < 0)
    {
      s.v = -perpendicular;
      s.p = relPos[1];
    }
    //! if meandiff > 0 the relPos[0] is desired position
    else
    {
      s.v = perpendicular;
      s.p = relPos[0];
    }
    //! Store the extension
    s.atom = atom;
    s.tabu = false;
    m_sites.push_back(s);
  } 
  else if(nb == 2) 
  {
    std::size_t idx[3] = {atom, b[0], b[1]};
    build_site(idx, coords::angle_type::from_deg(120), coords::angle_type::from_deg(180), true);
  }
}

void startopt::preoptimizers::Solvadd::build_site_group_16 (std::size_t atom)
{
  std::size_t nb(m_solvated_atoms.atom(atom).bonds().size());
  if (nb > 0u)
  {
    std::size_t ba(m_solvated_atoms.atom(atom).bonds(0u)), idx[3] = {atom, ba, 0u};
    if(nb == 1u)
    {
      std::size_t batomic(m_solvated_atoms.atom(ba).number()), nbb(m_solvated_atoms.atom(ba).bonds().size());
      if (nbb > 0u)
      {
        idx[2] = (m_solvated_atoms.atom(ba).bonds(0u) == atom) ? 
          m_solvated_atoms.atom(ba).bonds(1u) : m_solvated_atoms.atom(ba).bonds(0u);
        //! SP2 C or N --> Carbonyl O --> two
        if((nbb == 3u && batomic == 6u) || (nbb == 2u && batomic == 7u))
          build_multisite(idx, 2u, coords::angle_type::from_deg(120.0), coords::angle_type::from_deg(0.0), 180.0, false);
        //! SP3 C --> O-
        else if (nbb == 4u && batomic == 6u)
          build_multisite(idx, 3u, coords::angle_type::from_deg(109.5), coords::angle_type::from_deg(60.0), 120.0, false);
      }
    }
    //! O bound to two Atoms: tetraedric)
    else if (nb == 2u)
    {
      idx[2] = m_solvated_atoms.atom(atom).bonds(1u);
      build_multisite(idx, 2u, coords::angle_type::from_deg(109.5), coords::angle_type::from_deg(120.0), 120.0, true);
    }
  }
}

void startopt::preoptimizers::Solvadd::build_site_group_17 (std::size_t atom)
{
  if(m_solvated_atoms.atom(atom).bonds().size() > 1) return;
  std::size_t b(m_solvated_atoms.atom(atom).bonds(0u));
  if(m_solvated_atoms.atom(b).bonds().size() < 2u) return;
  std::size_t b0(m_solvated_atoms.atom(b).bonds(0u)), b1(m_solvated_atoms.atom(b).bonds(1u));
  std::size_t idx[3] = {atom, b, (b0 == atom ? b1 : b0)};
  build_multisite(idx, 3, coords::angle_type::from_deg(109.5), coords::angle_type::from_deg(60.0), 120.0, false);
}

void startopt::preoptimizers::Solvadd::build_site(std::size_t const index[3], 
  coords::angle_type const angle, coords::angle_type const dihedral, bool const tetraedric)
{
  solvadd::site s;
  s.atom = index[0];
  double const distance = Config::get().startopt.solvadd.defaultLenHB + Config::get().startopt.solvadd.water_bond;
  coords::Cartesian_Point const
    a(m_solvated_positions[index[2U]]), b(m_solvated_positions[index[1U]]), c(m_solvated_positions[index[0U]]), 
    ab(tetraedric ? c-b : b-a), cb(tetraedric ? a-c : b-c);
  s.p = appendNERF(c, ab, cb, distance, angle, dihedral);
  s.v = s.p - c;
  s.tabu = false;
  m_sites.push_back(s);
}

void startopt::preoptimizers::Solvadd::build_multisite(std::size_t const index[3], 
  std::size_t const n_sites, coords::angle_type const angle, 
  coords::angle_type const dihedral, coords::float_type const offset, bool const tetraedric)
{
  for (std::size_t i = 0; i < n_sites; ++i)
  {
    build_site(index, angle, dihedral + coords::angle_type::from_deg(i*offset), tetraedric);
  }
}

static coords::Cartesian_Point center_of_watermass 
(
  coords::Cartesian_Point const &o, 
  coords::Cartesian_Point const &h1, 
  coords::Cartesian_Point const &h2
)
{
  return coords::Cartesian_Point((o*15.9996+h1*1.00079+h2*1.00079)/18.00118);
}

bool startopt::preoptimizers::Solvadd::populate_site (solvadd::site const &s)
{
  using scon::randomized;
  //std::cout << "Populating site around " << solvated_coords.xyz(s.atom)+s.v << " box: " << m_solvated_cells.boxN(solvated_coords.xyz(s.atom)+s.v) << std::endl;
  cells_type::box_type const site_box(m_solvated_cells.box_of_point(m_solvated_positions[s.atom] + s.v));
  //for (auto const & atom_id : m_solvated_cells.box_of_element(12).adjacencies())
  //{
  //  atom_id 
  //}

  std::vector<std::size_t> const surrounding_atoms(site_box.adjacencies().begin(), site_box.adjacencies().end());
  coords::float_type const svl(len(s.v));
  //std::cout << "Sorroundings for populating:" << surrounding_atoms.size() << std::endl;
  coords::Cartesian_Point const b(-s.v);
  //static std::string n("x");
  for (std::size_t l(0u); l<5u; ++l)
  {
    solvadd::water w;
    if (m_solvated_atoms.atom(s.atom).number() == 1u)
    {
      w.o = m_solvated_positions[s.atom] + (s.v*(Config::get().startopt.solvadd.defaultLenHB + static_cast<coords::float_type>(l)*0.2)/svl);
    }
    else
    {
      w.h[0] = m_solvated_positions[s.atom] + (s.v*(Config::get().startopt.solvadd.defaultLenHB + static_cast<coords::float_type>(l)*0.2)/svl);
      w.o = w.h[0] + (s.v*(Config::get().startopt.solvadd.water_bond/svl));
    }
    for (std::size_t a(0u); a<360u; a+=20u)
    {
      if (m_solvated_atoms.atom(s.atom).number() == 1u)
      {

        w.h[0] = appendNERF(w.o, randomized<coords::Cartesian_Point>(), b, Config::get().startopt.solvadd.water_bond,
                                Config::get().startopt.solvadd.water_angle, coords::angle_type::from_deg(a));
        w.h[1] = appendNERF(w.o, -b, w.h[0]-w.o, Config::get().startopt.solvadd.water_bond, 
                                Config::get().startopt.solvadd.water_angle, coords::angle_type::from_deg(120.0));
      }
      else
      {
        w.h[1] = appendNERF(w.o, randomized<coords::Cartesian_Point>(), b, Config::get().startopt.solvadd.water_bond,
                                Config::get().startopt.solvadd.water_angle, coords::angle_type::from_deg(a));
      }
      //std::cout << "From " << m_solvated_atoms.size() << "\n";
      if (check_sterics(w, surrounding_atoms) && !check_out_of_boundary(center_of_watermass(w.o, w.h[0], w.h[1])) && w.check_geometry()) 
      {
        if (m_solvated_atoms.atom(s.atom).number() == 1u)
        {
          auto p = std::pair<std::size_t, std::size_t>(m_solvated_atoms.size() + 2, s.atom);
          m_interlinks.insert(p);
        }
        else
        {
          auto p = std::pair<std::size_t, std::size_t>(m_solvated_atoms.size(), s.atom);
          m_interlinks.insert(p);
        }
        add_water(w);
        //std::cout << "Added at " << scon::comma_delimeted(w.o, w.h[0], w.h[1]);
        return true;
      }
    }
  }
  return false;
}

bool startopt::preoptimizers::Solvadd::check_sterics (solvadd::water const &w, std::vector<std::size_t> const &atoms_around)  const
{
  //std::cout << "Sterics check with " << atoms_around.size() << " atoms around." << std::endl;
  // atoms_around
  for (auto const atom : atoms_around)
  {
    std::size_t const an(m_solvated_atoms.atom(atom).number());
    coords::float_type dist[3] = 
    {
      len(m_solvated_positions[atom] - w.h[0]), 
      len(m_solvated_positions[atom] - w.h[1]),
      len(m_solvated_positions[atom] - w.o)
    };
    //std::cout << "Distances to " << atom << ", " << m_solvated_positions[atom] 
    // << ": " << scon::comma_delimeted(dist[0], dist[1], dist[2]) << "\n";
    if (an == 1u && (dist[0] < 1.9 || dist[1] < 1.9 || dist[2] < 1.8)) 
      return false;
    else if (atomic::number_is_heteroatom(an) && (dist[0] < 1.7 || dist[1] < 1.7 || dist[2] < 2.4)) 
      return false;    
    else if (dist[0] < 1.9 || dist[1] < 1.9 || dist[2] < 2.2) 
      return false;
  }
  return true;
}

bool startopt::preoptimizers::Solvadd::check_out_of_boundary (coords::Cartesian_Point const & p) const
{
  switch (Config::get().startopt.solvadd.boundary)
  {
  case (config::startopt_conf::solvadd::boundary_types::BOX) :
    {
      if (p.x() > m_box_max.x() || p.y() > m_box_max.y() || p.z() > m_box_max.z()) return true;
      else if (p.x() < m_box_min.x() || p.y() < m_box_min.y() || p.z() < m_box_min.z()) return true;
      else return false;
    }
  case (config::startopt_conf::solvadd::boundary_types::SPHERE) :
    {
      if (len(p-m_init_center) > m_boundary) return true;
      else return false;
    }
  case (config::startopt_conf::solvadd::boundary_types::LAYER) :
    {
      cells_type::box_type const box(m_init_cells.box_of_point(p));
      for (auto atom : box.adjacencies())
      {
        if (len(coords.xyz(atom)-p) < m_boundary) return false;
      }
      return true;
    }
  }
  return true;
}

void startopt::preoptimizers::Solvadd::push_boundary ()
{
  m_boundary += 0.5;
  m_box_max += 0.25;
  m_box_min -= 0.25;
  if (Config::get().general.verbosity > 2U)
  {
    std::cout << "Pushing boundary to " << m_boundary << "\n";
  }
  m_init_cells.init(m_boundary);
}

void startopt::preoptimizers::Solvadd::add_water(solvadd::water const &w)
{
  // Set new energy
  coords::Atom h1(static_cast<std::size_t>(1u)), h2(static_cast<std::size_t>(1u)), o(static_cast<std::size_t>(8u));
  h1.set_energy_type(Config::get().startopt.solvadd.ffTypeHydrogen);
  h2.set_energy_type(Config::get().startopt.solvadd.ffTypeHydrogen);
  o.set_energy_type(Config::get().startopt.solvadd.ffTypeOxygen);
  // Set bonds
  std::size_t const atom_count(m_solvated_atoms.size());
  h1.bind_to(atom_count+2);
  h2.bind_to(atom_count+2);
  o.bind_to(atom_count);
  o.bind_to(atom_count+1);
  // Add atoms to atom list and representation
  m_solvated_atoms.add(h1);
  m_solvated_atoms.add(h2);
  m_solvated_atoms.add(o);
  m_solvated_positions.push_back(w.h[0u]);
  m_solvated_positions.push_back(w.h[1u]);
  m_solvated_positions.push_back(w.o);
  //std::cout << "SV " << scon::vector_delimeter('\n') << m_solvated_positions << "\n";
  m_solvated_cells.update();
  tabu_atoms.reserve(tabu_atoms.size()+3u);
  tabu_atoms.push_back(false);
  tabu_atoms.push_back(false);
  tabu_atoms.push_back(false);
}

void startopt::preoptimizers::Solvadd::populate_coords (std::size_t const added)
{
  solvated_coords.clear();
  //coords::Coordinates tc;
  if (m_solvated_atoms.size() != m_solvated_positions.size()) 
    throw std::logic_error("Atoms size does not equal representation size.");
  std::size_t const N(m_solvated_atoms.size()), M(N-added*3U);
  for (std::size_t i(0U); i<coords.size(); ++i)
  {
    m_solvated_atoms.atom(i).assign_to_system(0u);
    // or if i is not added this turn and we fix intermediate
    // insert coord position 
    // m_solvated_positions[i] = coords.xyz(i);
    // fix if we fix init and i is part of init 
    m_solvated_atoms.atom(i).fix(Config::get().startopt.solvadd.fix_initial);
  }
  for (std::size_t i(coords.size()); i<N; ++i)
  {
    m_solvated_atoms.atom(i).assign_to_system(1U);
    // fix if i is not added this turn and we fix intermediate
    m_solvated_atoms.atom(i).fix(Config::get().startopt.solvadd.fix_intermediate && i < M);
  }
  //if (Config::get().startopt.solvadd.intern_connect_waters)
  //  Config::set().coords.internal.connect.swap(m_interlinks);
  solvated_coords.init_in(m_solvated_atoms, coords::PES_Point(m_solvated_positions), true);
  //if (Config::get().startopt.solvadd.intern_connect_waters)
  //  Config::set().coords.internal.connect.swap(m_interlinks);
}

std::size_t startopt::preoptimizers::Solvadd::purge_coords (void)
{
  std::size_t const N(m_solvated_atoms.size());
  if (N != m_solvated_positions.size()) 
    throw std::logic_error("Atoms size does not equal representation size.");
  std::size_t removed_w(0U);
  if (N > coords.size()) 
  {
    for (std::size_t i(coords.size()); i < m_solvated_positions.size(); ++i)
    {
      std::size_t o(0U), h[2] = {0,0};
      if (m_solvated_atoms.atom(i).number() == 8U && i > 1 && 
        m_solvated_atoms.atom(i-1).number() == 1U && 
        m_solvated_atoms.atom(i-2).number() == 1U)
      {
        o = i;
        h[1] = i-1;
        h[0] = i-2;
      }
      else if (m_solvated_atoms.atom(i).number() == 1U && 
              (i+1) < m_solvated_atoms.size() && 
              m_solvated_atoms.atom(i+1).number() == 8U)
      {
        o = i+1;
        h[1] = i;
        h[0] = i-1;
      }
      else if (m_solvated_atoms.atom(i).number() == 1U && 
              (i+2) < m_solvated_atoms.size() && 
              m_solvated_atoms.atom(i+2).number() == 8U)
      {
        o = i+2;
        h[1] = i+1;
        h[0] = i;
      }
      else throw std::logic_error("Unexpected atomic number in Solvation process.");
      if (check_out_of_boundary(center_of_watermass(m_solvated_positions[o], m_solvated_positions[h[0]], m_solvated_positions[h[1]])))
      { // if center of mass is out of bounds we remove the water with those positions
        // swap with last water and remove last water to avoid resort
        if ((*(m_solvated_atoms.end()-1)).number() == 8U &&
            (*(m_solvated_atoms.end()-2)).number() == 1U &&
            (*(m_solvated_atoms.end()-3)).number() == 1U)
        {
          ++removed_w;
          // preserve positions for last water in selected water
          m_solvated_positions[o] = *(m_solvated_positions.end()-1);
          m_solvated_positions[h[1]] = *(m_solvated_positions.end()-2);
          m_solvated_positions[h[0]] = *(m_solvated_positions.end()-3);
          // remove last water 
          // (last water in fact survives, out of bounds water removed)
          m_solvated_positions.pop_back();
          m_solvated_positions.pop_back();
          m_solvated_positions.pop_back();
          m_solvated_atoms.pop_back();
          m_solvated_atoms.pop_back();
          m_solvated_atoms.pop_back();
        }
        else 
        {
          throw std::logic_error("Unexpected atomic numbers in Solvation process.");
        }
      }
    }
    std::size_t w((m_solvated_positions.size() - coords.size())/3);
    while (Config::get().startopt.solvadd.maxNumWater != 0 
      && w > Config::get().startopt.solvadd.maxNumWater)
    {
      ++removed_w;
      m_solvated_positions.pop_back();
      m_solvated_positions.pop_back();
      m_solvated_positions.pop_back();
      m_solvated_atoms.pop_back();
      m_solvated_atoms.pop_back();
      m_solvated_atoms.pop_back();
      --w;
    }
  }
  return removed_w;
}


/*

  GOSOL

*/



std::size_t startopt::preoptimizers::GOSol::f_solvate(std::size_t const max_num_w)
{
  startopt::preoptimizers::Solvadd sap(m_coords);
  std::size_t const current_w((m_coords.size() - m_num_init_atoms) / 3),
    missing_w(max_num_w - current_w);
  Config::set().startopt.solvadd.maxNumWater = std::min(missing_w, static_cast<std::size_t>(solvadd::wp_optimum(m_coords.size())));
  sap.generate(m_pes, current_w > 0 ? 1U : Config::get().startopt.number_of_structures / m_pes.size());
  //std::cout << "Coord fin: " << sap.m_coo();
  m_coords = sap.final_coords();
  m_pes = sap.PES();
  std::size_t const W = (m_coords.size() - m_num_init_atoms) / 3;
  std::stringstream sstrm;
  sstrm << W;
  std::ofstream gstream((Config::get().general.outputFilename + "_" + sstrm.str() + "_SA.xyz").c_str());
  for (auto const & pes : m_pes)
  {
    m_coords.set_pes(pes);
    gstream << coords::output::formats::tinker(m_coords);
  }
  return W;
}

void startopt::preoptimizers::GOSol::f_optimize(std::string const &suffix)
{
  m_coords.fix_all(false);
  if (!Config::get().coords.fixed.empty())
  {
    for (auto fix : Config::get().coords.fixed)
    {
      if (fix < m_coords.size()) m_coords.fix(fix);
    }
  }
  if (Config::get().startopt.solvadd.go_type == config::globopt_routine_type::BASINHOPPING)
  {
    optimization::global::optimizers::monteCarlo mc(m_coords, m_pes, suffix);
    mc.run(Config::get().optimization.global.iterations);
    m_pes = coords::Ensemble_PES(mc.accepted_minima.begin(), mc.accepted_minima.end());
    mc.write_range(suffix);
  }
  else if (Config::get().startopt.solvadd.go_type == config::globopt_routine_type::TABUSEARCH)
  {
    optimization::global::optimizers::tabuSearch gots(m_coords, m_pes, suffix);
    gots.run(Config::get().optimization.global.iterations);
    m_pes = coords::Ensemble_PES(gots.accepted_minima.begin(), gots.accepted_minima.end());
    gots.write_range(suffix);
  }
}

void startopt::preoptimizers::GOSol::run(std::size_t const num_water)
{
  std::size_t W(0u);
  while (W < num_water) 
  {
    W = f_solvate(num_water);
    std::ostringstream S;
    S << "_" << W;
    f_optimize(S.str());
  }
}

// ...

namespace
{

  startopt::solvadd::site site(coords::Representation_3D const &xyz, 
    std::size_t const index[3], coords::angle_type const angle, 
    coords::angle_type const dihedral, bool const tetraedric)
  {
    startopt::solvadd::site s;
    s.atom = index[0];
    double const distance = Config::get().startopt.solvadd.defaultLenHB + Config::get().startopt.solvadd.water_bond;
    coords::Cartesian_Point const
      a(xyz[index[2U]]), b(xyz[index[1U]]), c(xyz[index[0U]]),
      ab(tetraedric ? c - b : b - a), cb(tetraedric ? a - c : b - c);
    s.p = appendNERF(c, ab, cb, distance, angle, dihedral);
    s.v = s.p - c;
    s.tabu = false;
    return s;
  }

  void multi_sites(std::vector<startopt::solvadd::site> &sites, 
    coords::Representation_3D const &xyz,
    std::size_t const index[3], std::size_t const n_sites, 
    coords::angle_type const angle, coords::angle_type const dihedral, 
    coords::float_type const offset, bool const tetraedric)
  {
    for (std::size_t i = 0; i < n_sites; ++i)
    {
      sites.push_back(site(xyz, index, angle, 
        dihedral + coords::angle_type::from_deg(i*offset), tetraedric));
    }
  }

  void site_group_1(std::vector<startopt::solvadd::site> &sites,
    coords::Representation_3D const &xyz, 
    coords::Atoms const &atms, 
    std::size_t const i)
  {
    if (!atms.atom(i).bonds().empty())
    {
      auto b = atms.atom(i).bonds(0u);
      auto an = atms.atom(b).number();
      if (an == 7u || an == 8u || an == 15u || an == 16u)
      {
        startopt::solvadd::site s;
        s.v = xyz[i] - xyz[b];
        s.v *= Config::get().startopt.solvadd.defaultLenHB / len(s.v);
        s.p = xyz[i] + s.v;
        s.atom = i;
        s.donor = true;
        sites.push_back(s);
      }
    }
  }

  void site_group_15(std::vector<startopt::solvadd::site> &sites,
    coords::Representation_3D const &xyz,
    coords::Atoms const &atms,
    std::size_t const i)
  {
    config::startopt_conf::solvadd const & soc = Config::get().startopt.solvadd;
    coords::size_1d const & b(atms.atom(i).bonds());
    std::size_t nb(b.size());
    coords::Representation_3D const & p = xyz;
    if (nb == 3)
    {
      startopt::solvadd::site s;
      //! cross two Vectors in the plane of the three ligands for perpendicular vector
      auto perpendicular = normalized(cross(p[b[1]] - p[b[0]], p[b[2]] - p[b[0]]));
      //! set length to default hydrogen-bond length + water bond (reaching Oxygen)
      perpendicular *= soc.defaultLenHB + soc.water_bond;
      //! compute relative positions "Atom+p" and "Atom-p"
      coords::Cartesian_Point relPos[2] = { p[i] + perpendicular, p[i] - perpendicular };
      //! get distances between these relative positions and the bonded atoms A, B, C
      double dist1[3] = { distance(relPos[0], p[b[0]]),
        distance(relPos[0], p[b[1]]),
        distance(relPos[0], p[b[2]]) },
        dist2[3] = { distance(relPos[1], p[b[0]]),
        distance(relPos[1], p[b[1]]),
        distance(relPos[1], p[b[2]]) },
        meandist[2] = { (dist1[0] + dist1[1] + dist1[2]) / 3.0,
        (dist2[0] + dist2[1] + dist2[2]) / 3.0 },
        meandiff = meandist[0] - meandist[1],
        abmd = abs(meandiff);
      /*! if both rel posis have nearly the same mean distance to bonded atoms -> sp2
      * -> no electron-pair donation
      * -> no extension */
      if (abmd < 0.02) return;
      //! if meandiff is smaller than 0 we need the relPos[1] combination
      else if (meandiff < 0)
      {
        s.v = -perpendicular;
        s.p = relPos[1];
      }
      //! if meandiff > 0 the relPos[0] is desired position
      else
      {
        s.v = perpendicular;
        s.p = relPos[0];
      }
      //! Store the extension
      s.atom = i;
      s.tabu = false;
      sites.push_back(s);
    }
    else if (nb == 2)
    {
      std::size_t idx[3] = { i, b[0], b[1] };
      sites.push_back(site(xyz, idx, coords::angle_type::from_deg(120), coords::angle_type::from_deg(180), true));
    }
  }

  void site_group_16(std::vector<startopt::solvadd::site> &sites,
    coords::Representation_3D const &xyz,
    coords::Atoms const &atms, std::size_t const atom)
  {
    std::size_t nb(atms.atom(atom).bonds().size());
    if (nb > 0u)
    {
      std::size_t ba(atms.atom(atom).bonds(0u)), idx[3] = { atom, ba, 0u };
      if (nb == 1u)
      {
        std::size_t batomic(atms.atom(ba).number()), nbb(atms.atom(ba).bonds().size());
        if (nbb > 0u)
        {
          idx[2] = (atms.atom(ba).bonds(0u) == atom) ?
            atms.atom(ba).bonds(1u) : atms.atom(ba).bonds(0u);
          //! SP2 C or N --> Carbonyl O --> two
          if ((nbb == 3u && batomic == 6u) || (nbb == 2u && batomic == 7u))
            multi_sites(sites, xyz, idx, 2u, coords::angle_type::from_deg(120.0), 
              coords::angle_type::from_deg(0.0), 180.0, false);
          //! SP3 C --> O-
          else if (nbb == 4u && batomic == 6u)
            multi_sites(sites, xyz, idx, 3u, coords::angle_type::from_deg(109.5),
              coords::angle_type::from_deg(60.0), 120.0, false);
        }
      }
      //! O bound to two Atoms: tetraedric)
      else if (nb == 2u)
      {
        idx[2] = atms.atom(atom).bonds(1u);
        multi_sites(sites, xyz, idx, 2u, coords::angle_type::from_deg(109.5), 
          coords::angle_type::from_deg(120.0), 120.0, true);
      }
    }
  }

  void site_group_17(std::vector<startopt::solvadd::site> &sites,
    coords::Representation_3D const &xyz,
    coords::Atoms const &atms,
    std::size_t const atom)
  {
    if (atms.atom(atom).bonds().size() > 1) return;
    std::size_t b(atms.atom(atom).bonds(0u));
    if (atms.atom(b).bonds().size() < 2u) return;
    std::size_t b0(atms.atom(b).bonds(0u)), b1(atms.atom(b).bonds(1u));
    std::size_t idx[3] = { atom, b, (b0 == atom ? b1 : b0) };
    multi_sites(sites, xyz, idx, 3, coords::angle_type::from_deg(109.5), coords::angle_type::from_deg(60.0), 120.0, false);
  }

  bool has_enough_space(startopt::solvadd::water const &w,
    std::vector<std::size_t> const &atoms_around, 
    coords::Representation_3D const &xyz,
    coords::Atoms const &atms)
  {
    for (auto const a : atoms_around)
    {
      std::size_t const an(atms.atom(a).number());
      coords::float_type dist[3] =
      {
        len(xyz[a] - w.h[0]),
        len(xyz[a] - w.h[1]),
        len(xyz[a] - w.o)
      };
      if (an == 1u && (dist[0] < 1.7 || dist[1] < 1.7 || dist[2] < 1.6)) return false;
      else if (atomic::number_is_heteroatom(an) && (dist[0] < 1.5 || dist[1] < 1.5 || dist[2] < 2.2)) return false;
      else if (dist[0] < 1.8 || dist[1] < 1.8 || dist[2] < 2.0) return false;
    }
    return true;
  }
}

std::vector<startopt::solvadd::site> startopt::solvadd::build_hb_sites(
  coords::Representation_3D const &xyz, coords::Atoms const &atms,
   std::vector<bool> const &tabu)
{
  auto const n = xyz.size();
  if (n != atms.size())
  {
    throw std::logic_error("Atoms size does not match coordinates size in build_hb_sites.");
  }
  if (n != tabu.size())
  {
    throw std::logic_error("Tabu size does not match coordinates size in build_hb_sites.");
  }
  std::vector<startopt::solvadd::site> sites;
  for (std::size_t i{ 0 }; i < n; ++i)
  {
    if (tabu[i]) continue;
    switch (atms.atom(i).number())
    {
      case 1:
      {
        site_group_1(sites, xyz, atms, i);
        //build_site_group_1(i);
        break;
      }
      case 7:
      case 15:
      {
        site_group_15(sites, xyz, atms, i);
        //build_site_group_15(i);
        break;
      }
      case 8:
      case 16:
      {
        site_group_16(sites, xyz, atms, i);
        //build_site_group_16(i);
        break;
      }
      case 9:
      case 17:
      case 35:
      case 53:
      {
        site_group_17(sites, xyz, atms, i);
        //build_site_group_17(i);
        break;
      }
      default: break;
    }
  }
  return sites;
}

std::vector<startopt::solvadd::site> startopt::solvadd::accessible_sites(
  std::vector<startopt::solvadd::site> const & sites, coords::Representation_3D const & xyz)
{
  scon::linked::Cells < coords::float_type,
    coords::Cartesian_Point, coords::Representation_3D > cells{ xyz, 3.0, false, {}, 3.0 };
  std::vector<startopt::solvadd::site> ret;
  for (auto const &site : sites)
  {
    auto b = cells.box_of_point(site.p);
    bool accessible{ true };
    for (auto n : b.adjacencies())
    {
      if (n >= 0)
      {
        auto i = static_cast<std::size_t>(n);
        if (i != site.atom && len(xyz[i] - site.p) < 2.0)
        {
          accessible = false;
          break;
        }
      }
    }
    if (accessible)
    {
      ret.push_back(site);
    }
  }
  return ret;
}

std::pair<bool, startopt::solvadd::water> 
  startopt::solvadd::fit_water_into_site(startopt::solvadd::site const & s,
    std::vector<std::size_t> const &surrounding_atoms,
    coords::Representation_3D const & xyz,
    coords::Atoms const &atms)
{
  using scon::randomized;
  coords::float_type const svl(len(s.v));
  //std::cout << "Sorroundings for populating:" << surrounding_atoms.size() << std::endl;
  coords::Cartesian_Point const b(-s.v);
  //static std::string n("x");
  for (coords::float_type l(0.0); l < 1.2; l += coords::float_type(0.2))
  {
    startopt::solvadd::water w;
    if (s.donor)
    {
      w.o = xyz[s.atom] + (s.v*(Config::get().startopt.solvadd.defaultLenHB + l) / svl);
    }
    else
    {
      w.h[0] = xyz[s.atom] + (s.v*(Config::get().startopt.solvadd.defaultLenHB + l) / svl);
      w.o = w.h[0] + (s.v*(Config::get().startopt.solvadd.water_bond / svl));
    }
    for (std::size_t a(0u); a < 360u; a += 20u)
    {
      if (s.donor)
      {

        w.h[0] = appendNERF(w.o, randomized<coords::Cartesian_Point>(), b, Config::get().startopt.solvadd.water_bond,
          Config::get().startopt.solvadd.water_angle, coords::angle_type::from_deg(a));
        w.h[1] = appendNERF(w.o, -b, w.h[0] - w.o, Config::get().startopt.solvadd.water_bond,
          Config::get().startopt.solvadd.water_angle, coords::angle_type::from_deg(120.0));
      }
      else
      {
        w.h[1] = appendNERF(w.o, randomized<coords::Cartesian_Point>(), b, Config::get().startopt.solvadd.water_bond,
          Config::get().startopt.solvadd.water_angle, coords::angle_type::from_deg(a));
      }
      if (has_enough_space(w, surrounding_atoms, xyz, atms) && w.check_geometry())
      {
        return {true, w};
      }
    }
  }
  return{ false, water{} };
}

