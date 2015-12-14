#include <vector>
#include <string>
#include <stdexcept>
#include <random>
#include "configuration.h"
#include "atomic.h"
#include "scon_vect.h"
#include "evolution.h"
#include "coords.h"
#include "startopt_ringsearch.h"
#include "configuration.h"
#include "scon_utility.h"


startopt::ringsearch::ring::ring (std::size_t const _s) 
  : size(_s), atoms(_s), dihedrals(_s - 3), dihedral_direction(_s - 3)
{ }

startopt::ringsearch::ring::ring (std::size_t const _s, std::vector<std::size_t> _v) : 
  size(_s), atoms(_v.begin(), _v.begin()+size), dihedrals(size-3), dihedral_direction(size-3) 
{ }

startopt::ringsearch::ring::ring (ring const &r) : 
  size(r.size), atoms(r.atoms), dihedrals(r.dihedrals), 
  dihedral_direction(r.dihedral_direction) 
{ }

startopt::ringsearch::ring & startopt::ringsearch::ring::operator= (startopt::ringsearch::ring const & r) 
{
  if (size==r.size)
  {
    atoms = r.atoms;
    dihedrals = r.dihedrals;
    dihedral_direction = r.dihedral_direction;
  }
  return *this;
}

startopt::ringsearch::Search::Search (coords::Coordinates & coords) : 
  m_init_xyz(coords.xyz()),
  m_final_ensemble(),
  m_coord(coords),
  m_is_acceptor(coords.size(), false),
  m_rings(), 
  m_overlap(),
  m_ringcontainer(7u)
{
  std::size_t const n = m_coord.size();
  for (std::size_t i(0U); i < n; ++i)
  {
    std::size_t const atomic_number_i = m_coord.atoms(i).number();
    if (atomic_number_i == 1)
    {
      if (m_coord.atoms(i).bonds().size() == 1 &&
          atomic::number_is_heteroatom(m_coord.atoms(m_coord.atoms(i).bonds(0)).number()))
      {
        m_donors.push_back(i);
      }
    }
    else if (atomic::number_is_heteroatom(atomic_number_i))
    {
      m_is_acceptor[i] = true;
    }
  }
  find_rings();
}

startopt::ringsearch::Search::~Search()
{
  m_coord.set_xyz(m_init_xyz, true);
  m_coord.to_internal();
}


std::size_t startopt::ringsearch::Search::find_rings(void)
{
  std::size_t const n = m_coord.size();
  std::vector<bool> tabu_vector(n, false);
  std::size_t const m = m_donors.size();
  for (std::size_t i = 0; i<m; ++i)
  {
    m_ringcontainer[0u] = m_donors[i];
    find_route(m_donors[i], 0u, std::vector<bool>(n, false));
  }
  if (Config::get().general.verbosity > 19)
  {
    std::cout << "Found " << m_rings.size() << " ring-routes.\n";
    for (auto const & r : m_rings)
    {
      std::cout << "[" << r.size << "]: ";
      for (auto const & a : r.atoms)
      {
        std::cout << a << " ";
      }
      std::cout << "\n";
    }
  }
  find_torsions();
  find_overlap();

  return m_rings.size();
}

void startopt::ringsearch::Search::find_route(std::size_t const i, std::size_t const size, std::vector<bool> tabulist)
{
  if (size > 5) return;
  tabulist[i] = true;
  const std::size_t increased = size + 1;
  for (auto const bound : m_coord.atoms(i).bonds())
  {
    m_ringcontainer[increased] = bound;
    if (tabulist[bound]) continue;
    else if (m_is_acceptor[bound] && increased > 3)
    {
      m_rings.push_back(ring(increased + 1u, m_ringcontainer));
    }
    find_route(bound, increased, tabulist);
  }
}


std::size_t startopt::ringsearch::Search::find_torsion(double &direction, std::size_t const index_b, std::size_t const index_c) const
{
  const std::size_t M = m_coord.atoms().mains().size();
  for (std::size_t i = 0u; i<M; ++i)
  {
    std::size_t const index = m_coord.atoms().intern_of_main_idihedral(i);
    std::size_t const ib = m_coord.atoms(m_coord.atoms(index).ibond()).i_to_a();
    std::size_t const ia = m_coord.atoms(m_coord.atoms(index).iangle()).i_to_a();
    if (index_b == ib && index_c == ia)
    {
      direction = 1.0;
      return i;
    }
    else if (index_b == ia && index_c == ib)
    {
      direction = -1.0;
      return i;
    }
  }
  std::cout << "Indices " << index_b << " / " << index_c << "\n";
  throw std::runtime_error("Error: Could not match a ring bond to a main torsion of the system.");
  //return 0;
}

void startopt::ringsearch::Search::find_torsions(void)
{
  const std::size_t M = m_rings.size();
  for (std::size_t i = 0; i<M; ++i)
  {
    m_rings[i].dihedrals[0] = find_torsion(m_rings[i].dihedral_direction[0], m_rings[i].atoms[1], m_rings[i].atoms[2]);
    m_rings[i].dihedrals[1] = find_torsion(m_rings[i].dihedral_direction[1], m_rings[i].atoms[2], m_rings[i].atoms[3]);
    if (m_rings[i].size > 5)
    {
      m_rings[i].dihedrals[2] = find_torsion(m_rings[i].dihedral_direction[2], m_rings[i].atoms[3], m_rings[i].atoms[4]);
    }
    if (m_rings[i].size > 6)
    {
      m_rings[i].dihedrals[3] = find_torsion(m_rings[i].dihedral_direction[3], m_rings[i].atoms[4], m_rings[i].atoms[5]);
    }
  }
}

void startopt::ringsearch::Search::find_overlap(void)
{
  std::size_t const N(m_rings.size()), M = (N*N - N) / 2;
  m_overlap.resize(N, false);
  for (std::size_t i(0u), row(1u), col(0u); i<M; ++i)
  {
    for (auto dih_a : m_rings[row].dihedrals)
    {
      for (auto dih_b : m_rings[col].dihedrals)
      {
        std::size_t const rel_a_of_a(m_coord.atoms(dih_a).iangle()), rel_d_of_a(m_coord.atoms(dih_a).idihedral());
        std::size_t const rel_a_of_b(m_coord.atoms(dih_b).iangle()), rel_d_of_b(m_coord.atoms(dih_b).idihedral());
        if ((rel_a_of_a == rel_a_of_b && rel_d_of_a == rel_d_of_b)
            || (rel_a_of_a == rel_d_of_b && rel_d_of_a == rel_a_of_b))
        {
          m_overlap(row, col) = true;
        }
      }
    }
    ++col;
    if (col == row)
    {
      col = 0;
      ++row;
    }
  }
}

void startopt::ringsearch::Search::bias_ring(std::size_t const index, coords::float_type const force)
{
  if (index < m_rings.size())
  {
    config::biases::dihedral bias;
    bias.force = force;
    const std::size_t num_dihedrals = m_rings[index].size - 3, ring_dih_index = num_dihedrals - 2;
    for (std::size_t i = 0; i<num_dihedrals; ++i)
    {
      bias.ideal = RING_DIH[ring_dih_index][i];
      bias.a = m_rings[index].atoms[0 + i];
      bias.b = m_rings[index].atoms[1 + i];
      bias.c = m_rings[index].atoms[2 + i];
      bias.d = m_rings[index].atoms[3 + i];
      m_coord.potentials().add(bias);
    }
  }
}

void startopt::ringsearch::Search::bias_rings(std::vector<bool> const &close_it, coords::float_type const force)
{
  std::size_t const n = m_rings.size();
  m_coord.potentials().clear();
  if (close_it.size() == n)
  {
    for (std::size_t i(0u); i < n; ++i)
    {
      if (close_it[i])
      {
        bool can_be_closed(true);
        for (std::size_t j(0u); j < i; ++j)
        {
          if (close_it[j] && m_overlap(i, j) == true)
          {
            can_be_closed = false;
            break;
          }
        }
        if (can_be_closed) bias_ring(i, force);
      }
    }
  }
}

coords::float_type startopt::ringsearch::Search::energy(std::vector<bool> const &close_it, coords::float_type const force)
{
  // Set original coordinates
  m_coord.set_xyz(m_init_xyz, true);
  // Imply required bias potentials
  bias_rings(close_it, force);
  // Optimize with bias in place
  m_coord.o();
  // Remove bias
  //std::ofstream x("test.arc");
  //x << m_coord;
  m_coord.potentials().clear();
  // Relaxate
  return m_coord.o();
}

void startopt::ringsearch::Search::save_coords(std::size_t const i)
{
  m_coord.to_internal();
  m_final_ensemble.at(i) = m_coord.pes();
}

struct ring_energy
{
public:
  startopt::ringsearch::Search * p;
  ring_energy(startopt::ringsearch::Search & rso)
    : p(&rso)
  { }
  bool operator() (genetic::individuum<bool> & i, std::size_t const pop_index) const
  {
    auto E = p->energy(i.values(), Config::get().startopt.ringsearch.bias_force);
    //std::cout << scon::print_range(i.values(), " ") << ", " << pop_index << ", E = " << E << "\n";
    i.set_health(E);
    p->save_coords(pop_index);

    return p->coords().integrity();
  }
};

template<class G>
static typename G::options options_from_config()
{
  typename G::options ret_opt;
  ret_opt.crossing_chance = Config::get().optimization.global.evolution.chance_crossingover;
  ret_opt.mutation_chance_point = Config::get().optimization.global.evolution.chance_pointmutation;
  ret_opt.fitness_lower_bound = Config::get().optimization.global.selection.lin_rank_lower;
  ret_opt.fitness_upper_bound = Config::get().optimization.global.selection.lin_rank_upper;
  return ret_opt;
}

void startopt::ringsearch::Search::genetic_propagation(std::size_t const population, std::size_t const iterations)
{
  // individual holding a vector of bool, mutated by default bool mutator with float type from coords
  typedef genetic::individuum<bool, genetic::mutators::point_mutator<bool>, coords::float_type> ind_type;
  // rank based, exponential and truncated fitness
  typedef genetic::ranking::exponential<ind_type::float_type> fitness_type;
  // evolution type
  typedef genetic::evolution<ind_type, ring_energy, fitness_type> genetics_type;
  // Get options from config
  typedef genetics_type::options genetic_options_type;
  genetic_options_type opt_obj(options_from_config<genetics_type>());
  // Generate randomized initial population
  std::vector<ind_type> init_population;
  std::mt19937_64 rng = std::mt19937_64(std::random_device()());
  std::uniform_real_distribution<coords::float_type> p(0, 1);
  std::size_t const rn = m_rings.size();
  init_population.reserve(population);
  for (std::size_t i = 0u; i < population; ++i)
  {
    std::vector<bool> closure(rn, false);
    for (std::size_t j = 0u; j < rn; ++j)
    {
      if (p(rng) < Config::get().startopt.ringsearch.chance_close) closure[j] = true;
    }
    init_population.push_back(ind_type(closure));
  }
  // Resize coord saving
  m_final_ensemble.resize(population);
  // Create energy object
  ring_energy re_obj(*this);
  // Create genetics object
  genetics_type gen_obj(init_population, re_obj, opt_obj);
  // Propagte population
  gen_obj.propagate(iterations);
  
}

//
//
//startopt::ringsearch::search::search (coords::Coordinates const  & initial_coordinate) 
//  : coord_obj(initial_coordinate), N(coord_obj.size()), N_donor(0U), N_acceptor(0U), 
//  isAcceptor(N, false), isDonor(N, false), ringContainer(7U), overlap(N), initial_representation(N)
//{
//
//  // assign acceptor and donor checkvectors
//  for (std::size_t i(0U); i<N; ++i)
//  {
//    if (coord_obj.atoms(i).bonds().empty()) continue;
//    std::size_t const ia(coord_obj.atoms(i).number()), na(coord_obj.atoms(coord_obj.atoms(i).bonds(0U)).number());
//    if (ia == 1U && coord_obj.atoms(i).bonds().size() == 1 &&
//      (na == 7U || na == 8U || na == 15U || na == 16U))
//    {
//      ++N_donor;
//      isDonor[i] = true;
//    } 
//    else if (ia == 7U || ia == 8U || ia == 15U || ia == 16U
//      || ia == 9U || ia == 17U || ia == 35U || ia == 53U) 
//    {
//      ++N_acceptor;
//      isAcceptor[i] = true;
//    }
//  }
//  coord_obj.to_internal();
//  initial_representation.cartesian = coord_obj.positions;
//  initial_representation.intern = coord_obj.internals.coords;
//}
//
//std::size_t startopt::ringsearch::search::findRings (void)
//{
//  std::vector<bool> tabu(N);
//  tabu.assign(N, false);
//  for (std::size_t i=0; i<N; ++i)
//  {
//    if (isDonor[i])
//    {
//      ringContainer[0] = i;
//      routeSearch(i, 0, tabu);
//    }
//  }
//  findTorsions();
//  findOverlap();
//  return rings.size();
//}
//
//bool startopt::ringsearch::search::routeSearch (std::size_t i, std::size_t size, std::vector<bool> tabulist)
//{
//  if (size > 5) return false;
//  tabulist[i] = true;
//  const std::size_t M = coord_obj.topology[i].size(), increased = size+1;
//  std::size_t bound;
//  for (auto bound : coord_obj.atoms(i).bonds())
//  {
//    ringContainer[increased] = bound;
//    if (tabulist[bound]) continue;
//    else if (isAcceptor[bound])
//    {
//      // found 7-membered ring
//      if (increased == 6)
//      {
//        rings.push_back(ring(increased+1, ringContainer));
//        return true;
//      }
//      else if (increased == 5 || increased == 4)
//      {
//        ringContainer[increased] = bound;
//        rings.push_back(ring(increased+1, ringContainer));
//      }
//    }
//    routeSearch(bound, increased, tabulist);
//  }
//  return true;
//}
//
//void startopt::ringsearch::search::findTorsions (void)
//{
//  const std::size_t M = rings.size();
//  for (std::size_t i=0; i<M; ++i)
//  {
//    getTorsion(rings[i].dihedrals[0], rings[i].dihedral_direction[0], rings[i].atoms[1], rings[i].atoms[2]);
//    getTorsion(rings[i].dihedrals[1], rings[i].dihedral_direction[1], rings[i].atoms[2], rings[i].atoms[3]);
//    if (rings[i].size > 5) getTorsion(rings[i].dihedrals[2], rings[i].dihedral_direction[2], rings[i].atoms[3], rings[i].atoms[4]);
//    if (rings[i].size > 6) getTorsion(rings[i].dihedrals[3], rings[i].dihedral_direction[3], rings[i].atoms[4], rings[i].atoms[5]);
//  }
//  //coord_obj.write_zMatrix(cout);
//}
//
//double startopt::ringsearch::search::calcTorsion (const std::size_t A, const std::size_t B, const std::size_t C, const std::size_t D) const
//{
//  coords::Cartesian_Point temp[6] =
//  {
//    coord_obj.xyz(A) - coord_obj.xyz(B),
//    coord_obj.xyz(C) - coord_obj.xyz(B),
//    coord_obj.xyz(D) - coord_obj.xyz(C),
//    coords::Cartesian_Point(), 
//    coords::Cartesian_Point(),
//    coords::Cartesian_Point()
//  };
//  temp[3] = temp[1U].crossd(temp[0U]);
//  temp[4] = temp[1U].crossd(temp[2U]);
//  temp[5] = temp[3U].crossd(temp[4U]);
//  double torsion = temp[3].angle(temp[4]), norm = len(temp[1])*len(temp[5]);
//  return (abs(norm) > 1.0e-8 && (dot(temp[1], temp[5])/norm) < 0.0) ? torsion : -torsion;
//}
//
//void startopt::ringsearch::search::getTorsion (std::size_t &torsion_index, double &torsion_dir, const std::size_t indexB, const std::size_t indexC) const
//{
//  
//  const std::size_t M = coord_obj.atoms().mains();
//  if (M > 3)
//  {
//    std::size_t index, rel_indexA, rel_indexB;
//    std::size_t atomA, atomB;
//    for (std::size_t i=0; i<M; ++i)
//    {
//      index = coord_obj.atoms().intern_of_main_idihedral(i);
//      atomA = coord_obj.atoms(coord_obj.atoms(index).ibond()).i_to_a();
//      atomB = coord_obj.atoms(coord_obj.atoms(index).iangle()).i_to_a();
//      if (indexB == atomA && indexC == atomB)
//      {
//        torsion_index = i;
//        torsion_dir = 1.0;
//        return;
//      }
//      else if (indexB == atomB && indexC == atomA)
//      {
//        torsion_index = i;
//        torsion_dir = -1.0;
//        return;
//      }
//    }
//  }
//}
//
//void startopt::ringsearch::search::findOverlap (void)
//{
//  std::size_t const M(rings.size());
//  for (std::size_t a(0U); a<M; ++a)
//  {
//    for (std::size_t b(a+1U); b<M; ++b)
//    {
//      for (auto dih_a : rings[a].dihedrals)
//      {
//        for (auto dih_b : rings[b].dihedrals)
//        {
//          std::size_t const rel_a_of_a(coord_obj.atoms(dih_a).iangle()), rel_d_of_a(coord_obj.atoms(dih_a).idihedral()); 
//          std::size_t const rel_a_of_b(coord_obj.atoms(dih_b).iangle()), rel_d_of_b(coord_obj.atoms(dih_b).idihedral());
//          if ( (rel_a_of_a == rel_a_of_b && rel_d_of_a == rel_d_of_b) 
//            || (rel_a_of_a == rel_d_of_b && rel_d_of_a == rel_a_of_b) )
//          {
//            scon::sorted::insert_unique(overlap[a], b);
//            scon::sorted::insert_unique(overlap[b], a);
//          }
//        }
//      }
//    }
//  }
//}
//
//
//void startopt::ringsearch::search::closeRing (std::size_t index)
//{
//  if (index < rings.size())
//  {
//    //cout << "Closing ring " << index << std::endl;
//    const std::size_t num_dihedrals = rings[index].size-3, ring_dih_index = num_dihedrals-2;
//    std::size_t maintorsion_index, torsion_index;
//    double torsion_now, torsion_rotation;
//    for (std::size_t i=0; i<num_dihedrals; ++i)
//    {
//      maintorsion_index = rings[index].dihedrals[i];
//      torsion_index = coord_obj.internals.main_torsions[maintorsion_index];
//      torsion_now = calcTorsion(rings[index].atoms[0+i], rings[index].atoms[1+i], rings[index].atoms[2+i], rings[index].atoms[3+i]);
//      torsion_rotation = RING_DIH[ring_dih_index][i] - torsion_now;
//      if (torsion_rotation*torsion_now > 0.0 && rings[index].dihedral_direction[i] < 0.0) torsion_rotation *= -1.0;
//      coord_obj.internals.rotDihedral(torsion_index, torsion_rotation);
//    }
//  }
//}
//
//void ringsearch::search::opt_fixRing (std::size_t index)
//{
//  if (index < rings.size())
//  {
//    closeRing(index);
//    coord_obj.toCartesian();
//    coord_obj.fixed.assign(N, false);
//    for (std::size_t i=0; i<rings[index].size; ++i)
//    {
//      coord_obj.fixed[rings[index].atoms[i]] = true;
//    }
//    coord_obj.optimization(false);
//    coord_obj.toInternal();
//  }
//}
//
//void ringsearch::search::bias_ring (std::size_t index, const double force)
//{
//  if (index < rings.size())
//  {
//    dihBias biasPotential;
//    biasPotential.force = force;
//    biasPotential.type = bias_pot_types::QUADRATIC;
//    //cout << "Closing ring " << index << std::endl;
//    const std::size_t num_dihedrals = rings[index].size-3, ring_dih_index = num_dihedrals-2;
//    for (std::size_t i=0; i<num_dihedrals; ++i)
//    {
//      biasPotential.phi = RING_DIH[ring_dih_index][i];
//      biasPotential.a = rings[index].atoms[0+i];
//      biasPotential.b = rings[index].atoms[1+i];
//      biasPotential.c = rings[index].atoms[2+i];
//      biasPotential.d = rings[index].atoms[3+i];
//      coord_obj.potentials.dihedrals.push_back(biasPotential);
//    }
//  }
//}
