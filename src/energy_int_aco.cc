#include <sstream>
#include <cstddef>
#include "energy_int_aco.h"
#include "configuration.h"
#include "scon_utility.h"

::tinker::parameter::parameters energy::interfaces::aco::aco_ff::tp;

/*! Constructs a force-field energy interface
 *
 * Constructor for a force field energy interface.
 * Atom types are gathered and subsequently contracted
 * @todo: Describe this better once you understood it.
 *
 * @param cobj: Pointer to coordinates object for which energy interface will perform
 */
energy::interfaces::aco::aco_ff::aco_ff (coords::Coordinates *cobj) 
   : interface_base(cobj) 
{
  // tp are static tinker parameters envoked above 
  // (::tinker::parameter::parameters energy::interfaces::aco::aco_ff::tp;)

  interactions = true;
  if (!tp.valid())
  {
    // Here, we read the param file containing
    // the force field parameters
    // (at this point, we read ALL the ff-parameters,
    // even the ones we might not need
    tp.from_file(Config::get().get().general.paramFilename);
  }
  std::vector<std::size_t> types;
  for (auto atom : (*cobj).atoms())
  {
    scon::sorted::insert_unique(types, atom.energy_type());
  }
  cparams = tp.contract(types);
  refined.refine(*cobj, cparams);
}

energy::interfaces::aco::aco_ff::aco_ff (aco_ff const & rhs, 
  coords::Coordinates *cobj) : interface_base(cobj), 
  part_grad(rhs.part_grad), part_energy(rhs.part_energy), 
  cparams(rhs.cparams), refined(rhs.refined)
{
  interface_base::operator=(rhs);
}

energy::interfaces::aco::aco_ff::aco_ff (aco_ff && rhs, 
  coords::Coordinates *cobj) : interface_base(cobj), 
  part_grad(std::move(rhs.part_grad)), part_energy(std::move(rhs.part_energy)), 
  cparams(std::move(rhs.cparams)), refined(std::move(rhs.refined))
{
  interface_base::swap(rhs);
  
  /** 
   * This is necessary because the compiler-provided move-constructor for
   * tinker::refined is sometimes malfunctioning and containing a faulty pointer.
   * This is a (somewhat dirty) quickfix since I was to lazy to write
   * a proper, custom move constructor. If you read this message and you have 
   * nothing to do, you should probably do this :)
   */
  this->refined.setCoordsPointer(cobj);
}

void energy::interfaces::aco::aco_ff::swap (interface_base &rhs)
{
  swap(dynamic_cast<aco_ff&>(rhs));
}

void energy::interfaces::aco::aco_ff::swap (aco_ff &rhs)
{
  interface_base::swap(rhs);
  refined.swap_data(rhs.refined);
  cparams.swap(rhs.cparams);
  std::swap(part_energy, rhs.part_energy);
  for (std::size_t i(0u); i<part_grad.size(); ++i) part_grad[i].swap(rhs.part_grad[i]);
}

energy::interface_base * energy::interfaces::aco::aco_ff::clone (coords::Coordinates * coord_object) const
{
  aco_ff * tmp = new aco_ff(*this, coord_object);
  return tmp;
}

energy::interface_base * energy::interfaces::aco::aco_ff::move (coords::Coordinates * coord_object)
{
  aco_ff * tmp = new aco_ff(std::move(*this), coord_object);
  return tmp;
}

// initialize using coordinates pointer

// update structure (account for topology or rep change)
void energy::interfaces::aco::aco_ff::update (bool const skip_topology)
{
  if (!skip_topology) 
  {
    std::vector<std::size_t> types;
    for (auto atom : (*coords).atoms()) scon::sorted::insert_unique(types, atom.energy_type());
    cparams = tp.contract(types);
    refined.refine((*coords), cparams);
  }
  else 
  {
    refined.refine_nb((*coords));
  }
}

// Output functions
void energy::interfaces::aco::aco_ff::print_E (std::ostream &) const
{

}

void energy::interfaces::aco::aco_ff::print_E_head (std::ostream &S, bool const endline) const
{
  S << "Potentials\n";
  S << std::right << std::setw(24) << "B";
  S << std::right << std::setw(24) << "A";
  S << std::right << std::setw(24) << "U";
  S << std::right << std::setw(24) << "ID";
  S << std::right << std::setw(24) << "IT\n";
  S << std::right << std::setw(24) << "T";
  S << std::right << std::setw(24) << "V";
  S << std::right << std::setw(24) << "C";
  S << std::right << std::setw(24) << "SOLV";
  S << std::right << std::setw(24) << "-";
  S << std::right << std::setw(24) << "SUM\n";
  std::size_t const SSS(coords->subsystems().size());
  if (SSS > 1U)
  {
    S << std::right << std::setw(24) << "SubSystem IA:";
    for (std::size_t i(0U); i<SSS; ++i)
    {
      for (std::size_t j(0U); j<=i; ++j)
      {
        std::stringstream ss;
        ss << "[" << j << "<>" << i << "]";
        S << std::right << std::setw(24) << ss.str();
      }
    }
  }
  //std::size_t const IAS(coords->interactions().size());
  //if (IAS > 0U) 
  //{
  //  for (std::size_t i(0U); i<IAS; ++i)
  //  {
  //    std::stringstream ss;
  //    ss << "IA" << i;
  //    S << std::right << std::setw(24) << ss.str();
  //  }
  //}
  S << '\n';
  S << "Count\n";
  S << std::right << std::setw(24) << refined.bonds().size();
  S << std::right << std::setw(24) << refined.angles().size();
  S << std::right << std::setw(24) << refined.ureys().size();
  S << std::right << std::setw(24) << refined.impropers().size();
  S << std::right << std::setw(24) << refined.imptors().size() << '\n';
  S << std::right << std::setw(24) << refined.torsions().size();
  std::size_t ia(refined.ia_count());
  S << std::right << std::setw(24) << ia;
  if (Config::get().general.input == config::input_types::AMBER)
  {
    S << std::right << std::setw(24) << "see vdw";
  }
  else
  {
    S << std::right << std::setw(24) << (!cparams.charges().empty() ? ia : 0u);
  }
  S << std::right << std::setw(24) << "-";
  S << std::right << std::setw(24) << "-";
  S << std::right << std::setw(24) << "-";
  if (endline) S << '\n';
}

void energy::interfaces::aco::aco_ff::print_E_short (std::ostream &S, bool const endline) const
{

  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::BOND];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::ANGLE];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::UREY];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::IMPROPER];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::IMPTORSION] << '\n';
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::TORSION];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::VDW];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::CHARGE];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::SOLVATE];
  S << std::right << std::setw(24) << "-";
  S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << energy;
  std::size_t const IAS(coords->interactions().size());
  if (IAS > 0U) 
  {
    S << '\n' << std::right << std::setw(24) << "SubSystem IA:";
    for (std::size_t i(0U); i<IAS; ++i)
    {
      S << std::right << std::setw(24) << coords->interactions(i).energy;
    }
  }
  if (endline) S << '\n';
}

void energy::interfaces::aco::aco_ff::print_G_tinkerlike (std::ostream &S, bool const aggregate) const
{
  std::size_t const N(coords->size());
  if (!aggregate)
  {
    S << "Gradients\n";
    S << std::right << std::setw(24) << "B";
    S << std::right << std::setw(24) << "A";
    S << std::right << std::setw(24) << "U";
    S << std::right << std::setw(24) << "ID";
    S << std::right << std::setw(24) << "IT\n";
    S << std::right << std::setw(24) << "T";
    S << std::right << std::setw(24) << "V+C";
    S << std::right << std::setw(24) << "-";
    S << std::right << std::setw(24) << "-";
    S << std::right << std::setw(24) << "SUM\n";
    for (std::size_t i(0u); i<N; ++i)
    {
        S << i+1 << "-x\n";
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::BOND][i].x();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::ANGLE][i].x();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::UREY][i].x();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::IMPROPER][i].x();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::IMPTORSION][i].x() << '\n';
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::TORSION][i].x();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::VDWC][i].x();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << coords->g_xyz(i).x() << '\n';
        S << i+1 << "-y\n";
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::BOND][i].y();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::ANGLE][i].y();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::UREY][i].y();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::IMPROPER][i].y();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::IMPTORSION][i].y() << '\n';
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::TORSION][i].y();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::VDWC][i].y();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << coords->g_xyz(i).y() << '\n';
        S << i+1 << "-z\n";
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::BOND][i].z();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::ANGLE][i].z();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::UREY][i].z();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::IMPROPER][i].z();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::IMPTORSION][i].z() << '\n';
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::TORSION][i].z();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[types::VDWC][i].z();
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
        S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << coords->g_xyz(i).z() << '\n';
    }
  }
  else
  {
    coords::float_type p_b(0.0), p_a(0.0), p_u(0.0), p_im(0.0), p_it(0.0), p_t(0.0), p_vdwc(0.0), p_sum(0.0);
    for (std::size_t i(0u); i<N; ++i)
    {
      p_b += fabs(part_grad[types::BOND][i].x()) + fabs(part_grad[types::BOND][i].y()) + fabs(part_grad[types::BOND][i].z());
      p_a += fabs(part_grad[types::ANGLE][i].x()) + fabs(part_grad[types::ANGLE][i].y()) + fabs(part_grad[types::ANGLE][i].z());
      p_u += fabs(part_grad[types::UREY][i].x()) + fabs(part_grad[types::UREY][i].y()) + fabs(part_grad[types::UREY][i].z());
      p_im += fabs(part_grad[types::IMPROPER][i].x()) + fabs(part_grad[types::IMPROPER][i].y()) + fabs(part_grad[types::IMPROPER][i].z());
      p_it += fabs(part_grad[types::IMPTORSION][i].x()) + fabs(part_grad[types::IMPTORSION][i].y()) + fabs(part_grad[types::IMPTORSION][i].z());
      p_t += fabs(part_grad[types::TORSION][i].x()) + fabs(part_grad[types::TORSION][i].y()) + fabs(part_grad[types::TORSION][i].z());
      p_vdwc += fabs(part_grad[types::VDWC][i].x()) + fabs(part_grad[types::VDWC][i].y()) + fabs(part_grad[types::VDWC][i].z());
      p_sum += fabs(coords->g_xyz(i).x()) + fabs(coords->g_xyz(i).y()) + fabs(coords->g_xyz(i).z());
    }
    S << "Agg\n";
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_b;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_a;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_u;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_im;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_it << '\n';
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_t;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_vdwc;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_sum << '\n';
  }
}

void energy::interfaces::aco::aco_ff::pre    (void)
{
  for (auto & e : part_energy) e = 0.0;
  for (auto & g : part_grad) g.assign(coords->size(), coords::Cartesian_Point(0.0, 0.0, 0.0)); 
  std::array<coords::float_type, 3> za;
  za[0] = 0.0;
  za[1] = 0.0;
  za[2] = 0.0;
  std::array<std::array<coords::float_type, 3>, 3> zv;
  zv[0] = za;
  zv[1] = za;
  zv[2] = za;
  for (auto & v : part_virial) v = zv;
  integrity = true;
}

void energy::interfaces::aco::aco_ff::post   (void)
{
  energy = 0.0;
  for (auto const & e : part_energy) energy += e;
  coords->clear_g_xyz();
  for (auto const & g : part_grad) coords->sum_g_xyz(g);
  std::array<coords::float_type, 3> za;
  za[0] = 0.0;
  za[1] = 0.0;
  za[2] = 0.0;
  std::array<std::array<coords::float_type, 3>, 3> zv;
  zv[0] = za;
  zv[1] = za;
  zv[2] = za;
  for (auto const & v : part_virial)
  {
    zv[0][0] += v[0][0];
    zv[0][1] += v[0][1];
    zv[0][2] += v[0][2];
    zv[1][0] += v[1][0];
    zv[1][1] += v[1][1];
    zv[1][2] += v[1][2];
    zv[2][0] += v[2][0];
    zv[2][1] += v[2][1];
    zv[2][2] += v[2][2];
  }
  coords->set_virial(zv);
  //std::cout << coords->virial()[0][0] << "   " << coords->virial()[1][1] << "   " << coords->virial()[2][2] << std::endl;
}

void energy::interfaces::aco::aco_ff::to_stream (std::ostream &S) const 
{ 
  interface_base::to_stream(S);
  for (auto rep : part_grad) 
  {
    S << "Part grad:" << std::endl;
    for (auto grad : rep) { S << grad << std::endl; }
  }
  for (auto const & rep : part_virial)
  {
    S << "Part grad:" << std::endl;
    for (auto const & v : rep) 
    { 
      S << v[0] << "," << v[0] << "," << v[0] << std::endl; 
    }
  }
  for (auto en : part_energy) S << "Part energy: " << en << std::endl;
  S << "Params:" << std::endl;
  S << cparams;
  S << "Refined:" << std::endl;
  S << refined;
}