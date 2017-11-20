#include "energy_int_amoeba.h"



::tinker::parameter::parameters energy::interfaces::amoeba::amoeba_ff::tp;


energy::interfaces::amoeba::amoeba_ff::amoeba_ff(coords::Coordinates *cobj)
  : interface_base(cobj)
{
  if (!tp.valid()) tp.from_file(Config::get().get().general.paramFilename);
  //std::cout << tp << std::endl;
  std::vector<size_t> ntypes;
  for (auto atom : (*cobj).atoms()) scon::sorted::insert_unique(ntypes, atom.energy_type());
  cparams = tp.contract(ntypes);
  refined.refine(*cobj, cparams);
  parameters();
  Spackman_list_analytical1();

}



energy::interfaces::amoeba::amoeba_ff::amoeba_ff(amoeba_ff const & rhs, coords::Coordinates *cobj)
  : interface_base(cobj), part_grad(rhs.part_grad), part_energy(rhs.part_energy),
  cparams(rhs.cparams), refined(rhs.refined)
{
  interface_base::operator=(rhs);
}

energy::interfaces::amoeba::amoeba_ff::amoeba_ff(amoeba_ff && rhs, coords::Coordinates *cobj)
  : interface_base(cobj), part_grad(std::move(rhs.part_grad)), part_energy(std::move(rhs.part_energy)),
  cparams(std::move(rhs.cparams)), refined(std::move(rhs.refined))
{
  interface_base::swap(rhs);
}

void energy::interfaces::amoeba::amoeba_ff::swap(interface_base &rhs)
{
  swap(dynamic_cast<amoeba_ff&>(rhs));
}

void energy::interfaces::amoeba::amoeba_ff::swap(amoeba_ff &rhs)
{
  interface_base::swap(rhs);
  refined.swap_data(rhs.refined);
  std::swap(cparams, rhs.cparams);
  std::swap(part_energy, rhs.part_energy);
  for (size_t i(0u); i < part_grad.size(); ++i) part_grad[i].swap(rhs.part_grad[i]);
}

energy::interface_base * energy::interfaces::amoeba::amoeba_ff::clone(coords::Coordinates * coord_object) const
{
  amoeba_ff * tmp = new amoeba_ff(*this, coord_object);
  return tmp;
}

energy::interface_base * energy::interfaces::amoeba::amoeba_ff::move(coords::Coordinates * coord_object)
{
  amoeba_ff * tmp = new amoeba_ff(std::move(*this), coord_object);
  return tmp;
}

// initialize using coordinates pointer

// update structure (account for topology or rep change)
void energy::interfaces::amoeba::amoeba_ff::update(bool const skip_topology)
{
  if (!skip_topology)
  {
    std::vector<size_t> ntypes;
    for (auto atom : (*coords).atoms()) scon::sorted::insert_unique(ntypes, atom.energy_type());
    cparams = tp.contract(ntypes);
    refined.refine((*coords), cparams);
  }
  else refined.refine_nb((*coords));
}

// Output functions
void energy::interfaces::amoeba::amoeba_ff::print_E(std::ostream &) const
{

}

void energy::interfaces::amoeba::amoeba_ff::print_E_head(std::ostream &S, bool const endline) const
{
  S << "Potentials\n";
  S << std::right << std::setw(24) << "B";
  S << std::right << std::setw(24) << "A";
  S << std::right << std::setw(24) << "U";
  S << std::right << std::setw(24) << "T\n";
  S << std::right << std::setw(24) << "V";
  S << std::right << std::setw(24) << "MUL";
  S << std::right << std::setw(24) << "POL";
  S << std::right << std::setw(24) << "SR";
  S << std::right << std::setw(24) << "OOP";
  S << std::right << std::setw(24) << "-";
  S << std::right << std::setw(24) << "SUM\n";
  size_t const SSS(coords->subsystems().size());
  if (SSS > 1U)
  {
    S << std::right << std::setw(24) << "SubSystem IA:";
    for (size_t i(0U); i < SSS; ++i)
    {
      for (size_t j(0U); j <= i; ++j)
      {
        std::stringstream ss;
        ss << "[" << j << "<>" << i << "]";
        S << std::right << std::setw(24) << ss.str();
      }
    }
  }
  size_t const IAS(coords->interactions().size());
  if (IAS > 0U)
  {
    S << std::right << std::setw(24) << "SubSystem IA:";
    for (size_t i(0U); i < IAS; ++i)
    {
      S << std::right << std::setw(24) << i;
    }
  }
  S << '\n';
  S << "Count\n";
  S << std::right << std::setw(24) << refined.bonds().size();
  S << std::right << std::setw(24) << refined.angles().size();
  S << std::right << std::setw(24) << refined.ureys().size();
  S << std::right << std::setw(24) << refined.impropers().size();
  S << std::right << std::setw(24) << refined.imptors().size() << '\n';
  S << std::right << std::setw(24) << refined.torsions().size();
  size_t ia(refined.ia_count());
  S << std::right << std::setw(24) << ia;
  S << std::right << std::setw(24) << (!cparams.charges().empty() ? ia : 0u);
  S << std::right << std::setw(24) << "-";
  S << std::right << std::setw(24) << "-";
  S << std::right << std::setw(24) << "-";
  if (endline) S << '\n';
}

void energy::interfaces::amoeba::amoeba_ff::print_E_short(std::ostream &S, bool const endline) const
{

  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[BOND];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[ANGLE];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[UREY];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[TORSION] << '\n';
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[VDW];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[MULTIPOLE];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[POLARIZE];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[SHORTRANGE];
  S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[OPBEND];
  S << std::right << std::setw(24) << "-";
  S << std::right << std::setw(24) << "-";
  S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << energy;
  size_t const IAS(coords->interactions().size());
  if (IAS > 0U)
  {
    S << '\n' << std::right << std::setw(24) << "SubSystem IA:";
    for (size_t i(0U); i < IAS; ++i)
    {
      S << std::right << std::setw(24) << coords->interactions(i).energy;
    }
  }
  if (endline) S << '\n';
}

void energy::interfaces::amoeba::amoeba_ff::print_G_tinkerlike(std::ostream &S, bool const aggregate) const
{
  size_t const N(coords->size());
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
    for (size_t i(0u); i < N; ++i)
    {
      S << i + 1 << "-x\n";
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[BOND][i].x();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[ANGLE][i].x();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[UREY][i].x();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[TORSION][i].x() << '\n';
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[VDW][i].x();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[MULTIPOLE][i].x();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[POLARIZE][i].x();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[OPBEND][i].x();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << coords->g_xyz(i).x() << '\n';
      S << i + 1 << "-y\n";
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[BOND][i].y();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[ANGLE][i].y();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[UREY][i].y();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[TORSION][i].y() << '\n';
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[VDW][i].y();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[MULTIPOLE][i].y();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[POLARIZE][i].y();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[OPBEND][i].y();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << coords->g_xyz(i).y() << '\n';
      S << i + 1 << "-z\n";
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[BOND][i].z();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[ANGLE][i].z();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[UREY][i].z();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[TORSION][i].z() << '\n';
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[VDW][i].z();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[MULTIPOLE][i].z();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[POLARIZE][i].z();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << part_grad[OPBEND][i].z();
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
      S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << coords->g_xyz(i).z() << '\n';
    }
  }
  else
  {
    double p_b(0.0), p_a(0.0), p_u(0.0)/*, p_im(0.0), p_it(0.0)*/, p_t(0.0), p_vdw(0.0), p_mult(0.0), p_pol(0.0), p_oop(0.0), p_sum(0.0);
    for (size_t i(0u); i < N; ++i)
    {
      p_b += fabs(part_grad[BOND][i].x()) + fabs(part_grad[BOND][i].y()) + fabs(part_grad[BOND][i].z());
      p_a += fabs(part_grad[ANGLE][i].x()) + fabs(part_grad[ANGLE][i].y()) + fabs(part_grad[ANGLE][i].z());
      p_u += fabs(part_grad[UREY][i].x()) + fabs(part_grad[UREY][i].y()) + fabs(part_grad[UREY][i].z());
      p_t += fabs(part_grad[TORSION][i].x()) + fabs(part_grad[TORSION][i].y()) + fabs(part_grad[TORSION][i].z());
      p_vdw += fabs(part_grad[VDW][i].x()) + fabs(part_grad[VDW][i].y()) + fabs(part_grad[VDW][i].z());
      p_mult += fabs(part_grad[MULTIPOLE][i].x()) + fabs(part_grad[MULTIPOLE][i].y()) + fabs(part_grad[MULTIPOLE][i].z());
      p_pol += fabs(part_grad[POLARIZE][i].x()) + fabs(part_grad[POLARIZE][i].y()) + fabs(part_grad[POLARIZE][i].z());
      p_oop += fabs(part_grad[OPBEND][i].x()) + fabs(part_grad[OPBEND][i].y()) + fabs(part_grad[OPBEND][i].z());
      p_sum += fabs(coords->g_xyz(i).x()) + fabs(coords->g_xyz(i).y()) + fabs(coords->g_xyz(i).z());
    }
    S << "Agg\n";
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_b;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_a;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_u;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_t << '\n';
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_vdw;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_mult;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_pol;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_oop;
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << "-";
    S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << p_sum << '\n';
  }
}

void energy::interfaces::amoeba::amoeba_ff::pre(void)
{
  for (auto & e : part_energy) e = 0.0;
  for (auto & g : part_grad) g.assign(coords->size(), coords::Cartesian_Point());
  std::array<double, 3> za;
  za[0] = 0.0;
  za[1] = 0.0;
  za[2] = 0.0;
  std::array<std::array<double, 3>, 3> zv;
  zv[0] = za;
  zv[1] = za;
  zv[2] = za;
  for (auto & v : part_virial) v = zv;
  integrity = true;



}

void energy::interfaces::amoeba::amoeba_ff::post(void)
{
  energy = 0.0;
  for (auto const & e : part_energy) energy += e;
  coords->clear_g_xyz();
  for (auto const & g : part_grad) coords->sum_g_xyz(g);
  std::array<double, 3> za;
  za[0] = 0.0;
  za[1] = 0.0;
  za[2] = 0.0;
  std::array<std::array<double, 3>, 3> zv;
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
}

void energy::interfaces::amoeba::amoeba_ff::to_stream(std::ostream &S) const
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


