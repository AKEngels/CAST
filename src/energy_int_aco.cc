#include <sstream>
#include <cstddef>
#include "energy_int_aco.h"
#include "configuration.h"
#include "Scon/scon_utility.h"

::tinker::parameter::parameters energy::interfaces::aco::aco_ff::tp;
::tinker::parameter::parameters energy::interfaces::aco::aco_ff::cparams;

/*! Constructs a force-field energy interface
 *
 * Constructor for a force field energy interface.
 * Atom types are gathered and subsequently contracted
 * @todo: Describe this better once you understood it.
 *
 * @param cobj: Pointer to coordinates object for which energy interface will perform
 */
energy::interfaces::aco::aco_ff::aco_ff(coords::Coordinates* cobj)
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

  refined = ::tinker::refine::refined(*cobj, cparams);

  restrainInternals(*cobj, refined);

	double const min_cut = std::min({ Config::get().periodics.pb_box.x(), Config::get().periodics.pb_box.y(), Config::get().periodics.pb_box.z() }) / 2.0;
	if (Config::get().periodics.periodic && Config::get().energy.cutoff > min_cut)
	{
		std::cout << "\n!!! WARNING! Forcefield cutoff too big! Your cutoff should be smaller than " << min_cut << "! !!!\n\n";
	}
}


// NEW FOR FIXED INTERNALS: SET IDEAL VALUES APPROPRIATELY
void energy::interfaces::aco::restrainInternals(coords::Coordinates const& coords_in, ::tinker::refine::refined & refined, const float torsionalForce)
{
  coords::Coordinates const* coords = &coords_in;
  using scon::len;
  for (auto & bond : refined.set_bonds())
  {
    auto const bv(coords->xyz(bond.atoms[0]) - coords->xyz(bond.atoms[1])); // r_ij (i=1, j=2)
    auto const d = len(bv);
    bond.ideal = d;
    //bond.force = 10;
  }
  for (auto & angle : refined.set_angles())
  {
    auto const
      av1(coords->xyz(angle.atoms[0]) - coords->xyz(angle.atoms[1])),
      av2(coords->xyz(angle.atoms[2]) - coords->xyz(angle.atoms[1]));
    angle.ideal = scon::angle(av1, av2).degrees();
    //angle.force = 10.0;
  }
  for (auto & urey : refined.set_ureys())
  {
    coords::Cartesian_Point const bv =
      coords->xyz(urey.atoms[0]) - coords->xyz(urey.atoms[1]);
    coords::float_type const d = len(bv);
    urey.ideal = d;
  }
  for (auto & torsion : refined.set_torsions())
  {
    // Get bonding vectors
    coords::Cartesian_Point const b01 =
      coords->xyz(torsion.atoms[1]) - coords->xyz(torsion.atoms[0]);
    coords::Cartesian_Point const b12 =
      coords->xyz(torsion.atoms[2]) - coords->xyz(torsion.atoms[1]);
    coords::Cartesian_Point const b23 =
      coords->xyz(torsion.atoms[3]) - coords->xyz(torsion.atoms[2]);
    // Cross terms
    coords::Cartesian_Point const t = cross(b01, b12);
    coords::Cartesian_Point const u = cross(b12, b23);
    // Get length and variations
    coords::float_type const tl2 = dot(t, t);
    coords::float_type const ul2 = dot(u, u);
    // ...
    coords::float_type const tlul = sqrt(tl2*ul2);
    coords::float_type const r12 = len(b12);
    // scalar and length variations
    coords::float_type const cos_scalar0 = dot(t, u);
    coords::float_type const cos_scalar1 = tlul;
    // Get multiple sine and cosine values


    // cross of cross
    coords::Cartesian_Point const tu = cross(t, u);
    coords::float_type const sin_scalar0 = dot(b12, tu);
    coords::float_type const sin_scalar1 = r12 * tlul;
    coords::float_type const cos = cos_scalar0 / cos_scalar1;
    coords::float_type const sin = sin_scalar0 / sin_scalar1;


    torsion.p.number = 1;

    for (std::size_t j(0U); j < torsion.p.number; ++j)
    {
      torsion.p.order[j] = 1u;
      torsion.p.max_order = 1u;
      if (sin > 0.0)
        torsion.p.ideal[j] = std::acos(cos) * SCON_180PI;
      else
        torsion.p.ideal[j] = (SCON_2PI - std::acos(cos)) * SCON_180PI;
      //torsion.p.force[j] = torsionalForce*10;
    }
  }
  for (auto & imptor : refined.set_imptors())   //for every improper torsion
  {
    //energy calculation
    coords::Cartesian_Point const ba =
      coords->xyz(imptor.ligand[1]) - coords->xyz(imptor.ligand[0]);
    coords::Cartesian_Point const cb =
      coords->xyz(imptor.center) - coords->xyz(imptor.ligand[1]);
    coords::Cartesian_Point const dc =
      coords->xyz(imptor.twist) - coords->xyz(imptor.center);

    coords::Cartesian_Point const t = cross(ba, cb);
    coords::Cartesian_Point const u = cross(cb, dc);
    coords::float_type const tl2 = dot(t, t);
    coords::float_type const ul2 = dot(u, u);
    coords::float_type const tlul = sqrt(tl2*ul2);
    coords::float_type const r12 = len(cb);
    coords::Cartesian_Point const tu = cross(t, u);

    coords::float_type const cosine = dot(t, u) / tlul;
    coords::float_type const sine = dot(cb, tu) / (r12*tlul);

    if (sine > 0.0)
      imptor.p.ideal[0] = std::acos(cosine);
    else
      imptor.p.ideal[0] = SCON_2PI - std::acos(cosine);
    //imptor.p.force[0] = 10.0;
  }
  for (auto & improper : refined.set_impropers())
  {
    coords::Cartesian_Point const ba =
      coords->xyz(improper.ligand[0]) - coords->xyz(improper.center);
    coords::Cartesian_Point const cb =
      coords->xyz(improper.ligand[1]) - coords->xyz(improper.ligand[0]);
    coords::Cartesian_Point const dc =
      coords->xyz(improper.twist) - coords->xyz(improper.ligand[1]);
    coords::Cartesian_Point const t(cross(ba, cb));
    coords::Cartesian_Point const u(cross(cb, dc));
    coords::Cartesian_Point const tu(cross(t, u));
    coords::float_type const rt2(dot(t, t)), ru2(dot(u, u)), rtru = sqrt(rt2*ru2);
    coords::float_type const rcb(len(cb));
    coords::float_type const cosine(std::min(1.0, std::max(-1.0, (dot(t, u) / rtru))));
    coords::float_type const sine(dot(cb, tu) / (rcb*rtru));
    coords::float_type const angle(sine < 0.0 ?
      -acos(cosine) : acos(cosine));
    improper.p.ideal[0U] = angle;
    //improper.p.force[0] = 5.0;
  }
}


energy::interfaces::aco::aco_ff::aco_ff (aco_ff const & rhs, 
  coords::Coordinates *cobj) : interface_base(cobj), refined(rhs.refined),
  part_energy(rhs.part_energy), part_grad(rhs.part_grad) 
{
	interface_base::operator=(rhs);
}

energy::interfaces::aco::aco_ff::aco_ff(aco_ff&& rhs,
	coords::Coordinates* cobj) : interface_base(cobj), refined(std::move(rhs.refined)),
	part_energy(std::move(rhs.part_energy)), part_grad(std::move(rhs.part_grad))
{
  interface_base::swap(rhs);
}

void energy::interfaces::aco::aco_ff::swap(interface_base& rhs)
{
	swap(dynamic_cast<aco_ff&>(rhs));
}

void energy::interfaces::aco::aco_ff::swap(aco_ff& rhs)
{
	interface_base::swap(rhs);
	refined.swap_data(rhs.refined);
	std::swap(cparams, rhs.cparams);
	std::swap(part_energy, rhs.part_energy);
	for (std::size_t i(0u); i < part_grad.size(); ++i)
	{
		part_grad[i].swap(rhs.part_grad[i]);
	}
}

energy::interface_base* energy::interfaces::aco::aco_ff::clone(
	coords::Coordinates* coord_object) const
{
	aco_ff* tmp = new aco_ff(*this, coord_object);
	return tmp;
}

energy::interface_base* energy::interfaces::aco::aco_ff::move(
	coords::Coordinates* coord_object)
{
	aco_ff* tmp = new aco_ff(std::move(*this), coord_object);
	return tmp;
}

// initialize using coordinates pointer

// update structure (account for topology or rep change)
//
// In this function there is a pointer to coords. 
// As coords has a pointer to energy, this is recursive and should be removed in the future
void energy::interfaces::aco::aco_ff::update(bool const skip_topology)
{
  std::vector<std::size_t> types;
  for (auto && atom : (*coords).atoms())
  {
    scon::sorted::insert_unique(types, atom.energy_type());
  }
  
  if (!skip_topology) 
  {
    cparams = tp.contract(types);
    refined = ::tinker::refine::refined((*coords), cparams);
    restrainInternals(*coords, refined);
    purge_nb_at_same_molecule(*coords, refined);
  }
  else 
  {
    refined.refine_nb(*coords);
  }
}

std::vector<coords::float_type> energy::interfaces::aco::aco_ff::charges() const
{
	std::vector<coords::float_type> c;

	if (Config::get().general.input == config::input_types::AMBER || Config::get().general.chargefile)
	{
		c = Config::get().coords.amber_charges;  // get amber charges
		for (double& q : c) q = q / 18.2223;     // convert all charges to elementary units
	}

	else  // if no amber charges: get charges from charge parameters
	{
		if (tp.charges().empty())
		{
			throw std::runtime_error("No charges in parameters.");
		}
		for (auto&& atom : coords->atoms())
		{
			for (auto&& chg : tp.charges())
			{
				auto t_of_atom = tp.type(atom.energy_type(), tinker::CHARGE);
				if (chg.index == t_of_atom)
				{
					c.push_back(chg.c);
					break;
				}
			}
		}
	}

	if (c.size() != coords->size())  // check if correct number of parameters is found
	{
		std::cout << "found " << c.size() << " charges but system has " << coords->size() << " atoms\n";
		throw std::runtime_error("Didn't find all charges.");
	}
	return c;
}

// Output functions
void energy::interfaces::aco::aco_ff::print_E(std::ostream&) const
{

}

void energy::interfaces::aco::aco_ff::print_E_head(std::ostream& S, bool const endline) const
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
		for (std::size_t i(0U); i < SSS; ++i)
		{
			for (std::size_t j(0U); j <= i; ++j)
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

void energy::interfaces::aco::aco_ff::print_E_short(std::ostream& S, bool const endline) const
{
	S << "Energies\n";
	S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::BOND];
	S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::ANGLE];
	S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::UREY];
	S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::IMPROPER];
	S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::IMPTORSION] << '\n';
	S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::TORSION];
	S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::VDW];
	S << std::right << std::setw(24) << std::fixed << std::setprecision(8) << part_energy[types::CHARGE];
	S << std::right << std::setw(24) << "-";
	S << std::right << std::setw(24) << "-";
	S << std::right << std::setw(24) << std::fixed << std::setprecision(12) << energy;
	std::size_t const IAS(coords->interactions().size());
	if (IAS > 0U)
	{
		S << '\n' << std::right << std::setw(24) << "SubSystem IA:";
		for (std::size_t i(0U); i < IAS; ++i)
		{
			S << std::right << std::setw(24) << coords->interactions(i).energy;
		}
	}
	if (endline) S << '\n';
}

void energy::interfaces::aco::aco_ff::print_G_tinkerlike(std::ostream& S, bool const aggregate) const
{
	S << " Cartesian Gradient Breakdown over Individual Atoms :" << std::endl << std::endl;
	S << "  Type      Atom              dE/dX       dE/dY       dE/dZ          Norm" << std::endl << std::endl;
	for (std::size_t k = 0; k < coords->size(); ++k)
	{
		S << " Anlyt";
		S << std::right << std::setw(10) << k + 1U;
		S << "       ";
		S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).x();
		S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).y();
		S << std::right << std::fixed << std::setw(12) << std::setprecision(4) << coords->g_xyz(k).z();
		S << std::right << std::fixed << std::setw(12) << std::setprecision(4);
		S << std::sqrt(
			coords->g_xyz(k).x() * coords->g_xyz(k).x()
			+ coords->g_xyz(k).y() * coords->g_xyz(k).y()
			+ coords->g_xyz(k).z() * coords->g_xyz(k).z()) << "\n";
	}

	S << "\nSpecial stuff for ACO interface:\n";
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
		for (std::size_t i(0u); i < N; ++i)
		{
			S << i + 1 << "-x\n";
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
			S << i + 1 << "-y\n";
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
			S << i + 1 << "-z\n";
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
		for (std::size_t i(0u); i < N; ++i)
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

void energy::interfaces::aco::aco_ff::pre(void)
{
	for (auto& e : part_energy) e = 0.0;
	for (auto& g : part_grad) g.assign(coords->size(), coords::Cartesian_Point(0.0, 0.0, 0.0));
	std::array<coords::float_type, 3> za;
	za[0] = 0.0;
	za[1] = 0.0;
	za[2] = 0.0;
	std::array<std::array<coords::float_type, 3>, 3> zv;
	zv[0] = za;
	zv[1] = za;
	zv[2] = za;
	for (auto& v : part_virial) v = zv;
	integrity = true;
}

void energy::interfaces::aco::aco_ff::post(void)
{
	energy = 0.0;
	for (auto const& e : part_energy) energy += e;
	coords->clear_g_xyz();
	for (auto const& g : part_grad) coords->sum_g_xyz(g);
	std::array<coords::float_type, 3> za;
	za[0] = 0.0;
	za[1] = 0.0;
	za[2] = 0.0;
	std::array<std::array<coords::float_type, 3>, 3> zv;
	zv[0] = za;
	zv[1] = za;
	zv[2] = za;
	for (auto const& v : part_virial)
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

void energy::interfaces::aco::aco_ff::to_stream(std::ostream& S) const
{
	interface_base::to_stream(S);
	for (auto rep : part_grad)
	{
		S << "Part grad:" << std::endl;
		for (auto grad : rep) { S << grad << std::endl; }
	}
	for (auto const& rep : part_virial)
	{
		S << "Part grad:" << std::endl;
		for (auto const& v : rep)
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