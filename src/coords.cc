#include <cmath>
#include <stdexcept>
//#include <fftw3.h>
#include "atomic.h"
#include "coords.h"
#include "configuration.h"
#include "coords_io.h"
#include "lbfgs.h"
#include "Array.h"
#define SUPERPI 3.141592653589793238

#include "optimization_dimer.h"
std::ostream& coords::operator<< (std::ostream &stream, internal_relations const & inter)
{
  stream << "[Rel: B,A,D: " << inter.rel_b << ", " << inter.rel_a << ", " << inter.rel_d;
  stream << "]:[Atom: " << inter.atom_of_inter_index << ", Inter: " << inter.inter_of_atom_index << "]:[Follow: ";
  for (std::size_t i = 0; i < inter.follow_ups.size(); ++i)
  {
    stream << (i != 0 ? ", " : "") << inter.follow_ups[i];
  }
  stream << "]:[Dep: ";
  for (std::size_t i = 0; i < inter.dependant_torsions.size(); ++i)
  {
    stream << (i != 0 ? ", " : "") << inter.dependant_torsions[i];
  }
  stream << "] Main: " << inter.is_main_torsion << "\n";
  return stream;
}

std::ostream& coords::operator<< (std::ostream &stream, Atom const & atom)
{
  stream << atom.m_symbol << " (" << atom.m_number << ", " << atom.m_mass << "):[";
  for (std::size_t i = 0; i < atom.m_bonds.size(); ++i)
  {
    stream << (i != 0 ? ", " : "") << atom.m_bonds[i];
  }
  stream << "] Sys: " << atom.m_system << ", Type: " << atom.m_etype << ", Sub: " << atom.m_sub_id << '\n';
  stream << "[B,A,D: " << atom.rel_b << ", " << atom.rel_a << ", " << atom.rel_d;
  stream << "]:[Atom: " << atom.atom_of_inter_index << ", Inter: " << atom.inter_of_atom_index << "]:[Follow: ";
  for (std::size_t i = 0; i < atom.follow_ups.size(); ++i)
  {
    stream << (i != 0 ? ", " : "") << atom.follow_ups[i];
  }
  stream << "]:[Dep: ";
  for (std::size_t i = 0; i < atom.dependant_torsions.size(); ++i)
  {
    stream << (i != 0 ? ", " : "") << atom.dependant_torsions[i];
  }
  stream << "] Main: " << atom.is_main_torsion << '\n';
  return stream;
}

std::ostream& coords::operator<< (std::ostream &stream, Atoms const & atoms)
{
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    stream << atoms.atom(i);
  }
  return stream;
}

std::ostream& coords::operator<< (std::ostream &stream, Stereo const & stereo)
{
  for (auto const & center : stereo.centers())
  {
    stream << "Stereo, Atom: " << center.m_atom << ", Dir: " << center.m_dir << '\n';
  }
  return stream;
}


/* ######################################################

  class Stereo implementation

   ######  ######## ######## ########  ########  #######
  ##    ##    ##    ##       ##     ## ##       ##     ##
  ##          ##    ##       ##     ## ##       ##     ##
   ######     ##    ######   ########  ######   ##     ##
        ##    ##    ##       ##   ##   ##       ##     ##
  ##    ##    ##    ##       ##    ##  ##       ##     ##
   ######     ##    ######## ##     ## ########  #######

###################################################### */




coords::Stereo::Stereo(Atoms const &atoms, coords::Representation_3D const &xyz)
  : m_centers()
{
  using scon::cross;
  std::vector<Atom>::size_type const N = atoms.size();
  m_centers.clear();
  for (std::size_t i(0U); i < N; ++i)
  {
    if (atoms.atom(i).bonds().size() == 4)
    {
      pair stereo_pair;
      stereo_pair.m_atom = i;
      // Bond indices
      stereo_pair.bonds[0] = atoms.atom(i).bonds()[0];
      stereo_pair.bonds[1] = atoms.atom(i).bonds()[1];
      stereo_pair.bonds[2] = atoms.atom(i).bonds()[2];
      stereo_pair.bonds[3] = atoms.atom(i).bonds()[3];
      // Two hydrogens? 
      //if (atoms.atom(stereo_pair.bonds[0]).number() == 1 )
      // All of the residues need to differ if we have a stereo center
      if (
        !atoms.res_is_equal(stereo_pair.bonds[0], stereo_pair.bonds[1], i, i, 0) &&
        !atoms.res_is_equal(stereo_pair.bonds[0], stereo_pair.bonds[2], i, i, 0) &&
        !atoms.res_is_equal(stereo_pair.bonds[0], stereo_pair.bonds[3], i, i, 0) &&
        !atoms.res_is_equal(stereo_pair.bonds[1], stereo_pair.bonds[2], i, i, 0) &&
        !atoms.res_is_equal(stereo_pair.bonds[1], stereo_pair.bonds[3], i, i, 0) &&
        !atoms.res_is_equal(stereo_pair.bonds[2], stereo_pair.bonds[3], i, i, 0)
        )
      {
        coords::Cartesian_Point crossv(xyz[stereo_pair.bonds[1]] - xyz[stereo_pair.bonds[2]]);
        crossv = scon::cross(crossv, xyz[stereo_pair.bonds[3]] - xyz[stereo_pair.bonds[2]]);
        double const scalar(dot(xyz[stereo_pair.bonds[0]] - xyz[stereo_pair.bonds[2]], crossv));
        stereo_pair.m_dir = !(scalar < 0.0);
        m_centers.push_back(stereo_pair);
      }
    }
  }
}

void coords::Stereo::update(coords::Representation_3D const &xyz)
{
  for (auto & stereo_pair : m_centers)
  {
    coords::Cartesian_Point crossv(
      cross(xyz[stereo_pair.bonds[1]] - xyz[stereo_pair.bonds[2]],
        xyz[stereo_pair.bonds[3]] - xyz[stereo_pair.bonds[2]])
      );
    float_type const scalar(dot(xyz[stereo_pair.bonds[0]] -
      xyz[stereo_pair.bonds[2]], crossv));
    stereo_pair.m_dir = !(scalar < 0.0);
  }
}








/* ################################################################################################

 ######   #######   #######  ########  ########  #### ##    ##    ###    ######## ########  ######
##    ## ##     ## ##     ## ##     ## ##     ##  ##  ###   ##   ## ##      ##    ##       ##    ##
##       ##     ## ##     ## ##     ## ##     ##  ##  ####  ##  ##   ##     ##    ##       ##
##       ##     ## ##     ## ########  ##     ##  ##  ## ## ## ##     ##    ##    ######    ######
##       ##     ## ##     ## ##   ##   ##     ##  ##  ##  #### #########    ##    ##             ##
##    ## ##     ## ##     ## ##    ##  ##     ##  ##  ##   ### ##     ##    ##    ##       ##    ##
 ######   #######   #######  ##     ## ########  #### ##    ## ##     ##    ##    ########  ######

################################################################################################ */

/*! Constructor of empty Coordiantes object
 *
 * atoms, representations, stereo, potentials,
 * virial and stereo are empty. mult_struc_counter member
 * is set to zero (duh, it's an empty object).
 * However, keep in mind that the energy interface is already
 * initialized. And this will not change. The data for atoms etc.
 * will be filled in later using the function init_swap_in().
 * The energy interface and preinterface are created
 * according to the specifications from the global Config instance.
 *
 */
coords::Coordinates::Coordinates() :
  m_atoms(), m_representation(), m_stereo(),
  m_potentials(), m_virial(empty_virial()),
  m_interface(energy::new_interface(this)),
  m_preinterface(energy::pre_interface(this)),
  energy_valid(false),
  NEB_control(false),
  PathOpt_control(false),
  mult_struc_counter(0)
{
}

coords::Coordinates::Coordinates(Coordinates &&r) :
  m_atoms(std::move(r.m_atoms)),
  m_representation(std::move(r.m_representation)),
  m_stereo(std::move(r.m_stereo)),
  m_potentials(std::move(r.m_potentials)),
  m_virial(std::move(r.m_virial)),
  m_interface(r.m_interface->move(this)),
  m_preinterface(r.m_preinterface ? r.m_preinterface->move(this) : nullptr),
  energy_valid(r.energy_valid),
  NEB_control(r.NEB_control), PathOpt_control(r.PathOpt_control),
  mult_struc_counter(r.mult_struc_counter)
{
}

coords::Coordinates::Coordinates(Coordinates const &r) :
  m_atoms(r.m_atoms),
  m_representation(r.m_representation),
  m_stereo(r.m_stereo),
  m_potentials(r.m_potentials),
  m_virial(r.m_virial),
  m_interface(r.m_interface->clone(this)),
  m_preinterface(r.m_preinterface ? r.m_preinterface->clone(this) : nullptr),
  energy_valid(false),
  fep(r.fep),
  NEB_control(r.NEB_control),
  PathOpt_control(r.PathOpt_control),
  mult_struc_counter(r.mult_struc_counter)
{
}

coords::Coordinates& coords::Coordinates::operator=(Coordinates const & rhs)
{
  if (this != &rhs)
  {
    Coordinates tmp{ rhs };
    this->swap(tmp);
  }
  return *this;
}

coords::Coordinates& coords::Coordinates::operator=(Coordinates&& rhs)
{
  if (this != &rhs)
  {
    Coordinates tmp{ std::move(rhs) };
    this->swap(tmp);
  }
  return *this;
}

coords::Coordinates::~Coordinates()
{
  if (m_interface) delete m_interface;
  if (m_preinterface) delete m_preinterface;
}



void coords::Coordinates::set_fix(size_t const atom, bool const fix_it)
{
  fix(atom, fix_it);
}


void coords::Coordinates::init_swap_in(Atoms &a, PES_Point &p, bool const update)
{
  if (a.size() != p.size())
  {
	std::cout << "Can't swap in. Number of atoms "<< a.size() << " versus " << p.size() << ". Does not match.\n";
    throw std::logic_error("Initialization of Coordinates failed. Invalid sizes.");
  }
  energy_valid = false;
  a.swap(m_atoms);
  m_atoms.refine();
  p.swap(m_representation);
  Stereo(m_atoms, m_representation.structure.cartesian).swap(m_stereo);
  m_representation.ia_matrix.resize(subsystems().size());
  for (auto & interact : m_representation.ia_matrix)
  {
    interact.grad.resize(size());
    interact.energy = 0.0;
  }
  if (update)
  {
    std::size_t const N = m_atoms.size(), M = m_atoms.mains().size();
    m_representation.structure.cartesian.resize(N);
    m_representation.gradient.cartesian.resize(N);
    m_representation.structure.intern.resize(N);
    m_representation.gradient.intern.resize(N);
    m_representation.structure.main.resize(M);
    m_representation.gradient.main.resize(M);
    //m_atoms.c_to_i(m_representation); 
    m_interface->update(false);
    if (m_preinterface) m_preinterface->update(false);
  }
}


void coords::Coordinates::init_in(Atoms a, PES_Point p, bool const update)
{
  init_swap_in(a, p, update);
}

struct OCB
{
  coords::Coordinates * cp;
  OCB(coords::Coordinates & coordpointer) : cp(&coordpointer) {}
  float operator() (scon::vector<scon::c3<float>> const & v,
    scon::vector<scon::c3<float>> & g, std::size_t const S, bool & go_on)
  {
    cp->set_xyz(coords::Representation_3D(v.begin(), v.end()), false);
    float E = float(cp->g());
    go_on = cp->integrity();
    g.resize(cp->g_xyz().size());
    scon::explicit_transform(cp->g_xyz(), g);
    if (Config::get().general.verbosity >= 4)
    {
      std::cout << "Optimization: Energy of step " << S;
      std::cout << " is " << E << " integrity " << go_on << '\n';
    }
    return E;
  }
};

coords::float_type coords::Coordinates::lbfgs()
{
  using namespace  optimization::local;
  typedef coords::Container<scon::c3<float>> nc3_type;
  // Create optimizer
  auto optimizer = make_lbfgs(
    make_more_thuente(Coords_3d_float_callback(*this))
    );
  optimizer.ls.config.ignore_callback_stop = true;
  // Create Point
  using op_type = decltype(optimizer);
  op_type::point_type x(nc3_type(xyz().begin(), xyz().end()));
  // Optimize point
  optimizer.config.max_iterations =
    Config::get().optimization.local.bfgs.maxstep;
  optimizer.config.epsilon =
    (float)Config::get().optimization.local.bfgs.grad;
  optimizer(x);
  // Get structure, gradients and energy into coords
  m_representation.energy = optimizer.p().f;
  m_representation.structure.cartesian =
    coords::Representation_3D(optimizer.p().x.begin(), optimizer.p().x.end());
  //g();
  //std::cout << "Ene = " << m_representation.energy << "\n";
  m_representation.gradient.cartesian =
    coords::Gradients_3D(optimizer.p().g.begin(), optimizer.p().g.end());
  // Output
  if (Config::get().general.verbosity >= 4 ||
    (optimizer.state() < 0 && Config::get().general.verbosity >= 4))
  {
    std::cout << "Optimization done (status " << optimizer.state() <<
      "). Evaluations:" << optimizer.iter() << '\n';
  }
  if (Config::get().general.verbosity >= 4 && integrity())
  {
    std::cout << "Energy after optimization: \n";
    e_head_tostream_short(std::cout, energyinterface());
    e_tostream_short(std::cout, energyinterface());
  }
  // Return floating point
  return optimizer.p().f;
}

double coords::Coordinates::prelbfgs()
{
  using namespace  optimization::local;
  typedef coords::Container<scon::c3<float>> nc3_type;
  // Create optimizer
  auto optimizer = make_lbfgs(
    make_more_thuente(Coords_3d_float_pre_callback(*this))
    );
  // Create Point
  using op_type = decltype(optimizer);
  op_type::point_type x(nc3_type(xyz().begin(), xyz().end()));
  // Optimize point
  optimizer(x);
  m_representation.energy = x.f;
  m_representation.structure.cartesian =
    coords::Representation_3D(x.x.begin(), x.x.end());
  m_representation.gradient.cartesian =
    coords::Gradients_3D(x.g.begin(), x.g.end());
  if (Config::get().general.verbosity >= 4)
  {
    std::cout << "Optimization done (status " << optimizer.state() <<
      "). Evaluations:" << optimizer.iter() << '\n';
  }

  if (Config::get().general.verbosity >= 4 && m_interface->intact())
  {
    std::cout << "Energy after optimization: \n";
    e_head_tostream_short(std::cout, m_interface);
    e_tostream_short(std::cout, m_interface);
  }
  return x.f;
}

coords::Gradients_Main coords::Coordinates::dimermethod_dihedral
(std::vector<coords::Gradients_Main> const &D)
{
  if (this->atoms().mains().empty())
  {
    throw std::runtime_error("System does not contain any main dihedrals. Dimermethod cannot be applied.");
  }
  using Move_T = optimization::CG_DimerMover < Main_Callback >;
  Main_Callback C(*this);
  Move_T mover(0.1, 1000u);
  coords::Gradients_Main structure;
  std::size_t const N = m_representation.structure.main.size();
  structure.reserve(N);
  for (std::size_t i(0u); i < N; ++i)
  {
    structure.emplace_back(m_representation.structure.main[i].radians());
  }
  Move_T::minimum_type minimum(structure);
  minimum.directions = D;
  C = mover(minimum, C);
  return minimum.directions.back();
}


void coords::Coordinates::swap(Coordinates &rhs) // object swap
{
  m_atoms.swap(rhs.m_atoms);
  m_representation.swap(rhs.m_representation);
  m_stereo.swap(rhs.m_stereo);
  m_potentials.swap(rhs.m_potentials);
  m_virial.swap(rhs.m_virial);
  if (m_interface && rhs.m_interface) m_interface->swap(*rhs.m_interface);
  if (m_preinterface && rhs.m_preinterface) m_preinterface->swap(*rhs.m_preinterface);
  using std::swap;
  //m_sub_interaction.swap(rhs.m_sub_interaction);
  swap(energy_valid, rhs.energy_valid);
  swap(this->fep, rhs.fep);
  swap(this->hessian_def, rhs.hessian_def);
  swap(this->use_fep, rhs.use_fep);
  swap(this->mult_struc_counter, rhs.mult_struc_counter);
  swap(this->NEB_control, rhs.NEB_control);
  swap(this->orthogonalize, rhs.orthogonalize);
  swap(this->PathOpt_control, rhs.PathOpt_control);
}

coords::Cartesian_Point coords::Coordinates::center_of_mass() const
{
  coords::Cartesian_Point COM;
  std::vector<Atom>::size_type const N(m_atoms.size());
  coords::float_type M(0.0);
  for (std::vector<Atom>::size_type i(0U); i < N; ++i)
  {
    coords::float_type const mass(m_atoms.atom(i).mass());
    M += mass;
    COM += xyz(i)*mass;
  }
  COM /= M;
  return COM;
}


coords::Cartesian_Point coords::Coordinates::center_of_geometry() const
{
  std::size_t const N = xyz().size();
  coords::Cartesian_Point p(0);
  for (std::size_t i(0u); i < N; ++i)
  {
    p += xyz()[i];
  }
  p /= float_type(N);
  return p;
}


double coords::Coordinates::weight() const
{
  std::vector<Atom>::size_type const N(m_atoms.size());
  coords::float_type M(0.0);
  for (std::vector<Atom>::size_type i(0U); i < N; ++i)
  {
    M += m_atoms.atom(i).mass();
  }
  return M;
}


void coords::Coordinates::e_head_tostream_short(std::ostream &strm,
  energy::interface_base const * const ep) const
{
  if (ep) ep->print_E_head(strm);
  else m_interface->print_E_head(strm);
  if (!m_potentials.empty())
  {
    strm << "Bias Potentials: \n";
    strm << std::setw(24) << "DIH";
    strm << std::setw(24) << "ANG";
    strm << std::setw(24) << "DIST";
    strm << std::setw(24) << "SPHERICAL";
    strm << std::setw(24) << "CUBIC\n";
    strm << std::setw(24) << m_potentials.dihedrals().size();
    strm << std::setw(24) << m_potentials.angles().size();
    strm << std::setw(24) << m_potentials.distances().size();
    strm << std::setw(24) << m_potentials.sphericals().size();
    strm << std::setw(24) << m_potentials.cubic().size() << '\n';
  }
}


void coords::Coordinates::e_tostream_short(std::ostream &strm,
  energy::interface_base const * const ep) const
{
  if (ep) ep->print_E_short(strm);
  else if (m_interface) m_interface->print_E_short(strm);
  if (!m_potentials.empty())
  {
    strm << "Bias Energies: \n";
    strm << std::setw(24) << std::fixed << std::setprecision(8) << m_potentials.e_dihedral();
    strm << std::setw(24) << std::fixed << std::setprecision(8) << m_potentials.e_angle();
    strm << std::setw(24) << std::fixed << std::setprecision(8) << m_potentials.e_dist();
    strm << std::setw(24) << std::fixed << std::setprecision(8) << m_potentials.e_spherical();
    strm << std::setw(24) << std::fixed << std::setprecision(8) << m_potentials.e_cubic() << '\n';
  }
  strm << "Total energy: " << m_representation.energy << "\n";
  strm << '\n';
}


coords::Cartesian_Point coords::Coordinates::center_of_mass_mol(std::size_t index) const
{
  if (index >= m_atoms.molecules().size())
  {
    throw std::logic_error("Wrong molecule index for center_of_mass_mol() given.");
  }
  coords::Cartesian_Point COM = coords::Cartesian_Point();
  size_2d::size_type const N(m_atoms.molecules()[index].size());
  coords::float_type M = coords::float_type();
  for (std::vector<Atom>::size_type i(0U); i < N; ++i)
  {
    coords::float_type mass(m_atoms.atom(m_atoms.molecules(index, i)).mass());
    M += mass;
    COM += xyz(m_atoms.molecules(index, i))*mass;
  }
  return COM / M;
}


void coords::Coordinates::set_dih(size_type const int_index,
  coords::angle_type const target_angle,
  bool const move_dependants_along, bool const move_fixed_dih)
{
  if (int_index < size() && (!m_atoms.atom(int_index).ifix() || move_fixed_dih))
  {
    angle_type rot = target_angle - intern(int_index).azimuth();
    if (move_dependants_along)
    {
      std::size_t int_bond = m_atoms.atom(int_index).ibond();
      size_type const N = m_atoms.atom(int_bond).bound_internals().size();
      for (size_type i(0u); i < N; ++i)
      {
        m_representation.structure.intern[m_atoms.atom(int_bond).bound_internals()[i]].azimuth() += rot;
      }
    }
    else
    {
      m_representation.structure.intern[int_index].azimuth() = target_angle;
    }
  }
}


void coords::Coordinates::rotate_dih(size_type const int_index,
  coords::angle_type const rot_angle,
  bool const move_dependants_along, bool const move_fixed_dih)
{
  if (int_index < size() && (!m_atoms.atom(int_index).ifix() || move_fixed_dih))
  {
    if (move_dependants_along)
    {
      // If we rotate all we need to rotate every dihedral of every atom attached to the
      // atom to which the selected one is attached
      std::size_t int_bond = m_atoms.atom(int_index).ibond();
      size_type const N = m_atoms.atom(int_bond).bound_internals().size();
      for (size_type i(0u); i < N; ++i)
      {
        m_representation.structure.intern[m_atoms.atom(int_bond).bound_internals()[i]].azimuth() += rot_angle;
      }
    }
    else
    {
      m_representation.structure.intern[int_index].azimuth() += rot_angle;
    }
  }
}


void coords::Coordinates::rotate_main(size_type const main_index,
  coords::angle_type const rot_angle,
  bool const move_dependants_along, bool const move_fixed_dih)
{
  rotate_dih(m_atoms.intern_of_main_idihedral(main_index), rot_angle, move_dependants_along, move_fixed_dih);
}


void coords::Coordinates::set_all_main(coords::Representation_Main const & new_values,
  bool const apply_to_xyz, bool const move_dependants_along, bool const move_fixed_dih)
{
  size_type const N(main().size());
  if (new_values.size() != N)
  {
    throw std::logic_error("set_all_main called with wrong-sized vector");
  }
  for (size_type i(0u); i<N; ++i)
  {
    set_dih(m_atoms.intern_of_main_idihedral(i), new_values[i],
      move_dependants_along, move_fixed_dih);
    m_representation.structure.main[i] = new_values[i];
  }
  if (apply_to_xyz) to_xyz();
}

void coords::Coordinates::periodic_boxjump()
{
	std::size_t const N(molecules().size());
	Cartesian_Point const halfbox(Config::get().energy.pb_box / 2.0);
	for (std::size_t i = 0; i < N; ++i)
	{
		Cartesian_Point tmp_com(-center_of_mass_mol(i));
		bool move = false;
		if (std::abs(tmp_com.x() <= halfbox.x())) 
		{
			tmp_com.x() = 0;  
		}
		else
		{
			//std::cout << "x_vorher: " << -tmp_com.x() <<" Molekuel: "<<molecules(i)<< "\n";
			move = true;
			tmp_com.x() = tmp_com.x() / Config::get().energy.pb_box.x();
			int tmp_x = (int)tmp_com.x();
			if (tmp_x > 0)
			{
				tmp_com.x() = tmp_x;
			}
			else
			{
				tmp_com.x() = tmp_x + 1;
			}
		}

		if (std::abs(tmp_com.y() <= halfbox.y()))
		{
			tmp_com.y() = 0;
		}
		else
		{
			move = true;
			//std::cout << "y_vorher: " << -tmp_com.y() << " Molekuel: " << molecules(i) << "\n";
			tmp_com.y() = tmp_com.y() / Config::get().energy.pb_box.y();
			int tmp_y = (int)tmp_com.y();
			if (tmp_y > 0)
			{
				tmp_com.y() = tmp_y;
			}
			else
			{
				tmp_com.y() = tmp_y + 1;
			}
		}

		if (std::abs(tmp_com.z() <= halfbox.z()))
		{
			tmp_com.z() = 0;
		}
		else
		{
			move = true;
			//std::cout << "z_vorher: " << -tmp_com.z() << " Molekuel: " << molecules(i) << "\n";
			tmp_com.z() = tmp_com.z() / Config::get().energy.pb_box.z();
			int tmp_z = (int)tmp_com.z();
			if (tmp_z > 0)
			{
				tmp_com.z() = tmp_z;
			}
			else
			{
				tmp_com.z() = tmp_z + 1;
			}
		}
		tmp_com *= Config::get().energy.pb_box;
    for (auto const atom : molecules(i)) move_atom_by(atom, tmp_com, true);
	//if (move == true)
	//{
		//std::cout << "nachher " << center_of_mass_mol(i) << "\n";
	//}
  }
}


bool coords::Coordinates::validate_bonds()
{
  bool status = true;
  std::size_t const N(m_atoms.size());
  for (std::size_t i = 0; i < N; ++i)  // for every atom i
  {
    for (auto const & bound : m_atoms.atom(i).bonds())  // for every atom b that is bound to i
    {
      double const L(geometric_length(xyz(i) - xyz(bound)));
	  if (L < 0.3 || L > 5.0)  // test if bondlength between i and b is reasonable
	  {
		  status = false;  
		  if (i < bound)   // save all bonds with strange bondlengths in broken_bonds
		  {
			std::vector<float> bond;
			bond.push_back(i);
			bond.push_back(bound);
			bond.push_back(L);
			broken_bonds.push_back(bond);
		  }  
	  }
    }
  }
  return status;
}


void coords::Coordinates::set_pes(PES_Point const & pes_point,
  bool const overwrite_fixed /*= false*/)
{
  if (!m_representation.same_size(pes_point))
  {
    std::cout << "CS: " << (m_representation.structure.cartesian.size() ==
      pes_point.structure.cartesian.size()) <<
      "(" << m_representation.structure.cartesian.size() <<
      "," << pes_point.structure.cartesian.size() << ")\n";
    std::cout << "CG: " << (m_representation.gradient.cartesian.size() ==
      pes_point.gradient.cartesian.size()) <<
      "(" << m_representation.gradient.cartesian.size() <<
      "," << pes_point.gradient.cartesian.size() << ")\n";
    std::cout << "IS: " << (m_representation.structure.intern.size() ==
      pes_point.structure.intern.size()) <<
      "(" << m_representation.structure.intern.size() <<
      "," << pes_point.structure.intern.size() << ")\n";
    std::cout << "IG: " << (m_representation.gradient.intern.size() ==
      pes_point.gradient.intern.size()) <<
      "(" << m_representation.gradient.intern.size() <<
      "," << pes_point.gradient.intern.size() << ")\n";
    std::cout << "MS: " << (m_representation.structure.main.size() ==
      pes_point.structure.main.size()) <<
      "(" << m_representation.structure.main.size() <<
      "," << pes_point.structure.main.size() << ")\n";
    std::cout << "MG: " << (m_representation.gradient.main.size() ==
      pes_point.gradient.main.size()) <<
      "(" << m_representation.gradient.main.size() <<
      "," << pes_point.gradient.main.size() << ")\n";
    throw std::logic_error("Invalid PES size for coordinates.");
  }
  if (overwrite_fixed)
  {
    m_representation = pes_point;
    m_stereo.update(xyz());
  }
  else
  {
    std::size_t const N(size());
    m_representation.gradient.cartesian.resize(pes_point.gradient.cartesian.size());
    m_representation.structure.cartesian.resize(pes_point.structure.cartesian.size());
    for (std::size_t i(0U); i < N; ++i)
    {
      if (atoms(i).fixed())
      {
        m_representation.gradient.cartesian[i] = coords::cartesian_gradient_type();
      }
      else
      {
        m_representation.structure.cartesian[i] = pes_point.structure.cartesian[i];
        m_representation.gradient.cartesian[i] = pes_point.gradient.cartesian[i];
      }
      if (atoms(i).ifix())
      {
        m_representation.gradient.intern[i] = coords::internal_gradient_type();
      }
      else
      {
        m_representation.gradient.intern[i] = pes_point.gradient.intern[i];
        m_representation.structure.intern[i] = pes_point.structure.intern[i];
      }
    }
    std::size_t const MG(pes_point.gradient.main.size()),
      MS(pes_point.structure.main.size());
    m_representation.gradient.main.resize(MG);
    m_representation.structure.main.resize(MS);
    for (std::size_t i(0U); i < MG && i < MS; ++i)
    {
      std::size_t const ii = atoms().intern_of_main_idihedral(i);
      if (atoms(ii).ifix())
      {
        m_representation.gradient.main[i] = coords::main_gradient_type();
      }
      else
      {
        m_representation.gradient.main[i] = pes_point.gradient.main[i];
        m_representation.structure.main[i] = pes_point.structure.main[i];
      }
    }
    m_stereo.update(xyz());
  }
}


bool coords::Coordinates::check_superposition_xyz(Representation_3D const &a, Representation_3D const &b, double const x /*= 0.35*/) const
{
  std::size_t const N(a.size());
  if (size() != N || b.size() != N)
    throw std::logic_error("Evaluating the structural overlap between different sized structures.");
  for (std::size_t i(0u); i < N; ++i)
  {
    bool is_superposed(false);
    for (std::size_t j(0u); j < N; ++j)
    {
      if (j != i && m_atoms.atom(i).number() == m_atoms.atom(j).number())
      {
        if (scon::geometric_length(a[i] - b[j]) < x) is_superposed = true;
      }
    }
    if (!is_superposed) return false;
  }
  return true;
}

void coords::Coordinates::set_pes(PES_Point && pes_point,
  bool const overwrite_fixed /*= false*/)
{
  bool const valid_size = scon::equal_size_ranges(pes_point.structure.cartesian,
    pes_point.gradient.cartesian,
    pes_point.structure.intern,
    pes_point.gradient.intern,
    m_representation.structure.cartesian,
    m_representation.gradient.cartesian,
    m_representation.structure.intern,
    m_representation.gradient.intern) &&
    scon::equal_size_ranges(pes_point.structure.main,
      m_representation.structure.main,
      pes_point.gradient.main,
      m_representation.gradient.main);

  if (!valid_size)
  {
    throw std::logic_error("Invalid moved PES size for coordinates.");
  }
  m_representation.swap(pes_point);
  m_representation.resize(pes_point.structure.cartesian.size(),
    pes_point.structure.main.size());
  if (!overwrite_fixed)
  {
    std::size_t const N(m_representation.structure.cartesian.size());
    for (std::size_t i(0U); i < N; ++i)
    {
      if (atoms(i).fixed())
      {
        m_representation.structure.cartesian[i] = pes_point.structure.cartesian[i];
        m_representation.gradient.cartesian[i] = coords::cartesian_gradient_type();
      }
      if (atoms(i).ifix())
      {
        m_representation.structure.intern[i] = pes_point.structure.intern[i];
        m_representation.gradient.intern[i] = coords::internal_gradient_type();
      }
    }
    std::size_t const M(pes_point.structure.main.size());
    for (std::size_t i(0U); i < M; ++i)
    {
      std::size_t const ii = atoms().intern_of_main_idihedral(i);
      if (atoms(ii).ifix())
      {
        m_representation.structure.main[i] = pes_point.structure.main[i];
        m_representation.gradient.main[i] = coords::main_gradient_type();
      }
    }
  }
  m_stereo.update(xyz());
}

/* #############################################################################

class PES_Point and helper static function implementations

########  ########  ######          ########   #######  #### ##    ## ########
##     ## ##       ##    ##         ##     ## ##     ##  ##  ###   ##    ##
##     ## ##       ##               ##     ## ##     ##  ##  ####  ##    ##
########  ######    ######          ########  ##     ##  ##  ## ## ##    ##
##        ##             ##         ##        ##     ##  ##  ##  ####    ##
##        ##       ##    ##         ##        ##     ##  ##  ##   ###    ##
##        ########  ######  ####### ##         #######  #### ##    ##    ##

############################################################################# */

namespace
{
  bool out_of_angle_delta(coords::angle_type const & a,
    coords::angle_type const & l)
  {
    return (a > l || a < -l);
  }

}

bool coords::Coordinates::equal_structure(coords::PES_Point const & a, 
  coords::PES_Point const & b, 
  coords::main_type const md, 
  coords::internal_type const & id, 
  coords::Cartesian_Point const & cd) const
{
  using scon::operator-;
  auto const N = size();
  if (Config::get().coords.decouple_internals)
  {
    throw std::runtime_error("decouple_internals option "
      "may not be present to compare structure equality");
  }
  // Check main structure equality
  {
    bool equal = true;
    std::size_t i = 0;
    for (auto && d : (a.structure.main - b.structure.main))
    {
      // main only relevant if relative to atom (internal structure)
      auto in_of_main = m_atoms.intern_of_main_idihedral(i);
      auto rel_d = m_atoms.atom(in_of_main).idihedral();
      auto relevant = rel_d < N;
      if (relevant && (d > md || d < md))
      {
        equal = false;
        break;
      }
      ++i;
    }
    if (equal) return true;
  }
  // Check internal structure equality
  {
    bool equal = true;
    std::size_t i = 0;
    for (auto && d : (a.structure.intern - b.structure.intern))
    {
      auto rel_b = m_atoms.atom(i).ibond();
      auto rel_a = m_atoms.atom(i).iangle();
      auto rel_d = m_atoms.atom(i).idihedral();
      auto relevant = rel_b < N && rel_a < N && rel_d < N;
      if (relevant && (std::abs(d.radius()) > id.radius() ||
        out_of_angle_delta(d.inclination(), id.inclination()) ||
        out_of_angle_delta(d.azimuth(), id.azimuth())))
      {
        equal = false;
        break;
      }
      ++i;
    }
    if (equal) return true;
  }
  {
    bool equal = true;
    std::size_t i = 0;
    for (auto && r : (a.structure.cartesian - b.structure.cartesian))
    {
      if (std::abs(r.x()) > cd.x() ||
        std::abs(r.y()) > cd.y() ||
        std::abs(r.z()) > cd.z())
      {
        equal = false;
        break;
      }
      ++i;
    }
    if (equal) return true;
  }

  return false;
}

coords::float_type coords::Internal_Callback::operator()
(scon::vector<coords::float_type> const & v,
  scon::vector<coords::float_type>& g, std::size_t const S, bool & go_on)
{
  std::size_t i = 0;
  coords::Representation_Internal rin(v.size() / 3);
  for (auto & e : rin)
  {
    e.radius() = v[i++];
    e.inclination() = coords::angle_type(v[i++]);
    e.azimuth() = coords::angle_type(v[i++]);
  }
  cp->set_internal(rin);
  cp->to_xyz();
  float E = float(cp->g());
  cp->to_internal();
  go_on = cp->integrity();
  g.resize(v.size());
  i = 0;
  for (auto const & e : cp->g_intern())
  {
    g[i++] = e.x();
    g[i++] = e.y();
    g[i++] = e.z();
  }
  if (Config::get().general.verbosity >= 4)
  {
    std::cout << "Optimization: Energy of step " << S;
    std::cout << " is " << E << " integrity " << go_on << '\n';
  }
  return E;
}

coords::Gradients_Internal coords::Internal_Callback::from(
  coords::Representation_Internal const & v)
{
  std::size_t const N = v.size();
  coords::Gradients_Internal r;
  r.reserve(N);
  for (std::size_t i(0u); i < N; ++i)
  {
    r.emplace_back(v[i].radius(),
      v[i].inclination().radians(),
      v[i].azimuth().radians());
  }
  return r;
}
coords::Representation_Internal coords::Internal_Callback::to(
  scon::vector<scon::c3<coords::float_type>> const & v)
{
  std::size_t const N = v.size();
  coords::Representation_Internal r;
  r.reserve(N);
  for (std::size_t i(0u); i < N; ++i)
  {
    r.emplace_back(v[i].x(),
      scon::ang<coords::float_type>::from_rad(v[i].y()),
      scon::ang<coords::float_type>::from_rad(v[i].z()));
  }
  return r;
}

coords::float_type coords::Main_Callback::operator() (coords::Gradients_Main const & v,
  coords::Gradients_Main & g, std::size_t const S, bool & go_on)
{
  std::size_t const N = v.size();
  coords::Representation_Main r;
  r.reserve(N);
  for (std::size_t i(0u); i < N; ++i)
  {
    r.emplace_back(coords::main_type::from_rad(v[i]));
  }
  cp->set_all_main(r);
  auto E = cp->g();
  cp->to_internal();
  go_on = cp->integrity();
  g = cp->g_main();
  if (Config::get().general.verbosity >= 4)
  {
    std::cout << "Optimization: Energy of step " << S;
    std::cout << " is " << E << " integrity " << go_on << '\n';
  }
  return E;
}

void coords::cartesian_logfile_drain::operator() (coords::Representation_3D && xyz)
{
  if (!cp || !strm) return;
  // Save current state
  auto const tmp = (*cp).xyz();
  // plug snapshot into coords
  (*cp).set_xyz(std::move(xyz));
  // optimize snapshot
  if (opt) { (*cp).o(); }
  // Print to stream
  *strm << *cp;
  // reset state
  (*cp).set_xyz(std::move(tmp));
}

coords::offset_buffered_cartesian_logfile coords::make_buffered_cartesian_log(Coordinates &c,
  std::string file_suffix, std::size_t buffer_size,
  std::size_t log_offset, bool optimize)
{
  return scon::offset_call_buffer<coords::Representation_3D>(buffer_size, log_offset,
    cartesian_logfile_drain{ c, output::filename(file_suffix).c_str(), optimize });
}

float coords::Coords_3d_float_pre_callback::operator() (scon::vector<scon::c3<float>> const & v,
  scon::vector<scon::c3<float>> & g, std::size_t const S, bool & go_on)
{
  cp->set_xyz(coords::Representation_3D(v.begin(), v.end()));
  float E = float(cp->pg());
  go_on = cp->integrity();
  g = scon::vector<scon::c3<float>>(cp->g_xyz().begin(), cp->g_xyz().end());
  if (Config::get().general.verbosity >= 4)
    std::cout << "Optimization: Energy of step " <<
    S << " is " << E << " integrity " << go_on << '\n';
  return E;
}

coords::Representation_3D coords::Coords_3d_float_callback::to(scon::vector<scon::c3<float>> const & v)
{
  return coords::Representation_3D(v.begin(), v.end());
}

scon::vector<scon::c3<float>> coords::Coords_3d_float_callback::from(coords::Gradients_3D const & g)
{
  scon::vector<scon::c3<float>> r(g.size());
  std::transform(g.begin(), g.end(), r.begin(),
    [](coords::Cartesian_Point const &p) -> scon::c3<float>
  {
    return scon::c3<float>(static_cast<float>(p.x()),
      static_cast<float>(p.y()), static_cast<float>(p.z()));
  });
  return r;
}


float coords::Coords_3d_float_callback::operator() (scon::vector<scon::c3<float>> const & v,
  scon::vector<scon::c3<float>> & g, std::size_t const S, bool & go_on)
{
  cp->set_xyz(to(v), false);
  float E = float(cp->g());
  go_on = cp->integrity();
  g = from(cp->g_xyz());
  if (Config::get().general.verbosity >= 4)
  {
    std::cout << "Optimization: Energy of step " << S;
    std::cout << " is " << E << " integrity " << go_on << '\n';
    //std::cout << "totg " << scon::vector_delimeter('\n') << g << "\n";
  }
  return E;
}




//##################################################################################
//##################################################################################
//                    SMOOTH PARTICLE MESH EWALD
//##################################################################################
//##################################################################################

//void coords::Coordinates::pme_stuff(int atoms)
//{
//  pme.pmetemp.natoms = atoms;
//  double x, nenner;
//  std::vector<double> stupidarray(250);
//  std::vector<double> brecurs(6);
//  double eps = 1e-8;
//  double tempcoeff;
//  double rate, j, l, blow, bhigh;
//  //Define PME grid size. Must be even number with prime factores 2, 3 and 5
//  // get even numbers factorized by 2, 3 and 5
//  for (int i = 2; i < 600; i++)
//  {
//    int n = i;
//    // check if numnber is even
//    if (n % 2 != 0) continue;
//    // divide by 2 as long as possible
//    while (n % 2 == 0)
//    {
//      n = n / 2;
//    }
//    // get prime factors 3 and 5
//    for (int j = 3; j <= 5; j = j + 2)
//    {
//      // divide by 3 as long as possible, then by 5 as long as possible
//      while (n%j == 0)
//      {
//        n = n / j;
//      }
//    }
//    // if resulting number is 1 push number i into vector
//    if (n == 1) pme.pmetemp.gridnumbers.push_back(i);
//  }
//  //Get initial Ewald coefficient based on cutoff radius
//  rate = eps + 1.0;
//  pme.pmetemp.ewaldcoeff = 0.5;
//  j = 0;
//  // trial and error function for first initial guess
//  while (rate > eps){
//    j += 1;
//    pme.pmetemp.ewaldcoeff *= 2.0;
//    l = pme.pmetemp.ewaldcoeff * Config::get().energy.cutoff;
//    rate = erfc(l) / Config::get().energy.cutoff;
//  }
//  //run binary search according to Sander program to refine coefficient
//  blow = 0;
//  bhigh = pme.pmetemp.ewaldcoeff;
//  for (int i = 1; i <= 50; i++)
//  {
//    pme.pmetemp.ewaldcoeff = (blow + bhigh) / 2;
//    if (erfc(pme.pmetemp.ewaldcoeff*Config::get().energy.cutoff) / Config::get().energy.cutoff > eps)
//    {
//      blow = pme.pmetemp.ewaldcoeff;
//    }
//    else bhigh = pme.pmetemp.ewaldcoeff;
//  }
//  tempcoeff = pme.pmetemp.ewaldcoeff;
//  //Calculate the grid size from the periodic box dimensions
//  pme.pmetemp.treshold = 1e-8;
//  pme.pmetemp.dim1 = (Config::get().energy.pb_box.x() * 1.2 - pme.pmetemp.treshold) + 1;
//  pme.pmetemp.dim2 = (Config::get().energy.pb_box.y() * 1.2 - pme.pmetemp.treshold) + 1;
//  pme.pmetemp.dim3 = (Config::get().energy.pb_box.z() * 1.2 - pme.pmetemp.treshold) + 1;
//  // check that grid size is even-numbered for effective fftw
//  pme.pmetemp.nxpoints = pme.pmetemp.nypoints = pme.pmetemp.nzpoints = pme.pmetemp.fftgridmax;
//  // get grid size in each dimension
//  for (int i = pme.pmetemp.gridnumbers.size(); i >= 1; i--){
//    int temp;
//    temp = pme.pmetemp.gridnumbers[i];
//    if (temp < 250)
//    {
//      if (temp >= pme.pmetemp.dim1) pme.pmetemp.nxpoints = temp;
//      if (temp >= pme.pmetemp.dim2) pme.pmetemp.nypoints = temp;
//      if (temp >= pme.pmetemp.dim3) pme.pmetemp.nzpoints = temp;
//    }
//  }
//  // compute the moduli for the inverse discrete fourier transformation
//  // ########
//  // B-Spline Stuff
//  // ########
//  x = 0.0;
//  brecurs[1] = (1.0 - x);
//  brecurs[2] = x;
//  // recursive loop over the b-spline order 1 to x (here 5) to get spline coefficients
//  for (int k = 3; k <= pme.pmetemp.bsplineorder; k++)
//  {
//    nenner = 1.0 / (k - 1);
//    brecurs[k] = x * brecurs[k - 1] * nenner;
//    for (int i = 1; i <= (k - 2); i++)
//    {
//      brecurs[k - i] = ((x + i)*brecurs[k - i - 1] + (k - i - x)*brecurs[k - i]) * nenner;
//    }
//    brecurs[1] = (1.0 - x) * brecurs[1] * nenner;
//  }
//  // ########
//  // Fill spline array
//  // #######
//  for (int i = 0; i < 250; i++)
//  {
//    stupidarray[i] = 0.0;
//  }
//  // i = 2 till b-spline order + 1
//  for (int i = 2; i <= pme.pmetemp.bsplineorder + 1; i++)
//  {
//    stupidarray[i] = brecurs[i - 1];
//  }
//  // Calculate the modulus for the spline arrays
//  coords::Coordinates::pme_dftmodulus(stupidarray);
//  // Set chunks for parallel execution
//#ifdef _OPENMP
//  coords::Coordinates::roughgrid();
//#endif
//  // Allocate memory for needed arrays
//  pme.pmetemp.bscx.Allocate(4, pme.pmetemp.bsplineorder, pme.pmetemp.natoms);
//  pme.pmetemp.bscy.Allocate(4, pme.pmetemp.bsplineorder, pme.pmetemp.natoms);
//  pme.pmetemp.bscz.Allocate(4, pme.pmetemp.bsplineorder, pme.pmetemp.natoms);
//  size_t align = sizeof(Complex);
//  pme.pmetemp.charges.Allocate(pme.pmetemp.nxpoints, pme.pmetemp.nypoints, pme.pmetemp.nzpoints, align);
//  pme.pmetemp.fractionalcharges.Allocate(pme.pmetemp.nxpoints, pme.pmetemp.nypoints, pme.pmetemp.nzpoints);
//  pme.pmetemp.recivectors.Allocate(3, 3);
//  pme.pmetemp.initgrid.Allocate(3, pme.pmetemp.natoms);
//  pme.pmetemp.parallelpme.Allocate(pme.pmetemp.natoms, pme.pmetemp.rgridtotal);
//  pme.pmetemp.feinf.resize(pme.pmetemp.natoms);
//
//  // set up FFTW variables if OPENMP is used
//#ifdef _OPENMP
//  int threads;
//#pragma omp parallel
//  {
//    threads = omp_get_num_threads();
//  }
//  fftw_init_threads();
//  fftw_plan_with_nthreads(threads);
//#endif
//  // Assign memory for charge grid
//  pme.pmetemp.in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* (pme.pmetemp.nxpoints*pme.pmetemp.nypoints*pme.pmetemp.nzpoints));
//  // Set up FFTW plans
//  pme.pmetemp.forward = fftw_plan_dft_3d(pme.pmetemp.nxpoints, pme.pmetemp.nypoints, pme.pmetemp.nzpoints, pme.pmetemp.in, pme.pmetemp.in, FFTW_FORWARD, FFTW_PATIENT);
//  pme.pmetemp.backward = fftw_plan_dft_3d(pme.pmetemp.nxpoints, pme.pmetemp.nypoints, pme.pmetemp.nzpoints, pme.pmetemp.in, pme.pmetemp.in, FFTW_BACKWARD, FFTW_PATIENT);
//  // Calculate reciprocal lattice vectors
//  double boxvolume = Config::get().energy.pb_box.x() * Config::get().energy.pb_box.y() * Config::get().energy.pb_box.z();
//  double x1, x2, x3, y1, y2, y3, z1, z2, z3;
//  x1 = Config::get().energy.pb_box.x();
//  x2 = x3 = 0.0;
//  y1 = y3 = 0.0;
//  y2 = Config::get().energy.pb_box.y();
//  z1 = z2 = 0.0;
//  z3 = Config::get().energy.pb_box.z();
//  pme.pmetemp.recivectors(0, 0) = (y2*z3 - z2*y3) / boxvolume;
//  pme.pmetemp.recivectors(1, 0) = (y3*z1 - z3*y1) / boxvolume;
//  pme.pmetemp.recivectors(2, 0) = (y1*z2 - z1*y1) / boxvolume;
//  pme.pmetemp.recivectors(0, 1) = (z2*x3 - x2*z3) / boxvolume;
//  pme.pmetemp.recivectors(1, 1) = (z3*x1 - z1*x3) / boxvolume;
//  pme.pmetemp.recivectors(2, 1) = (z1*x2 - z2*x1) / boxvolume;
//  pme.pmetemp.recivectors(0, 2) = (x2*y3 - y2*x3) / boxvolume;
//  pme.pmetemp.recivectors(1, 2) = (x3*y1 - y3*x1) / boxvolume;
//  pme.pmetemp.recivectors(2, 2) = (x1*y2 - y1*x2) / boxvolume;
//  // set up charge array
//  pme.pmetemp.atomcharges.resize(pme.pmetemp.natoms);
//  std::ifstream ifs2;
//  ifs2.open(Config::get().general.paramFilename.c_str());
//  int iterator = 0;
//  char buffer[150];
//  std::string bufferc, tempname;
//  std::vector<std::string> tokens;
//  std::string bla;
//  bla = Config::get().general.paramFilename.substr(0, 5);
//  pme.pmetemp.elecfac = 1.0;
//  if (bla.substr(0, 3) == "cha") pme.pmetemp.elecfac = 332.0716;
//  else if (bla.substr(0, 3) == "opl") pme.pmetemp.elecfac = 332.06;
//  //else if (bla.substr(0, 3) == "amb") pme.pmetemp.elecfac = 1.0;
//  //else if (bla.substr(0, 3) == "amo") pme.pmetemp.elecfac
//}
//
//// get FEP data for the three different grids (IN, OUT, ALL) and set flags for total charge vector
//void coords::Coordinates::getfepinfo()
//{
//  std::ifstream ifs;
//  int iterator = 0;
//  char buffer[150];
//  std::string bufferc;
//  std::vector<std::string> tokens;
//
//  ifs.open(Config::get().general.inputFilename.c_str());
//  while (!ifs.eof())
//  {
//    ifs.getline(buffer, 150);
//    bufferc = buffer;
//    tokens.clear();
//    std::istringstream iss(bufferc);
//    std::copy(std::istream_iterator <std::string>(iss), std::istream_iterator <std::string>(), std::back_inserter <std::vector <std::string >>(tokens));
//    int g = tokens.size();
//    if (tokens.size() > 2)
//    {
//
//      if (tokens[g - 1] == "IN")
//      {
//        pme.pmetemp.fepi.push_back(atoi(tokens[0].c_str()));
//        pme.pmetemp.feinf[iterator].flag = 1;
//      }
//      else if (tokens[g - 1] == "OUT")
//      {
//        pme.pmetemp.fepo.push_back(atoi(tokens[0].c_str()));
//        pme.pmetemp.feinf[iterator].flag = 0;
//      }
//      else
//      {
//        pme.pmetemp.fepi.push_back(atoi(tokens[0].c_str()));
//        pme.pmetemp.fepo.push_back(atoi(tokens[0].c_str()));
//        pme.pmetemp.fepa.push_back(atoi(tokens[0].c_str()));
//        pme.pmetemp.feinf[iterator].flag = 2;
//      }
//      iterator += 1;
//    }
//
//  }
//}
//
//// Precalculate the DFT moduli
//void coords::Coordinates::pme_dftmodulus(std::vector<double> & stupidarray)
//{
//  double eps = 1e-7;
//  pme.pmetemp.moduli1.resize(250);
//  pme.pmetemp.moduli2.resize(250);
//  pme.pmetemp.moduli3.resize(250);
//  double sinsum, cossum, fact, sumarg, numbersum, zeta, tempsum1, tempsum2;
//  int acutoff, doubleorder, indexi, indexk, indexj;
//  fact = 2.0 * SUPERPI / pme.pmetemp.nxpoints;
//  for (int i = 0; i < pme.pmetemp.nxpoints; i++)
//  {
//    sinsum = cossum = 0.0;
//    for (int k = 0; k < pme.pmetemp.nxpoints; k++)
//    {
//      numbersum = double(i*k);
//      sumarg = fact * numbersum;
//      cossum += stupidarray[k] * cos(sumarg);
//      sinsum += stupidarray[k] * sin(sumarg);
//    }
//    pme.pmetemp.moduli1[i] = sinsum*sinsum + cossum*cossum;
//  }
//  // Coorection for the euler interpolation
//  double tresh = 1e-7;
//  if (pme.pmetemp.moduli1[0] < eps) pme.pmetemp.moduli1[0] = 0.5 * pme.pmetemp.moduli1[1];
//  for (int i = 1; i < pme.pmetemp.nxpoints - 1; i++)
//  {
//    if (pme.pmetemp.moduli1[i] < eps) pme.pmetemp.moduli1[i] = 0.5 * (pme.pmetemp.moduli1[i - 1] + pme.pmetemp.moduli1[i + 1]);
//  }
//  if (pme.pmetemp.moduli1[pme.pmetemp.nxpoints] < eps)  pme.pmetemp.moduli1[pme.pmetemp.nxpoints] = 0.5 * pme.pmetemp.moduli1[pme.pmetemp.nxpoints - 1];
//  // calculate some stupid factor someone called zeta...whatever this thing does....
//  acutoff = 50;
//  doubleorder = 2 * 5;
//  for (int i = 0; i < pme.pmetemp.nxpoints; i++)
//  {
//    indexk = i;
//    if (i > pme.pmetemp.nxpoints / 2) indexk -= pme.pmetemp.nxpoints;
//    if (indexk == 0) zeta = 1.0;
//    else
//    {
//      tempsum1 = tempsum2 = 1.0;
//      fact = SUPERPI * double(indexk) / double(pme.pmetemp.nxpoints);
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact + SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact - SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      zeta = tempsum2 / tempsum1;
//    }
//    pme.pmetemp.moduli1[i] *= (zeta*zeta);
//  }
//  // ###################
//  // SECOND MODULUS ARRAY
//  //####################
//  fact = 2.0 * SUPERPI / pme.pmetemp.nypoints;
//  for (int i = 0; i < pme.pmetemp.nypoints; i++)
//  {
//    sinsum = cossum = 0.0;
//    for (int k = 0; k < pme.pmetemp.nypoints; k++)
//    {
//      numbersum = double(i*k);
//      sumarg = fact * numbersum;
//      cossum += stupidarray[k] * cos(sumarg);
//      sinsum += stupidarray[k] * sin(sumarg);
//    }
//    pme.pmetemp.moduli2[i] = sinsum*sinsum + cossum*cossum;
//  }
//  // Coorection for the euler interpolation
//  if (pme.pmetemp.moduli2[0] < eps) pme.pmetemp.moduli2[0] = 0.5 * pme.pmetemp.moduli2[1];
//  for (int i = 1; i < pme.pmetemp.nypoints - 1; i++)
//  {
//    if (pme.pmetemp.moduli2[i] < eps) pme.pmetemp.moduli2[i] = 0.5 * (pme.pmetemp.moduli2[i - 1] + pme.pmetemp.moduli2[i + 1]);
//  }
//  if (pme.pmetemp.moduli2[pme.pmetemp.nypoints] < eps)  pme.pmetemp.moduli2[pme.pmetemp.nypoints] = 0.5 * pme.pmetemp.moduli2[pme.pmetemp.nypoints - 1];
//  // calculate some stupid factor someone called zeta...whatever this thing does....
//  acutoff = 50;
//  doubleorder = 2 * 5;
//  for (int i = 0; i < pme.pmetemp.nypoints; i++)
//  {
//    indexk = i;
//    if (i > pme.pmetemp.nypoints / 2) indexk -= pme.pmetemp.nypoints;
//    if (indexk == 0) zeta = 1.0;
//    else
//    {
//      tempsum1 = tempsum2 = 1.0;
//      fact = SUPERPI * double(indexk) / double(pme.pmetemp.nypoints);
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact + SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact - SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      zeta = tempsum2 / tempsum1;
//    }
//    pme.pmetemp.moduli2[i] *= (zeta*zeta);
//  }
//  // ###################
//  // THIRD MODULUS ARRAY
//  //####################
//  fact = 2.0 * SUPERPI / pme.pmetemp.nzpoints;
//  for (int i = 0; i < pme.pmetemp.nzpoints; i++)
//  {
//    sinsum = cossum = 0.0;
//    for (int k = 0; k < pme.pmetemp.nzpoints; k++)
//    {
//      numbersum = double(i*k);
//      sumarg = fact * numbersum;
//      cossum += stupidarray[k] * cos(sumarg);
//      sinsum += stupidarray[k] * sin(sumarg);
//    }
//    pme.pmetemp.moduli3[i] = sinsum*sinsum + cossum*cossum;
//  }
//  // Coorection for the euler interpolation
//  if (pme.pmetemp.moduli3[0] < eps) pme.pmetemp.moduli3[0] = 0.5 * pme.pmetemp.moduli3[1];
//  for (int i = 1; i < pme.pmetemp.nzpoints - 1; i++)
//  {
//    if (pme.pmetemp.moduli3[i] < eps) pme.pmetemp.moduli3[i] = 0.5 * (pme.pmetemp.moduli3[i - 1] + pme.pmetemp.moduli3[i + 1]);
//  }
//  if (pme.pmetemp.moduli3[pme.pmetemp.nzpoints] < eps)  pme.pmetemp.moduli3[pme.pmetemp.nzpoints] = 0.5 * pme.pmetemp.moduli3[pme.pmetemp.nzpoints - 1];
//  // calculate some stupid factor someone called zeta...whatever this thing does....
//  acutoff = 50;
//  doubleorder = 2 * 5;
//  for (int i = 0; i < pme.pmetemp.nzpoints; i++)
//  {
//    indexk = i;
//    if (i > pme.pmetemp.nzpoints / 2) indexk -= pme.pmetemp.nzpoints;
//    if (indexk == 0) zeta = 1.0;
//    else
//    {
//      tempsum1 = tempsum2 = 1.0;
//      fact = SUPERPI * double(indexk) / double(pme.pmetemp.nzpoints);
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact + SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      for (int j = 1; j <= acutoff; j++)
//      {
//        sumarg = fact / (fact - SUPERPI * double(j));
//        tempsum1 += pow(sumarg, 5.0);
//        tempsum2 += pow(sumarg, double(doubleorder));
//      }
//      zeta = tempsum2 / tempsum1;
//    }
//    pme.pmetemp.moduli3[i] *= (zeta*zeta);
//  }
//}
//
//#ifdef _OPENMP
//// calculate spatial sites in the box for partial parallelization
//void coords::Coordinates::roughgrid()
//{
//
//  pme.pmetemp.nrough1 = pme.pmetemp.nrough2 = pme.pmetemp.nrough3 = pme.pmetemp.rgridtotal = 1;
//  int threads;
//#pragma omp parallel
//  {
//    threads = omp_get_num_threads();
//  }
//  for (int i = 2; i <= 6; i++)
//  {
//    if (threads > pme.pmetemp.rgridtotal && (pme.pmetemp.nxpoints%i == 0))
//    {
//      pme.pmetemp.nrough1 = i;
//      pme.pmetemp.rgridtotal = pme.pmetemp.nrough1 * pme.pmetemp.nrough2 * pme.pmetemp.nrough3;
//    }
//    if (threads > pme.pmetemp.rgridtotal && (pme.pmetemp.nypoints%i == 0))
//    {
//      pme.pmetemp.nrough2 = i;
//      pme.pmetemp.rgridtotal = pme.pmetemp.nrough1 * pme.pmetemp.nrough2 * pme.pmetemp.nrough3;
//    }
//    if (threads > pme.pmetemp.rgridtotal && (pme.pmetemp.nzpoints%i == 0))
//    {
//      pme.pmetemp.nrough3 = i;
//      pme.pmetemp.rgridtotal = pme.pmetemp.nrough1 * pme.pmetemp.nrough2 * pme.pmetemp.nrough3;
//    }
//  }
//  // x,y,z-axis number of points per chunk
//  pme.pmetemp.rgrid1 = pme.pmetemp.nxpoints / pme.pmetemp.nrough1;
//  pme.pmetemp.rgrid2 = pme.pmetemp.nypoints / pme.pmetemp.nrough2;
//  pme.pmetemp.rgrid3 = pme.pmetemp.nzpoints / pme.pmetemp.nrough3;
//  // offset for bsplines and left and right points of central point
//  pme.pmetemp.roughleft = (pme.pmetemp.bsplineorder - 1) / 2;
//  pme.pmetemp.roughright = pme.pmetemp.bsplineorder - pme.pmetemp.roughleft - 1;
//  pme.pmetemp.bsoffset = (pme.pmetemp.bsplineorder + 1) / 2 + 1;
//
//}
//#endif

