#include <cmath>
#include <stdexcept>
#include "atomic.h"
#include "coords.h"
#include "configuration.h"
#include "coords_io.h"
#include "lbfgs.h"

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
  m_atoms.clear();
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
    if (m_preinterface) 
      m_preinterface->update(false);
  }
}

void coords::Coordinates::init_in(Atoms a, PES_Point p, bool const update)
{
  init_swap_in(a, p, update);
}

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
  m_representation.energy = optimizer.p().f;
  m_representation.structure.cartesian =
    coords::Representation_3D(optimizer.p().x.begin(), optimizer.p().x.end());
  m_representation.gradient.cartesian = coords::Gradients_3D(optimizer.p().g.begin(), optimizer.p().g.end());
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

void coords::Coordinates::h_tostream(std::ostream &S,
  energy::interface_base const * const ep) const
{
  std::vector<std::vector<double>> hess = m_representation.hessian;
  if (hess.size() == 0)
  {
    std::cout<<"ERROR: no hessian can be printed!\n";
  }
  else
  {
    S << "HESSIAN MATRIX";
    S << "\n\n            ";
    for (unsigned i=0; i<size(); i++)  // headline
    {
      S << "| X ("<< std::right << std::fixed << std::setw(6)<<i+1 <<") | Y ("<< std::right << std::fixed << std::setw(6)<<i+1 <<") | Z (" << std::right << std::fixed << std::setw(6) <<i+1<< ") ";
    }
    S<<"\n";
    for (unsigned i=0; i<(13*3*m_representation.size()+13); i++)  // second line
    {
      S << "-";
    }
    S<<"\n";
    for (unsigned i=0; i<m_representation.size()*3; i++)  // lines
    {
      if (i%3 == 0)
      {
        S << " X (" << std::right << std::fixed << std::setw(6) << i/3 + 1 << ") ";
      }
      else if (i%3 == 1)
      {
        S << " Y (" << std::right << std::fixed << std::setw(6) << i/3 + 1 << ") ";
      }
      else
      {
        S << " Z (" << std::right << std::fixed << std::setw(6) << i/3 + 1 << ") ";
      }
      for (unsigned j=0; j<3*m_representation.size();j++)  // columns
      {
        S <<"|"<<std::right<<std::fixed<<std::setw(11)<<std::setprecision(2)<<hess[i][j]<<" ";
      }
      S << "\n";
    }
  }
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
    coords::float_type mass(m_atoms.atom(m_atoms.atomOfMolecule(index, i)).mass());
    M += mass;
    COM += xyz(m_atoms.atomOfMolecule(index, i))*mass;
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
	for (size_type i(0u); i < N; ++i)
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
	Cartesian_Point const halfbox(Config::get().periodics.pb_box / 2.0);
	for (std::size_t i = 0; i < N; ++i)
	{
		Cartesian_Point tmp_com(-center_of_mass_mol(i));
		if (std::abs(tmp_com.x()) <= halfbox.x())
		{
			tmp_com.x() = 0;
		}
		else
		{
			tmp_com.x() = tmp_com.x() / Config::get().periodics.pb_box.x();
			int tmp_x = std::round(tmp_com.x());
			tmp_com.x() = tmp_x;
		}
		if (std::abs(tmp_com.y()) <= halfbox.y())
		{
			tmp_com.y() = 0;
		}
		else
		{
			tmp_com.y() = tmp_com.y() / Config::get().periodics.pb_box.y();
			int tmp_y = std::round(tmp_com.y());
			tmp_com.y() = tmp_y;
		}
		if (std::abs(tmp_com.z()) <= halfbox.z())
		{
			tmp_com.z() = 0;
		}
		else
		{
			tmp_com.z() = tmp_com.z() / Config::get().periodics.pb_box.z();
			int tmp_z = std::round(tmp_com.z());
			tmp_com.z() = tmp_z;
		}
		tmp_com *= Config::get().periodics.pb_box;
    for (auto const atom : molecule(i)) move_atom_by(atom, tmp_com, true);
  }
}

bool coords::Coordinates::validate_bonds()
{
  bool status = true;
  broken_bonds.clear();
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

bool coords::Coordinates::is_equal_structure(coords::PES_Point const & a, coords::PES_Point const & b) const
{
  using scon::operator-;
  config::coords::conditionsForStructuresToBeConsideredEqual equalityConditions = Config::get().coords.equals;
  coords::main_type const& md = equalityConditions.main;
  coords::internal_type const& id = equalityConditions.intern;
  coords::Cartesian_Point const& cd = equalityConditions.xyz;
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
    if (equal) 
      return true;
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
    if (equal) 
      return true;
  }
  // Check cartesian structure equality
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
    if (equal) 
      return true;
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
void coords::Coordinates::adapt_indexation(size_t no_dist, size_t no_angle, size_t no_dihedral,
  std::vector<std::vector<std::pair<std::vector<size_t>, double>>> const &reference,
  coords::Coordinates const *cPtr)
{
  size_t ibond = 0, iangle = 0, idihedral = 0;
  size_t N = atoms().size();
  for (size_t i = 0; i < N; ++i)
  {
    std::cout << atoms(i).a_to_i();
    atoms_changeable(i).set_i_of_a(i);
    std::cout << " changed to " << atoms(i).a_to_i() << " which is " << i << '\n';
  }
  for (size_t i = 0; i < N; ++i)
  {

    if (i == 0)
    {
      ibond = N;
      iangle = N + 1;
      idihedral = N + 2;
    }

    else if (i == 1)
    {
      ibond = 0;
      iangle = N;
      idihedral = N + 1;
    }

    else if (i == 2)
    {
      ibond = 1;
      iangle = 0;
      idihedral = N;
    }
    //shouldn't be needed anymore
    /*else if (i == no_dist)
    {
    if (reference[no_dist][0].first[0] == 1)
    {
    ibond = atoms(reference[no_dist][0].first[1] - 1).a_to_i();
    }
    else
    {
    ibond = atoms(reference[no_dist][0].first[0] - 1).a_to_i();
    }

    if (reference[no_dist][1].first[0] == 1)
    {
    iangle = atoms(reference[no_dist][1].first[2] - 1).a_to_i();
    }
    else
    {
    iangle = atoms(reference[no_dist][1].first[0] - 1).a_to_i();
    }

    if (reference[no_dist][2].first[0] == 1)
    {
    idihedral = atoms(reference[no_dist][2].first[3] - 1).a_to_i();
    }
    else
    {
    idihedral = atoms(reference[no_dist][2].first[0] - 1).a_to_i();
    }
    }
    else if (i == no_angle && i > 2)
    {
    if (reference[no_angle][0].first[0] == 2)
    {
    ibond = atoms(reference[no_angle][0].first[1] - 1).a_to_i();
    }
    else
    {
    ibond = atoms(reference[no_angle][0].first[0] - 1).a_to_i();
    }

    if (reference[no_angle][1].first[0] == 2)
    {
    iangle = atoms(reference[no_angle][1].first[2] - 1).a_to_i();
    }
    else
    {
    iangle = atoms(reference[no_angle][1].first[0] - 1).a_to_i();
    }

    if (reference[no_angle][2].first[0] == 2)
    {
    idihedral = atoms(reference[no_angle][2].first[3] - 1).a_to_i();
    }
    else
    {
    idihedral = atoms(reference[no_angle][2].first[0] - 1).a_to_i();
    }
    }
    else if (i == no_dihedral && i > 2)
    {
    if (reference[no_dihedral][0].first[0] == 3)
    {
    ibond = atoms(reference[no_dihedral][0].first[1] - 1).a_to_i();
    }
    else
    {
    ibond = atoms(reference[no_dihedral][0].first[0] - 1).a_to_i();
    }

    if (reference[no_dihedral][1].first[0] == 3)
    {
    iangle = atoms(reference[no_dihedral][1].first[2] - 1).a_to_i();
    }
    else
    {
    iangle = atoms(reference[no_dihedral][1].first[0] - 1).a_to_i();
    }

    if (reference[no_dihedral][2].first[0] == 3)
    {
    idihedral = atoms(reference[no_dihedral][2].first[3] - 1).a_to_i();
    }
    else
    {
    idihedral = atoms(reference[no_dihedral][2].first[0] - 1).a_to_i();
    }
    }*/
    else
    {
      if (reference[i][0].first[0] == i + 1)
      {
        ibond = atoms(reference[i][0].first[1] - 1).a_to_i();
      }
      else
      {
        ibond = atoms(reference[i][0].first[0] - 1).a_to_i();
      }

      if (reference[i][1].first[0] == i + 1)
      {
        iangle = atoms(reference[i][1].first[2] - 1).a_to_i();
      }
      else
      {
        iangle = atoms(reference[i][1].first[0] - 1).a_to_i();
      }

      if (reference[i][2].first[0] == i + 1)
      {
        idihedral = atoms(reference[i][2].first[3] - 1).a_to_i();
      }
      else
      {
        idihedral = atoms(reference[i][2].first[0] - 1).a_to_i();
      }
    }
    atoms_changeable(i).set_ibond(ibond);
    atoms_changeable(i).set_iangle(iangle);
    atoms_changeable(i).set_idihedral(idihedral);

    size_t temp = atoms(i).bonds().size();
    for (size_t j = 0; j < temp; ++j)
    {
      atoms_changeable(i).detach_from(atoms(i).bonds(0));
    }
    for (size_t j = 0; j < cPtr->atoms(i).bonds().size(); ++j)
    {
      atoms_changeable(i).bind_to(cPtr->atoms(i).bonds(j));
    }
  }
  for (size_t i = 0; i < N; ++i)
  {
    for (auto &a : atoms(i).bonds())
    {
      std::cout << a + 1 << ' ';
    }
    std::cout << '\n';
  }
}

std::vector<bool> const coords::Coordinates::terminal()
{
  size_t N = this->atoms().size();
  std::vector<bool> terminal;
  terminal.resize(N);
  for (size_t i = 0; i < N; ++i)
  {
    if (this->atoms(i).bonds().size() == 1 && this->atoms(i).bonds().size() != 0)
    {
      terminal[i] = true;
    }
    else
    {
      terminal[i] = false;
    }
  }
  return terminal;
}

//*gets terminal atoms and marks them with 1, ignores those atoms and gets atoms that become terminal in resulting replic, marks those with 2 and procedes analogously
//*through the rest of the structure, until either all atoms have been marked, or a core of (linked) rings is the result
std::vector<size_t> const coords::Coordinates::terminal_enum()
{
  std::vector<size_t> terminal_enum;
  std::vector<bool> terminal;
  coords::Coordinates red_replic = *this;
  coords::Coordinates red_replic_temp;
  size_t N = red_replic.atoms().size();
  size_t enum_index = 1;
  bool done;
  terminal_enum.resize(N);
  terminal = red_replic.terminal();
  for (size_t i = 0; i < N; ++i)
  {
    terminal_enum[i] = 0;
  }
  while (true)
  {
    red_replic_temp = red_replic;
    red_replic = red_replic.get_red_replic(terminal);
    for (size_t i = 0; i < N; ++i)
    {
      if (terminal[i])
      {
        terminal_enum[i] = enum_index;
      }
    }
    done = red_replic.atoms() == red_replic_temp.atoms();
    if (!done)
    {
      for (size_t i = 0; i < N; ++i)
      {
        for (size_t j = 0; j < N; ++j)
        {
          if (terminal[j] && red_replic.atoms(i).is_bound_to(j))
          {
            red_replic.atoms_changeable(i).detach_from(j);
          }
        }
      }
      terminal = red_replic.terminal();
      ++enum_index;
    }
    else
    {
      break;
    }
  }
  for (size_t i = 0; i < N; ++i)
  {
    if (terminal_enum[i] == 0)								//*terminal_enum-value for core atoms, as in rings and ring linking atoms, which will never become terminal, 
    {																					//*no matter how many times the replic is reduced. Core atoms are given the highest value, so as to be able
      terminal_enum[i] = enum_index;					//*to identify them more intuitively
    }
  }
  return terminal_enum;
}

coords::Coordinates coords::Coordinates::get_red_replic(std::vector<bool> criterion)
{
  coords::Coordinates replic = *this;
  size_t N = replic.atoms().size();
  coords::Coordinates red_replic;
  //red_replic.atoms().resize(N);
  for (size_t i = 0; i < N; ++i)
  {
    if (!criterion[i])
    {
      red_replic.atoms_changable().add(replic.atoms(i));
    }
    else
    {
      red_replic.atoms_changable().add(coords::Atom());
    }
  }
  return red_replic;
}