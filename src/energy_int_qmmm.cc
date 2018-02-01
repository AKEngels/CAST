#include <cstddef>

#include "energy_int_qmmm.h"
#include "scon_utility.h"

::tinker::parameter::parameters energy::interfaces::qmmm::QMMM::tp;

namespace
{

  std::vector<std::size_t> get_mm_atoms(std::size_t const num_atoms)
  {
    std::vector<std::size_t> mm_atoms;
    auto qm_size = Config::get().energy.qmmm.qmatoms.size();
    if (num_atoms > qm_size)
    {
      auto mm_size = num_atoms - qm_size;
      mm_atoms.reserve(mm_size);
      for (unsigned i = 0; i < num_atoms; i++)
      {
        if (!scon::sorted::exists(Config::get().energy.qmmm.qmatoms, i)) 
        { 
          mm_atoms.emplace_back(i); 
        }
      }
    }
    return mm_atoms;
  }

  std::vector<std::size_t> make_new_indices_qm(std::size_t const num_atoms)
  {
    std::vector<std::size_t> new_indices_qm;
    if (num_atoms > 0)
    {
      new_indices_qm.resize(num_atoms);
      std::size_t current_index = 0u;
      for (auto&& a : Config::get().energy.qmmm.qmatoms)
      {
        new_indices_qm.at(a) = current_index++;
      }
    }
    return new_indices_qm;
  }

  std::vector<std::size_t> make_new_indices_mm(std::size_t const num_atoms, 
    std::vector<std::size_t> const& mmi)
  {
    std::vector<std::size_t> new_indices_mm;
    new_indices_mm.resize(num_atoms);
    if (num_atoms > 0u)
    {
      std::size_t current_index = 0u;

      for (auto && a : mmi)
      {
        new_indices_mm.at(a) = current_index++;
      }
    }
    return new_indices_mm;
  }

  coords::Coordinates make_qm_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices)
  {
    auto tmp_i = Config::get().general.energy_interface;
    Config::set().general.energy_interface = Config::get().energy.qmmm.qminterface;
    coords::Coordinates new_qm_coords;
    if (cp->size() >= indices.size())
    {
      coords::Atoms new_qm_atoms;
      coords::PES_Point pes;
      pes.structure.cartesian.reserve(indices.size());
      for (auto && a : indices)
      {
        auto && ref_at = (*cp).atoms().atom(a);
        coords::Atom at{ (*cp).atoms().atom(a).number() };
        at.set_energy_type(ref_at.energy_type());
        auto bonds = ref_at.bonds();
        for (auto && b : bonds)
        {
          at.detach_from(b);
        }
        for (auto && b : bonds)
        {
          at.bind_to(new_indices.at(b));
        }
        new_qm_atoms.add(at);
        pes.structure.cartesian.push_back(cp->xyz(a));
      }
      new_qm_coords.init_swap_in(new_qm_atoms, pes);
    }
    Config::set().general.energy_interface = tmp_i;
    return new_qm_coords;
  }

  /**creates coordobject for MM interface
  @param cp: coordobj for whole system (QM + MM)
  @param indices: indizes of MM atoms
  @param new_indices: new indizes???*/
  coords::Coordinates make_aco_coords(coords::Coordinates const * cp,
    std::vector<std::size_t> const & indices, std::vector<std::size_t> const & new_indices)
  {
    auto tmp_i = Config::get().general.energy_interface;
    Config::set().general.energy_interface = Config::get().energy.qmmm.mminterface;
    coords::Coordinates new_aco_coords;
    if (cp->size() >= indices.size())
    {
      coords::Atoms new_aco_atoms;
      coords::PES_Point pes;
      pes.structure.cartesian.reserve(indices.size());
      for (auto && a : indices)
      {
        auto && ref_at = (*cp).atoms().atom(a);
        coords::Atom at{ (*cp).atoms().atom(a).number() };
        at.set_energy_type(ref_at.energy_type());
        auto bonds = ref_at.bonds();
        for (auto && b : bonds)
        {
          at.detach_from(b);
        }
        for (auto && b : bonds)
        {
          if (is_in(b,indices)) at.bind_to(new_indices[b]);  // only bind if bonding partner is also in ACO coords
        }
        new_aco_atoms.add(at);
        pes.structure.cartesian.push_back(cp->xyz(a));
      }
      new_aco_coords.init_swap_in(new_aco_atoms, pes);
    }
    Config::set().general.energy_interface = tmp_i;
    return new_aco_coords;
  }
}

energy::interfaces::qmmm::QMMM::QMMM(coords::Coordinates * cp) :
  interface_base(cp),
  qm_indices(Config::get().energy.qmmm.qmatoms),
  mm_indices(get_mm_atoms(cp->size())),
  new_indices_qm(make_new_indices_qm(cp->size())),
  new_indices_mm(make_new_indices_mm(cp->size(), mm_indices)),
  qmc(make_qm_coords(cp, qm_indices, new_indices_qm)),
  mmc(make_aco_coords(cp, mm_indices, new_indices_mm)),
  qm_energy(0.0), mm_energy(0.0), vdw_energy(0.0), bonded_energy(0.0)
{
  Config::set().energy.qmmm.mm_atoms_number = cp->size() - qm_indices.size();
  if (!tp.valid())
  {
    tp.from_file(Config::get().get().general.paramFilename);
  }
  std::vector<std::size_t> types;
  for (auto atom : (*cp).atoms())
  {
    scon::sorted::insert_unique(types, atom.energy_type());
  }
  cparams = tp.contract(types);
  torsionunit = cparams.torsionunit();
  prepare_bonded_qmmm();
}

energy::interfaces::qmmm::QMMM::QMMM(QMMM const & rhs, 
  coords::Coordinates *cobj) : interface_base(cobj), 
  cparams(rhs.cparams), distance(rhs.distance), 
  qm_indices(rhs.qm_indices), mm_indices(rhs.mm_indices),
  new_indices_qm(rhs.new_indices_qm), new_indices_mm(rhs.new_indices_mm),
  qmc(rhs.qmc), mmc(rhs.mmc), qm_charge_vector(rhs.qm_charge_vector), 
  mm_charge_vector(rhs.mm_charge_vector), vdw_energy(rhs.vdw_energy),  
  qm_energy(rhs.qm_energy), mm_energy(rhs.mm_energy), vdw_gradient(rhs.vdw_gradient),
  c_gradient(rhs.c_gradient), bonded_energy(rhs.bonded_energy), bonded_gradient(rhs.bonded_gradient)
{
  interface_base::operator=(rhs);
}

energy::interfaces::qmmm::QMMM::QMMM(QMMM&& rhs, coords::Coordinates *cobj) 
  : interface_base(cobj),
  cparams(std::move(rhs.cparams)), distance(std::move(rhs.distance)),
  qm_indices(std::move(rhs.qm_indices)), mm_indices(std::move(rhs.mm_indices)),
  new_indices_qm(std::move(rhs.new_indices_qm)), 
  new_indices_mm(std::move(rhs.new_indices_mm)),
  qmc(std::move(rhs.qmc)), mmc(std::move(rhs.mmc)), 
  qm_charge_vector(std::move(rhs.qm_charge_vector)),
  mm_charge_vector(std::move(rhs.mm_charge_vector)),
  vdw_energy(std::move(rhs.vdw_energy)),
  qm_energy(std::move(rhs.qm_energy)), mm_energy(std::move(rhs.mm_energy)),
  c_gradient(std::move(rhs.c_gradient)), vdw_gradient(std::move(rhs.vdw_gradient)), 
  bonded_energy(std::move(rhs.bonded_energy)), bonded_gradient(std::move(rhs.bonded_gradient))
{
  interface_base::operator=(rhs);
}

void energy::interfaces::qmmm::QMMM::find_parameters()
{
  for (auto &b : qmmm_bonds)  // find bond parameters
  {
    auto b_type_a = cparams.type(coords->atoms().atom(b.a).energy_type(), tinker::potential_keys::BOND);
    auto b_type_b = cparams.type(coords->atoms().atom(b.b).energy_type(), tinker::potential_keys::BOND);
    for (auto b_param : cparams.bonds())
    {
      if (b_param.index[0] == b_type_a && b_param.index[1] == b_type_b)
      {
        b.ideal = b_param.ideal;
        b.force = b_param.f;
      }
      else if (b_param.index[0] == b_type_b && b_param.index[1] == b_type_a)
      {
        b.ideal = b_param.ideal;
        b.force = b_param.f;
      }
    }
  }

  for (auto &a : qmmm_angles)  // find angle parameters
  {
    auto a_type_a = cparams.type(coords->atoms().atom(a.a).energy_type(), tinker::potential_keys::ANGLE);
    auto a_type_b = cparams.type(coords->atoms().atom(a.b).energy_type(), tinker::potential_keys::ANGLE);
    auto a_type_c = cparams.type(coords->atoms().atom(a.c).energy_type(), tinker::potential_keys::ANGLE);
    for (auto a_param : cparams.angles())
    {
      if (a_param.index[1] == a_type_c)
      {
        if (a_param.index[0] == a_type_a && a_param.index[2] == a_type_b)
        {
          a.ideal = a_param.ideal;
          a.force = a_param.f;
        }
        else if (a_param.index[0] == a_type_b && a_param.index[2] == a_type_b)
        {
          a.ideal = a_param.ideal;
          a.force = a_param.f;
        }
      }
    }
  }

  for (auto &d : qmmm_dihedrals)  // find parameters for dihedrals
  {
    auto d_type_a = cparams.type(coords->atoms().atom(d.a).energy_type(), tinker::potential_keys::TORSION);
    auto d_type_b = cparams.type(coords->atoms().atom(d.b).energy_type(), tinker::potential_keys::TORSION);
    auto d_type_c1 = cparams.type(coords->atoms().atom(d.c1).energy_type(), tinker::potential_keys::TORSION);
    auto d_type_c2 = cparams.type(coords->atoms().atom(d.c2).energy_type(), tinker::potential_keys::TORSION);
    for (auto d_param : cparams.torsions())
    {
      if (d_param.index[0] == d_type_a && d_param.index[1] == d_type_c1 && d_param.index[2] == d_type_c2 && d_param.index[3] == d_type_b)
      {
        d.max_order = d_param.max_order;
        d.number = d_param.number;
        d.orders = d_param.order;
        d.forces = d_param.force;
        d.ideals = d_param.ideal;
      }
      else if (d_param.index[0] == d_type_b && d_param.index[1] == d_type_c2 && d_param.index[2] == d_type_c1 && d_param.index[3] == d_type_a)
      {
        d.max_order = d_param.max_order;
        d.number = d_param.number;
        d.orders = d_param.order;
        d.forces = d_param.force;
        d.ideals = d_param.ideal;
      }
    }
  }
}

// calculate position of a link atom (see: doi 10.1002/jcc.20857)
coords::cartesian_type energy::interfaces::qmmm::QMMM::calc_position(bonded::LinkAtom link)
{
  coords::cartesian_type r_MM = coords->xyz(link.mm);
  coords::cartesian_type r_QM = coords->xyz(link.qm);
  double d_MM_QM = dist(r_MM, r_QM);
  std::cout << "MM atom: " << r_MM << " ; QM atom: " << r_QM << "\n";

  coords::cartesian_type pos = r_QM + ((r_MM - r_QM) / d_MM_QM) * link.d_L_QM;
  std::cout << "Pos: " << pos << "\n";
  return pos;
}

// creates link atom for every QM/MM bond
void energy::interfaces::qmmm::QMMM::create_link_atoms()
{
  for (auto b : qmmm_bonds)
  {
    //create link atom
    bonded::LinkAtom link;
    link.qm = b.b;
    link.mm = b.a;

    // determine equilibrium distance between link atom and QM atom from force field
    auto b_type_qm = cparams.type(coords->atoms().atom(b.b).energy_type(), tinker::potential_keys::BOND);
    auto b_type_L = cparams.type(85, tinker::potential_keys::BOND);
    for (auto b_param : cparams.bonds())
    {
      if (b_param.index[0] == b_type_qm && b_param.index[1] == b_type_L)  link.d_L_QM = b_param.ideal;
      else if (b_param.index[0] == b_type_L && b_param.index[1] == b_type_qm) link.d_L_QM = b_param.ideal;
    }
    if (link.d_L_QM == 0.0)  throw std::runtime_error("Determining position of link atom is not possible.\n");

    // calculate position of link atom
    link.position = calc_position(link);
    
    // add link atom to vector
    link_atoms.push_back(link);
  }
}

void energy::interfaces::qmmm::QMMM::find_bonds_etc()
{
  // find bonds between QM and MM region
  for (auto mma : mm_indices)
  {
    for (auto b : coords->atoms().atom(mma).bonds())
    {
      if (scon::sorted::exists(qm_indices, b))
      {
        if (Config::get().energy.qmmm.mminterface == config::interface_types::T::OPLSAA)
        {
          bonded::Bond bond(mma, b);
          qmmm_bonds.push_back(bond);
        }
        else
        {
          throw std::runtime_error("Breaking bonds is only possible with OPLSAA as MM interface.\n");
        }
      }
    }
  }

  // find angles between QM and MM region
  for (auto b : qmmm_bonds) // for every bond
  {
    for (auto p : coords->atoms().atom(b.a).bonds()) // atom a as angle center
    {
      if (b.b != p)
      {
        bonded::Angle angle(b.b, p, b.a);
        if (!bonded::is_in(angle, qmmm_angles))
        {
          qmmm_angles.push_back(angle);
        }
      }
    }
    for (auto p : coords->atoms().atom(b.b).bonds()) // atom b as angle center
    {
      if (b.a != p)
      {
        bonded::Angle angle(b.a, p, b.b);
        if (!bonded::is_in(angle, qmmm_angles))
        {
          qmmm_angles.push_back(angle);
        }
      }
    }
  }

  // find dihedrals between QM and MM region
  for (auto a : qmmm_angles)
  {
    for (auto p : coords->atoms().atom(a.a).bonds())   // expand angle at atom a
    {
      if (a.c != p) 
      {
        bonded::Dihedral dihed(p, a.b, a.a, a.c);
        if (!bonded::is_in(dihed, qmmm_dihedrals))
        {
          qmmm_dihedrals.push_back(dihed);
        }
      }
    }
    for (auto p : coords->atoms().atom(a.b).bonds())   // expand angle at atom b
    {
      if (a.c != p)
      {
        bonded::Dihedral dihed(p, a.a, a.b, a.c);
        if (!bonded::is_in(dihed, qmmm_dihedrals))
        {
          qmmm_dihedrals.push_back(dihed);
        }
      }
    }
  }
}

void energy::interfaces::qmmm::QMMM::prepare_bonded_qmmm()
{
  find_bonds_etc();    // find bonds, angles and dihedrals between QM and MM region
  find_parameters();   // find force field parameters for energy calculation
  create_link_atoms(); // create link atoms

  if (Config::get().general.verbosity > 4)  // Output
  {
    std::cout << "QM/MM-Bonds\n";
    for (auto b : qmmm_bonds)
    {
      std::cout << b.info() << "\n";
    }
    std::cout << "QM/MM-Angles\n";
    for (auto a : qmmm_angles)
    {
      std::cout << a.info() << "\n";
    }
    std::cout << "QM/MM-Dihedrals\n";
    for (auto a : qmmm_dihedrals)
    {
      std::cout << a.info() << "\n";
    }
  }
}

void energy::interfaces::qmmm::QMMM::update_representation()
{
  std::size_t qi = 0u;
  for (auto i : qm_indices)
  {
    qmc.move_atom_to(qi, coords->xyz()[i], true);
    ++qi;
  }
  std::size_t mi = 0u;
  for (auto j : mm_indices)
  {
    mmc.move_atom_to(mi, coords->xyz()[j], true);
    ++mi;
  }

}

// mol.in schreiben (for MOPAC, see http://openmopac.net/manual/QMMM.html)
void energy::interfaces::qmmm::QMMM::write_mol_in()
{
  auto elec_factor = 332.0;
  std::cout << std::setprecision(6);

  std::ofstream molstream{ "mol.in" };
  if (molstream)
  {
    auto const n_qm = qm_indices.size();
    molstream << '\n';
    molstream << n_qm << " 0\n";
    for (std::size_t i = 0; i < n_qm; ++i)
    {
      double qi{};
      for (std::size_t j = 0; j < mm_charge_vector.size(); ++j)
      {
        auto d = len(coords->xyz(qm_indices[i]) - coords->xyz(mm_indices[j]));
        qi += mm_charge_vector[j] / d;
      }
      qi *= elec_factor;
      molstream << "0 0 0 0 " << qi << "\n";
    }
    molstream.close();
  }
  else
  {
    throw std::runtime_error("Cannot write mol.in file.");
  }
}

// write gaussian inputfile
void energy::interfaces::qmmm::QMMM::write_gaussian_in(char calc_type)
{
  std::string id = qmc.energyinterface()->get_id();
  id += "_G_";
  std::string outstring(id);
  outstring.append(".gjf");

  std::ofstream out_file(outstring.c_str(), std::ios_base::out);

  if (out_file)
  {
    if (Config::get().energy.gaussian.link.length() != 0) { //if no link commands are issued the line wil be skipped

      std::istringstream iss(Config::get().energy.gaussian.link);

      std::vector<std::string> splitted_str(
        std::istream_iterator<std::string>{iss},
        std::istream_iterator<std::string>{}
      );

      for (std::size_t i = 0; i < splitted_str.size(); i++)
      {
        out_file << '%' << splitted_str[i] << '\n';
      }

    }
    out_file << "# " << Config::get().energy.gaussian.method << " " << Config::get().energy.gaussian.basisset << " " << Config::get().energy.gaussian.spec << " " << "Charge NoSymm ";

    switch (calc_type) {// to ensure the needed gaussian keywords are used in gausian inputfile for the specified calculation
    case 'g':
      out_file << " Force Prop=(Field,Read) Density";
      break;
    }

    out_file << '\n';
    out_file << '\n';
    out_file << Config::get().general.outputFilename;
    out_file << '\n';
    out_file << '\n';
    out_file << Config::get().energy.gaussian.charge << " ";
    out_file << Config::get().energy.gaussian.multipl;
    out_file << '\n';
    out_file << coords::output::formats::xyz(qmc);
    out_file << '\n';
    for (std::size_t j = 0; j < mm_charge_vector.size(); ++j)  // writing additional point charges (from MM atoms)
    {
      out_file << coords->xyz(mm_indices[j]).x() << " " << coords->xyz(mm_indices[j]).y() << " " << coords->xyz(mm_indices[j]).z() << " " << mm_charge_vector[j] << "\n";
    }
    out_file << '\n';
    if (Config::get().energy.gaussian.method == "DFTB=read")
    {
      std::vector<std::vector<std::string>> pairs = find_pairs(*coords);
      for (auto p : pairs)
      {
        std::string filename = p[0] + "-" + p[1] + ".skf";
        if (file_exists(filename) == false)
        {
          std::cout << "ERROR! Slater Koster file " << filename << " does not exist. Please download it from dftb.org and convert it with the task MODIFY_SK_FILES!\n";
          std::exit(0);
        }
        out_file << "@./" << filename << " /N\n";
      }
    }
    else if (Config::get().energy.gaussian.method == "DFTBA")
    {
      out_file << "@GAUSS_EXEDIR:dftba.prm\n";
    }
    if (calc_type == 'g')
    {
      out_file << '\n';
      for (std::size_t j = 0; j < mm_charge_vector.size(); ++j)  // writing points for electric field (positions of MM atoms)
      {
        out_file << coords->xyz(mm_indices[j]).x() << " " << coords->xyz(mm_indices[j]).y() << " " << coords->xyz(mm_indices[j]).z() << "\n";
      }
    }
    out_file.close();
  }
  else std::runtime_error("Writing Gaussian Inputfile failed.");
}

// write inputfiles for DFTB+ (real inputfile and file with MM charges)
void energy::interfaces::qmmm::QMMM::write_dftb_in(char calc_type)
{
  std::ofstream chargefile("charges.dat");
  for (int j = 0; j < mm_charge_vector.size(); j++)
  {
    chargefile << coords->xyz(mm_indices[j]).x() << " " << coords->xyz(mm_indices[j]).y() << " " << coords->xyz(mm_indices[j]).z() << "  " << mm_charge_vector[j] << "\n";
  }
  chargefile.close();

  // create a vector with all element symbols that are found in the structure
  // are needed for writing angular momenta into inputfile
  std::vector<std::string> elements;
  for (auto a : qmc.atoms())
  {
    if (is_in(a.symbol(), elements) == false)
    {
      elements.push_back(a.symbol());
    }
  }

  // create inputfile
  std::ofstream file("dftb_in.hsd");

  // write geometry
  file << "Geometry = GenFormat {\n";
  file << qmc.size()+link_atoms.size() << "  C\n";  // no supercells possible
  for (auto s : elements)
  {
    file << s << " ";
  }
  file << "\n";
  for (std::size_t i(0U); i < qmc.size(); ++i)  // atoms of QM coordobject
  {
    file << std::left << std::setw(5) << i + 1 << std::left << std::setw(5) << find_index(qmc.atoms(i).symbol(), elements) + 1;
    file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << qmc.xyz(i).x();
    file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << qmc.xyz(i).y();
    file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << qmc.xyz(i).z();
    file << '\n';
  }
  if (link_atoms.size() > 0)  // if link atoms: add them to inputfile
  {
    for (std::size_t i(0U); i < link_atoms.size(); ++i)
    {
      file << std::left << std::setw(5) << qmc.size() + i + 1 << std::left << std::setw(5) << find_index("H", elements) + 1;
      file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << link_atoms[i].position.x();
      file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << link_atoms[i].position.y();
      file << std::fixed << std::showpoint << std::right << std::setw(12) << std::setprecision(6) << link_atoms[i].position.z();
      file << '\n';
    }
  }
  file << "}\n\n";

  // write information that is needed for SCC calculation
  file << "Hamiltonian = DFTB {\n";
  file << "  SCC = Yes\n";
  file << "  SCCTolerance = " << std::scientific << Config::get().energy.dftb.scctol << "\n";
  file << "  MaxSCCIterations = " << Config::get().energy.dftb.max_steps << "\n";
  file << "  Charge = " << Config::get().energy.dftb.charge << "\n";
  file << "  SlaterKosterFiles = Type2FileNames {\n";
  file << "    Prefix = '" << Config::get().energy.dftb.sk_files << "'\n";
  file << "    Separator = '-'\n";
  file << "    Suffix = '.skf'\n";
  file << "  }\n";

  file << "  ElectricField = {\n";
  file << "    PointCharges = {\n";
  file << "      CoordsAndCharges [Angstrom] = DirectRead {\n";
  file << "        Records = " << Config::get().energy.qmmm.mm_atoms_number << "\n";
  file << "        File = 'charges.dat'\n";
  file << "      }\n";
  file << "    }\n";
  file << "  }\n";

  file << "  MaxAngularMomentum {\n";
  for (auto s : elements)
  {
    char angular_momentum = atomic::angular_momentum_by_symbol(s);
    if (angular_momentum == 'e')
    {
      std::cout << "Angular momentum for element " << s << " not defined. \n";
      std::cout << "Please go to file 'atomic.h' and define an angular momentum (s, p, d or f) in the array 'angular_momentum'.\n";
      std::cout << "Talk to a CAST developer if this is not possible for you.";
      std::exit(0);
    }
    file << "    " << s << " = " << angular_momentum << "\n";
  }
  file << "  }\n";
  file << "}\n\n";

  // which information will be saved after calculation?
  file << "Options {\n";
  file << "  WriteResultsTag = Yes\n";
  file << "  RestartFrequency = 0\n";   // does not work together with driver
  if (Config::get().energy.dftb.verbosity < 2) file << "  WriteDetailedOut = No\n";
  file << "}\n\n";

  // additional analysis that should be performed
  file << "Analysis = {\n";
  file << "  WriteBandOut = No\n";
  if (calc_type == 'g') file << "  CalculateForces = Yes\n";
  file << "}\n\n";

  // parser version (recommended so it is possible to use newer DFTB+ versions without adapting inputfile)
  file << "ParserOptions {\n";
  file << "  ParserVersion = 5\n";
  file << "}";
}

/**calculates energies and gradients
@paran if_gradient: true if gradients should be calculated, false if not*/
coords::float_type energy::interfaces::qmmm::QMMM::qmmm_calc(bool if_gradient)
{
  integrity = true;
  auto elec_factor = 332.0;
  
  mm_charge_vector = mmc.energyinterface()->charges();
  auto aco_p = dynamic_cast<energy::interfaces::aco::aco_ff const*>(mmc.energyinterface());
  if (aco_p)
  {
    elec_factor = aco_p->params().general().electric;
  }

  update_representation(); // update positions of QM and MM subsystem to those of coordinates object
  for (auto &l : link_atoms)
  {
    l.position = calc_position(l); // update positions of link atoms
    std::cout << "Link atom pos: " << l.position << "\n";
  }

  if (Config::get().energy.qmmm.qminterface == config::interface_types::T::MOPAC)
  {
    write_mol_in();
  }
  else if (Config::get().energy.qmmm.qminterface == config::interface_types::T::GAUSSIAN)
  {
    char calc_type = 'e';
    if (if_gradient == true) calc_type = 'g';
    write_gaussian_in(calc_type);
  }
  else if (Config::get().energy.qmmm.qminterface == config::interface_types::T::DFTB)
  {
    char calc_type = 'e';
    if (if_gradient == true) calc_type = 'g';
    write_dftb_in(calc_type);
  }
  else throw std::runtime_error("Chosen QM interface not implemented for QM/MM!");

  try {
    qm_energy = qmc.g();  // get energy for QM part and save gradients for QM part
    if (Config::get().energy.qmmm.qminterface == config::interface_types::T::GAUSSIAN ||
		    Config::get().energy.qmmm.qminterface == config::interface_types::T::DFTB)
    {   // electric field for QM and MM atoms (for GAUSSIAN) or coulomb gradients on MM atoms (for DFTB+)
      g_coul_mm = qmc.energyinterface()->get_g_coul_mm();  
    } 
  }
  catch(...)
  {
	std::cout << "QM programme failed. Treating structure as broken.\n";
    integrity = false;  // if QM programme fails: integrity is destroyed
  }
  if (Config::get().energy.qmmm.qminterface == config::interface_types::T::MOPAC && Config::get().energy.mopac.delete_input) std::remove("mol.in");
  
  ww_calc(if_gradient);  // calculate interactions between QM and MM part

  if (integrity == true)
  {
    if (if_gradient)  // if gradients should be calculated
    {
      mm_energy = mmc.g(); // get energy for MM part
      std::cout << "MM energy: " << mm_energy << "\n";

      // get gradients: QM + MM + vdW + Coulomb + bonded
	    auto new_grad = vdw_gradient + c_gradient + bonded_gradient;  // vdW + Coulomb + bonded
      auto g_qm = qmc.g_xyz(); // QM
      auto g_mm = mmc.g_xyz(); // MM

      std::cout << "QM-Grad: " << qmc.g_xyz(0) << "\n";
      std::cout << "bonded Grad: " << bonded_gradient[0] << "\n";
      std::cout << "vdW Grad: " << vdw_gradient << "\n";

      int counter = 0;
      for (auto&& qmi : qm_indices)
      {
        new_grad[qmi] += g_qm[new_indices_qm[qmi]];
        counter += 1;
        std::cout << "gradients of " << counter << " QM atoms read\n";
      }

      // calculate gradients from link atoms (see DOI 10.1002/(SICI)1096-987X(199703)18:4<463::AID-JCC2>3.0.CO;2-R)
      coords::Gradients_3D link_grads = qmc.energyinterface()->get_link_atom_grad();
      for (int j=0; j<link_atoms.size(); j++)
      {
        bonded::LinkAtom l = link_atoms[j];
        coords::r3 G_L = link_grads[j];
        std::cout << "Link atom: "<< G_L << "\n";

        double g = l.d_L_QM / dist(coords->xyz(l.mm), coords->xyz(l.qm));

        coords::r3 n = (coords->xyz(l.mm) - coords->xyz(l.qm)) / dist(coords->xyz(l.mm), coords->xyz(l.qm));
        std::cout << "Unit vector: " << n << "\n";

        double Fx_QM = g * scalar_product(G_L, n) *n.x() + (1 - g)*G_L.x();
        double Fy_QM = g * scalar_product(G_L, n) *n.y() + (1 - g)*G_L.y();
        double Fz_QM = g * scalar_product(G_L, n) *n.z() + (1 - g)*G_L.z();
        coords::r3 F_QM(Fx_QM, Fy_QM, Fz_QM);
        new_grad[l.qm] += F_QM;
        std::cout << "Gradient on QM atom: " << F_QM << "\n";

        double Fx_MM = g * G_L.x() - g * scalar_product(G_L, n) * n.x();
        double Fy_MM = g * G_L.y() - g * scalar_product(G_L, n) * n.y();
        double Fz_MM = g * G_L.z() - g * scalar_product(G_L, n) * n.z();
        coords::r3 F_MM(Fx_MM, Fy_MM, Fz_MM);
        new_grad[l.mm] += F_MM;
        std::cout << "Gradient on MM atom: " << F_MM << "\n";
      }

      for (auto && mmi : mm_indices)
      {
        new_grad[mmi] += g_mm[new_indices_mm[mmi]];
      }
      coords->swap_g_xyz(new_grad);
    }
    else  // only energy 
    {
      mm_energy = mmc.e();  // get energy for MM part
    }

    // energy = QM + MM + vdW + bonded (Coulomb is in QM energy)
    this->energy = qm_energy + mm_energy + vdw_energy + bonded_energy;
    if (check_bond_preservation() == false) integrity = false;
    else if (check_atom_dist() == false) integrity = false;
  }
  return energy;
}

/**calculate bonded energy and gradients*/
double energy::interfaces::qmmm::QMMM::calc_bonded(bool if_gradient)
{
  double E(0.0);
  if (if_gradient == false)  // only energy calculation
  {   // bonded gradient is given to functions because the function needs a parameter, is not used
    for (auto b : qmmm_bonds) E += b.calc_energy(coords, bonded_gradient);
    for (auto a : qmmm_angles) E += a.calc_energy(coords, bonded_gradient);
    for (auto d : qmmm_dihedrals) E += d.calc_energy(coords, torsionunit, bonded_gradient);
  }
  else   // gradient calculation
  {
    for (auto b : qmmm_bonds) E += b.calc_energy(coords, bonded_gradient, true);
    for (auto a : qmmm_angles) E += a.calc_energy(coords, bonded_gradient, true);
    for (auto d : qmmm_dihedrals) E += d.calc_energy(coords, torsionunit, bonded_gradient, true);
  }
  return E;
}

/**determines if a van der waals interaction between a QM and a MM atom should be calculated
(at least 3 bonds between those atoms)
@param qm: index of QM atom
@param mm: index of MM atom*/
bool energy::interfaces::qmmm::QMMM::calc_vdw(int qm, int mm)
{
  for (auto b : qmmm_bonds)
  {
    if (qm == b.b && mm == b.a) return false;
  }
  for (auto a : qmmm_angles)
  {
    if (qm == a.a && mm == a.b) return false;
    else if (qm == a.b && mm == a.a) return false;
  }
  return true;
}

/**calculates interaction between QM and MM part
energy is only vdW interactions
for MOPAC gradients are coulomb and vdW
for GAUSSIAN gradients are vdW and coulomb on MM atoms
@param if_gradient: true if gradients should be calculated, false if not*/
void energy::interfaces::qmmm::QMMM::ww_calc(bool if_gradient)
{
  // bonded interactions
  bonded_gradient.assign(coords->size(), coords::r3{});
  bonded_energy = calc_bonded(if_gradient);

  // preparation for calculation of non-bonded interactions
  auto elec_factor = 332.0;
  auto aco_p = dynamic_cast<energy::interfaces::aco::aco_ff const*>(mmc.energyinterface());
  if (aco_p)
  {
    elec_factor = aco_p->params().general().electric;
  }
  try { qm_charge_vector = qmc.energyinterface()->charges(); }
  catch (...) { integrity = false; }
  if (integrity == true)
  {
    c_gradient.assign(coords->size(), coords::r3{});
    vdw_gradient.assign(coords->size(), coords::r3{});

    std::size_t i2 = 0u;
    auto const & p = cparams.vdwc_matrices();
    auto const & p_vdw = p[5];
    vdw_energy = 0;

    for (auto i : qm_indices)  // for every QM atom
    {
      auto z = qmc.atoms(i2).energy_type();
      std::size_t i_vdw = cparams.type(z, tinker::potential_keys::VDW);
      coords::float_type charge_i = qm_charge_vector[i2];
      std::size_t j2 = 0u;
      for (auto j : mm_indices)  // for every MM atom
      {
        // preparation (parameters and so on)
        auto e_type = mmc.atoms(j2).energy_type();
        std::size_t j_vdw = cparams.type(e_type, tinker::potential_keys::VDW);
        auto const & p_ij = p_vdw(i_vdw, j_vdw);
        coords::float_type charge_j = mm_charge_vector[j2];
        auto r_ij = coords->xyz(j) - coords->xyz(i);
        coords::float_type d = len(r_ij);
        set_distance(d);
        coords::float_type b = (charge_i*charge_j) / d * elec_factor;
        auto R_r = std::pow(p_ij.R / d, 6);

        if (calc_vdw(i, j) == true)
        {
          std::cout << "calculate vdw energy between atoms " << i << " and " << j << "\n";
          // calculate vdW interaction
          if (cparams.general().radiustype.value ==
            ::tinker::parameter::radius_types::T::SIGMA)
          {
            vdw_energy += R_r * p_ij.E*(R_r - 1.0);
          }
          else if (cparams.general().radiustype.value ==
            ::tinker::parameter::radius_types::T::R_MIN)
          {
            vdw_energy += R_r * p_ij.E*(R_r - 2.0);
          }
          else
          {
            throw std::runtime_error("no valid radius_type");
          }
        }

        

        if (if_gradient)  // gradients
        {

          if (Config::get().energy.qmmm.qminterface == config::interface_types::T::MOPAC)
          {    // gradients of coulomb interaction (only for MOPAC here)
            coords::float_type db = b / d;
            auto c_gradient_ij = r_ij * db / d;
            c_gradient[i] += c_gradient_ij;
            c_gradient[j] -= c_gradient_ij;
          }

          if (calc_vdw(i, j) == true)
          {
            std::cout << "calculate vdw gradients between atoms " << i << " and " << j << "\n";
            // gradients of vdW interaction
            coords::float_type const V = p_ij.E*R_r;

            if (cparams.general().radiustype.value
              == ::tinker::parameter::radius_types::T::SIGMA)
            {
              auto vdw_r_grad_sigma = (V / d)*(6.0 - 12.0 * R_r);
              auto vdw_gradient_ij_sigma = (r_ij*vdw_r_grad_sigma) / d;
              vdw_gradient[i] -= vdw_gradient_ij_sigma;
              vdw_gradient[j] += vdw_gradient_ij_sigma;
              std::cout << "Gradient: " << vdw_gradient_ij_sigma << "\n";
            }
            else
            {
              auto vdw_r_grad_R_MIN = (V / d) * 12 * (1.0 - R_r);
              auto vdw_gradient_ij_R_MIN = (r_ij*vdw_r_grad_R_MIN) / d;
              vdw_gradient[i] -= vdw_gradient_ij_R_MIN;
              vdw_gradient[j] += vdw_gradient_ij_R_MIN;
              std::cout << "Gradient: " << vdw_gradient_ij_R_MIN << "\n";
            }
          }
        }
        ++j2;
      }
      ++i2;
    }

    if (Config::get().energy.qmmm.qminterface == config::interface_types::T::GAUSSIAN && if_gradient == true)
    {    // Coulomb gradients for GAUSSIAN (only on MM atoms)
      int j2 = 0;
      for (auto j : mm_indices)
      {   // additional force on MM atoms due to QM atoms (electrostatic interaction)
        double charge = mm_charge_vector[j2];
        coords::Cartesian_Point el_field = g_coul_mm[j2 + qm_indices.size()];
        double x = charge * el_field.x();
        double y = charge * el_field.y();
        double z = charge * el_field.z();
        coords::Cartesian_Point new_grad;
        new_grad.x() = -x;
        new_grad.y() = -y;
        new_grad.z() = -z;
        c_gradient[j] += new_grad;
        j2++;
      }
    }
    else if (Config::get().energy.qmmm.qminterface == config::interface_types::T::DFTB && if_gradient == true)
    {     // Coulomb gradients on MM atoms for DFTB+
      int j2 = 0;
      for (auto j : mm_indices)
      {
        c_gradient[j] += g_coul_mm[j2];
        j2++;
      }
    }
  }
}


energy::interface_base * energy::interfaces::qmmm::QMMM::clone(coords::Coordinates * c) const
{
  QMMM * tmp = new QMMM(*this, c);
  return tmp;
}

energy::interface_base * energy::interfaces::qmmm::QMMM::move(coords::Coordinates * c)
{
  QMMM * tmp = new QMMM(std::move(*this), c);
  return tmp;
}


void energy::interfaces::qmmm::QMMM::swap(interface_base& rhs)
{
  swap(dynamic_cast<QMMM&>(rhs));
}

void energy::interfaces::qmmm::QMMM::swap(QMMM& rhs)
{
  interface_base::swap(rhs);
  std::swap(cparams, rhs.cparams);
  qm_indices.swap(rhs.qm_indices);
  mm_indices.swap(rhs.mm_indices);
  new_indices_mm.swap(rhs.new_indices_mm);
  new_indices_qm.swap(rhs.new_indices_qm);
  qmc.swap(rhs.qmc);
  mmc.swap(rhs.mmc);
  qm_charge_vector.swap(rhs.qm_charge_vector);
  mm_charge_vector.swap(rhs.mm_charge_vector);
  std::swap(vdw_energy, rhs.vdw_energy);
  std::swap(qm_energy, rhs.qm_energy);
  std::swap(mm_energy, rhs.mm_energy);
  c_gradient.swap(rhs.c_gradient);
  vdw_gradient.swap(rhs.vdw_gradient);
}

// update structure (account for topology or rep change)
void energy::interfaces::qmmm::QMMM::update(bool const skip_topology)
{
  if (!skip_topology)
  {
    *this = QMMM(this->coords);
  }
  else
  {
    update_representation();
    qmc.energy_update(true);
    mmc.energy_update(true);
  }
}

coords::float_type energy::interfaces::qmmm::QMMM::g()
{
  return qmmm_calc(true);
}

coords::float_type energy::interfaces::qmmm::QMMM::e()
{
  return qmmm_calc(false);
}

coords::float_type energy::interfaces::qmmm::QMMM::h()
{
  throw std::runtime_error("no QMMM-function yet");
}

coords::float_type energy::interfaces::qmmm::QMMM::o()
{
  throw std::runtime_error("QMMM-cannot optimize");
}

std::vector<coords::float_type> energy::interfaces::qmmm::QMMM::charges() const
{
  std::vector<coords::float_type> v;
  v.assign(coords->size(), 0);
  std::size_t i = 0;
  std::size_t j = 0;
  if (qm_charge_vector.size() + mm_charge_vector.size() == coords->size())
  {
    for (auto && qmi : qm_indices)
    {
      v[qmi] = qm_charge_vector[i];
      i++;
    }

    for (auto && mmi : mm_indices)
    {
      v[mmi] = mm_charge_vector[j];
      j++;
    }

    if (v.size() != coords->size())
    {
      throw std::logic_error("Found " + std::to_string(v.size()) +
        " charges instead of " + std::to_string(coords->size()) + " charges.");
    }
  }
  else
  {
    throw std::logic_error("QM/MM-charges not yet calculated");
  }
  return v;
}

void energy::interfaces::qmmm::QMMM::print_E(std::ostream &) const
{
  throw std::runtime_error("no QMMM-function yet");
}

void energy::interfaces::qmmm::QMMM::print_E_head(std::ostream &S, bool const endline) const
{
  S << "QM-atoms: " << qm_indices.size() << '\n';
  S << "MM-atoms: " << mm_indices.size() << '\n';
  S << "Potentials\n";
  S << std::right << std::setw(24) << "QM";
  S << std::right << std::setw(24) << "MM";
  S << std::right << std::setw(24) << "VDW";
  S << std::right << std::setw(24) << "BONDED";
  S << std::right << std::setw(24) << "TOTAL";
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::QMMM::print_E_short(std::ostream &S, bool const endline) const
{
  S << '\n';
  S << std::right << std::setw(24) << qm_energy;
  S << std::right << std::setw(24) << mm_energy;
  S << std::right << std::setw(24) << vdw_energy;
  S << std::right << std::setw(24) << bonded_energy;
  S << std::right << std::setw(24) << energy;
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::QMMM::print_gnuplot(std::ostream &S, bool const endline) const
{
  
  S << std::right << std::setw(24) << distance;
  S << std::right << std::setw(24) << vdw_gradient + c_gradient;
  if (endline) S << '\n';
}

void energy::interfaces::qmmm::QMMM::print_G_tinkerlike(std::ostream &S, bool const endline) const
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
      + coords->g_xyz(k).z() * coords->g_xyz(k).z()) << std::endl;
  }
  
}

void energy::interfaces::qmmm::QMMM::to_stream(std::ostream &S) const
{
  S << '\n';
  interface_base::to_stream(S);
  throw std::runtime_error("no QMMM-function yet");
}

bool energy::interfaces::qmmm::QMMM::check_bond_preservation(void) const
{
  std::size_t const N(coords->size());
  for (std::size_t i(0U); i < N; ++i)
  { // cycle over all atoms i
    if (!coords->atoms(i).bonds().empty())
    {
      std::size_t const M(coords->atoms(i).bonds().size());
      for (std::size_t j(0U); j < M && coords->atoms(i).bonds(j) < i; ++j)
      { // cycle over all atoms bound to i
        double const L(geometric_length(coords->xyz(i) - coords->xyz(coords->atoms(i).bonds(j))));
        if (L > 2.2) return false;
      }
    }
  }
  return true;
}

bool energy::interfaces::qmmm::QMMM::check_atom_dist(void) const
{
  std::size_t const N(coords->size());
  for (std::size_t i(0U); i < N; ++i)
  {
    for (std::size_t j(0U); j < i; j++)
    {
      if (dist(coords->xyz(i), coords->xyz(j)) < 0.3)
      {
        return false;
      }
    }
  }
  return true;
}
