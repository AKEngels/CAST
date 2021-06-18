#include "pmf_ic_prep.h"

#include "pmf_interpolator_builder.h"

#include "InternalCoordinates/InternalCoordinateDecorator.h"
#include "InternalCoordinates/InternalCoordinateUtilities.h"
#include "InternalCoordinates/InternalCoordinates.h"

#include "ic_exec.h"

auto build_ic_system(coords::Coordinates& c) {
  std::vector<std::string> el_vec2;
  for (auto&& i : c.atoms())
  {
    el_vec2.emplace_back(i.symbol());
  }
  auto bonds = ic_util::bonds(el_vec2, c.xyz());

  struct graphinfo
  {
    std::size_t atom_serial;
    std::string atom_name, element;
    coords::Cartesian_Point cp;
  };
  std::vector<graphinfo> curGraphinfo;
  for (std::size_t i = 0u; i < c.xyz().size(); i++)
  {
    graphinfo tempinfo;
    tempinfo.cp = c.xyz(i);
    tempinfo.atom_serial = i + 1u;
    tempinfo.atom_name = c.atoms(i).symbol();
    tempinfo.element = tempinfo.atom_name;
    curGraphinfo.emplace_back(tempinfo);
  }

  return build_z_matrix_coords(ic_util::make_graph(bonds, curGraphinfo));
}

pmf_ic_prep::pmf_ic_prep(coords::Coordinates& c, coords::input::format& ci, std::string const& outfile, std::string const& splinefile) :
  coordobj(c), coord_input(&ci), outfilename(outfile), splinefilename(splinefile), dimension(Config::get().coords.umbrella.pmf_ic.indices_xi.size()),
  ic_system_(build_ic_system(c))
{
  std::tie(rc_, rc_index_) = ic_system_->get_coord_for_atom_indices(Config::get().coords.umbrella.pmf_ic.indices_xi[0]);
}

void pmf_ic_prep::run()
{
  calc_energies();
  write_to_file();
  if (dimension == 1) write_spline_1d();
  else write_spline_2d(); 
}

void pmf_ic_prep::calc_energies()
{
  // catch invalid low level interface
  if (Config::get().coords.umbrella.pmf_ic.LL_interface == config::interface_types::ILLEGAL) {
    throw std::runtime_error("Illegal low level interface!");
  }

  // create new coordinates object with low level interface
  Config::set().general.energy_interface = Config::get().coords.umbrella.pmf_ic.LL_interface;
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords_ll(ci->read(Config::get().general.inputFilename));

  for (auto const& pes : *coord_input)   // for every structure
  {
    coordobj.set_xyz(pes.structure.cartesian, true);
    coords_ll.set_xyz(pes.structure.cartesian, true);

    // calulate xi and z
    //auto xi = coords::bias::Potentials::calc_xi(coordobj.xyz(), Config::get().coords.umbrella.pmf_ic.indices_xi[0]);
    auto xi = rc_->val(coordobj.xyz()) * SCON_180PI;
    double xi_2;

    if (dimension == 1)  // in case of 1D
    {
      xis.emplace_back(xi);
    }
    if (dimension > 1)    // in case of 2D
    {
      xi_2 = coords::bias::Potentials::calc_xi(coordobj.xyz(), Config::get().coords.umbrella.pmf_ic.indices_xi[1]);
      xi_2d.emplace_back(std::make_pair( xi, xi_2 ));
    }

    // calculate energies
    auto E_hl = coordobj.g();
    E_HLs.emplace_back(E_hl);

    auto E_ll = coords_ll.g();
    E_LLs.emplace_back(E_ll);
    deltaEs.emplace_back(E_hl - E_ll);
    auto curr_grad = calc_gradient_difference(coordobj.xyz(), coordobj.g_xyz(), coords_ll.g_xyz());
    grad_Es.emplace_back(curr_grad);
    if (Config::get().general.verbosity > 3)
    {
      std::cout << xi << " , ";
      if (dimension > 1) std::cout << xi_2 << " , ";
      std::cout << E_hl << " , " << E_ll << " , " << curr_grad << "\n";
    }
  }
  if (Config::get().general.verbosity > 1) std::cout << "finished energy calculation\n";
}

void pmf_ic_prep::write_to_file()
{
  std::ofstream outfile(outfilename, std::ios_base::out);
  outfile.precision(10);
  outfile << "xi,";
  if (dimension > 1) outfile << "xi_2,";
  outfile<<"E_HL, E_LL, deltaE, dE_delta";
  for (auto i{ 0u }; i < E_HLs.size(); ++i)
  {
    if (dimension == 1) outfile << "\n" << xis[i] << ","<<E_HLs[i] << "," << E_LLs[i] << "," << deltaEs[i] << ',' << grad_Es[i];
    else outfile << "\n" << xi_2d[i].first << "," <<xi_2d[i].second << "," << E_HLs[i] << "," << E_LLs[i] << "," << deltaEs[i];
  }
  outfile.close();
}

void pmf_ic_prep::write_spline_1d()
{
  auto interpolator = pmf_ic::build_interpolator(xis, deltaEs, std::optional(grad_Es));

  // write spline to file
  std::ofstream splinefile(splinefilename, std::ios_base::out);
  splinefile.precision(10);
  splinefile << "xi,interpolated";   // headline
  auto const& start = Config::get().coords.umbrella.pmf_ic.ranges[0].start;
  auto const& stop = Config::get().coords.umbrella.pmf_ic.ranges[0].stop;
  auto const& step = Config::get().coords.umbrella.pmf_ic.ranges[0].step;
  for (auto xi{ start }; xi <= stop; xi += step)
  {
    auto y = interpolator->get_value(xi);
    splinefile << "\n" << xi << "," << y;
  }
  splinefile.close();
}

void pmf_ic_prep::write_spline_2d()
{
  auto interpolator = pmf_ic::build_interpolator(xi_2d, deltaEs);

  auto gpr = gpr::gpr_interpolator_2d(std::make_unique<gpr::SqExpKernel>(10), xi_2d, deltaEs);

  // write spline to file
  std::ofstream splinefile(splinefilename, std::ios_base::out);
  splinefile.precision(10);
  auto const& start1 = Config::get().coords.umbrella.pmf_ic.ranges[0].start;
  auto const& stop1 = Config::get().coords.umbrella.pmf_ic.ranges[0].stop;
  auto const& step1 = Config::get().coords.umbrella.pmf_ic.ranges[0].step;

  for (auto xi{ start1 }; xi <= stop1; xi += step1) {
    splinefile << "," << xi;  // headline
  }
  splinefile << ",xi_1";

  auto const& start2 = Config::get().coords.umbrella.pmf_ic.ranges[1].start;
  auto const& stop2 = Config::get().coords.umbrella.pmf_ic.ranges[1].stop;
  auto const& step2 = Config::get().coords.umbrella.pmf_ic.ranges[1].step;

  for (auto xi2{ start2 }; xi2 <= stop2; xi2 += step2)    // rows = xi_2
  {
    splinefile << "\n" << xi2;
    for (auto xi1{ start1 }; xi1 <= stop1; xi1 += step1)  // columns = xi_1
    {
      splinefile << "," << interpolator->get_value(xi1, xi2);
    }
  }
  splinefile << "\nxi_2";
}

double pmf_ic_prep::calc_gradient_difference(coords::Representation_3D const& xyz,
                                             coords::Gradients_3D const& grad_hl,
                                             coords::Gradients_3D const& grad_ll) {
  auto B = ic_system_->Bmat(xyz);
  auto G = ic_system_->Gmat(xyz);
  auto grad_vec_hl = scon::mathmatrix<double>::col_from_vec(ic_util::flatten_c3_vec(grad_hl));
  auto grad_vec_ll = scon::mathmatrix<double>::col_from_vec(ic_util::flatten_c3_vec(grad_ll));
  auto grad_vec = grad_vec_hl - grad_vec_ll;
  grad_vec /= energy::Hartree_Bohr2Kcal_MolAng;
  return G.solve(B * grad_vec)(rc_index_, 0);
}
