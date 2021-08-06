#include "pmf_ic_prep.h"

#include "pmf_interpolator_builder.h"

pmf_ic_prep::pmf_ic_prep(coords::Coordinates& c, coords::input::format& ci, std::string const& outfile, std::string const& splinefile) :
        pmf_ic_base(c, ci), outfilename(outfile), splinefilename(splinefile)
{
  if (Config::get().coords.umbrella.pmf_ic.reference_index >= coord_input->PES().size())
    throw std::runtime_error("Umbrella sampling reference index out of range. Check CAST.txt");
}

void pmf_ic_prep::run()
{
  calc_xis_zs_and_E_HLs();
  calc_E_LLs();
  calc_deltaEs();
  write_to_file();
  if (dimension == 1) write_spline_1d();
  else write_spline_2d(); 
}

void pmf_ic_base::calc_xis_zs_and_E_HLs()
{
  for (auto const& pes : *coord_input)   // for every structure
  {
    coordobj.set_xyz(pes.structure.cartesian, true);

    // calulate xi and z
    auto xi = coords::bias::Potentials::calc_xi(coordobj.xyz(), Config::get().coords.umbrella.pmf_ic.indices_xi[0]);
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

    // calculate high level energy
    auto E = coordobj.e();
    E_HLs.emplace_back(E);
    if (Config::get().general.verbosity > 3)
    {
      std::cout << xi << " , ";
      if (dimension > 1) std::cout << xi_2 << " , ";
      std::cout << E << "\n";
    }
  }
  if (Config::get().general.verbosity > 1) std::cout << "finished high level calculation\n";
}

void pmf_ic_base::calc_E_LLs()
{ 
  // catch unvalid low level interface
  if (Config::get().coords.umbrella.pmf_ic.LL_interface == config::interface_types::ILLEGAL) {
    throw std::runtime_error("Illegal low level interface!");
  }

  // create new coordinates object with low level interface
  Config::set().general.energy_interface = Config::get().coords.umbrella.pmf_ic.LL_interface;
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());
  coords::Coordinates coords(ci->read(Config::get().general.inputFilename));

  // calculate low level energies
  for (auto const& pes : *ci)   // for every structure
  {
    coords.set_xyz(pes.structure.cartesian, true);
    auto E = coords.e();
    E_LLs.emplace_back(E);
    if (Config::get().general.verbosity > 3) std::cout << E << "\n";
  }
  if (Config::get().general.verbosity > 1) std::cout << "finished low level calculation\n";
}

void pmf_ic_base::calc_deltaEs()
{
  auto ref_index = Config::get().coords.umbrella.pmf_ic.reference_index;
  auto ref_delta = E_LLs[ref_index] - E_HLs[ref_index];
  for (auto i{ 0u }; i < E_HLs.size(); ++i)
  {
    deltaEs.emplace_back(E_HLs[i] + ref_delta - E_LLs[i]);
  }
}

pmf_ic_base::pmf_ic_base(coords::Coordinates& coords, coords::input::format& ci)
: coordobj(coords), coord_input(&ci), dimension(Config::get().coords.umbrella.pmf_ic.indices_xi.size())
{
}

void pmf_ic_prep::write_to_file()
{
  std::ofstream outfile(outfilename, std::ios_base::out);
  outfile.precision(10);
  outfile << "xi,";
  if (dimension > 1) outfile << "xi_2,";
  outfile<<"E_HL, E_LL, deltaE";                            
  for (auto i{ 0u }; i < E_HLs.size(); ++i)
  {
    if (dimension == 1) outfile << "\n" << xis[i] << ","<<E_HLs[i] << "," << E_LLs[i] << "," << deltaEs[i];
    else outfile << "\n" << xi_2d[i].first << "," <<xi_2d[i].second << "," << E_HLs[i] << "," << E_LLs[i] << "," << deltaEs[i];
  }
  outfile.close();
}

void pmf_ic_prep::write_spline_1d()
{
  auto interpolator = pmf_ic::build_interpolator(xis, deltaEs);

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

  auto gpr = gpr::gpr_interpolator_2d(std::make_unique<gpr::SqExpCovariance>(10), xi_2d, deltaEs);

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

pmf_ic_test::pmf_ic_test(coords::Coordinates& coords, coords::input::format& ci) : pmf_ic_base(coords, ci), interpolator_(pmf_ic::load_interpolation())
{
}

void pmf_ic_test::run() {
  calc_xis_zs_and_E_HLs();
  calc_E_LLs();
  calc_deltaEs();
  calc_interpolation_errors();
}

void pmf_ic_test::calc_interpolation_errors() {
  std::vector<double> errors;
  std::cout << "Interpolation Errors:\n";
  for (std::size_t i=0; i<deltaEs.size(); ++i) {
    auto interpolated = [this, i](){
      if (dimension == 1) {
        auto xi = xis[i];
        std::cout << xi;
        auto const& interpolator = *std::get<std::unique_ptr<pmf_ic::Interpolator1DInterface>>(interpolator_);
        return interpolator.get_value(xi);
      } else {
        auto [xi1, xi2] = xi_2d[i];
        std::cout << '\n' << xi1 << ',' << xi2;
        auto const& interpolator = *std::get<std::unique_ptr<pmf_ic::Interpolator2DInterface>>(interpolator_);
        return interpolator.get_value(xi1, xi2);
      }
    }();
    auto curr_error = std::abs(interpolated - deltaEs[i]);
    std::cout << ',' << curr_error << '\n';
    errors.emplace_back(curr_error);
  }
  auto [min, max] = std::minmax_element(errors.begin(), errors.end());
  auto avg_error = std::transform_reduce(errors.begin(), errors.end(), 0., std::plus{}, [](auto x){return std::abs(x);}) / errors.size();
  auto rmsd = std::sqrt(std::transform_reduce(errors.begin(), errors.end(), 0., std::plus{}, [](auto x){return std::pow(x, 2);}) / errors.size());

  std::cout << "Min: " << *min << "; Max: " << *max << "; Average: " << avg_error << "; RMSD: " << rmsd << '\n';
}
