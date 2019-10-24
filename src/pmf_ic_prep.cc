#include "pmf_ic_prep.h"

pmf_ic_prep::pmf_ic_prep(coords::Coordinates& c, coords::input::format& ci, std::string const& outfile, std::string const& splinefile) :
  coordobj(c), coord_input(&ci), outfilename(outfile), splinefilename(splinefile) {}

void pmf_ic_prep::run()
{
  calc_xis_zs_and_E_HLs();
  calc_E_LLs();
  calc_deltaEs();
  write_to_file();
  write_spline();
}

void pmf_ic_prep::calc_xis_zs_and_E_HLs()
{
  for (auto const& pes : *coord_input)   // for every structure
  {
    coordobj.set_xyz(pes.structure.cartesian, true);
    auto xi = coords::bias::Potentials::calc_xi(coordobj.xyz());
    xis.emplace_back(xi);
    zs.emplace_back(mapping::xi_to_z(xi));
    auto E = coordobj.e();
    E_HLs.emplace_back(E);
    if (Config::get().general.verbosity > 3) std::cout << xi << " , " << E << "\n";
  }
  if (Config::get().general.verbosity > 1) std::cout << "finished high level calculation\n";
}

void pmf_ic_prep::calc_E_LLs()
{ 
  // set energy to low level method
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

void pmf_ic_prep::calc_deltaEs()
{
  for (auto i{ 0u }; i < xis.size(); ++i)
  {
    deltaEs.emplace_back(E_HLs[i] - E_LLs[i]);
  }
}

void pmf_ic_prep::write_to_file()
{
  std::ofstream outfile(outfilename, std::ios_base::out);
  outfile << "xi,z,E_HL,E_LL,deltaE";                            // headline
  for (auto i{ 0u }; i < xis.size(); ++i)
  {
    outfile << "\n" << xis[i] << "," << zs[i] << "," << E_HLs[i] << "," << E_LLs[i] << "," << deltaEs[i];
  }
  outfile.close();
}

void pmf_ic_prep::write_spline()
{
  Spline s(zs, deltaEs);   // create spline

  // write spline to file
  std::ofstream splinefile(splinefilename, std::ios_base::out);
  splinefile << "xi,spline";   // headline
  auto const& start = Config::get().coords.umbrella.pmf_ic.start;
  auto const& stop = Config::get().coords.umbrella.pmf_ic.stop;
  auto const& step = Config::get().coords.umbrella.pmf_ic.step;
  for (auto xi{ start }; xi <= stop; xi += step)
  {
    auto z = mapping::xi_to_z(xi);
    auto y = s.get_value(z);
    splinefile << "\n" << xi << "," << y;
  }
  splinefile.close();
}
