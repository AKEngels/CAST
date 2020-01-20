#include "md.h"
#include "md_thermostat.h"

// determine target temperature for heating
bool md::simulation::determine_current_desired_temperature(std::size_t const step, bool fep)
{
  if (Config::get().md.heat_steps.size() == 0)  // no temperature control
  {
    return false;
  }
  else  // temperature control
  {
    config::md_conf::config_heat last;
    last.raise = Config::get().md.T_init;
    if (fep == true)   // keep constant temperature in FEP calculation
    {
      desired_temp = Config::get().md.T_final;
      return true;
    }
    for (auto const& heatstep : Config::get().md.heat_steps)
    {
      if (heatstep.offset >= step)    // find first heatstep after current step
      {
        double const delta((heatstep.raise - last.raise) / static_cast<double>(heatstep.offset - last.offset));
        desired_temp += delta;  // adjust target temperature
        return true; // exit function
      }
      last = heatstep; // last heatstep before current step
    }
    if (step > Config::get().md.heat_steps[Config::get().md.heat_steps.size() - 1].offset)
    {             // after last heatstep: keep final temperature
      desired_temp = Config::get().md.T_final;
      return true;
    }
    else
    {
      std::cout << "This should not happen! Error in applying heating correction of thermostat.\n";
      throw std::exception();
    }
  }
}

// Nose-Hover thermostat. Variable names and implementation are identical to the book of
// Frenkel and Smit, Understanding Molecular Simulation, Appendix E
double md::simulation::nose_hoover_thermostat(void)
{
  return this->nose_hoover_thermostat_some_atoms(std::vector<size_t>(this->coordobj.size(), 1));
}

// Nose-Hover thermostat for inner atoms. Variable names and implementation are identical to the book of
// Frenkel and Smit, Understanding Molecular Simulation, Appendix E
// Note by Dustin okt2019: This implementation describes
// two Nose-Hoover Chains using Trotter Factorization of the Liouville Operator
double md::simulation::nose_hoover_thermostat_some_atoms(std::vector<size_t> active_atoms)
{
  md::nose_hoover_2chained & nht = this->thermostat.nht_2chained;
  double tempscale(0.0);
  std::size_t freedom_some = 3U * active_atoms.size();
  if (Config::get().periodics.periodic == true)
    freedom_some -= 3;
  else
    freedom_some -= 6;
  double const delt(Config::get().md.timeStep),
    d2(delt / 2.0), d4(d2 / 2.0), d8(d4 / 2.0),
    TR(desired_temp * md::R), fTR(freedom_some * TR);
  nht.G2 = (nht.Q1 * nht.v1 * nht.v1 - TR) / nht.Q2;
  nht.v2 += nht.G2 * d4;
  nht.v1 *= exp(-nht.v2 * d8);
  nht.G1 = (2.0 * E_kin - fTR) / nht.Q1;
  nht.v1 += nht.G1 * d4;
  nht.v1 *= exp(-nht.v2 * d8);
  nht.x1 += nht.v1 * d2;
  nht.x2 += nht.v2 * d2;
  tempscale = exp(-nht.v1 * d2);
  E_kin *= tempscale * tempscale;
  if (Config::get().general.verbosity >= 5u)
  {
    std::cout << "Nose-Hoover-Adjustment; Scaling factor: " << tempscale << '\n';
  }
  nht.v1 *= exp(-nht.v2 * d8);
  nht.G1 = (2.0 * E_kin - fTR) / nht.Q1;
  nht.v1 += nht.G1 * d4;
  nht.v1 *= exp(-nht.v2 * d8);
  nht.G2 = (nht.Q1 * nht.v1 * nht.v1 - TR) / nht.Q2;
  nht.v2 += nht.G2 * d4;
  return tempscale;
}

double md::simulation::nose_hoover_with_arbitrary_chain_length(std::vector<size_t> active_atoms, std::size_t const n_ys)
{
  using float_type = double;
  std::vector<float_type> omegas;
  md::nose_hoover_arbitrary_length & nht2 = this->thermostat.nht_v2;
  const std::size_t chainlength = nht2.chainlength;
  if (n_ys == 3u)
  {
    omegas = {1.0/(2.0-std::pow(2.,1./3.)),1-2.0/(2.0 - std::pow(2.,1. / 3.)),1.0 / (2.0 - std::pow(2.,1. / 3.)) };
  }
  else if (n_ys == 5u)
  {
    omegas =
    { 
      1.0 / (4.0 - std::pow(4.,1. / 3.)),
      1.0 / (4.0 - std::pow(4.,1. / 3.)),
      1 - 4.0 / (4.0 - std::pow(4.,1. / 3.)),
      1.0 / (4.0 - std::pow(4.,1. / 3.)),
      1.0 / (4.0 - std::pow(4.,1. / 3.)) 
    };
  }
  else
  {
    throw std::logic_error("Arbitrary chain-length Nose-Hoover thermostat is only implemented for n_ys=3 and n_ys=5. See also: DOI 10.1080/00268979600100761 . Exiting.");
  }
  float_type factor = 1.;
  constexpr float_type k_B = constants::gas_constant_R_kcal_per_mol_kelvin; // is per mole here, equal to gas constant R, as energy is in kcal/mol
  const float_type kBT = k_B * this->desired_temp;
  std::size_t freedom_some = 3U * active_atoms.size();
  if (Config::get().periodics.periodic == true)
    freedom_some -= 3u;
  else
    freedom_some -= 6u;
  nht2.forces.at(0u) = (2.*this->E_kin /*kcal per mole*/ - freedom_some * kBT /*kcal*/) / nht2.masses_param_Q.at(0);

  for (std::size_t iresn = 1u; iresn <= chainlength; iresn++) // FORTRAN stlye loops
  {
    for (std::size_t iyosh = 1u; iyosh <= n_ys; iyosh++) // FORTRAN stlye loops
    {
      const float_type wdti = omegas.at(iyosh - 1u) * dt / static_cast<float_type>(chainlength);
      const std::size_t currentArrayAccess = chainlength - 1u;
      nht2.velocities.at(currentArrayAccess) += nht2.forces.at(currentArrayAccess) * wdti/4.;
      for (std::size_t inos = 1u; inos <= chainlength - 1u; inos++)
      {
        const float_type AA = std::exp(-wdti/8. * nht2.velocities.at(chainlength - inos));
        nht2.velocities.at(chainlength - 1u - inos) *= AA*AA;
        nht2.velocities.at(chainlength - 1u - inos) += wdti / 4. * AA * nht2.forces.at(chainlength - 1u - inos);
      }
      float_type AA = std::exp((wdti/2. * -1.0 * nht2.velocities.at(0u)));
      factor *=AA;
      this->E_kin *= factor * factor;
      nht2.forces.at(0u) = (2.*this->E_kin - freedom_some * kBT) / nht2.masses_param_Q.at(0u);
      for (std::size_t inos = 1u; inos <= chainlength - 1u; inos++)
      {
        nht2.epsilons.at(inos -1u) = nht2.epsilons.at(inos - 1u) + nht2.velocities.at(inos - 1u) * wdti/2.;
      }
      for (std::size_t inos = 1u; inos <= chainlength - 1u; inos++)
      {
        AA = std::exp((wdti / 8. * -1.0 * nht2.velocities.at(inos)));
        nht2.velocities.at(inos - 1u) = nht2.velocities.at(inos - 1u) * AA * AA + wdti/4. * AA * nht2.forces.at(inos - 1u);
        nht2.forces.at(inos) = (nht2.masses_param_Q.at(inos - 1u)*nht2.velocities.at(inos - 1u)* nht2.velocities.at(inos - 1u)-kBT) / nht2.masses_param_Q.at(inos);
      }
      nht2.velocities.at(chainlength - 1u) = nht2.velocities.at(chainlength - 1u) + nht2.forces.at(chainlength - 1u) * wdti/4.;
    }
  }  
  return factor;
}

double md::simulation::tempcontrol(config::molecular_dynamics::thermostat_algorithms::T thermostat, bool half)
{
  using thermoalgo = config::molecular_dynamics::thermostat_algorithms;
  std::size_t const N = this->coordobj.size();  // total number of atoms
  if (Config::get().md.set_active_center == 1)
    updateEkin(inner_atoms);
  else if (Config::get().coords.fixed.size() != 0)
    updateEkin(movable_atoms);
  else
    updateEkin(range(N));
  double temp_after_scaling = 0.;
  const double instantaneous_temp_after_last_scaling = this->instantaneous_temp > 0. ? this->instantaneous_temp : Config::get().md.T_init;
  
  double scaling_factor = 1.;
  size_t dof = freedom;
  if (Config::get().md.set_active_center == 1)  // if biased potential
  {
    dof = 3u * inner_atoms.size();
    if (Config::get().periodics.periodic == true)
      dof -= 3;
    else
      dof -= 6;
  }
  else if (Config::get().coords.fixed.size() != 0)
  {
    dof = 3u * movable_atoms.size();
    if (Config::get().periodics.periodic == true)
      dof -= 3;
    else
      dof -= 6;
  }
  const double T_factor = (2.0 / (dof * md::R));
  double instantaneous_temp_before_scaling = this->E_kin * T_factor;
  if (Config::get().md.set_active_center == 1)  // if biased potential
  {
    double factor = 1.0;
    if (thermostat == thermoalgo::TWO_NOSE_HOOVER_CHAINS)
      factor = nose_hoover_thermostat_some_atoms(inner_atoms);     // calculate temperature scaling factor
    else if (thermostat == thermoalgo::ARBITRARY_CHAIN_LENGTH_NOSE_HOOVER)
      factor = nose_hoover_with_arbitrary_chain_length(inner_atoms);
    else if (thermostat == thermoalgo::HOOVER_EVANS)
    {
      factor = std::sqrt(instantaneous_temp_after_last_scaling / instantaneous_temp_before_scaling);
    }
    else if (thermostat == thermoalgo::BERENDSEN)
    {
      const double prefactor = (static_cast<double>(dof) - 1.) / static_cast<double>(dof); // See: Adv. Polym. Sci. (2005) 173:105–149 DOI:10.1007 / b99427 | Hünenberger | P.126
      double factor = std::sqrt(prefactor * desired_temp / instantaneous_temp_before_scaling);
      factor -= 1.0;
      factor *= dt * 1. / this->thermostat.berendsen_tB;
      factor += 1;
      factor = std::sqrt(factor);
      scaling_factor = factor;
    }
    else
    {
      const double prefactor = (static_cast<double>(dof) - 1.) / static_cast<double>(dof); // See: Adv. Polym. Sci. (2005) 173:105–149 DOI:10.1007 / b99427 | Hünenberger | P.126
      factor = std::sqrt(prefactor * desired_temp / instantaneous_temp_before_scaling);
    }
    scaling_factor = factor;
  }
  else if (Config::get().coords.fixed.size() != 0)
  {
    updateEkin(movable_atoms);
    double factor = 1.0;
    if (thermostat == thermoalgo::TWO_NOSE_HOOVER_CHAINS)
      factor = nose_hoover_thermostat_some_atoms(movable_atoms);     // calculate temperature scaling factor
    else if (thermostat == thermoalgo::ARBITRARY_CHAIN_LENGTH_NOSE_HOOVER)
      factor = nose_hoover_with_arbitrary_chain_length(movable_atoms);
    else if (thermostat == thermoalgo::HOOVER_EVANS)
    {
      factor = std::sqrt(instantaneous_temp_after_last_scaling / instantaneous_temp_before_scaling);
    }
    else if (thermostat == thermoalgo::BERENDSEN)
    {
      const double prefactor = (static_cast<double>(dof) - 1.) / static_cast<double>(dof); // See: Adv. Polym. Sci. (2005) 173:105–149 DOI:10.1007 / b99427 | Hünenberger | P.126
      double factor = std::sqrt(prefactor * desired_temp / instantaneous_temp_before_scaling);
      factor -= 1.0;
      factor *= dt * 1. / this->thermostat.berendsen_tB;
      factor += 1;
      factor = std::sqrt(factor);
      scaling_factor = factor;
    }
    else
    {
      const double prefactor = (static_cast<double>(dof) - 1.) / static_cast<double>(dof); // See: Adv. Polym. Sci. (2005) 173:105–149 DOI:10.1007 / b99427 | Hünenberger | P.126
      factor = std::sqrt(prefactor * desired_temp / instantaneous_temp_before_scaling);
    }
      
    scaling_factor = factor;
  }
  else if (thermostat == thermoalgo::TWO_NOSE_HOOVER_CHAINS)
  {
    const double factor = nose_hoover_thermostat();
    scaling_factor = factor;
  }
  else if (thermostat == thermoalgo::ARBITRARY_CHAIN_LENGTH_NOSE_HOOVER)
  {
    const double factor = nose_hoover_with_arbitrary_chain_length(range(N));;
    scaling_factor = factor;
  }
  else if (thermostat == thermoalgo::HOOVER_EVANS)
  {
    const double factor = std::sqrt(instantaneous_temp_after_last_scaling / instantaneous_temp_before_scaling);
    scaling_factor = factor;
  }
  else if (thermostat == thermoalgo::BERENDSEN)
  {
    const double prefactor = (static_cast<double>(dof) - 1.) / static_cast<double>(dof); // See: Adv. Polym. Sci. (2005) 173:105–149 DOI:10.1007 / b99427 | Hünenberger | P.126
    double factor = std::sqrt(prefactor * desired_temp / instantaneous_temp_before_scaling);
    factor -= 1.0;
    factor *= dt * 1./this->thermostat.berendsen_tB;
    factor += 1;
    factor = std::sqrt(factor);
    scaling_factor = factor;
  }
  else if (thermostat == thermoalgo::VELOCITY_RESCALING)
  {
    //const double prefactor = (static_cast<double>(dof) - 1.) / static_cast<double>(dof); // See: Adv. Polym. Sci. (2005) 173:105–149 DOI:10.1007 / b99427 | Hünenberger | P.126
    constexpr double prefactor = 1.;
    const double factor = std::sqrt(prefactor * desired_temp / instantaneous_temp_before_scaling);
    scaling_factor = factor;
  }
  else if (thermostat == thermoalgo::ANDERSEN)
  {
    // Keep this because we rescale atoms individually in this case!!
    scaling_factor = 1.;
    this->thermostat.normalDist = std::normal_distribution<double>{ 0, 1 }; // normal distribution with mean=0 and standard deviation=1
    for (auto i : movable_atoms)
    {
      double unifRandNum = this->thermostat.unifDist(this->thermostat.randomEngine);
      if (unifRandNum < this->thermostat.andersen_parameter * this->dt)
      {
        double newVelocity1 = this->thermostat.normalDist(this->thermostat.randomEngine) * std::sqrt(gasconstant_R_1 * desired_temp / M[i]);
        double newVelocity2 = this->thermostat.normalDist(this->thermostat.randomEngine) * std::sqrt(gasconstant_R_1 * desired_temp / M[i]);
        double newVelocity3 = this->thermostat.normalDist(this->thermostat.randomEngine) * std::sqrt(gasconstant_R_1 * desired_temp / M[i]);
        V[i] = scon::c3(newVelocity1, newVelocity2, newVelocity3);   // new velocity
      }
    }
  }
  for (auto i : movable_atoms)
  {
    V[i] *= scaling_factor;   // new velocities (for all atoms that have a velocity)
  }
  if (Config::get().md.set_active_center == 1)
    updateEkin(inner_atoms);
  else if (Config::get().coords.fixed.size() != 0)
    updateEkin(movable_atoms);
  else
    updateEkin(range(N));
  temp_after_scaling = this->E_kin * T_factor;
  if (Config::get().general.verbosity > 3 && half)
  {
    std::cout << "THERMOCONTROL - Half step: desired Temp: " << desired_temp << "K. current temp: " << temp_after_scaling << "K. Velocity Scaling Factor: " << scaling_factor << "\n";
  }
  else if (Config::get().general.verbosity > 3)
  {
    std::cout << "THERMOCONTROL - Full step: desired Temp: " << desired_temp << "K. current temp: " << temp_after_scaling << "K. Velocity Scaling Factor: " << scaling_factor << "\n";
  }
  return temp_after_scaling;
}
