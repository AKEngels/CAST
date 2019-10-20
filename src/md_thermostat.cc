#include "md.h"
#include "md_thermostat.h"
#pragma once

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
  return this->nose_hoover_thermostat_some_atoms(std::vector<size_t>(this->freedom, 1));
}

// Nose-Hover thermostat for inner atoms. Variable names and implementation are identical to the book of
// Frenkel and Smit, Understanding Molecular Simulation, Appendix E
// Note by Dustin okt2019: This implementation describes
// two Nose-Hoover Chains using Trotter Factorization of the Liouville Operator
double md::simulation::nose_hoover_thermostat_some_atoms(std::vector<size_t> active_atoms)
{
  double tempscale(0.0);
  int freedom_some = 3U * active_atoms.size();
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

double md::simulation::nose_hoover_with_arbitrary_chain_length(std::vector<size_t> active_atoms, std::size_t chainlength)
{
  double factor = 1.;
  md::nose_hoover_arbitrary_length& nht2 = this->nht2;
  return factor;
}

double md::simulation::tempcontrol(bool thermostat, bool half)
{
  std::size_t const N = this->coordobj.size();  // total number of atoms
  if (Config::get().md.set_active_center == 1)
    updateEkin(inner_atoms);
  else if (Config::get().coords.fixed.size() != 0)
    updateEkin(movable_atoms);
  else
    updateEkin(range(N));
  const double conversion_factor(2.0 / (freedom * md::R));     // factor for calculation of temperature from kinetic energy  
  double temp_after_scaling = 0.;
  const double instantaneous_temp_after_last_scaling = this->instantaneous_temp;
  double instantaneous_temp_before_scaling = this->E_kin * conversion_factor;
  double scaling_factor = 1.;

  if (Config::get().general.verbosity >= 5 && half)
  {
    std::cout << "Applying Nose-Hoover Thermostat halfstep.\n";
  }
  else if (Config::get().general.verbosity >= 5)
  {
    std::cout << "Applying Nose-Hoover Thermostat fullstep.\n";
  }
  size_t dof = freedom;
  if (Config::get().md.set_active_center == 1)  // if biased potential
  {
    dof = 3u * inner_atoms.size();
    if (Config::get().periodics.periodic == true)
      dof -= 3;
    else
      dof -= 6;
  }
  const double T_factor = (2.0 / (dof * md::R));
  if (Config::get().md.set_active_center == 1)  // if biased potential
  {
    double factor = 1.0;
    if (thermostat)
      factor = nose_hoover_thermostat_some_atoms(inner_atoms);     // calculate temperature scaling factor
    else
      factor = std::sqrt(desired_temp / instantaneous_temp_before_scaling);
    scaling_factor = factor;
  }
  else if (Config::get().coords.fixed.size() != 0)
  {
    updateEkin(movable_atoms);
    double factor = 1.0;
    if (thermostat)
      factor = nose_hoover_thermostat_some_atoms(movable_atoms);     // calculate temperature scaling factor
    else
      factor = std::sqrt(desired_temp / instantaneous_temp_before_scaling);
    scaling_factor = factor;
  }
  else if (thermostat)
  {
    const double factor = nose_hoover_thermostat();
    scaling_factor = factor;
  }
  //else if ("nose-hoover-chains")
  else
  {
    const double factor = std::sqrt(desired_temp / instantaneous_temp_before_scaling);
    scaling_factor = factor;
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
    std::cout << "Half step: desired Temp: " << desired_temp << "K. current temp: " << instantaneous_temp_after_last_scaling << "K. Velocity Scaling Factor: " << scaling_factor << "\n";
  }
  else if (Config::get().general.verbosity > 3)
  {
    std::cout << "Full step: desired Temp: " << desired_temp << "K. current temp: " << instantaneous_temp_after_last_scaling << "K. Velocity Scaling Factor: " << scaling_factor << "\n";
  }
  return temp_after_scaling;
}
