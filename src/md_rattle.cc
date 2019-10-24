#include "md.h"
#pragma once
// check if the two atom of a rattlepair are bonded with each other
void md::simulation::check_rattlepair_for_bond(config::md_conf::config_rattle::rattle_constraint_bond& rctemp)
{
  bool is_a_bond = false;
  for (auto bonding_partner : coordobj.atoms(rctemp.a).bonds())
  {
    if (bonding_partner == rctemp.b)
    {
      is_a_bond = true;
      break;
    }
  }
  if (is_a_bond == false)
  {
    throw std::runtime_error("No bond between " + std::to_string(rctemp.a + 1) + " and " + std::to_string(rctemp.b + 1) + ". It doesn't make sense to use a bonding parameter.");
  }
}

// set up constraints for H-X bonds if requested
// ideal bond lengths are taken from the foce field parameter file
// specified by RATTpar in the INPUTFILE
void md::simulation::rattlesetup(void)
{
  config::md_conf::config_rattle::rattle_constraint_bond rctemp;  // temporary rattlebond

  // temp vars and vectors;
  struct trat {
    std::size_t ia, ib;
    double ideal;
  };
  struct ratatoms {
    std::size_t ia, ga;
  };
  std::vector<trat> temprat;
  std::vector<ratatoms> ratoms;
  std::size_t tia = std::size_t(), tib = std::size_t();

  if (Config::get().md.rattle.use_paramfile)    //open rattle par file and look for atomtypes and distances
  {
    std::ifstream ifs;
    ifs.open(Config::get().general.paramFilename.c_str());
    if (ifs.good()) std::cout << "Opened file for RATTLE parameters successfully." << std::endl;
    if (!ifs.good()) {
      std::cout << "Couldn't open file for RATTLE parameters. Check your input" << std::endl;
      throw;
    }

    // more temporary variables
    ratatoms RTA;
    trat TRAT;
    char buffer[150];
    std::string bufferc;
    std::size_t found;
    std::vector < std::string > tokens;

    // read parameter file and extract ideal bond lengths
    while (!ifs.eof())
    {
      ifs.getline(buffer, 150);
      bufferc = buffer;
      found = bufferc.find("bond");
      if (found != std::string::npos)
      {
        tokens.clear();
        std::istringstream iss(bufferc);
        std::copy(std::istream_iterator <std::string>(iss), std::istream_iterator <std::string>(), std::back_inserter <std::vector <std::string >>(tokens));
        TRAT.ia = atoi(tokens[1].c_str());
        TRAT.ib = atoi(tokens[2].c_str());
        TRAT.ideal = atof(tokens[4].c_str());
        temprat.push_back(TRAT);
      }
    }// end of file check
    ifs.close();
    ifs.open(Config::get().general.paramFilename.c_str());
    // read parameter file and get atoms
    while (!ifs.eof())
    {
      ifs.getline(buffer, 150);
      bufferc = buffer;
      found = bufferc.find("atom");
      if (found != std::string::npos)
      {
        tokens.clear();
        std::istringstream iss(bufferc);
        std::copy(std::istream_iterator <std::string>(iss), std::istream_iterator <std::string>(), std::back_inserter <std::vector <std::string >>(tokens));
        RTA.ia = atoi(tokens[1].c_str());
        RTA.ga = atoi(tokens[2].c_str());
        ratoms.push_back(RTA);
      }
    }
    ifs.close();
  }

  // Generate vector with bonds which are to be constraint
  const std::size_t N = coordobj.size();

  if (Config::get().md.rattle.all == true)  // all H-bonds are constraint
  {
    for (std::size_t i = 0; i < N; i++) //loop over all atoms
    {
      if (coordobj.atoms(i).number() == 1) //check if atom is hydrogen
      {
        rctemp.a = i;                            // hydrogen atom
        rctemp.b = coordobj.atoms(i).bonds(0);   // bonding partner of hydrogen atom
        //loop over param vector
        for (unsigned j = 0; j < ratoms.size(); j++)
        {
          if (ratoms[j].ia == coordobj.atoms(rctemp.a).energy_type()) tia = ratoms[j].ga;
          if (ratoms[j].ia == coordobj.atoms(rctemp.b).energy_type()) tib = ratoms[j].ga;
        }
        for (unsigned k = 0; k < temprat.size(); k++)
        {
          // found matching parameters -> get ideal bond distance
          if ((tia == temprat[k].ia || tia == temprat[k].ib) && (tib == temprat[k].ia || tib == temprat[k].ib))
          {
            rctemp.len = temprat[k].ideal;
          }
        }
        rattle_bonds.push_back(rctemp);
      }
    }
  }
  else   // if MDrattle = 2 i.e. only specified distances are constrained
  {
    if (Config::get().md.rattle.use_paramfile == false && Config::get().md.rattle.dists.size() != Config::get().md.rattle.specified_rattle.size())
    {
      throw std::runtime_error("Wrong number of atom distances given for rattlepairs!");
    }

    for (auto i = 0u; i < Config::get().md.rattle.specified_rattle.size(); ++i)
    {
      rctemp = Config::get().md.rattle.specified_rattle[i];

      if (Config::get().md.rattle.use_paramfile)   // if distances are to be taken from parameterfile
      {
        // check if atoms share a bond
        check_rattlepair_for_bond(rctemp);

        //loop over param vector
        for (unsigned j = 0; j < ratoms.size(); j++)
        {
          if (ratoms[j].ia == coordobj.atoms(rctemp.a).energy_type()) tia = ratoms[j].ga;
          if (ratoms[j].ia == coordobj.atoms(rctemp.b).energy_type()) tib = ratoms[j].ga;
        }
        for (unsigned k = 0; k < temprat.size(); k++)
        {
          // found matching parameters -> get ideal bond distance
          if ((tia == temprat[k].ia || tia == temprat[k].ib) && (tib == temprat[k].ia || tib == temprat[k].ib))
          {
            rctemp.len = temprat[k].ideal;
            break;
          }
        }
      }

      else   // get rattledist from inputfile
      {
        rctemp.len = Config::get().md.rattle.dists[i];
      }
      rattle_bonds.push_back(rctemp);
    }
  }
}

//First part of the RATTLE algorithm to constrain H-X bonds ( half step)
void md::simulation::rattle_pre(void)
{
  std::size_t niter(0U);
  bool done = false;
  // main loop till convergence is reached
  do
  {
    niter += 1;
    done = true;
    // loop over all constraint bonds
    for (auto const& rattlebond : rattle_bonds)
    {
      // get bond distance
      coords::Cartesian_Point const d(coordobj.xyz(rattlebond.b) - coordobj.xyz(rattlebond.a));
      // difference of distance square to optimal distance square
      double const delta = rattlebond.len * rattlebond.len - dot(d, d);
      // if delta is greater than tolerance -> rattle
      if (std::abs(delta) > Config::get().md.rattle.tolerance)
      {
        done = false;
        coords::Cartesian_Point const d_old(P_old[rattlebond.b] - P_old[rattlebond.a]);
        double inv_ma = 1.0 / M[rattlebond.a];
        double inv_mb = 1.0 / M[rattlebond.b];
        //calculate lagrange multiplier
        coords::Cartesian_Point rattle(delta / (2.0 * (inv_ma + inv_mb) * dot(d, d_old)));
        rattle *= d_old;
        // update half step positions
        coordobj.move_atom_by(rattlebond.a, -(rattle * inv_ma));
        coordobj.move_atom_by(rattlebond.b, rattle * inv_mb);
        inv_ma /= Config::get().md.timeStep;
        inv_mb /= Config::get().md.timeStep;
        // update half step velocities
        V[rattlebond.a] -= rattle * inv_ma;
        V[rattlebond.b] += rattle * inv_mb;
      }
    }
  } while (niter < Config::get().md.rattle.num_iter && done == false);
}
//second part of RATTLE to constrain H.X bonds (full step)
void md::simulation::rattle_post(void)
{
  // initialize some local variables neede for the internal virial tensor
  std::size_t niter(0U);
  bool done = false;
  double tfact = 2.0 / (dt * 1.0e-3 * convert);
  double xvir, yvir, zvir;
  double vxx, vyx, vzx, vyy, vzy, vzz;
  std::array<std::array<double, 3>, 3> part_v;
  part_v[0][0] = part_v[0][1] = part_v[0][2] = part_v[1][0] = part_v[1][1] = part_v[1][2] = part_v[2][0] = part_v[2][1] = part_v[2][2] = 0.0;
  // main loop till convergence is reached
  do
  {
    niter += 1;
    done = true;
    // loop over all constraint bonds
    for (auto const& rattlebond : rattle_bonds)
    {
      // get bond distance
      coords::Cartesian_Point const d(coordobj.xyz(rattlebond.b) - coordobj.xyz(rattlebond.a));
      double inv_ma = 1.0 / M[rattlebond.a];
      double inv_mb = 1.0 / M[rattlebond.b];
      // calculate lagrange multiplier
      double const lagrange = -dot(V[rattlebond.b] - V[rattlebond.a], d) / (rattlebond.len * rattlebond.len * (inv_ma + inv_mb));
      if (std::fabs(lagrange) > Config::get().md.rattle.tolerance)
      {
        done = false;
        coords::Cartesian_Point rattle(d * lagrange);
        // update full step velocities
        V[rattlebond.a] -= rattle * inv_ma;
        V[rattlebond.b] += rattle * inv_mb;
        // compute RATTLE contributions to the internal virial tensor
        xvir = rattle.x() * tfact;
        yvir = rattle.y() * tfact;
        zvir = rattle.z() * tfact;
        vxx = d.x() * xvir;
        vyx = d.y() * xvir;
        vzx = d.z() * xvir;
        vyy = d.y() * yvir;
        vzy = d.z() * yvir;
        vzz = d.z() * zvir;
        part_v[0][0] -= vxx;
        part_v[1][0] -= vyx;
        part_v[2][0] -= vzx;
        part_v[0][1] -= vyx;
        part_v[1][1] -= vyy;
        part_v[2][1] -= vzy;
        part_v[0][2] -= vzx;
        part_v[1][2] -= vzy;
        part_v[2][2] -= vzz;
      }
    }
  } while (niter < Config::get().md.rattle.num_iter && done == false);
  // Add RATTLE contributions to the internal virial tensor
  coordobj.add_to_virial(part_v);
}

