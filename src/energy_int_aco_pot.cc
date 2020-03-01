/**
This file contains the calculation of energy and gradients for amber, oplsaa and charmm forcefield.
*/

#include <cmath>
#include <stddef.h>
#include <stdexcept>
#include <cstdlib>
#include "energy_int_aco.h"
#include "configuration.h"
#include "Scon/scon_utility.h"
#include "Scon/scon_mathmatrix.h"

/****************************************
*                                       *
*                                       *
*               Generics                *
*                                       *
*                                       *
*****************************************/

coords::float_type energy::interfaces::aco::aco_ff::e(void)
{
  pre();      // set everything to zero
  calc<0>();  // calc partial energies
  post();     // sum up partial energies
  return energy;
}

coords::float_type energy::interfaces::aco::aco_ff::g(void)
{
  pre();      // set everything to zero
  calc<1>();  // calc partial energies and gradients
  post();     // sum up partial energies
  return energy;
}

coords::float_type energy::interfaces::aco::aco_ff::h(void)
{
  constexpr double hesscut = 0.;
  // This is implemented in accordance to Ponder's tinker implementation, Version 8.7.2, file 'hessian.f'

  // zero out total number of indexed Hessian elements
  const std::size_t numberOfAtoms = this->coords->atoms().size();
  typedef scon::mathmatrix<coords::float_type> Matrix_Class;

  std::vector<std::size_t> hindex(numberOfAtoms,0u);
  std::vector<coords::float_type> h(numberOfAtoms,0.);
  Matrix_Class hinit(3u,numberOfAtoms,1.);
  Matrix_Class hstop(3u, numberOfAtoms, 0.);
  Matrix_Class hdiag(hstop);
  std::deque<bool> keep(numberOfAtoms,false);

  // Now we would compute induced dipoles at polaraizeable atoms for AMOEBE ff
  // skipping this...

  // Calculate reduced coordinates: I don't think we need this, skipping....

  Matrix_Class hessx(3u,numberOfAtoms,0.);
  Matrix_Class hessy(hessx),hessz(hessx);
  const double cutoff = 0.;

  pre();      // set everything to zero
  calc<2>();  // calc partial energies and gradients
  post();     // sum up partial energies

  // set diagonal elements:
  std::size_t nhess = 0u;
  for (std::size_t i = 0u; i < numberOfAtoms; ++i)
  {
    for (std::size_t j = 0u; j < 3u; ++j)
    {
      hdiag(j, i) += hessx(j, i);
      hdiag(j, i) += hessy(j, i);
      hdiag(j, i) += hessz(j, i);
    }
    // search each 3x3 block to see which blocks will be kept
    for (std::size_t k = i + 1u; k < numberOfAtoms; ++k)
    {
      if (!coords->atoms(k).fixed())
      {
        double hmax = 0.;
        for (std::size_t j = 0u; j < 3u; ++j)
        {
          hmax = std::max(hmax,hessx(j,k));
          hmax = std::max(hmax, hessy(j, k));
          hmax = std::max(hmax, hessz(j, k));
        }
        if (hmax >= hesscut)
          keep.at(k) = true;
      }
    }
    hinit(0u,i) = nhess + 1;
    for (std::size_t j = 2u; j <= 3u; ++j)
    {
      nhess += 1u;
      hindex[nhess] = 3*(i+1) + j - 3;
      h[nhess] = hessx(j-1u,i);
    }
    for (std::size_t k = i + 1u; k < numberOfAtoms; ++k)
    {
      if (keep.at(k))
      {
        for (std::size_t j = 1u; j <= 3u; ++j)
        {
          nhess += 1u;
          hindex[nhess] = 3 * (k + 1) + j - 3;
          h[nhess] = hessx(j - 1u, k);
        }
      }
    }
    hstop(0u,i) = nhess;
    hinit(1u,i) = nhess + 1u;
    nhess += 1u;
    hindex[nhess] = 3 * (i + 1u);
    h[nhess] = hessy(2u,i);
    for (std::size_t k = i + 1u; k < numberOfAtoms; ++k)
    {
      if (keep.at(k))
      {
        for (std::size_t j = 1u; j <= 3u; ++j)
        {
          nhess += 1u;
          hindex[nhess] = 3 * (k + 1) + j - 3;
          h[nhess] = hessy(j - 1u, k);
        }
      }
    }
    hstop(1u,i) = nhess;
    hinit(2u,i) = nhess + 1u;
    for (std::size_t k = i + 1u; k < numberOfAtoms; ++k)
    {
      if (keep.at(k))
      {
        for (std::size_t j = 1u; j <= 3u; ++j)
        {
          nhess += 1u;
          hindex[nhess] = 3 * (k + 1) + j - 3;
          h[nhess] = hessz(j - 1u, k);
        }
      }
    }
    // done, return h, hinit, hstop, hindex and hdiag
    // To-DO: Turn h into proper matrix
  }
  

  throw std::runtime_error("aco_ff doesn't provide a hessian matrix.");
}

coords::float_type energy::interfaces::aco::aco_ff::o(void)
{
  throw std::runtime_error("aco_ff doesn't provide any optimization routines.");
}

template<size_t DERIV>
void energy::interfaces::aco::aco_ff::calc(void)
{
#pragma omp parallel sections
  {
#pragma omp section
    part_energy[types::BOND] = f_12<DERIV>();
#pragma omp section
    part_energy[types::ANGLE] = f_13_a<DERIV>();
#pragma omp section
    part_energy[types::UREY] = f_13_u<DERIV>();
#pragma omp section
    part_energy[types::TORSION] = f_14<DERIV>();
#pragma omp section
    part_energy[types::IMPTORSION] = f_it<DERIV>();
#pragma omp section
    part_energy[types::IMPROPER] = f_imp<DERIV>();
  }

  // fill part_energy[CHARGE], part_energy[VDW] and part_grad[VDW], part_grad[CHARGE]
  if (cparams.radiustype() == ::tinker::parameter::radius_types::R_MIN)
  {
    g_nb< ::tinker::parameter::radius_types::R_MIN>();
  }
  else
  {
    g_nb< ::tinker::parameter::radius_types::SIGMA>();
  }

  if (Config::get().energy.qmmm.mm_charges.size() != 0)
  {
    calc_ext_charges_interaction(DERIV);   // adds to part_energy[EXTERNAL_CHARGES] and part_grad[EXTERNAL_CHARGES]
  }
}


/****************************************
*                                       *
*                                       *
*              Potentials               *
*                                       *
*                                       *
*****************************************/



/********************************
*                                *
*                                *
*  Bonding  Potential            *
*  Energy/Gradients/Hessians     *
*                                *
*                                *
*********************************/

namespace energy
{
  namespace interfaces
  {
    namespace aco
    {

      //! Bonding Energy
      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_12<0>(void)
      {
        using std::abs;
        coords::float_type E(0.0);
        for (auto const& bond : refined.bonds())
        {
          auto const bv(coords->xyz(bond.atoms[0]) - coords->xyz(bond.atoms[1]));
          auto const d = len(bv);
          auto const r = d - bond.ideal;
          E += bond.force * r * r;
          if (abs(r) > 0.5)
          {
            if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because length of bond " << bond << " is " << d << " but should be " << bond.ideal << "\n";
            integrity = false;
          }
        }
        return E;
      }

      //! Bonding Energy+Gradients
      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_12<1>(void)
      {
        using std::abs;
        using scon::len;
        coords::float_type E(0.0);
        for (auto const& bond : refined.bonds())
        {
          auto const bv(coords->xyz(bond.atoms[0]) - coords->xyz(bond.atoms[1])); // r_ij (i=1, j=2)
          auto const d = len(bv);
          auto const r = d - bond.ideal;
          auto dE = bond.force * r;
          E += dE * r;  // kcal/mol
          dE *= 2;  // kcal/(mol*Angstrom)  gradient without direction
          if (abs(d) > 0.0)
          {
            if (abs(r) > 0.5)
            {
              if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because length of bond " << bond << " is " << d << " but should be " << bond.ideal << "\n";
              integrity = false;
            }
            dE /= d;  // kcal/(mol*A^2)   gradient divided by distance because later it is multiplied with it again
            auto const gv = bv * dE;   // "force" on atom i due to atom j (kcal/(mol*A)), gradient with direction
            part_grad[BOND][bond.atoms[0]] += gv;   //(kcal/(mol*A))
            part_grad[BOND][bond.atoms[1]] -= gv;
            //increment internal virial tensor (no factor 1/2 because atoms i and j have the same contribution)
            auto const vxx = bv.x() * gv.x();
            auto const vyx = bv.y() * gv.x();
            auto const vzx = bv.z() * gv.x();
            auto const vyy = bv.y() * gv.y();
            auto const vzy = bv.z() * gv.y();
            auto const vzz = bv.z() * gv.z();
            part_virial[BOND][0][0] += vxx;
            part_virial[BOND][1][0] += vyx;
            part_virial[BOND][2][0] += vzx;
            part_virial[BOND][0][1] += vyx;
            part_virial[BOND][1][1] += vyy;
            part_virial[BOND][2][1] += vzy;
            part_virial[BOND][0][2] += vzx;
            part_virial[BOND][1][2] += vzy;
            part_virial[BOND][2][2] += vzz;
            /*std::cout << "Bond virial" << std::endl;
            std::cout << part_virial[BOND][0][0] << "   " << part_virial[BOND][1][0] << "   " << part_virial[BOND][2][0] << std::endl;
            std::cout << part_virial[BOND][0][1] << "   " << part_virial[BOND][1][1] << "   " << part_virial[BOND][2][1] << std::endl;
            std::cout << part_virial[BOND][0][2] << "   " << part_virial[BOND][1][2] << "   " << part_virial[BOND][2][2] << std::endl;*/
          }
          else
          {
            if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because of bond " << bond << "\n";
            integrity = false;
          }
        }
        return E;
      }

      //! Bonding Energy+Gradients+Hessians
      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_12<2>(void)
      {
        using std::abs;
        using scon::len;
        coords::float_type E(0.0);
        for (auto const& bond : refined.bonds())
        {
          auto const bv(coords->xyz(bond.atoms[0]) - coords->xyz(bond.atoms[1])); // r_ij (i=1, j=2)
          auto const d = len(bv);
          auto const r = d - bond.ideal;
          auto dE = bond.force * r;
          E += dE * r;  // kcal/mol
          dE *= 2.;  // kcal/(mol*Angstrom)  gradient without direction
          auto ddE = 2. * bond.force * r; // We left cbnd and qbnd terms
          if (abs(d) > 0.0)
          {
            if (abs(r) > 0.5)
            {
              if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because length of bond " << bond << " is " << d << " but should be " << bond.ideal << "\n";
              integrity = false;
            }
            dE /= d;  // kcal/(mol*A^2)   gradient divided by distance because later it is multiplied with it again
            auto const gv = bv * dE;   // "force" on atom i due to atom j (kcal/(mol*A)), gradient with direction
            part_grad[BOND][bond.atoms[0]] += gv;   //(kcal/(mol*A))
            part_grad[BOND][bond.atoms[1]] -= gv;
            //increment internal virial tensor (no factor 1/2 because atoms i and j have the same contribution)
            auto const vxx = bv.x() * gv.x();
            auto const vyx = bv.y() * gv.x();
            auto const vzx = bv.z() * gv.x();
            auto const vyy = bv.y() * gv.y();
            auto const vzy = bv.z() * gv.y();
            auto const vzz = bv.z() * gv.z();
            part_virial[BOND][0][0] += vxx;
            part_virial[BOND][1][0] += vyx;
            part_virial[BOND][2][0] += vzx;
            part_virial[BOND][0][1] += vyx;
            part_virial[BOND][1][1] += vyy;
            part_virial[BOND][2][1] += vzy;
            part_virial[BOND][0][2] += vzx;
            part_virial[BOND][1][2] += vzy;
            part_virial[BOND][2][2] += vzz;
            //
            // set the chain rule terms for the Hessian elements
            const double term = (ddE -dE) / (d*d);
            const double xab = (coords->xyz(bond.atoms[0u]).x() - coords->xyz(bond.atoms[1u]).x());
            const double yab = (coords->xyz(bond.atoms[0u]).y() - coords->xyz(bond.atoms[1u]).y());
            const double zab = (coords->xyz(bond.atoms[0u]).z() - coords->xyz(bond.atoms[1u]).z());
            const double termx = term * xab;
            const double termy = term * yab;
            const double termz = term * zab;
            std::array<std::array<coords::float_type,3u>,3u> d2e = std::array<std::array<coords::float_type, 3u>, 3u>();
            d2e.at(0).at(0) = termx * xab + dE;
            d2e.at(0).at(1) = termx * yab;
            d2e.at(0).at(2) = termx * zab;
            d2e.at(1).at(0) = d2e.at(0).at(0);
            d2e.at(1).at(1) = termy * yab + dE;
            d2e.at(1).at(2) = termy * zab;
            d2e.at(2).at(0) = d2e.at(0).at(2);
            d2e.at(2).at(1) = d2e.at(1).at(2);
            d2e.at(2).at(2) = termz * zab + dE;
            part_hessian[BOND].at(3 * bond.atoms[0u]).at( 3 * bond.atoms[0u])          += d2e.at(0).at(0); // dx dx same atom
            part_hessian[BOND].at(3 * bond.atoms[0u]).at(3 * bond.atoms[0u] + 1u)      += d2e.at(0).at(1); // dx dy same atom
            part_hessian[BOND].at(3 * bond.atoms[0u]).at(3 * bond.atoms[0u] + 2u)      += d2e.at(0).at(2); // dx dz same atom
            part_hessian[BOND].at(3 * bond.atoms[0u] + 1u).at(3 * bond.atoms[0u])      += d2e.at(1).at(0); // dy dx same atom
            part_hessian[BOND].at(3 * bond.atoms[0u] + 1u).at(3 * bond.atoms[0u] + 1u) += d2e.at(1).at(1); // dy dy same atom
            part_hessian[BOND].at(3 * bond.atoms[0u] + 1u).at(3 * bond.atoms[0u] + 2u) += d2e.at(1).at(2); // dy dz same atom
            part_hessian[BOND].at(3 * bond.atoms[0u] + 2u).at(3 * bond.atoms[0u])      += d2e.at(2).at(0); // dy dx same atom
            part_hessian[BOND].at(3 * bond.atoms[0u] + 2u).at(3 * bond.atoms[0u] + 1u) += d2e.at(2).at(1); // dy dy same atom
            part_hessian[BOND].at(3 * bond.atoms[0u] + 2u).at(3 * bond.atoms[0u] + 2u) += d2e.at(2).at(2); // dy dz same atom
            part_hessian[BOND].at(3 * bond.atoms[1u]).at(3 * bond.atoms[1u])           -= d2e.at(0).at(0); // dx dx same atom
            part_hessian[BOND].at(3 * bond.atoms[1u]).at(3 * bond.atoms[1u] + 1u)      -= d2e.at(0).at(1); // dx dy same atom
            part_hessian[BOND].at(3 * bond.atoms[1u]).at(3 * bond.atoms[1u] + 2u)      -= d2e.at(0).at(2); // dx dz same atom
            part_hessian[BOND].at(3 * bond.atoms[1u] + 1u).at(3 * bond.atoms[1u])      -= d2e.at(1).at(0); // dy dx same atom
            part_hessian[BOND].at(3 * bond.atoms[1u] + 1u).at(3 * bond.atoms[1u] + 1u) -= d2e.at(1).at(1); // dy dy same atom
            part_hessian[BOND].at(3 * bond.atoms[1u] + 1u).at(3 * bond.atoms[1u] + 2u) -= d2e.at(1).at(2); // dy dz same atom
            part_hessian[BOND].at(3 * bond.atoms[1u] + 2u).at(3 * bond.atoms[1u])      -= d2e.at(2).at(0); // dy dx same atom
            part_hessian[BOND].at(3 * bond.atoms[1u] + 2u).at(3 * bond.atoms[1u] + 1u) -= d2e.at(2).at(1); // dy dy same atom
            part_hessian[BOND].at(3 * bond.atoms[1u] + 2u).at(3 * bond.atoms[1u] + 2u) -= d2e.at(2).at(2); // dy dz same atom
          }
          else
          {
            if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because of bond " << bond << "\n";
            integrity = false;
          }
        }
        return E;
      }

      /*********************************
      *                                *
      *                                *
      *  Angle  Potential              *
      *  Energy/Gradients/Hessians     *
      *                                *
      *                                *
      *********************************/

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_13_a<0>(void)
      {
        using std::abs;
        coords::float_type E(0.0);
        for (auto const& angle : refined.angles())
        {
          auto const
            av1(coords->xyz(angle.atoms[0]) - coords->xyz(angle.atoms[1])),
            av2(coords->xyz(angle.atoms[2]) - coords->xyz(angle.atoms[1]));
          auto const d(scon::angle(av1, av2).degrees() - angle.ideal);
          auto const r(d * SCON_PI180);
          E += angle.force * r * r;
          if (abs(d) > 30.0)
          {
            if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because angle " << angle << " is " << scon::angle(av1, av2).degrees() << " but should be " << angle.ideal << "\n";
            integrity = false;
          }
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_13_a<1>(void)
      {
        using scon::angle;
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const& angle : refined.angles())
        {
          auto const
            av1(coords->xyz(angle.atoms[0]) - coords->xyz(angle.atoms[1])),
            av2(coords->xyz(angle.atoms[2]) - coords->xyz(angle.atoms[1]));
          auto const d = scon::angle(av1, av2).degrees() - angle.ideal;
          auto const r = d * SCON_PI180;
          /*std::cout << angle.atoms[0] << ", " << angle.atoms[1] << ", " << angle.atoms[2] << "\n";
          std::cout << "Delta deg: " << d << ", delta rad: " << r << ", ideal: " << angle.ideal << ", force = " << angle.force << '\n';
          std::cout << scon::angle(av1, av2).degrees() << "; " << angle.force*r*r << "\n";*/
          auto dE = angle.force * r;
          E += dE * r;
          //std::cout << "A " << angle.atoms[0] << "-" << angle.atoms[1] << "-" << angle.atoms[2] << ": ";

          coords::Cartesian_Point const cv(cross(av1, av2));
          coords::float_type const cvl(len(cv));
          if (abs(cvl) > 0.0)
          {
            if (abs(d) > 30.0)
            {
              if (Config::get().general.verbosity > 3) std::cout << "Integrity broke because angle " << angle << " is " << scon::angle(av1, av2).degrees() << " but should be " << angle.ideal << "\n";
              integrity = false;
            }
            dE *= 2.0 / cvl;
            //std::cout << dE << ", " << cvl << '\n';
            coords::Cartesian_Point const gv1(cross(av1, cv) * (dE / dot(av1, av1)));
            coords::Cartesian_Point const gv2(cross(cv, av2) * (dE / dot(av2, av2)));
            //std::cout << gv1 << ", " << gv2 << ", " << -(gv1 + gv2) << '\n';
            //std::cout << part_grad[ANGLE][angle.atoms[0]] << ", " << part_grad[ANGLE][angle.atoms[1]] << ", " << part_grad[ANGLE][angle.atoms[2]] << '\n';
            part_grad[ANGLE][angle.atoms[0]] += gv1;
            part_grad[ANGLE][angle.atoms[2]] += gv2;
            part_grad[ANGLE][angle.atoms[1]] += -(gv1 + gv2);
            //std::cout << part_grad[ANGLE][angle.atoms[0]] << ", " << part_grad[ANGLE][angle.atoms[1]] << ", " << part_grad[ANGLE][angle.atoms[2]] << '\n';
            //increment internal virial tensor
            coords::float_type const vxx = av1.x() * gv1.x() + av2.x() * gv2.x();
            coords::float_type const vyx = av1.y() * gv1.x() + av2.y() * gv2.x();
            coords::float_type const vzx = av1.z() * gv1.x() + av2.z() * gv2.x();
            coords::float_type const vyy = av1.y() * gv1.y() + av2.y() * gv2.y();
            coords::float_type const vzy = av1.z() * gv1.y() + av2.z() * gv2.y();
            coords::float_type const vzz = av1.z() * gv1.z() + av2.z() * gv2.z();
            part_virial[ANGLE][0][0] += vxx;
            part_virial[ANGLE][1][0] += vyx;
            part_virial[ANGLE][2][0] += vzx;
            part_virial[ANGLE][0][1] += vyx;
            part_virial[ANGLE][1][1] += vyy;
            part_virial[ANGLE][2][1] += vzy;
            part_virial[ANGLE][0][2] += vzx;
            part_virial[ANGLE][1][2] += vzy;
            part_virial[ANGLE][2][2] += vzz;
            /*std::cout << angle.atoms[0]+1 << "   " << angle.atoms[1]+1 << "   " << angle.atoms[2]+1 << std::endl;
            std::cout << "Angle virial" << std::endl;*/
            /* std::cout << part_virial[ANGLE][0][0] << "   " << part_virial[ANGLE][1][0] << "   " << part_virial[ANGLE][2][0] << std::endl;
             std::cout << part_virial[ANGLE][0][1] << "   " << part_virial[ANGLE][1][1] << "   " << part_virial[ANGLE][2][1] << std::endl;
             std::cout << part_virial[ANGLE][0][2] << "   " << part_virial[ANGLE][1][2] << "   " << part_virial[ANGLE][2][2] << std::endl;*/
             //std::cout << vxx << std::endl;
             //std::cout << vyx << std::endl; 
             //std::cout << vzx << std::endl; 
             //std::cout << vyy << std::endl; 
             //std::cout << vzy << std::endl; 
             //std::cout << vzz << std::endl; 
          }
          else
          {
            if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because of angle " << angle << "\n";
            integrity = false;
          }
        }
        /* std::cout << "Virial sum" << std::endl;
         std::cout << part_virial[ANGLE][0][0] + part_virial[BOND][0][0] << "   " << part_virial[ANGLE][1][0] + part_virial[BOND][1][0] << "   " << part_virial[ANGLE][2][0] + part_virial[BOND][2][0] << std::endl;
         std::cout << part_virial[ANGLE][0][1] + part_virial[BOND][0][1] << "   " << part_virial[ANGLE][1][1] + part_virial[BOND][1][1] << "   " << part_virial[ANGLE][2][1] + part_virial[BOND][2][1] << std::endl;
         std::cout << part_virial[ANGLE][0][2] + part_virial[BOND][0][2] << "   " << part_virial[ANGLE][1][2] + part_virial[BOND][1][2] << "   " << part_virial[ANGLE][2][2] + part_virial[BOND][2][2] << std::endl;*/
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_13_a<2>(void)
      {
        using scon::angle;
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const& angle : refined.angles())
        {
          auto const
            av1(coords->xyz(angle.atoms[0]) - coords->xyz(angle.atoms[1])),
            av2(coords->xyz(angle.atoms[2]) - coords->xyz(angle.atoms[1]));
          auto const d = scon::angle(av1, av2).degrees() - angle.ideal;
          auto const r = d * SCON_PI180; // This might be wrong, in the tinker code this is multiplied by 57.xx herefore SCON_180PI
          auto dE = angle.force * r;
          E += dE * r;
          coords::Cartesian_Point const cv(cross(av1, av2));
          coords::float_type const cvl(len(cv));
          if (abs(cvl) > 0.0)
          {
            if (abs(d) > 30.0)
            {
              if (Config::get().general.verbosity > 3) std::cout << "Integrity broke because angle " << angle << " is " << scon::angle(av1, av2).degrees() << " but should be " << angle.ideal << "\n";
              integrity = false;
            }
            dE *= 2.0 / cvl;
            coords::Cartesian_Point const gv1(cross(av1, cv) * (dE / dot(av1, av1)));
            coords::Cartesian_Point const gv2(cross(cv, av2) * (dE / dot(av2, av2)));
            part_grad[ANGLE][angle.atoms[0]] += gv1;
            part_grad[ANGLE][angle.atoms[2]] += gv2;
            part_grad[ANGLE][angle.atoms[1]] += -(gv1 + gv2);
            coords::float_type const vxx = av1.x() * gv1.x() + av2.x() * gv2.x();
            coords::float_type const vyx = av1.y() * gv1.x() + av2.y() * gv2.x();
            coords::float_type const vzx = av1.z() * gv1.x() + av2.z() * gv2.x();
            coords::float_type const vyy = av1.y() * gv1.y() + av2.y() * gv2.y();
            coords::float_type const vzy = av1.z() * gv1.y() + av2.z() * gv2.y();
            coords::float_type const vzz = av1.z() * gv1.z() + av2.z() * gv2.z();
            part_virial[ANGLE][0][0] += vxx;
            part_virial[ANGLE][1][0] += vyx;
            part_virial[ANGLE][2][0] += vzx;
            part_virial[ANGLE][0][1] += vyx;
            part_virial[ANGLE][1][1] += vyy;
            part_virial[ANGLE][2][1] += vzy;
            part_virial[ANGLE][0][2] += vzx;
            part_virial[ANGLE][1][2] += vzy;
            part_virial[ANGLE][2][2] += vzz;
            //
            // Begin Hessian (tinker v. 8.7.2. file eangle2.f)
            //
            const double xab = (coords->xyz(angle.atoms[0]) - coords->xyz(angle.atoms[1])).x();
            const double yab = (coords->xyz(angle.atoms[0]) - coords->xyz(angle.atoms[1])).y();
            const double zab = (coords->xyz(angle.atoms[0]) - coords->xyz(angle.atoms[1])).z();
            const double xcb = (coords->xyz(angle.atoms[2]) - coords->xyz(angle.atoms[1])).x();
            const double ycb = (coords->xyz(angle.atoms[2]) - coords->xyz(angle.atoms[1])).y();
            const double zcb = (coords->xyz(angle.atoms[2]) - coords->xyz(angle.atoms[1])).z();
            const double rab2 = xab * xab + yab * yab + zab * zab;
            const double rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;
            const double xp = ycb * zab - zcb * yab;
            const double yp = zcb * xab - xcb * zab;
            const double zp = xcb * yab - ycb * xab;
            double rp = sqrt(xp * xp + yp * yp + zp * zp);
            const double dot = xab * xcb + yab * ycb + zab * zcb;
            double cosine = dot / sqrt(rab2 * rcb2);
            cosine = std::min(1.0, std::max(-1.0, cosine));
            const double deddt = angle.force * d * SCON_PI180;
            const double d2eddt2 = angle.force * SCON_PI180 * SCON_PI180;
            // construct an orthogonal direction for linear angles
            bool linear = false;
            if (rp < 0.0001)
            {
              linear = true;
              double xp(0.),yp(0.),zp(0.);
              if (xab != 0.0  && yab != 0.0) 
              {
                xp = -yab;
                yp = xab ;
                zp = 0.0 ;
              }
              else if (xab == 0.0  && yab == 0.0) 
              {
                xp = 1.0;
                yp = 0.0;
                zp = 0.0;
              }
              else if (xab != 0.0  && yab == 0.0) 
              {
                xp = 0.0;
                yp = 1.0;
                zp = 0.0;
              }
              else if (xab == 0.0  && yab != 0.0)
              {
                xp = 1.0;
                yp = 0.0;
                zp = 0.0;
              }
              rp = sqrt(xp * xp + yp * yp + zp * zp);
            }
            //first derivatives of bond angle with respect to coordinates
            const double terma = -1.0 / (rab2 * rp)            ;
            const double termc = 1.0 / (rcb2 * rp)             ;
            const double ddtdxia = terma * (yab * zp - zab * yp) ;
            const double ddtdyia = terma * (zab * xp - xab * zp) ;
            const double ddtdzia = terma * (xab * yp - yab * xp) ;
            const double ddtdxic = termc * (ycb * zp - zcb * yp) ;
            const double ddtdyic = termc * (zcb * xp - xcb * zp) ;
            const double ddtdzic = termc * (xcb * yp - ycb * xp) ;
            const double ddtdxib = -ddtdxia - ddtdxic            ;
            const double ddtdyib = -ddtdyia - ddtdyic            ;
            const double ddtdzib = -ddtdzia - ddtdzic            ;
            //abbreviations used in defining chain rule terms
            const double xrab = 2.0 * xab / rab2               ;
            const double yrab = 2.0 * yab / rab2               ;
            const double zrab = 2.0 * zab / rab2               ;
            const double xrcb = 2.0 * xcb / rcb2               ;
            const double yrcb = 2.0 * ycb / rcb2               ;
            const double zrcb = 2.0 * zcb / rcb2               ;
            const double rp2 = 1.0 / (rp * rp)                 ;
            const double xabp = (yab * zp - zab * yp) * rp2      ;
            const double yabp = (zab * xp - xab * zp) * rp2      ;
            const double zabp = (xab * yp - yab * xp) * rp2      ;
            const double xcbp = (ycb * zp - zcb * yp) * rp2      ;
            const double ycbp = (zcb * xp - xcb * zp) * rp2      ;
            const double zcbp = (xcb * yp - ycb * xp) * rp2      ;
            //chain rule terms for second derivative components
            const double dxiaxia = terma * (xab * xcb - dot) + ddtdxia * (xcbp - xrab);
            const double dxiayia = terma * (zp + yab * xcb) + ddtdxia * (ycbp - yrab) ;
            const double dxiazia = terma * (zab * xcb - yp) + ddtdxia * (zcbp - zrab) ;
            const double dyiayia = terma * (yab * ycb - dot) + ddtdyia * (ycbp - yrab);
            const double dyiazia = terma * (xp + zab * ycb) + ddtdyia * (zcbp - zrab) ;
            const double dziazia = terma * (zab * zcb - dot) + ddtdzia * (zcbp - zrab);
            const double dxicxic = termc * (dot - xab * xcb) - ddtdxic * (xabp + xrcb);
            const double dxicyic = termc * (zp - ycb * xab) - ddtdxic * (yabp + yrcb) ;
            const double dxiczic = -termc * (yp + zcb * xab) - ddtdxic * (zabp + zrcb);
            const double dyicyic = termc * (dot - yab * ycb) - ddtdyic * (yabp + yrcb);
            const double dyiczic = termc * (xp - zcb * yab) - ddtdyic * (zabp + zrcb) ;
            const double dziczic = termc * (dot - zab * zcb) - ddtdzic * (zabp + zrcb);
            const double dxiaxic = terma * (yab * yab + zab * zab) - ddtdxia * xabp   ;
            const double dxiayic = -terma * xab * yab - ddtdxia * yabp                ;
            const double dxiazic = -terma * xab * zab - ddtdxia * zabp                ;
            const double dyiaxic = -terma * xab * yab - ddtdyia * xabp                ;
            const double dyiayic = terma * (xab * xab + zab * zab) - ddtdyia * yabp   ;
            const double dyiazic = -terma * yab * zab - ddtdyia * zabp                ;
            const double dziaxic = -terma * xab * zab - ddtdzia * xabp                ;
            const double dziayic = -terma * yab * zab - ddtdzia * yabp                ;
            const double dziazic = terma * (xab * xab + yab * yab) - ddtdzia * zabp   ;
            // get some second derivative chain rule terms by difference
            const double dxibxia = -dxiaxia - dxiaxic;
            const double dxibyia = -dxiayia - dyiaxic;
            const double dxibzia = -dxiazia - dziaxic;
            const double dyibxia = -dxiayia - dxiayic;
            const double dyibyia = -dyiayia - dyiayic;
            const double dyibzia = -dyiazia - dziayic;
            const double dzibxia = -dxiazia - dxiazic;
            const double dzibyia = -dyiazia - dyiazic;
            const double dzibzia = -dziazia - dziazic;
            const double dxibxic = -dxicxic - dxiaxic;
            const double dxibyic = -dxicyic - dxiayic;
            const double dxibzic = -dxiczic - dxiazic;
            const double dyibxic = -dxicyic - dyiaxic;
            const double dyibyic = -dyicyic - dyiayic;
            const double dyibzic = -dyiczic - dyiazic;
            const double dzibxic = -dxiczic - dziaxic;
            const double dzibyic = -dyiczic - dziayic;
            const double dzibzic = -dziczic - dziazic;
            const double dxibxib = -dxibxia - dxibxic;
            const double dxibyib = -dxibyia - dxibyic;
            const double dxibzib = -dxibzia - dxibzic;
            const double dyibyib = -dyibyia - dyibyic;
            const double dyibzib = -dyibzia - dyibzic;
            const double dzibzib = -dzibzia - dzibzic;
            //
            part_hessian[ANGLE].at(3 * angle.atoms[0u]     ).at(3 * angle.atoms[0u]     ) += deddt * dxiaxia + d2eddt2 * ddtdxia * ddtdxia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u]     ).at(3 * angle.atoms[0u] + 1u) += deddt * dxiayia + d2eddt2 * ddtdxia * ddtdyia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u]     ).at(3 * angle.atoms[0u] + 2u) += deddt * dxiazia + d2eddt2 * ddtdxia * ddtdzia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] +1u ).at(3 * angle.atoms[0u]     ) += deddt * dxiayia + d2eddt2 * ddtdyia * ddtdxia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] +1u ).at(3 * angle.atoms[0u] + 1u) += deddt * dyiayia + d2eddt2 * ddtdyia * ddtdyia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] +1u ).at(3 * angle.atoms[0u] + 2u) += deddt * dyiazia + d2eddt2 * ddtdyia * ddtdzia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 2u).at(3 * angle.atoms[0u] + 0u) += deddt * dxiazia + d2eddt2 * ddtdzia * ddtdxia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 2u).at(3 * angle.atoms[0u] + 1u) += deddt * dyiazia + d2eddt2 * ddtdzia * ddtdyia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 2u).at(3 * angle.atoms[0u] + 2u) += deddt * dziazia + d2eddt2 * ddtdzia * ddtdzia;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 0u).at(3 * angle.atoms[1u] + 0u) += deddt * dxibxia + d2eddt2 * ddtdxia * ddtdxib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 0u).at(3 * angle.atoms[1u] + 1u) += deddt * dyibxia + d2eddt2 * ddtdxia * ddtdyib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 0u).at(3 * angle.atoms[1u] + 2u) += deddt * dzibxia + d2eddt2 * ddtdxia * ddtdzib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 1u).at(3 * angle.atoms[1u] + 0u) += deddt * dxibyia + d2eddt2 * ddtdyia * ddtdxib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 1u).at(3 * angle.atoms[1u] + 1u) += deddt * dyibyia + d2eddt2 * ddtdyia * ddtdyib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 1u).at(3 * angle.atoms[1u] + 2u) += deddt * dzibyia + d2eddt2 * ddtdyia * ddtdzib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 2u).at(3 * angle.atoms[1u] + 0u) += deddt * dxibzia + d2eddt2 * ddtdzia * ddtdxib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 2u).at(3 * angle.atoms[1u] + 1u) += deddt * dyibzia + d2eddt2 * ddtdzia * ddtdyib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 2u).at(3 * angle.atoms[1u] + 2u) += deddt * dzibzia + d2eddt2 * ddtdzia * ddtdzib;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 0u).at(3 * angle.atoms[2u] + 0u) += deddt * dxibxic + d2eddt2 * ddtdxia * ddtdxic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 0u).at(3 * angle.atoms[2u] + 1u) += deddt * dyibxic + d2eddt2 * ddtdxia * ddtdyic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 0u).at(3 * angle.atoms[2u] + 2u) += deddt * dzibxic + d2eddt2 * ddtdxia * ddtdzic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 1u).at(3 * angle.atoms[2u] + 0u) += deddt * dxibyic + d2eddt2 * ddtdyia * ddtdxic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 1u).at(3 * angle.atoms[2u] + 1u) += deddt * dyibyic + d2eddt2 * ddtdyia * ddtdyic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 1u).at(3 * angle.atoms[2u] + 2u) += deddt * dzibyic + d2eddt2 * ddtdyia * ddtdzic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 2u).at(3 * angle.atoms[2u] + 0u) += deddt * dxibzic + d2eddt2 * ddtdzia * ddtdxic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 2u).at(3 * angle.atoms[2u] + 1u) += deddt * dyibzic + d2eddt2 * ddtdzia * ddtdyic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 2u).at(3 * angle.atoms[2u] + 2u) += deddt * dzibzic + d2eddt2 * ddtdzia * ddtdzic;
            //
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 0u).at(3 * angle.atoms[1u] + 0u) += deddt * dxibxib + d2eddt2 * ddtdxib * ddtdxib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 0u).at(3 * angle.atoms[1u] + 1u) += deddt * dxibyib + d2eddt2 * ddtdxib * ddtdyib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 0u).at(3 * angle.atoms[1u] + 2u) += deddt * dxibzib + d2eddt2 * ddtdxib * ddtdzib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 1u).at(3 * angle.atoms[1u] + 0u) += deddt * dxibyib + d2eddt2 * ddtdyib * ddtdxib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 1u).at(3 * angle.atoms[1u] + 1u) += deddt * dyibyib + d2eddt2 * ddtdyib * ddtdyib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 1u).at(3 * angle.atoms[1u] + 2u) += deddt * dyibzib + d2eddt2 * ddtdyib * ddtdzib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 2u).at(3 * angle.atoms[1u] + 0u) += deddt * dxibzib + d2eddt2 * ddtdzib * ddtdxib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 2u).at(3 * angle.atoms[1u] + 1u) += deddt * dyibzib + d2eddt2 * ddtdzib * ddtdyib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 2u).at(3 * angle.atoms[1u] + 2u) += deddt * dzibzib + d2eddt2 * ddtdzib * ddtdzib;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 0u).at(3 * angle.atoms[0u] + 0u) += deddt * dxibxia + d2eddt2 * ddtdxib * ddtdxia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 0u).at(3 * angle.atoms[0u] + 1u) += deddt * dxibyia + d2eddt2 * ddtdxib * ddtdyia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 0u).at(3 * angle.atoms[0u] + 2u) += deddt * dxibzia + d2eddt2 * ddtdxib * ddtdzia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 1u).at(3 * angle.atoms[0u] + 0u) += deddt * dyibxia + d2eddt2 * ddtdyib * ddtdxia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 1u).at(3 * angle.atoms[0u] + 1u) += deddt * dyibyia + d2eddt2 * ddtdyib * ddtdyia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 1u).at(3 * angle.atoms[0u] + 2u) += deddt * dyibzia + d2eddt2 * ddtdyib * ddtdzia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 2u).at(3 * angle.atoms[0u] + 0u) += deddt * dzibxia + d2eddt2 * ddtdzib * ddtdxia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 2u).at(3 * angle.atoms[0u] + 1u) += deddt * dzibyia + d2eddt2 * ddtdzib * ddtdyia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 2u).at(3 * angle.atoms[0u] + 2u) += deddt * dzibzia + d2eddt2 * ddtdzib * ddtdzia;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 0u).at(3 * angle.atoms[2u] + 0u) += deddt * dxibxic + d2eddt2 * ddtdxib * ddtdxic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 0u).at(3 * angle.atoms[2u] + 1u) += deddt * dxibyic + d2eddt2 * ddtdxib * ddtdyic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 0u).at(3 * angle.atoms[2u] + 2u) += deddt * dxibzic + d2eddt2 * ddtdxib * ddtdzic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 1u).at(3 * angle.atoms[2u] + 0u) += deddt * dyibxic + d2eddt2 * ddtdyib * ddtdxic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 1u).at(3 * angle.atoms[2u] + 1u) += deddt * dyibyic + d2eddt2 * ddtdyib * ddtdyic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 1u).at(3 * angle.atoms[2u] + 2u) += deddt * dyibzic + d2eddt2 * ddtdyib * ddtdzic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 2u).at(3 * angle.atoms[2u] + 0u) += deddt * dzibxic + d2eddt2 * ddtdzib * ddtdxic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 2u).at(3 * angle.atoms[2u] + 1u) += deddt * dzibyic + d2eddt2 * ddtdzib * ddtdyic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 2u).at(3 * angle.atoms[2u] + 2u) += deddt * dzibzic + d2eddt2 * ddtdzib * ddtdzic;
            //
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 0u).at(3 * angle.atoms[2u] + 0u) += deddt * dxicxic + d2eddt2 * ddtdxic * ddtdxic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 0u).at(3 * angle.atoms[2u] + 1u) += deddt * dxicyic + d2eddt2 * ddtdxic * ddtdyic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 0u).at(3 * angle.atoms[2u] + 2u) += deddt * dxiczic + d2eddt2 * ddtdxic * ddtdzic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 1u).at(3 * angle.atoms[2u] + 0u) += deddt * dxicyic + d2eddt2 * ddtdyic * ddtdxic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 1u).at(3 * angle.atoms[2u] + 1u) += deddt * dyicyic + d2eddt2 * ddtdyic * ddtdyic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 1u).at(3 * angle.atoms[2u] + 2u) += deddt * dyiczic + d2eddt2 * ddtdyic * ddtdzic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 2u).at(3 * angle.atoms[2u] + 0u) += deddt * dxiczic + d2eddt2 * ddtdzic * ddtdxic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 2u).at(3 * angle.atoms[2u] + 1u) += deddt * dyiczic + d2eddt2 * ddtdzic * ddtdyic;
            part_hessian[ANGLE].at(3 * angle.atoms[2u] + 2u).at(3 * angle.atoms[2u] + 2u) += deddt * dziczic + d2eddt2 * ddtdzic * ddtdzic;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 0u).at(3 * angle.atoms[1u] + 0u) += deddt * dxibxic + d2eddt2 * ddtdxic * ddtdxib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 0u).at(3 * angle.atoms[1u] + 1u) += deddt * dxibyic + d2eddt2 * ddtdxic * ddtdyib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 0u).at(3 * angle.atoms[1u] + 2u) += deddt * dxibzic + d2eddt2 * ddtdxic * ddtdzib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 1u).at(3 * angle.atoms[1u] + 0u) += deddt * dyibxic + d2eddt2 * ddtdyic * ddtdxib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 1u).at(3 * angle.atoms[1u] + 1u) += deddt * dyibyic + d2eddt2 * ddtdyic * ddtdyib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 1u).at(3 * angle.atoms[1u] + 2u) += deddt * dyibzic + d2eddt2 * ddtdyic * ddtdzib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 2u).at(3 * angle.atoms[1u] + 0u) += deddt * dzibxic + d2eddt2 * ddtdzic * ddtdxib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 2u).at(3 * angle.atoms[1u] + 1u) += deddt * dzibyic + d2eddt2 * ddtdzic * ddtdyib;
            part_hessian[ANGLE].at(3 * angle.atoms[1u] + 2u).at(3 * angle.atoms[1u] + 2u) += deddt * dzibzic + d2eddt2 * ddtdzic * ddtdzib;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 0u).at(3 * angle.atoms[0u] + 0u) += deddt * dxiaxic + d2eddt2 * ddtdxic * ddtdxia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 0u).at(3 * angle.atoms[0u] + 1u) += deddt * dxiayic + d2eddt2 * ddtdxic * ddtdyia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 0u).at(3 * angle.atoms[0u] + 2u) += deddt * dxiazic + d2eddt2 * ddtdxic * ddtdzia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 1u).at(3 * angle.atoms[0u] + 0u) += deddt * dyiaxic + d2eddt2 * ddtdyic * ddtdxia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 1u).at(3 * angle.atoms[0u] + 1u) += deddt * dyiayic + d2eddt2 * ddtdyic * ddtdyia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 1u).at(3 * angle.atoms[0u] + 2u) += deddt * dyiazic + d2eddt2 * ddtdyic * ddtdzia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 2u).at(3 * angle.atoms[0u] + 0u) += deddt * dziaxic + d2eddt2 * ddtdzic * ddtdxia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 2u).at(3 * angle.atoms[0u] + 1u) += deddt * dziayic + d2eddt2 * ddtdzic * ddtdyia;
            part_hessian[ANGLE].at(3 * angle.atoms[0u] + 2u).at(3 * angle.atoms[0u] + 2u) += deddt * dziazic + d2eddt2 * ddtdzic * ddtdzia;
            // end subrotuine eangle2a
            // Only AMOEBA FF ha in-plane angle parameters, we therefore can skip the rest of the tinker file and are done.
          }
          else
          {
            if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because of angle " << angle << "\n";
            integrity = false;
          }
        }
        return E;
      }


      /*!*******************************
      *                                *
      *                                *
      * CHARMM Urey-Bradley Potential  *
      * Energy/Gradients/Hessians      *
      *                                *
      *                                *
      *********************************/

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_13_u<0>(void)
      {
        using scon::len;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const& urey : refined.ureys())
        {
          coords::Cartesian_Point const bv =
            coords->xyz(urey.atoms[0]) - coords->xyz(urey.atoms[1]);
          coords::float_type const d = len(bv);
          coords::float_type const r = d - urey.ideal;
          E += urey.force * r * r;
          if (abs(d) < 1.0e-8)
          {
            if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because of urey " << urey << "\n";
            integrity = false;
          }
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_13_u<1>(void)
      {
        using scon::len;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const& urey : refined.ureys())
        {
          coords::Cartesian_Point const bv(coords->xyz(urey.atoms[0]) - coords->xyz(urey.atoms[1]));
          coords::float_type const d = len(bv);
          coords::float_type const r = d - urey.ideal;
          coords::float_type dE = urey.force * r;
          E += dE * r;
          dE *= 2;
          if (abs(d) > 0.0)
          {
            dE /= d;
            coords::Cartesian_Point const gv = bv * dE;
            part_grad[UREY][urey.atoms[0]] += gv;
            part_grad[UREY][urey.atoms[1]] -= gv;
            //increment internal virial tensor
            coords::float_type const vxx = bv.x() * gv.x();
            coords::float_type const vyx = bv.y() * gv.x();
            coords::float_type const vzx = bv.z() * gv.x();
            coords::float_type const vyy = bv.y() * gv.y();
            coords::float_type const vzy = bv.z() * gv.y();
            coords::float_type const vzz = bv.z() * gv.z();
            part_virial[UREY][0][0] += vxx;
            part_virial[UREY][1][0] += vyx;
            part_virial[UREY][2][0] += vzx;
            part_virial[UREY][0][1] += vyx;
            part_virial[UREY][1][1] += vyy;
            part_virial[UREY][2][1] += vzy;
            part_virial[UREY][0][2] += vzx;
            part_virial[UREY][1][2] += vzy;
            part_virial[UREY][2][2] += vzz;
          }
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_13_u<2>(void)
      {
        using scon::len;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const& urey : refined.ureys())
        {
          coords::Cartesian_Point const bv(coords->xyz(urey.atoms[0]) - coords->xyz(urey.atoms[1]));
          coords::float_type const d = len(bv);
          coords::float_type const r = d - urey.ideal;
          coords::float_type dE = urey.force * r;
          E += dE * r;
          dE *= 2;
          if (abs(d) > 0.0)
          {
            dE /= d;
            coords::Cartesian_Point const gv = bv * dE;
            part_grad[UREY][urey.atoms[0]] += gv;
            part_grad[UREY][urey.atoms[1]] -= gv;
            //increment internal virial tensor
            coords::float_type const vxx = bv.x() * gv.x();
            coords::float_type const vyx = bv.y() * gv.x();
            coords::float_type const vzx = bv.z() * gv.x();
            coords::float_type const vyy = bv.y() * gv.y();
            coords::float_type const vzy = bv.z() * gv.y();
            coords::float_type const vzz = bv.z() * gv.z();
            part_virial[UREY][0][0] += vxx;
            part_virial[UREY][1][0] += vyx;
            part_virial[UREY][2][0] += vzx;
            part_virial[UREY][0][1] += vyx;
            part_virial[UREY][1][1] += vyy;
            part_virial[UREY][2][1] += vzy;
            part_virial[UREY][0][2] += vzx;
            part_virial[UREY][1][2] += vzy;
            part_virial[UREY][2][2] += vzz;
            // Hessian
            const double d2eddt2 = 2.0 * urey.force;
            const double& de = dE;
            const double term = (d2eddt2 - de) / (d*d);
            const double xac = bv.x();
            const double yac = bv.y();
            const double zac = bv.z();
            const double termx = term * xac          ;
            const double termy = term * yac          ;
            const double termz = term * zac          ;
            std::array<std::array<coords::float_type, 3u>, 3u> d2e = std::array<std::array<coords::float_type, 3u>, 3u>();
            d2e.at(0).at(0) = termx * xac + dE;
            d2e.at(0).at(1) = termx * yac;
            d2e.at(0).at(2) = termx * zac;
            d2e.at(1).at(0) = d2e.at(0).at(0);
            d2e.at(1).at(1) = termy * yac + dE;
            d2e.at(1).at(2) = termy * zac;
            d2e.at(2).at(0) = d2e.at(0).at(2);
            d2e.at(2).at(1) = d2e.at(1).at(2);
            d2e.at(2).at(2) = termz * zac + dE;
            part_hessian[BOND].at(3 * urey.atoms[0u]).at(3 *      urey.atoms[0u]) += d2e.at(0).at(0); // dx dx same atom
            part_hessian[BOND].at(3 * urey.atoms[0u]).at(3 *      urey.atoms[0u] + 1u) += d2e.at(0).at(1); // dx dy same atom
            part_hessian[BOND].at(3 * urey.atoms[0u]).at(3 *      urey.atoms[0u] + 2u) += d2e.at(0).at(2); // dx dz same atom
            part_hessian[BOND].at(3 * urey.atoms[0u] + 1u).at(3 * urey.atoms[0u]) += d2e.at(1).at(0); // dy dx same atom
            part_hessian[BOND].at(3 * urey.atoms[0u] + 1u).at(3 * urey.atoms[0u] + 1u) += d2e.at(1).at(1); // dy dy same atom
            part_hessian[BOND].at(3 * urey.atoms[0u] + 1u).at(3 * urey.atoms[0u] + 2u) += d2e.at(1).at(2); // dy dz same atom
            part_hessian[BOND].at(3 * urey.atoms[0u] + 2u).at(3 * urey.atoms[0u]) += d2e.at(2).at(0); // dy dx same atom
            part_hessian[BOND].at(3 * urey.atoms[0u] + 2u).at(3 * urey.atoms[0u] + 1u) += d2e.at(2).at(1); // dy dy same atom
            part_hessian[BOND].at(3 * urey.atoms[0u] + 2u).at(3 * urey.atoms[0u] + 2u) += d2e.at(2).at(2); // dy dz same atom
            part_hessian[BOND].at(3 * urey.atoms[1u]).at(3 *      urey.atoms[1u]) -= d2e.at(0).at(0); // dx dx same atom
            part_hessian[BOND].at(3 * urey.atoms[1u]).at(3 *      urey.atoms[1u] + 1u) -= d2e.at(0).at(1); // dx dy same atom
            part_hessian[BOND].at(3 * urey.atoms[1u]).at(3 *      urey.atoms[1u] + 2u) -= d2e.at(0).at(2); // dx dz same atom
            part_hessian[BOND].at(3 * urey.atoms[1u] + 1u).at(3 * urey.atoms[1u]) -= d2e.at(1).at(0); // dy dx same atom
            part_hessian[BOND].at(3 * urey.atoms[1u] + 1u).at(3 * urey.atoms[1u] + 1u) -= d2e.at(1).at(1); // dy dy same atom
            part_hessian[BOND].at(3 * urey.atoms[1u] + 1u).at(3 * urey.atoms[1u] + 2u) -= d2e.at(1).at(2); // dy dz same atom
            part_hessian[BOND].at(3 * urey.atoms[1u] + 2u).at(3 * urey.atoms[1u]) -= d2e.at(2).at(0); // dy dx same atom
            part_hessian[BOND].at(3 * urey.atoms[1u] + 2u).at(3 * urey.atoms[1u] + 1u) -= d2e.at(2).at(1); // dy dy same atom
            part_hessian[BOND].at(3 * urey.atoms[1u] + 2u).at(3 * urey.atoms[1u] + 2u) -= d2e.at(2).at(2); // dy dz same atom
          }
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_14<0>(void)
      {
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const& torsion : refined.torsions())
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
          coords::float_type const tlul = sqrt(tl2 * ul2);
          coords::float_type const r12 = len(b12);
          // cross of cross
          coords::Cartesian_Point const tu = cross(t, u);
          // scalar and length variations
          coords::float_type const cos_scalar0 = dot(t, u);
          coords::float_type const cos_scalar1 = tlul;
          coords::float_type const sin_scalar0 = dot(b12, tu);
          coords::float_type const sin_scalar1 = r12 * tlul;
          // check whether 
          if (abs(cos_scalar1) < 1.0e-8 || abs(sin_scalar1) < 1.0e-8)
          {
            if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because of torsion " << torsion << "\n";
            integrity = false;
          }
          // Get multiple sine and cosine values
          coords::float_type cos[7], sin[7];
          cos[1] = cos_scalar0 / cos_scalar1;
          sin[1] = sin_scalar0 / sin_scalar1;
          for (std::size_t j(2U); j <= torsion.p.max_order; ++j)
          {
            std::size_t const k = j - 1;
            sin[j] = sin[k] * cos[1] + cos[k] * sin[1];
            cos[j] = cos[k] * cos[1] - sin[k] * sin[1];
          }

          // Coding note:
          // In Tinker-Code, v is the torsional force,
          // and c is the cosine-parameter.
          // Original tinker also supports sine-parameters, but 
          // They are always zero for the forcefields we use in CAST

          coords::float_type tE(0.0);
          //cout << "Number?:  " << torsions[i].paramPtr->n << std::endl;
          for (std::size_t j(0U); j < torsion.p.number; ++j)
          {
            coords::float_type const F = torsion.p.force[j] * cparams.torsionunit();
            std::size_t const k = torsion.p.order[j];
            coords::float_type const l = std::abs(torsion.p.ideal[j]) > 0.0 ? -1.0 : 1.0;
            tE += F * (1.0 + cos[k] * l);
          }
          E += tE;
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_14<1>(void)
      {
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const& torsion : refined.torsions())
        {
          E += f_14_1_calc_one_torsion(torsion);
        }
        return E;
      }

      coords::float_type energy::interfaces::aco::aco_ff::f_14_1_calc_one_torsion(::tinker::refine::types::torsion const& torsion)
      {
        coords::float_type E = 0.;
        coords::Cartesian_Point const b01 =
          coords->xyz(torsion.atoms[1]) - coords->xyz(torsion.atoms[0]);
        coords::Cartesian_Point const b12 =
          coords->xyz(torsion.atoms[2]) - coords->xyz(torsion.atoms[1]);
        coords::Cartesian_Point const b23 =
          coords->xyz(torsion.atoms[3]) - coords->xyz(torsion.atoms[2]);
        coords::Cartesian_Point const b02 =
          coords->xyz(torsion.atoms[2]) - coords->xyz(torsion.atoms[0]);
        coords::Cartesian_Point const b13 =
          coords->xyz(torsion.atoms[3]) - coords->xyz(torsion.atoms[1]);

        coords::Cartesian_Point const t = cross(b01, b12);
        coords::Cartesian_Point const u = cross(b12, b23);

        coords::float_type const tl2 = dot(t, t);
        coords::float_type const ul2 = dot(u, u);
        coords::float_type const tlul = sqrt(tl2 * ul2);
        coords::float_type const r12 = len(b12);

        coords::Cartesian_Point const tu = cross(t, u);

        coords::float_type const cos_scalar0 = dot(t, u);
        coords::float_type const cos_scalar1 = tlul;

        coords::float_type const sin_scalar0 = dot(b12, tu);
        coords::float_type const sin_scalar1 = r12 * tlul;

        if (abs(cos_scalar1) < 1.0e-8 || abs(sin_scalar1) < 1.0e-8)
        {
          if (Config::get().general.verbosity > 3) std::cout << "WARNING! Integrity broke because of torsion " << torsion << "\n";
          integrity = false;
        }

        coords::float_type cos[7], sin[7];
        cos[1] = cos_scalar0 / cos_scalar1;
        sin[1] = sin_scalar0 / sin_scalar1;

        for (std::size_t j(2U); j <= torsion.p.max_order; ++j)
        {
          std::size_t const k = j - 1;
          sin[j] = sin[k] * cos[1] + cos[k] * sin[1];
          cos[j] = cos[k] * cos[1] - sin[k] * sin[1];
        }

        coords::float_type tE(0.0), dE(0.0);
        //cout << "Number?:  " << torsions[i].paramPtr->n << std::endl;
        for (std::size_t j(0U); j < torsion.p.number; ++j)
        {
          coords::float_type const F = torsion.p.force[j] * cparams.torsionunit();
          std::size_t const k = torsion.p.order[j];
          coords::float_type const l = std::abs(torsion.p.ideal[j]) > 0.0 ? -1.0 : 1.0;
          tE += F * (1.0 + cos[k] * l);
          dE += -static_cast<coords::float_type>(k)* F* sin[k] * l;
        }
        E += tE;

        coords::Cartesian_Point const dt(cross(t, b12) * (dE / (tl2 * r12)));
        coords::Cartesian_Point const du(cross(u, b12) * (-dE / (ul2 * r12)));

        coords::Cartesian_Point const vir1 = cross(dt, b12);
        coords::Cartesian_Point const vir2 = cross(b02, dt) + cross(du, b23);
        coords::Cartesian_Point const vir3 = cross(dt, b01) + cross(b13, du);
        coords::Cartesian_Point const vir4 = cross(du, b12);

        part_grad[TORSION][torsion.atoms[0]] += vir1;
        part_grad[TORSION][torsion.atoms[1]] += vir2;
        part_grad[TORSION][torsion.atoms[2]] += vir3;
        part_grad[TORSION][torsion.atoms[3]] += vir4;

        //increment internal virial tensor
        coords::float_type const vxx = b12.x() * (vir3.x() + vir4.x()) -
          b01.x() * vir1.x() + b23.x() * vir4.x();
        coords::float_type const vyx = b12.y() * (vir3.x() + vir4.x()) -
          b01.y() * vir1.x() + b23.y() * vir4.x();
        coords::float_type const vzx = b12.z() * (vir3.x() + vir4.x()) -
          b01.z() * vir1.x() + b23.z() * vir4.x();
        coords::float_type const vyy = b12.y() * (vir3.y() + vir4.y()) -
          b01.y() * vir1.y() + b23.y() * vir4.y();
        coords::float_type const vzy = b12.z() * (vir3.y() + vir4.y()) -
          b01.z() * vir1.y() + b23.z() * vir4.y();
        coords::float_type const vzz = b12.z() * (vir3.z() + vir4.z()) -
          b01.z() * vir1.z() + b23.z() * vir4.z();
        part_virial[TORSION][0][0] += vxx;
        part_virial[TORSION][1][0] += vyx;
        part_virial[TORSION][2][0] += vzx;
        part_virial[TORSION][0][1] += vyx;
        part_virial[TORSION][1][1] += vyy;
        part_virial[TORSION][2][1] += vzy;
        part_virial[TORSION][0][2] += vzx;
        part_virial[TORSION][1][2] += vzy;
        part_virial[TORSION][2][2] += vzz;
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_14<2>(void)
      {
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt;
        using std::abs;
        coords::float_type E(0.0);
        std::vector < std::vector<double*>> hessx(4u,std::vector<double*>(unsigned int(coords->atoms().size()), nullptr));
        std::vector < std::vector<double*>> hessy(hessx), hessz(hessx);
        for (std::size_t i = 1u; i <= 3u; ++i)
        {
          for (std::size_t atom = 0u; atom < coords->atoms().size(); ++atom)
          {
            std::size_t param = atom * 3 + (i - 1u);
            hessx[i][atom] = &part_hessian[TORSION][atom * 3][atom * 3];
            hessy[i][atom] = &part_hessian[TORSION][atom * 3][atom * 3 + 1u];
            hessz[i][atom] = &part_hessian[TORSION][atom * 3][atom * 3 + 2u];
          }
        }

        for (auto const& torsion : refined.torsions())
        {
          E += f_14_1_calc_one_torsion(torsion);

          const double xia = coords->xyz(torsion.atoms[0]).x();
          const double yia = coords->xyz(torsion.atoms[0]).y();
          const double zia = coords->xyz(torsion.atoms[0]).z();
          const double xib = coords->xyz(torsion.atoms[1]).x();
          const double yib = coords->xyz(torsion.atoms[1]).y();
          const double zib = coords->xyz(torsion.atoms[1]).z();
          const double xic = coords->xyz(torsion.atoms[2]).x();
          const double yic = coords->xyz(torsion.atoms[2]).y();
          const double zic = coords->xyz(torsion.atoms[2]).z();
          const double xid = coords->xyz(torsion.atoms[3]).x();
          const double yid = coords->xyz(torsion.atoms[3]).y();
          const double zid = coords->xyz(torsion.atoms[3]).z();
          const double xba = xib - xia  ;
          const double   yba = yib - yia;
          const double   zba = zib - zia;
          const double   xcb = xic - xib;
          const double   ycb = yic - yib;
          const double   zcb = zic - zib;
          const double   xdc = xid - xic;
          const double   ydc = yid - yic;
          const double   zdc = zid - zic;
          const double xt = yba * zcb - ycb * zba            ;
          const double   yt = zba * xcb - zcb * xba          ;
          const double   zt = xba * ycb - xcb * yba          ;
          const double   xu = ycb * zdc - ydc * zcb          ;
          const double   yu = zcb * xdc - zdc * xcb          ;
          const double   zu = xcb * ydc - xdc * ycb          ;
          const double   xtu = yt * zu - yu * zt             ;
          const double   ytu = zt * xu - zu * xt             ;
          const double   ztu = xt * yu - xu * yt             ;
          const double   rt2 = xt * xt + yt * yt + zt * zt   ;
          const double   ru2 = xu * xu + yu * yu + zu * zu   ;
          const double   rtru = sqrt(rt2 * ru2)              ;
          if (rtru != 0.)
          {
            coords::float_type cos[7], sin[7];
            const double rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb)             ;
            cos[0] = (xt * xu + yt * yu + zt * zu) / rtru             ;
            const double& cosine = cos[0];
            sin[0] = (xcb * xtu + ycb * ytu + zcb * ztu) / (rcb * rtru) ;
            const double& sine = sin[0];
            cos[1] = cosine * cosine - sine * sine;
            const double& cosine2 = cos[1];
            sin[1] = 2. * cosine * sine;
            const double& sine2 = sin[1];
            cos[2] = cosine * cosine2 - sine * sine2;
            const double& cosine3 = cos[2];
            sin[2] = cosine * sine2 + sine * cosine2;
            const double sine3 = sin[2];
            cos[3] = cosine * cosine3 - sine * sine3;
            const double& cosine4 = cos[3];
            sin[3] = cosine * sine3 + sine * cosine3;
            const double& sine4 = sin[3];
            cos[4] = cosine * cosine4 - sine * sine4;
            const double cosine5 = cos[4];
            sin[4] = cosine * sine4 + sine * cosine4;
            const double& sine5 = sin[4];
            cos[5] = cosine * cosine5 - sine * sine5;
            const double& cosine6 = cos[5];
            sin[5] = cosine * sine5 + sine * cosine5;
            const double& sine6 = sin[5];

            for (std::size_t j(0U); j < torsion.p.number; ++j)
              {

                coords::float_type const F = torsion.p.force[j] * cparams.torsionunit();
                std::size_t const k = torsion.p.order[j];
                coords::float_type const l = std::abs(torsion.p.ideal[j]) > 0.0 ? -1.0 : 1.0;

                const double d2phi1 = -cos[k] * l;
                const double dphi1 = (sine*l);
                const double dedphi = cparams.torsionunit() * torsion.p.force[j] * dphi1;
                const double d2edphi2 = cparams.torsionunit() * torsion.p.force[j] * d2phi1;

                //abbreviations for first derivative chain rule terms
                const double xca = xic - xia;
                const double yca = yic - yia;
                const double zca = zic - zia;
                const double xdb = xid - xib;
                const double ydb = yid - yib;
                const double zdb = zid - zib;
                const double dphidxt = (yt * zcb - ycb * zt) / (rt2 * rcb) ;
                const double dphidyt = (zt * xcb - zcb * xt) / (rt2 * rcb) ;
                const double dphidzt = (xt * ycb - xcb * yt) / (rt2 * rcb) ;
                const double dphidxu = -(yu * zcb - ycb * zu) / (ru2 * rcb);
                const double dphidyu = -(zu * xcb - zcb * xu) / (ru2 * rcb);
                const double dphidzu = -(xu * ycb - xcb * yu) / (ru2 * rcb);

                //abbreviations for second derivative chain rule terms
                const double xycb2 = xcb * xcb + ycb * ycb              ;
                const double   xzcb2 = xcb * xcb + zcb * zcb            ;
                const double   yzcb2 = ycb * ycb + zcb * zcb            ;
                const double   rcbxt = -2.0 * rcb * dphidxt           ;
                const double   rcbyt = -2.0 * rcb * dphidyt           ;
                const double   rcbzt = -2.0 * rcb * dphidzt           ;
                const double   rcbt2 = rcb * rt2                        ;
                const double   rcbxu = 2.0 * rcb * dphidxu            ;
                const double   rcbyu = 2.0 * rcb * dphidyu            ;
                const double   rcbzu = 2.0 * rcb * dphidzu            ;
                const double   rcbu2 = rcb * ru2                        ;
                const double   dphidxibt = yca * dphidzt - zca * dphidyt;
                const double   dphidxibu = zdc * dphidyu - ydc * dphidzu;
                const double   dphidyibt = zca * dphidxt - xca * dphidzt;
                const double   dphidyibu = xdc * dphidzu - zdc * dphidxu;
                const double   dphidzibt = xca * dphidyt - yca * dphidxt;
                const double   dphidzibu = ydc * dphidxu - xdc * dphidyu;
                const double   dphidxict = zba * dphidyt - yba * dphidzt;
                const double   dphidxicu = ydb * dphidzu - zdb * dphidyu;
                const double   dphidyict = xba * dphidzt - zba * dphidxt;
                const double   dphidyicu = zdb * dphidxu - xdb * dphidzu;
                const double   dphidzict = yba * dphidxt - xba * dphidyt;
                const double   dphidzicu = xdb * dphidyu - ydb * dphidxu;

                //chain rule terms for first derivative components
                const double dphidxia = zcb * dphidyt - ycb * dphidzt;
                const double dphidyia = xcb * dphidzt - zcb * dphidxt;
                const double dphidzia = ycb * dphidxt - xcb * dphidyt;
                const double dphidxib = dphidxibt + dphidxibu        ;
                const double dphidyib = dphidyibt + dphidyibu        ;
                const double dphidzib = dphidzibt + dphidzibu        ;
                const double dphidxic = dphidxict + dphidxicu        ;
                const double dphidyic = dphidyict + dphidyicu        ;
                const double dphidzic = dphidzict + dphidzicu        ;
                const double dphidxid = zcb * dphidyu - ycb * dphidzu;
                const double dphidyid = xcb * dphidzu - zcb * dphidxu;
                const double dphidzid = ycb * dphidxu - xcb * dphidyu;

                //chain rule terms for second derivative components
                const double dxiaxia = rcbxt * dphidxia                                                      ;
                const double dxiayia = rcbxt * dphidyia - zcb * rcb / rt2                                    ;
                const double dxiazia = rcbxt * dphidzia + ycb * rcb / rt2                                    ;
                const double dxiaxic = rcbxt * dphidxict + xcb * xt / rcbt2                                  ;
                const double dxiayic = rcbxt * dphidyict - dphidzt - (xba * zcb * xcb + zba * yzcb2) / rcbt2 ;
                const double dxiazic = rcbxt * dphidzict + dphidyt + (xba * ycb * xcb + yba * yzcb2) / rcbt2 ;
                const double dxiaxid = 0.0                                                                 ;
                const double dxiayid = 0.0                                                                 ;
                const double dxiazid = 0.0                                                                 ;
                const double dyiayia = rcbyt * dphidyia                                                      ;
                const double dyiazia = rcbyt * dphidzia - xcb * rcb / rt2                                    ;
                const double dyiaxib = rcbyt * dphidxibt - dphidzt - (yca * zcb * ycb + zca * xzcb2) / rcbt2 ;
                const double dyiaxic = rcbyt * dphidxict + dphidzt + (yba * zcb * ycb + zba * xzcb2) / rcbt2 ;
                const double dyiayic = rcbyt * dphidyict + ycb * yt / rcbt2                                  ;
                const double dyiazic = rcbyt * dphidzict - dphidxt - (yba * xcb * ycb + xba * xzcb2) / rcbt2 ;
                const double dyiaxid = 0.0                                                                 ;
                const double dyiayid = 0.0                                                                 ;
                const double dyiazid = 0.0                                                                 ;
                const double dziazia = rcbzt * dphidzia                                                      ;
                const double dziaxib = rcbzt * dphidxibt + dphidyt + (zca * ycb * zcb + yca * xycb2) / rcbt2 ;
                const double dziayib = rcbzt * dphidyibt - dphidxt - (zca * xcb * zcb + xca * xycb2) / rcbt2 ;
                const double dziaxic = rcbzt * dphidxict - dphidyt - (zba * ycb * zcb + yba * xycb2) / rcbt2 ;
                const double dziayic = rcbzt * dphidyict + dphidxt + (zba * xcb * zcb + xba * xycb2) / rcbt2 ;
                const double dziazic = rcbzt * dphidzict + zcb * zt / rcbt2                                  ;
                const double dziaxid = 0.0                                                                 ;
                const double dziayid = 0.0                                                                 ;
                const double dziazid = 0.0                                                                 ;
                const double dxibxic = -xcb * dphidxib / (rcb * rcb) - (yca * (zba * xcb + yt) - zca * (yba * xcb - zt)) / rcbt2 - 2.0 * (yt * zba - yba * zt) * dphidxibt / rt2 - (zdc * (ydb * xcb + zu) - ydc * (zdb * xcb - yu)) / rcbu2 + 2.0 * (yu * zdb - ydb * zu) * dphidxibu / ru2;
                const double dxibyic = -ycb * dphidxib / (rcb * rcb) + dphidzt + dphidzu - (yca * (zba * ycb - xt) + zca * (xba * xcb + zcb * zba)) / rcbt2 - 2.0 * (zt * xba - zba * xt) * dphidxibt / rt2 + (zdc * (xdb * xcb + zcb * zdb) + ydc * (zdb * ycb + xu)) / rcbu2 + 2.0 * (zu * xdb - zdb * xu) * dphidxibu / ru2;
                const double dxibxid = rcbxu * dphidxibu + xcb * xu / rcbu2;
                const double dxibyid = rcbyu * dphidxibu - dphidzu - (ydc * zcb * ycb + zdc * xzcb2) / rcbu2;
                const double dxibzid = rcbzu * dphidxibu + dphidyu + (zdc * ycb * zcb + ydc * xycb2) / rcbu2;
                const double dyibzib = ycb * dphidzib / (rcb * rcb) - (xca * (xca * xcb + zcb * zca) + yca * (ycb * xca + zt)) / rcbt2 - 2.0 * (xt * zca - xca * zt) * dphidzibt / rt2 + (ydc * (xdc * ycb - zu) + xdc * (xdc * xcb + zcb * zdc)) / rcbu2 + 2.0 * (xu * zdc - xdc * zu) * dphidzibu / ru2;
                const double dyibxic = -xcb * dphidyib / (rcb * rcb) - dphidzt - dphidzu + (xca * (zba * xcb + yt) + zca * (zba * zcb + ycb * yba)) / rcbt2 - 2.0 * (yt * zba - yba * zt) * dphidyibt / rt2 - (zdc * (zdb * zcb + ycb * ydb) + xdc * (zdb * xcb - yu)) / rcbu2 + 2.0 * (yu * zdb - ydb * zu) * dphidyibu / ru2;
                const double dyibyic = -ycb * dphidyib / (rcb * rcb) - (zca * (xba * ycb + zt) - xca * (zba * ycb - xt)) / rcbt2 - 2.0 * (zt * xba - zba * xt) * dphidyibt / rt2 - (xdc * (zdb * ycb + xu) - zdc * (xdb * ycb - zu)) / rcbu2 + 2.0 * (zu * xdb - zdb * xu) * dphidyibu / ru2;
                const double dyibxid = rcbxu * dphidyibu + dphidzu + (xdc * zcb * xcb + zdc * yzcb2) / rcbu2;
                const double dyibyid = rcbyu * dphidyibu + ycb * yu / rcbu2;
                const double dyibzid = rcbzu * dphidyibu - dphidxu - (zdc * xcb * zcb + xdc * xycb2) / rcbu2;
                const double dzibxic = -xcb * dphidzib / (rcb * rcb) + dphidyt + dphidyu - (xca * (yba * xcb - zt) + yca * (zba * zcb + ycb * yba)) / rcbt2 - 2.0 * (yt * zba - yba * zt) * dphidzibt / rt2 + (ydc * (zdb * zcb + ycb * ydb) + xdc * (ydb * xcb + zu)) / rcbu2 + 2.0 * (yu * zdb - ydb * zu) * dphidzibu / ru2;
                const double dzibzic = -zcb * dphidzib / (rcb * rcb) - (xca * (yba * zcb + xt) - yca * (xba * zcb - yt)) / rcbt2 - 2.0 * (xt * yba - xba * yt) * dphidzibt / rt2 - (ydc * (xdb * zcb + yu) - xdc * (ydb * zcb - xu)) / rcbu2 + 2.0 * (xu * ydb - xdb * yu) * dphidzibu / ru2;
                const double dzibxid = rcbxu * dphidzibu - dphidyu - (xdc * ycb * xcb + ydc * yzcb2) / rcbu2;
                const double dzibyid = rcbyu * dphidzibu + dphidxu + (ydc * xcb * ycb + xdc * xzcb2) / rcbu2;
                const double dzibzid = rcbzu * dphidzibu + zcb * zu / rcbu2;
                const double dxicxid = rcbxu * dphidxicu - xcb * (zdb * ycb - ydb * zcb) / rcbu2;
                const double dxicyid = rcbyu * dphidxicu + dphidzu + (ydb * zcb * ycb + zdb * xzcb2) / rcbu2;
                const double dxiczid = rcbzu * dphidxicu - dphidyu - (zdb * ycb * zcb + ydb * xycb2) / rcbu2;
                const double dyicxid = rcbxu * dphidyicu - dphidzu - (xdb * zcb * xcb + zdb * yzcb2) / rcbu2;
                const double dyicyid = rcbyu * dphidyicu - ycb * (xdb * zcb - zdb * xcb) / rcbu2            ;
                const double dyiczid = rcbzu * dphidyicu + dphidxu + (zdb * xcb * zcb + xdb * xycb2) / rcbu2;
                const double dzicxid = rcbxu * dphidzicu + dphidyu + (xdb * ycb * xcb + ydb * yzcb2) / rcbu2;
                const double dzicyid = rcbyu * dphidzicu - dphidxu - (ydb * xcb * ycb + xdb * xzcb2) / rcbu2;
                const double dziczid = rcbzu * dphidzicu - zcb * (ydb * xcb - xdb * ycb) / rcbu2            ;
                const double dxidxid = rcbxu * dphidxid;
                const double dxidyid = rcbxu * dphidyid + zcb * rcb / ru2;
                const double dxidzid = rcbxu * dphidzid - ycb * rcb / ru2;
                const double dyidyid = rcbyu * dphidyid;
                const double dyidzid = rcbyu * dphidzid + xcb * rcb / ru2;
                const double dzidzid = rcbzu * dphidzid;

                // get some second derivative chain rule terms by difference
                const double dxiaxib = -dxiaxia - dxiaxic - dxiaxid;
                const double dxiayib = -dxiayia - dxiayic - dxiayid;
                const double dxiazib = -dxiazia - dxiazic - dxiazid;
                const double dyiayib = -dyiayia - dyiayic - dyiayid;
                const double dyiazib = -dyiazia - dyiazic - dyiazid;
                const double dziazib = -dziazia - dziazic - dziazid;
                const double dxibxib = -dxiaxib - dxibxic - dxibxid;
                const double dxibyib = -dyiaxib - dxibyic - dxibyid;
                const double dxibzib = -dxiazib - dzibxic - dzibxid;
                const double dxibzic = -dziaxib - dxibzib - dxibzid;
                const double dyibyib = -dyiayib - dyibyic - dyibyid;
                const double dyibzic = -dziayib - dyibzib - dyibzid;
                const double dzibzib = -dziazib - dzibzic - dzibzid;
                const double dzibyic = -dyiazib - dyibzib - dzibyid;
                const double dxicxic = -dxiaxic - dxibxic - dxicxid;
                const double dxicyic = -dyiaxic - dyibxic - dxicyid;
                const double dxiczic = -dziaxic - dzibxic - dxiczid;
                const double dyicyic = -dyiayic - dyibyic - dyicyid;
                const double dyiczic = -dziayic - dzibyic - dyiczid;
                const double dziczic = -dziazic - dzibzic - dziczid;

                // increment diagonal and off-diagonal Hessian elements
                const std::size_t ia = torsion.atoms[0];
                const std::size_t ib = torsion.atoms[1];
                const std::size_t ic = torsion.atoms[2];
                const std::size_t id = torsion.atoms[3];

                //else if (i.eq.ia) then
                *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxia + d2edphi2 * dphidxia * dphidxia;
                *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayia + d2edphi2 * dphidxia * dphidyia;
                *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazia + d2edphi2 * dphidxia * dphidzia;
                *hessx[2][ia] = *hessx[2][ia] + dedphi * dxiayia + d2edphi2 * dphidxia * dphidyia;
                *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayia + d2edphi2 * dphidyia * dphidyia;
                *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazia + d2edphi2 * dphidyia * dphidzia;
                *hessx[3][ia] = *hessx[3][ia] + dedphi * dxiazia + d2edphi2 * dphidxia * dphidzia;
                *hessy[3][ia] = *hessy[3][ia] + dedphi * dyiazia + d2edphi2 * dphidyia * dphidzia;
                *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazia + d2edphi2 * dphidzia * dphidzia;
                *hessx[1][ib] = *hessx[1][ib] + dedphi * dxiaxib + d2edphi2 * dphidxia * dphidxib;
                *hessy[1][ib] = *hessy[1][ib] + dedphi * dyiaxib + d2edphi2 * dphidyia * dphidxib;
                *hessz[1][ib] = *hessz[1][ib] + dedphi * dziaxib + d2edphi2 * dphidzia * dphidxib;
                *hessx[2][ib] = *hessx[2][ib] + dedphi * dxiayib + d2edphi2 * dphidxia * dphidyib;
                *hessy[2][ib] = *hessy[2][ib] + dedphi * dyiayib + d2edphi2 * dphidyia * dphidyib;
                *hessz[2][ib] = *hessz[2][ib] + dedphi * dziayib + d2edphi2 * dphidzia * dphidyib;
                *hessx[3][ib] = *hessx[3][ib] + dedphi * dxiazib + d2edphi2 * dphidxia * dphidzib;
                *hessy[3][ib] = *hessy[3][ib] + dedphi * dyiazib + d2edphi2 * dphidyia * dphidzib;
                *hessz[3][ib] = *hessz[3][ib] + dedphi * dziazib + d2edphi2 * dphidzia * dphidzib;
                *hessx[1][ic] = *hessx[1][ic] + dedphi * dxiaxic + d2edphi2 * dphidxia * dphidxic;
                *hessy[1][ic] = *hessy[1][ic] + dedphi * dyiaxic + d2edphi2 * dphidyia * dphidxic;
                *hessz[1][ic] = *hessz[1][ic] + dedphi * dziaxic + d2edphi2 * dphidzia * dphidxic;
                *hessx[2][ic] = *hessx[2][ic] + dedphi * dxiayic + d2edphi2 * dphidxia * dphidyic;
                *hessy[2][ic] = *hessy[2][ic] + dedphi * dyiayic + d2edphi2 * dphidyia * dphidyic;
                *hessz[2][ic] = *hessz[2][ic] + dedphi * dziayic + d2edphi2 * dphidzia * dphidyic;
                *hessx[3][ic] = *hessx[3][ic] + dedphi * dxiazic + d2edphi2 * dphidxia * dphidzic;
                *hessy[3][ic] = *hessy[3][ic] + dedphi * dyiazic + d2edphi2 * dphidyia * dphidzic;
                *hessz[3][ic] = *hessz[3][ic] + dedphi * dziazic + d2edphi2 * dphidzia * dphidzic;
                *hessx[1][id] = *hessx[1][id] + dedphi * dxiaxid + d2edphi2 * dphidxia * dphidxid;
                *hessy[1][id] = *hessy[1][id] + dedphi * dyiaxid + d2edphi2 * dphidyia * dphidxid;
                *hessz[1][id] = *hessz[1][id] + dedphi * dziaxid + d2edphi2 * dphidzia * dphidxid;
                *hessx[2][id] = *hessx[2][id] + dedphi * dxiayid + d2edphi2 * dphidxia * dphidyid;
                *hessy[2][id] = *hessy[2][id] + dedphi * dyiayid + d2edphi2 * dphidyia * dphidyid;
                *hessz[2][id] = *hessz[2][id] + dedphi * dziayid + d2edphi2 * dphidzia * dphidyid;
                *hessx[3][id] = *hessx[3][id] + dedphi * dxiazid + d2edphi2 * dphidxia * dphidzid;
                *hessy[3][id] = *hessy[3][id] + dedphi * dyiazid + d2edphi2 * dphidyia * dphidzid;
                *hessz[3][id] = *hessz[3][id] + dedphi * dziazid + d2edphi2 * dphidzia * dphidzid;
                //else if (i.eq.ib) then                                                           ;
                *hessx[1][ib] = *hessx[1][ib] + dedphi * dxibxib + d2edphi2 * dphidxib * dphidxib;
                *hessy[1][ib] = *hessy[1][ib] + dedphi * dxibyib + d2edphi2 * dphidxib * dphidyib;
                *hessz[1][ib] = *hessz[1][ib] + dedphi * dxibzib + d2edphi2 * dphidxib * dphidzib;
                *hessx[2][ib] = *hessx[2][ib] + dedphi * dxibyib + d2edphi2 * dphidxib * dphidyib;
                *hessy[2][ib] = *hessy[2][ib] + dedphi * dyibyib + d2edphi2 * dphidyib * dphidyib;
                *hessz[2][ib] = *hessz[2][ib] + dedphi * dyibzib + d2edphi2 * dphidyib * dphidzib;
                *hessx[3][ib] = *hessx[3][ib] + dedphi * dxibzib + d2edphi2 * dphidxib * dphidzib;
                *hessy[3][ib] = *hessy[3][ib] + dedphi * dyibzib + d2edphi2 * dphidyib * dphidzib;
                *hessz[3][ib] = *hessz[3][ib] + dedphi * dzibzib + d2edphi2 * dphidzib * dphidzib;
                *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxib + d2edphi2 * dphidxib * dphidxia;
                *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayib + d2edphi2 * dphidyib * dphidxia;
                *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazib + d2edphi2 * dphidzib * dphidxia;
                *hessx[2][ia] = *hessx[2][ia] + dedphi * dyiaxib + d2edphi2 * dphidxib * dphidyia;
                *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayib + d2edphi2 * dphidyib * dphidyia;
                *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazib + d2edphi2 * dphidzib * dphidyia;
                *hessx[3][ia] = *hessx[3][ia] + dedphi * dziaxib + d2edphi2 * dphidxib * dphidzia;
                *hessy[3][ia] = *hessy[3][ia] + dedphi * dziayib + d2edphi2 * dphidyib * dphidzia;
                *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazib + d2edphi2 * dphidzib * dphidzia;
                *hessx[1][ic] = *hessx[1][ic] + dedphi * dxibxic + d2edphi2 * dphidxib * dphidxic;
                *hessy[1][ic] = *hessy[1][ic] + dedphi * dyibxic + d2edphi2 * dphidyib * dphidxic;
                *hessz[1][ic] = *hessz[1][ic] + dedphi * dzibxic + d2edphi2 * dphidzib * dphidxic;
                *hessx[2][ic] = *hessx[2][ic] + dedphi * dxibyic + d2edphi2 * dphidxib * dphidyic;
                *hessy[2][ic] = *hessy[2][ic] + dedphi * dyibyic + d2edphi2 * dphidyib * dphidyic;
                *hessz[2][ic] = *hessz[2][ic] + dedphi * dzibyic + d2edphi2 * dphidzib * dphidyic;
                *hessx[3][ic] = *hessx[3][ic] + dedphi * dxibzic + d2edphi2 * dphidxib * dphidzic;
                *hessy[3][ic] = *hessy[3][ic] + dedphi * dyibzic + d2edphi2 * dphidyib * dphidzic;
                *hessz[3][ic] = *hessz[3][ic] + dedphi * dzibzic + d2edphi2 * dphidzib * dphidzic;
                *hessx[1][id] = *hessx[1][id] + dedphi * dxibxid + d2edphi2 * dphidxib * dphidxid;
                *hessy[1][id] = *hessy[1][id] + dedphi * dyibxid + d2edphi2 * dphidyib * dphidxid;
                *hessz[1][id] = *hessz[1][id] + dedphi * dzibxid + d2edphi2 * dphidzib * dphidxid;
                *hessx[2][id] = *hessx[2][id] + dedphi * dxibyid + d2edphi2 * dphidxib * dphidyid;
                *hessy[2][id] = *hessy[2][id] + dedphi * dyibyid + d2edphi2 * dphidyib * dphidyid;
                *hessz[2][id] = *hessz[2][id] + dedphi * dzibyid + d2edphi2 * dphidzib * dphidyid;
                *hessx[3][id] = *hessx[3][id] + dedphi * dxibzid + d2edphi2 * dphidxib * dphidzid;
                *hessy[3][id] = *hessy[3][id] + dedphi * dyibzid + d2edphi2 * dphidyib * dphidzid;
                *hessz[3][id] = *hessz[3][id] + dedphi * dzibzid + d2edphi2 * dphidzib * dphidzid;
                //else if (i.eq.ic) then                                                           ;
                *hessx[1][ic] = *hessx[1][ic] + dedphi * dxicxic + d2edphi2 * dphidxic * dphidxic;
                *hessy[1][ic] = *hessy[1][ic] + dedphi * dxicyic + d2edphi2 * dphidxic * dphidyic;
                *hessz[1][ic] = *hessz[1][ic] + dedphi * dxiczic + d2edphi2 * dphidxic * dphidzic;
                *hessx[2][ic] = *hessx[2][ic] + dedphi * dxicyic + d2edphi2 * dphidxic * dphidyic;
                *hessy[2][ic] = *hessy[2][ic] + dedphi * dyicyic + d2edphi2 * dphidyic * dphidyic;
                *hessz[2][ic] = *hessz[2][ic] + dedphi * dyiczic + d2edphi2 * dphidyic * dphidzic;
                *hessx[3][ic] = *hessx[3][ic] + dedphi * dxiczic + d2edphi2 * dphidxic * dphidzic;
                *hessy[3][ic] = *hessy[3][ic] + dedphi * dyiczic + d2edphi2 * dphidyic * dphidzic;
                *hessz[3][ic] = *hessz[3][ic] + dedphi * dziczic + d2edphi2 * dphidzic * dphidzic;
                *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxic + d2edphi2 * dphidxic * dphidxia;
                *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayic + d2edphi2 * dphidyic * dphidxia;
                *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazic + d2edphi2 * dphidzic * dphidxia;
                *hessx[2][ia] = *hessx[2][ia] + dedphi * dyiaxic + d2edphi2 * dphidxic * dphidyia;
                *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayic + d2edphi2 * dphidyic * dphidyia;
                *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazic + d2edphi2 * dphidzic * dphidyia;
                *hessx[3][ia] = *hessx[3][ia] + dedphi * dziaxic + d2edphi2 * dphidxic * dphidzia;
                *hessy[3][ia] = *hessy[3][ia] + dedphi * dziayic + d2edphi2 * dphidyic * dphidzia;
                *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazic + d2edphi2 * dphidzic * dphidzia;
                *hessx[1][ib] = *hessx[1][ib] + dedphi * dxibxic + d2edphi2 * dphidxic * dphidxib;
                *hessy[1][ib] = *hessy[1][ib] + dedphi * dxibyic + d2edphi2 * dphidyic * dphidxib;
                *hessz[1][ib] = *hessz[1][ib] + dedphi * dxibzic + d2edphi2 * dphidzic * dphidxib;
                *hessx[2][ib] = *hessx[2][ib] + dedphi * dyibxic + d2edphi2 * dphidxic * dphidyib;
                *hessy[2][ib] = *hessy[2][ib] + dedphi * dyibyic + d2edphi2 * dphidyic * dphidyib;
                *hessz[2][ib] = *hessz[2][ib] + dedphi * dyibzic + d2edphi2 * dphidzic * dphidyib;
                *hessx[3][ib] = *hessx[3][ib] + dedphi * dzibxic + d2edphi2 * dphidxic * dphidzib;
                *hessy[3][ib] = *hessy[3][ib] + dedphi * dzibyic + d2edphi2 * dphidyic * dphidzib;
                *hessz[3][ib] = *hessz[3][ib] + dedphi * dzibzic + d2edphi2 * dphidzic * dphidzib;
                *hessx[1][id] = *hessx[1][id] + dedphi * dxicxid + d2edphi2 * dphidxic * dphidxid;
                *hessy[1][id] = *hessy[1][id] + dedphi * dyicxid + d2edphi2 * dphidyic * dphidxid;
                *hessz[1][id] = *hessz[1][id] + dedphi * dzicxid + d2edphi2 * dphidzic * dphidxid;
                *hessx[2][id] = *hessx[2][id] + dedphi * dxicyid + d2edphi2 * dphidxic * dphidyid;
                *hessy[2][id] = *hessy[2][id] + dedphi * dyicyid + d2edphi2 * dphidyic * dphidyid;
                *hessz[2][id] = *hessz[2][id] + dedphi * dzicyid + d2edphi2 * dphidzic * dphidyid;
                *hessx[3][id] = *hessx[3][id] + dedphi * dxiczid + d2edphi2 * dphidxic * dphidzid;
                *hessy[3][id] = *hessy[3][id] + dedphi * dyiczid + d2edphi2 * dphidyic * dphidzid;
                *hessz[3][id] = *hessz[3][id] + dedphi * dziczid + d2edphi2 * dphidzic * dphidzid;
                //else if (i.eq.id) then                                                           ;
                *hessx[1][id] = *hessx[1][id] + dedphi * dxidxid + d2edphi2 * dphidxid * dphidxid;
                *hessy[1][id] = *hessy[1][id] + dedphi * dxidyid + d2edphi2 * dphidxid * dphidyid;
                *hessz[1][id] = *hessz[1][id] + dedphi * dxidzid + d2edphi2 * dphidxid * dphidzid;
                *hessx[2][id] = *hessx[2][id] + dedphi * dxidyid + d2edphi2 * dphidxid * dphidyid;
                *hessy[2][id] = *hessy[2][id] + dedphi * dyidyid + d2edphi2 * dphidyid * dphidyid;
                *hessz[2][id] = *hessz[2][id] + dedphi * dyidzid + d2edphi2 * dphidyid * dphidzid;
                *hessx[3][id] = *hessx[3][id] + dedphi * dxidzid + d2edphi2 * dphidxid * dphidzid;
                *hessy[3][id] = *hessy[3][id] + dedphi * dyidzid + d2edphi2 * dphidyid * dphidzid;
                *hessz[3][id] = *hessz[3][id] + dedphi * dzidzid + d2edphi2 * dphidzid * dphidzid;
                *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxid + d2edphi2 * dphidxid * dphidxia;
                *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayid + d2edphi2 * dphidyid * dphidxia;
                *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazid + d2edphi2 * dphidzid * dphidxia;
                *hessx[2][ia] = *hessx[2][ia] + dedphi * dyiaxid + d2edphi2 * dphidxid * dphidyia;
                *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayid + d2edphi2 * dphidyid * dphidyia;
                *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazid + d2edphi2 * dphidzid * dphidyia;
                *hessx[3][ia] = *hessx[3][ia] + dedphi * dziaxid + d2edphi2 * dphidxid * dphidzia;
                *hessy[3][ia] = *hessy[3][ia] + dedphi * dziayid + d2edphi2 * dphidyid * dphidzia;
                *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazid + d2edphi2 * dphidzid * dphidzia;
                *hessx[1][ib] = *hessx[1][ib] + dedphi * dxibxid + d2edphi2 * dphidxid * dphidxib;
                *hessy[1][ib] = *hessy[1][ib] + dedphi * dxibyid + d2edphi2 * dphidyid * dphidxib;
                *hessz[1][ib] = *hessz[1][ib] + dedphi * dxibzid + d2edphi2 * dphidzid * dphidxib;
                *hessx[2][ib] = *hessx[2][ib] + dedphi * dyibxid + d2edphi2 * dphidxid * dphidyib;
                *hessy[2][ib] = *hessy[2][ib] + dedphi * dyibyid + d2edphi2 * dphidyid * dphidyib;
                *hessz[2][ib] = *hessz[2][ib] + dedphi * dyibzid + d2edphi2 * dphidzid * dphidyib;
                *hessx[3][ib] = *hessx[3][ib] + dedphi * dzibxid + d2edphi2 * dphidxid * dphidzib;
                *hessy[3][ib] = *hessy[3][ib] + dedphi * dzibyid + d2edphi2 * dphidyid * dphidzib;
                *hessz[3][ib] = *hessz[3][ib] + dedphi * dzibzid + d2edphi2 * dphidzid * dphidzib;
                *hessx[1][ic] = *hessx[1][ic] + dedphi * dxicxid + d2edphi2 * dphidxid * dphidxic;
                *hessy[1][ic] = *hessy[1][ic] + dedphi * dxicyid + d2edphi2 * dphidyid * dphidxic;
                *hessz[1][ic] = *hessz[1][ic] + dedphi * dxiczid + d2edphi2 * dphidzid * dphidxic;
                *hessx[2][ic] = *hessx[2][ic] + dedphi * dyicxid + d2edphi2 * dphidxid * dphidyic;
                *hessy[2][ic] = *hessy[2][ic] + dedphi * dyicyid + d2edphi2 * dphidyid * dphidyic;
                *hessz[2][ic] = *hessz[2][ic] + dedphi * dyiczid + d2edphi2 * dphidzid * dphidyic;
                *hessx[3][ic] = *hessx[3][ic] + dedphi * dzicxid + d2edphi2 * dphidxid * dphidzic;
                *hessy[3][ic] = *hessy[3][ic] + dedphi * dzicyid + d2edphi2 * dphidyid * dphidzic;
                *hessz[3][ic] = *hessz[3][ic] + dedphi * dziczid + d2edphi2 * dphidzid * dphidzic;
              }
            }
          }
        
        return E;
      }
      /********************************
      *                                *
      *  Improper                      *
      *  Torsion  Potential            *
      *  Energy/Gradients/Hessians     *
      *                                *
      *                                *
      *********************************/

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_it<0>(void)
      {
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt;
        using std::abs;
        coords::float_type E(0.0);
        for (auto& imptor : refined.imptors())   //for every improper torsion
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
          coords::float_type const tlul = sqrt(tl2 * ul2);
          coords::float_type const r12 = len(cb);
          coords::Cartesian_Point const tu = cross(t, u);

          coords::float_type const cosine = dot(t, u) / tlul;
          coords::float_type const sine = dot(cb, tu) / (r12 * tlul);

          auto v2 = imptor.p.force[0];   // I don't know if this is correct (originally there was something more in this line)
          auto c2 = cos(imptor.p.ideal[0] * SCON_PI180);  // why index 0 everywhere?
          auto s2 = sin(imptor.p.ideal[0] * SCON_PI180);

          auto cosine2 = cosine * cosine - sine * sine;
          auto sine2 = 2.0 * cosine * sine;

          auto phi2 = 1.0 + (cosine2 * c2 + sine2 * s2);

          double E_add = 0.5 * (v2 * phi2);
          E += E_add;
        }
        return E;
      }

      coords::float_type energy::interfaces::aco::aco_ff::f_imptor_1_calc_one_torsion(::tinker::refine::types::imptor const& imptor)
      {
        double E = 0.;
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
        coords::float_type const tlul = sqrt(tl2 * ul2);
        coords::float_type const r12 = len(cb);
        coords::Cartesian_Point const tu = cross(t, u);
      
        coords::float_type const cosine = dot(t, u) / tlul;
        coords::float_type const sine = dot(cb, tu) / (r12 * tlul);
      
        auto v2 = imptor.p.force[0];   // I don't know if this is correct (originally there was something more in this line)
        auto c2 = cos(imptor.p.ideal[0] * SCON_PI180);  // why index 0 everywhere?
        auto s2 = sin(imptor.p.ideal[0] * SCON_PI180);
      
        auto cosine2 = cosine * cosine - sine * sine;
        auto sine2 = 2.0 * cosine * sine;
      
        auto phi2 = 1.0 + (cosine2 * c2 + sine2 * s2);
        auto dphi2 = 2.0 * (cosine2 * s2 - sine2 * c2);
      
        auto dedphi = 0.5 * (v2 * dphi2);
      
        double E_add = 0.5 * (v2 * phi2);
        E += E_add;
      
        //gradient calculation
        auto const ca(coords->xyz(imptor.center) - coords->xyz(imptor.ligand[0]));
        auto const db(coords->xyz(imptor.twist) - coords->xyz(imptor.ligand[1]));
      
        auto const dt = cross(t, cb) * (dedphi / (tl2 * r12));
        auto const du = cross(u, cb) * (-dedphi / (ul2 * r12));
      
        auto vir1 = cross(dt, ba) + cross(db, du);
        auto vir2 = cross(dt, cb);
        auto vir3 = cross(ca, dt) + cross(du, dc);
        auto vir4 = cross(du, cb);
      
        part_grad[IMPTORSION][imptor.ligand[0]] += vir2;
        part_grad[IMPTORSION][imptor.ligand[1]] += vir3;
        part_grad[IMPTORSION][imptor.center] += vir1;
        part_grad[IMPTORSION][imptor.twist] += vir4;
      
        //calculation of virial tensors (copied from function f_imp<1>)
        auto vxx = cb.x() * (vir3.x() + vir4.x()) - ba.x() * vir1.x() + dc.x() * vir4.x();
        auto vyx = cb.y() * (vir3.x() + vir4.x()) - ba.y() * vir1.x() + dc.y() * vir4.x();
        auto vzx = cb.z() * (vir3.x() + vir4.x()) - ba.z() * vir1.x() + dc.z() * vir4.x();
        auto vyy = cb.y() * (vir3.y() + vir4.y()) - ba.y() * vir1.y() + dc.y() * vir4.y();
        auto vzy = cb.z() * (vir3.y() + vir4.y()) - ba.z() * vir1.y() + dc.z() * vir4.y();
        auto vzz = cb.z() * (vir3.z() + vir4.z()) - ba.z() * vir1.z() + dc.z() * vir4.z();
      
        part_virial[IMPROPER][0][0] += vxx;
        part_virial[IMPROPER][1][0] += vyx;
        part_virial[IMPROPER][2][0] += vzx;
        part_virial[IMPROPER][0][1] += vyx;
        part_virial[IMPROPER][1][1] += vyy;
        part_virial[IMPROPER][2][1] += vzy;
        part_virial[IMPROPER][0][2] += vzx;
        part_virial[IMPROPER][1][2] += vzy;
        part_virial[IMPROPER][2][2] += vzz;
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_it<1>(void)
      {
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt;
        using std::abs;
        coords::float_type E(0.0);
        for (auto& imptor : refined.imptors())   //for every improper torsion
        {
          E += f_imptor_1_calc_one_torsion(imptor);
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_it<2>(void)
      {
        std::vector < std::vector<double*>> hessx(4u, std::vector<double*>(unsigned int(coords->atoms().size()), nullptr));
        std::vector < std::vector<double*>> hessy(hessx), hessz(hessx);
        for (std::size_t i = 1u; i <= 3u; ++i)
        {
          for (std::size_t atom = 0u; atom < coords->atoms().size(); ++atom)
          {
            std::size_t param = atom * 3 + (i - 1u);
            hessx[i][atom] = &part_hessian[TORSION][atom * 3][atom * 3];
            hessy[i][atom] = &part_hessian[TORSION][atom * 3][atom * 3 + 1u];
            hessz[i][atom] = &part_hessian[TORSION][atom * 3][atom * 3 + 2u];
          }
        }

        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt;
        using std::abs;
        coords::float_type E(0.0);
        const double itorunit = cparams.imptorunit();
        for (auto& imptor : refined.imptors())   //for every improper torsion
        {
          const double xia = coords->xyz(imptor.ligand[0]).x();
          const double yia = coords->xyz(imptor.ligand[0]).y();
          const double zia = coords->xyz(imptor.ligand[0]).z();
          const double xib = coords->xyz(imptor.ligand[1]).x();
          const double yib = coords->xyz(imptor.ligand[1]).y();
          const double zib = coords->xyz(imptor.ligand[1]).z();
          const double xic = coords->xyz(imptor.center).x();
          const double yic = coords->xyz(imptor.center).y();
          const double zic = coords->xyz(imptor.center).z();
          const double xid = coords->xyz(imptor.twist).x();
          const double yid = coords->xyz(imptor.twist).y();
          const double zid = coords->xyz(imptor.twist).z();
          const double xba = xib - xia                      ;
          const double yba = yib - yia                      ;
          const double zba = zib - zia                      ;
          const double xcb = xic - xib                      ;
          const double ycb = yic - yib                      ;
          const double zcb = zic - zib                      ;
          const double xdc = xid - xic                      ;
          const double ydc = yid - yic                      ;
          const double zdc = zid - zic                      ;
          const double xt = yba * zcb - ycb * zba           ;
          const double yt = zba * xcb - zcb * xba           ;
          const double zt = xba * ycb - xcb * yba           ;
          const double xu = ycb * zdc - ydc * zcb           ;
          const double yu = zcb * xdc - zdc * xcb           ;
          const double zu = xcb * ydc - xdc * ycb           ;
          const double xtu = yt * zu - yu * zt              ;
          const double ytu = zt * xu - zu * xt              ;
          const double ztu = xt * yu - xu * yt              ;
          const double rt2 = xt * xt + yt * yt + zt * zt    ;
          const double ru2 = xu * xu + yu * yu + zu * zu    ;
          const double rtru = sqrt(rt2 * ru2)               ;
          if (rtru != 0.)
          {
            coords::float_type cos[3], sin[3];
            const double rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb)            ;
            cos[0] = (xt * xu + yt * yu + zt * zu) / rtru            ;
            sin[0] = (xcb * xtu + ycb * ytu + zcb * ztu) / (rcb * rtru);
              //set the improper torsional parameters for this angle
            const double& cosine = cos[0];
            const double& sine = sin[0];
            cos[1] = cosine * cosine - sine * sine;
            const double& cosine2 = cos[1];
            sin[1] = 2. * cosine * sine;
            const double& sine2 = sin[1];
            cos[2] = cosine * cosine2 - sine * sine2;
            const double& cosine3 = cos[2];
            sin[2] = cosine * sine2 + sine * cosine2;
            const double sine3 = sin[2];


            // We only suppprt CHARMM / AMBER / OPLSAA, thus only order 2 params are important
            const double v1 = 0.  ;
            const double   c1 = 0.;
            const double   s1 = 0.;
            const double   v2 = imptor.p.force[0];
            const double   c2 = std::cos(imptor.p.ideal[0] * SCON_PI180);
            const double   s2 = std::sin(imptor.p.ideal[0] * SCON_PI180);
            const double   v3 = 0.;
            const double   c3 = 0.;
            const double   s3 = 0.;
            const double dphi1 = (cosine * s1 - sine * c1)              ;
            const double   dphi2 = 2.0 * (cosine2 * s2 - sine2 * c2)  ;
            const double   dphi3 = 3.0 * (cosine3 * s3 - sine3 * c3)  ;
            const double   d2phi1 = -(cosine * c1 + sine * s1)          ;
            const double   d2phi2 = -4.0 * (cosine2 * c2 + sine2 * s2);
            const double   d2phi3 = -9.0 * (cosine3 * c3 + sine3 * s3);

            // calculate the improper torsion master chain rule terms
            const double dedphi = itorunit * (v1 * dphi1 + v2 * dphi2 + v3 * dphi3);
            const double d2edphi2 = itorunit * (v1 * d2phi1 + v2 * d2phi2 + v3 * d2phi3);


            // chain rule terms for first derivative components
            const double xca = xic - xia;
            const double yca = yic - yia;
            const double zca = zic - zia;
            const double xdb = xid - xib;
            const double ydb = yid - yib;
            const double zdb = zid - zib;
            const double dphidxt = (yt * zcb - ycb * zt) / (rt2 * rcb)   ;
            const double   dphidyt = (zt * xcb - zcb * xt) / (rt2 * rcb) ;
            const double   dphidzt = (xt * ycb - xcb * yt) / (rt2 * rcb) ;
            const double   dphidxu = -(yu * zcb - ycb * zu) / (ru2 * rcb);
            const double   dphidyu = -(zu * xcb - zcb * xu) / (ru2 * rcb);
            const double   dphidzu = -(xu * ycb - xcb * yu) / (ru2 * rcb);

            // abbreviations for second derivative chain rule terms
            const double xycb2 = xcb * xcb + ycb * ycb            ;
            const double xzcb2 = xcb * xcb + zcb * zcb            ;
            const double yzcb2 = ycb * ycb + zcb * zcb            ;
            const double rcbxt = -2.0 * rcb * dphidxt           ;
            const double rcbyt = -2.0 * rcb * dphidyt           ;
            const double rcbzt = -2.0 * rcb * dphidzt           ;
            const double rcbt2 = rcb * rt2                        ;
            const double rcbxu = 2.0 * rcb * dphidxu            ;
            const double rcbyu = 2.0 * rcb * dphidyu            ;
            const double rcbzu = 2.0 * rcb * dphidzu            ;
            const double rcbu2 = rcb * ru2                        ;
            const double dphidxibt = yca * dphidzt - zca * dphidyt;
            const double dphidxibu = zdc * dphidyu - ydc * dphidzu;
            const double dphidyibt = zca * dphidxt - xca * dphidzt;
            const double dphidyibu = xdc * dphidzu - zdc * dphidxu;
            const double dphidzibt = xca * dphidyt - yca * dphidxt;
            const double dphidzibu = ydc * dphidxu - xdc * dphidyu;
            const double dphidxict = zba * dphidyt - yba * dphidzt;
            const double dphidxicu = ydb * dphidzu - zdb * dphidyu;
            const double dphidyict = xba * dphidzt - zba * dphidxt;
            const double dphidyicu = zdb * dphidxu - xdb * dphidzu;
            const double dphidzict = yba * dphidxt - xba * dphidyt;
            const double dphidzicu = xdb * dphidyu - ydb * dphidxu;

            // chain rule terms for first derivative components
            const double dphidxia = zcb * dphidyt - ycb * dphidzt  ;
            const double   dphidyia = xcb * dphidzt - zcb * dphidxt;
            const double   dphidzia = ycb * dphidxt - xcb * dphidyt;
            const double   dphidxib = dphidxibt + dphidxibu        ;
            const double   dphidyib = dphidyibt + dphidyibu        ;
            const double   dphidzib = dphidzibt + dphidzibu        ;
            const double   dphidxic = dphidxict + dphidxicu        ;
            const double   dphidyic = dphidyict + dphidyicu        ;
            const double   dphidzic = dphidzict + dphidzicu        ;
            const double   dphidxid = zcb * dphidyu - ycb * dphidzu;
            const double   dphidyid = xcb * dphidzu - zcb * dphidxu;
            const double   dphidzid = ycb * dphidxu - xcb * dphidyu;

            // chain rule terms for second derivative components
            const double dxiaxia = rcbxt * dphidxia                                                       ;
            const double   dxiayia = rcbxt * dphidyia - zcb * rcb / rt2                                   ;
            const double   dxiazia = rcbxt * dphidzia + ycb * rcb / rt2                                   ;
            const double   dxiaxic = rcbxt * dphidxict + xcb * xt / rcbt2                                 ;
            const double   dxiayic = rcbxt * dphidyict - dphidzt - (xba * zcb * xcb + zba * yzcb2) / rcbt2;
            const double   dxiazic = rcbxt * dphidzict + dphidyt + (xba * ycb * xcb + yba * yzcb2) / rcbt2;
            const double   dxiaxid = 0.0                                                                ;
            const double   dxiayid = 0.0                                                                ;
            const double   dxiazid = 0.0                                                                ;
            const double   dyiayia = rcbyt * dphidyia                                                     ;
            const double   dyiazia = rcbyt * dphidzia - xcb * rcb / rt2                                   ;
            const double   dyiaxib = rcbyt * dphidxibt - dphidzt - (yca * zcb * ycb + zca * xzcb2) / rcbt2;
            const double   dyiaxic = rcbyt * dphidxict + dphidzt + (yba * zcb * ycb + zba * xzcb2) / rcbt2;
            const double   dyiayic = rcbyt * dphidyict + ycb * yt / rcbt2                                 ;
            const double   dyiazic = rcbyt * dphidzict - dphidxt - (yba * xcb * ycb + xba * xzcb2) / rcbt2;
            const double   dyiaxid = 0.0                                                                ;
            const double   dyiayid = 0.0                                                                ;
            const double   dyiazid = 0.0                                                                ;
            const double   dziazia = rcbzt * dphidzia                                                     ;
            const double   dziaxib = rcbzt * dphidxibt + dphidyt + (zca * ycb * zcb + yca * xycb2) / rcbt2;
            const double   dziayib = rcbzt * dphidyibt - dphidxt - (zca * xcb * zcb + xca * xycb2) / rcbt2;
            const double   dziaxic = rcbzt * dphidxict - dphidyt - (zba * ycb * zcb + yba * xycb2) / rcbt2;
            const double   dziayic = rcbzt * dphidyict + dphidxt + (zba * xcb * zcb + xba * xycb2) / rcbt2;
            const double   dziazic = rcbzt * dphidzict + zcb * zt / rcbt2;
            const double   dziaxid = 0.0;
            const double   dziayid = 0.0;
            const double   dziazid = 0.0;
            const double   dxibxic = -xcb * dphidxib / (rcb * rcb) - (yca * (zba * xcb + yt) - zca * (yba * xcb - zt)) / rcbt2 - 2.0 * (yt * zba - yba * zt) * dphidxibt / rt2 - (zdc * (ydb * xcb + zu) - ydc * (zdb * xcb - yu)) / rcbu2 + 2.0 * (yu * zdb - ydb * zu) * dphidxibu / ru2;
            const double   dxibyic = -ycb * dphidxib / (rcb * rcb) + dphidzt + dphidzu - (yca * (zba * ycb - xt) + zca * (xba * xcb + zcb * zba)) / rcbt2 - 2.0 * (zt * xba - zba * xt) * dphidxibt / rt2 + (zdc * (xdb * xcb + zcb * zdb) + ydc * (zdb * ycb + xu)) / rcbu2 + 2.0 * (zu * xdb - zdb * xu) * dphidxibu / ru2;
            const double   dxibxid = rcbxu * dphidxibu + xcb * xu / rcbu2;
            const double   dxibyid = rcbyu * dphidxibu - dphidzu - (ydc * zcb * ycb + zdc * xzcb2) / rcbu2;
            const double   dxibzid = rcbzu * dphidxibu + dphidyu + (zdc * ycb * zcb + ydc * xycb2) / rcbu2;
            const double   dyibzib = ycb * dphidzib / (rcb * rcb) - (xca * (xca * xcb + zcb * zca) + yca * (ycb * xca + zt)) / rcbt2 - 2.0 * (xt * zca - xca * zt) * dphidzibt / rt2 + (ydc * (xdc * ycb - zu) + xdc * (xdc * xcb + zcb * zdc)) / rcbu2 + 2.0 * (xu * zdc - xdc * zu) * dphidzibu / ru2;
            const double   dyibxic = -xcb * dphidyib / (rcb * rcb) - dphidzt - dphidzu + (xca * (zba * xcb + yt) + zca * (zba * zcb + ycb * yba)) / rcbt2 - 2.0 * (yt * zba - yba * zt) * dphidyibt / rt2 - (zdc * (zdb * zcb + ycb * ydb) + xdc * (zdb * xcb - yu)) / rcbu2 + 2.0 * (yu * zdb - ydb * zu) * dphidyibu / ru2;
            const double   dyibyic = -ycb * dphidyib / (rcb * rcb) - (zca * (xba * ycb + zt) - xca * (zba * ycb - xt)) / rcbt2 - 2.0 * (zt * xba - zba * xt) * dphidyibt / rt2 - (xdc * (zdb * ycb + xu) - zdc * (xdb * ycb - zu)) / rcbu2 + 2.0 * (zu * xdb - zdb * xu) * dphidyibu / ru2;
            const double   dyibxid = rcbxu * dphidyibu + dphidzu + (xdc * zcb * xcb + zdc * yzcb2) / rcbu2;
            const double   dyibyid = rcbyu * dphidyibu + ycb * yu / rcbu2;
            const double   dyibzid = rcbzu * dphidyibu - dphidxu - (zdc * xcb * zcb + xdc * xycb2) / rcbu2;
            const double   dzibxic = -xcb * dphidzib / (rcb * rcb) + dphidyt + dphidyu - (xca * (yba * xcb - zt) + yca * (zba * zcb + ycb * yba)) / rcbt2 - 2.0 * (yt * zba - yba * zt) * dphidzibt / rt2 + (ydc * (zdb * zcb + ycb * ydb) + xdc * (ydb * xcb + zu)) / rcbu2 + 2.0 * (yu * zdb - ydb * zu) * dphidzibu / ru2;
            const double   dzibzic = -zcb * dphidzib / (rcb * rcb) - (xca * (yba * zcb + xt) - yca * (xba * zcb - yt)) / rcbt2 - 2.0 * (xt * yba - xba * yt) * dphidzibt / rt2 - (ydc * (xdb * zcb + yu) - xdc * (ydb * zcb - xu)) / rcbu2 + 2.0 * (xu * ydb - xdb * yu) * dphidzibu / ru2;
            const double   dzibxid = rcbxu * dphidzibu - dphidyu - (xdc * ycb * xcb + ydc * yzcb2) / rcbu2;
            const double   dzibyid = rcbyu * dphidzibu + dphidxu + (ydc * xcb * ycb + xdc * xzcb2) / rcbu2;
            const double   dzibzid = rcbzu * dphidzibu + zcb * zu / rcbu2                                 ;
            const double   dxicxid = rcbxu * dphidxicu - xcb * (zdb * ycb - ydb * zcb) / rcbu2            ;
            const double   dxicyid = rcbyu * dphidxicu + dphidzu + (ydb * zcb * ycb + zdb * xzcb2) / rcbu2;
            const double   dxiczid = rcbzu * dphidxicu - dphidyu - (zdb * ycb * zcb + ydb * xycb2) / rcbu2;
            const double   dyicxid = rcbxu * dphidyicu - dphidzu - (xdb * zcb * xcb + zdb * yzcb2) / rcbu2;
            const double   dyicyid = rcbyu * dphidyicu - ycb * (xdb * zcb - zdb * xcb) / rcbu2            ;
            const double   dyiczid = rcbzu * dphidyicu + dphidxu + (zdb * xcb * zcb + xdb * xycb2) / rcbu2;
            const double   dzicxid = rcbxu * dphidzicu + dphidyu + (xdb * ycb * xcb + ydb * yzcb2) / rcbu2;
            const double   dzicyid = rcbyu * dphidzicu - dphidxu - (ydb * xcb * ycb + xdb * xzcb2) / rcbu2;
            const double   dziczid = rcbzu * dphidzicu - zcb * (ydb * xcb - xdb * ycb) / rcbu2            ;
            const double   dxidxid = rcbxu * dphidxid                                                     ;
            const double   dxidyid = rcbxu * dphidyid + zcb * rcb / ru2                                   ;
            const double   dxidzid = rcbxu * dphidzid - ycb * rcb / ru2                                   ;
            const double   dyidyid = rcbyu * dphidyid                                                     ;
            const double   dyidzid = rcbyu * dphidzid + xcb * rcb / ru2                                   ;
            const double   dzidzid = rcbzu * dphidzid                                                     ;

            // get some second derivative chain rule terms by difference
            const double dxiaxib = -dxiaxia - dxiaxic - dxiaxid  ;
            const double   dxiayib = -dxiayia - dxiayic - dxiayid;
            const double   dxiazib = -dxiazia - dxiazic - dxiazid;
            const double   dyiayib = -dyiayia - dyiayic - dyiayid;
            const double   dyiazib = -dyiazia - dyiazic - dyiazid;
            const double   dziazib = -dziazia - dziazic - dziazid;
            const double   dxibxib = -dxiaxib - dxibxic - dxibxid;
            const double   dxibyib = -dyiaxib - dxibyic - dxibyid;
            const double   dxibzib = -dxiazib - dzibxic - dzibxid;
            const double   dxibzic = -dziaxib - dxibzib - dxibzid;
            const double   dyibyib = -dyiayib - dyibyic - dyibyid;
            const double   dyibzic = -dziayib - dyibzib - dyibzid;
            const double   dzibzib = -dziazib - dzibzic - dzibzid;
            const double   dzibyic = -dyiazib - dyibzib - dzibyid;
            const double   dxicxic = -dxiaxic - dxibxic - dxicxid;
            const double   dxicyic = -dyiaxic - dyibxic - dxicyid;
            const double   dxiczic = -dziaxic - dzibxic - dxiczid;
            const double   dyicyic = -dyiayic - dyibyic - dyicyid;
            const double   dyiczic = -dziayic - dzibyic - dyiczid;
            const double   dziczic = -dziazic - dzibzic - dziczid;

            // increment diagonal and off-diagonal Hessian elements
            const std::size_t ia = imptor.ligand[0];
            const std::size_t ib = imptor.ligand[1];
            const std::size_t ic = imptor.center;
            const std::size_t id = imptor.twist;
            //if (i .eq. ia) then;
            *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxia + d2edphi2 * dphidxia * dphidxia;
            *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayia + d2edphi2 * dphidxia * dphidyia;
            *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazia + d2edphi2 * dphidxia * dphidzia;
            *hessx[2][ia] = *hessx[2][ia] + dedphi * dxiayia + d2edphi2 * dphidxia * dphidyia;
            *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayia + d2edphi2 * dphidyia * dphidyia;
            *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazia + d2edphi2 * dphidyia * dphidzia;
            *hessx[3][ia] = *hessx[3][ia] + dedphi * dxiazia + d2edphi2 * dphidxia * dphidzia;
            *hessy[3][ia] = *hessy[3][ia] + dedphi * dyiazia + d2edphi2 * dphidyia * dphidzia;
            *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazia + d2edphi2 * dphidzia * dphidzia;
            *hessx[1][ib] = *hessx[1][ib] + dedphi * dxiaxib + d2edphi2 * dphidxia * dphidxib;
            *hessy[1][ib] = *hessy[1][ib] + dedphi * dyiaxib + d2edphi2 * dphidyia * dphidxib;
            *hessz[1][ib] = *hessz[1][ib] + dedphi * dziaxib + d2edphi2 * dphidzia * dphidxib;
            *hessx[2][ib] = *hessx[2][ib] + dedphi * dxiayib + d2edphi2 * dphidxia * dphidyib;
            *hessy[2][ib] = *hessy[2][ib] + dedphi * dyiayib + d2edphi2 * dphidyia * dphidyib;
            *hessz[2][ib] = *hessz[2][ib] + dedphi * dziayib + d2edphi2 * dphidzia * dphidyib;
            *hessx[3][ib] = *hessx[3][ib] + dedphi * dxiazib + d2edphi2 * dphidxia * dphidzib;
            *hessy[3][ib] = *hessy[3][ib] + dedphi * dyiazib + d2edphi2 * dphidyia * dphidzib;
            *hessz[3][ib] = *hessz[3][ib] + dedphi * dziazib + d2edphi2 * dphidzia * dphidzib;
            *hessx[1][ic] = *hessx[1][ic] + dedphi * dxiaxic + d2edphi2 * dphidxia * dphidxic;
            *hessy[1][ic] = *hessy[1][ic] + dedphi * dyiaxic + d2edphi2 * dphidyia * dphidxic;
            *hessz[1][ic] = *hessz[1][ic] + dedphi * dziaxic + d2edphi2 * dphidzia * dphidxic;
            *hessx[2][ic] = *hessx[2][ic] + dedphi * dxiayic + d2edphi2 * dphidxia * dphidyic;
            *hessy[2][ic] = *hessy[2][ic] + dedphi * dyiayic + d2edphi2 * dphidyia * dphidyic;
            *hessz[2][ic] = *hessz[2][ic] + dedphi * dziayic + d2edphi2 * dphidzia * dphidyic;
            *hessx[3][ic] = *hessx[3][ic] + dedphi * dxiazic + d2edphi2 * dphidxia * dphidzic;
            *hessy[3][ic] = *hessy[3][ic] + dedphi * dyiazic + d2edphi2 * dphidyia * dphidzic;
            *hessz[3][ic] = *hessz[3][ic] + dedphi * dziazic + d2edphi2 * dphidzia * dphidzic;
            *hessx[1][id] = *hessx[1][id] + dedphi * dxiaxid + d2edphi2 * dphidxia * dphidxid;
            *hessy[1][id] = *hessy[1][id] + dedphi * dyiaxid + d2edphi2 * dphidyia * dphidxid;
            *hessz[1][id] = *hessz[1][id] + dedphi * dziaxid + d2edphi2 * dphidzia * dphidxid;
            *hessx[2][id] = *hessx[2][id] + dedphi * dxiayid + d2edphi2 * dphidxia * dphidyid;
            *hessy[2][id] = *hessy[2][id] + dedphi * dyiayid + d2edphi2 * dphidyia * dphidyid;
            *hessz[2][id] = *hessz[2][id] + dedphi * dziayid + d2edphi2 * dphidzia * dphidyid;
            *hessx[3][id] = *hessx[3][id] + dedphi * dxiazid + d2edphi2 * dphidxia * dphidzid;
            *hessy[3][id] = *hessy[3][id] + dedphi * dyiazid + d2edphi2 * dphidyia * dphidzid;
            *hessz[3][id] = *hessz[3][id] + dedphi * dziazid + d2edphi2 * dphidzia * dphidzid;
            //else if (i .eq. ib) then;
            *hessx[1][ib] = *hessx[1][ib] + dedphi * dxibxib + d2edphi2 * dphidxib * dphidxib;
            *hessy[1][ib] = *hessy[1][ib] + dedphi * dxibyib + d2edphi2 * dphidxib * dphidyib;
            *hessz[1][ib] = *hessz[1][ib] + dedphi * dxibzib + d2edphi2 * dphidxib * dphidzib;
            *hessx[2][ib] = *hessx[2][ib] + dedphi * dxibyib + d2edphi2 * dphidxib * dphidyib;
            *hessy[2][ib] = *hessy[2][ib] + dedphi * dyibyib + d2edphi2 * dphidyib * dphidyib;
            *hessz[2][ib] = *hessz[2][ib] + dedphi * dyibzib + d2edphi2 * dphidyib * dphidzib;
            *hessx[3][ib] = *hessx[3][ib] + dedphi * dxibzib + d2edphi2 * dphidxib * dphidzib;
            *hessy[3][ib] = *hessy[3][ib] + dedphi * dyibzib + d2edphi2 * dphidyib * dphidzib;
            *hessz[3][ib] = *hessz[3][ib] + dedphi * dzibzib + d2edphi2 * dphidzib * dphidzib;
            *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxib + d2edphi2 * dphidxib * dphidxia;
            *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayib + d2edphi2 * dphidyib * dphidxia;
            *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazib + d2edphi2 * dphidzib * dphidxia;
            *hessx[2][ia] = *hessx[2][ia] + dedphi * dyiaxib + d2edphi2 * dphidxib * dphidyia;
            *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayib + d2edphi2 * dphidyib * dphidyia;
            *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazib + d2edphi2 * dphidzib * dphidyia;
            *hessx[3][ia] = *hessx[3][ia] + dedphi * dziaxib + d2edphi2 * dphidxib * dphidzia;
            *hessy[3][ia] = *hessy[3][ia] + dedphi * dziayib + d2edphi2 * dphidyib * dphidzia;
            *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazib + d2edphi2 * dphidzib * dphidzia;
            *hessx[1][ic] = *hessx[1][ic] + dedphi * dxibxic + d2edphi2 * dphidxib * dphidxic;
            *hessy[1][ic] = *hessy[1][ic] + dedphi * dyibxic + d2edphi2 * dphidyib * dphidxic;
            *hessz[1][ic] = *hessz[1][ic] + dedphi * dzibxic + d2edphi2 * dphidzib * dphidxic;
            *hessx[2][ic] = *hessx[2][ic] + dedphi * dxibyic + d2edphi2 * dphidxib * dphidyic;
            *hessy[2][ic] = *hessy[2][ic] + dedphi * dyibyic + d2edphi2 * dphidyib * dphidyic;
            *hessz[2][ic] = *hessz[2][ic] + dedphi * dzibyic + d2edphi2 * dphidzib * dphidyic;
            *hessx[3][ic] = *hessx[3][ic] + dedphi * dxibzic + d2edphi2 * dphidxib * dphidzic;
            *hessy[3][ic] = *hessy[3][ic] + dedphi * dyibzic + d2edphi2 * dphidyib * dphidzic;
            *hessz[3][ic] = *hessz[3][ic] + dedphi * dzibzic + d2edphi2 * dphidzib * dphidzic;
            *hessx[1][id] = *hessx[1][id] + dedphi * dxibxid + d2edphi2 * dphidxib * dphidxid;
            *hessy[1][id] = *hessy[1][id] + dedphi * dyibxid + d2edphi2 * dphidyib * dphidxid;
            *hessz[1][id] = *hessz[1][id] + dedphi * dzibxid + d2edphi2 * dphidzib * dphidxid;
            *hessx[2][id] = *hessx[2][id] + dedphi * dxibyid + d2edphi2 * dphidxib * dphidyid;
            *hessy[2][id] = *hessy[2][id] + dedphi * dyibyid + d2edphi2 * dphidyib * dphidyid;
            *hessz[2][id] = *hessz[2][id] + dedphi * dzibyid + d2edphi2 * dphidzib * dphidyid;
            *hessx[3][id] = *hessx[3][id] + dedphi * dxibzid + d2edphi2 * dphidxib * dphidzid;
            *hessy[3][id] = *hessy[3][id] + dedphi * dyibzid + d2edphi2 * dphidyib * dphidzid;
            *hessz[3][id] = *hessz[3][id] + dedphi * dzibzid + d2edphi2 * dphidzib * dphidzid;
            //else if (i .eq. ic) then;
            *hessx[1][ic] = *hessx[1][ic] + dedphi * dxicxic + d2edphi2 * dphidxic * dphidxic;
            *hessy[1][ic] = *hessy[1][ic] + dedphi * dxicyic + d2edphi2 * dphidxic * dphidyic;
            *hessz[1][ic] = *hessz[1][ic] + dedphi * dxiczic + d2edphi2 * dphidxic * dphidzic;
            *hessx[2][ic] = *hessx[2][ic] + dedphi * dxicyic + d2edphi2 * dphidxic * dphidyic;
            *hessy[2][ic] = *hessy[2][ic] + dedphi * dyicyic + d2edphi2 * dphidyic * dphidyic;
            *hessz[2][ic] = *hessz[2][ic] + dedphi * dyiczic + d2edphi2 * dphidyic * dphidzic;
            *hessx[3][ic] = *hessx[3][ic] + dedphi * dxiczic + d2edphi2 * dphidxic * dphidzic;
            *hessy[3][ic] = *hessy[3][ic] + dedphi * dyiczic + d2edphi2 * dphidyic * dphidzic;
            *hessz[3][ic] = *hessz[3][ic] + dedphi * dziczic + d2edphi2 * dphidzic * dphidzic;
            *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxic + d2edphi2 * dphidxic * dphidxia;
            *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayic + d2edphi2 * dphidyic * dphidxia;
            *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazic + d2edphi2 * dphidzic * dphidxia;
            *hessx[2][ia] = *hessx[2][ia] + dedphi * dyiaxic + d2edphi2 * dphidxic * dphidyia;
            *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayic + d2edphi2 * dphidyic * dphidyia;
            *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazic + d2edphi2 * dphidzic * dphidyia;
            *hessx[3][ia] = *hessx[3][ia] + dedphi * dziaxic + d2edphi2 * dphidxic * dphidzia;
            *hessy[3][ia] = *hessy[3][ia] + dedphi * dziayic + d2edphi2 * dphidyic * dphidzia;
            *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazic + d2edphi2 * dphidzic * dphidzia;
            *hessx[1][ib] = *hessx[1][ib] + dedphi * dxibxic + d2edphi2 * dphidxic * dphidxib;
            *hessy[1][ib] = *hessy[1][ib] + dedphi * dxibyic + d2edphi2 * dphidyic * dphidxib;
            *hessz[1][ib] = *hessz[1][ib] + dedphi * dxibzic + d2edphi2 * dphidzic * dphidxib;
            *hessx[2][ib] = *hessx[2][ib] + dedphi * dyibxic + d2edphi2 * dphidxic * dphidyib;
            *hessy[2][ib] = *hessy[2][ib] + dedphi * dyibyic + d2edphi2 * dphidyic * dphidyib;
            *hessz[2][ib] = *hessz[2][ib] + dedphi * dyibzic + d2edphi2 * dphidzic * dphidyib;
            *hessx[3][ib] = *hessx[3][ib] + dedphi * dzibxic + d2edphi2 * dphidxic * dphidzib;
            *hessy[3][ib] = *hessy[3][ib] + dedphi * dzibyic + d2edphi2 * dphidyic * dphidzib;
            *hessz[3][ib] = *hessz[3][ib] + dedphi * dzibzic + d2edphi2 * dphidzic * dphidzib;
            *hessx[1][id] = *hessx[1][id] + dedphi * dxicxid + d2edphi2 * dphidxic * dphidxid;
            *hessy[1][id] = *hessy[1][id] + dedphi * dyicxid + d2edphi2 * dphidyic * dphidxid;
            *hessz[1][id] = *hessz[1][id] + dedphi * dzicxid + d2edphi2 * dphidzic * dphidxid;
            *hessx[2][id] = *hessx[2][id] + dedphi * dxicyid + d2edphi2 * dphidxic * dphidyid;
            *hessy[2][id] = *hessy[2][id] + dedphi * dyicyid + d2edphi2 * dphidyic * dphidyid;
            *hessz[2][id] = *hessz[2][id] + dedphi * dzicyid + d2edphi2 * dphidzic * dphidyid;
            *hessx[3][id] = *hessx[3][id] + dedphi * dxiczid + d2edphi2 * dphidxic * dphidzid;
            *hessy[3][id] = *hessy[3][id] + dedphi * dyiczid + d2edphi2 * dphidyic * dphidzid;
            *hessz[3][id] = *hessz[3][id] + dedphi * dziczid + d2edphi2 * dphidzic * dphidzid;
            ////else if (i .eq. id) then;
            *hessx[1][id] = *hessx[1][id] + dedphi * dxidxid + d2edphi2 * dphidxid * dphidxid;
            *hessy[1][id] = *hessy[1][id] + dedphi * dxidyid + d2edphi2 * dphidxid * dphidyid;
            *hessz[1][id] = *hessz[1][id] + dedphi * dxidzid + d2edphi2 * dphidxid * dphidzid;
            *hessx[2][id] = *hessx[2][id] + dedphi * dxidyid + d2edphi2 * dphidxid * dphidyid;
            *hessy[2][id] = *hessy[2][id] + dedphi * dyidyid + d2edphi2 * dphidyid * dphidyid;
            *hessz[2][id] = *hessz[2][id] + dedphi * dyidzid + d2edphi2 * dphidyid * dphidzid;
            *hessx[3][id] = *hessx[3][id] + dedphi * dxidzid + d2edphi2 * dphidxid * dphidzid;
            *hessy[3][id] = *hessy[3][id] + dedphi * dyidzid + d2edphi2 * dphidyid * dphidzid;
            *hessz[3][id] = *hessz[3][id] + dedphi * dzidzid + d2edphi2 * dphidzid * dphidzid;
            *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxid + d2edphi2 * dphidxid * dphidxia;
            *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayid + d2edphi2 * dphidyid * dphidxia;
            *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazid + d2edphi2 * dphidzid * dphidxia;
            *hessx[2][ia] = *hessx[2][ia] + dedphi * dyiaxid + d2edphi2 * dphidxid * dphidyia;
            *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayid + d2edphi2 * dphidyid * dphidyia;
            *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazid + d2edphi2 * dphidzid * dphidyia;
            *hessx[3][ia] = *hessx[3][ia] + dedphi * dziaxid + d2edphi2 * dphidxid * dphidzia;
            *hessy[3][ia] = *hessy[3][ia] + dedphi * dziayid + d2edphi2 * dphidyid * dphidzia;
            *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazid + d2edphi2 * dphidzid * dphidzia;
            *hessx[1][ib] = *hessx[1][ib] + dedphi * dxibxid + d2edphi2 * dphidxid * dphidxib;
            *hessy[1][ib] = *hessy[1][ib] + dedphi * dxibyid + d2edphi2 * dphidyid * dphidxib;
            *hessz[1][ib] = *hessz[1][ib] + dedphi * dxibzid + d2edphi2 * dphidzid * dphidxib;
            *hessx[2][ib] = *hessx[2][ib] + dedphi * dyibxid + d2edphi2 * dphidxid * dphidyib;
            *hessy[2][ib] = *hessy[2][ib] + dedphi * dyibyid + d2edphi2 * dphidyid * dphidyib;
            *hessz[2][ib] = *hessz[2][ib] + dedphi * dyibzid + d2edphi2 * dphidzid * dphidyib;
            *hessx[3][ib] = *hessx[3][ib] + dedphi * dzibxid + d2edphi2 * dphidxid * dphidzib;
            *hessy[3][ib] = *hessy[3][ib] + dedphi * dzibyid + d2edphi2 * dphidyid * dphidzib;
            *hessz[3][ib] = *hessz[3][ib] + dedphi * dzibzid + d2edphi2 * dphidzid * dphidzib;
            *hessx[1][ic] = *hessx[1][ic] + dedphi * dxicxid + d2edphi2 * dphidxid * dphidxic;
            *hessy[1][ic] = *hessy[1][ic] + dedphi * dxicyid + d2edphi2 * dphidyid * dphidxic;
            *hessz[1][ic] = *hessz[1][ic] + dedphi * dxiczid + d2edphi2 * dphidzid * dphidxic;
            *hessx[2][ic] = *hessx[2][ic] + dedphi * dyicxid + d2edphi2 * dphidxid * dphidyic;
            *hessy[2][ic] = *hessy[2][ic] + dedphi * dyicyid + d2edphi2 * dphidyid * dphidyic;
            *hessz[2][ic] = *hessz[2][ic] + dedphi * dyiczid + d2edphi2 * dphidzid * dphidyic;
            *hessx[3][ic] = *hessx[3][ic] + dedphi * dzicxid + d2edphi2 * dphidxid * dphidzic;
            *hessy[3][ic] = *hessy[3][ic] + dedphi * dzicyid + d2edphi2 * dphidyid * dphidzic;
            *hessz[3][ic] = *hessz[3][ic] + dedphi * dziczid + d2edphi2 * dphidzid * dphidzic;
          }
          E += this->f_imptor_1_calc_one_torsion(imptor);
        }
        return E;
      }


      /********************************
      *                                *
      *  Improper                      *
      *  Dihedral Potential            *
      *  Energy/Gradients/Hessians     *
      *                                *
      *                                *
      *********************************/

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_imp<0>(void)
      {
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt; using std::acos;
        using std::abs;
        using std::min; using std::max;
        coords::float_type E(0.0);
        coords::Cartesian_Point vir1, vir2, vir3, vir4;
        for (auto& improper : refined.impropers())
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
          coords::float_type const rt2(dot(t, t)), ru2(dot(u, u)), rtru = sqrt(rt2 * ru2);
          if (abs(rtru) > 0.0)
          {
            coords::float_type const rcb(len(cb));
            coords::float_type const cosine(min(1.0, max(-1.0, (dot(t, u) / rtru))));
            coords::float_type const sine(dot(cb, tu) / (rcb * rtru));
            coords::float_type const angle(sine < 0.0 ?
              -acos(cosine) * SCON_180PI : acos(cosine) * SCON_180PI);
            coords::float_type da((abs(angle + improper.p.ideal) < abs(angle - improper.p.ideal)) ?
              angle + improper.p.ideal : angle - improper.p.ideal);
            if (da > 180.0) da -= 360.0;
            if (da < -180.0) da += 360.0;
            da *= SCON_PI180;
            coords::float_type dE = improper.p.force * da;
            E += dE * da;
          }
        }
        return E;
      }

      coords::float_type energy::interfaces::aco::aco_ff::f_improper_1_calc_one_torsion(::tinker::refine::types::improper const& improper)
      {
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt; using std::acos;
        using std::abs;
        using std::min; using std::max;
        coords::float_type vxx, vyx, vzx, vyy, vzy, vzz;
        coords::Cartesian_Point vir1, vir2, vir3, vir4;
        coords::float_type E = 0.;
        coords::Cartesian_Point const ba =
          coords->xyz(improper.ligand[0]) - coords->xyz(improper.center);
        coords::Cartesian_Point const cb =
          coords->xyz(improper.ligand[1]) - coords->xyz(improper.ligand[0]);
        coords::Cartesian_Point const dc =
          coords->xyz(improper.twist) - coords->xyz(improper.ligand[1]);
        coords::Cartesian_Point const t(cross(ba, cb));
        coords::Cartesian_Point const u(cross(cb, dc));
        coords::Cartesian_Point const tu(cross(t, u));
        coords::float_type const rt2(dot(t, t)), ru2(dot(u, u)), rtru = sqrt(rt2 * ru2);
        if (abs(rtru) > 0.0)
        {
          coords::float_type const rcb(len(cb));
          coords::float_type const cosine(min(1.0, max(-1.0, (dot(t, u) / rtru))));
          coords::float_type const sine(dot(cb, tu) / (rcb * rtru));
          coords::float_type const angle(sine < 0.0 ?
            -acos(cosine) * SCON_180PI : acos(cosine) * SCON_180PI);
          coords::float_type da((abs(angle + improper.p.ideal) < abs(angle - improper.p.ideal)) ?
            angle + improper.p.ideal : angle - improper.p.ideal);
          if (da > 180.0) da -= 360.0;
          if (da < -180.0) da += 360.0;
          da *= SCON_PI180;
          coords::float_type dE = improper.p.force * da;
          E += dE * da;
          dE *= 2.0;

          coords::Cartesian_Point const ca =
            coords->xyz(improper.ligand[1]) - coords->xyz(improper.center);
          coords::Cartesian_Point const db =
            coords->xyz(improper.twist) - coords->xyz(improper.ligand[0U]);
          coords::Cartesian_Point const dt(cross(t, cb) * (dE / (rt2 * rcb))),
            du(cross(u, cb) * (-dE / (ru2 * rcb)));

          vir1 = cross(dt, cb);
          vir2 = cross(ca, dt) + cross(du, dc);
          vir3 = cross(dt, ba) + cross(db, du);
          vir4 = cross(du, cb);

          part_grad[IMPROPER][improper.center] += vir1;
          part_grad[IMPROPER][improper.ligand[0]] += vir2;
          part_grad[IMPROPER][improper.ligand[1]] += vir3;
          part_grad[IMPROPER][improper.twist] += vir4;

          vxx = cb.x() * (vir3.x() + vir4.x()) - ba.x() * vir1.x() + dc.x() * vir4.x();
          vyx = cb.y() * (vir3.x() + vir4.x()) - ba.y() * vir1.x() + dc.y() * vir4.x();
          vzx = cb.z() * (vir3.x() + vir4.x()) - ba.z() * vir1.x() + dc.z() * vir4.x();
          vyy = cb.y() * (vir3.y() + vir4.y()) - ba.y() * vir1.y() + dc.y() * vir4.y();
          vzy = cb.z() * (vir3.y() + vir4.y()) - ba.z() * vir1.y() + dc.z() * vir4.y();
          vzz = cb.z() * (vir3.z() + vir4.z()) - ba.z() * vir1.z() + dc.z() * vir4.z();

          part_virial[IMPROPER][0][0] += vxx;
          part_virial[IMPROPER][1][0] += vyx;
          part_virial[IMPROPER][2][0] += vzx;
          part_virial[IMPROPER][0][1] += vyx;
          part_virial[IMPROPER][1][1] += vyy;
          part_virial[IMPROPER][2][1] += vzy;
          part_virial[IMPROPER][0][2] += vzx;
          part_virial[IMPROPER][1][2] += vzy;
          part_virial[IMPROPER][2][2] += vzz;
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_imp<1>(void)
      {
        coords::float_type E(0.0);
        for (auto& improper : refined.impropers())
        {
          E += f_improper_1_calc_one_torsion(improper);
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_imp<2>(void)
      {
        std::vector < std::vector<double*>> hessx(4u, std::vector<double*>(unsigned int(coords->atoms().size()), nullptr));
        std::vector < std::vector<double*>> hessy(hessx), hessz(hessx);
        for (std::size_t i = 1u; i <= 3u; ++i)
        {
          for (std::size_t atom = 0u; atom < coords->atoms().size(); ++atom)
          {
            std::size_t param = atom * 3 + (i - 1u);
            hessx[i][atom] = &part_hessian[TORSION][atom * 3][atom * 3];
            hessy[i][atom] = &part_hessian[TORSION][atom * 3][atom * 3 + 1u];
            hessz[i][atom] = &part_hessian[TORSION][atom * 3][atom * 3 + 2u];
          }
        }
        coords::float_type E(0.0);
        for (auto& improper : refined.impropers())
        {
          E += f_improper_1_calc_one_torsion(improper);
          const double xia = coords->xyz(improper.ligand[0]).x();
          const double yia = coords->xyz(improper.ligand[0]).y();
          const double zia = coords->xyz(improper.ligand[0]).z();
          const double xib = coords->xyz(improper.ligand[1]).x();
          const double yib = coords->xyz(improper.ligand[1]).y();
          const double zib = coords->xyz(improper.ligand[1]).z();
          const double xic = coords->xyz(improper.center).x();
          const double yic = coords->xyz(improper.center).y();
          const double zic = coords->xyz(improper.center).z();
          const double xid = coords->xyz(improper.twist).x();
          const double yid = coords->xyz(improper.twist).y();
          const double zid = coords->xyz(improper.twist).z();
          const double xba = xib - xia                  ;
          const double yba = yib - yia                  ;
          const double zba = zib - zia                  ;
          const double xcb = xic - xib                  ;
          const double ycb = yic - yib                  ;
          const double zcb = zic - zib                  ;
          const double xdc = xid - xic                  ;
          const double ydc = yid - yic                  ;
          const double zdc = zid - zic                  ;
          const double xt = yba * zcb - ycb * zba       ;
          const double yt = zba * xcb - zcb * xba       ;
          const double zt = xba * ycb - xcb * yba       ;
          const double xu = ycb * zdc - ydc * zcb       ;
          const double yu = zcb * xdc - zdc * xcb       ;
          const double zu = xcb * ydc - xdc * ycb       ;
          const double xtu = yt * zu - yu * zt          ;
          const double ytu = zt * xu - zu * xt          ;
          const double ztu = xt * yu - xu * yt          ;
          const double rt2 = xt * xt + yt * yt + zt * zt;
          const double ru2 = xu * xu + yu * yu + zu * zu;
          const double rtru = sqrt(rt2 * ru2);
          if (rtru != 0.0)
          {
            using std::max;
            using std::min;
            using std::abs;
            const double rcb = sqrt(xcb * xcb + ycb * ycb + zcb * zcb);
            double cosine = (xt * xu + yt * yu + zt * zu) / rtru;
            const double sine = (xcb * xtu + ycb * ytu + zcb * ztu) / (rcb * rtru);
            cosine = min(1.0, max(-1.0, cosine));
            const double radian = SCON_180PI;
            double angle = SCON_180PI * acos(cosine);
            if (sine < 0.0)  
              angle = -angle;
            // set the improper dihedral parameters for this angle
            double ideal = improper.p.ideal;
            const double force = improper.p.force;
            if (abs(angle + ideal) < abs(angle - ideal))
              ideal = -ideal;
            double dt = angle - ideal;
            while (dt > 180.0)
              dt = dt - 360.0;
            while (dt < - 180.0)
              dt = dt + 360.0;
            dt = dt / radian;
            //calculate the improper torsion master chain rule terms
            const double idihunit = 1.0;
            const double dedphi = 2.0 * idihunit * force * dt;
            const double d2edphi2 = 2.0 * idihunit * force;
            const double xca = xic - xia;
            const double yca = yic - yia;
            const double zca = zic - zia;
            const double xdb = xid - xib;
            const double ydb = yid - yib;
            const double zdb = zid - zib;
            const double dphidxt = (yt * zcb - ycb * zt) / (rt2 * rcb);
            const double dphidyt = (zt * xcb - zcb * xt) / (rt2 * rcb);
            const double dphidzt = (xt * ycb - xcb * yt) / (rt2 * rcb);
            const double dphidxu = -(yu * zcb - ycb * zu) / (ru2 * rcb);
            const double dphidyu = -(zu * xcb - zcb * xu) / (ru2 * rcb);
            const double dphidzu = -(xu * ycb - xcb * yu) / (ru2 * rcb);

            // abbreviations for second derivative chain rule terms
            const double xycb2 = xcb * xcb + ycb * ycb            ;
            const double xzcb2 = xcb * xcb + zcb * zcb            ;
            const double yzcb2 = ycb * ycb + zcb * zcb            ;
            const double rcbxt = -2.0 * rcb * dphidxt           ;
            const double rcbyt = -2.0 * rcb * dphidyt           ;
            const double rcbzt = -2.0 * rcb * dphidzt           ;
            const double rcbt2 = rcb * rt2                        ;
            const double rcbxu = 2.0 * rcb * dphidxu            ;
            const double rcbyu = 2.0 * rcb * dphidyu            ;
            const double rcbzu = 2.0 * rcb * dphidzu            ;
            const double rcbu2 = rcb * ru2                        ;
            const double dphidxibt = yca * dphidzt - zca * dphidyt;
            const double dphidxibu = zdc * dphidyu - ydc * dphidzu;
            const double dphidyibt = zca * dphidxt - xca * dphidzt;
            const double dphidyibu = xdc * dphidzu - zdc * dphidxu;
            const double dphidzibt = xca * dphidyt - yca * dphidxt;
            const double dphidzibu = ydc * dphidxu - xdc * dphidyu;
            const double dphidxict = zba * dphidyt - yba * dphidzt;
            const double dphidxicu = ydb * dphidzu - zdb * dphidyu;
            const double dphidyict = xba * dphidzt - zba * dphidxt;
            const double dphidyicu = zdb * dphidxu - xdb * dphidzu;
            const double dphidzict = yba * dphidxt - xba * dphidyt;
            const double dphidzicu = xdb * dphidyu - ydb * dphidxu;

            // chain rule terms for first derivative components
            const double dphidxia = zcb * dphidyt - ycb * dphidzt;
            const double dphidyia = xcb * dphidzt - zcb * dphidxt;
            const double dphidzia = ycb * dphidxt - xcb * dphidyt;
            const double dphidxib = dphidxibt + dphidxibu        ;
            const double dphidyib = dphidyibt + dphidyibu        ;
            const double dphidzib = dphidzibt + dphidzibu        ;
            const double dphidxic = dphidxict + dphidxicu        ;
            const double dphidyic = dphidyict + dphidyicu        ;
            const double dphidzic = dphidzict + dphidzicu        ;
            const double dphidxid = zcb * dphidyu - ycb * dphidzu;
            const double dphidyid = xcb * dphidzu - zcb * dphidxu;
            const double dphidzid = ycb * dphidxu - xcb * dphidyu;

            // chain rule terms for second derivative components
            const double dxiaxia = rcbxt * dphidxia;
            const double dxiayia = rcbxt * dphidyia - zcb * rcb / rt2;
            const double dxiazia = rcbxt * dphidzia + ycb * rcb / rt2;
            const double dxiaxic = rcbxt * dphidxict + xcb * xt / rcbt2;
            const double dxiayic = rcbxt * dphidyict - dphidzt - (xba * zcb * xcb + zba * yzcb2) / rcbt2;
            const double dxiazic = rcbxt * dphidzict + dphidyt + (xba * ycb * xcb + yba * yzcb2) / rcbt2;
            const double dxiaxid = 0.0;
            const double dxiayid = 0.0;
            const double dxiazid = 0.0;
            const double dyiayia = rcbyt * dphidyia;
            const double dyiazia = rcbyt * dphidzia - xcb * rcb / rt2;
            const double dyiaxib = rcbyt * dphidxibt - dphidzt - (yca * zcb * ycb + zca * xzcb2) / rcbt2;
            const double dyiaxic = rcbyt * dphidxict + dphidzt + (yba * zcb * ycb + zba * xzcb2) / rcbt2;
            const double dyiayic = rcbyt * dphidyict + ycb * yt / rcbt2;
            const double dyiazic = rcbyt * dphidzict - dphidxt - (yba * xcb * ycb + xba * xzcb2) / rcbt2;
            const double dyiaxid = 0.0;
            const double dyiayid = 0.0;
            const double dyiazid = 0.0;
            const double dziazia = rcbzt * dphidzia;
            const double dziaxib = rcbzt * dphidxibt + dphidyt + (zca * ycb * zcb + yca * xycb2) / rcbt2;
            const double dziayib = rcbzt * dphidyibt - dphidxt - (zca * xcb * zcb + xca * xycb2) / rcbt2;
            const double dziaxic = rcbzt * dphidxict - dphidyt - (zba * ycb * zcb + yba * xycb2) / rcbt2;
            const double dziayic = rcbzt * dphidyict + dphidxt + (zba * xcb * zcb + xba * xycb2) / rcbt2;
            const double dziazic = rcbzt * dphidzict + zcb * zt / rcbt2;
            const double dziaxid = 0.0;
            const double dziayid = 0.0;
            const double dziazid = 0.0;
            const double dxibxic = -xcb * dphidxib / (rcb * rcb) - (yca * (zba * xcb + yt) - zca * (yba * xcb - zt)) / rcbt2 - 2.0 * (yt * zba - yba * zt) * dphidxibt / rt2 - (zdc * (ydb * xcb + zu) - ydc * (zdb * xcb - yu)) / rcbu2 + 2.0 * (yu * zdb - ydb * zu) * dphidxibu / ru2;
            const double dxibyic = -ycb * dphidxib / (rcb * rcb) + dphidzt + dphidzu - (yca * (zba * ycb - xt) + zca * (xba * xcb + zcb * zba)) / rcbt2 - 2.0 * (zt * xba - zba * xt) * dphidxibt / rt2 + (zdc * (xdb * xcb + zcb * zdb) + ydc * (zdb * ycb + xu)) / rcbu2 + 2.0 * (zu * xdb - zdb * xu) * dphidxibu / ru2;
            const double dxibxid = rcbxu * dphidxibu + xcb * xu / rcbu2;
            const double dxibyid = rcbyu * dphidxibu - dphidzu - (ydc * zcb * ycb + zdc * xzcb2) / rcbu2;
            const double dxibzid = rcbzu * dphidxibu + dphidyu + (zdc * ycb * zcb + ydc * xycb2) / rcbu2;
            const double dyibzib = ycb * dphidzib / (rcb * rcb) - (xca * (xca * xcb + zcb * zca) + yca * (ycb * xca + zt)) / rcbt2 - 2.0 * (xt * zca - xca * zt) * dphidzibt / rt2 + (ydc * (xdc * ycb - zu) + xdc * (xdc * xcb + zcb * zdc)) / rcbu2 + 2.0 * (xu * zdc - xdc * zu) * dphidzibu / ru2;
            const double dyibxic = -xcb * dphidyib / (rcb * rcb) - dphidzt - dphidzu + (xca * (zba * xcb + yt) + zca * (zba * zcb + ycb * yba)) / rcbt2 - 2.0 * (yt * zba - yba * zt) * dphidyibt / rt2 - (zdc * (zdb * zcb + ycb * ydb) + xdc * (zdb * xcb - yu)) / rcbu2 + 2.0 * (yu * zdb - ydb * zu) * dphidyibu / ru2;
            const double dyibyic = -ycb * dphidyib / (rcb * rcb) - (zca * (xba * ycb + zt) - xca * (zba * ycb - xt)) / rcbt2 - 2.0 * (zt * xba - zba * xt) * dphidyibt / rt2 - (xdc * (zdb * ycb + xu) - zdc * (xdb * ycb - zu)) / rcbu2 + 2.0 * (zu * xdb - zdb * xu) * dphidyibu / ru2;
            const double dyibxid = rcbxu * dphidyibu + dphidzu + (xdc * zcb * xcb + zdc * yzcb2) / rcbu2;
            const double dyibyid = rcbyu * dphidyibu + ycb * yu / rcbu2;
            const double dyibzid = rcbzu * dphidyibu - dphidxu - (zdc * xcb * zcb + xdc * xycb2) / rcbu2;
            const double dzibxic = -xcb * dphidzib / (rcb * rcb) + dphidyt + dphidyu - (xca * (yba * xcb - zt) + yca * (zba * zcb + ycb * yba)) / rcbt2 - 2.0 * (yt * zba - yba * zt) * dphidzibt / rt2 + (ydc * (zdb * zcb + ycb * ydb) + xdc * (ydb * xcb + zu)) / rcbu2 + 2.0 * (yu * zdb - ydb * zu) * dphidzibu / ru2;
            const double dzibzic = -zcb * dphidzib / (rcb * rcb) - (xca * (yba * zcb + xt) - yca * (xba * zcb - yt)) / rcbt2 - 2.0 * (xt * yba - xba * yt) * dphidzibt / rt2 - (ydc * (xdb * zcb + yu) - xdc * (ydb * zcb - xu)) / rcbu2 + 2.0 * (xu * ydb - xdb * yu) * dphidzibu / ru2;
            const double dzibxid = rcbxu * dphidzibu - dphidyu - (xdc * ycb * xcb + ydc * yzcb2) / rcbu2;
            const double dzibyid = rcbyu * dphidzibu + dphidxu + (ydc * xcb * ycb + xdc * xzcb2) / rcbu2;
            const double dzibzid = rcbzu * dphidzibu + zcb * zu / rcbu2;
            const double dxicxid = rcbxu * dphidxicu - xcb * (zdb * ycb - ydb * zcb) / rcbu2;
            const double dxicyid = rcbyu * dphidxicu + dphidzu + (ydb * zcb * ycb + zdb * xzcb2) / rcbu2;
            const double dxiczid = rcbzu * dphidxicu - dphidyu - (zdb * ycb * zcb + ydb * xycb2) / rcbu2;
            const double dyicxid = rcbxu * dphidyicu - dphidzu - (xdb * zcb * xcb + zdb * yzcb2) / rcbu2;
            const double dyicyid = rcbyu * dphidyicu - ycb * (xdb * zcb - zdb * xcb) / rcbu2;
            const double dyiczid = rcbzu * dphidyicu + dphidxu + (zdb * xcb * zcb + xdb * xycb2) / rcbu2;
            const double dzicxid = rcbxu * dphidzicu + dphidyu + (xdb * ycb * xcb + ydb * yzcb2) / rcbu2;
            const double dzicyid = rcbyu * dphidzicu - dphidxu - (ydb * xcb * ycb + xdb * xzcb2) / rcbu2;
            const double dziczid = rcbzu * dphidzicu - zcb * (ydb * xcb - xdb * ycb) / rcbu2;
            const double dxidxid = rcbxu * dphidxid;
            const double dxidyid = rcbxu * dphidyid + zcb * rcb / ru2;
            const double dxidzid = rcbxu * dphidzid - ycb * rcb / ru2;
            const double dyidyid = rcbyu * dphidyid;
            const double dyidzid = rcbyu * dphidzid + xcb * rcb / ru2;
            const double dzidzid = rcbzu * dphidzid;

            // get some second derivative chain rule terms by difference
            const double dxiaxib = -dxiaxia - dxiaxic - dxiaxid;
            const double dxiayib = -dxiayia - dxiayic - dxiayid;
            const double dxiazib = -dxiazia - dxiazic - dxiazid;
            const double dyiayib = -dyiayia - dyiayic - dyiayid;
            const double dyiazib = -dyiazia - dyiazic - dyiazid;
            const double dziazib = -dziazia - dziazic - dziazid;
            const double dxibxib = -dxiaxib - dxibxic - dxibxid;
            const double dxibyib = -dyiaxib - dxibyic - dxibyid;
            const double dxibzib = -dxiazib - dzibxic - dzibxid;
            const double dxibzic = -dziaxib - dxibzib - dxibzid;
            const double dyibyib = -dyiayib - dyibyic - dyibyid;
            const double dyibzic = -dziayib - dyibzib - dyibzid;
            const double dzibzib = -dziazib - dzibzic - dzibzid;
            const double dzibyic = -dyiazib - dyibzib - dzibyid;
            const double dxicxic = -dxiaxic - dxibxic - dxicxid;
            const double dxicyic = -dyiaxic - dyibxic - dxicyid;
            const double dxiczic = -dziaxic - dzibxic - dxiczid;
            const double dyicyic = -dyiayic - dyibyic - dyicyid;
            const double dyiczic = -dziayic - dzibyic - dyiczid;
            const double dziczic = -dziazic - dzibzic - dziczid;

            // increment diagonal and off-diagonal Hessian elements
            const std::size_t ia = improper.center;
            const std::size_t ib = improper.ligand[0];
            const std::size_t ic = improper.ligand[1];
            const std::size_t id = improper.ligand[2];
            //if[i.eq.ia] then;
            *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxia + d2edphi2 * dphidxia * dphidxia;
            *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayia + d2edphi2 * dphidxia * dphidyia;
            *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazia + d2edphi2 * dphidxia * dphidzia;
            *hessx[2][ia] = *hessx[2][ia] + dedphi * dxiayia + d2edphi2 * dphidxia * dphidyia;
            *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayia + d2edphi2 * dphidyia * dphidyia;
            *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazia + d2edphi2 * dphidyia * dphidzia;
            *hessx[3][ia] = *hessx[3][ia] + dedphi * dxiazia + d2edphi2 * dphidxia * dphidzia;
            *hessy[3][ia] = *hessy[3][ia] + dedphi * dyiazia + d2edphi2 * dphidyia * dphidzia;
            *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazia + d2edphi2 * dphidzia * dphidzia;
            *hessx[1][ib] = *hessx[1][ib] + dedphi * dxiaxib + d2edphi2 * dphidxia * dphidxib;
            *hessy[1][ib] = *hessy[1][ib] + dedphi * dyiaxib + d2edphi2 * dphidyia * dphidxib;
            *hessz[1][ib] = *hessz[1][ib] + dedphi * dziaxib + d2edphi2 * dphidzia * dphidxib;
            *hessx[2][ib] = *hessx[2][ib] + dedphi * dxiayib + d2edphi2 * dphidxia * dphidyib;
            *hessy[2][ib] = *hessy[2][ib] + dedphi * dyiayib + d2edphi2 * dphidyia * dphidyib;
            *hessz[2][ib] = *hessz[2][ib] + dedphi * dziayib + d2edphi2 * dphidzia * dphidyib;
            *hessx[3][ib] = *hessx[3][ib] + dedphi * dxiazib + d2edphi2 * dphidxia * dphidzib;
            *hessy[3][ib] = *hessy[3][ib] + dedphi * dyiazib + d2edphi2 * dphidyia * dphidzib;
            *hessz[3][ib] = *hessz[3][ib] + dedphi * dziazib + d2edphi2 * dphidzia * dphidzib;
            *hessx[1][ic] = *hessx[1][ic] + dedphi * dxiaxic + d2edphi2 * dphidxia * dphidxic;
            *hessy[1][ic] = *hessy[1][ic] + dedphi * dyiaxic + d2edphi2 * dphidyia * dphidxic;
            *hessz[1][ic] = *hessz[1][ic] + dedphi * dziaxic + d2edphi2 * dphidzia * dphidxic;
            *hessx[2][ic] = *hessx[2][ic] + dedphi * dxiayic + d2edphi2 * dphidxia * dphidyic;
            *hessy[2][ic] = *hessy[2][ic] + dedphi * dyiayic + d2edphi2 * dphidyia * dphidyic;
            *hessz[2][ic] = *hessz[2][ic] + dedphi * dziayic + d2edphi2 * dphidzia * dphidyic;
            *hessx[3][ic] = *hessx[3][ic] + dedphi * dxiazic + d2edphi2 * dphidxia * dphidzic;
            *hessy[3][ic] = *hessy[3][ic] + dedphi * dyiazic + d2edphi2 * dphidyia * dphidzic;
            *hessz[3][ic] = *hessz[3][ic] + dedphi * dziazic + d2edphi2 * dphidzia * dphidzic;
            *hessx[1][id] = *hessx[1][id] + dedphi * dxiaxid + d2edphi2 * dphidxia * dphidxid;
            *hessy[1][id] = *hessy[1][id] + dedphi * dyiaxid + d2edphi2 * dphidyia * dphidxid;
            *hessz[1][id] = *hessz[1][id] + dedphi * dziaxid + d2edphi2 * dphidzia * dphidxid;
            *hessx[2][id] = *hessx[2][id] + dedphi * dxiayid + d2edphi2 * dphidxia * dphidyid;
            *hessy[2][id] = *hessy[2][id] + dedphi * dyiayid + d2edphi2 * dphidyia * dphidyid;
            *hessz[2][id] = *hessz[2][id] + dedphi * dziayid + d2edphi2 * dphidzia * dphidyid;
            *hessx[3][id] = *hessx[3][id] + dedphi * dxiazid + d2edphi2 * dphidxia * dphidzid;
            *hessy[3][id] = *hessy[3][id] + dedphi * dyiazid + d2edphi2 * dphidyia * dphidzid;
            *hessz[3][id] = *hessz[3][id] + dedphi * dziazid + d2edphi2 * dphidzia * dphidzid;
            //else if[i.eq.ib] then;
            *hessx[1][ib] = *hessx[1][ib] + dedphi * dxibxib + d2edphi2 * dphidxib * dphidxib;
            *hessy[1][ib] = *hessy[1][ib] + dedphi * dxibyib + d2edphi2 * dphidxib * dphidyib;
            *hessz[1][ib] = *hessz[1][ib] + dedphi * dxibzib + d2edphi2 * dphidxib * dphidzib;
            *hessx[2][ib] = *hessx[2][ib] + dedphi * dxibyib + d2edphi2 * dphidxib * dphidyib;
            *hessy[2][ib] = *hessy[2][ib] + dedphi * dyibyib + d2edphi2 * dphidyib * dphidyib;
            *hessz[2][ib] = *hessz[2][ib] + dedphi * dyibzib + d2edphi2 * dphidyib * dphidzib;
            *hessx[3][ib] = *hessx[3][ib] + dedphi * dxibzib + d2edphi2 * dphidxib * dphidzib;
            *hessy[3][ib] = *hessy[3][ib] + dedphi * dyibzib + d2edphi2 * dphidyib * dphidzib;
            *hessz[3][ib] = *hessz[3][ib] + dedphi * dzibzib + d2edphi2 * dphidzib * dphidzib;
            *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxib + d2edphi2 * dphidxib * dphidxia;
            *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayib + d2edphi2 * dphidyib * dphidxia;
            *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazib + d2edphi2 * dphidzib * dphidxia;
            *hessx[2][ia] = *hessx[2][ia] + dedphi * dyiaxib + d2edphi2 * dphidxib * dphidyia;
            *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayib + d2edphi2 * dphidyib * dphidyia;
            *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazib + d2edphi2 * dphidzib * dphidyia;
            *hessx[3][ia] = *hessx[3][ia] + dedphi * dziaxib + d2edphi2 * dphidxib * dphidzia;
            *hessy[3][ia] = *hessy[3][ia] + dedphi * dziayib + d2edphi2 * dphidyib * dphidzia;
            *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazib + d2edphi2 * dphidzib * dphidzia;
            *hessx[1][ic] = *hessx[1][ic] + dedphi * dxibxic + d2edphi2 * dphidxib * dphidxic;
            *hessy[1][ic] = *hessy[1][ic] + dedphi * dyibxic + d2edphi2 * dphidyib * dphidxic;
            *hessz[1][ic] = *hessz[1][ic] + dedphi * dzibxic + d2edphi2 * dphidzib * dphidxic;
            *hessx[2][ic] = *hessx[2][ic] + dedphi * dxibyic + d2edphi2 * dphidxib * dphidyic;
            *hessy[2][ic] = *hessy[2][ic] + dedphi * dyibyic + d2edphi2 * dphidyib * dphidyic;
            *hessz[2][ic] = *hessz[2][ic] + dedphi * dzibyic + d2edphi2 * dphidzib * dphidyic;
            *hessx[3][ic] = *hessx[3][ic] + dedphi * dxibzic + d2edphi2 * dphidxib * dphidzic;
            *hessy[3][ic] = *hessy[3][ic] + dedphi * dyibzic + d2edphi2 * dphidyib * dphidzic;
            *hessz[3][ic] = *hessz[3][ic] + dedphi * dzibzic + d2edphi2 * dphidzib * dphidzic;
            *hessx[1][id] = *hessx[1][id] + dedphi * dxibxid + d2edphi2 * dphidxib * dphidxid;
            *hessy[1][id] = *hessy[1][id] + dedphi * dyibxid + d2edphi2 * dphidyib * dphidxid;
            *hessz[1][id] = *hessz[1][id] + dedphi * dzibxid + d2edphi2 * dphidzib * dphidxid;
            *hessx[2][id] = *hessx[2][id] + dedphi * dxibyid + d2edphi2 * dphidxib * dphidyid;
            *hessy[2][id] = *hessy[2][id] + dedphi * dyibyid + d2edphi2 * dphidyib * dphidyid;
            *hessz[2][id] = *hessz[2][id] + dedphi * dzibyid + d2edphi2 * dphidzib * dphidyid;
            *hessx[3][id] = *hessx[3][id] + dedphi * dxibzid + d2edphi2 * dphidxib * dphidzid;
            *hessy[3][id] = *hessy[3][id] + dedphi * dyibzid + d2edphi2 * dphidyib * dphidzid;
            *hessz[3][id] = *hessz[3][id] + dedphi * dzibzid + d2edphi2 * dphidzib * dphidzid;
            //               else if [i .eq. ic] then;
            *hessx[1][ic] = *hessx[1][ic] + dedphi * dxicxic + d2edphi2 * dphidxic * dphidxic;
            *hessy[1][ic] = *hessy[1][ic] + dedphi * dxicyic + d2edphi2 * dphidxic * dphidyic;
            *hessz[1][ic] = *hessz[1][ic] + dedphi * dxiczic + d2edphi2 * dphidxic * dphidzic;
            *hessx[2][ic] = *hessx[2][ic] + dedphi * dxicyic + d2edphi2 * dphidxic * dphidyic;
            *hessy[2][ic] = *hessy[2][ic] + dedphi * dyicyic + d2edphi2 * dphidyic * dphidyic;
            *hessz[2][ic] = *hessz[2][ic] + dedphi * dyiczic + d2edphi2 * dphidyic * dphidzic;
            *hessx[3][ic] = *hessx[3][ic] + dedphi * dxiczic + d2edphi2 * dphidxic * dphidzic;
            *hessy[3][ic] = *hessy[3][ic] + dedphi * dyiczic + d2edphi2 * dphidyic * dphidzic;
            *hessz[3][ic] = *hessz[3][ic] + dedphi * dziczic + d2edphi2 * dphidzic * dphidzic;
            *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxic + d2edphi2 * dphidxic * dphidxia;
            *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayic + d2edphi2 * dphidyic * dphidxia;
            *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazic + d2edphi2 * dphidzic * dphidxia;
            *hessx[2][ia] = *hessx[2][ia] + dedphi * dyiaxic + d2edphi2 * dphidxic * dphidyia;
            *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayic + d2edphi2 * dphidyic * dphidyia;
            *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazic + d2edphi2 * dphidzic * dphidyia;
            *hessx[3][ia] = *hessx[3][ia] + dedphi * dziaxic + d2edphi2 * dphidxic * dphidzia;
            *hessy[3][ia] = *hessy[3][ia] + dedphi * dziayic + d2edphi2 * dphidyic * dphidzia;
            *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazic + d2edphi2 * dphidzic * dphidzia;
            *hessx[1][ib] = *hessx[1][ib] + dedphi * dxibxic + d2edphi2 * dphidxic * dphidxib;
            *hessy[1][ib] = *hessy[1][ib] + dedphi * dxibyic + d2edphi2 * dphidyic * dphidxib;
            *hessz[1][ib] = *hessz[1][ib] + dedphi * dxibzic + d2edphi2 * dphidzic * dphidxib;
            *hessx[2][ib] = *hessx[2][ib] + dedphi * dyibxic + d2edphi2 * dphidxic * dphidyib;
            *hessy[2][ib] = *hessy[2][ib] + dedphi * dyibyic + d2edphi2 * dphidyic * dphidyib;
            *hessz[2][ib] = *hessz[2][ib] + dedphi * dyibzic + d2edphi2 * dphidzic * dphidyib;
            *hessx[3][ib] = *hessx[3][ib] + dedphi * dzibxic + d2edphi2 * dphidxic * dphidzib;
            *hessy[3][ib] = *hessy[3][ib] + dedphi * dzibyic + d2edphi2 * dphidyic * dphidzib;
            *hessz[3][ib] = *hessz[3][ib] + dedphi * dzibzic + d2edphi2 * dphidzic * dphidzib;
            *hessx[1][id] = *hessx[1][id] + dedphi * dxicxid + d2edphi2 * dphidxic * dphidxid;
            *hessy[1][id] = *hessy[1][id] + dedphi * dyicxid + d2edphi2 * dphidyic * dphidxid;
            *hessz[1][id] = *hessz[1][id] + dedphi * dzicxid + d2edphi2 * dphidzic * dphidxid;
            *hessx[2][id] = *hessx[2][id] + dedphi * dxicyid + d2edphi2 * dphidxic * dphidyid;
            *hessy[2][id] = *hessy[2][id] + dedphi * dyicyid + d2edphi2 * dphidyic * dphidyid;
            *hessz[2][id] = *hessz[2][id] + dedphi * dzicyid + d2edphi2 * dphidzic * dphidyid;
            *hessx[3][id] = *hessx[3][id] + dedphi * dxiczid + d2edphi2 * dphidxic * dphidzid;
            *hessy[3][id] = *hessy[3][id] + dedphi * dyiczid + d2edphi2 * dphidyic * dphidzid;
            *hessz[3][id] = *hessz[3][id] + dedphi * dziczid + d2edphi2 * dphidzic * dphidzid;
            //else if[i.eq.id] then;
            *hessx[1][id] = *hessx[1][id] + dedphi * dxidxid + d2edphi2 * dphidxid * dphidxid;
            *hessy[1][id] = *hessy[1][id] + dedphi * dxidyid + d2edphi2 * dphidxid * dphidyid;
            *hessz[1][id] = *hessz[1][id] + dedphi * dxidzid + d2edphi2 * dphidxid * dphidzid;
            *hessx[2][id] = *hessx[2][id] + dedphi * dxidyid + d2edphi2 * dphidxid * dphidyid;
            *hessy[2][id] = *hessy[2][id] + dedphi * dyidyid + d2edphi2 * dphidyid * dphidyid;
            *hessz[2][id] = *hessz[2][id] + dedphi * dyidzid + d2edphi2 * dphidyid * dphidzid;
            *hessx[3][id] = *hessx[3][id] + dedphi * dxidzid + d2edphi2 * dphidxid * dphidzid;
            *hessy[3][id] = *hessy[3][id] + dedphi * dyidzid + d2edphi2 * dphidyid * dphidzid;
            *hessz[3][id] = *hessz[3][id] + dedphi * dzidzid + d2edphi2 * dphidzid * dphidzid;
            *hessx[1][ia] = *hessx[1][ia] + dedphi * dxiaxid + d2edphi2 * dphidxid * dphidxia;
            *hessy[1][ia] = *hessy[1][ia] + dedphi * dxiayid + d2edphi2 * dphidyid * dphidxia;
            *hessz[1][ia] = *hessz[1][ia] + dedphi * dxiazid + d2edphi2 * dphidzid * dphidxia;
            *hessx[2][ia] = *hessx[2][ia] + dedphi * dyiaxid + d2edphi2 * dphidxid * dphidyia;
            *hessy[2][ia] = *hessy[2][ia] + dedphi * dyiayid + d2edphi2 * dphidyid * dphidyia;
            *hessz[2][ia] = *hessz[2][ia] + dedphi * dyiazid + d2edphi2 * dphidzid * dphidyia;
            *hessx[3][ia] = *hessx[3][ia] + dedphi * dziaxid + d2edphi2 * dphidxid * dphidzia;
            *hessy[3][ia] = *hessy[3][ia] + dedphi * dziayid + d2edphi2 * dphidyid * dphidzia;
            *hessz[3][ia] = *hessz[3][ia] + dedphi * dziazid + d2edphi2 * dphidzid * dphidzia;
            *hessx[1][ib] = *hessx[1][ib] + dedphi * dxibxid + d2edphi2 * dphidxid * dphidxib;
            *hessy[1][ib] = *hessy[1][ib] + dedphi * dxibyid + d2edphi2 * dphidyid * dphidxib;
            *hessz[1][ib] = *hessz[1][ib] + dedphi * dxibzid + d2edphi2 * dphidzid * dphidxib;
            *hessx[2][ib] = *hessx[2][ib] + dedphi * dyibxid + d2edphi2 * dphidxid * dphidyib;
            *hessy[2][ib] = *hessy[2][ib] + dedphi * dyibyid + d2edphi2 * dphidyid * dphidyib;
            *hessz[2][ib] = *hessz[2][ib] + dedphi * dyibzid + d2edphi2 * dphidzid * dphidyib;
            *hessx[3][ib] = *hessx[3][ib] + dedphi * dzibxid + d2edphi2 * dphidxid * dphidzib;
            *hessy[3][ib] = *hessy[3][ib] + dedphi * dzibyid + d2edphi2 * dphidyid * dphidzib;
            *hessz[3][ib] = *hessz[3][ib] + dedphi * dzibzid + d2edphi2 * dphidzid * dphidzib;
            *hessx[1][ic] = *hessx[1][ic] + dedphi * dxicxid + d2edphi2 * dphidxid * dphidxic;
            *hessy[1][ic] = *hessy[1][ic] + dedphi * dxicyid + d2edphi2 * dphidyid * dphidxic;
            *hessz[1][ic] = *hessz[1][ic] + dedphi * dxiczid + d2edphi2 * dphidzid * dphidxic;
            *hessx[2][ic] = *hessx[2][ic] + dedphi * dyicxid + d2edphi2 * dphidxid * dphidyic;
            *hessy[2][ic] = *hessy[2][ic] + dedphi * dyicyid + d2edphi2 * dphidyid * dphidyic;
            *hessz[2][ic] = *hessz[2][ic] + dedphi * dyiczid + d2edphi2 * dphidzid * dphidyic;
            *hessx[3][ic] = *hessx[3][ic] + dedphi * dzicxid + d2edphi2 * dphidxid * dphidzic;
            *hessy[3][ic] = *hessy[3][ic] + dedphi * dzicyid + d2edphi2 * dphidyid * dphidzic;
            *hessz[3][ic] = *hessz[3][ic] + dedphi * dziczid + d2edphi2 * dphidzid * dphidzic;
          }

        }
        return E;
      }

      ////////////////////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////////////////////
      
      ////////////////////////////////////////////////////////////////////////////
      
      energy::interfaces::aco::nb_cutoff::nb_cutoff(
        coords::float_type const cutoffDistance, coords::float_type const switchDistance)
        : c(cutoffDistance), s(switchDistance), cc(c* c), ss(3.0 * s * s),
        cs((cc - s * s)* (cc - s * s)* (cc - s * s))
      { }

      /**the scaling factors are calculated according to https://doi.org/10.1002/jcc.540150702*/
      bool energy::interfaces::aco::nb_cutoff::factors(coords::float_type const rr,
        coords::float_type& r, coords::float_type& fQ, coords::float_type& fV)
      {
        using std::sqrt;
        using std::abs;
        r = sqrt(rr);
        if (r > c) return false;         // if distance bigger than cutoff -> return false
        coords::float_type const cr(cc - rr);
        fV = r < s ? 1.0 : (cr * cr * (cc + 2.0 * rr - ss)) / cs; // equation 4 in paper
        fQ = (1.0 - rr / cc) * (1.0 - rr / cc);              // equation 5 in paper
        return (abs(r) > 0.0);          // return true (distance smaller than cutoff)
      }

      void aco::aco_ff::calc_ext_charges_interaction(size_t deriv)
      {
        auto elec_factor = 332.0;  // factor for conversion of charge product into amber units
        grad_ext_charges.clear();  // reset vector for gradients of external charges
        coords::Cartesian_Point ext_grad;

        for (auto& c : Config::get().energy.qmmm.mm_charges) // loop over all external charges
        {
          ext_grad.x() = 0.0;  // set ext_grad for current charge to zero
          ext_grad.y() = 0.0;
          ext_grad.z() = 0.0;

          for (auto i = 0u; i < coords->size(); ++i)               // loop over all atoms
          {
            double atom_charge = charges()[i];
            double charge_product = c.scaled_charge * atom_charge * elec_factor;

            double dist_x = coords->xyz(i).x() - c.x;
            double dist_y = coords->xyz(i).y() - c.y;
            double dist_z = coords->xyz(i).z() - c.z;
            coords::Cartesian_Point vector{ dist_x,dist_y,dist_z }; // connection vector between charge and atom

            double dist = std::sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);  // distance or length of vector
            double inverse_dist = 1.0 / dist;  // get inverse distance

            if (deriv == 0u) part_energy[CHARGE] += eQ(charge_product, inverse_dist);  // energy calculation

            else if (deriv == 1u) // gradient calculation
            {
              coords::float_type dQ;
              part_energy[EXTERNAL_CHARGES] += gQ(charge_product, inverse_dist, dQ);

              coords::Cartesian_Point grad = (vector / dist) * dQ;     // dQ is a float, now the gradient gets a direction

              part_grad[EXTERNAL_CHARGES][i] += grad;  // gradient on atom
              ext_grad -= grad;              // gradient on external charge
            }
            else // if deriv > 1
            {
              throw std::runtime_error("External charge derivatives not implemented for derivatives > 1.\n");
            }
          }
          grad_ext_charges.push_back(ext_grad);  // add gradient on external charge to vector
        }
      }

      /**calculate coulomb potential;
      returns the energy
      @param C: product of the charges (in amber-units)
      @param ri: inverse distance between the two atoms */
      coords::float_type energy::interfaces::aco::aco_ff::eQ
      (coords::float_type const C, coords::float_type const ri)
      {
        return C * ri; //kcal/mol
      }

      /**calculate coulomb potential and gradient for FEP;
      returns the energy
      @param C: product of the charges (in amber-units)
      @param ri: inverse distance between the two atoms
      @param dQ: reference to variable that saves absolute value of gradient */
      coords::float_type energy::interfaces::aco::aco_ff::gQ
      (coords::float_type const C, coords::float_type const ri, coords::float_type& dQ)
      {
        coords::float_type const Q = C * ri; // Q = C/r [kcal/mol]
        dQ = -Q * ri; //[kcal/(mol*Angstrom)]
        return Q;
      }

      /**calculate coulomb potential and gradient for FEP;
      returns the energy
      @param C: product of the charges (in amber-units)
      @param ri: distance between the two atoms
      @param cout: lambda_el
      @param dQ: reference to variable that saves absolute value of gradient */
      coords::float_type energy::interfaces::aco::aco_ff::gQ_fep
      (coords::float_type const C, coords::float_type const ri,
        coords::float_type const c_out, coords::float_type& dQ)
      {
        using std::pow;
        coords::float_type const rmod = (1.0 - c_out) *
          Config::get().fep.cshift + ri * ri * ri * ri * ri * ri; //r^6 shifted
        coords::float_type const Q = c_out * C /
          pow(rmod, 0.16666666666666); //Q
        dQ = -c_out * C * pow(ri, 5.0) / pow(-Config::get().fep.cshift *
          c_out + Config::get().fep.cshift + std::pow(ri, 6.0), 1.16666666666666);  // derivative
        return Q;
      }


      /*
        lennard-jones potentials
      */

      /**calculate lenard-jones potential and gradient for charmm and amber forcefield (r_min-type);
      returns the energy
      @param E: epsilon-parameter
      @param R: r_min-parameter
      @param r: inverse distance 1/r between the two atoms*/
      template<> coords::float_type energy::interfaces::aco::aco_ff::eV
        < ::tinker::parameter::radius_types::R_MIN>
        (coords::float_type const E, coords::float_type const R, coords::float_type const r)
      {
        coords::float_type T = R * r;
        T = T * T * T; // T^3
        T = T * T; // T^6
        return E * T * (T - 2.0);
      }

      /**calculate lenard-jones potential and gradient for oplsaa forcefield (sigma-type);
      returns the energy
      @param E: 4 * epsilon-parameter
      @param R: sigma-parameter
      @param r: inverse distance 1/r between the two atoms*/
      template<> coords::float_type energy::interfaces::aco::aco_ff::eV
        < ::tinker::parameter::radius_types::SIGMA>
        (coords::float_type const E, coords::float_type const R, coords::float_type const r)
      {
        coords::float_type T = R * r;
        T = T * T * T; // T^3
        T = T * T; // T^6
        return E * T * (T - 1.0);
      }

      /**calculate lenard-jones potential and gradient for charmm and amber forcefield (r_min-type);
      returns the energy
      @param E: 4 * epsilon-parameter
      @param R: r_min-parameter
      @param r: inverse distance 1/r between the two atoms
      @param dV: reference to variable that saves absolute value of gradient*/
      template<> coords::float_type energy::interfaces::aco::aco_ff::gV
        < ::tinker::parameter::radius_types::R_MIN>
        (coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type& dV)
      {
        coords::float_type T = R * r;
        T = T * T * T; // T^3
        T = T * T; // T^6
        coords::float_type const V = E * T;
        dV = 12.0 * V * r * (1.0 - T);
        return V * (T - 2.0);
      }

      /**calculate lenard-jones potential and gradient for oplsaa-forcefield (sigma-type);
      returns the energy
      @param E: epsilon-parameter
      @param R: sigma-parameter
      @param r: inverse distance 1/r between the two atoms
      @param dV: reference to variable that saves absolute value of gradient*/
      template<> coords::float_type energy::interfaces::aco::aco_ff::gV
        < ::tinker::parameter::radius_types::SIGMA>
        (coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type& dV)
      {
        coords::float_type T = R * r;
        T = T * T * T; // T^3
        T = T * T; // T^6
        coords::float_type const V = E * T;    //[kcal/mol]
        dV = V * r * (6.0 - 12.0 * T);  // derivative  [kcal/(mol*Angstrom)]
        return V * (T - 1.0);  //potential
      }

      /**calculate lenard-jones potential and gradient for charmm and amber forcefield (r_min-type);
      returns the energy
      @param E: 4 * epsilon-parameter
      @param R: r_min-parameter
      @param r: distance between the two atoms
      @param vout: lambda_vdw
      @param dV: reference to variable that saves absolute value of gradient*/
      template<> coords::float_type energy::interfaces::aco::aco_ff::gV_fep
        < ::tinker::parameter::radius_types::R_MIN>
        (coords::float_type const E, coords::float_type const R, coords::float_type const r,
          coords::float_type const vout, coords::float_type& dV)
      {
        coords::float_type D6, D12, T2;
        coords::float_type T = R, D = r * r;
        T = T * T * T; // T^3
        T = T * T; // T^6
        T2 = T * T; // T^12
        D6 = D * D * D; //r^6 / D^3
        D6 = Config::get().fep.ljshift * (1 - vout) * (1 - vout) * T + D6; //r^6 shifted
        D12 = D6 * D6; //r^12 shifted
        coords::float_type V = vout * E * (T2 / D12 - 2 * T / D6);   //potential
        double numerator = T * (Config::get().fep.ljshift * (vout - 1) * (vout - 1) - 1) + r * r * r * r * r * r;
        double denominator = Config::get().fep.ljshift * (vout - 1) * (vout - 1) * T + r * r * r * r * r * r;
        dV = vout * E * 12.0 * T * r * r * r * r * r * numerator / (denominator * denominator * denominator);  //derivative
        return V;
      }

      /**calculate lenard-jones potential and gradient for oplsaa-forcefield (sigma-type);
      returns the energy
      @param E: epsilon-parameter
      @param R: sigma-parameter
      @param r: distance between the two atoms
      @param vout: lambda_vdw
      @param dV: reference to variable that saves absolute value of gradient*/
      template<> coords::float_type energy::interfaces::aco::aco_ff::gV_fep
        < ::tinker::parameter::radius_types::SIGMA>
        (coords::float_type const E, coords::float_type const R,
          coords::float_type const r, coords::float_type const vout,
          coords::float_type& dV)
      {
        coords::float_type T = R, D = r * r;
        T = T * T * T; // T^3
        T = T * T; // T^6
        D = D * D * D; // r^6
        D = Config::get().fep.ljshift * (1 - vout) * (1 - vout) * T + D; // r^6 shifted
        T /= D;
        coords::float_type V = vout * E * T;
        double numerator = R * R * R * R * R * R * (Config::get().fep.ljshift * (vout - 1) * (vout - 1) - 2) + r * r * r * r * r * r;
        dV = 6 * E * vout * R * R * R * R * R * R * r * r * r * r * r * numerator / (D * D * D);
        return V * (T - 1.0);
      }

      /**calculate non-bonding interactions and gradients between two atoms
      @param C: product of the charges
      @param E: epsilon-parameter or 4*epsilon-parameter
      @param R: sigma or Rmin parameter
      @param d: inverse distance between the two atoms
      @param e_c: reference to variable that saves coulomb-energy
      @param e_v: reference to variable that saves vdw-energy
      @param dE_c: reference to variable that saves coulomb gradient divided by distance
      @param dE_v: reference to variable that saves vdW gradient divided by distance*/
      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_QV
      (coords::float_type const C, coords::float_type const E,
        coords::float_type const R, coords::float_type const d,
        coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v)
      {
        coords::float_type dQ(0.0), dV(0.0);
        e_c += gQ(C, d, dQ);
        e_v += gV<RT>(E, R, d, dV);
        //division by distance because dQ and dV don't have a direction and get it by multiplying it with vector between atoms
        dE_c = dQ * d;
        dE_v = dV * d;
      }

      /**calculate non-bonding interactions and gradients between two atoms when one of them is IN or OUT (FEP)
     @param C: product of the charges
     @param E: epsilon-parameter or 4*epsilon-parameter
     @param R: sigma or Rmin parameter
     @param d: distance between the two atoms
     @param c_io: lambda_el
     @param v_io: lamda_vdw
     @param e_c: reference to variable that saves coulomb-energy
     @param e_v: reference to variable that saves vdw-energy
     @param dE_c: reference to variable that saves coulomb gradient divided by distance
     @param dE_v: reference to variable that saves vdW gradient divided by distance*/
      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_QV_fep
      (coords::float_type const C, coords::float_type const E,
        coords::float_type const R, coords::float_type const d,
        coords::float_type const c_io, coords::float_type const v_io,
        coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v)
      {
        coords::float_type dQ(0.0), dV(0.0);
        e_c += gQ_fep(C, d, c_io, dQ);
        e_v += gV_fep<RT>(E, R, d, v_io, dV);
        //division by distance because dQ and dV don't have a direction and get it by multiplying it with vector between atoms
        dE_c = dQ / d;
        dE_v = dV / d;
      }

      /**calculate non-bonding interactions and gradients between two atoms when a cutoff is applied
      @param C: product of the charges
      @param E: epsilon-parameter or 4*epsilon-parameter
      @param R: sigma or Rmin parameter
      @param d: inverse distance between the two atoms
      @param fQ: scaling factor for electrostatic interaction due to cutoff
      @param fV: scaling factor for vdw interaction due to cutoff
      @param e_c: reference to variable that saves coulomb-energy
      @param e_v: reference to variable that saves vdw-energy
      @param dE_c: reference to variable that saves coulomb gradient divided by distance
      @param dE_v: reference to variable that saves vdW gradient divided by distance*/
      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_QV_cutoff
      (coords::float_type const C, coords::float_type const E,
        coords::float_type const R, coords::float_type const d,
        coords::float_type const fQ, coords::float_type const fV,
        coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v)
      {
        auto const r = 1.0 / d;                              // distance between atoms
        double const& c = Config::get().energy.cutoff;       // cutoff distance
        double const& s = Config::get().energy.switchdist;   // distance where cutoff starts to kick in (only vdW)

        coords::float_type dQ(0.0), dV(0.0);
        auto const eQ = gQ(C, d, dQ);          // calculate coulomb energy and gradients without scaling
        auto const eV = gV<RT>(E, R, d, dV);   // calculate vdW energy and gradients without scaling
        e_c += eQ * fQ;                          // multiply coulomb energy with scaling factor
        e_v += eV * fV;                          // multiply coulomb energy with scaling factor

        auto dE_Q = dQ * fQ;                   // multiply coulomb gradients with scaling factor
        auto dE_V = dV * fV;                   // multiply vdW gradients with scaling factor

        if (r < c) {                 // for gradients: take into account that scaling factor is not constant but also changes with r
          dE_Q += eQ * (4 * r * (r * r - c * c)) / (c * c * c * c);
          if (r > s) dE_V += eV * (-12 * r * (c * c - r * r) * (r * r - s * s)) / ((c * c - s * s) * (c * c - s * s) * (c * c - s * s));
        }
        //division by distance because dE_Q and dE_V don't have a direction and get it by multiplying it with vector between atoms
        dE_c = dE_Q * d;
        dE_v = dE_V * d;
      }

      /**calculate non-bonding interactions and gradients between two atoms when a cutoff is applied
      and one of the atoms in IN or OUT in FEP
      @param C: product of the charges
      @param E: epsilon-parameter or 4*epsilon-parameter
      @param R: sigma or Rmin parameter
      @param d: distance between the two atoms
      @param c_out: lambda_el
      @param v_out: lamda_vdw
      @param fQ: scaling factor for electrostatic interaction due to cutoff
      @param fV: scaling factor for vdw interaction due to cutoff
      @param e_c: reference to variable that saves coulomb-energy
      @param e_v: reference to variable that saves vdw-energy
      @param dE_c: reference to variable that saves coulomb gradient divided by distance
      @param dE_v: reference to variable that saves vdW gradient divided by distance*/
      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_QV_fep_cutoff
      (coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const d,
        coords::float_type const c_out, coords::float_type const v_out, coords::float_type const fQ,
        coords::float_type const fV, coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v)
      {
        double const& c = Config::get().energy.cutoff;       // cutoff distance
        double const& s = Config::get().energy.switchdist;   // distance where cutoff starts to kick in (only vdW)

        coords::float_type dQ(0.0), dV(0.0);
        auto const eQ = gQ_fep(C, d, c_out, dQ);          // calculate coulomb energy and gradients without scaling
        auto const eV = gV_fep<RT>(E, R, d, v_out, dV);   // calculate vdW energy and gradients without scaling
        e_c += eQ * fQ;                          // multiply coulomb gradients with scaling factor
        e_v += eV * fV;                          // multiply vdW gradients with scaling factor

        auto dE_Q = dQ * fQ;                   // multiply coulomb gradients with scaling factor
        auto dE_V = dV * fV;                   // multiply vdW gradients with scaling factor

        if (d < c) {                 // for gradients: take into account that scaling factor is not constant but also changes with r
          dE_Q += eQ * (4 * d * (d * d - c * c)) / (c * c * c * c);
          if (d > s) dE_V += eV * (-12 * d * (c * c - d * d) * (d * d - s * s)) / ((c * c - s * s) * (c * c - s * s) * (c * c - s * s));
        }
        //division by distance because dQ and dV don't have a direction and get it by multiplying it with vector between atoms
        dE_c = dE_Q / d;
        dE_v = dE_V / d;
      }

      /**main function for calculating all non-bonding interactions*/
      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_nb(void)
      {
        part_energy[types::CHARGE] = 0.0;
        part_energy[types::VDW] = 0.0;
        part_grad[types::VDW].assign(part_grad[types::VDW].size(), coords::Cartesian_Point());
        part_grad[types::CHARGE].assign(part_grad[types::CHARGE].size(), coords::Cartesian_Point());

        coords->getFep().feptemp = energy::fepvect();
        for (auto& ia : coords->interactions()) ia.energy = 0.0;

        for (auto const& pairmatrix : refined.pair_matrices())
        {
          size_t const N(coords->interactions().size());
          for (size_t sub_ia_index(0u), row(0u), col(0u); sub_ia_index < N; ++sub_ia_index)
          {
            coords::float_type& e(coords->interactions(sub_ia_index).energy);
            coords::Representation_3D g_vdw(coords->interactions(sub_ia_index).grad);
            g_vdw.assign(coords->size(), coords::Cartesian_Point());
            coords::Representation_3D g_coul(coords->interactions(sub_ia_index).grad);
            g_coul.assign(coords->size(), coords::Cartesian_Point());
            std::vector< ::tinker::refine::types::nbpair> const& pl(pairmatrix.pair_matrix(sub_ia_index));
            scon::matrix< ::tinker::parameter::combi::vdwc, true> const& par(refined.vdwcm(pairmatrix.param_matrix_id));
            if (Config::get().md.fep)
            {
              if (Config::get().periodics.periodic)
              {
                if (coords->atoms().in_exists() && (coords->atoms().sub_in() == row || coords->atoms().sub_in() == col))
                  g_nb_QV_pairs_fep_io<RT, true, false>(e, g_vdw, g_coul, pl, par);
                else if (coords->atoms().out_exists() && (coords->atoms().sub_out() == row || coords->atoms().sub_out() == col))
                  g_nb_QV_pairs_fep_io<RT, true, true>(e, g_vdw, g_coul, pl, par);
                else
                  g_nb_QV_pairs_cutoff<RT, true>(e, g_vdw, g_coul, pl, par);
              }
              else  // no periodic boundaries
              {
                if (coords->atoms().in_exists() && (coords->atoms().sub_in() == row || coords->atoms().sub_in() == col))
                  g_nb_QV_pairs_fep_io<RT, false, false>(e, g_vdw, g_coul, pl, par);
                else if (coords->atoms().out_exists() && (coords->atoms().sub_out() == row || coords->atoms().sub_out() == col))
                  g_nb_QV_pairs_fep_io<RT, false, true>(e, g_vdw, g_coul, pl, par);
                else if (Config::get().energy.cutoff < 1000.0)
                  g_nb_QV_pairs_cutoff<RT, false>(e, g_vdw, g_coul, pl, par);
                else
                  g_nb_QV_pairs<RT>(e, g_vdw, g_coul, pl, par);
              }
            }
            else   // no fep
            {
              if (Config::get().periodics.periodic)
                g_nb_QV_pairs_cutoff<RT, true>(e, g_vdw, g_coul, pl, par);
              else if (Config::get().energy.cutoff < 1000.0)
                g_nb_QV_pairs_cutoff<RT, false>(e, g_vdw, g_coul, pl, par);
              else
                g_nb_QV_pairs<RT>(e, g_vdw, g_coul, pl, par);
            }
            if (col == row)
            {
              col = 0;
              ++row;
            }
            else ++col;

            part_grad[types::VDW] += g_vdw;
            part_grad[types::CHARGE] += g_coul;
          }
        }
        if (Config::get().md.fep)
        {
          coords->getFep().feptemp.dE = (coords->getFep().feptemp.e_c_l2 + coords->getFep().feptemp.e_vdw_l2) - (coords->getFep().feptemp.e_c_l1 + coords->getFep().feptemp.e_vdw_l1);
          coords->getFep().feptemp.dE_back = (coords->getFep().feptemp.e_c_l1 + coords->getFep().feptemp.e_vdw_l1) - (coords->getFep().feptemp.e_c_l0 + coords->getFep().feptemp.e_vdw_l0);
          coords->getFep().feptemp.dG = 0;
          coords->getFep().fepdata.push_back(coords->getFep().feptemp);
        }
      }

      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs
      (
        coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
        std::vector< ::tinker::refine::types::nbpair> const& pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const& params
      )
      {
        if (Config::get().general.single_charges) g_nb_QV_pairs_singleCharges<RT>(e_nb, grad_vdw, grad_coulomb, pairlist, params);
        else g_nb_QV_pairs_paramCharges<RT>(e_nb, grad_vdw, grad_coulomb, pairlist, params);
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC, bool ALCH_OUT>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io
      (
        coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
        std::vector< ::tinker::refine::types::nbpair> const& pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const& params
      )
      {
        if (Config::get().general.single_charges) g_nb_QV_pairs_fep_io_singleCharges<RT, PERIODIC, ALCH_OUT>(e_nb, grad_vdw, grad_coulomb, pairlist, params);
        else g_nb_QV_pairs_fep_io_paramCharges<RT, PERIODIC, ALCH_OUT>(e_nb, grad_vdw, grad_coulomb, pairlist, params);
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff
      (
        coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
        std::vector< ::tinker::refine::types::nbpair> const& pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const& params
      )
      {
        if (Config::get().general.single_charges) g_nb_QV_pairs_cutoff_singleCharges<RT, PERIODIC>(e_nb, grad_vdw, grad_coulomb, pairlist, params);
        else g_nb_QV_pairs_cutoff_paramCharges<RT, PERIODIC>(e_nb, grad_vdw, grad_coulomb, pairlist, params);
      }

      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_paramCharges
      (
        coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
        std::vector< ::tinker::refine::types::nbpair> const& pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const& params
      )
      {
        std::ptrdiff_t const M(pairlist.size());
        coords::float_type e_c(0.0), e_v(0.0);
#pragma omp parallel
        {
          coords::Representation_3D tmp_grad_vdw(grad_vdw.size());
          coords::Representation_3D tmp_grad_coul(grad_coulomb.size());
#pragma omp for reduction (+: e_c, e_v)
          for (std::ptrdiff_t i = 0; i < M; ++i)       // for every pair in pairlist
          {
            coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));
            coords::float_type const r = 1.0 / std::sqrt(dot(b, b));
            coords::float_type dE_c(0.0);
            coords::float_type dE_v(0.0);
            ::tinker::parameter::combi::vdwc const& p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));
            g_QV<RT>(p.C, p.E, p.R, r, e_c, e_v, dE_c, dE_v);
            auto grad_v = b * dE_v;
            auto grad_coul = b * dE_c;
            tmp_grad_vdw[pairlist[i].a] += grad_v;
            tmp_grad_vdw[pairlist[i].b] -= grad_v;
            tmp_grad_coul[pairlist[i].a] += grad_coul;
            tmp_grad_coul[pairlist[i].b] -= grad_coul;
          }
#pragma omp critical (nb_g_sum)
          {
            grad_vdw += tmp_grad_vdw;
            grad_coulomb += tmp_grad_coul;
          }
        }
        e_nb += e_c + e_v;

        part_energy[types::CHARGE] += e_c;
        part_energy[types::VDW] += e_v;
      }

      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_singleCharges
      (
        coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
        std::vector< ::tinker::refine::types::nbpair> const& pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const& params
      )
      {
        std::ptrdiff_t const M(pairlist.size());
        coords::float_type e_c(0.0), e_v(0.0);
#pragma omp parallel
        {
          coords::Representation_3D tmp_grad_vdw(grad_vdw.size());
          coords::Representation_3D tmp_grad_coul(grad_coulomb.size());
#pragma omp for reduction (+: e_c, e_v)
          for (std::ptrdiff_t i = 0; i < M; ++i)       // for every pair in pairlist
          {
            double ca = Config::get().coords.atom_charges[pairlist[i].a];
            double cb = Config::get().coords.atom_charges[pairlist[i].b];
            double current_c = ca * cb * cparams.general().electric;      // unit conversion
            if (refined.get_relation(pairlist[i].b, pairlist[i].a) == 3) current_c = current_c / cparams.general().chg_scale.value[3]; // 1,4 interactions are scaled down

            coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));
            coords::float_type const r = 1.0 / std::sqrt(dot(b, b));
            coords::float_type dE_c(0.0);
            coords::float_type dE_v(0.0);
            ::tinker::parameter::combi::vdwc const& p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));

            g_QV<RT>(current_c, p.E, p.R, r, e_c, e_v, dE_c, dE_v);  //calculate vdw and coulomb energy and gradients

            auto grad_v = b * dE_v;
            auto grad_coul = b * dE_c;
            tmp_grad_vdw[pairlist[i].a] += grad_v;
            tmp_grad_vdw[pairlist[i].b] -= grad_v;
            tmp_grad_coul[pairlist[i].a] += grad_coul;
            tmp_grad_coul[pairlist[i].b] -= grad_coul;
          }
#pragma omp critical (nb_g_sum)
          {
            grad_vdw += tmp_grad_vdw;
            grad_coulomb += tmp_grad_coul;
          }
        }
        e_nb += e_c + e_v;

        part_energy[types::CHARGE] += e_c;
        part_energy[types::VDW] += e_v;
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff_paramCharges
      (
        coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
        std::vector< ::tinker::refine::types::nbpair> const& pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const& params
      )
      {
        nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
        coords::float_type e_c(0.0), e_v(0.0);
        std::ptrdiff_t const M(pairlist.size());
#pragma omp parallel
        {
          coords::Representation_3D tmp_grad_vdw(grad_vdw.size());
          coords::Representation_3D tmp_grad_coul(grad_coulomb.size());
          coords::virial_t tempvir_vdw(coords::empty_virial());
          coords::virial_t tempvir_coul(coords::empty_virial());
#pragma omp for reduction (+: e_c, e_v)
          for (std::ptrdiff_t i = 0; i < M; ++i)  //for every pair in pairlist
          {
            coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));  //vector between the two atoms
            if (PERIODIC) boundary(b);  // for periodic boundaries: 
                              // if the absolute value of the distance in one of the coordinates is bigger than half the box size:
                              // subtract (or add) the box size
                              // => absolute value of the new box size is the smallest value between these atoms in any of the boxes
            coords::float_type const rr = dot(b, b);
            coords::float_type r(0.0), fQ(0.0), fV(0.0), dE_c(0.0), dE_v(0.0);
            if (!cutob.factors(rr, r, fQ, fV)) continue;   // cutoff applied? if yes: calculates scaling factors fQ (coulomb) and fV (vdW)
            r = 1.0 / r;
            ::tinker::parameter::combi::vdwc const& p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));   // get parameters for current pair
            g_QV_cutoff<RT>(p.C, p.E, p.R, r, fQ, fV, e_c, e_v, dE_c, dE_v);  //calculate vdw and coulomb energy and gradients
            auto const dist = b;
            // gradient dE/dr is getting a direction by muliplying it with vector between atoms
            auto grad_vdw = b * dE_v;
            auto grad_coul = b * dE_c;
            tmp_grad_vdw[pairlist[i].a] += grad_vdw;
            tmp_grad_vdw[pairlist[i].b] -= grad_vdw;
            tmp_grad_coul[pairlist[i].a] += grad_coul;
            tmp_grad_coul[pairlist[i].b] -= grad_coul;
            //Increment internal virial tensor
            // vdW
            coords::float_type vxx = grad_vdw.x() * dist.x();
            coords::float_type vyx = grad_vdw.x() * dist.y();
            coords::float_type vzx = grad_vdw.x() * dist.z();
            coords::float_type vyy = grad_vdw.y() * dist.y();
            coords::float_type vzy = grad_vdw.y() * dist.z();
            coords::float_type vzz = grad_vdw.z() * dist.z();
            tempvir_vdw[0][0] += vxx;
            tempvir_vdw[1][0] += vyx;
            tempvir_vdw[2][0] += vzx;
            tempvir_vdw[0][1] += vyx;
            tempvir_vdw[1][1] += vyy;
            tempvir_vdw[2][1] += vzy;
            tempvir_vdw[0][2] += vzx;
            tempvir_vdw[1][2] += vzy;
            tempvir_vdw[2][2] += vzz;
            // Coulomb
            vxx = grad_coul.x() * dist.x();
            vyx = grad_coul.x() * dist.y();
            vzx = grad_coul.x() * dist.z();
            vyy = grad_coul.y() * dist.y();
            vzy = grad_coul.y() * dist.z();
            vzz = grad_coul.z() * dist.z();
            tempvir_coul[0][0] += vxx;
            tempvir_coul[1][0] += vyx;
            tempvir_coul[2][0] += vzx;
            tempvir_coul[0][1] += vyx;
            tempvir_coul[1][1] += vyy;
            tempvir_coul[2][1] += vzy;
            tempvir_coul[0][2] += vzx;
            tempvir_coul[1][2] += vzy;
            tempvir_coul[2][2] += vzz;
          }
#pragma omp critical (nb_g_sum)
          {
            grad_vdw += tmp_grad_vdw;
            grad_coulomb += tmp_grad_coul;
            for (int i = 0; i <= 2; i++) {
              for (int k = 0; k <= 2; k++) {
                part_virial[VDW][i][k] += tempvir_vdw[i][k];
                part_virial[CHARGE][i][k] += tempvir_coul[i][k];
              }
            }
          }
        }
        e_nb += e_c + e_v;
        part_energy[types::CHARGE] += e_c;
        part_energy[types::VDW] += e_v;
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff_singleCharges
      (
        coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
        std::vector< ::tinker::refine::types::nbpair> const& pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const& params
      )
      {
        nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
        coords::float_type e_c(0.0), e_v(0.0);
        std::ptrdiff_t const M(pairlist.size());
#pragma omp parallel
        {
          coords::Representation_3D tmp_grad_vdw(grad_vdw.size());
          coords::Representation_3D tmp_grad_coul(grad_coulomb.size());
          coords::virial_t tempvir_vdw(coords::empty_virial());
          coords::virial_t tempvir_coul(coords::empty_virial());
#pragma omp for reduction (+: e_c, e_v)
          for (std::ptrdiff_t i = 0; i < M; ++i)  //for every pair in pairlist
          {
            coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));  //vector between the two atoms
            if (PERIODIC) boundary(b);  // for periodic boundaries: 
                              // if the absolute value of the distance in one of the coordinates is bigger than half the box size:
                              // subtract (or add) the box size
                              // => absolute value of the new box size is the smallest value between these atoms in any of the boxes
            coords::float_type const rr = dot(b, b);
            coords::float_type r(0.0), fQ(0.0), fV(0.0), dE_c(0.0), dE_v(0.0);
            if (!cutob.factors(rr, r, fQ, fV)) continue;   // cutoff applied? if yes: calculates scaling factors fQ (coulomb) and fV (vdW)
            r = 1.0 / r;
            double ca = Config::get().coords.atom_charges[pairlist[i].a];
            double cb = Config::get().coords.atom_charges[pairlist[i].b];
            double current_c = ca * cb * cparams.general().electric;  // unit conversion
            if (refined.get_relation(pairlist[i].b, pairlist[i].a) == 3) current_c = current_c / cparams.general().chg_scale.value[3]; // 1,4 interactions are scaled down
            ::tinker::parameter::combi::vdwc const& p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));   // get parameters for current pair

            g_QV_cutoff<RT>(current_c, p.E, p.R, r, fQ, fV, e_c, e_v, dE_c, dE_v);  //calculate vdw and coulomb energy and gradients

            auto const dist = b;
            // gradient dE/dr is getting a direction by muliplying it with vector between atoms
            auto grad_vdw = b * dE_v;
            auto grad_coul = b * dE_c;
            tmp_grad_vdw[pairlist[i].a] += grad_vdw;
            tmp_grad_vdw[pairlist[i].b] -= grad_vdw;
            tmp_grad_coul[pairlist[i].a] += grad_coul;
            tmp_grad_coul[pairlist[i].b] -= grad_coul;
            //Increment internal virial tensor
            // vdW
            coords::float_type vxx = grad_vdw.x() * dist.x();
            coords::float_type vyx = grad_vdw.x() * dist.y();
            coords::float_type vzx = grad_vdw.x() * dist.z();
            coords::float_type vyy = grad_vdw.y() * dist.y();
            coords::float_type vzy = grad_vdw.y() * dist.z();
            coords::float_type vzz = grad_vdw.z() * dist.z();
            tempvir_vdw[0][0] += vxx;
            tempvir_vdw[1][0] += vyx;
            tempvir_vdw[2][0] += vzx;
            tempvir_vdw[0][1] += vyx;
            tempvir_vdw[1][1] += vyy;
            tempvir_vdw[2][1] += vzy;
            tempvir_vdw[0][2] += vzx;
            tempvir_vdw[1][2] += vzy;
            tempvir_vdw[2][2] += vzz;
            // Coulomb
            vxx = grad_coul.x() * dist.x();
            vyx = grad_coul.x() * dist.y();
            vzx = grad_coul.x() * dist.z();
            vyy = grad_coul.y() * dist.y();
            vzy = grad_coul.y() * dist.z();
            vzz = grad_coul.z() * dist.z();
            tempvir_coul[0][0] += vxx;
            tempvir_coul[1][0] += vyx;
            tempvir_coul[2][0] += vzx;
            tempvir_coul[0][1] += vyx;
            tempvir_coul[1][1] += vyy;
            tempvir_coul[2][1] += vzy;
            tempvir_coul[0][2] += vzx;
            tempvir_coul[1][2] += vzy;
            tempvir_coul[2][2] += vzz;
          }
#pragma omp critical (nb_g_sum)
          {
            grad_vdw += tmp_grad_vdw;
            grad_coulomb += tmp_grad_coul;
            for (int i = 0; i <= 2; i++) {
              for (int k = 0; k <= 2; k++) {
                part_virial[VDW][i][k] += tempvir_vdw[i][k];
                part_virial[CHARGE][i][k] += tempvir_coul[i][k];
              }
            }
          }
        }
        e_nb += e_c + e_v;
        part_energy[types::CHARGE] += e_c;
        part_energy[types::VDW] += e_v;
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC, bool ALCH_OUT>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_paramCharges
      (
        coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
        std::vector< ::tinker::refine::types::nbpair> const& pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const& params
      )
      {
        nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
        coords::float_type e_c(0.0), e_v(0.0), e_c_l(0.0), e_vdw_l(0.0), e_c_dl(0.0), e_vdw_dl(0.0), e_c_ml(0.0), e_vdw_ml(0.0);
        fepvar const& fep = coords->getFep().window[coords->getFep().window[0].step];
        std::ptrdiff_t const M(pairlist.size());
#pragma omp parallel
        {
          coords::Representation_3D tmp_grad_vdw(grad_vdw.size());
          coords::Representation_3D tmp_grad_coul(grad_coulomb.size());
          coords::virial_t tempvir_vdw(coords::empty_virial());
          coords::virial_t tempvir_coul(coords::empty_virial());
#pragma omp for reduction (+: e_c, e_v, e_c_l, e_c_dl, e_vdw_l, e_vdw_dl, e_c_ml, e_vdw_ml)
          for (std::ptrdiff_t i = 0; i < M; ++i)      //for every pair in pairlist
          {
            coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));  //vector between atoms a and b
            if (PERIODIC) boundary(b);   // adjust vector to boundary conditions
            ::tinker::parameter::combi::vdwc const& p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));  // get parameters
            coords::float_type rr(dot(b, b)), dE_v(0.0), dE_c(0.0), Q(0.0), V(0.0);
            coords::Cartesian_Point dist;
            if (PERIODIC)
            {
              coords::float_type fQ(0.0), fV(0.0);
              coords::float_type r(0.0);
              if (!cutob.factors(rr, r, fQ, fV)) continue;
              g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.eout : fep.ein),
                (ALCH_OUT ? fep.vout : fep.vin), fQ, fV, Q, V, dE_c, dE_v);
              coords::float_type trash(0.0);
              g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.deout : fep.dein),
                (ALCH_OUT ? fep.dvout : fep.dvin), fQ, fV, e_c_dl, e_vdw_dl, trash, trash);
              g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.meout : fep.mein),
                (ALCH_OUT ? fep.mvout : fep.mvin), fQ, fV, e_c_ml, e_vdw_ml, trash, trash);
            }
            else    // not periodic
            {
              if (Config::get().energy.cutoff < 1000.0)
              {
                coords::float_type r(0.0);
                coords::float_type fQ(0.0), fV(0.0);
                if (cutob.factors(rr, r, fQ, fV))  //calculate r and see if r < cutoff
                {
                  g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.eout : fep.ein),
                    (ALCH_OUT ? fep.vout : fep.vin), fQ, fV, Q, V, dE_c, dE_v);  //calculate nb-energy(lambda)
                  coords::float_type trash(0.0);
                  g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.deout : fep.dein),
                    (ALCH_OUT ? fep.dvout : fep.dvin), fQ, fV, e_c_dl, e_vdw_dl, trash, trash);  //calculate nb-energy(lambda+dlambda)
                  g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.meout : fep.mein),
                    (ALCH_OUT ? fep.mvout : fep.mvin), fQ, fV, e_c_ml, e_vdw_ml, trash, trash);  //calculate nb-energy(lambda-dlambda)
                }
              }
              else  //if no cutoff
              {
                coords::float_type const r = sqrt(rr);
                g_QV_fep<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.eout : fep.ein),
                  (ALCH_OUT ? fep.vout : fep.vin), Q, V, dE_c, dE_v);  //calculate nb-energy(lambda)
                coords::float_type trash(0.0);
                g_QV_fep<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.deout : fep.dein),
                  (ALCH_OUT ? fep.dvout : fep.dvin), e_c_dl, e_vdw_dl, trash, trash); //calculate nb-energy(lambda+dlambda)
                g_QV_fep<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.meout : fep.mein),
                  (ALCH_OUT ? fep.mvout : fep.mvin), e_c_ml, e_vdw_ml, trash, trash); //calculate nb-energy(lambda-dlambda)
              }
            }
            dist = b;
            // gradient dE/dr is getting a direction by muliplying it with vector between atoms
            auto grad_vdw = b * dE_v;
            auto grad_coul = b * dE_c;
            tmp_grad_vdw[pairlist[i].a] += grad_vdw;
            tmp_grad_vdw[pairlist[i].b] -= grad_vdw;
            tmp_grad_coul[pairlist[i].a] += grad_coul;
            tmp_grad_coul[pairlist[i].b] -= grad_coul;
            e_c_l += Q;
            e_vdw_l += V;
            e_c += Q;
            e_v += V;
            //Increment internal virial tensor
            // vdW
            coords::float_type vxx = grad_vdw.x() * dist.x();
            coords::float_type vyx = grad_vdw.x() * dist.y();
            coords::float_type vzx = grad_vdw.x() * dist.z();
            coords::float_type vyy = grad_vdw.y() * dist.y();
            coords::float_type vzy = grad_vdw.y() * dist.z();
            coords::float_type vzz = grad_vdw.z() * dist.z();
            tempvir_vdw[0][0] += vxx;
            tempvir_vdw[1][0] += vyx;
            tempvir_vdw[2][0] += vzx;
            tempvir_vdw[0][1] += vyx;
            tempvir_vdw[1][1] += vyy;
            tempvir_vdw[2][1] += vzy;
            tempvir_vdw[0][2] += vzx;
            tempvir_vdw[1][2] += vzy;
            tempvir_vdw[2][2] += vzz;
            // Coulomb
            vxx = grad_coul.x() * dist.x();
            vyx = grad_coul.x() * dist.y();
            vzx = grad_coul.x() * dist.z();
            vyy = grad_coul.y() * dist.y();
            vzy = grad_coul.y() * dist.z();
            vzz = grad_coul.z() * dist.z();
            tempvir_coul[0][0] += vxx;
            tempvir_coul[1][0] += vyx;
            tempvir_coul[2][0] += vzx;
            tempvir_coul[0][1] += vyx;
            tempvir_coul[1][1] += vyy;
            tempvir_coul[2][1] += vzy;
            tempvir_coul[0][2] += vzx;
            tempvir_coul[1][2] += vzy;
            tempvir_coul[2][2] += vzz;
          }
#pragma omp critical (nb_g_sum)
          {
            grad_vdw += tmp_grad_vdw;
            grad_coulomb += tmp_grad_coul;
            for (int i = 0; i <= 2; i++) {
              for (int k = 0; k <= 2; k++) {
                part_virial[VDW][i][k] += tempvir_vdw[i][k];
                part_virial[CHARGE][i][k] += tempvir_coul[i][k];
              }
            }
          }
        }
        e_nb += e_c + e_v;
        coords->getFep().feptemp.e_c_l1 += e_c_l;    //lambda (Coulomb energy)
        coords->getFep().feptemp.e_c_l2 += e_c_dl;   //lambda + dlambda (Coulomb energy)
        coords->getFep().feptemp.e_c_l0 += e_c_ml;    //lambda - dlambda (Coulomb energy)
        coords->getFep().feptemp.e_vdw_l1 += e_vdw_l;  //lambda (vdW energy)
        coords->getFep().feptemp.e_vdw_l2 += e_vdw_dl;  //lambda + dlambda (vdW energy)
        coords->getFep().feptemp.e_vdw_l0 += e_vdw_ml;  //lambda - dlambda (vdW energy)
        part_energy[types::CHARGE] += e_c;
        part_energy[types::VDW] += e_v;
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC, bool ALCH_OUT>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_singleCharges
      (
        coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
        std::vector< ::tinker::refine::types::nbpair> const& pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const& params
      )
      {
        nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
        coords::float_type e_c(0.0), e_v(0.0), e_c_l(0.0), e_vdw_l(0.0), e_c_dl(0.0), e_vdw_dl(0.0), e_c_ml(0.0), e_vdw_ml(0.0);
        fepvar const& fep = coords->getFep().window[coords->getFep().window[0].step];
        std::ptrdiff_t const M(pairlist.size());
#pragma omp parallel
        {
          coords::Representation_3D tmp_grad_vdw(grad_vdw.size());
          coords::Representation_3D tmp_grad_coul(grad_coulomb.size());
          coords::virial_t tempvir_vdw(coords::empty_virial());
          coords::virial_t tempvir_coul(coords::empty_virial());
#pragma omp for reduction (+: e_c, e_v, e_c_l, e_c_dl, e_vdw_l, e_vdw_dl, e_c_ml, e_vdw_ml)
          for (std::ptrdiff_t i = 0; i < M; ++i)      //for every pair in pairlist
          {
            double ca = Config::get().coords.atom_charges[pairlist[i].a];
            double cb = Config::get().coords.atom_charges[pairlist[i].b];
            double current_c = ca * cb * cparams.general().electric;  // unit conversion
            if (refined.get_relation(pairlist[i].b, pairlist[i].a) == 3) current_c = current_c / cparams.general().chg_scale.value[3]; // 1,4 interactions are scaled down
            coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));  //vector between atoms a and b
            if (PERIODIC) boundary(b);   // adjust vector to boundary conditions
            ::tinker::parameter::combi::vdwc const& p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));  // get parameters
            coords::float_type rr(dot(b, b)), dE_v(0.0), dE_c(0.0), Q(0.0), V(0.0);
            coords::Cartesian_Point dist;
            if (PERIODIC)
            {
              coords::float_type fQ(0.0), fV(0.0);
              coords::float_type r(0.0);
              if (!cutob.factors(rr, r, fQ, fV)) continue;
              g_QV_fep_cutoff<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.eout : fep.ein),
                (ALCH_OUT ? fep.vout : fep.vin), fQ, fV, Q, V, dE_c, dE_v);
              coords::float_type trash(0.0);
              g_QV_fep_cutoff<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.deout : fep.dein),
                (ALCH_OUT ? fep.dvout : fep.dvin), fQ, fV, e_c_dl, e_vdw_dl, trash, trash);
              g_QV_fep_cutoff<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.meout : fep.mein),
                (ALCH_OUT ? fep.mvout : fep.mvin), fQ, fV, e_c_ml, e_vdw_ml, trash, trash);
            }
            else    // not periodic
            {
              if (Config::get().energy.cutoff < 1000.0)
              {
                coords::float_type r(0.0);
                coords::float_type fQ(0.0), fV(0.0);
                if (cutob.factors(rr, r, fQ, fV))  //calculate r and see if r < cutoff
                {
                  g_QV_fep_cutoff<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.eout : fep.ein),
                    (ALCH_OUT ? fep.vout : fep.vin), fQ, fV, Q, V, dE_c, dE_v);  //calculate nb-energy(lambda)
                  coords::float_type trash(0.0);
                  g_QV_fep_cutoff<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.deout : fep.dein),
                    (ALCH_OUT ? fep.dvout : fep.dvin), fQ, fV, e_c_dl, e_vdw_dl, trash, trash);  //calculate nb-energy(lambda+dlambda)
                  g_QV_fep_cutoff<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.meout : fep.mein),
                    (ALCH_OUT ? fep.mvout : fep.mvin), fQ, fV, e_c_ml, e_vdw_ml, trash, trash);  //calculate nb-energy(lambda-dlambda)
                }
              }
              else  //if no cutoff
              {
                coords::float_type const r = sqrt(rr);
                g_QV_fep<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.eout : fep.ein),
                  (ALCH_OUT ? fep.vout : fep.vin), Q, V, dE_c, dE_v);  //calculate nb-energy(lambda)
                coords::float_type trash(0.0);
                g_QV_fep<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.deout : fep.dein),
                  (ALCH_OUT ? fep.dvout : fep.dvin), e_c_dl, e_vdw_dl, trash, trash); //calculate nb-energy(lambda+dlambda)
                g_QV_fep<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.meout : fep.mein),
                  (ALCH_OUT ? fep.mvout : fep.mvin), e_c_ml, e_vdw_ml, trash, trash); //calculate nb-energy(lambda-dlambda)
              }
            }
            dist = b;
            // gradient dE/dr is getting a direction by muliplying it with vector between atoms
            auto grad_vdw = b * dE_v;
            auto grad_coul = b * dE_c;
            tmp_grad_vdw[pairlist[i].a] += grad_vdw;
            tmp_grad_vdw[pairlist[i].b] -= grad_vdw;
            tmp_grad_coul[pairlist[i].a] += grad_coul;
            tmp_grad_coul[pairlist[i].b] -= grad_coul;
            e_c_l += Q;
            e_vdw_l += V;
            e_c += Q;
            e_v += V;
            //Increment internal virial tensor
            // vdW
            coords::float_type vxx = grad_vdw.x() * dist.x();
            coords::float_type vyx = grad_vdw.x() * dist.y();
            coords::float_type vzx = grad_vdw.x() * dist.z();
            coords::float_type vyy = grad_vdw.y() * dist.y();
            coords::float_type vzy = grad_vdw.y() * dist.z();
            coords::float_type vzz = grad_vdw.z() * dist.z();
            tempvir_vdw[0][0] += vxx;
            tempvir_vdw[1][0] += vyx;
            tempvir_vdw[2][0] += vzx;
            tempvir_vdw[0][1] += vyx;
            tempvir_vdw[1][1] += vyy;
            tempvir_vdw[2][1] += vzy;
            tempvir_vdw[0][2] += vzx;
            tempvir_vdw[1][2] += vzy;
            tempvir_vdw[2][2] += vzz;
            // Coulomb
            vxx = grad_coul.x() * dist.x();
            vyx = grad_coul.x() * dist.y();
            vzx = grad_coul.x() * dist.z();
            vyy = grad_coul.y() * dist.y();
            vzy = grad_coul.y() * dist.z();
            vzz = grad_coul.z() * dist.z();
            tempvir_coul[0][0] += vxx;
            tempvir_coul[1][0] += vyx;
            tempvir_coul[2][0] += vzx;
            tempvir_coul[0][1] += vyx;
            tempvir_coul[1][1] += vyy;
            tempvir_coul[2][1] += vzy;
            tempvir_coul[0][2] += vzx;
            tempvir_coul[1][2] += vzy;
            tempvir_coul[2][2] += vzz;
          }
#pragma omp critical (nb_g_sum)
          {
            grad_vdw += tmp_grad_vdw;
            grad_coulomb += tmp_grad_coul;
            for (int i = 0; i <= 2; i++) {
              for (int k = 0; k <= 2; k++) {
                part_virial[VDW][i][k] += tempvir_vdw[i][k];
                part_virial[CHARGE][i][k] += tempvir_coul[i][k];
              }
            }
          }
        }
        e_nb += e_c + e_v;
        coords->getFep().feptemp.e_c_l1 += e_c_l;    //lambda (Coulomb energy)
        coords->getFep().feptemp.e_c_l2 += e_c_dl;   //lambda + dlambda (Coulomb energy)
        coords->getFep().feptemp.e_c_l0 += e_c_ml;    //lambda - dlambda (Coulomb energy)
        coords->getFep().feptemp.e_vdw_l1 += e_vdw_l;  //lambda (vdW energy)
        coords->getFep().feptemp.e_vdw_l2 += e_vdw_dl;  //lambda + dlambda (vdW energy)
        coords->getFep().feptemp.e_vdw_l0 += e_vdw_ml;  //lambda - dlambda (vdW energy)
        part_energy[types::CHARGE] += e_c;
        part_energy[types::VDW] += e_v;
      }

    }
  }
}

////////////////////////////////////////////////////
// And now this:
//
// if templates are defined in a different file than the one they are declared there must be a declaration for each type that the template is used with
// see https://stackoverflow.com/questions/115703/storing-c-template-function-definitions-in-a-cpp-file
////////////////////////////////////////////////////
template void energy::interfaces::aco::aco_ff::g_nb< ::tinker::parameter::radius_types::R_MIN >(void);
template void energy::interfaces::aco::aco_ff::g_nb< ::tinker::parameter::radius_types::SIGMA >(void);

template coords::float_type energy::interfaces::aco::aco_ff::f_12<0>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_12<1>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_13_a<0>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_13_a<1>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_13_u<0>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_13_u<1>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_14<0>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_14<1>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_it<0>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_it<1>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_imp<0>(void);
template coords::float_type energy::interfaces::aco::aco_ff::f_imp<1>(void);

template coords::float_type energy::interfaces::aco::aco_ff::eV< ::tinker::parameter::radius_types::R_MIN >
(coords::float_type const E, coords::float_type const R, coords::float_type const r);

template coords::float_type energy::interfaces::aco::aco_ff::eV< ::tinker::parameter::radius_types::SIGMA >
(coords::float_type const E, coords::float_type const R, coords::float_type const r);

template coords::float_type energy::interfaces::aco::aco_ff::gV< ::tinker::parameter::radius_types::R_MIN >
(coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type& dV);

template coords::float_type energy::interfaces::aco::aco_ff::gV< ::tinker::parameter::radius_types::SIGMA >
(coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type& dV);

template coords::float_type energy::interfaces::aco::aco_ff::gV_fep< ::tinker::parameter::radius_types::R_MIN >
(coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type const factor, coords::float_type& dV);

template coords::float_type energy::interfaces::aco::aco_ff::gV_fep< ::tinker::parameter::radius_types::SIGMA >
(coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type const factor, coords::float_type& dV);

template void energy::interfaces::aco::aco_ff::g_QV< ::tinker::parameter::radius_types::R_MIN >
(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const d,
  coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

template void energy::interfaces::aco::aco_ff::g_QV< ::tinker::parameter::radius_types::SIGMA >
(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const d,
  coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

template void energy::interfaces::aco::aco_ff::g_QV_fep< ::tinker::parameter::radius_types::R_MIN >
(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r,
  coords::float_type const c_out, coords::float_type const v_out,
  coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

template void energy::interfaces::aco::aco_ff::g_QV_fep< ::tinker::parameter::radius_types::SIGMA >
(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r,
  coords::float_type const c_out, coords::float_type const v_out,
  coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

template void energy::interfaces::aco::aco_ff::g_QV_cutoff< ::tinker::parameter::radius_types::R_MIN >
(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r,
  coords::float_type const fQ, coords::float_type const fV, coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

template void energy::interfaces::aco::aco_ff::g_QV_cutoff< ::tinker::parameter::radius_types::SIGMA >
(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r,
  coords::float_type const fQ, coords::float_type const fV, coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

template void energy::interfaces::aco::aco_ff::g_QV_fep_cutoff< ::tinker::parameter::radius_types::R_MIN >
(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r,
  coords::float_type const c_out, coords::float_type const v_out, coords::float_type const fQ, coords::float_type const fV,
  coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

template void energy::interfaces::aco::aco_ff::g_QV_fep_cutoff< ::tinker::parameter::radius_types::SIGMA >
(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r,
  coords::float_type const c_out, coords::float_type const v_out, coords::float_type const fQ, coords::float_type const fV,
  coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs< ::tinker::parameter::radius_types::R_MIN>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs< ::tinker::parameter::radius_types::SIGMA>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff< ::tinker::parameter::radius_types::R_MIN, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff< ::tinker::parameter::radius_types::SIGMA, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff< ::tinker::parameter::radius_types::R_MIN, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff< ::tinker::parameter::radius_types::SIGMA, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io< ::tinker::parameter::radius_types::R_MIN, true, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io< ::tinker::parameter::radius_types::SIGMA, true, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io< ::tinker::parameter::radius_types::R_MIN, false, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io< ::tinker::parameter::radius_types::SIGMA, false, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io< ::tinker::parameter::radius_types::R_MIN, true, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io< ::tinker::parameter::radius_types::SIGMA, true, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io< ::tinker::parameter::radius_types::R_MIN, false, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io< ::tinker::parameter::radius_types::SIGMA, false, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_singleCharges< ::tinker::parameter::radius_types::R_MIN, true, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_singleCharges< ::tinker::parameter::radius_types::SIGMA, true, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_singleCharges< ::tinker::parameter::radius_types::R_MIN, false, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_singleCharges< ::tinker::parameter::radius_types::SIGMA, false, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_singleCharges< ::tinker::parameter::radius_types::R_MIN, true, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_singleCharges< ::tinker::parameter::radius_types::SIGMA, true, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_singleCharges< ::tinker::parameter::radius_types::R_MIN, false, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_singleCharges< ::tinker::parameter::radius_types::SIGMA, false, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_paramCharges< ::tinker::parameter::radius_types::R_MIN, true, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_paramCharges< ::tinker::parameter::radius_types::SIGMA, true, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_paramCharges< ::tinker::parameter::radius_types::R_MIN, false, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_paramCharges< ::tinker::parameter::radius_types::SIGMA, false, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_paramCharges< ::tinker::parameter::radius_types::R_MIN, true, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_paramCharges< ::tinker::parameter::radius_types::SIGMA, true, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_paramCharges< ::tinker::parameter::radius_types::R_MIN, false, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io_paramCharges< ::tinker::parameter::radius_types::SIGMA, false, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff_singleCharges< ::tinker::parameter::radius_types::R_MIN, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff_singleCharges< ::tinker::parameter::radius_types::SIGMA, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff_singleCharges< ::tinker::parameter::radius_types::R_MIN, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff_singleCharges< ::tinker::parameter::radius_types::SIGMA, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff_paramCharges< ::tinker::parameter::radius_types::R_MIN, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff_paramCharges< ::tinker::parameter::radius_types::SIGMA, true>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff_paramCharges< ::tinker::parameter::radius_types::R_MIN, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff_paramCharges< ::tinker::parameter::radius_types::SIGMA, false>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_paramCharges< ::tinker::parameter::radius_types::R_MIN>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_paramCharges< ::tinker::parameter::radius_types::SIGMA>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_singleCharges< ::tinker::parameter::radius_types::R_MIN>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

template void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_singleCharges< ::tinker::parameter::radius_types::SIGMA>
(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb, std::vector< ::tinker::refine::types::nbpair> const& pairs,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);