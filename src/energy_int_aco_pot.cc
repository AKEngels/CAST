#include <cmath>
#include <stddef.h>
#include <stdexcept>
#include <cstdlib>
//#include <fftw3.h>
#include "energy_int_aco.h"
#include "configuration.h"
#include "scon_utility.h"
//#include "gbsa.h"

#define SUPERPI 3.141592653589793238
#define SQRTPI  sqrt(3.141592653589793238)
/****************************************
*                                       *
*                                       *
*               Generics                *
*                                       *
*                                       *
*****************************************/

coords::float_type energy::interfaces::aco::aco_ff::e (void)
{
  pre();
  // calc with derivates 0
  calc<0>();
  post();
  return energy;
}

coords::float_type energy::interfaces::aco::aco_ff::g (void)
{
  pre();
  // calc with derivatives 1
  calc<1>();
  post();
  return energy;
}

coords::float_type energy::interfaces::aco::aco_ff::h (void)
{
  pre();
  //... 
  post();
  return energy;
}

coords::float_type energy::interfaces::aco::aco_ff::o (void)
{
  throw std::runtime_error("aco_ff doesn't provide any optimization routines.");
}

template<size_t DERIV> 
void energy::interfaces::aco::aco_ff::calc (void)
{
#if defined (_OPENMP)
    #pragma omp parallel sections
    {
      #pragma omp section
        part_energy[types::BOND]       = f_12<DERIV>();
      #pragma omp section
        part_energy[types::ANGLE]      = f_13_a<DERIV>();
      #pragma omp section
        part_energy[types::UREY]       = f_13_u<DERIV>();
      #pragma omp section
        part_energy[types::TORSION]    = f_14<DERIV>();
      #pragma omp section
        part_energy[types::IMPTORSION] = f_it<DERIV>();
      #pragma omp section
        part_energy[types::IMPROPER]   = f_imp<DERIV>();
    }
#else
    part_energy[types::BOND]         = f_12<DERIV>();
    part_energy[types::ANGLE]        = f_13_a<DERIV>();
    part_energy[types::UREY]         = f_13_u<DERIV>();
    part_energy[types::TORSION]      = f_14<DERIV>();
    part_energy[types::IMPTORSION]   = f_it<DERIV>();
    part_energy[types::IMPROPER]     = f_imp<DERIV>();
#endif

	if (cparams.radiustype() == ::tinker::parameter::radius_types::R_MIN)
	{
		/*if (Config::get().energy.pme == true) g_nb_pme< ::tinker::parameter::radius_types::R_MIN>();
		else*/
        g_nb< ::tinker::parameter::radius_types::R_MIN>();
	}
	else
		g_nb< ::tinker::parameter::radius_types::SIGMA>();

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
        for (auto const & bond : refined.bonds())
        {
          auto const bv(coords->xyz(bond.atoms[0]) - coords->xyz(bond.atoms[1]));
          auto const d = len(bv);
          auto const r = d - bond.ideal;
          E += bond.force*r*r;
          if (abs(r) > 0.5) integrity = false;
        }
        return E;
      }

      //! Bonding Energy+Gradients
      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_12<1> (void)
      {
        using std::abs;
        using scon::len;
        coords::float_type E(0.0);
        for (auto const & bond : refined.bonds())
        {
          auto const bv(coords->xyz(bond.atoms[0]) - coords->xyz(bond.atoms[1])); // r_ij (i=1, j=2)
          auto const d = len(bv);
          auto const r = d - bond.ideal;
          auto dE = bond.force*r;
          E += dE*r;  // kcal/mol
          dE *= 2;  // kcal/(mol*Angstrom)
          if (abs(d) > 0.0) 
          {
            if (abs(r) > 0.5) integrity = false;
            dE /= d;  // kcal/(mol*A^2)
            auto const gv = bv*dE;   // "force" on atom i due to atom j (kcal/(mol*A))
            part_grad[BOND][bond.atoms[0]] += gv;
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
          else integrity = false;
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
      coords::float_type energy::interfaces::aco::aco_ff::f_13_a<0> (void)
      {
        using std::abs;
        coords::float_type E(0.0);
        for (auto const & angle : refined.angles())
        {
          auto const 
            av1(coords->xyz(angle.atoms[0]) - coords->xyz(angle.atoms[1])),
            av2(coords->xyz(angle.atoms[2]) - coords->xyz(angle.atoms[1]));
          auto const d(scon::angle(av1, av2).degrees() - angle.ideal);
          auto const r(d*SCON_PI180);
          E += angle.force*r*r;
          if (abs(d) > 20.0) integrity = false;
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_13_a<1> (void)
      {
        using scon::angle;
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const & angle : refined.angles())
        {
          auto const
            av1(coords->xyz(angle.atoms[0]) - coords->xyz(angle.atoms[1])),
            av2(coords->xyz(angle.atoms[2]) - coords->xyz(angle.atoms[1]));
          auto const d = scon::angle(av1, av2).degrees() - angle.ideal;
          auto const r = d*SCON_PI180;
          /*std::cout << angle.atoms[0] << ", " << angle.atoms[1] << ", " << angle.atoms[2] << "\n";
          std::cout << "Delta deg: " << d << ", delta rad: " << r << ", ideal: " << angle.ideal << ", force = " << angle.force << '\n';
          std::cout << scon::angle(av1, av2).degrees() << "; " << angle.force*r*r << "\n";*/
          auto dE = angle.force*r;
          E += dE*r;
          //std::cout << "A " << angle.atoms[0] << "-" << angle.atoms[1] << "-" << angle.atoms[2] << ": ";
          
          coords::Cartesian_Point const cv(cross(av1, av2));
          coords::float_type const cvl(len(cv));
          if (abs(cvl) > 0.0)
          {
            if (abs(d) > 20.0) integrity = false;
            dE *= 2.0/cvl;
            //std::cout << dE << ", " << cvl << '\n';
            coords::Cartesian_Point const gv1(cross(av1, cv) * (dE / dot(av1, av1)));
            coords::Cartesian_Point const gv2(cross(cv, av2) * (dE / dot(av2, av2)));
            //std::cout << gv1 << ", " << gv2 << ", " << -(gv1 + gv2) << '\n';
            //std::cout << part_grad[ANGLE][angle.atoms[0]] << ", " << part_grad[ANGLE][angle.atoms[1]] << ", " << part_grad[ANGLE][angle.atoms[2]] << '\n';
            part_grad[ANGLE][angle.atoms[0]] += gv1;
            part_grad[ANGLE][angle.atoms[2]] += gv2;
            part_grad[ANGLE][angle.atoms[1]] += -(gv1+gv2);
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
            integrity = false;
          }
        }
       /* std::cout << "Virial sum" << std::endl;
        std::cout << part_virial[ANGLE][0][0] + part_virial[BOND][0][0] << "   " << part_virial[ANGLE][1][0] + part_virial[BOND][1][0] << "   " << part_virial[ANGLE][2][0] + part_virial[BOND][2][0] << std::endl;
        std::cout << part_virial[ANGLE][0][1] + part_virial[BOND][0][1] << "   " << part_virial[ANGLE][1][1] + part_virial[BOND][1][1] << "   " << part_virial[ANGLE][2][1] + part_virial[BOND][2][1] << std::endl;
        std::cout << part_virial[ANGLE][0][2] + part_virial[BOND][0][2] << "   " << part_virial[ANGLE][1][2] + part_virial[BOND][1][2] << "   " << part_virial[ANGLE][2][2] + part_virial[BOND][2][2] << std::endl;*/
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
      coords::float_type energy::interfaces::aco::aco_ff::f_13_u<0> (void)
      {
        using scon::len;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const & urey : refined.ureys())
        {
          coords::Cartesian_Point const bv = 
            coords->xyz(urey.atoms[0]) - coords->xyz(urey.atoms[1]);
          coords::float_type const d = len(bv);
          coords::float_type const r = d - urey.ideal;
          E += urey.force*r*r;
          if (abs(d) < 1.0e-8) integrity = false;
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_13_u<1> (void)
      {
        using scon::len;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const & urey : refined.ureys())
        {
          coords::Cartesian_Point const bv(coords->xyz(urey.atoms[0]) - coords->xyz(urey.atoms[1]));
          coords::float_type const d = len(bv);
          coords::float_type const r = d - urey.ideal;
          coords::float_type dE = urey.force*r;
          E += dE*r;
          dE *= 2;
          if (abs(d) > 0.0)
          {
            dE /= d;
            coords::Cartesian_Point const gv = bv*dE;
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

      /********************************
      *                                *
      *                                *
      *  Torsion  Potential            *
      *  Energy/Gradients/Hessians     *
      *                                *
      *                                *
      *********************************/

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_14<0> (void)
      {
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const & torsion : refined.torsions())
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
          coords::float_type const tlul = sqrt(tl2*ul2);
          coords::float_type const r12 = len(b12);
          // cross of cross
          coords::Cartesian_Point const tu = cross(t, u);
          // scalar and length variations
          coords::float_type const cos_scalar0 = dot(t, u);
          coords::float_type const cos_scalar1 = tlul;
          coords::float_type const sin_scalar0 = dot(b12, tu);
          coords::float_type const sin_scalar1 = r12*tlul;
          // check whether 
          if (abs(cos_scalar1) < 1.0e-8 || abs(sin_scalar1) < 1.0e-8)
          {
            integrity = false;
          }
          // Get multiple sine and cosine values
          coords::float_type cos[7], sin[7];
          cos[1] = cos_scalar0/cos_scalar1;
          sin[1] = sin_scalar0/sin_scalar1;
          for (std::size_t j(2U); j <= torsion.p.max_order; ++j)
          {
            std::size_t const k = j - 1;
            sin[j] = sin[k]*cos[1] + cos[k]*sin[1];
            cos[j] = cos[k]*cos[1] - sin[k]*sin[1];
          }

          coords::float_type tE(0.0);
          //cout << "Number?:  " << torsions[i].paramPtr->n << std::endl;
          for (std::size_t j(0U); j < torsion.p.number; ++j)
          {
            coords::float_type const F = torsion.p.force[j]*cparams.torsionunit();
            std::size_t const k = torsion.p.order[j];
            coords::float_type const l = std::abs(torsion.p.ideal[j]) > 0.0 ? -1.0 : 1.0;
            tE += F * (1.0 + cos[k]*l);
          }
          E += tE;
        }
        return E;
      }

      template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_14<1> (void)
      {
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt;
        using std::abs;
        coords::float_type E(0.0);
        for (auto const & torsion : refined.torsions())
        {

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
          coords::float_type const tlul = sqrt(tl2*ul2);
          coords::float_type const r12 = len(b12);

          coords::Cartesian_Point const tu = cross(t, u);

          coords::float_type const cos_scalar0 = dot(t, u);
          coords::float_type const cos_scalar1 = tlul;

          coords::float_type const sin_scalar0 = dot(b12, tu);
          coords::float_type const sin_scalar1 = r12*tlul;

          if (abs(cos_scalar1) < 1.0e-8 || abs(sin_scalar1) < 1.0e-8)
          {
            integrity = false;
          }

          coords::float_type cos[7], sin[7];
          cos[1] = cos_scalar0/cos_scalar1;
          sin[1] = sin_scalar0/sin_scalar1;

          for (std::size_t j(2U); j <= torsion.p.max_order; ++j)
          {
            std::size_t const k = j - 1;
            sin[j] = sin[k]*cos[1] + cos[k]*sin[1];
            cos[j] = cos[k]*cos[1] - sin[k]*sin[1];
          }

          coords::float_type tE(0.0), dE(0.0);
          //cout << "Number?:  " << torsions[i].paramPtr->n << std::endl;
          for (std::size_t j(0U); j < torsion.p.number; ++j)
          {
            coords::float_type const F = torsion.p.force[j]*cparams.torsionunit();
            std::size_t const k = torsion.p.order[j];
            coords::float_type const l = std::abs(torsion.p.ideal[j]) > 0.0 ? -1.0 : 1.0;
            tE += F * (1.0 + cos[k]*l);
            dE += -static_cast<coords::float_type>(k) * F * sin[k] * l;
          }
          E += tE;

          coords::Cartesian_Point const dt(cross(t, b12) * (dE / (tl2*r12)));
          coords::Cartesian_Point const du(cross(u, b12) * (-dE / (ul2*r12)));

          coords::Cartesian_Point const vir1 = cross(dt, b12);
          coords::Cartesian_Point const vir2 = cross(b02, dt) + cross(du, b23);
          coords::Cartesian_Point const vir3 = cross(dt, b01) + cross(b13, du);
          coords::Cartesian_Point const vir4 = cross(du, b12);

          part_grad[TORSION][torsion.atoms[0]] += vir1;
          part_grad[TORSION][torsion.atoms[1]] += vir2;
          part_grad[TORSION][torsion.atoms[2]] += vir3;
          part_grad[TORSION][torsion.atoms[3]] += vir4;

          //increment internal virial tensor
          coords::float_type const vxx = b12.x()*(vir3.x() + vir4.x()) - 
            b01.x()*vir1.x() + b23.x()*vir4.x();
          coords::float_type const vyx = b12.y()*(vir3.x() + vir4.x()) - 
            b01.y()*vir1.x() + b23.y()*vir4.x();
          coords::float_type const vzx = b12.z()*(vir3.x() + vir4.x()) - 
            b01.z()*vir1.x() + b23.z()*vir4.x();
          coords::float_type const vyy = b12.y()*(vir3.y() + vir4.y()) - 
            b01.y()*vir1.y() + b23.y()*vir4.y();
          coords::float_type const vzy = b12.z()*(vir3.y() + vir4.y()) - 
            b01.z()*vir1.y() + b23.z()*vir4.y();
          coords::float_type const vzz = b12.z()*(vir3.z() + vir4.z()) - 
            b01.z()*vir1.z() + b23.z()*vir4.z();
          part_virial[TORSION][0][0] += vxx;
          part_virial[TORSION][1][0] += vyx;
          part_virial[TORSION][2][0] += vzx;
          part_virial[TORSION][0][1] += vyx;
          part_virial[TORSION][1][1] += vyy;
          part_virial[TORSION][2][1] += vzy;
          part_virial[TORSION][0][2] += vzx;
          part_virial[TORSION][1][2] += vzy;
          part_virial[TORSION][2][2] += vzz;

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
		  for (auto & imptor : refined.imptors())   //for every improper torsion
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
			  coords::float_type const tlul = sqrt(tl2*ul2);
			  coords::float_type const r12 = len(cb);
			  coords::Cartesian_Point const tu = cross(t, u);

			  coords::float_type const cosine = dot(t, u) / tlul;
			  coords::float_type const sine = dot(cb, tu) / (r12*tlul);

			  auto v2 = imptor.p.force[0];   // I don't know if this is correct (originally there was something more in this line)
			  auto c2 = cos(imptor.p.ideal[0] * SCON_PI180);  // why index 0 everywhere?
			  auto s2 = sin(imptor.p.ideal[0] * SCON_PI180);

			  auto cosine2 = cosine*cosine - sine*sine;
			  auto sine2 = 2.0 * cosine * sine;

			  auto phi2 = 1.0 + (cosine2*c2 + sine2*s2);
			  auto dphi2 = 2.0 * (cosine2*s2 - sine2*c2);

			  auto dedphi = 0.5 * (v2*dphi2);

			  double E_add = 0.5 * (v2*phi2);
			  E += E_add;
		  }
		  return E;
	  }

	  template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_it<1> (void)
      {
		  using scon::cross;
		  using scon::dot;
		  using scon::len;
		  using std::sqrt;
		  using std::abs;
		  coords::float_type E(0.0);
		  for (auto & imptor : refined.imptors())   //for every improper torsion
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
			  coords::float_type const tlul = sqrt(tl2*ul2);
			  coords::float_type const r12 = len(cb);
			  coords::Cartesian_Point const tu = cross(t, u);

			  coords::float_type const cosine = dot(t, u) / tlul;
			  coords::float_type const sine = dot(cb, tu) / (r12*tlul);

			  auto v2 = imptor.p.force[0];   // I don't know if this is correct (originally there was something more in this line)
			  auto c2 = cos(imptor.p.ideal[0] * SCON_PI180);  // why index 0 everywhere?
			  auto s2 = sin(imptor.p.ideal[0] * SCON_PI180);

			  auto cosine2 = cosine*cosine - sine*sine;
			  auto sine2 = 2.0 * cosine * sine;

			  auto phi2 = 1.0 + (cosine2*c2 + sine2*s2);
			  auto dphi2 = 2.0 * (cosine2*s2 - sine2*c2);

			  auto dedphi = 0.5 * (v2*dphi2);

			  double E_add = 0.5 * (v2*phi2);
			  E += E_add;

			  //gradient calculation
			  auto const ca(coords->xyz(imptor.center) - coords->xyz(imptor.ligand[0]));
			  auto const db(coords->xyz(imptor.twist) - coords->xyz(imptor.ligand[1]));

			  auto const dt = cross(t, cb) * (dedphi / (tl2*r12));
			  auto const du = cross(u, cb) * (-dedphi / (ul2*r12));

			  auto vir1 = cross(dt, ba) + cross(db, du);
			  auto vir2 = cross(dt, cb);
			  auto vir3 = cross(ca, dt) + cross(du, dc);
			  auto vir4 = cross(du, cb);

			  part_grad[IMPTORSION][imptor.ligand[0]] += vir2;
			  part_grad[IMPTORSION][imptor.ligand[1]] += vir3;
			  part_grad[IMPTORSION][imptor.center] += vir1;
			  part_grad[IMPTORSION][imptor.twist] += vir4;

			  //calculation of virial tensors (copied from function f_imp<1>)
			  auto vxx = cb.x()*(vir3.x() + vir4.x()) - ba.x()*vir1.x() + dc.x()*vir4.x();
			  auto vyx = cb.y()*(vir3.x() + vir4.x()) - ba.y()*vir1.x() + dc.y()*vir4.x();
			  auto vzx = cb.z()*(vir3.x() + vir4.x()) - ba.z()*vir1.x() + dc.z()*vir4.x();
			  auto vyy = cb.y()*(vir3.y() + vir4.y()) - ba.y()*vir1.y() + dc.y()*vir4.y();
			  auto vzy = cb.z()*(vir3.y() + vir4.y()) - ba.z()*vir1.y() + dc.z()*vir4.y();
			  auto vzz = cb.z()*(vir3.z() + vir4.z()) - ba.z()*vir1.z() + dc.z()*vir4.z();

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
		  for (auto & improper : refined.impropers())
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
			  coords::float_type const rt2(dot(t, t)), ru2(dot(u, u)), rtru = sqrt(rt2*ru2);
			  if (abs(rtru) > 0.0)
			  {
				  coords::float_type const rcb(len(cb));
				  coords::float_type const cosine(min(1.0, max(-1.0, (dot(t, u) / rtru))));
				  coords::float_type const sine(dot(cb, tu) / (rcb*rtru));
				  coords::float_type const angle(sine < 0.0 ?
					  -acos(cosine)*SCON_180PI : acos(cosine)*SCON_180PI);
				  coords::float_type da((abs(angle + improper.p.ideal[0U]) < abs(angle - improper.p.ideal[0U])) ?
					  angle + improper.p.ideal[0U] : angle - improper.p.ideal[0U]);
				  if (da > 180.0) da -= 360.0;
				  if (da < -180.0) da += 360.0;
				  da *= SCON_PI180;
				  coords::float_type dE = improper.p.force[0] * da;
				  E += dE*da;
			  }
		  }
		  return E;
	  }

	  template<>
      coords::float_type energy::interfaces::aco::aco_ff::f_imp<1> (void)
      {
        using scon::cross;
        using scon::dot;
        using scon::len;
        using std::sqrt; using std::acos;
        using std::abs;
        using std::min; using std::max;
        coords::float_type E(0.0);
        coords::float_type vxx, vyx, vzx, vyy, vzy, vzz;
        coords::Cartesian_Point vir1, vir2, vir3, vir4;
        for (auto & improper : refined.impropers())
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
          coords::float_type const rt2(dot(t, t)), ru2(dot(u, u)), rtru = sqrt(rt2*ru2);
          if (abs(rtru) > 0.0)
          {
            coords::float_type const rcb(len(cb));
            coords::float_type const cosine(min(1.0, max(-1.0, (dot(t, u)/rtru))));
            coords::float_type const sine(dot(cb, tu)/(rcb*rtru));
            coords::float_type const angle(sine < 0.0 ? 
              -acos(cosine)*SCON_180PI : acos(cosine)*SCON_180PI);
            coords::float_type da((abs(angle+improper.p.ideal[0U]) < abs(angle-improper.p.ideal[0U])) ? 
              angle+improper.p.ideal[0U] : angle-improper.p.ideal[0U]);
            if (da > 180.0) da -= 360.0;
            if (da < -180.0) da += 360.0;
            da *= SCON_PI180;
            coords::float_type dE = improper.p.force[0]*da;
            E += dE*da;
            dE *= 2.0;

            coords::Cartesian_Point const ca =
              coords->xyz(improper.ligand[1]) - coords->xyz(improper.center);
            coords::Cartesian_Point const db =
              coords->xyz(improper.twist) - coords->xyz(improper.ligand[0U]);
            coords::Cartesian_Point const dt(cross(t, cb)*(dE/(rt2*rcb))), 
              du(cross(u, cb)*(-dE/(ru2*rcb)));

            vir1 = cross(dt, cb);
            vir2 = cross(ca, dt) + cross(du, dc);
            vir3 = cross(dt, ba) + cross(db, du);
            vir4 = cross(du, cb);

            part_grad[IMPROPER][improper.center]  += vir1;
            part_grad[IMPROPER][improper.ligand[0]] += vir2;
            part_grad[IMPROPER][improper.ligand[1]] += vir3;
            part_grad[IMPROPER][improper.twist] += vir4;

            vxx = cb.x()*(vir3.x() + vir4.x()) - ba.x()*vir1.x() + dc.x()*vir4.x();
            vyx = cb.y()*(vir3.x() + vir4.x()) - ba.y()*vir1.x() + dc.y()*vir4.x();
            vzx = cb.z()*(vir3.x() + vir4.x()) - ba.z()*vir1.x() + dc.z()*vir4.x();
            vyy = cb.y()*(vir3.y() + vir4.y()) - ba.y()*vir1.y() + dc.y()*vir4.y();
            vzy = cb.z()*(vir3.y() + vir4.y()) - ba.z()*vir1.y() + dc.z()*vir4.y();
            vzz = cb.z()*(vir3.z() + vir4.z()) - ba.z()*vir1.z() + dc.z()*vir4.z();

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
        }
        return E;
      }


	  /********************************
	  *                                *
	  *   Implicit solvation Energy    *
	  *                                *
	  *                                *
	  *                                *
	  *                                *
	  *********************************/

	  template<>
	  coords::float_type energy::interfaces::aco::aco_ff::solv<0>()
	  {
		  //GB::born solvate(0);

		  //solvate.get_E();

		  //std::cout << '\n' << '\n' << "           Implicit solvation interactions: " << solvate.interactions << '\n' << '\n';

		  //return solvate.ES;
      return coords::float_type(0);
	  }

	  template<>
	  coords::float_type energy::interfaces::aco::aco_ff::solv<1>()
	  {
		  //GB::born solvate(1);

		  //solvate.get_E();

		  //if (Config::get().general.verbosity == 50) std::cout << "Solvation energy :\t" << solvate.ES << "\t";

		  //part_grad[SOLVATE] = solvate.dAES;
		  //part_virial[SOLVATE] = solvate.vir;

		  //return solvate.ES;
      return coords::float_type(0);
	  }


      /********************************
      *                                *
      *  Improper                      *
      *  Dihedral Potential            *
      *  Energy/Gradients/Hessians     *
      *                                *
      *                                *
      *********************************/


      energy::interfaces::aco::nb_cutoff::nb_cutoff (
        coords::float_type const ic, coords::float_type const is)
        : c(ic), s(is), cc(c*c), ss(3.0*s*s), 
        cs((cc-s*s)*(cc-s*s)*(cc-s*s)) 
      { }


      bool energy::interfaces::aco::nb_cutoff::factors (coords::float_type const rr, 
        coords::float_type & r, coords::float_type & fQ, coords::float_type & fV)
      {
        using std::sqrt;
        using std::abs;
        r = sqrt(rr);
        if (r > c) return false;  // if distance bigger than cutoff -> return false
        coords::float_type const cr(cc - rr);
        fV = r < s ? 1.0 : (cr*cr*(cc+2.0*rr-ss)) / cs;
        fQ = (1.0 - rr/cc);
        fQ *= fQ;
        return (abs(r) > 0.0);  // return true (always???)
      }


      void energy::interfaces::aco::aco_ff::boundary(coords::float_type & x, coords::float_type & y, coords::float_type & z) const
      {
        static coords::Cartesian_Point const halfbox(Config::get().energy.pb_box/2.0);
        if (x > halfbox.x())
        {
          x -= Config::get().energy.pb_box.x();
        }
        else if (x < -halfbox.x())
        {
          x += Config::get().energy.pb_box.x();
        }
		if (y > halfbox.y())
        {
          y -= Config::get().energy.pb_box.y();
        }
		else if (y < -halfbox.y())
        {
          y += Config::get().energy.pb_box.y();
        }
		if (z > halfbox.z())
        {
          z -= Config::get().energy.pb_box.z();
        }
		else if (z < -halfbox.z())
        {
          z += Config::get().energy.pb_box.z();
        }
      }


      // C = q1*q2, ri = 1/r, dQ = dQ/dr
      inline coords::float_type energy::interfaces::aco::aco_ff::eQ 
        (coords::float_type const C, coords::float_type const ri) const
      {
        return C*ri;
      }

        
      // C = q1*q2, ri = 1/r, dQ = dQ/dr
      inline coords::float_type energy::interfaces::aco::aco_ff::gQ 
        (coords::float_type const C, coords::float_type const ri, coords::float_type & dQ) const
      {
        coords::float_type const Q = C*ri; // Q = C/r
        dQ = -Q*ri; // dQ/dr = -C/r^2 (derivative)
        return Q;
      }


      inline coords::float_type energy::interfaces::aco::aco_ff::gQ_fep 
        (coords::float_type const C, coords::float_type const ri, 
        coords::float_type const c_out, coords::float_type & dQ) const
      {
        using std::pow;
        coords::float_type const rmod = (1.0 - c_out) * 
          Config::get().fep.cshift + ri*ri*ri*ri*ri*ri; //r^6 shifted
        coords::float_type const Q = c_out * C / 
          pow(rmod, 0.16666666666666); //Q
        dQ = - c_out * C *  pow(ri, 5.0) / pow(-Config::get().fep.cshift * 
          c_out + Config::get().fep.cshift + std::pow(ri,6.0), 1.16666666666666);  // dQ/dr
        dQ = dQ/pow(rmod, 0.16666666666666); // dQ/dr shifted
        return Q;
      }


      /*
        lennard-jones potentials
      */

      
      template<> inline coords::float_type energy::interfaces::aco::aco_ff::eV
        < ::tinker::parameter::radius_types::R_MIN> 
        (coords::float_type const E, coords::float_type const R, coords::float_type const r) const
      {
        coords::float_type T = R*r;
        T = T*T*T; // T^3
        T = T*T; // T^6
        return E*T*(T-2.0);
      }


      template<> inline coords::float_type energy::interfaces::aco::aco_ff::eV
        < ::tinker::parameter::radius_types::SIGMA> 
        (coords::float_type const E, coords::float_type const R, coords::float_type const r) const
      {
        coords::float_type T = R*r;
        T = T*T*T; // T^3
        T = T*T; // T^6
        return E*T*(T-1.0);
      }

      /**calculate lenard-jones potential and gradient for charmm and amber forcefield (r_min-type);
      returns the energy
      @param E: 4 * epsilon-parameter
      @param R: r_min-parameter
      @param r: inverse distance 1/r between the two atoms
      @param dV: reference to variable that saves gradient*/
      template<> inline coords::float_type energy::interfaces::aco::aco_ff::gV
        < ::tinker::parameter::radius_types::R_MIN> 
        (coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type &dV) const
      {
        coords::float_type T = R*r;
        T = T*T*T; // T^3
        T = T*T; // T^6
        coords::float_type const V = E*T;
        dV = 12.0*V*r*(1.0-T);
        return V*(T-2.0);
      }

      /**calculate lenard-jones potential and gradient for oplsaa-forcefield (sigma-type); 
      returns the energy
      @param E: epsilon-parameter
      @param R: sigma-parameter
      @param r: inverse distance 1/r between the two atoms
      @param dV: reference to variable that saves gradient*/
      template<> inline coords::float_type energy::interfaces::aco::aco_ff::gV
        < ::tinker::parameter::radius_types::SIGMA> 
        (coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type &dV) const
      {
        coords::float_type T = R*r;
        T = T*T*T; // T^3
        T = T*T; // T^6
        coords::float_type const V = E*T;
        dV = V*r*(6.0-12.0*T);  // derivative
        return V*(T-1.0);  //potential
      }


      template<> inline coords::float_type energy::interfaces::aco::aco_ff::gV_fep
        < ::tinker::parameter::radius_types::R_MIN> 
        (coords::float_type const E, coords::float_type const R, coords::float_type const r, 
          coords::float_type const vout, coords::float_type &dV) const
      {
        coords::float_type D6, D12, D13, T2;
        coords::float_type T = R, D = r*r;
        T = T*T*T; // T^3
        T = T*T; // T^6
        T2 = T*T; // T^12
        D6 = D*D*D; //r^6 / D^3
        D13 = D6*D6*r; //r^13
        D6 = Config::get().fep.ljshift * (1 - vout) * (1 - vout) * T + D6; //r^6 shifted
        D12 = D6 * D6; //r^12 shifted
        coords::float_type V = vout * E * (T2/D12 - 2*T/D6);
        dV = vout * E * 12.0 * (T * D6 - (T2))/D13;
        dV = dV/std::pow(D6, 0.16666666666666);
        return V;
      }
      // AMBER/CHARMM FEP vdW part
      template<> inline coords::float_type energy::interfaces::aco::aco_ff::gV_fep_cut
        < ::tinker::parameter::radius_types::R_MIN>
        (coords::float_type const E, coords::float_type const R, 
          coords::float_type const r, coords::float_type const vout, 
          coords::float_type const vout2, coords::float_type &dV, 
          coords::float_type &alche2, coords::float_type &fV, coords::float_type &fV2) const
      {
          coords::float_type A, B;
          coords::float_type K6, D6, T2;
          coords::float_type T = R, D = r*r, K = r*r;
          T = T*T*T; // T^3
          T = T*T; // T^6
          T2 = T*T; // T^12
          A = E*T2;
          B = E*T * 2;
          D += Config::get().fep.ljshift * (1 - vout); // r^2 shifted
          K += Config::get().fep.ljshift * (1 - vout2);
          D6 = D*D*D; // r^6 shifted
          K6 = K*K*K;
          coords::float_type V = A / (D6*D6) - B / D6;
          alche2 = (A / (K6*K6) - B / K6) * vout2 * fV;
          dV = -vout * ((12 * V + 6 * B / D6) / D * fV + V * fV2);
          V *= vout*fV;
          return V;
        }
      // OPLS-AA FEP vdW part
      template<> inline coords::float_type energy::interfaces::aco::aco_ff::gV_fep_cut
        < ::tinker::parameter::radius_types::SIGMA>
        (coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type const vout, 
          coords::float_type const vout2, coords::float_type &dV, coords::float_type &alche2, coords::float_type &fV, coords::float_type &fV2) const
      {
          coords::float_type K6, D6, T2, A, B;
          coords::float_type T = R, D = r*r, K = r*r;
          T = T*T*T; // T^3
          T = T*T; // T^6
          D += Config::get().fep.ljshift * (1 - vout); // r^2 shifted
          K += Config::get().fep.ljshift * (1 - vout2); //r^2 Dl shifted
          D6 = D*D*D; // r^6 shifted
          K6 = K*K*K;
          A = E*T2;
          B = E*T;
          coords::float_type V = A / (D6*D6) - B / D6;
          alche2 = (A / (K6*K6) - B / K6) * vout2 * fV;
          dV = -vout * ((12 * V + 6 * B / D6) / D * fV + V * fV2);
          V *= vout*fV;
          return V;
        }

        template<> inline coords::float_type energy::interfaces::aco::aco_ff::gV_fep
          < ::tinker::parameter::radius_types::SIGMA> 
        (coords::float_type const E, coords::float_type const R, 
          coords::float_type const r, coords::float_type const vout, 
          coords::float_type &dV) const
      {
        coords::float_type T = R, D = r*r;
        T = T*T*T; // T^3
        T = T*T; // T^6
        D = D*D*D; // r^6
        D = Config::get().fep.ljshift * (1-vout) * (1-vout) * T + D; // r^6 shifted
        T /= D;
        coords::float_type V = vout*E*T;
        dV = V*r*(6.0-12.0*T);
        return V*(T-1.0);
      }


      template< ::tinker::parameter::radius_types::T RT> 
      inline void energy::interfaces::aco::aco_ff::e_QV  
        (coords::float_type const C, coords::float_type const E, 
          coords::float_type const R, coords::float_type const d, 
        coords::float_type &e_c, coords::float_type &e_v) const
      {
        e_c += eQ(C, d);
        e_v += eV<RT>(E, R, d);
      }


      template< ::tinker::parameter::radius_types::T RT> 
      inline void energy::interfaces::aco::aco_ff::g_QV  
        (coords::float_type const C, coords::float_type const E, 
          coords::float_type const R, coords::float_type const d, 
        coords::float_type &e_c, coords::float_type &e_v, coords::float_type &dE) const
      {
        coords::float_type dQ(0.0), dV(0.0);
        e_c += gQ(C, d, dQ);
        e_v += gV<RT>(E, R, d, dV);
        dE = (dQ + dV)*d;
      }


      template< ::tinker::parameter::radius_types::T RT> 
      inline void energy::interfaces::aco::aco_ff::g_QV_fep  
        (coords::float_type const C, coords::float_type const E, 
          coords::float_type const R, coords::float_type const d, 
                coords::float_type const c_io, coords::float_type const v_io,
                coords::float_type &e_c, coords::float_type &e_v, coords::float_type &dE) const
      {
        coords::float_type dQ(0.0), dV(0.0);
        e_c += gQ_fep(C, d, c_io, dQ);
        e_v += gV_fep<RT>(E, R, d, v_io, dV);
        dE = (dQ + dV);
      }


      template< ::tinker::parameter::radius_types::T RT> 
      inline void energy::interfaces::aco::aco_ff::e_QV_cutoff  
        (coords::float_type const C, coords::float_type const E, 
          coords::float_type const R, coords::float_type const d, 
        coords::float_type const fQ, coords::float_type const fV,
        coords::float_type &e_c, coords::float_type &e_v) const
      {
        e_c += eQ(C, d)*fQ;
        e_v += eV<RT>(E, R, d)*fV;
      }


      template< ::tinker::parameter::radius_types::T RT> 
      inline void energy::interfaces::aco::aco_ff::g_QV_cutoff 
        (coords::float_type const C, coords::float_type const E, 
          coords::float_type const R, coords::float_type const d, 
        coords::float_type const fQ, coords::float_type const fV,
        coords::float_type &e_c, coords::float_type &e_v, coords::float_type &dE) const
      {
        coords::float_type dQ(0.0), dV(0.0);
        e_c += gQ(C, d, dQ)*fQ;
        e_v += gV<RT>(E, R, d, dV)*fV; 
        dE = (dQ*fQ+dV*fV)*d;
      }


      template< ::tinker::parameter::radius_types::T RT> 
      inline void energy::interfaces::aco::aco_ff::g_QV_fep_cutoff 
        (coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const d, 
        coords::float_type const c_out, coords::float_type const v_out, coords::float_type const fQ, 
        coords::float_type const fV, coords::float_type &e_c, coords::float_type &e_v, coords::float_type &dE) const
      {
        coords::float_type dQ, dV;
        e_c += gQ_fep(C, d, c_out, dQ)*fQ;
        e_v += gV_fep<RT>(E, R, d, v_out, dV)*fV;
        dE = (dQ*fQ + dV*fV);
      }


      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_nb (void)
      {

        part_energy[types::CHARGE] = 0.0;
        part_energy[types::VDW] = 0.0;
        part_energy[types::VDWC] = 0.0;
        part_grad[types::VDWC].assign(part_grad[types::VDWC].size(), coords::Cartesian_Point());

        coords->fep.feptemp = energy::fepvect();
        for (auto & ia : coords->interactions()) ia.energy = 0.0;

        for (auto const &pairmatrix : refined.pair_matrices())
        {
          
          size_t const N(coords->interactions().size());
          for (size_t sub_ia_index(0u), row(0u), col(0u); sub_ia_index<N; ++sub_ia_index)
          {
            coords::float_type & e(coords->interactions(sub_ia_index).energy);
            coords::Representation_3D & g(coords->interactions(sub_ia_index).grad);
            g.assign(coords->size(), coords::Cartesian_Point());
            std::vector< ::tinker::refine::types::nbpair> const & pl(pairmatrix.pair_matrix(sub_ia_index));
            scon::matrix< ::tinker::parameter::combi::vdwc, true> const & par(refined.vdwcm(pairmatrix.param_matrix_id));
            if (Config::get().md.fep)
            {
              if (Config::get().energy.periodic)
              { 
                  if (coords->atoms().in_exists() && (coords->atoms().sub_in() == row || coords->atoms().sub_in() == col))
                    g_nb_QV_pairs_fep_io<RT, true, false>(e, g, pl, par);
                  else if (coords->atoms().out_exists() && (coords->atoms().sub_out() == row || coords->atoms().sub_out() == col)) 
                    g_nb_QV_pairs_fep_io<RT, true, true>(e, g, pl, par);
                  else
                    g_nb_QV_pairs_cutoff<RT, true>(e, g, pl, par);
              }
              else  // no periodic boundaries
              {
                  if (coords->atoms().in_exists() && (coords->atoms().sub_in() == row || coords->atoms().sub_in() == col))
                    g_nb_QV_pairs_fep_io<RT, false, false>(e, g, pl, par);
                  else if (coords->atoms().out_exists() && (coords->atoms().sub_out() == row || coords->atoms().sub_out() == col))
                    g_nb_QV_pairs_fep_io<RT, false, true>(e, g, pl, par);
                  else if (Config::get().energy.cutoff < 1000.0)
                    g_nb_QV_pairs_cutoff<RT, false>(e, g, pl, par);
                  else
                    g_nb_QV_pairs<RT>(e, g, pl, par);
              }
            }
            else   // no fep
            {
              if (Config::get().energy.periodic)
                g_nb_QV_pairs_cutoff<RT, true>(e, g, pl, par);
              else if (Config::get().energy.cutoff < 1000.0)
                g_nb_QV_pairs_cutoff<RT, false>(e, g, pl, par);
              else
                g_nb_QV_pairs<RT>(e, g, pl, par);
            }
            if (col==row)
            {
              col = 0;
              ++row;
            }
            else ++col;
            part_grad[types::VDWC] += g;
          }
        }
        if (Config::get().md.fep)
        {
            coords->fep.feptemp.dE = (coords->fep.feptemp.e_c_l2 + coords->fep.feptemp.e_vdw_l2) - (coords->fep.feptemp.e_c_l1 + coords->fep.feptemp.e_vdw_l1);
            coords->fep.feptemp.dG = 0;
            coords->fep.fepdata.push_back(coords->fep.feptemp);
        }
      }


#ifndef _OPENMP

      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs 
      (
        coords::float_type &e_nb, coords::Representation_3D &grad_vector, 
        std::vector< ::tinker::refine::types::nbpair> const & pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params
      )
      {
        int number = 0;
        coords::float_type e_c(0.0), e_v(0.0);
		    //std::vector <coords::float_type> charge_vector(coords->size());
        std::cout << scon::c3_delimeter(' ');
        for (auto const & pair : pairlist)
        {
          number += 1;
          coords::Cartesian_Point b(coords->xyz(pair.a) - coords->xyz(pair.b));
          coords::float_type r = 1.0/len(b), dE(0.0);
          ::tinker::parameter::combi::vdwc const & p(params(refined.type(pair.a), refined.type(pair.b)));
          g_QV<RT>(p.C, p.E, p.R, r, e_c, e_v, dE);
          b *= dE;
		      //charge_vector[pair.a] +=eQ(p.C,r)/r;
		      //charge_vector[pair.b] +=eQ(p.C,r)/r;
          grad_vector[pair.a] += b;
          grad_vector[pair.b] -= b;
        }
        e_nb += e_c+e_v;
        //part_energy[types::VDWC] += e_c+e_v;
        part_energy[types::CHARGE] += e_c;
        part_energy[types::VDW] += e_v;
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff 
      (
        coords::float_type &e_nb, coords::Representation_3D &grad_vector, 
        std::vector< ::tinker::refine::types::nbpair> const & pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params
      )
      {
        nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
        coords::float_type e_c(0.0), e_v(0.0);
        //size_t ia(0u);
        for (auto const & pair : pairlist)
        {
          coords::Cartesian_Point b(coords->xyz(pair.a) - coords->xyz(pair.b));
          coords::Cartesian_Point dist;
          if (PERIODIC) boundary(b.x(), b.y(), b.z());
          coords::float_type const rr = dot(b, b);
          coords::float_type r(0.0), fQ(0.0), fV(0.0), dE(0.0);
          if (!cutob.factors(rr, r, fQ, fV))
          {
            continue;
          }
          //++ia;
          r = 1.0/r;
          ::tinker::parameter::combi::vdwc const & p(params(refined.type(pair.a), refined.type(pair.b)));
          g_QV_cutoff<RT>(p.C, p.E, p.R, r, fQ, fV, e_c, e_v, dE);
          dist = b;
          b *= dE;
          grad_vector[pair.a] += b;
          grad_vector[pair.b] -= b;
          //Increment internal virial tensor
          coords::float_type const vxx = b.x() * dist.x();
          coords::float_type const vyx = b.x() * dist.y();
          coords::float_type const vzx = b.x() * dist.z();
          coords::float_type const vyy = b.y() * dist.y();
          coords::float_type const vzy = b.y() * dist.z();
          coords::float_type const vzz = b.z() * dist.z();
          part_virial[VDWC][0][0] += vxx;
          part_virial[VDWC][1][0] += vyx;
          part_virial[VDWC][2][0] += vzx;
          part_virial[VDWC][0][1] += vyx;
          part_virial[VDWC][1][1] += vyy;
          part_virial[VDWC][2][1] += vzy;
          part_virial[VDWC][0][2] += vzx;
          part_virial[VDWC][1][2] += vzy;
          part_virial[VDWC][2][2] += vzz;
        }
        //std::cout << "IAs of " << &pairlist << " : " << ia << "\n";
        e_nb += e_c+e_v;
        //part_energy[types::VDWC] += e_c+e_v;
        part_energy[types::CHARGE] += e_c;
        part_energy[types::VDW] += e_v;
      }


	  template< ::tinker::parameter::radius_types::T RT, bool PERIODIC, bool ALCH_OUT>
	  void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io
	  (
		  coords::float_type &e_nb, coords::Representation_3D &grad_vector,
		  std::vector< ::tinker::refine::types::nbpair> const & pairlist,
		  scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params
	  )
	  {
		  nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
		  coords::float_type e_c(0.0), e_v(0.0), e_c_l(0.0), e_vdw_l(0.0), e_c_dl(0.0), e_vdw_dl(0.0);
		  fepvar const & fep = coords->fep.window[coords->fep.window[0].step];
		  std::ptrdiff_t const M(pairlist.size());
		  {
			  coords::Representation_3D tmp_grad(grad_vector.size());
			  coords::virial_t tempvir(coords::empty_virial());
			  for (std::ptrdiff_t i = 0; i<M; ++i)
			  {
				  coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));
				  if (PERIODIC) boundary(b.x(), b.y(), b.z());
				  ::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));
				  coords::float_type rr(dot(b, b)), dE, Q(0.0), V(0.0);
				  coords::Cartesian_Point dist;
				  if (PERIODIC)
				  {
					  coords::float_type fQ(0.0), fV(0.0);
					  coords::float_type r(0.0);
					  if (!cutob.factors(rr, r, fQ, fV)) continue;
					  g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout),
						  (ALCH_OUT ? fep.vin : fep.vout), fQ, fV, Q, V, dE);
					  coords::float_type trash(0.0);
					  g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout),
						  (ALCH_OUT ? fep.dvin : fep.dvout), fQ, fV, e_c_dl, e_vdw_dl, trash);
				  }
				  else
				  {
					  if (Config::get().energy.cutoff < 1000.0)
					  {
						  coords::float_type r(0.0);
						  coords::float_type fQ(0.0), fV(0.0);
						  if (!cutob.factors(rr, r, fQ, fV)) continue;
						  g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout),
							  (ALCH_OUT ? fep.vin : fep.vout), fQ, fV, Q, V, dE);
						  coords::float_type trash(0.0);
						  g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout),
							  (ALCH_OUT ? fep.dvin : fep.dvout), fQ, fV, e_c_dl, e_vdw_dl, trash);
					  }
					  else
					  {
						  coords::float_type const r = sqrt(rr);
						  g_QV_fep<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout),
							  (ALCH_OUT ? fep.vin : fep.vout), Q, V, dE);
						  coords::float_type trash(0.0);
						  g_QV_fep<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout),
							  (ALCH_OUT ? fep.dvin : fep.dvout), e_c_dl, e_vdw_dl, trash);
					  }
				  }
				  dist = b;
				  b *= dE;
				  e_c_l += Q;
				  e_vdw_l += V;
				  e_c += Q;
				  e_v += V;
				  tmp_grad[pairlist[i].a] += b;
				  tmp_grad[pairlist[i].b] -= b;
				  //Increment internal virial tensor
				  coords::float_type const vxx = b.x() * dist.x();
				  coords::float_type const vyx = b.x() * dist.y();
				  coords::float_type const vzx = b.x() * dist.z();
				  coords::float_type const vyy = b.y() * dist.y();
				  coords::float_type const vzy = b.y() * dist.z();
				  coords::float_type const vzz = b.z() * dist.z();
				  tempvir[0][0] += vxx;
				  tempvir[1][0] += vyx;
				  tempvir[2][0] += vzx;
				  tempvir[0][1] += vyx;
				  tempvir[1][1] += vyy;
				  tempvir[2][1] += vzy;
				  tempvir[0][2] += vzx;
				  tempvir[1][2] += vzy;
				  tempvir[2][2] += vzz;
			  }
			  {
				  grad_vector += tmp_grad;
				  for (int i = 0; i <= 2; i++) {
					  for (int k = 0; k <= 2; k++) {
						  part_virial[VDWC][i][k] += tempvir[i][k];
					  }
				  }
			  }
		  }
		  e_nb += e_c + e_v;
		  coords->fep.feptemp.e_c_l1 += e_c_l;
		  coords->fep.feptemp.e_c_l2 += e_c_dl;
		  coords->fep.feptemp.e_vdw_l1 += e_vdw_l;
		  coords->fep.feptemp.e_vdw_l2 += e_vdw_dl;
		  part_energy[types::CHARGE] += e_c;
		  part_energy[types::VDW] += e_v;
	  }


	  // #############################################################################################
	  //      PME CHARGE ENERGY (Third implementation, FEP included into reciprocal energy part)     #
	  // #############################################################################################
	  //template< ::tinker::parameter::radius_types::T RT>
	  //void energy::interfaces::aco::aco_ff::g_nb_pme(void)
	  //{
	  //  part_energy[types::CHARGE] = 0.0;
	  //  part_energy[types::VDW] = 0.0;
	  //  part_energy[types::VDWC] = 0.0;
	  //  part_grad[types::VDWC].assign(coords->size(), coords::Cartesian_Point());
	  //  // temp grad vector for reciprocal gradients
	  //  coords::Representation_3D g(coords->pme.pmetemp.natoms);
	  //  g.assign(coords->size(), coords::Cartesian_Point());
	  //  double e_corr(0.0), e_reci(0.0), e_dir(0.0), e_dir_scale(0.0), e_vdw(0.0);
	  //  // ################################################################
	  //  //                      FEP PART                                  #
	  //  //#################################################################
	  //  if (Config::get().md.fep)
	  //  {
	  //    double fepinreci(0.0), fepdinreci(0.0), fepoutreci(0.0), fepdoutreci(0.0), fepall(0.0);
	  //    coords->fep.feptemp = energy::fepvect();
	  //    // Calculate FEP-pme self interaction correction
	  //    pme_correction_fep(e_corr);
	  //    // ###################### IN ATOMS #############
	  //    // Get Spline and Grid stuff done
	  //    bsplinefepin();
	  //    setchargestofepindgrid();
	  //    // Forward fourier transformation of DIN charge grid  
	  //    fftw_execute(coords->pme.pmetemp.forward);
	  //    pme_reci_novir(fepdinreci);
	  //    setchargestofepingrid();
	  //    // Forward fourier transformation of IN charge grid  
	  //    fftw_execute(coords->pme.pmetemp.forward);
	  //    pme_reci(fepinreci);
	  //    // Backward fourier transformation of IN charge grid to get gradients
	  //    fftw_execute(coords->pme.pmetemp.backward);
	  //    // get reciprocal space gradients
	  //    pme_reci_fepingrad(g);
	  //    // ################## OUT ATOMS ################
	  //    bsplinefepout();
	  //    setchargestofepoutdgrid();
	  //    // Forward fourier transformation of DOUT charge grid
	  //    fftw_execute(coords->pme.pmetemp.forward);
	  //    pme_reci_novir(fepdoutreci);
	  //    setchargestofepoutgrid();
	  //    // Forward fourier transformation of OUT charge grid
	  //    fftw_execute(coords->pme.pmetemp.forward);
	  //    pme_reci(fepoutreci);
	  //    // Backward fourier transformation of OUT charge grid to get gradients
	  //    fftw_execute(coords->pme.pmetemp.backward);
	  //    pme_reci_fepoutgrad(g);
	  //    //################## ALL ATOM CORRECTION ###################
	  //    bsplinefepall();
	  //    setchargestoallgrid();
	  //    // Forward fourier transformation of ALL charge grid
	  //    fftw_execute(coords->pme.pmetemp.forward);
	  //    pme_reciall(fepall);
	  //    // Backward fourier transformation of ALL charge grid to get gradients
	  //    fftw_execute(coords->pme.pmetemp.backward);
	  //    pme_reci_fepallgrad(g);
	  //    // Correct charged cell
	  //    //pme_fepchargedcell();
	  //    // Add FEP contributions
	  //    coords->fep.feptemp.e_c_l2 += fepdinreci + fepdoutreci - fepall;
	  //    coords->fep.feptemp.e_c_l1 += fepinreci + fepoutreci - fepall;
	  //    // Add complete reciprocal energy
	  //    e_reci = fepinreci + fepoutreci - fepall;
	  //    // calculate scaled FEP contributions
	  //    pme_direct_scaled_fep(e_dir_scale, g);
	  //    // Calculate short range part of the electrostatic energy using standard CAST loops
	  //    for (auto & ia : coords->interactions()) ia.energy = 0.0;
	  //    for (auto const &pairmatrix : refined.pair_matrices())
	  //    {
	  //      std::size_t const N(coords->interactions().size());
	  //      for (std::size_t sub_ia_index(0u), row(0u), col(0u); sub_ia_index < N; ++sub_ia_index)
	  //      {
	  //        std::vector< ::tinker::refine::types::nbpair> const & pl(pairmatrix.pair_matrix(sub_ia_index));
	  //        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & par(refined.vdwcm(pairmatrix.param_matrix_id));



	  //        if (coords->atoms().in_exists() && (coords->atoms().sub_in() == row || coords->atoms().sub_in() == col))
	  //          pme_direct_fep<RT, false>(e_dir, g, pl, par, e_vdw);
	  //        else if (coords->atoms().out_exists() && (coords->atoms().sub_out() == row || coords->atoms().sub_out() == col))
	  //          pme_direct_fep<RT, true>(e_dir, g, pl, par, e_vdw);
	  //        else pme_direct<RT>(e_dir, g, pl, par, e_vdw);
	  //        // next matrix slot
	  //        if (col == row)
	  //        {
	  //          col = 0;
	  //          ++row;
	  //        }
	  //        else ++col;
	  //      }
	  //    }
	  //    // sum up FEP contributions
	  //    coords->fep.feptemp.dE = (coords->fep.feptemp.e_c_l2 + coords->fep.feptemp.e_vdw_l2) - (coords->fep.feptemp.e_c_l1 + coords->fep.feptemp.e_vdw_l1);
	  //    coords->fep.feptemp.dG = 0;
	  //    coords->fep.fepdata.push_back(coords->fep.feptemp);
	  //  }
	  //  // ###########################################################
	  //  //                      NORMAL PART (NO FEP)                 #
	  //  //############################################################
	  //  else
	  //  {
	  //    // Calculate self correction term
	  //    pme_correction(e_corr);
	  //    // Get spline coefficients, fill table array and set charges on the pme grid 
	  //    bsplinegeneration();
	  //    setchargestogrid();
	  //    // Forward fourier transformation of charge grid (complex)    
	  //    fftw_execute(coords->pme.pmetemp.forward);
	  //    // get recirpocal space energy and virial components
	  //    pme_reci(e_reci);
	  //    // Backward fourier transformation of charge grid (complex)
	  //    fftw_execute(coords->pme.pmetemp.backward);
	  //    // get reciprocal space gradients
	  //    pme_reci_grad(g);
	  //    // calculate direct scaled terms (1-2 and 1-3)
	  //    pme_direct_scaled(e_dir_scale, g);
	  //    for (auto & ia : coords->interactions()) ia.energy = 0.0;
	  //    // calculate direct energy via normal loop
	  //    for (auto const &pairmatrix : refined.pair_matrices())
	  //    {
	  //      std::size_t const N(coords->interactions().size());
	  //      for (std::size_t sub_ia_index(0u), row(0u), col(0u); sub_ia_index < N; ++sub_ia_index)
	  //      {
	  //        std::vector< ::tinker::refine::types::nbpair> const & pl(pairmatrix.pair_matrix(sub_ia_index));
	  //        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & par(refined.vdwcm(pairmatrix.param_matrix_id));
	  //        pme_direct<RT>(e_dir, g, pl, par, e_vdw);
	  //        if (col == row)
	  //        {
	  //          col = 0;
	  //          ++row;
	  //        }
	  //        else ++col;
	  //      }
	  //    }// end of matrix loop
	  //  }
	  //  part_grad[types::VDWC] += g;
	  //  part_energy[types::CHARGE] = e_corr + e_reci + e_dir + e_dir_scale;
	  //  part_energy[types::VDW] = e_vdw;
	  //}


   // template < ::tinker::parameter::radius_types::T RT, bool ALCH_OUT>
   // void energy::interfaces::aco::aco_ff::pme_direct_fep(double &e_ch,
   //   coords::Representation_3D &grad_vector,
   //   std::vector< ::tinker::refine::types::nbpair> const & pairlist,
   //   scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params, double & e_vdw)
   // {

   //   coords::float_type e_c(0.0), e_v(0.0), e_c_l(0.0), e_vdw_l(0.0), e_c_dl(0.0), e_vdw_dl(0.0);
   //   fepvar const & fep = coords->fep.window[coords->fep.window[0].step];

   //     coords::Representation_3D tmp_grad(grad_vector.size());
   //     coords::virial_t tempvir(coords::empty_virial());
   //     ptrdiff_t const M(pairlist.size());
   //     double dist, d, r, achg(0.0), distbuf, rbuf, factor, ewald, error;
   //     double energy, energyc2, den, cutoff, scale, e_v_temp, trash;
   //     coords::Cartesian_Point vir;
   //     coords::float_type dV, alche2;
   //     coords::float_type fQ(0.0), fV(0.0), fV2, switch2(0.0), cutoff2(0.0), switchfactor(0.0);
   //     switch2 = Config::get().energy.switchdist *  Config::get().energy.switchdist;
   //     cutoff2 = Config::get().energy.cutoff *  Config::get().energy.cutoff;
   //     double temp = cutoff2*cutoff2 - switch2*switch2;
   //     switchfactor = 1 / (temp*temp*temp);

   //     for (std::ptrdiff_t i = 0; i < M; ++i)
   //     {
   //       coords::Cartesian_Point b(coords->xyz(pairlist[i].b) - coords->xyz(pairlist[i].a));
   //       boundary(b.x(), b.y(), b.z());
   //       dist = dot(b, b);
   //       if (dist >= cutoff2) continue;
   //       r = sqrt(dist);
   //       d = 1 / r;
   //       rbuf = r + achg;
   //       vir = b;
   //       distbuf = rbuf*rbuf;
   //       factor = 332.0716 * coords->pme.pmetemp.atomcharges[pairlist[i].a] * coords->pme.pmetemp.atomcharges[pairlist[i].b];
   //       ewald = coords->pme.pmetemp.ewaldcoeff * r;
   //       error = erfc(ewald);
   //       energy = (factor / rbuf) * error;
   //       // FEP charge energies
   //       e_c += energy * (ALCH_OUT ? fep.ein : fep.eout);
   //       // increment FEP charge energies
   //       e_c_l += energy * (ALCH_OUT ? fep.ein : fep.eout);
   //       e_c_dl += energy * (ALCH_OUT ? fep.dein : fep.deout);
   //       // van der waals part
   //       ::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));
   //       fV = (dist > switch2) ? ((cutoff2 - dist) * (cutoff2 - dist) * (cutoff2 - 3 * switch2 + 2 * dist)) * switchfactor : 1.0;
   //       fV2 = (dist > switch2) ? 12 * switchfactor * (cutoff2 - dist)*(dist - switch2) : 0.0;
   //       e_v_temp = gV_fep_cut<RT>(p.E, p.R, r, (ALCH_OUT ? fep.vin : fep.vout), (ALCH_OUT ? fep.dvin : fep.dvout), dV, alche2, fV, fV2);
   //       e_v += e_v_temp;
   //       e_vdw_l += e_v_temp;
   //       e_vdw_dl += alche2;
   //       // gradients 
   //       den = -factor * (error / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r) * (ALCH_OUT ? fep.ein : fep.eout);
   //       b *= (dV + den*d);
   //       tmp_grad[pairlist[i].a] -= b;
   //       tmp_grad[pairlist[i].b] += b;
   //       // virial components
   //       coords::float_type const vxx = b.x() * vir.x();
   //       coords::float_type const vyx = b.x() * vir.y();
   //       coords::float_type const vzx = b.x() * vir.z();
   //       coords::float_type const vyy = b.y() * vir.y();
   //       coords::float_type const vzy = b.y() * vir.z();
   //       coords::float_type const vzz = b.z() * vir.z();
   //       tempvir[0][0] += vxx;
   //       tempvir[1][0] += vyx;
   //       tempvir[2][0] += vzx;
   //       tempvir[0][1] += vyx;
   //       tempvir[1][1] += vyy;
   //       tempvir[2][1] += vzy;
   //       tempvir[0][2] += vzx;
   //       tempvir[1][2] += vzy;
   //       tempvir[2][2] += vzz;
   //     }

   //       grad_vector += tmp_grad;
   //       for (int i = 0; i <= 2; i++){
   //         for (int k = 0; k <= 2; k++){
   //           part_virial[VDWC][i][k] += tempvir[i][k];
   //         }
   //       }
   //   coords->fep.feptemp.e_c_l1 += e_c_l;
   //   coords->fep.feptemp.e_c_l2 += e_c_dl;
   //   coords->fep.feptemp.e_vdw_l1 += e_vdw_l;
   //   coords->fep.feptemp.e_vdw_l2 += e_vdw_dl;
   //   e_ch += e_c;
   //   e_vdw += e_v;
   // }



   // void energy::interfaces::aco::aco_ff::pme_direct_scaled_fep(double & e_nb,
   //   coords::Representation_3D &grad_vector)
   // {
   //   double dist, r, achg, distbuf, rbuf, factor, ewald, error;
   //   double energy, den, cutoff, cutoff2, scale;
   //   cutoff = Config::get().energy.cutoff;
   //   cutoff2 = cutoff*cutoff;
   //   achg = 0.0;
   //   coords::Cartesian_Point vir;
   //   coords::float_type e_c(0.0), e_v(0.0), e_c_l(0.0), e_vdw_l(0.0), e_c_dl(0.0), e_vdw_dl(0.0);
   //   fepvar const & fep = coords->fep.window[coords->fep.window[0].step];
   //   // 1-2 scaled interactions 
   //   for (auto const & bond : refined.bonds())
   //   {
   //     scale = -1.0;
   //     coords::Cartesian_Point b((coords->xyz(bond.atoms[0])) - (coords->xyz(bond.atoms[1])));
   //     boundary(b.x(), b.y(), b.z());
   //     dist = b.x()*b.x() + b.y()*b.y() + b.z()*b.z();
   //     vir = b;
   //     r = sqrt(dist);
   //     rbuf = r + achg;
   //     distbuf = rbuf*rbuf;
   //     factor = 332.0716 * coords->pme.pmetemp.atomcharges[bond.atoms[0]] * coords->pme.pmetemp.atomcharges[bond.atoms[1]];
   //     ewald = coords->pme.pmetemp.ewaldcoeff * r;
   //     error = erfc(ewald);
   //     energy = (factor / rbuf) * (error + scale);
   //     // atom is in
   //     if (coords->pme.pmetemp.feinf[bond.atoms[0]].flag == 1 || coords->pme.pmetemp.feinf[bond.atoms[1]].flag == 1)
   //     {
   //       e_c_l += energy * fep.eout;
   //       e_c_dl += energy * fep.deout;
   //       e_nb += energy * fep.eout;
   //       den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r) * fep.eout;
   //       den /= r;
   //       b *= den;
   //       grad_vector[bond.atoms[0]] += b;
   //       grad_vector[bond.atoms[1]] -= b;;
   //       coords::float_type const vxx = b.x() * vir.x();
   //       coords::float_type const vyx = b.x() * vir.y();
   //       coords::float_type const vzx = b.x() * vir.z();
   //       coords::float_type const vyy = b.y() * vir.y();
   //       coords::float_type const vzy = b.y() * vir.z();
   //       coords::float_type const vzz = b.z() * vir.z();
   //       part_virial[VDWC][0][0] += vxx;
   //       part_virial[VDWC][1][0] += vyx;
   //       part_virial[VDWC][2][0] += vzx;
   //       part_virial[VDWC][0][1] += vyx;
   //       part_virial[VDWC][1][1] += vyy;
   //       part_virial[VDWC][2][1] += vzy;
   //       part_virial[VDWC][0][2] += vzx;
   //       part_virial[VDWC][1][2] += vzy;
   //       part_virial[VDWC][2][2] += vzz;
   //     }
   //     else if (coords->pme.pmetemp.feinf[bond.atoms[0]].flag == 0 || coords->pme.pmetemp.feinf[bond.atoms[1]].flag == 0)
   //     {
   //       e_c_l += energy * fep.ein;
   //       e_c_dl += energy * fep.dein;
   //       e_nb += energy * fep.ein;
   //       den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * 
   //         exp(-(ewald*ewald)) / r) * fep.ein;
   //       den /= r;
   //       b *= den;
   //       grad_vector[bond.atoms[0]] += b;
   //       grad_vector[bond.atoms[1]] -= b;;
   //       coords::float_type const vxx = b.x() * vir.x();
   //       coords::float_type const vyx = b.x() * vir.y();
   //       coords::float_type const vzx = b.x() * vir.z();
   //       coords::float_type const vyy = b.y() * vir.y();
   //       coords::float_type const vzy = b.y() * vir.z();
   //       coords::float_type const vzz = b.z() * vir.z();
   //       part_virial[VDWC][0][0] += vxx;
   //       part_virial[VDWC][1][0] += vyx;
   //       part_virial[VDWC][2][0] += vzx;
   //       part_virial[VDWC][0][1] += vyx;
   //       part_virial[VDWC][1][1] += vyy;
   //       part_virial[VDWC][2][1] += vzy;
   //       part_virial[VDWC][0][2] += vzx;
   //       part_virial[VDWC][1][2] += vzy;
   //       part_virial[VDWC][2][2] += vzz;
   //     }
   //     else
   //     {
   //       e_c_l += energy;
   //       e_c_dl += energy;
   //       e_nb += energy;
   //       den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * 
   //         exp(-(ewald*ewald)) / r);
   //       den /= r;
   //       b *= den;
   //       grad_vector[bond.atoms[0]] += b;
   //       grad_vector[bond.atoms[1]] -= b;;
   //       coords::float_type const vxx = b.x() * vir.x();
   //       coords::float_type const vyx = b.x() * vir.y();
   //       coords::float_type const vzx = b.x() * vir.z();
   //       coords::float_type const vyy = b.y() * vir.y();
   //       coords::float_type const vzy = b.y() * vir.z();
   //       coords::float_type const vzz = b.z() * vir.z();
   //       part_virial[VDWC][0][0] += vxx;
   //       part_virial[VDWC][1][0] += vyx;
   //       part_virial[VDWC][2][0] += vzx;
   //       part_virial[VDWC][0][1] += vyx;
   //       part_virial[VDWC][1][1] += vyy;
   //       part_virial[VDWC][2][1] += vzy;
   //       part_virial[VDWC][0][2] += vzx;
   //       part_virial[VDWC][1][2] += vzy;
   //       part_virial[VDWC][2][2] += vzz;
   //     }
   //   }
   //   // 1-3 scaled interactions
   //   for (auto const & angle : refined.angles())
   //   {
   //     scale = -1.0;
   //     coords::Cartesian_Point b((coords->xyz(angle.atoms[0])) - (coords->xyz(angle.atoms[2])));
   //     boundary(b.x(), b.y(), b.z());
   //     dist = b.x()*b.x() + b.y()*b.y() + b.z()*b.z();
   //     r = sqrt(dist);
   //     vir = b;
   //     rbuf = r + achg;
   //     distbuf = rbuf*rbuf;
   //     factor = 332.0716 * coords->pme.pmetemp.atomcharges[angle.atoms[2]] * coords->pme.pmetemp.atomcharges[angle.atoms[0]];
   //     ewald = coords->pme.pmetemp.ewaldcoeff * r;
   //     error = erfc(ewald);
   //     energy = (factor / rbuf) * (error + scale);
   //     if (coords->pme.pmetemp.feinf[angle.atoms[0]].flag == 1 || coords->pme.pmetemp.feinf[angle.atoms[1]].flag == 1 || coords->pme.pmetemp.feinf[angle.atoms[2]].flag == 1)
   //     {
   //       e_c_l += energy * fep.eout;
   //       e_c_dl += energy * fep.deout;
   //       e_nb += energy * fep.eout;
   //       den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r) * fep.eout;
   //       den /= r;
   //       b *= den;
   //       grad_vector[angle.atoms[0]] += b;
   //       grad_vector[angle.atoms[2]] -= b;;
   //       coords::float_type const vxx = b.x() * vir.x();
   //       coords::float_type const vyx = b.x() * vir.y();
   //       coords::float_type const vzx = b.x() * vir.z();
   //       coords::float_type const vyy = b.y() * vir.y();
   //       coords::float_type const vzy = b.y() * vir.z();
   //       coords::float_type const vzz = b.z() * vir.z();
   //       part_virial[VDWC][0][0] += vxx;
   //       part_virial[VDWC][1][0] += vyx;
   //       part_virial[VDWC][2][0] += vzx;
   //       part_virial[VDWC][0][1] += vyx;
   //       part_virial[VDWC][1][1] += vyy;
   //       part_virial[VDWC][2][1] += vzy;
   //       part_virial[VDWC][0][2] += vzx;
   //       part_virial[VDWC][1][2] += vzy;
   //       part_virial[VDWC][2][2] += vzz;
   //     }
   //     else if (coords->pme.pmetemp.feinf[angle.atoms[0]].flag == 0 || coords->pme.pmetemp.feinf[angle.atoms[1]].flag == 0 || coords->pme.pmetemp.feinf[angle.atoms[2]].flag == 0)
   //     {
   //       e_c_l += energy * fep.ein;
   //       e_c_dl += energy * fep.dein;
   //       e_nb += energy * fep.ein;
   //       den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r) * fep.ein;
   //       den /= r;
   //       b *= den;
   //       grad_vector[angle.atoms[0]] += b;
   //       grad_vector[angle.atoms[2]] -= b;;
   //       coords::float_type const vxx = b.x() * vir.x();
   //       coords::float_type const vyx = b.x() * vir.y();
   //       coords::float_type const vzx = b.x() * vir.z();
   //       coords::float_type const vyy = b.y() * vir.y();
   //       coords::float_type const vzy = b.y() * vir.z();
   //       coords::float_type const vzz = b.z() * vir.z();
   //       part_virial[VDWC][0][0] += vxx;
   //       part_virial[VDWC][1][0] += vyx;
   //       part_virial[VDWC][2][0] += vzx;
   //       part_virial[VDWC][0][1] += vyx;
   //       part_virial[VDWC][1][1] += vyy;
   //       part_virial[VDWC][2][1] += vzy;
   //       part_virial[VDWC][0][2] += vzx;
   //       part_virial[VDWC][1][2] += vzy;
   //       part_virial[VDWC][2][2] += vzz;
   //     }
   //     else
   //     {
   //       {
   //         e_c_l += energy;
   //         e_c_dl += energy;
   //         e_nb += energy;
   //         den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r);
   //         den /= r;
   //         b *= den;
   //         grad_vector[angle.atoms[0]] += b;
   //         grad_vector[angle.atoms[2]] -= b;;
   //         coords::float_type const vxx = b.x() * vir.x();
   //         coords::float_type const vyx = b.x() * vir.y();
   //         coords::float_type const vzx = b.x() * vir.z();
   //         coords::float_type const vyy = b.y() * vir.y();
   //         coords::float_type const vzy = b.y() * vir.z();
   //         coords::float_type const vzz = b.z() * vir.z();
   //         part_virial[VDWC][0][0] += vxx;
   //         part_virial[VDWC][1][0] += vyx;
   //         part_virial[VDWC][2][0] += vzx;
   //         part_virial[VDWC][0][1] += vyx;
   //         part_virial[VDWC][1][1] += vyy;
   //         part_virial[VDWC][2][1] += vzy;
   //         part_virial[VDWC][0][2] += vzx;
   //         part_virial[VDWC][1][2] += vzy;
   //         part_virial[VDWC][2][2] += vzz;
   //       }
   //     }
   //   }
   //   coords->fep.feptemp.e_c_l1 += e_c_l;
   //   coords->fep.feptemp.e_c_l2 += e_c_dl;

   // }

   // template < ::tinker::parameter::radius_types::T RT>
   // void energy::interfaces::aco::aco_ff::pme_direct(double &e_ch,
   //   coords::Representation_3D &grad_vector,
   //   std::vector< ::tinker::refine::types::nbpair> const & pairlist,
   //   scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params, double & e_vdw)
   // {

   //   coords::float_type e_c(0.0), e_v(0.0);

   //     nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
   //     coords::Representation_3D tmp_grad(grad_vector.size());
   //     coords::virial_t tempvir(coords::empty_virial());
   //     ptrdiff_t const M(pairlist.size());
   //     double dist, d, r, achg(0.0), distbuf, rbuf, factor, ewald, error;
   //     double energy, den, cutoff, cutoff2, scale;
   //     cutoff = Config::get().energy.cutoff;
   //     cutoff2 = cutoff*cutoff;
   //     coords::Cartesian_Point vir;
   //     coords::float_type dV;
   //     coords::float_type fQ(0.0), fV(0.0);

   //     for (ptrdiff_t i = 0; i < M; ++i)
   //     {
   //       coords::Cartesian_Point b(coords->xyz(pairlist[i].b) - coords->xyz(pairlist[i].a));
   //       boundary(b.x(), b.y(), b.z());
   //       dist = dot(b, b);
   //       if (dist >= cutoff2) continue;
   //       r = sqrt(dist);
   //       d = 1 / r;
   //       rbuf = r + achg;
   //       vir = b;
   //       distbuf = rbuf*rbuf;
   //       factor = 332.0716 * coords->pme.pmetemp.atomcharges[pairlist[i].a] * coords->pme.pmetemp.atomcharges[pairlist[i].b];
   //       ewald = coords->pme.pmetemp.ewaldcoeff * r;
   //       error = erfc(ewald);
   //       energy = (factor / rbuf) * error;
   //       e_c += energy;
   //       // van der waals part
   //       ::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));
   //       if (!cutob.factors(dist, r, fQ, fV)) continue;
   //       e_v += gV<RT>(p.E, p.R, d, dV)*fV;
   //       // gradients 
   //       den = -factor * (error / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r);
   //       b *= (dV*fV + den) * d;
   //       tmp_grad[pairlist[i].a] -= b;
   //       tmp_grad[pairlist[i].b] += b;
   //       // virial components
   //       coords::float_type const vxx = b.x() * vir.x();
   //       coords::float_type const vyx = b.x() * vir.y();
   //       coords::float_type const vzx = b.x() * vir.z();
   //       coords::float_type const vyy = b.y() * vir.y();
   //       coords::float_type const vzy = b.y() * vir.z();
   //       coords::float_type const vzz = b.z() * vir.z();
   //       tempvir[0][0] += vxx;
   //       tempvir[1][0] += vyx;
   //       tempvir[2][0] += vzx;
   //       tempvir[0][1] += vyx;
   //       tempvir[1][1] += vyy;
   //       tempvir[2][1] += vzy;
   //       tempvir[0][2] += vzx;
   //       tempvir[1][2] += vzy;
   //       tempvir[2][2] += vzz;
   //     }

   //       grad_vector += tmp_grad;
   //       for (int i = 0; i <= 2; i++){
   //         for (int k = 0; k <= 2; k++){
   //           part_virial[VDWC][i][k] += tempvir[i][k];
   //         }
   //       }

   //   coords->fep.feptemp.e_c_l1 += e_c;
   //   coords->fep.feptemp.e_c_l2 += e_c;
   //   coords->fep.feptemp.e_vdw_l1 += e_v;
   //   coords->fep.feptemp.e_vdw_l2 += e_v;
   //   e_ch += e_c;
   //   e_vdw += e_v;
   // }

#else

//// #############################################################################################
////      PME CHARGE ENERGY (Third implementation, FEP included into reciprocal energy part)     #
//// #############################################################################################
//template< ::tinker::parameter::radius_types::T RT>
//void energy::interfaces::aco::aco_ff::g_nb_pme(void)
//{
//  part_energy[types::CHARGE] = 0.0;
//  part_energy[types::VDW] = 0.0;
//  part_energy[types::VDWC] = 0.0;
//  part_grad[types::VDWC].assign(coords->size(), coords::Cartesian_Point());
//  // temp grad vector for reciprocal gradients
//  coords::Representation_3D g(coords->pme.pmetemp.natoms);
//  g.assign(coords->size(), coords::Cartesian_Point());
//  double e_corr(0.0), e_reci(0.0), e_dir(0.0), e_dir_scale(0.0), e_vdw(0.0);
//  // ################################################################
//  //                      FEP PART                                  #
//  //#################################################################
//  if (Config::get().md.fep)
//  {
//    double fepinreci(0.0), fepdinreci(0.0), fepoutreci(0.0), fepdoutreci(0.0), fepall(0.0);
//    coords->fep.feptemp = energy::fepvect();
//    // Calculate FEP-pme self interaction correction
//    pme_correction_fep(e_corr);
//    // ###################### IN ATOMS #############
//    // Get Spline and Grid stuff done
//    bsplinefepin();
//    generateintable();
//    setchargestofepindgrid();
//    // Forward fourier transformation of DIN charge grid  
//    fftw_execute(coords->pme.pmetemp.forward);
//    pme_reci_novir(fepdinreci);
//    setchargestofepingrid();
//    // Forward fourier transformation of IN charge grid  
//    fftw_execute(coords->pme.pmetemp.forward);
//    pme_reci(fepinreci);
//    // Backward fourier transformation of IN charge grid to get gradients
//    fftw_execute(coords->pme.pmetemp.backward);
//    // get reciprocal space gradients
//    pme_reci_fepingrad(g);
//    // ################## OUT ATOMS ################
//    bsplinefepout();
//    generateouttable();
//    setchargestofepoutdgrid();
//    // Forward fourier transformation of DOUT charge grid
//    fftw_execute(coords->pme.pmetemp.forward);
//    pme_reci_novir(fepdoutreci);
//    setchargestofepoutgrid();
//    // Forward fourier transformation of OUT charge grid
//    fftw_execute(coords->pme.pmetemp.forward);
//    pme_reci(fepoutreci);
//    // Backward fourier transformation of OUT charge grid to get gradients
//    fftw_execute(coords->pme.pmetemp.backward);
//    pme_reci_fepoutgrad(g);
//    //################## ALL ATOM CORRECTION ###################
//    bsplinefepall();
//    generatealltable();
//    setchargestoallgrid();
//    // Forward fourier transformation of ALL charge grid
//    fftw_execute(coords->pme.pmetemp.forward);
//    pme_reciall(fepall);
//    // Backward fourier transformation of ALL charge grid to get gradients
//    fftw_execute(coords->pme.pmetemp.backward);
//    pme_reci_fepallgrad(g);
//    // Add FEP contributions
//    coords->fep.feptemp.e_c_l2 += fepdinreci + fepdoutreci - fepall;
//    coords->fep.feptemp.e_c_l1 += fepinreci + fepoutreci - fepall;
//    // Add complete reciprocal energy
//    e_reci = fepinreci + fepoutreci - fepall;
//    // calculate scaled FEP contributions
//    pme_direct_scaled_fep(e_dir_scale, g);
//    // Calculate short range part of the electrostatic energy using standard CAST loops
//    for (auto & ia : coords->interactions()) ia.energy = 0.0;
//    for (auto const &pairmatrix : refined.pair_matrices())
//    {
//      std::size_t const N(coords->interactions().size());
//      for (std::size_t sub_ia_index(0u), row(0u), col(0u); sub_ia_index < N; ++sub_ia_index)
//      {
//        std::vector< ::tinker::refine::types::nbpair> const & pl(pairmatrix.pair_matrix(sub_ia_index));
//        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & par(refined.vdwcm(pairmatrix.param_matrix_id));
//
//        if (coords->atoms().in_exists() && (coords->atoms().sub_in() == row || coords->atoms().sub_in() == col))
//          pme_direct_fep<RT, false>(e_dir, g, pl, par, e_vdw);
//        else if (coords->atoms().out_exists() && (coords->atoms().sub_out() == row || coords->atoms().sub_out() == col))
//          pme_direct_fep<RT, true>(e_dir, g, pl, par, e_vdw);
//        else pme_direct<RT>(e_dir, g, pl, par, e_vdw);
//        // next matrix slot
//        if (col == row)
//        {
//          col = 0;
//          ++row;
//        }
//        else ++col;
//      }
//    }
//    // sum up FEP contributions
//    coords->fep.feptemp.dE = (coords->fep.feptemp.e_c_l2 + coords->fep.feptemp.e_vdw_l2) - (coords->fep.feptemp.e_c_l1 + coords->fep.feptemp.e_vdw_l1);
//    coords->fep.feptemp.dG = 0;
//    coords->fep.fepdata.push_back(coords->fep.feptemp);
//  }
//  // ###########################################################
//  //                      NORMAL PART (NO FEP)                 #
//  //############################################################
//  else
//  {
//    // Calculate self correction term
//    pme_correction(e_corr);
//    // Get spline coefficients, fill table array and set charges on the pme grid 
//    bsplinegeneration();
//    generatetable();
//    setchargestogrid();
//    // Forward fourier transformation of charge grid (complex)    
//    fftw_execute(coords->pme.pmetemp.forward);
//    // get recirpocal space energy and virial components
//    pme_reci(e_reci);
//    // Backward fourier transformation of charge grid (complex)
//    fftw_execute(coords->pme.pmetemp.backward);
//    // get reciprocal space gradients
//    pme_reci_grad(g);
//    // calculate direct scaled terms (1-2 and 1-3)
//    pme_direct_scaled(e_dir_scale, g);
//    // Standard loops for direct space term
//    for (auto & ia : coords->interactions()) ia.energy = 0.0;
//    // calculate direct energy via normal loop
//    for (auto const &pairmatrix : refined.pair_matrices())
//    {
//      std::size_t const N(coords->interactions().size());
//      for (std::size_t sub_ia_index(0u), row(0u), col(0u); sub_ia_index < N; ++sub_ia_index)
//      {
//        std::vector< ::tinker::refine::types::nbpair> const & pl(pairmatrix.pair_matrix(sub_ia_index));
//        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & par(refined.vdwcm(pairmatrix.param_matrix_id));
//        pme_direct<RT>(e_dir, g, pl, par, e_vdw);
//        if (col == row)
//        {
//          col = 0;
//          ++row;
//        }
//        else ++col;
//      }
//    }// end of matrix loop
//  }
//  part_grad[types::VDWC] += g;
//  part_energy[types::CHARGE] = e_corr + e_reci + e_dir + e_dir_scale;
//  part_energy[types::VDW] = e_vdw;
//}

    
      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs 
      (
        coords::float_type &e_nb, coords::Representation_3D &grad_vector, 
        std::vector< ::tinker::refine::types::nbpair> const & pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params
      )
      {
        std::ptrdiff_t const M(pairlist.size());
        coords::float_type e_c(0.0), e_v(0.0);
        #pragma omp parallel
        {
          coords::Representation_3D tmp_grad(grad_vector.size());
          #pragma omp for reduction (+: e_c, e_v)
          for (std::ptrdiff_t i=0; i<M ; ++i)       // for every pair in pairlist
          {
            double current_c;   // Q_a * Q_b from AMBER
            if (Config::get().general.input == config::input_types::AMBER)
            {    // calculate Q_a * Q_b from AMBER charges (better if this would be done while building up pairlist)
              double ca = Config::get().coords.charges[pairlist[i].a];
              double cb = Config::get().coords.charges[pairlist[i].b];
              current_c = ca * cb;
            }
            coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));
            coords::float_type const r = 1.0 / std::sqrt(dot(b, b));
            coords::float_type dE(0.0);
            ::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));
            if (Config::get().general.input == config::input_types::AMBER)
            {
              g_QV<RT>(current_c, p.E, p.R, r, e_c, e_v, dE);  //calculate vdw and coulomb energy and gradients
            }
            else
            {
              g_QV<RT>(p.C, p.E, p.R, r, e_c, e_v, dE);
            }
            b *= dE;
            tmp_grad[pairlist[i].a] += b;
            tmp_grad[pairlist[i].b] -= b;
          }
          #pragma omp critical (nb_g_sum)
          {
            grad_vector += tmp_grad;
          }
        }
        e_nb += e_c+e_v;
		
        //part_energy[types::VDWC] += e_c+e_v;
        part_energy[types::CHARGE] += e_c;
        part_energy[types::VDW] += e_v;
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_cutoff 
      (
        coords::float_type &e_nb, coords::Representation_3D &grad_vector, 
        std::vector< ::tinker::refine::types::nbpair> const & pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params
      )
      {
        nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
        coords::float_type e_c(0.0), e_v(0.0);
        std::ptrdiff_t const M(pairlist.size());
        #pragma omp parallel
        {
          coords::Representation_3D tmp_grad(grad_vector.size());
          coords::virial_t tempvir(coords::empty_virial());
          #pragma omp for reduction (+: e_c, e_v)
          for (std::ptrdiff_t i=0; i<M ; ++i)  //for every pair in pairlist
          {
            coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));  //vector between the two atoms
            if (PERIODIC) boundary(b.x(), b.y(), b.z());  // for periodic boundaries: 
			                        // if the absolute value of the distance in one of the coordinates is bigger than half the box size:
			                        // subtract (or add) the box size
			                        // => absolute value of the new box size is the smallest value between these atoms in any of the boxes
            coords::float_type const rr = dot(b, b);
            coords::float_type r(0.0), fQ(0.0), fV(0.0), dE(0.0);
            if(!cutob.factors(rr, r, fQ, fV)) continue;
            r = 1.0/r;
            double current_c;   // Q_a * Q_b from AMBER
            if (Config::get().general.input == config::input_types::AMBER)
            {    // calculate Q_a * Q_b from AMBER charges (better if this would be done while building up pairlist)
              double ca = Config::get().coords.charges[pairlist[i].a];
              double cb = Config::get().coords.charges[pairlist[i].b];
              current_c = ca * cb;
            }
            ::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a), 
              refined.type(pairlist[i].b)));   // get parameters for current pair
            if (Config::get().general.input == config::input_types::AMBER)
            {
              g_QV_cutoff<RT>(current_c, p.E, p.R, r, fQ, fV, e_c, e_v, dE);  //calculate vdw and coulomb energy and gradients
            }
            else
            {
              g_QV_cutoff<RT>(p.C, p.E, p.R, r, fQ, fV, e_c, e_v, dE);  //calculate vdw and coulomb energy and gradients
            }
            
            auto const dist = b;
            b *= dE;
            tmp_grad[pairlist[i].a] += b;
            tmp_grad[pairlist[i].b] -= b;
            //Increment internal virial tensor
            coords::float_type const vxx = b.x() * dist.x();
            coords::float_type const vyx = b.x() * dist.y();
            coords::float_type const vzx = b.x() * dist.z();
            coords::float_type const vyy = b.y() * dist.y();
            coords::float_type const vzy = b.y() * dist.z();
            coords::float_type const vzz = b.z() * dist.z();
            tempvir[0][0] += vxx;
            tempvir[1][0] += vyx;
            tempvir[2][0] += vzx;
            tempvir[0][1] += vyx;
            tempvir[1][1] += vyy;
            tempvir[2][1] += vzy;
            tempvir[0][2] += vzx;
            tempvir[1][2] += vzy;
            tempvir[2][2] += vzz;
          }
          #pragma omp critical (nb_g_sum)
          {
            grad_vector += tmp_grad;
            for (int i = 0; i <= 2; i++){
              for (int k = 0; k <= 2; k++){
                part_virial[VDWC][i][k] += tempvir[i][k];
              }
            }
          }
        }
        e_nb += e_c+e_v;
        //part_energy[types::VDWC] += e_c+e_v;
        part_energy[types::CHARGE] += e_c;
        part_energy[types::VDW] += e_v;
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC, bool ALCH_OUT>
      void energy::interfaces::aco::aco_ff::g_nb_QV_pairs_fep_io 
      (
        coords::float_type &e_nb, coords::Representation_3D &grad_vector, 
        std::vector< ::tinker::refine::types::nbpair> const & pairlist,
        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params
      )
      {
        nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
        coords::float_type e_c(0.0), e_v(0.0), e_c_l(0.0), e_vdw_l(0.0), e_c_dl(0.0), e_vdw_dl(0.0);
        fepvar const & fep = coords->fep.window[coords->fep.window[0].step];
        std::ptrdiff_t const M(pairlist.size());
        #pragma omp parallel
        {
          coords::Representation_3D tmp_grad(grad_vector.size());
          coords::virial_t tempvir(coords::empty_virial());
          #pragma omp for reduction (+: e_c, e_v, e_c_l, e_c_dl, e_vdw_l, e_vdw_dl)
          for (std::ptrdiff_t i=0; i<M ; ++i)      //for every pair in pairlist
          {
            double current_c;   // Q_a * Q_b from AMBER
            if (Config::get().general.input == config::input_types::AMBER)
            {    // calculate Q_a * Q_b from AMBER charges (better if this would be done while building up pairlist)
              double ca = Config::get().coords.charges[pairlist[i].a];
              double cb = Config::get().coords.charges[pairlist[i].b];
              current_c = ca * cb;
            }
            coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));  //vector between atoms a and b
            if (PERIODIC) boundary(b.x(), b.y(), b.z());   // adjust vector to boundary conditions
            ::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));  // get parameters
            coords::float_type rr(dot(b, b)), dE, Q(0.0), V(0.0);
            coords::Cartesian_Point dist;
            if (PERIODIC)
            {
              coords::float_type fQ(0.0), fV(0.0);
              coords::float_type r(0.0);
              if(!cutob.factors(rr, r, fQ, fV)) continue;
              if (Config::get().general.input == config::input_types::AMBER)
              {
                g_QV_fep_cutoff<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout),
                  (ALCH_OUT ? fep.vin : fep.vout), fQ, fV, Q, V, dE);
                coords::float_type trash(0.0);
                g_QV_fep_cutoff<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout),
                  (ALCH_OUT ? fep.dvin : fep.dvout), fQ, fV, e_c_dl, e_vdw_dl, trash);
              }
              else
              {
                g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout),
                  (ALCH_OUT ? fep.vin : fep.vout), fQ, fV, Q, V, dE);
                coords::float_type trash(0.0);
                g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout),
                  (ALCH_OUT ? fep.dvin : fep.dvout), fQ, fV, e_c_dl, e_vdw_dl, trash);
              }
            }
            else
            {
              if (Config::get().energy.cutoff < 1000.0)
              {
                coords::float_type r(0.0);
                coords::float_type fQ(0.0), fV(0.0);
                if (cutob.factors(rr, r, fQ, fV))  //calculate r and see if r < cutoff
                {
                  if (Config::get().general.input == config::input_types::AMBER)
                  {
                    g_QV_fep_cutoff<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout),
                      (ALCH_OUT ? fep.vin : fep.vout), fQ, fV, Q, V, dE);  //calculate nb-energy(lambda)
                    coords::float_type trash(0.0);
                    g_QV_fep_cutoff<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout),
                      (ALCH_OUT ? fep.dvin : fep.dvout), fQ, fV, e_c_dl, e_vdw_dl, trash);  //calculate nb-energy(lambda+dlambda)
                  }
                  else
                  {
                    g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout),
                      (ALCH_OUT ? fep.vin : fep.vout), fQ, fV, Q, V, dE);  //calculate nb-energy(lambda)
                    coords::float_type trash(0.0);
                    g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout),
                      (ALCH_OUT ? fep.dvin : fep.dvout), fQ, fV, e_c_dl, e_vdw_dl, trash);  //calculate nb-energy(lambda+dlambda)
                  }
                }
              }
              else  //if no cutoff
              {
                coords::float_type const r= sqrt(rr);
                if (Config::get().general.input == config::input_types::AMBER)
                {
                  g_QV_fep<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout),
                    (ALCH_OUT ? fep.vin : fep.vout), Q, V, dE);  //calculate nb-energy(lambda)
                  coords::float_type trash(0.0);
                  g_QV_fep<RT>(current_c, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout),
                    (ALCH_OUT ? fep.dvin : fep.dvout), e_c_dl, e_vdw_dl, trash); //calculate nb-energy(lambda+dlambda)
                }
                else
                {
                  g_QV_fep<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout),
                    (ALCH_OUT ? fep.vin : fep.vout), Q, V, dE);  //calculate nb-energy(lambda)
                  coords::float_type trash(0.0);
                  g_QV_fep<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout),
                    (ALCH_OUT ? fep.dvin : fep.dvout), e_c_dl, e_vdw_dl, trash); //calculate nb-energy(lambda+dlambda)
                }
              }
            }
            dist = b;
            b *= dE;
            e_c_l += Q;
            e_vdw_l += V;
            e_c += Q;
            e_v += V;
            tmp_grad[pairlist[i].a] += b;
            tmp_grad[pairlist[i].b] -= b;
            //Increment internal virial tensor
            coords::float_type const vxx = b.x() * dist.x();
            coords::float_type const vyx = b.x() * dist.y();
            coords::float_type const vzx = b.x() * dist.z();
            coords::float_type const vyy = b.y() * dist.y();
            coords::float_type const vzy = b.y() * dist.z();
            coords::float_type const vzz = b.z() * dist.z();
            tempvir[0][0] += vxx;
            tempvir[1][0] += vyx;
            tempvir[2][0] += vzx;
            tempvir[0][1] += vyx;
            tempvir[1][1] += vyy;
            tempvir[2][1] += vzy;
            tempvir[0][2] += vzx;
            tempvir[1][2] += vzy;
            tempvir[2][2] += vzz;
          }
          #pragma omp critical (nb_g_sum)
          {
            grad_vector += tmp_grad;
            for (int i = 0; i <= 2; i++){
               for (int k = 0; k <= 2; k++){
                  part_virial[VDWC][i][k] += tempvir[i][k];
               }
            }
          }
        }
        e_nb += e_c+e_v;
        coords->fep.feptemp.e_c_l1 += e_c_l;    //lambda (Coulomb energy)
        coords->fep.feptemp.e_c_l2 += e_c_dl;   //lambda + dlambda (Coulomb energy)
        coords->fep.feptemp.e_vdw_l1 += e_vdw_l;  //lambda (vdW energy)
        coords->fep.feptemp.e_vdw_l2 += e_vdw_dl;  //lambda + dlambda (vdW energy)
        part_energy[types::CHARGE] += e_c;  //gradients (coulomb)
        part_energy[types::VDW] += e_v;     //gradients (vdW)
      }
	  //      template < ::tinker::parameter::radius_types::T RT>
	  //      void energy::interfaces::aco::aco_ff::pme_direct(double &e_ch,
	  //        coords::Representation_3D &grad_vector,
	  //        std::vector< ::tinker::refine::types::nbpair> const & pairlist,
	  //        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params, double & e_vdw)
	  //      {
	  //       
	  //        coords::float_type e_c(0.0), e_v(0.0);
	  //
	  //#pragma omp parallel
	  //        {
	  //          nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
	  //          coords::Representation_3D tmp_grad(grad_vector.size());
	  //          coords::virial_t tempvir(coords::empty_virial());
	  //          ptrdiff_t const M(pairlist.size());
	  //          double dist, d, r, achg(0.0), distbuf, rbuf, factor, ewald, error;
	  //          double energy, den, cutoff, cutoff2, scale;
	  //          cutoff = Config::get().energy.cutoff;
	  //          cutoff2 = cutoff*cutoff;
	  //          coords::Cartesian_Point vir;
	  //          coords::float_type dV;
	  //          coords::float_type fQ(0.0), fV(0.0);
	  //
	  //          #pragma omp for reduction (+: e_c, e_v)
	  //          for (ptrdiff_t i = 0; i < M; ++i)
	  //          {
	  //            coords::Cartesian_Point b(coords->xyz(pairlist[i].b) - coords->xyz(pairlist[i].a));
	  //            boundary(b.x(), b.y(), b.z());
	  //            //dist =    b.scalar(b);
   //             dist = scon::dot(b, b);
	  //            if (dist >= cutoff2) continue;
	  //            r = sqrt(dist);
	  //            d = 1 / r;
	  //            rbuf = r + achg;
	  //            vir = b;
	  //            distbuf = rbuf*rbuf;
	  //            factor = 332.0716 * coords->pme.pmetemp.atomcharges[pairlist[i].a] * coords->pme.pmetemp.atomcharges[pairlist[i].b];
	  //            ewald = coords->pme.pmetemp.ewaldcoeff * r;
	  //            error = erfc(ewald);
	  //            energy = (factor / rbuf) * error;
	  //            e_c += energy;
	  //            // van der waals part
	  //            ::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));
	  //            if (!cutob.factors(dist, r, fQ, fV)) continue;
	  //            e_v += gV<RT>(p.E, p.R, d, dV)*fV;
	  //            // gradients 
	  //            den = -factor * (error / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r);
	  //            b *= (dV*fV + den) * d;
	  //            tmp_grad[pairlist[i].a] -= b;
	  //            tmp_grad[pairlist[i].b] += b;
	  //            // virial components
	  //            coords::float_type const vxx = b.x() * vir.x();
	  //            coords::float_type const vyx = b.x() * vir.y();
	  //            coords::float_type const vzx = b.x() * vir.z();
	  //            coords::float_type const vyy = b.y() * vir.y();
	  //            coords::float_type const vzy = b.y() * vir.z();
	  //            coords::float_type const vzz = b.z() * vir.z();
	  //            tempvir[0][0] += vxx;
	  //            tempvir[1][0] += vyx;
	  //            tempvir[2][0] += vzx;
	  //            tempvir[0][1] += vyx;
	  //            tempvir[1][1] += vyy;
	  //            tempvir[2][1] += vzy;
	  //            tempvir[0][2] += vzx;
	  //            tempvir[1][2] += vzy;
	  //            tempvir[2][2] += vzz;
	  //          }
	  //          #pragma omp critical (nb_g_sum)
	  //          {
	  //            grad_vector += tmp_grad;
	  //            for (int i = 0; i <= 2; i++){
	  //              for (int k = 0; k <= 2; k++){
	  //                part_virial[VDWC][i][k] += tempvir[i][k];
	  //              }
	  //            }
	  //          }// end of critical region
	  //        }// end of parallel region
	  //        coords->fep.feptemp.e_c_l1 += e_c;
	  //        coords->fep.feptemp.e_c_l2 += e_c;
	  //        coords->fep.feptemp.e_vdw_l1 += e_v;
	  //        coords->fep.feptemp.e_vdw_l2 += e_v;
	  //        e_ch += e_c;
	  //        e_vdw += e_v;
	  //      }
	  //
	  //
	  //      template < ::tinker::parameter::radius_types::T RT, bool ALCH_OUT>
	  //      void energy::interfaces::aco::aco_ff::pme_direct_fep(double &e_ch,
	  //        coords::Representation_3D &grad_vector,
	  //        std::vector< ::tinker::refine::types::nbpair> const & pairlist,
	  //        scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params, double & e_vdw)
	  //      {
	  //
	  //        coords::float_type e_c(0.0), e_v(0.0), e_c_l(0.0), e_vdw_l(0.0), e_c_dl(0.0), e_vdw_dl(0.0);
	  //        fepvar const & fep = coords->fep.window[coords->fep.window[0].step];
	  //
	  //#pragma omp parallel
	  //        {
	  //          coords::Representation_3D tmp_grad(grad_vector.size());
	  //          coords::virial_t tempvir(coords::empty_virial());
	  //          ptrdiff_t const M(pairlist.size());
	  //          double dist, d, r, achg(0.0), distbuf, rbuf, factor, ewald, error;
	  //          double energy, energyc2, den, cutoff, scale, e_v_temp, trash;
	  //          coords::Cartesian_Point vir;
	  //          coords::float_type dV, alche2;
	  //          coords::float_type fQ(0.0), fV(0.0), fV2, switch2(0.0), cutoff2(0.0), switchfactor(0.0);        
	  //          switch2 = Config::get().energy.switchdist *  Config::get().energy.switchdist;
	  //          cutoff2 = Config::get().energy.cutoff *  Config::get().energy.cutoff;
	  //          double temp = cutoff2*cutoff2 - switch2*switch2;
	  //          switchfactor = 1 / (temp*temp*temp);
	  //    
	  //#pragma omp for reduction (+: e_c, e_v, e_c_l, e_c_dl, e_vdw_l, e_vdw_dl)
	  //          for (ptrdiff_t i = 0; i < M; ++i)
	  //          {
	  //            coords::Cartesian_Point b(coords->xyz(pairlist[i].b) - coords->xyz(pairlist[i].a));
	  //            boundary(b.x(), b.y(), b.z());
	  //            //dist = b.scalar(b);
   //             dist = scon::dot(b, b);
	  //            if (dist >= cutoff2) continue;
	  //            r = sqrt(dist);
	  //            d = 1 / r;
	  //            rbuf = r + achg;
	  //            vir = b;
	  //            distbuf = rbuf*rbuf;
	  //            factor = 332.0716 * coords->pme.pmetemp.atomcharges[pairlist[i].a] * coords->pme.pmetemp.atomcharges[pairlist[i].b];
	  //            ewald = coords->pme.pmetemp.ewaldcoeff * r;
	  //            error = erfc(ewald);
	  //            energy = (factor / rbuf) * error;
	  //            // FEP charge energies
	  //            e_c += energy * (ALCH_OUT ? fep.ein : fep.eout);
	  //            // increment FEP charge energies
	  //            e_c_l += energy * (ALCH_OUT ? fep.ein : fep.eout);
	  //            e_c_dl += energy * (ALCH_OUT ? fep.dein : fep.deout);
	  //            // van der waals part
	  //            ::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));
	  //            fV = (dist > switch2) ? ((cutoff2 - dist) * (cutoff2 - dist) * (cutoff2 - 3 * switch2 + 2 * dist)) * switchfactor : 1.0;
	  //            fV2 = (dist > switch2) ? 12 * switchfactor * (cutoff2 - dist)*(dist - switch2) : 0.0;
	  //            e_v_temp = gV_fep_cut<RT>(p.E, p.R, r, (ALCH_OUT ? fep.vin : fep.vout), (ALCH_OUT ? fep.dvin : fep.dvout), dV, alche2, fV, fV2);
	  //            e_v += e_v_temp;
	  //            e_vdw_l += e_v_temp;
	  //            e_vdw_dl += alche2;
	  //            // gradients 
	  //            den = -factor * (error / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r) * (ALCH_OUT ? fep.ein : fep.eout);
	  //            b *= (dV + den*d);
	  //            tmp_grad[pairlist[i].a] -= b;
	  //            tmp_grad[pairlist[i].b] += b;
	  //            // virial components
	  //            coords::float_type const vxx = b.x() * vir.x();
	  //            coords::float_type const vyx = b.x() * vir.y();
	  //            coords::float_type const vzx = b.x() * vir.z();
	  //            coords::float_type const vyy = b.y() * vir.y();
	  //            coords::float_type const vzy = b.y() * vir.z();
	  //            coords::float_type const vzz = b.z() * vir.z();
	  //            tempvir[0][0] += vxx;
	  //            tempvir[1][0] += vyx;
	  //            tempvir[2][0] += vzx;
	  //            tempvir[0][1] += vyx;
	  //            tempvir[1][1] += vyy;
	  //            tempvir[2][1] += vzy;
	  //            tempvir[0][2] += vzx;
	  //            tempvir[1][2] += vzy;
	  //            tempvir[2][2] += vzz;
	  //          }
	  //#pragma omp critical (nb_g_sum)
	  //          {
	  //            grad_vector += tmp_grad;
	  //            for (int i = 0; i <= 2; i++){
	  //              for (int k = 0; k <= 2; k++){
	  //                part_virial[VDWC][i][k] += tempvir[i][k];
	  //              }
	  //            }
	  //          }// end of critical region
	  //        }// end of parallel region
	  //        coords->fep.feptemp.e_c_l1 += e_c_l;
	  //        coords->fep.feptemp.e_c_l2 += e_c_dl;
	  //        coords->fep.feptemp.e_vdw_l1 += e_vdw_l;
	  //        coords->fep.feptemp.e_vdw_l2 += e_vdw_dl;
	  //        e_ch += e_c;
	  //        e_vdw += e_v;
	  //      }
	  //
	  // 
	  //
	  //  void energy::interfaces::aco::aco_ff::pme_direct_scaled_fep(double & e_nb,
	  //  coords::Representation_3D &grad_vector)
	  //{
	  //  double dist, r, achg, distbuf, rbuf, factor, ewald, error;
	  //  double energy, den, cutoff, cutoff2, scale;
	  //  cutoff = Config::get().energy.cutoff;
	  //  cutoff2 = cutoff*cutoff;
	  //  achg = 0.0;
	  //  coords::Cartesian_Point vir;
	  //  coords::float_type e_c(0.0), e_v(0.0), e_c_l(0.0), e_vdw_l(0.0), e_c_dl(0.0), e_vdw_dl(0.0);
	  //  fepvar const & fep = coords->fep.window[coords->fep.window[0].step];
	  //  // 1-2 scaled interactions 
	  //  for (auto const & bond : refined.bonds())
	  //  {
	  //    scale = -1.0;
	  //    coords::Cartesian_Point b((coords->xyz(bond.atoms[0])) - (coords->xyz(bond.atoms[1])));
	  //    boundary(b.x(), b.y(), b.z());
	  //    dist = b.x()*b.x() + b.y()*b.y() + b.z()*b.z();
	  //    vir = b;
	  //    r = sqrt(dist);
	  //    rbuf = r + achg;
	  //    distbuf = rbuf*rbuf;
	  //    factor = 332.0716 * coords->pme.pmetemp.atomcharges[bond.atoms[0]] * coords->pme.pmetemp.atomcharges[bond.atoms[1]];
	  //    ewald = coords->pme.pmetemp.ewaldcoeff * r;
	  //    error = erfc(ewald);
	  //    energy = (factor / rbuf) * (error + scale);
	  //    // atom is in
	  //    if (coords->pme.pmetemp.feinf[bond.atoms[0]].flag == 1 || coords->pme.pmetemp.feinf[bond.atoms[1]].flag == 1)
	  //    {
	  //      e_c_l += energy * fep.eout;
	  //      e_c_dl += energy * fep.deout;
	  //      e_nb += energy * fep.eout;
	  //      den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r) * fep.eout;
	  //      den /= r;
	  //      b *= den;
	  //      grad_vector[bond.atoms[0]] += b;
	  //      grad_vector[bond.atoms[1]] -= b;;
	  //      coords::float_type const vxx = b.x() * vir.x();
	  //      coords::float_type const vyx = b.x() * vir.y();
	  //      coords::float_type const vzx = b.x() * vir.z();
	  //      coords::float_type const vyy = b.y() * vir.y();
	  //      coords::float_type const vzy = b.y() * vir.z();
	  //      coords::float_type const vzz = b.z() * vir.z();
	  //      part_virial[VDWC][0][0] += vxx;
	  //      part_virial[VDWC][1][0] += vyx;
	  //      part_virial[VDWC][2][0] += vzx;
	  //      part_virial[VDWC][0][1] += vyx;
	  //      part_virial[VDWC][1][1] += vyy;
	  //      part_virial[VDWC][2][1] += vzy;
	  //      part_virial[VDWC][0][2] += vzx;
	  //      part_virial[VDWC][1][2] += vzy;
	  //      part_virial[VDWC][2][2] += vzz;
	  //    }
	  //    else if (coords->pme.pmetemp.feinf[bond.atoms[0]].flag == 0 || coords->pme.pmetemp.feinf[bond.atoms[1]].flag == 0)
	  //    {
	  //      e_c_l += energy * fep.ein;
	  //      e_c_dl += energy * fep.dein;
	  //      e_nb += energy * fep.ein;
	  //      den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r) * fep.ein;
	  //      den /= r;
	  //      b *= den;
	  //      grad_vector[bond.atoms[0]] += b;
	  //      grad_vector[bond.atoms[1]] -= b;;
	  //      coords::float_type const vxx = b.x() * vir.x();
	  //      coords::float_type const vyx = b.x() * vir.y();
	  //      coords::float_type const vzx = b.x() * vir.z();
	  //      coords::float_type const vyy = b.y() * vir.y();
	  //      coords::float_type const vzy = b.y() * vir.z();
	  //      coords::float_type const vzz = b.z() * vir.z();
	  //      part_virial[VDWC][0][0] += vxx;
	  //      part_virial[VDWC][1][0] += vyx;
	  //      part_virial[VDWC][2][0] += vzx;
	  //      part_virial[VDWC][0][1] += vyx;
	  //      part_virial[VDWC][1][1] += vyy;
	  //      part_virial[VDWC][2][1] += vzy;
	  //      part_virial[VDWC][0][2] += vzx;
	  //      part_virial[VDWC][1][2] += vzy;
	  //      part_virial[VDWC][2][2] += vzz;
	  //    }
	  //    else
	  //    {
	  //        e_c_l += energy;
	  //        e_c_dl += energy;
	  //        e_nb += energy;
	  //        den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r);
	  //        den /= r;
	  //        b *= den;
	  //        grad_vector[bond.atoms[0]] += b;
	  //        grad_vector[bond.atoms[1]] -= b;;
	  //        coords::float_type const vxx = b.x() * vir.x();
	  //        coords::float_type const vyx = b.x() * vir.y();
	  //        coords::float_type const vzx = b.x() * vir.z();
	  //        coords::float_type const vyy = b.y() * vir.y();
	  //        coords::float_type const vzy = b.y() * vir.z();
	  //        coords::float_type const vzz = b.z() * vir.z();
	  //        part_virial[VDWC][0][0] += vxx;
	  //        part_virial[VDWC][1][0] += vyx;
	  //        part_virial[VDWC][2][0] += vzx;
	  //        part_virial[VDWC][0][1] += vyx;
	  //        part_virial[VDWC][1][1] += vyy;
	  //        part_virial[VDWC][2][1] += vzy;
	  //        part_virial[VDWC][0][2] += vzx;
	  //        part_virial[VDWC][1][2] += vzy;
	  //        part_virial[VDWC][2][2] += vzz;
	  //    }
	  //  }
	  //  // 1-3 scaled interactions
	  //  for (auto const & angle : refined.angles())
	  //  {
	  //    scale = -1.0;
	  //    coords::Cartesian_Point b((coords->xyz(angle.atoms[0])) - (coords->xyz(angle.atoms[2])));
	  //    boundary(b.x(), b.y(), b.z());
	  //    dist = b.x()*b.x() + b.y()*b.y() + b.z()*b.z();
	  //    r = sqrt(dist);
	  //    vir = b;
	  //    rbuf = r + achg;
	  //    distbuf = rbuf*rbuf;
	  //    factor = 332.0716 * coords->pme.pmetemp.atomcharges[angle.atoms[2]] * coords->pme.pmetemp.atomcharges[angle.atoms[0]];
	  //    ewald = coords->pme.pmetemp.ewaldcoeff * r;
	  //    error = erfc(ewald);
	  //    energy = (factor / rbuf) * (error + scale);
	  //    if (coords->pme.pmetemp.feinf[angle.atoms[0]].flag == 1 || coords->pme.pmetemp.feinf[angle.atoms[1]].flag == 1 ||  coords->pme.pmetemp.feinf[angle.atoms[2]].flag == 1)
	  //    {
	  //      e_c_l += energy * fep.eout;
	  //      e_c_dl += energy * fep.deout;
	  //      e_nb += energy * fep.eout;
	  //      den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r) * fep.eout;
	  //      den /= r;
	  //      b *= den;
	  //      grad_vector[angle.atoms[0]] += b;
	  //      grad_vector[angle.atoms[2]] -= b;;
	  //      coords::float_type const vxx = b.x() * vir.x();
	  //      coords::float_type const vyx = b.x() * vir.y();
	  //      coords::float_type const vzx = b.x() * vir.z();
	  //      coords::float_type const vyy = b.y() * vir.y();
	  //      coords::float_type const vzy = b.y() * vir.z();
	  //      coords::float_type const vzz = b.z() * vir.z();
	  //      part_virial[VDWC][0][0] += vxx;
	  //      part_virial[VDWC][1][0] += vyx;
	  //      part_virial[VDWC][2][0] += vzx;
	  //      part_virial[VDWC][0][1] += vyx;
	  //      part_virial[VDWC][1][1] += vyy;
	  //      part_virial[VDWC][2][1] += vzy;
	  //      part_virial[VDWC][0][2] += vzx;
	  //      part_virial[VDWC][1][2] += vzy;
	  //      part_virial[VDWC][2][2] += vzz;
	  //    }
	  //    else if (coords->pme.pmetemp.feinf[angle.atoms[0]].flag == 0 || coords->pme.pmetemp.feinf[angle.atoms[1]].flag == 0 || coords->pme.pmetemp.feinf[angle.atoms[2]].flag == 0)
	  //    {
	  //      e_c_l += energy * fep.ein;
	  //      e_c_dl += energy * fep.dein;
	  //      e_nb += energy * fep.ein;
	  //      den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r) * fep.ein;
	  //      den /= r;
	  //      b *= den;
	  //      grad_vector[angle.atoms[0]] += b;
	  //      grad_vector[angle.atoms[2]] -= b;;
	  //      coords::float_type const vxx = b.x() * vir.x();
	  //      coords::float_type const vyx = b.x() * vir.y();
	  //      coords::float_type const vzx = b.x() * vir.z();
	  //      coords::float_type const vyy = b.y() * vir.y();
	  //      coords::float_type const vzy = b.y() * vir.z();
	  //      coords::float_type const vzz = b.z() * vir.z();
	  //      part_virial[VDWC][0][0] += vxx;
	  //      part_virial[VDWC][1][0] += vyx;
	  //      part_virial[VDWC][2][0] += vzx;
	  //      part_virial[VDWC][0][1] += vyx;
	  //      part_virial[VDWC][1][1] += vyy;
	  //      part_virial[VDWC][2][1] += vzy;
	  //      part_virial[VDWC][0][2] += vzx;
	  //      part_virial[VDWC][1][2] += vzy;
	  //      part_virial[VDWC][2][2] += vzz;
	  //    }
	  //    else
	  //    {
	  //      {
	  //        e_c_l += energy;
	  //        e_c_dl += energy;
	  //        e_nb += energy;
	  //        den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r);
	  //        den /= r;
	  //        b *= den;
	  //        grad_vector[angle.atoms[0]] += b;
	  //        grad_vector[angle.atoms[2]] -= b;;
	  //        coords::float_type const vxx = b.x() * vir.x();
	  //        coords::float_type const vyx = b.x() * vir.y();
	  //        coords::float_type const vzx = b.x() * vir.z();
	  //        coords::float_type const vyy = b.y() * vir.y();
	  //        coords::float_type const vzy = b.y() * vir.z();
	  //        coords::float_type const vzz = b.z() * vir.z();
	  //        part_virial[VDWC][0][0] += vxx;
	  //        part_virial[VDWC][1][0] += vyx;
	  //        part_virial[VDWC][2][0] += vzx;
	  //        part_virial[VDWC][0][1] += vyx;
	  //        part_virial[VDWC][1][1] += vyy;
	  //        part_virial[VDWC][2][1] += vzy;
	  //        part_virial[VDWC][0][2] += vzx;
	  //        part_virial[VDWC][1][2] += vzy;
	  //        part_virial[VDWC][2][2] += vzz;
	  //      }
	  //      /* den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r);
	  //       den /= r;
	  //       b *= den;
	  //       e_nb += energy;
	  //       grad_vector[angle.atoms[0]] += b;
	  //       grad_vector[angle.atoms[2]] -= b;
	  //       coords::float_type const vxx = b.x() * vir.x();
	  //       coords::float_type const vyx = b.x() * vir.y();
	  //       coords::float_type const vzx = b.x() * vir.z();
	  //       coords::float_type const vyy = b.y() * vir.y();
	  //       coords::float_type const vzy = b.y() * vir.z();
	  //       coords::float_type const vzz = b.z() * vir.z();
	  //       part_virial[VDWC][0][0] += vxx;
	  //       part_virial[VDWC][1][0] += vyx;
	  //       part_virial[VDWC][2][0] += vzx;
	  //       part_virial[VDWC][0][1] += vyx;
	  //       part_virial[VDWC][1][1] += vyy;
	  //       part_virial[VDWC][2][1] += vzy;
	  //       part_virial[VDWC][0][2] += vzx;
	  //       part_virial[VDWC][1][2] += vzy;
	  //       part_virial[VDWC][2][2] += vzz;*/
	  //    }
	  //  }
	  //  coords->fep.feptemp.e_c_l1 += e_c_l;
	  //  coords->fep.feptemp.e_c_l2 += e_c_dl;
	  //
	  //}

	
#endif

    }
  }
}
