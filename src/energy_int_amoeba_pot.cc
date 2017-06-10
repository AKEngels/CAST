#include <cmath>
#include <stddef.h>
#include <stdexcept>
#include <cstdlib>
#include "energy_int_amoeba.h"
#include "scon_angle.h"
#include "configuration.h"
#include "scon_utility.h"
#include "scon_c3.h"
#include <algorithm>
#include "math.h"


#ifdef _MSC_VER
#pragma warning(disable: 4996)
#endif

/****************************************
*                                       *
*                                       *
*               Generics                *
*                                       *
*                                       *
*****************************************/

double energy::interfaces::amoeba::amoeba_ff::e(void)
{

  pre();


  if (Config::get().energy.spackman.on)

  {
    Spackman_mol();
    Spackman_vec();	
	Spackman1(); 

  }



  multipole_sites();


  rot_matrix(coords->xyz());


  refine_pair_lists();
  e_ind();



  e_perm();
  calc<0>();
  post();
  return energy;
}







double energy::interfaces::amoeba::amoeba_ff::g(void)
{
  pre();
  if (Config::get().energy.spackman.on)
  {
    Spackman_mol();
    Spackman_vec();
	Spackman1();
    
    if (Config::get().energy.spackman.interp)
      SpackmanGrad_3();
    else Spackman_GRAD();

  }
  multipole_sites();
  rot_matrix(coords->xyz());
  refine_pair_lists();
  e_ind();
  e_perm();
  calc<1>();
  post();
  return energy;
}

double energy::interfaces::amoeba::amoeba_ff::h(void)
{
  pre();
  //... 
  post();
  return energy;
}

double energy::interfaces::amoeba::amoeba_ff::o(void)
{
  throw std::runtime_error("amoeba_ff doesn't provide any optimization routines.");
}

template<size_t DERIV>
void energy::interfaces::amoeba::amoeba_ff::calc(void)
{

  //#if defined (_OPENMP)
  //    #pragma omp parallel sections
  //    {
  //      #pragma omp section
  //       part_energy[BOND]       = f_12<DERIV>();
  //      #pragma omp section
  //        part_energy[ANGLE]      = f_13_a<DERIV>();
  //      #pragma omp section
  //        part_energy[UREY]       = f_13_u<DERIV>();
  //      #pragma omp section
  //        part_energy[TORSION]    = f_14<DERIV>();
  //      #pragma omp section
  //        part_energy[TORSION]    = f_14p();
  //      #pragma omp section 
  //        part_energy[IMPTORSION] = g_it();
  //      #pragma omp section
  //        part_energy[IMPTORSION] = f_it<DERIV>();
  //      #pragma omp section 
  //        part_energy[IMPROPER]   = g_imp();
  //
  //    }
  //#else
  part_energy[BOND] = f_12<DERIV>();
  part_energy[ANGLE] = f_13_a<DERIV>();
  part_energy[UREY] = f_13_u<DERIV>();
  part_energy[TORSION] = f_14<DERIV>();
  part_energy[MULTIPOLE] = em;
  part_energy[POLARIZE] = ep;
  part_energy[OPBEND] = f_oop<DERIV>();
  if (Config::get().energy.spackman.on)
  {
	
    if (Config::get().energy.spackman.interp) part_energy[SHORTRANGE] = Spackman_energy_analytical();


    else part_energy[SHORTRANGE] = Spackman_Energy();

  }

  //part_energy[types::IMPTORSION]   = f_it<DERIV>();
    //part_energy[types::IMPROPER]     = g_imp();
//#endif
  if (cparams.radiustype() == ::tinker::parameter::radius_types::R_MIN)
  {
    g_nb< ::tinker::parameter::radius_types::R_MIN>();
  }
  else
  {
    g_nb< ::tinker::parameter::radius_types::SIGMA>();
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
    namespace amoeba
    {

      //! Bonding Energy
      template<>
      double energy::interfaces::amoeba::amoeba_ff::f_12<0>(void)
      {
        double E(0.0);
        for (auto const & bond : refined.bonds())
        {
          coords::Cartesian_Point const bv(coords->xyz(bond.atoms[0]) - coords->xyz(bond.atoms[1]));
          double d = len(bv);
          double r = d - bond.ideal;
          E += bond.force*r*r*(1.0 - 2.55*(r)+3.793125*(r*r));
          if (std::abs(r) > 1.0) integrity = false;
        }
        return E;
      }

      //! Bonding Energy+Gradients
      template<>
      double energy::interfaces::amoeba::amoeba_ff::f_12<1>(void)
      {
        double E(0.0);
        for (auto bond : refined.bonds())
        {
          coords::Cartesian_Point const bv(coords->xyz(bond.atoms[0]) - coords->xyz(bond.atoms[1]));
          double d = len(bv);
          double r = d - bond.ideal;

          E += bond.force*r*r*(1.0 - 2.55*(r)+3.793125*(r*r));
          double dE = bond.force*r*((2.0 + (3.0*(-2.55))*r) + ((4.0*3.793125)*(r*r)));
          if (std::abs(d) > 1.0e-12)
          {
            if (std::abs(r) > 1.0) integrity = false;
            dE /= d;
            auto gv = bv*dE;
            part_grad[BOND][bond.atoms[0]] += gv;
            part_grad[BOND][bond.atoms[1]] -= gv;
            //increment internal virial tensor
            double vxx = bv.x() * gv.x();
            double vyx = bv.y() * gv.x();
            double vzx = bv.z() * gv.x();
            double vyy = bv.y() * gv.y();
            double vzy = bv.z() * gv.y();
            double vzz = bv.z() * gv.z();

            part_virial[BOND][0][0] += vxx;
            part_virial[BOND][1][0] += vyx;
            part_virial[BOND][2][0] += vzx;
            part_virial[BOND][0][1] += vyx;
            part_virial[BOND][1][1] += vyy;
            part_virial[BOND][2][1] += vzy;
            part_virial[BOND][0][2] += vzx;
            part_virial[BOND][1][2] += vzy;
            part_virial[BOND][2][2] += vzz;
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
      double energy::interfaces::amoeba::amoeba_ff::f_13_a<0>(void)
      {
        double E(0.0);
        for (auto const & ang13 : refined.angles())
        {
          coords::Cartesian_Point const
            av1(coords->xyz(ang13.atoms[0]) - coords->xyz(ang13.atoms[1])),
            av2(coords->xyz(ang13.atoms[2]) - coords->xyz(ang13.atoms[1]));
          coords::float_type const d(scon::angle(av1, av2).degrees() - ang13.ideal);
          coords::float_type const r(d*SCON_PI180);
          E += (SCON_PI180*SCON_PI180)*ang13.force*d*d*(1.0 - 0.014*d + 5.60E-05*(d*d) - 7.0E-07*(d*d*d) + 2.2E-08*(d*d*d*d));

          if (std::abs(d) > 30.0) integrity = false;
        }
        return E;
      }

      template<>
      double energy::interfaces::amoeba::amoeba_ff::f_13_a<1>(void)
      {
        double E(0.0);
        double vxx, vyx, vzx, vyy, vzy, vzz;
        for (auto angle : refined.angles())
        {
          coords::Cartesian_Point const
            av1(coords->xyz(angle.atoms[0]) - coords->xyz(angle.atoms[1])),
            av2(coords->xyz(angle.atoms[2]) - coords->xyz(angle.atoms[1]));
          double d = scon::angle(av1, av2).degrees();
          d -= angle.ideal;
          E += (SCON_PI180*SCON_PI180)*angle.force*d*d*(1.0 - 0.014*d + 5.60E-05*(d*d) - 7.0E-07*(d*d*d) + 2.2E-08*(d*d*d*d));
          double dE = (SCON_PI180*SCON_PI180)*angle.force*d*SCON_180PI*(2.0 + (3.0*(-0.014*d)) + 4.0*5.60E-05*(d*d) + 5.0*-7.0E-07*(d*d*d) + 6.0*2.2E-08*(d*d*d*d));
          coords::Cartesian_Point cv(av1);
          cv = cross(cv, av2);
          double const cvl(len(cv));
          if (std::abs(cvl) > 0.0)
          {
            if (std::abs(d) > 30.0) integrity = false;
            dE *= 1.0 / cvl;
            auto gv1 = cross(av1, cv);
            auto gv2 = cross(cv, av2);
            gv1 *= dE / dot(av1, av1);
            gv2 *= dE / dot(av2, av2);
            part_grad[ANGLE][angle.atoms[0]] += gv1;
            part_grad[ANGLE][angle.atoms[2]] += gv2;
            part_grad[ANGLE][angle.atoms[1]] += -(gv2 + gv1);
            //increment internal virial tensor
            vxx = av1.x() * gv1.x() + av2.x() * gv2.x();
            vyx = av1.y() * gv1.x() + av2.y() * gv2.x();
            vzx = av1.z() * gv1.x() + av2.z() * gv2.x();
            vyy = av1.y() * gv1.y() + av2.y() * gv2.y();
            vzy = av1.z() * gv1.y() + av2.z() * gv2.y();
            vzz = av1.z() * gv1.z() + av2.z() * gv2.z();
            part_virial[ANGLE][0][0] += vxx;
            part_virial[ANGLE][1][0] += vyx;
            part_virial[ANGLE][2][0] += vzx;
            part_virial[ANGLE][0][1] += vyx;
            part_virial[ANGLE][1][1] += vyy;
            part_virial[ANGLE][2][1] += vzy;
            part_virial[ANGLE][0][2] += vzx;
            part_virial[ANGLE][1][2] += vzy;
            part_virial[ANGLE][2][2] += vzz;
          }
          else
          {
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
      double energy::interfaces::amoeba::amoeba_ff::f_13_u<0>(void)
      {
        double E(0.0);
        for (auto urey : refined.ureys())
        {
          coords::Cartesian_Point const bv(coords->xyz(urey.atoms[0]) - coords->xyz(urey.atoms[1]));
          coords::float_type const d = len(bv);
          coords::float_type const r = d - urey.ideal;
          E += urey.force*r*r;
          if (fabs(d) < 1.0e-8) integrity = false;
        }
        return E;
      }

      template<>
      double energy::interfaces::amoeba::amoeba_ff::f_13_u<1>(void)
      {
        double E(0.0);
        double vxx, vyx, vzx, vyy, vzy, vzz;
        for (auto urey : refined.ureys())
        {
          coords::Cartesian_Point const bv(coords->xyz(urey.atoms[0]) - coords->xyz(urey.atoms[1]));
          coords::float_type const d = len(bv);
          coords::float_type const r = d - urey.ideal;
          double dE = urey.force*r;
          E += dE*r;
          dE *= 2;
          if (fabs(d) > 0)
          {
            dE /= d;
            auto gv = bv*dE;
            part_grad[UREY][urey.atoms[0]] += gv;
            part_grad[UREY][urey.atoms[1]] -= gv;
            //increment internal virial tensor
            vxx = bv.x() * gv.x();
            vyx = bv.y() * gv.x();
            vzx = bv.z() * gv.x();
            vyy = bv.y() * gv.y();
            vzy = bv.z() * gv.y();
            vzz = bv.z() * gv.z();
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
      double energy::interfaces::amoeba::amoeba_ff::f_14<0>(void)
      {
        coords::float_type E(0.0);
        for (auto const & torsion : refined.torsions())
        {
          // Get bonding vectors
          coords::Cartesian_Point const b01(coords->xyz(torsion.atoms[1]) - coords->xyz(torsion.atoms[0]));
          coords::Cartesian_Point const b12(coords->xyz(torsion.atoms[2]) - coords->xyz(torsion.atoms[1]));
          coords::Cartesian_Point const b23(coords->xyz(torsion.atoms[3]) - coords->xyz(torsion.atoms[2]));
          // Cross terms
          coords::Cartesian_Point const t = cross(b01, b12);
          coords::Cartesian_Point const u = cross(b12, b23);
          // Get length and variations
          coords::float_type const tl2 = dot(t, t);
          coords::float_type const ul2 = dot(u, u);
          // ...
          coords::float_type const tlul = std::sqrt(tl2*ul2);
          coords::float_type const r12 = len(b12);
          // cross of cross
          coords::Cartesian_Point const tu = cross(t, u);
          // scalar and length variations
          coords::float_type const cos_scalar0 = dot(t, u);
          coords::float_type const cos_scalar1 = tlul;
          coords::float_type const sin_scalar0 = dot(b12, tu);
          coords::float_type const sin_scalar1 = r12*tlul;
          // check whether 
          if (std::abs(cos_scalar1) < 1.0e-8 || std::abs(sin_scalar1) < 1.0e-8)
          {
            integrity = false;
            continue;
          }
          // Get multiple sine and cosine values
          coords::float_type cos[7], sin[7];
          cos[1] = cos_scalar0 / cos_scalar1;
          sin[1] = sin_scalar0 / sin_scalar1;
          for (size_t j(2U); j <= torsion.p.max_order; ++j)
          {
            size_t const k = j - 1;
            sin[j] = sin[k] * cos[1] + cos[k] * sin[1];
            cos[j] = cos[k] * cos[1] - sin[k] * sin[1];
          }

          coords::float_type tE(0.0);
          //cout << "Number?:  " << torsions[i].paramPtr->n << std::endl;
          for (size_t j(0U); j < torsion.p.number; ++j)
          {
            coords::float_type const F = torsion.p.force[j] * cparams.torsionunit();
            size_t const k = torsion.p.order[j];
            coords::float_type const l = std::abs(torsion.p.ideal[j]) > 0.0 ? -1.0 : 1.0;
            tE += F * (1.0 + cos[k] * l);
          }
          E += tE;
        }
        return E;
      }

      template<>
      double energy::interfaces::amoeba::amoeba_ff::f_14<1>(void)
      {
        double E(0.0);
        double vxx, vyx, vzx, vyy, vzy, vzz;
        coords::Cartesian_Point vir1, vir2, vir3, vir4;
        for (auto torsion : refined.torsions())
        {

          coords::Cartesian_Point b01(coords->xyz(torsion.atoms[1]) - coords->xyz(torsion.atoms[0]));
          coords::Cartesian_Point b12(coords->xyz(torsion.atoms[2]) - coords->xyz(torsion.atoms[1]));
          coords::Cartesian_Point b23(coords->xyz(torsion.atoms[3]) - coords->xyz(torsion.atoms[2]));
          coords::Cartesian_Point b02(coords->xyz(torsion.atoms[2]) - coords->xyz(torsion.atoms[0]));
          coords::Cartesian_Point b13(coords->xyz(torsion.atoms[3]) - coords->xyz(torsion.atoms[1]));

          coords::Cartesian_Point t = cross(b01, b12);
          coords::Cartesian_Point u = cross(b12, b23);
          double tl2 = dot(t, t);
          double ul2 = dot(u, u);
          double tlul = std::sqrt(tl2*ul2);
          double r12 = len(b12);
          coords::Cartesian_Point tu = cross(t, u);

          double cos_scalar0 = dot(t, u);
          double cos_scalar1 = tlul;

          double sin_scalar0 = dot(b12, tu);
          double sin_scalar1 = r12*tlul;

          if (std::abs(cos_scalar1) < 0.0000001 || std::abs(sin_scalar1) < 0.0000001) continue;
          double cos[7], sin[7];
          cos[1] = cos_scalar0 / cos_scalar1;
          sin[1] = sin_scalar0 / sin_scalar1;

          for (size_t j(2U); j <= torsion.p.max_order; ++j)
          {
            size_t k = j - 1;
            sin[j] = sin[k] * cos[1] + cos[k] * sin[1];
            cos[j] = cos[k] * cos[1] - sin[k] * sin[1];
          }

          double tE(0.0), dE(0.0);
          //cout << "Number?:  " << torsions[i].paramPtr->n << std::endl;
          for (size_t j(0U); j < torsion.p.number; ++j)
          {
            double F = torsion.p.force[j] * cparams.torsionunit();
            size_t k = torsion.p.order[j];
            double l = std::abs(torsion.p.ideal[j]) > 0.0 ? -1.0 : 1.0;
            tE += F * (1.0 + static_cast<double>(cos[k] * l));
            dE += -static_cast<double>(k * F * sin[k] * l);
          }
          E += tE;

          auto dt = cross(t, b12);
          auto du = cross(u, b12);
          dt *= dE / (tl2*r12);
          du *= -dE / (ul2*r12);

          //part_grad[TORSION][torsion.atoms[0]] += cross(dt, b12);
          //part_grad[TORSION][torsion.atoms[1]] += cross(b02, dt) + cross(du, b23);
          //part_grad[TORSION][torsion.atoms[2]] += cross(dt, b01) + cross(b13, du);
          //part_grad[TORSION][torsion.atoms[3]] += cross(du, b12);
          vir1 = cross(dt, b12);
          vir2 = cross(b02, dt) + cross(du, b23);
          vir3 = cross(dt, b01) + cross(b13, du);
          vir4 = cross(du, b12);

          part_grad[TORSION][torsion.atoms[0]] += vir1;
          part_grad[TORSION][torsion.atoms[1]] += vir2;
          part_grad[TORSION][torsion.atoms[2]] += vir3;
          part_grad[TORSION][torsion.atoms[3]] += vir4;

          //increment internal virial tensor
          vxx = b12.x()*(vir3.x() + vir4.x()) - b01.x()*vir1.x() + b23.x()*vir4.x();
          vyx = b12.y()*(vir3.x() + vir4.x()) - b01.y()*vir1.x() + b23.y()*vir4.x();
          vzx = b12.z()*(vir3.x() + vir4.x()) - b01.z()*vir1.x() + b23.z()*vir4.x();
          vyy = b12.y()*(vir3.y() + vir4.y()) - b01.y()*vir1.y() + b23.y()*vir4.y();
          vzy = b12.z()*(vir3.y() + vir4.y()) - b01.z()*vir1.y() + b23.z()*vir4.y();
          vzz = b12.z()*(vir3.z() + vir4.z()) - b01.z()*vir1.z() + b23.z()*vir4.z();
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






      //template<>
      //double energy::interfaces::amoeba::amoeba_ff::f_strbnd<0>(void)
      //{
      //	double eba, e, dr1, dr2;
      //	double force1, force2;
      //	double xia, yia, zia, xib, yib, zib, xic, yic, zic, xab, yab, zab, xcb, ycb, zcb, xp, yp, zp, dot, rp;
      //	size_t i, ia, ib, ic;
      //	double dt;
      //	double radian = 180 / PI, angleunit = (PI / 180)*(PI / 180), cosi;
      //	double stbnunit = 1.74532925199432955E-002;
      //	double term1, term1t, term2t;
      //	double term2, termr;
      //	double rab2, rcb2, rab, rcb;
      //	double energy = 0.0;
      //	eba = 0.0;
      //	const scon::nv3d &positions = coords->xyz();
      //	scon::v3d bv, b, gv;
      //
      //	for (auto strbnds : refined.strbends())
      //
      //	{
      //
      //
      //		force1 = strbnds.force;
      //		force2 = strbnds.ff;
      //		/*i = */
      //		//i = isb1[n + 1];
      //		i = strbnds.atoms[0];
      //		if (i == 0)continue;
      //		ia = strbnds.atoms[0];
      //		ib = strbnds.atoms[1];
      //		ic = strbnds.atoms[2];
      //
      //		xia = positions[ia - 1].x();
      //		yia = positions[ia - 1].y();
      //		zia = positions[ia - 1].z();
      //
      //		xib = positions[ib - 1].x();
      //		yib = positions[ib - 1].y();
      //		zib = positions[ib - 1].z();
      //
      //		xic = positions[ic - 1].x();
      //		yic = positions[ic - 1].y();
      //		zic = positions[ic - 1].z();
      //
      //		xab = xia - xib;
      //		yab = yia - yib;
      //		zab = zia - zib;
      //		xcb = xic - xib;
      //		ycb = yic - yib;
      //		zcb = zic - zib;
      //
      //		rab2 = xab*xab + yab*yab + zab*zab;
      //		rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
      //		if (rab2 != 0.0 && rcb2 != 0.0){
      //			rab = sqrt(rab2);
      //			rcb = sqrt(rcb2);
      //
      //			xp = ycb*zab - zcb*yab;
      //			yp = zcb*xab - xcb*zab;
      //			zp = xcb*yab - ycb*xab;
      //
      //			rp = sqrt(xp*xp + yp*yp + zp*zp);
      //			rp = max(rp, 0.001);
      //			dot = xab*xcb + yab*ycb + zab*zcb;
      //			cosi = dot / (rab*rcb);
      //			cosi = acos(cosi) * 180 / PI;
      //
      //			//dt = cosi - anat[i - 1];
      //			dt = cosi - 
      //
      //			term1 = -radian / (rab2*rp);
      //			term2 = radian / (rcb2*rp);
      //
      //	/*		ddtdxia = term1 * (yab*zp - zab*yp);
      //			ddtdyia = term1 * (zab*xp - xab*zp);
      //			ddtdzia = term1 * (xab*yp - yab*xp);
      //
      //			ddtdxic = term2 * (ycb*zp - zcb*yp);
      //			ddtdyic = term2 * (zcb*xp - xcb*zp);
      //			ddtdzic = term2 * (xcb*yp - ycb*xp);*/
      //
      //
      //
      //			dr1 = rab - bonds[isb2[i] - 1].relValue;
      //
      //
      //
      //			term1 = 1.0 / rab;
      //			dr2 = rcb - bonds[isb3[i] - 1].relValue;
      //			term2 = 1.0 / rcb;
      //
      //			/*ddrdxia = term1 * xab;
      //			ddrdyia = term1 * yab;
      //			ddrdzia = term1 * zab;
      //
      //			ddrdxic = term2 * xcb;
      //			ddrdyic = term2 * ycb;
      //			ddrdzic = term2 * zcb;*/
      //
      //
      //
      //
      //			term1 = force1*stbnunit;
      //			term2 = force2*stbnunit;
      //			// 	cout << term1 << endl;
      //			termr = term1*dr1 + term2*dr2;
      //			term1t = term1 * dt;
      //			term2t = term2 * dt;
      //
      //			e = termr*dt;
      //
      //
      //			energy += e;
      //		}
      //	}
      //}


      template<>
      double energy::interfaces::amoeba::amoeba_ff::f_oop<0>(void)
      {

        double eopb;
        size_t ia, ib, ic, id;
        double force;
        double xia, yia, zia, xib, yib, zib, xic, yic, zic, xid, yid, zid;
        double xab, yab, zab, xcb, ycb, zcb, xdb, ydb, zdb, xad, yad, zad, xcd, ycd, zcd;
        double ee, e, dt, dt2, dt3, dt4, rdb2, bkk2, cc, rad2, rcd2, dot, cosi;
        double energytemp(0.0);

        auto const &positions = coords->xyz();
        coords::Cartesian_Point bv, b;


        double E(0.0);
        for (auto opbends : refined.opbends())

        {


          eopb = 0.0;


          ia = opbends.atoms[0];
          ib = opbends.atoms[1];
          ic = opbends.atoms[2];
          id = opbends.atoms[3];

          if (id == 0) return energytemp = 0.0;
          force = opbends.force;


          xia = positions[ia].x();
          yia = positions[ia].y();
          zia = positions[ia].z();
          xib = positions[ib].x();
          yib = positions[ib].y();
          zib = positions[ib].z();

          xic = positions[ic].x();
          yic = positions[ic].y();
          zic = positions[ic].z();

          xid = positions[id].x();
          yid = positions[id].y();
          zid = positions[id].z();


          xab = xia - xib;
          yab = yia - yib;
          zab = zia - zib;
          xcb = xic - xib;
          ycb = yic - yib;
          zcb = zic - zib;
          xdb = xid - xib;
          ydb = yid - yib;
          zdb = zid - zib;
          xad = xia - xid;
          yad = yia - yid;
          zad = zia - zid;
          xcd = xic - xid;
          ycd = yic - yid;
          zcd = zic - zid;


          rad2 = xad*xad + yad*yad + zad*zad;
          rcd2 = xcd*xcd + ycd*ycd + zcd*zcd;
          dot = xad*xcd + yad*ycd + zad*zcd;
          cc = rad2*rcd2 - dot*dot;


          ee = xdb*(yab*zcb - zab*ycb) + ydb*(zab*xcb - xab*zcb) + zdb*(xab*ycb - yab*xcb);
          rdb2 = xdb*xdb + ydb*ydb + zdb*zdb;
          if (rdb2 != 0.0 && cc != 0.0) {
            bkk2 = rdb2 - ee*ee / cc;
            cosi = sqrt(bkk2 / rdb2);
            cosi = acos(cosi) * 180 / SCON_PI;
            dt = cosi;

            dt2 = dt*dt;
            dt3 = dt2 * dt;
            dt4 = dt2 * dt2;
            e = 3.04617419E-004*force*dt2*(1.0 - 0.014*dt + 0.000056*dt2 - 0.0000007*dt3 + 0.000000022*dt4);

            eopb += e;
            energytemp = eopb;


          }



        }





        return E;
      }

      template<>
      double energy::interfaces::amoeba::amoeba_ff::f_oop<1>(void)
      {

        double eopb;
        size_t ia, ib, ic, id;
        double force;
        double xia, yia, zia, xib, yib, zib, xic, yic, zic, xid, yid, zid;
        double xab, yab, zab, xcb, ycb, zcb, xdb, ydb, zdb, xad, yad, zad, xcd, ycd, zcd;
        double ee, e, dt, dt2, dt3, dt4, rdb2, bkk2, cc, rad2, rcd2, dot, cosi;

        double energytemp(0.0);

        double deddt, dedcos, deedxia, deedyia, deedzia, deedxic, deedyic, deedzic, deedxid, deedyid, deedzid, dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, dedzic, dedxid, dedyid, dedzid;
        double term, dccdxia, dccdyia, dccdzia, dccdxic, dccdyic, dccdzic, dccdxid, dccdyid, dccdzid;
        auto const & positions = coords->xyz();
        coords::Cartesian_Point bv, b, gv;
        double radian = 180 / SCON_PI;

        double E(0.0);
        for (auto opbends : refined.opbends())

        {




          eopb = 0.0;


          ia = opbends.atoms[0];
          ib = opbends.atoms[1];
          ic = opbends.atoms[2];
          id = opbends.atoms[3];
          if (id == 0) return energytemp = 0.0;
          force = opbends.force;


          xia = positions[ia].x();
          yia = positions[ia].y();
          zia = positions[ia].z();
          //std::cout << ia << "   " << ib << "  " << ic << "  " << id << '\n';
          xib = positions[ib].x();
          yib = positions[ib].y();
          zib = positions[ib].z();

          xic = positions[ic].x();
          yic = positions[ic].y();
          zic = positions[ic].z();

          xid = positions[id].x();
          yid = positions[id].y();
          zid = positions[id].z();


          xab = xia - xib;
          yab = yia - yib;
          zab = zia - zib;
          xcb = xic - xib;
          ycb = yic - yib;
          zcb = zic - zib;
          xdb = xid - xib;
          ydb = yid - yib;
          zdb = zid - zib;
          xad = xia - xid;
          yad = yia - yid;
          zad = zia - zid;
          xcd = xic - xid;
          ycd = yic - yid;
          zcd = zic - zid;


          rad2 = xad*xad + yad*yad + zad*zad;
          rcd2 = xcd*xcd + ycd*ycd + zcd*zcd;
          dot = xad*xcd + yad*ycd + zad*zcd;
          cc = rad2*rcd2 - dot*dot;




          ee = xdb*(yab*zcb - zab*ycb) + ydb*(zab*xcb - xab*zcb) + zdb*(xab*ycb - yab*xcb);
          rdb2 = xdb*xdb + ydb*ydb + zdb*zdb;
          if (rdb2 != 0.0 && cc != 0.0) {
            bkk2 = rdb2 - ee*ee / cc;
            cosi = sqrt(bkk2 / rdb2);
            cosi = acos(cosi) * 180 / SCON_PI;
            dt = cosi;

            dt2 = dt*dt;
            dt3 = dt2 * dt;
            dt4 = dt2 * dt2;
            e = 3.04617419E-004*force*dt2*(1.0 - 0.014*dt + 0.000056*dt2 - 0.0000007*dt3 + 0.000000022*dt4);
            deddt = 3.04617419E-004*force*dt*radian*(2.0 + 3.0*-0.014*dt + 4.0*0.000056*dt2 + 5.0*-0.0000007*dt3 + 6.0*0.000000022*dt4);
            double sign = 1.0;
            if (ee <= 0.0) sign = -1.0;
            dedcos = -deddt * sign / sqrt(cc*bkk2);

            term = ee / cc;
            dccdxia = (xad*rcd2 - xcd*dot) * term;
            dccdyia = (yad*rcd2 - ycd*dot) * term;
            dccdzia = (zad*rcd2 - zcd*dot) * term;
            dccdxic = (xcd*rad2 - xad*dot) * term;
            dccdyic = (ycd*rad2 - yad*dot) * term;
            dccdzic = (zcd*rad2 - zad*dot) * term;
            dccdxid = -dccdxia - dccdxic;
            dccdyid = -dccdyia - dccdyic;
            dccdzid = -dccdzia - dccdzic;

            term = ee / rdb2;
            deedxia = ydb*zcb - zdb*ycb;
            deedyia = zdb*xcb - xdb*zcb;
            deedzia = xdb*ycb - ydb*xcb;
            deedxic = yab*zdb - zab*ydb;
            deedyic = zab*xdb - xab*zdb;
            deedzic = xab*ydb - yab*xdb;
            deedxid = ycb*zab - zcb*yab + xdb*term;
            deedyid = zcb*xab - xcb*zab + ydb*term;
            deedzid = xcb*yab - ycb*xab + zdb*term;

            dedxia = dedcos * (dccdxia + deedxia);
            dedyia = dedcos * (dccdyia + deedyia);
            dedzia = dedcos * (dccdzia + deedzia);
            dedxic = dedcos * (dccdxic + deedxic);
            dedyic = dedcos * (dccdyic + deedyic);
            dedzic = dedcos * (dccdzic + deedzic);
            dedxid = dedcos * (dccdxid + deedxid);
            dedyid = dedcos * (dccdyid + deedyid);
            dedzid = dedcos * (dccdzid + deedzid);
            dedxib = -dedxia - dedxic - dedxid;
            dedyib = -dedyia - dedyic - dedyid;
            dedzib = -dedzia - dedzic - dedzid;

            eopb += e;
            energytemp = eopb;

            gv.x() = dedxia;
            gv.y() = dedyia;
            gv.z() = dedzia;

            part_grad[OPBEND][ia] += gv;

            gv.x() = dedxic;
            gv.y() = dedyic;
            gv.z() = dedzic;

            part_grad[OPBEND][ic] += gv;

            gv.x() = dedxid;
            gv.y() = dedyid;
            gv.z() = dedzid;

            part_grad[OPBEND][id] += gv;

            gv.x() = dedxib;
            gv.y() = dedyib;
            gv.z() = dedzib;

            part_grad[OPBEND][ib] += gv;


          }



        }





        return E;

      }


      /********************************
      *                                *
      *                        *
      *  VDW         *
      *  Energy/Gradients/Hessians     *
      *                                *
      *                                *
      *********************************/




      energy::interfaces::amoeba::nb_cutoff::nb_cutoff(double const ic, double const is)
        : c(ic), s(is), cc(c*c), ss(3.0*s*s), cs((cc - s*s)*(cc - s*s)*(cc - s*s)) { }

      bool energy::interfaces::amoeba::nb_cutoff::factors(double const rr, double & r, double & fQ, double & fV)
      {
        r = sqrt(rr);
        if (r > c) return false;
        double cr = cc - rr;
        fV = r < s ? 1.0 : (cr*cr*(cc + 2.0*rr - ss)) / cs;
        fQ = (1.0 - rr / cc);
        fQ *= fQ;
        return (std::abs(r) > 0.0);
      }


      void energy::interfaces::amoeba::amoeba_ff::boundary(double & x, double & y, double & z) const
      {
        static coords::Cartesian_Point const halfbox(Config::get().energy.pb_box / 2.0);

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





      /*
      buffered 14-7 potentials
      */

      template<> inline double energy::interfaces::amoeba::amoeba_ff::eV
        < ::tinker::parameter::radius_types::R_MIN>
        (double const E, double const R, double const r) const
      {
        double  T, T7, rho, V, rv7, R6, R7, rv;

        rv = R;
        rv7 = rv*rv*rv*rv*rv*rv*rv;
        R6 = r*r*r*r*r*r;
        R7 = R6*r;
        rho = R7 + 0.12 * rv7;
        T = 1.07 / (r + 0.07*rv);
        T7 = T*T*T*T*T*T*T;
        V = E * T7 * rv7* (1.12*rv7 / rho - 2.0);

        return V;


      }

      template<> inline double energy::interfaces::amoeba::amoeba_ff::eV
        < ::tinker::parameter::radius_types::SIGMA>
        (double const E, double const R, double const r) const
      {
        double  T, T7, rho, V, rv7, R6, R7, rv;

        rv = R;
        rv7 = rv*rv*rv*rv*rv*rv*rv;
        R6 = r*r*r*r*r*r;
        R7 = R6*r;
        rho = R7 + 0.12 * rv7;
        T = 1.07 / (r + 0.07*rv);
        T7 = T*T*T*T*T*T*T;
        V = E * T7 * rv7* (1.12*rv7 / rho - 1.0);

        return V;
      }

      template<> inline double energy::interfaces::amoeba::amoeba_ff::gV
        < ::tinker::parameter::radius_types::R_MIN>
        (double const E, double const R, double const r, double &dV) const
      {

        double  T, T7, rho, V, rv7, gt, dt, R6, R7, rv;

        rv = R;
        rv7 = rv*rv*rv*rv*rv*rv*rv;
        R6 = r*r*r*r*r*r;
        R7 = R6*r;
        rho = R7 + 0.12 * rv7;
        T = 1.07 / (r + 0.07*rv);
        T7 = T*T*T*T*T*T*T;
        V = E * T7 * rv7* (1.12*rv7 / rho - 2.0);
        dt = T / (1.07);
        gt = E * T7 *R6 * 1.12 * (rv7 / rho) * (rv7 / rho);
        dV = -7.0 * (dt*V + gt);
        dV /= r;
        return V;

      }

      template<> inline double energy::interfaces::amoeba::amoeba_ff::gV
        < ::tinker::parameter::radius_types::SIGMA>
        (double const E, double const R, double const r, double &dV) const
      {
        double  T, T7, rho, V, rv7, gt, dt, R6, R7, rv;

        rv = R;
        rv7 = rv*rv*rv*rv*rv*rv*rv;
        R6 = r*r*r*r*r*r;
        R7 = R6*r;
        rho = R7 + 0.12 * rv7;
        T = 1.07 / (r + 0.07*rv);
        T7 = T*T*T*T*T*T*T;
        V = E * T7 * rv7* (1.12*rv7 / rho - 1.0);
        dt = T / (1.07);
        gt = E * T7 *R6 * 1.12 * (rv7 / rho) * (rv7 / rho);
        dV = -7.0 * (dt*V + gt);
        dV /= r;

        return V;
      }


      template<> inline double energy::interfaces::amoeba::amoeba_ff::gV_fep
        < ::tinker::parameter::radius_types::R_MIN>
        (double const E, double const R, double const r, double const vout, double &dV) const
      {
        double D6, D12, D13, T2;
        double T = R, D = r*r;
        T = T*T*T; // T^3
        T = T*T; // T^6
        T2 = T*T; // T^12
        D6 = D*D*D; //r^6 / D^3
        D13 = D6*D6*r; //r^13
        D6 = Config::get().fep.ljshift * (1 - vout) * (1 - vout) * T + D6; //r^6 shifted
        D12 = D6 * D6; //r^12 shifted
        double V = vout * E * (T2 / D12 - 2 * T / D6);
        dV = vout * E * 12.0 * (T * D6 - (T2)) / D13;
        dV = dV / std::pow(D6, 0.16666666666666);
        return V;
      }

      template<> inline double energy::interfaces::amoeba::amoeba_ff::gV_fep
        < ::tinker::parameter::radius_types::SIGMA>
        (double const E, double const R, double const r, double const vout, double &dV) const
      {
        double T = R, D = r*r;
        T = T*T*T; // T^3
        T = T*T; // T^6
        D = D*D*D; // r^6
        double const vd = 0;
        D = Config::get().fep.ljshift * vd * vd * T + D; // r^6 shifted
        T /= D;
        double V = vout*E*T;
        dV = V*r*(6.0 - 12.0*T);
        return V*(T - 1.0);
      }

      template< ::tinker::parameter::radius_types::T RT>
      inline void energy::interfaces::amoeba::amoeba_ff::e_QV
        (double const C, double const E, double const R, double const d,
          double &e_c, double &e_v) const
      {
        //e_c += eQ(C, d);
        e_v += eV<RT>(E, R, d);
      }


      template< ::tinker::parameter::radius_types::T RT>
      inline void energy::interfaces::amoeba::amoeba_ff::g_QV
        (double const E, double const R, double const d,
          double &e_v, double &dE) const
      {
        double  dV(0.0);
        //e_c += gQ(C, d, dQ);
        e_v += gV<RT>(E, R, d, dV);
        //dE = (dQ + dV)*d;
        dE = dV;
      }

      template< ::tinker::parameter::radius_types::T RT>
      inline void energy::interfaces::amoeba::amoeba_ff::g_QV_fep
        (double const E, double const R, double const d,
          double const v_io,
          double &e_v, double &dE) const
      {
        double  dV(0.0);
        //e_c += gQ_fep(C, d, c_io, dQ);
        e_v += gV_fep<RT>(E, R, d, v_io, dV);
        dE = dV;
      }

      template< ::tinker::parameter::radius_types::T RT>
      inline void energy::interfaces::amoeba::amoeba_ff::e_QV_cutoff
        (double const C, double const E, double const R, double const d,
          double const fQ, double const fV,
          double &e_c, double &e_v) const
      {
        //e_c += eQ(C, d)*fQ;
        e_v += eV<RT>(E, R, d)*fV;
      }


      template< ::tinker::parameter::radius_types::T RT>
      inline void energy::interfaces::amoeba::amoeba_ff::g_QV_cutoff
        (double const E, double const R, double const d,
          double const fV,
          double &e_v, double &dE) const
      {
        double  dV(0.0);
        //e_c += gQ(C, d, dQ)*fQ;
        e_v += gV<RT>(E, R, d, dV)*fV;
        dE = dV*fV*d;
      }

      template< ::tinker::parameter::radius_types::T RT>
      inline void energy::interfaces::amoeba::amoeba_ff::g_QV_fep_cutoff
        (double const E, double const R, double const d,
          double const v_out, double const fV, double &e_v, double &dE) const
      {
        double/* dQ,*/ dV;
        //e_c += gQ_fep(C, d, c_out, dQ)*fQ;
        e_v += gV_fep<RT>(E, R, d, v_out, dV)*fV;
        dE = (/*dQ*fQ*/ dV*fV);
      }


      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::amoeba::amoeba_ff::g_nb(void)
      {


        part_energy[CHARGE] = 0.0;
        part_energy[VDW] = 0.0;
        part_energy[VDWC] = 0.0;
        part_grad[VDWC].assign(part_grad[VDWC].size(), coords::Cartesian_Point());

        coords->fep.feptemp = energy::fepvect();
        for (auto & ia : coords->interactions()) ia.energy = 0.0;

        for (auto const &pairmatrix : refined.pair_matrices())
        {

          size_t const N(coords->interactions().size());
          for (size_t sub_ia_index(0u), row(0u), col(0u); sub_ia_index < N; ++sub_ia_index)
          {
            coords::float_type & e(coords->interactions(sub_ia_index).energy);
            coords::Representation_3D & g(coords->interactions(sub_ia_index).grad);
            g.assign(coords->size(), coords::Cartesian_Point());
            std::vector< ::tinker::refine::types::nbpair> const & pl(pairmatrix.pair_matrix(sub_ia_index));
            scon::matrix< ::tinker::parameter::combi::vdwc, true> const & par(refined.vdwcm(pairmatrix.param_matrix_id));
            refine_vdw_h_bonds(pl, par);
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
              else
              {
                if (coords->atoms().in_exists() && (coords->atoms().sub_in() == row || coords->atoms().sub_in() == col))
                  g_nb_QV_pairs_fep_io<RT, false, false>(e, g, pl, par);
                else if (coords->atoms().out_exists() && (coords->atoms().sub_out() == row || coords->atoms().sub_out() == col))
                  g_nb_QV_pairs_fep_io<RT, false, true>(e, g, pl, par);
                else if (Config::get().energy.cutoff < 1000.0)
                  g_nb_QV_pairs_cutoff<RT, false>(e, g, pl, par);
                else
                  g_nb_QV_pairs<RT>(e, g, pl, par);
              } // end of no periodic else
            }
            else
            {
              if (Config::get().energy.periodic)
                g_nb_QV_pairs_cutoff<RT, true>(e, g, pl, par);
              else if (Config::get().energy.cutoff < 1000.0)
                g_nb_QV_pairs_cutoff<RT, false>(e, g, pl, par);
              else
                g_nb_QV_pairs<RT>(e, g, pl, par);
            }
            if (col == row)
            {
              col = 0;
              ++row;
            }
            else ++col;
            part_grad[VDW] += g;
          }
        }
        if (Config::get().md.fep)
        {
          // std::cout << coords->fep.feptemp.e_c_l2 << "  " << coords->fep.feptemp.e_vdw_l2 << "   " << coords->fep.feptemp.e_c_l1 << "   " << coords->fep.feptemp.e_vdw_l1 << std::endl;
            coords->fep.feptemp.dE = (coords->fep.feptemp.e_c_l2 + coords->fep.feptemp.e_vdw_l2) - (coords->fep.feptemp.e_c_l1 + coords->fep.feptemp.e_vdw_l1);
            coords->fep.feptemp.dG = 0;
            coords->fep.fepdata.push_back(coords->fep.feptemp);
        }
      }

      template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
      void energy::interfaces::amoeba::amoeba_ff::g_nb_QV_pairs_fep_switch_periodic
        (coords::float_type &e_nb, coords::Representation_3D &grad_vector,
          std::vector< ::tinker::refine::types::nbpair> const & pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const & parameters)
      {

      }

      //#ifndef _OPENMP

      template< ::tinker::parameter::radius_types::T RT>
      void energy::interfaces::amoeba::amoeba_ff::g_nb_QV_pairs
        (
          coords::float_type &e_nb, coords::Representation_3D &grad_vector,
          std::vector< ::tinker::refine::types::nbpair> const & pairlist,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params
          )
      {
        coords::float_type  e_v(0.0);
        for (auto const & pair : pairlist)
        {
          ::tinker::parameter::combi::vdwc const & p(params(refined.type(pair.a), refined.type(pair.b)));
          coords::Cartesian_Point b(vdwnew[pair.a] - vdwnew[pair.b]);



          coords::float_type r = len(b), dE(0.0);
          //std::cout << r << '\n';
          // std::cout << "Atoms: " << std::endl;
          // std::cout << pair.a + 1 << "  " << pair.b + 1 << std::endl;

          g_QV<RT>(p.E, p.R, r, e_v, dE);
          b *= dE;
          //charge_vector[pair.a] +=eQ(p.C,r)/r;
          //charge_vector[pair.b] +=eQ(p.C,r)/r;
          grad_vector[pair.a] += b;
          grad_vector[pair.b] -= b;
        }
        //std::cout << "e_nb=" << e_c + e_v << std::endl;
        e_nb += e_v;
        //part_energy[types::VDWC] += e_c+e_v;
        //part_energy[types::CHARGE] += e_c;
        part_energy[VDW] += e_v;
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC>
      void energy::interfaces::amoeba::amoeba_ff::g_nb_QV_pairs_cutoff
        (
          coords::float_type &e_nb, coords::Representation_3D &grad_vector,
          std::vector< ::tinker::refine::types::nbpair> const & pairlist,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params
          )
      {
        nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
        coords::float_type  e_v(0.0);
        //size_t ia(0u);
        for (auto const & pair : pairlist)
        {
          coords::Cartesian_Point b(coords->xyz(pair.a) - coords->xyz(pair.b));
          coords::Cartesian_Point dist;
          ::tinker::parameter::combi::vdwc const & p(params(refined.type(pair.a), refined.type(pair.b)));

          if (PERIODIC) boundary(b.x(), b.y(), b.z());
          coords::float_type const rr = dot(b, b);
          coords::float_type r(0.0), fQ(0.0), fV(0.0), dE(0.0);
          if (!cutob.factors(rr, r, fQ, fV))
          {
            continue;
          }
          //++ia;
          r = r;


          g_QV_cutoff<RT>(p.E, p.R, r, fV, e_v, dE);
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
        e_nb += e_v;
        //part_energy[types::VDWC] += e_c+e_v;
        //part_energy[types::CHARGE] += e_c;
        part_energy[VDW] += e_v;
      }

      template< ::tinker::parameter::radius_types::T RT, bool PERIODIC, bool ALCH_OUT>
      void energy::interfaces::amoeba::amoeba_ff::g_nb_QV_pairs_fep_io
        (
          double &e_nb, coords::Representation_3D &grad_vector,
          std::vector< ::tinker::refine::types::nbpair> const & pairlist,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params
          )
      {
        nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
        double e_c(0.0), e_v(0.0), e_c_l(0.0), e_vdw_l(0.0), e_c_dl(0.0), e_vdw_dl(0.0);
        double vxx, vyx, vzx, vyy, vzy, vzz;
        fepvar const & fep = coords->fep.window[coords->fep.window[0].step];
        for (auto const & pair : pairlist)
        {
          coords::Cartesian_Point b(coords->xyz(pair.a) - coords->xyz(pair.b));
          coords::Cartesian_Point dist;


          if (PERIODIC) boundary(b.x(), b.y(), b.z());
          double const rr(dot(b, b));
          double dE(0.0), Q(0.0), V(0.0);
          ::tinker::parameter::combi::vdwc const & p(params(refined.type(pair.a), refined.type(pair.b)));

          if (PERIODIC)
          {
            double fQ(0.0), fV(0.0);
            double r(0.0);
            if (!cutob.factors(rr, r, fQ, fV)) continue;
            g_QV_fep_cutoff<RT>(p.E, p.R, r, (ALCH_OUT ? fep.vin : fep.vout), fV, V, dE);
            double trash(0.0);
            g_QV_fep_cutoff<RT>(p.E, p.R, r, (ALCH_OUT ? fep.dvin : fep.dvout), fV, e_vdw_dl, trash);
          }
          else
          {
            if (Config::get().energy.cutoff < 1000.0) {
              double fV(0.0);
              double r(0.0);
              g_QV_fep_cutoff<RT>(p.E, p.R, r, (ALCH_OUT ? fep.vin : fep.vout), fV, V, dE);
              double trash(0.0);
              g_QV_fep_cutoff<RT>(p.E, p.R, r, (ALCH_OUT ? fep.dvin : fep.dvout), fV, e_vdw_dl, trash);
            }
            else {
              double const r = sqrt(rr);
              g_QV_fep<RT>(p.E, p.R, r, (ALCH_OUT ? fep.vin : fep.vout), V, dE);
              double trash(0.0);
              g_QV_fep<RT>(p.E, p.R, r, (ALCH_OUT ? fep.dvin : fep.dvout), e_vdw_dl, trash);
            }
          }
          dist = b;
          b *= dE;
          e_c_l += Q;
          e_vdw_l += V;
          e_c += Q;
          e_v += V;
          grad_vector[pair.a] += b;
          grad_vector[pair.b] -= b;
          //Increment internal virial tensor
          vxx = b.x() * dist.x();
          vyx = b.x() * dist.y();
          vzx = b.x() * dist.z();
          vyy = b.y() * dist.y();
          vzy = b.y() * dist.z();
          vzz = b.z() * dist.z();
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
        e_nb += e_c + e_v;
        coords->fep.feptemp.e_c_l1 += e_c_l;
        coords->fep.feptemp.e_c_l2 += e_c_dl;
        coords->fep.feptemp.e_vdw_l1 += e_vdw_l;
        coords->fep.feptemp.e_vdw_l2 += e_vdw_dl;
        //part_energy[types::VDWC] += e_c+e_v;
        //part_energy[types::CHARGE] += e_c;
        part_energy[VDW] += e_v;
      }

      //#else

      // #################### OPENMP  START ####################


      //double energy::interfaces::amoeba::amoeba_ff::f_14p(void)
      //{
      //	double start = omp_get_wtime();
      //	std::cout << "New almost parallel Torsion function!" << std::endl;
      //	std::cout << "Number of refined torsions:  " << refined.torsions().size() << std::endl;
      //	double E(0.0);
      //#pragma omp parallel 
      //	{
      //		double vxx, vyx, vzx, vyy, vzy, vzz;
      //		coords::Cartesian_Point vir1, vir2, vir3, vir4;
      //#pragma omp for reduction(+: E)
      //		//for (auto torsion : refined.torsions())
      //		for (int i = 0; i < refined.torsions().size(); i++)
      //		{
      //
      //			coords::Cartesian_Point b01(coords->xyz(refined.torsions()[i].atoms[1]) - coords->xyz(refined.torsions()[i].atoms[0]));
      //			coords::Cartesian_Point b12(coords->xyz(refined.torsions()[i].atoms[2]) - coords->xyz(refined.torsions()[i].atoms[1]));
      //			coords::Cartesian_Point b23(coords->xyz(refined.torsions()[i].atoms[3]) - coords->xyz(refined.torsions()[i].atoms[2]));
      //			coords::Cartesian_Point b02(coords->xyz(refined.torsions()[i].atoms[2]) - coords->xyz(refined.torsions()[i].atoms[0]));
      //			coords::Cartesian_Point b13(coords->xyz(refined.torsions()[i].atoms[3]) - coords->xyz(refined.torsions()[i].atoms[1]));
      //
      //			coords::Cartesian_Point t = b01.crossd(b12);
      //			coords::Cartesian_Point u = b12.crossd(b23);
      //			double tl2 = t.scalar(t);
      //			double ul2 = u.scalar(u);
      //			double tlul = sqrt(tl2*ul2);
      //			double r12 = b12.len();
      //			coords::Cartesian_Point tu = t.crossd(u);
      //
      //			double cos_scalar0 = t.scalar(u);
      //			double cos_scalar1 = tlul;
      //
      //			double sin_scalar0 = b12.scalar(tu);
      //			double sin_scalar1 = r12*tlul;
      //
      //			if (fabs(cos_scalar1) < 0.0000001 || fabs(sin_scalar1) < 0.0000001) continue;
      //			double cos[7], sin[7];
      //			cos[1] = cos_scalar0 / cos_scalar1;
      //			sin[1] = sin_scalar0 / sin_scalar1;
      //
      //			for (size_t j(2U); j <= refined.torsions()[i].p.max_order; ++j)
      //			{
      //				size_t k = j - 1;
      //				sin[j] = sin[k] * cos[1] + cos[k] * sin[1];
      //				cos[j] = cos[k] * cos[1] - sin[k] * sin[1];
      //			}
      //
      //			double tE(0.0), dE(0.0);
      //			//cout << "Number?:  " << torsions[i].paramPtr->n << std::endl;
      //			for (size_t j(0U); j < refined.torsions()[i].p.number; ++j)
      //			{
      //				double F = refined.torsions()[i].p.force[j] * cparams.torsionunit();
      //				size_t k = refined.torsions()[i].p.order[j];
      //				double l = fabs(refined.torsions()[i].p.ideal[j]) > 0.0 ? -1.0 : 1.0;
      //				tE += F * (1.0 + static_cast<double>(cos[k] * l));
      //				dE += -static_cast<double>(k * F * sin[k] * l);
      //			}
      //			E += tE;
      //			coords::Cartesian_Point dt(t);
      //			coords::Cartesian_Point du(u);
      //			dt.cross(b12);
      //			du.cross(b12);
      //			dt *= dE / (tl2*r12);
      //			du *= -dE / (ul2*r12);
      //
      //			vir1 = dt.crossd(b12);
      //			vir2 = b02.crossd(dt) + du.crossd(b23);
      //			vir3 = dt.crossd(b01) + b13.crossd(du);
      //			vir4 = du.crossd(b12);
      //
      //			part_grad[TORSION][refined.torsions()[i].atoms[0]] += vir1;
      //			part_grad[TORSION][refined.torsions()[i].atoms[1]] += vir2;
      //			part_grad[TORSION][refined.torsions()[i].atoms[2]] += vir3;
      //			part_grad[TORSION][refined.torsions()[i].atoms[3]] += vir4;
      //
      //			//increment internal virial tensor
      //			vxx = b12.x()*(vir3.x() + vir4.x()) - b01.x()*vir1.x() + b23.x()*vir4.x();
      //			vyx = b12.y()*(vir3.x() + vir4.x()) - b01.y()*vir1.x() + b23.y()*vir4.x();
      //			vzx = b12.z()*(vir3.x() + vir4.x()) - b01.z()*vir1.x() + b23.z()*vir4.x();
      //			vyy = b12.y()*(vir3.y() + vir4.y()) - b01.y()*vir1.y() + b23.y()*vir4.y();
      //			vzy = b12.z()*(vir3.y() + vir4.y()) - b01.z()*vir1.y() + b23.z()*vir4.y();
      //			vzz = b12.z()*(vir3.z() + vir4.z()) - b01.z()*vir1.z() + b23.z()*vir4.z();
      //			part_virial[TORSION][0][0] += vxx;
      //			part_virial[TORSION][1][0] += vyx;
      //			part_virial[TORSION][2][0] += vzx;
      //			part_virial[TORSION][0][1] += vyx;
      //			part_virial[TORSION][1][1] += vyy;
      //			part_virial[TORSION][2][1] += vzy;
      //			part_virial[TORSION][0][2] += vzx;
      //			part_virial[TORSION][1][2] += vzy;
      //			part_virial[TORSION][2][2] += vzz;
      //		}
      //	}//parallel region ends
      //	double end = omp_get_wtime();
      //	std::cout << "Time:  " << end - start << std::endl;
      //	return E;
      //}
      //
      //
      //
      //template< ::tinker::parameter::radius_types::T RT>
      //void energy::interfaces::amoeba::amoeba_ff::g_nb_QV_pairs
      //(
      //double &e_nb, coords::Representation_3D &grad_vector,
      //std::vector< ::tinker::refine::types::nbpair> const & pairlist,
      //scon::matrix< ::tinker::parameter::combi::vdwc> const & params
      //)
      //{
      //	ptrdiff_t const M(pairlist.size());
      //	double e_c(0.0), e_v(0.0);
      //#pragma omp parallel
      //	{
      //		coords::Representation_3D tmp_grad(grad_vector.size());
      //#pragma omp for reduction (+: e_c, e_v)
      //		for (ptrdiff_t i = 0; i<M; ++i)
      //		{
      //			coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));
      //			double r = 1.0 / std::sqrt(b.scalar(b)), dE(0.0);
      //			::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));
      //			g_QV<RT>(p.C, p.E, p.R, r, e_c, e_v, dE);
      //			b *= dE;
      //			tmp_grad[pairlist[i].a] += b;
      //			tmp_grad[pairlist[i].b] -= b;
      //		}
      //#pragma omp critical (nb_g_sum)
      //		{
      //			grad_vector += tmp_grad;
      //		}
      //	}
      //	e_nb += e_c + e_v;
      //
      //	//part_energy[types::VDWC] += e_c+e_v;
      //	part_energy[types::CHARGE] += e_c;
      //	part_energy[types::VDW] += e_v;
      //}
      //
      //template< ::tinker::parameter::radius_types::T RT, bool PERIODIC>
      //void energy::interfaces::amoeba::amoeba_ff::g_nb_QV_pairs_cutoff
      //(
      //double &e_nb, coords::Representation_3D &grad_vector,
      //std::vector< ::tinker::refine::types::nbpair> const & pairlist,
      //scon::matrix< ::tinker::parameter::combi::vdwc> const & params
      //)
      //{
      //	nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
      //	double e_c(0.0), e_v(0.0);
      //	ptrdiff_t const M(pairlist.size());
      //#pragma omp parallel
      //	{
      //		coords::Representation_3D tmp_grad(grad_vector.size());
      //		std::array<std::array<double, 3, 3> tempvir(coords::empty_virial());
      //#pragma omp for reduction (+: e_c, e_v)
      //		for (ptrdiff_t i = 0; i<M; ++i)
      //		{
      //			coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));
      //			coords::Cartesian_Point dist;
      //			double vxx, vyx, vzx, vyy, vzy, vzz;
      //			if (PERIODIC) boundary(b.x(), b.y(), b.z());
      //			double rr = b.scalar(b), r(0.0), fQ(0.0), fV(0.0), dE(0.0);
      //			if (!cutob.factors(rr, r, fQ, fV)) continue;
      //			r = 1.0 / r;
      //			::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a),
      //				refined.type(pairlist[i].b)));
      //			g_QV_cutoff<RT>(p.C, p.E, p.R, r, fQ, fV, e_c, e_v, dE);
      //			dist = b;
      //			b *= dE;
      //			tmp_grad[pairlist[i].a] += b;
      //			tmp_grad[pairlist[i].b] -= b;
      //			//Increment internal virial tensor
      //			vxx = b.x() * dist.x();
      //			vyx = b.x() * dist.y();
      //			vzx = b.x() * dist.z();
      //			vyy = b.y() * dist.y();
      //			vzy = b.y() * dist.z();
      //			vzz = b.z() * dist.z();
      //			tempvir[0][0] += vxx;
      //			tempvir[1][0] += vyx;
      //			tempvir[2][0] += vzx;
      //			tempvir[0][1] += vyx;
      //			tempvir[1][1] += vyy;
      //			tempvir[2][1] += vzy;
      //			tempvir[0][2] += vzx;
      //			tempvir[1][2] += vzy;
      //			tempvir[2][2] += vzz;
      //		}
      //#pragma omp critical (nb_g_sum)
      //		{
      //			grad_vector += tmp_grad;
      //			for (int i = 0; i <= 2; i++){
      //				for (int k = 0; k <= 2; k++){
      //					part_virial[VDWC][i][k] += tempvir[i][k];
      //				}
      //			}
      //		}
      //	}
      //	e_nb += e_c + e_v;
      //	//part_energy[types::VDWC] += e_c+e_v;
      //	part_energy[types::CHARGE] += e_c;
      //	part_energy[types::VDW] += e_v;
      //}
      //
      //template< ::tinker::parameter::radius_types::T RT, bool PERIODIC, bool ALCH_OUT>
      //void energy::interfaces::amoeba::amoeba_ff::g_nb_QV_pairs_fep_io
      //(
      //double &e_nb, coords::Representation_3D &grad_vector,
      //std::vector< ::tinker::refine::types::nbpair> const & pairlist,
      //scon::matrix< ::tinker::parameter::combi::vdwc> const & params
      //)
      //{
      //	nb_cutoff cutob(Config::get().energy.cutoff, Config::get().energy.switchdist);
      //	double e_c(0.0), e_v(0.0), e_c_l(0.0), e_vdw_l(0.0), e_c_dl(0.0), e_vdw_dl(0.0);
      //	fepvar const & fep = coords->fep.window[coords->fep.window[0].step];
      //	ptrdiff_t const M(pairlist.size());
      //#pragma omp parallel
      //	{
      //		coords::Representation_3D tmp_grad(grad_vector.size());
      //		std::array<std::array<double, 3, 3> tempvir(coords::empty_virial());
      //#pragma omp for reduction (+: e_c, e_v, e_c_l, e_c_dl, e_vdw_l, e_vdw_dl)
      //		for (ptrdiff_t i = 0; i<M; ++i)
      //		{
      //			coords::Cartesian_Point b(coords->xyz(pairlist[i].a) - coords->xyz(pairlist[i].b));
      //			if (PERIODIC) boundary(b.x(), b.y(), b.z());
      //			::tinker::parameter::combi::vdwc const & p(params(refined.type(pairlist[i].a), refined.type(pairlist[i].b)));
      //			double rr(b.scalar(b)), dE, Q(0.0), V(0.0);
      //			double vxx, vyx, vzx, vyy, vzy, vzz;
      //			coords::Cartesian_Point dist;
      //			if (PERIODIC)
      //			{
      //				double fQ(0.0), fV(0.0);
      //				double r(0.0);
      //				if (!cutob.factors(rr, r, fQ, fV)) continue;
      //				g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout), (ALCH_OUT ? fep.vin : fep.vout), fQ, fV, Q, V, dE);
      //				double trash(0.0);
      //				g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout), (ALCH_OUT ? fep.dvin : fep.dvout), fQ, fV, e_c_dl, e_vdw_dl, trash);
      //			}
      //			else
      //			{
      //				if (Config::get().energy.cutoff < 1000.0){
      //					double r(0.0);
      //					double fQ(0.0), fV(0.0);
      //					if (!cutob.factors(rr, r, fQ, fV)) continue;
      //					g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout), (ALCH_OUT ? fep.vin : fep.vout), fQ, fV, Q, V, dE);
      //					double trash(0.0);
      //					g_QV_fep_cutoff<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout), (ALCH_OUT ? fep.dvin : fep.dvout), fQ, fV, e_c_dl, e_vdw_dl, trash);
      //				}
      //				else{
      //					double const r = sqrt(rr);
      //					g_QV_fep<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.ein : fep.eout), (ALCH_OUT ? fep.vin : fep.vout), Q, V, dE);
      //					double trash(0.0);
      //					g_QV_fep<RT>(p.C, p.E, p.R, r, (ALCH_OUT ? fep.dein : fep.deout), (ALCH_OUT ? fep.dvin : fep.dvout), e_c_dl, e_vdw_dl, trash);
      //				}
      //			}
      //			dist = b;
      //			b *= dE;
      //			e_c_l += Q;
      //			e_vdw_l += V;
      //			e_c += Q;
      //			e_v += V;
      //			tmp_grad[pairlist[i].a] += b;
      //			tmp_grad[pairlist[i].b] -= b;
      //			//Increment internal virial tensor
      //			vxx = b.x() * dist.x();
      //			vyx = b.x() * dist.y();
      //			vzx = b.x() * dist.z();
      //			vyy = b.y() * dist.y();
      //			vzy = b.y() * dist.z();
      //			vzz = b.z() * dist.z();
      //			tempvir[0][0] += vxx;
      //			tempvir[1][0] += vyx;
      //			tempvir[2][0] += vzx;
      //			tempvir[0][1] += vyx;
      //			tempvir[1][1] += vyy;
      //			tempvir[2][1] += vzy;
      //			tempvir[0][2] += vzx;
      //			tempvir[1][2] += vzy;
      //			tempvir[2][2] += vzz;
      //		}
      //#pragma omp critical (nb_g_sum)
      //		{
      //			grad_vector += tmp_grad;
      //			for (int i = 0; i <= 2; i++){
      //				for (int k = 0; k <= 2; k++){
      //					part_virial[VDWC][i][k] += tempvir[i][k];
      //				}
      //			}
      //		}
      //	}
      //	e_nb += e_c + e_v;
      //	coords->fep.feptemp.e_c_l1 += e_c_l;
      //	coords->fep.feptemp.e_c_l2 += e_c_dl;
      //	coords->fep.feptemp.e_vdw_l1 += e_vdw_l;
      //	coords->fep.feptemp.e_vdw_l2 += e_vdw_dl;
      //	//part_energy[types::VDWC] += e_c+e_v;
      //	part_energy[types::CHARGE] += e_c;
      //	part_energy[types::VDW] += e_v;
      //}
      //
      //#endif


    }
  }
}


void energy::interfaces::amoeba::amoeba_ff::e_ind(void)
{
  std::vector <std::vector <double> > field, fieldp, udir, uold, uoldp, udirp;
  size_t i, k, j;
  size_t ii, kk;
  size_t iter, maxiter;
  double xr, yr, zr, r2, r;
  double ci, dix, diy, diz;
  double qixx, qixy, qixz;
  double qiyy, qiyz, qizz;
  double ck, dkx, dky, dkz;
  double qkxx, qkxy, qkxz;
  double qkyy, qkyz, qkzz;
  double duix, duiy, duiz;
  double dukx, duky, dukz;
  double pukx, puky, pukz;
  double puix, puiy, puiz;
  double duir, puir;
  double dukr, pukr;
  double scale3, scale5, scale7;
  double damp, expdamp;
  double rr3, rr5, rr7, dir;
  double qix, qiy, qiz, qir, dkr, qkx, qky, qkz, qkr;
  std::vector <double> fid, fkd, fip, fkp;
  bool done;
  double eps, epsold, epsd, epsp;
  double debye = 4.803210, poleps = 0.000001;
  size_t n, alloc(multipole_sites());
  pscale2 = 0.0; pscale3 = 0.0; pscale4 = 1.0; pscale5 = 1.0;
  mscale2 = 0.0; mscale3 = 0.0; mscale4 = 0.4; mscale5 = 0.8;
  dscale1 = 0.0; dscale2 = 1.0; dscale3 = 1.0; dscale4 = 1.0;
  uscale1 = 1.0; uscale2 = 1.0; uscale3 = 1.0; uscale4 = 1.0;
  p4scale = 1.0;
  double cutoff(Config::get().energy.cutoff);

  if (Config::get().energy.periodic == true){
  cutoff = len(Config::get().energy.pb_box);
  }
  else cutoff = Config::get().energy.cutoff;

  double pdi, pti, pgamma;


  //for (auto axes : refined.multipole_vecs())
  //{
  //	for (auto mult : axes)
  //	{
  //	


  //	}
  //}




  //    npole=npole-1;
  auto const & positions = coords->xyz();
  coords::Cartesian_Point bv, b, gv;
  uind.resize(5);
  uinp.resize(5);
  pscale.resize(alloc + 100);
  dscale.resize(alloc + 100);
  fid.resize(5);
  fkd.resize(5);
  fip.resize(5);
  fkp.resize(5);


  for (n = 0; n < uind.size(); n++) {
    uind[n].resize(alloc + 1);
  }
  for (n = 0; n < uinp.size(); n++) {
    uinp[n].resize(alloc + 1);
  }

  field.resize(5);
  fieldp.resize(5);
  udir.resize(5);
  udirp.resize(5);
  uold.resize(5);
  uoldp.resize(5);




  for (n = 0; n < field.size(); n++) {
    field[n].resize(alloc + 1);
  }
  for (n = 0; n < fieldp.size(); n++) {
    fieldp[n].resize(alloc + 1);
  }
  for (n = 0; n < udir.size(); n++) {
    udir[n].resize(alloc + 1);
  }
  for (n = 0; n < udirp.size(); n++) {
    udirp[n].resize(alloc + 1);
  }
  for (n = 0; n < uold.size(); n++) {
    uold[n].resize(alloc + 1);
  }
  for (n = 0; n < uoldp.size(); n++) {
    uoldp[n].resize(alloc + 1);
  }



  // ! zero out induced dipoles


  for (i = 1; i <= alloc; i++)
  {
    if (alloc == 0) break;
    for (j = 1; j <= 3; j++)
    {

      uind[j][i] = 0.0;
      uinp[j][i] = 0.0;
      field[j][i] = 0.0;
      fieldp[j][i] = 0.0;
    }
  }

  for (i = 1; i <= alloc - 1; i++) {

    if (alloc == 0) break;
    ii = ipole[i] + 1;
    if (ii == 0) continue;

    pdi = pdamp[i];

    pti = thole[i];

    ci = rp[0][i];
    dix = rp[1][i];

    diy = rp[2][i];
    diz = rp[3][i];
    qixx = rp[4][i];
    qixy = rp[5][i];
    qixz = rp[6][i];
    qiyy = rp[8][i];
    qiyz = rp[9][i];
    qizz = rp[12][i];


    for (j = i + 1; j <= alloc; j++) {

      pscale[ipole[j] + 1] = 1.0;
      dscale[ipole[j] + 1] = 1.0;

    }
    //std::cout << "n13.size" << n13[ii] << '\n';
    for (j = 1; j <= n13[ii]; j++) {
      pscale[i13[j][ii]] = pscale3;
      //std::cout << "i13 " << i13[j][ii] << " " << pscale3<< '\n';

    }
    for (j = 1; j <= n14[ii]; j++) {
      pscale[i14[j][ii]] = pscale4;
      for (k = 1; k <= np11[ii]; k++) {
        if (i14[j][ii] == ip11[k][ii]) pscale[i14[j][ii]] = pscale4*p4scale;

      }
    }

    for (j = 1; j <= n15[ii]; j++) {
      pscale[i15[j][ii]] = pscale5;

    }
    for (j = 1; j <= np11[ii]; j++) {
      dscale[ip11[j][ii]] = dscale1;

    }

    for (j = 1; j <= np12[ii]; j++) {
      dscale[ip12[j][ii]] = dscale2;
    }
    for (j = 1; j <= np13[ii]; j++) {
      dscale[ip13[j][ii]] = dscale3;
    }
    for (j = 1; j <= np14[ii]; j++) {
      dscale[ip14[j][ii]] = dscale4;
    }



    //std::cout << "TEST6_2\n";
    for (k = i + 1; k <= alloc; k++) {
      kk = ipole[k] + 1;
      if (kk == 0)continue;
      //std::cout << "TEST7\n";
      xr = positions[kk - 1].x() - positions[ii - 1].x();
      yr = positions[kk - 1].y() - positions[ii - 1].y();
      zr = positions[kk - 1].z() - positions[ii - 1].z();

      if( Config::get().energy.periodic == true){
      boundary(xr, yr, zr);
      }
      /*	std::cout << "TEST8\n";*/
      r2 = xr*xr + yr*yr + zr*zr;
      r = sqrt(r2);


      if (r < cutoff) {

        for (size_t ll = 0; ll < 13; ll++) {
          //if (abs(rp[ll][k]) == -NaN) rp[ll][k] = 0.0;
        }
        ck = rp[0][k];
        // 	cout << ck << endl;
        dkx = rp[1][k];

        dky = rp[2][k];
        dkz = rp[3][k];
        qkxx = rp[4][k];
        qkxy = rp[5][k];
        qkxz = rp[6][k];
        qkyy = rp[8][k];
        qkyz = rp[9][k];
        qkzz = rp[12][k];

        /*	     std::cout << ci << '\n';
        std::cout << diy << '\n';
        std::cout << dix << '\n';
        std::cout << diz << '\n';
        std::cout << qixy <<'\n';
        std::cout << qixx << '\n';*/

        scale3 = 1.0;
        scale5 = 1.0;
        scale7 = 1.0;

        damp = pdi * pdamp[k];

        /*	std::cout << "damp  " << damp << '\n';*/

        if (damp != 0.0) {
          pgamma = std::min(pti, thole[k]);

          damp = -pgamma * ((r / damp)*(r / damp)*(r / damp));

          if (damp > -50.0) {
            expdamp = exp(damp);
            scale3 = 1.0 - expdamp;
            scale5 = 1.0 - expdamp*(1.0 - damp);
            scale7 = 1.0 - expdamp*(1.0 - damp + 0.6*(damp*damp));
          }

        }

        rr3 = scale3 / (r*r2);

        rr5 = 3.0 * scale5 / (r*r2*r2);

        rr7 = 15.0 * scale7 / (r*r2*r2*r2);

        dir = dix*xr + diy*yr + diz*zr;

        qix = qixx*xr + qixy*yr + qixz*zr;
        qiy = qixy*xr + qiyy*yr + qiyz*zr;
        qiz = qixz*xr + qiyz*yr + qizz*zr;



        qir = qix*xr + qiy*yr + qiz*zr;
        dkr = dkx*xr + dky*yr + dkz*zr;
        qkx = qkxx*xr + qkxy*yr + qkxz*zr;
        qky = qkxy*xr + qkyy*yr + qkyz*zr;
        qkz = qkxz*xr + qkyz*yr + qkzz*zr;
        qkr = qkx*xr + qky*yr + qkz*zr;


        fid[0] = -xr*(rr3*ck - rr5*dkr + rr7*qkr) - rr3*dkx + 2.0*rr5*qkx;
        fid[1] = -yr*(rr3*ck - rr5*dkr + rr7*qkr) - rr3*dky + 2.0*rr5*qky;
        fid[2] = -zr*(rr3*ck - rr5*dkr + rr7*qkr) - rr3*dkz + 2.0*rr5*qkz;
        fkd[0] = xr*(rr3*ci + rr5*dir + rr7*qir) - rr3*dix - 2.0*rr5*qix;
        fkd[1] = yr*(rr3*ci + rr5*dir + rr7*qir) - rr3*diy - 2.0*rr5*qiy;
        fkd[2] = zr*(rr3*ci + rr5*dir + rr7*qir) - rr3*diz - 2.0*rr5*qiz;



        for (j = 1; j <= 3; j++) {
          field[j][i] = field[j][i] + fid[j - 1] * dscale[kk];
          field[j][k] = field[j][k] + fkd[j - 1] * dscale[kk];
          fieldp[j][i] = fieldp[j][i] + fid[j - 1] * pscale[kk];
          fieldp[j][k] = fieldp[j][k] + fkd[j - 1] * pscale[kk];
          //std::cout << dscale[kk] << '\n';

        }



      }



    }


  }

  for (i = 1; i <= alloc; i++) {
    if (alloc == 0) break;
    for (j = 1; j <= 3; j++) {
      udir[j][i] = polarity[i] * field[j][i];
      udirp[j][i] = polarity[i] * fieldp[j][i];
      uind[j][i] = udir[j][i];
      uinp[j][i] = udirp[j][i];


    }
  }

  done = false;
  maxiter = 10000;
  iter = 0;
  eps = 100.0;



  while (!done) {
    if (alloc == 0) break;

    // cout << "Hallo" << endl;
    for (i = 1; i <= alloc; i++) {
      for (j = 1; j <= 3; j++) {

        field[j][i] = 0.0;
        fieldp[j][i] = 0.0;
      }
    }
    for (i = 1; i <= alloc - 1; i++) {

      ii = ipole[i] + 1;
      if (ii == 0)continue;
      pdi = pdamp[i];
      pti = thole[i];
      duix = uind[1][i];
      duiy = uind[2][i];
      duiz = uind[3][i];
      puix = uinp[1][i];
      puiy = uinp[2][i];
      puiz = uinp[3][i];


      for (j = i + 1; j <= alloc; j++) {
        dscale[ipole[j - 1] + 1] = 1.0;
      }

      for (j = 1; j <= np11[ii]; j++) {
        dscale[ip11[j][ii]] = uscale1;
        //std::cout << ip11[j][ii] << "  " << uscale1 << '\n';
      }
      for (j = 1; j <= np12[ii]; j++) {
        dscale[ip12[j][ii]] = uscale2;
      }
      for (j = 1; j <= np13[ii]; j++) {
        dscale[ip13[j][ii]] = uscale3;
      }
      for (j = 1; j <= np14[ii]; j++) {
        dscale[ip14[j][ii]] = uscale4;
      }

      for (k = i + 1; k <= alloc; k++) {
        kk = ipole[k] + 1;


        if (kk == 0) continue;
        //              bool proceed =true;
        //         //      
        //             if(proceed) {
        xr = positions[kk - 1].x() - positions[ii - 1].x();
        yr = positions[kk - 1].y() - positions[ii - 1].y();
        zr = positions[kk - 1].z() - positions[ii - 1].z();

        /*if (configuration::e.periodic == true){
        boundary(xr, yr, zr);
        }
        */
        r2 = xr*xr + yr*yr + zr*zr;
        r = sqrt(r2);

        if (r < cutoff) {


          dukx = uind[1][k];
          duky = uind[2][k];
          dukz = uind[3][k];
          pukx = uinp[1][k];
          puky = uinp[2][k];
          pukz = uinp[3][k];


          scale3 = dscale[kk];
          scale5 = dscale[kk];

          damp = pdi * pdamp[k];
          if (damp != 0.0) {
            pgamma = std::min(pti, thole[k]);

            damp = -pgamma * ((r / damp)*(r / damp)*(r / damp));

            if (damp > -50.0) {
              expdamp = exp(damp);
              scale3 = scale3 * (1.0 - expdamp);
              scale5 = scale5 * (1.0 - expdamp*(1.0 - damp));

            }

          }

          rr3 = scale3 / (r*r2);

          rr5 = 3.0 * scale5 / (r*r2*r2);



          duir = xr*duix + yr*duiy + zr*duiz;
          dukr = xr*dukx + yr*duky + zr*dukz;
          puir = xr*puix + yr*puiy + zr*puiz;
          pukr = xr*pukx + yr*puky + zr*pukz;
          fid[0] = -rr3*dukx + rr5*dukr*xr;
          fid[1] = -rr3*duky + rr5*dukr*yr;
          fid[2] = -rr3*dukz + rr5*dukr*zr;
          fkd[0] = -rr3*duix + rr5*duir*xr;
          fkd[1] = -rr3*duiy + rr5*duir*yr;
          fkd[2] = -rr3*duiz + rr5*duir*zr;
          fip[0] = -rr3*pukx + rr5*pukr*xr;
          fip[1] = -rr3*puky + rr5*pukr*yr;
          fip[2] = -rr3*pukz + rr5*pukr*zr;
          fkp[0] = -rr3*puix + rr5*puir*xr;
          fkp[1] = -rr3*puiy + rr5*puir*yr;
          fkp[2] = -rr3*puiz + rr5*puir*zr;

          for (j = 1; j <= 3; j++) {

            field[j][i] = field[j][i] + fid[j - 1];
            field[j][k] = field[j][k] + fkd[j - 1];
            fieldp[j][i] = fieldp[j][i] + fip[j - 1];
            fieldp[j][k] = fieldp[j][k] + fkp[j - 1];

          }
        }









        //            }
      }
    }

    iter++;
    epsold = eps;
    epsd = 0.0;
    epsp = 0.0;

    for (i = 1; i <= alloc; i++) {
      if (alloc == 0) break;
      for (j = 1; j <= 3; j++) {
        uold[j][i] = uind[j][i];
        uoldp[j][i] = uinp[j][i];
        uind[j][i] = udir[j][i] + polarity[i] * field[j][i];
        uinp[j][i] = udirp[j][i] + polarity[i] * fieldp[j][i];
        uind[j][i] = uold[j][i] + 0.55*(uind[j][i] - uold[j][i]);
        uinp[j][i] = uoldp[j][i] + 0.55*(uinp[j][i] - uoldp[j][i]);
        epsd = epsd + (uind[j][i] - uold[j][i])*(uind[j][i] - uold[j][i]);
        epsp = epsp + (uinp[j][i] - uoldp[j][i])*(uinp[j][i] - uoldp[j][i]);




      }
    }
    eps = std::max(epsd, epsp);

    eps = debye * sqrt(eps / (alloc));

    // cout << npole << endl;
    if (eps < poleps) done = true;
    if (eps > epsold) done = true;
    if (iter >= maxiter) done = true;
  }
  if (eps > poleps) std::cout << "induced dipoles may not converged\n";

  //for (i = 1; i <= alloc; i++){
  //	if (alloc == 0) break;
  //	for (j = 1; j <= 3; j++){
  //		std::cout << field[j][i] << '\n';
  //		std::cout << fieldp[j][i] << '\n';
  //	}
  //}







  //


}
void energy::interfaces::amoeba::amoeba_ff::e_perm(void)
{


  size_t i, j, k;
  size_t ii, kk;
  size_t ix, iy, iz, kx, ky, kz;

  double e(0.0), ei(0.0), f(0.0);
  double damp(0.0), expdamp(0.0);
  double tempmulti(0.0), temppol(0.0);
  double pdi(0.0), pti(0.0), pgamma(0.0);
  double scale3(0.0), scale3i(0.0), scale7(0.0), scale7i(0.0), scale5(0.0), scale5i(0.0);
  double temp3(0.0), temp5(0.0), temp7(0.0);
  double psc3(0.0), psc5(0.0), psc7(0.0), dsc3(0.0), dsc5(0.0), dsc7(0.0);
  double xr(0.0), yr(0.0), zr(0.0);
  double r1(0.0), r2(0.0), rr1(0.0), rr(0.0), rr3(0.0), rr5(0.0), rr7(0.0), rr9(0.0), rr11(0.0);
  double ci(0.0), ck(0.0);
  double cutoff(0.0), dd(0.0), cc(0.0), fQ(0.0);


  if (Config::get().energy.periodic == true){
  cutoff = len(Config::get().energy.pb_box);
  }
  else  cutoff = Config::get().energy.cutoff;
  cc = cutoff *cutoff;

  std::vector <double> di(3), qi(9), dk(3), qk(9);
  std::vector <double> frcxi(4), frcxk(4), frcyi(4), frcyk(4), frczi(4), frczk(4);
  std::vector <double> fridmp(4), findmp(4), ftm2(4), ftm2i(4), ttm2(4), ttm3(4), ttm2i(4), ttm3i(4);
  std::vector <double> dixdk(3), dkxui(3), dixukp(3), dkxuip(3), uixqkr(3), ukxqir(3), uixqkrp(3), ukxqirp(3), qiuk(3), qkui(3), qiukp(3), qkuip(3), rxqiuk(3), rxqqkui(3), rxqiukp(3), rxkuip(3), qidk(3), qkdi(3), qir(3), qkr(3), qiqkr(3), qkqir(3);
  std::vector <double> qixqk(3), rxqir(3), dixr(3), dkxr(3), dixqkr(3), rxqkr(3), qkrxqir(3), rxqikr(3), ryqkir(3), rxqidk(3), rxqkdi(3), ddsc3(3), ddsc5(3), ddsc7(3);
  std::vector <double> dixuk(3), rxqkir(3), dkxqir(3), rxqkui(3), rxqkuip(3);
  std::vector <double> sc(11), sci(9), scip(9), gli(8), glip(8), gf(8), gfi(8), gti(8);
  std::vector <double> gl(9);
  double dielec = 332.063714000;
  dem.resize(5);
  dep.resize(5);
  const size_t N = alloc_glob;
  auto const &positions = coords->xyz();
  coords::Cartesian_Point bv, b, gv_multipole, gv_polarization;
  //scon::v3d &gradients = part_grad[energy::interfaces::amoeba::types::mpp];

  pscale.resize(N + 1);
  dscale.resize(N + 1);
  mscale.resize(N + 1);
  uscale.resize(N + 1);
  //  cout << "chekc8" << endl;
  f = dielec / 1.0;
  for (i = 0; i < dem.size(); i++) {
    dem[i].resize(N + 1);
  }
  for (i = 0; i < dep.size(); i++) {
    dep[i].resize(N + 1);
  }


  em = 0.0;
  ep = 0.0;
  //   npole=npole-1;


  for (i = 1; i <= N; i++) {
    for (j = 1; j <= 3; j++) {
      dem[j][i] = 0.0;
      dep[j][i] = 0.0;
    }
  }

  for (i = 1; i <= N; i++) {
    mscale[i] = 1.0;
    pscale[i] = 1.0;
    dscale[i] = 1.0;
    uscale[i] = 1.0;
  }
  //std::cout << "ANFANG  " << N << '\n';
  for (i = 1; i <= N - 1; i++) {
    if (N == 0) break;

    ii = ipole[i] + 1;
    if (ii == 0)continue;
    iz = zaxis[i];
    ix = xaxis[i];
    iy = yaxis[i];

    pdi = pdamp[i];
    pti = thole[i];
    ci = rp[0][i];
    di[0] = rp[1][i];
    di[1] = rp[2][i];
    di[2] = rp[3][i];
    qi[0] = rp[4][i];
    qi[1] = rp[5][i];
    qi[2] = rp[6][i];
    qi[3] = rp[7][i];
    qi[4] = rp[8][i];
    qi[5] = rp[9][i];
    qi[6] = rp[10][i];
    qi[7] = rp[11][i];
    qi[8] = rp[12][i];



    for (j = 0; j < coords->atoms(ii - 1).bonds().size(); j++) {
      pscale[coords->atoms(ii - 1).bonds()[j] + 1] = pscale2;
      mscale[coords->atoms(ii - 1).bonds()[j] + 1] = mscale2;

    }
    for (j = 1; j <= n13[ii]; j++) {
      pscale[i13[j][ii]] = pscale3;
      mscale[i13[j][ii]] = mscale3;

    }
    for (j = 1; j <= n14[ii]; j++) {
      pscale[i14[j][ii]] = pscale4;
      mscale[i14[j][ii]] = mscale4;

      for (k = 1; k <= np11[ii]; k++) {
        if (i14[j][ii] == ip11[k][ii]) pscale[i14[j][ii]] = pscale4*p4scale;

      }
    }

    for (j = 1; j <= n15[ii]; j++) {
      pscale[i15[j][ii]] = pscale5;
      mscale[i15[j][ii]] = mscale5;
    }
    for (j = 1; j <= np11[ii]; j++) {
      dscale[ip11[j][ii]] = dscale1;
      uscale[ip11[j][ii]] = uscale1;
    }
    for (j = 1; j <= np12[ii]; j++) {
      dscale[ip12[j][ii]] = dscale2;
      uscale[ip12[j][ii]] = uscale2;
    }
    for (j = 1; j <= np13[ii]; j++) {
      dscale[ip13[j][ii]] = dscale3;
      uscale[ip13[j][ii]] = uscale3;
    }
    for (j = 1; j <= np14[ii]; j++) {

      dscale[ip14[j][ii]] = dscale4;
      uscale[ip14[j][ii]] = uscale4;
    }
    //    cout << "chekc8" << endl;
    for (k = i + 1; k <= N; k++) {
      kk = ipole[k] + 1;
      if (kk == 0)continue;
      kz = zaxis[k];
      kx = xaxis[k];
      ky = yaxis[k];
      xr = positions[kk - 1].x() - positions[ii - 1].x();
      yr = positions[kk - 1].y() - positions[ii - 1].y();
      zr = positions[kk - 1].z() - positions[ii - 1].z();

      if (Config::get().energy.periodic == true){
      boundary(xr, yr, zr);
      }

      r2 = xr*xr + yr*yr + zr*zr;
      r1 = sqrt(r2);


      if (r1 < cutoff) {

        dd = r2;
        fQ = (1 - dd / cc);
        fQ *= fQ;
        //cout << cutoff << "   " << fQ << endl;

        ck = rp[0][k];
        dk[0] = rp[1][k];
        dk[1] = rp[2][k];
        dk[2] = rp[3][k];
        qk[0] = rp[4][k];
        qk[1] = rp[5][k];
        qk[2] = rp[6][k];
        qk[3] = rp[7][k];
        qk[4] = rp[8][k];
        qk[5] = rp[9][k];
        qk[6] = rp[10][k];
        qk[7] = rp[11][k];
        qk[8] = rp[12][k];


        rr1 = 1.0 / r1;
        rr3 = rr1 / r2;
        rr5 = 3.0 *rr3 / r2;
        rr7 = 5.0 *rr5 / r2;
        rr9 = 7.0 *rr7 / r2;
        rr11 = 9.0 *rr9 / r2;
        scale3 = 1.0;
        scale5 = 1.0;
        scale7 = 1.0;

        //}

        for (j = 1; j <= 3; j++) {
          ddsc3[j - 1] = 0.0;
          ddsc5[j - 1] = 0.0;
          ddsc7[j - 1] = 0.0;
        }

        damp = pdi * pdamp[k];

        if (damp != 0.0) {
          pgamma = std::min(pti, thole[k]);

          damp = -pgamma * ((r1 / damp)*(r1 / damp)*(r1 / damp));

          if (damp > -50.0) {
            expdamp = exp(damp);
            scale3 = 1.0 - expdamp;
            scale5 = 1.0 - (1.0 - damp)*expdamp;
            scale7 = 1.0 - (1.0 - damp + 0.6*(damp*damp))*expdamp;
            temp3 = -3.0 *damp *expdamp / r2;
            temp5 = -damp;
            temp7 = -0.2 - 0.6*damp;

            ddsc3[0] = temp3*xr;
            ddsc3[1] = temp3*yr;
            ddsc3[2] = temp3*zr;
            ddsc5[0] = temp5 * ddsc3[0];
            ddsc5[1] = temp5 * ddsc3[1];
            ddsc5[2] = temp5 * ddsc3[2];
            ddsc7[0] = temp7 * ddsc5[0];
            ddsc7[1] = temp7 * ddsc5[1];
            ddsc7[2] = temp7 * ddsc5[2];

          }
        }
        scale3i = scale3 *uscale[kk];
        scale5i = scale5 * uscale[kk];
        scale7i = scale7 * uscale[kk];

        dsc3 = scale3 *dscale[kk];
        dsc5 = scale5 *dscale[kk];
        dsc7 = scale7 *dscale[kk];

        psc3 = scale3 *pscale[kk];
        psc5 = scale5 *pscale[kk];
        psc7 = scale7 *pscale[kk];

        //!construction of necessary auxilliary vectors	

        dixdk[0] = di[1] * dk[2] - di[2] * dk[1];
        dixdk[1] = di[2] * dk[0] - di[0] * dk[2];
        dixdk[2] = di[0] * dk[1] - di[1] * dk[0];
        //std::cout << dixdk[0] << '\n';
        dixuk[0] = di[1] * uind[3][k] - di[2] * uind[2][k];
        dixuk[1] = di[2] * uind[1][k] - di[0] * uind[3][k];
        dixuk[2] = di[0] * uind[2][k] - di[1] * uind[1][k];

        dkxui[0] = dk[1] * uind[3][i] - dk[2] * uind[2][i];
        dkxui[1] = dk[2] * uind[1][i] - dk[0] * uind[3][i];
        dkxui[2] = dk[0] * uind[2][i] - dk[1] * uind[1][i];

        dixukp[0] = di[1] * uinp[3][k] - di[2] * uinp[2][k];
        dixukp[1] = di[2] * uinp[1][k] - di[0] * uinp[3][k];
        dixukp[2] = di[0] * uinp[2][k] - di[1] * uinp[1][k];

        dkxuip[0] = dk[1] * uinp[3][i] - dk[2] * uinp[2][i];
        dkxuip[1] = dk[2] * uinp[1][i] - dk[0] * uinp[3][i];
        dkxuip[2] = dk[0] * uinp[2][i] - dk[1] * uinp[1][i];
        /*std::cout << dixdk[0] << '\n';*/
        //

        dixr[0] = di[1] * zr - di[2] * yr;
        dixr[1] = di[2] * xr - di[0] * zr;
        dixr[2] = di[0] * yr - di[1] * xr;

        dkxr[0] = dk[1] * zr - dk[2] * yr;
        dkxr[1] = dk[2] * xr - dk[0] * zr;
        dkxr[2] = dk[0] * yr - dk[1] * xr;

        qir[0] = qi[0] * xr + qi[3] * yr + qi[6] * zr;
        qir[1] = qi[1] * xr + qi[4] * yr + qi[7] * zr;
        qir[2] = qi[2] * xr + qi[5] * yr + qi[8] * zr;

        qkr[0] = qk[0] * xr + qk[3] * yr + qk[6] * zr;
        qkr[1] = qk[1] * xr + qk[4] * yr + qk[7] * zr;
        qkr[2] = qk[2] * xr + qk[5] * yr + qk[8] * zr;

        qiqkr[0] = qi[0] * qkr[0] + qi[3] * qkr[1] + qi[6] * qkr[2];
        qiqkr[1] = qi[1] * qkr[0] + qi[4] * qkr[1] + qi[7] * qkr[2];
        qiqkr[2] = qi[2] * qkr[0] + qi[5] * qkr[1] + qi[8] * qkr[2];

        qkqir[0] = qk[0] * qir[0] + qk[3] * qir[1] + qk[6] * qir[2];
        qkqir[1] = qk[1] * qir[0] + qk[4] * qir[1] + qk[7] * qir[2];
        qkqir[2] = qk[2] * qir[0] + qk[5] * qir[1] + qk[8] * qir[2];

        qixqk[0] = qi[1] * qk[2] + qi[4] * qk[5] + qi[7] * qk[8] - qi[2] * qk[1] - qi[5] * qk[4] - qi[8] * qk[7];
        qixqk[1] = qi[2] * qk[0] + qi[5] * qk[3] + qi[8] * qk[6] - qi[0] * qk[2] - qi[3] * qk[5] - qi[6] * qk[8];
        qixqk[2] = qi[0] * qk[1] + qi[3] * qk[4] + qi[6] * qk[7] - qi[1] * qk[0] - qi[4] * qk[3] - qi[7] * qk[6];

        rxqir[0] = yr*qir[2] - zr*qir[1];
        rxqir[1] = zr*qir[0] - xr*qir[2];
        rxqir[2] = xr*qir[1] - yr*qir[0];

        rxqkr[0] = yr*qkr[2] - zr*qkr[1];
        rxqkr[1] = zr*qkr[0] - xr*qkr[2];
        rxqkr[2] = xr*qkr[1] - yr*qkr[0];

        rxqikr[0] = yr*qiqkr[2] - zr*qiqkr[1];
        rxqikr[1] = zr*qiqkr[0] - xr*qiqkr[2];
        rxqikr[2] = xr*qiqkr[1] - yr*qiqkr[0];

        rxqkir[0] = yr*qkqir[2] - zr*qkqir[1];
        rxqkir[1] = zr*qkqir[0] - xr*qkqir[2];
        rxqkir[2] = xr*qkqir[1] - yr*qkqir[0];

        qkrxqir[0] = qkr[1] * qir[2] - qkr[2] * qir[1];
        qkrxqir[1] = qkr[2] * qir[0] - qkr[0] * qir[2];
        qkrxqir[2] = qkr[0] * qir[1] - qkr[1] * qir[0];

        qidk[0] = qi[0] * dk[0] + qi[3] * dk[1] + qi[6] * dk[2];
        qidk[1] = qi[1] * dk[0] + qi[4] * dk[1] + qi[7] * dk[2];
        qidk[2] = qi[2] * dk[0] + qi[5] * dk[1] + qi[8] * dk[2];

        qkdi[0] = qk[0] * di[0] + qk[3] * di[1] + qk[6] * di[2];
        qkdi[1] = qk[1] * di[0] + qk[4] * di[1] + qk[7] * di[2];
        qkdi[2] = qk[2] * di[0] + qk[5] * di[1] + qk[8] * di[2];

        qiuk[0] = qi[0] * uind[1][k] + qi[3] * uind[2][k] + qi[6] * uind[3][k];
        qiuk[1] = qi[1] * uind[1][k] + qi[4] * uind[2][k] + qi[7] * uind[3][k];
        qiuk[2] = qi[2] * uind[1][k] + qi[5] * uind[2][k] + qi[8] * uind[3][k];

        qkui[0] = qk[0] * uind[1][i] + qk[3] * uind[2][i] + qk[6] * uind[3][i];
        qkui[1] = qk[1] * uind[1][i] + qk[4] * uind[2][i] + qk[7] * uind[3][i];
        qkui[2] = qk[2] * uind[1][i] + qk[5] * uind[2][i] + qk[8] * uind[3][i];

        qiukp[0] = qi[0] * uinp[1][k] + qi[3] * uinp[2][k] + qi[6] * uinp[3][k];
        qiukp[1] = qi[1] * uinp[1][k] + qi[4] * uinp[2][k] + qi[7] * uinp[3][k];
        qiukp[2] = qi[2] * uinp[1][k] + qi[5] * uinp[2][k] + qi[8] * uinp[3][k];

        qkuip[0] = qk[0] * uinp[1][i] + qk[3] * uinp[2][i] + qk[6] * uinp[3][i];
        qkuip[1] = qk[1] * uinp[1][i] + qk[4] * uinp[2][i] + qk[7] * uinp[3][i];
        qkuip[2] = qk[2] * uinp[1][i] + qk[5] * uinp[2][i] + qk[8] * uinp[3][i];
        // 	cout << uind[1][k] << endl;	

        dixqkr[0] = di[1] * qkr[2] - di[2] * qkr[1];
        dixqkr[1] = di[2] * qkr[0] - di[0] * qkr[2];
        dixqkr[2] = di[0] * qkr[1] - di[1] * qkr[0];

        dkxqir[0] = dk[1] * qir[2] - dk[2] * qir[1];
        dkxqir[1] = dk[2] * qir[0] - dk[0] * qir[2];
        dkxqir[2] = dk[0] * qir[1] - dk[1] * qir[0];

        uixqkr[0] = uind[2][i] * qkr[2] - uind[3][i] * qkr[1];
        uixqkr[1] = uind[3][i] * qkr[0] - uind[1][i] * qkr[2];
        uixqkr[2] = uind[1][i] * qkr[1] - uind[2][i] * qkr[0];

        ukxqir[0] = uind[2][k] * qir[2] - uind[3][k] * qir[1];
        ukxqir[1] = uind[3][k] * qir[0] - uind[1][k] * qir[2];
        ukxqir[2] = uind[1][k] * qir[1] - uind[2][k] * qir[0];

        uixqkrp[0] = uinp[2][i] * qkr[2] - uinp[3][i] * qkr[1];
        uixqkrp[1] = uinp[3][i] * qkr[0] - uinp[1][i] * qkr[2];
        uixqkrp[2] = uinp[1][i] * qkr[1] - uinp[2][i] * qkr[0];

        ukxqirp[0] = uinp[2][k] * qir[2] - uinp[3][k] * qir[1];
        ukxqirp[1] = uinp[3][k] * qir[0] - uinp[1][k] * qir[2];
        ukxqirp[2] = uinp[1][k] * qir[1] - uinp[2][k] * qir[0];

        rxqidk[0] = yr*qidk[2] - zr*qidk[1];
        rxqidk[1] = zr*qidk[0] - xr*qidk[2];
        rxqidk[2] = xr*qidk[1] - yr*qidk[0];

        rxqkdi[0] = yr*qkdi[2] - zr*qkdi[1];
        rxqkdi[1] = zr*qkdi[0] - xr*qkdi[2];
        rxqkdi[2] = xr*qkdi[1] - yr*qkdi[0];


        rxqiuk[0] = yr*qiuk[2] - zr*qiuk[1];
        rxqiuk[1] = zr*qiuk[0] - xr*qiuk[2];
        rxqiuk[2] = xr*qiuk[1] - yr*qiuk[0];

        rxqkui[0] = yr*qkui[2] - zr*qkui[1];
        rxqkui[1] = zr*qkui[0] - xr*qkui[2];
        rxqkui[2] = xr*qkui[1] - yr*qkui[0];

        rxqiukp[0] = yr*qiukp[2] - zr*qiukp[1];
        rxqiukp[1] = zr*qiukp[0] - xr*qiukp[2];
        rxqiukp[2] = xr*qiukp[1] - yr*qiukp[0];

        rxqkuip[0] = yr*qkuip[2] - zr*qkuip[1];
        rxqkuip[1] = zr*qkuip[0] - xr*qkuip[2];
        rxqkuip[2] = xr*qkuip[1] - yr*qkuip[0];

        //! calculation of scalarproducts for permanent components

        sc[2] = di[0] * dk[0] + di[1] * dk[1] + di[2] * dk[2];
        sc[3] = di[0] * xr + di[1] * yr + di[2] * zr;
        sc[4] = dk[0] * xr + dk[1] * yr + dk[2] * zr;
        sc[5] = qir[0] * xr + qir[1] * yr + qir[2] * zr;
        sc[6] = qkr[0] * xr + qkr[1] * yr + qkr[2] * zr;
        sc[7] = qir[0] * dk[0] + qir[1] * dk[1] + qir[2] * dk[2];
        sc[8] = qkr[0] * di[0] + qkr[1] * di[1] + qkr[2] * di[2];
        sc[9] = qir[0] * qkr[0] + qir[1] * qkr[1] + qir[2] * qkr[2];
        sc[10] = qi[0] * qk[0] + qi[1] * qk[1] + qi[2] * qk[2] + qi[3] * qk[3] + qi[4] * qk[4] + qi[5] * qk[5] + qi[6] * qk[6] + qi[7] * qk[7] + qi[8] * qk[8];

        //! calculation of the scalproducts for induced components

        sci[1] = uind[1][i] * dk[0] + uind[2][i] * dk[1] + uind[3][i] * dk[2] + di[0] * uind[1][k] + di[1] * uind[2][k] + di[2] * uind[3][k];
        sci[2] = uind[1][i] * uind[1][k] + uind[2][i] * uind[2][k] + uind[3][i] * uind[3][k];
        sci[3] = uind[1][i] * xr + uind[2][i] * yr + uind[3][i] * zr;
        sci[4] = uind[1][k] * xr + uind[2][k] * yr + uind[3][k] * zr;
        sci[7] = qir[0] * uind[1][k] + qir[1] * uind[2][k] + qir[2] * uind[3][k];
        sci[8] = qkr[0] * uind[1][i] + qkr[1] * uind[2][i] + qkr[2] * uind[3][i];

        scip[1] = uinp[1][i] * dk[0] + uinp[2][i] * dk[1] + uinp[3][i] * dk[2] + di[0] * uinp[1][k] + di[1] * uinp[2][k] + di[2] * uinp[3][k];
        scip[2] = uind[1][i] * uinp[1][k] + uind[2][i] * uinp[2][k] + uind[3][i] * uinp[3][k] + uinp[1][i] * uind[1][k] + uinp[2][i] * uind[2][k] + uinp[3][i] * uind[3][k];
        scip[3] = uinp[1][i] * xr + uinp[2][i] * yr + uinp[3][i] * zr;
        scip[4] = uinp[1][k] * xr + uinp[2][k] * yr + uinp[3][k] * zr;
        scip[7] = qir[0] * uinp[1][k] + qir[1] * uinp[2][k] + qir[2] * uinp[3][k];
        scip[8] = qkr[0] * uinp[1][i] + qkr[1] * uinp[2][i] + qkr[2] * uinp[3][i];

        //! calculation of the gl functions for permanent moments

        gl[0] = ci*ck;
        gl[1] = ck*sc[3] - ci*sc[4];
        gl[2] = ci*sc[6] + ck*sc[5] - sc[3] * sc[4];
        gl[3] = sc[3] * sc[6] - sc[4] * sc[5];
        gl[4] = sc[5] * sc[6];
        gl[5] = -4.0 * sc[9];
        gl[6] = sc[2];
        gl[7] = 2.0 * (sc[7] - sc[8]);
        gl[8] = 2.0 *sc[10];

        //! calculate the gl function for induced moments

        gli[1] = ck*sci[3] - ci*sci[4];
        gli[2] = -sc[3] * sci[4] - sci[3] * sc[4];
        gli[3] = sci[3] * sc[6] - sci[4] * sc[5];
        gli[6] = sci[1];
        gli[7] = 2.0 * (sci[7] - sci[8]);
        glip[1] = ck*scip[3] - ci*scip[4];
        glip[2] = -sc[3] * scip[4] - scip[3] * sc[4];
        glip[3] = scip[3] * sc[6] - scip[4] * sc[5];
        glip[6] = scip[1];
        glip[7] = 2.0 * (scip[7] - scip[8]);

        //! compute the energy contribution
        	if (Config::get().energy.periodic == true){
        e = rr1*gl[0] + rr3*(gl[1] + gl[6]) + rr5*(gl[2] + gl[7] + gl[8]) + rr7*(gl[3] + gl[5]) + rr9*gl[4];
        ei = 0.50*(rr3*(gli[1] + gli[6])*psc3 + rr5*(gli[2] + gli[7])*psc5 + rr7*gli[3] * psc7);
        e = f*mscale[kk] * e * fQ;
        ei = f *ei * fQ;
        }
        else{
        e = rr1*gl[0] + rr3*(gl[1] + gl[6]) + rr5*(gl[2] + gl[7] + gl[8]) + rr7*(gl[3] + gl[5]) + rr9*gl[4];
        ei = 0.50*(rr3*(gli[1] + gli[6])*psc3 + rr5*(gli[2] + gli[7])*psc5 + rr7*gli[3] * psc7);
        e = f*mscale[kk] * e;
        ei = f *ei;
        }


        em = em + e;
        ep = ep + ei;




        tempmulti = em;
        temppol = ep;
        //std::cout << "Test1\n";

              //!intermediate variables for the permanent components

        gf[1] = rr3*gl[0] + rr5*(gl[1] + gl[6]) + rr7*(gl[2] + gl[7] + gl[8]) + rr9*(gl[3] + gl[5]) + rr11*gl[4];
        gf[2] = -ck*rr3 + sc[4] * rr5 - sc[6] * rr7;
        gf[3] = ci*rr3 + sc[3] * rr5 + sc[5] * rr7;
        gf[4] = 2.0 *rr5;
        gf[5] = 2.0 * (-ck*rr5 + sc[4] * rr7 - sc[6] * rr9);
        gf[6] = 2.0* (-ci*rr5 - sc[3] * rr7 - sc[5] * rr9);
        gf[7] = 4.0 *rr7;

        //!intermediate variables for the induced components

        gfi[1] = 0.5 * rr5 * ((gli[1] + gli[6])*psc3 + (glip[1] + glip[6])*dsc3 + scip[2] * scale3i) + 0.5*rr7*((gli[7] + gli[2])*psc5 + (glip[7] + glip[2])*dsc5 - (sci[3] * scip[4] + scip[3] * sci[4])*scale5i) + 0.5*rr9*(gli[3] * psc7 + glip[3] * dsc7);
        gfi[2] = -rr3*ck + rr5*sc[4] - rr7*sc[6];
        gfi[3] = rr3*ci + rr5*sc[3] + rr7*sc[5];
        gfi[4] = 2.0 * rr5;
        gfi[5] = rr7 * (sci[4] * psc7 + scip[4] * dsc7);
        gfi[6] = -rr7 * (sci[3] * psc7 + scip[3] * dsc7);

        //! get the permanent force components

        ftm2[1] = gf[1] * xr + gf[2] * di[0] + gf[3] * dk[0] + gf[4] * (qkdi[0] - qidk[0]) + gf[5] * qir[0] + gf[6] * qkr[0] + gf[7] * (qiqkr[0] + qkqir[0]);
        ftm2[2] = gf[1] * yr + gf[2] * di[1] + gf[3] * dk[1] + gf[4] * (qkdi[1] - qidk[1]) + gf[5] * qir[1] + gf[6] * qkr[1] + gf[7] * (qiqkr[1] + qkqir[1]);
        ftm2[3] = gf[1] * zr + gf[2] * di[2] + gf[3] * dk[2] + gf[4] * (qkdi[2] - qidk[2]) + gf[5] * qir[2] + gf[6] * qkr[2] + gf[7] * (qiqkr[2] + qkqir[2]);


        //std::cout << "Test2\n";
        //! get the induced force components

        ftm2i[1] = gfi[1] * xr + 0.5*(-rr3*ck*(uind[1][i] * psc3 + uinp[1][i] * dsc3) + rr5*sc[4] * (uind[1][i] * psc5 + uinp[1][i] * dsc5) - rr7*sc[6] * (uind[1][i] * psc7 + uinp[1][i] * dsc7))
          + (rr3*ci*(uind[1][k] * psc3 + uinp[1][k] * dsc3) + rr5*sc[3] * (uind[1][k] * psc5 + uinp[1][k] * dsc5) + rr7*sc[5] * (uind[1][k] * psc7 + uinp[1][k] * dsc7))*0.5
          + rr5*scale5i*(sci[4] * uinp[1][i] + scip[4] * uind[1][i] + sci[3] * uinp[1][k] + scip[3] * uind[1][k])*0.5
          + 0.5*(sci[4] * psc5 + scip[4] * dsc5)*rr5*di[0] + 0.5*(sci[3] * psc5 + scip[3] * dsc5)*rr5*dk[0] + 0.5* gfi[4] * ((qkui[0] - qiuk[0])*psc5 + (qkuip[0] - qiukp[0])*dsc5)
          + gfi[5] * qir[0] + gfi[6] * qkr[0];

        ftm2i[2] = gfi[1] * yr + 0.5*(-rr3*ck*(uind[2][i] * psc3 + uinp[2][i] * dsc3) + rr5*sc[4] * (uind[2][i] * psc5 + uinp[2][i] * dsc5) - rr7*sc[6] * (uind[2][i] * psc7 + uinp[2][i] * dsc7))
          + (rr3*ci*(uind[2][k] * psc3 + uinp[2][k] * dsc3) + rr5*sc[3] * (uind[2][k] * psc5 + uinp[2][k] * dsc5) + rr7*sc[5] * (uind[2][k] * psc7 + uinp[2][k] * dsc7))*0.5
          + rr5*scale5i*(sci[4] * uinp[2][i] + scip[4] * uind[2][i] + sci[3] * uinp[2][k] + scip[3] * uind[2][k])*0.5
          + 0.5*(sci[4] * psc5 + scip[4] * dsc5)*rr5*di[1] + 0.5*(sci[3] * psc5 + scip[3] * dsc5)*rr5*dk[1] + 0.5* gfi[4] * ((qkui[1] - qiuk[1])*psc5 + (qkuip[1] - qiukp[1])*dsc5)
          + gfi[5] * qir[1] + gfi[6] * qkr[1];
        ftm2i[3] = gfi[1] * zr + 0.5*(-rr3*ck*(uind[3][i] * psc3 + uinp[3][i] * dsc3) + rr5*sc[4] * (uind[3][i] * psc5 + uinp[3][i] * dsc5) - rr7*sc[6] * (uind[3][i] * psc7 + uinp[3][i] * dsc7))
          + (rr3*ci*(uind[3][k] * psc3 + uinp[3][k] * dsc3) + rr5*sc[3] * (uind[3][k] * psc5 + uinp[3][k] * dsc5) + rr7*sc[5] * (uind[3][k] * psc7 + uinp[3][k] * dsc7))*0.5
          + rr5*scale5i*(sci[4] * uinp[3][i] + scip[4] * uind[3][i] + sci[3] * uinp[3][k] + scip[3] * uind[3][k])*0.5
          + 0.5*(sci[4] * psc5 + scip[4] * dsc5)*rr5*di[2] + 0.5*(sci[3] * psc5 + scip[3] * dsc5)*rr5*dk[2] + 0.5* gfi[4] * ((qkui[2] - qiuk[2])*psc5 + (qkuip[2] - qiukp[2])*dsc5)
          + gfi[5] * qir[2] + gfi[6] * qkr[2];


        //! account for part. excluded induced interactions

        temp3 = 0.5 * rr3 * ((gli[1] + gli[6])*pscale[kk] + (glip[1] + glip[6]) * dscale[kk]);
        temp5 = 0.5 * rr5 * ((gli[2] + gli[7])*pscale[kk] + (glip[2] + glip[7]) * dscale[kk]);
        temp7 = 0.5 * rr7 * (gli[3] * pscale[kk] + glip[3] * dscale[kk]);

        fridmp[1] = temp3 * ddsc3[0] + temp5 * ddsc5[0] + temp7 * ddsc7[0];
        fridmp[2] = temp3 * ddsc3[1] + temp5 * ddsc5[1] + temp7 * ddsc7[1];
        fridmp[3] = temp3 * ddsc3[2] + temp5 * ddsc5[2] + temp7 * ddsc7[2];

        //! find some scaling for induced-induced force

        temp3 = 0.5 * rr3 * uscale[kk] * scip[2];
        temp5 = -0.5 * rr5 * uscale[kk] * (sci[3] * scip[4] + scip[3] * sci[4]);

        findmp[1] = temp3 * ddsc3[0] + temp5 * ddsc5[0];
        findmp[2] = temp3 * ddsc3[1] + temp5 * ddsc5[1];
        findmp[3] = temp3 * ddsc3[2] + temp5 * ddsc5[2];


        //std::cout << "Test3\n";


        //! modifiy induced force for partially excluded interactions

        ftm2i[1] = ftm2i[1] - fridmp[1] - findmp[1];
        ftm2i[2] = ftm2i[2] - fridmp[2] - findmp[2];
        ftm2i[3] = ftm2i[3] - fridmp[3] - findmp[3];


        //!intermediate terms for induced torque on multipoles

        gti[2] = 0.5 * rr5 * (sci[4] * psc5 + scip[4] * dsc5);
        gti[3] = 0.5 * rr5 * (sci[3] * psc5 + scip[3] * dsc5);
        gti[4] = gfi[4];
        gti[5] = gfi[5];
        gti[6] = gfi[6];

        //! permanent torque components

        ttm2[1] = -rr3*dixdk[0] + gf[2] * dixr[0] - gf[5] * rxqir[0] + gf[4] * (dixqkr[0] + dkxqir[0] + rxqidk[0] - 2.0 * qixqk[0]) - gf[7] * (rxqikr[0] + qkrxqir[0]);
        ttm2[2] = -rr3*dixdk[1] + gf[2] * dixr[1] - gf[5] * rxqir[1] + gf[4] * (dixqkr[1] + dkxqir[1] + rxqidk[1] - 2.0 * qixqk[1]) - gf[7] * (rxqikr[1] + qkrxqir[1]);
        ttm2[3] = -rr3*dixdk[2] + gf[2] * dixr[2] - gf[5] * rxqir[2] + gf[4] * (dixqkr[2] + dkxqir[2] + rxqidk[2] - 2.0 * qixqk[2]) - gf[7] * (rxqikr[2] + qkrxqir[2]);


        ttm3[1] = rr3*dixdk[0] + gf[3] * dkxr[0] - gf[6] * rxqkr[0] - gf[4] * (dixqkr[0] + dkxqir[0] + rxqkdi[0] - 2.0*qixqk[0]) - gf[7] * (rxqkir[0] - qkrxqir[0]);
        ttm3[2] = rr3*dixdk[1] + gf[3] * dkxr[1] - gf[6] * rxqkr[1] - gf[4] * (dixqkr[1] + dkxqir[1] + rxqkdi[1] - 2.0*qixqk[1]) - gf[7] * (rxqkir[1] - qkrxqir[1]);
        ttm3[3] = rr3*dixdk[2] + gf[3] * dkxr[2] - gf[6] * rxqkr[2] - gf[4] * (dixqkr[2] + dkxqir[2] + rxqkdi[2] - 2.0*qixqk[2]) - gf[7] * (rxqkir[2] - qkrxqir[2]);

        //! induced torque components

        ttm2i[1] = -rr3* (dixuk[0] * psc3 + dixukp[0] * dsc3) * 0.5 + gti[2] * dixr[0] + gti[4] * ((ukxqir[0] + rxqiuk[0]) * psc5 + (ukxqirp[0] + rxqiukp[0]) * dsc5) * 0.5 - gti[5] * rxqir[0];
        ttm2i[2] = -rr3* (dixuk[1] * psc3 + dixukp[1] * dsc3) * 0.5 + gti[2] * dixr[1] + gti[4] * ((ukxqir[1] + rxqiuk[1]) * psc5 + (ukxqirp[1] + rxqiukp[1]) * dsc5) * 0.5 - gti[5] * rxqir[1];
        ttm2i[3] = -rr3* (dixuk[2] * psc3 + dixukp[2] * dsc3) * 0.5 + gti[2] * dixr[2] + gti[4] * ((ukxqir[2] + rxqiuk[2]) * psc5 + (ukxqirp[2] + rxqiukp[2]) * dsc5) * 0.5 - gti[5] * rxqir[2];

        ttm3i[1] = -rr3 * (dkxui[0] * psc3 + dkxuip[0] * dsc3) * 0.5 + gti[3] * dkxr[0] - gti[4] * ((uixqkr[0] + rxqkui[0]) * psc5 + (uixqkrp[0] + rxqkuip[0]) * dsc5) *  0.5 - gti[6] * rxqkr[0];
        ttm3i[2] = -rr3 * (dkxui[1] * psc3 + dkxuip[1] * dsc3) * 0.5 + gti[3] * dkxr[1] - gti[4] * ((uixqkr[1] + rxqkui[1]) * psc5 + (uixqkrp[1] + rxqkuip[1]) * dsc5) *  0.5 - gti[6] * rxqkr[1];
        ttm3i[3] = -rr3 * (dkxui[2] * psc3 + dkxuip[2] * dsc3) * 0.5 + gti[3] * dkxr[2] - gti[4] * ((uixqkr[2] + rxqkui[2]) * psc5 + (uixqkrp[2] + rxqkuip[2]) * dsc5) *  0.5 - gti[6] * rxqkr[2];

        //! handle the case were scaling is used

        for (j = 1; j <= 3; j++) {
          ftm2[j] = f * ftm2[j] * mscale[kk];
          ftm2i[j] = f * ftm2i[j];
          ttm2[j] = f * ttm2[j] * mscale[kk];
          ttm2i[j] = f * ttm2i[j];
          ttm3[j] = f * ttm3[j] * mscale[kk];
          ttm3i[j] = f * ttm3i[j];

        }

        //! increment gradient due to force and torque on first sites
      /*	vector <double> test;
        test.reserve(n_atom + 100);*/

        dem[1][ii] += ftm2[1];
        dem[2][ii] += ftm2[2];
        dem[3][ii] += ftm2[3];


        dep[1][ii] += ftm2i[1];
        dep[2][ii] += ftm2i[2];
        dep[3][ii] += ftm2i[3];



        //std::cout << "Test4\n";

        // 	cout << ftm2[2] << endl;

        //! calling torque: convert single site torque to force

        std::vector <double> frcx(4), frcy(4), frcz(4);
        ptrdiff_t ia, ic, ib, id;
        std::vector<double> u(4), v(4), w(4), r(4), s(4), t(4), t1(4), t2(4);
        std::vector<double> uv(4), uw(4), vw(4), ur(4), us(4), vs(4), ws(4);
        double du(0.0), dv(0.0), dw(0.0), usiz(0.0), vsiz(0.0), wsiz(0.0), rsiz(0.0), ssiz(0.0);
        double t1siz(0.0), t2siz(0.0), uvsiz(0.0), uwsiz(0.0), vwsiz(0.0), ursiz(0.0), ussiz(0.0), vssiz(0.0), wssiz(0.0);
        double uvcos(0.0), uwcos(0.0), vwcos(0.0), urcos(0.0), uscos(0.0), vscos(0.0), wscos(0.0);
        double ut1cos(0.0), ut2cos(0.0), uvsin(0.0), uwsin(0.0), vwsin(0.0), ursin(0.0), ussin(0.0), vssin(0.0), wssin(0.0), ut1sin(0.0), ut2sin(0.0);
        double dphidu(0.0), dphidw(0.0), dphidv(0.0), dphids(0.0), dphidr(0.0);



        for (j = 1; j <= 3; j++) {

          frcz[j] = 0.0;
          frcx[j] = 0.0;
          frcy[j] = 0.0;
        }
        //std::cout << "Test5\n";
        ia = zaxis[i] + 1;
        ib = ipole[i] + 1;
        if (ib == 0)continue;
        ic = xaxis[i] + 1;
        id = yaxis[i] + 1;


        if (axistype[i] == 0) continue;

        u[1] = positions[ia - 1].x() - positions[ib - 1].x();
        u[2] = positions[ia - 1].y() - positions[ib - 1].y();
        u[3] = positions[ia - 1].z() - positions[ib - 1].z();


        if (axistype[i] != 1) {
          v[1] = positions[ic - 1].x() - positions[ib - 1].x();
          v[2] = positions[ic - 1].y() - positions[ib - 1].y();
          v[3] = positions[ic - 1].z() - positions[ib - 1].z();

        }
        if (axistype[i] == 4 || axistype[i] == 5) {
          w[1] = positions[id - 1].x() - positions[ib - 1].x();
          w[2] = positions[id - 1].y() - positions[ib - 1].y();
          w[3] = positions[id - 1].z() - positions[ib - 1].z();
        }
        else {
          w[1] = u[2] * v[3] - u[3] * v[2];
          w[2] = u[3] * v[1] - u[1] * v[3];
          w[3] = u[1] * v[2] - u[2] * v[1];
          // 	  cout<< w[1] << endl;
        }

        usiz = sqrt(u[1] * u[1] + u[2] * u[2] + u[3] * u[3]);
        vsiz = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
        wsiz = sqrt(w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
        if (usiz == 0.0) usiz = 1.0;
        if (vsiz == 0.0) vsiz = 1.0;
        if (wsiz == 0.0) wsiz = 1.0;
        // 	cout<< usiz << endl;
        for (j = 1; j <= 3; j++) {
          u[j] = u[j] / usiz;
          v[j] = v[j] / vsiz;
          w[j] = w[j] / wsiz;

        }

        if (axistype[i] == 4) {
          r[1] = v[1] + w[1];
          r[2] = v[2] + w[2];
          r[3] = v[3] + w[3];

          s[1] = u[2] * r[3] - u[3] * r[2];
          s[2] = u[3] * r[1] - u[1] * r[3];
          s[3] = u[1] * r[2] - u[2] * r[1];

          rsiz = sqrt(r[1] * r[1] + r[2] * r[2] + r[3] * r[3]);
          ssiz = sqrt(s[1] * s[1] + s[2] * s[2] + s[3] * s[3]);
          if (rsiz == 0.0) rsiz = 1.0;
          if (ssiz == 0.0) ssiz = 1.0;


          for (j = 1; j <= 3; j++) {
            r[j] = r[j] / rsiz;
            s[j] = s[j] / ssiz;
          }
        }
        //! find the perpendicularand angle for each pair of axes    

        uv[1] = v[2] * u[3] - v[3] * u[2];
        uv[2] = v[3] * u[1] - v[1] * u[3];
        uv[3] = v[1] * u[2] - v[2] * u[1];

        uw[1] = w[2] * u[3] - w[3] * u[2];
        uw[2] = w[3] * u[1] - w[1] * u[3];
        uw[3] = w[1] * u[2] - w[2] * u[1];

        vw[1] = w[2] * v[3] - w[3] * v[2];
        vw[2] = w[3] * v[1] - w[1] * v[3];
        vw[3] = w[1] * v[2] - w[2] * v[1];

        uvsiz = sqrt(uv[1] * uv[1] + uv[2] * uv[2] + uv[3] * uv[3]);
        uwsiz = sqrt(uw[1] * uw[1] + uw[2] * uw[2] + uw[3] * uw[3]);
        vwsiz = sqrt(vw[1] * vw[1] + vw[2] * vw[2] + vw[3] * vw[3]);
        if (uvsiz == 0.0) uvsiz = 1.0;
        if (uwsiz == 0.0) uwsiz = 1.0;
        if (vwsiz == 0.0) vwsiz = 1.0;

        for (j = 1; j <= 3; j++) {
          uv[j] = uv[j] / uvsiz;
          uw[j] = uw[j] / uwsiz;
          vw[j] = vw[j] / vwsiz;
        }


        if (axistype[i] == 4) {
          ur[1] = r[2] * u[3] - r[3] * u[2];
          ur[2] = r[3] * u[1] - r[1] * u[3];
          ur[3] = r[1] * u[2] - r[2] * u[1];
          us[1] = s[2] * u[3] - s[3] * u[2];
          us[2] = s[3] * u[1] - s[1] * u[3];
          us[3] = s[1] * u[2] - s[2] * u[1];
          vs[1] = s[2] * v[3] - s[3] * v[2];
          vs[2] = s[3] * v[1] - s[1] * v[3];
          vs[3] = s[1] * v[2] - s[2] * v[1];
          ws[1] = s[2] * w[3] - s[3] * w[2];
          ws[2] = s[3] * w[1] - s[1] * w[3];
          ws[3] = s[1] * w[2] - s[2] * w[1];

          ursiz = sqrt(ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]);
          ussiz = sqrt(us[1] * us[1] + us[2] * us[2] + us[3] * us[3]);
          vssiz = sqrt(vs[1] * vs[1] + vs[2] * vs[2] + vs[3] * vs[3]);
          wssiz = sqrt(ws[1] * ws[1] + ws[2] * ws[2] + ws[3] * ws[3]);

          if (ursiz == 0.0) ursiz = 1.0;
          if (ussiz == 0.0) ussiz = 1.0;
          if (vssiz == 0.0) vssiz = 1.0;
          if (wssiz == 0.0) wssiz = 1.0;
          for (j = 1; j <= 3; j++) {
            ur[j] = ur[j] / ursiz;
            us[j] = us[j] / ussiz;
            vs[j] = vs[j] / vssiz;
            ws[j] = ws[j] / wssiz;
          }
        }

        uvcos = u[1] * v[1] + u[2] * v[2] + u[3] * v[3];
        uvsin = sqrt(1.0 - uvcos*uvcos);
        uwcos = u[1] * w[1] + u[2] * w[2] + u[3] * w[3];
        uwsin = sqrt(1.0 - uwcos*uwcos);
        vwcos = v[1] * w[1] + v[2] * w[2] + v[3] * w[3];
        vwsin = sqrt(1.0 - vwcos*vwcos);
        // 	cout << u[1] << endl;
        if (axistype[i] == 4) {
          urcos = u[1] * r[1] + u[2] * r[2] + u[3] * r[3];
          ursin = sqrt(1.0 - urcos*urcos);
          uscos = u[1] * s[1] + u[2] * s[2] + u[3] * s[3];
          ussin = sqrt(1.0 - uscos*uscos);
          vscos = v[1] * s[1] + v[2] * s[2] + v[3] * s[3];
          vssin = sqrt(1.0 - vscos*vscos);
          wscos = w[1] * s[1] + w[2] * s[2] + w[3] * s[3];
          wssin = sqrt(1.0 - wscos*wscos);
        }

        //!compute the projection of v and w onto the ru-plane

        if (axistype[i] == 4) {
          for (j = 1; j <= 3; j++) {
            t1[j] = v[j] - s[j] * vscos;
            t2[j] = w[j] - s[j] * wscos;
          }
          t1siz = sqrt(t1[1] * t1[1] + t1[2] * t1[2] + t1[3] * t1[3]);
          t2siz = sqrt(t2[1] * t2[1] + t2[2] * t2[2] + t2[3] * t2[3]);
          if (t1siz == 0.0) t1siz = 1.0;
          if (t2siz == 0.0) t2siz = 1.0;
          for (j = 1; j <= 3; j++) {
            t1[j] = t1[j] / t1siz;
            t2[j] = t2[j] / t2siz;
          }
          ut1cos = u[1] * t1[1] + u[2] * t1[2] + u[3] * t1[3];
          ut1sin = sqrt(1.0 - ut1cos*ut1cos);
          ut2cos = u[1] * t2[1] + u[2] * t2[2] + u[3] * t2[3];
          ut2sin = sqrt(1.0 - ut2cos*ut2cos);
        }

        //! negative of dot product of torque with unit vectors gives
        //! result of the infinitesimal rotation along the rot vectors

        dphidu = -ttm2[1] * u[1] - ttm2[2] * u[2] - ttm2[3] * u[3];
        dphidv = -ttm2[1] * v[1] - ttm2[2] * v[2] - ttm2[3] * v[3];
        dphidw = -ttm2[1] * w[1] - ttm2[2] * w[2] - ttm2[3] * w[3];

        if (axistype[i] == 4) {
          dphidr = ttm2[1] * r[1] - ttm2[2] * r[2] - ttm2[3] * r[3];
          dphids = ttm2[1] * s[1] - ttm2[2] * s[2] - ttm2[3] * s[3];
        }
        //! force distribution for the z-only local coordinate frame	
        if (axistype[i] == 1) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + uw[j] * dphidw / usiz;
            dem[j][ia] = dem[j][ia] + du;
            dem[j][ib] = dem[j][ib] - du;
          }
          //! force distribution for the z-thenx local coordinate method    
        }
        else if (axistype[i] == 2) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + uw[j] * dphidw / usiz;
            dv = -uv[j] * dphidu / (vsiz*uvsin);
            dem[j][ia] = dem[j][ia] + du;
            dem[j][ic] = dem[j][ic] + dv;
            dem[j][ib] = dem[j][ib] - du - dv;
            frcz[j] = frcz[j] + du;
            frcx[j] = frcx[j] + dv;
            // 	  cout << dem[j][ia] << endl;


          }
        }
        else if (axistype[i] == 3) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + 0.5*uw[j] * dphidw / usiz;
            dv = -uv[j] * dphidu / (vsiz*uvsin) + 0.5*vw[j] * dphidw / vsiz;
            dem[j][ia] = dem[j][ia] + du;
            dem[j][ic] = dem[j][ic] + dv;
            dem[j][ib] = dem[j][ib] - du - dv;
            frcz[j] = frcz[j] + du;
            frcx[j] = frcx[j] + dv;


          }
        }
        else if (axistype[i] == 4) {
          for (j = 1; j <= 3; j++) {
            du = ur[j] * dphidr / (usiz*ursin) + us[j] * dphids / usiz;
            dv = (vssin*s[j] - vscos*t1[j])*dphidu / (vsiz*(ut1sin + ut2sin));
            dw = (wssin*s[j] - wscos*t2[j]) *dphidu / (wsiz*(ut1sin + ut2sin));
            dem[j][ia] = dem[j][ia] + du;
            dem[j][ic] = dem[j][ic] + dv;
            dem[j][id] = dem[j][id] + dw;
            dem[j][ib] = dem[j][ib] - du - dv - dw;
            frcz[j] = frcz[j] + du;
            frcx[j] = frcx[j] + dv;
            frcy[j] = frcy[j] + dw;
          }
        }
        else if (axistype[i] == 5) {
          for (j = 1; j <= 3; j++) {
            du = uw[j] * dphidw / (usiz*uwsin) + uv[j] * dphidv / (usiz*uvsin) - uw[j] * dphidu / (usiz * uwsin) - uv[j] * dphidu / (usiz*uvsin);
            dv = vw[j] * dphidw / (vsiz*vwsin) - uv[j] * dphidu / (vsiz*uvsin) - vw[j] * dphidv / (vsiz * vwsin) + uv[j] * dphidv / (vsiz*uvsin);
            dw = -uw[j] * dphidu / (wsiz*uwsin) - vw[j] * dphidv / (wsiz * vwsin) + uw[j] * dphidw / (wsiz*uwsin) + vw[j] * dphidw / (wsiz*vwsin);
            du = du / 3.0;
            dv = dv / 3.0;
            dw = dw / 3.0;

            dem[j][ia] = dem[j][ia] + du;
            dem[j][ic] = dem[j][ic] + dv;
            dem[j][id] = dem[j][id] + dw;
            dem[j][ib] = dem[j][ib] - du - dv - dw;
            frcz[j] = frcz[j] + du;
            frcx[j] = frcx[j] + dv;
            frcy[j] = frcy[j] + dw;
          }
        }
        //! negative of dot product of torque with unit vectors gives
        //! result of infinitesimal rotation along these vectors

        dphidu = -ttm2i[1] * u[1] - ttm2i[2] * u[2] - ttm2i[3] * u[3];
        dphidv = -ttm2i[1] * v[1] - ttm2i[2] * v[2] - ttm2i[3] * v[3];
        dphidw = -ttm2i[1] * w[1] - ttm2i[2] * w[2] - ttm2i[3] * w[3];

        if (axistype[i] == 4) {
          dphidr = ttm2i[1] * r[1] - ttm2i[2] * r[2] - ttm2i[3] * r[3];
          dphids = ttm2i[1] * s[1] - ttm2i[2] * s[2] - ttm2i[3] * s[3];

        }

        //! force distribution for the z-only local coordinate frame	
        if (axistype[i] == 1) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + uw[j] * dphidw / usiz;
            dep[j][ia] = dep[j][ia] + du;
            dep[j][ib] = dep[j][ib] - du;
            frcz[j] = frcz[j] + du;
          }
          //! force distribution for the z-thenx local coordinate method    
        }
        else if (axistype[i] == 2) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + uw[j] * dphidw / usiz;
            dv = -uv[j] * dphidu / (vsiz*uvsin);

            dep[j][ia] = dep[j][ia] + du;
            dep[j][ic] = dep[j][ic] + dv;
            dep[j][ib] = dep[j][ib] - du - dv;
            // 	cout<<dep[j][ic] << endl;
            frcz[j] = frcz[j] + du;
            frcx[j] = frcx[j] + dv;


          }
        }
        else if (axistype[i] == 3) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + 0.5*uw[j] * dphidw / usiz;
            dv = -uv[j] * dphidu / (vsiz*uvsin) + 0.5*vw[j] * dphidw / vsiz;

            dep[j][ia] = dep[j][ia] + du;
            dep[j][ic] = dep[j][ic] + dv;
            dep[j][ib] = dep[j][ib] - du - dv;

            frcz[j] = frcz[j] + du;
            frcx[j] = frcx[j] + dv;


          }
        }
        else if (axistype[i] == 4) {
          for (j = 1; j <= 3; j++) {
            du = ur[j] * dphidr / (usiz*ursin) + us[j] * dphids / usiz;
            dv = (vssin*s[j] - vscos*t1[j])*dphidu / (vsiz*(ut1sin + ut2sin));
            dw = (wssin*s[j] - wscos*t2[j]) *dphidu / (wsiz*(ut1sin + ut2sin));
            dep[j][ia] = dep[j][ia] + du;
            dep[j][ic] = dep[j][ic] + dv;
            dep[j][id] = dep[j][id] + dw;
            dep[j][ib] = dep[j][ib] - du - dv - dw;
            frcz[j] = frcz[j] + du;
            frcx[j] = frcx[j] + dv;
            frcy[j] = frcy[j] + dw;
          }
        }
        else if (axistype[i] == 5) {
          for (j = 1; j <= 3; j++) {
            du = uw[j] * dphidw / (usiz*uwsin) + uv[j] * dphidv / (usiz*uvsin) - uw[j] * dphidu / (usiz * uwsin) - uv[j] * dphidu / (usiz*uvsin);
            dv = vw[j] * dphidw / (vsiz*vwsin) - uv[j] * dphidu / (vsiz*uvsin) - vw[j] * dphidv / (vsiz * vwsin) + uv[j] * dphidv / (vsiz*uvsin);
            dw = -uw[j] * dphidu / (wsiz*uwsin) - vw[j] * dphidv / (wsiz * vwsin) + uw[j] * dphidw / (wsiz*uwsin) + vw[j] * dphidw / (wsiz*vwsin);
            du = du / 3.0;
            dv = dv / 3.0;
            dw = dw / 3.0;

            dep[j][ia] = dep[j][ia] + du;
            dep[j][ic] = dep[j][ic] + dv;
            dep[j][id] = dep[j][id] + dw;
            dep[j][ib] = dep[j][ib] - du - dv - dw;
            frcz[j] = frcz[j] + du;
            frcx[j] = frcx[j] + dv;
            frcy[j] = frcy[j] + dw;
          }
        }


        //! increment gradients due to force and torque on the second sites

        dem[1][kk] -= ftm2[1];
        dem[2][kk] -= ftm2[2];
        dem[3][kk] -= ftm2[3];

        dep[1][kk] -= ftm2i[1];
        dep[2][kk] -= ftm2i[2];
        dep[3][kk] -= ftm2i[3];



        //!"""""""""""""""""""222222222222222222222222222222"""""""""""""""""""""""""""""""""
        //!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        //!'''''''''''''''''''#############################""""""""""""""""""""""""""""""""""
        //! calling torque: convert single site torque to force
        // 	ptrdiff_t ia,ic,ib,id;
        // 	vector<double> u(4),v(4),w(4), r(4),s(4),t(4),t1(4),t2(4);
        // 	vector<double> uv(4), uw(4), vw(4), ur(4), us(4), vs(4), ws(4);
        // 	double du, dv, dw, random, usiz, vsiz, wsiz, rsiz, ssiz;
        // 	double t1siz, t2siz, uvsiz, uwsiz, vwsiz, ursiz, ussiz, vssiz, wssiz;
        // 	double uvcos, uwcos, vwcos, urcos, uscos, vscos, wscos;
        // 	double ut1cos, ut2cos, uvsin, uwsin, vwsin, ursin, ussin, vssin, wssin,ut1sin, ut2sin;
        // 	double dphidu,dphidw,dphidv, dphids, dphidr;
        // 	ptrdiff_t tempi;
        // 	tempi=i;
        // 	i=k;

        for (j = 1; j <= 3; j++) {

          frczk[j] = 0.0;
          frcxk[j] = 0.0;
          frcyk[j] = 0.0;
        }

        ia = zaxis[k] + 1;
        ib = ipole[k] + 1;
        if (ib == 0)continue;
        ic = xaxis[k] + 1;
        id = yaxis[k] + 1;
        // 	string axetype;


        if (axistype[k] == 0) continue;

        u[1] = positions[ia - 1].x() - positions[ib - 1].x();
        u[2] = positions[ia - 1].y() - positions[ib - 1].y();
        u[3] = positions[ia - 1].z() - positions[ib - 1].z();


        if (axistype[k] != 1) {
          v[1] = positions[ic - 1].x() - positions[ib - 1].x();
          v[2] = positions[ic - 1].y() - positions[ib - 1].y();
          v[3] = positions[ic - 1].z() - positions[ib - 1].z();

        }
        if (axistype[k] == 4 || axistype[k] == 5) {
          w[1] = positions[id - 1].x() - positions[ib - 1].x();
          w[2] = positions[id - 1].y() - positions[ib - 1].y();
          w[3] = positions[id - 1].z() - positions[ib - 1].z();
        }
        else {
          w[1] = u[2] * v[3] - u[3] * v[2];
          w[2] = u[3] * v[1] - u[1] * v[3];
          w[3] = u[1] * v[2] - u[2] * v[1];

        }
        usiz = sqrt(u[1] * u[1] + u[2] * u[2] + u[3] * u[3]);
        vsiz = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
        wsiz = sqrt(w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
        if (usiz == 0.0) usiz = 1.0;
        if (vsiz == 0.0) vsiz = 1.0;
        if (wsiz == 0.0) wsiz = 1.0;

        for (j = 1; j <= 3; j++) {
          u[j] = u[j] / usiz;
          v[j] = v[j] / vsiz;
          w[j] = w[j] / wsiz;

        }

        if (axistype[k] == 4) {
          r[1] = v[1] + w[1];
          r[2] = v[2] + w[2];
          r[3] = v[3] + w[3];

          s[1] = u[2] * r[3] - u[3] * r[2];
          s[2] = u[3] * r[1] - u[1] * r[3];
          s[3] = u[1] * r[2] - u[2] * r[1];

          rsiz = sqrt(r[1] * r[1] + r[2] * r[2] + r[3] * r[3]);
          ssiz = sqrt(s[1] * s[1] + s[2] * s[2] + s[3] * s[3]);

          if (rsiz == 0.0) rsiz = 1.0;
          if (ssiz == 0.0) ssiz = 1.0;

          for (j = 1; j <= 3; j++) {
            r[j] = r[j] / rsiz;
            s[j] = s[j] / ssiz;
          }
        }
        //! find the perpendicularand angle for each pair of axes       
        uv[1] = v[2] * u[3] - v[3] * u[2];
        uv[2] = v[3] * u[1] - v[1] * u[3];
        uv[3] = v[1] * u[2] - v[2] * u[1];

        uw[1] = w[2] * u[3] - w[3] * u[2];
        uw[2] = w[3] * u[1] - w[1] * u[3];
        uw[3] = w[1] * u[2] - w[2] * u[1];

        vw[1] = w[2] * v[3] - w[3] * v[2];
        vw[2] = w[3] * v[1] - w[1] * v[3];
        vw[3] = w[1] * v[2] - w[2] * v[1];

        uvsiz = sqrt(uv[1] * uv[1] + uv[2] * uv[2] + uv[3] * uv[3]);
        uwsiz = sqrt(uw[1] * uw[1] + uw[2] * uw[2] + uw[3] * uw[3]);
        vwsiz = sqrt(vw[1] * vw[1] + vw[2] * vw[2] + vw[3] * vw[3]);

        if (uvsiz == 0.0) uvsiz = 1.0;
        if (uwsiz == 0.0) uwsiz = 1.0;
        if (vwsiz == 0.0) vwsiz = 1.0;

        for (j = 1; j <= 3; j++) {
          uv[j] = uv[j] / uvsiz;
          uw[j] = uw[j] / uwsiz;
          vw[j] = vw[j] / vwsiz;
        }


        if (axistype[k] == 4) {
          ur[1] = r[2] * u[3] - r[3] * u[2];
          ur[2] = r[3] * u[1] - r[1] * u[3];
          ur[3] = r[1] * u[2] - r[2] * u[1];
          us[1] = s[2] * u[3] - s[3] * u[2];
          us[2] = s[3] * u[1] - s[1] * u[3];
          us[3] = s[1] * u[2] - s[2] * u[1];
          vs[1] = s[2] * v[3] - s[3] * v[2];
          vs[2] = s[3] * v[1] - s[1] * v[3];
          vs[3] = s[1] * v[2] - s[2] * v[1];
          ws[1] = s[2] * w[3] - s[3] * w[2];
          ws[2] = s[3] * w[1] - s[1] * w[3];
          ws[3] = s[1] * w[2] - s[2] * w[1];

          ursiz = sqrt(ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]);
          ussiz = sqrt(us[1] * us[1] + us[2] * us[2] + us[3] * us[3]);
          vssiz = sqrt(vs[1] * vs[1] + vs[2] * vs[2] + vs[3] * vs[3]);
          wssiz = sqrt(ws[1] * ws[1] + ws[2] * ws[2] + ws[3] * ws[3]);

          if (ursiz == 0.0) ursiz = 1.0;
          if (ussiz == 0.0) ussiz = 1.0;
          if (vssiz == 0.0) vssiz = 1.0;
          if (wssiz == 0.0) wssiz = 1.0;

          for (j = 1; j <= 3; j++) {
            ur[j] = ur[j] / ursiz;
            us[j] = us[j] / ussiz;
            vs[j] = vs[j] / vssiz;
            ws[j] = ws[j] / wssiz;
          }
        }

        uvcos = u[1] * v[1] + u[2] * v[2] + u[3] * v[3];
        uvsin = sqrt(1.0 - uvcos*uvcos);
        uwcos = u[1] * w[1] + u[2] * w[2] + u[3] * w[3];
        uwsin = sqrt(1.0 - uwcos*uwcos);
        vwcos = v[1] * w[1] + v[2] * w[2] + v[3] * w[3];
        vwsin = sqrt(1.0 - vwcos*vwcos);

        if (axistype[k] == 4) {
          urcos = u[1] * r[1] + u[2] * r[2] + u[3] * r[3];
          ursin = sqrt(1.0 - urcos*urcos);
          uscos = u[1] * s[1] + u[2] * s[2] + u[3] * s[3];
          ussin = sqrt(1.0 - uscos*uscos);
          vscos = v[1] * s[1] + v[2] * s[2] + v[3] * s[3];
          vssin = sqrt(1.0 - vscos*vscos);
          wscos = w[1] * s[1] + w[2] * s[2] + w[3] * w[3];
          wssin = sqrt(1.0 - wscos*wscos);
        }

        //!compute the projection of v and w onto the ru-plane

        if (axistype[k] == 4) {
          for (j = 1; j <= 3; j++) {
            t1[j] = v[j] - s[j] * vscos;
            t2[j] = w[j] - s[j] * wscos;
          }
          t1siz = sqrt(t1[1] * t1[1] + t1[2] * t1[2] + t1[3] * t1[3]);
          t2siz = sqrt(t2[1] * t2[1] + t2[2] * t2[2] + t2[3] * t2[3]);

          if (t1siz == 0.0) t1siz = 1.0;
          if (t2siz == 0.0) t2siz = 1.0;

          for (j = 1; j <= 3; j++) {
            t1[j] = t1[j] / t1siz;
            t2[j] = t2[j] / t2siz;
          }
          ut1cos = u[1] * t1[1] + u[2] * t1[2] + u[3] * t1[3];
          ut1sin = sqrt(1.0 - ut1cos*ut1cos);
          ut2cos = u[1] * t2[1] + u[2] * t2[2] + u[3] * t2[3];
          ut2sin = sqrt(1.0 - ut2cos*ut2cos);
        }

        //! negative of dot product of torque with unit vectors gives
        //! result of the infinitesimal rotation along the rot vectors

        dphidu = -ttm3[1] * u[1] - ttm3[2] * u[2] - ttm3[3] * u[3];
        dphidv = -ttm3[1] * v[1] - ttm3[2] * v[2] - ttm3[3] * v[3];
        dphidw = -ttm3[1] * w[1] - ttm3[2] * w[2] - ttm3[3] * w[3];

        if (axistype[k] == 4) {
          dphidr = ttm3[1] * r[1] - ttm3[2] * r[2] - ttm3[3] * r[3];
          dphids = ttm3[1] * s[1] - ttm3[2] * s[2] - ttm3[3] * s[3];
        }
        //! force distribution for the z-only local coordinate frame	
        if (axistype[k] == 1) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + uw[j] * dphidw / usiz;
            dem[j][ia] = dem[j][ia] + du;
            dem[j][ib] = dem[j][ib] - du;
          }
          //! force distribution for the z-thenx local coordinate method    
        }
        else if (axistype[k] == 2) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + uw[j] * dphidw / usiz;
            dv = -uv[j] * dphidu / (vsiz*uvsin);
            dem[j][ia] = dem[j][ia] + du;
            dem[j][ic] = dem[j][ic] + dv;
            dem[j][ib] = dem[j][ib] - du - dv;
            frczk[j] = frczk[j] + du;
            frcxk[j] = frcxk[j] + dv;


          }
        }
        else if (axistype[k] == 3) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + 0.5*uw[j] * dphidw / usiz;
            dv = -uv[j] * dphidu / (vsiz*uvsin) + 0.5*vw[j] * dphidw / vsiz;
            dem[j][ia] = dem[j][ia] + du;
            dem[j][ic] = dem[j][ic] + dv;
            dem[j][ib] = dem[j][ib] - du - dv;
            frczk[j] = frczk[j] + du;
            frcxk[j] = frcxk[j] + dv;


          }
        }
        else if (axistype[k] == 4) {
          for (j = 1; j <= 3; j++) {
            du = ur[j] * dphidr / (usiz*ursin) + us[j] * dphids / usiz;
            dv = (vssin*s[j] - vscos*t1[j])*dphidu / (vsiz*(ut1sin + ut2sin));
            dw = (wssin*s[j] - wscos*t2[j]) *dphidu / (wsiz*(ut1sin + ut2sin));
            dem[j][ia] = dem[j][ia] + du;
            dem[j][ic] = dem[j][ic] + dv;
            dem[j][id] = dem[j][id] + dw;
            dem[j][ib] = dem[j][ib] - du - dv - dw;
            frczk[j] = frczk[j] + du;
            frcxk[j] = frcxk[j] + dv;
            frcyk[j] = frcyk[j] + dw;
          }
        }
        else if (axistype[k] == 5) {
          for (j = 1; j <= 3; j++) {
            du = uw[j] * dphidw / (usiz*uwsin) + uv[j] * dphidv / (usiz*uvsin) - uw[j] * dphidu / (usiz * uwsin) - uv[j] * dphidu / (usiz*uvsin);
            dv = vw[j] * dphidw / (vsiz*vwsin) - uv[j] * dphidu / (vsiz*uvsin) - vw[j] * dphidv / (vsiz * vwsin) + uv[j] * dphidv / (vsiz*uvsin);
            dw = -uw[j] * dphidu / (wsiz*uwsin) - vw[j] * dphidv / (wsiz * vwsin) + uw[j] * dphidw / (wsiz*uwsin) + vw[j] * dphidw / (wsiz*vwsin);
            du = du / 3.0;
            dv = dv / 3.0;
            dw = dw / 3.0;

            dem[j][ia] = dem[j][ia] + du;
            dem[j][ic] = dem[j][ic] + dv;
            dem[j][id] = dem[j][id] + dw;
            dem[j][ib] = dem[j][ib] - du - dv - dw;
            frczk[j] = frczk[j] + du;
            frcxk[j] = frcxk[j] + dv;
            frcyk[j] = frcyk[j] + dw;
          }
        }
        //! negative of dot product of torque with unit vectors gives
        //! result of infinitesimal rotation along these vectors

        dphidu = -ttm3i[1] * u[1] - ttm3i[2] * u[2] - ttm3i[3] * u[3];
        dphidv = -ttm3i[1] * v[1] - ttm3i[2] * v[2] - ttm3i[3] * v[3];
        dphidw = -ttm3i[1] * w[1] - ttm3i[2] * w[2] - ttm3i[3] * w[3];

        if (axistype[k] == 4) {
          dphidr = ttm3i[1] * r[1] - ttm3i[2] * r[2] - ttm3i[3] * r[3];
          dphids = ttm3i[1] * s[1] - ttm3i[2] * s[2] - ttm3i[3] * s[3];

        }

        //! force distribution for the z-only local coordinate frame	
        if (axistype[k] == 1) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + uw[j] * dphidw / usiz;
            dep[j][ia] = dep[j][ia] + du;
            dep[j][ib] = dep[j][ib] - du;
          }
          //! force distribution for the z-thenx local coordinate method    
        }
        else if (axistype[k] == 2) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + uw[j] * dphidw / usiz;
            dv = -uv[j] * dphidu / (vsiz*uvsin);
            dep[j][ia] = dep[j][ia] + du;
            dep[j][ic] = dep[j][ic] + dv;
            dep[j][ib] = dep[j][ib] - du - dv;
            frczk[j] = frczk[j] + du;
            frcxk[j] = frcxk[j] + dv;


          }
        }
        else if (axistype[k] == 3) {
          for (j = 1; j <= 3; j++) {
            du = uv[j] * dphidv / (usiz*uvsin) + 0.5*uw[j] * dphidw / usiz;
            dv = -uv[j] * dphidu / (vsiz*uvsin) + 0.5*vw[j] * dphidw / vsiz;
            dep[j][ia] = dep[j][ia] + du;
            dep[j][ic] = dep[j][ic] + dv;
            dep[j][ib] = dep[j][ib] - du - dv;
            frczk[j] = frczk[j] + du;
            frcxk[j] = frcxk[j] + dv;


          }
        }
        else if (axistype[k] == 4) {
          for (j = 1; j <= 3; j++) {
            du = ur[j] * dphidr / (usiz*ursin) + us[j] * dphids / usiz;
            dv = (vssin*s[j] - vscos*t1[j])*dphidu / (vsiz*(ut1sin + ut2sin));
            dw = (wssin*s[j] - wscos*t2[j]) *dphidu / (wsiz*(ut1sin + ut2sin));
            dep[j][ia] = dep[j][ia] + du;
            dep[j][ic] = dep[j][ic] + dv;
            dep[j][id] = dep[j][id] + dw;
            dep[j][ib] = dep[j][ib] - du - dv - dw;
            frczk[j] = frczk[j] + du;
            frcxk[j] = frcxk[j] + dv;
            frcyk[j] = frcyk[j] + dw;
          }
        }
        else if (axistype[k] == 5) {
          for (j = 1; j <= 3; j++) {
            du = uw[j] * dphidw / (usiz*uwsin) + uv[j] * dphidv / (usiz*uvsin) - uw[j] * dphidu / (usiz * uwsin) - uv[j] * dphidu / (usiz*uvsin);
            dv = vw[j] * dphidw / (vsiz*vwsin) - uv[j] * dphidu / (vsiz*uvsin) - vw[j] * dphidv / (vsiz * vwsin) + uv[j] * dphidv / (vsiz*uvsin);
            dw = -uw[j] * dphidu / (wsiz*uwsin) - vw[j] * dphidv / (wsiz * vwsin) + uw[j] * dphidw / (wsiz*uwsin) + vw[j] * dphidw / (wsiz*vwsin);
            du = du / 3.0;
            dv = dv / 3.0;
            dw = dw / 3.0;

            dep[j][ia] = dep[j][ia] + du;
            dep[j][ic] = dep[j][ic] + dv;
            dep[j][id] = dep[j][id] + dw;
            dep[j][ib] = dep[j][ib] - du - dv - dw;
            frczk[j] = frczk[j] + du;
            frcxk[j] = frcxk[j] + dv;
            frcyk[j] = frcyk[j] + dw;
          }
        }


      }
    }


    for (j = 0; j < coords->atoms(ii - 1).bonds().size(); j++) {
      pscale[coords->atoms(ii - 1).bonds()[j] + 1] = 1.0;
      mscale[coords->atoms(ii - 1).bonds()[j] + 1] = 1.0;

    }
    for (j = 1; j <= n13[ii]; j++) {
      pscale[i13[j][ii]] = 1.0;
      mscale[i13[j][ii]] = 1.0;

    }
    for (j = 1; j <= n14[ii]; j++) {
      pscale[i14[j][ii]] = 1.0;
      mscale[i14[j][ii]] = 1.0;


    }

    for (j = 1; j <= n15[ii]; j++) {
      pscale[i15[j][ii]] = 1.0;
      mscale[i15[j][ii]] = 1.0;

    }
    for (j = 1; j <= np11[ii]; j++) {
      dscale[ip11[j][ii]] = 1.0;
      uscale[ip11[j][ii]] = 1.0;
    }
    for (j = 1; j <= np12[ii]; j++) {
      dscale[ip12[j][ii]] = 1.0;
      uscale[ip12[j][ii]] = 1.0;
    }
    for (j = 1; j <= np13[ii]; j++) {
      dscale[ip13[j][ii]] = 1.0;
      uscale[ip13[j][ii]] = 1.0;
    }
    for (j = 1; j <= np14[ii]; j++) {

      dscale[ip14[j][ii]] = 1.0;
      uscale[ip14[j][ii]] = 1.0;
    }





  }

  for (size_t n = 1; n <= N; n++)
  {
    gv_multipole.x() = dem[1][n];
    gv_multipole.y() = dem[2][n];
    gv_multipole.z() = dem[3][n];
    gv_polarization.x() = dep[1][n];
    gv_polarization.y() = dep[2][n];
    gv_polarization.z() = dep[3][n];

    part_grad[MULTIPOLE][n - 1] = gv_multipole;
    part_grad[POLARIZE][n - 1] = gv_polarization;


  }





}

size_t energy::interfaces::amoeba::amoeba_ff::multipole_sites(void)
{

  for (auto axes : refined.multipole_vecs())
  {
    for (auto mult : axes)
    {
      alloc_glob = mult.npole;

    }
  }
  quadro.clear();
  dipole.clear();
  charges.clear();
  ipole.clear();
  pdamp.clear();
  axistype.clear();
  thole.clear();
  polarity.clear();
  xaxis.clear();
  zaxis.clear();
  yaxis.clear();
  quadro.resize(alloc_glob + 1);
  dipole.resize(alloc_glob + 1);
  charges.resize(alloc_glob + 1);
  ipole.resize(alloc_glob + 1);
  pdamp.resize(alloc_glob + 1);
  thole.resize(alloc_glob + 1);
  polarity.resize(alloc_glob + 1);
  axistype.resize(alloc_glob + 1);
  xaxis.resize(alloc_glob + 1);
  yaxis.resize(alloc_glob + 1);
  zaxis.resize(alloc_glob + 1);
  plrgrp.resize(8u);
  size_t m(0U);
  for (auto axes : refined.multipole_vecs())
  {
    for (auto mult : axes)
    {
      m++;
      charges[m] = (mult.p_rot.charge);
      dipole[m].push_back(mult.p_rot.dipole.x());
      dipole[m].push_back(mult.p_rot.dipole.y());
      dipole[m].push_back(mult.p_rot.dipole.z());
      quadro[m].push_back(mult.p_rot.quadrupole[0]);
      quadro[m].push_back(mult.p_rot.quadrupole[1]);
      quadro[m].push_back(mult.p_rot.quadrupole[2]);
      quadro[m].push_back(mult.p_rot.quadrupole[3]);
      quadro[m].push_back(mult.p_rot.quadrupole[4]);
      quadro[m].push_back(mult.p_rot.quadrupole[5]);
      quadro[m].push_back(mult.p_rot.quadrupole[6]);
      quadro[m].push_back(mult.p_rot.quadrupole[7]);
      quadro[m].push_back(mult.p_rot.quadrupole[8]);
      ipole[m] = mult.center;
      axistype[m] = mult.p_rot.axt;
      xaxis[m] = mult.axes.x();
      yaxis[m] = mult.axes.y();
      zaxis[m] = mult.axes.z();




    }

  }

  m = 0;

  for (auto polaxes : refined.polarize_vecs())
  {
    for (auto pol : polaxes)
    {
      m++;
      pdamp[m] = pol.pdamp;
      thole[m] = pol.ff;
      polarity[m] = pol.force;
      plrgrp[0].push_back(pol.atoms[0]);
      plrgrp[1].push_back(pol.atoms[1]);
      plrgrp[2].push_back(pol.atoms[2]);
    }

  }

  //for (auto polaxes : cparams.polarizes())
  //{
  //	std::cout << polaxes.bond_index[0] << '\n';
  //}


  return alloc_glob;






}
void energy::interfaces::amoeba::amoeba_ff::refine_vdw_h_bonds(std::vector< ::tinker::refine::types::nbpair> const & pairlist,
  scon::matrix< ::tinker::parameter::combi::vdwc, true> const & params)
{



  std::vector<double> reduce(coords->xyz().size());


  for (auto const & pair : pairlist)
  {

    ::tinker::parameter::combi::vdwc const & p(params(refined.type(pair.a), refined.type(pair.b)));

    if (p.RR != 0.0)
    {

      reduce[pair.a] = p.RR;

    }
    else
    {
      reduce[pair.a] = 1.0;
      reduce[pair.b] = 1.0;

    }
  }
  coords::Cartesian_Point intxyz;
  coords::Representation_3D int2xyz(coords->xyz().size());
  vdwnew = coords->xyz();
  for (auto bnd : refined.bonds())
  {
    intxyz.x() = reduce[bnd.atoms[0]] * (coords->xyz(bnd.atoms[0]).x() - coords->xyz(bnd.atoms[1]).x()) + coords->xyz(bnd.atoms[1]).x();
    intxyz.y() = reduce[bnd.atoms[0]] * (coords->xyz(bnd.atoms[0]).y() - coords->xyz(bnd.atoms[1]).y()) + coords->xyz(bnd.atoms[1]).y();
    intxyz.z() = reduce[bnd.atoms[0]] * (coords->xyz(bnd.atoms[0]).z() - coords->xyz(bnd.atoms[1]).z()) + coords->xyz(bnd.atoms[1]).z();
    int2xyz[bnd.atoms[0]] = intxyz;
    intxyz.x() = reduce[bnd.atoms[1]] * (coords->xyz(bnd.atoms[1]).x() - coords->xyz(bnd.atoms[0]).x()) + coords->xyz(bnd.atoms[0]).x();
    intxyz.y() = reduce[bnd.atoms[1]] * (coords->xyz(bnd.atoms[1]).y() - coords->xyz(bnd.atoms[0]).y()) + coords->xyz(bnd.atoms[0]).y();
    intxyz.z() = reduce[bnd.atoms[1]] * (coords->xyz(bnd.atoms[1]).z() - coords->xyz(bnd.atoms[0]).z()) + coords->xyz(bnd.atoms[0]).z();
    int2xyz[bnd.atoms[1]] = intxyz;

  }
  vdwnew = (int2xyz);

}

void energy::interfaces::amoeba::amoeba_ff::refine_pair_lists(void)
{

  size_t n_atom, i(0U), jj(0U), kk(0U);
  n_atom = (coords->size());
  bool control_seq(false);

  //binary.atoms[0];
  n13.resize(alloc_glob + 50);

  i13.resize(alloc_glob + 50);
  i12.resize(refined.bonds().size());

  for (size_t n = 0; n < i13.size(); n++) {
    i13[n].resize(alloc_glob + 50);
  }


  //POLARIZATION GROUP VARIABLES

  size_t it(0U);
  size_t jt(0U);
  bool done(false);
  size_t start(0U), stop(0U), nlist(0U), nkeep(0U), maxval(0U);

  np11.resize(alloc_glob + 50);
  ip11.resize(alloc_glob + 50);
  list.resize(alloc_glob + 50);
  mask.resize(alloc_glob + 50);
  keep.resize(alloc_glob + 50);



  //!1-3 attached atom list for multipoles
  for (i = 1; i <= n_atom; i++) {
    n13[i] = 0;
    /* cout << C->topology[i-1].V[1]+1<<endl;*/
    for (size_t j = 0; j < coords->atoms(i - 1).bonds().size(); j++) {
      jj = coords->atoms(i - 1).bonds()[j] + 1;
      //	std::cout << i12[i].size() << '\n';

      for (size_t k = 0; k < coords->atoms(jj - 1).bonds().size(); k++) {

        kk = coords->atoms(jj - 1).bonds()[k] + 1;
        //std::cout << kk << '\n';
        if (kk == i)continue;

        for (size_t m = 0; m < coords->atoms(i - 1).bonds().size(); m++) {

          if (kk == coords->atoms(i - 1).bonds()[m] + 1) break;
        }

        n13[i] = n13[i] + 1;

        i13[n13[i]][i] = kk;
        //std::cout << " n13 " << n13[i] << " i13  " << i13[n13[i]][i]<<'\n';

      }
    }
  }

  //!1-4 attached atom list for multipoles


  n14.resize(alloc_glob + 50);
  i14.resize(alloc_glob + 50);


  for (size_t n = 0; n < i14.size(); n++) {
    i14[n].resize(alloc_glob + 50);
  }



  for (i = 1; i <= n_atom; i++) {
    n14[i] = 0;
    for (size_t j = 1; j <= n13[i]; j++) {
      jj = i13[j][i];

      for (size_t k = 0; k < coords->atoms(jj - 1).bonds().size(); k++) {
        //       cout << k<<"k"<<endl;
        control_seq = false;
        kk = coords->atoms(jj - 1).bonds()[k] + 1;
        // 	cout << kk<< endl;
        if (kk == i) continue;
        for (size_t m = 0; m < coords->atoms(i - 1).bonds().size(); m++) {

          if (kk == coords->atoms(i - 1).bonds()[m] + 1) control_seq = true;
        }
        if (control_seq == true) continue;

        for (size_t m = 1; m <= n13[i]; m++) {
          if (kk == i13[m][i]) control_seq = true;
        }
        if (control_seq == true) continue;

        n14[i] = n14[i] + 1;

        i14[n14[i]][i] = kk;

        /*zwei:
          continue;*/
      }
    }
  }


  //!1-5 attached atom list for multipoles

  n15.resize(alloc_glob + 50);
  i15.resize(alloc_glob + 50);

  for (size_t n = 0; n < i15.size(); n++) {
    i15[n].resize(alloc_glob + 50);
  }


  control_seq = false;
  for (i = 1; i <= n_atom; i++) {
    n15[i] = 0;
    for (size_t j = 1; j <= n14[i]; j++) {
      jj = i14[j][i];
      for (size_t k = 0; k < coords->atoms(jj - 1).bonds().size(); k++) {
        kk = coords->atoms(jj - 1).bonds()[k] + 1;
        control_seq = false;
        if (kk == i) continue;
        for (size_t m = 0; m < coords->atoms(i - 1).bonds().size(); m++) {

          if (kk == coords->atoms(i - 1).bonds()[m] + 1) control_seq = true;
        }
        if (control_seq == true) continue;

        for (size_t m = 1; m <= n13[i]; m++) {
          if (kk == i13[m][i]) control_seq = true;

        }
        if (control_seq == true) continue;

        for (size_t m = 1; m <= n14[i]; m++) {
          if (kk == i14[m][i]) control_seq = true;
        }
        if (control_seq == true) continue;

        n15[i] = n15[i] + 1;
        i15[n15[i]][i] = kk;


      }
    }
  }


  //!polargrp




  for (size_t n = 0; n < ip11.size(); n++) {
    ip11[n].resize(alloc_glob + 50);
  }


  //! assign direct list
  maxval = 8;

  for (i = 1; i <= n_atom; i++) {
    np11[i] = 1;
    ip11[1][i] = i;
    it = coords->atoms(i - 1).energy_type();


    for (size_t j = 0; j < coords->atoms(i - 1).bonds().size(); j++) {
      jj = coords->atoms(i - 1).bonds()[j] + 1;
      jt = coords->atoms(jj - 1).energy_type();
      /*std::cout << "jj " << jj << " jt" << jt << '\n';*/

      for (size_t k = 0; k < 3; k++) {

        //kk = plrgrp[k][it-1];

        //std::cout << plrgrp[k][i-1] << '\n';
        if (plrgrp[k][i - 1] == 0) break;
        if (plrgrp[k][i - 1] == jt) {
          np11[i]++;
          //std::cout << np11[i] << endl;
          if (np11[i] <= 100) {
            ip11[np11[i]][i] = jj;


          }
          else { std::cout << "TOO MANY ATOMS IN POLARIZATION GROUP\n"; }

        }
      }

    }
  }




  //!find any group members

  for (i = 1; i <= n_atom; i++) {
    list[i] = 0;
  }

  for (i = 1; i <= n_atom; i++) {

    done = false;
    start = 1;
    stop = np11[i];
    //std::cout << np11[i] << '\n';
    for (size_t j = start; j <= stop; j++) {
      jj = ip11[j][i];

      if (jj < i) {
        done = true;
        np11[i] = np11[jj];

        for (size_t k = 1; k <= np11[i]; k++) {
          ip11[k][i] = ip11[k][jj];
          //std::cout << ip11[k][i] << '\n';
        }
      }
      else { list[jj] = i; }
    }
    while (!done) {
      done = true;
      for (size_t j = start; j <= stop; j++) {
        jj = ip11[j][i];
        //         cout << j<< endl;
        for (size_t k = 1; k <= np11[jj]; k++) {
          kk = ip11[k][jj];

          if (list[kk] != i) {
            np11[i]++;
            if (np11[i] <= 100) {
              ip11[np11[i]][i] = kk;

            }
            list[kk] = i;
          }
        }
      }

      if (np11[i] != stop) {
        done = false;
        start = stop + 1;
        stop = np11[i];
      }
    }
  }






  //! assign 1-2 relations 
  np12.resize(alloc_glob + 50);
  ip12.resize(alloc_glob + 50);
  for (size_t n = 0; n < ip12.size(); n++) {
    ip12[n].resize(alloc_glob);
  }


  for (i = 1; i <= n_atom; i++) {
    mask[i] = 0;
  }
  for (i = 1; i <= n_atom; i++) {
    nlist = 0;
    for (size_t j = 1; j <= np11[i]; j++) {
      jj = ip11[j][i];
      //cout << j<< endl;
      nlist++;

      list[nlist] = jj;
      mask[jj] = i;
    }
    nkeep = 0;
    for (size_t j = 1; j <= nlist; j++) {
      jj = list[j];
      //cout << nlist << endl;
      for (size_t k = 0; k < coords->atoms(jj - 1).bonds().size(); k++) {
        kk = coords->atoms(jj - 1).bonds()[k] + 1;

        //std::cout << kk << '\n';
        if (mask[kk] != i) {
          nkeep++;
          keep[nkeep] = kk;
        }
      }
    }
    nlist = 0;

    for (size_t j = 1; j <= nkeep; j++) {
      jj = keep[j];

      for (size_t k = 1; k <= np11[jj]; k++) {
        kk = ip11[k][jj];
        nlist++;
        list[nlist] = kk;
      }
    }

    if (nlist <= 50) {
      np12[i] = nlist;
      for (size_t j = 1; j <= nlist; j++) {
        ip12[j][i] = list[j];
        //std::cout << ip12[j][i] << '\n';

      }


    }
    else { break; }
  }
  //!1-3relations

  np13.resize(alloc_glob + 50);
  ip13.resize(alloc_glob + 50);

  for (size_t n = 0; n < ip13.size(); n++) {
    ip13[n].resize(alloc_glob);
  }

  //!!!!!

  for (i = 1; i <= n_atom; i++) {
    mask[i] = 0;
  }
  for (i = 1; i <= n_atom; i++) {
    for (size_t j = 1; j <= np11[i]; j++) {
      jj = ip11[j][i];
      mask[jj] = i;
    }
    for (size_t j = 1; j <= np12[i]; j++) {
      jj = ip12[j][i];
      mask[jj] = i;
    }


    nlist = 0;
    for (size_t j = 1; j <= np12[i]; j++) {
      kk = ip12[j][jj];
      for (size_t k = 1; k <= np12[jj]; k++) {
        if (mask[kk] != i) {
          nlist++;
          list[nlist] = kk;
        }
      }
    }

    if (nlist <= 50) {
      np13[i] = nlist;
      for (size_t j = 1; j <= nlist; j++) {
        ip13[j][i] = list[j];

      }
    }
    else { break; }
  }


  //!1-4relations

  np14.resize(alloc_glob + 50);
  ip14.resize(alloc_glob + 50);

  for (size_t n = 0; n < ip14.size(); n++) {
    ip14[n].resize(alloc_glob + 50);
  }


  for (i = 1; i <= n_atom; i++) {
    mask[i] = 0;
  }
  for (i = 1; i <= n_atom; i++) {
    for (size_t j = 1; j <= np11[i]; j++) {
      jj = ip11[j][i];
      mask[jj] = i;
    }
    for (size_t j = 1; j <= np12[i]; j++) {
      jj = ip12[j][i];
      mask[jj] = i;
    }
    for (size_t j = 1; j <= np13[i]; j++) {
      jj = ip13[j][i];
      mask[jj] = i;
    }
    nlist = 0;
    for (size_t j = 1; j <= np13[i]; j++) {
      kk = ip13[j][jj];
      for (size_t k = 1; k <= np12[jj]; k++) {
        if (mask[kk] != i) {
          nlist++;
          list[nlist] = kk;
        }
      }
    }
    if (nlist <= 50) {
      np14[i] = nlist;
      for (size_t j = 1; j <= nlist; j++) {
        ip14[j][i] = list[j];
      }
    }
    else { break; }
  }



}
void energy::interfaces::amoeba::amoeba_ff::rot_matrix(coords::Representation_3D const &pos)
{
  scon::c3<scon::c3<double>> am;
  std::vector <std::vector<double> > m2, r2;
  std::array<std::array<double, 3>, 3> am2;
  size_t n(1U);





  rp.resize(14);

  for (size_t i = 0; i < rp.size(); i++) {
    rp[i].resize(alloc_glob + 1);
  }

  m2.resize(14);
  for (size_t i = 0; i < m2.size(); i++) {
    m2[i].resize(alloc_glob + 1);
  }
  r2.resize(14);
  for (size_t i = 0; i < r2.size(); i++) {
    r2[i].resize(alloc_glob + 1);
  }

  am2[0][0] = 1.0;
  am2[1][0] = 0.0;
  am2[2][0] = 0.0;
  am2[0][2] = 0.0;
  am2[1][2] = 0.0;
  am2[2][2] = 1.0;

  using tinker::parameter::multipole;

  for (auto axes : refined.multipole_vecs())
  {
    for (auto mult : axes)
    {

      //std::cout <<pos[mult.center].x() << "   " << pos[mult.center].y() << "   " <<pos[mult.center].z() << '\n';
      //std::cout << mult.p_rot.axt << '\n';

      switch (mult.p_rot.axt)
      {

        case multipole::axtype::NONE: // 0
        {
          am.x().x() = 1.0;
          am.z().z() = 1.0;
          break;
        }

        case multipole::axtype::Z_AXIS: // 1
        {
          am.z() = normalized(pos[mult.axes.z()] - pos[mult.center]);
          coords::Cartesian_Point d = scon::randomized<coords::Cartesian_Point>();
          d -= am.z() * dot(d, am.z());
          am.x() = normalized(d);
          break;
        }

        case multipole::axtype::Z_THEN_X: // 2
        {
          am.z() = normalized(pos[mult.axes.z()] - pos[mult.center]);
          coords::Cartesian_Point d(pos[mult.axes.x()] - pos[mult.center]);
          d -= am.z()*dot(d, am.z());
          am.x() = normalized(d);
          break;
        }

        case multipole::axtype::BISECTOR: // 3
        {
          coords::Cartesian_Point const d1(normalized(pos[mult.axes.z()] - pos[mult.center]));
          coords::Cartesian_Point const d2(normalized(pos[mult.axes.x()] - pos[mult.center]));
          am.z() = normalized(d1 + d2);
          am.x() = normalized(d2 - (am.z() * scon::dot(d2, am.z())));
          break;
        }

        case multipole::axtype::Z_BISECTOR: // 4
        {
          am.z() = normalized(pos[mult.axes.z()] - pos[mult.center]);
          coords::Cartesian_Point const d1(normalized(pos[mult.axes.x()] - pos[mult.center]));
          coords::Cartesian_Point const d2(normalized(pos[mult.axes.y()] - pos[mult.center]));
          coords::Cartesian_Point d(normalized(d1 + d2));
          d -= am.z()*dot(d, am.z());
          am.x() = normalized(d);
          break;
        }

        case multipole::axtype::THREEFOLD: // 5
        {
          coords::Cartesian_Point const d1(normalized(pos[mult.axes.z()] - pos[mult.center]));
          coords::Cartesian_Point const d2(normalized(pos[mult.axes.x()] - pos[mult.center]));
          coords::Cartesian_Point const d3(normalized(pos[mult.axes.y()] - pos[mult.center]));
          am.z() = normalized(d1 + d2 + d3);
          am.x() = normalized(d2 - am.z()*scon::dot(d2, am.z()));
          break;

        }

      }

      am.y() = cross(am.z(), am.x());

      am2[0][0] = am.x().x();
      am2[0][1] = am.y().x();
      am2[0][2] = am.z().x();
      am2[1][0] = am.x().y();
      am2[1][1] = am.y().y();
      am2[1][2] = am.z().y();
      am2[2][0] = am.x().z();
      am2[2][1] = am.y().z();
      am2[2][2] = am.z().z();



      //ROTATION of MONOPOLES


      rp[0][n] = charges[n];



      //ROTATION of DIPOLES


      for (size_t i = 1; i < 4; i++)
      {
        rp[i][n] = 0.0;
        for (size_t j = 0; j < 3; j++)
        {
          rp[i][n] = rp[i][n] + dipole[n][j] * am2[i - 1][j];

        }
      }





      //ROTATION OF QUADRUPOLES	


      size_t k = 0;

      for (size_t i = 1; i < 4; i++) {
        for (size_t j = 1; j < 4; j++) {
          m2[i][j] = quadro[n][k];
          r2[i][j] = 0.0;

          k++;



        }
      }



      for (size_t i = 1; i < 4; i++) {
        for (size_t j = 1; j < 4; j++) {

          if (j < i) {
            r2[i][j] = r2[j][i];

          }
          else {
            for (k = 1; k < 4; k++) {

              for (size_t m = 1; m < 4; m++) {

                r2[i][j] = r2[i][j] + (am2[i - 1][k - 1] * am2[j - 1][m - 1] * m2[k][m]);
                //std::cout << am2[i - 1][k - 1] << '\n';
              }
            }
          }
        }
      }
      k = 4;
      for (size_t i = 1; i < 4; i++) {
        for (size_t j = 1; j < 4; j++) {
          rp[k][n] = r2[i][j];
          k++;
        }
      }

      //}
      n++;

    }
  }



}

//ENERGY and GRADIENT FUNCTIONS for short-range correction
//********************************************************
//********************************************************


//DEFINITION OF MOLECULAR UNITS

void energy::interfaces::amoeba::amoeba_ff::Spackman_mol() {
  const size_t N = (coords->size());;



  /* alloc = 25;*/

  mol.resize(N + 1);

  size_t i, j, k, ii, mi, mj, mk;
  size_t nmol;

  nmol = 0;

  for (i = 0; i < N; i++) {
    mol[i] = 0;
  }

  for (i = 0; i < N; i++) {

    if (mol[i] == 0) {
      nmol++;
      mol[i] = nmol;
    }
    mi = mol[i];

    for (ii = 0; ii < coords->atoms(i).bonds().size(); ii++) {

      j = coords->atoms(i).bonds()[ii];
      mj = mol[j];

      if (mj == 0) {
        mol[j] = mi;
      }
      else if (mi < mj) {
        nmol--;
        for (k = 0; k < N; k++) {
          mk = mol[k];
          if (mk == mj) {
            mol[k] = mi;
          }
          else if (mk > mj) {
            mol[k] = mk - 1;
          }

        }
      }
      else if (mi > mj) {
        nmol--;
        for (k = 0; k < N; k++) {
          mk = mol[k];
          if (mk == mi) {
            mol[k] = mj;
          }
          else if (mk > mi) {
            mol[k] = mk - 1;

          }
        }
        mi = mj;
      }
    }



  }




}


//DEFINITION OF MOLECULAR INTERACTIONVECTOR

void energy::interfaces::amoeba::amoeba_ff::Spackman_vec() {


  double dist_3;
  auto const &positions = coords->xyz();
  coords::Cartesian_Point bv, gv;
  spack_list spacktemp;
  vec_spack.clear();
  const size_t N = coords->xyz().size();


  cutoff_spackman = 10.0;
  for (size_t i = 0; i < N; i++) {

    for (size_t ii = i + 1; ii < N; ii++) {

      if (i == ii) continue;
      if (mol[i] == mol[ii]) continue;
      bv = positions[i] - positions[ii];

      dist_3 = len(bv);




      if (dist_3 > cutoff_spackman) continue;
      spacktemp.atom[0] = i;
      spacktemp.atom[1] = ii;

      vec_spack.push_back(spacktemp);

    }
  }
}

//FUNCTION FOR CALCULATING THE FACULTY

size_t energy::interfaces::amoeba::amoeba_ff::fak_iter1(size_t n)
{
  size_t k(0U);
  fak = 1;
  for (k = 1; k <= n; k++)
    fak *= k;

  return fak;
}


//FUNCTION FOR CALCULATING ESSENTIALS FOR SHORT_RANGE CORRECTION

double energy::interfaces::amoeba::amoeba_ff::fff_f1(size_t n, double s, double z)
{
  if (n == 2) { fff = 2.0*z / ((s*s + z*z)*(s*s + z*z)); }
  if (n == 3) { fff = 2.0*(3.0*(z*z) - (s*s)) / ((s*s + z*z)*(s*s + z*z)*(s*s + z*z)); }
  if (n == 4) { fff = 24.0*z*(z*z - s*s) / ((s*s + z*z)*(s*s + z*z)*(s*s + z*z)*(s*s + z*z)); }
  if (n == 5) { fff = 24.0*((5.0*(z*z*z*z)) - (10.0*((z*s)*(z*s))) + (s*s*s*s)) / ((s*s + z*z)*(s*s + z*z)*(s*s + z*z)*(s*s + z*z)*(s*s + z*z)); }
  if (n == 6) { fff = 240.0*z*((s*s) - (3.0*(z*z))*((3.0*(s*s)) - z*z)) / ((s*s + z*z)*(s*s + z*z)*(s*s + z*z)*(s*s + z*z)*(s*s + z*z)*(s*s + z*z)); }

  return fff;
}

void energy::interfaces::amoeba::amoeba_ff::parameters()
{


  char control;
  nijx.clear();
  zetax.clear();
  char buffer[200];
  std::ifstream nijf("spackman.prm", std::ios::in);


  while (!nijf.eof())
  {

    nijf.getline(buffer, 200);
    sscanf(buffer, "%c", &control);

    k kbuffer;
    hh hbuffer;
	oo oobuffer;
    nn nnbuffer;


	std::string forcefield;







    if (control == 'n') {
      std::vector <ptrdiff_t> tempn(19);
      sscanf(buffer, "%*s %td  %td  %td  %td  %td  %td  %td  %td  %td  %td  %td  %td  %td  %td  %td  %td  %td  %td  %td",
        &tempn[0], &tempn[1], &tempn[2], &tempn[3], &tempn[4], &tempn[5], &tempn[6], &tempn[7], &tempn[8], &tempn[9], &tempn[10], &tempn[11], &tempn[12], &tempn[13], &tempn[14], &tempn[15], &tempn[16], &tempn[17], &tempn[18]);
      nijx.push_back(tempn);
    }


    else if (control == 'z') {
      std::vector <float> tempx(19);

      sscanf(buffer, "%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
        &tempx[0], &tempx[1], &tempx[2], &tempx[3], &tempx[4], &tempx[5], &tempx[6], &tempx[7], &tempx[8], &tempx[9], &tempx[10], &tempx[11], &tempx[12], &tempx[13], &tempx[14], &tempx[15], &tempx[16], &tempx[17], &tempx[18]);


      zetax.push_back(tempx);

    }

    else if (control == 'c') {
      std::vector <float> tempc(19);
      sscanf(buffer, "%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
        &tempc[0], &tempc[1], &tempc[2], &tempc[3], &tempc[4], &tempc[5], &tempc[6], &tempc[7], &tempc[8], &tempc[9], &tempc[10], &tempc[11], &tempc[12], &tempc[13], &tempc[14], &tempc[15], &tempc[16], &tempc[17], &tempc[18]);
      //if (!cf.eof())
      cs.push_back(tempc);

    }

    else if (control == 'o') {
      std::vector <ptrdiff_t> tempo(3);
      sscanf(buffer, "%*s  %td  %td  %td ",
        &tempo[0], &tempo[1], &tempo[2]);
      occ.push_back(tempo);

    }


    else if (control == '6') {

      sscanf(buffer, "%*s %lf ", &kbuffer.temp5);
      kappan.push_back(kbuffer.temp5);

    }
	else if (control == '1') {

		sscanf(buffer, "%*s %lf ", &hbuffer.temp6);
		kappan.push_back(hbuffer.temp6);

	}
	else if (control == '5') {

		sscanf(buffer, "%*s %lf ", &nnbuffer.temp7);
		kappan.push_back(nnbuffer.temp7);
	}
	else if (control == '8') {

		sscanf(buffer, "%*s %lf ", &oobuffer.temp8);
		kappan.push_back(oobuffer.temp8);

	}
    else if (control == 'e') break;
  }






}
void energy::interfaces::amoeba::amoeba_ff::Spackman_list_analytical1() {




  //size_t i,j,i1,j1,j2,k,it,ia,n,kk;
  //		  ptrdiff_t nbas=19;
  //		  ptrdiff_t nat=10;
  //		  ptrdiff_t nref=7;
  //double bohr=0.52917720859;
  //double pi=3.141592653589793238;
  //double coulomb=332.063709;
  //		  std::vector <std::vector<double> > coef;
  //		  std::vector <std::vector<std::vector<double> > > pij;
  //		  std::vector <double> sscat;
  //		  double aLimes,bLimes,fak,smax;
  //double sum1, sum2, tnm,del,ddel,sum,kap,f,argBessel0,e,q,qsum=0,esum=0;
  //double zz,pp;
  //		  //ptrdiff_t n_atom=24;
  //		 // vector <double> kappa/*atomic(24)*/;
  //std::vector <size_t>  atomic;		
  //		  double distx=0, disty=0, distz=0, dist=0, dist_3=0;
  //std::vector <double> kappa;
  //size_t n_atom=2;

  //const scon::c3<double> &positions = coords->xyz();




  //double dextemp1;
  //Coord::vect gradient;
  //std::ofstream fout1("grad1.in");
  //std::ofstream fout2("grad2.in");
  //std::ofstream fout3("grad3.in");
  //std::ofstream fout4("grad4.in");
  //std::ofstream fout5("grad5.in");
  //std::ofstream fout6("grad6.in");
  /* size_t kappamax=2;  */
  double dist;
  int cut = 11002;
  //double cut=(10.0*1000)+1001;  

  exa11.resize(cut);
  exa22.resize(cut);
  exa33.resize(cut);
  eveca1.resize(cut);
  dex11.resize(cut);
  dex22.resize(cut);
  dex33.resize(cut);

  exa44.resize(cut);
  exa55.resize(cut);
  exa66.resize(cut);
  dex44.resize(cut);
  dex55.resize(cut);
  dex66.resize(cut);

  exa77.resize(cut);
  exa88.resize(cut);
  exa99.resize(cut);
  exa1010.resize(cut);
  dex77.resize(cut);
  dex88.resize(cut);
  dex99.resize(cut);
  dex1010.resize(cut);


  std::ifstream in1("HH_GRAD.in", std::ios::in);
  std::ifstream in2("CC_GRAD.in", std::ios::in);
  std::ifstream in3("CH_GRAD.in", std::ios::in);
  std::ifstream in4("HH_EN.in", std::ios::in);
  std::ifstream in5("CC_EN.in", std::ios::in);
  std::ifstream in6("CH_EN.in", std::ios::in);

  std::ifstream in7("OH_GRAD.in", std::ios::in);
  std::ifstream in8("OC_GRAD.in", std::ios::in);
  std::ifstream in9("OO_GRAD.in", std::ios::in);
  std::ifstream in10("OH_EN.in", std::ios::in);
  std::ifstream in11("OC_EN.in", std::ios::in);
  std::ifstream in12("OO_EN.in", std::ios::in);

  std::ifstream in13("NH_GRAD.in", std::ios::in);
  std::ifstream in14("NC_GRAD.in", std::ios::in);
  std::ifstream in15("NN_GRAD.in", std::ios::in);
  std::ifstream in16("NO_GRAD.in", std::ios::in);
  std::ifstream in17("NH_EN.in", std::ios::in);
  std::ifstream in18("NC_EN.in", std::ios::in);
  std::ifstream in19("NN_EN.in", std::ios::in);
  std::ifstream in20("NO_EN.in", std::ios::in);

  char buffer[200];
  double temp;

 size_t i = 0;
  while (!in1.eof())
  {
    in1.getline(buffer, 200);
    sscanf(buffer, "%lf", &temp);
    exa11[i] = temp;
    i++;
  }
  i = 0;
  while (!in2.eof())
  {
    in2.getline(buffer, 200);
    sscanf(buffer, "%lf", &temp);
    exa22[i] = temp;
    i++;
  }
  i = 0;
  while (!in3.eof())
  {
    in3.getline(buffer, 200);
    sscanf(buffer, "%lf", &temp);
    exa33[i] = temp;
    i++;
  }
  i = 0;
  while (!in4.eof())
  {
    in4.getline(buffer, 200);
    sscanf(buffer, "%lf", &temp);
    dex11[i] = temp;

    i++;
  }
  i = 0;
  while (!in5.eof())
  {
    in5.getline(buffer, 200);
    sscanf(buffer, "%lf", &temp);
    dex22[i] = temp;
    i++;
  }
  i = 0;
  while (!in6.eof())
  {
    in6.getline(buffer, 200);
    sscanf(buffer, "%lf", &temp);
    dex33[i] = temp;
    i++;
  }
  i = 0;
  while (!in7.eof())
  {
	  in7.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  exa44[i] = temp;
	  i++;
  }
  i = 0;
  while (!in8.eof())
  {
	  in8.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  exa55[i] = temp;
	  i++;
  }
  i = 0;
  while (!in9.eof())
  {
	  in9.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  exa66[i] = temp;
	  i++;
  }
  i = 0;
  while (!in10.eof())
  {
	  in10.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  dex44[i] = temp;

	  i++;
  }
  i = 0;
  while (!in11.eof())
  {
	  in11.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  dex55[i] = temp;

	  i++;
  }
  i = 0;
  while (!in12.eof())
  {
	  in12.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  dex66[i] = temp;

	  i++;
  }
  i = 0;



  while (!in13.eof())
  {
	  in13.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  exa77[i] = temp;
	  i++;
  }
  i = 0;
  while (!in14.eof())
  {
	  in14.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  exa88[i] = temp;
	  i++;
  }
  i = 0;
  while (!in15.eof())
  {
	  in15.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  exa99[i] = temp;
	  i++;
  }
  i = 0;
  while (!in16.eof())
  {
	  in16.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  exa1010[i] = temp;
	  i++;
  }
  i = 0;


  while (!in17.eof())
  {
	  in17.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  dex77[i] = temp;

	  i++;
  }
  i = 0;
  while (!in18.eof())
  {
	  in18.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  dex88[i] = temp;

	  i++;
  }
  i = 0;
  while (!in19.eof())
  {
	  in19.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  dex99[i] = temp;

	  i++;
  }
  i = 0;
  while (!in20.eof())
  {
	  in20.getline(buffer, 200);
	  sscanf(buffer, "%lf", &temp);
	  dex1010[i] = temp;

	  i++;
  }





  for (i = 1; i <= 11002; i++) {
    dist = i*0.001;
    eveca1[i-1] = dist;

  }

  // 	//	!Initialisierung	

  //
  // //!Initialsisierung	
  // atomic.clear();
  // kappa.clear();
  // pij.resize(20);
  // for(i=0;i<pij.size();i++){
  //   pij[i].resize(20);
  //   for(j=0;j<pij[i].size();j++){
  //     pij[i][j].resize(20);

  //   }
  // }




  // coef.resize(400);
  // for(i=0;i<coef.size();i++){
  //   coef[i].resize(400);
  // }
  // fscat.resize(n_atom+1000);
  // for(i=0;i<fscat.size();i++){
  //   fscat[i].resize(n_atom+1000);
  // }
  // sscat.resize(730);
  // kappa.resize(n_atom+1);
  // energy_short =4.0;





  // smax=10.0;
  // aLimes=0.0;
  // bLimes=smax*4*pi;





  // //! calculate the coefficients for scattering vector



  // for (i1=0;i1<nat;i1++){

  //   if (i1 == 9) continue;
  //   if((i1 > 3 ) && (i1 < 5)) continue;

  //   for(j=0;j<nbas;j++){

  //     if(nijx[i1][j] == 0) continue;
  //     coef[i1][j]=(pow(2.0*zetax[i1][j]/bohr,(int)nijx[i1][j])*sqrt(2.0*zetax[i1][j]/bohr)/sqrt(double(fak_iter1(2*nijx[i1][j]))));

  //   }
  // } 

  // //!atomic density matrix

  // for (i1=0;i1<nat;i1++){
  //   if (i1 == 9) continue;
  //   if((i1 > 3) && (i1 < 5)) continue;
  //   //!1s	
  //   for(i=0;i<7;i++){
  //     for(j=i;j<7;j++){
  //       pij[i][j][i1]+=occ[i1][0]*cs[i1][i]*cs[i1][j];



  //     }
  //   }
  //   //!2s				
  //   for(i=7;i<14;i++){
  //     for(j=i;j<14;j++){
  //       pij[i][j][i1]+=occ[i1][1]*cs[i1][i]*cs[i1][j];


  //     }
  //   }	

  //   for(i=14;i<19;i++){
  //     for(j=i;j<19;j++){
  //       pij[i][j][i1]+=occ[i1][2]*cs[i1][i]*cs[i1][j];
  //     }
  //   }
  // }



  // //!get integration points

  // kp=1;
  // sscat[kp]=0.50*(aLimes+bLimes);
  // for (i=2;i<=7;i++){
  //   it=(size_t)pow(3.0,double(i-2));
  //   tnm=it;
  //   del=(bLimes-aLimes)/(3.0*tnm);
  //   ddel=del+del;
  //   kp++;
  //   sscat[kp]=aLimes+0.5*del;
  //   for (j=1;j<=it;j++){
  //     kp++;
  //     sscat[kp]=sscat[kp-1]+ddel;
  //     if (j == it) continue;
  //     kp++;

  //     sscat[kp]=sscat[kp-1]+del;
  //   }
  // }









  //
  // for( size_t nn=0;nn<=2;nn++){

  //   /* n=0;*/  


  //   for (size_t ii=1; ii<cut; ii++){  
  //     kappa.clear();
  //     atomic.clear();
  //     dextemp1=0.0;
  //     esum=0.0;
  //     qsum=0.0;
  //   //  !Definition of atomic parameters
  //   //  !H-H
  //     if (nn == 0 ){ for (i=1;i<2;i++){
  //       atomic.push_back(1);
  //       kappa.push_back(kappan[1]); 
  //     }
  //     for (i=1;i<2;i++){
  //       		
  //       atomic.push_back(7);
  //       kappa.push_back(1.0);
  //     }


  //     }
  //   //  !C-C
  //     else if (nn == 1 ){ for (i=1;i<2;i++){
  //       atomic.push_back(7);
  //       kappa.push_back(1.0);
  //     }
  //     for (i=1;i<2;i++){
  //       	
  //       atomic.push_back(8);
  //       kappa.push_back(1.0);
  //     }



  //     }




  //   //  !c-H
  //     else if (nn == 2 ){ for (i=1;i<2;i++){
  //       atomic.push_back(7);
  //       kappa.push_back(1.0); 
  //     }
  //     for (i=1;i<2;i++){
  //   
  //       atomic.push_back(6);
  //       kappa.push_back(kappan[0]);
  //     }


  //     }	

  //

  // //   //// !definition of further interactions here

  // //   //  else if (nn == 3 ){ for (i=1;i<2;i++){
  // //   //  				   atomic.push_back(X);
  // //   //  				   kappa.push_back(X); 
  // //   //  				  }
  // //   //  			for (i=1;i<2;i++){
  // //   //  				  // cout<<"Hallo3"<<endl;
  // //   //  		                 atomic.push_back(X);
  // //   //  				 kappa.push_back(X);
  // //   //  				}
  // //   //// !......

  //     for( ia=0; ia <2 ; ia++){
  //       i1=atomic[ia];
  //       kap=1.0;
  //       //!Core
  //       if( (i1 == 1) && (kappa[ia] != 0.0)) kap=kappa[ia]; 			
  //       for (i=0;i<this->kp;i++){
  //         sum1=0.0;
  //         for(j1=0;j1<7;j1++){
  //           for(j2=j1;j2<7;j2++){
  //             zz=(zetax[i1-1][j1]+zetax[i1-1][j2])/bohr;

  //             k=nijx[i1-1][j1]+nijx[i1-1][j2];
  //             /* kk++;*/		

  //             pp=pij[j1][j2][i1-1];

  //             if(j1 != j2) {pp=2.0*pp;}
  //             sum1+= pp*coef[i1-1][j1]*coef[i1-1][j2]*(fff_f1(k,(this->sscat[i+1]/kap),zz)); //!function

  //           }
  //         }
  //         fscat[ia][i]=sum1;

  //       }

  //       if (i1 == 1) continue;


  //       kap=1.0;
  //       if (kappa[ia] != 0.0) kap= kappa[ia];
  //       //!Valence			   
  //       for (i=0;i<this->kp;i++){
  //         sum2=0.0;
  //         for(j1=7;j1<nbas;j1++){
  //           for(j2=j1;j2<nbas;j2++){
  //             zz=(zetax[i1-1][j1]+zetax[i1-1][j2])/bohr;
  //             k=nijx[i1-1][j1]+nijx[i1-1][j2];
  //             pp=pij[j1][j2][i1-1];
  //             if(j1 != j2) {pp=2.0*pp;}
  //             sum2+= pp*coef[i1-1][j1]*coef[i1-1][j2]*(fff_f1(k,(this->sscat[i+1]/kap),zz)); //!function
  //           }
  //         }
  //         fscat[ia][i]+=sum2;	

  //       }

  //     }



  //     dist=(ii*0.001);
  //     eveca1[ii]=((ii)*0.001);

  //     n_atom=1;


  //     for(size_t i=1; i<2;i++){
  //       kk=2;

  //       k=1;

  //       // 		
  //       f=fscat[i-1][k-1]*fscat[kk-1][k-1]-atomic[i-1]*fscat[kk-1][k-1]-atomic[kk-1]*fscat[i-1][k-1];

  //	
  //       argBessel0=sscat[k]*dist;

  //      
  //       f=f*((argBessel0*cos(argBessel0)-sin(argBessel0))/(sscat[k]*(dist*dist)));
  //       e=(bLimes-aLimes)*f;

  //       for(j=2; j <=nref; j++){
  //         it=(size_t)pow(3.0,double(j-2));
  //         // cout<<"Test2"<<endl;
  //         tnm=it;

  //         sum=0.0;

  //         for(n=1;n<=it;n++){
  //           k++;
  //           f=fscat[i-1][k-1]*fscat[kk-1][k-1]-atomic[i-1]*fscat[kk-1][k-1]-atomic[kk-1]*fscat[i-1][k-1];
  //           argBessel0=sscat[k]*dist;

  //           f=f*((argBessel0*cos(argBessel0)-sin(argBessel0))/(sscat[k]*(dist*dist)));

  //           sum=sum+f;

  //           k++;
  //           f=fscat[i-1][k-1]*fscat[kk-1][k-1]-atomic[i-1]*fscat[kk-1][k-1]-atomic[kk-1]*fscat[i-1][k-1];
  //           argBessel0=sscat[k]*dist;
  //           f=f*((argBessel0*cos(argBessel0)-sin(argBessel0))/(sscat[k]*(dist*dist)));

  //           sum=sum+f;

  //         }
  //         e=(e+(bLimes-aLimes)*sum/tnm)/3.0;
  //         //cout<<e<<endl;
  //       }
  //       q=-double(atomic[i-1]*atomic[kk-1])/(dist*dist);

  //       esum=esum+e;

  //       qsum=qsum+q;   
  //       // cout<<esum<<endl;  



  //       e=(qsum+((2.0*esum)/pi))*coulomb;






  //       if ( nn == 0) {
  //         exa11[ii]=e;
  //         fout1<<std::setprecision(16)<<exa11[ii]<<std::endl;
  //       }

  //       else if ( nn == 1) {
  //         exa22[ii]=e;
  //	  fout2<<std::setprecision(16)<<exa22[ii]<<std::endl;
  //       }
  //       else if ( nn == 2 ) {
  //         exa33[ii]=e;
  //	  fout3<<std::setprecision(16)<<exa33[ii]<<std::endl;
  //       }
  //  }
  //  
  ////  }
  //  

  //  n_atom=1;
  //        dist=(ii*0.001);
  //     eveca1[ii]=((ii)*0.001);

  //	  for(size_t i=1; i<2;i++){
  //       kk=2;

  //       k=1;

  //       // 		
  //       f=fscat[i-1][k-1]*fscat[kk-1][k-1]-atomic[i-1]*fscat[kk-1][k-1]-atomic[kk-1]*fscat[i-1][k-1];

  //	
  //       argBessel0=sscat[k]*dist;

  //      
  //      f=f*sin(argBessel0)/argBessel0;
  //       e=(bLimes-aLimes)*f;

  //       for(j=2; j <=nref; j++){
  //         it=(size_t)pow(3.0,double(j-2));
  //         // cout<<"Test2"<<endl;
  //         tnm=it;

  //         sum=0.0;

  //         for(n=1;n<=it;n++){
  //           k++;
  //           f=fscat[i-1][k-1]*fscat[kk-1][k-1]-atomic[i-1]*fscat[kk-1][k-1]-atomic[kk-1]*fscat[i-1][k-1];
  //           argBessel0=sscat[k]*dist;

  //           f=f*sin(argBessel0)/argBessel0;

  //           sum=sum+f;

  //           k++;
  //           f=fscat[i-1][k-1]*fscat[kk-1][k-1]-atomic[i-1]*fscat[kk-1][k-1]-atomic[kk-1]*fscat[i-1][k-1];
  //           argBessel0=sscat[k]*dist;
  //            f=f*sin(argBessel0)/argBessel0;

  //           sum=sum+f;

  //         }
  //         e=(e+(bLimes-aLimes)*sum/tnm)/3.0;
  //         //cout<<e<<endl;
  //       }
  //       q=double(atomic[i-1]*atomic[kk-1])/dist;
  //	
  //       esum=esum+e;

  //       qsum=qsum+q;   
  //       // cout<<esum<<endl;  

  //	  

  //       e=(qsum+((2.0*esum)/pi))*coulomb;

  //	  if ( nn == 0) {
  //         dex11[ii]=e;
  //        fout4<<std::setprecision(16)<<dex11[ii]<<std::endl;
  //       }

  //       else if ( nn == 1) {
  //         dex22[ii]=e;
  //       fout5<<std::setprecision(16)<<dex22[ii]<<std::endl;
  //       }
  //       else if ( nn == 2 ) {
  //         dex33[ii]=e;
  //	fout6<<std::setprecision(16)<<dex33[ii]<<std::endl;
  //       }
  //  
  //	  }

  //       }
  //	
  //	

  //   }

  SPACKrefine = true;
}


void energy::interfaces::amoeba::amoeba_ff::Spackman1() {




  const size_t
    n_atom = coords->xyz().size();


  size_t i, j, i1, j1, j2, k, it, ia;
  nbas = 19;
  nref = 7;
  nat = 10;
  double bohr = 0.52917720859;
  double pi = 3.141592653589793238;
  /* double coulomb=332.063709;
  double distx=0, disty=0, distz=0, dist_2=0, dist=0, dist_3=0;*/
  double  sum2, tnm, del, ddel, kap, sum1;
  double zz, pp;
  coords::Cartesian_Point bv, b, gv;
  //!temporary zeta




  //!Initialsisierung	
  atomic.clear();
  kappa.clear();
  pij.resize(200);
  for (i = 0; i < pij.size(); i++) {
    pij[i].resize(200);
    for (j = 0; j < pij[i].size(); j++) {
      pij[i][j].resize(200);

    }
  }




  coef.resize(400);
  for (i = 0; i < coef.size(); i++) {
    coef[i].resize(400);
  }
  fscat.resize(n_atom + 1000);
  for (i = 0; i < fscat.size(); i++) {
    fscat[i].resize(n_atom + 1000);
  }
  sscat.resize(730);
  kappa.resize(n_atom + 1);





  smax = 10.0;
  aLimes = 0.0;
  bLimes = smax * 4 * pi;





  //! calculate the coefficients for scattering vector



  for (i1 = 0; i1 < nat; i1++) {

    if (i1 == 9) continue;
    if ((i1 > 3) && (i1 < 5)) continue;

    for (j = 0; j < nbas; j++) {

      if (nijx[i1][j] == 0) continue;
      coef[i1][j] = (pow(2.0*zetax[i1][j] / bohr, (int)nijx[i1][j])*sqrt(2.0*zetax[i1][j] / bohr) / sqrt(double(fak_iter1(2 * nijx[i1][j]))));

    }
  }
  //!atomic density matrix

  for (i1 = 0; i1 < nat; i1++) {
    if (i1 == 9) continue;
    if ((i1 > 3) && (i1 < 5)) continue;
    //!1s	
    for (i = 0; i < 7; i++) {
      for (j = i; j < 7; j++) {
        pij[i][j][i1] += occ[i1][0] * cs[i1][i] * cs[i1][j];



      }
    }
    //!2s				
    for (i = 7; i < 14; i++) {
      for (j = i; j < 14; j++) {
        pij[i][j][i1] += occ[i1][1] * cs[i1][i] * cs[i1][j];


      }
    }

    for (i = 14; i < 19; i++) {
      for (j = i; j < 19; j++) {
        pij[i][j][i1] += occ[i1][2] * cs[i1][i] * cs[i1][j];
      }
    }
  }

  //!get integration points
  kp = 1;
  sscat[kp] = 0.50*(aLimes + bLimes);
  for (i = 2; i <= 7; i++) {
    it = (size_t)pow(3.0, double(i - 2));
    tnm = double(it);
    del = (bLimes - aLimes) / (3.0*tnm);
    ddel = del + del;
    kp++;
    sscat[kp] = aLimes + 0.5*del;
    for (j = 1; j <= it; j++) {
      kp++;
      sscat[kp] = sscat[kp - 1] + ddel;
      if (j == it) continue;
      kp++;

      sscat[kp] = sscat[kp - 1] + del;
    }
  }


  //! scattering factors for all atoms


  for (i = 0; i < n_atom; i++) {
    atomic.push_back(coords->atoms(i).number());
    if (coords->atoms(i).symbol() == "h" || coords->atoms(i).symbol() == "H") kappa[i] = kappan[1];
    else if (coords->atoms(i).symbol() == "c" || coords->atoms(i).symbol() == "C") kappa[i] = kappan[0];
	else if (coords->atoms(i).symbol() == "n" || coords->atoms(i).symbol() == "N") kappa[i] = kappan[2];
	else if (coords->atoms(i).symbol() == "o" || coords->atoms(i).symbol() == "O") kappa[i] = kappan[3];
  }
  
  //!getting atomic charges
  for (i = 0; i < n_atom; i++) {



    if (coords->atoms(i).symbol() == "c" || coords->atoms(i).symbol() == "C") { atomic.push_back(6); }

    else if (coords->atoms(i).symbol() == "o" || coords->atoms(i).symbol() == "O") { atomic.push_back(8); }
    else if (coords->atoms(i).symbol() == "n" || coords->atoms(i).symbol() == "N") { atomic.push_back(5); }
    else if (coords->atoms(i).symbol() == "s" || coords->atoms(i).symbol() == "S") { atomic.push_back(16); }
    else if (coords->atoms(i).symbol() == "h" || coords->atoms(i).symbol() == "H") { atomic.push_back(1); }
    else if (coords->atoms(i).symbol() == "f" || coords->atoms(i).symbol() == "F") { atomic.push_back(7); }
    else if (coords->atoms(i).symbol() == "p" || coords->atoms(i).symbol() == "P") { atomic.push_back(15); }

    //!....

    else { atomic.push_back(1); }


  }


  for (ia = 0; ia < n_atom; ia++) {
    i1 = atomic[ia];

    kap = 1.0;
    //!Core

    if ((i1 == 1) && (kappa[ia] != 0.0)) kap = kappa[ia];

    for (i = 0; i < kp; i++) {
      sum1 = 0.0;
      for (j1 = 0; j1 < 7; j1++) {
        for (j2 = j1; j2 < 7; j2++) {
          zz = (zetax[i1 - 1][j1] + zetax[i1 - 1][j2]) / bohr;

          k = nijx[i1 - 1][j1] + nijx[i1 - 1][j2];
          /* kk++;*/

          pp = pij[j1][j2][i1 - 1];

          if (j1 != j2) { pp = 2.0*pp; }
          sum1 += pp*coef[i1 - 1][j1] * coef[i1 - 1][j2] * (fff_f1(k, (sscat[i + 1] / kap), zz)); //!function

        }
      }
      fscat[ia][i] = sum1;


    }

    if (i1 == 1) continue;


    kap = 1.0;
    if (kappa[ia] != 0.0) kap = kappa[ia];
    //!Valence			   
    for (i = 0; i < kp; i++) {
      sum2 = 0.0;
      for (j1 = 7; j1 < nbas; j1++) {
        for (j2 = j1; j2 < nbas; j2++) {
          zz = (zetax[i1 - 1][j1] + zetax[i1 - 1][j2]) / bohr;
          k = nijx[i1 - 1][j1] + nijx[i1 - 1][j2];
          //cout<<zz<<endl;
          pp = pij[j1][j2][i1 - 1];
          if (j1 != j2) { pp = 2.0*pp; }
          sum2 += pp*coef[i1 - 1][j1] * coef[i1 - 1][j2] * (fff_f1(k, (sscat[i + 1] / kap), zz)); //!function
        }
      }
      fscat[ia][i] += sum2;

    }

  }


}
//
//
double energy::interfaces::amoeba::amoeba_ff::Spackman_Energy() {

  double r_cutoff;
  r_cutoff = Config::get().energy.spackman.cut;
  const size_t
    n_atom = coords->xyz().size();


  size_t j, k, it, n;
  nbas = 19;
  nref = 7;
  nat = 10;
  double coulomb = 332.063709;
  double  dist_2 = 0;
  double  tnm, sum, f, argBessel0, en, q, qsum = 0, esum = 0;
  auto const &positions = coords->xyz();
  coords::Cartesian_Point bv, b, gv;

  for (size_t i = 0; i < vec_spack.size(); i++) {

    bv = positions[vec_spack[i].atom[0]] - positions[vec_spack[i].atom[1]];

    dist_2 = len(bv);


    //	
    //	
    if (dist_2 > r_cutoff) continue;
    //  if(dist_2 <= 3.0) continue;
    k = 1;
    f = fscat[vec_spack[i].atom[0]][k - 1] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[0]] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[1]] * fscat[vec_spack[i].atom[0]][k - 1];

    argBessel0 = sscat[k] * dist_2;

    f = f*std::sin(argBessel0) / argBessel0;
    en = (bLimes - aLimes)*f;

    for (j = 2; j <= nref; j++) {
      it = (size_t)pow(3.0, double(j - 2));

      tnm = double(it);

      sum = 0.0;

      for (n = 1; n <= it; n++) {
        k++;
        f = fscat[vec_spack[i].atom[0]][k - 1] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[0]] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[1]] * fscat[vec_spack[i].atom[0]][k - 1];
        argBessel0 = sscat[k] * dist_2;
        f = f*std::sin(argBessel0) / argBessel0;
        sum += f;
        k++;
        f = fscat[vec_spack[i].atom[0]][k - 1] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[0]] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[1]] * fscat[vec_spack[i].atom[0]][k - 1];
        argBessel0 = sscat[k] * dist_2;
        f = f*std::sin(argBessel0) / argBessel0;

        sum += f;

      }
      en = (en + (bLimes - aLimes)*sum / tnm) / 3.0;

    }
    q = atomic[vec_spack[i].atom[0]] * atomic[vec_spack[i].atom[1]] / dist_2;

    esum += en;
    qsum += q;
  }


  en = ((qsum + ((2.0*esum) / SCON_PI))*coulomb);
  // 	 cout<<end-start<<"SP-energy"<<endl;
  energy_short = en;



  //	
  return energy_short;
}
void energy::interfaces::amoeba::amoeba_ff::Spackman_GRAD() {








  size_t j, k, it, n;


  double coulomb = 332.063709;

  double tnm, sum, f, argBessel0, e, q, qsum = 0, esum = 0;

  double dist_2, xgrad, ygrad, zgrad, distx, disty, distz;


  //		
  auto const &positions = coords->xyz();
  coords::Cartesian_Point bv, b, gv;




  smax = 10.0;
  aLimes = 0.0;
  bLimes = smax * 4 * SCON_PI;




  double r_cutoff;
  r_cutoff = Config::get().energy.spackman.cut;

  xgrad = 0.0;
  ygrad = 0.0;
  zgrad = 0.0;
  for (size_t i = 0; i < vec_spack.size(); i++) {

    bv = positions[vec_spack[i].atom[0]] - positions[vec_spack[i].atom[1]];
    distx = positions[vec_spack[i].atom[0]].x() - positions[vec_spack[i].atom[1]].x();
    disty = positions[vec_spack[i].atom[0]].y() - positions[vec_spack[i].atom[1]].y();
    distz = positions[vec_spack[i].atom[0]].z() - positions[vec_spack[i].atom[1]].z();
    dist_2 = len(bv);



    //	
    //	
    if (dist_2 > r_cutoff) continue;
    //  if(dist_2 <= 3.0) continue;
    k = 1;
    f = fscat[vec_spack[i].atom[0]][k - 1] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[0]] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[1]] * fscat[vec_spack[i].atom[0]][k - 1];

    argBessel0 = sscat[k] * dist_2;

    f = f*((argBessel0*std::cos(argBessel0) - std::sin(argBessel0)) / (sscat[k] * (dist_2*dist_2)));
    e = (bLimes - aLimes)*f;

    for (j = 2; j <= nref; j++) {
      it = (size_t)pow(3.0, double(j - 2));

      tnm = double(it);

      sum = 0.0;

      for (n = 1; n <= it; n++) {
        k++;
        f = fscat[vec_spack[i].atom[0]][k - 1] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[0]] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[1]] * fscat[vec_spack[i].atom[0]][k - 1];
        argBessel0 = sscat[k] * dist_2;
        f = f*((argBessel0*std::cos(argBessel0) - std::sin(argBessel0)) / (sscat[k] * (dist_2*dist_2)));
        sum += f;
        k++;
        f = fscat[vec_spack[i].atom[0]][k - 1] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[0]] * fscat[vec_spack[i].atom[1]][k - 1] - atomic[vec_spack[i].atom[1]] * fscat[vec_spack[i].atom[0]][k - 1];
        argBessel0 = sscat[k] * dist_2;
        f = f*((argBessel0*std::cos(argBessel0) - std::sin(argBessel0)) / (sscat[k] * (dist_2*dist_2)));

        sum += f;

      }
      e = (e + (bLimes - aLimes)*sum / tnm) / 3.0;

    }
    q = -double(atomic[vec_spack[i].atom[0]] * atomic[vec_spack[i].atom[1]]) / (dist_2*dist_2);

    esum = e;

    qsum = q;
    e = ((qsum + ((2.0*esum) / SCON_PI))*coulomb);


    xgrad = e * (distx / dist_2);
    ygrad = e * (disty / dist_2);
    zgrad = e * (distz / dist_2);

    gv.x() = xgrad;
    gv.y() = ygrad;
    gv.z() = zgrad;

    part_grad[SHORTRANGE][vec_spack[i].atom[0]] += gv;
    part_grad[SHORTRANGE][vec_spack[i].atom[1]] -= gv;

  }


}
//!calculation of analytical gradients for Spackman correction using generated list from Spackman_list	  
void energy::interfaces::amoeba::amoeba_ff::SpackmanGrad_3()
{
  size_t n(0), contr(0);
  double distx(0.0), disty(0.0), distz(0.0), dist_3(0.0),
	  xx_in(0.0), fac_x(0.0), fac_y(0.0), fac_z(0.0),
	  xgrad(0.0), ygrad(0.0), zgrad(0.0), y(0.0);
  std::vector <double> dist_2;
  auto const &positions = coords->xyz();
  coords::Cartesian_Point bv, b, gv;
  //scon::nv3d &gradients = part_grad[energy::interfaces::amoeba::types::short_range_an];
  //	    double start = omp_get_wtime();
  //
  //	    
  //
  //
  //
  xgrad = 0.0;
  ygrad = 0.0;
  zgrad = 0.0;
  //!initialize Spline routine
  //!Spline_interp == cubic spline
  //!Poly_interp == polynominal spline
  //!Linear_interp == linear spline	    
  ///CC-interaction
  Linear_interp_sorted myfunc11(eveca1, exa22);
  ///HH-interaction
  Linear_interp_sorted myfunc22(eveca1, exa11);
  ///CH-interaction
  Linear_interp_sorted myfunc33(eveca1, exa33);
  ///OH-interaction
  Linear_interp_sorted myfunc44(eveca1, exa44);
  ///OO-interaction
  Linear_interp_sorted myfunc55(eveca1, exa55);
  ///OC-interaction
  Linear_interp_sorted myfunc66(eveca1, exa66);
  ///NH-interaction
  Linear_interp_sorted myfunc77(eveca1, exa77);
  ///NN-interaction
  Linear_interp_sorted myfunc88(eveca1, exa88);
  ///NC-interaction
  Linear_interp_sorted myfunc99(eveca1, exa99);
  ///NO-interaction
  Linear_interp_sorted myfunc1010(eveca1, exa1010);
  //
  //
  ////!loop over monomer interactions 	    
  //	    
  for (size_t i = 0; i < vec_spack.size(); i++) {
    //	    
    //	    
    //	    
    //
    ////!definition of interaction type depending on atomtype parameters

    if (coords->atoms(vec_spack[i].atom[0]).symbol() == "C" && coords->atoms(vec_spack[i].atom[1]).symbol() == "C") contr = 1;
    if (coords->atoms(vec_spack[i].atom[0]).symbol() == "H" && coords->atoms(vec_spack[i].atom[1]).symbol() == "H") contr = 2;
    if (coords->atoms(vec_spack[i].atom[0]).symbol() == "C" && coords->atoms(vec_spack[i].atom[1]).symbol() == "H") contr = 3;
    if (coords->atoms(vec_spack[i].atom[0]).symbol() == "H" && coords->atoms(vec_spack[i].atom[1]).symbol() == "C") contr = 3;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "O" && coords->atoms(vec_spack[i].atom[1]).symbol() == "H") contr = 4;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "H" && coords->atoms(vec_spack[i].atom[1]).symbol() == "O") contr = 4;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "O" && coords->atoms(vec_spack[i].atom[1]).symbol() == "O") contr = 5;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "O" && coords->atoms(vec_spack[i].atom[1]).symbol() == "C") contr = 6;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "C" && coords->atoms(vec_spack[i].atom[1]).symbol() == "O") contr = 6;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "N" && coords->atoms(vec_spack[i].atom[1]).symbol() == "H") contr = 7;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "H" && coords->atoms(vec_spack[i].atom[1]).symbol() == "N") contr = 7;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "N" && coords->atoms(vec_spack[i].atom[1]).symbol() == "N") contr = 8;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "N" && coords->atoms(vec_spack[i].atom[1]).symbol() == "C") contr = 9;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "C" && coords->atoms(vec_spack[i].atom[1]).symbol() == "N") contr = 9;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "N" && coords->atoms(vec_spack[i].atom[1]).symbol() == "O") contr = 10;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "O" && coords->atoms(vec_spack[i].atom[1]).symbol() == "N") contr = 10;


    distx = positions[vec_spack[i].atom[0]].x() - positions[vec_spack[i].atom[1]].x();
    disty = positions[vec_spack[i].atom[0]].y() - positions[vec_spack[i].atom[1]].y();
    distz = positions[vec_spack[i].atom[0]].z() - positions[vec_spack[i].atom[1]].z();
    bv = positions[vec_spack[i].atom[0]] - positions[vec_spack[i].atom[1]];

    dist_3 = len(bv);

    //


    if (dist_3 <= 10.0) {
      //!for next entries in analytical list rounded up	    
      /*dist_4=dist_3*10000.00;

      ptrdiff_t l = dist_4/10;

      ptrdiff_t t= dist_4%10;

      if ( t >= 5) l++ ;*/

      xx_in = dist_3;
      fac_x = distx / xx_in;
      fac_y = disty / xx_in;
      fac_z = distz / xx_in;
      n = 1000;
      if (contr == 1) {


        y = myfunc11.interpolate(xx_in);

        xgrad = (fac_x*y);
        ygrad = (fac_y*y);
        zgrad = (fac_z*y);

      }

      else if (contr == 2) {

        y = myfunc22.interpolate(xx_in);
        xgrad = (fac_x*y);
        ygrad = (fac_y*y);
        zgrad = (fac_z*y);

      }

      else if (contr == 3) {

        y = myfunc33.interpolate(xx_in);
        xgrad = (fac_x*y);
        ygrad = (fac_y*y);
        zgrad = (fac_z*y);

      }
	  else if (contr == 4) {

		  y = myfunc44.interpolate(xx_in);

		  xgrad = (fac_x*y);
		  ygrad = (fac_y*y);
		  zgrad = (fac_z*y);

	  }
	  else if (contr == 5) {

		  y = myfunc55.interpolate(xx_in);

		  xgrad = (fac_x*y);
		  ygrad = (fac_y*y);
		  zgrad = (fac_z*y);

	  }
	  else if (contr == 6) {

		  y = myfunc66.interpolate(xx_in);

		  xgrad = (fac_x*y);
		  ygrad = (fac_y*y);
		  zgrad = (fac_z*y);

	  }
	  else if (contr == 7) {

		  y = myfunc77.interpolate(xx_in);

		  xgrad = (fac_x*y);
		  ygrad = (fac_y*y);
		  zgrad = (fac_z*y);

	  }
	  else if (contr == 8) {

		  y = myfunc88.interpolate(xx_in);

		  xgrad = (fac_x*y);
		  ygrad = (fac_y*y);
		  zgrad = (fac_z*y);

	  }
	  else if (contr == 9) {

		  y = myfunc99.interpolate(xx_in);

		  xgrad = (fac_x*y);
		  ygrad = (fac_y*y);
		  zgrad = (fac_z*y);

	  }
	  else if (contr == 10) {

		  y = myfunc1010.interpolate(xx_in);

		  xgrad = (fac_x*y);
		  ygrad = (fac_y*y);
		  zgrad = (fac_z*y);

	  }

    }
    else  continue;


    gv.x() = xgrad;
    gv.y() = ygrad;
    gv.z() = zgrad;

    part_grad[SHORTRANGE][vec_spack[i].atom[0]] += gv;
    part_grad[SHORTRANGE][vec_spack[i].atom[1]] -= gv;

  }

}
//
double energy::interfaces::amoeba::amoeba_ff::Spackman_energy_analytical()
{
  size_t n(0), contr(0);
  double distx = 0, disty = 0, distz = 0, dist_3 = 0,xx_in(0.0),
	  y(0.0), energytemp(0.0);
  std::vector <double> dist_2;
  auto const &positions = coords->xyz();
  coords::Cartesian_Point bv;

  //!initialize Spline routine
  //!Spline_interp == cubic spline
  //!Poly_interp == polynominal spline
  //!Linear_interp == linear spline
  ///CC-interaction
  Linear_interp_sorted myfunc11(eveca1, dex22);
  ///HH-interaction
  Linear_interp_sorted myfunc22(eveca1, dex11);
  ///CH-interaction
  Linear_interp_sorted myfunc33(eveca1, dex33);
  ///OH-interaction
  Linear_interp_sorted myfunc44(eveca1, dex44);
  ///OO-interaction
  Linear_interp_sorted myfunc55(eveca1, dex55);
  ///OC-interaction
  Linear_interp_sorted myfunc66(eveca1, dex66);
  ///NH-interaction
  Linear_interp_sorted myfunc77(eveca1, dex77);
  ///NN-interaction
  Linear_interp_sorted myfunc88(eveca1, dex88);
  ///NC-interaction
  Linear_interp_sorted myfunc99(eveca1, dex99);
  ///NO-interaction
  Linear_interp_sorted myfunc1010(eveca1, dex1010);
  //
  //
  ////!loop over monomer interactions 	    
  //	
  for (size_t i = 0; i < vec_spack.size(); i++) {
    //	    
    //	    
    //	    
    //
    ////!definition of interaction type depending on atomtype parameters	 
    if (coords->atoms(vec_spack[i].atom[0]).symbol() == "C" && coords->atoms(vec_spack[i].atom[1]).symbol() == "C") contr = 1;
    if (coords->atoms(vec_spack[i].atom[0]).symbol() == "H" && coords->atoms(vec_spack[i].atom[1]).symbol() == "H") contr = 2;
    if (coords->atoms(vec_spack[i].atom[0]).symbol() == "C" && coords->atoms(vec_spack[i].atom[1]).symbol() == "H") contr = 3;
    if (coords->atoms(vec_spack[i].atom[0]).symbol() == "H" && coords->atoms(vec_spack[i].atom[1]).symbol() == "C") contr = 3;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "O" && coords->atoms(vec_spack[i].atom[1]).symbol() == "H") contr = 4;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "H" && coords->atoms(vec_spack[i].atom[1]).symbol() == "O") contr = 4;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "O" && coords->atoms(vec_spack[i].atom[1]).symbol() == "O") contr = 5;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "O" && coords->atoms(vec_spack[i].atom[1]).symbol() == "C") contr = 6;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "C" && coords->atoms(vec_spack[i].atom[1]).symbol() == "O") contr = 6;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "N" && coords->atoms(vec_spack[i].atom[1]).symbol() == "H") contr = 7;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "H" && coords->atoms(vec_spack[i].atom[1]).symbol() == "N") contr = 7;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "N" && coords->atoms(vec_spack[i].atom[1]).symbol() == "N") contr = 8;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "N" && coords->atoms(vec_spack[i].atom[1]).symbol() == "C") contr = 9;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "C" && coords->atoms(vec_spack[i].atom[1]).symbol() == "N") contr = 9;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "N" && coords->atoms(vec_spack[i].atom[1]).symbol() == "O") contr = 10;
	if (coords->atoms(vec_spack[i].atom[0]).symbol() == "O" && coords->atoms(vec_spack[i].atom[1]).symbol() == "N") contr = 10;

    distx = positions[vec_spack[i].atom[0]].x() - positions[vec_spack[i].atom[1]].x();
    disty = positions[vec_spack[i].atom[0]].y() - positions[vec_spack[i].atom[1]].y();
    distz = positions[vec_spack[i].atom[0]].z() - positions[vec_spack[i].atom[1]].z();
    bv = positions[vec_spack[i].atom[0]] - positions[vec_spack[i].atom[1]];

    dist_3 = scon::len(bv);
    //
    if (dist_3 <= 10.0) {

      xx_in = dist_3;
      n = 1000;
      if (contr == 1) {
        y = myfunc11.interpolate(xx_in);

      }

      else if (contr == 2) {
        y = myfunc22.interpolate(xx_in);


      }

      else if (contr == 3) {

        y = myfunc33.interpolate(xx_in);

      }
	  else if (contr == 4) {

		  y = myfunc44.interpolate(xx_in);


	  }
	  else if (contr == 5) {

		  y = myfunc55.interpolate(xx_in);


	  }
	  else if (contr == 6) {

		  y = myfunc66.interpolate(xx_in);


	  }
	  else if (contr == 7) {

		  y = myfunc77.interpolate(xx_in);


	  }
	  else if (contr == 8) {

		  y = myfunc88.interpolate(xx_in);


	  }
	  else if (contr == 9) {

		  y = myfunc99.interpolate(xx_in);


	  }
	  else if (contr == 10) {

		  y = myfunc1010.interpolate(xx_in);


	  }
    }

    else {
      y = 0.0;

    }

    energytemp += y;




  }
  return energytemp;
}
