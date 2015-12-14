//#include <cmath>
//#include <stddef.h>
//#include <stdexcept>
//#include <cstdlib>
//#include "energy_int_aco.h"
//#include "configuration.h"
//#include "scon_utility.h"
//
//#if defined _OPENMP
//#include <omp.h>
//#endif
//
//#define SUPERPI 3.141592653589793238
//#define SQRTPI  sqrt(3.141592653589793238)
//
//
//// Direct energy for 1-2 and 1-3 interactions
//void energy::interfaces::aco::aco_ff::pme_direct_scaled(double & e_nb,
//  coords::Representation_3D &grad_vector)
//{
//  double dist, r, achg, distbuf, rbuf, factor, ewald, error;
//  double energy, den, cutoff, cutoff2, scale;
//  cutoff = Config::get().energy.cutoff;
//  cutoff2 = cutoff*cutoff;
//  achg = 0.0;
//  coords::Cartesian_Point vir;
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
//    den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r);
//    den /= r;
//    b *= den;
//    e_nb += energy;
//    grad_vector[bond.atoms[0]] += b;
//    grad_vector[bond.atoms[1]] -= b;;
//    coords::float_type const vxx = b.x() * vir.x();
//    coords::float_type const vyx = b.x() * vir.y();
//    coords::float_type const vzx = b.x() * vir.z();
//    coords::float_type const vyy = b.y() * vir.y();
//    coords::float_type const vzy = b.y() * vir.z();
//    coords::float_type const vzz = b.z() * vir.z();
//    part_virial[VDWC][0][0] += vxx;
//    part_virial[VDWC][1][0] += vyx;
//    part_virial[VDWC][2][0] += vzx;
//    part_virial[VDWC][0][1] += vyx;
//    part_virial[VDWC][1][1] += vyy;
//    part_virial[VDWC][2][1] += vzy;
//    part_virial[VDWC][0][2] += vzx;
//    part_virial[VDWC][1][2] += vzy;
//    part_virial[VDWC][2][2] += vzz;
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
//    den = -factor * ((error + scale) / distbuf + (2.0*coords->pme.pmetemp.ewaldcoeff / SQRTPI) * exp(-(ewald*ewald)) / r);
//    den /= r;
//    b *= den;
//    e_nb += energy;
//    grad_vector[angle.atoms[0]] += b;
//    grad_vector[angle.atoms[2]] -= b;
//    coords::float_type const vxx = b.x() * vir.x();
//    coords::float_type const vyx = b.x() * vir.y();
//    coords::float_type const vzx = b.x() * vir.z();
//    coords::float_type const vyy = b.y() * vir.y();
//    coords::float_type const vzy = b.y() * vir.z();
//    coords::float_type const vzz = b.z() * vir.z();
//    part_virial[VDWC][0][0] += vxx;
//    part_virial[VDWC][1][0] += vyx;
//    part_virial[VDWC][2][0] += vzx;
//    part_virial[VDWC][0][1] += vyx;
//    part_virial[VDWC][1][1] += vyy;
//    part_virial[VDWC][2][1] += vzy;
//    part_virial[VDWC][0][2] += vzx;
//    part_virial[VDWC][1][2] += vzy;
//    part_virial[VDWC][2][2] += vzz;
//  }
//
//}
//
////// Preliminary correction term for charged simulation cells (NOT WORKING!) See GROMACS source code!
////void energy::interfaces::aco::aco_ff::pme_chargedcell()
////{
////  double q(0.0);
////  std::cout << "Just a simple test" << std::endl;
////  for (int i = 0; i < coords->pme.pmetemp.natoms; i++)
////  {
////    q += coords->pme.pmetemp.atomcharges[i];
////  }
////  std::cout << "Total absolute charge of the system:  " << q << std::endl;
////  q *= 332.0716;
////  q = abs(q*q);
////
////  q = q / (8 * Config::get().energy.pb_box.x() * Config::get().energy.pb_box.y() * Config::get().energy.pb_box.z() * coords->pme.pmetemp.ewaldcoeff * coords->pme.pmetemp.ewaldcoeff);
////  std::cout << "Corrected charge value:  " << q << std::endl;
////}
//
//// Reciprocal space energy 
//void energy::interfaces::aco::aco_ff::pme_reci(double & energy)
//{
//  double f, pifact, volumestuff, height1, height2, height3, square;
//  double tempterm, expo, nenner, chargeval, e_rec(0.0), e_temp(0.0), virial;
//  int gridpoints, planepoints, tempx, tempy, tempz, i, j, k;
//  int k1, k2, k3, m1, m2, m3, b1, b2, b3;
//  f = 0.5 * 332.0716 / 1;
//  int stuffx = coords->pme.pmetemp.nxpoints;
//  int stuffy = coords->pme.pmetemp.nypoints;
//  int stuffz = coords->pme.pmetemp.nzpoints;
//  gridpoints = coords->pme.pmetemp.nxpoints * coords->pme.pmetemp.nypoints * coords->pme.pmetemp.nzpoints;
//  pifact = (SUPERPI / coords->pme.pmetemp.ewaldcoeff)*(SUPERPI / coords->pme.pmetemp.ewaldcoeff);
//  volumestuff = SUPERPI * Config::get().energy.pb_box.x() * Config::get().energy.pb_box.y() * Config::get().energy.pb_box.z();
//  planepoints = coords->pme.pmetemp.nxpoints * coords->pme.pmetemp.nypoints;
//  tempx = (coords->pme.pmetemp.nxpoints + 1) / 2;
//  tempy = (coords->pme.pmetemp.nypoints + 1) / 2;
//  tempz = (coords->pme.pmetemp.nzpoints + 1) / 2;
//  //loop over all gridpoints
//  for (i = 0; i < gridpoints - 1; i++)
//  {
//    k3 = (i + 1) / planepoints + 1;
//    j = (i + 1) - (k3 - 1)*planepoints;
//    k2 = j / coords->pme.pmetemp.nxpoints + 1;
//    k1 = j - (k2 - 1)*coords->pme.pmetemp.nxpoints + 1;
//    m1 = k1 - 1;
//    m2 = k2 - 1;
//    m3 = k3 - 1;
//    if (k1 > tempx) m1 -= coords->pme.pmetemp.nxpoints;
//    if (k2 > tempy) m2 -= coords->pme.pmetemp.nypoints;
//    if (k3 > tempz) m3 -= coords->pme.pmetemp.nzpoints;
//    height1 = coords->pme.pmetemp.recivectors(0, 0) * m1 + coords->pme.pmetemp.recivectors(0, 1) * m2 + coords->pme.pmetemp.recivectors(0, 2) * m3;
//    height2 = coords->pme.pmetemp.recivectors(1, 0) * m1 + coords->pme.pmetemp.recivectors(1, 1) * m2 + coords->pme.pmetemp.recivectors(1, 2) * m3;
//    height3 = coords->pme.pmetemp.recivectors(2, 0) * m1 + coords->pme.pmetemp.recivectors(2, 1) * m2 + coords->pme.pmetemp.recivectors(2, 2) * m3;
//    square = height1*height1 + height2*height2 + height3*height3;
//    tempterm = -pifact * square;
//    expo = 0.0;
//    if (tempterm > -50)
//    {
//      nenner = volumestuff * square * coords->pme.pmetemp.moduli1[k1 - 1] * coords->pme.pmetemp.moduli2[k2 - 1] * coords->pme.pmetemp.moduli1[k3 - 1];
//      expo = exp(tempterm) / nenner;
//      // check for octahedral box
//      if (Config::get().energy.periodic == false)
//      {
//          expo = expo * (1.0 - cos(SUPERPI*Config::get().energy.pb_box.x()*sqrt(square)));
//      }
//      else if (Config::get().energy.pb_box.x() != Config::get().energy.pb_box.y() || Config::get().energy.pb_box.x() != Config::get().energy.pb_box.z() || Config::get().energy.pb_box.y() != Config::get().energy.pb_box.z())
//      {
//         if ((m1 + m2 + m3) % 2 != 0) expo = 0.0;
//      }
//      chargeval = coords->pme.pmetemp.in[((k1 - 1)*stuffy + (k2 - 1))*stuffz + (k3 - 1)][0] * coords->pme.pmetemp.in[((k1 - 1)*stuffy + (k2 - 1))*stuffz + (k3 - 1)][0] + coords->pme.pmetemp.in[((k1 - 1)*stuffy + (k2 - 1))*stuffz + (k3 - 1)][1] * coords->pme.pmetemp.in[((k1 - 1)*stuffy + (k2 - 1))*stuffz + (k3 - 1)][1];
//      e_temp = f * expo * chargeval;
//      energy += e_temp;
//      // get terms for virial tensor
//      virial = (2.0 / square) * (1.0 - tempterm) * e_temp;
//      part_virial[VDWC][0][0] += height1*height1*virial - e_temp;
//      part_virial[VDWC][1][0] += height1*height2*virial;
//      part_virial[VDWC][2][0] += height1*height3*virial;
//      part_virial[VDWC][0][1] += height2*height1*virial;
//      part_virial[VDWC][1][1] += height2*height2*virial - e_temp;
//      part_virial[VDWC][2][1] += height2*height3*virial;
//      part_virial[VDWC][0][2] += height3*height1*virial;
//      part_virial[VDWC][1][2] += height3*height2*virial;
//      part_virial[VDWC][2][2] += height3*height3*virial - e_temp;
//    }// end of if clause
//    coords->pme.pmetemp.in[((k1 - 1)*stuffy + (k2 - 1))*stuffz + (k3 - 1)][0] = expo * coords->pme.pmetemp.in[((k1 - 1)*stuffy + (k2 - 1))*stuffz + (k3 - 1)][0];
//    coords->pme.pmetemp.in[((k1 - 1)*stuffy + (k2 - 1))*stuffz + (k3 - 1)][1] = expo * coords->pme.pmetemp.in[((k1 - 1)*stuffy + (k2 - 1))*stuffz + (k3 - 1)][1];
//  }// end of gridpoint loop 
//}
//
//
//// reciprocal space gradients
//void energy::interfaces::aco::aco_ff::pme_reci_grad(coords::Representation_3D &grad)
//{
//  double factor = 332.0716;
//  double nx, ny, nz, xgrid, ygrid, zgrid, chargefactor;
//  double dex, dey, dez, xg, yg, zg, thet1, thet2, thet3;
//  double dtet3, dtet2, dtet1, tcharge;
//  int atom, k, j, i;
//  coords::Cartesian_Point tempgrad, difgrad;
//  nx = coords->pme.pmetemp.nxpoints;
//  ny = coords->pme.pmetemp.nypoints;
//  nz = coords->pme.pmetemp.nzpoints;
//  for (int u = 0; u < coords->pme.pmetemp.natoms; u++)
//  {
//    atom = u;
//    xgrid = coords->pme.pmetemp.initgrid(0, atom);
//    ygrid = coords->pme.pmetemp.initgrid(1, atom);
//    zgrid = coords->pme.pmetemp.initgrid(2, atom);
//    chargefactor = factor * coords->pme.pmetemp.atomcharges[atom];
//    tempgrad.x() = tempgrad.y() = tempgrad.z() = 0.0;
//    dex = dey = dez = 0.0;
//    zg = zgrid;
//    for (int i1 = 1; i1 <= coords->pme.pmetemp.bsplineorder; i1++)
//    {
//      zg += 1;
//      k = zg + 1 + (coords->pme.pmetemp.nzpoints - copysign(coords->pme.pmetemp.nzpoints, zg)) / 2;
//      thet3 = coords->pme.pmetemp.bscz(0, i1 - 1, atom);
//      dtet3 = nz * coords->pme.pmetemp.bscz(1, i1 - 1, atom);
//      yg = ygrid;
//      for (int k1 = 1; k1 <= coords->pme.pmetemp.bsplineorder; k1++)
//      {
//        yg += 1;
//        j = yg + 1 + (coords->pme.pmetemp.nypoints - copysign(coords->pme.pmetemp.nypoints, yg)) / 2;
//        thet2 = coords->pme.pmetemp.bscy(0, k1 - 1, atom);
//        dtet2 = ny * coords->pme.pmetemp.bscy(1, k1 - 1, atom);
//        xg = xgrid;
//        for (int j1 = 1; j1 <= coords->pme.pmetemp.bsplineorder; j1++)
//        {
//          xg += 1;
//          i = xg + 1 + (coords->pme.pmetemp.nxpoints - copysign(coords->pme.pmetemp.nxpoints, xg)) / 2;
//          thet1 = coords->pme.pmetemp.bscx(0, j1 - 1, atom);
//          dtet1 = nx * coords->pme.pmetemp.bscx(1, j1 - 1, atom);
//          tcharge = coords->pme.pmetemp.in[((i - 1)*coords->pme.pmetemp.nxpoints + (j - 1))*coords->pme.pmetemp.nxpoints + (k - 1)][0];
//          tempgrad.x() += tcharge*dtet1*thet2*thet3;
//          tempgrad.y() += tcharge*dtet2*thet1*thet3;
//          tempgrad.z() += tcharge*dtet3*thet1*thet2;
//        }
//      }
//    }// end of for
//    difgrad.x() = (chargefactor*(coords->pme.pmetemp.recivectors(0, 0) * tempgrad.x() + coords->pme.pmetemp.recivectors(0, 1) * tempgrad.y() + coords->pme.pmetemp.recivectors(0, 2) * tempgrad.z()));
//    difgrad.y() = (chargefactor*(coords->pme.pmetemp.recivectors(1, 0) * tempgrad.x() + coords->pme.pmetemp.recivectors(1, 1) * tempgrad.y() + coords->pme.pmetemp.recivectors(1, 2) * tempgrad.z()));
//    difgrad.z() = (chargefactor*(coords->pme.pmetemp.recivectors(2, 0) * tempgrad.x() + coords->pme.pmetemp.recivectors(2, 1) * tempgrad.y() + coords->pme.pmetemp.recivectors(2, 2) * tempgrad.z()));
//    grad[u] += difgrad;
//  }// end of atoms for loop
//
//}
//
//
//// Self interaction correction term
//void energy::interfaces::aco::aco_ff::pme_correction(double & e_corr)
//{
//  double constant = 332.0716;
//  double selffactor;
//  e_corr = 0.0;
//  selffactor = -constant * coords->pme.pmetemp.ewaldcoeff / SQRTPI;
//  // get charges for each atom. -> Precalculating this stuff should be possible
//  for (int i = 0; i < coords->pme.pmetemp.natoms; i++)
//  {
//    std::size_t const type = cparams.indextype(tinker::CHARGE) == tinker::GROUP ? cparams.group_by_type(refined.type(i)) : refined.type(i);
//    double const charge_i = cparams.charges()[type - 1].c;
//    coords->pme.pmetemp.atomcharges[i] = charge_i;
//    e_corr += selffactor * charge_i * charge_i;
//  }
//
//}
//
//// Generate bsplines for the atomic sites
//void energy::interfaces::aco::aco_ff::bsplinegeneration()
//{
//  double tresh = coords->pme.pmetemp.treshold;
//  double x, y, z, sr, orr;
//  int o;
//  coords::Cartesian_Point tempcor;
//  // get spline coeffcicients for all atoms via loop over single atoms
//  for (int i = 0; i < coords->pme.pmetemp.natoms; i++)
//  {
//    tempcor = coords->xyz(i);
//    //generate spline for atomic site i in x direction
//    sr = tempcor.x() * coords->pme.pmetemp.recivectors(0, 0) + tempcor.y() * coords->pme.pmetemp.recivectors(1, 0) + tempcor.z() * coords->pme.pmetemp.recivectors(2, 0);
//    orr = coords->pme.pmetemp.nxpoints * (sr - lround(sr) + 0.5);
//    o = int(orr - tresh);
//    sr = orr - double(o);
//    coords->pme.pmetemp.initgrid(0, i) = o - coords->pme.pmetemp.bsplineorder;
//    singlespline(sr, i, 1);
//    //generate spline for atomic site i in y direction
//    sr = tempcor.x() * coords->pme.pmetemp.recivectors(0, 1) + tempcor.y() * coords->pme.pmetemp.recivectors(1, 1) + tempcor.z() * coords->pme.pmetemp.recivectors(2, 1);
//    orr = coords->pme.pmetemp.nypoints * (sr - lround(sr) + 0.5);
//    o = int(orr - tresh);
//    sr = orr - double(o);
//    coords->pme.pmetemp.initgrid(1, i) = o - coords->pme.pmetemp.bsplineorder;
//    singlespline(sr, i, 2);
//    //generate spline for atomic site i in z direction
//    sr = tempcor.x() * coords->pme.pmetemp.recivectors(0, 2) + tempcor.y() * coords->pme.pmetemp.recivectors(1, 2) + tempcor.z() * coords->pme.pmetemp.recivectors(2, 2);
//    orr = coords->pme.pmetemp.nzpoints * (sr - lround(sr) + 0.5);
//    o = int(orr - tresh);
//    sr = orr - double(o);
//    coords->pme.pmetemp.initgrid(2, i) = o - coords->pme.pmetemp.bsplineorder;
//    singlespline(sr, i, 3);
//
//  }
//}
//
//// set up splines at atomic sites in each dimension up to the bspline order
//void energy::interfaces::aco::aco_ff::singlespline(double term, int atomnum, int splinedirection)
//{
//  int level = 2;
//  int k;
//  double nenner;
//  double temp[11][11];
//  // initialize first and second order terms
//  temp[2][2] = term;
//  temp[2][1] = 1.0 - term;
//  // get 3rd order terms
//  temp[3][3] = 0.5 * term * temp[2][2];
//  temp[3][2] = 0.5 * ((1.0 + term)*temp[2][1] + (2.0 - term)*temp[2][2]);
//  temp[3][1] = 0.5 * (1.0 - term) * temp[2][1];
//  // get all the higher orders
//  for (int i = 4; i <= coords->pme.pmetemp.bsplineorder; i++)
//  {
//    k = i - 1;
//    nenner = 1.0 / double(k);
//    temp[i][i] = nenner * term *temp[k][k];
//    for (int j = 1; j <= (i - 2); j++)
//    {
//      temp[i][i - j] = nenner *((term + double(j))*temp[k][i - j - 1] + (double(i - j) - term) * temp[k][i - j]);
//    }
//    temp[i][1] = nenner * (1.0 - term) * temp[k][1];
//  }
//  // get first derivative coefficients
//  k = coords->pme.pmetemp.bsplineorder - 1;
//  temp[k][coords->pme.pmetemp.bsplineorder] = temp[k][coords->pme.pmetemp.bsplineorder - 1];
//  for (int i = coords->pme.pmetemp.bsplineorder - 1; i >= 2; i--)
//  {
//    temp[k][i] = temp[k][i - 1] - temp[k][i];
//  }
//  temp[k][1] = -temp[k][1];
//  // temporary -> final array copying
//  for (int i = 1; i <= coords->pme.pmetemp.bsplineorder; i++)
//  {
//    for (int k = 1; k <= level; k++)
//    {
//      //std::cout << "Final values" << std::endl;
//      if (splinedirection == 1)      coords->pme.pmetemp.bscx(k - 1, i - 1, atomnum) = temp[coords->pme.pmetemp.bsplineorder - k + 1][i];
//      else if (splinedirection == 2) coords->pme.pmetemp.bscy(k - 1, i - 1, atomnum) = temp[coords->pme.pmetemp.bsplineorder - k + 1][i];
//      else if (splinedirection == 3) coords->pme.pmetemp.bscz(k - 1, i - 1, atomnum) = temp[coords->pme.pmetemp.bsplineorder - k + 1][i];
//    }
//  }
//}
//
//// get the spatial chunk for each atomic site and store it in the pmetable
//void energy::interfaces::aco::aco_ff::generatetable()
//{
//  int bounds1[6], bounds2[6], points[3];
//  std::vector<int> pointid;
//  pointid.resize(3);
//  bool xn, yn, zn, xp, yp, zp, xm, ym, zm;
//  int p;
//  // zero out array
//  for (int k = 0; k < coords->pme.pmetemp.rgridtotal; k++)
//  {
//    for (int i = 0; i < coords->pme.pmetemp.natoms; i++)
//    {
//      coords->pme.pmetemp.parallelpme(i, k) = 0;
//    }
//  }
//  // loop over all atomic sites, put site onto respective chunk on pme grid
//  for (int i = 0; i < coords->pme.pmetemp.natoms; i++)
//  {
//    points[0] = coords->pme.pmetemp.initgrid(0, i) + coords->pme.pmetemp.bsoffset;
//    points[1] = coords->pme.pmetemp.initgrid(1, i) + coords->pme.pmetemp.bsoffset;
//    points[2] = coords->pme.pmetemp.initgrid(2, i) + coords->pme.pmetemp.bsoffset;
//    // x direction
//    if (points[0] < 1)
//    {
//      points[0] = points[0] % coords->pme.pmetemp.nxpoints + coords->pme.pmetemp.nxpoints;
//    }
//    else if (points[0] > coords->pme.pmetemp.nxpoints)
//    {
//      points[0] = points[0] % coords->pme.pmetemp.nxpoints;
//    }
//    // y direction
//    if (points[1] < 1)
//    {
//      points[1] = points[1] % coords->pme.pmetemp.nypoints + coords->pme.pmetemp.nypoints;
//    }
//    else if (points[1] > coords->pme.pmetemp.nypoints)
//    {
//      points[1] = points[1] % coords->pme.pmetemp.nypoints;
//    }
//    // z direction
//    if (points[2] < 1)
//    {
//      points[2] = points[2] % coords->pme.pmetemp.nzpoints + coords->pme.pmetemp.nzpoints;
//    }
//    else if (points[2] > coords->pme.pmetemp.nzpoints)
//    {
//      points[2] = points[2] % coords->pme.pmetemp.nzpoints;
//    }
//    bounds1[0] = points[0] - coords->pme.pmetemp.roughleft;
//    bounds1[1] = points[0] + coords->pme.pmetemp.roughright;
//    bounds1[2] = points[1] - coords->pme.pmetemp.roughleft;
//    bounds1[3] = points[1] + coords->pme.pmetemp.roughright;
//    bounds1[4] = points[2] - coords->pme.pmetemp.roughleft;
//    bounds1[5] = points[2] + coords->pme.pmetemp.roughright;
//    pointid[0] = (points[0] - 1) / coords->pme.pmetemp.rgrid1 + 1;
//    pointid[1] = (points[1] - 1) / coords->pme.pmetemp.rgrid2 + 1;
//    pointid[2] = (points[2] - 1) / coords->pme.pmetemp.rgrid3 + 1;
//    bounds2[0] = (pointid[0] - 1)* coords->pme.pmetemp.rgrid1 + 1;
//    bounds2[1] = bounds2[0] + coords->pme.pmetemp.rgrid1 - 1;
//    bounds2[2] = (pointid[1] - 1)* coords->pme.pmetemp.rgrid2 + 1;
//    bounds2[3] = bounds2[2] + coords->pme.pmetemp.rgrid2 - 1;
//    bounds2[4] = (pointid[2] - 1)* coords->pme.pmetemp.rgrid3 + 1;
//    bounds2[5] = bounds2[4] + coords->pme.pmetemp.rgrid3 - 1;
//    // get the chunk for the current atomic site
//    p = (pointid[2] - 1)*coords->pme.pmetemp.nrough1*coords->pme.pmetemp.nrough2 + (pointid[1] - 1)*coords->pme.pmetemp.nrough1 + pointid[0];
//    coords->pme.pmetemp.parallelpme(i, p - 1) = 1;
//    // set bool flags for each chunk
//    xn = (bounds1[0] < bounds2[0]) ? 1 : 0;
//    yn = (bounds1[2] < bounds2[2]) ? 1 : 0;
//    zn = (bounds1[4] < bounds2[4]) ? 1 : 0;
//    xp = (bounds1[1] > bounds2[1]) ? 1 : 0;
//    yp = (bounds1[3] > bounds2[3]) ? 1 : 0;
//    zp = (bounds1[5] > bounds2[5]) ? 1 : 0;;
//    xm = (xn == false && xp == false) ? 1 : 0;
//    ym = (yn == false && yp == false) ? 1 : 0;
//    zm = (zn == false && zp == false) ? 1 : 0;
//    // if part of central chunk, go to next atomic site
//    if (xm == true && ym == true && zm == true) continue;
//    xm = (xn == false || xp == false) ? 1 : 0;
//    ym = (yn == false || yp == false) ? 1 : 0;
//    zm = (zn == false || zp == false) ? 1 : 0;
//    // if not part of central chunk check for overlap with surrounding 26 chunks
//    if (xm == true && ym == true && zn == true) generatesinglesite(i, pointid, 0, 0, -1);
//    if (xm == true && ym == true && zp == true) generatesinglesite(i, pointid, 0, 0, 1);
//    if (xm == true && yn == true && zm == true) generatesinglesite(i, pointid, 0, -1, 0);
//    if (xm == true && yp == true && zm == true) generatesinglesite(i, pointid, 0, 1, 0);
//    if (xn == true && ym == true && zm == true) generatesinglesite(i, pointid, -1, 0, 0);
//    if (xp == true && ym == true && zm == true) generatesinglesite(i, pointid, 1, 0, 0);
//    if (xm == true && yn == true && zn == true) generatesinglesite(i, pointid, 0, -1, -1);
//    if (xm == true && yn == true && zp == true) generatesinglesite(i, pointid, 0, -1, 1);
//    if (xm == true && yp == true && zn == true) generatesinglesite(i, pointid, 0, 1, -1);
//    if (xm == true && yp == true && zp == true) generatesinglesite(i, pointid, 0, 1, 1);
//    if (xn == true && ym == true && zn == true) generatesinglesite(i, pointid, -1, 0, -1);
//    if (xn == true && ym == true && zp == true) generatesinglesite(i, pointid, -1, 0, 1);
//    if (xp == true && ym == true && zn == true) generatesinglesite(i, pointid, 1, 0, -1);
//    if (xp == true && ym == true && zp == true) generatesinglesite(i, pointid, 1, 0, 1);
//    if (xn == true && yn == true && zm == true) generatesinglesite(i, pointid, -1, -1, 0);
//    if (xn == true && yp == true && zm == true) generatesinglesite(i, pointid, -1, 1, 0);
//    if (xp == true && yn == true && zm == true) generatesinglesite(i, pointid, 1, -1, 0);
//    if (xp == true && yp == true && zm == true) generatesinglesite(i, pointid, 1, 1, 0);
//    if (xn == true && yn == true && zn == true) generatesinglesite(i, pointid, -1, -1, -1);
//    if (xn == true && yn == true && zp == true) generatesinglesite(i, pointid, -1, -1, 1);
//    if (xn == true && yp == true && zn == true) generatesinglesite(i, pointid, -1, 1, -1);
//    if (xp == true && yn == true && zn == true) generatesinglesite(i, pointid, 1, -1, -1);
//    if (xn == true && yp == true && zp == true) generatesinglesite(i, pointid, -1, 1, 1);
//    if (xp == true && yn == true && zp == true) generatesinglesite(i, pointid, 1, -1, 1);
//    if (xp == true && yp == true && zn == true) generatesinglesite(i, pointid, 1, 1, -1);
//    if (xp == true && yp == true && zp == true) generatesinglesite(i, pointid, 1, 1, 1);
//  }
//}
//
//// correction for the spatial sites
//void energy::interfaces::aco::aco_ff::generatesinglesite(int & i, std::vector<int> & pointid, const int idix, const int idiy, const int idiz)
//{
//  int temp[3], q;
//  temp[0] = pointid[0] + idix;
//  if (temp[0] < 1) temp[0] = coords->pme.pmetemp.nrough1;
//  if (temp[0] > coords->pme.pmetemp.nrough1) temp[0] = 1;
//  temp[1] = pointid[1] + idiy;
//  if (temp[1] < 1) temp[1] = coords->pme.pmetemp.nrough2;
//  if (temp[1] > coords->pme.pmetemp.nrough2) temp[1] = 1;
//  temp[2] = pointid[2] + idiz;
//  if (temp[2] < 1) temp[2] = coords->pme.pmetemp.nrough3;
//  if (temp[2] > coords->pme.pmetemp.nrough3) temp[2] = 1;
//  q = (temp[2] - 1)*coords->pme.pmetemp.nrough1*coords->pme.pmetemp.nrough2 + (temp[1] - 1)*coords->pme.pmetemp.nrough1 + temp[0];
//  coords->pme.pmetemp.parallelpme(i, q - 1) = 1;
//}
//
//
//void energy::interfaces::aco::aco_ff::sitecorrection(int & offset, int & npoints, int & nspat, int & ab1, int & ab2, int & cb1, int & cb2)
//{
//  offset = 0;
//  if (nspat != 1)
//  {
//    if (ab1 < cb1 || ab2 > cb2)
//    {
//      if (ab1 < 1 || ab2 > npoints)
//      {
//        if (cb1 == 1)
//        {
//          offset = 1 - ab1;
//          ab1 = 1;
//        }
//        else if (cb2 == npoints)
//        {
//          ab2 = npoints;
//          ab1 = ab1 + npoints;
//        }
//      }
//      else
//      {
//        if (cb1 > ab1)
//        {
//          offset = cb1 - ab1;
//          ab1 = cb1;
//        }
//        else ab2 = cb2;
//      }
//    }
//  }
//  offset = offset + 1 - ab1;
//}
//
//
//#ifndef _OPENMP
//
//// put charges on PME grid, serial version
//void energy::interfaces::aco::aco_ff::setchargestogrid()
//{
//  int chargei[3], temp[3];
//  int chargebounds[6], atombounds[6];
//  int offx(0.0), offy(0.0), offz(0.0), m;
//  double x0, y0, z0, term;
//  int i, j, k, i1, j1, k1, spat, atom;
//  // set all grid values to zero
//  for (int i = 0; i < coords->pme.pmetemp.nxpoints; i++)
//  {
//    for (int j = 0; j < coords->pme.pmetemp.nypoints; j++)
//    {
//      for (int k = 0; k < coords->pme.pmetemp.nzpoints; k++)
//      {
//        coords->pme.pmetemp.charges(i, j, k) = 0.0;
//        coords->pme.pmetemp.in[k + coords->pme.pmetemp.nxpoints * (j + coords->pme.pmetemp.nxpoints * i)][0] = 0.0;
//        coords->pme.pmetemp.in[k + coords->pme.pmetemp.nxpoints * (j + coords->pme.pmetemp.nxpoints * i)][1] = 0.0;
//      }
//    }
//  }
//  // loop over all atoms
//  for (atom = 0; atom < coords->pme.pmetemp.natoms; atom++)
//  {
//    // loop over all Spline coeffcients
//    for (i1 = 1; i1 <= coords->pme.pmetemp.bsplineorder; i1++)
//    {
//      k = i1 + coords->pme.pmetemp.initgrid(2, atom) + 1;
//      m = i1;
//      if (k < 1) k = k + coords->pme.pmetemp.nzpoints;
//      x0 = coords->pme.pmetemp.bscz(0, m - 1, atom) * coords->pme.pmetemp.atomcharges[atom];
//      for (j1 = 1; j1 <= coords->pme.pmetemp.bsplineorder; j1++)
//      {
//        j = j1 + coords->pme.pmetemp.initgrid(1, atom) + 1;
//        m = j1;
//        if (j < 1) j = j + coords->pme.pmetemp.nypoints;
//        y0 = coords->pme.pmetemp.bscy(0, m - 1, atom);
//        term = x0 * y0;
//        for (k1 = 1; k1 <= coords->pme.pmetemp.bsplineorder; k1++)
//        {
//          i = k1 + coords->pme.pmetemp.initgrid(0, atom) + 1;
//          m = k1;
//          if (i < 1) i = i + coords->pme.pmetemp.nxpoints;
//          z0 = coords->pme.pmetemp.bscx(0, m - 1, atom);
//          coords->pme.pmetemp.charges(i-1, j-1, k-1) += term*z0;
//          coords->pme.pmetemp.in[(k-1) + coords->pme.pmetemp.nxpoints * ((j-1 )+ coords->pme.pmetemp.nxpoints * (i-1))][0] += term*z0;
//        }
//      }
//    }
//  }// end of all atoms loop;
//}
//
//
//
//#else
//
////put atomic charges on the pme grid, parallel version
//void energy::interfaces::aco::aco_ff::setchargestogrid()
//{
//  int chargei[3], temp[3];
//  int chargebounds[6], atombounds[6];
//  int offx(0.0), offy(0.0), offz(0.0), m;
//  double x0, y0, z0, term;
//  int i, j, k, i1, j1, k1, spat, atom;
//  // set all grid values to zero
//  for (int i = 0; i < coords->pme.pmetemp.nxpoints; i++)
//  {
//    for (int j = 0; j < coords->pme.pmetemp.nypoints; j++)
//    {
//      for (int k = 0; k < coords->pme.pmetemp.nzpoints; k++)
//      {
//        coords->pme.pmetemp.charges(i, j, k) = 0.0;
//        coords->pme.pmetemp.in[k + coords->pme.pmetemp.nxpoints * (j + coords->pme.pmetemp.nxpoints * i)][0] = 0.0;
//        coords->pme.pmetemp.in[k + coords->pme.pmetemp.nxpoints * (j + coords->pme.pmetemp.nxpoints * i)][1] = 0.0;
//      }
//    }
//  }
////  std::cout << coords->pme.pmetemp.rgridtotal << std::endl;
//#pragma omp parallel private (i, j, k, m, i1, j1, k1, spat, atom, chargei, temp, chargebounds, atombounds, offx, offy, offz, x0, y0, z0, term) 
//  {
//    // loop over all spatial sites
//    #pragma omp for schedule(static,1)
//    for (spat = 0; spat < coords->pme.pmetemp.rgridtotal; spat++)
//    {
//      chargei[0] = (spat) % coords->pme.pmetemp.nrough1;
//      chargei[1] = ((spat - chargei[0]) / coords->pme.pmetemp.nrough1) % coords->pme.pmetemp.nrough2;
//      chargei[2] = (spat / (coords->pme.pmetemp.nrough1*coords->pme.pmetemp.nrough2)) % coords->pme.pmetemp.nrough3;
//      chargebounds[0] = chargei[0] * coords->pme.pmetemp.rgrid1 + 1;
//      chargebounds[1] = chargebounds[0] + coords->pme.pmetemp.rgrid1 - 1;
//      chargebounds[2] = chargei[1] * coords->pme.pmetemp.rgrid2 + 1;
//      chargebounds[3] = chargebounds[2] + coords->pme.pmetemp.rgrid2 - 1;
//      chargebounds[4] = chargei[2] * coords->pme.pmetemp.rgrid3 + 1;
//      chargebounds[5] = chargebounds[4] + coords->pme.pmetemp.rgrid3 - 1;
//      // loop over all atoms
//      for (atom = 0; atom < coords->pme.pmetemp.natoms; atom++)
//      {
//        if (coords->pme.pmetemp.parallelpme(atom, spat) == 1)
//        {
//          temp[0] = coords->pme.pmetemp.initgrid(0, atom) + coords->pme.pmetemp.bsoffset;
//          temp[1] = coords->pme.pmetemp.initgrid(1, atom) + coords->pme.pmetemp.bsoffset;
//          temp[2] = coords->pme.pmetemp.initgrid(2, atom) + coords->pme.pmetemp.bsoffset;
//          atombounds[0] = temp[0] - coords->pme.pmetemp.roughleft;
//          atombounds[1] = temp[0] + coords->pme.pmetemp.roughright;
//          atombounds[2] = temp[1] - coords->pme.pmetemp.roughleft;
//          atombounds[3] = temp[1] + coords->pme.pmetemp.roughright;
//          atombounds[4] = temp[2] - coords->pme.pmetemp.roughleft;
//          atombounds[5] = temp[2] + coords->pme.pmetemp.roughright;
//          sitecorrection(offx, coords->pme.pmetemp.nxpoints, coords->pme.pmetemp.nrough1, atombounds[0], atombounds[1], chargebounds[0], chargebounds[1]);
//          sitecorrection(offy, coords->pme.pmetemp.nypoints, coords->pme.pmetemp.nrough2, atombounds[2], atombounds[3], chargebounds[2], chargebounds[3]);
//          sitecorrection(offz, coords->pme.pmetemp.nzpoints, coords->pme.pmetemp.nrough3, atombounds[4], atombounds[5], chargebounds[4], chargebounds[5]);
//          for (i1 = atombounds[4]; i1 <= atombounds[5]; i1++)
//          {
//            k = i1;
//            m = k + offz;
//            if (k < 1) k = k + coords->pme.pmetemp.nzpoints;
//            x0 = coords->pme.pmetemp.bscz(0, m - 1, atom) * coords->pme.pmetemp.atomcharges[atom];
//            for (j1 = atombounds[2]; j1 <= atombounds[3]; j1++)
//            {
//              j = j1;
//              m = j + offy;
//              if (j < 1) j = j + coords->pme.pmetemp.nypoints;
//              y0 = coords->pme.pmetemp.bscy(0, m - 1, atom);
//              term = x0 * y0;
//              for (k1 = atombounds[0]; k1 <= atombounds[1]; k1++)
//              {
//                i = k1;
//                m = i + offx;
//                if (i < 1) i = i + coords->pme.pmetemp.nxpoints;
//                z0 = coords->pme.pmetemp.bscx(0, m - 1, atom);
//                coords->pme.pmetemp.charges(i - 1, j - 1, k - 1) = coords->pme.pmetemp.charges(i - 1, j - 1, k - 1) + term*z0;
//                coords->pme.pmetemp.in[(k-1) + coords->pme.pmetemp.nxpoints * ((j-1) + coords->pme.pmetemp.nxpoints * (i-1))][0] += term*z0;
//              }
//            }
//          }
//        } // end of if clause
//      }// end of all atoms loop
//    }
//  }//end of parallel section
//}
//
//
//#endif
