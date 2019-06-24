#include "coords.h"
#include "scon_angle.h"
#include"helperfunctions.h"


#include <algorithm>
#include <cmath>

/* ######################################################


  ########  ####    ###     ######
  ##     ##  ##    ## ##   ##    ##
  ##     ##  ##   ##   ##  ##
  ########   ##  ##     ##  ######
  ##     ##  ##  #########       ##
  ##     ##  ##  ##     ## ##    ##
  ########  #### ##     ##  ######


###################################################### */

void coords::bias::Potentials::append_config()
{
  m_distances.insert(m_distances.end(),
    Config::get().coords.bias.distance.begin(),
    Config::get().coords.bias.distance.end());
  m_angles.insert(m_angles.end(),
    Config::get().coords.bias.angle.begin(),
    Config::get().coords.bias.angle.end());
  m_dihedrals.insert(m_dihedrals.end(),
    Config::get().coords.bias.dihedral.begin(),
    Config::get().coords.bias.dihedral.end());
  m_angles.insert(m_angles.end(),
    Config::get().coords.bias.angle.begin(),
    Config::get().coords.bias.angle.end());
  m_utors.insert(m_utors.end(),
    Config::get().coords.bias.utors.begin(),
    Config::get().coords.bias.utors.end());
  m_ucombs.insert(m_ucombs.end(),
    Config::get().coords.bias.ucombs.begin(),
    Config::get().coords.bias.ucombs.end());
  m_thresh.insert(m_thresh.end(),
    Config::get().coords.bias.threshold.begin(),
    Config::get().coords.bias.threshold.end());
  m_threshBottom.insert(m_threshBottom.end(),
    Config::get().coords.bias.thresholdBottom.begin(),
    Config::get().coords.bias.thresholdBottom.end());
}

void coords::bias::Potentials::swap(Potentials & rhs)
{
  std::swap(b, rhs.b);
  std::swap(a, rhs.a);
  std::swap(d, rhs.d);
  std::swap(s, rhs.s);
  std::swap(c, rhs.c);
  std::swap(u, rhs.u);
  std::swap(thr, rhs.thr);
  m_dihedrals.swap(rhs.m_dihedrals);
  m_angles.swap(rhs.m_angles);
  m_distances.swap(rhs.m_distances);
  m_spherical.swap(rhs.m_spherical);
  m_cubic.swap(rhs.m_cubic);
  m_utors.swap(rhs.m_utors);
  m_uangles.swap(rhs.m_uangles);
  m_udist.swap(rhs.m_udist);
  m_ucombs.swap(rhs.m_ucombs);
  m_thresh.swap(rhs.m_thresh);
  m_threshBottom.swap(rhs.m_threshBottom);
}

coords::bias::Potentials::Potentials()
  : b{}, a{}, d{}, s{}, c{}, u{},
  m_dihedrals{Config::get().coords.bias.dihedral},
  m_angles{Config::get().coords.bias.angle},
  m_distances{Config::get().coords.bias.distance},
  m_spherical{Config::get().coords.bias.spherical},
  m_cubic{Config::get().coords.bias.cubic},
  m_thresh{Config::get().coords.bias.threshold},
  m_utors{Config::get().coords.bias.utors},
  m_uangles{ Config::get().coords.bias.uangles },
  m_udist{Config::get().coords.bias.udist},
  m_threshBottom {Config::get().coords.bias.thresholdBottom},
  m_ucombs{ Config::get().coords.bias.ucombs }
{ }

bool coords::bias::Potentials::empty() const
{
  return scon::empty(m_dihedrals, m_angles, m_distances,
    m_spherical, m_cubic, m_utors, m_udist, m_thresh, m_threshBottom, m_ucombs);
}

double coords::bias::Potentials::apply(Representation_3D const & xyz,
  Gradients_3D & g_xyz, Cartesian_Point maxPos, Cartesian_Point minPos, Cartesian_Point const & center)
{
  if (!m_dihedrals.empty()) 
    d = dih(xyz, g_xyz);
  if (!m_angles.empty()) 
    a = ang(xyz, g_xyz);
  if (!m_distances.empty()) 
    b = dist(xyz, g_xyz);
  if (!m_spherical.empty()) 
    s = spherical(xyz, g_xyz, center);
  if (!m_cubic.empty()) 
    c = cubic(xyz, g_xyz, center);
  if( !m_thresh.empty())
    thr = thresh(xyz, g_xyz, maxPos);
  if (Config::set().coords.umbrella.use_comb && !m_ucombs.empty())
    u = umbrellacomb(xyz, g_xyz);
  return b + a + d + s + c + u;
  if (!m_threshBottom.empty())
    thrB = thresh_bottom(xyz, g_xyz, minPos);
  return b + a + d + s + c;
}

/**apply umbrella potentials and save data for 'umbrella.txt' into uout*/
void coords::bias::Potentials::umbrellaapply(Representation_3D const & xyz,
  Gradients_3D & g_xyz, std::vector<double> &uout) 
{
  if (!m_utors.empty()) // torsion restraints
    umbrelladih(xyz, g_xyz, uout);
  if (!m_uangles.empty())   // angle restraints
    umbrellaang(xyz, g_xyz, uout);
  if (!m_udist.empty()) // distance restraints
    umbrelladist(xyz, g_xyz, uout);
  if (!m_ucombs.empty()) // restraints of combined distances
    umbrellacomb(xyz, g_xyz, uout);
}

double coords::bias::Potentials::dih(Representation_3D const &positions,
  Gradients_3D & gradients)
{
  float_type dih_bias_energy = 0;
  for (auto &dih : m_dihedrals)
  {
    // distances
    auto const ba = positions[dih.b] - positions[dih.a];
    auto const cb = positions[dih.c] - positions[dih.b];
    auto const dc = positions[dih.d] - positions[dih.c];
    // cross products
    auto const t = cross(ba, cb);
    auto const u = cross(cb, dc);
    auto const tu = cross(t, u);
    // lengths
    auto const rt2 = dot(t, t);
    auto const ru2 = dot(u, u);
    auto const  rtru = std::sqrt(rt2*ru2);
    // check for length of rtru
    if (std::abs(rtru) > 0)
    {
      auto const rcb = len(cb);
      // sine and cosine of current angle:
      auto const cosine = dot(t, u) / rtru;
      auto const sine = dot(cb, tu) / (rcb*rtru);
      // sine and cosine of "ideal":
      auto const c1 = cos(dih.ideal);
      auto const s1 = sin(dih.ideal);
      // energy of this bias
      auto const e = dih.force * (1 - (cosine*c1 + sine*s1));
      dih_bias_energy += e;
      // set value
      dih.value = coords::angle_type::from_rad(sine > 0 ?
        std::acos(cosine) : -std::acos(cosine));
      // derivative (f*cos(p0-p))' = -f*sin(p0-p)
      auto const dphi = -dih.force * (cosine*s1 - sine*c1);
      // per atom derivatives
      auto const dt = cross(t, cb)*(dphi / (rt2*rcb));
      auto const du = cross(cb, u)*(dphi / (ru2*rcb));
      auto const d_a = cross(dt, cb);
      auto const d_b = cross(positions[dih.c] - positions[dih.a], dt) + cross(du, dc);
      auto const d_c = cross(dt, ba) + cross(positions[dih.d] - positions[dih.b], du);
      auto const d_d = cross(du, cb);
      // add to gradient
      gradients[dih.a] += d_a;
      gradients[dih.b] += d_b;
      gradients[dih.c] += d_c;
      gradients[dih.d] += d_d;
    } // if rtru != 0
  } // for all m_dihedrals
  return dih_bias_energy;
}

void coords::bias::Potentials::umbrelladih(Representation_3D const &positions,
  Gradients_3D & gradients, std::vector<double> &uout) const
{
  using std::abs;
  for (auto const &dih : m_utors)
  {
    float_type dE(0.0), diff(0.0);
    Cartesian_Point const b01(positions[dih.index[1]] - positions[dih.index[0]]);
    Cartesian_Point const b12(positions[dih.index[2]] - positions[dih.index[1]]);
    Cartesian_Point const b23(positions[dih.index[3]] - positions[dih.index[2]]);
    Cartesian_Point const b02(positions[dih.index[2]] - positions[dih.index[0]]);
    Cartesian_Point const b13(positions[dih.index[3]] - positions[dih.index[1]]);
    Cartesian_Point t(cross(b01, b12));
    Cartesian_Point u(cross(b12, b23));
    float_type const tl2(dot(t, t));
    float_type const ul2(dot(u, u));
    float_type const r12(geometric_length(b12));
    Cartesian_Point const tu(cross(t, u));
    float_type torsion = angle(t, u).degrees();
    float_type const norm = r12*geometric_length(tu);
    torsion = (abs(norm) > float_type(0) && (dot(b12, tu) / norm) < 0.0) ? -torsion : torsion;
    // Apply half harmonic bias potential according to torsion value
    if (dih.angle > 0) {
      if (torsion > 0) {
        diff = torsion - dih.angle;
        if (diff < -180) {
          diff = diff + 360;
        }
        else if (diff > 180) {
          diff = diff - 360;
        }
        dE = dih.force * diff * SCON_180PI;
      }
      else {
        diff = torsion - dih.angle;
        if (diff < -180) {
          diff = diff + 360;
        }
        else if (diff > 180) {
          diff = diff - 360;
        }
        dE = dih.force * diff * SCON_180PI;
      }
    }// end of angle > 0
    else if (dih.angle < 0) {
      if (torsion < 0) {
        diff = torsion - dih.angle;
        if (diff < -180) {
          diff = diff + 360;
        }
        else if (diff > 180) {
          diff = diff - 360;
        }
        dE = dih.force * diff * SCON_180PI;
      }
      else {
        diff = torsion - dih.angle;
        if (diff > 180) {
          diff = diff - 360;
        }
        else if (diff < -180) {
          diff = diff + 360;
        }
        dE = dih.force * diff * SCON_180PI;
      }
    }// end of angle < 0
    else if (dih.angle == 0) {
      diff = torsion;
      dE = dih.force * diff * SCON_180PI;
    }
    // Derivatives
    t = cross(t, b12);
    t *= dE / (tl2*r12);
    u = cross(u, b12);
    u *= -dE / (ul2*r12);
    // fill udatacontainer
    double diff_tors = std::abs(torsion - dih.angle);
    double diff_tors_m = std::abs(torsion - 360 - dih.angle);
    double diff_tors_p = std::abs(torsion + 360 - dih.angle);
    if (diff_tors_m < diff_tors) uout.push_back(torsion - 360);
    else if (diff_tors_p < diff_tors) uout.push_back(torsion + 360);
    else uout.push_back(torsion);

    if (Config::get().general.verbosity >= 4)
    {
      std::cout << "U W: " << torsion << ", AW: " << diff << ", SOLL: " << dih.angle;
      std::cout << "; FORCE: " << dih.force << ", DE: " << dE << ", PW: " << uout.back() << std::endl;
      std::cout << dih.index[0] << "   " << dih.index[1] << "   " << dih.index[2] << "   " << dih.index[3] << std::endl;
    }
    gradients[dih.index[0]] += cross(t, b12);
    gradients[dih.index[1]] += cross(b02, t) + cross(u, b23);
    gradients[dih.index[2]] += cross(t, b01) + cross(b13, u);
    gradients[dih.index[3]] += cross(u, b12);
  }
}

void coords::bias::Potentials::umbrellaang(Representation_3D const & xyz, Gradients_3D & g_xyz, std::vector<double>& uout) const
{
  for (auto const &angle : m_uangles)
  {
    // calculate current angle
    auto pos_a = xyz[angle.index[0]];
    auto pos_b = xyz[angle.index[1]];
    auto pos_c = xyz[angle.index[2]];
    auto vec1 = pos_a - pos_b;
    auto vec2 = pos_c - pos_b;
    auto current_angle = scon::angle(vec1, vec2).degrees();
    uout.push_back(current_angle);        // fill uout

    // calculate gradients
    auto diff = (current_angle - angle.angle)*SCON_PI180;     // difference to ideal angle in rad
    auto d1 = geometric_length(vec1);
    auto d2 = geometric_length(vec2);
    auto scalar_product = scon::dot(vec1, vec2);
    auto prefactor = angle.force * diff * (-1.0) / (std::sqrt(1 - std::cos(current_angle*SCON_PI180) * std::cos(current_angle*SCON_PI180)));

    auto grad_a_x = prefactor * (((pos_b.x() - pos_a.x())*scalar_product) / (d2 * std::pow(d1, 1.5)) + (pos_c.x() - pos_b.x()) / (d1*d2));
    auto grad_a_y = prefactor * (((pos_b.y() - pos_a.y())*scalar_product) / (d2 * std::pow(d1, 1.5)) + (pos_c.y() - pos_b.y()) / (d1*d2));
    auto grad_a_z = prefactor * (((pos_b.z() - pos_a.z())*scalar_product) / (d2 * std::pow(d1, 1.5)) + (pos_c.z() - pos_b.z()) / (d1*d2));

    auto grad_b_x = prefactor * ((pos_a.x() - pos_b.x()) *scalar_product / (d2*std::pow(d1, 1.5))
      + (pos_c.x() - pos_b.x())*scalar_product / (d1*std::pow(d2, 1.5)) + (2 * pos_b.x() - pos_a.x() - pos_c.x()) / (d1*d2));
    auto grad_b_y = prefactor * ((pos_a.y() - pos_b.y()) *scalar_product / (d2*std::pow(d1, 1.5))
      + (pos_c.y() - pos_b.y())*scalar_product / (d1*std::pow(d2, 1.5)) + (2 * pos_b.y() - pos_a.y() - pos_c.y()) / (d1*d2));
    auto grad_b_z = prefactor * ((pos_a.z() - pos_b.z()) *scalar_product / (d2*std::pow(d1, 1.5))
      + (pos_c.z() - pos_b.z())*scalar_product / (d1*std::pow(d2, 1.5)) + (2 * pos_b.z() - pos_a.z() - pos_c.z()) / (d1*d2));

    auto grad_c_x = prefactor * ((pos_a.x() - pos_b.x()) / (d1*d2) + (pos_b.x() - pos_c.x()) * scalar_product / (d1*std::pow(d2, 1.5)));
    auto grad_c_y = prefactor * ((pos_a.y() - pos_b.y()) / (d1*d2) + (pos_b.y() - pos_c.y()) * scalar_product / (d1*std::pow(d2, 1.5)));
    auto grad_c_z = prefactor * ((pos_a.z() - pos_b.z()) / (d1*d2) + (pos_b.z() - pos_c.z()) * scalar_product / (d1*std::pow(d2, 1.5)));

    if (Config::get().general.verbosity > 4)
    {
      std::cout << "bias gradients\n";
      std::cout << grad_a_x << " , " << grad_a_y << " , " << grad_a_z << "\n";
      std::cout << grad_b_x << " , " << grad_b_y << " , " << grad_b_z << "\n";
      std::cout << grad_c_x << " , " << grad_c_y << " , " << grad_c_z << "\n";
    }

    // add bias gradients
    g_xyz[angle.index[0]] += coords::r3(grad_a_x, grad_a_y, grad_a_z);
    g_xyz[angle.index[1]] += coords::r3(grad_b_x, grad_b_y, grad_b_z);
    g_xyz[angle.index[2]] += coords::r3(grad_c_x, grad_c_y, grad_c_z);
  }
}

void coords::bias::Potentials::umbrelladist(Representation_3D const &positions,
  Gradients_3D & gradients, std::vector<double> &uout) const
{
  for (auto const &dist : m_udist)
  {
    coords::Cartesian_Point bv(positions[dist.index[0]] - positions[dist.index[1]]);
    float_type md = geometric_length(bv);
    uout.push_back(md);                       // fill value for restraint into uout
    float_type diff(md - dist.dist);
    float_type dE = dist.force * diff / md;
    Cartesian_Point gv = bv*dE;
    gradients[dist.index[0]] += gv;           //apply half harmonic potential on gradients
    gradients[dist.index[1]] -= gv;
  }
}

void coords::bias::Potentials::umbrellacomb(Representation_3D const &positions,
  Gradients_3D & gradients, std::vector<double> &uout)
{
  for (auto &comb : set_ucombs())   // for every restraint combination
  {
    // raise force constant in the first half of equilibration
    if (is_smaller_than(comb.force_current, comb.force_final)){
      comb.force_current += (comb.force_final * 2) / Config::get().md.usequil;
    }

    // calculate value for 'umbrella.txt' file and save it
    double reactioncoord{ 0.0 };     // total reaction coordinate
    for (auto &d : comb.dists)   // for every distance in combination
    {
      coords::Cartesian_Point vec(positions[d.index1] - positions[d.index2]); // vector between atoms
      double distance = geometric_length(vec);                                // distance between atoms
      reactioncoord += distance * d.factor;                                   // add to reaction coordinate
    }
    uout.push_back(reactioncoord);    // fill value for restraint into uout

		// Console output
		if (Config::get().general.verbosity > 3)
		{
			std::cout << "Force constant: " << comb.force_current << ", reaction coordinate: " << reactioncoord << "\n";
		}
    
    // apply bias potential on gradients
    float_type diff(reactioncoord - comb.value);        // difference to desired value (restraint)
    for (auto &d : comb.dists)   // for every distance in combination
    {
      coords::Cartesian_Point vec(positions[d.index1] - positions[d.index2]); // vector between atoms (r1-r2)
      double distance = geometric_length(vec);                                // distance between atoms

      gradients[d.index1].x() += comb.force_current * diff * d.factor * (1.0 / distance) * vec.x();
      gradients[d.index1].y() += comb.force_current * diff * d.factor * (1.0 / distance) * vec.y();
      gradients[d.index1].z() += comb.force_current * diff * d.factor * (1.0 / distance) * vec.z();

      gradients[d.index2].x() -= comb.force_current * diff * d.factor * (1.0 / distance) * vec.x();
      gradients[d.index2].y() -= comb.force_current * diff * d.factor * (1.0 / distance) * vec.y();
      gradients[d.index2].z() -= comb.force_current * diff * d.factor * (1.0 / distance) * vec.z();
    }
  }
}

double coords::bias::Potentials::umbrellacomb(Representation_3D const &positions, Gradients_3D & gradients)
{
  double bias_energy{ 0.0 };
  for (auto &comb : set_ucombs())   // for every restraint combination
  {
    // calculate value for reactioncoord
    double reactioncoord{ 0.0 };     // total reaction coordinate
    for (auto &d : comb.dists)   // for every distance in combination
    {
      coords::Cartesian_Point vec(positions[d.index1] - positions[d.index2]); // vector between atoms
      double distance = geometric_length(vec);                                // distance between atoms
      reactioncoord += distance * d.factor;                                   // add to reaction coordinate
    }

    // apply bias potential on gradients
    float_type diff(reactioncoord - comb.value);        // difference to desired value (restraint)
    for (auto &d : comb.dists)   // for every distance in combination
    {
      coords::Cartesian_Point vec(positions[d.index1] - positions[d.index2]); // vector between atoms (r1-r2)
      double distance = geometric_length(vec);                                // distance between atoms

      gradients[d.index1].x() += comb.force_final * diff * d.factor * (1.0 / distance) * vec.x();
      gradients[d.index1].y() += comb.force_final * diff * d.factor * (1.0 / distance) * vec.y();
      gradients[d.index1].z() += comb.force_final * diff * d.factor * (1.0 / distance) * vec.z();

      gradients[d.index2].x() -= comb.force_final * diff * d.factor * (1.0 / distance) * vec.x();
      gradients[d.index2].y() -= comb.force_final * diff * d.factor * (1.0 / distance) * vec.y();
      gradients[d.index2].z() -= comb.force_final * diff * d.factor * (1.0 / distance) * vec.z();
    }
    bias_energy += 0.5 * comb.force_final * diff * diff;
  }
  return bias_energy;
}

double coords::bias::Potentials::dist(Representation_3D const &positions, Gradients_3D &gradients)
{
  double E(0.0);
  for (auto &distance : m_distances) 
  {
    auto bv = positions[distance.a] - positions[distance.b];
    auto md = len(bv);
    distance.value = md;
    auto diff = md - distance.ideal;
    //apply half harmonic potential
    auto e = distance.force * diff;
    E += e*diff;
    auto dE = e * 2;
    auto gv = bv*(dE / md);
    gradients[distance.a] += gv;
    gradients[distance.b] -= gv;
  }
  return E;
}

double coords::bias::Potentials::ang(Representation_3D const &positions, Gradients_3D &gradients)
{
	double E(0.0);
	for (auto &angle : m_angles)
	{
    // calculate energy
    auto pos_a = positions[angle.a];
    auto pos_b = positions[angle.b];
    auto pos_c = positions[angle.c];
		auto vec1 = pos_a - pos_b;
		auto vec2 = pos_c - pos_b;

		angle.value = scon::angle(vec1, vec2).degrees();        // get current angle in degrees
		auto diff = (angle.value - angle.ideal)*SCON_PI180;     // difference to ideal angle in rad
		E += angle.force * diff * diff;                         // apply half harmonic potential

    // calculate gradients
    auto d1 = geometric_length(vec1);
    auto d2 = geometric_length(vec2);
    auto scalar_product = scon::dot(vec1, vec2);
    auto prefactor = 2 * angle.force * diff * (-1.0) / (std::sqrt(1 - std::cos(angle.value*SCON_PI180) * std::cos(angle.value*SCON_PI180)));  

    auto grad_a_x = prefactor * (((pos_b.x() - pos_a.x())*scalar_product) / (d2 * std::pow(d1, 1.5)) + (pos_c.x() - pos_b.x()) / (d1*d2));
    auto grad_a_y = prefactor * (((pos_b.y() - pos_a.y())*scalar_product) / (d2 * std::pow(d1, 1.5)) + (pos_c.y() - pos_b.y()) / (d1*d2));
    auto grad_a_z = prefactor * (((pos_b.z() - pos_a.z())*scalar_product) / (d2 * std::pow(d1, 1.5)) + (pos_c.z() - pos_b.z()) / (d1*d2));

    auto grad_b_x = prefactor * ( (pos_a.x() - pos_b.x()) *scalar_product / (d2*std::pow(d1, 1.5)) 
      + (pos_c.x() - pos_b.x())*scalar_product / (d1*std::pow(d2, 1.5)) + (2 * pos_b.x() - pos_a.x() - pos_c.x()) / (d1*d2));
    auto grad_b_y = prefactor * ((pos_a.y() - pos_b.y()) *scalar_product / (d2*std::pow(d1, 1.5))
      + (pos_c.y() - pos_b.y())*scalar_product / (d1*std::pow(d2, 1.5)) + (2 * pos_b.y() - pos_a.y() - pos_c.y()) / (d1*d2));
    auto grad_b_z = prefactor * ((pos_a.z() - pos_b.z()) *scalar_product / (d2*std::pow(d1, 1.5))
      + (pos_c.z() - pos_b.z())*scalar_product / (d1*std::pow(d2, 1.5)) + (2 * pos_b.z() - pos_a.z() - pos_c.z()) / (d1*d2));

    auto grad_c_x = prefactor * ((pos_a.x() - pos_b.x()) / (d1*d2) + (pos_b.x() - pos_c.x()) * scalar_product / (d1*std::pow(d2, 1.5)));
    auto grad_c_y = prefactor * ((pos_a.y() - pos_b.y()) / (d1*d2) + (pos_b.y() - pos_c.y()) * scalar_product / (d1*std::pow(d2, 1.5)));
    auto grad_c_z = prefactor * ((pos_a.z() - pos_b.z()) / (d1*d2) + (pos_b.z() - pos_c.z()) * scalar_product / (d1*std::pow(d2, 1.5)));

    if (Config::get().general.verbosity > 4)
    {
      std::cout << "bias gradients\n";
      std::cout << grad_a_x << " , " << grad_a_y << " , " << grad_a_z << "\n";
      std::cout << grad_b_x << " , " << grad_b_y << " , " << grad_b_z << "\n";
      std::cout << grad_c_x << " , " << grad_c_y << " , " << grad_c_z << "\n";
    }
    
    // add bias gradients
    gradients[angle.a] += coords::r3(grad_a_x, grad_a_y, grad_a_z);
    gradients[angle.b] += coords::r3(grad_b_x, grad_b_y, grad_b_z);
    gradients[angle.c] += coords::r3(grad_c_x, grad_c_y, grad_c_z);
	}
  return E;
}

coords::float_type coords::bias::Potentials::spherical(Representation_3D const &positions,
  Gradients_3D & gradients, Cartesian_Point const & center)
{
  using std::max;
  using std::pow;
  coords::float_type energy_total(0.0);
  std::size_t const numberOfAtoms(positions.size());
  for (auto const & orb : m_spherical)
  {
    if (Config::get().general.verbosity >= 4)
    {
      std::cout << "Applying spherical bias with radius " << orb.radius;
      std::cout << " around " << center << " (exponent = " << orb.exponent << ")\n";
    }
    for (std::size_t i = 0; i < numberOfAtoms; ++i)
    {
      Cartesian_Point sphericalBiasGradient = positions[i] - center;
      coords::float_type const distanceFromCenterOfBias = geometric_length(sphericalBiasGradient);
      if (distanceFromCenterOfBias < orb.radius)
      {
        continue;
      }
      coords::float_type const delta = distanceFromCenterOfBias - orb.radius;
      sphericalBiasGradient *= delta / distanceFromCenterOfBias;
      coords::float_type const exponent = max(1.0, orb.exponent);
      coords::float_type const energy_current = orb.force * pow(delta, exponent);
      energy_total += energy_current;
      coords::float_type dE = orb.force * exponent * pow(delta, exponent - 1.0);
      sphericalBiasGradient *= dE;
      if (Config::get().general.verbosity >= 4)
      {
        std::cout << "Spherical bias gradient of atom " << i + 1 << " at ";
        std::cout << positions[i] << " (which is " << distanceFromCenterOfBias << " from " << center;
        std::cout << " making delta = " << delta << ") is " << sphericalBiasGradient << " with energy " << energy_current << '\n';
      }
      gradients[i] += sphericalBiasGradient;
    }
    if (Config::get().general.verbosity >= 4)
    {
      std::cout << "Spherical Bias Energy = " << energy_total << '\n';
    }
  }
  return energy_total;
}

double coords::bias::Potentials::cubic(Representation_3D const &positions,
  Gradients_3D & gradients, Cartesian_Point const & center)
{
  double totalEnergy(0.0);
  using scon::abs;
  std::size_t const N(positions.size());
  for (auto const & box : m_cubic)
  {
    Cartesian_Point const halfdim(box.dim / 2.);
    if (Config::get().general.verbosity >= 4)
    {
      std::cout << "Applying cubic boundary with side length " << box.dim;
      std::cout << " and halfdim " << halfdim << " around " << center;
      std::cout << " (exponent = " << box.exponent << ")\n";
    }
    for (std::size_t i = 0; i < N; ++i)
    {
      Cartesian_Point const db = positions[i] - center;
      Cartesian_Point const delta = scon::abs(db) - scon::abs(halfdim);
      double const expo = std::max(1.0, box.exponent);
      double const fexpo = box.force * expo;
      double currentEnergy = 0.;
      Cartesian_Point de;
      if (delta.x() > 0.)
      {
        currentEnergy += box.force * std::pow(delta.x(), expo);
        de.x() = (db.x() < 0. ? -fexpo : fexpo) * std::pow(delta.x(), expo - 1.);
        if (Config::get().general.verbosity >= 4)
        {
          std::cout << positions[i].x() << " x out of " << center.x();
          std::cout << " +/- " << halfdim.x() << " by " << delta.x() << " (" << db.x() << ") ";
        }
      }
      if (delta.y() > 0.)
      {
        currentEnergy += box.force * std::pow(delta.y(), expo);
        de.y() = (db.y() < 0. ? -fexpo : fexpo) * std::pow(delta.y(), expo - 1.);
        if (Config::get().general.verbosity >= 4)
        {
          std::cout << positions[i].y() << " y out of " << center.y();
          std::cout << " +/- " << halfdim.y() << " by " << delta.y() << " (" << db.y() << ") ";
        }
      }
      if (delta.z() > 0.)
      {
        currentEnergy += box.force * std::pow(delta.z(), expo);
        de.z() = (db.z() < 0. ? -fexpo : fexpo) * std::pow(delta.z(), expo - 1.);
        if (Config::get().general.verbosity >= 4)
        {
          std::cout << positions[i].z() << " z out of " << center.z();
          std::cout << " +/- " << halfdim.z() << " by " << delta.z() << " (" << db.z() << ") ";
        }
      }
      Cartesian_Point const g = delta*de;
      totalEnergy += currentEnergy;
      if (Config::get().general.verbosity >= 4 && currentEnergy > 0.)
      {
        std::cout << " E = " << currentEnergy << " and grad " << g << '\n';
      }
      gradients[i] += g;
    }
  }
  return totalEnergy;
}

double coords::bias::Potentials::thresh(Representation_3D const &positions, Gradients_3D &gradients, Cartesian_Point maxPos)//special potential for special task layerdeposiotion , who is a very special task. SPECIAL!
{
  double E(0.0);
  std::size_t const N(positions.size());

  for (auto &thresholdstr : m_thresh)
  {
    for (std::size_t i = 0u; i < N; ++i)
    {
      switch(Config::get().layd.laydaxis)
      {
      case 'x':
      {
        if (positions[i].x() > maxPos.x() + thresholdstr.th_dist)
        {
          double force = thresholdstr.forceconstant * (positions[i].x() - (maxPos.x() + thresholdstr.th_dist));
          gradients[i] += force;
        }
        break;
      }
      case 'y':
      {
        if (positions[i].y() > maxPos.y() + thresholdstr.th_dist)
        {
          double force = thresholdstr.forceconstant * (positions[i].y() - (maxPos.y() + thresholdstr.th_dist));
          gradients[i] += force;
        }
        break;
      }
      case 'z':
      {
        if (positions[i].z() > maxPos.z() + thresholdstr.th_dist)
        {
          double force = thresholdstr.forceconstant * (positions[i].z() - (maxPos.z() + thresholdstr.th_dist));
          gradients[i] += force;
        }
        break;
      }
      default:
      {
        throw std::runtime_error("Entered invalid dimension. Only x,y and z possible.");
        break;
      }//default end
      }
    }
  }
  return E;
}

double coords::bias::Potentials::thresh_bottom(Representation_3D const &positions, Gradients_3D &gradients, Cartesian_Point minPos)//special potential for layerdeposition to prevent molecules leaving layerto far
{
  double E(0.0);
  std::size_t const N(positions.size());

  for (auto &thresholdstr : m_thresh)
  {
    for (std::size_t i = 0u; i < N; ++i)
    {
      switch (Config::get().layd.laydaxis)
      {
      case 'x':
      {
        if (positions[i].x() < minPos.x() - thresholdstr.th_dist)
        {
          double force = thresholdstr.forceconstant * (positions[i].x() - (minPos.x() - thresholdstr.th_dist));
          gradients[i] += force;
        }
        break;
      }
      case 'y':
      {
        if (positions[i].y() < minPos.y() - thresholdstr.th_dist)
        {
          double force = thresholdstr.forceconstant * (positions[i].y() - (minPos.y() - thresholdstr.th_dist));
          gradients[i] += force;
        }
        break;
      }
      case 'z':
      {
        if (positions[i].z() < minPos.z() - thresholdstr.th_dist)
        {
          double force = thresholdstr.forceconstant * (positions[i].z() - (minPos.z() - thresholdstr.th_dist));
          gradients[i] += force;
        }
        break;
      }
      default:
      {
        throw std::runtime_error("Entered invalid dimension. Only x,y and z possible.");
        break;
      }//default end
      }
    }
  }
  return E;
}
