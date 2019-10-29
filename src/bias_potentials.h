#pragma once
#include "coords_rep.h"
#include "spline.h"

/* ######################################################


########  ####    ###     ######
##     ##  ##    ## ##   ##    ##
##     ##  ##   ##   ##  ##
########   ##  ##     ##  ######
##     ##  ##  #########       ##
##     ##  ##  ##     ## ##    ##
########  #### ##     ##  ######


###################################################### */

/**namespace that contains only bias potentials*/
namespace coords::bias
{
  /**struct for every kind of bias potentials*/
  struct Potentials
  {
    /**constructor, fetches biases from config object*/
    Potentials();

    /**are there any biases?*/
    bool empty() const;
    /**are there any umbrella restraints?*/
    bool uempty() const { return m_utors.empty() && m_udist.empty() && m_ucombs.empty() && m_uangles.empty(); }
    /**returns bias energy*/
    double energy() const { return b + a + d + s + c; }

    /**clear all biases*/
    void clear()
    {
      b = a = d = s = c = thr = u = 0.0;
      scon::clear(m_dihedrals, m_angles, m_distances,
        m_spherical, m_cubic, m_utors, m_udist, m_ucombs, m_thresh, m_uangles);
    }

    /**energy for distance biases*/
    double e_dist() const { return b; }
    /**energy for angle biases*/
    double e_angle() const { return a; }
    /**energy for dihedral biases*/
    double e_dihedral() const { return d; }
    /**energy for spherical biases*/
    double e_spherical() const { return s; }
    /**energy for cubic biases*/
    double e_cubic() const { return c; }
    /**energy for threshold biases*/
    double e_thresh() const { return thr; }
    /**energy for umbrella combination biases*/
    double e_ucomb() const { return u; }

    /**add a new dihedral bias*/
    void add(config::biases::dihedral const& new_d) { m_dihedrals.push_back(new_d); }
    /**add a new angle bias*/
    void add(config::biases::angle const& new_a) { m_angles.push_back(new_a); }
    /**add a new distance bias*/
    void add(config::biases::distance const& new_d) { m_distances.push_back(new_d); }
    /**add a new spherical bias*/
    void add(config::biases::spherical const& new_d) { m_spherical.push_back(new_d); }
    /**add a new cubic bias*/
    void add(config::biases::cubic const& new_d) { m_cubic.push_back(new_d); }
    /**add a new threshold bias*/
    void add(config::biases::thresholdstr const& new_thr) { m_thresh.push_back(new_thr); }

    /**return all dihedral biases*/
    std::vector<config::biases::dihedral> const& dihedrals() const { return m_dihedrals; }
    /**return all angle biases*/
    std::vector<config::biases::angle> const& angles() const { return m_angles; }
    /**return all distance biases*/
    std::vector<config::biases::distance> const& distances() const { return m_distances; }
    /**return all spherical biases*/
    std::vector<config::biases::spherical> const& sphericals() const { return m_spherical; }
    /**return all cubic biases*/
    std::vector<config::biases::cubic> const& cubic() const { return m_cubic; }
    /**return all threshold biases*/
    std::vector<config::biases::thresholdstr> const& thresholds() const { return m_thresh; }
    /**return all umbrella combination biases*/
    std::vector<config::coords::umbrellas::umbrella_comb> const& ucombs() const { return m_ucombs; }
    /**function to change umbrella combinations
    necessary for raising force constant during equilibration*/
    std::vector<config::coords::umbrellas::umbrella_comb>& set_ucombs() { return m_ucombs; }

    /**apply "normal" biases (dihedrals, angles, distances, spherical, cubic, threshold, if desired umbrella_combs)
    returns bias energy*/
    double apply(Representation_3D const& xyz, Representation_3D& g_xyz,
      Cartesian_Point maxPos, Cartesian_Point minPos, Cartesian_Point const& center = Cartesian_Point());
    /**apply umbrella biases, i.e. apply bias gradients and fill values for reaction coordinate into uout
    @param xyz: cartesian coordinates of molecule
    @param g_xyz: cartesian gradients of molecule (are changed according to bias)
    @param iout: vector of values for umbrella reaction coordinate that are later written into 'umbrella.txt'*/
    void umbrellaapply(Representation_3D const& xyz, Representation_3D& g_xyz, std::vector<double>& uout, Spline const& s);

    /**calculate size of a torsion in degrees
    @param xyz: cartesian coordinates of molecule
    @param dih: atom indices of torsion*/
    static double calc_tors(Representation_3D const& xyz, std::vector<std::size_t> const& dih);

    /**calculate angle (in degree)
    @param xyz: cartesian coordinates of molecule
    @param ang: atom indices of angle*/
    static double calc_angle(Representation_3D const& xyz, std::vector<std::size_t> const& ang);

    /**calculate distance
    @param xyz: cartesian coordinates of molecule
    @param dist: atom indices of distance*/
    static double calc_dist(Representation_3D const& xyz, std::vector<std::size_t> const& dist);

    /**calculate xi value in PMF_IC_PREP
    @param xyz: cartesian coordinates of molecule
    @param indices: indices that define the reaction coordinate xi*/
    static double calc_xi(Representation_3D const& xyz, std::vector<std::size_t> const& indices);

    /**get biases from config*/
    void append_config();

    void swap(Potentials& rhs);

  private:

    //energies of biases
    /**energy of distance bias*/
    double b;
    /**energy of angle bias*/
    double a;
    /**energy of dihedral bias*/
    double d;
    /**energy of spherical bias*/
    double s;
    /**energy of cubic bias*/
    double c;
    /**energy of threshold*/
    double thr;
    /**energy of bottom-threshold*/
    double thrB;
    /**energy of umbrella combination bias*/
    double u;

    //biases

    /**distance biases*/
    std::vector<config::biases::distance>  m_distances;
    /**angle biases*/
    std::vector<config::biases::angle>     m_angles;
    /**dihedral biases*/
    std::vector<config::biases::dihedral>  m_dihedrals;
    /**spherical biases*/
    std::vector<config::biases::spherical> m_spherical;
    /**cubic biases*/
    std::vector<config::biases::cubic>     m_cubic;
    /**threshold biases*/
    std::vector<config::biases::thresholdstr>  m_thresh;
    std::vector<config::biases::thresholdstr>  m_threshBottom;
    // umbrella biases
    /**distance biases for umbrella*/
    std::vector<config::coords::umbrellas::umbrella_dist> m_udist;
    /**angle biases for umbrella*/
    std::vector<config::coords::umbrellas::umbrella_angle> m_uangles;
    /**dihedral biases for umbrella*/
    std::vector<config::coords::umbrellas::umbrella_tor> m_utors;
    /**linear combinations biases for umbrella (can also be used as a "normal" bias)*/
    std::vector<config::coords::umbrellas::umbrella_comb> m_ucombs;

    // applying biases

    /**function to apply bias potential on dihedral*/
    double dih(Representation_3D const& xyz, Gradients_3D& g_xyz);
    /**function to apply bias potential on distance*/
    double dist(Representation_3D const& xyz, Gradients_3D& g_xyz);
    /**function to apply bias potential on angle*/
    double ang(Representation_3D const& xyz, Gradients_3D& g_xyz);
    /**function to apply spherical potential*/
    double spherical(Representation_3D const& xyz, Gradients_3D& g_xyz,
      Cartesian_Point const& center = Cartesian_Point());
    /**function to apply cubic potential*/
    double cubic(Representation_3D const& xyz, Gradients_3D& g_xyz,
      Cartesian_Point const& center = Cartesian_Point());
    /**function to apply a restraint on an umbrella dihedral and saves the values for the 'umbrella.txt' file
    @param xyz: coordinates of system
    @param g_xyz: cartesian gradients of system
    @param uout: vector with the real values for the restraint coordinate*/
    void umbrelladih(Representation_3D const& xyz, Gradients_3D& g_xyz, std::vector<double>& uout) const;
    /**function to apply a restraint on an umbrella angle and saves the values for the 'umbrella.txt' file
    @param xyz: coordinates of system
    @param g_xyz: cartesian gradients of system
    @param uout: vector with the real values for the restraint coordinate*/
    void umbrellaang(Representation_3D const& xyz, Gradients_3D& g_xyz, std::vector<double>& uout) const;
    /**function to apply a restraint on an umbrella distance and saves the values for the 'umbrella.txt' file
    @param xyz: coordinates of system
    @param g_xyz: cartesian gradients of system
    @param uout: vector with the real values for the restraint coordinate*/
    void umbrelladist(Representation_3D const& xyz, Gradients_3D& g_xyz, std::vector<double>& uout)  const;
    /**function to apply a restraint on an umbrella combination of distances and saves the values for the 'umbrella.txt' file
    @param xyz: coordinates of system
    @param g_xyz: cartesian gradients of system
    @param uout: vector with the real values for the restraint coordinate*/
    void umbrellacomb(Representation_3D const& xyz, Gradients_3D& g_xyz, std::vector<double>& uout);
    /**function to apply a umbrella combination biases on energy and gradients outside of umbrella sampling
    @param xyz: coordinates of system
    @param g_xyz: cartesian gradients of system
    returns the additional bias energy*/
    double umbrellacomb(Representation_3D const& xyz, Gradients_3D& g_xyz);
    /**function to apply threshold potential*/
    double thresh(Representation_3D const& xyz, Gradients_3D& g_xyz, Cartesian_Point maxPos);
    double thresh_bottom(Representation_3D const& xyz, Gradients_3D& g_xyz, Cartesian_Point minPos);

    /**function to apply spline on potential gradients (calculates prefactor and calls appy_spline)
      @param s: spline function
      @param xyz: coordinates of system
      @param g_xyz: cartesian gradients of system (are changed in function)*/
    void pmf_ic_spline(Spline const& s, Representation_3D const& xyz, Gradients_3D& g_xyz);
    /**applies spline gradient for different kinds of xi (only on one dimension)
    @param prefactor: derivative of spline by xi
    @param xyz: coordinates of system
    @param indices: indices that define xi
    @param g_xyz; cartesian gradients of system (are changed in this function)*/
    void apply_spline_1d(double prefactor, Representation_3D const& xyz, std::vector<std::size_t> const& indices, Gradients_3D& g_xyz);
    /**function to apply a spline on a torsion
    @param prefactor: derivative of spline by xi
    @param xyz: coordinates of system
    @param indices: indices that define xi
    @param g_xyz; cartesian gradients of system (are changed in this function)*/
    void apply_spline_on_torsion(double prefactor, Representation_3D const& xyz, std::vector<std::size_t> const& indices, Gradients_3D& g_xyz);
    /**function to apply a spline on an angle
    @param prefactor: derivative of spline by xi
    @param xyz: coordinates of system
    @param indices: indices that define xi
    @param g_xyz; cartesian gradients of system (are changed in this function)*/
    void apply_spline_on_angle(double prefactor, Representation_3D const& xyz, std::vector<std::size_t> const& indices, Gradients_3D& g_xyz);
    /**function to apply a spline on a distance
    @param prefactor: derivative of spline by xi
    @param xyz: coordinates of system
    @param indices: indices that define xi
    @param g_xyz; cartesian gradients of system (are changed in this function)*/
    void apply_spline_on_distance(double prefactor, Representation_3D const& xyz, std::vector<std::size_t> const& indices, Gradients_3D& g_xyz);
  };
}
