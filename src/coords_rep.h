#if !defined(COORDS_REP_H)
#define COORDS_REP_H

#include <iostream>
#include <vector>
#include <cstddef>
#include <utility>

#include "ls.h"
#include "scon_angle.h"
#include "scon_spherical.h"
#include "scon_vect.h"
#include "scon_utility.h"
#include "scon_matrix.h"

#include "ls.h"
#include "lbfgs.h"

#include "representation.h"

/**struct for a point charge*/
struct PointCharge
{
  /**position*/
	double x, y, z;
  /**charge*/
	double charge;

  /**function to set position*/
	void set_xyz(double m_x, double m_y, double m_z)
	{
		x = m_x;
		y = m_y;
		z = m_z;
	}
};

/**overloaded output operator for PointCharge*/
inline std::ostream& operator<< (std::ostream& stream, const PointCharge& c)
{
	stream << c.x << ", " << c.y << ", " << c.z << ", charge: " << c.charge;
	return stream;
}

namespace coords
{

  class Coordinates;

  template<class T> using Container = scon::vector < T >;

  // which floating point type do we use?
  using float_type = double;
  // size type vectors: 1d and 2d
  using size_1d = Container<std::size_t>;
  using size_2d = Container< Container<std::size_t> >;

  // Angles
  using angle_type = scon::ang<float_type>;

  using r3 = scon::c3<float_type>;
  using s3 = scon::sphericals<float_type>;


  /* #############################################

  Coordinate representations
  (intern/cartesian/dihedral)

  ########  ######## ########  ########
  ##     ## ##       ##     ## ##     ##
  ##     ## ##       ##     ## ##     ##
  ########  ######   ########  ########
  ##   ##   ##       ##        ##   ##
  ##    ##  ##       ##        ##    ##  ###
  ##     ## ######## ##        ##     ## ###

  ############################################# */

  // 1d, 3d and N*3d points

  template<class P, class G>
  using State = optimization::State < P, G >;

  using Cartesian = State < Container<r3>, Container<r3> >;
  using Internal = State < Container<s3>, Container<r3> >;
  using Main_Dihedral = State < Container<angle_type>, Container<float_type> >;


  /* #############################################

  Coordinate fixations
  (intern/cartesian/dihedral)

  ######## #### ##     ##
  ##        ##   ##   ##
  ##        ##    ## ##
  ######    ##     ###
  ##        ##    ## ##
  ##        ##   ##   ##
  ##       #### ##     ##

  ############################################# */

  using fix3 = scon::c3<bool>;



  /* #############################################

  ######   #######  ##     ## ########  ####
  ##    ## ##     ## ###   ### ##     ##  ##
  ##       ##     ## #### #### ##     ##  ##
  ##       ##     ## ## ### ## ########   ##
  ##       ##     ## ##     ## ##     ##  ##
  ##    ## ##     ## ##     ## ##     ##  ##
  ######   #######  ##     ## ########  ####

  ############################################# */

  // Cartesian
  using Cartesian_Point = r3;
  using cartesian_type = r3;
  using cartesian_gradient_type = r3;
  using Representation_3D = Container<Cartesian_Point>;
  using Gradients_3D = Container<cartesian_gradient_type>;
  using Fixations_3D = std::vector<fix3>;
  // Internal
  using internal_type = scon::sphericals<float_type>;
  using internal_gradient_type = scon::c3<float_type>;
  using Representation_Internal = Container<internal_type>;
  using Gradients_Internal = Container<internal_gradient_type>;
  using Fixations_Internal = std::vector<fix3>;
  // Mains
  using main_type = internal_type::angle_type;
  using main_gradient_type = float_type;
  using Representation_Main = Container< main_type >;
  using Gradients_Main = Container< main_gradient_type >;
  using Fixations_Main = std::vector< bool >;

  using Representation_1D = Container<float_type>;

  // subsystem interaction energy and gradients
  struct sub_ia
  {
    coords::Representation_3D grad;
    double energy;
  };
  // interaction matrix
  using sub_ia_matrix_t = scon::matrix<sub_ia, true>;

  // template<class 

  template<class Rep3D, class RepInt, class RepMain>
  struct RepType
  {

    using size_type = std::size_t;

    Rep3D cartesian;
    RepInt intern;
    RepMain main;

    RepType()
      : cartesian(), intern(), main()
    { }

    RepType(size_type const n, size_type const m = 0)
      : cartesian(n), intern(n), main(m)
    { }

    RepType(Representation_3D const & xyz)
      : cartesian(xyz), intern(xyz.size()), main()
    { }

    RepType(Representation_3D && xyz)
      : cartesian(std::forward<Representation_3D>(xyz)),
      intern(cartesian.size()), main()
    { }

    RepType(Rep3D const & xyz, RepInt const & inter, RepMain const & maindih)
      : cartesian(xyz), intern(inter), main(maindih)
    { }

    RepType(Rep3D && xyz, RepInt && inter, RepMain && maindih)
      : cartesian(std::forward<Representation_3D>(xyz)),
      intern(std::forward<Representation_3D>(inter)),
      main(maindih)
    { }

    size_type size() const { return cartesian.size(); }

    void resize(size_type const n, size_type const mains = 0)
    {
      cartesian.resize(n);
      intern.resize(n);
      main.resize(mains);
    }

    void clear()
    {
      //scon::clear(cartesian, intern, main);
      cartesian.clear();
      intern.clear();
      main.clear();
    }

    void swap(RepType &r)
    {
      cartesian.swap(r.cartesian);
      intern.swap(r.intern);
      main.swap(r.main);
    }
    bool empty(void) {
      return cartesian.empty();
    }
    Cartesian_Point cartesian_vector(size_type a, size_type b) const
    {
      return (cartesian[a] - cartesian[b]);
    }

    friend inline std::ostream& operator<< (std::ostream &stream, RepType<Rep3D, RepInt, RepMain> const & rep)
    {
      stream << "C:";
      for (auto const & e : rep.cartesian) stream << e << '\n';
      stream << "I:";
      for (auto const & e : rep.intern) stream << e << '\n';
      stream << "M:";
      for (auto const & e : rep.main) stream << e << '\n';
      return stream;
    }

  };

  //template<class 

  //using 

  using Representation_Dual = RepType<Representation_3D, Representation_Internal, Representation_Main>;
  using Gradients_Dual = RepType<Gradients_3D, Gradients_Internal, Gradients_Main>;

  template<class Rep3D, class RepInt, class RepMain>
  inline void swap(RepType<Rep3D, RepInt, RepMain> &a, RepType<Rep3D, RepInt, RepMain> &b) { a.swap(b); }

  using Ensemble_3d = Container<Representation_3D>;
  using Ensemble_Dual = Container<Representation_Dual>;


  /* #############################################################################

  class PES_Point

  ########  ########  ######          ########   #######  #### ##    ## ########
  ##     ## ##       ##    ##         ##     ## ##     ##  ##  ###   ##    ##
  ##     ## ##       ##               ##     ## ##     ##  ##  ####  ##    ##
  ########  ######    ######          ########  ##     ##  ##  ## ## ##    ##
  ##        ##             ##         ##        ##     ##  ##  ##  ####    ##
  ##        ##       ##    ##         ##        ##     ##  ##  ##   ###    ##
  ##        ########  ######  ####### ##         #######  #### ##    ##    ##

  ############################################################################# */

  struct PES_Point
  {
    using size_type = Representation_Dual::size_type;
    Representation_Dual structure;
    Gradients_Dual gradient;
    /**hessian matrix: outer vector is line, inner vector is column
    the order of the values is like this: 
              x(atom1)  y(atom1)  z(atom1)  x(atom2) ...
    x(atom1) .............................................
    y(atom1) .............................................
    z(atom1) .............................................
    x(atom2) .............................................*/
    std::vector<std::vector<double>> hessian;
    sub_ia_matrix_t ia_matrix;
    float_type energy;
    bool integrity;
    PES_Point(void)
      : structure(), gradient(), energy(), integrity(false)
    { }
    PES_Point(size_type n)
      : structure(n, 0), gradient(n, 0), energy(), integrity(false)
    { }
    PES_Point(Representation_3D xyz, bool const structure_integrity = true) :
      structure(xyz, Representation_Internal(xyz.size()), Representation_Main()),
      gradient(Gradients_3D(xyz.size()), Gradients_Internal(xyz.size()), Gradients_Main()),
      energy(0.0), integrity(structure_integrity)
    { }
    PES_Point(double E,
      Representation_3D xyz, Gradients_3D xyz_force,
      Representation_Internal inter, Gradients_Internal inter_force,
      Representation_Main maindih, Gradients_Main maindih_force,
      bool const structure_integrity = true) :
      structure(xyz, inter, maindih), gradient(xyz_force, inter_force, maindih_force),
      energy(E), integrity(structure_integrity)
    { }
    void swap(PES_Point &rhs)
    {
      using std::swap;
      structure.swap(rhs.structure);
      gradient.swap(rhs.gradient);
      ia_matrix.swap(rhs.ia_matrix);
      swap(energy, rhs.energy);
      swap(integrity, rhs.integrity);
    }
    //bool equal_compare(PES_Point const &) const;
    bool empty(void)
    {
      return structure.empty();
    }
    size_type size(void) const { return structure.size(); }
    void resize(size_type atoms, size_type const mains = 0)
    {
      structure.resize(atoms, mains);
      gradient.resize(atoms, mains);
    }
    bool same_size(PES_Point const &rhs) const
    {
      return structure.cartesian.size() == rhs.structure.cartesian.size()
        && gradient.cartesian.size() == rhs.gradient.cartesian.size()
        && structure.intern.size() == rhs.structure.intern.size()
        && gradient.intern.size() == rhs.gradient.intern.size()
        && structure.main.size() == rhs.structure.main.size()
        && gradient.main.size() == rhs.gradient.main.size();
    }
    void clear(void)
    {
      structure.clear();
      gradient.clear();
      energy = 0.0;
    }
    friend std::ostream& operator<< (std::ostream &stream, PES_Point const & rep)
    {
      stream << "S:" << rep.structure;
      stream << "G:" << rep.gradient;
      stream << "E:" << rep.energy;
      stream << "I:" << rep.integrity;
      return stream;
    }
  };

  inline void swap(PES_Point &a, PES_Point &b) { a.swap(b); }
  //inline bool operator== (const PES_Point& lhs, const PES_Point& rhs) { return lhs.equal_compare(rhs); }
  //inline bool operator!= (const PES_Point& lhs, const PES_Point& rhs) { return !operator==(lhs, rhs); }
  inline bool operator<  (const PES_Point& lhs, const PES_Point& rhs) { return lhs.energy < rhs.energy; }
  inline bool operator>(const PES_Point& lhs, const PES_Point& rhs) { return  operator< (rhs, lhs); }
  inline bool operator<= (const PES_Point& lhs, const PES_Point& rhs) { return !operator> (lhs, rhs); }
  inline bool operator>= (const PES_Point& lhs, const PES_Point& rhs) { return !operator< (lhs, rhs); }

  // ...
  using Ensemble_PES = Container<PES_Point>;

}

#endif
