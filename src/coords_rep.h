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

namespace coords
{

  class Coordinates;

  template<class T> using Container = scon::vector < T >;

  // which floating point type do we use?
  typedef double float_type;
  // size type vectors: 1d and 2d
  typedef Container<std::size_t> size_1d;
  typedef Container< Container<std::size_t> > size_2d;

  // Angles
  typedef scon::ang<float_type> angle_type;

  typedef scon::c3<float_type> r3;
  typedef scon::sphericals<float_type> s3;



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

  typedef scon::c3<bool> fix3;



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
  typedef r3                                      Cartesian_Point;
  typedef r3                                      cartesian_type;
  typedef r3                                      cartesian_gradient_type;
  typedef Container<Cartesian_Point>              Representation_3D;
  typedef Container<cartesian_gradient_type>      Gradients_3D;
  typedef std::vector<fix3>                       Fixations_3D;
  // Internal
  typedef scon::sphericals<float_type>            internal_type;
  typedef scon::c3<float_type>                    internal_gradient_type;
  typedef Container<internal_type>                Representation_Internal;
  typedef Container<internal_gradient_type>       Gradients_Internal;
  typedef std::vector<fix3>                       Fixations_Internal;
  // Mains
  typedef internal_type::angle_type               main_type;
  typedef float_type                              main_gradient_type;
  typedef Container< main_type >                  Representation_Main;
  typedef Container< main_gradient_type >         Gradients_Main;
  typedef std::vector< bool >                     Fixations_Main;

  typedef Container<float_type>                   Representation_1D;

  // subsystem interaction energy and gradients
  struct sub_ia
  {
    coords::Representation_3D grad;
    double energy;
  };
  // interaction matrix
  typedef scon::matrix<sub_ia, true> sub_ia_matrix_t;

  // template<class 

  template<class Rep3D, class RepInt, class RepMain>
  struct RepType
  {

    typedef std::size_t size_type;

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

  typedef RepType<Representation_3D, Representation_Internal, Representation_Main> Representation_Dual;
  typedef RepType<Gradients_3D, Gradients_Internal, Gradients_Main> Gradients_Dual;

  template<class Rep3D, class RepInt, class RepMain>
  inline void swap(RepType<Rep3D, RepInt, RepMain> &a, RepType<Rep3D, RepInt, RepMain> &b) { a.swap(b); }

  typedef Container<Representation_3D> Ensemble_3d;
  typedef Container<Representation_Dual> Ensemble_Dual;


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
    typedef Representation_Dual::size_type size_type;
    Representation_Dual structure;
    Gradients_Dual gradient;
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
  typedef Container<PES_Point> Ensemble_PES;

}

#endif
