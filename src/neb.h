#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <cstdlib>
#include "atomic.h"
#include "configuration.h"
#include "coords.h"
#include "coords_io.h"
#include "scon_vect.h"
#include "tinker_refine.h"
#include "optimization_global.h"
#include <math.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "interpolation.h"
#include "nr3.h"

using namespace scon;
class neb
{
public:
  bool ts, reversed, ClimbingImage;

  neb(void);

  coords::Coordinates *cPtr;
  coords::Atoms atoms_neb;
  neb(coords::Coordinates *c);
  optimization::global::optimizers::monteCarlo* mc;

  std::vector <coords::Representation_3D> imagi, imagederiv, tau, image_ini;
  coords::Representation_3D ts_pathstruc;
  std::vector <scon::c3 <float> > tau_int;
  std::vector <double> energies, ts_energies, min_energies,energies_NEB;
  double grad_v, grad_v_temp;
  size_t N, num_images, global_imagex;

  // IDPP start
  coords::Representation_3D start_structure, final_structure;

  void idpp_prep();

  template<typename ContainerIt, typename Result>
  void concatenate(ContainerIt, ContainerIt, Result);

  template<typename T>
  using dist_T = std::vector<std::vector<T>>;

  template<typename T>
  dist_T<T> euclid_dist(coords::Representation_3D const&);

  std::vector<coords::Representation_3D> bond_dir(coords::Representation_3D const&);

  template<typename T, typename Euclid>
  std::vector<T> interpol_dist(Euclid&, Euclid&, size_t const);

  template<typename T, typename Euclid>
  coords::Representation_3D idpp_gradients(std::vector<coords::Representation_3D> const&,
	  Euclid&, Euclid&, Euclid&, size_t const);

private:
	std::vector<coords::Representation_3D> bond_st;
	coords::Representation_3D Fidpp;
	std::vector<std::vector<double>> eu_st, eu_fi;
	// IDPP end

public:

  void run(ptrdiff_t &count);
  void preprocess(ptrdiff_t &count);
  void preprocess(ptrdiff_t &image, ptrdiff_t &count, const coords::Representation_3D &start, const coords::Representation_3D &fi, const std::vector <double> &ts_energy, const std::vector <double> &min_energy, bool reverse, const coords::Representation_3D &ts_path);
  void preprocess(ptrdiff_t &file, ptrdiff_t &image, ptrdiff_t &count, const coords::Representation_3D &start, const coords::Representation_3D &fi, bool reverse);
  void initial(void);
  void final(void);
  void initial(const coords::Representation_3D &start);
  void final(const coords::Representation_3D &fi);
  void create();
  void get_energies(void);
  void calc_tau(void);
  void opt_io(ptrdiff_t &count);
  void print(std::string const&, std::vector <coords::Representation_3D > &, ptrdiff_t &);
  void print_rev(std::string const &, std::vector <coords::Representation_3D > &, ptrdiff_t &);
  void printmono(std::string const &, coords::Representation_3D &print, ptrdiff_t &);
  double g_new();
  double g_new_maxflux();
  double g_int(std::vector <scon::c3 <float> >  tx);
  void calc_shift(void);
  double dot_uneq(coords::Representation_3D const &a, coords::Representation_3D const &b);

  //Julians implementation
  void run(ptrdiff_t &count, std::vector<size_t>& image_remember, std::vector<std::vector<size_t> >& atoms_remember);

  void internal_execute(std::vector <coords::Representation_3D> &input, std::vector<std::vector<size_t> >& atoms_remember);
  void biggest(const std::vector<std::vector<std::vector<double> > >& dist, const std::vector<std::vector<std::vector<double> > >& anglevec, const std::vector<std::vector<std::vector<double> > >& dihedralvec, std::vector<std::vector<std::vector<double> > >& comp,
    std::vector<std::vector<std::vector<size_t> > >& which_img, std::vector<std::vector<size_t> >& i_remember,
    std::vector<std::vector<size_t> >& j_remember, std::vector <std::vector <size_t> > & i_remember_all, std::vector<std::vector<size_t> > & j_remember_all);
  //void very_biggest(std::vector<std::vector<std::vector<double> > >& comp, double* biggest_dist, double* biggest_angle, double* biggest_dihedral, std::vector<size_t>& i_remember, std::vector<size_t>& j_remember);
  void double_or_not(std::vector<std::vector<std::vector<size_t> > >& involved_bonds, const std::vector<size_t>& pivot, const size_t& a, const size_t& b, std::vector<std::vector<double> >& vec);
  void Sort(const size_t& imgs, std::vector<std::string>& name, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y, std::vector<std::vector<double> >& z, std::vector<std::vector<size_t> >& bonds);
  void bond_Sort(std::vector<size_t>& bonds, const size_t& start, const size_t& length);
  void bond_Sort(size_t* a, size_t* b, size_t* c, size_t* d);
  void bond_Sort(size_t* a, size_t* b, size_t* c);
  void bond_Sort(size_t* a, size_t* b);
  void hetero_Sort(const size_t& imgs, const size_t& start, const size_t& length, std::vector<std::string>& name, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y, std::vector<std::vector<double> >& z, std::vector<std::vector<size_t> >& bonds);
  bool dihedral_or_not(const int& first_atom, const std::vector<int>& dist, const std::vector<int>& ang, const int& last_atom);
  void get_values(std::vector<double>& x_val, std::vector<double>& y_val, std::vector<double>& z_val, const size_t& atom_iter, std::vector<std::vector<size_t> >& which_bonds, std::vector<std::vector<double> >* dist, std::vector<std::vector<double> >* anglevec, std::vector<std::vector<double> >* dihedralvec, std::vector<std::vector<std::vector<std::vector<size_t> > > >* involved_bonds, std::vector<std::vector<size_t> >& bonds);
  void get_values(std::vector<double>& x_val, std::vector<double>& y_val, std::vector<double>& z_val, std::vector<std::vector<double> >* dist, std::vector<std::vector<double> >* anglevec, std::vector<std::vector<double> >* dihedralvec, std::vector<std::vector<std::vector<std::vector<size_t> > > >& involved_bonds, std::vector<std::vector<size_t> >& which_bonds);
  size_t factorial(const size_t& fac);
  double inline distance(const double& x1, const double& x2, const double& y1, const double& y2, const double& z1, const double& z2);
  double inline angle(const double& x1, const double& x2, const double& x3, const double& y1, const double& y2, const double& y3, const double& z1, const double& z2, const double& z3);
  double inline dihedral(const double& x1, const double& x2, const double& x3, const double& x4, const double& y1, const double& y2, const double& y3, const double& y4, const double& z1, const double& z2, const double& z3, const double& z4);
  double inline dihedral_same_atom(const double& x1, const double& x2, const double& x3, const double& x4, const double& y1, const double& y2, const double& y3, const double& y4, const double& z1, const double& z2, const double& z3, const double& z4);
  void no_torsion(const size_t& image_remember, const std::vector<int>& atoms_remember);
  void execute_fix(const std::vector<size_t>& atoms_remember);
  void execute_defix(const std::vector<size_t>& atoms_remember);
  void defix_all(void);
  void opt_internals(ptrdiff_t &count, const std::vector<std::vector<size_t> >& atoms_remember);

  double lbfgs();
  double lbfgs_int(std::vector <scon::c3 <float> > t);
  double lbfgs_maxflux();

  struct GradCallBack
  {
    neb * p;
    GradCallBack(neb & nebo) : p(&nebo) {}
    float operator() (scon::vector< scon::c3<float> > const &x,
      scon::vector< scon::c3<float> > &g, size_t const, bool & go_on)
    {
      using ot = scon::vector< scon::c3<float> >;

      using ct = coords::Representation_3D;
      //p->cPtr->set_xyz(std::move(scon::explicit_transform<ct>(x)));
	  p->images_initial = scon::explicit_transform<ct>(x);
      float E = static_cast<float>(p->g_new());
      go_on = p->cPtr->integrity();
      g = scon::explicit_transform<ot>(p->grad_tot);
      return E;
    }
  };

  struct GradCallBackMaxFlux
  {
	  neb * p;
	  GradCallBackMaxFlux(neb & nebo) : p(&nebo) {}
	  float operator() (scon::vector< scon::c3<float> > const &x,
		  scon::vector< scon::c3<float> > &g, size_t const, bool & go_on)
	  {
		  using ot = scon::vector< scon::c3<float> >;

		  using ct = coords::Representation_3D;
		  //p->cPtr->set_xyz(std::move(scon::explicit_transform<ct>(x)));
		  p->images_initial = scon::explicit_transform<ct>(x);
		  float E = static_cast<float>(p->g_new_maxflux());
		  go_on = p->cPtr->integrity();
		  g = scon::explicit_transform<ot>(p->grad_tot);
		  return E;
	  }
  };

  struct GradCallBack_int
  {
    neb * p;
    GradCallBack_int(neb & nebo) : p(&nebo) {}
    float operator() (scon::vector< scon::c3<float> > const &x,
      scon::vector< scon::c3<float> > &g, size_t const, bool & go_on)
    {
      using ot = scon::vector< scon::c3<float> >;

      using ct = coords::Representation_3D;
      p->cPtr->set_xyz(std::move(scon::explicit_transform<ct>(x)));
      float E = static_cast<float>(p->g_int(p->tau_int));
      go_on = p->cPtr->integrity();
      g = scon::explicit_transform<ot>(p->cPtr->g_xyz());
      return E;
    }
  };
  

private:
  double springconstant;
  double tauderiv;
  double EnergyPml, EnergyPpl;
  const double _KT_;

  ptrdiff_t natoms;
  ptrdiff_t nvar;
  ptrdiff_t CIMaximum;

  coords::Representation_3D images_initial, grad_tot;
  coords::Representation_3D images, tempimage_final, tempimage_ini, Rm1, Rp1, Fvertical, Fpar;

};

enum { DIST = 0, ANGLE, DIHEDRAL }; //Juli