#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <cstdlib>
#include "error.h"
#include "atomic.h"
#include "configuration.h"
#include "coords.h"
#include "coords_io.h"
#include "scon_vect.h"

using namespace scon;
class neb
{
public:
	bool controlrun,ts,reversed,ClimbingImage;

	neb(void);

	coords::Coordinates *cPtr;
	//coords::Coordinates const * coords;
	coords::Atoms atoms_neb;
	coords::Coordinates coords_neb;
	neb(coords::Coordinates *c) ;
	optimization::global::optimizers::monteCarlo* mc;
	neb(const std::size_t);

	~neb(void);

	
		
		
	coords::Representation_3D images_initial;
  coords::Representation_3D images,tempimage_final,tempimage_ini,Rm1,Rp1,Fvertical;
	std::vector <coords::Representation_3D> imagi,imagederiv,tau,image_ini,Fpar,imagi_test;
    coords::Representation_3D ts_pathstruc;
    std::vector <coords::Representation_3D > NEB_force;
	std::vector <double> energies,ts_energies,min_energies;
	double grad_v, grad_v_temp;
    size_t N,num_images,global_imagex;
		

	void run (ptrdiff_t &count);
	void preprocess(ptrdiff_t &count, ptrdiff_t const &image);
    void preprocess(ptrdiff_t &image, ptrdiff_t &count, const coords::Representation_3D &start, const coords::Representation_3D &fi, const std::vector <double> &ts_energy, const std::vector <double> &min_energy, bool reverse, const coords::Representation_3D &ts_path);
	void preprocess(ptrdiff_t &file, ptrdiff_t &image,ptrdiff_t &count, const coords::Representation_3D &start, const coords::Representation_3D &fi, bool reverse);
	void initial(void);
	void final(void);
	void initial(const coords::Representation_3D &start);
	void final(const coords::Representation_3D &fi);
	void create(ptrdiff_t &count);
	void get_energies_grad(void);
	void calc_tau(void);
	void opt_mep(ptrdiff_t &count);
	void print (std::string const&,std::vector <coords::Representation_3D > &,ptrdiff_t &);
	void print_rev (std::string const &,std::vector <coords::Representation_3D > &,ptrdiff_t &);
	void printmono (std::string const &,coords::Representation_3D &print,ptrdiff_t &);
	double g_new(ptrdiff_t im);
	void calc_shift(void);
	
	//Julians implementation
	void run(ptrdiff_t &count, std::vector<size_t>& image_remember, std::vector<std::vector<int> >& atoms_remember);

	void internal_execute(std::vector <coords::Representation_3D> &input, std::vector<std::vector<int> >& atoms_remember);
	void biggest(const std::vector<std::vector<std::vector<double> > >& dist, const std::vector<std::vector<std::vector<double> > >& anglevec, const std::vector<std::vector<std::vector<double> > >& dihedralvec, std::vector<std::vector<std::vector<double> > >& comp,
		std::vector<std::vector<std::vector<int> > >& which_img, std::vector<std::vector<size_t> >& i_remember,
		std::vector<std::vector<size_t> >& j_remember, std::vector <std::vector <size_t> > & i_remember_all, std::vector<std::vector<size_t> > & j_remember_all);
	void very_biggest(std::vector<std::vector<std::vector<double> > >& comp, double* biggest_dist, double* biggest_angle, double* biggest_dihedral, std::vector<size_t>& i_remember, std::vector<size_t>& j_remember);
	void double_or_not(std::vector<std::vector<std::vector<int> > >& involved_bonds, const std::vector<int>& pivot, const size_t& a, const size_t& b, std::vector<std::vector<double> >& vec);
	void Sort(const size_t& imgs, std::vector<std::string>& name, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y, std::vector<std::vector<double> >& z, std::vector<std::vector<int> >& bonds);
	void bond_Sort(std::vector<int>& bonds, const size_t& start, const size_t& length);
	void bond_Sort(int* a, int* b, int* c, int* d);
	void bond_Sort(int* a, int* b, int* c);
	void bond_Sort(int* a, int* b);
	void hetero_Sort(const size_t& imgs, const size_t& start, const size_t& length, std::vector<std::string>& name, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y, std::vector<std::vector<double> >& z, std::vector<std::vector<int> >& bonds);
	bool dihedral_or_not(const int& first_atom, const std::vector<int>& dist, const std::vector<int>& ang, const int& last_atom);
	void get_values(std::vector<double>& x_val, std::vector<double>& y_val, std::vector<double>& z_val, const size_t& atom_iter, std::vector<std::vector<int> >& which_bonds, std::vector<std::vector<double> >* dist, std::vector<std::vector<double> >* anglevec, std::vector<std::vector<double> >* dihedralvec, std::vector<std::vector<std::vector<std::vector<int> > > >* involved_bonds, std::vector<std::vector<int> >& bonds);
	void get_values(std::vector<double>& x_val, std::vector<double>& y_val, std::vector<double>& z_val, std::vector<std::vector<double> >* dist, std::vector<std::vector<double> >* anglevec, std::vector<std::vector<double> >* dihedralvec, std::vector<std::vector<std::vector<std::vector<int> > > >& involved_bonds, std::vector<std::vector<int> >& which_bonds);
	size_t factorial(const size_t& fac);
	double inline distance(const double& x1, const double& x2, const double& y1, const double& y2, const double& z1, const double& z2);
	double inline angle(const double& x1, const double& x2, const double& x3, const double& y1, const double& y2, const double& y3, const double& z1, const double& z2, const double& z3);
	double inline dihedral(const double& x1, const double& x2, const double& x3, const double& x4, const double& y1, const double& y2, const double& y3, const double& y4, const double& z1, const double& z2, const double& z3, const double& z4);
	double inline dihedral_same_atom(const double& x1, const double& x2, const double& x3, const double& x4, const double& y1, const double& y2, const double& y3, const double& y4, const double& z1, const double& z2, const double& z3, const double& z4);
	void no_torsion(const size_t& image_remember, const std::vector<int>& atoms_remember);
	void execute_fix(const std::vector<int>& atoms_remember);
	void execute_defix(const std::vector<int>& atoms_remember);
	void defix_all(void);
	void opt_internals(ptrdiff_t &count, const std::vector<std::vector<int> >& atoms_remember);

	
	double lbfgs(ptrdiff_t imagex);

	struct GradCallBack
	{
    neb * p;
    GradCallBack (neb & nebo) : p(&nebo) {}
    float operator() (scon::vector< scon::c3<float> > const &x, 
      scon::vector< scon::c3<float> > &g, size_t const, bool & go_on)
	  {
      using ot = scon::vector< scon::c3<float> >;
      using ct = coords::Representation_3D;
		  p->cPtr->set_xyz(std::move(scon::explicit_transform<ct>(x)));
		  float E = static_cast<float>(p->g_new(p->global_imagex));
		  go_on = p->cPtr->integrity();	
		  g = scon::explicit_transform<ot>(p->cPtr->g_xyz());
      return E;
	  }
	};



private:
	double springconstant;
	double tauderiv;
	double EnergyPml,EnergyPpl;
		 
	std::size_t natoms;
  std::size_t nvar;
  std::size_t CIMaximum;


};

enum{ DIST = 0, ANGLE, DIHEDRAL }; //Juli