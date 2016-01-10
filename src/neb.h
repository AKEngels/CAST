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
#include "tinker_refine.h"
#include <math.h>
#include <algorithm>
#include<cmath>
#include<fstream>
#include<iomanip>
#include "nr3.h"

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

	static ::tinker::parameter::parameters tp;
	::tinker::parameter::parameters cparams;
	::tinker::refine::refined refined;
		
		
	coords::Representation_3D images_initial;
	coords::Representation_3D images,tempimage_final,tempimage_ini,Rm1,Rp1,Fvertical;
	std::vector <coords::Representation_3D> imagi,imagederiv,tau,image_ini,Fpar,imagi_test;
    coords::Representation_3D ts_pathstruc;
    std::vector <coords::Representation_3D > NEB_force;
	std::vector <scon::c3 <float> > tau_int;
	std::vector <double> energies,ts_energies,min_energies,energies_NEB;
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
	double g_int(std::vector <scon::c3 <float> >  tx);
	void calc_shift(void);
	
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
	
	double lbfgs(ptrdiff_t imagex);
	double lbfgs_int(std::vector <scon::c3 <float> > t);

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
	struct Base_interp


	{
		int n, mm, jsav, cor, dj;
		const Doub *xx, *yy;
		Base_interp(VecDoub_I &x, const Doub *y, int m)
			: n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
			dj = MIN(1, (int)pow((Doub)n, 0.25));
		}




		Doub interp(Doub x) {
			int j1o = cor ? hunt(x) : locate(x);
			return rawinterp(j1o, x);
		}

		int locate(const Doub x)

		{
			int ju, jm, j1;
			if (n < 2 || mm < 2 || mm > n) throw ("locate size error");
			bool ascnd = (xx[n - 1] >= xx[0]);
			j1 = 0;
			ju = n - 1;
			while (ju - j1 > 1) {
				jm = (ju + j1) >> 1;
				if ((x >= xx[jm]) == ascnd)
					j1 = jm;
				else
					ju = jm;
			}
			cor = abs(j1 - jsav) > dj ? 0 : 1;
			jsav = j1;
			return MAX(0, MIN(n - mm, j1 - ((mm - 2) >> 1)));
		}

		int hunt(const Doub x)
		{
			int j1 = jsav, jm, ju, inc = 1;
			if (n < 2 || mm<2 || mm > n) throw ("hunt size error");
			bool ascnd = (xx[n - 1] >= xx[0]);
			if (j1< 0 || j1>n - 1) {
				j1 = 0;
				ju = n - 1;
			}
			else {
				if ((x >= xx[j1]) == ascnd) {
					for (;;) {
						ju = j1 + inc;
						if (ju >= n - 1) { ju = n - 1; break; }
						else if ((x < xx[ju]) == ascnd) break;
						else {
							j1 = ju;
							inc += inc;
						}
					}
				}
				else {
					ju = j1;
					for (;;) {
						j1 = j1 - inc;
						if (j1 <= 0) { j1 = 0; break; }
						else if ((x >= xx[j1]) == ascnd) break;
						else {
							ju = j1;
							inc += inc;
						}
					}
				}
			}
			while (ju - j1 > 1) {
				jm = (ju + j1) >> 1;
				if ((x >= xx[jm]) == ascnd)
					j1 = jm;
				else
					ju = jm;
			}
			cor = abs(j1 - jsav) < dj ? 0 : 1;
			jsav = j1;
			return MAX(0, MIN(n - mm, j1 - ((mm - 2) >> 1)));
		}

		Doub virtual rawinterp(int j1o, Doub x) = 0;


	};

	struct Poly_interp : Base_interp
	{

		Doub dy;
		Poly_interp(VecDoub_I &xv, VecDoub_I &yv, int m)
			: Base_interp(xv, &yv[0], m), dy(0.) {}




		Doub rawinterp(int j1, Doub x)


		{






			int i, m, ns = 0;
			Doub y, den, dif, dift, ho, hp, w;
			const Doub *xa = &xx[j1], *ya = &yy[j1];
			VecDoub c(mm), d(mm);
			dif = abs(x - xa[0]);
			for (i = 0; i < mm; i++) {
				if ((dift = abs(x - xa[i])) < dif) {
					ns = i;
					dif = dift;
				}
				c[i] = ya[i];
				d[i] = ya[i];
			}
			y = ya[ns--];
			for (m = 1; m < mm; m++) {
				for (i = 0; i < mm - m; i++) {
					ho = xa[i] - x;
					hp = xa[i + m] - x;
					w = c[i + 1] - d[i];
					if ((den = ho - hp) == 0.0) throw ("Poly_interp error");
					den = w / den;
					d[i] = hp*den;
					c[i] = ho*den;
				}
				y += (dy = (2 * (ns + 1) < (mm - m) ? c[ns + 1] : d[ns--]));
			}
			return y;

		}

	};


	struct Spline_interp : Base_interp
	{
		VecDoub y2;

		Spline_interp(VecDoub_I &xv, VecDoub_I &yv, Doub yp1 = 1.e99, Doub ypn = 1.e99)
			: Base_interp(xv, &yv[0], 2), y2(xv.size())
		{
			sety2(&xv[0], &yv[0], yp1, ypn);
		}

		Spline_interp(VecDoub_I &xv, const Doub *yv, Doub yp1 = 1.e99, Doub ypn = 1.e99)
			: Base_interp(xv, yv, 2), y2(xv.size())
		{
			sety2(&xv[0], yv, yp1, ypn);
		}







		void sety2(const Doub *xv, const Doub *yv, Doub yp1, Doub ypn)
		{
			Int i, k;
			Doub p, qn, sig, un;
			Int n = y2.size();
			VecDoub u(n - 1);
			if (yp1 > 0.99e99)
				y2[0] = u[0] = 0.0;
			else {
				y2[0] = -0.5;
				u[0] = (3.0 / (xv[1] - xv[0]))*((yv[1] - yv[0]) / (xv[1] - xv[0]) - yp1);
			}
			for (i = 1; i < n - 1; i++) {
				sig = (xv[i] - xv[i - 1]) / (xv[i + 1] - xv[i - 1]);
				p = sig*y2[i - 1] + 2.0;
				y2[i] = (sig - 1.0) / p;
				u[i] = (yv[i + 1] - yv[i]) / (xv[i + 1] - xv[i]) - (yv[i] - yv[i - 1]) / (xv[i] - xv[i - 1]);
				u[i] = (6.0*u[i] / (xv[i + 1] - xv[i - 1]) - sig*u[i - 1]) / p;
			}
			if (ypn > 0.99e99)
				qn = un = 0.0;
			else {
				qn = 0.5;
				un = (3.0 / (xv[n - 1] - xv[n - 2]))*(ypn - (yv[n - 1] - yv[n - 2]) / (xv[n - 1] - xv[n - 2]));
			}
			y2[n - 1] = (un - qn*u[n - 2]) / (qn*y2[n - 2] + 1.0);
			for (k = n - 2; k >= 0; k--)
				y2[k] = y2[k] * y2[k + 1] + u[k];
		}
		Doub rawinterp(Int jl, Doub x)
		{
			Int klo = jl, khi = jl + 1;
			Doub y, h, b, a;
			h = xx[khi] - xx[klo];
			if (h == 0.0) throw("Bad input to routine splint");
			a = (xx[khi] - x) / h;
			b = (x - xx[klo]) / h;
			y = a*yy[klo] + b*yy[khi] + ((a*a*a - a)*y2[klo]
				+ (b*b*b - b)*y2[khi])*(h*h) / 6.0;
			return y;
		}

	};








private:
	double springconstant;
	double tauderiv;
	double EnergyPml, EnergyPpl;

	ptrdiff_t natoms;
	ptrdiff_t nvar;
	ptrdiff_t CIMaximum;


};

enum{ DIST = 0, ANGLE, DIHEDRAL }; //Juli