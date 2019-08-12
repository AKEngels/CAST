/**
* PATHOPT class - Global reaction pathway determination
* via constraint minimization on "n-1" dimensional hyperplanes
*
*  Version 1.1
*
*  @ Daniel Bellinger
*
*/
#pragma once 

#include <vector>
#include "coords.h"
#include <ostream>
#include <memory>
#include <string>
#include "Scon/scon_vect.h"
#include "coords_rep.h"
#include "neb.h"
#include <cstdlib>
#include <iomanip>



class pathx
{
private:
	double randvec[3];
	double maxvar, mcstepsize;
	ptrdiff_t mciteration;
	double STARTENERGY, ENDENERGY, MCEN;
	std::vector<double>  MCcoordglob, MCcoordlast, MCcoordin;
	const double _KT_;
	std::vector <std::vector <double> > global_path_minima_energy;
	std::vector < std::vector <coords::Representation_3D> >  global_path_minima, global_path_minima_temp;

	struct perp_point
	{
		coords::PES_Point pes;
		std::vector<coords::Representation_Main> main_direction;
		std::vector<coords::Representation_Internal> intern_direction;
		std::vector<coords::Representation_3D> xyz_direction;
		size_t visited;
		perp_point(void) : visited() {}
		perp_point(coords::PES_Point const& pes_point)
			: pes(pes_point), visited()
		{}
		void swap(perp_point& rhs)
		{
			pes.swap(rhs.pes);
			main_direction.swap(rhs.main_direction);
			intern_direction.swap(rhs.intern_direction);
			xyz_direction.swap(rhs.xyz_direction);
			std::swap(visited, rhs.visited);
		}
		operator coords::PES_Point() const { return pes; }
	};

	struct GradCallBack
	{
		pathx* po;
		GradCallBack(pathx& path_object)
			: po(&path_object)
		{ }

		float operator() (scon::vector<scon::c3<float>> const& x,
			scon::vector<scon::c3<float>>& g, size_t const, bool& go_on)
		{
			using ot = scon::vector < scon::c3<float> >;
			using ct = coords::Representation_3D;
			po->cPtr->set_xyz(std::move(scon::explicit_transform<ct>(x)));
			auto const E = static_cast<float>(po->g_new(po->global_imagex));
			go_on = po->cPtr->integrity();
			g = std::move(scon::explicit_transform<ot>(po->cPtr->g_xyz()));
			return E;
		}
	};

public:

	pathx(neb* NEB, coords::Coordinates* c);
	coords::Representation_3D images;
	neb* N;
	neb* NEB;
	coords::Coordinates* c;
	coords::Coordinates* cPtr;
	coords::Representation_3D minima;
	std::vector <coords::Representation_3D> global_maxima, global_minima;
	ptrdiff_t global_image, counter, global_imagex;
	void pathx_ini();
	void MCM_PO(ptrdiff_t opt);
	void MC_PO(ptrdiff_t opt);
	bool testcoord(coords::Representation_3D& coords);
	void randvect();
	void proof_connect();
	double g_new(ptrdiff_t im);
	double lbfgs(ptrdiff_t imagex);
	void printmono(std::string const& name, coords::Representation_3D& print, ptrdiff_t& count);
	void move_main(perp_point& direction);
	std::vector<coords::Representation_Main> move_main_directions;

};