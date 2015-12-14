#include "error.h"
#include <vector>
#include "coords.h"
#include <ostream>
#include <memory>
#include <string>
#include "scon_vect.h"
#include <cstdlib>
#include <iomanip>

using namespace scon;

class path_perp
	{

   

	private:
   
	void MCM_NEB (ptrdiff_t opt);

	void randvect ();

	double randvec[3];

	ptrdiff_t nvar,natom;

	
	double STARTENERGY,ENDENERGY;
	std::vector<double>  MCcoordglob, MCcoordlast, MCcoordin;
	std::vector <std::vector <double> > path_ready, global_path_minima_energy;
	std::vector < std::vector <coords::Representation_3D> >  global_path_minima,global_path_minima_temp;
	coords::Representation_3D image1,image2,tau;

	std::vector <std::vector <size_t> > segindex;
	//vect<double> randvec_perp;
	bool testcoord (coords::Representation_3D &coords);
	bool file_dimer;
	bool eins_drei,controlseg;
	void calc_tau (void);
	double lbfgs(ptrdiff_t im);
	void final(void);
	void initial(void);
	void printmono(std::string const &name, coords::Representation_3D &print,ptrdiff_t &count);

	struct pathready {

		size_t index;

	};

	std::vector <pathready> pathready_index;

	public:

    double g_new(ptrdiff_t im);

	void pathx_ini (void);
    std::string id;

	coords::Coordinates *c;
	coords::Coordinates *cPtr;

	ptrdiff_t counter,global_imagex;

	path_perp ( coords::Coordinates *c);
	//pathx (NEB::pathway *NEB, coordinates *c, coordinates_internal *C_I);
	path_perp (void);
	path_perp&  operator=    (path_perp);
	coords::Representation_3D minima;
	double MCEN, Maxvar, MCSTEPSIZE;
	size_t count,mciteration;
	ptrdiff_t global_image,global_vibrate,mcit,temp_vibrate,global_run,total_struc_num;
	double globalenergy;



	};