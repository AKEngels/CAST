#include "error.h"
#include <vector>
#include "coords.h"
#include <ostream>
#include <memory>
#include <string>
#include "scon_vect.h"
#include "coords_rep.h"
#include "neb.h"
#include <cstdlib>
#include <iomanip>



class pathx
{



private:

  void MCM_NEB(ptrdiff_t opt);

  void randvect();


  double randvec[3];


  ptrdiff_t nvar, natom;

  double Maxvar, MCSTEPSIZE;
  double STARTENERGY, ENDENERGY;
  std::vector<double>  MCcoordglob, MCcoordlast, MCcoordin;
  std::vector <ptrdiff_t> plane1, plane2, plane3;
  std::vector <double> ndist1, ndist2, ndist3;
  std::vector <std::vector <double> > path_ready, global_path_minima_energy;
  std::vector <coords::Representation_3D > PS_MIN_A, PS_MIN_B, PS_TS;
  std::vector <coords::Representation_3D > PS_MIN_A_TEMP, PS_MIN_B_TEMP, PS_TS_TEMP;
  std::vector < std::vector <coords::Representation_3D> >  global_path_minima, global_path_minima_temp;


  std::vector <std::vector <size_t> > segindex;
  std::vector <bool> tabu_entries1, tabu_entries3;
  //vect<double> randvec_perp;
  bool testcoord(coords::Representation_3D &coords);
  bool file_dimer;
  bool eins_drei, controlseg;

  void proof_connect();
  double g_new(ptrdiff_t im);
  double lbfgs();
  void printmono(std::string const &name, coords::Representation_3D &print, ptrdiff_t &count);
 
  




  struct perp_point
  {
	  coords::PES_Point pes;
	  std::vector<coords::Representation_Main> main_direction;
	  std::vector<coords::Representation_Internal> intern_direction;
	  std::vector<coords::Representation_3D> xyz_direction;
	  size_t visited;
	  perp_point(void) : visited() {}
	  perp_point(coords::PES_Point const &pes_point)
		  : pes(pes_point), visited()
	  {}
	  void swap(perp_point &rhs)
	  {
		  pes.swap(rhs.pes);
		  main_direction.swap(rhs.main_direction);
		  intern_direction.swap(rhs.intern_direction);
		  xyz_direction.swap(rhs.xyz_direction);
		  std::swap(visited, rhs.visited);
	  }
	  operator coords::PES_Point() const { return pes; }
  };


  void move_main(perp_point & direction);
  std::vector<coords::Representation_Main> move_main_directions;











  //void proof_segments(void);
  //bool proof_list(std::vector <double> &, size_t,std::vector <double> e1,std::vector <double> e2, std::vector <double> e3,size_t);
  //void call(void);
  //void prepare_tinker(void);
  //int call_tinker(void);

  //   void print_tinkerInput (void);
  //   void read_tinkerOutput (bool const grad = true, bool const hess = false, bool const opt = false);

  struct pathready {


    size_t index;

  };


  struct GradCallBack
  {
    GradCallBack(pathx & path_object)
      : po(&path_object)
    { }
    pathx * po;
    float operator() (scon::vector<scon::c3<float>> const &x,
      scon::vector<scon::c3<float>> &g, size_t const, bool & go_on)
    {
      using c3fv = scon::vector < scon::c3<float> > ;
      po->cPtr->set_xyz(std::move(scon::explicit_transform<coords::Representation_3D>(x)));
      auto const E = static_cast<float>(po->g_new(po->global_imagex));
      go_on = po->cPtr->integrity();
      g = std::move(scon::explicit_transform<c3fv>(x));
      return E;
    }
  };

  std::vector <pathready> pathready_index;


public:
  void pathx_ini();
  std::string id;
  coords::Representation_3D images;
  neb *N;
  neb *NEB;
  coords::Coordinates *c;
  coords::Coordinates *cPtr;

  ptrdiff_t counter, global_imagex;

  pathx(neb *NEB, coords::Coordinates *c);
  coords::Representation_3D minima;
  std::vector <coords::Representation_3D> global_maxima, global_minima, global_maxima2;
  double MCEN;
  std::vector <double>  energy1, energy2, energy3;
  ptrdiff_t global_image, global_vibrate, mcit, temp_vibrate, global_run, total_struc_num;
  double globalenergy;
  std::vector <double> rmsd_list;
  std::string outstring;
  std::string keystring;
  std::string logstring;

  std::vector <coords::Representation_3D > perp;






};