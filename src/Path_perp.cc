#include <cstdlib>
#include <iomanip>
#include <stdio.h>

#include "scon_vect.h"
#include "configuration.h"
#include "coords.h"
#include "Path_perp.h"
#include "scon_utility.h"
#include "ls.h"
#include "scon_vect.h"
#if defined(_WIN32)

	#include <direct.h>
	#include "sys/stat.h"
#else
	#include <unistd.h>
	#include <dirent.h>
    #include <sys/stat.h>
    
#endif


path_perp::path_perp(coords::Coordinates *c)
{
	cPtr = c;
	file_dimer = false;
}

void path_perp::initial(void)

{

 /* images_initial = cPtr->xyz();
  imagi[0] = cPtr->xyz();*/
  std::ifstream initial (Config::get().general.inputFilename);
  std::string buffer;
  //getline(final,buffer);
  getline(initial,buffer);
  size_t number;
  char atom[2];
  
  for (size_t i=0;i<cPtr->size();i++)
  {
    getline(initial,buffer);
    std::istringstream line_coord(buffer);
    line_coord >>number>>atom >> image1[i].x() >> image1[i].y() >> image1[i].z() ;
   	
  }

 
}

void path_perp::final(void)
{
  std::ifstream final (Config::get().neb.FINAL_STRUCTURE.c_str());
  std::string buffer;
  getline(final,buffer);
  //getline(final,buffer);
  size_t number;
  char atom[2];

  for (size_t i=0;i<cPtr->size();i++)
  {
    getline(final,buffer);
    std::istringstream line_coord(buffer);
    line_coord >>number>>atom >> image2[i].x() >> image2[i].y() >> image2[i].z() ;
  }
 
}

void path_perp::calc_tau (void)
{

	double abso(0.0);
	tau.resize(cPtr->size());
	   for (size_t j=0; j<cPtr->size(); j++){
          tau[j].x() = (image2[j].x() -image1[j].x());
          tau[j].y() = (image2[j].y() -image1[j].y());
          tau[j].z() = (image2[j].z() -image1[j].z());						  
		  abso = 0.0;
		  if (tau[j].x() - tau[j].x() != 0) tau[j].x() = 0.0;
		  if (tau[j].y() - tau[j].y() != 0) tau[j].y() = 0.0;
		  if (tau[j].z() - tau[j].z() != 0) tau[j].z() = 0.0;
		  abso = len(tau[j]);
		  if (abso != 0.0) normalize(tau[j]);

        }	
}



void path_perp::pathx_ini()
{
  
  std::cout << "**************INITIALIZATION OF PERPENDICULAR SEARCH*************"<<'\n';
  image1.resize(cPtr->size());
  image2.resize(cPtr->size());
  initial();
  final();
  calc_tau();
  STARTENERGY=cPtr->g();
  counter=0; 
  cPtr->set_xyz(image1);	

  MCM_NEB(0);
  //proof_connect();
  
}




//BASIN HOPPING
void path_perp::MCM_NEB(ptrdiff_t opt)
{
  Maxvar=Config::get().neb.VARIATION;
  std::vector <double> MCpmin_vec;
  double MCmin,MCpmin,MCgmin,factor;
  double boltzman=0.0, kt=1/(0.0019872966*Config::get().neb.TEMPERATURE), trial=(double)rand()/(double)RAND_MAX;
  ptrdiff_t mciteration=Config::get().neb.MCITERATION;
  ptrdiff_t nancounter{ 0 }, nbad{ 0 }, status{ 0 };
  bool  l_disp=false;
  MCSTEPSIZE=Config::get().neb.MCSTEPSIZE;
  counter =0;
  std::string out_opt("PATHOPT_PERP_SAMPLE_STRUCTURES");
  coords::Representation_3D positions;

  std::ofstream output2("PATHOPT_PERP_SAMPLE_ENERGIES.dat",std::ios::app);
  coords::Representation_3D tempstart,tempstart2,coord_in,coord_glob,coord_last;

  global_path_minima.resize(2);
  global_path_minima_temp.resize(2);
  global_path_minima_energy.resize(2);
  
  	tempstart = image1;
	tempstart2 = image2;
  

  
  for (size_t t=0; t<Config::get().neb.GLOBALITERATION;++t){

	 
 for (size_t i = 0; i < global_path_minima.size();++i)
  {
	   global_path_minima[i].resize(mciteration*(t+1));
	   global_path_minima_temp[i].resize(mciteration*(t+1));
	   global_path_minima_energy[i].resize(mciteration*(t+1));

  }
  //INITIALIZE COORDS


	coord_in = cPtr->xyz();
	coord_glob = cPtr->xyz();
  //CALCULATE SP ENERGY
  MCmin=cPtr->g();



  MCgmin=MCmin;
  MCpmin=MCmin;


  //ITERATIONS

  for (ptrdiff_t mcstep=0; mcstep<=mciteration;mcstep++){

	  coord_last = cPtr->xyz();

	  MCpmin_vec.push_back(MCmin);
	
	  positions = cPtr->xyz();
	  for(size_t j=0; j < cPtr->size();j++)
	  {
		randvect();
		factor = MCSTEPSIZE * (double)rand()/(double)RAND_MAX;
		double rv_p[3];
		coords::Cartesian_Point RV;
		double abs=0.0; 
		rv_p[0] = tau[j].y() * this->cPtr->g_xyz(j).z() - tau[j].z() *this->cPtr->g_xyz(j).y();
		rv_p[1] = tau[j].z() * this->cPtr->g_xyz(j).x() - tau[j].x() *this->cPtr->g_xyz(j).z();
		rv_p[2] = tau[j].x() * this->cPtr->g_xyz(j).y() - tau[j].y() *this->cPtr->g_xyz(j).x();
		abs += rv_p[0]*rv_p[0];
		abs += rv_p[1]*rv_p[1];
		abs += rv_p[2]*rv_p[2];
		abs= sqrt(abs);
		if(abs == 0.0) abs=1.0;
		rv_p[0] /= abs;
		rv_p[1] /= abs;
		rv_p[2] /= abs;
		RV.x() = (factor*rv_p[0]);
		RV.y() = (factor*rv_p[1]);
		RV.z() = (factor*rv_p[2]);	
		positions[j] += RV;
 
    }

	
	cPtr->set_xyz(positions);

	//this->cPtr->biascontrol=true;
	this->cPtr->NEB_control=false;
	
	//this->cPtr->g2(N->tau,opt,mcstep,this->N->image_ini);

    MCmin=lbfgs(opt);

	 
	 	

	if (MCmin != MCmin)
	{

		nancounter++;
		//std::cout << "***size of counter***:  "<<nancounter << lineend;
		status = 0;
		nbad++;
		cPtr->set_xyz(coord_in);
		if (nancounter>10) break;
	}
	//TEST FOR LOW ENERGIES THAT ARE NOT REASONABLE
	if (MCmin< (-1.5*std::abs(STARTENERGY)))
	{
		MCmin = 1000000.00;
		nbad++;
	}


	//METROPOLIS MONTE CARLO CRITERIUM AND TEST FOR IDENTICAL MINIMA
	for (auto mp : MCpmin_vec) { if (abs(MCmin - mp) < 0.001) status = 2; nbad = 0; }


	if (status == 2)
	{
		nbad = 0;
		//!same minimum
		MCpmin_vec[mcstep] = MCmin;
	}
	else if (MCmin < MCpmin_vec[mcstep])
	{
		nbad = 0;
		status = 1; //accepted as next minimum
		MCpmin_vec[mcstep] = MCmin;
	}
	else
	{
		nbad = 0;
		boltzman = exp(-kt*(MCmin - MCpmin_vec[mcstep]));
		trial = (double)rand() / (double)RAND_MAX;
		if (boltzman < trial)
		{
			status = 0;
		}
		else
		{
			status = 1;
			MCpmin_vec[mcstep] = MCmin;
		}
	}
	if (testcoord(coord_in))
	{
		l_disp = false;
	}
	else
	{
		//std::cout << "DISPLACEMENT TOO BIG" << lineend;
		l_disp = true;
		status = 0;
	}

	if (std::abs(MCmin - STARTENERGY) > Config::get().neb.PO_ENERGY_RANGE || std::abs(MCmin) == INFINITY)
	{
		MCmin = 1000000.00;
		status = 0;
		nbad++;
	}
	//SAVE ENERGY OF NEXT GLOBAL MINIMUM
	if ((MCmin < MCgmin) && status == 1)
	{
		coord_glob = cPtr->xyz();
		MCgmin = MCmin;
	}
	//RESTORE GLOBAL MINIMUM AFTER 3 BAD ITERATIONS
	if (nbad>3)
	{
		nbad = 0;
		cPtr->set_xyz(coord_glob);
	}
	//RESTORE COORDS FROM PREVIOUS ITERATION
	else if (status != 1)
	{
		cPtr->set_xyz(coord_last);
	}

    //std::cout << "ITERATION: " << mcstep << "    Current ENERGY: "<<MCmin <<"      " << " global MINIMUM" << MCgmin;
    //if(status == 0) std::cout << "REJECTED STRUCTURE\n";
    else if (status == 1)  {
      //std::cout << "     ACCEPT(B:T)" << boltzman << ":" << trial<<'\n';
      MCEN=MCmin;
     // writetinker(output,opt);
	  output2 <<mcstep <<"    "<< opt <<"    "<< MCEN << '\n';
	 
	  counter++;
	 
		global_path_minima_energy[opt][counter]=MCEN;
     for (size_t i=0; i< cPtr->size();i++){
		  global_path_minima[opt][counter].push_back(cPtr->xyz(i));

	 }
     
	 printmono(out_opt,global_path_minima[opt][counter],opt);

    }
    else if (status == 2) /*std::cout << "SAME STRUCTURE\n";*/
    status=0;
	

  }
  
  }





}





//GENERATES 3D-UNIT VECTOR
//G: MARSAGLIA, ANN. MATH. STAT., 43, 645 (1972)
void path_perp::randvect(void)
{
  double x(0.0), y(0.0), s(0.0);
  s=2.0;
  while(s>=1.0){
    x=2.0 * (double)rand()/(double)RAND_MAX -1.0;
    y=2.0 * (double)rand()/(double)RAND_MAX -1.0;
    s= x*x + y*y;
  }
  randvec[2]=1.0-2.0*s;
  s=2.0 * sqrt(1.0-s);
  randvec[1]=s*y;
  randvec[0]=s*x;



}



bool path_perp::testcoord(coords::Representation_3D &coords)
{

	for (size_t i =0; i<cPtr->size();i++){

    if(abs(cPtr->xyz(i).x() - coords[i].x())   > Maxvar) return false;
	if(abs(cPtr->xyz(i).y() - coords[i].y())   > Maxvar) return false;
	if(abs(cPtr->xyz(i).z() - coords[i].z())   > Maxvar) return false;
	
	}
  
  
  return true;


}




struct GradCallBack
{
  path_perp * p;
  GradCallBack(path_perp & pp) : p(&pp) { }
  float operator() (scon::vector< scon::c3 <float> > const &x,
    scon::vector< scon::c3 <float> > &g, size_t const, bool & go_on)
  {
      p->cPtr->set_xyz(coords::Representation_3D(x.begin(), x.end()));
      auto const E = (float)p->g_new();
      go_on = p->cPtr->integrity();
      p->cPtr->g_xyz();
      g = scon::vector< scon::c3 <float> >(p->cPtr->g_xyz().begin(), p->cPtr->g_xyz().end());
      return E;
  }
};



double path_perp::lbfgs (ptrdiff_t imagex)
{
  using namespace optimization::local;
  global_imagex=imagex;
  auto optimizer = make_lbfgs(make_more_thuente(GradCallBack(*this)));
  using op_type = decltype(optimizer);
  using vec_type = op_type::rep_type;
  op_type::point_type point(vec_type(cPtr->xyz().begin(), cPtr->xyz().end()));
  optimizer(point);
  if (Config::get().general.verbosity > 3 && optimizer.state() != status::SUCCESS)
  {
    std::cout << "Optimization returned " << optimizer.state() << '\n';
  }
  if (Config::get().general.verbosity > 3)
  {
    std::cout << "Optimization done. Evaluations:" << optimizer.iter() << '\n';
  }
  return point.f;
}



double path_perp::g_new()
{

	double cosi,energytemp,pot(0.0),dpot(0.0);
	coords::Cartesian_Point rv_p;
	energytemp=cPtr->g();
	if(Config::get().neb.OPTMODE == "PROJECTED")
	{
	for (size_t i=0;i<cPtr->size();i++)
	{
    auto const l = scon::geometric_length(tau[i]);
    if (l > 0.0)
	  {
      cosi = (cPtr->g_xyz(i).x()*tau[i].x() + cPtr->g_xyz(i).y()*tau[i].y() + cPtr->g_xyz(i).z()*tau[i].z()) / (l*l);
	  }else
	  {
		   cosi= 0.0;
	  }
	  rv_p.x()=  cPtr->g_xyz(i).x() - (cosi * tau[i].x());
	  rv_p.y() = cPtr->g_xyz(i).y() - (cosi * tau[i].y());
	  rv_p.z() = cPtr->g_xyz(i).z() - (cosi * tau[i].z());
	  
	  cPtr->update_g_xyz(i,rv_p);

	}
	}else if(Config::get().neb.OPTMODE == "BIAS")
	{
	for (size_t i=0;i<cPtr->size();i++)
	{

		 if(scon::geometric_length(tau[i]) > 0.0)
	  {
		cosi=(cPtr->xyz(i).x()*tau[i].x()+cPtr->xyz(i).y()*tau[i].y()+cPtr->xyz(i).z()*tau[i].z()-tau[i].x()*image1[i].x()-tau[i].y()*image1[i].y()-tau[i].z()*image1[i].z())/sqrt(tau[i].x()*tau[i].x()+tau[i].y()*tau[i].y()+tau[i].z()*tau[i].z());
		 }else
		 {
			 cosi=0.0;
		 }
		cosi= cosi-0.0;
		pot += cosi*cosi;
		dpot=cosi*2*Config::get().neb.BIASCONSTANT;

		rv_p.x()= cPtr->g_xyz(i).x() + dpot*cPtr->xyz(i).x();
		rv_p.y()= cPtr->g_xyz(i).y() + dpot*cPtr->xyz(i).y();
		rv_p.z()= cPtr->g_xyz(i).z() + dpot*cPtr->xyz(i).z();

		cPtr->update_g_xyz(i,rv_p);

	}
	
	cosi=0.0;
	}
  return energytemp;
}

void path_perp::printmono(std::string const &name, coords::Representation_3D &print,ptrdiff_t &count)
{

  std::ofstream out (name.c_str(),std::ios::app);
  std::string temp;


  out << "     " << cPtr->size() << "  global counter:  "<<count <<'\n';
  for(size_t j = 0; j <cPtr->size() ;j++){

    out << std::right << std::setw(6) << j+1;
    out << std::left << "  " << std::setw(12)<<cPtr->atoms(j).symbol().substr(0U, 2U);	
    out << std::fixed << std::right << std::setprecision(6) << std::setw(13) << print[j].x();
    out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[j].y();
    out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[j].z();
    out << std::right << std::setw(6)<<cPtr->atoms(j).energy_type();
    size_t const bSize(cPtr->atoms(j).bonds().size());
    for (size_t n(0U); n<bSize; ++n)
    {
      out << std::right << std::setw(6) << cPtr->atoms(j).bonds()[n]+1U;
    }
    out << '\n';

  }

  out.close();


}