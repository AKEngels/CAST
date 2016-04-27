#include "scon_vect.h"
#include "configuration.h"
#include "coords.h"
//#include "neb.h"
#include <cstdlib>
#include <iomanip>
#include "pathopt.h"
#include <stdio.h>



#include "scon_utility.h"
#include "ls.h"
#include "lbfgs.h"
#if defined(_WIN32)

#include <direct.h>
#include "sys/stat.h"
#else
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

#endif


pathx::pathx(neb *NEB, coords::Coordinates *c)
{
  N = NEB;
  cPtr = c;
  global_image = 0;
}

void pathx::pathx_ini()
{

  std::cout << "**************INITIALIZATION OF PATHOPT*************\n";
  natom = cPtr->size();
  nvar = 3 * natom;
  global_minima.resize(this->N->num_images);
  STARTENERGY = N->energies[0];
  ENDENERGY = N->energies[N->num_images - 1];
  ptrdiff_t temp_image;
  temp_image = N->num_images;
  counter = 0;
  for (ptrdiff_t i = 1; i < temp_image - 1; i++)
  {
    N->num_images = temp_image;
    cPtr->set_xyz(N->imagi[i]);
    MCM_NEB(i);
  }
  proof_connect();
}




//BASIN HOPPING
void pathx::MCM_NEB(ptrdiff_t opt)
{
  Maxvar = Config::get().neb.VARIATION;
  double MCmin, MCgmin, factor,MCpmin;
  std::vector <double> MCpmin_vec;
  double boltzman = 0.0, kt = 1 / (0.0019872966*Config::get().neb.TEMPERATURE), trial = (double)rand() / (double)RAND_MAX;
  ptrdiff_t mciteration = Config::get().neb.MCITERATION;
  ptrdiff_t nancounter = 0, nbad = 0, status;
  bool  l_disp = false;
  MCSTEPSIZE = Config::get().neb.MCSTEPSIZE;
  std::string output, output3, output4;

  global_image = 0;
  counter = 0;
  coords::Representation_3D positions;

  std::ofstream output2("PATHOPT_BASIN_ENERGIES", std::ios::app);
  double sum, sum2;
  coords::Representation_3D  coord_in, coord_glob, coord_last;

  global_path_minima.resize(N->num_images+natom);
  global_path_minima_temp.resize(N->num_images+natom);
  global_path_minima_energy.resize(N->num_images+natom);





  for (size_t t = 0; t < Config::get().neb.GLOBALITERATION; t++) {


    for (size_t i = 0; i < global_path_minima.size(); i++)
    {
      global_path_minima[i].resize(mciteration*(t + 1));
      global_path_minima_temp[i].resize(mciteration*(t + 1));
      global_path_minima_energy[i].resize(mciteration*(t + 1));

    }
    //INITIALIZE COORDS


    coord_in = cPtr->xyz();
    coord_glob = cPtr->xyz();



	//CALCULATE SP ENERGY
    MCmin = cPtr->g();
    MCgmin = MCmin;
	MCpmin = MCmin;


    //ITERATIONS

    for (ptrdiff_t mcstep = 0; mcstep <= mciteration; mcstep++) {

      bool move_control(false);
      coord_last = cPtr->xyz();
	  MCpmin_vec.push_back(MCmin);

      sum = 0.0;
	  sum2 = 0.0;
	  positions.clear();
      positions.resize(natom);
	

	  move_control = Config::get().neb.MIXED_MOVE;
	  if (move_control && (mcstep % 10 == 0 && mcstep > 1)) 
		  {

			  for (ptrdiff_t j = 0; j < natom; j++) {
				  randvect();
				  factor = MCSTEPSIZE * (double)rand() / (double)RAND_MAX;
				  double rv_p[3];
				  coords::Cartesian_Point RV;
				  double abs = 0.0;

				  rv_p[0] = N->tau[opt][j].y() * randvec[2] - N->tau[opt][j].z() * randvec[1];
				  rv_p[1] = N->tau[opt][j].z() * randvec[0] - N->tau[opt][j].x() * randvec[2];
				  rv_p[2] = N->tau[opt][j].x() * randvec[1] - N->tau[opt][j].y() * randvec[0];
				  abs += rv_p[0] * rv_p[0];
				  abs += rv_p[1] * rv_p[1];
				  abs += rv_p[2] * rv_p[2];

				  abs = sqrt(abs);
				  if (abs == 0.0) abs = 1.0;
				  rv_p[0] /= abs;
				  rv_p[1] /= abs;
				  rv_p[2] /= abs;

				  RV.x() = rv_p[0];
				  RV.y() = rv_p[1];
				  RV.z() = rv_p[2];

				  positions[j] = RV;

			  }
			  coords::Coordinates temp_perp(*cPtr);
			  temp_perp.set_xyz(positions);
			  perp_point test_perp;
			  test_perp = perp_point(temp_perp.pes());
			  move_main(test_perp);
			  positions = cPtr->xyz();
			  //cPtr->o();
			  positions = cPtr->xyz();	  
	  }
	  else
	  {
		  positions = cPtr->xyz();
        for (ptrdiff_t j = 0; j < natom; j++) {
          randvect();
          factor = MCSTEPSIZE * (double)rand() / (double)RAND_MAX;
          double rv_p[3];
          coords::Cartesian_Point RV;
          double abs = 0.0;

          rv_p[0] = N->tau[opt][j].y() * randvec[2] - N->tau[opt][j].z() * randvec[1];
          rv_p[1] = N->tau[opt][j].z() * randvec[0] - N->tau[opt][j].x() * randvec[2];
          rv_p[2] = N->tau[opt][j].x() * randvec[1] - N->tau[opt][j].y() * randvec[0];
          abs += rv_p[0] * rv_p[0];
          abs += rv_p[1] * rv_p[1];
          abs += rv_p[2] * rv_p[2];

          abs = sqrt(abs);
          if (abs == 0.0) abs = 1.0;
          rv_p[0] /= abs;
          rv_p[1] /= abs;
          rv_p[2] /= abs;

          RV.x() = (factor*rv_p[0]);
          RV.y() = (factor*rv_p[1]);
          RV.z() = (factor*rv_p[2]);

          positions[j] += RV;

        }
      }
 
      cPtr->set_xyz(positions);

      this->cPtr->mult_struc_counter = 0;
      //this->cPtr->biascontrol=true;
      this->cPtr->NEB_control = false;

      //this->cPtr->g2(N->tau,opt,mcstep,this->N->image_ini);
	  MCmin = lbfgs(opt);
  
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
	for (auto mp : MCpmin_vec) { if (abs(MCmin - mp) < 0.001 ) status = 2; nbad = 0; }


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
	//if(status == 0) std::cout << "REJECTED STRUCTURE" << lineend;
	else if (status == 1) 
	{
		MCEN = MCmin;
		global_image = opt;
		std::ostringstream struc_opt;
		struc_opt << "PATHOPT_STRUCTURES_" << opt << ".xyz";
		output2 << mcstep << "    " << opt << "    " << MCEN << "\n";
		counter++;
		global_path_minima_energy[opt][counter] = MCEN;
		for (ptrdiff_t i = 0; i<natom; i++)
		{
			global_path_minima[opt][counter].push_back(cPtr->xyz(i));
		}
	printmono(struc_opt.str(), global_path_minima[opt][counter], counter);
	}
	else if (status == 2) status = 0;
    }
  }
}



void pathx::proof_connect()
{


  natom = cPtr->size();
  STARTENERGY = N->energies[0];
  ENDENERGY = N->energies[N->num_images - 1];
  ptrdiff_t temp_image;
  temp_image = N->num_images;
  ptrdiff_t image = Config::get().neb.CONNECT_NEB_NUMBER;

  std::vector < std::vector <std::vector  <double> > > RMSD;
  std::vector < std::vector <ptrdiff_t> >PARTNER;
  std::vector < std::vector <ptrdiff_t> >  TABU;
  coords::Cartesian_Point cog1, cog2, tempcoord1, tempcoord2;
  coords::Coordinates  coord1, coord2;
  PARTNER.resize(N->num_images);
  TABU.resize(N->num_images);
  RMSD.resize(N->num_images);

  for (size_t l = 0; l < RMSD.size(); l++)RMSD[l].resize(Config::get().neb.MCITERATION*(Config::get().neb.GLOBALITERATION));
  for (size_t l = 0; l < PARTNER.size(); l++)PARTNER[l].resize(Config::get().neb.MCITERATION*(Config::get().neb.GLOBALITERATION));
  for (size_t l = 0; l < TABU.size(); l++)TABU[l].resize(Config::get().neb.MCITERATION*(Config::get().neb.GLOBALITERATION));


  bool reverse(false);
  coords::Representation_3D tempstart, tempstart2;
  tempstart = N->imagi[0];
  tempstart2 = N->imagi[N->num_images-1];

  coord1.clear();
  coord2.clear();

  coord1.init_in(cPtr->atoms(), coords::PES_Point(N->imagi[0]), true);
  coord2.init_in(cPtr->atoms(), coords::PES_Point(N->imagi[0]), true);






  for (ptrdiff_t i = 1; i < temp_image - 1; i++)
  {
    N->num_images = temp_image;
    coord2.set_xyz(N->imagi[0]);
    cog2 = coord2.center_of_geometry();

    for (size_t j = 0; j < global_path_minima[i].size(); j++)
    {
      if (global_path_minima[i][j].size() == 0) continue;

      //coord1.init_in(cPtr->atoms(), coords::PES_Point(global_path_minima[i][j]), true);
      coord1.set_xyz(global_path_minima[i][j]);

      cog1 = coord1.center_of_geometry();

      for (size_t l = 0; l < global_path_minima[i][j].size(); l++)
      {
        tempcoord2 = coord1.xyz(l);
        global_path_minima_temp[i][j].push_back(scon::cog_align(tempcoord2, cog1, cog2));
      }

      global_path_minima[i][j] = global_path_minima_temp[i][j];

    }




  }



  for (ptrdiff_t i = 1; i < temp_image - 1; i++)
  {
    N->num_images = temp_image;

    for (size_t j = 0; j < global_path_minima[i].size(); j++)
    {
      for (size_t n = 0; n < global_path_minima[i].size(); n++)
      {
        if (global_path_minima[i][j].empty())continue;
        auto rmsd1 = root_mean_square_deviation(global_path_minima[i][j], global_path_minima[i + 1][n]);
        if (rmsd1 == 0.0 || rmsd1 != rmsd1)continue;
        RMSD[i][j].push_back(rmsd1);

      }
    }

  }


  for (size_t mm = 1; mm < 4; mm++) {


    std::ostringstream NEB1;

    NEB1 << mm;

#if defined(_MSC_VER)
    _mkdir(NEB1.str().c_str());
    _chdir(NEB1.str().c_str());
#else
    mkdir(NEB1.str().c_str(), 0777);
    chdir(NEB1.str().c_str());
#endif





    for (ptrdiff_t i = 1; i < temp_image - 1; i++)
    {


      for (size_t j = 0; j < global_path_minima[i].size(); j++)
      {
        if (RMSD[i][j].empty()) continue;
        size_t index = 0;
        double min = abs(RMSD[i][j][0]);
        for (size_t p = 0; p < RMSD[i][j].size(); p++) {



          if (abs(RMSD[i][j][p]) < min) {
            index = p;
            min = abs(RMSD[i][j][p]);
          }


        }

        PARTNER[i][j] = (index + (mm - 1));
        RMSD[i][j].erase(RMSD[i][j].begin() + index);


      }
    }



    ptrdiff_t i = 1, jj;
    ptrdiff_t tempcount = 0;
    bool tempproff = false;
    /*ptrdiff_t j=2;*/
    if (Config::get().neb.NEB_CONN == true)
    {

      for (ptrdiff_t j = 0; j < ptrdiff_t(global_path_minima[1].size()); j++)

      {
        if (global_path_minima[1][j].empty()) continue;

        reverse = false;

        N->preprocess(j, image, j, global_path_minima[1][j], tempstart, reverse);
        reverse = true;
        tempcount = 0;
        tempproff = false;
        jj = j;
        for (i = 1; i < temp_image - 2; i++)
        {


          if (global_path_minima[i - tempcount][jj].empty()) continue;
          if (global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]].empty()) continue;
          //if(global_path_minima[(i-tempcount)+1][PARTNER[(i-tempcount)+1][jj]].empty()){/*tempcount++;*/ tempproff=true;continue;}
          N->preprocess(j, image, j, global_path_minima[i - tempcount][jj], global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]], reverse);
          jj = PARTNER[i - tempcount][jj];

          tempproff = false;

        }
        /*if(global_path_minima[i-tempcount][jj].empty()) continue ;*/
        /*if(tempproff==true)continue;*/

        reverse = true;
        /*if(tempcount < temp_image-2) N->preprocess(j,image,j,global_path_minima[tempcount+1][j],tempstart2,reverse);*/
        /*if(tempproff=true) {N->preprocess(j,image,i,global_path_minima[tempcount+1][j],tempstart2,reverse);}*/
		if (global_path_minima[i - tempcount][jj].size() == 0)continue;
        N->preprocess(j, image, j, global_path_minima[i - tempcount][jj], tempstart2, reverse);
      }
    }
    else
    {
      for (ptrdiff_t j = 0; j < ptrdiff_t(global_path_minima[1].size()); j++)

      {
        if (global_path_minima[1][j].empty()) continue;

        reverse = false;


        std::ostringstream en, img;
        en << "ENERGY_" << j << ".dat";
        img << "PATH_" << j << ".xyz";
        std::ofstream energy(en.str().c_str(), std::ios::app);

        printmono(img.str().c_str(), tempstart, j);
        cPtr->set_xyz(tempstart);
        energy << cPtr->g() << '\n';
        printmono(img.str().c_str(), global_path_minima[1][j], j);
        cPtr->set_xyz(global_path_minima[1][j]);
        energy << cPtr->g() << '\n';

        tempcount = 0;
        tempproff = false;
        jj = j;
        for (i = 1; i < temp_image - 2; i++)
        {
          if (global_path_minima[i - tempcount][jj].empty()) continue;
          if (global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]].empty()) continue;
          // if(global_path_minima[(i-tempcount)+1][PARTNER[(i-tempcount)+1][jj]].empty()){/*tempcount++;*/ tempproff=true;continue;}
          printmono(img.str().c_str(), global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]], j);
          cPtr->set_xyz(global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]]);
          energy << cPtr->g() << '\n';
          jj = PARTNER[i - tempcount][jj];

          tempproff = false;

        }
        //if(global_path_minima[i-tempcount][jj].empty()) continue ;
      //if(tempproff==true)continue;

        printmono(img.str().c_str(), tempstart2, j);
        cPtr->set_xyz(tempstart2);
        energy << cPtr->g() << '\n';
      }







    }


#if defined (_MSC_VER)
    _chdir("..");
#else
    chdir("..");
#endif
  }






}


//GENERATES 3D-UNIT VECTOR
//G: MARSAGLIA, ANN. MATH. STAT., 43, 645 (1972)
void pathx::randvect(void)
{
  double x(0.0), y(0.0), s(0.0);
  s = 2.0;
  while (s >= 1.0) {
    x = 2.0 * (double)rand() / (double)RAND_MAX - 1.0;
    y = 2.0 * (double)rand() / (double)RAND_MAX - 1.0;
    s = x*x + y*y;
  }
  randvec[2] = 1.0 - 2.0*s;
  s = 2.0 * sqrt(1.0 - s);
  randvec[1] = s*y;
  randvec[0] = s*x;



}



bool pathx::testcoord(coords::Representation_3D &coords)
{

  for (size_t i = 0; i<cPtr->size(); i++) {

    if (abs(cPtr->xyz(i).x() - coords[i].x()) > Maxvar) return false;
    if (abs(cPtr->xyz(i).y() - coords[i].y()) > Maxvar) return false;
    if (abs(cPtr->xyz(i).z() - coords[i].z()) > Maxvar) return false;

  }


  return true;


}








double pathx::lbfgs(ptrdiff_t imagex)
{

	using namespace  optimization::local;
	//typedef coords::Container<scon::c3<float>> nc3_type;
	// Create optimizer
	global_imagex = imagex;
	coords::Cartesian_Point g;
	auto optimizer = make_lbfgs(
		make_more_thuente(GradCallBack(*this))
	);
	//optimizer.ls.config.ignore_callback_stop = true;
	// Create Point
	using op_type = decltype(optimizer);
	op_type::point_type x(scon::explicit_transform<op_type::rep_type>(cPtr->xyz()));
	// Optimize point
	optimizer.config.max_iterations = Config::get().optimization.local.bfgs.maxstep;
	//optimizer.config.max_iterations = Config::get().neb.LBFGS_IT;
	optimizer.config.epsilon = (float)Config::get().optimization.local.bfgs.grad;
	optimizer(x);

	if (Config::get().general.verbosity > 4)
	{
		std::cout << "Optimization done (status " << optimizer.state() << "). Evaluations:" << optimizer.iter() << '\n';
	}


 


  return optimizer.p().f;
}



double pathx::g_new(ptrdiff_t im)
{

  double cosi, energytemp, pot(0.0), dpot(0.0);
  coords::Cartesian_Point rv_p;
  energytemp = cPtr->g();
  if (Config::get().neb.OPTMODE == "PROJECTED")
  {
    for (size_t i = 0; i < cPtr->size(); i++)
    {

      auto L = scon::geometric_length(N->tau[im][i]);
	

      if (L != 0.0)
      {
        cosi = (cPtr->g_xyz(i).x()*N->tau[im][i].x() + cPtr->g_xyz(i).y()*N->tau[im][i].y() + cPtr->g_xyz(i).z()*N->tau[im][i].z()) / (L * L);
	
      }
      else
      {
        cosi = 0.0;
      }
      rv_p.x() = cPtr->g_xyz(i).x() - (cosi * N->tau[im][i].x());
      rv_p.y() = cPtr->g_xyz(i).y() - (cosi * N->tau[im][i].y());
      rv_p.z() = cPtr->g_xyz(i).z() - (cosi * N->tau[im][i].z());

      cPtr->update_g_xyz(i, rv_p);

    }
  }
  else if (Config::get().neb.OPTMODE == "BIAS")
  {
    for (size_t i = 0; i < cPtr->size(); i++)
    {

      if (scon::geometric_length(N->tau[im][i]) != 0.0)
      {
        cosi = (cPtr->xyz(i).x()*N->tau[im][i].x() + cPtr->xyz(i).y()*N->tau[im][i].y() + cPtr->xyz(i).z()*N->tau[im][i].z() - N->tau[im][i].x()*N->image_ini[im][i].x() - N->tau[im][i].y()*N->image_ini[im][i].y() - N->tau[im][i].z()*N->image_ini[im][i].z()) / sqrt(N->tau[im][i].x()*N->tau[im][i].x() + N->tau[im][i].y()*N->tau[im][i].y() + N->tau[im][i].z()*N->tau[im][i].z());
      }
      else
      {
        cosi = 0.0;
      }
      cosi = cosi - 0.0;
      pot += cosi*cosi;
      dpot = cosi * 2 * Config::get().neb.BIASCONSTANT;

      rv_p.x() = cPtr->g_xyz(i).x() + dpot*cPtr->xyz(i).x();
      rv_p.y() = cPtr->g_xyz(i).y() + dpot*cPtr->xyz(i).y();
      rv_p.z() = cPtr->g_xyz(i).z() + dpot*cPtr->xyz(i).z();

      cPtr->update_g_xyz(i, rv_p);

    }

    cosi = 0.0;
  }
  return energytemp;
}


void pathx::printmono(std::string const &name, coords::Representation_3D &print, ptrdiff_t &count)
{

  std::ofstream out(name.c_str(), std::ios::app);
  std::string temp;


  out << "     " << natom << "  global counter:  " << count << '\n';
  for (ptrdiff_t j = 0; j < natom; j++) {

    out << std::right << std::setw(6) << j + 1;
    out << std::left << "  " << std::setw(12) << cPtr->atoms(j).symbol().substr(0U, 2U);
    out << std::fixed << std::right << std::setprecision(6) << std::setw(13) << print[j].x();
    out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[j].y();
    out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[j].z();
    out << std::right << std::setw(6) << cPtr->atoms(j).energy_type();
    size_t const bSize(cPtr->atoms(j).bonds().size());
    for (size_t n(0U); n < bSize; ++n)
    {
      out << std::right << std::setw(6) << cPtr->atoms(j).bonds()[n] + 1U;
    }
    out << '\n';

  }

  out.close();


}

void pathx::move_main(perp_point & direction)
{

  size_t const NUM(cPtr->main().size());

  coords::Representation_Main main_tors(NUM), bla(NUM);

  // choose number of torsions to modify
  coords::float_type const NUM_MAINS(static_cast<coords::float_type>(NUM));
  size_t const NUM_MOD(std::min((static_cast<size_t>(-std::log(scon::rand<coords::float_type>(0.0, 1.0))) + 1U), NUM));
  // apply those torsions
  if (Config::get().general.verbosity > 9U) std::cout << "Changing " << NUM_MOD << " of " << NUM << " mains.\n";
  for (size_t i(0U); i < NUM_MOD; ++i)
  {
    size_t const K(static_cast<size_t>(NUM_MAINS*scon::rand<coords::float_type>(0.0, 1.0)));
    coords::float_type const F = scon::rand<coords::float_type>(-Config::get().optimization.global.montecarlo.dihedral_max_rot,
      Config::get().optimization.global.montecarlo.dihedral_max_rot);
    if (Config::get().general.verbosity > 9U) std::cout << "Changing main " << K << " by " << F << '\n';
    main_tors[K] = coords::main_type::from_deg(F);
  }
  // orthogonalize movement to main directions
  for (auto perp_direction : direction.main_direction)
  {

    perp_direction.resize(main_tors.size());
    scon::orthogonalize_to_normal(main_tors, perp_direction);
  }
  for (size_t i = 0; i < NUM; ++i)
  {
    cPtr->rotate_main(i, main_tors[i]);
  }
  //main_tors.print(cout);
  cPtr->to_xyz();
}