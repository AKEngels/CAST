#include <cstdlib>
#include <iomanip>
#include "pathopt.h"
#include <stdio.h>
#include "scon_utility.h"
#include "ls.h"
#include "scon_vect.h"
#include "configuration.h"
#include "coords.h"
#include "lbfgs.h"
#if defined(_WIN32)
#include <direct.h>
#include "sys/stat.h"
#else
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

#endif
#if defined(_MSC_VER)
#include <process.h>
#define pid_func _getpid
#else 
#include <unistd.h>
#define pid_func getpid
#endif

/** 
* CONSTRUCTOR OF PATHX-CLASS AND INITIALIZATION OF GLOBAL VARIABLES
*/
pathx::pathx(neb *NEB, coords::Coordinates *c) :_KT_(1 / (0.0019872966*Config::get().neb.TEMPERATURE))
{
  N = NEB;
  cPtr = c;
  global_image = 0;
  mcstepsize = Config::get().neb.MCSTEPSIZE;
  maxvar = Config::get().neb.VARIATION;
  mciteration = Config::get().neb.MCITERATION;
  STARTENERGY = N->energies[0];
  ENDENERGY = N->energies[N->num_images - 1];
}

/**
* INITIALIZATION FUNCTION FOR PATHOPT RUN
*/
void pathx::pathx_ini()
{
  std::cout << "**************INITIALIZATION OF PATHOPT*************\n";
  global_minima.resize(this->N->num_images);
  ptrdiff_t temp_image;
  temp_image = N->num_images;
  for (ptrdiff_t i = 1; i < temp_image - 1; i++)
  {
	   
    N->num_images = temp_image;
	/**
	* initialize global optimization by using NEB starting structure 
	*/
    cPtr->set_xyz(N->imagi[i]);
	/**
	* performing global optimization on i-th n-1 dimensional hyperplane 
	*/

	if (Config::get().neb.MCM_OPT)
	{
		MCM_PO(i);
	}
	else
	{
		MC_PO(i);
	}
		
	
  }
  /**
  * Connection procedure on the total number of obtained minima
  * after optimizaiton on n-1 subspace
  */
  if (Config::get().neb.CONN)
  {
	  proof_connect();
  }
 
}

/**
* BASIN HOPPING BY USING MODIFIED RANDOM JUMPS AND CONSTRAINT MINIMIZATION 
*/
void pathx::MCM_PO(ptrdiff_t opt)
{
  
  double MCmin{ 0.0 }, MCgmin{ 0.0 }, factor{ 0.0 };
  std::vector <double> MCpmin_vec;
  /**
  * initialize Boltzman and trial number generation
  */
  double boltzman { 0.0 }, trial = (double)rand() / (double)RAND_MAX;
  ptrdiff_t nancounter (0), nbad (0), status (0);
  bool  nanstatus(false);
  global_image = 0;
  counter = 0;
  coords::Representation_3D positions;
  std::ostringstream basinout;
  basinout << "PATHOPT_BASIN_ENERGIES_" << cPtr->mult_struc_counter << ".dat";
  std::ofstream output(basinout.str(), std::ios::app);
  coords::Representation_3D  coord_in, coord_glob, coord_last;
  global_path_minima.resize(N->num_images);
  global_path_minima_temp.resize(N->num_images);
  global_path_minima_energy.resize(N->num_images);

  /**
  * global iterations for multiple runs
  */
  for (size_t t = 0; t < Config::get().neb.GLOBALITERATION; t++) 
  {
	for (size_t i = 0; i < global_path_minima.size(); i++)
    {
		global_path_minima[i].resize(mciteration*(t + 1));
		global_path_minima_temp[i].resize(mciteration*(t + 1));
		global_path_minima_energy[i].resize(mciteration*(t + 1));
    }
    /**
	* initialize coords
	*/
    coord_in = cPtr->xyz();
    coord_glob = cPtr->xyz();
	/** 
	* calculate SP-energies and setting start global minimum
	*/
    MCmin = cPtr->g();
    MCgmin = MCmin;
    /**
	* MCM iterations
	*/
    for (ptrdiff_t mcstep = 0; mcstep <= mciteration; mcstep++)
	{
      coord_last = cPtr->xyz();
	  MCpmin_vec.push_back(MCmin);
	  positions.clear();
      positions.resize(cPtr->size());
	  nanstatus = false;
	  /**
	  * Decision which jump strategy is used
	  * 1. possibility --> MIXED MOVE / rotation of main dihedrals
	  * 2. possibility --> STANDARD MOVE / Cartesian steps
	  */
	  if (Config::get().neb.MIXED_MOVE && (mcstep % 10 == 0 && mcstep > 1))
		  {
			for (size_t j = 0; j < cPtr->size(); j++)
			{
				randvect();
				factor = mcstepsize * (double)rand() / (double)RAND_MAX;
				scon::c3 <double> rv_p{ 3 };
				coords::Cartesian_Point RV;
				double abs = 0.0;
				rv_p.x() = N->tau[opt][j].y() * randvec[2] - N->tau[opt][j].z() * randvec[1];
				rv_p.y() = N->tau[opt][j].z() * randvec[0] - N->tau[opt][j].x() * randvec[2];
				rv_p.z() = N->tau[opt][j].x() * randvec[1] - N->tau[opt][j].y() * randvec[0];
				if (rv_p.x() - rv_p.x() != 0) rv_p.x() = 0.0;
				if (rv_p.y() - rv_p.y() != 0) rv_p.y() = 0.0;
				if (rv_p.z() - rv_p.z() != 0) rv_p.z() = 0.0;
				abs = len(rv_p);
				if (abs != 0.0) normalize(rv_p);
				RV.x() = rv_p.x();
				RV.y() = rv_p.y();
				RV.z() = rv_p.z();
				positions[j] = RV;
			}
			coords::Coordinates temp_perp(*cPtr);
			temp_perp.set_xyz(positions);
			perp_point test_perp;
			test_perp = perp_point(temp_perp.pes());
			move_main(test_perp);
			positions = cPtr->xyz();
	  }
	  else
	  {
			positions = cPtr->xyz();
			for (size_t j = 0; j < cPtr->size(); j++) 
			{
				randvect();
				factor = mcstepsize * (double)rand() / (double)RAND_MAX;
				scon::c3 <double> rv_p{ 3 };
				coords::Cartesian_Point RV;
				double abs = 0.0;
				rv_p.x() = N->tau[opt][j].y() * randvec[2] - N->tau[opt][j].z() * randvec[1];
				rv_p.y() = N->tau[opt][j].z() * randvec[0] - N->tau[opt][j].x() * randvec[2];
				rv_p.z() = N->tau[opt][j].x() * randvec[1] - N->tau[opt][j].y() * randvec[0];
				if (rv_p.x() - rv_p.x() != 0) rv_p.x() = 0.0;
				if (rv_p.y() - rv_p.y() != 0) rv_p.y() = 0.0;
				if (rv_p.z() - rv_p.z() != 0) rv_p.z() = 0.0;
				abs = len(rv_p);
				if (abs != 0.0) normalize(rv_p);
				RV.x() = (factor*rv_p.x());
				RV.y() = (factor*rv_p.y());
				RV.z() = (factor*rv_p.z());
				positions[j] += RV;
			}
      }
      cPtr->set_xyz(positions);
      this->cPtr->NEB_control = false;
	  /**
	  * Optimization using projection or biased gradients
	  */

		  MCmin = lbfgs(opt);

	  /**
	  * MCM Criteria for accepting new minimum
	  */
	  if(Config::get().general.verbosity > 4) std::cout << "MCM energy of step (not proofed) " << mcstep << " is " << MCmin << '\n';
	  if (MCmin != MCmin) 
	  {
		nancounter++;
		nanstatus = true;
		status = 0;
		nbad++;
		cPtr->set_xyz(coord_in);
		if (nancounter>(mcstep/2)) break;
	  }
	  /**
	  * test for low energies that are not reasonable
	  */
	  if (MCmin< (-2.5*std::abs(STARTENERGY))) 
	  {
		MCmin = 1000000.00;
		nbad++;
	  }
	  /**
	  * Metropolis Monte Carlo criterium and test for identical minima
	  */
	  for (auto mp : MCpmin_vec) { if (abs(MCmin - mp) < 0.001 ) status = 2; nbad = 0; }
      /**
	  * testing for identical minima
	  */
	  if (status == 2) 
	  {
		nbad = 0;
	    ///same minimum
		MCpmin_vec[mcstep] = MCmin;
	  }
	 else if (MCmin < MCpmin_vec[mcstep])
	 {
		nbad = 0;
		status = 1; ///accepted as next minimum
		MCpmin_vec[mcstep] = MCmin;
	 }
	 else 
     /**
	 * Metropolis Monte Carlo criterium
	 */
	 {
		nbad = 0;
		boltzman = exp(-_KT_*(MCmin - MCpmin_vec[mcstep]));
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
	 /// test if the displacement is in an acceptable range 
	 if (!testcoord(coord_in)) 
	 {
		status = 0;
	 }
	 ///apply energy range criterium
	 if (std::abs(MCmin - STARTENERGY) > Config::get().neb.PO_ENERGY_RANGE || std::abs(MCmin) == INFINITY) 
	 {
		status = 0;
		nbad++;
	 }
	 ///save energy of next global minimum
	 if ((MCmin < MCgmin) && status == 1) 
	 {
		coord_glob = cPtr->xyz();
		MCgmin = MCmin;
	}
	///restore global minimum after three bad iterations
	if (nbad>3)
	{
		nbad = 0;
		cPtr->set_xyz(coord_glob);
	}
	///restore coords from previous iteration
	else if (status != 1)
	{
		cPtr->set_xyz(coord_last);
	}
	///saving the accepted minima
	else if (status == 1 && nanstatus == false) 
	{
		MCEN = MCmin;
		global_image = opt;
		std::ostringstream struc_opt;
		struc_opt << "PATHOPT_STRUCTURES_" << cPtr->mult_struc_counter << "_" << opt << ".arc";
		output << mcstep << "    " << opt << "    " << std::right << std::fixed << std::setprecision(6) << MCEN << "\n";
		counter++;
		global_path_minima_energy[opt][counter] = MCEN;
		for (size_t i = 0; i<cPtr->size(); i++)
		{
			global_path_minima[opt][counter].push_back(cPtr->xyz(i));
		}
	printmono(struc_opt.str(), global_path_minima[opt][counter], counter);
	}
	else if (status == 2) status = 0;
    }
  }
}

void pathx::MC_PO(ptrdiff_t opt)
{

	double MCmin{ 0.0 }, MCgmin{ 0.0 }, factor{ 0.0 };
	
	/**
	* initialize Boltzman and trial number generation
	*/
	double boltzman{ 0.0 }, trial = (double)rand() / (double)RAND_MAX, start_image_energy{ 0.0 };
	ptrdiff_t nancounter(0), nbad(0), status(0),same_counter(0);
	bool  nanstatus(false);
	global_image = 0;
	counter = 0;
	coords::Representation_3D positions;
	std::ostringstream basinout;
	basinout << "PATHOPT_MC_ENERGIES_" << cPtr->mult_struc_counter << ".dat";
	std::ofstream output(basinout.str(), std::ios::app);
	coords::Representation_3D  coord_in, coord_glob, coord_last;
	global_path_minima.resize(N->num_images);
	global_path_minima_temp.resize(N->num_images);
	global_path_minima_energy.resize(N->num_images);

	/**
	* global iterations for multiple runs
	*/
	for (size_t t = 0; t < Config::get().neb.GLOBALITERATION; t++)
	{
		std::cout << "global iterator: " << t << "\n";
		srand((unsigned int)time(NULL) + pid_func());
		for (size_t i = 0; i < global_path_minima.size(); i++)
		{
			global_path_minima[i].resize(mciteration*(t + 1));
			global_path_minima_temp[i].resize(mciteration*(t + 1));
			global_path_minima_energy[i].resize(mciteration*(t + 1));
		}
		/**
		* initialize coords
		*/
		cPtr->set_xyz(N->imagi[opt]);
		coord_in = cPtr->xyz();
		coord_glob = cPtr->xyz();
		/**
		* calculate SP-energies and setting start global minimum
		*/
		MCmin = cPtr->g();
		MCgmin = MCmin;
		start_image_energy = MCmin;
		std::vector <double> MCpmin_vec;
		MCpmin_vec.clear();
		/**
		* MCM iterations
		*/
		for (ptrdiff_t mcstep = 0; mcstep <= mciteration; mcstep++)
		{
			coord_last = cPtr->xyz();
			MCpmin_vec.push_back(MCmin);
			positions.clear();
			positions.resize(cPtr->size());
			nanstatus = false;
			status = 0;
			if (same_counter > 10) same_counter = 0;
			/**
			* Decision which jump strategy is used
			* 1. possibility --> MIXED MOVE / rotation of main dihedrals
			* 2. possibility --> STANDARD MOVE / Cartesian steps
			*/
			if (Config::get().neb.MIXED_MOVE && (mcstep % 10 == 0 && mcstep > 1))
			{
				for (size_t j = 0; j < cPtr->size(); j++)
				{
					randvect();
					factor = mcstepsize * (double)rand() / (double)RAND_MAX;
					scon::c3 <double> rv_p{ 3 };
					coords::Cartesian_Point RV;
					double abs = 0.0;
					rv_p.x() = N->tau[opt][j].y() * randvec[2] - N->tau[opt][j].z() * randvec[1];
					rv_p.y() = N->tau[opt][j].z() * randvec[0] - N->tau[opt][j].x() * randvec[2];
					rv_p.z() = N->tau[opt][j].x() * randvec[1] - N->tau[opt][j].y() * randvec[0];
					if (rv_p.x() - rv_p.x() != 0) rv_p.x() = 0.0;
					if (rv_p.y() - rv_p.y() != 0) rv_p.y() = 0.0;
					if (rv_p.z() - rv_p.z() != 0) rv_p.z() = 0.0;
					abs = len(rv_p);
					if (abs != 0.0) normalize(rv_p);
					RV.x() = rv_p.x();
					RV.y() = rv_p.y();
					RV.z() = rv_p.z();
					positions[j] = RV;
				}
				coords::Coordinates temp_perp(*cPtr);
				temp_perp.set_xyz(positions);
				perp_point test_perp;
				test_perp = perp_point(temp_perp.pes());
				move_main(test_perp);
				positions = cPtr->xyz();
			}
			else
			{
				positions = cPtr->xyz();
				for (size_t j = 0; j < cPtr->size(); j++)
				{
					randvect();
					factor = mcstepsize * (double)rand() / (double)RAND_MAX;
					scon::c3 <double> rv_p{ 3 };
					coords::Cartesian_Point RV;
					double abs = 0.0;
					rv_p.x() = N->tau[opt][j].y() * randvec[2] - N->tau[opt][j].z() * randvec[1];
					rv_p.y() = N->tau[opt][j].z() * randvec[0] - N->tau[opt][j].x() * randvec[2];
					rv_p.z() = N->tau[opt][j].x() * randvec[1] - N->tau[opt][j].y() * randvec[0];
					if (rv_p.x() - rv_p.x() != 0) rv_p.x() = 0.0;
					if (rv_p.y() - rv_p.y() != 0) rv_p.y() = 0.0;
					if (rv_p.z() - rv_p.z() != 0) rv_p.z() = 0.0;
					abs = len(rv_p);
					if (abs != 0.0) normalize(rv_p);
					RV.x() = (factor*rv_p.x());
					RV.y() = (factor*rv_p.y());
					RV.z() = (factor*rv_p.z());
					positions[j] += RV;
				}
			}
			cPtr->set_xyz(positions);
			this->cPtr->NEB_control = false;
			/**
			* Optimization using projection or biased gradients
			*/

			MCmin = cPtr->g();
			/**
			* MCM Criteria for accepting new minimum
			*/
			if (Config::get().general.verbosity > 5) std::cout << "MCM energy of step (not proofed) " << mcstep << " is " << MCmin << '\n';
			if (MCmin != MCmin)
			{
				nancounter++;
				nanstatus = true;
				status = 0;
				nbad++;
				cPtr->set_xyz(coord_in);
				if (nancounter > (mcstep / 2)) break;
			}

			/**
			* Metropolis Monte Carlo criterium and test for identical minima
			*/
			for (auto mp : MCpmin_vec) { if (abs(MCmin - mp) < 0.00001) status = 2; nbad = 0; }
			/**
			* testing for identical minima
			*/
			if (status == 2)
			{
				nbad++;
				///same minimum
				same_counter++;
				if (Config::get().general.verbosity > 4)
				{
					std::cout << "same Minimum \n";
				}
			}
			else if (MCmin < MCpmin_vec[mcstep])
			{
				nbad = 0;
				status = 1; ///accepted as next minimum
				MCpmin_vec[mcstep] = MCmin;
				coord_glob = cPtr->xyz();
				if (Config::get().general.verbosity > 4)
				{
					std::cout << "Accepted due to lower energy criterium \n";
				}
			}
			else
			{
				nbad = 0;
				boltzman = exp(-_KT_*(MCmin - start_image_energy));
				trial = (double)rand() / (double)RAND_MAX;
				if (boltzman < trial)
				{
					status = 0;
					if (Config::get().general.verbosity > 4)
					{
						std::cout << "Rejected due to Metropolis criterium \n";
					}
					nbad++;
				}
				else
				{
					status = 1;
					MCpmin_vec[mcstep] = MCmin;
					nbad = 0;
					if (Config::get().general.verbosity > 4)
					{
						std::cout << "Accepted due to Metropolis criterium \n";
					}
				}
			}
			/// test if the displacement is in an acceptable range 
			if (!testcoord(coord_in))
			{
				status = 0;
			}
			///apply energy range criterium
			if (std::abs(MCmin - STARTENERGY) > Config::get().neb.PO_ENERGY_RANGE || std::abs(MCmin) == INFINITY)
			{
				status = 0;
				nbad++;
				if (Config::get().general.verbosity > 4)
				{
					std::cout << "Out of Energy Range: " << MCmin << "\n";
				}
			}
			///save energy of next global minimum
			if ((MCmin < MCgmin) && status == 1)
			{
				coord_glob = cPtr->xyz();
				MCgmin = MCmin;
				if (Config::get().general.verbosity > 4)
				{
					std::cout << "New global minimum: " << MCgmin << "\n";
				}
			}
			///restore global minimum after five bad iterations
			if (nbad > 3 || same_counter >= 10)
			{
				nbad = 0;
				cPtr->set_xyz(coord_glob);
				if (Config::get().general.verbosity > 4)
				{
					std::cout << "Restoring old coords due to bad_iterator \n";
				}
			}
			///restore coords from previous iteration
			else if (status != 1)
			{
				cPtr->set_xyz(coord_last);
				if (Config::get().general.verbosity > 4)
				{
					std::cout << "Restoring old coords \n";
				}
			}
			///saving the accepted minima
			else if (status == 1 && nanstatus == false)
			{
				
				if (Config::get().general.verbosity > 4)
				{
					std::cout << "structure is saved \n";
				}
				MCEN = MCmin;
				global_image = opt;
				std::ostringstream struc_opt;
				struc_opt << "PATHOPT_STRUCTURES_" << cPtr->mult_struc_counter << "_" << opt << ".arc";
				output << mcstep << "    " << opt << "    " << std::right << std::fixed << std::setprecision(6) << MCEN << "\n";
				counter++;
				global_path_minima_energy[opt][counter] = MCEN;
				for (size_t i = 0; i < cPtr->size(); i++)
				{
					global_path_minima[opt][counter].push_back(cPtr->xyz(i));
				}
				if (counter % Config::get().neb.MCM_SAVEITER == 0 && mcstep > 1)
				{
					printmono(struc_opt.str(), global_path_minima[opt][counter], counter);
				}
			}
			else if (status == 2) status = 0;
		}
	}
}

/*
* CONNECTION PROCEDURE FOR THE FOUND MINIMA ON N-HYPERPLANES
*/

void pathx::proof_connect()
{

	ptrdiff_t temp_image{ ptrdiff_t(N->num_images) };
	ptrdiff_t image{ ptrdiff_t(Config::get().neb.CONNECT_NEB_NUMBER) };
	std::vector < std::vector <std::vector  <double> > > RMSD{ N->num_images };
	std::vector < std::vector <ptrdiff_t> >PARTNER{ N->num_images };
	std::vector < std::vector <ptrdiff_t> >  TABU{ N->num_images };
	coords::Cartesian_Point cog1, cog2, tempcoord1, tempcoord2;
	coords::Coordinates  coord1, coord2;
	for (size_t l = 0; l < RMSD.size(); l++)RMSD[l].resize(Config::get().neb.MCITERATION*(Config::get().neb.GLOBALITERATION));
    for (size_t l = 0; l < PARTNER.size(); l++)PARTNER[l].resize(Config::get().neb.MCITERATION*(Config::get().neb.GLOBALITERATION));
    for (size_t l = 0; l < TABU.size(); l++)TABU[l].resize(Config::get().neb.MCITERATION*(Config::get().neb.GLOBALITERATION));
	bool reverse(false);
	coords::Representation_3D tempstart{ N->imagi[0] }, tempstart2{ N->imagi[N->num_images - 1] };
	coord1.init_in(cPtr->atoms(), coords::PES_Point(N->imagi[0]), true);
	coord2.init_in(cPtr->atoms(), N->imagi[0], true);
	/// center of geometry 
	for (ptrdiff_t i = 1; i < temp_image - 1; i++)
	{
		N->num_images = temp_image;
		cog2 = coord2.center_of_geometry();
		for (size_t j = 0; j < global_path_minima[i].size(); j++)
		{
			if (global_path_minima[i][j].size() == 0) continue;
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
	 ///caclulating the rmsd value starting from hyperplane n to n+1
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
	 /// loop over the first next up to the third nearest neighbors
	 size_t arrhenius_counter(0U);
	 double arrhenius(0.0);
	 std::ofstream arrhenius_file("arrhenius_global.dat", std::ios::app);
	 for (size_t mm = 1; mm < 4; mm++) 
	 {
		std::ostringstream NEB1;
		NEB1 << mm;

#if defined(_MSC_VER)
    _mkdir(NEB1.str().c_str());
    _chdir(NEB1.str().c_str());
#else
    mkdir(NEB1.str().c_str(), 0777);
    int res = chdir(NEB1.str().c_str());
	if (res != 0) std::cout<<"Something went wrong!\n";
#endif




	    /// assigning the partners with respect to the nearest neighbors on the next hyperplane
		for (ptrdiff_t i = 1; i < temp_image - 1; i++)	
		{
			for (size_t j = 0; j < global_path_minima[i].size(); j++)
			{
				if (RMSD[i][j].empty()) continue;
				size_t index = 0;
				double min = abs(RMSD[i][j][0]);
				for (size_t p = 0; p < RMSD[i][j].size(); p++)
				{
					if (abs(RMSD[i][j][p]) < min)
					{
						index = p;
						min = abs(RMSD[i][j][p]);
					}
				}
				PARTNER[i][j] = (index + (mm - 1));
				RMSD[i][j].erase(RMSD[i][j].begin() + index);
			}	
		}
		ptrdiff_t i(0U), jj(0U);
		ptrdiff_t tempcount (0U);
		///Connecting I/O 1: VIA NEB 2: DIRECT 
		if (Config::get().neb.NEB_CONN == true)
		{
			for (ptrdiff_t j = 0; j < ptrdiff_t(global_path_minima[1].size()); j++)
			{
				if (global_path_minima[1][j].empty()) continue;
				reverse = false;
				N->preprocess(j, image, j, global_path_minima[1][j], tempstart, reverse);
				reverse = true;
				tempcount = 0;
				jj = j;
				for (i = 1; i < temp_image - 2; i++)
				{
					if (global_path_minima[i - tempcount][jj].empty()) continue;
					if (global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]].empty()) continue;
					N->preprocess(j, image, j, global_path_minima[i - tempcount][jj], global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]], reverse);
					jj = PARTNER[i - tempcount][jj];
				}       
				reverse = true;
				if (global_path_minima[i - tempcount][jj].size() == 0)continue;
				N->preprocess(j, image, j, global_path_minima[i - tempcount][jj], tempstart2, reverse);
			}
			}
			else
			{
				for (ptrdiff_t j = 0; j < ptrdiff_t(global_path_minima[1].size()); j++)
				{
					std::vector <double> energy_connect;
					if (global_path_minima[1][j].empty()) continue;
					reverse = false;
					std::ostringstream en, img;
					en << "ENERGY_" <<cPtr->mult_struc_counter<<"_"<< j << ".dat";
					img << "PATH_" << cPtr->mult_struc_counter << "_"<< j << ".arc";
					std::ofstream energy(en.str().c_str(), std::ios::app);
					printmono(img.str().c_str(), tempstart, j);
					cPtr->set_xyz(tempstart);
					energy_connect.push_back(cPtr->g());
					energy << std::right << std::fixed << std::setprecision(6) << energy_connect[0] << '\n';
					printmono(img.str().c_str(), global_path_minima[1][j], j);
					cPtr->set_xyz(global_path_minima[1][j]);
					energy_connect.push_back(cPtr->g());
					energy << std::right << std::fixed << std::setprecision(6) << energy_connect[1] << '\n';
					tempcount = 0;
					jj = j;

					//During the path sempling each path is written in an array and printed out later
					std::vector<coords::Representation_3D>* final_path = nullptr;
					if (Config::get().neb.MAXFLUX_PATHOPT) {
						final_path = new std::vector<coords::Representation_3D>;
						final_path->resize(temp_image);
						final_path->at(0) = tempstart;
						final_path->at(1) = global_path_minima[1][j];
					}

					for (i = 1; i < temp_image - 2; i++)
					{
						if (global_path_minima[i - tempcount][jj].empty()) continue;
						if (global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]].empty()) continue;
						printmono(img.str().c_str(), global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]], j);
						
						if (Config::get().neb.MAXFLUX_PATHOPT) {
							final_path->at(i + 1) = global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]];
						}
						
						cPtr->set_xyz(global_path_minima[(i - tempcount) + 1][PARTNER[i - tempcount][jj]]);
						energy_connect.push_back(cPtr->g());
						energy << std::right << std::fixed << std::setprecision(6) << energy_connect[i+1] << '\n';
						jj = PARTNER[i - tempcount][jj];
					}
					printmono(img.str().c_str(), tempstart2, j);

					if (Config::get().neb.MAXFLUX_PATHOPT) {
						auto kill = false;
						final_path->at(final_path->size() - 1U) = tempstart2;
						for (auto && path : *(final_path)) {
							if (path.empty()) kill = true;
						}
						if (!kill) {
							N->preprocess(*(final_path), j);
						}
						/*std::for_each(final_path->begin(), final_path->end(), [&](auto && path) {
						printmono(std::string("MAXFLUX_PATH_" + std::to_string(j) + ".arc").c_str(), path, j);
						});*/
						delete final_path;
						final_path = nullptr;
					}

					cPtr->set_xyz(tempstart2);
					energy_connect.push_back(cPtr->g());
					energy << std::right << std::fixed << std::setprecision(6) << energy_connect.back() << '\n';
					energy << "Barrier_rel to start: "<< *max_element(energy_connect.begin(), energy_connect.end()) - energy_connect[0]<<'\n';
					energy << "RATE via Arrhenius: " << exp(-((*max_element(energy_connect.begin(), energy_connect.end()) - energy_connect[0]) / _KT_)) << '\n';
					arrhenius += exp(-((*max_element(energy_connect.begin(), energy_connect.end()) - energy_connect[0]) / _KT_));
					arrhenius_counter++;
					arrhenius_file << exp(-((*max_element(energy_connect.begin(), energy_connect.end()) - energy_connect[0]) / _KT_)) << '\n';
				}
			}
#if defined (_MSC_VER)
    _chdir("..");
#else
    res = chdir("..");
	if (res != 0) std::cout<<"Something went wrong!\n";
#endif
  }
  arrhenius_file << arrhenius/arrhenius_counter << '\n';
}


/**
 *GENERATES 3D-UNIT VECTOR
 *G: MARSAGLIA, ANN. MATH. STAT., 43, 645 (1972)
 */
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

/** 
* FUNCTION FOR TESTING THE COORDINATE VARIATION WITH RESPECT TO MAXVARIATION LIMIT 
*/
bool pathx::testcoord(coords::Representation_3D &coords)
{

  for (size_t i = 0; i<cPtr->size(); i++) {

    if (abs(cPtr->xyz(i).x() - coords[i].x()) > maxvar) return false;
    if (abs(cPtr->xyz(i).y() - coords[i].y()) > maxvar) return false;
    if (abs(cPtr->xyz(i).z() - coords[i].z()) > maxvar) return false;

  }
  return true;
}

/** 
* LBFGS WITH PATHOPT CALLBACK FOR GRADIENTS 
*/
double pathx::lbfgs(ptrdiff_t imagex)
{

	using namespace  optimization::local;

	/**
	* Create optimizer
	*/
	global_imagex = imagex;
	coords::Cartesian_Point g;
	auto optimizer = make_lbfgs(
		make_more_thuente(GradCallBack(*this))
	);
	/** 
	* Create Point 
	*/
	using op_type = decltype(optimizer);
	op_type::point_type x(scon::explicit_transform<op_type::rep_type>(cPtr->xyz()));
	/** 
	* Optimize point 
	*/
	optimizer.config.max_iterations = Config::get().optimization.local.bfgs.maxstep;
	optimizer.config.epsilon = (float)Config::get().optimization.local.bfgs.grad;
	optimizer(x);

	if (Config::get().general.verbosity > 3)
	{
		std::cout << "Optimization done (status " << optimizer.state() << "). Evaluations:" << optimizer.iter() << '\n';
	}

  return optimizer.p().f;
}


/** 
* MODIFICATION OF GRADIENT CALCULATION 
*/
double pathx::g_new(ptrdiff_t im)
{
  double cosi(0.0), energytemp(0.0), pot(0.0), dpot(0.0);
  coords::Cartesian_Point rv_p, rv_n;
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

  }	  
  return energytemp;
}

/**
* PRINT FCT FOR SINGLE STRUCTURE
*/
void pathx::printmono(std::string const &name, coords::Representation_3D &print, ptrdiff_t &count)
{

  std::ofstream out(name.c_str(), std::ios::app);
  std::string temp;
  out << "     " << cPtr->size() << "  global counter:  " << count << '\n';
  for (size_t j = 0; j < cPtr->size(); j++) {

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

/*
* FUNCTION FOR CHANGING RANDOM DIHEDRAL ALONG GIVEN DIRECTION
*/

void pathx::move_main(perp_point & direction)
{
  size_t const NUM(cPtr->main().size());
  coords::Representation_Main main_tors(NUM), bla(NUM);
  /// choose number of torsions to modify
  coords::float_type const NUM_MAINS(static_cast<coords::float_type>(NUM));
  size_t const NUM_MOD(std::min((static_cast<size_t>(-std::log(scon::rand<coords::float_type>(0.0, 1.0))) + 1U), NUM));
  /// apply those torsions
  if (Config::get().general.verbosity > 4) std::cout << "Changing " << NUM_MOD << " of " << NUM << " mains.\n";
  for (size_t i(0U); i < NUM_MOD; ++i)
  {
    size_t const K(static_cast<size_t>(NUM_MAINS*scon::rand<coords::float_type>(0.0, 1.0)));
    coords::float_type const F = scon::rand<coords::float_type>(-Config::get().optimization.global.montecarlo.dihedral_max_rot,
      Config::get().optimization.global.montecarlo.dihedral_max_rot);
    if (Config::get().general.verbosity > 4) std::cout << "Changing main " << K << " by " << F << '\n';
    main_tors[K] = coords::main_type::from_deg(F);
  }
  /// orthogonalize movement to main directions
  for (auto perp_direction : direction.main_direction)
  {
    perp_direction.resize(main_tors.size());
    scon::orthogonalize_to_normal(main_tors, perp_direction);
  }
  for (size_t i = 0; i < NUM; ++i)
  {
    cPtr->rotate_main(i, main_tors[i]);
  }
  /// main_tors.print(cout);
  cPtr->to_xyz();
}