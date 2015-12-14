#include "neb.h"
#include <vector>
#include <string>
#include <cmath>
#include "scon_vect.h"
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "ls.h"
#include "lbfgs.h"
#include "coords.h"
#include"optimization_global.h"




::tinker::parameter::parameters neb::tp;


neb::neb(coords::Coordinates *cptr)
{
	cPtr = cptr;
	reversed = false;


	if (!tp.valid()) tp.from_file(Config::get().get().general.paramFilename);
	std::vector<size_t> types;
	for (auto atom : (*cptr).atoms()) scon::sorted::insert_unique(types, atom.energy_type());
	cparams = tp.contract(types);
	refined.refine(*cptr, cparams);
}

neb::~neb(void)
{
}
void neb::preprocess (ptrdiff_t &count, ptrdiff_t const &image)
{
  std::vector<size_t> image_remember;
  std::vector<std::vector<int> > atoms_remember;
  N = cPtr->size();
  CIMaximum=0;
  num_images=image;
  ts=false;
  ClimbingImage=true;
  reversed=true;
  imagi.resize(num_images+N);
  image_ini.resize(num_images+N);
  imagederiv.resize(num_images+N);
  tau.resize(num_images+N);
  images_initial.resize(N);
  images.resize(N);
  tempimage_final.resize(num_images+N);
  tempimage_ini.resize(num_images+N);
  energies.resize(num_images+N);
  springconstant = Config::get().neb.SPRINGCONSTANT;
  initial();
  final();
  create(count);
  if(Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
  else run(count);
  controlrun =true;
}

void neb::preprocess(ptrdiff_t &image, ptrdiff_t &count, const coords::Representation_3D &start, const coords::Representation_3D &fi, const std::vector <double> &ts_energy, const std::vector <double> &min_energy,bool reverse,const coords::Representation_3D &ts_path)
{


  std::vector<size_t> image_remember;
  std::vector<std::vector<int> > atoms_remember;
  N = cPtr->size();
  CIMaximum=0;
  ts=true;
  reversed=reverse;
  num_images=image;
  imagi.clear();
  image_ini.clear();
  imagederiv.clear();
  energies.clear();
  tau.clear();
  images.clear();
  ClimbingImage=true;
  ts_energies.resize(ts_energy.size());
  min_energies.resize(min_energy.size());
  imagi.resize(start.size()+num_images);
  ts_pathstruc.resize(start.size()+num_images);
  image_ini.resize(start.size()+num_images);
  imagederiv.resize(start.size()+num_images);
  images_initial.resize(N);
  images.resize(N);
  tempimage_final.resize(N+num_images);
  tempimage_ini.resize(N+num_images);
  energies.resize(start.size()+num_images);
  ts_energies=ts_energy;
  min_energies=min_energy;
  ts_pathstruc=ts_path;
  springconstant = Config::get().neb.SPRINGCONSTANT;
  initial(start);
  final(fi);
  create(count);
  if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
  else run(count);
  controlrun =true;





}

void neb::preprocess(ptrdiff_t &file, ptrdiff_t &image,ptrdiff_t &count, const coords::Representation_3D &start, const coords::Representation_3D &fi, bool reverse)
{


  std::vector<size_t> image_remember;
  std::vector<std::vector<int> > atoms_remember;
  N = cPtr->size();
  this->cPtr->mult_struc_counter=file;
  CIMaximum=0;
  ts=false;
  reversed=reverse;
  num_images=image;
  imagi.clear();
  image_ini.clear();
  imagederiv.clear();
  energies.clear();
  tau.clear();
  images.clear();
  ClimbingImage=true;
 
  imagi.resize(start.size()+num_images);
  ts_pathstruc.resize(start.size()+num_images);
  image_ini.resize(start.size()+num_images);
  imagederiv.resize(start.size()+num_images);
  images_initial.resize(N);
  images.resize(N);
  tempimage_final.resize(N+num_images);
  tempimage_ini.resize(N+num_images);
  energies.resize(start.size()+num_images);
  springconstant = Config::get().neb.SPRINGCONSTANT;
  initial(start);
  final(fi);
  create(count);
  if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
  else run(count);
  controlrun =true;





}


void neb::initial(void)

{

  images_initial = cPtr->xyz();
  imagi[0] = cPtr->xyz();
  std::ifstream final (Config::get().neb.START_STRUCTURE.c_str());
  std::string buffer;
  //getline(final,buffer);
  getline(final,buffer);
  size_t number;
  char atom[2];
  
  /*for (size_t i=0;i<N;i++)
  {
    getline(final,buffer);
    std::istringstream line_coord(buffer);
    line_coord >>number>>atom >> images[i].x() >> images[i].y() >> images[i].z() ;
    imagi[0].push_back(images[i]);
	
  }*/

  images_initial = imagi[0];
 
}

void neb::final(void)
{
  std::ifstream final (Config::get().neb.FINAL_STRUCTURE.c_str());
  std::string buffer;
  getline(final,buffer);
  //getline(final,buffer);
  size_t number;
  char atom[2];
  coords::Cartesian_Point cog1,cog2,tempcoord;

  for (size_t i=0;i<N;i++)
  {
    getline(final,buffer);
    std::istringstream line_coord(buffer);
    line_coord >>number>>atom >> images[i].x() >> images[i].y() >> images[i].z() ;
    imagi[num_images-1].push_back(images[i]);	
  }
 
}

void neb::initial(const coords::Representation_3D &start)

{


  for(size_t i=0;i<cPtr->size();i++)
  {
    imagi[0].push_back(start[i]);
  }

}

void neb::final(const coords::Representation_3D &fi)
{
  for (size_t i=0;i<cPtr->size();i++)
  {
    imagi[num_images-1].push_back(fi[i]);
  }
}
void neb::create(ptrdiff_t &count)
{

  double diff;
  tempimage_ini = imagi[0];
  tempimage_final = imagi[num_images-1];
  //std::string name;
  std::ostringstream name,name2,name3;
  // off << tempimage_final << endl;
  name <<"IMAGES_INI"<<cPtr->mult_struc_counter<<".xyz";



  cPtr->set_xyz(imagi[0]);
  cPtr->to_internal();
  tempimage_ini=cPtr->xyz();

  cPtr->set_xyz(imagi[num_images-1]);
  cPtr->to_internal();
  tempimage_final=cPtr->xyz();

  for (size_t j=1; j<(num_images-1);j++){

    diff=(double)j/num_images;

    for(size_t i=0; i<this->cPtr->size();i++){

      images[i].x()= tempimage_ini[i].x() + diff * (tempimage_final[i].x()- tempimage_ini[i].x());
      images[i].y()= tempimage_ini[i].y() + diff * (tempimage_final[i].y()- tempimage_ini[i].y());
      images[i].z()= tempimage_ini[i].z() + diff * (tempimage_final[i].z()- tempimage_ini[i].z());


      imagi[j].push_back(images[i]);
      image_ini[j].push_back(images[i]);

    }	
	cPtr->set_xyz(imagi[j]);

	for (auto const & bond : refined.bonds())
	{
		coords::Cartesian_Point const bv(cPtr->xyz(bond.atoms[0]) - cPtr->xyz(bond.atoms[1]));
		coords::float_type const d = len(bv);
		//std::cout << "dist " << d << " ideal " << bond.ideal << lineend;
		if (d < bond.ideal && (1.0 - d / bond.ideal) > 0.10)   std::cout << "WARNING_1-bond_length not reasonable" << lineend;
		if (d > bond.ideal && abs(1.0 - (d / bond.ideal)) > 0.10) std::cout <<  "WARNING_1-bond_length not reasonable" << lineend;




	}

	for (auto const & angle : refined.angles())
	{
		coords::Cartesian_Point const
			av1(cPtr->xyz(angle.atoms[0]) - cPtr->xyz(angle.atoms[1])),
			av2(cPtr->xyz(angle.atoms[2]) - cPtr->xyz(angle.atoms[1]));
		coords::float_type const d(scon::angle(av1, av2).degrees());
		//std::cout << "angle " << d << " ideal " << angle.ideal << lineend;
		if (d < angle.ideal && (1.0 - d / angle.ideal) > 0.10) std::cout << "WARNING_1_dihedral not reasonable" << lineend;
		if (d > angle.ideal && abs(1.0 - (d / angle.ideal)) > 0.10) std::cout << "WARNING_2_dihedral not reasonable" << lineend;
	}
  }
  
  if(Config::get().neb.INT_PATH) calc_shift();
  print(name.str(),imagi,count);

}





void neb::run (ptrdiff_t &count)
{

  get_energies_grad();
  calc_tau();
  opt_mep(count);


}
void neb::run(ptrdiff_t &count, std::vector<size_t>& , std::vector<std::vector<int> >& atoms_remember)
{

	get_energies_grad();
	calc_tau();
	opt_mep(count);

	internal_execute(imagi, atoms_remember);

	//no_torsion(image_remember[DIHEDRAL],atoms_remember[DIHEDRAL]);

	opt_internals(count, atoms_remember);

	//internal_execute(imagi,image_remember,atoms_remember);

}



void neb::get_energies_grad (void)
{

  for(size_t i = 0; i<num_images;i++)
  {

    cPtr->set_xyz(imagi[i]);

    energies[i] = cPtr->g();
    if(i>0){
      if (energies[i] > energies[i-1]) CIMaximum = i;}
    for (size_t k=0 ; k < N;k++){

      images[k] = cPtr->g_xyz(k);


      imagederiv[i].push_back(images[k]);
    }



  }

  // cout << "maximum energy image is: " << CIMaximum+1 << endl;
}

void neb::calc_tau (void)
{
    double Em1,Ep1,Emax,Emin,abso(0.0);
	/*ofstream off ("energies");*/


  for(size_t i =1; i<num_images-1;i++){
	 
	  if (Config::get().neb.TAU == false)
	  {
		  for (size_t j = 0; j < N; j++) {
			  images[j].x() = imagi[i][j].x() - imagi[i - 1][j].x() / abs(imagi[i][j].x() - imagi[i - 1][j].x()) + (imagi[i + 1][j].x() - imagi[i][j].x()) / abs(imagi[i + 1][j].x() - imagi[i][j].x());
			  images[j].y() = imagi[i][j].y() - imagi[i - 1][j].y() / abs(imagi[i][j].y() - imagi[i - 1][j].y()) + (imagi[i + 1][j].y() - imagi[i][j].y()) / abs(imagi[i + 1][j].y() - imagi[i][j].y());
			  images[j].z() = imagi[i][j].z() - imagi[i - 1][j].z() / abs(imagi[i][j].z() - imagi[i - 1][j].z()) + (imagi[i + 1][j].z() - imagi[i][j].z()) / abs(imagi[i + 1][j].z() - imagi[i][j].z());
			  abso = 0.0;
			  tau[i].push_back(images[j]);
			  abso += tau[i][j].x()*tau[i][j].x();
			  abso += tau[i][j].y()*tau[i][j].y();
			  abso += tau[i][j].z()*tau[i][j].z();
			  abso = sqrt(abso);
			  if (abso == 0.0) abso = 1.0;
			  tau[i][j].x() = tau[i][j].x() / abso;
			  tau[i][j].y() = tau[i][j].y() / abso;
			  tau[i][j].z() = tau[i][j].z() / abso;
		  }
	  }
	  else
	  {

		  EnergyPml = 0.0;
		  EnergyPpl = 0.0;

		  if (energies[i - 1] > energies[i]) EnergyPml = energies[i - 1];
		  else EnergyPml = energies[i];
		  if (energies[i + 1] > energies[i]) EnergyPpl = energies[i + 1];
		  else EnergyPpl = energies[i];

		  if (EnergyPml != EnergyPml) {
			  if (EnergyPml > EnergyPpl) {

				  for (size_t j = 0; j < N; j++) {

					  abso = 0.0;
					  images[j].x() = (imagi[i][j].x() - imagi[i - 1][j].x());
					  images[j].y() = (imagi[i][j].y() - imagi[i - 1][j].y());
					  images[j].z() = (imagi[i][j].z() - imagi[i - 1][j].z());
					  tau[i].push_back(images[j]);
					  abso += tau[i][j].x()*tau[i][j].x();
					  abso += tau[i][j].y()*tau[i][j].y();
					  abso += tau[i][j].z()*tau[i][j].z();
					  abso = sqrt(abso);
					  if (abso == 0.0) abso = 1.0;
					  tau[i][j].x() = tau[i][j].x() / abso;
					  tau[i][j].y() = tau[i][j].y() / abso;
					  tau[i][j].z() = tau[i][j].z() / abso;

				  }
			  }
			  else {
				  for (size_t j = 0; j < N; j++) {

					  abso = 0.0;

					  images[j].x() = (imagi[i + 1][j].x() - imagi[i][j].x());
					  images[j].y() = (imagi[i + 1][j].y() - imagi[i][j].y());
					  images[j].z() = (imagi[i + 1][j].z() - imagi[i][j].z());
					  tau[i].push_back(images[j]);
					  //tau[i].norm();
					  abso += tau[i][j].x()*tau[i][j].x();
					  abso += tau[i][j].y()*tau[i][j].y();
					  abso += tau[i][j].z()*tau[i][j].z();
					  abso = sqrt(abso);
					  if (abso == 0.0) abso = 1.0;
					  tau[i][j].x() = tau[i][j].x() / abso;
					  tau[i][j].y() = tau[i][j].y() / abso;
					  tau[i][j].z() = tau[i][j].z() / abso;

				  }
			  }
		  }
		  else {
			  Em1 = energies[i - 1] - energies[i];
			  Ep1 = energies[i + 1] - energies[i];

			  Emin = std::min(abs(Ep1), abs(Em1));
			  Emax = std::max(abs(Ep1), abs(Em1));

			  if (Em1 > Ep1) {
				  for (size_t j = 0; j < N; j++) {
					  abso = 0.0;
					  images[j].x() = (imagi[i + 1][j].x() - imagi[i][j].x()) * Emin + (imagi[i][j].x() - imagi[i - 1][j].x()) * Emax;
					  images[j].y() = (imagi[i + 1][j].y() - imagi[i][j].y()) * Emin + (imagi[i][j].y() - imagi[i - 1][j].y()) * Emax;
					  images[j].z() = (imagi[i + 1][j].z() - imagi[i][j].z()) * Emin + (imagi[i][j].z() - imagi[i - 1][j].z()) * Emax;
					  tau[i].push_back(images[j]);
					  //tau[i].norm();
					  abso += tau[i][j].x()*tau[i][j].x();
					  abso += tau[i][j].y()*tau[i][j].y();
					  abso += tau[i][j].z()*tau[i][j].z();
					  abso = sqrt(abso);
					  if (abso == 0.0) abso = 1.0;

					  tau[i][j].x() = tau[i][j].x() / abso;
					  tau[i][j].y() = tau[i][j].y() / abso;
					  tau[i][j].z() = tau[i][j].z() / abso;


				  }

			  }


			  else {
				  for (size_t j = 0; j < N; j++) {

					  abso = 0.0;
					  images[j].x() = (imagi[i + 1][j].x() - imagi[i][j].x()) * Emax + (imagi[i][j].x() - imagi[i - 1][j].x()) * Emin;
					  images[j].y() = (imagi[i + 1][j].y() - imagi[i][j].y()) * Emax + (imagi[i][j].y() - imagi[i - 1][j].y()) * Emin;
					  images[j].z() = (imagi[i + 1][j].z() - imagi[i][j].z()) * Emax + (imagi[i][j].z() - imagi[i - 1][j].z()) * Emin;
					  tau[i].push_back(images[j]);
					  abso += tau[i][j].x()*tau[i][j].x();
					  abso += tau[i][j].y()*tau[i][j].y();
					  abso += tau[i][j].z()*tau[i][j].z();
					  abso = sqrt(abso);
					  if (abso == 0.0) abso = 1.0;
					  tau[i][j].x() = tau[i][j].x() / abso;
					  tau[i][j].y() = tau[i][j].y() / abso;
					  tau[i][j].z() = tau[i][j].z() / abso;
				  }
			  }

		  }
	  }
	
  }
}

void neb::opt_mep(ptrdiff_t &count)
{



	

	std::ostringstream energies, name, energies_ini;
	energies << "ENERGIES_COMPLETE_" << this->cPtr->mult_struc_counter;
	energies_ini << "ENERGIES_COMPLETE_LIN_" << this->cPtr->mult_struc_counter;
	std::string temp;
	temp = energies.str();

	std::fstream off(temp.c_str(), std::ios::app);
	temp = energies_ini.str();
	std::fstream off2(temp.c_str(), std::ios::app);
	energies_NEB.resize(num_images);
	grad_v = 1.0;
	grad_v_temp = 2.0;
	int maxit(0U);
	energies_NEB[0] = this->energies[0];
	energies_NEB[num_images] = this->energies[num_images];

	while (std::abs(grad_v - grad_v_temp) > Config::get().neb.NEB_RMSD && maxit < Config::get().neb.NEB_INT_IT) {
		
		grad_v_temp = grad_v;
		maxit++;
		calc_tau();

		if (reversed == true)
		{
		

			for (size_t imagecount = 1; imagecount < num_images - 1; imagecount++)
			{

				
				cPtr->set_xyz(imagi[imagecount]);
			
				energies_NEB[imagecount]=lbfgs(imagecount);
				
				imagi[imagecount] = cPtr->xyz();

				grad_v += grad_v;

			}



		}
		else{


			for (ptrdiff_t imagecount = num_images - 2; imagecount >= 1; imagecount--)
			{


				

				cPtr->set_xyz(imagi[imagecount]);

				energies_NEB[imagecount] = lbfgs(imagecount);

				imagi[imagecount] = cPtr->xyz();

				grad_v += grad_v;


			}

		

		
		}

	
		grad_v /= (num_images - 2);

	}



	if (reversed == true)
	{
		name << "IMAGES_FINAL" << this->cPtr->mult_struc_counter << ".xyz";
		if (ts == true)
		{
			off << "TS          " << ts_energies[count] << lineend;
			off << "MIN         " << min_energies[count] << lineend;
			printmono(name.str(), imagi[0], count);
		}
		off << "ENERGIE:    " << this->energies[0] << "   START" << lineend;
		off2 << "ENERGIE:    " << this->energies[0] << "   START" << lineend;
	
		for (size_t imagecount = 1; imagecount < num_images - 1; imagecount++)
		{

			cPtr->set_xyz(imagi[imagecount]);





			off << "ENERGIE:    " << energies_NEB[imagecount]  << "     IMAGE:   " << imagecount << "   GLOBAL-COUNT:  " << count << lineend;
			off2 << "ENERGIE:    " << this->energies[imagecount] << "     IMAGE:   " << imagecount << "   GLOBAL-COUNT:  " << count << lineend;


		}

		off << "ENERGIE:    " << this->energies[num_images - 1] << "   FINAL" << lineend;
		off2 << "ENERGIE:    " << this->energies[num_images - 1] << "   FINAL" << lineend;
		print(name.str(), imagi, count);


	}
	else{

		name << "IMAGES_FINAL" << this->cPtr->mult_struc_counter << ".xyz";

		off << "ENERGIE:    " << this->energies[num_images - 1] << "   FINAL" << lineend;
		off2 << "ENERGIE:    " << this->energies[num_images - 1] << "   FINAL" << lineend;
		for (ptrdiff_t imagecount = num_images - 2; imagecount >= 1; imagecount--)
		{



			cPtr->set_xyz(imagi[imagecount]);

			



			off << "ENERGIE:    " << energies_NEB[imagecount]  << "     IMAGE:   " << imagecount << "   GLOBAL-COUNT:  " << count << lineend;
			off2 << "ENERGIE:    " << this->energies[imagecount] << "     IMAGE:   " << imagecount << "   GLOBAL-COUNT:  " << count << lineend;

		}

		print_rev(name.str(), imagi, count);


		if (ts == true)
		{
			off << "MIN         " << min_energies[count] << lineend;
			off << "TS          " << ts_energies[count] << lineend;

			printmono(name.str(), ts_pathstruc, count);
		}
	}



}

void neb::print(std::string const &name, std::vector <coords::Representation_3D> &print,ptrdiff_t &count)
{

  std::ofstream out (name.c_str(),std::ios::app);
  std::string temp;

  for (size_t i(0U); i< num_images;i++){
    out << "     " << N <<" IMAGE: " << i+1 << "  global counter:  "<<count <<lineend;
    for(size_t j = 0; j <N ;j++){

      out << std::right << std::setw(6) << j+1;
      out << std::left << "  " << std::setw(12)<<cPtr->atoms(j).symbol().substr(0U, 2U);	
      out << std::fixed << std::right << std::setprecision(6) << std::setw(13) << print[i][j].x();
      out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[i][j].y();
      out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[i][j].z();
      out << std::right << std::setw(6)<<cPtr->atoms(j).energy_type();
      size_t const bSize(cPtr->atoms(j).bonds().size());
      for (size_t n(0U); n<bSize; ++n)
      {
        out << std::right << std::setw(6) << cPtr->atoms(j).bonds()[n]+1U;
      }
      out << lineend;
    }
  }
  out.close();


}

void neb::print_rev(std::string const &name, std::vector <coords::Representation_3D> &print,ptrdiff_t &count)
{

  std::ofstream out (name.c_str(),std::ios::app);
  std::string temp;

  for (ptrdiff_t i =num_images-1; i >= 0;i--){
    out << "     " << N <<" IMAGE: " << i+1 << "  global counter:  "<<count <<lineend;
    for(size_t j = 0; j <N ;j++){

      out << std::right << std::setw(6) << j+1;
      out << std::left << "  " << std::setw(12)<<cPtr->atoms(j).symbol().substr(0U, 2U);	
      out << std::fixed << std::right << std::setprecision(6) << std::setw(13) << print[i][j].x();
      out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[i][j].y();
      out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[i][j].z();
      out << std::right << std::setw(6)<<cPtr->atoms(j).energy_type();
      size_t const bSize(cPtr->atoms(j).bonds().size());
      for (size_t n(0U); n<bSize; ++n)
      {
        out << std::right << std::setw(6) << cPtr->atoms(j).bonds()[n]+1U;
      }
      out << lineend;
    }
  }
  out.close();


}


void neb::printmono(std::string const &name, coords::Representation_3D &print,ptrdiff_t &count)
{

  std::ofstream out (name.c_str(),std::ios::app);
  std::string temp;


  out << "     " << N <<" IMAGE: " << "  global counter:  "<<count <<lineend;
  for(size_t j = 0; j <N ;j++){

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
    out << lineend;

  }

  out.close();


}




double neb::lbfgs (ptrdiff_t imagex)
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
  //optimizer.config.max_iterations = Config::get().optimization.local.bfgs.maxstep;
  optimizer.config.max_iterations = Config::get().neb.LBFGS_IT;
  optimizer.config.epsilon = (float)Config::get().optimization.local.bfgs.grad;
  optimizer(x);
  

 


  if (Config::get().general.verbosity > 4)
  {
    std::cout << "Optimization done (status " << optimizer.state() << "). Evaluations:" << optimizer.iter() << lineend;
  }

  //cPtr->set_xyz(scon::explicit_transform<coords::Representation_3D>(x.x));

  return cPtr->g();
}



double neb::g_new(ptrdiff_t im)
{

  //coords::Gradients_3D g;
	std::vector <coords::Representation_3D> temp(N);
	double Rp1mag, Rm1mag, energytemp;
	energytemp = this->cPtr->g();

	/* std::cout << "g vor neb: " << cPtr->g_xyz();
	std::cout<<energytemp<<lineend;*/
	Rp1.resize(num_images + cPtr->size());
	Rm1.resize(num_images + cPtr->size());
	Fvertical.resize(num_images + cPtr->size() );
	Fpar.resize(num_images + cPtr->size());
	//g.resize(cPtr->size());

	for (size_t i = 0; i<Fpar.size(); i++){
		Fpar[i].resize(cPtr->size());
	}

	/*for (size_t im = 1; im < num - 1; im++)
	{*/
	if (ClimbingImage == true && num_images == CIMaximum){
		double magni = 0.0;

		magni = dot(cPtr->g_xyz(), tau[im]);

		if (len(tau[im]) == 0.0){ magni = 0.0; }
		else{ magni = magni / len(tau[im]); }
		for (size_t i = 0; i < cPtr->size(); i++)
    {
			cPtr->update_g_xyz(i, cPtr->g_xyz(i) - tau[im][i] * magni*2.0);
		}



	}
	else{

		tauderiv = dot(cPtr->g_xyz(), tau[im]);
		if (len(tau[im]) == 0.0){ tauderiv = 0.0; }
		else{ tauderiv /= len(tau[im]); }

		for (size_t j = 0; j < cPtr->size(); j++){
		
			Rm1[j].x() = imagi[im - 1][j].x() - imagi[im][j].x();
			Rm1[j].y() = imagi[im - 1][j].y() - imagi[im][j].y();
			Rm1[j].z() = imagi[im - 1][j].z() - imagi[im][j].z();

			Rp1[j].x() = imagi[im + 1][j].x() - imagi[im][j].x();
			Rp1[j].y() = imagi[im + 1][j].y() - imagi[im][j].y();
			Rp1[j].z() = imagi[im + 1][j].z() - imagi[im][j].z();
		}
		Rm1mag = len(imagi[im - 1]) - len(imagi[im]);
		Rp1mag = len(imagi[im + 1]) - len(imagi[im]);
	
	
		for (size_t i = 0; i < cPtr->size(); i++){
			Fvertical[i].x() = cPtr->g_xyz(i).x() - tauderiv * tau[im][i].x();
			Fvertical[i].y() = cPtr->g_xyz(i).y() - tauderiv * tau[im][i].y();
			Fvertical[i].z() = cPtr->g_xyz(i).z() - tauderiv * tau[im][i].z();

			Fpar[im][i].x() = springconstant * (Rp1mag - Rm1mag) * tau[im][i].x();
			Fpar[im][i].y() = springconstant * (Rp1mag - Rm1mag) * tau[im][i].y();
			Fpar[im][i].z() = springconstant * (Rp1mag - Rm1mag) * tau[im][i].z();

		}

	
	}

	
	double fparvx(0), fparvy(0), fparvz(0);
	for (size_t j = 0; j < cPtr->size(); j++){

    auto const g = Fvertical[j] + Fpar[im][j];
	
		cPtr->update_g_xyz(j, g);
		fparvx += g.x();
		fparvy += g.y();
		fparvz += g.z();
	}
	//}
	grad_v = (fparvx + fparvy + fparvz) / cPtr->size();

	return energytemp;
}

void neb::calc_shift(void)
{
	double diff, gridp;
	VecDoub posx(num_images), posy(num_images), posz(num_images), gridx(num_images), gridy(num_images), gridz(num_images), shiftx(cPtr->size()), shifty(cPtr->size()), shiftz(cPtr->size());
	Doub x, y, z, pathlenx(0.0), pathleny(0.0), pathlenz(0.0);
	double distx(0.0), disty(0.0), distz(0.0);
	ptrdiff_t laf(0), counter(0);
	std::vector <double> temp;
	tempimage_ini = imagi[0];
	tempimage_final = imagi[num_images - 1];
	//std::string name;
	std::ostringstream name, name2, name3;
	std::vector <std::vector < std::vector < scon::c3<float> > > > position;
	scon::c3 <float> pos3;
	std::vector <coords::Cartesian_Point> interpol_position;
	coords::Representation_3D temp_int;
	position.resize(this->cPtr->size());




	for (size_t i = 0; i < this->cPtr->size(); i++){

		for (size_t j = 0; j < (num_images); j++){

			diff = (double)j / num_images;

			image_ini[j][i];
			posx[j] = imagi[j][i].x();
			posy[j] = imagi[j][i].y();
			posz[j] = imagi[j][i].z();
			gridx[j] = Doub(j);
			gridy[j] = Doub(j);
			gridz[j] = Doub(j);
			//std::cout << gridx[j] << "   " << posx[j] << endl;
		}
		//std::cout << lineend;

		Spline_interp splinex(gridx, posx);
		Spline_interp spliney(gridy, posy);
		Spline_interp splinez(gridz, posz);
		position[i].resize(num_images);

		for (size_t k = 0; k < num_images; k++)

		{
			
			laf = 0;
			gridp = k;
			size_t jj, it;
			Doub x, tnm, sumx, sumy, sumz, del, sx, sy, sz, skx0, sky0, skz0, skx, sky, skz;
			skx0 = abs((splinex.interp(k + 0.00001) - splinex.interp(k)) / 0.00001);
			skx = abs((splinex.interp((k + 1) + 0.00001) - splinex.interp(k + 1)) / 0.00001);
			sky0 = abs((spliney.interp(k + 0.00001) - spliney.interp(k)) / 0.00001);
			sky = abs((spliney.interp((k + 1) + 0.00001) - spliney.interp(k + 1)) / 0.00001);
			skz0 = abs((splinez.interp(k + 0.00001) - splinez.interp(k)) / 0.00001);
			skz = abs((splinez.interp((k + 1) + 0.00001) - splinez.interp(k + 1)) / 0.00001);
			//s = 0.5*((k + 1) - k)*(splinex.interp(k + 1) + splinex.interp(k));
			sx = 0.5*((k + 1) - k)*(skx + skx0);
			sy = 0.5*((k + 1) - k)*(sky + sky0);
			sz = 0.5*((k + 1) - k)*(skz + skz0);

			pathlenx = 0.0;
			pathleny = 0.0;
			pathlenz = 0.0;
			for (size_t m = 1; m < 8; m++)
			{
				for (it = 1, jj = 1; jj < m - 1; jj++) it <<= 1;

				tnm = it;
				if (tnm == 0)tnm = 1;
				del = ((k + 1) - k) / tnm;
				x = k + 0.5*del;
				/* std::cout << "tnm  " << x<< lineend;*/
				for (sumx = 0.0, sumy = 0.0, sumz = 0.0, jj = 0; jj < it; jj++, x += del)
				{
					skx = abs((splinex.interp(x + 0.00001) - splinex.interp(x)) / 0.00001);
					sky = abs((spliney.interp(x + 0.00001) - spliney.interp(x)) / 0.00001);
					skz = abs((splinez.interp(x + 0.00001) - splinez.interp(x)) / 0.00001);

					sumx += skx;
					sumy += sky;
					sumz += skz;

				}
				sx = 0.5*(sx + ((k + 1) - k)*sumx / tnm);
				sy = 0.5*(sy + ((k + 1) - k)*sumy / tnm);
				sz = 0.5*(sz + ((k + 1) - k)*sumz / tnm);



			}
			//std::cout << "i  " << i << "k  " << k << "  Wertx:     " << abs(sx) << "  Werty:     " << abs(sy) << "  Wertz:     " << abs(sz) << lineend;
			pathlenx += sx;
			pathleny += sy;
			pathlenz += sz;
			/*  std::cout << "total_length  " << pathlen/(num_images-2 )<< lineend;
			std::cout << "new shifted coordinate   "<<imagi[k][i].x() + (pathlen / (num_images - 2)) << lineend;*/
			std::cout << "K " << k << lineend;
			
			while (gridp <= (k+1))
			{
			
			 laf++;
			 gridp += 0.5;
			 //std::cout << gridp << endl;
			 x = splinex.interp(gridp);
			 y = spliney.interp(gridp);
			 z = splinez.interp(gridp);
			 distx += x;
			 disty += y;
			 distz += z;
			 pos3.x() = x;
			 pos3.y() = y;
			 pos3.z() = z;
		     
			 position[i][k].push_back(pos3);
			}
			std::cout << lineend;
			//std::cout << laf << lineend;

			distx /= laf;
			disty /= laf;
			distz /= laf;
			//std::cout << distx << "  " << disty << "   " << distz << lineend;
		}
		shiftx[i] = pathlenx / (num_images - 2);
		shifty[i] = pathleny / (num_images - 2);
		shiftz[i] = pathlenz / (num_images - 2);

		//cout << " shift x :  " << shiftx[i] << " shift y :  " << shifty[i] << " shift x :  " << shiftz[i] << lineend;
	}

	imagi_test.resize(num_images);
	//name2 << "BETA.xyz";
	//for (size_t j = 1; j < (num_images - 1); j++){



	//	for (size_t i = 0; i < this->cPtr->size(); i++){

	//		imagi[j][i].x() += shiftx[i];
	//		imagi[j][i].y() += shifty[i];
	//		imagi[j][i].z() += shiftz[i];


	//		imagi_test[j].push_back(imagi[j][i]);
	//		//image_ini[j].push_back(images[i]);


	//	}



	//	

	//	//out << "     " << N << " IMAGE_FINAL: " << "  global counter:  " << num_images << lineend;
	//	//for (size_t j = 0; j <N; j++) {

	//	//	out << std::right << std::setw(6) << j + 1;
	//	//	out << std::left << "  " << std::setw(12) << cPtr->atoms(j).symbol().substr(0U, 2U);
	//	//	out << std::fixed << std::right << std::setprecision(6) << std::setw(13) << imagi[num_images][j].x();
	//	//	out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << imagi[num_images][j].y();
	//	//	out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << imagi[num_images][j].z();
	//	//	out << std::right << std::setw(6) << cPtr->atoms(j).energy_type();
	//	//	size_t const bSize(cPtr->atoms(j).bonds().size());
	//	//	for (size_t n(0U); n<bSize; ++n)
	//	//	{
	//	//		out << std::right << std::setw(6) << cPtr->atoms(j).bonds()[n] + 1U;
	//	//	}
	//	//	out << lineend;

	//	//}
	//



	//	//printmono(name2.str(), imagi_test[j], counter);

	//	//cout << imagi_test[j] << endl;


	//}

	std::ofstream out("INTERPOL_preopt.arc", std::ios::app), out2("INTERPOL_opt.arc", std::ios::app),out3("ENERGY_INT_PREOPT.dat", std::ios::app), out4("ENERGY_INT_OPT.dat", std::ios::app);
	interpol_position.resize(N);
	tau_int.resize(N);
	interpol_position.resize(N);
	temp_int.resize(N);
	
	for (ptrdiff_t i = 1; i < num_images - 1; i++)
	{

		for (ptrdiff_t k = 0; k < laf; k++) {
			out << "     " << N << " IMAGE: " << k << "  global counter:  " << i << lineend;
			temp_int.clear();
			for (size_t j = 0; j < N; j++) {




				double abso(0.0);
				/*ofstream off ("energies");*/



				tau_int[j].x() = position[j][i][k].x() - position[j][i - 1][k].x() / abs(position[j][i][k].x() - position[j][i - 1][k].x()) + (position[j][i + 1][k].x() - position[j][i][k].x()) / abs(position[j][i + 1][k].x() - position[j][i][k].x());
				tau_int[j].y() = position[j][i][k].y() - position[j][i - 1][k].y() / abs(position[j][i][k].y() - position[j][i - 1][k].y()) + (position[j][i + 1][k].y() - position[j][i][k].y()) / abs(position[j][i + 1][k].y() - position[j][i][k].y());
				tau_int[j].z() = position[j][i][k].z() - position[j][i - 1][k].z() / abs(position[j][i][k].z() - position[j][i - 1][k].z()) + (position[j][i + 1][k].z() - position[j][i][k].z()) / abs(position[j][i + 1][k].z() - position[j][i][k].z());
				abso = 0.0;
				if (abs(scon::len(tau_int[j])) != (abs(scon::len(tau_int[j]))))
				{

					tau_int[j].x() = 1.0;
					tau_int[j].y() = 1.0;
					tau_int[j].z() = 1.0;

				}
				abso += tau_int[j].x()*tau_int[j].x();
				abso += tau_int[j].y()*tau_int[j].y();
				abso += tau_int[j].z()*tau_int[j].z();
				abso = sqrt(abso);
			
				
				if (abso == 0.0) abso = 1.0;
				tau_int[j].x() = tau_int[j].x() / abso;
				tau_int[j].y() = tau_int[j].y() / abso;
				tau_int[j].z() = tau_int[j].z() / abso;
				interpol_position[j].x() = position[j][i][k].x();
				interpol_position[j].y() = position[j][i][k].y();
				interpol_position[j].z() = position[j][i][k].z();
				temp_int.push_back(interpol_position[j]);






				out << std::right << std::setw(6) << j + 1;
				out << std::left << "  " << std::setw(12) << cPtr->atoms(j).symbol().substr(0U, 2U);
				out << std::fixed << std::right << std::setprecision(6) << std::setw(13) << position[j][i][k].x();
				out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << position[j][i][k].y();
				out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << position[j][i][k].z();
				out << std::right << std::setw(6) << cPtr->atoms(j).energy_type();
				size_t const bSize(cPtr->atoms(j).bonds().size());
				for (size_t n(0U); n < bSize; ++n)
				{
					out << std::right << std::setw(6) << cPtr->atoms(j).bonds()[n] + 1U;
				}
				out << lineend;


			}
			cPtr->set_xyz(temp_int);
			out3 << cPtr->g() << lineend;
			out4 << lbfgs_int(tau_int) << lineend;

			out2 << "     " << N << lineend;

			for (size_t j = 0; j < N; j++) {


				out2 << std::right << std::setw(6) << j + 1;
				out2 << std::left << "  " << std::setw(12) << cPtr->atoms(j).symbol().substr(0U, 2U);
				out2 << std::fixed << std::right << std::setprecision(6) << std::setw(13) << cPtr->xyz(j).x();
				out2 << std::fixed << std::right << std::setprecision(6) << std::setw(12) << cPtr->xyz(j).y();
				out2 << std::fixed << std::right << std::setprecision(6) << std::setw(12) << cPtr->xyz(j).z();
				out2 << std::right << std::setw(6) << cPtr->atoms(j).energy_type();
				size_t const bSize(cPtr->atoms(j).bonds().size());
				for (size_t n(0U); n < bSize; ++n)
				{
					out2 << std::right << std::setw(6) << cPtr->atoms(j).bonds()[n] + 1U;
				}
				out2 << lineend;


			}


		}
	}
	out2.close();
	out.close();
}



double neb::lbfgs_int(std::vector <scon::c3 <float> > t)
{


	using namespace  optimization::local;
	//typedef coords::Container<scon::c3<float>> nc3_type;
	// Create optimizer
	
	coords::Cartesian_Point g;
	auto optimizer = make_lbfgs(
		make_more_thuente(GradCallBack_int(*this))
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
		std::cout << "Optimization done (status " << optimizer.state() << "). Evaluations:" << optimizer.iter() << lineend;
	}

	//cPtr->set_xyz(scon::explicit_transform<coords::Representation_3D>(x.x));

	return cPtr->g();
}

double neb::g_int(std::vector<scon::c3 <float> > t)
{

	double cosi, energytemp, pot(0.0), dpot(0.0);
	coords::Cartesian_Point rv_p;
	energytemp = cPtr->g();
	if (Config::get().neb.OPTMODE == "PROJECTED")
	{
		for (size_t i = 0; i < cPtr->size(); i++)
		{

			auto L = scon::geometric_length(t[i]);
			if (L != 0.0)
			{
				cosi = (cPtr->g_xyz(i).x()*t[i].x() + cPtr->g_xyz(i).y()*t[i].y() + cPtr->g_xyz(i).z()*t[i].z()) / (L * L);
			}
			else
			{
				cosi = 0.0;
			}
			rv_p.x() = cPtr->g_xyz(i).x() - (cosi * t[i].x());
			rv_p.y() = cPtr->g_xyz(i).y() - (cosi * t[i].y());
			rv_p.z() = cPtr->g_xyz(i).z() - (cosi * t[i].z());

			cPtr->update_g_xyz(i, rv_p);









		}
	}
	else if (Config::get().neb.OPTMODE == "BIAS")
	{
		for (size_t i = 0; i < cPtr->size(); i++)
		{

			if (scon::geometric_length(t[i]) != 0.0)
			{
				cosi = (cPtr->xyz(i).x()*t[i].x() + cPtr->xyz(i).y()*t[i].y() + cPtr->xyz(i).z()*t[i].z() - t[i].x()*t[i].x() - t[i].y()*t[i].y() - t[i].z()*t[i].z()) / sqrt(t[i].x()*t[i].x() + t[i].y()*t[i].y() + t[i].z()*t[i].z());
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



//}

/***********************************************************************
************************************************************************
****                                                                ****
****    *******               **                                    ****
****    *******               **                                    ****
****         **               **        **                          ****
****         **               **                                    ****
****         **   **    **    **        **    **** **   ** ****     ****
****         **   **    **    **        **   **  ****   ****  **    ****
****    ******     **  **     **  **    **   **  ****   **    **    ****
****    *****       ****       ****     **    **** **   **    **    ****
****                                                                ****
************************************************************************
***********************************************************************/

void neb::opt_internals(ptrdiff_t &count, const std::vector<std::vector<int> >& atoms_remember){

	coords::Ensemble_PES *ep = new coords::Ensemble_PES;
	std::ostringstream name;
	
	

	//Just to get xyz-data beginning by the origin
	//
	for (size_t i = 0U; i<num_images; i++){
		cPtr->set_xyz(imagi[i]);
		cPtr->to_internal();
		cPtr->to_xyz();
		imagi[i] = cPtr->xyz();
	}
	//Get output for comparison
	name << "Formatted.xyz";
	print(name.str(), imagi, count);
	name.str("");
	//std::cout <<"SIZE_ALL "<< atoms_remember[1].size() << lineend;
	//Loop to optimize every image
	for (size_t i = 1U; i<num_images - 1; i++){

		//reverse every fixation
		defix_all();

		cPtr->set_xyz(imagi[i]);

		//fill ep with current values
		ep->clear();
		ep->push_back(imagi[i]);

		//fix atoms for this and last image to guarantee no change in dihedrals during the change
		//std::cout << "SIZE_B " << atoms_remember[i].size() << lineend;
		execute_fix(atoms_remember[i]);
		execute_fix(atoms_remember[i - 1]);

		//Just to show which atoms are fixed
		for (size_t j = 0U; j<N; j++){
			std::cout << j + 1 << ": " << (cPtr->atoms(j).fixed() ? "fixed" : "not fixed") << std::endl;
		}

		//creating a new class for Monte Carlo optimization.
		
    name << "_Image_" << i + 1;

		mc = new optimization::global::optimizers::monteCarlo(*cPtr, *ep, name.str(), true);
		mc->run(Config::get().optimization.global.iterations);
		//write range
		mc->write_range(name.str());
		name.str("");
		//get garbage lost
		delete mc;
		mc = nullptr;
	}

	//here too
	delete ep;
	ep = nullptr;

}

void neb::internal_execute(std::vector <coords::Representation_3D> &input, std::vector<std::vector<int> >& atoms_remember){

	size_t i = 0U, j = 0U, imgs = num_images, atom_iter = 0U;

	std::vector<std::string> Atom_name_input, distance_string_swap, angle_string_swap, dihedral_string_swap;
	std::vector<std::vector<std::string> > DAD_names;
	std::vector<std::vector<std::vector<std::string> > > name;

	std::vector<double> x_input, y_input, z_input;
	std::vector<std::vector<double> > imgs_x_input, imgs_y_input, imgs_z_input, dist, anglevec, dihedralvec;
	std::vector<std::vector<std::vector<double> > > dist_all, anglevec_all, dihedralvec_all, compare;

	std::vector<int> int_swap;
	std::vector<std::vector<int> > bonds, which_bonds;
	std::vector<std::vector<std::vector<int> > > which_img;
	std::vector<std::vector<std::vector<std::vector<int> > > > involved_bonds;

	std::vector<std::vector<size_t> > i_remember, j_remember, i_remember_all, j_remember_all;

	std::ofstream myfile;

	std::ostringstream string_swap;
    
	//clearing variables used to carry results outside the function
	atoms_remember.clear();

	//writing input in own variables
	for (i = 0U; i<N; i++){
		int_swap.clear();

		Atom_name_input.push_back(cPtr->atoms(i).symbol().substr(0U, 2U));
		x_input.push_back(input[0][i].x());
		y_input.push_back(input[0][i].y());
		z_input.push_back(input[0][i].z());

		int_swap.push_back(i + 1);
		for (j = 0U; j<cPtr->atoms(i).bonds().size(); j++)
		{
			int_swap.push_back(cPtr->atoms(i).bonds()[j] + 1U);
		}
		bonds.push_back(int_swap);

	}

	imgs_x_input.push_back(x_input);
	imgs_y_input.push_back(y_input);
	imgs_z_input.push_back(z_input);

	//because variable bonds won't change during images just the x,y and z values are needed from now on
	for (i = 1U; i<imgs; i++){
		x_input.clear();
		y_input.clear();
		z_input.clear();

		for (j = 0U; j<N; j++){
			x_input.push_back(input[i][j].x());
			y_input.push_back(input[i][j].y());
			z_input.push_back(input[i][j].z());
		}
		imgs_x_input.push_back(x_input);
		imgs_y_input.push_back(y_input);
		imgs_z_input.push_back(z_input);
	}

	//A sorting was implimented in the beginning but the need was bypassed and thrown out for more performance
	//
	//Sort(imgs,Atom_name_input, imgs_x_input, imgs_y_input,imgs_z_input,bonds);

	//array which is worked with is created
	//getting the backbone of the molecule on which to move hand over hand along the backbone
	for (i = 0U; i<N; i++){
		if (bonds[i].size() > 2){
			which_bonds.push_back(bonds[i]); atom_iter++;
		}
	}

	//std::cout << atom_iter << std::endl;

	//now the first values can be obtained and it has been taken care that there is no need to
	//use such a expensive function again
	get_values(imgs_x_input[0], imgs_y_input[0], imgs_z_input[0], atom_iter, which_bonds, &dist, &anglevec, &dihedralvec, &involved_bonds, bonds);

	//The name of the atoms are obtained. Just for output
	//
	for (i = 0U; i < involved_bonds[DIST].size(); i++){
		for (j = 0U; j<involved_bonds[DIST][i].size(); j++){
			double_or_not(involved_bonds[DIST], involved_bonds[DIST][i][j], i, j, dist);
		}
		for (j = 0U; j<involved_bonds[ANGLE][i].size(); j++){
			double_or_not(involved_bonds[ANGLE], involved_bonds[ANGLE][i][j], i, j, anglevec);
		}
		for (j = 0U; j<involved_bonds[DIHEDRAL][i].size(); j++){
			double_or_not(involved_bonds[DIHEDRAL], involved_bonds[DIHEDRAL][i][j], i, j, dihedralvec);
		}
	}

	for (i = 0U; i < involved_bonds[DIST].size(); i++){

		distance_string_swap.clear();
		angle_string_swap.clear();
		dihedral_string_swap.clear();

		DAD_names.clear();
		for (j = 0U; j<involved_bonds[DIST][i].size(); j++){
			string_swap << involved_bonds[DIST][i][j][0] << Atom_name_input[involved_bonds[DIST][i][j][0] - 1] << " ";
			string_swap << involved_bonds[DIST][i][j][1] << Atom_name_input[involved_bonds[DIST][i][j][1] - 1];
			distance_string_swap.push_back(string_swap.str()); string_swap.str("");
		}
		for (j = 0U; j<involved_bonds[ANGLE][i].size(); j++){
			string_swap << involved_bonds[ANGLE][i][j][0] << Atom_name_input[involved_bonds[ANGLE][i][j][0] - 1] << " ";
			string_swap << involved_bonds[ANGLE][i][j][1] << Atom_name_input[involved_bonds[ANGLE][i][j][1] - 1] << " ";
			string_swap << involved_bonds[ANGLE][i][j][2] << Atom_name_input[involved_bonds[ANGLE][i][j][2] - 1];
			angle_string_swap.push_back(string_swap.str()); string_swap.str("");
		}
		for (j = 0U; j<involved_bonds[DIHEDRAL][i].size(); j++){
			string_swap << involved_bonds[DIHEDRAL][i][j][0] << Atom_name_input[involved_bonds[DIHEDRAL][i][j][0] - 1] << " ";
			string_swap << involved_bonds[DIHEDRAL][i][j][1] << Atom_name_input[involved_bonds[DIHEDRAL][i][j][1] - 1] << " ";
			string_swap << involved_bonds[DIHEDRAL][i][j][2] << Atom_name_input[involved_bonds[DIHEDRAL][i][j][2] - 1] << " ";
			string_swap << involved_bonds[DIHEDRAL][i][j][3] << Atom_name_input[involved_bonds[DIHEDRAL][i][j][3] - 1];
			dihedral_string_swap.push_back(string_swap.str()); string_swap.str("");
		}
		DAD_names.push_back(distance_string_swap);
		DAD_names.push_back(angle_string_swap);
		DAD_names.push_back(dihedral_string_swap);

		name.push_back(DAD_names);
	}
	//end

	dist_all.push_back(dist);
	anglevec_all.push_back(anglevec);
	dihedralvec_all.push_back(dihedralvec);

	dist.clear();
	anglevec.clear();
	dihedralvec.clear();

	//values for the other images is obtained
	for (i = 1U; i<imgs; i++){

		get_values(imgs_x_input[i], imgs_y_input[i], imgs_z_input[i], &dist, &anglevec, &dihedralvec, involved_bonds, bonds);

		dist_all.push_back(dist);
		anglevec_all.push_back(anglevec);
		dihedralvec_all.push_back(dihedralvec);

		dist.clear();
		anglevec.clear();
		dihedralvec.clear();

	}

	//to find out where the biggest change in dihedrals during change of images is
	biggest(dist_all, anglevec_all, dihedralvec_all,
		compare, which_img, i_remember, j_remember, i_remember_all, j_remember_all);

	atoms_remember.resize(num_images);
	for (size_t li(0U); li < num_images - 1; li++)
	{
		if (Config::get().neb.NUMBER_OF_DIHEDRALS > i_remember_all[li].size()) Config::set().neb.NUMBER_OF_DIHEDRALS = i_remember_all[li].size();
		for (size_t lk(i_remember_all[li].size()-Config::get().neb.NUMBER_OF_DIHEDRALS ); lk < i_remember_all[li].size(); lk++)
		{
			
			//std::cout << i_remember_all[li][lk] << "   " << j_remember_all[li][lk] << lineend;
		
			atoms_remember[li].push_back(involved_bonds[DIHEDRAL][i_remember_all[li][lk]][j_remember_all[li][lk]][0]);
			atoms_remember[li].push_back(involved_bonds[DIHEDRAL][i_remember_all[li][lk]][j_remember_all[li][lk]][1]);
			atoms_remember[li].push_back(involved_bonds[DIHEDRAL][i_remember_all[li][lk]][j_remember_all[li][lk]][2]);
			atoms_remember[li].push_back(involved_bonds[DIHEDRAL][i_remember_all[li][lk]][j_remember_all[li][lk]][3]);
			atoms_remember.push_back(involved_bonds[DIST][i_remember[DIST][li]][j_remember[DIST][li]]);
			//atoms_remember[li].push_back(involved_bonds[DIHEDRAL][i_remember_all[li][lk]][j_remember_all[li][lk]]);
		}
	}

	for (size_t lo(0U); lo <num_images-1; lo++)
	{
		for (size_t lp(0U); lp < atoms_remember[lo].size(); lp++)
		{
			//std::cout << "ATOMS_REM " <<"lo "<<lo<<"  lp "<<lp <<"  "<< atoms_remember[lo][lp] << lineend;
		}
	}
	//for (i = 0U; i<i_remember[DIST].size(); i++){
	//	std::cout << "SIZE#" << i_remember[DIST].size()<<lineend;
	//	std::cout << i_remember[DIST][i] << " " << j_remember[DIST][i] << std::endl;
	//	atoms_remember.push_back(involved_bonds[DIST][i_remember[DIST][i]][j_remember[DIST][i]]);
	//	//atoms_remember.push_back(involved_bonds[DIHEDRAL][i_remember[DIHEDRAL][i]][j_remember[DIHEDRAL][i]]);
	//}

	//write output
	//
	myfile.open("Internals.txt");


	for (size_t k = 0U; k<num_images; k++){


		for (size_t i = 0U; i<dist_all[k].size(); i++){
			if (dist_all[k][i].size()>anglevec_all[k][i].size() && dist_all[k][i].size()>dihedralvec_all[k][i].size())for (size_t j = 0U; j<dist_all[k][i].size(); j++){
				if (dist_all[k][i].size()>j){
					myfile << std::right << std::setw(10) << name[i][DIST][j] << " ";
					myfile << std::left << std::setw(10) << dist_all[k][i][j];
				}
				else myfile << std::right << std::setw(21) << " ";
				if (anglevec_all[k][i].size()>j){
					myfile << std::right << std::setw(10) << name[i][ANGLE][j] << " ";
					myfile << std::left << std::setw(10) << anglevec_all[k][i][j];
				}
				else myfile << std::right << std::setw(21) << " ";
				if (dihedralvec_all[k][i].size()>j){
					myfile << std::right << std::setw(10) << name[i][DIHEDRAL][j] << " ";
					myfile << std::left << std::setw(10) << dihedralvec_all[k][i][j];
				}
				else myfile << std::right << std::setw(21) << " ";
				myfile << std::endl;
			}
			else if (anglevec_all[k][i].size()>dist_all[k][i].size() && anglevec_all[k][i].size()>dihedralvec_all[k][i].size())for (size_t j = 0U; j<anglevec_all[k][i].size(); j++){
				if (dist_all[k][i].size()>j){

					myfile << std::right << std::setw(10) << name[i][DIST][j] << " ";
					myfile << std::left << std::setw(10) << dist_all[k][i][j];
				}
				else myfile << std::right << std::setw(21) << " ";
				if (anglevec_all[k][i].size()>j){
					myfile << std::right << std::setw(10) << name[i][ANGLE][j] << " ";
					myfile << std::left << std::setw(10) << anglevec_all[k][i][j];
				}
				else myfile << std::right << std::setw(21) << " ";
				if (dihedralvec_all[k][i].size()>j){
					myfile << std::right << std::setw(10) << name[i][DIHEDRAL][j] << " ";
					myfile << std::left << std::setw(10) << dihedralvec_all[k][i][j];
				}
				else myfile << std::right << std::setw(21) << " ";
				myfile << std::endl;
			}
			else for (size_t j = 0U; j<dihedralvec[i].size(); j++){
				if (dist[i].size()>j){
					myfile << std::right << std::setw(10) << name[i][DIST][j] << " ";
					myfile << std::left << std::setw(10) << dist_all[k][i][j];
				}
				else myfile << std::right << std::setw(21) << " ";
				if (anglevec[i].size()>j){
					myfile << std::right << std::setw(10) << name[i][ANGLE][j] << " ";
					myfile << std::left << std::setw(10) << anglevec_all[k][i][j];
				}
				else myfile << std::right << std::setw(21) << " ";
				if (dihedralvec[i].size()>j){
					myfile << std::right << std::setw(10) << name[i][DIHEDRAL][j] << " ";
					myfile << std::left << std::setw(10) << dihedralvec_all[k][i][j];
				}
				else myfile << std::right << std::setw(21) << " ";
				myfile << std::endl;
			}
			myfile << std::endl;
		}
	}
	myfile.close();
	//output written

	//get fixations in an output
	myfile.open("Fixation.txt");

	for (i = 0U; i<atoms_remember.size(); i++){
		myfile << "Image " << i + 1 << ": ";
		for (j = 0U; j<atoms_remember[i].size(); j++){
			myfile << atoms_remember[i][j] << " ";
		}
		myfile << std::endl;
	}

	myfile.close();

	//get lost of garbage
	imgs_x_input.clear();
	imgs_y_input.clear();
	imgs_z_input.clear();
	x_input.clear();
	y_input.clear();
	z_input.clear();
	bonds.clear();
	int_swap.clear();
	Atom_name_input.clear();
	which_bonds.clear();
	dist.clear();
	anglevec.clear();
	dihedralvec.clear();
	dist_all.clear();
	anglevec_all.clear();
	dihedralvec_all.clear();
	compare.clear();
	which_img.clear();
	name.clear();
	DAD_names.clear();
	distance_string_swap.clear();
	angle_string_swap.clear();
	dihedral_string_swap.clear();
	involved_bonds.clear();
	i_remember.clear();
	j_remember.clear();

}

void neb::no_torsion(const size_t& image_remember, const std::vector<int>& atoms_remember){
	for (size_t i = 0U; i<atoms_remember.size(); i++){
		imagi[image_remember][atoms_remember[i] - 1].x() = imagi[image_remember - 1][atoms_remember[i] - 1].x();
		imagi[image_remember][atoms_remember[i] - 1].y() = imagi[image_remember - 1][atoms_remember[i] - 1].y();
		imagi[image_remember][atoms_remember[i] - 1].z() = imagi[image_remember - 1][atoms_remember[i] - 1].z();
	}
}

//void execute_fix
//Function to fixate atoms
//

void neb::execute_fix(const std::vector<int>& atoms_remember){
	//std::cout <<"SIZE_exe "<< atoms_remember.size() << endl;
	for (size_t i = 0U; i<atoms_remember.size(); i++){
		cPtr->set_fix(atoms_remember[i] - 1,true);
		//std::cout <<"REM "<< atoms_remember[i] - 1 << lineend;
		//scon::insert_unique_sorted(Config::set().coords.fixed, atoms_remember[i] - 1);
	}
}//end execute_fix

void neb::execute_defix(const std::vector<int>& atoms_remember){
	for (size_t i = 0U; i<atoms_remember.size(); i++){
		cPtr->set_fix(atoms_remember[i] - 1, false);
		
	}

}

//void defix_all
//Function to loosen all atoms
void neb::defix_all(){
	for (size_t i = 0U; i<N; i++){
		cPtr->set_fix(i, false);
	}
}//end defix_all

//void biggest
//Function to get the biggest change during two images
//
void neb::biggest(const std::vector<std::vector<std::vector<double> > >& dist, const std::vector<std::vector<std::vector<double> > >& anglevec,
	const std::vector<std::vector<std::vector<double> > >& dihedralvec, std::vector<std::vector<std::vector<double> > >& , std::vector<std::vector<std::vector<int> > >& , std::vector<std::vector<size_t> >& i_remember,
	std::vector<std::vector<size_t> >& j_remember, std::vector<std::vector<size_t> > & i_remember_all, std::vector<std::vector<size_t> > & j_remember_all){

	size_t i = 0U, j = 0U, k = 0U;
	double double_swap;

	std::vector<double> comp_value_swap;

	std::vector<int> bonds_value_swap;
	std::vector<std::vector<int> > bonds_atom_swap;

	std::ofstream myfile;
	std::ostringstream name;

	std::vector<size_t> i_temp, j_temp, i_remembered, j_remembered;
	std::vector< std::vector <size_t> > i_remember_two, j_remember_two;

	i_remember.clear();
	j_remember.clear();
	i_remember_all.clear();

	i_remember_two.clear();
	j_remember_two.clear();
	i_remember.resize(DIHEDRAL + 1);
	j_remember.resize(DIHEDRAL + 1);

	i_remembered.clear();
	j_remembered.clear();
	i_remembered.resize(num_images - 1);
	j_remembered.resize(num_images - 1);

	i_remember_two.resize(num_images - 1);
	j_remember_two.resize(num_images - 1);

	//start a loop for every image
	for (i = 0U; i<num_images - 1; i++){
		comp_value_swap.clear();
		i_temp.clear();
		j_temp.clear();

		//Fill the compare array with the first values to initialize the comparison
		for (j = 0U; j<dist[i].size(); j++){
			for (k = 0U; k<dist[i][j].size(); k++){
				comp_value_swap.push_back(fabs(dist[i][j][k] - dist[i + 1][j][k]));
				i_temp.push_back(j);
				j_temp.push_back(k);
			}
		}
		//Compare the change of every value to get the biggest
		double_swap = comp_value_swap[0];
		for (j = 1U; j < comp_value_swap.size(); j++){
			if (double_swap<comp_value_swap[j]){
				double_swap = comp_value_swap[j];
				i_remembered[i] = i_temp[j];
				j_remembered[i] = j_temp[j];
			}
		}
	}
	i_remember[DIST] = i_remembered;
	j_remember[DIST] = j_remembered;

	for (i = 0U; i<num_images - 1; i++){
		comp_value_swap.clear();
		i_temp.clear();
		j_temp.clear();

		for (j = 0U; j<anglevec[i].size(); j++){
			for (k = 0U; k<anglevec[i][j].size(); k++){
				comp_value_swap.push_back(fabs(anglevec[i][j][k] - anglevec[i + 1][j][k]));
				i_temp.push_back(j);
				j_temp.push_back(k);
			}
		}

		double_swap = comp_value_swap[0];
		for (j = 1U; j < comp_value_swap.size(); j++){
			if (double_swap<comp_value_swap[j]){
				double_swap = comp_value_swap[j];
				i_remembered[i] = i_temp[j];
				j_remembered[i] = j_temp[j];
			}
		}
	}
	i_remember[ANGLE] = i_remembered;
	j_remember[ANGLE] = j_remembered;

	myfile.open("Logfile.txt");

	for (i = 0U; i<num_images - 1; i++){
		comp_value_swap.clear();
		i_temp.clear();
		j_temp.clear();

		for (j = 0U; j<dihedralvec[i].size(); j++){
			for (k = 0U; k<dihedralvec[i][j].size(); k++){
				if ((dihedralvec[i][j][k]<0.0&&dihedralvec[i + 1][j][k]<0.0) || (dihedralvec[i][j][k]>0.0&&dihedralvec[i + 1][j][k]>0.0)){
					comp_value_swap.push_back(fabs(fabs(dihedralvec[i][j][k]) - fabs(dihedralvec[i + 1][j][k])));
					i_temp.push_back(j);
					j_temp.push_back(k);
				}
				else{
					if (fabs(dihedralvec[i][j][k]) + fabs(dihedralvec[i + 1][j][k])>180.0){
						comp_value_swap.push_back(360.0 - fabs(dihedralvec[i][j][k]) - fabs(dihedralvec[i + 1][j][k]));
						i_temp.push_back(j);
						j_temp.push_back(k);
					}
					else{
						comp_value_swap.push_back(fabs(fabs(dihedralvec[i][j][k]) + fabs(dihedralvec[i + 1][j][k])));
						i_temp.push_back(j);
						j_temp.push_back(k);
					}
				}
			}
		}

		for (j = 0U; j<comp_value_swap.size(); j++){
			if ((int)comp_value_swap[j]>0) myfile << "Image " << i + 1 << " to " << i + 2 << ": " << (int)comp_value_swap[j] << std::endl;
		}

		double_swap = comp_value_swap[0];
		bool swapped(true);
		
		while(swapped == true){
			
			swapped = false;
			for (j = 0U; j < comp_value_swap.size(); j++){

				if (comp_value_swap[j] > comp_value_swap[j + 1]){
					swapped = true;
					int temp_i, temp_j;
					double temp_d;

					temp_i = i_temp[j];
					i_temp[j] = i_temp[j + 1];
					i_temp[j + 1] = temp_i;

					temp_j = j_temp[j];
					j_temp[j] = j_temp[j + 1];
					j_temp[j + 1] = temp_j;

					temp_d = comp_value_swap[j];
					comp_value_swap[j] = comp_value_swap[j + 1];
					comp_value_swap[j + 1] = temp_d;
					//newnn = j + 1;
					/*std::cout << " IMG " << i << "DOUBLE SWAP " << double_swap << " Index i " << i_temp[j] << " Index j " << j_temp[j] << lineend;
					i_remember_two[i].push_back(i_temp[j]);
					j_remember_two[i].push_back(j_temp[j]);
					i_remembered[i] = i_temp[j];
					j_remembered[i] = j_temp[j];*/
				}
				//nn = newnn;
			}
		
			
		} 

			for (j = 0U; j < comp_value_swap.size(); j++){

				//std::cout << "COMP_swap " << comp_value_swap[j] << " Index i " << i_temp[j] << " Index j " << j_temp[j] << lineend;
				i_remember_two[i].push_back(i_temp[j]);
				j_remember_two[i].push_back(j_temp[j]);

			}
			
		
	}
	myfile.close();
	i_remember[DIHEDRAL] = i_remembered;
	j_remember[DIHEDRAL] = j_remembered;
	i_remember_all.resize(i_remember_two.size());
	j_remember_all.resize(j_remember_two.size());
	i_remember_all = i_remember_two;
	j_remember_all = j_remember_two;
	i_remembered.clear();
	j_remembered.clear();
}

//void double_or_not
//Function to find out if a distance, angle or dihedral is double
//
void neb::double_or_not(std::vector<std::vector<std::vector<int> > >& involved_bonds, const std::vector<int>& pivot, const size_t& a, const size_t& b, std::vector<std::vector<double> >& vec){

	bool all_the_same;

	size_t i, j, k, l;
	std::vector<bool> its_double;

	its_double.resize(pivot.size());

	//check every whatever if it is double by compareing the atoms involved. Afterwards killing all double values
	for (i = 0U; i < involved_bonds.size(); i++){
		for (j = 0U; j < involved_bonds[i].size(); j++){
			if ((j<(involved_bonds[i].size() - 1)) && (i == a && j == b)) j++;
			else if ((j == (involved_bonds[i].size() - 1)) && (i == a && j == b)) break;
			for (k = 0U; k<pivot.size(); k++){
				its_double[k] = false;
			}
			for (k = 0U; k<involved_bonds[i][j].size(); k++){
				for (l = 0U; l<pivot.size(); l++){
					if (involved_bonds[i][j][k] == pivot[l]){
						its_double[k] = true;
						break;
					}
				}
			}
			all_the_same = true;
			if (pivot.size()>0) for (k = 0U; k<pivot.size() - 1; k++){
				if (its_double[k] != its_double[k + 1]) all_the_same = false;
			}
			if (all_the_same&&its_double[0]){
				vec[i].erase(vec[i].begin() + j);
				involved_bonds[i].erase(involved_bonds[i].begin() + j);
			}
		}
	}
}

bool neb::dihedral_or_not(const int& first_atom, const std::vector<int>& dist, const std::vector<int>& ang, const int& last_atom){

	if ((first_atom == dist[1] || first_atom == dist[2] || first_atom == dist[3] || first_atom == dist[4]) && (last_atom == ang[1] || last_atom == ang[2] || last_atom == ang[3] || last_atom == ang[4])) return true;
	else if ((first_atom == ang[1] || first_atom == ang[2] || first_atom == ang[3] || first_atom == ang[4]) && (last_atom == dist[1] || last_atom == dist[2] || last_atom == dist[3] || last_atom == dist[4])) return true;
	else return false;

}

void neb::get_values(std::vector<double>& x_val, std::vector<double>& y_val, std::vector<double>& z_val, const size_t& atom_iter, std::vector<std::vector<int> >& which_bonds, std::vector<std::vector<double> >* dist, std::vector<std::vector<double> >* anglevec, std::vector<std::vector<double> >* dihedralvec, std::vector<std::vector<std::vector<std::vector<int> > > >* involved_bonds, std::vector<std::vector<int> >& bonds){

	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	size_t i, j, k, l;
	int a, b, c, d;

	std::vector<int> int_swap;
	std::vector<std::vector<int> > dist_int_swap, angle_int_swap, dihedral_int_swap;
	std::vector<std::vector<std::vector<int> > > dist_int, angle_int, dihedral_int;

	std::vector<double> dist_swap, angle_swap, dihedral_swap;

	//search down the backbone
	for (i = 0U; i<atom_iter; i++){

		dist_swap.clear();
		angle_swap.clear();
		dihedral_swap.clear();

		dist_int_swap.clear();
		angle_int_swap.clear();
		dihedral_int_swap.clear();

		//start atom
		x1 = x_val[which_bonds[i][0] - 1];
		y1 = y_val[which_bonds[i][0] - 1];
		z1 = z_val[which_bonds[i][0] - 1];

		for (j = 0U; j<which_bonds[i].size() - 1; j++){

			//first neighbor to get
			x2 = x_val[which_bonds[i][j + 1] - 1];
			y2 = y_val[which_bonds[i][j + 1] - 1];
			z2 = z_val[which_bonds[i][j + 1] - 1];

			//getting distance
			dist_swap.push_back(distance(x1, x2, y1, y2, z1, z2));

			a = which_bonds[i][0];
			b = which_bonds[i][j + 1];


			//safe the bonds involved
			int_swap.push_back(a);
			int_swap.push_back(b);

			dist_int_swap.push_back(int_swap);
			int_swap.clear();

			for (k = 0U; k<(which_bonds[i].size() - (j + 2)); k++){

				//second neighbor
				x3 = x_val[which_bonds[i][(j + 1) + (k + 1)] - 1];
				y3 = y_val[which_bonds[i][(j + 1) + (k + 1)] - 1];
				z3 = z_val[which_bonds[i][(j + 1) + (k + 1)] - 1];

				//getting angle
				angle_swap.push_back(angle(x1, x2, x3, y1, y2, y3, z1, z2, z3));

				//safe the bonds involved
				a = which_bonds[i][0];
				b = which_bonds[i][j + 1];
				c = which_bonds[i][(j + 1) + (k + 1)];

				int_swap.push_back(a);
				int_swap.push_back(b);
				int_swap.push_back(c);

				angle_int_swap.push_back(int_swap);
				int_swap.clear();

				//search on the start atom for dihedrals
				for (l = 0U; l<which_bonds[i].size() - 1; l++){

					//of cause it shoudn't be an atom which takes part either way
					if (which_bonds[i][l + 1] != which_bonds[i][(j + 1) + (k + 1)] && which_bonds[i][l + 1] != which_bonds[i][j + 1] && which_bonds[i][l + 1] != which_bonds[i][0]){

						//get the values
						x4 = x_val[which_bonds[i][l + 1] - 1];
						y4 = y_val[which_bonds[i][l + 1] - 1];
						z4 = z_val[which_bonds[i][l + 1] - 1];

						//same as above
						a = which_bonds[i][0];
						b = which_bonds[i][j + 1];
						c = which_bonds[i][(j + 1) + (k + 1)];
						d = which_bonds[i][l + 1];

						//calculating diheral angle
						dihedral_swap.push_back(dihedral_same_atom(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4));

						int_swap.push_back(a);
						int_swap.push_back(b);
						int_swap.push_back(c);
						int_swap.push_back(d);

						dihedral_int_swap.push_back(int_swap);
						int_swap.clear();
					}
				}

				//Search on the "distance" atom for more dihedrals (same procedure as on the other dihedral search)
				//of course first check if it is not terminal
				if (bonds[which_bonds[i][j + 1] - 1].size()>1){
					for (l = 0U; l<(bonds[which_bonds[i][j + 1] - 1].size() - 1); l++){

						if (bonds[which_bonds[i][j + 1] - 1][l + 1] != which_bonds[i][0] && bonds[which_bonds[i][j + 1] - 1][l + 1] != bonds[i][(j + 1) + (k + 1)]){

							x4 = x_val[bonds[which_bonds[i][j + 1] - 1][l + 1] - 1];
							y4 = y_val[bonds[which_bonds[i][j + 1] - 1][l + 1] - 1];
							z4 = z_val[bonds[which_bonds[i][j + 1] - 1][l + 1] - 1];

							dihedral_swap.push_back(dihedral(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4));

							a = which_bonds[i][0];
							b = which_bonds[i][j + 1];
							c = which_bonds[i][(j + 1) + (k + 1)];
							d = bonds[which_bonds[i][j + 1] - 1][l + 1];

							int_swap.push_back(a);
							int_swap.push_back(b);
							int_swap.push_back(c);
							int_swap.push_back(d);

							dihedral_int_swap.push_back(int_swap);
							int_swap.clear();
						}
					}
				}


				//Search on the "angle" atom for more dihedrals (same procedure as on the other dihedral search)
				//of cause first check if it is not terminal
				if (bonds[which_bonds[i][(j + 1) + (k + 1)] - 1].size()>1){
					for (l = 0U; l<(bonds[which_bonds[i][(j + 1) + (k + 1)] - 1].size() - 1); l++){

						if (bonds[which_bonds[i][(j + 1) + (k + 1)] - 1][l + 1] != which_bonds[i][(j + 1) + (k + 1)] && bonds[which_bonds[i][(j + 1) + (k + 1)] - 1][l + 1] != which_bonds[i][0] && bonds[which_bonds[i][(j + 1) + (k + 1)] - 1][l + 1] != which_bonds[i][j + 1]){

							x4 = x_val[bonds[which_bonds[i][(j + 1) + (k + 1)] - 1][l + 1] - 1];
							y4 = y_val[bonds[which_bonds[i][(j + 1) + (k + 1)] - 1][l + 1] - 1];
							z4 = z_val[bonds[which_bonds[i][(j + 1) + (k + 1)] - 1][l + 1] - 1];

							a = which_bonds[i][0];
							b = which_bonds[i][j + 1];
							c = which_bonds[i][(j + 1) + (k + 1)];
							d = bonds[which_bonds[i][(j + 1) + (k + 1)] - 1][l + 1];

							dihedral_swap.push_back(dihedral(x1, x3, x2, x4, y1, y3, y2, y4, z1, z3, z2, z4));

							int_swap.push_back(a);
							int_swap.push_back(c);
							int_swap.push_back(b);
							int_swap.push_back(d);

							dihedral_int_swap.push_back(int_swap);
							int_swap.clear();
						}
					}
				}
			}
		}
		//saving values
		dist->push_back(dist_swap);
		anglevec->push_back(angle_swap);
		dihedralvec->push_back(dihedral_swap);
		dist_int.push_back(dist_int_swap);
		angle_int.push_back(angle_int_swap);
		dihedral_int.push_back(dihedral_int_swap);
	}
	//clean up
	involved_bonds->push_back(dist_int);
	involved_bonds->push_back(angle_int);
	involved_bonds->push_back(dihedral_int);
}//end get_values

//void get_values
//getting values without searching by just using safed atoms from the former used get_values-function
//
void neb::get_values(std::vector<double>& x_val, std::vector<double>& y_val, std::vector<double>& z_val, std::vector<std::vector<double> >* dist, std::vector<std::vector<double> >* anglevec, std::vector<std::vector<double> >* dihedralvec, std::vector<std::vector<std::vector<std::vector<int> > > >& involved_bonds, std::vector<std::vector<int> >& which_bonds){

	std::vector<double> dist_swap, angle_swap, dihedral_swap;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	size_t i, j;
	int a, b;
	bool same_atom;


	for (i = 0U; i<involved_bonds[DIST].size(); i++){
		for (j = 0U; j<involved_bonds[DIST][i].size(); j++){
			x1 = x_val[involved_bonds[DIST][i][j][0] - 1];
			y1 = y_val[involved_bonds[DIST][i][j][0] - 1];
			z1 = z_val[involved_bonds[DIST][i][j][0] - 1];

			x2 = x_val[involved_bonds[DIST][i][j][1] - 1];
			y2 = y_val[involved_bonds[DIST][i][j][1] - 1];
			z2 = z_val[involved_bonds[DIST][i][j][1] - 1];

			dist_swap.push_back(distance(x1, x2, y1, y2, z1, z2));
		}
		for (j = 0U; j<involved_bonds[ANGLE][i].size(); j++){

			x1 = x_val[involved_bonds[ANGLE][i][j][0] - 1];
			y1 = y_val[involved_bonds[ANGLE][i][j][0] - 1];
			z1 = z_val[involved_bonds[ANGLE][i][j][0] - 1];

			x2 = x_val[involved_bonds[ANGLE][i][j][1] - 1];
			y2 = y_val[involved_bonds[ANGLE][i][j][1] - 1];
			z2 = z_val[involved_bonds[ANGLE][i][j][1] - 1];

			x3 = x_val[involved_bonds[ANGLE][i][j][2] - 1];
			y3 = y_val[involved_bonds[ANGLE][i][j][2] - 1];
			z3 = z_val[involved_bonds[ANGLE][i][j][2] - 1];

			angle_swap.push_back(angle(x1, x2, x3, y1, y2, y3, z1, z2, z3));
		}
		for (j = 0U; j<involved_bonds[DIHEDRAL][i].size(); j++){

			x1 = x_val[involved_bonds[DIHEDRAL][i][j][0] - 1];
			y1 = y_val[involved_bonds[DIHEDRAL][i][j][0] - 1];
			z1 = z_val[involved_bonds[DIHEDRAL][i][j][0] - 1];

			x2 = x_val[involved_bonds[DIHEDRAL][i][j][1] - 1];
			y2 = y_val[involved_bonds[DIHEDRAL][i][j][1] - 1];
			z2 = z_val[involved_bonds[DIHEDRAL][i][j][1] - 1];

			x3 = x_val[involved_bonds[DIHEDRAL][i][j][2] - 1];
			y3 = y_val[involved_bonds[DIHEDRAL][i][j][2] - 1];
			z3 = z_val[involved_bonds[DIHEDRAL][i][j][2] - 1];

			x4 = x_val[involved_bonds[DIHEDRAL][i][j][3] - 1];
			y4 = y_val[involved_bonds[DIHEDRAL][i][j][3] - 1];
			z4 = z_val[involved_bonds[DIHEDRAL][i][j][3] - 1];

			//there is a need to look if the atom sits on the start atom if it is so
			//the transfer parameters must be donated by another order
			a = (int)involved_bonds[DIHEDRAL][i][j][0] - 1;
			b = (int)involved_bonds[DIHEDRAL][i][j][3];
			same_atom = false;
			for (size_t k = 1U; k<which_bonds[a].size(); k++){
				if (which_bonds[a][k] == b) same_atom = true;
			}
			if (same_atom){
				dihedral_swap.push_back(dihedral_same_atom(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4));
			}
			else{
				dihedral_swap.push_back(dihedral(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4));
			}
		}

		dist->push_back(dist_swap);
		anglevec->push_back(angle_swap);
		dihedralvec->push_back(dihedral_swap);
		dist_swap.clear();
		angle_swap.clear();
		dihedral_swap.clear();
	}
}

void neb::Sort(const size_t& imgs, std::vector<std::string>& name, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y, std::vector<std::vector<double> >& z, std::vector<std::vector<int> >& bonds){

	size_t i = 0U, j = 0U, k = 0U, iter = 0U, lenght, start;
	std::string string_swap;
	double double_swap;
	std::vector<int> vec_swap, remember_C, remember;

	for (i = 0U; i<bonds.size(); i++){

		if (name[i] == "C" && i != iter){
			string_swap = name[iter];
			name[iter] = name[i];
			name[i] = string_swap;

			for (size_t n = 0U; n<imgs; n++){

				double_swap = x[n][iter];
				x[n][iter] = x[n][i];
				x[n][i] = double_swap;
				double_swap = y[n][iter];
				y[n][iter] = y[n][i];
				y[n][i] = double_swap;
				double_swap = z[n][iter];
				z[n][iter] = z[n][i];
				z[n][i] = double_swap;

			}

			remember.push_back(iter + 1);
			remember_C.push_back(i + 1);
			vec_swap.clear();
			vec_swap = bonds[iter];
			bonds[iter] = bonds[i];
			bonds[i] = vec_swap;

			iter++;
		}
		else if (name[i] == "C" && i == iter) iter++;
	}

	for (i = 0U; i<bonds.size(); i++){
		for (j = 0U; j<bonds[i].size(); j++){
			for (k = 0U; k<remember.size(); k++){
				if (bonds[i][j] == remember[k]) bonds[i][j] = remember_C[k];
				else if (bonds[i][j] == remember_C[k]) bonds[i][j] = remember[k];
			}
		}
	}
	start = iter;
	lenght = bonds.size() - 1U;
	hetero_Sort(imgs, start, lenght, name, x, y, z, bonds);

	for (i = 0U; i<bonds.size(); i++){
		lenght = bonds[i].size() - 1U;
		start = 1;
		bond_Sort(bonds[i], start, lenght);
	}
}

void neb::bond_Sort(std::vector<int>& bonds, const size_t& start, const size_t& length){

	size_t i = start, j = length;
	int int_swap;
	int mid = bonds[(start + length) / 2];

	while (i <= j){
		while (bonds[i] < mid) i++;
		while (bonds[j] > mid) j--;
		if (i <= j){
			int_swap = bonds[i];
			bonds[i] = bonds[j];
			bonds[j] = int_swap;

			i++;
			j--;
		}
	};

	if (start < j) bond_Sort(bonds, start, j);
	if (i < length) bond_Sort(bonds, i, length);

}

void neb::bond_Sort(int* a, int* b, int* c, int* d){

	int int_swap;

	bool again = false;

	if (*a>*b){
		int_swap = *a;
		*a = *b;
		*b = int_swap;
		again = true;
	}
	else if (*b>*c){
		int_swap = *b;
		*b = *c;
		*c = int_swap;
		again = true;
	}
	else if (*c>*d){
		int_swap = *c;
		*c = *d;
		*d = int_swap;
		again = true;
	}
	if (again) bond_Sort(a, b, c, d);
}

void neb::bond_Sort(int* a, int* b, int* c){

	int int_swap;

	bool again = false;

	if (*a>*b){
		int_swap = *a;
		*a = *b;
		*b = int_swap;
		again = true;
	}
	else if (*b>*c){
		int_swap = *b;
		*b = *c;
		*c = int_swap;
		again = true;
	}
	if (again) bond_Sort(a, b, c);
}

void neb::bond_Sort(int* a, int* b){

	int int_swap;

	if (*a>*b){
		int_swap = *a;
		*a = *b;
		*b = int_swap;
	}
}

void neb::hetero_Sort(const size_t& imgs, const size_t& start, const size_t& length, std::vector<std::string>& name, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y, std::vector<std::vector<double> >& z, std::vector<std::vector<int> >& bonds){

	size_t a = 0U, b = 0U, c = 0U, i = start, j = length;
	double double_swap;
	std::vector<int> vector_swap, remember_i, remember_j;
	std::string string_swap;
	std::string mid = name[(start + length) / 2];

	while (i <= j){

		if (name[i] < "H"){
			while (name[i] < mid) i++;
			while (name[j] > mid) j--;
		}
		else if (name[i] >= "H"){
			while (name[i] > mid) i++;
			while (name[j] < mid) j--;
		}
		if (i <= j){
			string_swap = name[i];
			name[i] = name[j];
			name[j] = string_swap;

			for (size_t n = 0U; n<imgs; n++){

				double_swap = x[n][i];
				x[n][i] = x[n][j];
				x[n][j] = double_swap;
				double_swap = y[n][i];
				y[n][i] = y[n][j];
				y[n][j] = double_swap;
				double_swap = z[n][i];
				z[n][i] = z[n][j];
				z[n][j] = double_swap;

			}

			remember_i.push_back(i + 1);
			remember_j.push_back(j + 1);
			vector_swap.clear();
			vector_swap = bonds[i];
			bonds[i] = bonds[j];
			bonds[j] = vector_swap;

			i++;
			j--;
		}
	};

	for (a = 0U; a<bonds.size(); a++){
		for (b = 0U; b<bonds[a].size(); b++){
			for (c = 0U; c<remember_i.size(); c++){
				if (bonds[a][b] == remember_i[c]) bonds[a][b] = remember_j[c];
				else if (bonds[a][b] == remember_j[c]) bonds[a][b] = remember_i[c];
			}
		}
	}

	if (start < j) hetero_Sort(imgs, start, j, name, x, y, z, bonds);
	if (i < length) hetero_Sort(imgs, i, length, name, x, y, z, bonds);

}

size_t neb::factorial(const size_t& fac) {

	if (fac <= 1) {
		return 1;
	}

	return factorial(fac - 1) * fac;
}

//double distance
//calculating distance
//
double inline neb::distance(const double& x1, const double& x2, const double& y1, const double& y2, const double& z1, const double& z2){

	double x, y, z;

	x = x2 - x1;
	y = y2 - y1;
	z = z2 - z1;

	return sqrt(x*x + y*y + z*z);

}//end distance

//double angle
//calculating the vectors from 1 to 2 and from 1 to three and calculating the angle
//
double inline neb::angle(const double& x1, const double& x2, const double& x3, const double& y1, const double& y2, const double& y3, const double& z1, const double& z2, const double& z3){

	double re, x, y, z, a, b, c;

	x = x2 - x1;
	y = y2 - y1;
	z = z2 - z1;

	a = x3 - x1;
	b = y3 - y1;
	c = z3 - z1;

	re = (x*a + y*b + z*c) / (sqrt(x*x + y*y + z*z)*sqrt(a*a + b*b + c*c));
	re = acos(re);

	return re*(180.0f / 3.1415926535f);

}//end angle

//double dihedral
//generating the normal vectors of both planes mounted by 1,2 and 3 and 1,2 and 4 to get their angle
//which is equal to the dihedral angle and looking whether the point is +/- rotated
double inline neb::dihedral(const double& x1, const double& x2, const double& x3, const double& x4, const double& y1, const double& y2, const double& y3, const double& y4, const double& z1, const double& z2, const double& z3, const double& z4){

	double re, x, y, z, a, b, c, a_norm, b_norm, c_norm, length_normal, length_reflection, lambda, x_1, y_1, z_1, x_2, y_2, z_2, d;

	a = (y4 - y2)*(z1 - z2) - (z4 - z2)*(y1 - y2);
	b = (z4 - z2)*(x1 - x2) - (x4 - x2)*(z1 - z2);
	c = (x4 - x2)*(y1 - y2) - (y4 - y2)*(x1 - x2);

	x = (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1);
	y = (z2 - z1)*(x3 - x1) - (x2 - x1)*(z3 - z1);
	z = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1);

	//to see if the dihedral has a negative or positive torsion we need to reflect it
	//get d of plain_1 (n_1*x+n_2*y+n_3*z+d=0)
	d = a*x2 + b*y2 + c*z2;

	//normalize normal vector
	a_norm = a / (sqrt(a*a + b*b + c*c));
	b_norm = b / (sqrt(a*a + b*b + c*c));
	c_norm = c / (sqrt(a*a + b*b + c*c));

	//get lambda
	lambda = (a*x3 + b*y3 + c*z3 - d) / ((-1)*(a*a_norm) - (b*b_norm) - (c*c_norm));

	//reflect x3,y3,z3
	x_1 = x3 + lambda*a_norm;
	y_1 = y3 + lambda*b_norm;
	z_1 = z3 + lambda*c_norm;

	x_2 = x_1 - x3;
	y_2 = y_1 - y3;
	z_2 = z_1 - z3;

	x_1 += x_2;
	y_1 += y_2;
	z_1 += z_2;

	//construct a point lying on the "normal" side of the plain
	//and check whether distance to the normal point is farther or
	//its reflection

	length_normal = sqrt(((x1 + a_norm) - x3)*((x1 + a_norm) - x3) + ((y1 + b_norm) - y3)*((y1 + b_norm) - y3) + ((z1 + c_norm) - z3)*((z1 + c_norm) - z3));
	length_reflection = sqrt(((x1 + a_norm) - x_1)*((x1 + a_norm) - x_1) + ((y1 + b_norm) - y_1)*((y1 + b_norm) - y_1) + ((z1 + c_norm) - z_1)*((z1 + c_norm) - z_1));

	re = (x*a + y*b + z*c) / (sqrt(x*x + y*y + z*z)*sqrt(a*a + b*b + c*c));

	re = acos(re);

	if (length_normal<length_reflection){
		return (re*(180.0f / 3.1415926535f));
	}
	else{
		return -(re*(180.0f / 3.1415926535f));
	}

}//end dihedral

//double dihedral_same_atom
//same as above just for dihedrals where both points are at the start atom
//
double inline neb::dihedral_same_atom(const double& x1, const double& x2, const double& x3, const double& x4, const double& y1, const double& y2, const double& y3, const double& y4, const double& z1, const double& z2, const double& z3, const double& z4){

	double re, x, y, z, a, b, c, a_norm, b_norm, c_norm, length_normal, length_reflection, lambda, x_1, y_1, z_1, x_2, y_2, z_2, d;

	//first normal vector
	a = (y4 - y1)*(z2 - z1) - (z4 - z1)*(y2 - y1);
	b = (z4 - z1)*(x2 - x1) - (x4 - x1)*(z2 - z1);
	c = (x4 - x1)*(y2 - y1) - (y4 - y1)*(x2 - x1);

	//secund normal vector
	x = (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1);
	y = (z2 - z1)*(x3 - x1) - (x2 - x1)*(z3 - z1);
	z = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1);

	//to see if the dihedral has a negative or positive torsion we need to reflect it
	//get d of plain_1 (n_1*x+n_2*y+n_3*z+d=0)
	d = a*x1 + b*y1 + c*z1;

	//normalize normal vector
	a_norm = a / (sqrt(a*a + b*b + c*c));
	b_norm = b / (sqrt(a*a + b*b + c*c));
	c_norm = c / (sqrt(a*a + b*b + c*c));

	//get lambda
	lambda = (a*x3 + b*y3 + c*z3 - d) / ((-1)*(a*a_norm) - (b*b_norm) - (c*c_norm));

	//reflect x3,y3,z3
	x_1 = x3 + lambda*a_norm;
	y_1 = y3 + lambda*b_norm;
	z_1 = z3 + lambda*c_norm;

	x_2 = x_1 - x3;
	y_2 = y_1 - y3;
	z_2 = z_1 - z3;

	x_1 += x_2;
	y_1 += y_2;
	z_1 += z_2;

	//construct a point lying on the "normal" side of the plain
	//and check whether distance to the normal point is farther or
	//its reflection

	length_normal = sqrt(((x1 + a_norm) - x3)*((x1 + a_norm) - x3) + ((y1 + b_norm) - y3)*((y1 + b_norm) - y3) + ((z1 + c_norm) - z3)*((z1 + c_norm) - z3));
	length_reflection = sqrt(((x1 + a_norm) - x_1)*((x1 + a_norm) - x_1) + ((y1 + b_norm) - y_1)*((y1 + b_norm) - y_1) + ((z1 + c_norm) - z_1)*((z1 + c_norm) - z_1));

	re = (x*a + y*b + z*c) / (sqrt(x*x + y*y + z*z)*sqrt(a*a + b*b + c*c));

	re = acos(re);

	if (length_normal<length_reflection){
		return (re*(180.0f / 3.1415926535f));
	}
	else{
		return -(re*(180.0f / 3.1415926535f));
	}

}//end dihedral_same_atom



/***************************************
****************************************
****                                ****
****    *******                *    ****
****    **                     *    ****
****    **                     *    ****
****    *******   * ***     ** *    ****
****    **        **   *   *   *    ****
****    **        *    *   *   *    ****
****    *******   *    *    ** *    ****
****                                ****
****************************************
***************************************/

