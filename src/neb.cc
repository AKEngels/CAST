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
#include <iterator>  
#include <algorithm>  
#include <utility>
#include <unordered_set>
#include "ls.h"
#include "lbfgs.h"
#include "coords.h"
#include "optimization_global.h"
#include "matop.h"
#include "scon_mathmatrix.h"


/**
* NEB constructor
*/

neb::neb(coords::Coordinates *cptr):_KT_(1 / (0.0019872966*Config::get().neb.TEMPERATURE))
{
  cPtr = cptr;
  reversed = true;
  CIMaximum = 0;
  num_images = Config::get().neb.IMAGES;
  ClimbingImage = Config::get().neb.CLIMBING;
  springconstant = Config::get().neb.SPRINGCONSTANT;
}

/**
* INITIALIZING function I. standard NEB call from main with Inputstructures defined in the INPUTFILE
*/

void neb::preprocess(ptrdiff_t &count)
{
  std::vector<size_t> image_remember;
  std::vector<std::vector<size_t> > atoms_remember;
  N = cPtr->size();
  ts = false;
  N = cPtr->size();
  imagi.clear();
  image_ini.clear();
  energies.clear();
  tau.clear();
  images.clear();
  images_initial.clear();
  ts_pathstruc.resize(num_images);
  image_ini.resize(num_images);
  images.resize(N);
  tau.resize(num_images);
  tempimage_final.resize(num_images);
  tempimage_ini.resize(num_images);
  energies.resize(num_images);
  imagi.resize(num_images);
  initial();
  final();
  if (Config::get().neb.IDPP)
  {
    if (!Config::get().neb.INTERNAL_INTERPOLATION)
    {
      idpp_prep();
      create_cartesian_interpolation();
      if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
      else run(count);
    }
    else
    {
      idpp_prep();
      create_internal_interpolation(imagi);
      if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
      else run(count);
    }
  }
  else
  {
    if (!Config::get().neb.INTERNAL_INTERPOLATION)
    {
      create_cartesian_interpolation();
      if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
      else run(count);
    }
    else
    {
      create_internal_interpolation(imagi);
      if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
      else run(count);
    }
  }
}

/**
* INITIALIZING function II. Pathopt NEB call with structures defined in Pathopt / forward ordered
*/

void neb::preprocess(ptrdiff_t &image, ptrdiff_t &count, const coords::Representation_3D &start, const coords::Representation_3D &fi, const std::vector <double> &ts_energy, const std::vector <double> &min_energy, bool reverse, const coords::Representation_3D &ts_path)
{
  std::vector<size_t> image_remember;
  std::vector<std::vector<size_t> > atoms_remember;
  N = cPtr->size();
  ts = true;
  reversed = reverse;
  num_images = image;
  imagi.clear();
  image_ini.clear();
  energies.clear();
  tau.clear();
  images.clear();
  ts_energies.resize(ts_energy.size());
  min_energies.resize(min_energy.size());
  imagi.resize(num_images);
  ts_pathstruc.resize(num_images);
  image_ini.resize(num_images);
  images_initial.clear();
  images.resize(N);
  energies.resize(num_images);
  tau.resize(num_images);
  ts_energies = ts_energy;
  min_energies = min_energy;
  ts_pathstruc = ts_path;
  initial(start);
  final(fi);
  if (!Config::get().neb.INTERNAL_INTERPOLATION)
  {
    create_cartesian_interpolation();
    images.clear();
    if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
    else run(count);
  }
  else
  {
    create_internal_interpolation(imagi);
    images.clear();
    if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
    else run(count);
  }
  images.clear();
  if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
  else run(count);
}

/**
* INITIALIZING function III. Pathopt NEB call with structures defined in Pathopt / reverse ordered
*/

void neb::preprocess(ptrdiff_t &file, ptrdiff_t &image, ptrdiff_t &count, const coords::Representation_3D &start, const coords::Representation_3D &fi, bool reverse)
{
  std::vector<size_t> image_remember;
  std::vector<std::vector<size_t> > atoms_remember;
  N = cPtr->size();
  this->cPtr->mult_struc_counter = file;
  CIMaximum = 0;
  ts = false;
  reversed = reverse;
  num_images = image;
  imagi.clear();
  image_ini.clear();
  energies.clear();
  tau.clear();
  images.clear();
  imagi.resize(num_images);
  ts_pathstruc.resize(num_images);
  image_ini.resize(num_images);
  images_initial.clear();
  tempimage_final.resize(N);
  tempimage_ini.resize(N);
  images.resize(N);
  tau.resize(num_images);
  energies.resize(num_images);
  initial(start);
  final(fi);
  if (!Config::get().neb.INTERNAL_INTERPOLATION)
  {
    create_cartesian_interpolation();
    if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
    else run(count);
  }
  else
  {
    create_internal_interpolation(imagi);
    if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
    else run(count);
  }
}

/**
* defining start structures from INPUT
*/
void neb::preprocess(std::vector<coords::Representation_3D> & ini_path, ptrdiff_t &count) {
	std::vector<size_t> image_remember;
	std::vector<std::vector<size_t> > atoms_remember;
	coords::Representation_3D ini_path_x;
	N = cPtr->size();
	ts = false;
	num_images = ini_path.size();
	imagi.clear();
	image_ini.clear();
	energies.clear();
	tau.clear();
	images.clear();
	images_initial.clear();
	ts_pathstruc.resize(num_images);
	image_ini.resize(num_images);
	images.resize(N);
	tau.resize(num_images);
	tempimage_final.resize(num_images);
	tempimage_ini.resize(num_images);
	energies.resize(num_images);
	imagi.resize(num_images);
	ini_path_x.resize(N);
	ini_path_x=ini_path[0];
	initial(ini_path_x);
	ini_path_x = ini_path[num_images-1];
	final(ini_path_x);
	create_ini_path(ini_path);
	images.clear();
	if (Config::get().neb.CONSTRAINT_GLOBAL)run(count, image_remember, atoms_remember);
	else run(count);

}

void neb::initial(void)

{
  imagi[0] = cPtr->xyz();
}

/**
* Reading second structure for double ended optimization from special NEB INPUT definition
*/
void neb::final(void)
{
  std::ifstream final(Config::get().neb.FINAL_STRUCTURE.c_str());
  if (!final)
  {
    throw std::runtime_error("Final NEB Structure specified was not found. Aborting.");
  }

  std::string buffer;
  getline(final, buffer);
  size_t number;
  char atom[3];
  for (size_t i = 0; i < N; i++)
  {
    getline(final, buffer);
    std::istringstream line_coord(buffer);
    line_coord >> number >> atom >> images[i].x() >> images[i].y() >> images[i].z();
    imagi[num_images - 1].push_back(images[i]);
  }

  
}

/**
* Reading second structure for double ended optimization from special NEB INPUT definition
*/
void neb::final_align(void)
{
  std::ifstream final(Config::get().neb.FINAL_STRUCTURE.c_str());
  if (!final)
  {
    throw std::runtime_error("Final NEB Structure specified was not found. Aborting.");
  }
  std::string buffer;
  getline(final, buffer);
  size_t number;
  char atom[3];
  coords::Coordinates final_struct, initial_struct;
  for (size_t i = 0; i < N; i++)
  {
    getline(final, buffer);
    std::istringstream line_coord(buffer);
    line_coord >> number >> atom >> images[i].x() >> images[i].y() >> images[i].z();
  }
  initial_struct.init_in(cPtr->atoms(), coords::PES_Point(imagi[0]), true);
  final_struct.init_in(cPtr->atoms(), coords::PES_Point(images), true);
  //std::cout << initial_struct.size() << '\n' << final_struct.size() << '\n';
  for (size_t i = 0; i < this->cPtr->size(); ++i)
  {
    scon::align_cog(final_struct.xyz(i), final_struct.center_of_geometry(), initial_struct.center_of_geometry());
  }
  matop::align::kabschAlignment(final_struct, initial_struct, false);
  imagi[num_images - 1] = final_struct.xyz();
}

/**
* defining start structure from PATHOPT call
*/

void neb::initial(const coords::Representation_3D &start)

{
  for (size_t i = 0; i < cPtr->size(); i++)
  {
    imagi[0].push_back(start[i]);
  }
}

/**
* defining end structure from Pathopt call
*/

void neb::final(const coords::Representation_3D &fi)
{
  for (size_t i = 0; i < cPtr->size(); i++)
  {
    imagi[num_images - 1].push_back(fi[i]);
  }
}

/**
* linear interpolation
*/
void neb::create_cartesian_interpolation()
{

  tempimage_ini = imagi[0];
  tempimage_final = imagi[num_images - 1];
  std::ostringstream name;
  name << "IMAGES_INI" << cPtr->mult_struc_counter << ".dat";
  cPtr->set_xyz(imagi[0]);
  tempimage_ini = cPtr->xyz();
  cPtr->set_xyz(imagi[num_images - 1]);
  tempimage_final = cPtr->xyz();
  for (size_t j = 1; j < (num_images - 1); j++) {

    double diff = (double)j / num_images;

    for (size_t i = 0; i < this->cPtr->size(); i++) {

      images[i].x() = tempimage_ini[i].x() + diff * (tempimage_final[i].x() - tempimage_ini[i].x());
      images[i].y() = tempimage_ini[i].y() + diff * (tempimage_final[i].y() - tempimage_ini[i].y());
      images[i].z() = tempimage_ini[i].z() + diff * (tempimage_final[i].z() - tempimage_ini[i].z());
	  

      imagi[j].push_back(images[i]);
      image_ini[j].push_back(images[i]);
      images_initial.push_back(images[i]);

    }
	
	
  }

  std::ostringstream na;
  na << "IMAGES_START" << this->cPtr->mult_struc_counter << ".arc";
  std::ptrdiff_t s{ 0 };
  print(na.str(), imagi, s);

  if (Config::get().neb.INT_PATH) calc_shift();
}


void neb::create_ini_path(const std::vector<coords::Representation_3D> &ini)
{
	for (size_t j = 1; j < (num_images - 1); j++)
	{
		for (size_t i = 0; i < this->cPtr->size(); i++) 
		{
			images[i].x() = ini[j][i].x();
			images[i].y() = ini[j][i].y();
			images[i].z() = ini[j][i].z();

			imagi[j].push_back(images[i]);
			image_ini[j].push_back(images[i]);
			images_initial.push_back(images[i]);
		}

	}
	std::ostringstream na;
	na << "IMAGES_START" << this->cPtr->mult_struc_counter << ".arc";
	std::ptrdiff_t s{ 0 };
	print(na.str(), imagi, s);

}

/**
* NEB run function for standard execution
*/

void neb::run(ptrdiff_t &count)
{

  calc_tau();
  for (auto tx : tau[1]) tau[0].push_back(tx);
  for (auto te : tau[num_images - 2])tau[num_images - 1].push_back(te);
  if (Config::get().neb.MAXFLUX) lbfgs_maxflux();
  else lbfgs();
  tau[0].clear();
  tau[num_images - 1].clear();
  for (auto tx : tau[1]) tau[0].push_back(tx);
  for (auto te : tau[num_images - 2])tau[num_images - 1].push_back(te);
  opt_io(count);

}
void neb::run(ptrdiff_t &count, std::vector<size_t>&, std::vector<std::vector<size_t> >& atoms_remember)
{
  lbfgs();
  opt_io(count);
  internal_execute(imagi, atoms_remember);
  //no_torsion(image_remember[DIHEDRAL],atoms_remember[DIHEDRAL]);
  opt_internals(count, atoms_remember);
  //internal_execute(imagi,image_remember,atoms_remember);
  for (size_t t = 0; t < num_images; t++) tau[t].clear();
  calc_tau();
}

/**
* Calculation of energies which stem from the actual optimization step
*/

void neb::get_energies(void)
{

  for (size_t i = 0; i < num_images; i++)
  {
    cPtr->set_xyz(imagi[i]);
    energies[i] = cPtr->g();
  }

  if (Config::get().general.verbosity > 4) std::cout << "maximum energy image is: " << CIMaximum << "\n";
}

/**
* Calculation of the band defining vectors tau --> standard / improved tangent estimate
*/

void neb::calc_tau(void)
{

  double Em1{ 0.0 }, Ep1{ 0.0 }, Emax{ 0.0 }, Emin{ 0.0 };
  get_energies();
  for (size_t i = 1; i < num_images - 1; i++) {
	if (energies[i] > energies[i - 1]) CIMaximum = i;
	/// standard tangent estimate
    if (Config::get().neb.TAU == false)
    {
      for (size_t j = 0; j < N; j++)
	  {
        images[j].x() = imagi[i][j].x() - imagi[i - 1][j].x() / abs(imagi[i][j].x() - imagi[i - 1][j].x()) + (imagi[i + 1][j].x() - imagi[i][j].x()) / abs(imagi[i + 1][j].x() - imagi[i][j].x());
        images[j].y() = imagi[i][j].y() - imagi[i - 1][j].y() / abs(imagi[i][j].y() - imagi[i - 1][j].y()) + (imagi[i + 1][j].y() - imagi[i][j].y()) / abs(imagi[i + 1][j].y() - imagi[i][j].y());
        images[j].z() = imagi[i][j].z() - imagi[i - 1][j].z() / abs(imagi[i][j].z() - imagi[i - 1][j].z()) + (imagi[i + 1][j].z() - imagi[i][j].z()) / abs(imagi[i + 1][j].z() - imagi[i][j].z());
        if (images[j].x() - images[j].x() != 0) images[j].x() = 0.0;
        if (images[j].y() - images[j].y() != 0) images[j].y() = 0.0;
        if (images[j].z() - images[j].z() != 0) images[j].z() = 0.0;
		tau[i].push_back(images[j]);		
      }
    }
	/// improved tangent estimate
    else
    {

      EnergyPml = 0.0;
      EnergyPpl = 0.0;

      if (energies[i - 1] > energies[i]) EnergyPml = energies[i - 1];
      else EnergyPml = energies[i];
      if (energies[i + 1] > energies[i]) EnergyPpl = energies[i + 1];
      else EnergyPpl = energies[i];

      if (EnergyPml != EnergyPml) 
	  {
        if (EnergyPml > EnergyPpl) 
		{

          for (size_t j = 0; j < N; j++) {
            images[j].x() = (imagi[i][j].x() - imagi[i - 1][j].x());
            images[j].y() = (imagi[i][j].y() - imagi[i - 1][j].y());
            images[j].z() = (imagi[i][j].z() - imagi[i - 1][j].z());
            if (images[j].x() != images[j].x()) images[j].x() = 0.0;
            if (images[j].y() != images[j].y()) images[j].y() = 0.0;
            if (images[j].z() != images[j].z()) images[j].z() = 0.0;
            tau[i].push_back(images[j]);
          }
        }
        else {
          for (size_t j = 0; j < N; j++) {

          

            images[j].x() = (imagi[i + 1][j].x() - imagi[i][j].x());
            images[j].y() = (imagi[i + 1][j].y() - imagi[i][j].y());
            images[j].z() = (imagi[i + 1][j].z() - imagi[i][j].z());
            if (images[j].x() - images[j].x() != 0) images[j].x() = 0.0;
            if (images[j].y() - images[j].y() != 0) images[j].y() = 0.0;
            if (images[j].z() - images[j].z() != 0) images[j].z() = 0.0;
            tau[i].push_back(images[j]);

          }
        }
      }
      else 
	  {
        Em1 = energies[i - 1] - energies[i];
        Ep1 = energies[i + 1] - energies[i];

        Emin = std::min(abs(Ep1), abs(Em1));
        Emax = std::max(abs(Ep1), abs(Em1));

        if (Em1 > Ep1) 
		{
          for (size_t j = 0; j < N; j++) {
            images[j].x() = (imagi[i + 1][j].x() - imagi[i][j].x()) * Emin + (imagi[i][j].x() - imagi[i - 1][j].x()) * Emax;
            images[j].y() = (imagi[i + 1][j].y() - imagi[i][j].y()) * Emin + (imagi[i][j].y() - imagi[i - 1][j].y()) * Emax;
            images[j].z() = (imagi[i + 1][j].z() - imagi[i][j].z()) * Emin + (imagi[i][j].z() - imagi[i - 1][j].z()) * Emax;
			if (images[j].x() - images[j].x() != 0) images[j].x() = 0.0;
			if (images[j].y() - images[j].y() != 0) images[j].y() = 0.0;
			if (images[j].z() - images[j].z() != 0) images[j].z() = 0.0;
			tau[i].push_back(images[j]);
          }

        }


        else {
          for (size_t j = 0; j < N; j++) {

            
            images[j].x() = (imagi[i + 1][j].x() - imagi[i][j].x()) * Emax + (imagi[i][j].x() - imagi[i - 1][j].x()) * Emin;
            images[j].y() = (imagi[i + 1][j].y() - imagi[i][j].y()) * Emax + (imagi[i][j].y() - imagi[i - 1][j].y()) * Emin;
            images[j].z() = (imagi[i + 1][j].z() - imagi[i][j].z()) * Emax + (imagi[i][j].z() - imagi[i - 1][j].z()) * Emin;
            if (images[j].x() - images[j].x() != 0) images[j].x() = 0.0;
            if (images[j].y() - images[j].y() != 0) images[j].y() = 0.0;
            if (images[j].z() - images[j].z() != 0) images[j].z() = 0.0;
            tau[i].push_back(images[j]);				
          }
        }

      }

    }
	for (size_t j = 0; j < N; j++) 
	{
		if (len(tau[i][j]) != 0.0)
		{
			tau[i][j].x() /= len(tau[i][j]);
			tau[i][j].y() /= len(tau[i][j]);
			tau[i][j].z() /= len(tau[i][j]);
		}
	}

  }

}

/**
* IDPP start  
*/
void neb::idpp_prep()
{
  start_structure = cPtr->xyz();
  final_structure = imagi.back();
  bond_st = bond_dir(start_structure);
  eu_st = euclid_dist<double>(start_structure);
  eu_fi = euclid_dist<double>(final_structure);
}

template<typename ContainerIt, typename Result>
void neb::concatenate(ContainerIt st, ContainerIt fi, Result res)
{
  while (st != fi)
  {
    res = std::move(st->begin(), st->end(), res);
    ++st;
  }
}

template<typename T>
neb::dist_T<T> neb::euclid_dist(coords::Representation_3D const& struc)
{
  dist_T<T> result;
  std::vector<T> dist;
  for (auto& i : struc)
  {
    dist.clear();
    for (auto& s : struc)
    {
      auto f{ (i - s)*(i - s) };
      auto d(std::sqrt(f.x() + f.y() + f.z()));
      if (d != 0.0) dist.push_back(d);
      else continue;
    }
    result.push_back(dist);
  }
  return result;
}

std::vector<coords::Representation_3D> neb::bond_dir(coords::Representation_3D const& struc)
{
  coords::Representation_3D direction;
  std::vector<coords::Representation_3D> all_directions;
  for (auto& i : struc)
  {
    direction.clear();
    for (auto& s : struc)
    {
      auto f{ s - i };
      auto d{ f.x() + f.y() + f.z() };
      if (d != 0.0) direction.push_back(f);
      else continue;
    }
    all_directions.push_back(direction);
  }
  return all_directions;
}

template<typename T, typename Euclid>
std::vector<T> neb::interpol_dist(Euclid& st, Euclid& fi, size_t const img_no)
{
  /*if (!std::all_of(st.cbegin(), st.cend(),
  [&](auto const& x) {return x.size() == st.front().size(); }))
  throw ("Size error!");

  if (!std::all_of(fi.cbegin(), fi.cend(),
  [&](auto const& x) {return x.size() == fi.front().size(); }))
  throw ("Size error!");*/

  std::vector<T> di_st, di_fi;
  di_st.reserve(st.size()*(st.size() - 1));
  di_fi.reserve(fi.size()*(fi.size() - 1));
  concatenate(st.begin(), st.end(), std::back_inserter(di_st));
  concatenate(fi.begin(), fi.end(), std::back_inserter(di_fi));
  auto st_it{ di_st.cbegin() };
  auto fi_it{ di_fi.cbegin() };
  std::vector<T> result;
  for (; st_it != di_st.cend() && fi_it != di_fi.cend(); ++st_it, ++fi_it)
  {
    auto s{ *st_it + ((double)img_no / (num_images - 1)) * (*fi_it - *st_it) };
    result.push_back(s);
  }
  return result;
}

template<typename T, typename Euclid>
coords::Representation_3D neb::idpp_gradients(std::vector<coords::Representation_3D> const& bond_st,
	Euclid& start, Euclid& st, Euclid& fi, size_t const im_no)
{
  auto interp_di{ interpol_dist<double, dist_T<double>>(start, fi, im_no) };
  coords::Representation_3D all_grad;
  auto bond_it{ bond_st.cbegin() };
  auto euclid_it{ st.cbegin() };
  auto interp_it{ interp_di.cbegin() };
  std::vector<scon::c3<coords::float_type>> dummy;
  for (; bond_it != bond_st.cend() && euclid_it != st.cend()
    && interp_it != interp_di.cend(); ++bond_it, ++euclid_it, ++interp_it)
  {
    dummy.clear();
    scon::c3<coords::float_type> dummy_sum{ 0.0, 0.0, 0.0 };
    std::vector<scon::c3<coords::float_type>> bonds{ *bond_it };
    std::vector<T> eu_dist{ *euclid_it };
    std::vector<T> interp_dist{ *interp_it };
    auto b_it{ bonds.cbegin() };
    auto eu_it{ eu_dist.cbegin() };
    auto ip_it{ interp_dist.cbegin() };
    for (; b_it != bonds.cend() && eu_it != eu_dist.cend()
      && ip_it != interp_dist.cend(); ++b_it, ++eu_it, ++ip_it)
    {
      auto s{ (*b_it / *eu_it) * (((*ip_it - 2.0 * *eu_it)
        * (*ip_it - *eu_it)) / std::pow(*eu_it, 5.0)) };
      dummy.push_back(s);
    }
    for (auto i : dummy) { dummy_sum += i; }
    scon::c3<coords::float_type> dummy_sum_mult{ (-2.0)*dummy_sum.x(),
      (-2.0)*dummy_sum.y(), (-2.0)*dummy_sum.z() };
    all_grad.push_back(dummy_sum_mult);
  }
  return all_grad;
}
//IDPP end  

/**
* I/O of optimized structures and energies
*/
void neb::opt_io(ptrdiff_t &count)
{

  std::ostringstream energies_out, name;
  if (Config::get().general.task == 8) energies_out << "ENERGIES_COMPLETE_" << this->cPtr->mult_struc_counter << ".dat";
  else energies_out << "ENERGIES_COMPLETE_" << this->cPtr->mult_struc_counter <<"_"<< count << ".dat";
  std::string temp{ energies_out.str() };
  std::fstream off(temp.c_str(), std::ios::app);
  if (reversed == true)
  {
	  if (Config::get().general.task == 8) name << "IMAGES_FINAL" << this->cPtr->mult_struc_counter << ".arc";
	  else name << "IMAGES_FINAL" << this->cPtr->mult_struc_counter << "_" << count << ".arc";
    if (ts == true)
    {
      off << "TS          " << std::right << std::fixed << std::setprecision(6) << ts_energies[count] << '\n';
      off << "MIN         " << std::right << std::fixed << std::setprecision(6) << min_energies[count] << '\n';
      printmono(name.str(), imagi[0], count);
    }
    off << "ENERGIE:    " << std::right << std::fixed << std::setprecision(6) << energies[0] << "   START\n";
    for (size_t imagecount = 1; imagecount < num_images - 1; imagecount++)
    {
      off << "ENERGIE:    " << std::right << std::fixed << std::setprecision(6) << energies[imagecount] << "     IMAGE:   " << imagecount << "   GLOBAL-COUNT:  " << count << '\n';
    }
    off << "ENERGIE:    " << std::right << std::fixed << std::setprecision(6) <<energies[num_images - 1] << "   FINAL\n";
    print(name.str(), imagi, count);
  }
  else
  {
	  if (Config::get().general.task == 8) name << "IMAGES_FINAL" << this->cPtr->mult_struc_counter << ".arc";
	  else name << "IMAGES_FINAL" << this->cPtr->mult_struc_counter << "_" << count << ".arc";
    off << "ENERGIE:    " << std::right << std::fixed << std::setprecision(6) << energies[num_images - 1] << "   FINAL\n";
    for (ptrdiff_t imagecount = num_images - 2; imagecount >= 1; imagecount--)
    {
      cPtr->set_xyz(imagi[imagecount]);
      off << "ENERGIE:    " << std::right << std::fixed << std::setprecision(6) << energies_NEB[imagecount] << "     IMAGE:   " << imagecount << "   GLOBAL-COUNT:  " << count << '\n';
    }

    print_rev(name.str(), imagi, count);
    if (ts == true)
    {
      off << "MIN         " << std::right << std::fixed << std::setprecision(6) << min_energies[count] << '\n';
      off << "TS          " << std::right << std::fixed << std::setprecision(6) << ts_energies[count] << '\n';
      printmono(name.str(), ts_pathstruc, count);
    }
  }
}

void neb::print(std::string const &name, std::vector <coords::Representation_3D> &print, ptrdiff_t const &count)
{

  std::ofstream out(name.c_str(), std::ios::app);
  std::string temp;

  for (size_t i(0U); i < num_images; i++) {
    out << "     " << N << " IMAGE: " << i + 1 << "  global counter:  " << count << '\n';
    for (size_t j = 0; j < N; j++) {

      out << std::right << std::setw(6) << j + 1;
      out << std::left << "  " << std::setw(12) << cPtr->atoms(j).symbol().substr(0U, 2U);
      out << std::fixed << std::right << std::setprecision(6) << std::setw(13) << print[i][j].x();
      out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[i][j].y();
      out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[i][j].z();
      out << std::right << std::setw(6) << cPtr->atoms(j).energy_type();
      size_t const bSize(cPtr->atoms(j).bonds().size());
      for (size_t n(0U); n < bSize; ++n)
      {
        out << std::right << std::setw(6) << cPtr->atoms(j).bonds()[n] + 1U;
      }
      out << '\n';
    }
  }
  out.close();


}

void neb::print_rev(std::string const &name, std::vector <coords::Representation_3D> &print, ptrdiff_t &count)
{

  std::ofstream out(name.c_str(), std::ios::app);
  std::string temp;

  for (ptrdiff_t i = num_images - 1; i >= 0; i--) {
    out << "     " << N << " IMAGE: " << i + 1 << "  global counter:  " << count << '\n';
    for (size_t j = 0; j < N; j++) {

      out << std::right << std::setw(6) << j + 1;
      out << std::left << "  " << std::setw(12) << cPtr->atoms(j).symbol().substr(0U, 2U);
      out << std::fixed << std::right << std::setprecision(6) << std::setw(13) << print[i][j].x();
      out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[i][j].y();
      out << std::fixed << std::right << std::setprecision(6) << std::setw(12) << print[i][j].z();
      out << std::right << std::setw(6) << cPtr->atoms(j).energy_type();
      size_t const bSize(cPtr->atoms(j).bonds().size());
      for (size_t n(0U); n < bSize; ++n)
      {
        out << std::right << std::setw(6) << cPtr->atoms(j).bonds()[n] + 1U;
      }
      out << '\n';
    }
  }
  out.close();


}


void neb::printmono(std::string const &name, coords::Representation_3D &print, ptrdiff_t const &count)
{

  std::ofstream out(name.c_str(), std::ios::app);
  std::string temp;


  out << "     " << N << " IMAGE: " << "  global counter:  " << count << '\n';
  for (size_t j = 0; j < N; j++) {

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




double neb::lbfgs()
{

  using namespace  optimization::local;
  energies_NEB.resize(num_images);
  //typedef coords::Container<scon::c3<float>> nc3_type;
  // Create optimizer
  coords::Cartesian_Point g;
  auto optimizer = make_lbfgs(
    make_more_thuente(GradCallBack(*this))
  );
  //optimizer.ls.config.ignore_callback_stop = true;
  // Create Point
  using op_type = decltype(optimizer);
  op_type::point_type x(scon::explicit_transform<op_type::rep_type>(images_initial));
  // Optimize point
  optimizer.config.max_iterations = Config::get().optimization.local.bfgs.maxstep;
  //optimizer.config.max_iterations = Config::get().neb.LBFGS_IT;
  optimizer.config.epsilon = (float)Config::get().optimization.local.bfgs.grad;
  optimizer(x);

  if (Config::get().general.verbosity > 4)
  {
    std::cout << "Optimization done (status " << optimizer.state() << "). Evaluations:" << optimizer.iter() << '\n';
  }

  return cPtr->g();
}

double neb::lbfgs_maxflux()
{

	using namespace  optimization::local;
	energies_NEB.resize(num_images);
	//typedef coords::Container<scon::c3<float>> nc3_type;
	// Create optimizer
	coords::Cartesian_Point g;
	auto optimizer = make_lbfgs(
		make_more_thuente(GradCallBackMaxFlux(*this))
	);
	//optimizer.ls.config.ignore_callback_stop = true;
	// Create Point
	using op_type = decltype(optimizer);
	op_type::point_type x(scon::explicit_transform<op_type::rep_type>(images_initial));
	// Optimize point
	optimizer.config.max_iterations = Config::get().optimization.local.bfgs.maxstep;
	//optimizer.config.max_iterations = Config::get().neb.LBFGS_IT;
	optimizer.config.epsilon = (float)Config::get().optimization.local.bfgs.grad;
	optimizer(x);

	if (Config::get().general.verbosity > 4)
	{
		std::cout << "Optimization done (status " << optimizer.state() << "). Evaluations:" << optimizer.iter() << '\n';
	}

	return cPtr->g();
}

double neb::g_new()
{
  double Rp1mag{ 0.0 }, Rm1mag{ 0.0 }, energytemp{ 0.0 };
  Rp1.resize(cPtr->size());
  Rm1.resize(cPtr->size());
  Fvertical.resize(cPtr->size());
  Fpar.resize(cPtr->size());
  if (Config::get().neb.IDPP) Fidpp.resize(cPtr->size());
  grad_tot.clear();
 
  calc_tau();
  


  for (size_t im = 1; im < num_images - 1; im++)
  {
	
    imagi[im].clear();
    for (size_t kk = (im - 1)*cPtr->size(); kk < im * cPtr->size(); kk++)  imagi[im].push_back(images_initial[kk]);

    cPtr->set_xyz(imagi[im]);
    energies_NEB[im] = energytemp = this->cPtr->g();
    energytemp += energytemp;
    if (ClimbingImage == true && num_images == static_cast<decltype(num_images)>(CIMaximum))
    {
      double magni = 0.0;
      magni = dot_uneq(cPtr->g_xyz(), tau[im]);
      if (len(tau[im]) == 0.0) { magni = 0.0; }
      else { magni = magni / len(tau[im]); }
      for (size_t i = 0; i < cPtr->size(); i++) cPtr->update_g_xyz(i, cPtr->g_xyz(i) - tau[im][i] * magni*2.0);
    }
    else
    {

      tauderiv = dot_uneq(cPtr->g_xyz(), tau[im]);

	  if (len(tau[im]) == 0.0) { tauderiv = 0.0; }
      else { tauderiv /= len(tau[im]); }
      if (tauderiv != tauderiv) tauderiv = 0.0;
      for (size_t j = 0; j < cPtr->size(); j++)
      {

        Rm1[j].x() = imagi[im - 1][j].x() - imagi[im][j].x();
        Rm1[j].y() = imagi[im - 1][j].y() - imagi[im][j].y();
        Rm1[j].z() = imagi[im - 1][j].z() - imagi[im][j].z();

        Rp1[j].x() = imagi[im + 1][j].x() - imagi[im][j].x();
        Rp1[j].y() = imagi[im + 1][j].y() - imagi[im][j].y();
        Rp1[j].z() = imagi[im + 1][j].z() - imagi[im][j].z();
		
      }
      Rm1mag = len(imagi[im - 1]) - len(imagi[im]);
      Rp1mag = len(imagi[im + 1]) - len(imagi[im]);
      if (Rm1mag != Rm1mag) Rm1mag = 0.0;
      if (Rp1mag != Rp1mag) Rp1mag = 0.0;
      for (size_t i = 0; i < cPtr->size(); i++)
      {
        Fvertical[i].x() = cPtr->g_xyz(i).x() - tauderiv * tau[im][i].x();
        Fvertical[i].y() = cPtr->g_xyz(i).y() - tauderiv * tau[im][i].y();
        Fvertical[i].z() = cPtr->g_xyz(i).z() - tauderiv * tau[im][i].z();

        Fpar[i].x() = springconstant * (Rp1mag - Rm1mag) * tau[im][i].x();
        Fpar[i].y() = springconstant * (Rp1mag - Rm1mag) * tau[im][i].y();
        Fpar[i].z() = springconstant * (Rp1mag - Rm1mag) * tau[im][i].z();

      }
	  if (Config::get().neb.IDPP)
	  {
		  auto bond_dummy{ bond_dir(imagi[im]) };
		  auto eu_dummy{ euclid_dist<double>(imagi[im]) };
		  Fidpp = idpp_gradients<double, dist_T<double>>(bond_dummy, eu_st, eu_dummy, eu_fi, im);
	  }
    }

    for (size_t j = 0; j < cPtr->size(); j++)
    {
      if (Config::get().neb.IDPP)
      {
        auto const g = Fvertical[j] + Fpar[j] + Fidpp[j];
        cPtr->update_g_xyz(j, g);
        grad_tot.push_back(g);
      }
      else
      {
        auto const g = Fvertical[j] + Fpar[j];
        cPtr->update_g_xyz(j, g);
        grad_tot.push_back(g);
		
      }
	 
    }
	
  }
  
  return energytemp;
}

double neb::g_new_maxflux()
{
	double Rp1mag{ 0.0 }, Rm1mag{ 0.0 }, energytemp{ 0.0 }, cosi2{ 0.0 }, kappa{ 0.0 };
	coords::Cartesian_Point rv_p, rv_n;
	Rp1.resize(cPtr->size());
	Rm1.resize(cPtr->size());
	Fvertical.resize(cPtr->size());
	Fpar.resize(cPtr->size());
	grad_tot.clear();
	calc_tau();
	for (size_t im = 1; im < num_images - 1; im++)
	{
		imagi[im].clear();
		for (size_t kk = (im - 1)*cPtr->size(); kk < im * cPtr->size(); kk++)  imagi[im].push_back(images_initial[kk]);

		cPtr->set_xyz(imagi[im]);
		energies_NEB[im] = energytemp = this->cPtr->g();
		energytemp += energytemp;
		if (ClimbingImage == true && num_images == static_cast<decltype(num_images)>(CIMaximum))
		{
			double magni = 0.0;
			magni = dot_uneq(cPtr->g_xyz(), tau[im]);

			if (len(tau[im]) == 0.0) { magni = 0.0; }
			else { magni = magni / len(tau[im]); }
			for (size_t i = 0; i < cPtr->size(); i++) cPtr->update_g_xyz(i, cPtr->g_xyz(i) - tau[im][i] * magni*2.0);
		}
		else
		{
			tauderiv = dot_uneq(cPtr->g_xyz(), tau[im]);
			if (len(tau[im]) == 0.0) { tauderiv = 0.0; }
			else { tauderiv /= len(tau[im]); }
			if (tauderiv != tauderiv) tauderiv = 0.0;
			for (size_t j = 0; j < cPtr->size(); j++)
			{

				Rm1[j].x() = imagi[im - 1][j].x() - imagi[im][j].x();
				Rm1[j].y() = imagi[im - 1][j].y() - imagi[im][j].y();
				Rm1[j].z() = imagi[im - 1][j].z() - imagi[im][j].z();

				Rp1[j].x() = imagi[im + 1][j].x() - imagi[im][j].x();
				Rp1[j].y() = imagi[im + 1][j].y() - imagi[im][j].y();
				Rp1[j].z() = imagi[im + 1][j].z() - imagi[im][j].z();
			}
			Rm1mag = len(imagi[im - 1]) - len(imagi[im]);
			Rp1mag = len(imagi[im + 1]) - len(imagi[im]);
			if (Rm1mag != Rm1mag) Rm1mag = 0.0;
			if (Rp1mag != Rp1mag) Rp1mag = 0.0;
			for (size_t i = 0; i < cPtr->size(); i++)
			{
				Fvertical[i].x() = cPtr->g_xyz(i).x() - tauderiv * tau[im][i].x();
				Fvertical[i].y() = cPtr->g_xyz(i).y() - tauderiv * tau[im][i].y();
				Fvertical[i].z() = cPtr->g_xyz(i).z() - tauderiv * tau[im][i].z();

				Fpar[i].x() = springconstant * (Rp1mag - Rm1mag) * tau[im][i].x();
				Fpar[i].y() = springconstant * (Rp1mag - Rm1mag) * tau[im][i].y();
				Fpar[i].z() = springconstant * (Rp1mag - Rm1mag) * tau[im][i].z();
				
			}
		
		}

		
		for (size_t i = 0; i < cPtr->size(); i++)
		{
			auto L = scon::geometric_length(tau[im][i]);
			if (L != 0.0)
			{
				cosi2 = (cPtr->xyz(i).x()*tau[im][i].x() + cPtr->xyz(i).y()*tau[im][i].y() + cPtr->xyz(i).z()*tau[im][i].z()) / (L * L);
			}
			else
			{
				cosi2 = 0.0;
			}
			if (cosi2 != cosi2)cosi2 = 0.0;
			rv_n.x() = cPtr->xyz(i).x() - (cosi2 * tau[im][i].x());
			rv_n.y() = cPtr->xyz(i).y() - (cosi2 * tau[im][i].y());
			rv_n.z() = cPtr->xyz(i).z() - (cosi2 * tau[im][i].z());

			kappa = acos(dot(tau[im - 1][i], tau[im + 1][i])) / (len(imagi[im][i] - imagi[im - 1][i]) + len(imagi[im + 1][i] - imagi[im][i]));
			if (kappa != kappa)kappa = 0.0;
			rv_p.x() =  (kappa / _KT_ )*rv_n.x();
			rv_p.y() =  (kappa / _KT_ )*rv_n.y();
			rv_p.z() =  (kappa / _KT_ )*rv_n.z();

			auto const g = Fpar[i] + Fvertical[i] - rv_p;
			cPtr->update_g_xyz(i, g);
			grad_tot.push_back(g);
		}
		
	}

	return energytemp;
}

void neb::calc_shift(void)
{
  std::ptrdiff_t laf{ 0 };
  double diff{ 0.0 }, gridp{ 0.0 },
	  x{ 0.0 }, y{ 0.0 }, z{ 0.0 },
	  distx{ 0.0 }, disty{ 0.0 }, distz{ 0.0 };
  std::vector<double> posx(num_images), posy(num_images), posz(num_images), gridx(num_images), gridy(num_images),
	  gridz(num_images), shiftx(cPtr->size()), shifty(cPtr->size()), shiftz(cPtr->size());
  std::vector <std::vector < std::vector < scon::c3<double> > > > position{ N };
  scon::c3 <double> pos3(0.0);
  std::vector <coords::Cartesian_Point> interpol_position{ N };
  coords::Representation_3D temp_int;

  for (std::size_t i = 0; i < this->cPtr->size(); i++) {

    for (std::size_t j = 0; j < (num_images); j++) {

      diff = (double)j / num_images;

      image_ini[j][i];
      posx[j] = imagi[j][i].x();
      posy[j] = imagi[j][i].y();
      posz[j] = imagi[j][i].z();
      gridx[j] = double(j);
      gridy[j] = double(j);
      gridz[j] = double(j);
      //std::cout << gridx[j] << "   " << posx[j] << endl;
    }
    //std::cout << '\n';

   Lagrange_interp splinex(gridx, posx);
   Lagrange_interp spliney(gridy, posy);
   Lagrange_interp splinez(gridz, posz);
    position[i].resize(num_images);

    for (size_t k = 0; k < num_images; k++)

    {


      //int jj { 0 }, it{ 0 };
      //double  tnm { 0.0 }, sumx { 0.0 }, sumy{ 0.0 }, sumz{ 0.0 }, del{ 0.0 };
      //double skx0 { abs((splinex.interp(k + 0.00001) - splinex.interp(k)) / 0.00001) };
      //double skx { abs((splinex.interp((k + 1) + 0.00001) - splinex.interp(k + 1)) / 0.00001) };
      //double sky0 { abs((spliney.interp(k + 0.00001) - spliney.interp(k)) / 0.00001) };
      //double sky { abs((spliney.interp((k + 1) + 0.00001) - spliney.interp(k + 1)) / 0.00001) };
      //double skz0 { abs((splinez.interp(k + 0.00001) - splinez.interp(k)) / 0.00001) };
      //double skz { abs((splinez.interp((k + 1) + 0.00001) - splinez.interp(k + 1)) / 0.00001) };
          //double s { 0.5*((k + 1) - k)*(splinex.interp(k + 1) + splinex.interp(k)) };
      //double sx { 0.5*((k + 1) - k)*(skx + skx0) };
      //double sy { 0.5*((k + 1) - k)*(sky + sky0) };
      //double sz { 0.5*((k + 1) - k)*(skz + skz0) };

      //   pathlenx = 0.0;
      //   pathleny = 0.0;
      //   pathlenz = 0.0;
      //for (int m = 1; m < 8; m++)
      //{
      //	for (it = 1, jj = 1; jj < m - 1; jj++) it <<= 1;

       //	tnm = it;
       //	if (tnm == 0)tnm = 1;
       //	del = ((k + 1) - k) / tnm;
       //	x = k + 0.5*del;
       //  /* std::cout << "tnm  " << x<< '\n';*/
       //  for (sumx = 0.0, sumy = 0.0, sumz = 0.0, jj = 0; jj < it; jj++, x += del)
       //  {
       //    skx = abs((splinex.interp(x + 0.00001) - splinex.interp(x)) / 0.00001);
       //    sky = abs((spliney.interp(x + 0.00001) - spliney.interp(x)) / 0.00001);
       //    skz = abs((splinez.interp(x + 0.00001) - splinez.interp(x)) / 0.00001);

       //    sumx += skx;
       //    sumy += sky;
       //    sumz += skz;

       //  }
       //  sx = 0.5*(sx + ((k + 1) - k)*sumx / tnm);
       //  sy = 0.5*(sy + ((k + 1) - k)*sumy / tnm);
       //  sz = 0.5*(sz + ((k + 1) - k)*sumz / tnm);
       //}
       //pathlenx += sx;
       //pathleny += sy;
       //pathlenz += sz;
      laf = 0;
      gridp = double(k);
      while (gridp <= (k + 1))
      {

        laf++;
        gridp += Config::get().neb.INT_IT;
        x = splinex.interpolate(gridp);
        y = spliney.interpolate(gridp);
        z = splinez.interpolate(gridp);
        distx += x;
        disty += y;
        distz += z;
        pos3.x() = x;
        pos3.y() = y;
        pos3.z() = z;

        position[i][k].push_back(pos3);
      }
      distx /= laf;
      disty /= laf;
      distz /= laf;
    }
    //shiftx[i] = pathlenx / (num_images - 2);
    //shifty[i] = pathleny / (num_images - 2);
    //shiftz[i] = pathlenz / (num_images - 2);

  }

  std::ofstream out("INTERPOL_preopt.arc", std::ios::app), out2("INTERPOL_opt.arc", std::ios::app),
    out3("ENERGY_INT_PREOPT.dat", std::ios::app), out4("ENERGY_INT_OPT.dat", std::ios::app);
  interpol_position.resize(N);
  tau_int.resize(N);
  interpol_position.resize(N);
  temp_int.resize(N);

  for (std::ptrdiff_t i = 1; i < std::ptrdiff_t(num_images - 1); i++)
  {

    for (std::ptrdiff_t k = 0; k < laf; k++) {
      out << "     " << N << " IMAGE: " << k << "  global counter:  " << i << '\n';
      temp_int.clear();
      for (size_t j = 0; j < N; j++) {




        float abso(0.0);
        /*ofstream off ("energies");*/



        tau_int[j].x() = float(position[j][i][k].x() - position[j][i - 1][k].x() / abs(position[j][i][k].x() - position[j][i - 1][k].x()) 
			+ (position[j][i + 1][k].x() - position[j][i][k].x()) / abs(position[j][i + 1][k].x() - position[j][i][k].x()));
        tau_int[j].y() = float(position[j][i][k].y() - position[j][i - 1][k].y() / abs(position[j][i][k].y() - position[j][i - 1][k].y())
			+ (position[j][i + 1][k].y() - position[j][i][k].y()) / abs(position[j][i + 1][k].y() - position[j][i][k].y()));
        tau_int[j].z() = float(position[j][i][k].z() - position[j][i - 1][k].z() / abs(position[j][i][k].z() - position[j][i - 1][k].z())
			+ (position[j][i + 1][k].z() - position[j][i][k].z()) / abs(position[j][i + 1][k].z() - position[j][i][k].z()));
        abso = 0.0;
        if (tau_int[j].x() - tau_int[j].x() != 0) tau_int[j].x() = 0.0;
        if (tau_int[j].y() - tau_int[j].y() != 0) tau_int[j].y() = 0.0;
        if (tau_int[j].z() - tau_int[j].z() != 0) tau_int[j].z() = 0.0;
        abso = len(tau_int[j]);
        if (abso != 0.0) normalize(tau_int[j]);
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
        out << '\n';


      }
      cPtr->set_xyz(temp_int);
      out3 << cPtr->g() << '\n';
      out4 << lbfgs_int(tau_int) << '\n';

      out2 << "     " << N << '\n';

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
        out2 << '\n';


      }


    }
  }
  out2.close();
  out.close();
}



double neb::lbfgs_int(std::vector <scon::c3 <float> >)
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
    std::cout << "Optimization done (status " << optimizer.state() << "). Evaluations:" << optimizer.iter() << '\n';
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

double neb::dot_uneq (coords::Representation_3D const &a, coords::Representation_3D const &b)
{
	double temp(0.0);
	for (size_t i=0; i < a.size(); i++)
	{
		temp += a[i].x()*b[i].x();
		temp += a[i].y()*b[i].y();
		temp += a[i].z()*b[i].z();
	}
	return temp;
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

void neb::opt_internals(std::ptrdiff_t &count, const std::vector<std::vector<size_t> >& atoms_remember) {

  coords::Ensemble_PES *ep = new coords::Ensemble_PES;
  std::ostringstream name;



  //Just to get xyz-data beginning by the origin
  //
  for (size_t i = 0U; i < num_images; i++) {
    cPtr->set_xyz(imagi[i]);
    cPtr->to_internal();
    cPtr->to_xyz();

    imagi[i] = cPtr->xyz();
  }
  //Get output for comparison
  name << "Formatted.xyz";
  print(name.str(), imagi, count);
  name.str("");
  //std::cout <<"SIZE_ALL "<< atoms_remember[1].size() << '\n';
  //Loop to optimize every image
  for (size_t i = 1U; i < num_images - 1; i++) {

    //reverse every fixation
    defix_all();

    cPtr->set_xyz(imagi[i]);

    //fill ep with current values
    ep->clear();
    ep->push_back(imagi[i]);

    //fix atoms for this and last image to guarantee no change in dihedrals during the change
    //std::cout << "SIZE_B " << atoms_remember[i].size() << '\n';
    execute_fix(atoms_remember[i]);
    execute_fix(atoms_remember[i - 1]);

    //Just to show which atoms are fixed
    for (size_t j = 0U; j < N; j++) {
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

void neb::internal_execute(std::vector <coords::Representation_3D> &input, std::vector<std::vector<size_t> >& atoms_remember) {

  size_t i = 0U, j = 0U, imgs = num_images, atom_iter = 0U;

  std::vector<std::string> Atom_name_input, distance_string_swap, angle_string_swap, dihedral_string_swap;
  std::vector<std::vector<std::string> > DAD_names;
  std::vector<std::vector<std::vector<std::string> > > name;

  std::vector<double> x_input, y_input, z_input;
  std::vector<std::vector<double> > imgs_x_input, imgs_y_input, imgs_z_input, dist, anglevec, dihedralvec;
  std::vector<std::vector<std::vector<double> > > dist_all, anglevec_all, dihedralvec_all, compare;

  std::vector<size_t> int_swap;
  std::vector<std::vector<size_t> > bonds, which_bonds;
  std::vector<std::vector<std::vector<size_t> > > which_img;
  std::vector<std::vector<std::vector<std::vector<size_t> > > > involved_bonds;

  std::vector<std::vector<size_t> > i_remember, j_remember, i_remember_all, j_remember_all;

  std::ofstream myfile;

  std::ostringstream string_swap;

  //clearing variables used to carry results outside the function
  atoms_remember.clear();

  //writing input in own variables
  for (i = 0U; i < N; i++) {
    int_swap.clear();

    Atom_name_input.push_back(cPtr->atoms(i).symbol().substr(0U, 2U));
    x_input.push_back(input[0][i].x());
    y_input.push_back(input[0][i].y());
    z_input.push_back(input[0][i].z());

    int_swap.push_back(i + 1);
    for (j = 0U; j < cPtr->atoms(i).bonds().size(); j++)
    {
      int_swap.push_back(cPtr->atoms(i).bonds()[j] + 1U);
    }
    bonds.push_back(int_swap);

  }

  imgs_x_input.push_back(x_input);
  imgs_y_input.push_back(y_input);
  imgs_z_input.push_back(z_input);

  //because variable bonds won't change during images just the x,y and z values are needed from now on
  for (i = 1U; i < imgs; i++) {
    x_input.clear();
    y_input.clear();
    z_input.clear();

    for (j = 0U; j < N; j++) {
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
  for (i = 0U; i < N; i++) {
    if (bonds[i].size() > 2) {
      which_bonds.push_back(bonds[i]); atom_iter++;
    }
  }

  //std::cout << atom_iter << std::endl;

  //now the first values can be obtained and it has been taken care that there is no need to
  //use such a expensive function again
  get_values(imgs_x_input[0], imgs_y_input[0], imgs_z_input[0], atom_iter, which_bonds, &dist, &anglevec, &dihedralvec, &involved_bonds, bonds);

  //The name of the atoms are obtained. Just for output
  //
  for (i = 0U; i < involved_bonds[DIST].size(); i++) {
    for (j = 0U; j < involved_bonds[DIST][i].size(); j++) {
      double_or_not(involved_bonds[DIST], involved_bonds[DIST][i][j], i, j, dist);
    }
    for (j = 0U; j < involved_bonds[ANGLE][i].size(); j++) {
      double_or_not(involved_bonds[ANGLE], involved_bonds[ANGLE][i][j], i, j, anglevec);
    }
    for (j = 0U; j < involved_bonds[DIHEDRAL][i].size(); j++) {
      double_or_not(involved_bonds[DIHEDRAL], involved_bonds[DIHEDRAL][i][j], i, j, dihedralvec);
    }
  }

  for (i = 0U; i < involved_bonds[DIST].size(); i++) {

    distance_string_swap.clear();
    angle_string_swap.clear();
    dihedral_string_swap.clear();

    DAD_names.clear();
    for (j = 0U; j < involved_bonds[DIST][i].size(); j++) {
      string_swap << involved_bonds[DIST][i][j][0] << Atom_name_input[involved_bonds[DIST][i][j][0] - 1] << " ";
      string_swap << involved_bonds[DIST][i][j][1] << Atom_name_input[involved_bonds[DIST][i][j][1] - 1];
      distance_string_swap.push_back(string_swap.str()); string_swap.str("");
    }
    for (j = 0U; j < involved_bonds[ANGLE][i].size(); j++) {
      string_swap << involved_bonds[ANGLE][i][j][0] << Atom_name_input[involved_bonds[ANGLE][i][j][0] - 1] << " ";
      string_swap << involved_bonds[ANGLE][i][j][1] << Atom_name_input[involved_bonds[ANGLE][i][j][1] - 1] << " ";
      string_swap << involved_bonds[ANGLE][i][j][2] << Atom_name_input[involved_bonds[ANGLE][i][j][2] - 1];
      angle_string_swap.push_back(string_swap.str()); string_swap.str("");
    }
    for (j = 0U; j < involved_bonds[DIHEDRAL][i].size(); j++) {
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
  for (i = 1U; i < imgs; i++) {

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
    for (size_t lk(i_remember_all[li].size() - Config::get().neb.NUMBER_OF_DIHEDRALS); lk < i_remember_all[li].size(); lk++)
    {

      //std::cout << i_remember_all[li][lk] << "   " << j_remember_all[li][lk] << '\n';

      atoms_remember[li].push_back(involved_bonds[DIHEDRAL][i_remember_all[li][lk]][j_remember_all[li][lk]][0]);
      atoms_remember[li].push_back(involved_bonds[DIHEDRAL][i_remember_all[li][lk]][j_remember_all[li][lk]][1]);
      atoms_remember[li].push_back(involved_bonds[DIHEDRAL][i_remember_all[li][lk]][j_remember_all[li][lk]][2]);
      atoms_remember[li].push_back(involved_bonds[DIHEDRAL][i_remember_all[li][lk]][j_remember_all[li][lk]][3]);
      atoms_remember.push_back(involved_bonds[DIST][i_remember[DIST][li]][j_remember[DIST][li]]);
      //atoms_remember[li].push_back(involved_bonds[DIHEDRAL][i_remember_all[li][lk]][j_remember_all[li][lk]]);
    }
  }

  for (size_t lo(0U); lo < num_images - 1; lo++)
  {
    for (size_t lp(0U); lp < atoms_remember[lo].size(); lp++)
    {
      //std::cout << "ATOMS_REM " <<"lo "<<lo<<"  lp "<<lp <<"  "<< atoms_remember[lo][lp] << '\n';
    }
  }
  //for (i = 0U; i<i_remember[DIST].size(); i++){
  //	std::cout << "SIZE#" << i_remember[DIST].size()<<'\n';
  //	std::cout << i_remember[DIST][i] << " " << j_remember[DIST][i] << std::endl;
  //	atoms_remember.push_back(involved_bonds[DIST][i_remember[DIST][i]][j_remember[DIST][i]]);
  //	//atoms_remember.push_back(involved_bonds[DIHEDRAL][i_remember[DIHEDRAL][i]][j_remember[DIHEDRAL][i]]);
  //}

  //write output
  //
  myfile.open("Internals.txt");


  for (size_t k = 0U; k < num_images; k++) {


    for (i = 0U; i < dist_all[k].size(); i++) {
      if (dist_all[k][i].size() > anglevec_all[k][i].size() && dist_all[k][i].size() > dihedralvec_all[k][i].size())for (j = 0U; j < dist_all[k][i].size(); j++) {
        if (dist_all[k][i].size() > j) {
          myfile << std::right << std::setw(10) << name[i][DIST][j] << " ";
          myfile << std::left << std::setw(10) << dist_all[k][i][j];
        }
        else myfile << std::right << std::setw(21) << " ";
        if (anglevec_all[k][i].size() > j) {
          myfile << std::right << std::setw(10) << name[i][ANGLE][j] << " ";
          myfile << std::left << std::setw(10) << anglevec_all[k][i][j];
        }
        else myfile << std::right << std::setw(21) << " ";
        if (dihedralvec_all[k][i].size() > j) {
          myfile << std::right << std::setw(10) << name[i][DIHEDRAL][j] << " ";
          myfile << std::left << std::setw(10) << dihedralvec_all[k][i][j];
        }
        else myfile << std::right << std::setw(21) << " ";
        myfile << std::endl;
      }
      else if (anglevec_all[k][i].size() > dist_all[k][i].size() && anglevec_all[k][i].size() > dihedralvec_all[k][i].size())for (j = 0U; j < anglevec_all[k][i].size(); j++) {
        if (dist_all[k][i].size() > j) {

          myfile << std::right << std::setw(10) << name[i][DIST][j] << " ";
          myfile << std::left << std::setw(10) << dist_all[k][i][j];
        }
        else myfile << std::right << std::setw(21) << " ";
        if (anglevec_all[k][i].size() > j) {
          myfile << std::right << std::setw(10) << name[i][ANGLE][j] << " ";
          myfile << std::left << std::setw(10) << anglevec_all[k][i][j];
        }
        else myfile << std::right << std::setw(21) << " ";
        if (dihedralvec_all[k][i].size() > j) {
          myfile << std::right << std::setw(10) << name[i][DIHEDRAL][j] << " ";
          myfile << std::left << std::setw(10) << dihedralvec_all[k][i][j];
        }
        else myfile << std::right << std::setw(21) << " ";
        myfile << std::endl;
      }
      else for (j = 0U; j < dihedralvec[i].size(); j++) {
        if (dist[i].size() > j) {
          myfile << std::right << std::setw(10) << name[i][DIST][j] << " ";
          myfile << std::left << std::setw(10) << dist_all[k][i][j];
        }
        else myfile << std::right << std::setw(21) << " ";
        if (anglevec[i].size() > j) {
          myfile << std::right << std::setw(10) << name[i][ANGLE][j] << " ";
          myfile << std::left << std::setw(10) << anglevec_all[k][i][j];
        }
        else myfile << std::right << std::setw(21) << " ";
        if (dihedralvec[i].size() > j) {
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

  for (i = 0U; i < atoms_remember.size(); i++) {
    myfile << "Image " << i + 1 << ": ";
    for (j = 0U; j < atoms_remember[i].size(); j++) {
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

void neb::no_torsion(const size_t& image_remember, const std::vector<int>& atoms_remember) {
  for (size_t i = 0U; i < atoms_remember.size(); i++) {
    imagi[image_remember][atoms_remember[i] - 1].x() = imagi[image_remember - 1][atoms_remember[i] - 1].x();
    imagi[image_remember][atoms_remember[i] - 1].y() = imagi[image_remember - 1][atoms_remember[i] - 1].y();
    imagi[image_remember][atoms_remember[i] - 1].z() = imagi[image_remember - 1][atoms_remember[i] - 1].z();
  }
}

//void execute_fix
//Function to fixate atoms
//

void neb::execute_fix(const std::vector<size_t>& atoms_remember) {
  //std::cout <<"SIZE_exe "<< atoms_remember.size() << endl;
  for (size_t i = 0U; i < atoms_remember.size(); i++) {
    cPtr->set_fix(atoms_remember[i] - 1, true);
    //std::cout <<"REM "<< atoms_remember[i] - 1 << '\n';
    //scon::insert_unique_sorted(Config::set().coords.fixed, atoms_remember[i] - 1);
  }
}//end execute_fix

void neb::execute_defix(const std::vector<size_t>& atoms_remember) {
  for (size_t i = 0U; i < atoms_remember.size(); i++) {
    cPtr->set_fix(atoms_remember[i] - 1, false);

  }

}

//void defix_all
//Function to loosen all atoms
void neb::defix_all() {
  for (size_t i = 0U; i < N; i++) {
    cPtr->set_fix(i, false);
  }
}//end defix_all

//void biggest
//Function to get the biggest change during two images
//
void neb::biggest(const std::vector<std::vector<std::vector<double> > >& dist, const std::vector<std::vector<std::vector<double> > >& anglevec,
  const std::vector<std::vector<std::vector<double> > >& dihedralvec, std::vector<std::vector<std::vector<double> > >&, std::vector<std::vector<std::vector<size_t> > >&, std::vector<std::vector<size_t> >& i_remember,
  std::vector<std::vector<size_t> >& j_remember, std::vector<std::vector<size_t> > & i_remember_all, std::vector<std::vector<size_t> > & j_remember_all) {

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
  for (i = 0U; i < num_images - 1; i++) {
    comp_value_swap.clear();
    i_temp.clear();
    j_temp.clear();

    //Fill the compare array with the first values to initialize the comparison
    for (j = 0U; j < dist[i].size(); j++) {
      for (k = 0U; k < dist[i][j].size(); k++) {
        comp_value_swap.push_back(fabs(dist[i][j][k] - dist[i + 1][j][k]));
        i_temp.push_back(j);
        j_temp.push_back(k);
      }
    }
    //Compare the change of every value to get the biggest
    double_swap = comp_value_swap[0];
    for (j = 1U; j < comp_value_swap.size(); j++) {
      if (double_swap < comp_value_swap[j]) {
        double_swap = comp_value_swap[j];
        i_remembered[i] = i_temp[j];
        j_remembered[i] = j_temp[j];
      }
    }
  }
  i_remember[DIST] = i_remembered;
  j_remember[DIST] = j_remembered;

  for (i = 0U; i < num_images - 1; i++) {
    comp_value_swap.clear();
    i_temp.clear();
    j_temp.clear();

    for (j = 0U; j < anglevec[i].size(); j++) {
      for (k = 0U; k < anglevec[i][j].size(); k++) {
        comp_value_swap.push_back(fabs(anglevec[i][j][k] - anglevec[i + 1][j][k]));
        i_temp.push_back(j);
        j_temp.push_back(k);
      }
    }

    double_swap = comp_value_swap[0];
    for (j = 1U; j < comp_value_swap.size(); j++) {
      if (double_swap < comp_value_swap[j]) {
        double_swap = comp_value_swap[j];
        i_remembered[i] = i_temp[j];
        j_remembered[i] = j_temp[j];
      }
    }
  }
  i_remember[ANGLE] = i_remembered;
  j_remember[ANGLE] = j_remembered;

  myfile.open("Logfile.txt");

  for (i = 0U; i < num_images - 1; i++) {
    comp_value_swap.clear();
    i_temp.clear();
    j_temp.clear();

    for (j = 0U; j < dihedralvec[i].size(); j++) {
      for (k = 0U; k < dihedralvec[i][j].size(); k++) {
        if ((dihedralvec[i][j][k] < 0.0&&dihedralvec[i + 1][j][k] < 0.0) || (dihedralvec[i][j][k] > 0.0&&dihedralvec[i + 1][j][k] > 0.0)) {
          comp_value_swap.push_back(fabs(fabs(dihedralvec[i][j][k]) - fabs(dihedralvec[i + 1][j][k])));
          i_temp.push_back(j);
          j_temp.push_back(k);
        }
        else {
          if (fabs(dihedralvec[i][j][k]) + fabs(dihedralvec[i + 1][j][k]) > 180.0) {
            comp_value_swap.push_back(360.0 - fabs(dihedralvec[i][j][k]) - fabs(dihedralvec[i + 1][j][k]));
            i_temp.push_back(j);
            j_temp.push_back(k);
          }
          else {
            comp_value_swap.push_back(fabs(fabs(dihedralvec[i][j][k]) + fabs(dihedralvec[i + 1][j][k])));
            i_temp.push_back(j);
            j_temp.push_back(k);
          }
        }
      }
    }

    for (j = 0U; j < comp_value_swap.size(); j++) {
      if ((int)comp_value_swap[j] > 0) myfile << "Image " << i + 1 << " to " << i + 2 << ": " << (int)comp_value_swap[j] << std::endl;
    }

    double_swap = comp_value_swap[0];
    bool swapped(true);

    while (swapped == true) {

      swapped = false;
      for (j = 0U; j < comp_value_swap.size(); j++) {

        if (comp_value_swap[j] > comp_value_swap[j + 1]) {
          swapped = true;
          size_t temp_i, temp_j;
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
          /*std::cout << " IMG " << i << "DOUBLE SWAP " << double_swap << " Index i " << i_temp[j] << " Index j " << j_temp[j] << '\n';
          i_remember_two[i].push_back(i_temp[j]);
          j_remember_two[i].push_back(j_temp[j]);
          i_remembered[i] = i_temp[j];
          j_remembered[i] = j_temp[j];*/
        }
        //nn = newnn;
      }


    }

    for (j = 0U; j < comp_value_swap.size(); j++) {

      //std::cout << "COMP_swap " << comp_value_swap[j] << " Index i " << i_temp[j] << " Index j " << j_temp[j] << '\n';
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
void neb::double_or_not(std::vector<std::vector<std::vector<size_t> > >& involved_bonds, const std::vector<size_t>& pivot, const size_t& a, const size_t& b, std::vector<std::vector<double> >& vec) {

  bool all_the_same;

  size_t i, j, k, l;
  std::vector<bool> its_double;

  its_double.resize(pivot.size());

  //check every whatever if it is double by compareing the atoms involved. Afterwards killing all double values
  for (i = 0U; i < involved_bonds.size(); i++) {
    for (j = 0U; j < involved_bonds[i].size(); j++) {
      if ((j < (involved_bonds[i].size() - 1)) && (i == a && j == b)) j++;
      else if ((j == (involved_bonds[i].size() - 1)) && (i == a && j == b)) break;
      for (k = 0U; k < pivot.size(); k++) {
        its_double[k] = false;
      }
      for (k = 0U; k < involved_bonds[i][j].size(); k++) {
        for (l = 0U; l < pivot.size(); l++) {
          if (involved_bonds[i][j][k] == pivot[l]) {
            its_double[k] = true;
            break;
          }
        }
      }
      all_the_same = true;
      if (pivot.size() > 0) for (k = 0U; k < pivot.size() - 1; k++) {
        if (its_double[k] != its_double[k + 1]) all_the_same = false;
      }
      if (all_the_same&&its_double[0]) {
        vec[i].erase(vec[i].begin() + j);
        involved_bonds[i].erase(involved_bonds[i].begin() + j);
      }
    }
  }
}

bool neb::dihedral_or_not(const int& first_atom, const std::vector<int>& dist, const std::vector<int>& ang, const int& last_atom) {

  if ((first_atom == dist[1] || first_atom == dist[2] || first_atom == dist[3] || first_atom == dist[4]) && (last_atom == ang[1] || last_atom == ang[2] || last_atom == ang[3] || last_atom == ang[4])) return true;
  else if ((first_atom == ang[1] || first_atom == ang[2] || first_atom == ang[3] || first_atom == ang[4]) && (last_atom == dist[1] || last_atom == dist[2] || last_atom == dist[3] || last_atom == dist[4])) return true;
  else return false;

}

void neb::get_values(std::vector<double>& x_val, std::vector<double>& y_val, std::vector<double>& z_val, const size_t& atom_iter, std::vector<std::vector<size_t> >& which_bonds, std::vector<std::vector<double> >* dist, std::vector<std::vector<double> >* anglevec, std::vector<std::vector<double> >* dihedralvec, std::vector<std::vector<std::vector<std::vector<size_t> > > >* involved_bonds, std::vector<std::vector<size_t> >& bonds) {

  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  size_t i, j, k, l;
  size_t a, b, c, d;

  std::vector<size_t> int_swap;
  std::vector<std::vector<size_t> > dist_int_swap, angle_int_swap, dihedral_int_swap;
  std::vector<std::vector<std::vector<size_t> > > dist_int, angle_int, dihedral_int;

  std::vector<double> dist_swap, angle_swap, dihedral_swap;

  //search down the backbone
  for (i = 0U; i < atom_iter; i++) {

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

    for (j = 0U; j < which_bonds[i].size() - 1; j++) {

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

      for (k = 0U; k < (which_bonds[i].size() - (j + 2)); k++) {

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
        for (l = 0U; l < which_bonds[i].size() - 1; l++) {

          //of cause it shoudn't be an atom which takes part either way
          if (which_bonds[i][l + 1] != which_bonds[i][(j + 1) + (k + 1)] && which_bonds[i][l + 1] != which_bonds[i][j + 1] && which_bonds[i][l + 1] != which_bonds[i][0]) {

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
        if (bonds[which_bonds[i][j + 1] - 1].size() > 1) {
          for (l = 0U; l < (bonds[which_bonds[i][j + 1] - 1].size() - 1); l++) {

            if ((bonds[which_bonds[i][j + 1] - 1][l + 1] != which_bonds[i][0]) && (bonds[which_bonds[i][j + 1] - 1][l + 1] != which_bonds[i][(j + 1) + (k + 1)])) {

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
        if (bonds[which_bonds[i][(j + 1) + (k + 1)] - 1].size() > 1) {
          for (l = 0U; l < (bonds[which_bonds[i][(j + 1) + (k + 1)] - 1].size() - 1); l++) {

            if (bonds[which_bonds[i][(j + 1) + (k + 1)] - 1][l + 1] != which_bonds[i][(j + 1) + (k + 1)] && bonds[which_bonds[i][(j + 1) + (k + 1)] - 1][l + 1] != which_bonds[i][0] && bonds[which_bonds[i][(j + 1) + (k + 1)] - 1][l + 1] != which_bonds[i][j + 1]) {

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
void neb::get_values(std::vector<double>& x_val, std::vector<double>& y_val, std::vector<double>& z_val, std::vector<std::vector<double> >* dist, std::vector<std::vector<double> >* anglevec, std::vector<std::vector<double> >* dihedralvec, std::vector<std::vector<std::vector<std::vector<size_t> > > >& involved_bonds, std::vector<std::vector<size_t> >& which_bonds) {

  std::vector<double> dist_swap, angle_swap, dihedral_swap;
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  size_t i, j;
  size_t a, b;
  bool same_atom;


  for (i = 0U; i < involved_bonds[DIST].size(); i++) {
    for (j = 0U; j < involved_bonds[DIST][i].size(); j++) {
      x1 = x_val[involved_bonds[DIST][i][j][0] - 1];
      y1 = y_val[involved_bonds[DIST][i][j][0] - 1];
      z1 = z_val[involved_bonds[DIST][i][j][0] - 1];

      x2 = x_val[involved_bonds[DIST][i][j][1] - 1];
      y2 = y_val[involved_bonds[DIST][i][j][1] - 1];
      z2 = z_val[involved_bonds[DIST][i][j][1] - 1];

      dist_swap.push_back(distance(x1, x2, y1, y2, z1, z2));
    }
    for (j = 0U; j < involved_bonds[ANGLE][i].size(); j++) {

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
    for (j = 0U; j < involved_bonds[DIHEDRAL][i].size(); j++) {

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
      for (size_t k = 1U; k < which_bonds[a].size(); k++) {
        if (which_bonds[a][k] == b) same_atom = true;
      }
      if (same_atom) {
        dihedral_swap.push_back(dihedral_same_atom(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4));
      }
      else {
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

void neb::Sort(const size_t& imgs, std::vector<std::string>& name, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y, std::vector<std::vector<double> >& z, std::vector<std::vector<size_t> >& bonds) {

  size_t i = 0U, j = 0U, k = 0U, iter = 0U, lenght, start;
  std::string string_swap;
  double double_swap;
  std::vector<size_t> vec_swap, remember_C, remember;

  for (i = 0U; i < bonds.size(); i++) {

    if (name[i] == "C" && i != iter) {
      string_swap = name[iter];
      name[iter] = name[i];
      name[i] = string_swap;

      for (size_t n = 0U; n < imgs; n++) {

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

  for (i = 0U; i < bonds.size(); i++) {
    for (j = 0U; j < bonds[i].size(); j++) {
      for (k = 0U; k < remember.size(); k++) {
        if (bonds[i][j] == remember[k]) bonds[i][j] = remember_C[k];
        else if (bonds[i][j] == remember_C[k]) bonds[i][j] = remember[k];
      }
    }
  }
  start = iter;
  lenght = bonds.size() - 1U;
  hetero_Sort(imgs, start, lenght, name, x, y, z, bonds);

  for (i = 0U; i < bonds.size(); i++) {
    lenght = bonds[i].size() - 1U;
    start = 1;
    bond_Sort(bonds[i], start, lenght);
  }
}

void neb::bond_Sort(std::vector<size_t>& bonds, const size_t& start, const size_t& length) {

  size_t i = start, j = length;
  size_t int_swap;
  size_t mid = bonds[(start + length) / 2];

  while (i <= j) {
    while (bonds[i] < mid) i++;
    while (bonds[j] > mid) j--;
    if (i <= j) {
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

void neb::bond_Sort(size_t* a, size_t* b, size_t* c, size_t* d) {

  size_t int_swap;

  bool again = false;

  if (*a > *b) {
    int_swap = *a;
    *a = *b;
    *b = int_swap;
    again = true;
  }
  else if (*b > *c) {
    int_swap = *b;
    *b = *c;
    *c = int_swap;
    again = true;
  }
  else if (*c > *d) {
    int_swap = *c;
    *c = *d;
    *d = int_swap;
    again = true;
  }
  if (again) bond_Sort(a, b, c, d);
}

void neb::bond_Sort(size_t* a, size_t* b, size_t* c) {

  size_t int_swap;

  bool again = false;

  if (*a > *b) {
    int_swap = *a;
    *a = *b;
    *b = int_swap;
    again = true;
  }
  else if (*b > *c) {
    int_swap = *b;
    *b = *c;
    *c = int_swap;
    again = true;
  }
  if (again) bond_Sort(a, b, c);
}

void neb::bond_Sort(size_t* a, size_t* b) {

  size_t int_swap;

  if (*a > *b) {
    int_swap = *a;
    *a = *b;
    *b = int_swap;
  }
}

void neb::hetero_Sort(const size_t& imgs, const size_t& start, const size_t& length, std::vector<std::string>& name, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y, std::vector<std::vector<double> >& z, std::vector<std::vector<size_t> >& bonds) {

  size_t a = 0U, b = 0U, c = 0U, i = start, j = length;
  double double_swap;
  std::vector<size_t> vector_swap, remember_i, remember_j;
  std::string string_swap;
  std::string mid = name[(start + length) / 2];

  while (i <= j) {

    if (name[i] < "H") {
      while (name[i] < mid) i++;
      while (name[j] > mid) j--;
    }
    else if (name[i] >= "H") {
      while (name[i] > mid) i++;
      while (name[j] < mid) j--;
    }
    if (i <= j) {
      string_swap = name[i];
      name[i] = name[j];
      name[j] = string_swap;

      for (size_t n = 0U; n < imgs; n++) {

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

  for (a = 0U; a < bonds.size(); a++) {
    for (b = 0U; b < bonds[a].size(); b++) {
      for (c = 0U; c < remember_i.size(); c++) {
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
double inline neb::distance(const double& x1, const double& x2, const double& y1, const double& y2, const double& z1, const double& z2) {

  double x, y, z;

  x = x2 - x1;
  y = y2 - y1;
  z = z2 - z1;

  return sqrt(x*x + y*y + z*z);

}//end distance

//double angle
//calculating the vectors from 1 to 2 and from 1 to three and calculating the angle
//
double inline neb::angle(const double& x1, const double& x2, const double& x3, const double& y1, const double& y2, const double& y3, const double& z1, const double& z2, const double& z3) {

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
double inline neb::dihedral(const double& x1, const double& x2, const double& x3, const double& x4, const double& y1, const double& y2, const double& y3, const double& y4, const double& z1, const double& z2, const double& z3, const double& z4) {

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
  if (abs(re) > 1.) re = (re < 0. ? -1. : 1.);
  re = acos(re);

  if (length_normal < length_reflection) {
    return (re*(180.0f / 3.1415926535f));
  }
  else {
    return -(re*(180.0f / 3.1415926535f));
  }

}//end dihedral

//double dihedral_same_atom
//same as above just for dihedrals where both points are at the start atom
//
double inline neb::dihedral_same_atom(const double& x1, const double& x2, const double& x3, const double& x4, const double& y1, const double& y2, const double& y3, const double& y4, const double& z1, const double& z2, const double& z3, const double& z4) {

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

  if (length_normal < length_reflection) {
    return (re*(180.0f / 3.1415926535f));
  }
  else {
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


void neb::create_internal_interpolation(std::vector <coords::Representation_3D> &input)
{
  size_t i = 0U, j = 0U, imgs = num_images, N_main = 0U;

  std::vector<double> x_input, y_input, z_input;
  std::vector<std::vector<double>> imgs_x_input, imgs_y_input,
    imgs_z_input, dist, anglevec,
    dihedralvec;
  std::vector<std::vector<std::vector<double>>> dist_all, anglevec_all,
    dihedralvec_all;

  std::vector<size_t> int_swap;
  std::vector<std::vector<size_t>> bonds, which_bonds;
  std::vector<std::vector<std::vector<size_t>>> which_img;
  std::vector<std::vector<std::vector<std::vector<size_t>>>> involved_bonds;
  std::vector<std::vector<std::vector<std::pair<std::vector<size_t>, double>>>> Z_matrices;
  std::vector<size_t> backbone_indeces;
  size_t backbone_iterator = 1;

  for (i = 0U; i < N; i++)
  {
    int_swap.clear();

    x_input.push_back(input[0][i].x());
    y_input.push_back(input[0][i].y());
    z_input.push_back(input[0][i].z());

    int_swap.push_back(i + 1);
    for (j = 0U; j < cPtr->atoms(i).bonds().size(); j++)
    {
      int_swap.push_back(cPtr->atoms(i).bonds()[j] + 1U);
    }
    bonds.push_back(int_swap);
  }

  imgs_x_input.push_back(x_input);
  imgs_y_input.push_back(y_input);
  imgs_z_input.push_back(z_input);

  x_input.clear();
  y_input.clear();
  z_input.clear();

  for (i = 0U; i < N; i++)
  {
    x_input.push_back(input[imgs - 1][i].x());
    y_input.push_back(input[imgs - 1][i].y());
    z_input.push_back(input[imgs - 1][i].z());
  }

  imgs_x_input.push_back(x_input);
  imgs_y_input.push_back(y_input);
  imgs_z_input.push_back(z_input);

  //creates backbone enumeration (0 for terminal atoms)
  for (i = 0U; i < N; i++)
  {
    if (bonds[i].size() > 2)
    {
      which_bonds.push_back(bonds[i]);
      ++N_main;
      backbone_indeces.push_back(backbone_iterator++);
    }
    else
    {
      backbone_indeces.push_back(0);
    }
  }
  get_values(imgs_x_input[0], imgs_y_input[0], imgs_z_input[0],
    N_main, which_bonds, &dist, &anglevec, &dihedralvec,
    &involved_bonds, bonds);

  dihedralvec_all.resize(imgs);
  for (size_t i = 0; i < imgs; ++i)
  {
    dist_all.push_back(dist);
    anglevec_all.push_back(anglevec);
  }
  dihedralvec_all[0] = dihedralvec;

  dist.clear();
  anglevec.clear();
  dihedralvec.clear();

  get_values(imgs_x_input[1], imgs_y_input[1], imgs_z_input[1],
    &dist, &anglevec, &dihedralvec, involved_bonds,
    bonds);

  dihedralvec_all[imgs - 1] = dihedralvec;

  dihedralvec.clear();

  std::vector<std::vector<std::vector<std::pair<std::vector<size_t>, double>>>> redundant_dists,
    redundant_angles,
    redundant_dihedrals;
  std::vector<size_t> temporary;

  redundant_dists.resize(imgs);
  redundant_angles.resize(imgs);
  redundant_dihedrals.resize(imgs);

  redundant_dists[0].resize(N_main);
  redundant_angles[0].resize(N_main);
  redundant_dihedrals[0].resize(N_main);

  redundant_dists[imgs - 1].resize(N_main);
  redundant_angles[imgs - 1].resize(N_main);
  redundant_dihedrals[imgs - 1].resize(N_main);


  /*for (auto &a : involved_bonds[DIHEDRAL])
  {
  for (auto &b : a)
  {
  for (auto &c : b)
  {
  std::cout << c << ' ';
  }
  std::cout << '\n';
  }
  }
  system("pause");*/


  for (size_t j = 0; j < N_main; ++j)
  {
    for (size_t k = 0; k < involved_bonds[DIST][j].size(); ++k)
    {
      redundant_dists[0][j].push_back(std::make_pair(involved_bonds[DIST][j][k],
        dist_all[0][j][k]));
      redundant_dists[imgs - 1][j].push_back(std::make_pair(involved_bonds[DIST][j][k],
        dist_all[imgs - 1][j][k]));
    }
    for (size_t k = 0; k < involved_bonds[ANGLE][j].size(); ++k)
    {
      temporary = { involved_bonds[ANGLE][j][k][1],
        involved_bonds[ANGLE][j][k][0],
        involved_bonds[ANGLE][j][k][2] };
      redundant_angles[0][j].push_back(std::make_pair(temporary,
        anglevec_all[0][j][k]));
      redundant_angles[imgs - 1][j].push_back(std::make_pair(temporary,
        anglevec_all[imgs - 1][j][k]));
    }
    for (size_t k = 0; k < involved_bonds[DIHEDRAL][j].size(); ++k)
    {
      if (cPtr->atoms(involved_bonds[DIHEDRAL][j][k][0] - 1).is_bound_to(involved_bonds[DIHEDRAL][j][k][1] - 1)
        && cPtr->atoms(involved_bonds[DIHEDRAL][j][k][0] - 1).is_bound_to(involved_bonds[DIHEDRAL][j][k][2] - 1)
        && cPtr->atoms(involved_bonds[DIHEDRAL][j][k][0] - 1).is_bound_to(involved_bonds[DIHEDRAL][j][k][3] - 1)) continue;
      temporary = { involved_bonds[DIHEDRAL][j][k][2],
        involved_bonds[DIHEDRAL][j][k][0],
        involved_bonds[DIHEDRAL][j][k][1],
        involved_bonds[DIHEDRAL][j][k][3] };
      redundant_dihedrals[0][j].push_back(std::make_pair(temporary,
        dihedralvec_all[0][j][k]));
      redundant_dihedrals[imgs - 1][j].push_back(std::make_pair(temporary,
        dihedralvec_all[imgs - 1][j][k]));
    }
  }


  /*for (auto &a : redundant_dihedrals)
  {
  for (auto &b : a)
  {
  for (auto &c : b)
  {
  for (auto &d : c.first)
  {
  std::cout << d << ' ';
  }
  std::cout << c.second << '\n';
  }
  }
  }
  system("pause");*/


  std::vector<std::vector<scon::sphericals<float_type>>> Z_matrices_s3;
  coords::PES_Point point;
  coords::Representation_Internal Rep;
  Rep.resize(N);
  Z_matrices.resize(imgs);
  Z_matrices_s3.resize(imgs);

  //Z_matrices[0] = redundant_to_Z_backbone(redundant_dists[0],
  //  redundant_angles[0],
  //  redundant_dihedrals[0],
  //  backbone_indeces);
  //Z_matrices[imgs - 1] = redundant_to_Z_backbone(redundant_dists[imgs - 1],
  //  redundant_angles[imgs - 1],
  //  redundant_dihedrals[imgs - 1],
  //  backbone_indeces);



  //for implementation of NeRF algorithm in CAST
  double constexpr pi = 3.1415926535897932384626433832795;
  coords::s3 mega_temp;
  coords::Coordinates coords_initial, coords_fi;
  coords_initial.init_in(cPtr->atoms(), coords::PES_Point(imagi[0]), true);
  coords_fi.init_in(cPtr->atoms(), coords::PES_Point(imagi[imgs - 1]), true);
  std::tuple<coords::Coordinates, std::vector<size_t>> ini = coords_initial.to_internal();
  std::tuple<coords::Coordinates, std::vector<size_t>> final = coords_fi.to_internal();

  coords::Coordinates coords_ini, coords_final;
  std::vector<size_t> new_order_ini, new_order_final;
  std::tie(coords_ini, new_order_ini) = ini;
  std::tie(coords_final, new_order_final) = final;

  coords::output::formats::tinker tinker_ini_writer(coords_ini);
  tinker_ini_writer.to_stream(std::cout);
  coords::output::formats::zmatrix intern_ini_writer(coords_ini);
  intern_ini_writer.to_stream(std::cout);

  std::cout << coords_final << std::endl;
  coords::output::formats::zmatrix intern_final_writer(coords_final);
  intern_final_writer.to_stream(std::cout);



  Z_matrices[0].clear();
  Z_matrices[0].reserve(N);
  Z_matrices[imgs - 1].clear();
  Z_matrices[imgs - 1].reserve(N);
  for (size_t j = 0; j < N; ++j) {
    Z_matrices[0].emplace_back(j);
    Z_matrices[0][j].resize(3);

    Z_matrices[0][j][0].first = std::vector<std::size_t>{ j + 1, coords_ini.atoms_changeable(j).ibond() + 1 };
    Z_matrices[0][j][0].second = coords_ini.intern(j).radius();
    Z_matrices[0][j][1].first = std::vector<std::size_t>{ j + 1, coords_ini.atoms_changeable(j).ibond() + 1, coords_ini.atoms_changeable(j).iangle() + 1 };
    Z_matrices[0][j][1].second = coords_ini.intern(j).inclination().degrees();
    Z_matrices[0][j][2].first = std::vector<std::size_t>{ j + 1, coords_ini.atoms_changeable(j).ibond() + 1, coords_ini.atoms_changeable(j).iangle() + 1 , coords_ini.atoms_changeable(j).idihedral() + 1 };
    Z_matrices[0][j][2].second = coords_ini.intern(j).azimuth().degrees();

    Z_matrices[imgs - 1].emplace_back(j);
    Z_matrices[imgs - 1][j].resize(3);

    Z_matrices[imgs - 1][j][0].first = std::vector<std::size_t>{ j + 1, coords_final.atoms_changeable(j).ibond() + 1 };
    Z_matrices[imgs - 1][j][0].second = coords_final.intern(j).radius();
    Z_matrices[imgs - 1][j][1].first = std::vector<std::size_t>{ j + 1, coords_final.atoms_changeable(j).ibond() + 1, coords_final.atoms_changeable(j).iangle() + 1 };
    Z_matrices[imgs - 1][j][1].second = coords_final.intern(j).inclination().degrees();
    Z_matrices[imgs - 1][j][2].first = std::vector<std::size_t>{ j + 1, coords_final.atoms_changeable(j).ibond() + 1, coords_final.atoms_changeable(j).iangle() + 1 , coords_final.atoms_changeable(j).idihedral() + 1 };
    Z_matrices[imgs - 1][j][2].second = coords_final.intern(j).azimuth().degrees();
  }



  coords::Cartesian_Point x_ref(1, 0, 0);
  coords::Cartesian_Point y_ref(0, 1, 0);
  coords::Cartesian_Point z_ref(0, 0, 1);
  coords::Cartesian_Point atom_A(-Z_matrices[0][2][0].second + Z_matrices[0][1][0].second * std::sin((pi / 180) * (Z_matrices[0][2][1].second - 90)),
    Z_matrices[0][1][0].second * std::cos(Z_matrices[0][2][1].second - 90),
    0);
  coords::Cartesian_Point atom_B(-Z_matrices[0][2][0].second,
    0,
    0);
  coords::Cartesian_Point atom_C(0, 0, 0);

  double const incr(1. / imgs);

  for (size_t i = 1; i < (imgs - 1); ++i)
  {
    Z_matrices[i].resize(N);
    for (size_t j = 0; j < N; ++j)
    {

      auto const start_radius = coords_final.intern(j).radius();
      auto const start_inclination = coords_final.intern(j).inclination().degrees();
      auto const start_azimuth = coords_final.intern(j).azimuth().degrees();

      auto const total_change_radius = start_radius - coords_ini.intern(j).radius();
      auto const total_change_inclination = start_inclination - coords_ini.intern(j).inclination().degrees();
      auto const total_change_azimuth = start_azimuth - coords_ini.intern(j).azimuth().degrees();

      auto const change_radius = total_change_radius * static_cast<double>(i) / static_cast<double>(imgs);
      auto const change_inclination = total_change_inclination * static_cast<double>(i) / static_cast<double>(imgs);
      auto const change_azimuth = total_change_azimuth * static_cast<double>(i) / static_cast<double>(imgs);

      Z_matrices[i][j].resize(3);
/*
      if (j == 0)
      {*/

        Z_matrices[i][j][0].first = Z_matrices[0][j][0].first;
        Z_matrices[i][j][1].first = Z_matrices[0][j][1].first;
        Z_matrices[i][j][2].first = Z_matrices[0][j][2].first;

        Z_matrices[i][j][0].second = start_radius + change_radius;
        Z_matrices[i][j][1].second = start_inclination + change_inclination;
        Z_matrices[i][j][2].second = start_azimuth + change_azimuth;
     /* }
      else if (j == 1)
      {

        change = i * incr * total_change;
        Z_matrices[i][j][0] = std::make_pair(Z_matrices[0][j][0].first,
          Z_matrices[0][j][0].second + change);

        mega_temp = spherical(atom_B, atom_A, x_ref - atom_A, y_ref - atom_A);

        temporary = { 1, 0, N };
        Z_matrices[i][j][1].first = temporary;

        temporary = { 1, 0, N, N + 1 };
        Z_matrices[i][j][2].first = temporary;

        Z_matrices[i][j][1].second = (double)coords_ini.intern(j).inclination()
          + (double)(coords_final.intern(j).inclination()
            - coords_ini.intern(j).inclination()) * incr * (double)i;
        Z_matrices[i][j][2].second = (double)coords_ini.intern(j).azimuth()
          + (double)(coords_final.intern(j).azimuth()
            - coords_ini.intern(j).azimuth())			* incr * (double)i;
      }
      else if (j == 2)
      {

        total_change = Z_matrices[imgs - 1][j][0].second
          - Z_matrices[0][j][0].second;
        change = i * incr * total_change;
        Z_matrices[i][j][0] = std::make_pair(Z_matrices[0][j][0].first,
          Z_matrices[0][j][0].second + change);

        total_change = Z_matrices[imgs - 1][j][1].second
          - Z_matrices[0][j][1].second;
        if (abs(total_change) >= (180.))
          total_change = total_change - (total_change < 0 ? -360. : 360.);
        change = i * incr * total_change;
        Z_matrices[i][j][1] = std::make_pair(Z_matrices[0][j][1].first,
          Z_matrices[0][j][1].second + change);

        mega_temp = spherical(atom_C, atom_B, atom_A - atom_B, x_ref - atom_B);

        temporary = { 2, 1, 0, N };
        Z_matrices[i][j][2].first = temporary;

        Z_matrices[i][j][2].second = (double)coords_ini.intern(j).azimuth()
          + (double)(coords_final.intern(j).azimuth()
            - coords_ini.intern(j).azimuth()) * incr * (double)i;
      }
      else
      {

        total_change = Z_matrices[imgs - 1][j][0].second
          - Z_matrices[0][j][0].second;
        change = i * incr * total_change;
        Z_matrices[i][j][0] = std::make_pair(Z_matrices[0][j][0].first,
          Z_matrices[0][j][0].second + change);

        total_change = Z_matrices[imgs - 1][j][1].second
          - Z_matrices[0][j][1].second;
        if (abs(total_change) >= (180.))
          total_change = total_change - (total_change < 0 ? -360. : 360.);
        change = i * incr * total_change;
        Z_matrices[i][j][1] = std::make_pair(Z_matrices[0][j][1].first,
          Z_matrices[0][j][1].second + change);

        total_change = Z_matrices[imgs - 1][j][2].second
          - Z_matrices[0][j][2].second;
        if (abs(total_change) >= (180.))
          total_change = total_change - (total_change < 0 ? -360. : 360.);
        change = i * incr * total_change;
        Z_matrices[i][j][2] = std::make_pair(Z_matrices[0][j][2].first,
          Z_matrices[0][j][2].second + change);
      }*/
    }
    Z_matrices_s3[i].resize(N);
    for (size_t j = 0; j < N; ++j)
    {
      Z_matrices_s3[i][j].radius() = Z_matrices[i][j][0].second;
      Z_matrices_s3[i][j].inclination() = (coords::angle_type)Z_matrices[i][j][1].second;
      Z_matrices_s3[i][j].azimuth() = (coords::angle_type)Z_matrices[i][j][2].second;
    }
  }

  coords::Coordinates coords;
  coords = *cPtr;

  /*write_gzmat("Cyclohexan_Zmat3.gzmat", Z_matrices[3], coords_ini);*/

  /*size_t no_dist = Z_matrices[0][N - 1][2].first[1],
    no_angle = Z_matrices[0][N - 1][2].first[2],
    no_dihedral = Z_matrices[0][N - 1][2].first[3];


  coords.adapt_indexation(no_dist, no_angle, no_dihedral,
    Z_matrices[0], cPtr);*/
 

  //system("pause");

  //new coords object for saving newly generated structures
  std::unique_ptr<coords::input::format> format_ptr(coords::input::new_format());
  coords::Coordinates structure;

    ////converts Z-matrix to cartesian structure using OpenBabel
    //write_gzmat("NEB_" + to_string(i) + ".gzmat", Z_matrices[i], coords);

    //std::string command = "obabel -i gzmat NEB_"
    //  + to_string(i)
    //  + ".gzmat -o txyz -O NEB_coordinates_cartesian_"
    //  + to_string(i);
    //system(command.c_str());

    ////saves new structure
    //structure = format_ptr->read("NEB_coordinates_cartesian_"
    //  + to_string(i));
    //structure.adapt_indexation(no_dist, no_angle, no_dihedral,
    //  Z_matrices[0], cPtr);

    //prepares parameters for NEB MEP finding


    /*images[j].x() = structure.xyz(j).x();
    images[j].y() = structure.xyz(j).y();
    images[j].z() = structure.xyz(j).z();
    */

  //converting Z-matrix to cartesian structure, NeRF

  for (size_t i = 1; i < (imgs - 1); ++i)
  { 
    coords::Representation_3D current_images(N);
    coords::Representation_3D CartesianStructure(N);
    
    for (size_t j = 0; j < N; ++j){
        if (j == 0) {
          CartesianStructure[j] = coords::Cartesian_Point(0.,0.,0.);

        }
        else if (j == 1) {
          CartesianStructure[j] = coords::Cartesian_Point(
            Z_matrices[i][j][0].second,
            0.,
            0.
          );

        }
        else if (j == 2) {
          auto const & ABC = Z_matrices[i][j][1].first;
          auto const & A = Z_matrices[i][ABC[0] - 1];
          auto const & B = Z_matrices[i][ABC[1] - 1];
          auto const & C = Z_matrices[i][ABC[2] - 1];

          auto r1 = B[0].second;
          auto r2 = C[0].second;
          auto theta = scon::ang<coords::float_type>::from_deg(C[1].second).radians();

          auto x = r2 * cos(pi - theta);
          auto y = r2 * sin(pi - theta);
          
          CartesianStructure[j] = coords::Cartesian_Point(
            (r1 + x), y, 0.
          );

         }
        else {
          std::unordered_set<std::size_t> check_list;
          auto is_not_in_vec = [&](std::vector<std::size_t> const & vec) -> std::size_t {
            for (auto const & v : vec) {
              if (check_list.insert(v).second) {
                return v;
              }
            }
            return 0;
          };

          auto const & D = is_not_in_vec({ j + 1 })-1;
          auto const & C = is_not_in_vec(Z_matrices[i][j][0].first)-1;
          auto const & B = is_not_in_vec(Z_matrices[i][j][1].first)-1;
          auto const & A = is_not_in_vec(Z_matrices[i][j][2].first)-1;

          auto const & DD = Z_matrices[i][D];

          auto r = DD[0].second;
          auto theta = scon::ang<coords::float_type>::from_deg(DD[1].second).radians();
          auto phi = scon::ang<coords::float_type>::from_deg(DD[2].second).radians();

          auto x = r * cos(phi) * sin(theta);
          auto y = r * sin(phi) * sin(theta);
          auto z = r * cos(theta);

          Eigen::Vector3d D2(z, x, y);
          Eigen::Vector3d Dvec;

          coords::Cartesian_Point BA = normalized(CartesianStructure[B] - CartesianStructure[A]);
          coords::Cartesian_Point BC = normalized(CartesianStructure[B] - CartesianStructure[C]);
          coords::Cartesian_Point N = normalized(cross(BA, BC)); 
          coords::Cartesian_Point NcrBC = normalized(cross(N, BC));

          Eigen::Matrix3d M;
          M << BC.x(), NcrBC.x(), N.x(),
               BC.y(), NcrBC.y(), N.y(),
               BC.z(), NcrBC.z(), N.z();

          Dvec = (M*D2);

          CartesianStructure[j] = coords::Cartesian_Point(Dvec(0), Dvec(1), Dvec(2)) + CartesianStructure[C];

        }

      current_images[j].x() = CartesianStructure[j].x();
      current_images[j].y() = CartesianStructure[j].y();
      current_images[j].z() = CartesianStructure[j].z();

    }
    // Nach Konvertierung in xyz müssen die Images wieder in die alte Reihenfolge gebracht werden:

    for (auto const & n : new_order_ini) {
      images[new_order_ini[n]].x() = current_images[n].x();
      images[new_order_ini[n]].y() = current_images[n].y();
      images[new_order_ini[n]].z() = current_images[n].z();
    }

    for (auto const & j : images) {
      imagi[i].emplace_back(j);
      image_ini[i].emplace_back(j);
      images_initial.emplace_back(j);
    }

  }
  

  coords::Coordinates new_coords = coords;
  new_coords.set_xyz(imagi[1], false);
  coords::output::formats::tinker tinker_writer(new_coords);
  tinker_writer.to_stream(std::cout);



  /*printmono("Cyclohexan_Image2.xyz", imagi[1], 1);*/
  //print("Trideca_All_Images_z_to_xyz.xyz", imagi, 1);
  /*printmono("Cyclohexan_test_image2.xyz", images, 1);*/

  std::ostringstream names;
  names << "IMAGES_INI" << cPtr->mult_struc_counter << ".dat";

  std::ostringstream na;
  na << "IMAGES_START" << this->cPtr->mult_struc_counter << ".arc";
  std::ptrdiff_t s{ 0 };
  print(na.str(), imagi, s);
  if (Config::get().neb.INT_PATH) calc_shift();
}

//creates Z-matrix from redundant internal coords that are only defined for nonterminal atoms

std::vector<std::vector<std::pair<std::vector<size_t>, double>>> neb::redundant_to_Z_backbone(std::vector<std::vector<std::pair<std::vector<size_t>, double>>> &redundant_dists,
  std::vector<std::vector<std::pair<std::vector<size_t>, double>>> &redundant_angles,
  std::vector<std::vector<std::pair<std::vector<size_t>, double>>> &redundant_dihedrals,
  std::vector<size_t> backbone_indeces)
{
  std::vector<std::vector<std::pair<std::vector<size_t>, double>>> Z_matrix;
  std::vector<std::string> atom_labels;
  std::vector<size_t> terminal_enum;																//*terminality indeces
  std::pair<std::vector<size_t>, double> temp, dihedral_vect;				//*for saving different index combinations which define the same geometrical instance
  std::vector<std::pair<std::vector<size_t>, double>> unique_dists,
    unique_angles,
    unique_dihedrals;
  std::vector<bool> defined;																				//for checking if there is an angle or dihedral defined for the current atom
  bool doesnt_exist, linked = false;																//*for checking if parameter to be defined has already been defined for another atom
                                                                    //* and if the atoms are linked properly
  size_t no_dist = 0, no_angle = 0, no_dihedral = 0;								//for saving the indeces of the atoms that don't get any coords, no angle and no dihedral and no dihedral
  size_t N = cPtr->atoms().size();
  if (!(redundant_dists.size() == redundant_angles.size()
    && redundant_angles.size() == redundant_dihedrals.size()))
  {
    throw std::logic_error("Wrong definition of parameter and index vectors.");
  }
  size_t N_main = redundant_dists.size();
  size_t term_deepness;
  Z_matrix.resize(N + 1);
  unique_dists.resize(N);
  unique_angles.resize(N);
  unique_dihedrals.resize(N);
  defined.resize(N);
  for (size_t i = 0; i < N; ++i)
  {
    defined[i] = false;
  }

  terminal_enum = cPtr->terminal_enum();
  term_deepness = *std::max_element(terminal_enum.begin(),
    terminal_enum.end());


  /*for (auto &a : terminal_enum)
  {
  std::cout << a << '\n';
  }
  system("pause");*/


  //parameter reduction
  //finding and treating first three atoms in Z-matrix differently
  for (size_t j = 0; j < N; ++j)
  {
    if (!(terminal_enum[j] == term_deepness
      || N == 3)) continue;
    if (N == 3)
    {
      for (size_t k = 0; k < N; ++k)
      {
        if (terminal_enum[k] == term_deepness)
        {
          no_angle = k;
          no_dist = cPtr->atoms(no_angle).bonds(0);
          no_dihedral = cPtr->atoms(no_angle).bonds(1);
          break;
        }
      }
    }
    else
    {
      no_dist = j;
      for (size_t k = 0; k < cPtr->atoms(no_dist).bonds().size(); ++k)
      {
        if (terminal_enum[cPtr->atoms(no_dist).bonds(k)] == 1) continue;
        no_angle = cPtr->atoms(no_dist).bonds(k);
        break;
      }
      for (size_t k = 0; k < cPtr->atoms(no_angle).bonds().size(); ++k)
      {
        if ((terminal_enum[cPtr->atoms(no_angle).bonds(k)] == 1 && term_deepness > 3)
          || cPtr->atoms(no_angle).bonds(k) == no_dist) continue;
        no_dihedral = cPtr->atoms(no_angle).bonds(k);
        break;
      }
    }
    defined[no_dist] = true;
    temp.first = { no_angle + 1,
      no_dist + 1 };
    for (size_t k = 0; k < redundant_dists[backbone_indeces[no_angle] - 1].size(); ++k)
    {
      if (find(redundant_dists[backbone_indeces[no_angle] - 1][k].first.begin(),
        redundant_dists[backbone_indeces[no_angle] - 1][k].first.end(),
        no_dist + 1) == redundant_dists[backbone_indeces[no_angle] - 1][k].first.end()) continue;
      temp.second = redundant_dists[backbone_indeces[no_angle] - 1][k].second;
      unique_dists[no_angle] = temp;
      defined[no_angle] = true;
      break;
    }
    temp.first = { no_dihedral + 1,
      no_angle + 1 };
    if (terminal_enum[no_dihedral] > 1)
    {
      for (size_t k = 0; k < redundant_dists[backbone_indeces[no_dihedral] - 1].size(); ++k)
      {
        if (find(redundant_dists[backbone_indeces[no_dihedral] - 1][k].first.begin(),
          redundant_dists[backbone_indeces[no_dihedral] - 1][k].first.end(),
          no_angle + 1) == redundant_dists[backbone_indeces[no_dihedral] - 1][k].first.end()) continue;
        temp.second = redundant_dists[backbone_indeces[no_dihedral] - 1][k].second;
        unique_dists[no_dihedral] = temp;
        defined[no_dihedral] = true;
        break;
      }
    }
    else
    {
      for (size_t k = 0; k < redundant_dists[backbone_indeces[no_angle] - 1].size(); ++k)
      {
        if (find(redundant_dists[backbone_indeces[no_angle] - 1][k].first.begin(),
          redundant_dists[backbone_indeces[no_angle] - 1][k].first.end(),
          no_dihedral + 1) == redundant_dists[backbone_indeces[no_angle] - 1][k].first.end()) continue;
        temp.second = redundant_dists[backbone_indeces[no_angle] - 1][k].second;
        unique_dists[no_dihedral] = temp;
        defined[no_dihedral] = true;
        break;
      }
    }
    temp.first = { no_dihedral + 1,
      no_angle + 1,
      no_dist + 1 };
    for (size_t k = 0; k < redundant_angles[backbone_indeces[no_angle] - 1].size(); ++k)
    {
      if (find(redundant_angles[backbone_indeces[no_angle] - 1][k].first.begin(),
        redundant_angles[backbone_indeces[no_angle] - 1][k].first.end(),
        no_dist + 1) == redundant_angles[backbone_indeces[no_angle] - 1][k].first.end()
        || find(redundant_angles[backbone_indeces[no_angle] - 1][k].first.begin(),
          redundant_angles[backbone_indeces[no_angle] - 1][k].first.end(),
          no_dihedral + 1) == redundant_angles[backbone_indeces[no_angle] - 1][k].first.end()) continue;
      temp.second = redundant_angles[backbone_indeces[no_angle] - 1][k].second;
      unique_angles[no_dihedral] = temp;
      break;
    }
    break;
  }

  //gets distances
  for (size_t i = 0; i < term_deepness; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      if (defined[j]) continue;
      if (i + 1 == 1
        && 1 == terminal_enum[j])
      {
        if (cPtr->atoms(j).bonds().size() != 1)
          throw std::logic_error("Terminal atom "
            + to_string(j + 1)
            + "not actually terminal.");
        for (size_t k = 0;
          k < redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1].size();
          ++k)
        {
          if (redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k].first[1] != j + 1) continue;
          temp.first = { redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k].first[1],
            redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k].first[0] };
          temp.second = redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k].second;
          doesnt_exist = find(unique_dists.begin(),
            unique_dists.end(),
            redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k]) == unique_dists.end()
            && find(unique_dists.begin(),
              unique_dists.end(),
              temp) == unique_dists.end();
          linked = cPtr->atoms(j).is_bound_to(redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k].first[0] - 1);
          if (!(doesnt_exist && linked)) continue;
          unique_dists[j] = temp;
          defined[j] = true;
          break;
        }
      }
      else if (terminal_enum[j] == i + 1)
      {
        if (backbone_indeces[j] == 0)
          throw std::logic_error("Backbone atom "
            + to_string(j + 1)
            + " not actually part of backbone.");
        for (size_t k = 0; k < redundant_dists[backbone_indeces[j] - 1].size(); ++k)
        {
          if (terminal_enum[redundant_dists[backbone_indeces[j] - 1][k].first[1] - 1] < terminal_enum[j]) continue;
          temp.first = { redundant_dists[backbone_indeces[j] - 1][k].first[1],
            redundant_dists[backbone_indeces[j] - 1][k].first[0] };
          temp.second = redundant_dists[backbone_indeces[j] - 1][k].second;
          doesnt_exist = find(unique_dists.begin(),
            unique_dists.end(),
            redundant_dists[backbone_indeces[j] - 1][k]) == unique_dists.end()
            && find(unique_dists.begin(),
              unique_dists.end(),
              temp) == unique_dists.end();
          linked = cPtr->atoms(j).is_bound_to(redundant_dists[backbone_indeces[j] - 1][k].first[1] - 1);
          if (!(doesnt_exist && linked)) continue;
          unique_dists[j] = redundant_dists[backbone_indeces[j] - 1][k];
          defined[j] = true;
          break;
        }
      }
      else if (terminal_enum[j] < i + 1 && !defined[j])
      {
        for (size_t k = 0;
          k < redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1].size();
          ++k)
        {
          if (redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k].first[1] != j + 1) continue;
          temp.first = { redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k].first[1],
            redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k].first[0] };
          temp.second = redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k].second;
          doesnt_exist = find(unique_dists.begin(),
            unique_dists.end(),
            redundant_dists[cPtr->atoms(j).bonds(0)][0]) == unique_dists.end()
            && find(unique_dists.begin(),
              unique_dists.end(),
              temp) == unique_dists.end();
          linked = cPtr->atoms(j).is_bound_to(redundant_dists[backbone_indeces[cPtr->atoms(j).bonds(0)] - 1][k].first[0] - 1);
          if (!(doesnt_exist && linked)) continue;
          unique_dists[j] = temp;
          defined[j] = true;
          break;
        }
      }
    }
  }
  //checks, if distances have been defined for every atom
  for (size_t i = 0; i < N; ++i)
  {
    if (!defined[i])
      throw std::runtime_error("No bond parameter could be defined for atom "
        + to_string(i + 1));
    defined[i] = false;
  }
  //gets angles
  defined[no_dist] = true;
  defined[no_angle] = true;
  defined[no_dihedral] = true;
  for (size_t i = 0; i < term_deepness; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      if (defined[j]) continue;
      if (terminal_enum[j] == i + 1)
      {
        for (size_t k = 0;
          k < cPtr->atoms(j).bonds().size() && !defined[j];
          ++k)
        {
          if (terminal_enum[cPtr->atoms(j).bonds(k)] < terminal_enum[j]) continue;
          for (size_t l = 0;
            l < redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1].size();
            ++l)
          {
            if (redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[0] == j + 1
              && redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[2] != j + 1
              && (terminal_enum[redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[2] - 1] > 2 || term_deepness == 2))
            {
              temp.first = { redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[2],
                redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[1],
                j + 1, };
              temp.second = redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].second;
            }
            else if (redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[0] != j + 1
              && redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[2] == j + 1
              && (terminal_enum[redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[0] - 1] > 2 || term_deepness == 2))
            {
              temp.first = { j + 1,
                redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[1],
                redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[0] };
              temp.second = redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].second;
            }
            else
            {
              continue;
            }
            doesnt_exist = find(unique_angles.begin(),
              unique_angles.end(),
              redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l]) == unique_angles.end()
              && find(unique_angles.begin(),
                unique_angles.end(),
                temp) == unique_angles.end();
            linked = cPtr->atoms(cPtr->atoms(j).bonds(k)).is_bound_to(redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[0] - 1)
              && cPtr->atoms(cPtr->atoms(j).bonds(k)).is_bound_to(redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[2] - 1);
            if (!(doesnt_exist && linked)) continue;
            unique_angles[j] = redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l];
            defined[j] = true;
            break;
          }
        }
      }
      else if (terminal_enum[j] < i + 1 && !defined[j])
      {
        for (size_t k = 0; k <
          cPtr->atoms(j).bonds().size() && !defined[j];
          ++k)
        {
          for (size_t l = 0;
            l < redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1].size();
            ++l)
          {
            if (redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[0] == j + 1
              && redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[2] != j + 1)
            {
              temp.first = { redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[2],
                redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[1],
                j + 1 };
              temp.second = redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].second;
            }
            else if (redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[0] != j + 1
              && redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[2] == j + 1)
            {
              temp.first = { j + 1,
                redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[1],
                redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[0] };
              temp.second = redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].second;
            }
            else continue;
            doesnt_exist = find(unique_angles.begin(),
              unique_angles.end(),
              redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l]) == unique_angles.end()
              && find(unique_angles.begin(),
                unique_angles.end(),
                temp) == unique_angles.end();
            linked = cPtr->atoms(cPtr->atoms(j).bonds(k)).is_bound_to(redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[0] - 1)
              && cPtr->atoms(cPtr->atoms(j).bonds(k)).is_bound_to(redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l].first[2] - 1);
            if (!(doesnt_exist && linked)) continue;
            unique_angles[j] = redundant_angles[backbone_indeces[cPtr->atoms(j).bonds(k)] - 1][l];
            defined[j] = true;
            break;
          }
        }
      }
    }
  }
  //checks, if angles have been defined for every atom
  for (size_t i = 0; i < N; ++i)
  {
    if (!defined[i])
      throw std::runtime_error("No angle parameter could be defined for atom "
        + to_string(i + 1));
    defined[i] = false;
  }
  //gets torsions
  defined[no_dist] = true;
  defined[no_angle] = true;
  defined[no_dihedral] = true;
  for (size_t i = 0; i < term_deepness; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      if (defined[j]) continue;
      if (terminal_enum[j] == i + 1
        && i + 1 == 1)
      {
        for (size_t k = 0;
          k < cPtr->atoms(j).bonds().size() && !defined[j];
          ++k)
        {
          if (terminal_enum[cPtr->atoms(j).bonds(k)] == 1) continue;
          for (size_t m = 0;
            m < cPtr->atoms(cPtr->atoms(j).bonds(k)).bonds().size() && !defined[j];
            ++m)
          {
            if (terminal_enum[cPtr->atoms(cPtr->atoms(j).bonds(k)).bonds(m)] == 1) continue;
            for (size_t o = 0;
              o < cPtr->atoms(cPtr->atoms(cPtr->atoms(j).bonds(k)).bonds(m)).bonds().size() && !defined[j];
              ++o)
            {
              if (terminal_enum[cPtr->atoms(cPtr->atoms(cPtr->atoms(j).bonds(k)).bonds(m)).bonds(o)] == 1) continue;
              for (size_t p = 0;
                p < redundant_dihedrals[backbone_indeces[cPtr->atoms(cPtr->atoms(cPtr->atoms(j).bonds(k)).bonds(m)).bonds(o)] - 1].size();
                ++p)
              {
                dihedral_vect.first = redundant_dihedrals[backbone_indeces[cPtr->atoms(cPtr->atoms(cPtr->atoms(j).bonds(k)).bonds(m)).bonds(o)] - 1][p].first;
                if (!(terminal_enum[dihedral_vect.first[0] - 1] != 1
                  || terminal_enum[dihedral_vect.first[3] - 1] != 1)) continue;
                dihedral_vect.second = redundant_dihedrals[backbone_indeces[cPtr->atoms(cPtr->atoms(cPtr->atoms(j).bonds(k)).bonds(m)).bonds(o)] - 1][p].second;
                temp.first = { dihedral_vect.first[3],
                  dihedral_vect.first[2],
                  dihedral_vect.first[1],
                  dihedral_vect.first[0] };
                temp.second = redundant_dihedrals[backbone_indeces[cPtr->atoms(cPtr->atoms(cPtr->atoms(j).bonds(k)).bonds(m)).bonds(o)] - 1][p].second;
                doesnt_exist = find(unique_dihedrals.begin(),
                  unique_dihedrals.end(),
                  dihedral_vect) == unique_dihedrals.end()
                  && find(unique_dihedrals.begin(),
                    unique_dihedrals.end(),
                    temp) == unique_dihedrals.end();
                linked = cPtr->atoms(dihedral_vect.first[0] - 1).is_bound_to(dihedral_vect.first[1] - 1)
                  && cPtr->atoms(dihedral_vect.first[1] - 1).is_bound_to(dihedral_vect.first[2] - 1)
                  && cPtr->atoms(dihedral_vect.first[2] - 1).is_bound_to(dihedral_vect.first[3] - 1);
                if (!(doesnt_exist && linked)) continue;
                if (dihedral_vect.first[3] == j + 1
                  && dihedral_vect.first[0] != dihedral_vect.first[2]
                  && dihedral_vect.first[0] != dihedral_vect.first[3]
                  && dihedral_vect.first[1] != dihedral_vect.first[3])
                {
                  unique_dihedrals[j] = temp;
                }
                else if (dihedral_vect.first[0] == j + 1
                  && dihedral_vect.first[0] != dihedral_vect.first[2]
                  && dihedral_vect.first[0] != dihedral_vect.first[3]
                  && dihedral_vect.first[1] != dihedral_vect.first[3])
                {
                  unique_dihedrals[j] = dihedral_vect;
                }
                else continue;
                defined[j] = true;
                break;
              }
            }
          }
        }
      }
      else if (i + 1 == terminal_enum[j])
      {
        for (auto &a : redundant_dihedrals)
        {
          if (defined[j]) break;
          for (auto &dih : a)
          {
            if (terminal_enum[dih.first[0] - 1] == 1
              || terminal_enum[dih.first[3] - 1] == 1) continue;
            if ((dih.first[0] != j + 1)
              && (dih.first[3] != j + 1)) continue;
            temp.first = { dih.first[3],
              dih.first[2],
              dih.first[1],
              dih.first[0] };
            temp.second = dih.second;
            doesnt_exist = find(unique_dihedrals.begin(),
              unique_dihedrals.end(),
              dih) == unique_dihedrals.end()
              && find(unique_dihedrals.begin(),
                unique_dihedrals.end(),
                temp) == unique_dihedrals.end();
            linked = cPtr->atoms(dih.first[0] - 1).is_bound_to(dih.first[1] - 1)
              && cPtr->atoms(dih.first[1] - 1).is_bound_to(dih.first[2] - 1)
              && cPtr->atoms(dih.first[2] - 1).is_bound_to(dih.first[3] - 1);
            if (!(doesnt_exist && linked)) continue;
            unique_dihedrals[j] = dih;
            defined[j] = true;
            break;
          }
        }
      }
      else if (i + 1 > terminal_enum[j] && !defined[j])
      {
        for (auto &a : redundant_dihedrals)
        {
          if (defined[j]) break;
          for (auto &dih : a)
          {
            if ((dih.first[0] != j + 1)
              && (dih.first[3] != j + 1)) continue;
            temp.first = { dih.first[3],
              dih.first[2],
              dih.first[1],
              dih.first[0] };
            temp.second = dih.second;
            doesnt_exist = find(unique_dihedrals.begin(),
              unique_dihedrals.end(),
              dih) == unique_dihedrals.end()
              && find(unique_dihedrals.begin(),
                unique_dihedrals.end(),
                temp) == unique_dihedrals.end();
            linked = cPtr->atoms(dih.first[0] - 1).is_bound_to(dih.first[1] - 1)
              && cPtr->atoms(dih.first[1] - 1).is_bound_to(dih.first[2] - 1)
              && cPtr->atoms(dih.first[2] - 1).is_bound_to(dih.first[3] - 1);
            if (!(doesnt_exist && linked)) continue;
            unique_dihedrals[j] = dih;
            defined[j] = true;
            break;
          }
        }
      }
    }
  }
  //checks, if torsions have been defined for every atom
  for (size_t i = 0; i < N; ++i)
  {
    if (!defined[i])
      throw std::runtime_error("No torsion parameter could be defined for atom "
        + to_string(i + 1));
  }
  //end of param reduction
  //getting actual coords

  //Needed for NeRF algorithm in i_to_c, but doesn't work properly, yet.
  double constexpr pi = 3.1415926535897932384626433832795;
  coords::Cartesian_Point x_ref(1, 0, 0);
  coords::Cartesian_Point y_ref(0, 1, 0);
  coords::Cartesian_Point z_ref(0, 0, 1);

  //checks, if first three atoms for Z-matrix are bound to each other properly
  if (!((unique_dists[no_angle].first[0] - 1 == no_dist || unique_dists[no_angle].first[1] - 1 == no_dist)
    && (unique_dists[no_dihedral].first[0] - 1 == no_angle || unique_dists[no_dihedral].first[1] - 1 == no_angle)
    && (unique_angles[no_dihedral].first[1] - 1 == no_angle && (unique_angles[no_dihedral].first[0] - 1 == no_dist
      || unique_angles[no_dihedral].first[2] - 1 == no_dist)))) throw std::runtime_error("Z-Matrix couldn't be calculated in NEB.");

  coords::Cartesian_Point atom_A(-unique_dists[no_dihedral].second + unique_dists[no_angle].second * std::sin((pi / 180) * (unique_angles[no_dihedral].second - 90)),
    unique_dists[no_angle].second * std::cos(unique_angles[no_dihedral].second - 90),
    0);
  coords::Cartesian_Point atom_B(-unique_dists[no_dihedral].second,
    0,
    0);
  coords::Cartesian_Point atom_C(0, 0, 0);
  coords::s3 Rep;
  std::vector<size_t> giga_temp;

  for (size_t i = 0; i < N; ++i)
  {
    if (i == no_dist)
    {
      Z_matrix[i].resize(3);
      giga_temp = { 1, N + 1 };
      Z_matrix[i][0].first = giga_temp;
      giga_temp = { 1, N + 1, N + 2 };
      Z_matrix[i][1].first = giga_temp;
      giga_temp = { 1, N + 1, N + 2, N + 3 };
      Z_matrix[i][2].first = giga_temp;
      Rep = scon::spherical(atom_A, x_ref, y_ref - x_ref, z_ref - x_ref);
      Z_matrix[i][0].second = (double)Rep.radius();
      Z_matrix[i][1].second = (double)Rep.inclination();
      Z_matrix[i][2].second = (double)Rep.azimuth();
    }
    else if (i == no_angle)
    {
      Z_matrix[i].resize(3);
      Z_matrix[i][0] = unique_dists[i];
      giga_temp = { 2, 1, N + 1 };
      Z_matrix[i][1].first = giga_temp;
      giga_temp = { 2, 1, N + 1, N + 2 };
      Z_matrix[i][2].first = giga_temp;
      Rep = scon::spherical(atom_B, atom_A, x_ref - atom_A, y_ref - atom_A);
      Z_matrix[i][1].second = (double)Rep.inclination();
      Z_matrix[i][2].second = (double)Rep.azimuth();
    }
    else if (i == no_dihedral)
    {
      Z_matrix[i].resize(3);
      Z_matrix[i][0] = unique_dists[i];
      Z_matrix[i][1] = unique_angles[i];
      giga_temp = { 3, 2, 1, N + 1 };
      Z_matrix[i][2].first = giga_temp;
      Rep = scon::spherical(atom_C, atom_B, atom_A - atom_B, x_ref - atom_B);
      Z_matrix[i][2].second = (double)Rep.azimuth();
    }
    else
    {
      Z_matrix[i].push_back(unique_dists[i]);
      Z_matrix[i].push_back(unique_angles[i]);
      Z_matrix[i].push_back(unique_dihedrals[i]);
    }
  }
  //switching places in Z-matrix so that atom coords with undefined 
  //internal coords are located at the start of the matrix
  std::vector<std::pair<std::vector<size_t>, double>> tump;
  std::vector<size_t> switch_rememberer_pre = { no_dist,
    no_angle,
    no_dihedral };
  size_t ultra_temp = 0;
  tump = Z_matrix[0];
  Z_matrix[0] = Z_matrix[no_dist];
  Z_matrix[no_dist] = tump;
  if (no_angle == 0)
  {
    ultra_temp = no_angle;
    no_angle = no_dist;
    no_dist = ultra_temp;
  }
  else if (no_dihedral == 0)
  {
    ultra_temp = no_dihedral;
    no_dihedral = no_dist;
    no_dist = ultra_temp;
  }
  tump = Z_matrix[1];
  Z_matrix[1] = Z_matrix[no_angle];
  Z_matrix[no_angle] = tump;
  if (no_dihedral == 1)
  {
    ultra_temp = no_dihedral;
    no_dihedral = no_angle;
    no_angle = ultra_temp;
  }
  tump = Z_matrix[2];
  Z_matrix[2] = Z_matrix[no_dihedral];
  Z_matrix[no_dihedral] = tump;
  std::vector<size_t> switch_rememberer_post = { no_dist,
    no_angle,
    no_dihedral };

  //for future reference
  Z_matrix[N].push_back(std::make_pair(switch_rememberer_post, 0));

  //switching indexation in Z-matrix, since atom numbers have been changed
  for (size_t i = 1; i < N; ++i)
  {
    for (auto &cdist : Z_matrix[i][0].first)
    {
      if (cdist == switch_rememberer_post[0] + 1)
      {
        cdist = 1;
      }
      else if (cdist == switch_rememberer_post[1] + 1)
      {
        cdist = 2;
      }
      else if (cdist == switch_rememberer_post[2] + 1)
      {
        cdist = 3;
      }
      else if (cdist == 1)
      {
        cdist = switch_rememberer_post[0] + 1;
      }
      else if (cdist == 2)
      {
        cdist = switch_rememberer_post[1] + 1;
      }
      else if (cdist == 3)
      {
        cdist = switch_rememberer_post[2] + 1;
      }
    }
    if (i == 1) continue;
    for (auto &cangle : Z_matrix[i][1].first)
    {
      if (cangle == switch_rememberer_post[0] + 1)
      {
        cangle = 1;
      }
      else if (cangle == switch_rememberer_post[1] + 1)
      {
        cangle = 2;
      }
      else if (cangle == switch_rememberer_post[2] + 1)
      {
        cangle = 3;
      }
      else if (cangle == 1)
      {
        cangle = switch_rememberer_post[0] + 1;
      }
      else if (cangle == 2)
      {
        cangle = switch_rememberer_post[1] + 1;
      }
      else if (cangle == 3)
      {
        cangle = switch_rememberer_post[2] + 1;
      }
    }
    if (i == 2) continue;
    for (auto &cdih : Z_matrix[i][2].first)
    {
      if (cdih == switch_rememberer_post[0] + 1)
      {
        cdih = 1;
      }
      else if (cdih == switch_rememberer_post[1] + 1)
      {
        cdih = 2;
      }
      else if (cdih == switch_rememberer_post[2] + 1)
      {
        cdih = 3;
      }
      else if (cdih == 1)
      {
        cdih = switch_rememberer_post[0] + 1;
      }
      else if (cdih == 2)
      {
        cdih = switch_rememberer_post[1] + 1;
      }
      else if (cdih == 3)
      {
        cdih = switch_rememberer_post[2] + 1;
      }
    }
  }
  return Z_matrix;
}

void neb::write_gzmat(std::string const &filename,
  std::vector<std::vector<std::pair<std::vector<size_t>, double>>> const &Z_matrix,
  coords::Coordinates const &coords) const
{
  size_t N = cPtr->atoms().size();
  std::ofstream stream;
  stream.open(filename);
  stream << "#Put Keywords Here, check Charge and Multiplicity.\n\n \n\n0  1\n";
  for (size_t i = 0; i < N; ++i)
  {
    if (i == 0)
    {
      stream << coords.atoms(i).symbol() << '\n';
    }
    else if (i == 1)
    {
      stream << coords.atoms(i).symbol() << "  " << coords.atoms(i).ibond() + 1 << "  ";
      stream << 'r' << i + 1 << '\n';
    }
    else if (i == 2)
    {
      stream << coords.atoms(i).symbol() << "  " << coords.atoms(i).ibond() + 1 << "  ";
      stream << 'r' << i + 1 << "  " << coords.atoms(i).iangle() + 1 << "  ";
      stream << 'a' << i + 1 << '\n';
    }
    else
    {
      stream << coords.atoms(i).symbol() << "  " << coords.atoms(i).ibond() + 1 << "  ";
      stream << 'r' << i + 1 << "  " << coords.atoms(i).iangle() + 1 << "  ";
      stream << 'a' << i + 1 << "  " << coords.atoms(i).idihedral() + 1 << "  ";
      stream << 'd' << i + 1 << '\n';
    }
  }
  stream << "Variables:\n";
  double temp = 0;
  for (size_t i = 1; i < N; ++i)
  {
    stream << 'r' << i + 1 << "= " << Z_matrix[i][0].second << '\n';
    if (i == 1) continue;
    stream << 'a' << i + 1 << "= " << Z_matrix[i][1].second << '\n';
    if (i == 2) continue;
    if (Z_matrix[i][2].second < 0.) temp = Z_matrix[i][2].second + 360.;
    else temp = Z_matrix[i][2].second;
    stream << 'd' << i + 1 << "= " << temp << '\n';
  }
}