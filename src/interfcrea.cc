#include<string>
#include<memory>
#include<iomanip>
#include <iostream>
#include "interfcrea.h"
#include "coords_io.h"
#include "coords.h"


//Function to read a second structure and add it to the structure read by standart input

coords::Coordinates interface_creation(char iaxis, double idist, coords::Coordinates norm_coord, coords::Coordinates add_coords)//norm_coord is coord object from read standart input in main.cc
{

  std::ofstream new_interface(Config::set().general.outputFilename, std::ios_base::out);

  std::size_t origN = norm_coord.xyz().size();
  std::size_t Nges = origN + add_coords.xyz().size();//Number of atoms in new structure file

  double displ; //displacement of atoms from add_coords in new interface


  if (iaxis == 'x') //Check fom which direction the other structure shall come
  {
    double maxx(0.);
    for (std::size_t i = 0u; i < norm_coord.xyz().size(); i++) //loop to find biggest x value in standart input
    {
      if (norm_coord.xyz(i).x() > maxx) maxx = norm_coord.xyz(i).x();    
    } 
    displ = maxx + idist;

    new_interface << Nges << '\n';

    for (std::size_t i = 0u; i < norm_coord.xyz().size(); i++) //loop for writing standard input into new interface structure file
    {
      new_interface << std::right << std::fixed << std::setprecision(7) << std::setw(4) << i + 1 << std::setw(6) << norm_coord.atoms(i).symbol() 
        << std::setw(13) << norm_coord.xyz(i).x() << std::setw(13) << norm_coord.xyz(i).y() << std::setw(13) << norm_coord.xyz(i).z() 
        << std::setw(8) << norm_coord.atoms(i).energy_type();

      for (std::size_t j = 0u; j < norm_coord.atoms(i).bonds().size(); j++) // loop for writing bonding information
      {
        new_interface << std::right << std::setw(7) << norm_coord.atoms(i).bonds(j) + 1;
      }
      new_interface << '\n';
    }// end standart input writing loop

    for (std::size_t i = 0u; i < add_coords.xyz().size(); i++)//loop for writing second structure
    {
      new_interface << std::right << std::fixed << std::setprecision(7) << std::setw(4) << i + origN + 1 << std::setw(6) << add_coords.atoms(i).symbol() 
        << std::setw(13) << add_coords.xyz(i).x() + displ << std::setw(13) << add_coords.xyz(i).y() << std::setw(13) << add_coords.xyz(i).z() 
        << std::setw(8) << add_coords.atoms(i).energy_type();

      for (std::size_t j = 0u; j < add_coords.atoms(i).bonds().size(); j++) // loop for writing bonding information
      {
        new_interface << std::right << std::setw(7) << add_coords.atoms(i).bonds(j) + origN + 1;
      }
      new_interface << '\n';
    }

  }
  else if (iaxis == 'y')
  {
    double maxy(0.);
    for (std::size_t i = 0u; i < norm_coord.xyz().size(); i++) //loop to find biggest x value in standart input
    {
      if (norm_coord.xyz(i).y() > maxy) maxy = norm_coord.xyz(i).y();
    }
    displ = maxy + idist;

    new_interface << Nges << '\n';

    for (std::size_t i = 0u; i < norm_coord.xyz().size(); i++) //loop for writing standard input into new interface structure file
    {
      new_interface << std::right << std::fixed << :: std::setprecision(7) << std::setw(4) << i + 1 << std::setw(6) << norm_coord.atoms(i).symbol() 
        << std::setw(13) << norm_coord.xyz(i).x() << std::setw(13) << norm_coord.xyz(i).y() << std::setw(13) << norm_coord.xyz(i).z() 
        << std::setw(8) << norm_coord.atoms(i).energy_type();

      for (std::size_t j = 0u; j < norm_coord.atoms(i).bonds().size(); j++) // loop for writing bonding information
      {
        new_interface << std::right << std::setw(7) << norm_coord.atoms(i).bonds(j) + 1;
      }
      new_interface << '\n';
    }// end standart input writing loop

    for (std::size_t i = 0u; i < add_coords.xyz().size(); i++)//loop for writing second structure
    {
      new_interface << std::right << std::fixed << ::std::setprecision(7) << std::setw(4) << i + origN + 1 << std::setw(6) << add_coords.atoms(i).symbol() 
        << std::setw(13) << add_coords.xyz(i).x() << std::setw(13) << add_coords.xyz(i).y() + displ << std::setw(13) << add_coords.xyz(i).z() 
        << std::setw(8) << add_coords.atoms(i).energy_type();

      for (std::size_t j = 0u; j < add_coords.atoms(i).bonds().size(); j++) // loop for writing bonding information
      {
        new_interface << std::right << std::setw(7) << add_coords.atoms(i).bonds(j) + origN + 1;
      }
      new_interface << '\n';
    }

  }
  else if (iaxis == 'z')
  {
    double maxz(0.);
    for (std::size_t i = 0u; i < norm_coord.xyz().size(); i++) //loop to find biggest x value in standart input
    {
      if (norm_coord.xyz(i).z() > maxz) maxz = norm_coord.xyz(i).z();
    }
    displ = maxz + idist;

    new_interface << Nges << '\n';

    for (std::size_t i = 0u; i < norm_coord.xyz().size(); i++) //loop for writing standard input into new interface structure file
    {
      new_interface << std::right << std::fixed << ::std::setprecision(7) << std::setw(4) << i + 1 << std::setw(6) << norm_coord.atoms(i).symbol() 
        << std::setw(13) << norm_coord.xyz(i).x() << std::setw(13) << norm_coord.xyz(i).y() << std::setw(13) << norm_coord.xyz(i).z() 
        << std::setw(8) << norm_coord.atoms(i).energy_type();

      for (std::size_t j = 0u; j < norm_coord.atoms(i).bonds().size(); j++) // loop for writing bonding information
      {
        new_interface << std::right << std::setw(7) << norm_coord.atoms(i).bonds(j) + 1;
      }
      new_interface << '\n';
    }// end standart input writing loop

    for (std::size_t i = 0u; i < add_coords.xyz().size(); i++)//loop for writing second structure
    {
      new_interface << std::right << std::fixed << ::std::setprecision(7) << std::setw(4) << i + origN + 1 << std::setw(6) << add_coords.atoms(i).symbol() 
        << std::setw(13) << add_coords.xyz(i).x() << std::setw(13) << add_coords.xyz(i).y() << std::setw(13) << add_coords.xyz(i).z() + displ 
        << std::setw(8) << add_coords.atoms(i).energy_type();

      for (std::size_t j = 0u; j < add_coords.atoms(i).bonds().size(); j++) // loop for writing bonding information
      {
        new_interface << std::right << std::setw(7) << add_coords.atoms(i).bonds(j) + origN + 1;
      }
      new_interface << '\n';
    }

  }
  else
  {
    throw std::runtime_error("Entered invalid dimension. Only x,y and z possible.");
  }

  new_interface.close();

  std::unique_ptr<coords::input::format> new_strukt_uptr(coords::input::new_interf_format());
  coords::Coordinates new_interf(new_strukt_uptr->read(Config::set().general.outputFilename));

  return new_interf;

}