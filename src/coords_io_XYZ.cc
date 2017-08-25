/**
CAST 3
Purpose: Reading structures from XYZ-files
no atom types are assigned
bonds are created by distance criterion (1.2 times sum of covalent radii)

@author Susanne Sauer
@version 1.0
*/
#pragma once
#include "coords_io.h"
#include "helperfunctions.h"

/**function that reads the structure
@ param file: name of the xyz-file
@ return: Coordinates object that is created out of file*/
coords::Coordinates coords::input::formats::xyz::read(std::string file)
{
    std::cout<<"WARNING!!!\n";
    std::cout<<"You are reading the structure from yxz-file.\n";
    std::cout<<"No atom types are read so you can't use the structure with a forcefield interface!\n";

    Coordinates coord_object;
    std::ifstream config_file_stream(file.c_str(), std::ios_base::in);  // read file to ifstream
    
    std::string line, element;  // a few variables
    double x,y,z;
    Representation_3D positions;

    std::getline(config_file_stream, line);  // first line: number of atoms
    double N = std::stoi(line); 
    
    std::getline(config_file_stream, line);  // discard second line (comment)
    
    while (config_file_stream >> element >> x >> y >> z)  // for every line
    {
        // create atom
        Atom current(element);  
        current.set_energy_type(0);  // because of this no forcefield interfaces are available with XYZ inputfile
        atoms.add(current);
        
        // create position
        position.x() = x;
        position.y() = y;
        position.z() = z;
        positions.push_back(position);
    }

    input_ensemble.push_back(positions);  // fill the positions into PES_Point
    coords::PES_Point pes(input_ensemble[0u]);

    // loop over all atompairs and bind them if they fulfill distance criterion 
    // i.e. the distance is smaller than 1.2 * sum of covalent radiuses
    for (unsigned i=0; i<N; i++)
    {
        for(unsigned j=0; j<i; j++)
        {
            double d = dist(positions[i], positions[j]);
            double d_max = 1.2*(atoms.atom(i).cov_radius() + atoms.atom(j).cov_radius());
            if (d < d_max)
            {
                atoms.atom(i).bind_to(j);
                atoms.atom(j).bind_to(i);
            }
        }
    }

    coord_object.init_swap_in(atoms, pes);  // fill atoms and positions into coord_object

    for (auto & p : input_ensemble)  // do some important stuff (see coords_io_AMBER.cc)
    {
      p.gradient.cartesian.resize(p.structure.cartesian.size());
      coord_object.set_xyz(p.structure.cartesian);
      coord_object.to_internal_light();
      p = coord_object.pes();
    }
    
    return coord_object;
}