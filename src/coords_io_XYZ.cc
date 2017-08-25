/**
CAST 3
Purpose: Reading from XYZ-files

@author Susanne Sauer
@version 1.0
*/
#pragma once
#include "coords_io.h"
#include "helperfunctions.h"


coords::Coordinates coords::input::formats::xyz::read(std::string file)
{
    std::cout<<"WARNING!!!\n";
    std::cout<<"You are reading the structure from yxz-file.\n";
    std::cout<<"No atom types are read so you can't use the structure with a forcefield interface!\n";

    Coordinates coord_object;
    double N = 3; // number of atoms

    Atom o("O");
    Atom h1("H");
    Atom h2("H");
    
    // because of this no forcefield interfaces are available with XYZ inputfile
    o.set_energy_type(0);
    h1.set_energy_type(0);
    h2.set_energy_type(0);

    atoms.add(o);
    atoms.add(h1);
    atoms.add(h2);
    
    Representation_3D positions;

    position.x() = -7.734581;
    position.y() = 3.971126;
    position.z() = 1.221690;
    positions.push_back(position);

    position.x() = -6.923463;
    position.y() = 3.427050;
    position.z() = 1.081216;
    positions.push_back(position);

    position.x() = -8.470271;
    position.y() = 3.365787;
    position.z() = 0.832401;
    positions.push_back(position);

    input_ensemble.push_back(positions);
    coords::PES_Point x(input_ensemble[0u]);

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

    coord_object.init_swap_in(atoms, x);

    for (auto & p : input_ensemble)
    {
      p.gradient.cartesian.resize(p.structure.cartesian.size());
      coord_object.set_xyz(p.structure.cartesian);
      coord_object.to_internal_light();
      p = coord_object.pes();
    }
    
    return coord_object;
}