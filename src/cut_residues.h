#pragma once 
#include<string>
#include "coords.h"
#include "coords_io.h"
#include "helperfunctions.h"

/**takes a coordinates object, cuts out residues around a given reaction site
and writes these residues into a new tinkerstructure
indizes of 'reactive atoms' and a cutoff radius have to be given by inpufile
furthermore it writes all atoms of protein backbone into new file "fix.txt"
@param c: coordinates object
@param s: stream where the tinkerstucture is written to*/
void cut_residues(coords::Coordinates c, std::ostream & s);