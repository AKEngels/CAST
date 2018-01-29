﻿/**
CAST 3
modify_sk.h
Purpose: modify Slater-Koster files from dftb.org format to gaussian format

@author Susanne Sauer
@version 1.0
*/

#include "helperfunctions.h"
#include "atomic.h"

/**finds all pairs of elements of a given structure
returns a vector of these pairs
every pair is a vector of the two element symbols*/
std::vector<std::vector<std::string>> find_pairs(coords::Coordinates);

/**modifies the slater koster file for a given element pair 
(details see here: http://gaussian.com/dftb/, Modifying Slater-Koster Files)
@param pair: vector of the two element symbols
returns true if modifying was successfull,
returns false if not (e.g. if slater koster file for given element pair is not found by CAST)*/
bool modify_file(std::vector<std::string> pair);