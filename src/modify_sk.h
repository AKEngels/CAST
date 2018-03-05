/**
CAST 3
modify_sk.h
Purpose: do things with Slater-Koster files
    - modify Slater-Koster files from dftb.org format to gaussian format
    - read angular momenta from Slater-Koster files

@author Susanne Sauer
@version 1.0
*/
#pragma once

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

/**returns highest angular momentum for an element (from slater-koster files)
@param s: element symbol*/
char angular_momentum_by_symbol(std::string s);