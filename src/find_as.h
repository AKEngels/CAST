/**
CAST 3
find_as.h
Purpose: perform task FIND_AS

@author Susanne Sauer
@version 1.0
*/
#pragma once

#include "coords.h"

/**main function for task FIND_AS
@param coords: Coordinates object
@param filename: name of the outputfile*/
void find_as(coords::Coordinates const & coords, std::string const& filename);