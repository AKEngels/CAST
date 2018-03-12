#pragma once

#include <iostream>
#include "coords.h"
#include "coords_io.h"
#include "ic_core.h"

class ic_testing
{
public: 
	coords::Coordinates* cPtr;
	coords::Representation_3D inp_struc_cartesian;
	coords::Representation_Internal inp_struc_internal;
	coords::Representation_Main inp_struc_main;

	void ic_execution(coords::DL_Coordinates & coords);
};