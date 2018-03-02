#pragma once

#include <iostream>
#include "coords.h"
#include "coords_io.h"

class ic_testing
{
public: 
	coords::Coordinates* cPtr;
	coords::Representation_3D inp_struc_cartesian;
	coords::Representation_Internal inp_struc_internal;
	coords::Representation_Main inp_struc_main;

	void ic_execution(std::shared_ptr<coords::input::format>, coords::Coordinates &);
};