#pragma once
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "configuration.h"
#include "coords_io.h"
#include "energy_int_gaussian.h"


namespace couplings {

	struct coupling {
	public:
		//MO, excitation energies and dipolemoments
		std::vector <double> c_occMO, c_virtMO, c_excitE;
		coords::Representation_3D  c_ex_ex_trans, c_gz_ex_trans;
		std::vector <int> c_state_i, c_state_j, c_gz_i_state;

		void kopplung();

		void INDO(coords::Coordinates, std::string, std::string, std::string);
		void ZINDO(coords::Coordinates, std::string, std::string, std::string);
		void write();

		//kopplungen
		std::vector <double> V_el, V_hole, V_ex, V_ct, V_rek;//actual variables for couplings
		std::vector <int> pSC_homo_1, pSC_homo_2, nSC_homo_1, nSC_homo_2, hetero_pSC, hetero_nSC;//vectors to keep track between which monomers the coupling is

	};
}