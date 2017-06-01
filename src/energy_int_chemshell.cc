#include "energy_int_chemshell.h"

std::string energy::interfaces::chemshell::sysCallInterface::create_pdb() {
	auto inp_file_name = tmp_file_name;
	auto xyz_file_name = inp_file_name;
	xyz_file_name.append(".xyz");
	
	std::ofstream xyz_file(xyz_file_name);
	xyz_file << coords->xyz();

	std::stringstream ss;
	ss << "babel -ixyz " << xyz_file_name << " -opdb " << inp_file_name << ".pdb";

	scon::system_call(ss.str());
}