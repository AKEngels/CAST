#include "modify_sk.h"

std::vector<std::vector<std::string>> find_pairs(coords::Coordinates coordobj)
{
	std::vector<std::vector<std::string>> pairs;
	for (auto i = 0u; i < coordobj.atoms().size(); i++)
	{
		for (auto j = 0u; j <= i; j++)
		{
			std::string s1 = coordobj.atoms(i).symbol();
			std::string s2 = coordobj.atoms(j).symbol();
			std::vector<std::string> pair;
			pair.push_back(s1);
			pair.push_back(s2);
			if (!is_in(pair, pairs))
			{
				pairs.push_back(pair);
				std::vector<std::string> pair_r;
				pair_r.push_back(s2);
				pair_r.push_back(s1);
				if (!is_in(pair_r, pairs))
				{
					pairs.push_back(pair_r);
				}
			}
		}
	}
	return pairs;
}

bool modify_file(std::vector<std::string> pair)
{
	std::string filename = pair[0] + "-" + pair[1] + ".skf";

	// look if slater koster file exists
	if (file_exists(filename) == false)
	{
		std::cout << "WARNING! Slater Koster file " << filename << " does not exist. Please download it from dftb.org!\n";
		return false;
	}

	// get atomic numbers for elements
	int atom_number1 = atomic::atomic_number_by_symbol(pair[0]);
	int atom_number2 = atomic::atomic_number_by_symbol(pair[1]);

	// read slater koster file
	std::ifstream file_stream(filename.c_str(), std::ios_base::in);
	std::string line;

	std::getline(file_stream, line); // get first line
	std::vector<std::string> linevec = split(line, ',');
	double starting_number = std::stod(linevec[0]); // get first number out of first line
	bool search = true; //still search for first line containing grid points
	int begin, end;
	int counter = 1;
	std::vector<std::string> stringvec;
	while (file_stream)
	{
		std::getline(file_stream, line);
		stringvec.push_back(line);  // save every line from the second on in a vector (only first line is modified)
		counter += 1;
		if (search && is_in('*', line))  // where is the first grid point?
		{
			begin = counter;
			search = false;
		}
		if (search == false && line == "Spline")  // where is the end of the grid points?
		{
			end = counter;
		}
	}

	// test if line before spline is empty and if yes, remove it
	for (auto i = 0u; i < stringvec.size(); ++i)
	{
		if (stringvec[i] == "Spline")
		{
			if (stringvec[i - 1] == "")
			{
				end = end - 1;                           // decrement number of grid points
				stringvec.erase(stringvec.begin() + i - 1);  // delete line from vector
			}
			break;
		}
	}

	// write all lines from the second on in a new string
	std::string new_string = "";
	for (auto line : stringvec) new_string += "\n" + line;

	// construct new first line
	int second_number = end - begin;  // second number = number of lines with grid points
	std::string new_first_line = std::to_string(starting_number) + "  " + std::to_string(second_number) + "  " + std::to_string(atom_number1) + " " + std::to_string(atom_number2);

	// delete old file
	std::remove(filename.c_str());

	// put new first line and the rest of the old file into a new file
	std::string new_filestring = new_first_line + new_string;
	std::ofstream new_file(filename.c_str());
	new_file << new_filestring;

	if (Config::get().general.verbosity > 2) std::cout << "Successfully converted " << filename << "!\n";
	return true;
}

/**returns highest angular momentum for an element (from slater-koster file for homonuclear atom pair)
@param s: element symbol*/
char angular_momentum_by_symbol(std::string s)
{
	std::string filename = Config::get().energy.dftb.sk_files + "/" + s + "-" + s + ".skf";

	std::string line, angulars;
	int begin, end;

	if (file_exists(filename) == true)
	{
		std::ifstream file_stream(filename.c_str(), std::ios_base::in);
		while (!file_stream.eof())
		{
			std::getline(file_stream, line);
			if (line.find("<Shells>") != std::string::npos)
			{
				begin = int(line.find("<Shells>")) + 8;
				end = int(line.find("</Shells>"));
				angulars = line.substr(begin, end - begin);

				if (is_in('f', angulars) == true) return 'f';
				else if (is_in('d', angulars) == true) return 'd';
				else if (is_in('p', angulars) == true) return 'p';
				else if (is_in('s', angulars) == true) return 's';
				else throw std::runtime_error("Something went wrong. No angular momentum for element " + s + " found.\n");
			}
		}
		throw std::runtime_error("Something went wrong. No angular momentum for element " + s + " found.\n");
	}
	else throw std::runtime_error("No Slater Koster file for element " + s + " found.\n");
}

double get_zeta()
{
	std::string filename = Config::get().energy.dftb.sk_files + "/dftb3.info";
	if (file_exists(filename) == true)
	{
		std::ifstream file_stream(filename.c_str(), std::ios_base::in);
		std::string line;
		std::getline(file_stream, line);
		std::vector<std::string> linevec = split(line, ' ', true);
		return std::stod(linevec[1]);
	}
	else throw std::runtime_error("No dftb3.info file found.\n");
}

double hubbard_deriv_by_symbol(std::string s)
{
	std::string filename = Config::get().energy.dftb.sk_files + "/dftb3.info";
	if (file_exists(filename) == true)
	{
		std::ifstream file_stream(filename.c_str(), std::ios_base::in);
		std::string line;

		while (!file_stream.eof())
		{
			std::getline(file_stream, line);
			std::vector<std::string> linevec = split(line, ' ', true);
			if (linevec[0] == s) return std::stod(linevec[1]);
		}
		throw std::runtime_error("Something went wrong. No hubbard derivative for element " + s + " found.\n");
	}
	else throw std::runtime_error("No dftb3.info file found.\n");
}