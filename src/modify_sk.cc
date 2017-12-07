#include "modify_sk.h"

std::vector<std::vector<std::string>> find_pairs(coords::Coordinates coordobj)
{
  std::vector<std::vector<std::string>> pairs;
  for (int i = 0; i < coordobj.atoms().size(); i++)
  {
    for (int j = 0; j < i; j++)
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
    std::cout << "Slater Koster file " << filename << " does not exist. Please download it from dftb.org!\n";
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

  int second_number, counter = 1;
  std::string new_string = "";
  while (file_stream)
  {
    std::getline(file_stream, line);
    new_string += "\n"+line;  // write every line from the second on in a new string (only first line is modified)
    counter += 1;
    if (line == "Spline")  // get second number for new first line
    {
      if (pair[0] == pair[1]) second_number = counter - 3;  
      else second_number = counter - 2;
    }
  }

  // construct new first line
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