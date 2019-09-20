#include "find_as.h"
#include "coords_io.h"
#include "helperfunctions.h"

void find_as(coords::Coordinates const& coords, std::string const& filename)
{
  // some stuff before really doing something
  std::ofstream outputfile(filename);
  auto counter{ 0u };
  coords::Atoms atoms = coords.atoms();

  // find aminoacids
  coords::AtomtypeFinder atf(atoms);
  auto amino_acids = atf.get_aminoacids();   
  for (auto& as : amino_acids)
  {
    // find name for aminoacid
    as.determine_aminoacid(atoms);       
    if (as.get_res_name() != "XXX") as.correct_residue_names(atoms);

    // write into file
    counter++;
    if (Config::get().stuff.find_as_md_regions) outputfile << "MDregion     ";
    outputfile << as.get_res_name() << "_" << counter << "    ";
    std::vector<std::size_t> vec(as.get_indices());               // in order to use for_each we have to copy the vector as get_indices() is const
    std::for_each(vec.begin(), vec.end(), [](std::size_t& i) {return ++i; }); // increment each element of vector -> conversion to tinkernumbering
    outputfile << vec_to_string(vec, ",") << "\n";
  }

  // find molecules that don't consist of aminoacids
  for (auto m : coords.molecules()) {
    for (auto a : m) {
      if (atf.recognized_atom(a) == false)   // for every molecule that contains atoms that are not recognized
      {
        // write into file
        counter++;
        if (Config::get().stuff.find_as_md_regions) outputfile << "MDregion     ";
        outputfile << coords.molecule_name(m) << "_" << counter << "    ";
        std::for_each(m.begin(), m.end(), [](std::size_t& i) {return ++i; }); // increment each element of vector -> conversion to tinkernumbering
        outputfile << vec_to_string(m, ",") << "\n";
        break;
      }
    }
  }
}
