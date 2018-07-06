#ifndef ENERGY_INT_PSI_H
#define ENERGY_INT_PSI_H

#include<vector>
#include<utility>
#include<algorithm>
#include<sstream>

#include"coords_io.h"
#include"energy.h"
#include"coords.h"

namespace energy{
  namespace interfaces{
    namespace psi4{
      using coords::float_type;
      class sysCallInterface final: public interface_base{
      public:
        sysCallInterface(coords::Coordinates* coords_ptr)
          : interface_base(coords_ptr),
          tmp_file_name(create_random_file_name(Config::get().general.outputFilename))
        {
          optimizer = true;
        }
        sysCallInterface(sysCallInterface const& other, coords::Coordinates* coord)
          : interface_base(coord),
          tmp_file_name(create_random_file_name(Config::get().general.outputFilename))
        {
          interface_base::operator=(other);
          optimizer = true;
        }

        void swap(interface_base & other) final override;
        interface_base * clone(coords::Coordinates* coord_object) const final override;
        interface_base * move(coords::Coordinates* coord_object) final override;

        void update(bool const) final override{
          /*std::stringstream ss;
          ss << "No need for an update function. \n";
          ss << std::boolalpha <<  "Skip Topology set to " << skip_topology << "\n";
          throw std::runtime_error(ss.str());*/
        };

				/**calculates energy*/
        float_type e(void) final override;
				/**calculates gradients*/
        float_type g(void) final override;
				/**should calculate hessian, in the moment does nothing*/
        float_type h(void) final override;
				/**performs an optimization*/
        float_type o(void) final override;

				/**should print total energy, but I think it prints partial energies as well as print_E_short()*/
        void print_E(std::ostream&) const final override;
				/**prints headline for partial energies*/
        void print_E_head(std::ostream&, bool const endline=true) const final override;
				/**prints partial energies*/
        void print_E_short(std::ostream&, bool const endline=true) const final override;
				/**does nothing*/
        void to_stream(std::ostream&) const final override;

				/**reads mulliken charges from output file*/
				std::vector<float_type> charges() const override;
				/**calculates gradients on external charges
				uses coulomb potential between external charge and mulliken charges of atoms*/
				std::vector<coords::Cartesian_Point> get_g_ext_chg() const override;

      private:

        /**possible calculation types*/
        enum class Calc {
          energy,
          gradient,
          optimize
        };

				/**writes input file
				@param kind: calculation type (e.g. energy or gradient)*/
				void write_input(Calc kind = Calc::energy) const;
				/**writes everything in input file that doesn't depend on calculation type
        @param os: stream where external charges are written to (directed into PSI4 inputfile)*/
        void write_head(std::ostream&) const;
				/**writes molecule to input file
        @param os: stream where external charges are written to (directed into PSI4 inputfile)*/
        void write_molecule(std::ostream&) const;
				/**writes external charges into input file and positions of external charges into file 'grid.dat' for reading external field
        (needed for QM/MM) 
        @param os: stream where external charges are written to (directed into PSI4 inputfile)*/
				void write_ext_charges(std::ostream & os) const;
				/**writes everything except for head in inputfile for energy calculation
        @param os: stream where external charges are written to (directed into PSI4 inputfile)*/
        void write_energy_input(std::ostream&) const;
				/**writes everything except for head in inputfile for gradients calculation
        @param os: stream where external charges are written to (directed into PSI4 inputfile)*/
        void write_gradients_input(std::ostream&) const;
				/**writes everything except for head in inputfile for optimization
        @param os: stream where external charges are written to (directed into PSI4 inputfile)*/
        void write_optimize_input(std::ostream&) const;

				/**makes system call to psi4, throws error if this fails 3 times*/
        void make_call()const;

        std::vector<std::string> parse_specific_position(std::istream& is, std::string const& delim, int space) const;
        std::vector<std::string> get_last_gradients() const;
				/**reads final geometry from output file of optimization*/
        coords::Representation_3D get_final_geometry() const;

				/**reads and calculates the energy from outputfile*/
        coords::float_type parse_energy();
        std::pair<coords::float_type, coords::Representation_3D> parse_gradients();
        std::tuple<coords::float_type, coords::Representation_3D, coords::Representation_3D>
        parse_geometry_and_gradients();

        template<typename StrCont>
        coords::Representation_3D extract_Rep3D(StrCont && lines)const;

				/**randomly created name of input and output file*/
        std::string tmp_file_name;

				/**partial energies (in hartree!!!)*/
        std::vector<std::pair<std::string, float_type>> energies;//<- energy in Hartree!
      };
    }
  }
}

template<typename StrCont>
inline coords::Representation_3D energy::interfaces::psi4::sysCallInterface::extract_Rep3D(StrCont && lines)const{
  coords::Representation_3D ret;
  for(auto && line: lines){
    std::istringstream iss{line};
    std::vector<std::string> words{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
    ret.emplace_back(coords::Cartesian_Point(
      std::stod(words.at(1)), std::stod(words.at(2)), std::stod(words.at(3))
    ));
  }
  return ret;
}

#endif
