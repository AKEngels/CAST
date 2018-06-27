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

        void update(bool const skip_topology = false) final override{
          /*std::stringstream ss;
          ss << "No need for an update function. \n";
          ss << std::boolalpha <<  "Skip Topology set to " << skip_topology << "\n";
          throw std::runtime_error(ss.str());*/

        }

        float_type e(void) final override;
        float_type g(void) final override;
        float_type h(void) final override;
        float_type o(void) final override;

        void print_E(std::ostream&) const final override;
        void print_E_head(std::ostream&, bool const endline=true) const final override;
        void print_E_short(std::ostream&, bool const endline=true) const final override;
        void to_stream(std::ostream&) const final override;

				std::vector<float_type> charges() const override;
        std::vector<coords::Cartesian_Point> get_g_ext_chg() const override{
          throw std::runtime_error("Function not implemented for psi4 interface.\n");
        }

      private:
        enum class Calc {
          energy,
          gradient,
          optimize
        };
        void write_head(std::ostream&) const;
        void write_input(Calc kind = Calc::energy) const;
        void write_molecule(std::ostream&) const;
				void write_ext_charges(std::ostream & os) const;
        void write_energy_input(std::ostream&) const;
        void write_gradients_input(std::ostream&) const;
        void write_optimize_input(std::ostream&) const;

        void make_call()const;

        std::vector<std::string> parse_specific_position(std::istream& is, std::string const& delim, int space) const;
        std::vector<std::string> get_last_gradients() const;
        coords::Representation_3D get_final_geometry() const;

        coords::float_type parse_energy();
        std::pair<coords::float_type, coords::Representation_3D> parse_gradients();
        std::tuple<coords::float_type, coords::Representation_3D, coords::Representation_3D>
        parse_geometry_and_gradients();

        template<typename StrCont>
        coords::Representation_3D extract_Rep3D(StrCont && lines)const;

        std::string tmp_file_name;

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
