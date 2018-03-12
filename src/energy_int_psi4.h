#ifndef ENERGY_INT_PSI_H
#define ENERGY_INT_PSI_H

#include<vector>
#include<utility>

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

        void swap(interface_base & other) final;
        interface_base * clone(coords::Coordinates* coord_object) const final;
        interface_base * move(coords::Coordinates* coord_object) final;

        void update(bool const skip_topology = false) final;

        float_type e(void) final;
        float_type g(void) final;
        float_type h(void) final;
        float_type o(void) final;

        void print_E(std::ostream&) const final;
        void print_E_head(std::ostream&, bool const endline=true) const final;
        void print_E_short(std::ostream&, bool const endline=true) const final;
        void print_G_tinkerlike(std::ostream&, bool const aggregate = false) const final;
        void to_stream(std::ostream&) const final;

        std::vector<float_type> charges() const override{
          throw std::runtime_error("TODO: Implement charge getter for AMOEBA.\n");
        }
        std::vector<coords::Cartesian_Point> get_g_coul_mm() const override{
          throw std::runtime_error("Function not implemented for psi4 interface.\n");
        }
        std::string get_id() const override{
          throw std::runtime_error("Function not implemented for psi4 interface.\n");
        }
        /**overwritten function, should not be called*/
        coords::Gradients_3D get_link_atom_grad() const override
        {
          throw std::runtime_error("function not implemented\n");
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
        void write_energy_input(std::ostream&) const;
        void write_gradients_input(std::ostream&) const;
        void write_optimize_input(std::ostream&) const;

        void make_call()const;

        std::vector<std::string> parse_specific_position(std::istream& is, std::string const& delim, int space) const;
        std::vector<std::string> get_last_gradients() const;
        coords::Representation_3D get_final_geometry() const;

        coords::float_type parse_energy()const;
        std::pair<coords::float_type, coords::Representation_3D> parse_gradients() const;
        std::tuple<coords::float_type, coords::Representation_3D, coords::Representation_3D>
        parse_geometry_and_gradients() const;

        coords::Representation_3D extract_Rep3D(std::vector<std::string> const& lines)const;

        std::string tmp_file_name;
      };
    }
  }
}

#endif
