#pragma once
#include "coords.h"
#include <vector>
#include "coords_atoms.h"
#include "energy_int_aco.h"
#include "energy_int_mopac.h"
#include "tinker_parameters.h"

namespace energy
{
  namespace interfaces
  {
    namespace qmmm
    {
      class QMMM
        : public interface_base
      {

        static ::tinker::parameter::parameters tp;
        ::tinker::parameter::parameters cparams;

      public:


        QMMM(coords::Coordinates*);
        QMMM(QMMM const&, coords::Coordinates*);
        QMMM(QMMM&&, coords::Coordinates*);

        double distance;
        interface_base * clone(coords::Coordinates*) const;
        interface_base * move(coords::Coordinates*);

        void swap(interface_base &);
        void swap(QMMM &);

        // update interface information from coordinates
        void update(bool const skip_topology = false);
        
        // Energy function
        coords::float_type e() override;
        // Energy+Gradient function
        coords::float_type g() override;
        // Energy+Gradient+Hessian function
        coords::float_type h() override;
        // Optimization in the intface or interfaced program
        coords::float_type o() override;

        // Return charges (QM und MM charges?)
        std::vector<coords::float_type> charges() const override;
        void print_E(std::ostream&) const override;
        void print_E_head(std::ostream&, 
          bool const endline = true) const override;
        void print_gnuplot(std::ostream&, 
          bool const endline = true) const;
        void print_E_short(std::ostream&, 
          bool const endline = true) const override;
        void print_G_tinkerlike(std::ostream&, 
          bool const aggregate = false) const override;
        void to_stream(std::ostream&) const;

        void set_distance(double d)
        {
          distance = d;
        }

      private:

        void ww_calc(bool);
        coords::float_type qmmm_calc(bool);

        std::vector<std::size_t> qm_indices;
        std::vector<std::size_t> mm_indices;
  
        std::vector<std::size_t> new_indices_qm;
        std::vector<std::size_t> new_indices_mm;
       
        coords::Coordinates qmc;
        coords::Coordinates mmc;
        
        std::vector<double> qm_charge_vector;
        std::vector<double> mm_charge_vector;

        //coords::float_type c_energy;
        coords::float_type vdw_energy;
        coords::float_type qm_energy;
        coords::float_type mm_energy;

        coords::Gradients_3D c_gradient;
        coords::Gradients_3D vdw_gradient;
               
        
        
        

        
      };
    }
  }
}
