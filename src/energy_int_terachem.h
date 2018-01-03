/*DW                            *
    TeraChem MPI Interface
*                            DW*/


#pragma once 

#if defined(USE_MPI)

#include <array> 
#include <vector>
#include <string>
#include <iostream>
#include "mpi_cast.h"
#include "energy.h"
#include "scon_vect.h"


namespace energy 
{

  namespace interfaces
  {

    namespace terachem
    {

      class mpiInterface 
        : public energy::interface_base
      {
        // MPI Connector object
        static cmpi::CMPI mpo;
        // Check whether option init is done
        static bool       option_init_done;

      public:

        mpiInterface (coords::Coordinates *cobj);
        ~mpiInterface (void);

        interface_base * clone (coords::Coordinates * coord_object) const;
        interface_base * move (coords::Coordinates * coord_object);

        void update (bool const skip_topology = false);

        /*
         Energy class functions that need to be overloaded
        */

        virtual void swap (interface_base &other);

        void swap (mpiInterface &other);

        // Energy function
        coords::float_type e (void);
        // Energy+Gradient function
        coords::float_type g (void);
        // Energy+Gradient+Hessian function
        coords::float_type h (void);
        // Optimization in the interface(d program)
        coords::float_type o (void);

        /**overwritten function, should not be called*/
        std::vector<coords::float_type> charges() const override
        {
          throw std::runtime_error("TODO: Implement charge getter for TeraChem.\n");
        }
        /**overwritten function, should not be called*/
        std::vector<coords::Cartesian_Point> get_el_field() const override
        {
          throw std::runtime_error("TODO: Implement electric field.\n");
        }
        /**overwritten function*/
        std::string get_id() const override { return "bullshit"; }

        // Output functions
        void print_E (std::ostream&) const;
        void print_E_head (std::ostream&, bool const endline = true) const;
        void print_E_short (std::ostream&, bool const endline = true) const;
        void print_G_tinkerlike (std::ostream&, bool const aggregate = false) const;
        void to_stream (std::ostream&) const;

      private:

        mpiInterface (mpiInterface const & rhs, coords::Coordinates *cobj);
        mpiInterface (mpiInterface && rhs, coords::Coordinates *cobj);

        // initialize MPI
        static void init_terachem (void);
        // MPI send calculation request
        void   mpi_send_data (const int);
        // MPI recieve calculated data
        void mpi_recv_energy    (void);
        void mpi_recv_gradients (void);
        void mpi_recv_positions (void);

        std::vector<double> qm_population_charges;

      };

    }

  }

}

#endif
