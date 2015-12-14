#pragma once
#include "energy.h"
#include "tinker_parameters.h"
#include "tinker_refine.h"
#include "coords.h"

namespace energy
{
  namespace interfaces
  {
    namespace aco
    {
      

      class nb_cutoff
      {
      public:
        nb_cutoff (coords::float_type const ic, coords::float_type const is);
        inline bool factors (coords::float_type const rr, coords::float_type & r, coords::float_type & fQ, coords::float_type & fV);
      private:
        coords::float_type const c, s, cc, ss, cs;
        nb_cutoff& operator=(nb_cutoff const&);
      };

      class aco_ff 
        : public interface_base
      {
      public:
        aco_ff (coords::Coordinates *cobj);

        interface_base * clone (coords::Coordinates * coord_object) const;
        interface_base * move (coords::Coordinates * coord_object);
        // initialize using coordinates pointer

        // update structure (account for topology or rep change)
        void update (bool const skip_topology = false);

        // Energy function
        coords::float_type e (void);
        // Energy+Gradient function
        coords::float_type g (void);
        // Energy+Gradient+Hessian function
        coords::float_type h (void);
        // Optimization in the intface or interfaced program
        coords::float_type o (void);

        // Virial Tensor

        std::array<std::array<coords::float_type, 3>, 3> virial();

        // Output functions
        void print_E (std::ostream&) const;
        void print_E_head (std::ostream&, bool const endline = true) const;
        void print_E_short (std::ostream&, bool const endline = true) const;
        void print_G_tinkerlike (std::ostream&, bool const aggregate = false) const;
        void to_stream (std::ostream&) const;
        void swap (interface_base&);
        void swap (aco_ff&);

      private:
        
        aco_ff (aco_ff const & rhs, coords::Coordinates *cobj);
        aco_ff (aco_ff && rhs, coords::Coordinates *cobj);


        enum types { BOND = 0, ANGLE, IMPROPER, IMPTORSION, TORSION, MULTIPOLE, CHARGE, 
          POLARIZE, VDW, UREY, STRBEND, OPBEND, VDWC, SOLVATE, TYPENUM};

        // Parameters
        static ::tinker::parameter::parameters tp;
        // Partial Gradients / Energies / Virials
        std::array<coords::Representation_3D, TYPENUM> part_grad;
        std::array<coords::float_type, TYPENUM>  part_energy;
        std::array<std::array<std::array<coords::float_type, 3>, 3>, TYPENUM> part_virial;
        
        ::tinker::parameter::parameters cparams;
        ::tinker::refine::refined refined;

        // Functions for periodic boundary conditions
        void	boundary (coords::float_type&, coords::float_type&, coords::float_type&) const;
        void	boxjump (void);
        inline std::ptrdiff_t sign (coords::float_type const) const;
 
        // Gradient functions (energy functions left since single point is 
        coords::float_type g_it       (void);
        coords::float_type g_imp      (void);
        // selection of the correct nonbonded function
        template< ::tinker::parameter::radius_types::T RADIUS_TYPE > void   g_nb (void);
        // nonbonded function including PME
        template< ::tinker::parameter::radius_types::T RADIUS_TYPE > void   g_nb_pme(void);
        template< ::tinker::parameter::radius_types::T RADIUS_TYPE > void   g_pme_QV_fep(void);
        template< ::tinker::parameter::radius_types::T RADIUS_TYPE > void   g_pme_QV_(void);

        void pre    (void);
        void post   (void);

        template<std::size_t DERIV> void calc (void);

        ::tinker::parameter::combi::vdwc const & vdwcpar (std::size_t const ta, std::size_t const tb, 
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const &pm)
        {
          return pm(cparams.contract_type(ta), cparams.contract_type(tb));
        }

        // 12-interactions 
        template<std::size_t T_DERV> coords::float_type f_12 (void);
        // 13-angles
        template<std::size_t T_DERV> coords::float_type f_13_a (void);
        template<std::size_t T_DERV> coords::float_type f_13_u (void);

        template<std::size_t T_DERV> coords::float_type f_14 (void);
        template<std::size_t T_DERV> coords::float_type f_it (void);
        template<std::size_t T_DERV> coords::float_type f_imp (void);

        // implicit solvation energy
        template<std::size_t T_DERV> coords::float_type solv();

        // charge energy
        inline coords::float_type eQ (coords::float_type const C, coords::float_type const r) const;
        // charge gradients
        inline coords::float_type gQ (coords::float_type const C, coords::float_type const r, coords::float_type &dQ) const;
        // charge gradients fep version
        inline coords::float_type gQ_fep_cut(coords::float_type const C, coords::float_type const ri, coords::float_type const c_out, coords::float_type const dc_out, coords::float_type & Q2, coords::float_type &dQ, coords::float_type & fQ);
        inline coords::float_type gQ_fep(coords::float_type const C, coords::float_type const r, coords::float_type const c_out, coords::float_type &dQ) const;
        // vdw energy
        template< ::tinker::parameter::radius_types::T T_RADIUSTYPE > inline coords::float_type eV 
          (coords::float_type const E, coords::float_type const R, coords::float_type const r) const;
        // vdw gradients
        template< ::tinker::parameter::radius_types::T T_RADIUSTYPE > inline coords::float_type gV 
          (coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type &dV) const;
        template< ::tinker::parameter::radius_types::T T_RADIUSTYPE > inline coords::float_type gV_fep 
          (coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type const factor, coords::float_type &dV) const;
        template< ::tinker::parameter::radius_types::T T_RADIUSTYPE > inline coords::float_type gV_fep_cut
          (coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type const factor, coords::float_type const factor2, coords::float_type &dV, coords::float_type &alche2, coords::float_type &fV, coords::float_type &fV2) const;
        // charge+vdw (no cutoff, no fep, no periodics)

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE> 
        inline void e_QV  (coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r, 
          coords::float_type &e_c, coords::float_type &e_v) const;
        // charge+vdw gradients
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE> 
        inline void g_QV  (coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r, 
          coords::float_type &e_c, coords::float_type &e_v, coords::float_type &dE) const;

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
        inline void g_QV_fep  (coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r, 
          coords::float_type const c_out, coords::float_type const v_out,
          coords::float_type &e_c, coords::float_type &e_v, coords::float_type &dE) const;

        // charge+vdw (cutoff, no fep, no periodics)

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE> 
        inline void e_QV_cutoff  (coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r, 
          coords::float_type const fQ, coords::float_type const fV,coords::float_type &e_c, coords::float_type &e_v) const;
        // charge+vdw gradients
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE> 
        inline void g_QV_cutoff  (coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r, 
          coords::float_type const fQ, coords::float_type const fV, coords::float_type &e_c, coords::float_type &e_v, coords::float_type &dE) const;
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE> 
        inline void g_QV_fep_cutoff  (coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r, 
          coords::float_type const c_out, coords::float_type const v_out, coords::float_type const fQ, coords::float_type const fV, coords::float_type &e_c, coords::float_type &e_v, coords::float_type &dE) const;

        // pairs functions
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE> 
        void g_nb_QV_pairs (coords::float_type &e_nb, coords::Representation_3D &grad_vector, 
          std::vector< ::tinker::refine::types::nbpair> const & pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const & parameters);

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool PERIODIC> 
        void g_nb_QV_pairs_cutoff (coords::float_type &e_nb, coords::Representation_3D &grad_vector, 
          std::vector< ::tinker::refine::types::nbpair> const & pairs, 
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const & parameters);
        
        //PME functions
        template < ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
        void pme_direct(double &, coords::Representation_3D &grad_vector,
          std::vector< ::tinker::refine::types::nbpair> const & pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const & parameters, double &);

        template < ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool ALCH_OUT>
        void pme_direct_fep(double &, coords::Representation_3D &grad_vector,
          std::vector< ::tinker::refine::types::nbpair> const & pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const & parameters, double &);


        void pme_direct_scaled(double &, coords::Representation_3D &);
        void pme_direct_scaled_fep(double &, coords::Representation_3D &);
        void pme_reci(double &);
        void pme_reciall(double &);
        void pme_reci_novir(double &);
        void pme_reci_grad(coords::Representation_3D &);
        void pme_reci_fepingrad(coords::Representation_3D &);
        void pme_reci_fepoutgrad(coords::Representation_3D &);
        void pme_reci_fepallgrad(coords::Representation_3D &);
        void pme_chargedcell();
        void pme_correction(double &);
        void pme_correction_fep(double &);
        void bsplinegeneration();
        void bsplinefepin();
        void bsplinefepout();
        void bsplinefepall();
        void fepinsplines();
        void fepoutsplines();
        void fepallsplines();
        void singlespline(double, int, int);
        void generatetable();
        void generateintable();
        void generateouttable();
        void generatealltable();
        void generatesinglesite(int &, std::vector<int> &, const int, const int, const int);
        void generateinsite(int &, std::vector<int> &, const int, const int, const int);
        void generateoutsite(int &, std::vector<int> &, const int, const int, const int);
        void generateallsite(int &, std::vector<int> &, const int, const int, const int);
        void setchargestogrid();
        void setchargestofepingrid();
        void setchargestofepoutgrid();
        void setchargestofepindgrid();
        void setchargestofepoutdgrid();
        void setchargestoallgrid();
        void sitecorrection(int &, int &, int &, int &, int &, int &, int &);

        // FEP pair functions
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE> 
        void g_nb_QV_pairs_fep_switch_periodic (coords::float_type &e_nb, coords::Representation_3D &grad_vector, 
          std::vector< ::tinker::refine::types::nbpair> const & pairs, 
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const & parameters);

        //template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool PERIODIC> 
        //void g_nb_QV_pairs_fep (coords::float_type &e_nb, coords::Representation_3D &grad_vector, 
        //  std::vector< ::tinker::refine::types::nbpair> const & pairs, 
        //  scon::tsmatrix< ::tinker::parameter::combi::vdwc> const & parameters);

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool PERIODIC, bool IS_OUT> 
        void g_nb_QV_pairs_fep_io(coords::float_type &e_nb, coords::Representation_3D &grad_vector,
          std::vector< ::tinker::refine::types::nbpair> const & pairs, 
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const & parameters);

      };
    }
  }
}