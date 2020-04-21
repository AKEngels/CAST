#pragma once
#include "energy.h"
#include "tinker_parameters.h"
#include "tinker_refine.h"
#include "coords.h"
#define private public //for testing reasons

namespace energy
{
  namespace interfaces
  {
    /**namespace for Amber, Charmm and Oplsaa*/
    namespace aco
    {

      /**class for non-bonded cutoff*/
      class nb_cutoff
      {
      public:
        /**create cutoff object
        @param ic: cutoff-distance
        @param is: switchdist-distance*/
        nb_cutoff(coords::float_type const cutoffDistance, coords::float_type const switchDistance);
        /**test if distance of an atom pair is smaller than cutoff
        @param rr: scalar product of the vector between two atoms
        @param r: reference to distance between two atoms (is calculated during function)
        @param fQ: scaling factor for coulomb energy
        @param fV: scaling factor for vdw energy between switchdist and cutoff*/
        bool factors(coords::float_type const rr, coords::float_type& r, coords::float_type& fQ, coords::float_type& fV);
      private:
        /**cutoff distance*/
        coords::float_type const c;
        /**switchdist*/
        coords::float_type const s;
        /**scalar product of c*c */
        coords::float_type const cc;
        /** 3.0*s*s */
        coords::float_type const ss;
        /** (cc-s*s)*(cc-s*s)*(cc-s*s) */
        coords::float_type const cs;
      };

      /**class for amber, charmm and oplsaa forcefield functions*/
      class aco_ff
        : public interface_base
      {
      public:
        aco_ff(coords::Coordinates* cobj);

        interface_base* clone(coords::Coordinates* coord_object) const;
        interface_base* move(coords::Coordinates* coord_object);
        // initialize using coordinates pointer

        /**update structure (account for topology or rep change, or just non-bonded list)*/
        void update(bool const skip_topology = false);

        /**Energy function*/
        coords::float_type e(void);
        /**Energy+Gradient function*/
        coords::float_type g(void);
        /** Energy+Gradient+Hessian function
        at the moment does nothing because Hessians are not implemented yet*/
        coords::float_type h(void);
        /** Optimization in the intface or interfaced program*/
        coords::float_type o(void);

        /**get charges*/
        std::vector<coords::float_type> charges() const override;
        /**function to get coulomb gradients on external charges*/
        coords::Gradients_3D get_g_ext_chg() const override {
          return grad_ext_charges;
        }

        // Virial Tensor
        std::array<std::array<coords::float_type, 3>, 3> virial();

        // Output functions
        void print_E(std::ostream&) const;
        void print_E_head(std::ostream&, bool const endline = true) const;
        void print_E_short(std::ostream&, bool const endline = true) const;
        void print_G_tinkerlike(std::ostream&, bool const aggregate = false) const final override;
        void to_stream(std::ostream&) const;
        void swap(interface_base&);
        void swap(aco_ff&);


        ::tinker::parameter::parameters const& params() const
        {
          return cparams;
        }
#ifdef GOOGLE_MOCK
      public:
#else
      protected:
#endif
        aco_ff(aco_ff const& rhs, coords::Coordinates* cobj);
        aco_ff(aco_ff&& rhs, coords::Coordinates* cobj);


        enum types {
          BOND = 0, ANGLE, IMPROPER, IMPTORSION, TORSION, MULTIPOLE, CHARGE,
          POLARIZE, VDW, UREY, STRBEND, OPBEND, EXTERNAL_CHARGES, TYPENUM
        };

        /** uncontracted arameters */
        static ::tinker::parameter::parameters tp;
        /** contracted parameters */
        ::tinker::parameter::parameters cparams;
        /** refined parameters */
        ::tinker::refine::refined refined;
        /**main function for calculating all non-bonding interactions:
        - selection of the correct nonbonded function
        - fills part_energy[types::CHARGE], part_energy[types::VDW] and part_grad[types::VDW], part_grad[types::CHARGE]
        */
        template< ::tinker::parameter::radius_types::T RADIUS_TYPE > void   g_nb(void);

        /** Partial Energies for every atom */
        std::array<coords::float_type, TYPENUM>  part_energy;
        /** Partial Gradients for every atom */
        std::array<coords::Representation_3D, TYPENUM> part_grad;
        /** Partial Virials (are they implemented completely???)
        formula: (1/2)* sum_(i,j) r_ij * f_ij   (see http://www.strodel.info/index_files/lecture/MDthermostats_handout.pdf)
        with r_ij = vector from i to j and f_ij = force on i due to j
        unit: kcal/mol*/
        std::array<std::array<std::array<coords::float_type, 3>, 3>, TYPENUM> part_virial;

        void pre(void);
        void post(void);

        /**calculate energy and gradients
        DERIV = 0 -> energy
        DERIV = 1 -> energy and gradients*/
        template<std::size_t DERIV> void calc(void);

        ::tinker::parameter::combi::vdwc const& vdwcpar(std::size_t const ta, std::size_t const tb,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const& pm)
        {
          return pm(cparams.contract_type(ta), cparams.contract_type(tb));
        }

        /** 12-interactions (bonds)*/
        template<std::size_t T_DERV> coords::float_type f_12(void);
        /** 13-interactions (angles) */
        template<std::size_t T_DERV> coords::float_type f_13_a(void);
        /** 13-interactions (urey-bradley) */
        template<std::size_t T_DERV> coords::float_type f_13_u(void);

        /**14-interactions (torsions)*/
        template<std::size_t T_DERV> coords::float_type f_14(void);
        /**imptorsions interactions*/
        template<std::size_t T_DERV> coords::float_type f_it(void);
        /**impropers interactions*/
        template<std::size_t T_DERV> coords::float_type f_imp(void);
        /**calculates interaction to external charges*/
        void calc_ext_charges_interaction(size_t);

        /**gradients on external charges*/
        coords::Gradients_3D grad_ext_charges;

        /**charge energy */
        static coords::float_type eQ(coords::float_type const C, coords::float_type const r);
        /** charge gradients */
        static coords::float_type gQ(coords::float_type const C, coords::float_type const r, coords::float_type& dQ);
        /** charge gradients fep version */
        static coords::float_type gQ_fep(coords::float_type const C, coords::float_type const r, coords::float_type const c_out, coords::float_type& dQ);
        /** vdw energy */
        template< ::tinker::parameter::radius_types::T T_RADIUSTYPE >
        static coords::float_type eV(coords::float_type const E, coords::float_type const R, coords::float_type const r);
        /** vdw gradients */
        template< ::tinker::parameter::radius_types::T T_RADIUSTYPE >
        static coords::float_type gV(coords::float_type const E, coords::float_type const R, coords::float_type const r, coords::float_type& dV);
        /** vdw gradients FEP version*/
        template< ::tinker::parameter::radius_types::T T_RADIUSTYPE >
        static coords::float_type gV_fep(coords::float_type const E, coords::float_type const R,
          coords::float_type const r, coords::float_type const factor, coords::float_type& dV);

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
        static void g_QV(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r,
          coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

        /** charge+vdw gradients fep version (no cutoff, no periodics) */
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
        static void g_QV_fep(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r,
          coords::float_type const c_out, coords::float_type const v_out,
          coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

        /** charge+vdw gradients (cutoff, no fep, no periodics)
        @param C: product of charges of the two atoms
        @param E: epsilon value for lenard-jones
        @param R: equilibrium distance
        @param r: current distance
        @param fQ: scaling factor for coulomb interaction
        @param fV: scaling factor for vdw interaction
        @param e_c: reference to coulomb energy
        @param e_v: reference to vdw energy
        @param dE_c: reference to coulomb gradient
        @param dE_v: reference to vdW gradient*/
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
        static void g_QV_cutoff(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r,
          coords::float_type const fQ, coords::float_type const fV, coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);
        /** charge+vdw gradients fep version (cutoff, no periodics) */
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
        static void g_QV_fep_cutoff(coords::float_type const C, coords::float_type const E, coords::float_type const R, coords::float_type const r,
          coords::float_type const c_out, coords::float_type const v_out, coords::float_type const fQ, coords::float_type const fV,
          coords::float_type& e_c, coords::float_type& e_v, coords::float_type& dE_c, coords::float_type& dE_v);

        /** gradient function for non-bonded pairs */
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
        void g_nb_QV_pairs(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
          std::vector< ::tinker::refine::types::nbpair> const& pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

        /** gradient function for non-bonded pairs with cutoff (with and without periodics)*/
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool PERIODIC>
        void g_nb_QV_pairs_cutoff(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
          std::vector< ::tinker::refine::types::nbpair> const& pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

        /** FEP gradient function for non-bonded pairs
        depends on the fact if the current atom is appearing or disappearing
        @param T_RADIUS_TYPE: r_min or sigma
        @param PERIODIC: periodic boundaries true/false
        @param IS_OUT: true if disappearing atom, false if appearing atom
        */
        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool PERIODIC, bool IS_OUT>
        void g_nb_QV_pairs_fep_io(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
          std::vector< ::tinker::refine::types::nbpair> const& pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

        // in the following there are some functions for non-bonded pairs that are called by the above three functions
        // depending whether the charges are taken from the forcefield (_paramCharges) or from atom charges in config vector (_singleCharges)
        // parameters are the same as explained above

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
        void g_nb_QV_pairs_paramCharges(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
          std::vector< ::tinker::refine::types::nbpair> const& pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
        void g_nb_QV_pairs_singleCharges(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
          std::vector< ::tinker::refine::types::nbpair> const& pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool PERIODIC>
        void g_nb_QV_pairs_cutoff_paramCharges(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
          std::vector< ::tinker::refine::types::nbpair> const& pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool PERIODIC>
        void g_nb_QV_pairs_cutoff_singleCharges(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
          std::vector< ::tinker::refine::types::nbpair> const& pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool PERIODIC, bool IS_OUT>
        void g_nb_QV_pairs_fep_io_paramCharges(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
          std::vector< ::tinker::refine::types::nbpair> const& pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);

        template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool PERIODIC, bool IS_OUT>
        void g_nb_QV_pairs_fep_io_singleCharges(coords::float_type& e_nb, coords::Representation_3D& grad_vdw, coords::Representation_3D& grad_coulomb,
          std::vector< ::tinker::refine::types::nbpair> const& pairs,
          scon::matrix< ::tinker::parameter::combi::vdwc, true> const& parameters);
      };

      /**???*/
      void restrainInternals(coords::Coordinates const& coords_in, ::tinker::refine::refined& refined);
    }
  }
}
