#pragma once
#include "energy.h"
#include "tinker_parameters.h"
#include "tinker_refine.h"
#include "coords.h"
#include "interpolation.h"

namespace energy
{
  namespace interfaces
  {
    namespace amoeba
    {
      
	 
      class nb_cutoff
      {
      public:
        nb_cutoff (double const ic, double const is);
        inline bool factors (double const rr, double & r, double & fQ, double & fV);
      private:
        double const c, s, cc, ss, cs;
       
      };

      class amoeba_ff
        : public interface_base
      {
      public:
        amoeba_ff (coords::Coordinates *cobj);
	/*	amoeba_ff(coords::Coordinates const &ref);*/
		size_t alloc_glob;
		interface_base * clone(coords::Coordinates * coord_object) const;
        interface_base * move (coords::Coordinates * coord_object);
		
        // initialize using coordinates pointer

        // update structure (account for topology or rep change)
        void update (bool const skip_topology = false);

        // Energy function
        double e (void);
        // Energy+Gradient function
        double g (void);
        // Energy+Gradient+Hessian function
        double h (void);
        // Optimization in the intface or interfaced program
        double o (void);

        /**overwritten function, should not be called*/
        std::vector<coords::float_type> charges() const override
        {
          throw std::runtime_error("TODO: Implement charge getter for AMOEBA.\n");
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
        void swap (interface_base&);
        void swap (amoeba_ff&);



		//!Definition of SPACKMAN-Variables and Functions
		//***********************************************
		//***********************************************
		struct atom
		{
			char symbol[3];
			ptrdiff_t index, group;
			size_t nBonds, atomNumber;
			double mass, charge;

		};
		struct bond
		{
			ptrdiff_t index[2];
			double force, relativeValue;
		};
		struct n
		{
			ptrdiff_t temp1[19];
		};
		struct z
		{
			float temp2[19];
		};
		struct c
		{
			float temp3[19];
		};
		struct o
		{
			ptrdiff_t temp4[3];
		};
		struct k
		{
			double temp5;
		};
		struct hh
		{
			double temp6;
		};
		struct nn
		{
			double temp7;
		};
		struct oo
		{
			double temp8;
		};


		ptrdiff_t MON;
		size_t kp;
		float nonbond_cutoff;
		double cutoff_spackman;
		double energy_short;
		double energy_short_analytical;
		double energy_total;
		double fff;
		bool SPACKrefine;
		std::vector<double> dex11, dex22, dex33, exa11, exa22, exa33;
		std::vector<double> dex44, dex55, dex66, exa44, exa55, exa66;
		std::vector<double> dex77, dex88, dex99, dex1010, exa77, exa88, exa99, exa1010;
		std::vector<double> dex0, xvec1, evec1, eveca1;

		std::vector<atom>     atoms;
		std::vector<bond>     bonds;
		std::vector<n>        nij1;
		std::vector<z>        zeta;
		std::vector <std::vector <float> > zetax, cs;
		std::vector <double> kappan;
		std::vector <std::vector <ptrdiff_t> > nijx, occ;


		std::vector <std::vector<float> > zetan;
		size_t nbas;
		size_t nat;
		size_t nref;

		std::vector <std::vector<double> > coef;
		std::vector <std::vector<std::vector<double> > > pij;
		std::vector <double> sscat;
		double aLimes, bLimes, smax;
		size_t fak;
		std::vector <double> kappa;
		std::vector <size_t>  atomic;
		std::vector <std::vector<double> > fscat;


		std::vector<size_t> mol;
		struct spack_list {

			/* size_t atom1;
			size_t atom2;*/
			size_t atom[2];


		};
		std::vector <spack_list> vec_spack;
		size_t count;


		struct ccoord{
			double x;
			double y;
			double z;
		};
		std::vector < ccoord > ccenter;

		void Spackman_mol(void);
		void parameters();

		void SpackmanGrad_3(void);
		void Spackman_GRAD(void);

		void Spackman_vec(void);
		void Spackman_list(void);
		void Spackman_list_analytical1(void);
		double fvalue_f(ptrdiff_t, ptrdiff_t);
		double fvalue_f_2(ptrdiff_t, double, ptrdiff_t);
		double fff_f1(size_t, double, double);
		size_t fak_iter1(size_t);

		void Spackman1(void);
		double Spackman_Energy(void);
		double Spackman_energy_analytical(void);

      private:
        


        amoeba_ff (amoeba_ff const & rhs, coords::Coordinates *cobj);
        amoeba_ff (amoeba_ff && rhs, coords::Coordinates *cobj);


        enum types {TYPENUM=15, BOND=0, ANGLE, IMPROPER, IMPTORSION, TORSION, MULTIPOLE, CHARGE, POLARIZE, VDW, UREY, STRBEND, OPBEND, VDWC, INDUCE, SHORTRANGE };

        std::array<coords::Representation_3D, TYPENUM> part_grad;
        std::array<double, TYPENUM>  part_energy;
        std::array<std::array<std::array<double, 3>, 3>, TYPENUM> part_virial;

		

        static ::tinker::parameter::parameters tp;
        ::tinker::parameter::parameters cparams;
        ::tinker::refine::refined refined;
		::tinker::parameter::multipole multi;
		::tinker::refine::types::multipole types;
		::tinker::refine::types::binary_quadratic binary;
		
		
        // Functions for periodic boundary conditions
        void	boundary (double&, double&, double&) const;
        inline ptrdiff_t sign (double const) const;

        // Gradient functions (energy functions left since single point is 
        double g_it       (void);
        double g_imp      (void);
        // selection of the correct nonbonded function
        template< ::tinker::parameter::radius_types::T RADIUS_TYPE > void   g_nb (void);

        void pre    (void);
        void post   (void);
		

		template<size_t DERIV> void calc(void);

		::tinker::parameter::combi::vdwc const & vdwcpar(size_t const ta, size_t const tb,
			scon::matrix< ::tinker::parameter::combi::vdwc, true> const &pm)
		{
			return pm(cparams.contract_type(ta), cparams.contract_type(tb));
		}

        // 12-interactions 
        template<size_t T_DERV> double f_12 (void);
        // 13-angles
        template<size_t T_DERV> double f_13_a (void);
        template<size_t T_DERV> double f_13_u (void);
		template<size_t T_DERV> double f_strbnd (void);

        template<size_t T_DERV> double f_14 (void);
        template<size_t T_DERV> double f_it (void);
        template<size_t T_DERV> double f_imp (void);
		template<size_t T_DERV> double f_oop(void);


        // multipole gradients and energies
        void e_perm (void);



        void e_ind (void);
		inline size_t multipole_sites(void);
		void rot_matrix(coords::Representation_3D const &pos);
		//std::std::vector <double> ci, dx, dy, dz, qxx, qyy, qzz, qxy, qxz, qyx, qyz, qzx, qzy;
		std::vector <std::vector <double> > quadro, rp, dipole;
		std::vector <double> m_charges;
		std::vector <std::vector <double> > dem, dep;
		double em, ep;
	
		//multipole scalingfactors
		double pscale2, pscale3, pscale4, pscale5;
		double mscale2, mscale3, mscale4, mscale5;
		double dscale1, dscale2, dscale3, dscale4;
		double uscale1, uscale2, uscale3, uscale4;
		double p4scale;

		//Polarization GROUP params

		std::vector <size_t> mask, list, keep;


        

       //multipole rotation into coordinate framework

		std::vector < std::vector <double> > uind, uinp;
		std::vector < size_t > ipole,xaxis,zaxis,yaxis,axistype;
		std::vector < std::vector < size_t >  > plrgrp;
		std::vector <double> pdamp, thole, pscale, dscale, polarity, mscale, uscale;
		
		//pair lists for multipoles


		void refine_pair_lists(void);
		void refine_vdw_h_bonds(std::vector< ::tinker::refine::types::nbpair> const & pairs,
			scon::matrix< ::tinker::parameter::combi::vdwc, true> const & parameters);
	
		coords::Representation_3D vdwnew;
		std::vector< std::vector<size_t> > i12, i13, i14, i15, ip11, ip13, ip12, ip14;
		std::vector <size_t> n13, n14, n15, np11, np12, np13, np14;





		// charge energy
		inline double eQ(double const C, double const r) const;
		// charge gradients
		inline double gQ(double const C, double const r, double &dQ) const;
		// charge gradients fep version
		inline double gQ_fep(double const C, double const r, double const c_out, double &dQ) const;
		// vdw energy
		template< ::tinker::parameter::radius_types::T T_RADIUSTYPE > inline double eV
			(double const E, double const R, double const r) const;
		// vdw gradients
		template< ::tinker::parameter::radius_types::T T_RADIUSTYPE > inline double gV
			(double const E, double const R, double const r, double &dV) const;
		template< ::tinker::parameter::radius_types::T T_RADIUSTYPE > inline double gV_fep
			(double const E, double const R, double const r, double const factor, double &dV) const;

		// charge+vdw (no cutoff, no fep, no periodics)

		template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
		inline void e_QV(double const C, double const E, double const R, double const r,
			double &e_c, double &e_v) const;
		// charge+vdw gradients
		template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
		inline void g_QV( double const E, double const R, double const r,
			 double &e_v, double &dE) const;

		template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
		inline void g_QV_fep(double const E, double const R, double const r,
		 double const v_out, double &e_v, double &dE) const;

		// charge+vdw (cutoff, no fep, no periodics)

		template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
		inline void e_QV_cutoff(double const C, double const E, double const R, double const r,
			double const fQ, double const fV, double &e_c, double &e_v) const;
		// charge+vdw gradients
		template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
		inline void g_QV_cutoff( double const E, double const R, double const r,
			 double const fV, double &e_v, double &dE) const;
		template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
		inline void g_QV_fep_cutoff( double const E, double const R, double const r,
			double const v_out, double const fV, double &e_v, double &dE) const;

		// pairs functions

		template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
		void g_nb_QV_pairs(coords::float_type &e_nb, coords::Representation_3D &grad_vector,
			std::vector< ::tinker::refine::types::nbpair> const & pairs,
			scon::matrix< ::tinker::parameter::combi::vdwc, true> const & parameters);

		template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE, bool PERIODIC>
		void g_nb_QV_pairs_cutoff(coords::float_type &e_nb, coords::Representation_3D &grad_vector,
			std::vector< ::tinker::refine::types::nbpair> const & pairs,
			scon::matrix< ::tinker::parameter::combi::vdwc, true> const & parameters);

		// FEP pair functions

		template< ::tinker::parameter::radius_types::T T_RADIUS_TYPE>
		void g_nb_QV_pairs_fep_switch_periodic(coords::float_type &e_nb, coords::Representation_3D &grad_vector,
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