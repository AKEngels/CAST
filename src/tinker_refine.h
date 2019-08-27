#pragma once
#include <array>
#include <vector>
#include <ostream>
#include "tinker_parameters.h"
#include "Scon/scon_matrix.h"
#include <boost/optional.hpp>

namespace tinker
{
  namespace refine
  {
    namespace types
    {
      struct binary_quadratic
      {
        double force, ideal;
        std::array<std::size_t, 2> atoms;
        binary_quadratic(void) : force(), ideal(), atoms() {}
      };

			std::ostream& operator<< (std::ostream& stream, binary_quadratic const& bq);

			struct ternary_quadratic
			{
				double force, ideal;
				std::array<std::size_t, 3> atoms;
				ternary_quadratic(void) : force(), ideal(), atoms() {}
			};

			std::ostream& operator<< (std::ostream& stream, ternary_quadratic const& bq);

			struct torsion
			{
				std::array<std::size_t, 4> atoms;
				::tinker::parameter::torsion p;
				torsion(void) : atoms(), p() {}
			};

			std::ostream& operator<< (std::ostream& stream, torsion const& bq);

			struct improper
			{
				std::size_t center, ligand[2], twist;
				::tinker::parameter::improper p;
				improper(std::size_t c, std::size_t t, std::size_t l1, std::size_t l2, ::tinker::parameter::improper par, double sym = 1.0)
					: center(c), ligand(), twist(t), p(par)
				{
					ligand[0] = l1;
					ligand[1] = l2;
					for (auto& f : p.force) f /= sym;
				}
				improper(void)
					: center(), ligand(), twist(), p() {}
			};

			std::ostream& operator<< (std::ostream& stream, improper const& bq);

			struct imptor
			{
				std::size_t center, ligand[2], twist;
				::tinker::parameter::imptor p;
				imptor(std::size_t c, std::size_t t, std::size_t l1, std::size_t l2, ::tinker::parameter::imptor par, double sym = 1.0)
					: center(c), ligand(), twist(t), p(par)
				{
					ligand[0] = l1;
					ligand[1] = l2;
					for (auto& f : p.force) f /= sym;
				}
				imptor(void)
					: center(), ligand(), twist(), p() {}
			};

			std::ostream& operator<< (std::ostream& stream, imptor const& bq);

			struct multipole
			{
				//std::array<std::size_t, 4> atoms;
				scon::c3<std::size_t> axes;
				std::size_t center, npole;
				::tinker::parameter::multipole p_rot;
				::tinker::parameter::multipole const* p_nonrot;
				//multipole (T axistype, 
				multipole(void)
					: axes(), center(), npole(), p_rot(), p_nonrot(nullptr) {}
				multipole(::tinker::parameter::multipole const& mp,
					std::size_t const center_index, scon::c3<std::size_t> const axis_indices)
					: axes(axis_indices), center(center_index), npole(), p_rot(), p_nonrot(&mp) {}
				scon::m3<double> rotation_matrix(coords::Representation_3D const&) const;
				void rotate_to_frame(coords::Representation_3D const&);
			};

			std::ostream& operator<< (std::ostream& stream, multipole const& bq);

			struct opbend
			{
				double force;
				std::array<std::size_t, 4> atoms;
				opbend(void) : force(), atoms() {}
			};

			std::ostream& operator<< (std::ostream& stream, opbend const& bq);

			struct polarize
			{
				double force, ff, pdamp;
				std::array<std::size_t, 3u> atoms;
				std::size_t center;
				scon::c3<double> thole, polarity;
				scon::c3<std::size_t> p_group;
				polarize(void) : force(), ff(), pdamp(), atoms(), center(), thole(), polarity(), p_group() {}
			};

			std::ostream& operator<< (std::ostream& stream, polarize const& bq);

			struct strbend
			{
				double force, ff;
				std::array<std::size_t, 3u> atoms;
				strbend(void) : force(), ff(), atoms() {}
			};

			std::ostream& operator<< (std::ostream& stream, strbend const& bq);

			// pair of nonbonded interactions, only having indices
			struct nbpair
			{
				std::size_t a, b;
				nbpair() : a(), b() {}
				nbpair(std::size_t ia, std::size_t ib) : a(ia), b(ib) {}
				void swap_them() { std::swap(a, b); }
			};

			std::ostream& operator<< (std::ostream& stream, nbpair const& bq);

      struct nbpm
      {
        typedef std::vector<types::nbpair> vector_pairs;
        scon::matrix<vector_pairs, true> pair_matrix;
        std::size_t param_matrix_id;
        nbpm(void) : pair_matrix(), param_matrix_id() {}
        nbpm(std::size_t matrix_size, std::size_t param_id)
          : pair_matrix(matrix_size), param_matrix_id(param_id) {}
      };

			std::ostream& operator<< (std::ostream& stream, nbpm const& bq);

    }
    // Typedefs
    typedef std::vector<types::binary_quadratic>   vector_biquad;
    typedef std::vector<types::ternary_quadratic>  vector_triquad;
    typedef std::vector<types::torsion>            vector_tors;
    typedef std::vector<types::imptor>             vector_imptors;
    typedef std::vector<types::improper>           vector_improper;
    typedef std::vector<types::multipole>          vector_multipole;
    typedef std::vector<types::opbend>             vector_opbend;
    typedef std::vector<types::polarize>           vector_polarize;
    typedef std::vector<types::strbend>            vector_strbend;
    typedef std::vector<types::nbpair>             vector_pairs_1d;
    typedef scon::matrix<vector_pairs_1d, true>     pairs_matrix;
    typedef std::vector<std::size_t>               vector_size_1d;
    typedef std::vector<vector_size_1d>            vector_size_2d;

    class refined
    {
    public:
      enum rel { R11, R12, R13, R14, R15, R1N };

      refined(){};
      refined(coords::Coordinates const & cobj, tinker::parameter::parameters const & pobj);
      

      void clear(void);
      void swap_data(refined&);



      vector_biquad const &   bonds(void) const { return m_bonds; }
      vector_triquad const &  angles(void) const { return m_angles; }
      vector_improper const & impropers(void) const { return m_impropers; }
      vector_imptors const &  imptors(void) const { return m_imptors; }
      vector_opbend const & opbends(void) const { return m_opbends; }
      vector_strbend const & strbends(void) const { return m_strbends; }
      vector_tors const &     torsions(void) const { return m_torsions; }
      vector_biquad const &   ureys(void) const { return m_ureys; }
      std::vector<types::nbpm> const & pair_matrices(void) const { return m_pair_matrices; };
      vector_size_1d const &   remove_relations(std::size_t const index) { return m_removes[index]; }
      std::vector<vector_multipole> const & multipole_vecs(void) const { return m_multipole_vec; };
      std::vector<vector_polarize> const & polarize_vecs(void) const { return m_polarize_vec; };

      vector_biquad  &                    set_bonds(void) { return m_bonds; }
      vector_triquad  &                   set_angles(void) { return m_angles; }
      vector_improper  &                  set_impropers(void) { return m_impropers; }
      vector_imptors  &                   set_imptors(void) { return m_imptors; }
      vector_opbend  &                    set_opbends(void) { return m_opbends; }
      vector_strbend  &                   set_strbends(void) { return m_strbends; }
      vector_tors  &                      set_torsions(void) { return m_torsions; }
      vector_biquad  &                    set_ureys(void) { return m_ureys; }
      std::vector<types::nbpm>  &         set_pair_matrices(void)  { return m_pair_matrices; };
      vector_size_1d  &                   set_remove_relations(std::size_t const index) { return m_removes[index]; }
      std::vector<vector_multipole>  &    set_multipole_vecs(void) { return m_multipole_vec; };
      std::vector<vector_polarize>  &     set_polarize_vecs(void) { return m_polarize_vec; };


      scon::matrix<parameter::combi::vdwc, true> const & vdwcm(std::size_t index) const { return m_vdwc_matrices[index]; }
      std::array<scon::matrix<parameter::combi::vdwc, true>, 6u> const & vdwcm(void) const { return m_vdwc_matrices; }

			std::size_t const& type(std::size_t const index) const { return m_red_types[index]; }

			std::size_t ia_count(void) const
			{
				std::size_t ia(0u);
				for (auto const& nbpmatrix : m_pair_matrices)
				{
					for (auto const& pairvec : nbpmatrix.pair_matrix)
					{
						ia += pairvec.size();
					}
				}
				return ia;
			}

      void refine_nb(coords::Coordinates const & cobj);

    private:
      //
      ::tinker::parameter::parameters m_cparams;
      //
      bool has_tighter_relation(std::size_t const atom, std::size_t const related, std::size_t const relation_to_check);
      void remove_loose_relations(std::size_t const atom, std::size_t const related, std::size_t const relation_to_check);
      void add_relation(tinker::parameter::parameters const & pobj, std::size_t const atom, std::size_t const related, std::size_t const relation);
      template<rel RELATION> void build_pairs_direct(coords::Coordinates const & coords);
      template<rel RELATION> bool add_pair(coords::Coordinates const & coords,std::size_t const row, std::size_t const col, std::array<std::size_t, 5u> const & to_matrix_id);

      //void refine_mp(coords::Coordinates const & coords, tinker::parameter::parameters const & params);
      //void refine_pol(void);



      // Standard Types
      vector_triquad                                                m_angles;
      vector_biquad                                                 m_bonds;
      vector_improper                                               m_impropers;
      vector_imptors                                                m_imptors;
      vector_tors                                                   m_torsions;
      vector_biquad                                                 m_ureys;
      vector_opbend                                                 m_opbends;
      vector_strbend                                                m_strbends;
      std::vector<types::nbpm>                                      m_pair_matrices;
      std::vector<vector_multipole>                                 m_multipole_vec;
      std::vector<vector_polarize>                                  m_polarize_vec;
      vector_pairs_1d                                               m_1d_nbpairs;
      std::vector<vector_pairs_1d>                                  m_1d_nbpairs_vec;
      // Refined relations 
      // 5u -> (-11- -12- -13- -14- -15-)
      std::array<vector_size_2d, 5u>                                m_relations;
      // removed relations
      vector_size_2d                                                m_removes;
      vector_size_1d                                                m_red_types;
      // Refined vdw matrices
      // 6u -> (-11- -12- -13- -14- -15- -1n-)
      std::array<scon::matrix<parameter::combi::vdwc, true>, 6u>    m_vdwc_matrices;
    };

    //Re-Implemented as Free Function, much better code design! Dustin 07/2019
    //tinker::refine::refined refine(coords::Coordinates const & coords, tinker::parameter::parameters const &);
    tinker::refine::vector_imptors find_imptors(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t center, std::size_t a, std::size_t b, std::size_t c);
    tinker::refine::vector_improper find_improper(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t center, std::size_t a, std::size_t b, std::size_t c);
    boost::optional<types::binary_quadratic> find_bond(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b);
    boost::optional<types::ternary_quadratic> find_angle(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b, std::size_t c);

    boost::optional<tinker::refine::types::torsion> find_torsion(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b, std::size_t c, std::size_t d);

    boost::optional<tinker::refine::types::binary_quadratic> find_urey(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b, std::size_t c);
    boost::optional<tinker::refine::types::opbend> find_opbend(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b, std::size_t c, std::size_t d);
    boost::optional<tinker::refine::types::strbend> find_strbend(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b, std::size_t c);
    
    tinker::refine::vector_polarize  refine_pol(coords::Coordinates const & coords, tinker::parameter::parameters const & params);
    tinker::refine::vector_multipole refine_mp(coords::Coordinates const & coords, tinker::parameter::parameters const & params);

    std::pair<tinker::refine::vector_imptors, tinker::refine::vector_improper> refine_imptor_improper_of_atom(coords::Coordinates const & cobj, tinker::parameter::parameters const & pobj, std::size_t atom);
    
    // Output of refined class
    std::ostream& operator<< (std::ostream &stream, refined const &ref);

    void purge_nb_at_same_molecule(coords::Coordinates const & coords, tinker::refine::refined & params);

  }
}
