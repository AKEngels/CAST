#pragma once

#include <array>
#include <vector>
#include <string>
#include <stddef.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "scon_angle.h"
#include "scon_matrix.h"
#include "scon_tsqcmatrix.h"
#include "scon_vect.h"
#include "coords.h"

namespace tinker
{

  enum potential_keys { KEYS=13, ANGLE=0, BOND, CHARGE, IMPROPER, IMPTORS, MULTIPOLE, 
                      OPBEND, POLARIZE, STRBEND, TORSION, UREYBRAD, VDW, VDW14 };

  enum index_types { GROUP, TYPE };

  namespace parameter
  {

    struct forcefield_types 
    { 
      enum value_type { OPLSAA, CHARMM22, AMBER, AMOEBA } value; 
      forcefield_types (void) : value(AMBER) {}
      void parse_line (std::string const & line);
    };
    std::ostream& operator<< (std::ostream & stream, forcefield_types const & t);
    struct vdw_types          
    { 
      enum { LENNARD_JONES, BUFFERED_14_7, BUCKINGHAM, MM3_HBOND, GAUSSIAN } value;
      vdw_types (void) : value(LENNARD_JONES) {}
      void parse_line (std::string const & line);
    };
    std::ostream& operator<< (std::ostream & stream, vdw_types const & t);
    struct average_rules          
    { 
      enum { GEOMETRIC, ARITHMETIC, CUBIC_MEAN, HARMONIC, HHG } value;
      average_rules (void) : value(ARITHMETIC) {}
      void parse_line (std::string const & line);
      inline double process (double x) const
      {
        using std::abs;
        switch (value)
        {
        case CUBIC_MEAN:
        case ARITHMETIC:
        case HARMONIC:
        case HHG:
          {
            return x;
          }
        default:
          { // GEOMETRIC
            return abs(x);
          }
        }
      }
      inline double process (double a, double b) const
      {
        using std::sqrt;
        using std::pow;
        using std::abs;
        if (a==0.0&&b==0.0) return 0.0;
        switch (value)
        {
        case CUBIC_MEAN:
          {
						   auto a2 = a*a;
						   auto b2 = b*b;
						   return(a2*a + b2*b) / (a2 + b2);
          }
        case ARITHMETIC:
          {
            return (a+b)/2;
          }
        case HARMONIC:
          {
            return 2.0*(a*b)/(a+b);
          }
        case HHG:
          {
            double const s(sqrt(a)+sqrt(b)), s2(s*s);
            if(abs(s2) > 0.0) return 4.0*a*b/s2;
            else return 0.0;
          }
        default: 
          { // GEOMETRIC
            return sqrt(a*b);
          }
        }
      }
    };
    std::ostream& operator<< (std::ostream & stream, average_rules const & t);
    struct radius_types
    {
      enum T { R_MIN, SIGMA } value;
      radius_types (void) : value(R_MIN) {}
      void parse_line (std::string const & line);
    };
    std::ostream& operator<< (std::ostream & stream, radius_types const & t);
    struct radius_sizes
    {
      enum { RADIUS, DIAMETER } value;
      radius_sizes (void) : value(RADIUS) {}
      void parse_line (std::string const & line);
      inline double process (double x) const
      {
        switch (value)
        {
        case DIAMETER:
          return x;
        default:
          return x+x;
        }
      }
    };
    std::ostream& operator<< (std::ostream & stream, radius_sizes const & t);
    struct scales
    { // 0 = -11-, 1 = -12-, 2 = -13-, 3 = -14-, 4 = -15-
      std::array<double, 5u> value;
      std::array<bool, 5u> use;
      scales (void) : value(), use() 
      { 
        value[0u]=value[1u]=value[2u]=value[3u]=value[4u]=1.0; 
        use[0u]=use[1u]=use[2u]=use[3u]=use[4u]=true; 
      }
      void parse_line (std::string const & line);
      void to_stream (std::ostream & stream, std::string type) const;
      bool required (std::size_t i) const { return (i > 4U || (use[i] && (std::fabs(value[i]) > 0.0))); }
      bool scaled (std::size_t i) const { return (i < 5U && (std::fabs(value[i]-1.0) > 0.0)); }
      bool used (std::size_t const i) const { return (i > 4u || use[i]); }
      bool sized (std::size_t const i) const { return (fabs(value[i]-1.0) > 0.0); }
      double factor (std::size_t const i) const 
      {
        return i > 4U ? 1.0 : (std::fabs(value[i]-1.0) > 0.0) ? 1.0/value[i] : 1.0;
      }
    };

    struct factors
    {
      std::array<double, 4u> value;
      std::size_t order;
      factors (void) : value(), order() { }
      void parse_line (std::string const & line);
      void to_stream (std::ostream & stream, std::string type) const;
    };

    struct index
    {
      index_types value;
      index (void) : value(GROUP) {}
      void parse_line (std::string const & line);
    };
    std::ostream& operator<< (std::ostream & stream, index const & t);

    struct global
    {

      enum itsequence { CENTER, TWIST, LIG1, LIG2 };

      static const std::size_t improper_sequence[4], imptors_sequence[4];

      forcefield_types forcefield;
      vdw_types vdwtype;
      average_rules  radiusrule, epsilonrule;
      radius_types radiustype;
      radius_sizes radiussize;
      scales vdw_scale, chg_scale, mpole_scale, polar_scale, direct_scale, mutual_scale;
      factors bond_factor, angle_factor, opbend_factor;

      //index_types angleindex, bondindex, chargeindex, 
      //improperindex, imptorsindex, multipoleindex, 
      //opbendindex, polarizeindex, strbendindex,
      //torsionindex, ureyindex, vdwindex;

      index indices[KEYS];

      double angleunit;
      double bondunit;
      double torsionunit;
      double dielectric;
      double electric;

      global (void) 
        : angleunit(SCON_PI180*SCON_PI180), bondunit(1.0), torsionunit(1.0), dielectric(1.0), electric(332.06)
      {
        vdw_scale.use[0u] = vdw_scale.use[1u] = vdw_scale.use[2u] = false; // do not use vdw 11,12,13 per default
        vdw_scale.value[0u] = vdw_scale.value[1u] = vdw_scale.value[2u] = 0.0; // do not use vdw 11,12,13 per default
        chg_scale.use[0u] = chg_scale.use[1u] = chg_scale.use[2u] = false; // do not use charge 11,12,13 per default
        chg_scale.value[0u] = chg_scale.value[1u] = chg_scale.value[2u] = 0.0; // do not use charge 11,12,13 per default
        epsilonrule.value = average_rules::GEOMETRIC;
        indices[CHARGE].value = index_types::TYPE;
        indices[MULTIPOLE].value = index_types::TYPE;
        indices[POLARIZE].value = index_types::TYPE;
      }

      void swap (global&);

    };

    std::ostream& operator<< (std::ostream & stream, global const & g);

    

    // atoms

    struct atom
    {
      struct typeless
      { 
        bool operator() (atom const &a, atom const &b) { return a.type < b.type; }
      };
      double mass;
      std::size_t type, group, atomic, bonds;
      std::string symbol, description;
      atom (void) : mass(), type(), group(), atomic(), bonds() {}
      atom (std::string const &);
      std::size_t req_mem (void) const;
    };
    std::ostream& operator<< (std::ostream & stream, atom const & a);

    // angles

    struct angle
    {
      double f, ideal;
      std::array<std::size_t, 3u> index;
      angle (std::string const &);
      bool check (std::size_t a, std::size_t b, std::size_t c) const;
    };
    std::ostream& operator<< (std::ostream & stream, angle const & a);

    // urey bradleys

    struct ureybrad
    {
      double f, ideal;
      std::array<std::size_t, 3u> index;
      ureybrad (std::string const &);
      bool check (std::size_t a, std::size_t b, std::size_t c) const;
    };
    std::ostream& operator<< (std::ostream & stream, ureybrad const & a);

    // bonds

    struct bond
    {
      double f, ideal;
      std::array<std::size_t, 2u> index;
      bond (std::string const &);
      bool check (std::size_t a, std::size_t b) const;
    };
    std::ostream& operator<< (std::ostream & stream, bond const & a);

    // charges

    struct charge
    {
      double c;
      std::size_t index;
      charge (std::string const &);
      bool check (std::size_t a) const;
    };
    std::ostream& operator<< (std::ostream & stream, charge const & a);

    // multipoles

    struct multipole
    {
		enum axtype { NONE, Z_AXIS, Z_THEN_X, BISECTOR, Z_BISECTOR, THREEFOLD };
      double charge;
      scon::c3<double> dipole;
      scon::tsqcmatrix<double,4> quadrupole;
      std::array<std::size_t, 4u> index;
      axtype axt;
      bool has_4th_index;
      multipole (void) 
        : charge(), dipole(), quadrupole(), index(), axt(Z_THEN_X), has_4th_index(false) 
      {
      }
      multipole (std::string const &);
      void parse_line_2 (std::string const & line);
      void parse_line_3 (std::string const & line);
      void parse_line_4 (std::string const & line);
      void parse_line_5 (std::string const & line);
      void multiply (double const quadro_factor, double const dip_factor);
    };
    std::ostream& operator<< (std::ostream & stream, multipole const & a);

    // out of plane bendung

    struct opbend
    {
      double f;
      std::array<std::size_t, 4u> index;
      opbend (std::string const &);
	  bool check(std::size_t a, std::size_t b) const;
    };
    std::ostream& operator<< (std::ostream & stream, opbend const & a);

    // polarization

    struct polarize
    {
      double p, pp;
      size_t index;
	  std::array<size_t, 3u> bonded;
	  std::array<size_t, 3u> not_contracted_bond;
	  polarize(std::string const &);
    };
    std::ostream& operator<< (std::ostream & stream, polarize const & a);

    // stretch-bend

    struct strbend
    {
      double f,ff;
      std::array<std::size_t, 3u> index;
      strbend (std::string const &);
	  std::size_t check(std::size_t a, std::size_t b, std::size_t c) const;
    };

    std::ostream& operator<< (std::ostream & stream, strbend const & a);

    // torsion

    struct torsion
    {
      std::array<double, 4u> force, ideal;
      std::array<std::size_t, 4u> index, order;
      std::size_t number, max_order;
      torsion (std::string const &);
      torsion (void) : force(), ideal(), index(), order(), number(), max_order() {}
      void to_stream (std::ostream & stream, std::string type) const;
      bool empty (void) const;
      std::size_t check (std::size_t a, std::size_t b, std::size_t c, std::size_t d) const;
      std::size_t check2 (std::size_t a, std::size_t b, std::size_t c, std::size_t d) const;
    };

    // improper torsion data

    struct imptor
    { 
      std::array<double, 4u> force, ideal;
      std::array<std::size_t, 4u> order;
      std::size_t number, max_order, center;
      std::array<std::size_t, 3u> ligand;
      imptor (void) : force(), ideal(), order(), number(), max_order(), center(), ligand() {}
      imptor (std::string const &);
      void to_stream (std::ostream & stream, std::string type) const;
      bool empty (void) const;
      bool check_lig (std::size_t lig1, std::size_t lig2, std::size_t lig3) const
      {
        return (ligand[0] == lig1 && ligand[1] == lig2 && ligand[2] == lig3);
      }
    };

    struct improper
    { 
      std::array<double, 4u> force, ideal;
      std::array<std::size_t, 4u> order;
      std::size_t number, max_order, center;
      std::array<std::size_t, 3u> ligand;
      improper (void) : force(), ideal(), order(), number(), max_order(), center(), ligand() {}
      improper (std::string const &);
      void to_stream (std::ostream & stream, std::string type) const;
      bool empty (void) const;
      bool check_lig (std::size_t lig1, std::size_t lig2, std::size_t lig3) const
      {
        return (ligand[0] == lig1 && ligand[1] == lig2 && ligand[2] == lig3);
      }
    };

    // VdW data

    struct vdw
    {
      double r, e, rf;
      std::size_t index;
	    std::array<size_t, 2u> indices;
      vdw (std::string const &);
      bool check (std::size_t a) const;
	  bool check(size_t a, size_t b) const;
    };
    std::ostream& operator<< (std::ostream & stream, vdw const & a);

    namespace combi
    {
      struct vdwc
      {
        double C, E, R;
		double RR;
        vdwc (void) : C(), E(), R(),RR() {}
      };
      std::ostream& operator<< (std::ostream &stream, vdwc const &);

      struct vdw
      {
        double E, R;
      };
      std::ostream& operator<< (std::ostream &stream, vdw const &);
    }

    class parameters
    {

    public:
      typedef std::array<scon::matrix<combi::vdw, true>, 6u> vdw_matrices_t;
      typedef std::array<scon::matrix<combi::vdwc, true>, 6u> vdwc_matrices_t;

      parameters (void) : contracted(false), m_valid(false) { }

      bool valid (void) const { return m_valid; }

      void from_file (std::string const & filename);
    
      parameters contract (std::vector<std::size_t> actual_types) const;

      index_types indextype (potential_keys key) const { return m_general.indices[key].value; }
      radius_types::T radiustype (void) const { return m_general.radiustype.value; }

      std::size_t req_mem (void) const;

      std::size_t group_by_type (std::size_t const i) const { return m_atoms[i-1].group; }

      std::size_t type (std::size_t const i, potential_keys key = BOND, bool const skip_contraction = false) const
      {
        //contracted_type()
        if (contracted)
        {
          if (skip_contraction)
          {
            if (m_general.indices[key].value == index_types::GROUP)
            {
              return m_uncontracted_atoms[i-1].group;
            }
            else
              return i;
          }
          if (m_general.indices[key].value == index_types::GROUP && (contraction_map[i]-1) < m_atoms.size())
            return m_atoms[contraction_map[i]-1].group;
          else 
            return contraction_map[i];
        }
        else
        {
          if (m_general.indices[key].value == TYPE) 
            return i;
          else
            return m_atoms[i-1].group;
        }
      }

      std::size_t contract_type (std::size_t const type) const { return contraction_map[type]; }
      bool has_vdw14s (void) const { return !m_vdw14s.empty(); }
      bool vdwc_used (std::size_t const i) const { return (m_general.chg_scale.used(i) || m_general.vdw_scale.used(i)); }
      bool vdwc_scaled (std::size_t const i) const { return (m_general.chg_scale.sized(i) || m_general.vdw_scale.sized(i)); }
      bool vdwc_same (std::size_t const a, std::size_t const b) const 
      {  
        bool none_is_scaled_14(m_vdw14s.empty() || b == a || (a != 3 && b != 3 ));
        return (none_is_scaled_14 && (m_general.chg_scale.value[a] == m_general.chg_scale.value[b]) &&
                (m_general.vdw_scale.value[a] == m_general.vdw_scale.value[b]));
      }

      vdwc_matrices_t vdwc_matrices (void) const;
    
      std::vector<angle> const &      angles (void) const { return m_angles; }
      std::vector<bond> const &       bonds (void) const { return m_bonds; }
      std::vector<charge> const &     charges (void) const { return m_charges; }
      std::vector<improper> const &   impropers (void) const { return m_impropers; }
      std::vector<imptor> const &     imptors (void) const { return m_imptors; }
      std::vector<multipole> const &  multipoles (void) const { return m_multipoles; }
      std::vector<opbend> const &     opbends (void) const { return m_opbends; }
      std::vector<polarize> const &   polarizes (void) const { return m_polarizes; }
      std::vector<strbend> const &    strbends (void) const { return m_strbends; }
      std::vector<torsion> const &    torsions (void) const { return m_torsions; }
      std::vector<ureybrad> const &   ureybrads (void) const { return m_ureybrads; }
      std::vector<vdw> const &        vdws (void) const { return m_vdws; }
      std::vector<vdw> const &        vdw14s (void) const { return m_vdw14s; }
	  std::vector<combi::vdwc> const & vdwsc(void) const { return m_vdwsc; }

      double torsionunit (void) const { return m_general.torsionunit; }
      double angleunit (void) const { return m_general.angleunit; }
      double bondunit (void) const { return m_general.bondunit; }

      void swap (parameters&);

    private:



      // get combination of vdw and charge parameters for given types
      combi::vdwc    vdwc_combi (std::size_t const a, std::size_t const b, bool const use14vdw=false) const;
      // find vdw parameters of given type
      vdw const &    find_vdw (std::size_t type) const;
      vdw const &    find_vdw14 (std::size_t type) const;
      // find charge parameters of given type
      charge const & find_chg (std::size_t type) const;
      // parse parameter file line
      void parse_lines (std::vector<std::string> const & line);
      // get contracted type/group
      std::size_t contracted_type (std::size_t const i, potential_keys key = BOND)
      {
        if (m_general.indices[key].value == index_types::GROUP) return group_contraction[i];
        else return contraction_map[i];
      }



      // general stuff
      bool                    contracted, m_valid;
      global                  m_general;
      // reduction info
      std::vector<std::size_t>     m_type_of_group;
      std::vector<std::size_t>     m_reduced_types, m_reduced_groups, contraction_map, group_contraction;
      // atoms
      std::vector<atom>       m_atoms;
      // atoms
      std::vector<atom>       m_uncontracted_atoms;
      // potentials
      std::vector<angle>      m_angles;
      std::vector<bond>       m_bonds;
      std::vector<charge>     m_charges;
      std::vector<improper>   m_impropers;
      std::vector<imptor>     m_imptors;
      std::vector<multipole>  m_multipoles;
      std::vector<opbend>     m_opbends;
      std::vector<polarize>   m_polarizes;
      std::vector<strbend>    m_strbends;
      std::vector<torsion>    m_torsions;
      std::vector<ureybrad>   m_ureybrads;
      std::vector<vdw>        m_vdws;
      std::vector<vdw>        m_vdw14s;
	  std::vector<combi::vdwc> m_vdwsc;

      friend std::ostream& operator<< (std::ostream & stream, parameters const & p);

    };

  }

}
