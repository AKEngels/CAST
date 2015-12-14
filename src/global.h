/*
  Forward declarations.
*/
#pragma once 

#include <iostream>
#include <cstdlib>
#include <cmath>

#if !defined(USE_CUDA) && defined COMPILE_WITH_CUDA
#define USE_CUDA double
#endif

#if defined(_MSC_VER)
#include <yvals.h>
#endif

#if __cplusplus > 199711L
#define DCPP11_SUPPORT 1
#elif defined(_MSC_VER) && defined(_HAS_CPP0X)
#define DCPP11_SUPPORT 1
#endif

#if defined(DCPP11_SUPPORT)
#define DNULL nullptr
#else
#define DNULL NULL
#endif

#if defined(_MSC_VER) && defined(_M_X64)
#define COMPILEX64 
#elif !defined(_MSC_VER) && (defined(__LP64__) || defined(__x86_64))
#define COMPILEX64
#else
#endif

#if defined(TERACHEM_MPI)
#define USE_MPI
#endif

static const double
  PI = 3.1415926535897932384626433,
  RATIOPI180 = 0.01745329251994329576923690768489, 
  RATIO180PI = 57.295779513082320876798154814105,
  D_RAND_MAX = static_cast<double>(RAND_MAX);

static char const lineend = '\n';

//struct lineend { std::nullptr_t x; };
//std::ostream & operator<< (std::ostream& strm, lineend const &l) { strm.put('\n'); }

inline double d_rand (void) { return static_cast<double>(std::rand())/D_RAND_MAX; }
inline double angle_to_fullcircle (double const angle) { return (angle >= 360.0 ? (angle-360.0) : (angle < 0.0 ? (angle+360.0) : angle)); }
//inline double angle_to_range (double const angle) { return (angle > 180.0 ? (angle-360.0) : (angle < -180.0 ? (angle+360.0) : angle)); }
inline double angle_to_range(double const value) { return (value >= 180.0 || value <= -180.0) ? (value - std::floor((value+180.0) / 360.0)*360.0) : value; }

inline double ranged_angle_rotation(double const from, double const to)
{
  return angle_to_range(angle_to_fullcircle(from) - angle_to_fullcircle(to));
}

template <class T>
inline std::size_t num_digits (T number)
{
  std::size_t digits(0);
  while (number) {
    number /= 10;
    ++digits;
  }
  return digits;
}

namespace input
{
  class format;
  namespace formats
  {
    class tinker;
    class tinker_multi;
  }
}

namespace linked
{
  class environment_iterator;
}
class LinkedCells;

namespace md 
{ 
  class simulation; 
}

namespace startopt 
{ 
  class Preoptimizer; 
  namespace preoptimizers
  {
    class Solvadd;
    class R_evolution;
    //class Fold;
  }
}

namespace solvadd 
{ 
  struct configurables;
  class solvation; 
}

namespace rdf 
{ 
  class distribution; 
}

namespace dimermethod 
{ 
  class dimer; 
}

namespace energy 
{ 
  class interface_base;
  namespace interfaces 
  { 
    namespace aco
    {
      class aco_ff;
    }
    namespace mopac
    {
      class sysCallInterface;
    }
    //namespace tinker
    //{
    //  class tinkerCaller;
    //}
    namespace terachem
    {
      class mpiInterface;
    }
	  namespace amoeba
	  {
		  class amoeba_ff;
	  }
  } 
}

// ###################################### CAST 3 Stuff

namespace coords
{
  namespace input
  {
    class format;
    namespace formats
    {
      class tinker;
      class amber;
    }
  }
  namespace output
  {
    class format;
    namespace formats
    {
      class tinker;
      class zmatrix;
    }
  }
  class Coordinates;
}

namespace tinker
{
  namespace parameter
  {
    namespace combi
    {
      struct vdw;
      struct vdwc;
    }
    struct atom;
    struct angle;
    struct bond;
    struct charge;
    struct imptor;
    struct improper;
    struct multipole;
    struct opbend;
    struct polarize;
    struct strbend;
    struct torsion;
    struct ureybrad;
    struct vdw;
    class parameters;
  }
  namespace refine
  {
    namespace types
    {
      struct binary_quadratic;
      struct ternary_quadratic;
      struct torsion;
      struct improper;
      struct imptor;
      struct multipole;
      struct opbend;
      struct polarize;
      struct strbend;
      struct nbpair;
    }
    class refined;
  }
  
}

namespace optimization
{
  namespace global
  {
    struct Tabu_Point;
    struct Tabu_List;
    class optimizer;
    namespace optimizers
    {
      class monteCarlo;
      class tabuSearch;
    }
  }
  namespace local
  {
    namespace linesearch
    {
      template<class callback>
      struct more_thuente;
      template<class callback>
      struct none;
    }
    template<class LineSearch, class Log> class lbfgs;
    template<class LineSearch, class Log> class conjugate_gradient;
  }
}
