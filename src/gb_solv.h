#ifndef GB_SOLV_HEADER
#define GB_SOLV_HEADER

#include <vector>

#include "atomic.h"
#include "coords.h"
#include "tinker_parameters.h"
#include "configuration.h"

#ifndef GBSA_PI
#define GBSA_PI 3.14159265358979323846264338327950288419716939937510
#endif

namespace gbsa
{

  enum METHODS  { VAC = -1, STILL = 0, HCT, OBC, GRYCUK, ACE, ONION, METHODNUM };
  enum SURFACES { TINKER = 0, SASASTILL, GAUSS, SURFACESNUM };
  enum RADIUS   { STD = 0, VDW };


  // GAUSS berechnet die exakte Oberfläche von sich überschneidenden
  // Kugeln nach dem Theorem von Gauss-Bonnet.
  //
  //     literature references :
  //
  //     T.J.Richmond, "Solvent Accessible Surface Area and
  //     Excluded Volume in Proteins", Journal of Molecular Biology,
  //     178, 63 - 89 (1984)
  //
  //     L.Wesson and D.Eisenberg, "Atomic Solvation Parameters
  //     Applied to Molecular Dynamics of Proteins in Solution",
  //     Protein Science, 1, 227 - 235 (1992)
  //

  class gborn
  {

    coords::float_type probe;

    std::vector<coords::float_type> born_radius;
    std::vector<coords::float_type> solvation_radius;
    std::vector<coords::float_type> solvent_access;

    tinker::parameter::parameters const * param_ptr;
    coords::Coordinates const * coord_ptr;

    coords::float_type solv_rad (size_t const atom_index) const;

    void solv_radii();

    void sa_tinker();

    void sa_still();


  public:

    gborn(coords::Coordinates const &co, tinker::parameter::parameters const &pp)
      : probe(1.4), param_ptr(&pp), coord_ptr(&co)
    { }

  };



}

#endif
