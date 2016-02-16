/** **********************************************************
** Class declarations for the subroutine StartOpt           **
** StartOpt is a class for searching a molecule for         **
** startingstructures.                                      **
** RingSearch will search the molecule for 5, 6, and 7      **
** membered rings                                           **
** Fold will search the program for secondary structure     **
** elements                                                 **
********************************************************DW **/

#pragma once

#include "coords.h"

namespace startopt
{
  
  class Preoptimizer
  {

  protected:

    coords::Ensemble_PES m_ensemble;
    coords::Coordinates const & coords;
    coords::Coordinates m_final_coords;
    Preoptimizer& operator= (Preoptimizer const &);

  public:

    Preoptimizer (coords::Coordinates const & co) :
      m_ensemble(), coords(co), m_final_coords(co)
    {}
    virtual ~Preoptimizer() { }
    virtual void generate (coords::Ensemble_PES const & init_ensemble, std::size_t const multiplier) = 0;
    coords::Coordinates & final_coords (void) { return m_final_coords; }
    coords::Ensemble_PES const & PES (void) const { return m_ensemble; }

  };

  void apply(coords::Coordinates & c, coords::Ensemble_PES & e);

} /* / namespace_startopt */
