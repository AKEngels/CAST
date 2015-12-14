#pragma once

#include <vector>

#include "global.h"
#include "configuration.h"              
#include "coords.h"



namespace dimermethod
{
  class kaestner_sherwood_dih
  {
  public:
    struct dimer_dih
    {
      double e0, e1, C;
      coords::Representation_Main tau, x0, x1;
      coords::Gradients_Main f0, f1, f_rot, f_trans;
      coords::Coordinates & coordobject;
      dimer_dih(coords::Coordinates &cobj, coords::Representation_Main const &, coords::Gradients_Main const &);
      void span (std::vector<coords::Representation_Main> const * rep_ptr = nullptr);
      void update (bool const update_x0 = false, bool const update_x1 = true);
      void rotate(coords::Representation_Main const &O, double const phi);
      void invert (void);
      void move_to(coords::Representation_Main const &new_x0);
      void advance_along_axis (coords::float_type const d);
      size_t linesearch (void);
      size_t size (void) const { return x0.size(); }
    private:
      void update_force(coords::Representation_Main const &dihedrals, double &energy, coords::Gradients_Main &forces);
      dimer_dih& operator= (dimer_dih const &);
    };
  private:
    dimer_dih m_dimer;
    void update_force(coords::Representation_Main const &dihedrals, double &energy, coords::Gradients_Main &forces);
    void rotate (void);
  public:
    kaestner_sherwood_dih(coords::Coordinates &coordinates_object, std::vector<coords::Representation_Main> const * rep_ptr = nullptr, double const x0_precalced = true);
    bool rotate_to_min_curve_cg (void);
    void translate_across_transition_cg (void);
    coords::Representation_Main translate_across_transition_lbfgs(void);
  };

}
