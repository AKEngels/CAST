#pragma once
#include <vector>

using float_type = double;

namespace md
{
  /** Nose-Hover thermostat. Variable names and implementation are identical to the book of
  Frenkel and Smit, Understanding Molecular Simulation, Appendix E */
  struct nose_hoover_2chained
  {
    double v1, v2, x1, x2;
    double Q1, Q2, G1, G2;
    nose_hoover_2chained(void) :
      v1(0.0), v2(0.0), x1(0.0), x2(0.0),
      Q1(0.1), Q2(0.1), G1(0.0), G2(0.0)
    { }

    void setQ1(double Q1inp) { Q1 = Q1inp; }
    void setQ2(double Q2inp) { Q2 = Q2inp; }
  };

  class nose_hoover_arbitrary_length
  {
    public:
    nose_hoover_arbitrary_length::nose_hoover_arbitrary_length() : chainlength(5u)
    {
      epsilons = std::vector<float_type>(chainlength,0.);
      velocities = std::vector<float_type>(chainlength, 0.);
      forces = std::vector<float_type>(chainlength, 0.);
      masses_param_Q = std::vector<float_type>(chainlength, 1.);
      // According to https://doi.org/10.1063/1.463940 M=5 and Q=1 converges fastes to the canonical distribution.
    }
    nose_hoover_arbitrary_length::nose_hoover_arbitrary_length(std::size_t chainlength_) : chainlength(chainlength_)
    {
      epsilons = std::vector<float_type>(chainlength, 0.);
      velocities = std::vector<float_type>(chainlength, 0.);
      forces = std::vector<float_type>(chainlength, 0.);
      masses_param_Q = std::vector<float_type>(chainlength, 1.);
    }
    nose_hoover_arbitrary_length::nose_hoover_arbitrary_length(std::vector<float_type> const& masses_param_Q_) : chainlength(masses_param_Q_.size())
    {
      epsilons = std::vector<float_type>(chainlength, 0.);
      velocities = std::vector<float_type>(chainlength, 0.);
      forces = std::vector<float_type>(chainlength, 0.);
      masses_param_Q = std::vector<float_type>(masses_param_Q_);
    }
    const std::size_t chainlength;
    std::vector<float_type> epsilons; // "positions" of each thermostat
    std::vector<float_type> velocities; // "velocities" of each thermostat
    std::vector<float_type> forces; // "G" or forces acting on each thermostat
    std::vector<float_type> masses_param_Q; // "Q" or masses of each thermostat
  };
}