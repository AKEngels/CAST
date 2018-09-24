#pragma once

#define _USE_MATH_DEFINES

#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<string>
#include<random>
#include<vector>
#include<math.h>

#include "configuration.h"


//namespace for ecitonDiffusion with dimer approach
namespace exciD {

  const double h_quer = (4.135667662e-15)/ (2 * M_PI);//eV * s
  const double boltzmann_const = 8.6173303e-5; //  gauss units

  struct Couplings
  {
    std::size_t monA;
    std::size_t monB;
    double coupling;
    coords::Cartesian_Point position;

    /**default constructor*/
    Couplings() : monA(), monB(), coupling() {}

    /** another constructor taking members as parameters*/
    Couplings(std::size_t monAp, std::size_t monBp, double couplingp) :
      monA(monAp), monB(monBp), coupling(couplingp) {}

    Couplings(std::size_t monAp, std::size_t monBp, double couplingp, coords::Cartesian_Point pos) :
      monA(monAp), monB(monBp), coupling(couplingp), position(pos) {}

  };

  struct Partners
  {
    std::size_t partnerIndex;
    std::vector<std::size_t> connect;
    double avgCoup;

    Partners() : partnerIndex(), connect(),avgCoup() {}

    Partners(std::size_t partnerIndexp, std::vector<std::size_t> connectp) :
      partnerIndex(partnerIndexp), connect(connectp) {}

    Partners(std::size_t partnerIndexp, std::vector<std::size_t> connectp, double avgCoupp) :
      partnerIndex(partnerIndexp), connect(connectp), avgCoup(avgCoupp) {}

    Partners(std::size_t partnerIndexp) :
      partnerIndex(partnerIndexp) {}
  };

  struct Exciton
  {
    std::size_t location;
    char state;

    Exciton() : location(), state('e') {} //state of exciton is by default set to e for Exciton in constructor, is changed to different states when exciton changes (vgl excitonbreakup task)

    Exciton(std::size_t locationp, char statep) :
      location(locationp), state(statep) {}

    Exciton(std::size_t locationp) :
      location(locationp), state('e') {} //state of exciton is by default set to e for Exciton in constructor, is changed to different states when exciton changes (vgl excitonbreakup task)
  };
 
  void dimexc(std::string, std::string, double, double);
  coords::Cartesian_Point avgDimCoM(coords::Cartesian_Point, coords::Cartesian_Point);
  double length(coords::Cartesian_Point, coords::Cartesian_Point);
  double marcus(double, double, double);
  coords::Cartesian_Point structCenter(std::vector <exciD::Couplings>);
  coords::Cartesian_Point min(std::vector<exciD::Couplings>);
}