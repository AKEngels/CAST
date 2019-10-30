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

  const double h_quer = (4.135667662e-15) / (2 * M_PI);//eV * s
  const double boltzmann_const = 8.6173303e-5; //  gauss units

  struct Couplings
  {
    std::size_t monA;
    std::size_t monB;
    double coupling;
    double seccoupling;
    coords::Cartesian_Point position;

    /**default constructor*/
    Couplings() : monA(), monB(), coupling(), seccoupling() {}

    /** another constructor taking members as parameters*/
    Couplings(std::size_t monAp, std::size_t monBp, double couplingp) :
      monA(monAp), monB(monBp), coupling(couplingp), seccoupling(0.) {}

    Couplings(std::size_t monAp, std::size_t monBp, double couplingp, double seccouplingp) :
      monA(monAp), monB(monBp), coupling(couplingp), seccoupling(seccouplingp) {}

    Couplings(std::size_t monAp, std::size_t monBp, double couplingp, double seccouplingp, coords::Cartesian_Point pos) :
      monA(monAp), monB(monBp), coupling(couplingp), seccoupling(seccouplingp), position(pos) {}

  };

  struct Partners
  {
    std::size_t partnerIndex;
    std::vector<std::size_t> connect;
    double avgCoup, avgsecCoup;

    Partners() : partnerIndex(), connect(), avgCoup(), avgsecCoup() {}

    Partners(std::size_t partnerIndexp, std::vector<std::size_t> connectp) :
      partnerIndex(partnerIndexp), connect(connectp), avgCoup(), avgsecCoup() {}

    Partners(std::size_t partnerIndexp, std::vector<std::size_t> connectp, double avgCoupp, double avgsecCoupp) :
      partnerIndex(partnerIndexp), connect(connectp), avgCoup(avgCoupp), avgsecCoup(avgsecCoupp) {}

    Partners(std::size_t partnerIndexp) :
      partnerIndex(partnerIndexp), connect(), avgCoup(), avgsecCoup() {}
  };

  struct Exciton
  {
    std::size_t location;
    char state;
    std::size_t h_location, location_lastS, h_location_lastS;

    Exciton() : location(), state('e'), h_location(), location_lastS(), h_location_lastS() {} //state of exciton is by default set to e for Exciton in constructor, is changed to different states when exciton changes (vgl excitonbreakup task)

    Exciton(std::size_t locationp, char statep) :
      location(locationp), state(statep), h_location(), location_lastS(), h_location_lastS() {}

    Exciton(std::size_t locationp) :
      location(locationp), state('e'), h_location(), location_lastS(), h_location_lastS() {} //state of exciton is by default set to e for Exciton in constructor, is changed to different states when exciton changes (vgl excitonbreakup task)
  };

  void dimexc(std::string, std::string, std::size_t, int, char, double, std::size_t);
  coords::Cartesian_Point avgDimCoM(coords::Cartesian_Point, coords::Cartesian_Point);
  double length(coords::Cartesian_Point, coords::Cartesian_Point);
  double marcus(double, double, double);
  double coulomb(coords::Cartesian_Point, coords::Cartesian_Point, double);
  coords::Cartesian_Point structCenter(coords::Representation_3D);
  coords::Cartesian_Point min(coords::Representation_3D);
  coords::Cartesian_Point max(coords::Representation_3D);
}