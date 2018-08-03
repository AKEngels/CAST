#pragma once
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<string>
#include<random>
#include<vector>

#include "configuration.h"


//namespace for ecitonDiffusion with dimer approach
namespace exciD {

  struct Couplings
  {
    std::size_t monA;
    std::size_t monB;
    std::vector<double> coupling;

    /**default constructor*/
    Couplings() : monA(), monB(), coupling() {}

    /** another constructor taking members as parameters*/
    Couplings(std::size_t monAp, std::size_t monBp, std::vector<double> couplingp) :
      monA(monAp), monB(monBp), coupling(couplingp) {}

  };
 
  void dimexc(std::string, std::string);
}