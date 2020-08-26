#pragma once
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<limits>
#include<string>
#include<random>
#include<vector>

#include "constants.h"
#include "helperfunctions.h"

//Original code by Charlotte. German comments are from original code, english comments were inserted during implementation into CAST.
//Arrays were replaced by vectors.



namespace XB
{

  inline double rate(double const coupling, double const deltaG, double const reorganisation, double const temperatureInK = 298.)
  {
    constexpr double pi = constants::pi;
    constexpr double h_quer = constants::h_bar_gaussian_units;
    constexpr double boltzmann_constant_kb = constants::boltzmann_constant_kb_gaussian_units; //  in gauß einheiten // Dustin July19: is in eV/K
    const double prefactor = (coupling * coupling) / h_quer * sqrt(pi / (reorganisation * boltzmann_constant_kb * temperatureInK));
    const double exponential_part = std::exp(-(reorganisation + deltaG) * (reorganisation + deltaG) / (4. * boltzmann_constant_kb * temperatureInK * reorganisation));
    const double l = prefactor * exponential_part;
    return l;
  }


  class ExcitonBreakup
  {
  public:
    //ExcitonBreakup(std::string filename);

    ExcitonBreakup(std::string masscenters, std::string nscpairrates, std::string pscpairexrates, std::string pscpairchrates, std::string pnscpairrates)//LEGACY
      : totalNumberOfMonomers(0u), reorganisationsenergie_exciton(Config::get().exbreak.ReorgE_exc), reorganisationsenergie_ladung(Config::get().exbreak.ReorgE_ch),
      fullerenreorganisationsenergie(Config::get().exbreak.ReorgE_nSC), ct_reorganisation(Config::get().exbreak.ReorgE_ct), chargetransfertriebkraft(Config::get().exbreak.ct_triebkraft),
      rekombinationstriebkraft(Config::get().exbreak.rek_triebkraft), rek_reorganisation(Config::get().exbreak.ReorgE_rek), oszillatorstrength(Config::get().exbreak.oscillatorstrength),
      wellenzahl(Config::get().exbreak.wellenzahl), k_rad(wellenzahl* wellenzahl* oszillatorstrength), 
      startingPscaling(Config::get().exbreak.startingPscaling), nbrStatingpoins(Config::get().exbreak.nbrStatingpoins),//added options for startingPoint-algorithm
      avg_position_total__x(0.), avg_position_total__y(0.), avg_position_total__z(0.), numberOf_p_SC(0u), numberOf_n_SC(0u), numberOfStartingPoints(0u + 1u),
      avg_position_p_sc__x(0.), avg_position_p_sc__y(0.), avg_position_p_sc__z(0.), avg_position_n_sc__x(0.), avg_position_n_sc__y(0.), avg_position_n_sc__z(0.)
    {
      this->read(Config::get().exbreak.pscnumber, Config::get().exbreak.nscnumber, masscenters, nscpairrates, pscpairexrates, pscpairchrates, pnscpairrates);
    };

    void runAndWrite(char direction)
    {
      this->writeAuxFiles(direction);
      this->run(direction);
      this->analyseResults();
    }

    std::vector <std::size_t> calculateStartingpoints(char direction) const;

    void run(
      char direction,
      std::size_t numberOfRunsPerStartingPoint = 100u,
      //std::size_t const maxNumSteps = 0u,
      double excitonicDrivingForce_GaussianSigma = 0.0338987,
      double chargecarrierDrivingForce_GaussianSigma = 0.068584577);

    void analyseResults(std::size_t numberOfRunsPerStartingPoint = 100u) const;

    void writeAuxFiles(char direction) const;

    std::size_t getTotalNumberOfMonomers() const
    {
      return this->totalNumberOfMonomers;
    }

#ifndef GOOGLE_MOCK
  private:
    ExcitonBreakup();
#endif
    void read(std::size_t numberOf_p_SC_, std::size_t numberOf_n_SC_, std::string masscenters, std::string nscpairrates, std::string pscpairexrates, std::string pscpairchrates, std::string pnscpairrates);

    void processAfterFilereading();

    double evaluateCoulomb(std::size_t particle1_iterator, std::size_t particle2_iterator, double e_relative) const;

    // Definition der length-berechnung
    inline double distance(std::size_t particle1_iterator, std::size_t particle2_iterator) const
    {
      std::size_t const& p = particle1_iterator;
      std::size_t const& q = particle2_iterator;
      std::vector <double> const& arr1 = this->x;
      std::vector <double> const& arr2 = this->y;
      std::vector <double> const& arr3 = this->z;
      const double l = std::sqrt((arr1[p] - arr1[q]) * (arr1[p] - arr1[q]) + (arr2[p] - arr2[q]) * (arr2[p] - arr2[q]) + (arr3[p] - arr3[q]) * (arr3[p] - arr3[q]));
      return l;
    }

    std::size_t totalNumberOfMonomers;

    const double reorganisationsenergie_exciton;
    const double reorganisationsenergie_ladung;
    const double fullerenreorganisationsenergie;
    const double ct_reorganisation;
    const double chargetransfertriebkraft;
    const double rekombinationstriebkraft;
    const double rek_reorganisation;
    const double oszillatorstrength;
    const double wellenzahl;
    const double k_rad;
    
    const double startingPscaling;
    const double nbrStatingpoins;
    double avg_position_total__x, avg_position_total__y, avg_position_total__z;
    
    std::size_t numberOf_p_SC, numberOf_n_SC;
    std::size_t numberOfStartingPoints; // Number of starting points

    double avg_position_p_sc__x, avg_position_p_sc__y, avg_position_p_sc__z, avg_position_n_sc__x, avg_position_n_sc__y, avg_position_n_sc__z;

    std::size_t numberOfExcitonPairs, numberOfNSemiconductorHomopairs, numberOfHeteroDimers;
    std::vector <std::size_t> numberOfPartnerPerMonomer; // How many partners does one specific monomer hvae?
    std::vector <std::vector<std::size_t>> partner; // Matrices for accessing the partners
    std::vector <std::vector<double>> coupling_exciton;
    std::vector <std::vector<double>> coupling_ladung;
    std::vector <std::vector<double>> coupling_ct;
    std::vector <std::vector<double>> coupling_rek;
    std::vector <std::vector<double>> coupling_fulleren;
    std::vector <double> x, y, z; // Coordinates of the mass-points of each monomer
    std::vector <std::size_t> startpunkt; // Iterator-numbers of the starting points

    struct results {
      std::vector <std::size_t> ex_diss, ch_diss, rek, trapping, radiativ;
      std::vector <std::vector<double>> vel_ex, vel_ch, zeit_ex, zeit_ch;
      results(std::size_t index, std::size_t numberOfRunsPerStartingPoint)
        : ex_diss(index + 1), ch_diss(index + 1), rek(index + 1), trapping(index + 1), radiativ(index + 1),
        vel_ex(index + 1, std::vector <double>(numberOfRunsPerStartingPoint)),
        vel_ch(index + 1, std::vector <double>(numberOfRunsPerStartingPoint)),
        zeit_ex(index + 1, std::vector <double>(numberOfRunsPerStartingPoint)),
        zeit_ch(index + 1, std::vector <double>(numberOfRunsPerStartingPoint))
      {
        for (std::size_t i = 0u; i < (index + 1); i++) //initializing the vectors with 0
        {
          this->ex_diss[i] = 0;
          this->ch_diss[i] = 0;
          this->rek[i] = 0;
          this->radiativ[i] = 0;
          this->trapping[i] = 0;
        }
      }
      results() {}
    };
    results m_results;
  };

}

