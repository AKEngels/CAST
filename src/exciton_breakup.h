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

  inline double rate(double coupling, double deltaG, double reorganisation)
  {
    constexpr double pi = constants::pi;
    constexpr double h_quer = constants::h_quer;
    constexpr double boltzmann_constant_kb = constants::boltzmann_constant_kb; //  in gauß einheiten // Dustin July19: is in eV/K
    const double l = (coupling*coupling) / h_quer * sqrt(pi / (reorganisation*boltzmann_constant_kb * 298.))*exp(-(reorganisation + deltaG)*(reorganisation + deltaG) / (4. * boltzmann_constant_kb * 298. * reorganisation));
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
      wellenzahl(Config::get().exbreak.wellenzahl), k_rad(wellenzahl * wellenzahl*oszillatorstrength), numberOfRunsPerStartingPoint(101u),
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
#ifndef GOOGLE_MOCK
  private:
    ExcitonBreakup();
#endif
    void read(std::size_t numberOf_p_SC_, std::size_t numberOf_n_SC_, std::string masscenters, std::string nscpairrates, std::string pscpairexrates, std::string pscpairchrates, std::string pnscpairrates)
    {
      this->numberOf_p_SC = numberOf_p_SC_;
      this->numberOf_n_SC = numberOf_n_SC_;
      /////////////////////////////////// INPUT-READING
      std::ifstream com_file;
      com_file.open(masscenters);

      com_file >> totalNumberOfMonomers;
      std::cout << "Read number of Monomers: " << totalNumberOfMonomers << std::endl;
      if (totalNumberOfMonomers == numberOf_p_SC + numberOf_n_SC) //test if correct number of molecules was given
      {
        std::cout << "Number of monomers is correct, proceeding." << std::endl;
      }
      else //totalNumberOfMonomers != numberOf_p_SC + numberOf_n_SC
      {
        std::cout << "Wrong number of monomers detected!" << std::endl;
        throw std::logic_error("Wrong numbers of p- and n-type molecules in inputfile. Expected: " + std::to_string(totalNumberOfMonomers) + " | Is: " + std::to_string(numberOf_p_SC + numberOf_n_SC));
      }
      com_file >> skipline;
      x = std::vector <double>(totalNumberOfMonomers + 1);
      y = std::vector <double>(totalNumberOfMonomers + 1);
      z = std::vector <double>(totalNumberOfMonomers + 1);

      for (std::size_t i = 1u; i < (totalNumberOfMonomers + 1u); i++)
      {
        std::string zeile;
        com_file >> zeile >> x[i] >> y[i] >> z[i]; //reading and saving balance points for molecules
      }
      ///////////////////////////////////
      std::ifstream exciton;
      exciton.open(pscpairexrates);
      numberOfExcitonPairs = 0u;
      std::string zeile;
      while (getline(exciton, zeile))
      { //counting of excitonpairs in homodimers
        numberOfExcitonPairs++;
      }
      exciton.close();
      ///////////////////////////////////

      coupling_exciton = std::vector <std::vector<double>>(totalNumberOfMonomers + 1, std::vector <double>(totalNumberOfMonomers + 1)); //in original code 2d-arrays were used
      coupling_ladung = std::vector <std::vector<double>>(totalNumberOfMonomers + 1, std::vector <double>(totalNumberOfMonomers + 1));
      coupling_ct = std::vector <std::vector<double>>(totalNumberOfMonomers + 1, std::vector <double>(totalNumberOfMonomers + 1));
      coupling_rek = std::vector <std::vector<double>>(totalNumberOfMonomers + 1, std::vector <double>(totalNumberOfMonomers + 1));
      coupling_fulleren = std::vector <std::vector<double>>(totalNumberOfMonomers + 1, std::vector <double>(totalNumberOfMonomers + 1));

      for (std::size_t i = 1u; i < (totalNumberOfMonomers + 1u); i++) //initializaton of the 2d vectors with 0 in all places
      {
        for (std::size_t j = 1u; j < (totalNumberOfMonomers + 1u); j++) {
          coupling_exciton[i][j] = 0;
          coupling_ladung[i][j] = 0;
          coupling_ct[i][j] = 0;
          coupling_rek[i][j] = 0;
          coupling_fulleren[i][j] = 0;
        }
      }

      std::vector <int> exciton_1(numberOfExcitonPairs + 1), exciton_2(numberOfExcitonPairs + 1); //vectors for exciton pairs

      exciton.open(pscpairexrates);
      for (std::size_t i = 1u; i < (numberOfExcitonPairs + 1u); i++)
      {
        exciton >> exciton_1[i] >> exciton_2[i];
        exciton >> coupling_exciton[exciton_1[i]][exciton_2[i]];
        coupling_exciton[exciton_2[i]][exciton_1[i]] = coupling_exciton[exciton_1[i]][exciton_2[i]];
      }

      exciton.close();
      numberOfPartnerPerMonomer = std::vector <size_t>(numberOf_p_SC + numberOf_n_SC + 1);
      for (std::size_t i = 1u; i < (numberOf_p_SC + numberOf_n_SC + 1); i++)
      {
        numberOfPartnerPerMonomer[i] = 0u;
      } //inizialisation of all elements with 0

      for (std::size_t i = 1u; i < (numberOfExcitonPairs + 1); i++) // counting homodimer partners for j
      {
        for (std::size_t j = 1u; j < (numberOf_p_SC + 1); j++)
        {
          if ((exciton_1[i] == j) || (exciton_2[i] == j))
          {
            numberOfPartnerPerMonomer[j]++;
          }
        }
      }
      ////////////////////////////////////
      exciton.open(pscpairchrates);
      for (std::size_t i = 1u; i < (numberOfExcitonPairs + 1); i++)
      {
        std::size_t h = 0u;
        std::size_t j = 0u;
        exciton >> j >> h; //j=exciton_1[i], h=exciton_2[i]
        exciton >> coupling_ladung[j][h];
        coupling_ladung[h][j] = coupling_ladung[j][h];
      }
      exciton.close();

      numberOfHeteroDimers = 0u;
      exciton.open(pnscpairrates);
      while (getline(exciton, zeile)) //counting of heterodimers
      {
        numberOfHeteroDimers++;
      }
      exciton.close();

      std::vector <int> hetero_1(numberOfHeteroDimers + 1), hetero_2(numberOfHeteroDimers + 1);

      exciton.open(pnscpairrates);
      for (std::size_t i = 1; i < (numberOfHeteroDimers + 1); i++)
      {
        exciton >> hetero_1[i] >> hetero_2[i];
        exciton >> coupling_ct[hetero_1[i]][hetero_2[i]] >> coupling_rek[hetero_1[i]][hetero_2[i]];
        coupling_rek[hetero_2[i]][hetero_1[i]] = coupling_rek[hetero_1[i]][hetero_2[i]];
        coupling_ct[hetero_2[i]][hetero_1[i]] = coupling_ct[hetero_1[i]][hetero_2[i]];
      }
      for (std::size_t i = 1u; i < (numberOfHeteroDimers + 1); i++) //counting of heteropartners for j and adding to known number of partners
      {
        for (std::size_t j = 1; j < (totalNumberOfMonomers + 1); j++)
        {
          if ((hetero_1[i] == j) || (hetero_2[i] == j))
          {
            numberOfPartnerPerMonomer[j]++;
          }
        }
      }
      exciton.close();

      ////////////////////////////////////////
      exciton.open(nscpairrates);
      numberOfNSemiconductorHomopairs = 0u;
      while (getline(exciton, zeile)) //counting fullerene homopairs
      {
        numberOfNSemiconductorHomopairs++;
      }
      exciton.close();

      std::cout << "Number of n-semiconductor pairs " << numberOfNSemiconductorHomopairs << std::endl;
      std::vector <int> fulleren_1(numberOfNSemiconductorHomopairs + 1), fulleren_2(numberOfNSemiconductorHomopairs + 1), test(2000);

      exciton.open(nscpairrates);
      for (std::size_t i = 1; i < (numberOfNSemiconductorHomopairs + 1); i++)
      {
        exciton >> fulleren_1[i] >> fulleren_2[i];
        test[i] = fulleren_2[i];
        exciton >> coupling_fulleren[fulleren_1[i]][fulleren_2[i]];
        coupling_fulleren[fulleren_2[i]][fulleren_1[i]] = coupling_fulleren[fulleren_1[i]][fulleren_2[i]];
      }

      exciton.close();

      for (std::size_t i = 1; i < (numberOfNSemiconductorHomopairs + 1); i++) //counting of fullerenhomopartners and adding to known partners
      {
        for (std::size_t j = 1; j < (totalNumberOfMonomers + 1); j++)
        {
          if ((fulleren_1[i] == j) || (fulleren_2[i] == j))
          {
            numberOfPartnerPerMonomer[j]++;
          }
        }
      }

      partner = std::vector<std::vector<std::size_t>>(totalNumberOfMonomers + 1, std::vector<std::size_t>());//2D-vector with variing length for second vector

      for (std::size_t i = 1; i < (totalNumberOfMonomers + 1); i++) //dynamic allocation for length of 2nd vector
      {
        partner[i].resize(numberOfPartnerPerMonomer[i] + 1);
      }

      for (std::size_t i = 1; i < (totalNumberOfMonomers + 1); i++) //initializing all elements of vector partner with 0
      {
        for (std::size_t j = 1; j < (numberOfPartnerPerMonomer[i] + 1); j++) {
          partner[i][j] = 0;
        }
      }

      for (std::size_t i = 1; i < (totalNumberOfMonomers + 1); i++)
      {
        std::size_t j = 1; //j is here the number of the partner to particle i 

        for (std::size_t h = 1; h < (numberOfExcitonPairs + 1); h++)  //e = number of exciton-pairs [homodimer-pairs?]
        {
          if (exciton_1[h] == i)
          {
            partner[i][j] = exciton_2[h];
            j++; //j is always incremented when an element is added to the list of partners of particle i
          }
          if (exciton_2[h] == i)
          {
            partner[i][j] = exciton_1[h];
            j++;
          }
        }

        for (std::size_t h = 1; h < (numberOfHeteroDimers + 1); h++) //het = number of heterodimer-pairs
        {
          if (hetero_1[h] == i)
          {
            partner[i][j] = hetero_2[h];
            j++;
          }
          if (hetero_2[h] == i)
          {
            partner[i][j] = hetero_1[h];
            j++;
          }
        }

        for (std::size_t h = 1; h < (numberOfNSemiconductorHomopairs + 1); h++) //full = number of fullerene homopairs
        {
          if (fulleren_1[h] == i)
          {
            partner[i][j] = fulleren_2[h];
            j++;
          }
          if (fulleren_2[h] == i)
          {
            partner[i][j] = fulleren_1[h];
            j++; // since the 2nd dimension length of vector partner was set to numberOfPartnerPerMonomer[i] in thes logic construction j must always end up to be equal to numberOfPartnerPerMonomer[i]
          }
        }
        if (numberOfPartnerPerMonomer[i] != j - 1) //Sanity Check
        {
          std::cout << "Error with number of partners for monomer " << i << std::endl;
          throw std::runtime_error("Error with number of partners for monomer " + std::to_string(i) + ". Aborting.");
        }
      }

      // INPUT-END


      //////////////////////////////////////////////////////////////////////
      this->processAfterFilereading();
    }

    void writeAuxFiles(char direction) const
    {
      std::ofstream kopplung;
      kopplung.open("partner.txt"); // Writing file partner.txt
      for (std::size_t i = 1u; i < (totalNumberOfMonomers + 1); i++)
      {
        kopplung << std::setw(6) << i << std::setw(6) << numberOfPartnerPerMonomer[i]; //writes the indices of the molecules and the ammount of partners they posess
        for (std::size_t j = 1; j < (numberOfPartnerPerMonomer[i] + 1); j++)
        {
          kopplung << std::setw(6) << partner[i][j]; //writes the indices of the partners j
        }
        kopplung << '\n';
      }
      kopplung.close();

      kopplung.open("couplings.txt"); // Writing file couplings.txt
      for (std::size_t i = 1u; i < (numberOf_p_SC + 1); i++)
      {
        kopplung << std::setw(6) << i << std::setw(6) << numberOfPartnerPerMonomer[i]; //writes the indices of the molecules and the ammount of partners they posess

        for (std::size_t j = 1; j < (numberOfPartnerPerMonomer[i] + 1); j++)
        {
          kopplung << std::setw(6) << partner[i][j]; //writes the indices of the partners j
          if (partner[i][j] < (numberOf_p_SC + 1))
          {
            kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_exciton[i][partner[i][j]]; //writes the exciton-coupling between i and j
          }
          else if (partner[i][j] > numberOf_p_SC)
          {
            kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_ct[i][partner[i][j]];     //writes charge-transfer-coupling between i and j 
          }
        }
        kopplung << '\n';
      }

      for (std::size_t i = (numberOf_p_SC + 1); i < (totalNumberOfMonomers + 1); i++)
      {
        kopplung << std::setw(6) << i << std::setw(6) << numberOfPartnerPerMonomer[i]; //writes the indices of the molecules and the ammount of partners they posess
        for (std::size_t j = 1; j < (numberOfPartnerPerMonomer[i] + 1); j++)
        {
          kopplung << std::setw(6) << partner[i][j]; //writes the indices of the partners j
          if (partner[i][j] < (numberOf_p_SC))
          {
            kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_rek[i][partner[i][j]]; //writes the recombination coupling between i and j
          }
          else if (partner[i][j] > numberOf_p_SC)
          {
            kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_fulleren[i][partner[i][j]]; // writes some coupling regarding fullerens?
          }
        }
        kopplung << '\n';
      }
      kopplung.close();

      std::ofstream interface;
      interface.open("masspoints_general.xyz"); //writing out average balance points for all groupings of monomers
      interface << "4" << '\n' << '\n';
      interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << avg_position_p_sc__x << std::setw(12) << std::setprecision(6) << std::fixed << avg_position_p_sc__y << std::setw(12) << std::setprecision(6) << std::fixed << avg_position_p_sc__z << '\n';
      interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << avg_position_n_sc__x << std::setw(12) << std::setprecision(6) << std::fixed << avg_position_n_sc__y << std::setw(12) << std::setprecision(6) << std::fixed << avg_position_n_sc__z << '\n';
      interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << avg_position_total__x << std::setw(12) << std::setprecision(6) << std::fixed << avg_position_total__y << std::setw(12) << std::setprecision(6) << std::fixed << avg_position_total__z << '\n';
      interface.close();


      interface.open("startingpoints.xyz");
      std::vector<size_t> startpunkte_tmp;
      std::size_t numberStartpoints = 0u;
      calculateStartingpoints(direction, numberStartpoints, startpunkte_tmp);

      interface << numberStartpoints << '\n' << '\n'; //writes the number of startingponts
      for (std::size_t i = 1; i < (numberStartpoints + 1); i++) 	//writes coordinates of startingpoints
      {
        interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x[startpunkte_tmp[i]];
        interface << std::setw(12) << std::setprecision(6) << std::fixed << y[startpunkte_tmp[i]];
        interface << std::setw(12) << std::setprecision(6) << std::fixed << z[startpunkte_tmp[i]] << '\n';
      }
      interface.close();
    }

    void processAfterFilereading()
    {
      // Startpunkte bestimmen ##################################################################################################################


      for (std::size_t i = 1u; i < (numberOf_p_SC + 1); i++)
      {
        avg_position_p_sc__x += (x[i] / numberOf_p_SC);
        avg_position_p_sc__y += (y[i] / numberOf_p_SC);
        avg_position_p_sc__z += (z[i] / numberOf_p_SC);
      }

      for (std::size_t i = (numberOf_p_SC + 1); i < (totalNumberOfMonomers + 1); i++)     //using fact, that fullerens always have larger indices than other monomers
      {
        avg_position_n_sc__x += (x[i] / numberOf_n_SC);
        avg_position_n_sc__y += (y[i] / numberOf_n_SC);
        avg_position_n_sc__z += (z[i] / numberOf_n_SC);
      }

      this->avg_position_total__x = (avg_position_p_sc__x + avg_position_n_sc__x) / 2.;
      this->avg_position_total__y = (avg_position_p_sc__y + avg_position_n_sc__y) / 2.;
      this->avg_position_total__z = (avg_position_p_sc__z + avg_position_n_sc__z) / 2.;
    }

    void calculateStartingpoints(char direction, std::size_t& numPoints, std::vector <std::size_t>& vecOfStartingPoints, double procentualDist2Interf = 0.85) const
    {
      std::size_t& index = numPoints;
      index = 0u; // Yes, this is correct...
      double max = std::numeric_limits<double>::lowest();
      double min = std::numeric_limits<double>::max();
      std::vector <std::size_t>& returner = vecOfStartingPoints;
      vecOfStartingPoints = std::vector <std::size_t>(this->numberOf_p_SC + 1);

      switch (direction)
      { //different cases for the possible planes of the interface
      case 'x':
        for (std::size_t i = 1u; i < (numberOf_p_SC + 1u); i++) //determining the maximal distance to interface
        {
          if (x[i] > max) {
            max = x[i];
          }
          if (x[i] < min) {
            min = x[i];
          }
        }


        for (std::size_t i = 1u; i < (numberOf_p_SC + 1u); i++)  //determining the necessary number of starting points? 
        {
          double comparison = max;
          if (avg_position_p_sc__x < avg_position_n_sc__x)
            comparison = min;
          if (std::abs(x[i] - avg_position_total__x) > std::abs(procentualDist2Interf*(comparison - avg_position_total__x)))
          {
            index++;
            vecOfStartingPoints[index] = i;
          }
        }
        break;

      case 'y':
        for (std::size_t i = 1u; i < (numberOf_p_SC + 1); i++)  //determining the maximal distance to interface
        {
          if (y[i] > max) {
            max = y[i];
          }
          if (y[i] < min) {
            min = y[i];
          }
        }
        for (std::size_t i = 1; i < (numberOf_p_SC + 1); i++) //determining the necessary number of starting points? 
        {
          double comparison = max;
          if (avg_position_p_sc__y < avg_position_n_sc__y)
            comparison = min;
          if (std::abs(y[i] - avg_position_total__y) > std::abs(procentualDist2Interf*(comparison - avg_position_total__y)))
          {
            index++;
            vecOfStartingPoints[index] = i;
          }
        }
        break;

      case 'z':
        for (std::size_t i = 1u; i < (numberOf_p_SC + 1); i++) //determining the maximal distance to interace
        {
          if (z[i] > max) {
            max = z[i];
          }
          if (z[i] < min) {
            min = z[i];
          }
        }
        for (std::size_t i = 1; i < (numberOf_p_SC + 1); i++) //determining the necessary number of starting points? 
        {
          double comparison = max;
          if (avg_position_p_sc__z < avg_position_n_sc__z)
            comparison = min;
          if (std::abs(z[i] - avg_position_total__z) > std::abs(procentualDist2Interf*(comparison - avg_position_total__z)))
          {
            index++;
            vecOfStartingPoints[index] = i;
          }
        }
        break;
      }
    }

    void run(char direction, double const excitonicDrivingForce_GaussianSigma = 0.0338987, double const chargecarrierDrivingForce_GaussianSigma = 0.068584577) // hier neue standardabweichung eintragen
    {
      calculateStartingpoints(direction, this->numberOfStartingPoints, this->startpunkt);

      // ################################################################################## Beginn der Simulation ##############################################################################
      // Variablen


      // Schrittanzahl pro MC-Simulation
      std::size_t const numberOfSteps = 2 * (numberOfExcitonPairs + numberOfNSemiconductorHomopairs) + 400u; // MAGIC NUMBER?
      std::cout << "Number of steps for MC-Simulation: " << numberOfSteps << "." << std::endl;

      // ###################################################################################

      m_results = results(numberOfStartingPoints, numberOfRunsPerStartingPoint);
      std::vector <std::vector<double>> & vel_ex = m_results.vel_ex;
      std::vector <std::vector<double>> & vel_ch = m_results.vel_ch;
      std::vector <std::vector<double>> & zeit_ex = m_results.zeit_ex;
      std::vector <std::vector<double>> & zeit_ch = m_results.zeit_ch;
      std::vector <std::size_t>& ex_diss = m_results.ex_diss;
      std::vector <std::size_t>& ch_diss = m_results.ch_diss;
      std::vector <std::size_t>& rek = m_results.rek;
      std::vector <std::size_t>& trapping = m_results.trapping;
      std::vector <std::size_t>& radiativ = m_results.radiativ;

      std::vector <std::vector<char>> zustand(numberOfStartingPoints + 1, std::vector <char>(numberOfRunsPerStartingPoint));
      std::vector <std::size_t> punkt(numberOfSteps + 1), punkt_ladung(numberOfSteps + 1);

      for (std::size_t i = 1; i < (numberOfStartingPoints + 1); i++) //initializing the vectors with 0
      {
        for (std::size_t j = 1; j < numberOfRunsPerStartingPoint; j++)
        {
          vel_ex[i][j] = 0u;
          vel_ch[i][j] = 0u;
          zeit_ex[i][j] = 0u;
          zeit_ch[i][j] = 0u;
          zustand[i][j] = 'e';
        }
      }

      // k: index für startpunkte
      // j: index für durchläufe
      // i: index für schritt
      for (std::size_t i = 0u; i < (numberOfSteps + 1); i++)
      {
        punkt[i] = 0u;
      }

      std::ofstream run;
      if (Config::get().general.verbosity >= 4u)
        run.open("debug_xb_log.txt");

      std::random_device rd;
      for (std::size_t k = 1; k < (numberOfStartingPoints + 1); k++) // schleife über startpunkte "index durch 1 vertauscht"
      {
        run << "Startingpoint(k)-Iterator is " << k << std::endl;

        for (std::size_t j = 1; j < numberOfRunsPerStartingPoint; j++)   // schleife über durchläufe für den gleichen startpunkt " 101 durch 11 vertauscht"
        {
          double zeit(0.), zeit_1(0.), zeit_2(0.);

          punkt[0] = startpunkt[k];

          for (std::size_t i = 1; i < (numberOfSteps + 1); i++)
          {
            std::mt19937 engine(rd());
            if (zustand[k][j] == 'c')
            {
              double r_sum = 0.;
              // site energies berechnen
              //########## raten addieren für monomere ##########

              std::normal_distribution<double> distribution0(0.0, chargecarrierDrivingForce_GaussianSigma); 
              const double zufall1 = distribution0(engine); //generating normal-distributed random number

              std::vector<double> raten(numberOfPartnerPerMonomer[punkt_ladung[i - 1]] + 1);
              for (std::size_t h = 0; h < (numberOfPartnerPerMonomer[punkt_ladung[i - 1]] + 1); h++)
              {
                const double zufall = distribution0(engine);
                if (partner[punkt_ladung[i - 1]][h] < (numberOf_p_SC + 1))
                {
                  //const double coulombenergy = coulomb(x, y, z, punkt[i - 1], partner[punkt_ladung[i - 1]][h], 3.4088) - coulomb(x, y, z, punkt[i - 1], punkt_ladung[i - 1], 3.4088);
                  const double coulombenergy = evaluateCoulomb(punkt[i - 1], partner[punkt_ladung[i - 1]][h], 3.4088) - evaluateCoulomb(punkt[i - 1], punkt_ladung[i - 1], 3.4088);
                  r_sum = r_sum + rate(coupling_ladung[punkt_ladung[i - 1]][partner[punkt_ladung[i - 1]][h]], ((zufall - zufall1) + coulombenergy), reorganisationsenergie_ladung);

                  raten[h] = r_sum;
                }
                if ((partner[punkt_ladung[i - 1]][h] > (numberOf_p_SC)) && (partner[punkt_ladung[i - 1]][h] == punkt[i - 1]))
                {
                  // coulomb energie berechnen	   
                  const double coulombenergy = evaluateCoulomb(punkt_ladung[i - 1], partner[punkt_ladung[i - 1]][h], 1);
                  r_sum = r_sum + rate(coupling_rek[punkt_ladung[i - 1]][partner[punkt_ladung[i - 1]][h]], (zufall - zufall1) + rekombinationstriebkraft - coulombenergy, rek_reorganisation);
                  raten[h] = r_sum;
                }
                //ACHTUNG: hier Korrektur ///////////////////////////////////////////////////
                else if ((partner[punkt_ladung[i - 1]][h] > (numberOf_p_SC)) && (partner[punkt_ladung[i - 1]][h] != punkt[i - 1]))
                {
                  r_sum = r_sum;
                  raten[h] = 0;
                }
              }

              // hier raten für fullerene addieren 
              double r_sum_n_sc = 0.;
              const double zufall2 = distribution0(engine);
              std::vector <double> raten_fulleren(numberOfPartnerPerMonomer[punkt[i - 1]] + 1);
              for (std::size_t h = 1; h < (numberOfPartnerPerMonomer[punkt[i - 1]] + 1); h++)
              {
                if (partner[punkt[i - 1]][h] > (numberOf_p_SC))
                {
                  const double zufall = distribution0(engine);
                  const double coulombenergy = evaluateCoulomb(punkt_ladung[i - 1], partner[punkt[i - 1]][h], 3.4088) - evaluateCoulomb(punkt[i - 1], punkt_ladung[i - 1], 3.4088);
                  r_sum_n_sc = r_sum_n_sc + rate(coupling_fulleren[punkt[i - 1]][partner[punkt[i - 1]][h]], ((zufall - zufall2) + coulombenergy), fullerenreorganisationsenergie);

                  raten_fulleren[h] = r_sum_n_sc;
                }
                if ((partner[punkt[i - 1]][h] < (numberOf_p_SC + 1)) && (partner[punkt[i - 1]][h] == punkt_ladung[i - 1]))
                {
                  const double zufall = distribution0(engine);

                  // coulomb energie berechnen
                  const double coulombenergy = evaluateCoulomb(punkt[i - 1], partner[punkt[i - 1]][h], 1);

                  r_sum_n_sc = r_sum_n_sc + rate(coupling_rek[punkt[i - 1]][partner[punkt[i - 1]][h]], (zufall - zufall2) + rekombinationstriebkraft - coulombenergy, rek_reorganisation);

                  raten_fulleren[h] = r_sum_n_sc;
                }
                //ACHTUNG: hier Korrektur ////////////////////////////////////////////////////////
                else if ((partner[punkt[i - 1]][h] < (numberOf_p_SC + 1)) && (partner[punkt[i - 1]][h] != punkt_ladung[i - 1]))
                {
                  r_sum_n_sc = r_sum_n_sc;
                  raten_fulleren[h] = 0;
                }
                //////////////////////////////////////////////////////////////////////////////////
              }

              // hüpfendes teilchen bestimmen
              if ((1 / r_sum - zeit_1) < (1 / r_sum_n_sc - zeit_2))
              {
                run << "P-SC hopps first." << std::endl;

                //Update der Zeiten
                if ((1. / r_sum - zeit_1) > 0.)
                {
                  zeit = zeit + (1. / r_sum - zeit_1);
                  zeit_2 = zeit_2 + (1 / r_sum - zeit_1);
                  zeit_1 = 0.;
                }
                else if ((1. / r_sum - zeit_1) < 0.)
                {
                  zeit = zeit;
                  zeit_2 = zeit_2;
                  zeit_1 = 0.;
                }
                else
                {
                  run << "ERROR!" << std::endl;
                  throw std::runtime_error("Critical Error in Exciton Breakup Task, Aborting!");
                }
                // monomerhüpfen ausführen
                std::uniform_real_distribution<double> distribution1(0, 1);
                const double zufall = distribution1(engine);
                const double r_i = zufall * r_sum;
                for (std::size_t g = 1; g < (numberOfPartnerPerMonomer[punkt_ladung[i - 1]] + 1); g++)
                {
                  if ((raten[g] > r_i) && (partner[punkt_ladung[i - 1]][g] < (numberOf_p_SC + 1)))
                  {
                    run << "Chargetransport" << std::endl;
                    punkt_ladung[i] = partner[punkt_ladung[i - 1]][g];
                    punkt[i] = punkt[i - 1];

                    // Abbruchkriterium für Ladungstrennung
                    constexpr double distanceCriterion = 0.75;
                    switch (direction)
                    {
                    case 'x':
                      if (((x[punkt_ladung[i]]) - avg_position_total__x) > (distanceCriterion*(x[startpunkt[k]] - avg_position_total__x)))
                      {
                        ch_diss[k]++;
                        zeit_ch[k][j] = zeit - zeit_ex[k][j];
                        vel_ch[k][j] = ((x[punkt_ladung[i]]) - avg_position_total__x) / zeit_ch[k][j];

                        run << "Charges separated" << std::endl;
                        zustand[k][j] = 's';
                      }
                      break;
                    case 'y':
                      if (((y[punkt_ladung[i]]) - avg_position_total__y) > (distanceCriterion*(y[startpunkt[k]] - avg_position_total__y)))
                      {
                        ch_diss[k]++;
                        zeit_ch[k][j] = zeit - zeit_ex[k][j];
                        vel_ch[k][j] = ((y[punkt_ladung[i]]) - avg_position_total__y) / zeit_ch[k][j];

                        run << "Charges separated" << std::endl;
                        zustand[k][j] = 's';
                      }
                      break;
                    case 'z':
                      if (((z[punkt_ladung[i]]) - avg_position_total__z) > (distanceCriterion*(z[startpunkt[k]] - avg_position_total__z)))
                      {
                        ch_diss[k]++;
                        zeit_ch[k][j] = zeit - zeit_ex[k][j];
                        vel_ch[k][j] = ((z[punkt_ladung[i]]) - avg_position_total__z) / zeit_ch[k][j];

                        run << "Charges separated" << std::endl;
                        zustand[k][j] = 's';
                      }
                      break;
                    }

                    //#########################################################################################################################
                    run << "old p-SC " << std::setw(5) << punkt_ladung[i - 1] << std::endl;
                    run << "new p-SC " << std::setw(5) << punkt_ladung[i] << std::endl;
                    run << "Coupling " << std::setw(12) << std::setprecision(6) << std::fixed << coupling_ladung[punkt_ladung[i - 1]][punkt_ladung[i]] << std::endl;
                    run << "n-SC " << std::setw(5) << punkt[i] << std::setw(5) << punkt[i - 1] << std::endl;
                    break;
                  }
                  else if ((raten[g] > r_i) && (partner[punkt_ladung[i - 1]][g] > (numberOf_p_SC)))
                  {
                    run << "Recombination" << std::endl;
                    zustand[k][j] = 't';
                    rek[k]++;
                    break;
                  }
                  else if (g == (numberOfPartnerPerMonomer[punkt_ladung[i - 1]]))
                  {
                    run << "WARING: ERROR during p-semiconductor chargetransfer." << std::endl;
                    throw std::runtime_error("WARING: ERROR during p-semiconductor chargetransfer. Aborting");
                  }
                }
              } // endes des hüpfenden p-halbleiters

              else if ((1 / r_sum - zeit_1) > (1 / r_sum_n_sc - zeit_2))
              {
                run << "N-SC hopped first." << std::endl;
                if ((1 / r_sum_n_sc - zeit_2) > 0)
                {
                  zeit = zeit + (1 / r_sum_n_sc - zeit_2);
                  zeit_1 = zeit_1 + (1 / r_sum_n_sc - zeit_2);
                  zeit_2 = 0;
                }
                else if ((1 / r_sum_n_sc - zeit_2) < 0)
                {
                  zeit = zeit;
                  zeit_1 = zeit_1;
                  zeit_2 = 0;
                }
                else
                {
                  run << "ERROR!" << std::endl;
                  throw std::runtime_error("Critical Error in Exciton Breakup Task, Aborting!");
                }

                // fullerenhüpfen ausführen
                std::uniform_real_distribution<double> distribution1(0, 1);
                const double zufall = distribution1(engine);
                const double r_i = zufall * r_sum_n_sc;

                for (std::size_t g = 1u; g < ((numberOfPartnerPerMonomer[punkt[i - 1]]) + 1); g++)
                {
                  if ((raten_fulleren[g] > r_i) && ((partner[punkt[i - 1]][g]) > numberOf_p_SC))
                  {
                    run << "Chargetransfer in n-SC phase" << std::endl;

                    punkt[i] = partner[punkt[i - 1]][g];
                    punkt_ladung[i] = punkt_ladung[i - 1];
                    run << "old n-SC " << std::setw(5) << punkt[i - 1] << std::endl;
                    run << "new n-SC " << std::setw(5) << punkt[i] << std::endl;
                    run << "p-SC " << std::setw(5) << punkt_ladung[i] << std::endl;
                    run << "Coupling " << std::setw(12) << std::setprecision(6) << coupling_fulleren[punkt[i]][punkt[i - 1]] << std::endl;
                    break;
                  }
                  else if ((raten_fulleren[g] > r_i) && ((partner[punkt[i - 1]][g]) < (numberOf_p_SC + 1)))
                  {

                    run << "Recombination." << std::endl;
                    zustand[k][j] = 't';
                    rek[k]++;
                    break;
                  }
                  else if (g == (numberOfPartnerPerMonomer[punkt[i - 1]]))
                  {
                    run << "WARNING: ERROR during fullerene chargetransport." << std::endl;
                    throw std::runtime_error("Critical Error in Exciton Breakup Task, Aborting!");
                  }
                }
              } // endes des hüpfenden fullerens
            } // ende des 'c'-zustands
          //__________________________________________________________________________________________________________

            else if (zustand[k][j] == 'e')
            {
              // site energies berechnen
              double r_summe = 0.;
              std::normal_distribution<double> distribution0(0.0, excitonicDrivingForce_GaussianSigma); 
              const double zufall1 = distribution0(engine); //generating an normal-distributed random number

              std::vector <double> raten(numberOfPartnerPerMonomer[punkt[i - 1]] + 1);

              for (std::size_t h = 1; h < (numberOfPartnerPerMonomer[punkt[i - 1]] + 1); h++)
              {
                // Jump to p SC
                if (partner[punkt[i - 1]][h] < (numberOf_p_SC + 1))
                {
                  const double zufall = distribution0(engine);// generatinjg a second normal distributed random number
                  const double testrate = rate(coupling_exciton[punkt[i - 1]][partner[punkt[i - 1]][h]], (zufall - zufall1), reorganisationsenergie_exciton);
                  r_summe += testrate;
                  raten[h] = r_summe;
                  run << "A: " << punkt[i - 1] << "   B: " << partner[punkt[i - 1]][h] << "rate   " << testrate << std::endl;
                }
                // Jump tp n SC
                else if (partner[punkt[i - 1]][h] > (numberOf_p_SC))
                {
                  const double zufall = distribution0(engine);
                  // coulomb energie berechnen
                  const double coulombenergy = evaluateCoulomb(punkt[i - 1], partner[punkt[i - 1]][h], 1);
                  const double testrate = rate(coupling_ct[punkt[i - 1]][partner[punkt[i - 1]][h]], (zufall - zufall1) + chargetransfertriebkraft + coulombenergy, ct_reorganisation);
                  r_summe += testrate;
                  raten[h] = r_summe;
                  run << "coulomb  " << coulombenergy << "  rate   " << testrate << std::endl;
                }
              } // end of h

              // fluoreszenz dazuaddieren
              r_summe = r_summe + k_rad;

              // schritt bestimmen
              std::uniform_real_distribution<double> distribution1(0, 1);
              double zufall = distribution1(engine);
              const double r_i = zufall * r_summe;
              zeit = zeit + 1 / r_summe;

              //falls trapping
              zufall = distribution1(engine);
              if (zufall*(900e-1 + 1 / r_summe) > (900e-1))
              {
                run << "Exciton trapped!" << std::endl;
                trapping[k]++;
                zustand[k][j] = 't';
                break;
              }

              for (std::size_t g = 1u; g < (numberOfPartnerPerMonomer[punkt[i - 1]] + 1); g++)
              {
                if (raten[g] > r_i)
                {
                  punkt[i] = partner[punkt[i - 1]][g];

                  if (punkt[i] < (numberOf_p_SC + 1))
                  {
                    run << "hopped to " << punkt[i] << std::endl;
                  }
                  else if (punkt[i] > numberOf_p_SC)
                  {
                    punkt[i] = partner[punkt[i - 1]][g];
                    punkt_ladung[i] = punkt[i - 1];
                    zustand[k][j] = 'c';
                    run << "Chargeseparation." << std::endl;

                    vel_ex[k][j] = distance(punkt[0], punkt[i]) / zeit;
                    //run << "Exzitonspeed " << vel_ex[k][j] * 1e-9 << std::endl;
                    ex_diss[k]++;
                  }

                  break;
                }
                else if (raten[numberOfPartnerPerMonomer[punkt[i - 1]]] < r_i)
                {
                  run << "radiating decay." << std::endl;
                  radiativ[k]++;
                  zustand[k][j] = 't';
                  break;
                }
              }
            } // end of 'e'-zustand
          //____________________________________________________________________________________________

            else if (zustand[k][j] == 't')
            {
              run << "BROKEN!" << std::endl;
              break;
            } // end of 't'-zustand
          //__________________________________________________________________________________________________

            else if (zustand[k][j] == 's')
            {
              run << "SUCCESS!" << std::endl;
              break;
            } //end of 's'-zustand 
          //_________________________________________________________________________________________________

            else
            {
              run << "Warning. State undefined!" << std::endl;
              throw std::runtime_error("Warning: Undefined State encountered. Aborting.");
            } //end of undefined zustand
          //___________________________________________________________________________________________________

          } // ende über schleife i
        } // Ende über Schleife über durchläufe für den gleichen startpunkt j
      } // Ende über Schleife der Startpunkte k

      run.close();

    }

    void analyseResults() const
    {
      // Auswertung: Prozentsätze
      std::ofstream auswertung;
      auswertung.open("evaluation.txt");
      auswertung << std::setw(4) << "k" << std::setw(4) << "IX" << std::setw(9) << "Ex_Diss." << std::setw(9) << "Ch_Diss" << std::setw(9) << "Rek." << std::setw(9) << "Trapp." << std::setw(9) << "Fluor." << '\n';

      double mittel_ex = 0;
      double mittel_ch = 0;
      double mittel_rek = 0;
      double mittel_trapp = 0;
      double mittel_rad = 0;
      std::vector <double> mittel_ex_vel(numberOfStartingPoints + 1), mittel_ch_vel(numberOfStartingPoints + 1), standard_ex(numberOfStartingPoints + 1), standard_ch(numberOfStartingPoints + 1);
      std::vector <std::size_t> const& ex_diss = m_results.ex_diss;
      std::vector <std::size_t> const& ch_diss = m_results.ch_diss;
      std::vector <std::size_t> const& rek = m_results.rek;
      std::vector <std::size_t> const& trapping = m_results.trapping;
      std::vector <std::size_t> const& radiativ = m_results.radiativ;

      std::vector <std::vector<double>> const & vel_ex = m_results.vel_ex;
      std::vector <std::vector<double>> const & vel_ch = m_results.vel_ch;
      std::vector <std::vector<double>> const & zeit_ex = m_results.zeit_ex;
      std::vector <std::vector<double>> const & zeit_ch = m_results.zeit_ch;

      for (std::size_t k = 1u; k < (numberOfStartingPoints + 1); k++)
      {
        mittel_ex = mittel_ex + (ex_diss[k] * 1.0);
        mittel_ch = mittel_ch + (ch_diss[k] * 1.0);
        mittel_rek = mittel_rek + (rek[k] * 1.0);
        mittel_trapp = mittel_trapp + (trapping[k] * 1.0);
        mittel_rad = mittel_rad + (radiativ[k] * 1.0);
        double const ex_diss_efficiency = ex_diss[k];
        double const ch_diss_efficiency = ch_diss[k];
        double const rek_efficiency = rek[k];
        double const trapp_efficiency = trapping[k];
        double const rad_efficiency = radiativ[k];
        auswertung << std::setw(4) << k << std::setw(4) << startpunkt[k] << std::setw(9) << std::setprecision(5) << ex_diss_efficiency;
        auswertung << std::setw(9) << std::setprecision(5) << ch_diss_efficiency;
        auswertung << std::setw(9) << std::setprecision(5) << rek_efficiency << std::setw(9) << std::setprecision(5) << trapp_efficiency << std::setw(9) << std::setprecision(5) << rad_efficiency << '\n';
      }

      auswertung << std::setw(9) << "Average " << std::setw(9) << std::setprecision(5) << std::fixed << mittel_ex / numberOfStartingPoints << std::setw(9) << std::setprecision(5) << std::fixed << mittel_ch / numberOfStartingPoints;
      auswertung << std::setw(9) << std::setprecision(5) << std::fixed << mittel_rek / numberOfStartingPoints << std::setw(9) << std::setprecision(5) << std::fixed << mittel_trapp / numberOfStartingPoints;
      auswertung << std::setw(9) << std::setprecision(5) << std::fixed << mittel_rad / numberOfStartingPoints << '\n';
      auswertung << "Velocities" << '\n';
      auswertung << std::setw(4) << "k" << std::setw(5) << "IX" << std::setw(11) << "Ex_vel" << std::setw(11) << "Ex_s_dev" << std::setw(11) << "Ch_vel" << std::setw(11) << "Ch_s_dev" << '\n';

      // mittlere Geschwindigkeit
      for (std::size_t k = 1; k < (numberOfStartingPoints + 1); k++)
      {
        mittel_ch_vel[k] = 0;
        mittel_ex_vel[k] = 0;

        // Mittelwert berechnen
        for (std::size_t j = 1; j < numberOfRunsPerStartingPoint; j++)
        {
          if (vel_ch[k][j] > 0.0001) {
            mittel_ch_vel[k] = mittel_ch_vel[k] + vel_ch[k][j];
          }
          if (vel_ex[k][j] > 0.0001) {
            mittel_ex_vel[k] = mittel_ex_vel[k] + vel_ex[k][j];
          }
        }
        if (ch_diss[k] > 0) {
          mittel_ch_vel[k] = mittel_ch_vel[k] / (ch_diss[k] * 1.0);
        }
        else if (ch_diss[k] == 0) {
          mittel_ch_vel[k] = 0;
        }
        if (ex_diss[k] > 0) {
          mittel_ex_vel[k] = mittel_ex_vel[k] / (ex_diss[k] * 1.0);
        }
        else if (ex_diss[k] == 0) {
          mittel_ex_vel[k] = 0;
        }

        // Standardabweichunb berechnen
        standard_ex[k] = 0;
        standard_ch[k] = 0;
        for (std::size_t j = 1; j < numberOfRunsPerStartingPoint; j++)
        {
          if (vel_ch[k][j] > 0.0001) {
            standard_ch[k] = standard_ch[k] + (vel_ch[k][j] - mittel_ch_vel[k])*(vel_ch[k][j] - mittel_ch_vel[k]);
          }
          if (vel_ex[k][j] > 0.0001) {
            standard_ex[k] = standard_ex[k] + (vel_ex[k][j] - mittel_ex_vel[k])*(vel_ex[k][j] - mittel_ex_vel[k]);
          }
        }
        if (ch_diss[k] > 1) {
          standard_ch[k] = sqrt(standard_ch[k] / (ch_diss[k] - 1));
        }
        else if (ch_diss[k] < 2) {
          standard_ch[k] = 0;
        }
        if (ex_diss[k] > 1) {
          standard_ex[k] = sqrt(standard_ex[k] / (ex_diss[k] - 1));
        }
        else if (ex_diss[k] < 2) {
          standard_ex[k] = 0;
        }
        auswertung << std::setw(4) << k << std::setw(5) << startpunkt[k] << std::setw(11) << std::setprecision(5) << std::fixed << mittel_ex_vel[k] * 1e-9;
        auswertung << std::setw(11) << std::setprecision(5) << std::fixed << standard_ex[k] * 1e-9;
        auswertung << std::setw(11) << std::setprecision(5) << std::fixed << mittel_ch_vel[k] * 1e-9;
        auswertung << std::setw(11) << std::setprecision(5) << std::fixed << standard_ch[k] * 1e-9 << '\n';
      }
      double mittelwert_geschw_exciton = 0;
      double mittelwert_geschw_ladung = 0;
      for (std::size_t k = 1; k < (numberOfStartingPoints + 1); k++)
      {
        mittelwert_geschw_exciton = mittelwert_geschw_exciton + mittel_ex_vel[k];
        mittelwert_geschw_ladung = mittelwert_geschw_ladung + mittel_ch_vel[k];
      }
      auswertung << std::left << std::setw(7) << " Average    " << std::left << std::setw(22) << std::setprecision(5) << std::fixed << mittelwert_geschw_exciton / numberOfStartingPoints * 1e-9;
      auswertung << std::left << std::setw(9) << std::setprecision(5) << std::fixed << mittelwert_geschw_ladung / numberOfStartingPoints * 1e-9 << '\n';

      // Verteilung Ladungen und Exzitonengeschwindigkeiten
      std::ofstream exciton_verteilung;
      exciton_verteilung.open("exciton_distribution.txt");
      for (std::size_t i = 1; i < 21; i++) {
        std::size_t zahl = 0u;
        for (std::size_t k = 1; k < (numberOfStartingPoints + 1); k++) {
          for (std::size_t j = 1; j < numberOfRunsPerStartingPoint; j++) {
            if ((vel_ex[k][j] > (i * 50 * 1e9))) {
              zahl++;
            }
          }
        }
        exciton_verteilung << std::setw(9) << std::setprecision(5) << i * 50 << std::setw(9) << zahl / numberOfStartingPoints << '\n';
      }

      exciton_verteilung.close();
      exciton_verteilung.open("charge_distribution.txt");
      for (std::size_t i = 1; i < 21; i++) {
        std::size_t zahl = 0u;
        for (std::size_t k = 1; k < (numberOfStartingPoints + 1); k++) {
          for (std::size_t j = 1; j < numberOfRunsPerStartingPoint; j++) {
            if ((vel_ch[k][j] > (i * 50 * 1e9))) {
              zahl++;
            }
          }
        }
        exciton_verteilung << std::setw(9) << std::setprecision(5) << i * 50 << std::setw(9) << zahl / numberOfStartingPoints << '\n';
      }
      exciton_verteilung.close();
    }

    double evaluateCoulomb(std::size_t particle1_iterator, std::size_t particle2_iterator, double e_relative) const
    {
      const double l = this->distance(particle1_iterator, particle2_iterator);
      constexpr double pi = constants::pi;
      constexpr double epsilon_0 = constants::epsilon_0;
      constexpr double elementar = 1.60217662e-19;
      const double c = -elementar / (4. * pi*epsilon_0*e_relative*l*1e-10);
      return c;
    }

    // Definition der length-berechnung
    inline double distance(std::size_t particle1_iterator, std::size_t particle2_iterator) const
    {
      std::size_t const& p = particle1_iterator;
      std::size_t const& q = particle2_iterator;
      std::vector <double> const& arr1 = x;
      std::vector <double> const& arr2 = y;
      std::vector <double> const& arr3 = z;
      const double l = std::sqrt((arr1[p] - arr1[q])*(arr1[p] - arr1[q]) + (arr2[p] - arr2[q])*(arr2[p] - arr2[q]) + (arr3[p] - arr3[q])*(arr3[p] - arr3[q]));
      return l;
    }

    std::size_t totalNumberOfMonomers;
    std::size_t numberOf_p_SC, numberOf_n_SC;
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

    double avg_position_total__x, avg_position_total__y, avg_position_total__z;
    std::size_t numberOfStartingPoints; // Number of starting points
    double avg_position_p_sc__x, avg_position_p_sc__y, avg_position_p_sc__z, avg_position_n_sc__x, avg_position_n_sc__y, avg_position_n_sc__z;



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
    const std::size_t numberOfRunsPerStartingPoint;
  };

}

