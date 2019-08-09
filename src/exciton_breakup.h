#pragma once
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<limits>
#include<string>
#include<random>
#include<vector>

#include "helperfunctions.h"

//Original code by Charlotte. German comments are from original code, english comments were inserted during implementation into CAST.
//Arrays were replaced by vectors.



namespace XB
{

  // Definition der length-berechnung
  double length(std::vector <double> const& arr1, std::vector <double> const& arr2, std::vector <double> const& arr3, int p, int q)
  {
    const double l = sqrt((arr1[p] - arr1[q])*(arr1[p] - arr1[q]) + (arr2[p] - arr2[q])*(arr2[p] - arr2[q]) + (arr3[p] - arr3[q])*(arr3[p] - arr3[q]));
    return l;
  }

  double rate(double coupling, double energy, double reorganisation)
  {
    constexpr double pi = 3.141592654;
    constexpr double h_quer = 1 / (2 * pi)*4.135667662e-15;
    constexpr double boltzmann_konstante = 8.6173303e-5; //  in gauß einheiten // Dustin July19: is in eV/K
    const double l = (coupling*coupling) / h_quer * sqrt(pi / (reorganisation*boltzmann_konstante * 298))*exp(-(reorganisation + energy)*(reorganisation + energy) / (4 * boltzmann_konstante * 298 * reorganisation));
    return l;
  }

  double coulomb(std::vector<double> const& arr1, std::vector<double> const& arr2 , std::vector<double> const& arr3, int p, int q, double e_relative)
  {
    const double l = sqrt((arr1[p] - arr1[q])*(arr1[p] - arr1[q]) + (arr2[p] - arr2[q])*(arr2[p] - arr2[q]) + (arr3[p] - arr3[q])*(arr3[p] - arr3[q]));
    constexpr double pi = 3.141592654;
    constexpr double e_0 = 8.854187e-12;
    constexpr double elementar = 1.60217662e-19;
    const double c = -elementar / (4 * pi*e_0*e_relative*l*1e-10);
    return c;
  }

  class ExcitonBreakup
  {
  public:
    //ExcitonBreakup(std::string filename);
    //ExcitonBreakup();
    ExcitonBreakup(std::string masscenters, std::string nscpairrates, std::string pscpairexrates, std::string pscpairchrates, std::string pnscpairrates)//LEGACY
      : gesamtanzahl(0u), reorganisationsenergie_exciton(Config::get().exbreak.ReorgE_exc), reorganisationsenergie_ladung(Config::get().exbreak.ReorgE_ch),
        fullerenreorganisationsenergie(Config::get().exbreak.ReorgE_nSC), ct_reorganisation(Config::get().exbreak.ReorgE_ct), chargetransfertriebkraft(Config::get().exbreak.ct_triebkraft),
        rekombinationstriebkraft(Config::get().exbreak.rek_triebkraft), rek_reorganisation(Config::get().exbreak.ReorgE_rek), oszillatorstrength(Config::get().exbreak.oscillatorstrength),
        wellenzahl(Config::get().exbreak.wellenzahl), k_rad(wellenzahl * wellenzahl*oszillatorstrength), procentualDist2Interf(0.5), vectorlen(101u),
        x_mittel(0.), y_mittel(0.), z_mittel(0.)
    {
      this->read(Config::get().exbreak.pscnumber, Config::get().exbreak.nscnumber, masscenters, nscpairrates, pscpairexrates, pscpairchrates, pnscpairrates);
    };
    
    void runAndWrite(char direction)
    {
      this->processAndWriteAuxFiles(direction);
      this->run(direction);
      this->analyseResults();
    }
  private:
    //void read(std::string filename);
    void read(std::size_t pscanzahl, std::size_t nscanzahl, std::string masscenters, std::string nscpairrates, std::string pscpairexrates, std::string pscpairchrates, std::string pnscpairrates)
    {
      /////////////////////////////////// INPUT-READING
      std::ifstream schwerpunkt;
      schwerpunkt.open(masscenters);
      
      schwerpunkt >> gesamtanzahl;
      std::cout << "Read number of Monomers: " << gesamtanzahl << std::endl;
      if (gesamtanzahl == pscanzahl + nscanzahl) //test if correct number of molecules was given
      {
        std::cout << "Number of monomers is correct, proceeding." << std::endl;
      }
      else //gesamtanzahl != pscanzahl + nscanzahl
      {
        std::cout << "Wrong number of monomers detected!" << std::endl;
        throw std::logic_error("Wrong numbers of p- and n-type molecules in inputfile. Expected: " + std::to_string(gesamtanzahl) + " | Is: " + std::to_string(pscanzahl + nscanzahl));
      }
      schwerpunkt >> skipline;
      x = std::vector <double>(gesamtanzahl + 1);
      y = std::vector <double>(gesamtanzahl + 1); 
      z = std::vector <double> (gesamtanzahl + 1);

      for (std::size_t i = 1u; i < (gesamtanzahl + 1u); i++)
      {
        std::string zeile;
        schwerpunkt >> zeile >> x[i] >> y[i] >> z[i]; //reading and saving balance points for molecules
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

      coupling_exciton = std::vector <std::vector<double>>(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1)); //in original code 2d-arrays were used
      coupling_ladung = std::vector <std::vector<double>>(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1));
      coupling_ct = std::vector <std::vector<double>>(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1));
      coupling_rek = std::vector <std::vector<double>>(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1));
      coupling_fulleren = std::vector <std::vector<double>>(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1));

      for (std::size_t i = 1u; i < (gesamtanzahl + 1u); i++) //initializaton of the 2d vectors with 0 in all places
      {
        for (std::size_t j = 1u; j < (gesamtanzahl + 1u); j++) {
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
      partneranzahl = std::vector <size_t>(pscanzahl + nscanzahl + 1);
      for (std::size_t i = 1u; i < (pscanzahl + nscanzahl + 1); i++)
      {
        partneranzahl[i] = 0u;
      } //inizialisation of all elements with 0

      for (std::size_t i = 1u; i < (numberOfExcitonPairs + 1); i++) // counting homodimer partners for j
      {
        for (std::size_t j = 1u; j < (pscanzahl + 1); j++)
        {
          if ((exciton_1[i] == j) || (exciton_2[i] == j))
          {
            partneranzahl[j]++;
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
        for (std::size_t j = 1; j < (gesamtanzahl + 1); j++)
        {
          if ((hetero_1[i] == j) || (hetero_2[i] == j))
          {
            partneranzahl[j]++;
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
        for (std::size_t j = 1; j < (gesamtanzahl + 1); j++)
        {
          if ((fulleren_1[i] == j) || (fulleren_2[i] == j))
          {
            partneranzahl[j]++;
          }
        }
      }

      partner = std::vector<std::vector<std::size_t>>(gesamtanzahl + 1, std::vector<std::size_t>());//2D-vector with variing length for second vector

      for (std::size_t i = 1; i < (gesamtanzahl + 1); i++) //dynamic allocation for length of 2nd vector
      {
        partner[i].resize(partneranzahl[i] + 1);
      }

      for (std::size_t i = 1; i < (gesamtanzahl + 1); i++) //initializing all elements of vector partner with 0
      {
        for (std::size_t j = 1; j < (partneranzahl[i] + 1); j++) {
          partner[i][j] = 0;
        }
      }

      for (std::size_t i = 1; i < (gesamtanzahl + 1); i++)
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
            j++; // since the 2nd dimension length of vector partner was set to partneranzahl[i] in thes logic construction j must always end up to be equal to partneranzahl[i]
          }
        }
        if (partneranzahl[i] != j - 1) //Sanity Check
        {
          std::cout << "Error with number of partners for monomer " << i << std::endl;
          throw std::runtime_error("Error with number of partners for monomer " + std::to_string(i) + ". Aborting.");
        }
      }

      // INPUT-END
    }
    //void write() const;

    void processAndWriteAuxFiles(char direction)
    {
      std::ofstream kopplung;
      kopplung.open("partner.txt"); // Writing file partner.txt
      for (std::size_t i = 1u; i < (gesamtanzahl + 1); i++)
      {
        kopplung << std::setw(6) << i << std::setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess
        for (std::size_t j = 1; j < (partneranzahl[i] + 1); j++)
        {
          kopplung << std::setw(6) << partner[i][j]; //writes the indices of the partners j
        }
        kopplung << '\n';
      }
      kopplung.close();

      kopplung.open("couplings.txt"); // Writing file couplings.txt
      for (std::size_t i = 1u; i < (pscanzahl + 1); i++)
      {
        kopplung << std::setw(6) << i << std::setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess

        for (std::size_t j = 1; j < (partneranzahl[i] + 1); j++)
        {
          kopplung << std::setw(6) << partner[i][j]; //writes the indices of the partners j
          if (partner[i][j] < (pscanzahl + 1))
          {
            kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_exciton[i][partner[i][j]]; //writes the exciton-coupling between i and j
          }
          else if (partner[i][j] > pscanzahl)
          {
            kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_ct[i][partner[i][j]];     //writes charge-transfer-coupling between i and j 
          }
        }
        kopplung << '\n';
      }

      for (std::size_t i = (pscanzahl + 1); i < (gesamtanzahl + 1); i++)
      {
        kopplung << std::setw(6) << i << std::setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess
        for (std::size_t j = 1; j < (partneranzahl[i] + 1); j++)
        {
          kopplung << std::setw(6) << partner[i][j]; //writes the indices of the partners j
          if (partner[i][j] < (pscanzahl))
          {
            kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_rek[i][partner[i][j]]; //writes the recombination coupling between i and j
          }
          else if (partner[i][j] > pscanzahl)
          {
            kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_fulleren[i][partner[i][j]]; // writes some coupling regarding fullerens?
          }
        }
        kopplung << '\n';
      }
      kopplung.close();

      // Startpunkte bestimmen ##################################################################################################################
      double x_monomer(0.), y_monomer(0.), z_monomer(0.), x_fulleren(0.), y_fulleren(0.), z_fulleren(0.),
        x_gesamt(0.), y_gesamt(0.), z_gesamt(0.);

      for (std::size_t i = 1u; i < (pscanzahl + 1); i++)
      {
        x_monomer += (x[i] / pscanzahl);
        y_monomer += (y[i] / pscanzahl);
        z_monomer += (z[i] / pscanzahl);
      }

      for (std::size_t i = (pscanzahl + 1); i < (gesamtanzahl + 1); i++)     //using fact, that fullerens always have larger indices than other monomers
      {
        x_fulleren += (x[i] / nscanzahl);
        y_fulleren += (y[i] / nscanzahl);
        z_fulleren += (z[i] / nscanzahl);
      }

      for (std::size_t i = 1u; i < (gesamtanzahl + 1); i++)
      {
        x_gesamt += (x[i] / gesamtanzahl);
        y_gesamt += (y[i] / gesamtanzahl);
        z_gesamt += (z[i] / gesamtanzahl);
      }

      this->x_mittel = (x_monomer + x_fulleren) / 2;
      this->y_mittel = (y_monomer + y_fulleren) / 2;
      this->z_mittel = (z_monomer + z_fulleren) / 2;

      std::ofstream interface;
      interface.open("massponts_general.xyz"); //writing out average balance points for all groupings of monomers
      interface << "4" << '\n' << '\n';
      interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x_monomer << std::setw(12) << std::setprecision(6) << std::fixed << y_monomer << std::setw(12) << std::setprecision(6) << std::fixed << z_monomer << '\n';
      interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x_fulleren << std::setw(12) << std::setprecision(6) << std::fixed << y_fulleren << std::setw(12) << std::setprecision(6) << std::fixed << z_fulleren << '\n';
      interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x_gesamt << std::setw(12) << std::setprecision(6) << std::fixed << y_gesamt << std::setw(12) << std::setprecision(6) << std::fixed << z_gesamt << '\n';
      interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x_mittel << std::setw(12) << std::setprecision(6) << std::fixed << y_mittel << std::setw(12) << std::setprecision(6) << std::fixed << z_mittel << '\n';
      interface.close();


      index = 0.;
      double max = 0.;
      startpunkt = std::vector <std::size_t>(pscanzahl + 1);

      switch (direction)
      { //different cases for the possible planes of the interface
      case 'x':
        for (std::size_t i = 1u; i < (pscanzahl + 1); i++) //determining the maximal distance to interface
        {
          if (x[i] > max) {
            max = x[i];
          }
        }
        std::cout << "Maxdistance is " << std::setw(12) << std::setprecision(6) << std::fixed << max << std::endl;

        for (std::size_t i = 1; i < (pscanzahl + 1); i++)  //determining the necessary number of starting points? 
        {
          if ((x[i] - x_mittel) > (procentualDist2Interf*(max - x_mittel))) {
            index++;
            startpunkt[index] = i;
          }
        }
        break;

      case 'y':
        for (std::size_t i = 1u; i < (pscanzahl + 1); i++)  //determining the maximal distance to interface
        {
          if (y[i] > max) {
            max = y[i];
          }
        }
        std::cout << "Maxdistance is " << std::setw(12) << std::setprecision(6) << std::fixed << max << std::endl;

        for (std::size_t i = 1; i < (pscanzahl + 1); i++) //determining the necessary number of starting points? 
        {
          if ((y[i] - y_mittel) > (procentualDist2Interf*(max - y_mittel))) {
            index++;
            startpunkt[index] = i;

          }
        }
        break;

      case 'z':
        for (std::size_t i = 1u; i < (pscanzahl + 1); i++) //determining the maximal distance to interace
        {
          if (z[i] > max) {
            max = z[i];
          }
        }
        std::cout << "Maxdistance is " << std::setw(12) << std::setprecision(6) << std::fixed << max << std::endl;

        for (std::size_t i = 1; i < (pscanzahl + 1); i++) //determining the necessary number of starting points? 
        {
          if ((z[i] - z_mittel) > (procentualDist2Interf*(max - z_mittel)))
          {
            index++;
            startpunkt[index] = i;

          }
        }
        break;
      }

      interface.open("startingpoint.xyz");
      interface << index << '\n' << '\n'; //writes the number of startingponts
      for (std::size_t i = 1; i < (index + 1); i++) 	//writes coordinates of startingpoints
      {
        interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x[startpunkt[i]];
        interface << std::setw(12) << std::setprecision(6) << std::fixed << y[startpunkt[i]];
        interface << std::setw(12) << std::setprecision(6) << std::fixed << z[startpunkt[i]] << '\n';
      }
      interface.close();

    }

    void run(char direction)
    {

      this->processAndWriteAuxFiles(direction);

      // ################################################################################## Beginn der Simulation ##############################################################################
      // Variablen


      // Schrittanzahl pro MC-Simulation
      std::size_t const schritt = 2 * (numberOfExcitonPairs + numberOfNSemiconductorHomopairs) + 400u; // MAGIC NUMBER?
      std::cout << "Number of stepps for MC-Simulation: " << schritt << "." << std::endl;

      // ###################################################################################
      

      std::vector <std::vector<double>> & vel_ex  = m_results.vel_ex;
      std::vector <std::vector<double>> & vel_ch = m_results.vel_ch;
      std::vector <std::vector<double>> & zeit_ex = m_results.zeit_ex;
      std::vector <std::vector<double>> & zeit_ch = m_results.zeit_ch;


      m_results = results(index, vectorlen);
      std::vector <std::size_t>& ex_diss = m_results.ex_diss;
      std::vector <std::size_t>& ch_diss = m_results.ch_diss;
      std::vector <std::size_t>& rek = m_results.rek;
      std::vector <std::size_t>& trapping = m_results.trapping;
      std::vector <std::size_t>& radiativ = m_results.radiativ;

      std::vector <std::vector<char>> zustand(index + 1, std::vector <char>(vectorlen));
      std::vector <int> punkt(schritt + 1), punkt_ladung(schritt + 1);

      for (std::size_t i = 1; i < (index + 1); i++) //initializing the vectors with 0
      {
        for (std::size_t j = 1; j < vectorlen; j++)
        {
          vel_ex[i][j] = 0;
          vel_ch[i][j] = 0;
          zeit_ex[i][j] = 0;
          zeit_ch[i][j] = 0;
          zustand[i][j] = 'e';
        }
      }

      // k: index für startpunkte
      // j: index für durchläufe
      // i: index für schritt
      for (std::size_t i = 0u; i < (schritt + 1); i++)
      {
        punkt[i] = 0;
      }

      std::ofstream run;
      if (Config::get().general.verbosity >= 4u)
        run.open("debug_xb_log.txt");

      for (std::size_t k = 1; k < (index + 1); k++) // schleife über startpunkte "index durch 1 vertauscht"
      {
        run << "k ist " << k << std::endl;

        for (std::size_t j = 1; j < vectorlen; j++)   // schleife über durchläufe für den gleichen startpunkt " 101 durch 11 vertauscht"
        {
          double zeit(0.), zeit_1(0.), zeit_2(0.);

          punkt[0] = startpunkt[k];

          for (std::size_t i = 1; i < (schritt + 1); i++)
          {

            if (zustand[k][j] == 'c')
            {
              double r_summe = 0.;
              // site energies berechnen
              //########## raten addieren für monomere ##########

              std::random_device rd;
              std::default_random_engine engine(rd());
              std::normal_distribution<double> distribution0(0.0, 0.068584577); // TODO: hier neue standardabweichung eintragen
              double zufall1 = distribution0(engine); //generating normal-distributed random number

              std::vector<double> raten(partneranzahl[punkt_ladung[i - 1]] + 1);
              for (std::size_t h = 0; h < (partneranzahl[punkt_ladung[i - 1]] + 1); h++)
              {
                const double zufall = distribution0(engine);
                if (partner[punkt_ladung[i - 1]][h] < (pscanzahl + 1))
                {
                  const double coulombenergy = coulomb(x, y, z, punkt[i - 1], partner[punkt_ladung[i - 1]][h], 3.4088) - coulomb(x, y, z, punkt[i - 1], punkt_ladung[i - 1], 3.4088);
                  r_summe = r_summe + rate(coupling_ladung[punkt_ladung[i - 1]][partner[punkt_ladung[i - 1]][h]], ((zufall - zufall1) + coulombenergy), reorganisationsenergie_ladung);

                  raten[h] = r_summe;
                }
                if ((partner[punkt_ladung[i - 1]][h] > (pscanzahl)) && (partner[punkt_ladung[i - 1]][h] == punkt[i - 1]))
                {
                  // coulomb energie berechnen	   
                  const double coulombenergy = coulomb(x, y, z, punkt_ladung[i - 1], partner[punkt_ladung[i - 1]][h], 1);
                  r_summe = r_summe + rate(coupling_rek[punkt_ladung[i - 1]][partner[punkt_ladung[i - 1]][h]], (zufall - zufall1) + rekombinationstriebkraft - coulombenergy, rek_reorganisation);
                  raten[h] = r_summe;
                }
                //ACHTUNG: hier Korrektur ///////////////////////////////////////////////////
                else if ((partner[punkt_ladung[i - 1]][h] > (pscanzahl)) && (partner[punkt_ladung[i - 1]][h] != punkt[i - 1]))
                {
                  r_summe = r_summe;
                  raten[h] = 0;
                }
              }

              // hier raten für fullerene addieren 
              double r_summe_fulleren = 0.;
              zufall1 = distribution0(engine);
              std::vector <double> raten_fulleren(partneranzahl[punkt[i - 1]] + 1);
              for (std::size_t h = 1; h < (partneranzahl[punkt[i - 1]] + 1); h++)
              {
                if (partner[punkt[i - 1]][h] > (pscanzahl))
                {
                  const double zufall = distribution0(engine);
                  const double coulombenergy = coulomb(x, y, z, punkt_ladung[i - 1], partner[punkt[i - 1]][h], 3.4088) - coulomb(x, y, z, punkt[i - 1], punkt_ladung[i - 1], 3.4088);
                  r_summe_fulleren = r_summe_fulleren + rate(coupling_fulleren[punkt[i - 1]][partner[punkt[i - 1]][h]], ((zufall - zufall1) + coulombenergy), fullerenreorganisationsenergie);

                  raten_fulleren[h] = r_summe_fulleren;
                }
                if ((partner[punkt[i - 1]][h] < (pscanzahl + 1)) && (partner[punkt[i - 1]][h] == punkt_ladung[i - 1]))
                {
                  const double zufall = distribution0(engine);

                  // coulomb energie berechnen
                  const double coulombenergy = coulomb(x, y, z, punkt[i - 1], partner[punkt[i - 1]][h], 1);

                  r_summe_fulleren = r_summe_fulleren + rate(coupling_rek[punkt[i - 1]][partner[punkt[i - 1]][h]], (zufall - zufall1) + rekombinationstriebkraft - coulombenergy, rek_reorganisation);

                  raten_fulleren[h] = r_summe_fulleren;
                }
                //ACHTUNG: hier Korrektur ////////////////////////////////////////////////////////
                else if ((partner[punkt[i - 1]][h] < (pscanzahl + 1)) && (partner[punkt[i - 1]][h] != punkt_ladung[i - 1]))
                {
                  r_summe_fulleren = r_summe_fulleren;
                  raten_fulleren[h] = 0;
                }
                //////////////////////////////////////////////////////////////////////////////////
              }

              // hüpfendes teilchen bestimmen
              if ((1 / r_summe - zeit_1) < (1 / r_summe_fulleren - zeit_2))
              {
                run << "Monomer hopps first." << std::endl;

                //Update der Zeiten
                if ((1. / r_summe - zeit_1) > 0.)
                {
                  zeit = zeit + (1. / r_summe - zeit_1);
                  zeit_2 = zeit_2 + (1 / r_summe - zeit_1);
                  zeit_1 = 0.;
                }
                else if ((1. / r_summe - zeit_1) < 0.)
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
                const double r_i = zufall * r_summe;
                for (std::size_t g = 1; g < (partneranzahl[punkt_ladung[i - 1]] + 1); g++)
                {
                  if ((raten[g] > r_i) && (partner[punkt_ladung[i - 1]][g] < (pscanzahl + 1)))
                  {
                    run << "Chargetransport" << std::endl;
                    punkt_ladung[i] = partner[punkt_ladung[i - 1]][g];
                    punkt[i] = punkt[i - 1];

                    // Abbruchkriterium für Ladungstrennung
                    switch (direction)
                    {
                    case 'x':
                      if (((x[punkt_ladung[i]]) - x_mittel) > (0.75*(x[startpunkt[k]] - x_mittel)))
                      {
                        ch_diss[k]++;
                        zeit_ch[k][j] = zeit - zeit_ex[k][j];
                        vel_ch[k][j] = ((x[punkt_ladung[i]]) - x_mittel) / zeit_ch[k][j];

                        run << "Charges separated" << std::endl;
                        zustand[k][j] = 's';
                      }
                      break;
                    case 'y':
                      if (((y[punkt_ladung[i]]) - y_mittel) > (0.75*(y[startpunkt[k]] - y_mittel)))
                      {
                        ch_diss[k]++;
                        zeit_ch[k][j] = zeit - zeit_ex[k][j];
                        vel_ch[k][j] = ((y[punkt_ladung[i]]) - y_mittel) / zeit_ch[k][j];

                        run << "Charges separated" << std::endl;
                        zustand[k][j] = 's';
                      }
                      break;
                    case 'z':
                      if (((z[punkt_ladung[i]]) - z_mittel) > (0.6*(z[startpunkt[k]] - z_mittel)))
                      {
                        ch_diss[k]++;
                        zeit_ch[k][j] = zeit - zeit_ex[k][j];
                        vel_ch[k][j] = ((z[punkt_ladung[i]]) - z_mittel) / zeit_ch[k][j];

                        run << "Charges separated" << std::endl;
                        zustand[k][j] = 's';
                      }
                      break;
                    }

                    //#########################################################################################################################
                    run << "old Monomer " << std::setw(5) << punkt_ladung[i - 1] << std::endl;
                    run << "new Monomer " << std::setw(5) << punkt_ladung[i] << std::endl;
                    run << "Coupling " << std::setw(12) << std::setprecision(6) << std::fixed << coupling_ladung[punkt_ladung[i - 1]][punkt_ladung[i]] << std::endl;
                    run << "Fulleren " << std::setw(5) << punkt[i] << std::setw(5) << punkt[i - 1] << std::endl;
                    break;
                  }
                  else if ((raten[g] > r_i) && (partner[punkt_ladung[i - 1]][g] > (pscanzahl)))
                  {
                    run << "Rekombination" << std::endl;
                    zustand[k][j] = 't';
                    rek[k]++;
                    break;
                  }
                  else if (g == (partneranzahl[punkt_ladung[i - 1]]))
                  {
                    run << "WARING: ERROR during p-semiconductor chargetransfer." << std::endl;
                    throw std::runtime_error("WARING: ERROR during p-semiconductor chargetransfer. Aborting");
                  }
                }
              } // endes des hüpfenden p-halbleiters

              else if ((1 / r_summe - zeit_1) > (1 / r_summe_fulleren - zeit_2))
              {
                run << "Fulleren hopped first." << std::endl;
                if ((1 / r_summe_fulleren - zeit_2) > 0)
                {
                  zeit = zeit + (1 / r_summe_fulleren - zeit_2);
                  zeit_1 = zeit_1 + (1 / r_summe_fulleren - zeit_2);
                  zeit_2 = 0;
                }
                else if ((1 / r_summe_fulleren - zeit_2) < 0)
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
                const double r_i = zufall * r_summe_fulleren;

                for (std::size_t g = 1u; g < ((partneranzahl[punkt[i - 1]]) + 1); g++)
                {
                  if ((raten_fulleren[g] > r_i) && ((partner[punkt[i - 1]][g]) > pscanzahl))
                  {
                    run << "Chargetransfer in Fullerenephase" << std::endl;

                    punkt[i] = partner[punkt[i - 1]][g];
                    punkt_ladung[i] = punkt_ladung[i - 1];
                    run << "old Fulleren " << std::setw(5) << punkt[i - 1] << std::endl;
                    run << "new Fulleren " << std::setw(5) << punkt[i] << std::endl;
                    run << "Monomer " << std::setw(5) << punkt_ladung[i] << std::endl;
                    run << "Coupling " << std::setw(12) << std::setprecision(6) << coupling_fulleren[punkt[i]][punkt[i - 1]] << std::endl;
                    break;
                  }
                  else if ((raten_fulleren[g] > r_i) && ((partner[punkt[i - 1]][g]) < (pscanzahl + 1)))
                  {

                    run << "Rekombination." << std::endl;
                    zustand[k][j] = 't';
                    rek[k]++;
                    break;
                  }
                  else if (g == (partneranzahl[punkt[i - 1]]))
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
              std::random_device rd;
              std::default_random_engine engine(rd());
              std::normal_distribution<double> distribution0(0.0, 0.0338987); // hier neue standardabweichung eintragen
              const double zufall1 = distribution0(engine); //generating an normal-distributed random number

              std::vector <double> raten(partneranzahl[punkt[i - 1]] + 1);

              for (std::size_t h = 1; h < (partneranzahl[punkt[i - 1]] + 1); h++)
              {
                if (partner[punkt[i - 1]][h] < (pscanzahl + 1))
                {
                  const double zufall = distribution0(engine);// generatinjg a second normal distributed random number
                  const double testrate = rate(coupling_exciton[punkt[i - 1]][partner[punkt[i - 1]][h]], (zufall - zufall1), reorganisationsenergie_exciton);
                  r_summe += testrate;
                  raten[h] = r_summe;
                  run << "A: " << punkt[i - 1] << "   B: " << partner[punkt[i - 1]][h] << "rate   " << testrate << std::endl;
                }

                else if (partner[punkt[i - 1]][h] > (pscanzahl))
                {
                  const double zufall = distribution0(engine);
                  // coulomb energie berechnen
                  const double coulombenergy = coulomb(x, y, z, punkt[i - 1], partner[punkt[i - 1]][h], 1);
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

              for (std::size_t g = 1u; g < (partneranzahl[punkt[i - 1]] + 1); g++)
              {
                if (raten[g] > r_i)
                {
                  punkt[i] = partner[punkt[i - 1]][g];

                  if (punkt[i] < (pscanzahl + 1))
                  {
                    run << "hopped to " << punkt[i] << std::endl;
                  }
                  else if (punkt[i] > pscanzahl)
                  {
                    punkt[i] = partner[punkt[i - 1]][g];
                    punkt_ladung[i] = punkt[i - 1];
                    zustand[k][j] = 'c';
                    run << "Chargeseparation." << std::endl;

                    vel_ex[k][j] = length(x, y, z, punkt[0], punkt[i]) / zeit;
                    //run << "Exzitonspeed " << vel_ex[k][j] * 1e-9 << std::endl;
                    ex_diss[k]++;
                  }

                  break;
                }
                else if (raten[partneranzahl[punkt[i - 1]]] < r_i)
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
      std::vector <double> mittel_ex_vel(index + 1), mittel_ch_vel(index + 1), standard_ex(index + 1), standard_ch(index + 1);
      std::vector <std::size_t> const& ex_diss = m_results.ex_diss;
      std::vector <std::size_t> const& ch_diss = m_results.ch_diss;
      std::vector <std::size_t> const& rek = m_results.rek;
      std::vector <std::size_t> const& trapping = m_results.trapping;
      std::vector <std::size_t> const& radiativ = m_results.radiativ;

      std::vector <std::vector<double>> const & vel_ex = m_results.vel_ex;
      std::vector <std::vector<double>> const & vel_ch = m_results.vel_ch;
      std::vector <std::vector<double>> const & zeit_ex = m_results.zeit_ex;
      std::vector <std::vector<double>> const & zeit_ch = m_results.zeit_ch;

      for (std::size_t k = 1u; k < (index + 1); k++)
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

      auswertung << std::setw(9) << "Average " << std::setw(9) << std::setprecision(5) << std::fixed << mittel_ex / index << std::setw(9) << std::setprecision(5) << std::fixed << mittel_ch / index;
      auswertung << std::setw(9) << std::setprecision(5) << std::fixed << mittel_rek / index << std::setw(9) << std::setprecision(5) << std::fixed << mittel_trapp / index;
      auswertung << std::setw(9) << std::setprecision(5) << std::fixed << mittel_rad / index << '\n';
      auswertung << "Velocities" << '\n';
      auswertung << std::setw(4) << "k" << std::setw(5) << "IX" << std::setw(11) << "Ex_vel" << std::setw(11) << "Ex_s_dev" << std::setw(11) << "Ch_vel" << std::setw(11) << "Ch_s_dev" << '\n';

      // mittlere Geschwindigkeit
      for (std::size_t k = 1; k < (index + 1); k++)
      {
        mittel_ch_vel[k] = 0;
        mittel_ex_vel[k] = 0;

        // Mittelwert berechnen
        for (std::size_t j = 1; j < vectorlen; j++)
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
        for (std::size_t j = 1; j < vectorlen; j++)
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
      for (std::size_t k = 1; k < (index + 1); k++)
      {
        mittelwert_geschw_exciton = mittelwert_geschw_exciton + mittel_ex_vel[k];
        mittelwert_geschw_ladung = mittelwert_geschw_ladung + mittel_ch_vel[k];
      }
      auswertung << std::left << std::setw(7) << " Average    " << std::left << std::setw(22) << std::setprecision(5) << std::fixed << mittelwert_geschw_exciton / index * 1e-9;
      auswertung << std::left << std::setw(9) << std::setprecision(5) << std::fixed << mittelwert_geschw_ladung / index * 1e-9 << '\n';

      // Verteilung Ladungen und Exzitonengeschwindigkeiten
      std::ofstream exciton_verteilung;
      exciton_verteilung.open("exciton_distribution.txt");
      for (std::size_t i = 1; i < 21; i++) {
        std::size_t zahl = 0u;
        for (std::size_t k = 1; k < (index + 1); k++) {
          for (std::size_t j = 1; j < 101; j++) {
            if ((vel_ex[k][j] > (i * 50 * 1e9))) {
              zahl++;
            }
          }
        }
        exciton_verteilung << std::setw(9) << std::setprecision(5) << i * 50 << std::setw(9) << zahl / index << '\n';
      }

      exciton_verteilung.close();
      exciton_verteilung.open("charge_distribution.txt");
      for (std::size_t i = 1; i < 21; i++) {
        std::size_t zahl = 0u;
        for (std::size_t k = 1; k < (index + 1); k++) {
          for (std::size_t j = 1; j < 101; j++) {
            if ((vel_ch[k][j] > (i * 50 * 1e9))) {
              zahl++;
            }
          }
        }
        exciton_verteilung << std::setw(9) << std::setprecision(5) << i * 50 << std::setw(9) << zahl / index << '\n';
      }
      exciton_verteilung.close();
    }



    std::size_t gesamtanzahl;
    std::size_t pscanzahl, nscanzahl;
    std::size_t numberOfExcitonPairs, numberOfNSemiconductorHomopairs, numberOfHeteroDimers;
    std::vector <std::size_t> partneranzahl;
    std::vector<std::vector<std::size_t>> partner;
    std::vector <std::vector<double>> coupling_exciton;
    std::vector <std::vector<double>> coupling_ladung;
    std::vector <std::vector<double>> coupling_ct;
    std::vector <std::vector<double>> coupling_rek;
    std::vector <std::vector<double>> coupling_fulleren;
    std::vector <double> x, y, z;
    std::vector <std::size_t> startpunkt;

    double x_mittel, y_mittel, z_mittel;
    double index;

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
    const double procentualDist2Interf; //0.85
    const std::size_t vectorlen;

    struct results {
      std::vector <std::size_t> ex_diss, ch_diss, rek, trapping, radiativ;
      std::vector <std::vector<double>> vel_ex, vel_ch, zeit_ex, zeit_ch;
      results(std::size_t index, std::size_t vectorlen)
        : ex_diss(index + 1), ch_diss(index + 1), rek(index + 1), trapping(index + 1), radiativ(index + 1),
          vel_ex(index + 1, std::vector <double>(vectorlen)),
          vel_ch(index + 1, std::vector <double>(vectorlen)),
          zeit_ex(index + 1, std::vector <double>(vectorlen)),
          zeit_ch(index + 1, std::vector <double>(vectorlen))
      {
        for (std::size_t i = 1; i < (index + 1); i++) //initializing the vectors with 0
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



  int exciton_breakup(int pscanzahl, int nscanzahl, char ebene, std::string masscenters, std::string nscpairrates,
    std::string pscpairexrates, std::string pscpairchrates, std::string pnscpairrates)
  {
    ////////////////////////////////////// hier Parameter angeben ///////////////////////////////////////////////
    const double reorganisationsenergie_exciton = Config::get().exbreak.ReorgE_exc; //0.561; // SCS-CC2 wert: vorher 0.561 // vermutlich in eV
    const double reorganisationsenergie_ladung = Config::get().exbreak.ReorgE_ch;//0.194;
    const double fullerenreorganisationsenergie = Config::get().exbreak.ReorgE_nSC;//0.178;
    const double ct_reorganisation = Config::get().exbreak.ReorgE_ct;//0.156;
    const double chargetransfertriebkraft = Config::get().exbreak.ct_triebkraft;//1.550;
    const double rekombinationstriebkraft = Config::get().exbreak.rek_triebkraft;//-4.913;
    const double rek_reorganisation = Config::get().exbreak.ReorgE_rek;//0.184;
    const double oszillatorstrength = Config::get().exbreak.oscillatorstrength;//0.0852;
    const double wellenzahl = Config::get().exbreak.wellenzahl;//28514.91;
    const double k_rad = wellenzahl * wellenzahl*oszillatorstrength; // fluoreszenz
    const double procentualDist2Interf = 0.5; //0.85

    /////////////////////////////////// INPUT-READING
    std::ifstream schwerpunkt;
    schwerpunkt.open(masscenters);
    std::size_t gesamtanzahl = 0u;
    schwerpunkt >> gesamtanzahl;
    std::cout << "Read number of Monomers: " << gesamtanzahl << std::endl;
    if (gesamtanzahl == pscanzahl + nscanzahl) //test if correct number of molecules was given
    {
      std::cout << "Number of monomers is correct, proceeding." << std::endl;
    }
    else //gesamtanzahl != pscanzahl + nscanzahl
    {
      std::cout << "Wrong number of monomers detected!" << std::endl;
      throw std::logic_error("Wrong numbers of p- and n-type molecules in inputfile. Expected: " + std::to_string(gesamtanzahl) + " | Is: " + std::to_string(pscanzahl + nscanzahl));
      return -1;
    }
    schwerpunkt >> skipline;
    std::vector <double> x(gesamtanzahl + 1), y(gesamtanzahl + 1), z(gesamtanzahl + 1);

    for (std::size_t i = 1u; i < (gesamtanzahl + 1u); i++)
    {
      std::string zeile;
      schwerpunkt >> zeile >> x[i] >> y[i] >> z[i]; //reading and saving balance points for molecules
    }
    ///////////////////////////////////
    std::ifstream exciton;
    exciton.open(pscpairexrates);
    std::size_t numberOfExcitonPairs = 0u;
    std::string zeile;
    while (getline(exciton, zeile))
    { //counting of excitonpairs in homodimers
      numberOfExcitonPairs++;
    }
    exciton.close();
    ///////////////////////////////////

    std::vector <std::vector<double>> coupling_exciton(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1)); //in original code 2d-arrays were used
    std::vector <std::vector<double>> coupling_ladung(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1));
    std::vector <std::vector<double>> coupling_ct(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1));
    std::vector <std::vector<double>> coupling_rek(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1));
    std::vector <std::vector<double>> coupling_fulleren(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1));

    for (std::size_t i = 1u; i < (gesamtanzahl + 1u); i++) //initializaton of the 2d vectors with 0 in all places
    {
      for (std::size_t j = 1u; j < (gesamtanzahl + 1u); j++) {
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
    std::vector <int> partneranzahl(pscanzahl + nscanzahl + 1);
    for (std::size_t i = 1u; i < (pscanzahl + nscanzahl + 1); i++)
    {
      partneranzahl[i] = 0;
    } //inizialisation of all elements with 0

    for (std::size_t i = 1u; i < (numberOfExcitonPairs + 1); i++) // counting homodimer partners for j
    {
      for (std::size_t j = 1u; j < (pscanzahl + 1); j++)
      {
        if ((exciton_1[i] == j) || (exciton_2[i] == j))
        {
          partneranzahl[j]++;
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

    size_t numberOfHeteroDimers = 0u;
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
      for (std::size_t j = 1; j < (gesamtanzahl + 1); j++)
      {
        if ((hetero_1[i] == j) || (hetero_2[i] == j))
        {
          partneranzahl[j]++;
        }
      }
    }
    exciton.close();

    ////////////////////////////////////////
    exciton.open(nscpairrates);
    std::size_t numberOfNSemiconductorHomopairs = 0u;
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
      for (std::size_t j = 1; j < (gesamtanzahl + 1); j++)
      {
        if ((fulleren_1[i] == j) || (fulleren_2[i] == j))
        {
          partneranzahl[j]++;
        }
      }
    }

    std::vector<std::vector<int>> partner(gesamtanzahl + 1, std::vector<int>());//2D-vector with variing length for second vector

    for (std::size_t i = 1; i < (gesamtanzahl + 1); i++) //dynamic allocation for length of 2nd vector
    {
      partner[i].resize(partneranzahl[i] + 1);
    }

    for (std::size_t i = 1; i < (gesamtanzahl + 1); i++) //initializing all elements of vector partner with 0
    {
      for (std::size_t j = 1; j < (partneranzahl[i] + 1); j++) {
        partner[i][j] = 0;
      }
    }

    for (std::size_t i = 1; i < (gesamtanzahl + 1); i++)
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
          j++; // since the 2nd dimension length of vector partner was set to partneranzahl[i] in thes logic construction j must always end up to be equal to partneranzahl[i]
        }
      }
      if (partneranzahl[i] != j - 1) //Sanity Check
      {
        std::cout << "Error with number of partners for monomer " << i << std::endl;
        throw std::runtime_error("Error with number of partners for monomer " + std::to_string(i) + ". Aborting.");
      }
    }

    // INPUT-END

    std::ofstream kopplung;
    kopplung.open("partner.txt"); // Writing file partner.txt
    for (std::size_t i = 1u; i < (gesamtanzahl + 1); i++)
    {
      kopplung << std::setw(6) << i << std::setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess
      for (std::size_t j = 1; j < (partneranzahl[i] + 1); j++)
      {
        kopplung << std::setw(6) << partner[i][j]; //writes the indices of the partners j
      }
      kopplung << '\n';
    }
    kopplung.close();

    kopplung.open("couplings.txt"); // Writing file couplings.txt
    for (std::size_t i = 1u; i < (pscanzahl + 1); i++)
    {
      kopplung << std::setw(6) << i << std::setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess

      for (std::size_t j = 1; j < (partneranzahl[i] + 1); j++)
      {
        kopplung << std::setw(6) << partner[i][j]; //writes the indices of the partners j
        if (partner[i][j] < (pscanzahl + 1))
        {
          kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_exciton[i][partner[i][j]]; //writes the exciton-coupling between i and j
        }
        else if (partner[i][j] > pscanzahl)
        {
          kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_ct[i][partner[i][j]];     //writes charge-transfer-coupling between i and j 
        }
      }
      kopplung << '\n';
    }

    for (std::size_t i = (pscanzahl + 1); i < (gesamtanzahl + 1); i++)
    {
      kopplung << std::setw(6) << i << std::setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess
      for (std::size_t j = 1; j < (partneranzahl[i] + 1); j++)
      {
        kopplung << std::setw(6) << partner[i][j]; //writes the indices of the partners j
        if (partner[i][j] < (pscanzahl))
        {
          kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_rek[i][partner[i][j]]; //writes the recombination coupling between i and j
        }
        else if (partner[i][j] > pscanzahl)
        {
          kopplung << std::setw(12) << std::setprecision(6) << std::fixed << coupling_fulleren[i][partner[i][j]]; // writes some coupling regarding fullerens?
        }
      }
      kopplung << '\n';
    }
    kopplung.close();

    // Startpunkte bestimmen ##################################################################################################################
    double x_monomer(0.), y_monomer(0.), z_monomer(0.), x_fulleren(0.), y_fulleren(0.), z_fulleren(0.),
      x_gesamt(0.), y_gesamt(0.), z_gesamt(0.), x_mittel(0.), y_mittel(0.), z_mittel(0.);

    for (std::size_t i = 1u; i < (pscanzahl + 1); i++)
    {
      x_monomer += (x[i] / pscanzahl);
      y_monomer += (y[i] / pscanzahl);
      z_monomer += (z[i] / pscanzahl);
    }

    for (std::size_t i = (pscanzahl + 1); i < (gesamtanzahl + 1); i++)     //using fact, that fullerens always have larger indices than other monomers
    {
      x_fulleren += (x[i] / nscanzahl);
      y_fulleren += (y[i] / nscanzahl);
      z_fulleren += (z[i] / nscanzahl);
    }

    for (std::size_t i = 1u; i < (gesamtanzahl + 1); i++)
    {
      x_gesamt += (x[i] / gesamtanzahl);
      y_gesamt += (y[i] / gesamtanzahl);
      z_gesamt += (z[i] / gesamtanzahl);
    }

    x_mittel = (x_monomer + x_fulleren) / 2;
    y_mittel = (y_monomer + y_fulleren) / 2;
    z_mittel = (z_monomer + z_fulleren) / 2;

    std::ofstream interface;
    interface.open("massponts_general.xyz"); //writing out average balance points for all groupings of monomers
    interface << "4" << '\n' << '\n';
    interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x_monomer << std::setw(12) << std::setprecision(6) << std::fixed << y_monomer << std::setw(12) << std::setprecision(6) << std::fixed << z_monomer << '\n';
    interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x_fulleren << std::setw(12) << std::setprecision(6) << std::fixed << y_fulleren << std::setw(12) << std::setprecision(6) << std::fixed << z_fulleren << '\n';
    interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x_gesamt << std::setw(12) << std::setprecision(6) << std::fixed << y_gesamt << std::setw(12) << std::setprecision(6) << std::fixed << z_gesamt << '\n';
    interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x_mittel << std::setw(12) << std::setprecision(6) << std::fixed << y_mittel << std::setw(12) << std::setprecision(6) << std::fixed << z_mittel << '\n';
    interface.close();

    double index = 0.;
    double max = 0.;
    std::vector <int> startpunkt(pscanzahl + 1);

    switch (ebene) 
    { //different cases for the possible planes of the interface
    case 'x':
      for (std::size_t i = 1u; i < (pscanzahl + 1); i++) //determining the maximal distance to interface
      {
        if (x[i] > max) {
          max = x[i];
        }
      }
      std::cout << "Maxdistance is " << std::setw(12) << std::setprecision(6) << std::fixed << max << std::endl;

      for (std::size_t i = 1; i < (pscanzahl + 1); i++)  //determining the necessary number of starting points? 
      {
        if ((x[i] - x_mittel) > (procentualDist2Interf*(max - x_mittel))) {
          index++;
          startpunkt[index] = i;
        }
      }
      break;

    case 'y':
      for (std::size_t i = 1u; i < (pscanzahl + 1); i++)  //determining the maximal distance to interface
      {
        if (y[i] > max) {
          max = y[i];
        }
      }
      std::cout << "Maxdistance is " << std::setw(12) << std::setprecision(6) << std::fixed << max << std::endl;

      for (std::size_t i = 1; i < (pscanzahl + 1); i++) //determining the necessary number of starting points? 
      {
        if ((y[i] - y_mittel) > (procentualDist2Interf*(max - y_mittel))) {
          index++;
          startpunkt[index] = i;

        }
      }
      break;

    case 'z':
      for (std::size_t i = 1u; i < (pscanzahl + 1); i++) //determining the maximal distance to interace
      {
        if (z[i] > max) {
          max = z[i];
        }
      }
      std::cout << "Maxdistance is " << std::setw(12) << std::setprecision(6) << std::fixed << max << std::endl;

      for (std::size_t i = 1; i < (pscanzahl + 1); i++) //determining the necessary number of starting points? 
      {
        if ((z[i] - z_mittel) > (procentualDist2Interf*(max - z_mittel)))
        {
          index++;
          startpunkt[index] = i;

        }
      }
      break;
    }


    interface.open("startingpoint.xyz");
    interface << index << '\n' << '\n'; //writes the number of startingponts
    for (std::size_t i = 1; i < (index + 1); i++) 	//writes coordinates of startingpoints
    {
      interface << std::setw(5) << "X" << std::setw(12) << std::setprecision(6) << std::fixed << x[startpunkt[i]];
      interface << std::setw(12) << std::setprecision(6) << std::fixed << y[startpunkt[i]];
      interface << std::setw(12) << std::setprecision(6) << std::fixed << z[startpunkt[i]] << '\n';
    }
    interface.close();

    // ################################################################################## Beginn der Simulation ##############################################################################
    // Variablen


    // Schrittanzahl pro MC-Simulation
    std::size_t schritt = 2 * (numberOfExcitonPairs + numberOfNSemiconductorHomopairs) + 400u; // MAGIC NUMBER?
    std::cout << "Number of stepps for MC-Simulation: " << schritt << "." << std::endl;

    // ###################################################################################
    constexpr std::size_t vectorlen = 101u;

    std::vector <std::vector<double>> vel_ex(index + 1, std::vector <double>(vectorlen)),
      vel_ch(index + 1, std::vector <double>(vectorlen)),
      zeit_ex(index + 1, std::vector <double>(vectorlen)),
      zeit_ch(index + 1, std::vector <double>(vectorlen));
    std::vector <int> ex_diss(index + 1), ch_diss(index + 1), rek(index + 1), trapping(index + 1), radiativ(index + 1);
    std::vector <std::vector<char>> zustand(index + 1, std::vector <char>(vectorlen));
    std::vector <int> punkt(schritt + 1), punkt_ladung(schritt + 1);

    for (std::size_t i = 1; i < (index + 1); i++) //initializing the vectors with 0
    {
      ex_diss[i] = 0;
      ch_diss[i] = 0;
      rek[i] = 0;
      radiativ[i] = 0;
      trapping[i] = 0;
      for (std::size_t j = 1; j < vectorlen; j++)
      {
        vel_ex[i][j] = 0;
        vel_ch[i][j] = 0;
        zeit_ex[i][j] = 0;
        zeit_ch[i][j] = 0;
        zustand[i][j] = 'e';
      }
    }

    // k: index für startpunkte
    // j: index für durchläufe
    // i: index für schritt
    for (std::size_t i = 0u; i < (schritt + 1); i++)
    {
      punkt[i] = 0;
    }

    std::ofstream run;
    if (Config::get().general.verbosity >= 4u)
      run.open("run.txt");

    for (std::size_t k = 1; k < (index + 1); k++) // schleife über startpunkte "index durch 1 vertauscht"
    {
      run << "k ist " << k << std::endl;

      for (std::size_t j = 1; j < vectorlen; j++)   // schleife über durchläufe für den gleichen startpunkt " 101 durch 11 vertauscht"
      {
        double zeit(0.), zeit_1(0.), zeit_2(0.);

        punkt[0] = startpunkt[k];

        for (std::size_t i = 1; i < (schritt + 1); i++)
        {

          if (zustand[k][j] == 'c')
          {
            double r_summe = 0.;
            // site energies berechnen
            //########## raten addieren für monomere ##########
            
            std::random_device rd;
            std::default_random_engine engine(rd());
            std::normal_distribution<double> distribution0(0.0, 0.068584577); // TODO: hier neue standardabweichung eintragen
            double zufall1 = distribution0(engine); //generating normal-distributed random number

            std::vector<double> raten(partneranzahl[punkt_ladung[i - 1]] + 1);
            for (std::size_t h = 0; h < (partneranzahl[punkt_ladung[i - 1]] + 1); h++)
            {
              const double zufall = distribution0(engine);
              if (partner[punkt_ladung[i - 1]][h] < (pscanzahl + 1))
              {
                const double coulombenergy = coulomb(x, y, z, punkt[i - 1], partner[punkt_ladung[i - 1]][h], 3.4088) - coulomb(x, y, z, punkt[i - 1], punkt_ladung[i - 1], 3.4088);
                r_summe = r_summe + rate(coupling_ladung[punkt_ladung[i - 1]][partner[punkt_ladung[i - 1]][h]], ((zufall - zufall1) + coulombenergy), reorganisationsenergie_ladung);

                raten[h] = r_summe;
              }
              if ((partner[punkt_ladung[i - 1]][h] > (pscanzahl)) && (partner[punkt_ladung[i - 1]][h] == punkt[i - 1]))
              {
                // coulomb energie berechnen	   
                const double coulombenergy = coulomb(x, y, z, punkt_ladung[i - 1], partner[punkt_ladung[i - 1]][h], 1);
                r_summe = r_summe + rate(coupling_rek[punkt_ladung[i - 1]][partner[punkt_ladung[i - 1]][h]], (zufall - zufall1) + rekombinationstriebkraft - coulombenergy, rek_reorganisation);
                raten[h] = r_summe;
              }
              //ACHTUNG: hier Korrektur ///////////////////////////////////////////////////
              else if ((partner[punkt_ladung[i - 1]][h] > (pscanzahl)) && (partner[punkt_ladung[i - 1]][h] != punkt[i - 1]))
              {
                r_summe = r_summe;
                raten[h] = 0;
              }
            }

            // hier raten für fullerene addieren 
            double r_summe_fulleren = 0.;
            zufall1 = distribution0(engine);
            std::vector <double> raten_fulleren(partneranzahl[punkt[i - 1]] + 1);
            for (std::size_t h = 1; h < (partneranzahl[punkt[i - 1]] + 1); h++)
            {
              if (partner[punkt[i - 1]][h] > (pscanzahl))
              {
                const double zufall = distribution0(engine);
                const double coulombenergy = coulomb(x, y, z, punkt_ladung[i - 1], partner[punkt[i - 1]][h], 3.4088) - coulomb(x, y, z, punkt[i - 1], punkt_ladung[i - 1], 3.4088);
                r_summe_fulleren = r_summe_fulleren + rate(coupling_fulleren[punkt[i - 1]][partner[punkt[i - 1]][h]], ((zufall - zufall1) + coulombenergy), fullerenreorganisationsenergie);

                raten_fulleren[h] = r_summe_fulleren;
              }
              if ((partner[punkt[i - 1]][h] < (pscanzahl + 1)) && (partner[punkt[i - 1]][h] == punkt_ladung[i - 1]))
              {
                const double zufall = distribution0(engine);

                // coulomb energie berechnen
                const double coulombenergy = coulomb(x, y, z, punkt[i - 1], partner[punkt[i - 1]][h], 1);

                r_summe_fulleren = r_summe_fulleren + rate(coupling_rek[punkt[i - 1]][partner[punkt[i - 1]][h]], (zufall - zufall1) + rekombinationstriebkraft - coulombenergy, rek_reorganisation);

                raten_fulleren[h] = r_summe_fulleren;
              }
              //ACHTUNG: hier Korrektur ////////////////////////////////////////////////////////
              else if ((partner[punkt[i - 1]][h] < (pscanzahl + 1)) && (partner[punkt[i - 1]][h] != punkt_ladung[i - 1]))
              {
                r_summe_fulleren = r_summe_fulleren;
                raten_fulleren[h] = 0;
              }
              //////////////////////////////////////////////////////////////////////////////////
            }

            // hüpfendes teilchen bestimmen
            if ((1 / r_summe - zeit_1) < (1 / r_summe_fulleren - zeit_2))
            {
              run << "Monomer hopps first." << std::endl;

              //Update der Zeiten
              if ((1. / r_summe - zeit_1) > 0.)
              {
                zeit = zeit + (1. / r_summe - zeit_1);
                zeit_2 = zeit_2 + (1 / r_summe - zeit_1);
                zeit_1 = 0.;
              }
              else if ((1. / r_summe - zeit_1) < 0.)
              {
                zeit = zeit;
                zeit_2 = zeit_2;
                zeit_1 = 0.;
              }
              else 
              {
                run << "ERROR!" << std::endl;
                throw std::runtime_error("Critical Error in Exciton Breakup Task, Aborting!");
                return 0;
              }
              // monomerhüpfen ausführen
              std::uniform_real_distribution<double> distribution1(0, 1);
              const double zufall = distribution1(engine);
              const double r_i = zufall * r_summe;
              for (std::size_t g = 1; g < (partneranzahl[punkt_ladung[i - 1]] + 1); g++)
              {
                if ((raten[g] > r_i) && (partner[punkt_ladung[i - 1]][g] < (pscanzahl + 1)))
                {
                  run << "Chargetransport" << std::endl;
                  punkt_ladung[i] = partner[punkt_ladung[i - 1]][g];
                  punkt[i] = punkt[i - 1];

                  // Abbruchkriterium für Ladungstrennung
                  switch (ebene)
                  {
                  case 'x':
                    if (((x[punkt_ladung[i]]) - x_mittel) > (0.75*(x[startpunkt[k]] - x_mittel)))
                    {
                      ch_diss[k]++;
                      zeit_ch[k][j] = zeit - zeit_ex[k][j];
                      vel_ch[k][j] = ((x[punkt_ladung[i]]) - x_mittel) / zeit_ch[k][j];

                      run << "Charges separated" << std::endl;
                      zustand[k][j] = 's';
                    }
                    break;
                  case 'y':
                    if (((y[punkt_ladung[i]]) - y_mittel) > (0.75*(y[startpunkt[k]] - y_mittel)))
                    {
                      ch_diss[k]++;
                      zeit_ch[k][j] = zeit - zeit_ex[k][j];
                      vel_ch[k][j] = ((y[punkt_ladung[i]]) - y_mittel) / zeit_ch[k][j];

                      run << "Charges separated" << std::endl;
                      zustand[k][j] = 's';
                    }
                    break;
                  case 'z':
                    if (((z[punkt_ladung[i]]) - z_mittel) > (0.6*(z[startpunkt[k]] - z_mittel)))
                    {
                      ch_diss[k]++;
                      zeit_ch[k][j] = zeit - zeit_ex[k][j];
                      vel_ch[k][j] = ((z[punkt_ladung[i]]) - z_mittel) / zeit_ch[k][j];

                      run << "Charges separated" << std::endl;
                      zustand[k][j] = 's';
                    }
                    break;
                  }

                  //#########################################################################################################################
                  run << "old Monomer " << std::setw(5) << punkt_ladung[i - 1] << std::endl;
                  run << "new Monomer " << std::setw(5) << punkt_ladung[i] << std::endl;
                  run << "Coupling " << std::setw(12) << std::setprecision(6) << std::fixed << coupling_ladung[punkt_ladung[i - 1]][punkt_ladung[i]] << std::endl;
                  run << "Fulleren " << std::setw(5) << punkt[i] << std::setw(5) << punkt[i - 1] << std::endl;
                  break;
                }
                else if ((raten[g] > r_i) && (partner[punkt_ladung[i - 1]][g] > (pscanzahl)))
                {
                  run << "Rekombination" << std::endl;
                  zustand[k][j] = 't';
                  rek[k]++;
                  break;
                }
                else if (g == (partneranzahl[punkt_ladung[i - 1]]))
                {
                  run << "WARING: ERROR during p-semiconductor chargetransfer." << std::endl;
                  throw std::runtime_error("WARING: ERROR during p-semiconductor chargetransfer. Aborting");
                  return 0;
                }
              }
            } // endes des hüpfenden p-halbleiters

            else if ((1 / r_summe - zeit_1) > (1 / r_summe_fulleren - zeit_2))
            {
              run << "Fulleren hopped first." << std::endl;
              if ((1 / r_summe_fulleren - zeit_2) > 0)
              {
                zeit = zeit + (1 / r_summe_fulleren - zeit_2);
                zeit_1 = zeit_1 + (1 / r_summe_fulleren - zeit_2);
                zeit_2 = 0;
              }
              else if ((1 / r_summe_fulleren - zeit_2) < 0)
              {
                zeit = zeit;
                zeit_1 = zeit_1;
                zeit_2 = 0;
              }
              else {
                run << "ERROR!" << std::endl;
                throw std::runtime_error("Critical Error in Exciton Breakup Task, Aborting!");
                return 0;
              }

              // fullerenhüpfen ausführen
              std::uniform_real_distribution<double> distribution1(0, 1);
              const double zufall = distribution1(engine);
              const double r_i = zufall * r_summe_fulleren;

              for (std::size_t g = 1u; g < ((partneranzahl[punkt[i - 1]]) + 1); g++)
              {
                if ((raten_fulleren[g] > r_i) && ((partner[punkt[i - 1]][g]) > pscanzahl))
                {
                  run << "Chargetransfer in Fullerenephase" << std::endl;

                  punkt[i] = partner[punkt[i - 1]][g];
                  punkt_ladung[i] = punkt_ladung[i - 1];
                  run << "old Fulleren " << std::setw(5) << punkt[i - 1] << std::endl;
                  run << "new Fulleren " << std::setw(5) << punkt[i] << std::endl;
                  run << "Monomer " << std::setw(5) << punkt_ladung[i] << std::endl;
                  run << "Coupling " << std::setw(12) << std::setprecision(6) << coupling_fulleren[punkt[i]][punkt[i - 1]] << std::endl;
                  break;
                }
                else if ((raten_fulleren[g] > r_i) && ((partner[punkt[i - 1]][g]) < (pscanzahl + 1)))
                {

                  run << "Rekombination." << std::endl;
                  zustand[k][j] = 't';
                  rek[k]++;
                  break;
                }
                else if (g == (partneranzahl[punkt[i - 1]]))
                {
                  run << "WARNING: ERROR during fullerene chargetransport." << std::endl;
                  throw std::runtime_error("Critical Error in Exciton Breakup Task, Aborting!");
                  return 0;
                }
              }
            } // endes des hüpfenden fullerens
          } // ende des 'c'-zustands
        //__________________________________________________________________________________________________________

          else if (zustand[k][j] == 'e')
          {
            // site energies berechnen
            double r_summe = 0.;
            std::random_device rd;
            std::default_random_engine engine(rd());
            std::normal_distribution<double> distribution0(0.0, 0.0338987); // hier neue standardabweichung eintragen
            const double zufall1 = distribution0(engine); //generating an normal-distributed random number

            std::vector <double> raten(partneranzahl[punkt[i - 1]] + 1);

            for (std::size_t h = 1; h < (partneranzahl[punkt[i - 1]] + 1); h++)
            {
              if (partner[punkt[i - 1]][h] < (pscanzahl + 1))
              {
                const double zufall = distribution0(engine);// generatinjg a second normal distributed random number
                const double testrate = rate(coupling_exciton[punkt[i - 1]][partner[punkt[i - 1]][h]], (zufall - zufall1), reorganisationsenergie_exciton);
                r_summe += testrate;
                raten[h] = r_summe;
                run << "A: " << punkt[i - 1] << "   B: " << partner[punkt[i - 1]][h] << "rate   " << testrate << std::endl;
              }

              else if (partner[punkt[i - 1]][h] > (pscanzahl))
              {
                const double zufall = distribution0(engine);
                // coulomb energie berechnen
                const double coulombenergy = coulomb(x, y, z, punkt[i - 1], partner[punkt[i - 1]][h], 1);
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

            for (std::size_t g = 1u; g < (partneranzahl[punkt[i - 1]] + 1); g++)
            {
              if (raten[g] > r_i)
              {
                punkt[i] = partner[punkt[i - 1]][g];

                if (punkt[i] < (pscanzahl + 1))
                {
                  run << "hopped to " << punkt[i] << std::endl;
                }
                else if (punkt[i] > pscanzahl)
                {
                  punkt[i] = partner[punkt[i - 1]][g];
                  punkt_ladung[i] = punkt[i - 1];
                  zustand[k][j] = 'c';
                  run << "Chargeseparation." << std::endl;

                  vel_ex[k][j] = length(x, y, z, punkt[0], punkt[i]) / zeit;
                  //run << "Exzitonspeed " << vel_ex[k][j] * 1e-9 << std::endl;
                  ex_diss[k]++;
                }

                break;
              }
              else if (raten[partneranzahl[punkt[i - 1]]] < r_i)
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
            return 0;
          } //end of undefined zustand
        //___________________________________________________________________________________________________

        } // ende über schleife i
      } // Ende über Schleife über durchläufe für den gleichen startpunkt j
    } // Ende über Schleife der Startpunkte k

    run.close();

    // Auswertung: Prozentsätze
    std::ofstream auswertung;
    auswertung.open("evaluation.txt");
    auswertung << std::setw(4) << "k" << std::setw(4) << "IX" << std::setw(9) << "Ex_Diss." << std::setw(9) << "Ch_Diss" << std::setw(9) << "Rek." << std::setw(9) << "Trapp." << std::setw(9) << "Fluor." << '\n';

    //double ex_diss_efficiency, ch_diss_efficiency, rek_efficiency, trapp_efficiency, rad_efficiency;
    double mittel_ex = 0;
    double mittel_ch = 0;
    double mittel_rek = 0;
    double mittel_trapp = 0;
    double mittel_rad = 0;
    std::vector <double> mittel_ex_vel(index + 1), mittel_ch_vel(index + 1), standard_ex(index + 1), standard_ch(index + 1);
    //int zahl;

    for (std::size_t k = 1u; k < (index + 1); k++) 
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

    auswertung << std::setw(9) << "Average " << std::setw(9) << std::setprecision(5) << std::fixed << mittel_ex / index << std::setw(9) << std::setprecision(5) << std::fixed << mittel_ch / index;
    auswertung << std::setw(9) << std::setprecision(5) << std::fixed << mittel_rek / index << std::setw(9) << std::setprecision(5) << std::fixed << mittel_trapp / index;
    auswertung << std::setw(9) << std::setprecision(5) << std::fixed << mittel_rad / index << '\n';
    auswertung << "Velocities" << '\n';
    auswertung << std::setw(4) << "k" << std::setw(5) << "IX" << std::setw(11) << "Ex_vel" << std::setw(11) << "Ex_s_dev" << std::setw(11) << "Ch_vel" << std::setw(11) << "Ch_s_dev" << '\n';

    // mittlere Geschwindigkeit
    for (std::size_t k = 1; k < (index + 1); k++) 
    {
      mittel_ch_vel[k] = 0;
      mittel_ex_vel[k] = 0;

      // Mittelwert berechnen
      for (std::size_t j = 1; j < vectorlen; j++) 
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
      for (std::size_t j = 1; j < vectorlen; j++) 
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
    for (std::size_t k = 1; k < (index + 1); k++) 
    {
      mittelwert_geschw_exciton = mittelwert_geschw_exciton + mittel_ex_vel[k];
      mittelwert_geschw_ladung = mittelwert_geschw_ladung + mittel_ch_vel[k];
    }
    auswertung << std::left << std::setw(7) << " Average    " << std::left << std::setw(22) << std::setprecision(5) << std::fixed << mittelwert_geschw_exciton / index * 1e-9;
    auswertung << std::left << std::setw(9) << std::setprecision(5) << std::fixed << mittelwert_geschw_ladung / index * 1e-9 << '\n';

    // Verteilung Ladungen und Exzitonengeschwindigkeiten
    std::ofstream exciton_verteilung;
    exciton_verteilung.open("exciton_distribution.txt");
    for (std::size_t i = 1; i < 21; i++) {
      std::size_t zahl = 0u;
      for (std::size_t k = 1; k < (index + 1); k++) {
        for (std::size_t j = 1; j < 101; j++) {
          if ((vel_ex[k][j] > (i * 50 * 1e9))) {
            zahl++;
          }
        }
      }
      exciton_verteilung << std::setw(9) << std::setprecision(5) << i * 50 << std::setw(9) << zahl / index << '\n';
    }

    exciton_verteilung.close();
    exciton_verteilung.open("charge_distribution.txt");
    for (std::size_t i = 1; i < 21; i++) {
      std::size_t zahl = 0u;
      for (std::size_t k = 1; k < (index + 1); k++) {
        for (std::size_t j = 1; j < 101; j++) {
          if ((vel_ch[k][j] > (i * 50 * 1e9))) {
            zahl++;
          }
        }
      }
      exciton_verteilung << std::setw(9) << std::setprecision(5) << i * 50 << std::setw(9) << zahl / index << '\n';
    }
    exciton_verteilung.close();

    return 0;
  }

}