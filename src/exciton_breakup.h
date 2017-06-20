#pragma once
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<limits>
#include<string>
#include<random>
#include<vector>

//Original code by Charlotte. German comments are from original code, english comments were inserted during implementation into CAST.
//Arrays were replaced by vectors.



istream& skipline(istream& in)
{
  return in.ignore(numeric_limits < streamsize >::max(), '\n');
}

// Definition der length-berechnung
double length(std::vector <double> arr1, std::vector <double> arr2, std::vector <double> arr3, int p, int q)
{
  double l = 1;
  l = sqrt((arr1[p] - arr1[q])*(arr1[p] - arr1[q]) + (arr2[p] - arr2[q])*(arr2[p] - arr2[q]) + (arr3[p] - arr3[q])*(arr3[p] - arr3[q]));
  return l;
}

double rate(double coupling, double energy, double reorganisation)
{
  double l = 1;
  double pi = 3.141592654;
  double h_quer = 1 / (2 * pi)*4.135667662e-15;
  double boltzmann_konstante = 8.6173303e-5; //  in gauß einheiten
  l = (coupling*coupling) / h_quer*sqrt(pi / (reorganisation*boltzmann_konstante * 298))*exp(-(reorganisation + energy)*(reorganisation + energy) / (4 * boltzmann_konstante * 298 * reorganisation));
  return l;
}

double coulomb(std::vector<double> arr1, std::vector<double> arr2, std::vector<double> arr3, int p, int q, double e_relative) {
  double l = 1;
  l = sqrt((arr1[p] - arr1[q])*(arr1[p] - arr1[q]) + (arr2[p] - arr2[q])*(arr2[p] - arr2[q]) + (arr3[p] - arr3[q])*(arr3[p] - arr3[q]));
  double pi = 3.141592654;
  double e_0 = 8.854187e-12;
  double elementar = 1.60217662e-19;
  double c = 1;
  c = -elementar / (4 * pi*e_0*e_relative*l*1e-10);
  return c;
}

int exciton_breakup(int pscanzahl, int nscanzahl, char ebene, std::string masscenters, std::string nscpairrates,
  std::string pscpairexrates, std::string pscpairchrates, std::string pnscpairrates)
{
  std::string zeile;
  int i, j, h, index, k, g;
  int gesamtanzahl;
  int e = 0;
  int het = 0;
  int full = 0;
  double max, zufall;



  ////////////////////////////////////// hier Parameter angeben ///////////////////////////////////////////////
  double reorganisationsenergie_exciton = Config::get().exbreak.ReorgE_exc; //0.561; // SCS-CC2 wert: vorher 0.561
  double reorganisationsenergie_ladung = Config::get().exbreak.ReorgE_ch;//0.194;
  double fullerenreorganisationsenergie = Config::get().exbreak.ReorgE_nSC;//0.178;
  double ct_reorganisation = Config::get().exbreak.ReorgE_ct;//0.156;
  double chargetransfertriebkraft = Config::get().exbreak.ct_triebkraft;//1.550;
  double rekombinationstriebkraft = Config::get().exbreak.rek_triebkraft;//-4.913;
  double rek_reorganisation = Config::get().exbreak.ReorgE_rek;//0.184;
  double oszillatorstrength = Config::get().exbreak.oscillatorstrength;//0.0852;
  double wellenzahl = Config::get().exbreak.wellenzahl;//28514.91;

  ////////////////////////////////////
  double pi = 3.141592654;
  double h_quer = 1 / (2 * pi)*4.135667662e-15;
  double boltzmann_konstante = 8.6173303e-5; //  in gauß einheiten
  double k_rad = wellenzahl*wellenzahl*oszillatorstrength; // fluoreszenz

  /////////////////////////////////// INPUT-READING
  ifstream schwerpunkt;
  schwerpunkt.open(masscenters);
  schwerpunkt >> gesamtanzahl;
  std::cout << "Read number of Monomers: " << gesamtanzahl << std::endl;
  if (gesamtanzahl == pscanzahl + nscanzahl) //test if correct number of molecules was given
  {  
    std::cout << "OK!" << std::endl;
  }
  else //gesamtanzahl != pscanzahl + nscanzahl
  { 
    std::cout << "WRONG!" << std::endl;
    return 0;
  }
  schwerpunkt >> skipline;
  std::vector <double> x(gesamtanzahl + 1), y(gesamtanzahl + 1), z(gesamtanzahl + 1);

  for (i = 1; i < (gesamtanzahl + 1); i++) 
  {
    schwerpunkt >> zeile >> x[i] >> y[i] >> z[i]; //reading and saving balance points for molecules
  }
  ///////////////////////////////////
  ifstream exciton;
  exciton.open(pscpairexrates);
  while (getline(exciton, zeile)) 
  { //counting of excitonpairs in homodimers
    e++;
  }
  exciton.close();

  std::vector <std::vector<double>> coupling_exciton(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1)), //in original code 2d-arrays were used
    coupling_ladung(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1)),
    coupling_ct(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1)),
    coupling_rek(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1)),
    coupling_fulleren(gesamtanzahl + 1, std::vector <double>(gesamtanzahl + 1));
  for (i = 1; i < (gesamtanzahl + 1); i++) //initializaton of the 2d vectors with 0 in all places
  { 
    for (j = 1; j < (gesamtanzahl + 1); j++) {
      coupling_exciton[i][j] = 0;
      coupling_ladung[i][j] = 0;
      coupling_ct[i][j] = 0;
      coupling_rek[i][j] = 0;
      coupling_fulleren[i][j] = 0;
    }
  }

  std::vector <int> exciton_1(e + 1), exciton_2(e + 1); //vectors for exciton pairs

  exciton.open(pscpairexrates);
  for (i = 1; i < (e + 1); i++) 
  {
    exciton >> exciton_1[i] >> exciton_2[i];
    exciton >> coupling_exciton[exciton_1[i]][exciton_2[i]];
    coupling_exciton[exciton_2[i]][exciton_1[i]] = coupling_exciton[exciton_1[i]][exciton_2[i]];
  }

  exciton.close();
  std::vector <int> partneranzahl(pscanzahl + nscanzahl + 1);
  for (i = 1; i < (pscanzahl + nscanzahl + 1); i++) 
  {
    partneranzahl[i] = 0;
  } //inizialisation of all elements with 0

  for (i = 1; i < (e + 1); i++) // counting homodimer partners for j
  {	
    for (j = 1; j < (pscanzahl + 1); j++) 
    {
      if ((exciton_1[i] == j) || (exciton_2[i] == j)) 
      {
        partneranzahl[j]++;
      }
    }
  }
  ////////////////////////////////////
  exciton.open(pscpairchrates);
  for (i = 1; i < (e + 1); i++) 
  {
    exciton >> j >> h; //j=exciton_1[i], h=exciton_2[i]
    exciton >> coupling_ladung[j][h];
    coupling_ladung[h][j] = coupling_ladung[j][h];
  }
  exciton.close();

  exciton.open(pnscpairrates);
  while (getline(exciton, zeile)) //counting of heterodimers
  { 
    het++;
  }
  exciton.close();

  std::vector <int> hetero_1(het + 1), hetero_2(het + 1);

  exciton.open(pnscpairrates);
  for (i = 1; i < (het + 1); i++) 
  {
    exciton >> hetero_1[i] >> hetero_2[i];
    exciton >> coupling_ct[hetero_1[i]][hetero_2[i]] >> coupling_rek[hetero_1[i]][hetero_2[i]];
    coupling_rek[hetero_2[i]][hetero_1[i]] = coupling_rek[hetero_1[i]][hetero_2[i]];
    coupling_ct[hetero_2[i]][hetero_1[i]] = coupling_ct[hetero_1[i]][hetero_2[i]];
  }
  for (i = 1; i < (het + 1); i++) //counting of heteropartners for j and adding to known number of partners
  {				
    for (j = 1; j < (gesamtanzahl + 1); j++) 
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
  while (getline(exciton, zeile)) //counting fullerene homopairs
  { 
    full++;
  }
  exciton.close();

  std::cout << "Number of n-semiconductor pairs " << full << std::endl;
  std::vector <int> fulleren_1(full + 1), fulleren_2(full + 1), test(2000);

  exciton.open(nscpairrates);
  for (i = 1; i < (full + 1); i++) 
  {
    exciton >> fulleren_1[i] >> fulleren_2[i];
    test[i] = fulleren_2[i];
    exciton >> coupling_fulleren[fulleren_1[i]][fulleren_2[i]];
    coupling_fulleren[fulleren_2[i]][fulleren_1[i]] = coupling_fulleren[fulleren_1[i]][fulleren_2[i]];
  }

  exciton.close();

  for (i = 1; i < (full + 1); i++) //counting of fullerenhomopartners and adding to known partners
  {				
    for (j = 1; j < (gesamtanzahl + 1); j++) 
    {
      if ((fulleren_1[i] == j) || (fulleren_2[i] == j)) 
      {
        partneranzahl[j]++;
      }
    }
  }

  std::vector<std::vector<int>> partner(gesamtanzahl + 1, std::vector<int>());//2D-vector with variing length for second vector

  for (i = 1; i < (gesamtanzahl + 1); i++) //dynamic allocation for length of 2nd vector
  { 
    partner[i].resize(partneranzahl[i] + 1);
  }

  for (i = 1; i < (gesamtanzahl + 1); i++) //initializing all elements of vector partner with 0
  { 
    for (j = 1; j < (partneranzahl[i] + 1); j++) {
      partner[i][j] = 0;
    }
  }

  for (i = 1; i < (gesamtanzahl + 1); i++) 
  {
    j = 1; //j is here the number of the partner to particle i 

    for (h = 1; h < (e + 1); h++)  //e = number of exciton-pairs [homodimer-pairs?]
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

    for (h = 1; h < (het + 1); h++) //het = number of heterodimer-pairs
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

    for (h = 1; h < (full + 1); h++) //full = number of fullerene homopairs
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
    if (partneranzahl[i] != j - 1) //after implementation into CAST replace by throw
    { 
      cout << "Error with number of partners for monomer " << i << std::endl;
    }
  }

  // INPUT-END

  ofstream kopplung;
  kopplung.open("partner.txt");
  for (i = 1; i < (gesamtanzahl + 1); i++) 
  {
    kopplung << setw(6) << i << setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess
    for (j = 1; j < (partneranzahl[i] + 1); j++) 
    {
      kopplung << setw(6) << partner[i][j]; //writes the indices of the partners j
    }
    kopplung << '\n';
  }
  kopplung.close();

  kopplung.open("couplings.txt");
  for (i = 1; i < (pscanzahl + 1); i++) 
  {
    kopplung << setw(6) << i << setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess

    for (j = 1; j < (partneranzahl[i] + 1); j++) 
    {
      kopplung << setw(6) << partner[i][j]; //writes the indices of the partners j
      if (partner[i][j] < (pscanzahl + 1)) 
      {
        kopplung << setw(12) << setprecision(6) << fixed << coupling_exciton[i][partner[i][j]]; //writes the exciton-coupling between i and j
      }
      else if (partner[i][j] > pscanzahl)
      {
        kopplung << setw(12) << setprecision(6) << fixed << coupling_ct[i][partner[i][j]];     //writes charge-transfer-coupling between i and j 
      }
    }
    kopplung << '\n';
  }

  for (i = (pscanzahl + 1); i < (gesamtanzahl + 1); i++)
  {
    kopplung << setw(6) << i << setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess
    for (j = 1; j < (partneranzahl[i] + 1); j++) 
    {
      kopplung << setw(6) << partner[i][j]; //writes the indices of the partners j
      if (partner[i][j] < (pscanzahl)) 
      {
        kopplung << setw(12) << setprecision(6) << fixed << coupling_rek[i][partner[i][j]]; //writes the recombination coupling between i and j
      }
      else if (partner[i][j] > pscanzahl)
      {
        kopplung << setw(12) << setprecision(6) << fixed << coupling_fulleren[i][partner[i][j]]; // writes some coupling regarding fullerens?
      }
    }
    kopplung << '\n';
  }
  kopplung.close();

  // Startpunkte bestimmen ##################################################################################################################
  double x_monomer(0.), y_monomer(0.), z_monomer(0.), x_fulleren(0.), y_fulleren(0.), z_fulleren(0.),
    x_gesamt(0.), y_gesamt(0.), z_gesamt(0.), x_mittel(0.), y_mittel(0.), z_mittel(0.);

  for (i = 1; i < (pscanzahl + 1); i++)
  {
    x_monomer += (x[i] / pscanzahl);
    y_monomer += (y[i] / pscanzahl);
    z_monomer += (z[i] / pscanzahl);
  }

  for (i = (pscanzahl + 1); i < (gesamtanzahl + 1); i++)     //using fact, that fullerens always have larger indices than other monomers
  { 
    x_fulleren += (x[i] / nscanzahl);
    y_fulleren += (y[i] / nscanzahl);
    z_fulleren += (z[i] / nscanzahl);
  }

  for (i = 1; i < (gesamtanzahl + 1); i++)
  {
    x_gesamt += (x[i] / gesamtanzahl);
    y_gesamt += (y[i] / gesamtanzahl);
    z_gesamt += (z[i] / gesamtanzahl);
  }

  x_mittel = (x_monomer + x_fulleren) / 2;
  y_mittel = (y_monomer + y_fulleren) / 2;
  z_mittel = (z_monomer + z_fulleren) / 2;

  ofstream interface;
  interface.open("massponts_general.xyz"); //writing out average balance points for all groupings of monomers
  interface << "4" << '\n' << '\n';
  interface << setw(5) << "X" << setw(12) << setprecision(6) << fixed << x_monomer << setw(12) << setprecision(6) << fixed << y_monomer << setw(12) << setprecision(6) << fixed << z_monomer << '\n';
  interface << setw(5) << "X" << setw(12) << setprecision(6) << fixed << x_fulleren << setw(12) << setprecision(6) << fixed << y_fulleren << setw(12) << setprecision(6) << fixed << z_fulleren << '\n';
  interface << setw(5) << "X" << setw(12) << setprecision(6) << fixed << x_gesamt << setw(12) << setprecision(6) << fixed << y_gesamt << setw(12) << setprecision(6) << fixed << z_gesamt << '\n';
  interface << setw(5) << "X" << setw(12) << setprecision(6) << fixed << x_mittel << setw(12) << setprecision(6) << fixed << y_mittel << setw(12) << setprecision(6) << fixed << z_mittel << '\n';
  interface.close();

  index = 0;
  max = 0;
  std::vector <int> startpunkt(pscanzahl + 1);

  switch (ebene) { //different cases for the possible planes of the interface

  case 'x':
    for (i = 1; i < (pscanzahl + 1); i++) //determining the maximal distance to interace
    { 
      if (x[i] > max) {
        max = x[i];
      }
    }
    std::cout << "Maxdistance is " << setw(12) << setprecision(6) << fixed << max << std::endl;

    for (i = 1; i < (pscanzahl + 1); i++)  //determining the necessary number of starting points? 
    {
      if ((x[i] - x_mittel) > (0.85*(max - x_mittel))) {
        index++;
        startpunkt[index] = i;
      }
    }
    break;

  case 'y':
    for (i = 1; i < (pscanzahl + 1); i++)  //determining the maximal distance to interface
    { 
      if (y[i] > max) {
        max = y[i];
      }
    }
    std::cout << "Maxdistance is " << setw(12) << setprecision(6) << fixed << max << std::endl;

    for (i = 1; i < (pscanzahl + 1); i++) //determining the necessary number of starting points? 
    { 
      if ((y[i] - y_mittel) > (0.85*(max - y_mittel))) {
        index++;
        startpunkt[index] = i;

      }
    }
    break;

  case 'z':
    for (i = 1; i < (pscanzahl + 1); i++) //determining the maximal distance to interace
    { 
      if (z[i] > max) {
        max = z[i];
      }
    }
    std::cout << "Maxdistance is " << setw(12) << setprecision(6) << fixed << max << std::endl;

    for (i = 1; i < (pscanzahl + 1); i++) //determining the necessary number of starting points? 
    { 
      if ((z[i] - z_mittel) > (0.85*(max - z_mittel))) 
      {
        index++;
        startpunkt[index] = i;

      }
    }
    break;
  }


  interface.open("startingpoint.xyz");
  interface << index << '\n' << '\n'; //writes the number of startingponts
  for (i = 1; i < (index + 1); i++) 	//writes coordinates of startingpoints
  {			
    interface << setw(5) << "X" << setw(12) << setprecision(6) << fixed << x[startpunkt[i]];
    interface << setw(12) << setprecision(6) << fixed << y[startpunkt[i]];
    interface << setw(12) << setprecision(6) << fixed << z[startpunkt[i]] << '\n';
  }
  interface.close();

  // ################################################################################## Beginn der Simulation ##############################################################################
  // Variablen
  double zeit(0.), zeit_1(0.), zeit_2(0.);

  // Schrittanzahl pro MC-Simulation
  int schritt = 2 * (e + full) + 400;
  std::cout << "Number of stepps " << schritt << std::endl;

  // ###################################################################################

  std::vector <std::vector<double>> vel_ex(index + 1, std::vector <double>(101)),
    vel_ch(index + 1, std::vector <double>(101)),
    zeit_ex(index + 1, std::vector <double>(101)),
    zeit_ch(index + 1, std::vector <double>(101));
  std::vector <int> ex_diss(index + 1), ch_diss(index + 1), rek(index + 1), trapping(index + 1), radiativ(index + 1);
  std::vector <std::vector<char>> zustand(index + 1, std::vector <char>(101));
  std::vector <int> punkt(schritt + 1), punkt_ladung(schritt + 1);
  double r_summe, r_i, zufall1, coulombenergy, dwelltime, r_summe_fulleren;

  for (i = 1; i < (index + 1); i++) //initializing the vectors with 0
  { 
    ex_diss[i] = 0;
    ch_diss[i] = 0;
    rek[i] = 0;
    radiativ[i] = 0;
    trapping[i] = 0;
    for (j = 1; j < 101; j++)
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
  for (i = 0; i < (schritt + 1); i++) 
  {
    punkt[i] = 0;
  }

  ofstream run;
  run.open("run.txt");

  for (k = 1; k < (index + 1); k++) // schleife über startpunkte "index durch 1 vertauscht"
  { 
    run << "k ist " << k << std::endl;
    ofstream wtf;
    wtf.open("wtf.txt");
    wtf << i << "  " << startpunkt[k] << std::endl;
    wtf.close();
    for (j = 1; j < 101; j++)   // schleife über durchläufe für den gleichen startpunkt " 101 durch 11 vertauscht"
    { 
      zeit = 0;
      zeit_1 = 0;
      zeit_2 = 0;
      punkt[0] = startpunkt[k];

      for (i = 1; i < (schritt + 1); i++)
      {

        if (zustand[k][j] == 'c')
        {
          // site energies berechnen
          //########## raten addieren für monomere ##########
          r_summe = 0;
          random_device rd;
          default_random_engine engine(rd());
          normal_distribution<double> distribution0(0.0, 0.068584577); // hier neue standardabweichung eintragen
          zufall1 = distribution0(engine); //generating normal-distributed random number

          std::vector<double> raten(partneranzahl[punkt_ladung[i - 1]] + 1);
          for (h = 0; h < (partneranzahl[punkt_ladung[i - 1]] + 1); h++)
          {
            if (partner[punkt_ladung[i - 1]][h] < (pscanzahl + 1))
            {
              zufall = distribution0(engine);
              coulombenergy = coulomb(x, y, z, punkt[i - 1], partner[punkt_ladung[i - 1]][h], 3.4088) - coulomb(x, y, z, punkt[i - 1], punkt_ladung[i - 1], 3.4088);
              r_summe = r_summe + rate(coupling_ladung[punkt_ladung[i - 1]][partner[punkt_ladung[i - 1]][h]], ((zufall - zufall1) + coulombenergy), reorganisationsenergie_ladung);

              raten[h] = r_summe;
            }
            if ((partner[punkt_ladung[i - 1]][h] > (pscanzahl)) && (partner[punkt_ladung[i - 1]][h] == punkt[i - 1]))
            {
              zufall = distribution0(engine);

              // coulomb energie berechnen	   
              coulombenergy = coulomb(x, y, z, punkt_ladung[i - 1], partner[punkt_ladung[i - 1]][h], 1);

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
          r_summe_fulleren = 0;
          zufall1 = distribution0(engine);
          std::vector <double> raten_fulleren(partneranzahl[punkt[i - 1]] + 1);
          for (h = 1; h < (partneranzahl[punkt[i - 1]] + 1); h++)
          {
            if (partner[punkt[i - 1]][h] > (pscanzahl))
            {
              zufall = distribution0(engine);
              coulombenergy = coulomb(x, y, z, punkt_ladung[i - 1], partner[punkt[i - 1]][h], 3.4088) - coulomb(x, y, z, punkt[i - 1], punkt_ladung[i - 1], 3.4088);
              r_summe_fulleren = r_summe_fulleren + rate(coupling_fulleren[punkt[i - 1]][partner[punkt[i - 1]][h]], ((zufall - zufall1) + coulombenergy), fullerenreorganisationsenergie);

              raten_fulleren[h] = r_summe_fulleren;
            }
            if ((partner[punkt[i - 1]][h] < (pscanzahl + 1)) && (partner[punkt[i - 1]][h] == punkt_ladung[i - 1]))
            {
              zufall = distribution0(engine);

              // coulomb energie berechnen
              coulombenergy = coulomb(x, y, z, punkt[i - 1], partner[punkt[i - 1]][h], 1);

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
            if ((1 / r_summe - zeit_1) > 0)
            {
              zeit = zeit + (1 / r_summe - zeit_1);
              zeit_2 = zeit_2 + (1 / r_summe - zeit_1);
              zeit_1 = 0;
            }
            else if ((1 / r_summe - zeit_1) < 0)
            {
              zeit = zeit;
              zeit_2 = zeit_2;
              zeit_1 = 0;
            }
            else {

              run << "ERROR!" << std::endl;
              return 0;
            }
            // monomerhüpfen ausführen
            uniform_real_distribution<double> distribution1(0, 1);
            zufall = distribution1(engine);
            r_i = zufall*r_summe;
            for (g = 1; g < (partneranzahl[punkt_ladung[i - 1]] + 1); g++)
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
                run << "old Monomer " << setw(5) << punkt_ladung[i - 1] << std::endl;
                run << "new Monomer " << setw(5) << punkt_ladung[i] << std::endl;
                run << "Coupling " << setw(12) << setprecision(6) << fixed << coupling_ladung[punkt_ladung[i - 1]][punkt_ladung[i]] << std::endl;
                run << "Fulleren " << setw(5) << punkt[i] << setw(5) << punkt[i - 1] << std::endl;
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
              return 0;
            }

            // fullerenhüpfen ausführen
            uniform_real_distribution<double> distribution1(0, 1);
            zufall = distribution1(engine);
            r_i = zufall*r_summe_fulleren;

            for (g = 1; g < ((partneranzahl[punkt[i - 1]]) + 1); g++)
            {
              if ((raten_fulleren[g] > r_i) && ((partner[punkt[i - 1]][g]) > pscanzahl))
              {
                run << "Chargetransfer in Fullerenephase" << std::endl;

                punkt[i] = partner[punkt[i - 1]][g];
                punkt_ladung[i] = punkt_ladung[i - 1];
                run << "old Fulleren " << setw(5) << punkt[i - 1] << std::endl;
                run << "new Fulleren " << setw(5) << punkt[i] << std::endl;
                run << "Monomer " << setw(5) << punkt_ladung[i] << std::endl;
                run << "Coupling " << setw(12) << setprecision(6) << coupling_fulleren[punkt[i]][punkt[i - 1]] << std::endl;
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
                return 0;
              }
            }
          } // endes des hüpfenden fullerens
        } // ende des 'c'-zustands
      //__________________________________________________________________________________________________________

        else if (zustand[k][j] == 'e')
        {
          // site energies berechnen
          r_summe = 0;
          random_device rd;
          default_random_engine engine(rd());
          normal_distribution<double> distribution0(0.0, 0.0338987); // hier neue standardabweichung eintragen
          zufall1 = distribution0(engine); //generating an normal-distributed random number

          std::vector <double> raten(partneranzahl[punkt[i - 1]] + 1);

          for (h = 1; h < (partneranzahl[punkt[i - 1]] + 1); h++)
          {
            if (partner[punkt[i - 1]][h] < (pscanzahl + 1))
            {
              zufall = distribution0(engine);// generatinjg a second normal distributed random number

              r_summe += rate(coupling_exciton[punkt[i - 1]][partner[punkt[i - 1]][h]], (zufall - zufall1), reorganisationsenergie_exciton);
              raten[h] = r_summe;
            }

            else if (partner[punkt[i - 1]][h] > (pscanzahl))
            {
              zufall = distribution0(engine);
              // coulomb energie berechnen
              coulombenergy = coulomb(x, y, z, punkt[i - 1], partner[punkt[i - 1]][h], 1);
              r_summe += rate(coupling_ct[punkt[i - 1]][partner[punkt[i - 1]][h]], (zufall - zufall1) + chargetransfertriebkraft + coulombenergy, ct_reorganisation);
              raten[h] = r_summe;
            }
          } // end of h

          // fluoreszenz dazuaddieren
          r_summe = r_summe + k_rad;
          dwelltime = 1 / r_summe;

          // schritt bestimmen
          uniform_real_distribution<double> distribution1(0, 1);
          zufall = distribution1(engine);
          r_i = zufall*r_summe;
          zeit = zeit + 1 / r_summe;

          //falls trapping
          zufall = distribution1(engine);
          if (zufall*(900e-1 + 1 / r_summe) > (900e-1))
          {
            run << "Exziton trapped!" << std::endl;
            trapping[k]++;
            zustand[k][j] = 't';
            break;
          }

          for (g = 1; g < (partneranzahl[punkt[i - 1]] + 1); g++)
          {
            ofstream wtf;
            wtf.open("wtf.txt");
            wtf << g << "  " << raten[g] << "  " << r_i << std::endl;
            wtf.close();
            if (raten[g] > r_i)
            {
              punkt[i] = partner[punkt[i - 1]][g];
              run << punkt[i] << std::endl;

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
                run << "Exzitonspeed " << vel_ex[k][j] * 1e-9 << std::endl;
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
          return 0;
        } //end of undefined zustand
      //___________________________________________________________________________________________________

      } // ende über schleife i
    } // Ende über Schleife über durchläufe für den gleichen startpunkt j
  } // Ende über Schleife der Startpunkte k

  run.close();

  // Auswertung: Prozentsätze
  ofstream auswertung;
  auswertung.open("evaluation.txt");
  auswertung << setw(4) << "k" << setw(4) << "IX" << setw(9) << "Ex_Diss." << setw(9) << "Ch_Diss" << setw(9) << "Rek." << setw(9) << "Trapp." << setw(9) << "Fluor." << '\n';

  double ex_diss_efficiency, ch_diss_efficiency, rek_efficiency, trapp_efficiency, rad_efficiency;
  double mittel_ex = 0;
  double mittel_ch = 0;
  double mittel_rek = 0;
  double mittel_trapp = 0;
  double mittel_rad = 0;
  std::vector <double> mittel_ex_vel(index + 1), mittel_ch_vel(index + 1), standard_ex(index + 1), standard_ch(index + 1);
  int zahl;

  for (k = 1; k < (index + 1); k++) {
    mittel_ex = mittel_ex + (ex_diss[k] * 1.0);
    mittel_ch = mittel_ch + (ch_diss[k] * 1.0);
    mittel_rek = mittel_rek + (rek[k] * 1.0);
    mittel_trapp = mittel_trapp + (trapping[k] * 1.0);
    mittel_rad = mittel_rad + (radiativ[k] * 1.0);
    ex_diss_efficiency = ex_diss[k];
    ch_diss_efficiency = ch_diss[k];
    rek_efficiency = rek[k];
    trapp_efficiency = trapping[k];
    rad_efficiency = radiativ[k];
    auswertung << setw(4) << k << setw(4) << startpunkt[k] << setw(9) << setprecision(5) << ex_diss_efficiency;
    auswertung << setw(9) << setprecision(5) << ch_diss_efficiency;
    auswertung << setw(9) << setprecision(5) << rek_efficiency << setw(9) << setprecision(5) << trapp_efficiency << setw(9) << setprecision(5) << rad_efficiency << '\n';
  }

  auswertung << setw(9) << "Average " << setw(9) << setprecision(5) << fixed << mittel_ex / index << setw(9) << setprecision(5) << fixed << mittel_ch / index;
  auswertung << setw(9) << setprecision(5) << fixed << mittel_rek / index << setw(9) << setprecision(5) << fixed << mittel_trapp / index;
  auswertung << setw(9) << setprecision(5) << fixed << mittel_rad / index << '\n';
  auswertung << "Velocities" << '\n';
  auswertung << setw(4) << "k" << setw(5) << "IX" << setw(11) << "Ex_vel" << setw(11) << "Ex_s_dev" << setw(11) << "Ch_vel" << setw(11) << "Ch_s_dev" << '\n';

  // mittlere Geschwindigkeit
  for (k = 1; k < (index + 1); k++) {
    mittel_ch_vel[k] = 0;
    mittel_ex_vel[k] = 0;

    // Mittelwert berechnen
    for (j = 1; j < 101; j++) {
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
    for (j = 1; j < 101; j++) {
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
    auswertung << setw(4) << k << setw(5) << startpunkt[k] << setw(11) << setprecision(5) << fixed << mittel_ex_vel[k] * 1e-9;
    auswertung << setw(11) << setprecision(5) << fixed << standard_ex[k] * 1e-9;
    auswertung << setw(11) << setprecision(5) << fixed << mittel_ch_vel[k] * 1e-9;
    auswertung << setw(11) << setprecision(5) << fixed << standard_ch[k] * 1e-9 << '\n';
  }
  double mittelwert_geschw_exciton = 0;
  double mittelwert_geschw_ladung = 0;
  for (k = 1; k < (index + 1); k++) {
    mittelwert_geschw_exciton = mittelwert_geschw_exciton + mittel_ex_vel[k];
    mittelwert_geschw_ladung = mittelwert_geschw_ladung + mittel_ch_vel[k];
  }
  auswertung << left << setw(7) << " Average    " << left << setw(22) << setprecision(5) << fixed << mittelwert_geschw_exciton / index*1e-9;
  auswertung << left << setw(9) << setprecision(5) << fixed << mittelwert_geschw_ladung / index*1e-9 << '\n';

  // Verteilung Ladungen und Exzitonengeschwindigkeiten
  ofstream exciton_verteilung;
  exciton_verteilung.open("exciton_distribution.txt");
  for (i = 1; i < 21; i++) {
    zahl = 0;
    for (k = 1; k < (index + 1); k++) {
      for (j = 1; j < 101; j++) {
        if ((vel_ex[k][j] > (i * 50 * 1e9))) {
          zahl++;
        }
      }
    }
    exciton_verteilung << setw(9) << setprecision(5) << i * 50 << setw(9) << zahl / index << '\n';
  }

  exciton_verteilung.close();
  exciton_verteilung.open("charge_distribution.txt");
  for (i = 1; i < 21; i++) {
    zahl = 0;
    for (k = 1; k < (index + 1); k++) {
      for (j = 1; j < 101; j++) {
        if ((vel_ch[k][j] > (i * 50 * 1e9))) {
          zahl++;
        }
      }
    }
    exciton_verteilung << setw(9) << setprecision(5) << i * 50 << setw(9) << zahl / index << '\n';
  }
  exciton_verteilung.close();

  return 0;
}

