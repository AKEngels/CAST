#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<limits>
#include<string>
#include<random>
#include<vector>

//Original code by Charlotte. German comments are from original code, english comments were inserted during implementation into CAST.
//In original code loops begun at 1 and ran to x+1. Also many arrays were replaced by vectors.

using namespace std;

istream& skipline(istream& in )
{
    return in.ignore( numeric_limits < streamsize >::max(), '\n' );
}

// Definition der length-berechnung
double length(vector <double> arr1, vector <double> arr2, vector <double> arr3, int p, int q)
{
  double l=1;
  l=sqrt((arr1[p]-arr1[q])*(arr1[p]-arr1[q])+(arr2[p]-arr2[q])*(arr2[p]-arr2[q])+(arr3[p]-arr3[q])*(arr3[p]-arr3[q]));  
  return l;
}

double rate(double coupling, double energy, double reorganisation)
{
  double l=1;
  double pi=3.141592654;
  double h_quer=1/(2*pi)*4.135667662e-15;
  double boltzmann_konstante=8.6173303e-5; //  in gauß einheiten
  l=(coupling*coupling)/h_quer*sqrt(pi/(reorganisation*boltzmann_konstante*298))*exp(-(reorganisation+energy)*(reorganisation+energy)/(4*boltzmann_konstante*298*reorganisation));
  return l;
}

double coulomb(vector<double> arr1, vector<double> arr2, vector<double> arr3, int p, int q, double e_relative){
  double l=1;
  l=sqrt((arr1[p]-arr1[q])*(arr1[p]-arr1[q])+(arr2[p]-arr2[q])*(arr2[p]-arr2[q])+(arr3[p]-arr3[q])*(arr3[p]-arr3[q]));  
  double pi=3.141592654;
  double e_0=8.854187e-12;
  double elementar=1.60217662e-19;
  double c=1;
  c=-elementar/(4*pi*e_0*e_relative*l*1e-10);
  return c;
}

int main()
{
	string zeile;
	int i, j, h, index, k, g;
	int monomer_anzahl, monomer_atom, full_anzahl, gesamtanzahl;
	char ebene;
	int e = 0;
	int het = 0;
	int full = 0;
	double max, zufall;
	//////////////////// Definition der Konstanten  
	double reorganisationsenergie_exciton, reorganisationsenergie_ladung, fullerenreorganisationsenergie, trappingrate, chargetransfertriebkraft, rekombinationstriebkraft, temperatur, ct_reorganisation, rek_reorganisation;
	double oszillatorstrength, wellenzahl;
	////////////////////////////////////// hier konstanten angeben ///////////////////////////////////////////////
	reorganisationsenergie_exciton = 0.561; // SCS-CC2 wert: vorher 0.561
	reorganisationsenergie_ladung = 0.194;
	fullerenreorganisationsenergie = 0.178;
	ct_reorganisation = 0.156;
	trappingrate = 0.0000005;
	chargetransfertriebkraft = 1.550;
	rekombinationstriebkraft = -4.913;
	rek_reorganisation = 0.184;
	temperatur = 298;
	oszillatorstrength = 0.0852;
	wellenzahl = 28514.91;
	////////////////////////////////////
	double pi = 3.141592654;
	double h_quer = 1 / (2 * pi)*4.135667662e-15;
	double boltzmann_konstante = 8.6173303e-5; //  in gauß einheiten
	double k_rad = wellenzahl*wellenzahl*oszillatorstrength; // fluoreszenz
	/////////////////////////////////// INPUT-MANUELL
	std::cout << "Anzahl der Monomere" << endl;
	cin >> monomer_anzahl;
	std::cout << "Anzahl der Monomeratome" << endl;
	cin >> monomer_atom;
	std::cout << "Anzahl der Fullerene" << endl;
	cin >> full_anzahl;
	std::cout << "Ebene " << endl;
	cin >> ebene;
	/////////////////////////////////// INPUT-READING
	ifstream schwerpunkt;
	schwerpunkt.open("schwerpunkt.xyz");
	schwerpunkt >> gesamtanzahl;
	std::cout << "Eingelesen " << gesamtanzahl << endl;
	if (gesamtanzahl == monomer_anzahl + full_anzahl) {  //test if correct number of molecules was given
		std::cout << "OK!" << endl;
	}
	else { //gesamtanzahl != monomer_anzahl + full_anzahl
		std::cout << "FALSCH!" << endl;
		return 0;
	}
	schwerpunkt >> skipline;
	std::vector <double> x(gesamtanzahl), y(gesamtanzahl), z(gesamtanzahl);
	for (i = 0; i < (gesamtanzahl); i++) {
		schwerpunkt >> zeile >> x[i] >> y[i] >> z[i]; //reading and saving balance point-coordinates for molecules
	}
	///////////////////////////////////
	ifstream exciton;
	exciton.open("homodimer_exciton.txt");
	while (getline(exciton, zeile)) { //counting of excitonpairs in homodimers
		e++;
	}
	exciton.close();
	//cout << e << " Exzitonenpaare." << endl;
	std::vector <vector<double>> coupling_exciton(gesamtanzahl, std::vector <double>(gesamtanzahl)), //in original code 2d-arrays were used
		coupling_ladung(gesamtanzahl, std::vector <double>(gesamtanzahl)),							 
		coupling_ct(gesamtanzahl, std::vector <double>(gesamtanzahl)),
		coupling_rek(gesamtanzahl, std::vector <double>(gesamtanzahl)),
		coupling_fulleren(gesamtanzahl, std::vector <double>(gesamtanzahl));
	for (i = 0; i < (gesamtanzahl); i++) { //initializaton of the 2d vectors with 0 in all places
		for (j = 0; j < (gesamtanzahl); j++) {
			coupling_exciton[i][j] = 0;
			coupling_ladung[i][j] = 0;
			coupling_ct[i][j] = 0;
			coupling_rek[i][j] = 0;
			coupling_fulleren[i][j] = 0;
		}
	}

	std::vector <int> exciton_1(e), exciton_2(e); //vectors for exciton pairs

	exciton.open("homodimer_exciton.txt");
	for (i = 0; i < e; i++) {
		exciton >> exciton_1[i] >> exciton_2[i]; 
		exciton >> coupling_exciton[exciton_1[i]][exciton_2[i]];
		coupling_exciton[exciton_2[i]][exciton_1[i]] = coupling_exciton[exciton_1[i]][exciton_2[i]];
	}
	exciton.close();
	std::vector <int> partneranzahl(monomer_anzahl + full_anzahl);
	for (i = 0; i < (monomer_anzahl + full_anzahl); i++) {
		partneranzahl[i] = 0;
	} //inizialisation of all elements with 0

	for (i = 0; i < e; i++) {	// counting homodimer partners for j
		for (j = 0; j < monomer_anzahl; j++) {
			if ((exciton_1[i] == j) || (exciton_2[i] == j)) {
				partneranzahl[j]++;
			}
		}
	}
	////////////////////////////////////
	exciton.open("homodimer_ladung.txt");
	for (i = 0; i < e; i++) {
		exciton >> j >> h; //j=exciton_1[i], h=exciton_2[i]
		exciton >> coupling_ladung[j][h];
		coupling_ladung[h][j] = coupling_ladung[j][h];
	}
	exciton.close();

	exciton.open("heterodimer.txt");
	while (getline(exciton, zeile)) { //counting of heterodimers
		het++;
	}
	exciton.close();

	std::vector <int> hetero_1(het), hetero_2(het);

	exciton.open("heterodimer.txt");
	for (i = 0; i < het; i++) {
		exciton >> hetero_1[i] >> hetero_2[i];
		exciton >> coupling_ct[hetero_1[i]][hetero_2[i]] >> coupling_rek[hetero_1[i]][hetero_2[i]];
		coupling_rek[hetero_2[i]][hetero_1[i]] = coupling_rek[hetero_1[i]][hetero_2[i]];
		coupling_ct[hetero_2[i]][hetero_1[i]] = coupling_ct[hetero_1[i]][hetero_2[i]];
	}
	for (i = 0; i < het; i++) {				//counting of heteropartners for j and adding to known number of partners
		for (j = 0; j < gesamtanzahl; j++) {
			if ((hetero_1[i] == j) || (hetero_2[i] == j)) {
				partneranzahl[j]++;
			}
		}
	}
	exciton.close();

	////////////////////////////////////////
	exciton.open("fulleren.txt");
	while (getline(exciton, zeile)) { //counting fullerene homopairs
		full++;
	}
	exciton.close();

	std::cout << "Anzahl an Fullerenpaaren " << full << endl;
	std::vector <int> fulleren_1(full), fulleren_2(full), test(full);

	exciton.open("fulleren.txt");
	for (i = 0; i < full; i++) {
		exciton >> fulleren_1[i] >> fulleren_2[i];
		test[i] = fulleren_2[i];
		exciton >> coupling_fulleren[fulleren_1[i]][fulleren_2[i]];
		coupling_fulleren[fulleren_2[i]][fulleren_1[i]] = coupling_fulleren[fulleren_1[i]][fulleren_2[i]];
	}
	//for (i = 0; i < full; i++) { //??? why? above test[i] = fulleren_2[i] no changes to fulleren_2 appear in between
	//	fulleren_2[i] = test[i];
	//}
	exciton.close();

	for (i = 0; i < full; i++) {				//counting of fullerenhomopartners and adding to known partners
		for (j = 0; j < gesamtanzahl; j++) {
			if ((fulleren_1[i] == j) || (fulleren_2[i] == j)) {
				partneranzahl[j]++;
			}
		}
	}

	std::vector<vector<int>> partner(gesamtanzahl, std::vector<int>());//2D-vector with variing length for second vector

for (i=0;i<gesamtanzahl;i++){ //dynamic allocation for length of 2nd vector
  partner[i].resize(partneranzahl[i]+1);
}

for (i=0;i<gesamtanzahl;i++){ //initializing all elements of vector partner with 0
  for (j=0;j<partneranzahl[i];j++){
    partner[i][j]=0;
  }
}

for (i=0;i<gesamtanzahl;i++){
  j=0; //j is here the number of the partner to particle i 
	
  for (h=0;h<e;h++){ //e = number of exciton-pairs [homodimer-pairs?]
		if (exciton_1[h]==i){
		 partner[i][j]=exciton_2[h];
		 j++; //j is always incremented when an element is added to the list of partners of particle i
		}
		if (exciton_2[h]==i){
		 partner[i][j]=exciton_1[h];
		 j++;
		}
	}

	for (h=0;h<het;h++){ //het = number of heterodimer-pairs
	 if (hetero_1[h]==i){
      partner[i][j]=hetero_2[h];
      j++;
	 }
	 if (hetero_2[h]==i){
      partner[i][j]=hetero_1[h];
      j++;
     }
	}

   for (h=0;h<full;h++){ //full = number of fullerene homopairs
    if (fulleren_1[h]==i){
      partner[i][j]=fulleren_2[h];
//      cout << fulleren_2[h] << endl;      
      j++;
    }
    if (fulleren_2[h]==i){
      partner[i][j]=fulleren_1[h];
//      cout << fulleren_1[h] << endl;
      j++; // since the 2nd dimension length of vector partner was set to partneranzahl[i] in thes logic construction j must always end up to be equal to partneranzahl[i]
    }
  }  
   if (partneranzahl[i] != j ) { //after implementation into CAST replace by throw
	   cout << "Fehler bei Partneranzahl für Monomer:" << i << endl;
  }
}

// INPUT-END

ofstream kopplung;
kopplung.open("partner.txt");
for (i=0;i<gesamtanzahl;i++){
  kopplung << setw(6) << i << setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess
  for (j=0;j<partneranzahl[i];j++){
    kopplung << setw(6) << partner[i][j]; //writes the indices of the partners j
  }
  kopplung << '\n';
}
kopplung.close();

kopplung.open("couplings.txt");
for (i=0;i<monomer_anzahl;i++){
  kopplung << setw(6) << i << setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess

  for (j=0;j<partneranzahl[i];j++){
    kopplung << setw(6) << partner[i][j]; //writes the indices of the partners j
    if (partner[i][j] < monomer_anzahl){
      kopplung << setw(12) << setprecision(6) << fixed << coupling_exciton[i] [partner[i][j]]; //writes the exciton-coupling between i and j
    }
    else if (partner[i][j]>monomer_anzahl){
      kopplung << setw(12) << setprecision(6) << fixed << coupling_ct[i][ partner[i][j] ];     //writes charge-transfer-coupling between i and j 
    }
  }
  kopplung << '\n';  
}

for (i=monomer_anzahl;i<gesamtanzahl;i++){
  kopplung << setw(6) << i << setw(6) << partneranzahl[i]; //writes the indices of the molecules and the ammount of partners they posess
  for (j=0;j<partneranzahl[i];j++){
    kopplung << setw(6) << partner[i][j]; //writes the indices of the partners j
    if (partner[i][j]<(monomer_anzahl)){
      kopplung << setw(12) << setprecision(6) << fixed << coupling_rek[i][ partner[i][j] ]; //writes the recombination coupling between i and j
    }
    else if (partner[i][j]> monomer_anzahl){
      kopplung << setw(12) << setprecision(6) << fixed << coupling_fulleren[i][ partner[i][j] ]; // writes some coupling regarding fullerens?
    }
  }
  kopplung << '\n';
}
kopplung.close();

// Startpunkte bestimmen ##################################################################################################################
double x_monomer(0.), y_monomer(0.), z_monomer(0.), x_fulleren(0.), y_fulleren(0.), z_fulleren(0.), 
	   x_gesamt(0.), y_gesamt(0.), z_gesamt(0.), x_mittel(0.), y_mittel(0.), z_mittel(0.);

for (i=0;i<monomer_anzahl;i++){
  x_monomer+=(x[i]/monomer_anzahl); 
  y_monomer+=(y[i]/monomer_anzahl);
  z_monomer+=(z[i]/monomer_anzahl);
}

for (i=monomer_anzahl;i<gesamtanzahl;i++){ //using fact, that fullerens always have larger indices than other monomers
  x_fulleren+=(x[i]/full_anzahl);
  y_fulleren+=(y[i]/full_anzahl);
  z_fulleren+=(z[i]/full_anzahl);
}

for (i=0;i<gesamtanzahl;i++){
  x_gesamt+=(x[i]/gesamtanzahl);
  y_gesamt+=(y[i]/gesamtanzahl);
  z_gesamt+=(z[i]/gesamtanzahl);
}

x_mittel=(x_monomer+x_fulleren)/2;
y_mittel=(y_monomer+y_fulleren)/2;
z_mittel=(z_monomer+z_fulleren)/2;

ofstream interface;
interface.open("schwerpunkt_allgemein.xyz"); //writing out average balance points for all groupings of monomers
interface << "4" << '\n' << '\n';
interface << setw(5) << "X" << setw(12) << setprecision(6) << fixed << x_monomer << setw(12) << setprecision(6) << fixed << y_monomer << setw(12) << setprecision(6) << fixed << z_monomer << '\n';
interface << setw(5) << "X" << setw(12) << setprecision(6) << fixed << x_fulleren << setw(12) << setprecision(6) << fixed << y_fulleren << setw(12) << setprecision(6) << fixed << z_fulleren << '\n';
interface << setw(5) << "X" << setw(12) << setprecision(6) << fixed << x_gesamt << setw(12) << setprecision(6) << fixed << y_gesamt << setw(12) << setprecision(6) << fixed << z_gesamt << '\n';
interface << setw(5) << "X" << setw(12) << setprecision(6) << fixed << x_mittel << setw(12) << setprecision(6) << fixed << y_mittel << setw(12) << setprecision(6) << fixed << z_mittel << '\n';
interface.close();

index=0;
max=0;
std::vector <int> startpunkt(monomer_anzahl);

switch (ebene){ //different cases for the possible planes of the interface

  case 'x':
    for (i=0;i<monomer_anzahl;i++){ //determining the maximal distance to interace
      if (x[i]>max){
	max=x[i];
      }
    }
	std::cout << "Maxabstand ist " << setw(12) << setprecision(6) << fixed << max << endl;

    for (i=0;i<monomer_anzahl;i++){ //determining the necessary number of starting points? 
      if ((x[i]-x_mittel)>(0.85*(max-x_mittel))){
	startpunkt[index]=i;
	index++;
      }
    }
    break;

  case 'y':
    for (i=0;i<monomer_anzahl;i++){ //determining the maximal distance to interace
      if (y[i]>max){
	max=y[i];
      }
    }
	std::cout << "Maxabstand ist " << setw(12) << setprecision(6) << fixed << max << endl;
	std::cout << "Kriterium ist " << setw(12) << setprecision(6) << fixed << (max-y_mittel) << endl; //why only for y-plane?

    for (i=0;i<monomer_anzahl;i++){ //determining the necessary number of starting points? 
      if ((y[i]-y_mittel)>(0.85*(max-y_mittel))){ 
	startpunkt[index]=i;
	index++;
      }
    }    
    break;

  case 'z':
    for (i=0;i<monomer_anzahl;i++){ //determining the maximal distance to interace
      if (z[i]>max){
	max=z[i];
      }
    }
	std::cout << "Maxabstand ist " << setw(12) << setprecision(6) << fixed << max << endl;

    for (i=0;i<monomer_anzahl;i++){ //determining the necessary number of starting points? 
      if ((z[i]-z_mittel)>(0.85*(max-z_mittel))){
	startpunkt[index]=i;
	index++;
      }
    }    
    break;
}


interface.open("startpunkt.xyz");
interface << index << '\n' << '\n'; //writes the number of startingponts
for (i=0;i<index;i++){				//writes coordinates of startingpoints
  interface << setw(5) << "X" << setw(12) << setprecision(6) << fixed << x[ startpunkt[i] ];
  interface << setw(12) << setprecision(6) << fixed << y[ startpunkt[i] ];  
  interface << setw(12) << setprecision(6) << fixed << z[ startpunkt[i] ] << '\n';    
}
interface.close();

// ################################################################################## Beginn der Simulation ##############################################################################
// Variablen
double zeit(0.), zeit_1(0.), zeit_2(0.);

// Schrittanzahl pro MC-Simulation
int schritt=2*(e+full)+400;
std::cout << "Anzahl der Schritte " << schritt << endl;

// ###################################################################################

std::vector <vector<double>> vel_ex (index, std::vector <double> (101)), 
							vel_ch(index, std::vector <double>(101)),
							zeit_ex(index, std::vector <double>(101)),
							zeit_ch(index, std::vector <double>(101));
std::vector <int> ex_diss(index), ch_diss(index), rek(index), trapping(index), radiativ(index);
std::vector <vector<char>> zustand(index, std::vector <char> (101));
std::vector <int> punkt(schritt+1), punkt_ladung(schritt);
double r_summe, r_i, zufall1, coulombenergy, dwelltime, r_summe_fulleren;


for (i=0;i<index;i++){ //initializing the vectors with 0
    ex_diss[i]=0;
    ch_diss[i]=0;
    rek[i]=0;
    radiativ[i]=0;
    trapping[i]=0;  
  for (j=0;j<100;j++){
    vel_ex[i][j]=0;
    vel_ch[i][j]=0;
    zeit_ex[i][j]=0;
    zeit_ch[i][j]=0;
    zustand[i][j]='e';
  }
}

// k: index für startpunkte
// j: index für durchläufe
// i: index für schritt
for (i=0;i<schritt;i++){
  punkt[i]=0;
}

ofstream run; //outputs to control calculation
	run.open("run.txt");
ofstream raten_out;
	raten_out.open("raten.txt");

for (k=0;k<index;k++){ // schleife über startpunkte "index durch 1 vertauscht"
	run << "k ist " << k << endl;

	

for (j=0;j<100;j++){ // schleife über durchläufe für den gleichen startpunkt " 101 durch 11 vertauscht"
//run << "j ist " << j << endl;
zeit=0;
zeit_1=0;
zeit_2=0;

punkt[0]=startpunkt[k];

//cout << "Punkt 0 " << punkt[0] << endl;
   
int vschritt = 0;                                 //=schritt-1;//additional variable to access the preceding step set to maximum value for first step after that i-1

  for (i=0;i<schritt;i++){ 
//run << "i ist " << i << endl;

//______________________________________________________________________________________________________________________________________________________________________________________________

    if (zustand[k][j]=='c'){
      // site energies berechnen
      //############################################################################################################ raten addieren für monomere ##############################################################
      r_summe=0;
      random_device rd;
      default_random_engine engine(rd() );
      normal_distribution<double> distribution0(0.0,0.068584577); //##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! hier neue standardabweichung eintragen
      zufall1 = distribution0(engine); //generating normal-distributed random number

      std::vector<double> raten( partneranzahl [punkt_ladung[vschritt]]);
      for (h=0;h<partneranzahl [punkt_ladung[vschritt]];h++){
	if (partner [punkt_ladung[vschritt]] [h]<(monomer_anzahl+1)){
          zufall=distribution0(engine);	 
	  coulombenergy=coulomb( x, y, z, punkt[vschritt], partner [punkt_ladung[vschritt] ][h],3.4088)-coulomb( x, y, z, punkt[vschritt], punkt_ladung[vschritt], 3.4088);
	  r_summe=r_summe+rate(coupling_ladung[ punkt_ladung[vschritt] ][ partner[ punkt_ladung[vschritt] ][h] ], ((zufall-zufall1)+coulombenergy), reorganisationsenergie_ladung);
//	  cout << "LADUNGEN" << endl;
//	  cout << "Rate " << setw(12) << setprecision(3) << scientific << rate(coupling_ladung[ punkt_ladung[vschritt] ][ partner[ punkt_ladung[vschritt] ][h] ], (zufall-zufall1), reorganisationsenergie_ladung) << endl;
//	  cout << "Energie " << setw(12) << setprecision(6) << fixed << (zufall-zufall1)+coulombenergy << endl;
//	  cout << "Kopplung " << setw(12) << setprecision(6) << fixed << coupling_ladung[ punkt_ladung[vschritt] ][ partner[ punkt_ladung[vschritt] ][h] ] << endl;
	  raten[h]=r_summe;
	}
	if ((partner[ punkt_ladung[vschritt] ][h]>(monomer_anzahl))&&(partner[ punkt_ladung[vschritt] ][h]==punkt[vschritt])){
          zufall=distribution0(engine);	 
	  // coulomb energie berechnen	   
	  coulombenergy=coulomb( x, y, z, punkt_ladung[vschritt], partner[ punkt_ladung[vschritt] ][h], 1);
//	  cout << "LADUNGEN-HETERODIMER" << endl;	  
	  r_summe=r_summe+rate(coupling_rek[ punkt_ladung[vschritt] ][ partner[ punkt_ladung[vschritt] ][h] ], (zufall-zufall1)+rekombinationstriebkraft-coulombenergy, rek_reorganisation);
//	  cout << "Rate " << setw(12) << setprecision(3) << scientific  << rate(coupling_rek[ punkt_ladung[vschritt] ][ partner[ punkt_ladung[vschritt] ][h] ], (zufall-zufall1)+rekombinationstriebkraft-coulombenergy, rek_reorganisation) << endl;
	  raten[h]=r_summe;
	}
	//ACHTUNG: hier Korrektur ///////////////////////////////////////////////////
	else if ((partner[ punkt_ladung[vschritt] ][h]>(monomer_anzahl))&&(partner[ punkt_ladung[vschritt] ][h]!=punkt[vschritt])){
	  r_summe=r_summe;
	  raten[h]=0;
	}
	/////////////////////////////////////////////////////////////////////////////
      }
      // #######################################################################################################################################################################################################
      // hier raten für fullerene addieren #####################################################################################################################################################################
      r_summe_fulleren=0;
      zufall1=distribution0(engine);  
      std::vector <double> raten_fulleren (partneranzahl [punkt[vschritt]]);
      for (h=0;h<partneranzahl [punkt[vschritt]];h++){
	if (partner [punkt[vschritt]][h] > (monomer_anzahl)){
          zufall=distribution0(engine);	
	  coulombenergy=coulomb( x, y, z, punkt_ladung[vschritt], partner[ punkt[vschritt] ][h],3.4088)-coulomb( x, y, z, punkt[vschritt], punkt_ladung[vschritt], 3.4088);
	  r_summe_fulleren=r_summe_fulleren+rate(coupling_fulleren[ punkt[vschritt] ][ partner[ punkt[vschritt] ][h] ], ((zufall-zufall1)+coulombenergy), fullerenreorganisationsenergie);
//	  cout << "FULLERENE" << endl;
//	  cout << "Rate " << setw(12) << setprecision(3) << scientific << rate(coupling_fulleren[ punkt[vschritt] ][ partner[ punkt[vschritt] ][h] ], (zufall-zufall1)+coulombenergy, fullerenreorganisationsenergie) << endl;	
//	  cout << "Energie " << setw(12) << setprecision(6) << fixed << (zufall-zufall1)+coulombenergy << endl;
//	  cout << "Kopplung " << setw(12) << setprecision(6) << fixed << coupling_fulleren[ punkt[vschritt] ][ partner[ punkt[vschritt] ][h] ] << endl;
	  raten_fulleren[h]=r_summe_fulleren;
	}
	if ((partner[ punkt[vschritt] ][h]<(monomer_anzahl+1))&&(partner[ punkt[vschritt] ][h]==punkt_ladung[vschritt])){
          zufall=distribution0(engine);	 
	  // coulomb energie berechnen
	  coulombenergy=coulomb( x, y, z, punkt[vschritt], partner[ punkt[vschritt] ][h], 1);
//	  cout << "FULLEREN-HETERODIMER" << endl;
	  r_summe_fulleren=r_summe_fulleren+rate(coupling_rek[ punkt[vschritt] ][ partner[ punkt[vschritt] ][h] ], (zufall-zufall1)+rekombinationstriebkraft-coulombenergy, rek_reorganisation);
//	  cout << "Rate " << setw(12) << setprecision(3) << scientific << rate(coupling_rek[ punkt[vschritt] ][ partner[ punkt[vschritt] ][h] ], (zufall-zufall1)+rekombinationstriebkraft-coulombenergy, rek_reorganisation) << endl;
	  raten_fulleren[h]=r_summe_fulleren;
	}
	//ACHTUNG: hier Korrektur ////////////////////////////////////////////////////////
	else if ((partner[ punkt[vschritt] ][h]<(monomer_anzahl+1))&&(partner[ punkt[vschritt] ][h]!=punkt_ladung[vschritt])){
	  r_summe_fulleren=r_summe_fulleren;
	  raten_fulleren[h]=0;
	}
	//////////////////////////////////////////////////////////////////////////////////
      }
//      cout << "R_Summe_fulleren ist " << setw(12) << setprecision(3) << scientific << r_summe_fulleren << endl;
//      cout << "Dwelltime Monomer " << setw(12) << setprecision(4) << fixed << 1/r_summe*1e12 << endl;
//      cout << "Dwelltime Fulleren " << setw(12) << setprecision(4) << fixed << 1/r_summe_fulleren*1e12 << endl;  
//      cout << "Zeit Monomer " << setw(12) << setprecision(4) << fixed << zeit_1*1e12 << endl;
//      cout << "Zeit Fulleren " << setw(12) << setprecision(4) << fixed << zeit_2*1e12 << endl;
      // #######################################################################################################################################################################################################     
      // hüpfendes teilchen bestimmen
      if ((1/r_summe-zeit_1)<(1/r_summe_fulleren-zeit_2)){
		 // std::cout << "Monomer hüpft zuerst." << endl;
		  run << "Monomer hüpft zuerst." << endl;
	//Update der Zeiten
	if ((1/r_summe-zeit_1)>0){
	  zeit=zeit+(1/r_summe-zeit_1);
	  zeit_2=zeit_2+(1/r_summe-zeit_1);
	  zeit_1=0;
	}
	else if ((1/r_summe-zeit_1)<0){
	  zeit=zeit;
	  zeit_2=zeit_2;
	  zeit_1=0;
	}
	else {
		//std::cout << "FEHLER!" << endl;
		run << "FEHLER!" << endl;
	  return 0;
	}
	// monomerhüpfen ausführen
        uniform_real_distribution<double> distribution1(0,1); 
        zufall=distribution1(engine);
        r_i=zufall*r_summe;	
	for (g=0;g< partneranzahl [punkt_ladung[vschritt]];g++){
	  if ((raten[g]>r_i)&&(partner [punkt_ladung[vschritt]][g]<(monomer_anzahl))){
		  std::cout << "Ladungstransport" << endl;
	    punkt_ladung[i]=partner [punkt_ladung[vschritt]][g];
	    punkt[i]=punkt[vschritt];
	    // Abbruchkriterium für Ladungstrennung
	    switch(ebene){
	      case 'x':
		if (((x[ punkt_ladung[i] ])-x_mittel)>(0.75*(x[ startpunkt[k] ]-x_mittel))){
		  ch_diss[k]++;
		  zeit_ch[k][j]=zeit-zeit_ex[k][j];
		  vel_ch[k][j]=((x[ punkt_ladung[i] ])-x_mittel)/zeit_ch[k][j];			  
		 // std::cout << "Ladungen getrennt" << endl;
		  run << "Ladungen getrennt" << endl;
		  zustand[k][j]='s';
		}
		break;
	      case 'y':
		if (((y[ punkt_ladung[i] ])-y_mittel)>(0.75*(y[ startpunkt[k] ]-y_mittel))){
		  ch_diss[k]++;
		  zeit_ch[k][j]=zeit-zeit_ex[k][j];
		  vel_ch[k][j]=((y[ punkt_ladung[i] ])-y_mittel)/zeit_ch[k][j];
		  //std::cout << "Ladungen getrennt" << endl;
		  run << "Ladungen getrennt" << endl;
		  zustand[k][j]='s';
		}
		break;
	      case 'z':
		if (((z[ punkt_ladung[i] ])-z_mittel)>(0.6*(z[ startpunkt[k] ]-z_mittel))){
		  ch_diss[k]++;
		  zeit_ch[k][j]=zeit-zeit_ex[k][j];
		  vel_ch[k][j]=((z[ punkt_ladung[i] ])-z_mittel)/zeit_ch[k][j];		  
		  //std::cout << "Ladungen getrennt" << endl;
		  run << "Ladungen getrennt" << endl; 
		  zustand[k][j]='s';
		}		
		break;
	    }
	    //#########################################################################################################################
		std::cout << "altes Monomer " << setw(5) << punkt_ladung[vschritt] << endl;
		std::cout << "neues Monomer " << setw(5) << punkt_ladung[i] << endl;
		std::cout << "Kopplung " << setw(12) << setprecision(6) << fixed << coupling_ladung[ punkt_ladung[vschritt] ][ punkt_ladung[i] ] << endl;
		std::cout << "Fulleren " << setw(5) << punkt[i] << setw(5) << punkt[vschritt] << endl;
	    break;
	  }
	  else if ((raten[g]>r_i)&&(partner[ punkt_ladung[vschritt] ][g]>(monomer_anzahl))){
		  //std::cout << "Rekombination" << endl;
		  run << "Rekombination" << endl;
	    zustand[k][j]='t';
	    rek[k]++;
	    break;
	  }
	  else if (g==(partneranzahl[ punkt_ladung[vschritt] ])) {
		  //std::cout << "ACHTUNG: FEHLER im Ladungstransport des p-Halbleiters." << endl;
		  run << "ACHTUNG: FEHLER im Ladungstransport des p-Halbleiters." << endl;
	    return 0;
	  }
	}
      } // endes des hüpfenden p-halbleiters
      else if ((1/r_summe-zeit_1)>(1/r_summe_fulleren-zeit_2)){
		  //std::cout << "Fulleren hüpft zuerst." << endl;
		  run << "Fulleren hüpft zuerst." << endl;
	if ((1/r_summe_fulleren-zeit_2)>0){
	  zeit=zeit+(1/r_summe_fulleren-zeit_2);
	  zeit_1=zeit_1+(1/r_summe_fulleren-zeit_2);
	  zeit_2=0;
	}
	else if ((1/r_summe_fulleren-zeit_2)<0){
	  zeit=zeit;
	  zeit_1=zeit_1;
	  zeit_2=0;
	}
	else {
		//std::cout << "FEHLER!" << endl;
		run << "FEHLER!" << endl;
	  return 0;
	}
	// fullerenhüpfen ausführen
        uniform_real_distribution<double> distribution1(0,1); 
        zufall=distribution1(engine);
        r_i=zufall*r_summe_fulleren;
//	cout << "R_i ist " << setw(12) << setprecision(6) << scientific  << r_i << endl;
//	cout << "Partneranzahl " << partneranzahl[ punkt[vschritt] ] << endl;
	for (g=0;g<(partneranzahl [punkt[vschritt]]);g++){
//	  cout << "g anfangs " << g << endl;
//	  cout << "Partner " << partner[ punkt[vschritt] ][g] << endl;
	  if ((raten_fulleren[g]>r_i)&&((partner [punkt[vschritt]][g])>monomer_anzahl)){
		  std::cout << "Ladungstranport im Fullerenphase" << endl;
//	    cout << "Raten_fulleren ist " << setw(12) << setprecision(6) << scientific << raten_fulleren[g] << endl; 	    
	    punkt[i]=partner [punkt[vschritt]][g];
	    punkt_ladung[i]=punkt_ladung[vschritt]; 
		std::cout << "altes Fulleren " << setw(5) << punkt[vschritt] << endl;
		std::cout << "neues Fulleren " << setw(5) << punkt[i] << endl;
		std::cout << "Monomer " << setw(5) << punkt_ladung[i] << endl;
		std::cout << "Kopplung " << setw(12) << setprecision(6) << coupling_fulleren[ punkt[i] ][ punkt[vschritt] ] << endl;
	    break;
	  }
	  else if ((raten_fulleren[g]>r_i)&&((partner[ punkt[vschritt] ][g])<(monomer_anzahl+1))){
//	    cout << "g ist " << g << endl;
//	    cout << "Raten_fulleren ist " << setw(12) << setprecision(6) << scientific << raten_fulleren[g] << endl; 
		  //std::cout << "Rekombination." << endl;
		  run << "Rekombination." << endl;
	    zustand[k][j]='t';
	    rek[k]++;
	    break;
	  }
	  else if (g==(partneranzahl[ punkt[vschritt] ])) {
		  //std::cout << "ACHTUNG: FEHLER im Ladungstransport des Fullerens." << endl;
		  run << "ACHTUNG: FEHLER im Ladungstransport des Fullerens." << endl;
	    return 0;
	  }
	}
      } // endes des hüpfenden fullerens
    } // ende des 'c'-zustands
	//__________________________________________________________________________________________________________

    else if (zustand[k][j]=='e'){
      // site energies berechnen
      r_summe=0;
      random_device rd;
      default_random_engine engine(rd() );
      normal_distribution<double> distribution0(0.0,0.0338987); //##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! hier neue standardabweichung eintragen
      zufall1=distribution0(engine); //generating an normal-distributed random number

      std::vector <double> raten (partneranzahl [punkt[vschritt]]+1);

//      cout << "neuer Zustand " << endl;
//      cout << "k_rad " << setw(12) << setprecision(5) << scientific << k_rad << endl;

      for (h=0;h<(partneranzahl [punkt[vschritt]]);h++)
	  {
		if (partner [punkt[vschritt]][h]<(monomer_anzahl)){ //interaction with a monomer
           zufall=distribution0(engine);// generatinjg a second normal distributed random number

//			 cout << "Excitonrate " << setw(12) << setprecision(5) << scientific << rate(coupling_exciton[ punkt[vschritt] ][ partner[ punkt[vschritt] ][h] ], (zufall-zufall1), reorganisationsenergie_exciton) << endl; 
		   r_summe+= rate(coupling_exciton[ punkt[vschritt]] [partner[punkt[vschritt]] [h] ], (zufall-zufall1), reorganisationsenergie_exciton);
		   raten[h]= r_summe;

raten_out << h << "	" << setw(5) << raten[h] << "	" << zufall << endl;
    	}

		else if (partner [punkt[vschritt]] [h]>(monomer_anzahl))
		{ //interaction with a fulleren
          zufall=distribution0(engine);	 
		// coulomb energie berechnen
		  coulombenergy=coulomb( x, y, z, punkt[vschritt], partner[ punkt[vschritt] ][h], 1);
		   r_summe+= rate(coupling_ct[ punkt[vschritt]] [partner [punkt[vschritt]] [h]], (zufall-zufall1)+chargetransfertriebkraft+coulombenergy, ct_reorganisation);
		   raten[h]= r_summe;

raten_out << h << "	" << setw(5) << raten[h] << "	" << zufall << endl;
		}
      } // end of h

raten_out << "r_summe: " << r_summe << endl << "-----------------------------------------------------" << endl;

      // fluoreszenz dazuaddieren
      r_summe=r_summe+k_rad;
      dwelltime=1/r_summe;
      // schritt bestimmen
      uniform_real_distribution<double> distribution1(0,1); 
      zufall=distribution1(engine);

raten_out << zufall << endl << "-------------------------------------------------------------------" << endl;

      r_i=zufall*r_summe;
      zeit=zeit+1/r_summe;
//      cout << "R_Summe ist " << setw(12) << setprecision(5) << scientific << r_summe << endl;
//      cout << "R_i ist " << setw(12) << setprecision(5) << scientific << r_i << endl;
//      cout << "Zeitschritt " << setw(12) << setprecision(4) << fixed << 1/r_summe*1e12 << endl;
      //falls trapping
      zufall=distribution1(engine);      
      if (zufall*(900e-1+1/r_summe)>(900e-1)){
		 // std::cout << "Exziton getrappt!" << endl;
		  run << "Exziton getrappt!" << endl;
	trapping[k]++;
	zustand[k][j]='t';
	break;
      }
//      cout << "Zeit ist " << setw(12) << setprecision(3) << zeit*1e10 << endl;

  for (g=0;g<((partneranzahl [punkt[vschritt]]));g++){// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	if (raten[g]>r_i){
	  punkt[i]= partner [punkt[vschritt]] [g]; //partner [punkt[vschritt]] doesn't reach needed values &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	  run << "	neuer Punkt " << punkt[i] << endl;
//	  if (punkt[i]<(monomer_anzahl+1)){
//	    cout << "Exzitonentransport." << endl;
//	    cout << "alter Punkt " << punkt[vschritt] <<endl;
//	    cout << "neuer Punkt " << punkt[i] << endl;
//	    cout << "Kopplung " << coupling_exciton[ punkt[vschritt] ][ punkt[i] ] << endl;	    
//	    cout << "Raten_g " << setw(12) << setprecision(5) << scientific << raten[g] << endl;
//	    cout << "Raten_g-1 " << setw(12) << setprecision(5) << scientific << raten[g-1] << endl;
//	    cout << "r_i " << setw(12) << setprecision(5) << scientific << r_i << endl;
//	  }
	  /*else*/ if (punkt[i] > monomer_anzahl) {//PROBLEM!!!!!!!!!!!!!!!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	  //	    punkt[i]=partner [punkt[vschritt]] [g];	
		  punkt_ladung[i] = punkt[vschritt];
		  zustand[k][j] = 'c';
		  //std::cout << "Ladungstrennung." << endl;
		  run << "Ladungstrennung." << endl;
		  //	    cout << "Kopplung " << coupling_ct[ punkt[vschritt] ][ punkt[i] ] << endl;
		  //	    cout << "Coulomb-Energie " << setw(12) << setprecision(6) << fixed <<  coulomb( x, y, z, punkt[vschritt], punkt[i], 1) << endl;
		  //	    cout << "Länge " << setw(12) << setprecision(6) << fixed << length(x,y,z,punkt[vschritt], punkt[i]) << endl;
		  vel_ex[k][j] = length(x, y, z, punkt[0], punkt[i]) / zeit;
		  //std::cout << "Exzitonengeschw. " << vel_ex[k][j]*1e-9 << endl;
		  run << "Exzitonengeschw. " << vel_ex[k][j] * 1e-9 << endl;
		  ex_diss[k]++;


		  break;// put break statement into if-statement so g can get >0 §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
	  }

	}
	else if (raten [partneranzahl [punkt[vschritt]]] < r_i){
		//std::cout << "strahlender Zerfall." << endl;
		run << "strahlender Zerfall." << endl;
	  radiativ[k]++;
	  zustand[k][j]='t';
	  break;
	}
      }
    } // end of 'e'-zustand
	//____________________________________________________________________________________________

    else if (zustand[k][j]=='t'){
		//std::cout << "KAPUTT!" << endl;
		run << "KAPUTT!" << endl;
      break;
    } // end of 't'-zustand
	//__________________________________________________________________________________________________

    else if (zustand[k][j]=='s'){
		//std::cout << "ERFOLG!" << endl;
		run << "ERFOLG!" << endl;
      break;
    } //end of 's'-zustand 
	//_________________________________________________________________________________________________

    else {
		//std::cout << "ACHTUNG. Zustand undefiniert!" << endl;
		run << "ACHTUNG. Zustand undefiniert!" << endl;
      return 0;
    } //end of undefined zustand
	//___________________________________________________________________________________________________

	vschritt = i; //saving number of previous schritt for use in next iteration

  } // ende über schleife i

} // Ende über Schleife über durchläufe für den gleichen startpunkt j



} // Ende über Schleife der Startpunkte k

run.close();
raten_out.close();
// Auswertung: Prozentsätze
ofstream auswertung;
auswertung.open("auswertung.txt");
auswertung << setw(4) << "k" << setw(5) << "IX" << setw(12) << "Ex_Diss." << setw(12) << "Ch_Diss" << setw(12) << "Rek." << setw(12) << "Trapp." << setw(12) << "Fluor." << '\n';
double ex_diss_efficiency, ch_diss_efficiency, rek_efficiency, trapp_efficiency, rad_efficiency;
double mittel_ex=0;
double mittel_ch=0;
double mittel_rek=0;
double mittel_trapp=0;
double mittel_rad=0;
std::vector <double> mittel_ex_vel (index), mittel_ch_vel (index), standard_ex (index), standard_ch(index);
int zahl;
for (k=0;k<index;k++){
  mittel_ex=mittel_ex+(ex_diss[k]*1.0);
  mittel_ch=mittel_ch+(ch_diss[k]*1.0);
  mittel_rek=mittel_rek+(rek[k]*1.0);
  mittel_trapp=mittel_trapp+(trapping[k]*1.0);
  mittel_rad=mittel_rad+(radiativ[k]*1.0);
  ex_diss_efficiency=ex_diss[k];
  ch_diss_efficiency=ch_diss[k];  
  rek_efficiency=rek[k];
  trapp_efficiency=trapping[k];
  rad_efficiency=radiativ[k];
  auswertung << setw(4) << k << setw(5) << startpunkt[k] << setw(12) << setprecision(5) << ex_diss_efficiency;
  auswertung << setw(12) << setprecision(5) << ch_diss_efficiency;  
  auswertung << setw(12) << setprecision(5) << rek_efficiency << setw(12) << setprecision(5) << trapp_efficiency << setw(12) << setprecision(5) << rad_efficiency << '\n';  
}
auswertung << setw(9) << "Mittel " << setw(12) << setprecision(5) << fixed << mittel_ex/index << setw(12) << setprecision(5) << fixed << mittel_ch/index;
auswertung << setw(12) << setprecision(5) << fixed << mittel_rek/index << setw(12) << setprecision(5) << fixed << mittel_trapp/index;
auswertung << setw(12) << setprecision(5) << fixed << mittel_rad/index << '\n';
auswertung << "GESCHWINDIGKEITEN" << '\n';
auswertung << setw(4) << "k" << setw(5) << "IX" << setw(12) << "Ex_vel" << setw(12) << "Ex_s_dev" << setw(12) << "Ch_vel" << setw(12) << "Ch_s_dev" << '\n';
// mittlere Geschwindigkeit
for (k=0;k<index;k++){
  mittel_ch_vel[k]=0;
  mittel_ex_vel[k]=0;
  // Mittelwert berechnen
  for (j=0;j<100;j++){
    if (vel_ch[k][j]>0.0001){
      mittel_ch_vel[k]=mittel_ch_vel[k]+vel_ch[k][j];
    }
    if (vel_ex[k][j]>0.0001){
      mittel_ex_vel[k]=mittel_ex_vel[k]+vel_ex[k][j];
    } 
  }
  if (ch_diss[k]>0){
    mittel_ch_vel[k]=mittel_ch_vel[k]/(ch_diss[k]*1.0);
  }
  else if (ch_diss[k]==0){
    mittel_ch_vel[k]=0;
  }
  if (ex_diss[k]>0){
    mittel_ex_vel[k]=mittel_ex_vel[k]/(ex_diss[k]*1.0);
  }
  else if (ex_diss[k]==0){
    mittel_ex_vel[k]=0;
  }
 // Standardabweichunb berechnen
 standard_ex[k]=0;
 standard_ch[k]=0;
  for (j=0;j<100;j++){
    if (vel_ch[k][j]>0.0001){
      standard_ch[k]=standard_ch[k]+(vel_ch[k][j]-mittel_ch_vel[k])*(vel_ch[k][j]-mittel_ch_vel[k]);
    }
    if (vel_ex[k][j]>0.0001){
      standard_ex[k]=standard_ex[k]+(vel_ex[k][j]-mittel_ex_vel[k])*(vel_ex[k][j]-mittel_ex_vel[k]);
    } 
  }
  if (ch_diss[k]>1){
    standard_ch[k]=sqrt(standard_ch[k]/(ch_diss[k]-1));
  }
  else if (ch_diss[k]<2){
    standard_ch[k]=0;    
  }
  if (ex_diss[k]>1){
    standard_ex[k]=sqrt(standard_ex[k]/(ex_diss[k]-1));
  }
  else if (ex_diss[k]<2){
    standard_ex[k]=0;    
  }
  auswertung << setw(4) << k << setw(5) << startpunkt[k] << setw(12) << setprecision(5) << fixed << mittel_ex_vel[k]*1e-9;
  auswertung << setw(12) << setprecision(5) << fixed << standard_ex[k]*1e-9;
  auswertung << setw(12) << setprecision(5) << fixed << mittel_ch_vel[k]*1e-9;
  auswertung << setw(12) << setprecision(5) << fixed << standard_ch[k]*1e-9 << '\n';
}
double mittelwert_geschw_exciton=0;
double mittelwert_geschw_ladung=0;
for (k=0;k<index;k++){
  mittelwert_geschw_exciton=mittelwert_geschw_exciton+mittel_ex_vel[k];
  mittelwert_geschw_ladung=mittelwert_geschw_ladung+mittel_ch_vel[k];
}
auswertung << setw(9) << "Mittel " << setw(12) << setprecision(5) << fixed << mittelwert_geschw_exciton/index*1e-9;
auswertung << setw(12) << setprecision(5) << fixed << mittelwert_geschw_ladung/index*1e-9 << '\n';
// Verteilung Ladungen und Exzitonengeschwindigkeiten
ofstream exciton_verteilung;
exciton_verteilung.open("exciton_verteilung.txt");
for (i=0;i<20;i++){
  zahl=0;
  for (k=0;k<index;k++){
    for (j=0;j<100;j++){
      if ((vel_ex[k][j]>(i*50*1e9))){
	zahl++;
      }
    }
  }
  exciton_verteilung << setw(12) << setprecision(5) << i*50 << setw(12) << zahl/index << '\n';
}
exciton_verteilung.close();
exciton_verteilung.open("ladung_verteilung.txt");
for (i=0;i<20;i++){
  zahl=0;
  for (k=0;k<index;k++){
    for (j=0;j<100;j++){
      if ((vel_ch[k][j]>(i*50*1e9))){
	zahl++;
      }
    }
  }
  exciton_verteilung << setw(12) << setprecision(5) << i*50 << setw(12) << zahl/index << '\n';
}
exciton_verteilung.close();
}

