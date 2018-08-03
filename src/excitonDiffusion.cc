#include "excitonDiffusion.h"

void exciD::dimexc(std::string masscenters, std::string couplings) {

  double reorganisationsenergie_exciton = Config::get().exbreak.ReorgE_exc;//noch extra variablen in config.h und config.cc einfügen
  double oszillatorstrength = Config::get().exbreak.oscillatorstrength;
  double wellenzahl = Config::get().exbreak.wellenzahl;
  double krad = wellenzahl*wellenzahl * oszillatorstrength; // fluoreszenz

  

  std::ifstream comf;
  std::ifstream coupf;
  double numbermon;
  std::vector<coords::Cartesian_Point> com;
  
  comf.open(masscenters);
  comf >> numbermon;
  comf.ignore(256, '\n');//ignore all till newline

  coords::Cartesian_Point tmp;
  while (!comf.eof())
  {
    comf >> tmp;
    com.push_back(tmp);
  }

  if(numbermon != com.size())
  {throw std::logic_error("Unclear number of centers of mass.");
  }

  coupf.open(couplings);

  std::size_t tmpA, tmpB;
  double tmpC;
  std::vector <double> tmpD;
  std::vector<exciD::Couplings> excCoup;

 while (!coupf.eof())
  {
    coupf >> tmpA >> tmpB >> tmpC;
    tmpD.push_back(tmpC);
    exciD::Couplings tmpE(tmpA, tmpB, tmpD);
    excCoup.push_back(tmpE);
  }

 for (std::size_t i = 0u; i < excCoup.size(); i++)
 {
   std::cout << excCoup[i].monA << " " << excCoup[i].monB << " " << excCoup[i].coupling[0] << '\n';
 }
}