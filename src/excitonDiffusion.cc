#include "excitonDiffusion.h"

coords::Cartesian_Point exciD::avgDimCoM(coords::Cartesian_Point posA, coords::Cartesian_Point posB)//funtion to calculate average of two cartesian coordinates, used for average com of a dimer
{
  double x = (posA.x() + posB.x()) / 2;
  double y = (posA.y() + posB.y()) / 2;
  double z = (posA.z() + posB.z()) / 2;

  coords::Cartesian_Point ret;
  ret.x() = x;
  ret.y() = y;
  ret.z() = z;

  return ret;
}

double exciD::length(coords::Cartesian_Point pointA, coords::Cartesian_Point pointB)//function to calculate the length of a vector between two cartesian points
{
  double len = sqrt((pointA.x() - pointB.x())*(pointA.x() - pointB.x()) + (pointA.y() - pointB.y())*(pointA.y() - pointB.y()) + (pointA.z() - pointB.z())*(pointA.z() - pointB.z()));
  return len;
}

double exciD::marcus(double coupling, double drivingF, double reorganisation)
{
  double marc;
  marc = (coupling*coupling) / exciD::h_quer * sqrt(M_PI / (reorganisation*exciD::boltzmann_const * 298))*exp(-(reorganisation + drivingF)*(reorganisation + drivingF) / (4 * exciD::boltzmann_const * 298 * reorganisation));//replace 298 by user defined temperatuer?
  return marc;
}

double exciD::coulomb(coords::Cartesian_Point aktPos, coords::Cartesian_Point targetPos, double e_relative)
{
  double e_0 = 8.854187e-12;
  double elementar = 1.60217662e-19;
  double coulomb = -elementar / (4 * M_PI*e_0*e_relative*length(aktPos, targetPos)*1e-10);
  return coulomb;
}

coords::Cartesian_Point exciD::structCenter(std::vector<exciD::Couplings> excCoup)//calculate average position of the given positions (more than two)
{
  coords::Cartesian_Point avg(0.0, 0.0, 0.0);
  for (std::size_t i = 0u; i < excCoup.size(); i++)
  {
    avg += excCoup[i].position;
  }
  avg /= excCoup.size();
  return avg;
}

coords::Cartesian_Point exciD::min(std::vector<exciD::Couplings> coords)
{
  coords::Cartesian_Point min;

  min = coords[0].position;//initialize with coordinates of first dimer

  for (std::size_t i = 1u; i < coords.size(); i++)
  {
    if (coords[i].position.x() < min.x())
    {
      min.x() = coords[i].position.x();
    }

    if (coords[i].position.y() < min.y())
    {
      min.y() = coords[i].position.y();
    }

    if (coords[i].position.z() < min.z())
    {
      min.z() = coords[i].position.z();
    }
  }
  return min;
}

void exciD::dimexc(std::string masscenters, std::string couplings, int pscnumber, int nscnumber) {
  try {

    double reorganisationsenergie_exciton = Config::get().exbreak.ReorgE_exc;//noch extra variablen in config.h und config.cc einfügen
    double reorganisationsenergie_charge = Config::get().exbreak.ReorgE_ch;
    double reorganisationsenergie_nSC = Config::get().exbreak.ReorgE_nSC;
    double reorganisationsenergie_ct = Config::get().exbreak.ReorgE_ct;
    double reorganisationsenergie_rek = Config::get().exbreak.ReorgE_rek;
    double triebkraft_ct = Config::get().exbreak.ct_triebkraft;
    double triebkraft_rek = Config::get().exbreak.rek_triebkraft;
    double oszillatorstrength = Config::get().exbreak.oscillatorstrength;
    double wellenzahl = Config::get().exbreak.wellenzahl;
    double k_rad = wellenzahl * wellenzahl * oszillatorstrength; // fluoreszenz

    char plane = 'z';//don't forget to replace by userinput

    std::ifstream comf;
    std::ifstream coupf;
    std::size_t numbermon;
    std::vector<std::size_t> startPind, viablePartners, h_viablePartners;
    std::vector<exciD::Partners> partnerConnections, h_partnerConnections;
    coords::Representation_3D com;
    coords::Cartesian_Point avg;
    exciD::Exciton excPos;
 

    comf.open(masscenters);
    comf >> numbermon;
    //comf.ignore(256, '\n');//ignore all till newline 

    coords::Cartesian_Point tmp;
    int unnecessaryindex;
    std::string line;

    //read centers of mass
    while (!comf.eof())
    {
      std::getline(comf, line);
      std::istringstream iss(line);
    
      if (line.size() != 0)
      {
        iss >> unnecessaryindex >> tmp;
        com.push_back(tmp);
      }
    }

    //check if the number of masscenters matches the expected value
      std::cout << numbermon << " " << com.size() << '\n';    
    if (numbermon != com.size())
    {
      throw std::logic_error("Unclear number of centers of mass.");
    }

    coupf.open(couplings);

    std::string tmpcount;
    std::size_t tmpA, tmpB;
    double tmpC, tmpD;
    exciD::Couplings tmpE;
    std::vector<exciD::Couplings> excCoup;
    double avgCoup;

    //read monomer indices of dimers and corresponding exciton-coupling and save in vector of couplings
    while (!coupf.eof())
    {
      std::getline(coupf, line);
      std::istringstream iss(line);
      std::size_t numberofelements(0);
      while (!iss.eof())
      {
        iss >> tmpcount;
        numberofelements++;
      }
      std::istringstream iss1(line);
      if (numberofelements < 4)
      {
        iss1 >> tmpA >> tmpB >> tmpC;
        tmpD = 0.0;
      }
      else if (numberofelements < 6)
      {
        iss1 >> tmpA >> tmpB >> tmpC >> tmpD;
      }
      else
      {
        throw std::logic_error("Unexpected number of elements in couplingsfile.");
      }
      exciD::Couplings tmpE(tmpA, tmpB, tmpC, tmpD);
      excCoup.push_back(tmpE);
    }
    std::cout << "Number of dimer-pairs " << excCoup.size() << '\n';

    //calculate average centers of mass for relevant dimers and add to struct
    for (std::size_t i = 0u; i < excCoup.size(); i++)
    {
      excCoup[i].position = (exciD::avgDimCoM(com[excCoup[i].monA], com[excCoup[i].monB]));
    }

    avg = exciD::structCenter(excCoup);

    coords::Cartesian_Point minV = min(excCoup);

    std::random_device rd; //prepare rng
    std::default_random_engine engine(rd());
    std::normal_distribution<double> distributionN(0.0, 0.068584577);
    std::uniform_real_distribution<double> distributionR(0, 1); //beispiel für rng var: double rng = distributionN(engine);

    //choose startingpoints
    switch (plane)
    {
    case'x':

      for (std::size_t i = 0u; i < excCoup.size(); i++)
      {
        if (excCoup[i].position.x() - avg.x() < 0.85 * (minV.x() - avg.x()))
        {
          startPind.push_back(i);
        }
      }

      break;

    case'y':

      for (std::size_t i = 0u; i < excCoup.size(); i++)
      {
        if (excCoup[i].position.y() - avg.y() < 0.85 * (minV.y() - avg.y()))
        {
          startPind.push_back(i);
        }
      }

      break;



    case'z':
      for (std::size_t i = 0u; i < excCoup.size(); i++)
      {
        if (excCoup[i].position.z() - avg.z() < 0.85 * (minV.z() - avg.z()))
        {
          startPind.push_back(i);
        }
      }
      break;

    }

    //loop for writing starting points
       /* for (std::size_t i = 0u; i < startPind.size(); i++)
        {
          std::cout << "Startingpoint " << i << ": " << startPind[i] << " Monomer A " << excCoup[startPind[i]].monA << " Monomer B " << excCoup[startPind[i]].monB << '\n';
        }*/

    std::vector <int> trapped(startPind.size(), 0);//for counting the trapped excitons unable to reach the interface from each startingpoint
    std::vector <int> ex_diss(startPind.size(), 0);
    std::vector <int> radiating(startPind.size(), 0);
    std::vector <int> rekombined(startPind.size(), 0);
    std::vector <int> ch_separation(startPind.size(), 0);

    std::vector <std::vector<double>> time_ch(startPind.size(), std::vector<double> (100,0.)), //vectors to keep time/velocities for different startingpoints and tries
                                      time_ex(startPind.size(), std::vector<double> (100,0.)),
                                      vel_ch(startPind.size(), std::vector<double> (100,0.)), 
                                      vel_ex(startPind.size(), std::vector<double> (100,0.));

    double time(0.0), time_p(0.0), time_n(0.0);

    //loop over all startingpoints 
    for (std::size_t i = 0u; i < startPind.size(); i++)
    {


      //loop ensures to start 100 times from every startingpoint
      for (std::size_t j = 0u; j < 5; j++)//don't forget to set to 100 when all works fine
      {

        excPos.location = startPind[i];//startPind[i] is the startingpoint for the actual simulation, excPos is used to keep track of the position of the exciton during simulation

        std::cout << "Startingpoint " << i << ": " << startPind[i] << " Monomer A: " << excCoup[startPind[i]].monA << " Monomer B: " << excCoup[startPind[i]].monB << " Coupling: " << excCoup[excPos.location].coupling << '\n';

        for (std::size_t h = 0u; h < 5; h++)//steps in each try for testing hardcode
        {
          std::cout << "Exciton Position: " << excPos.location << '\n';

          //loop over all dimerpairs for viable partners
          for (std::size_t k = 0; k < excCoup.size(); k++)
          {

            //in this if-cause the viable partners to the actual location are determined and the dimer pairs necessary to calculate the average couplings between the dimers are determined.
            if (excPos.location != k)//skip k if k is the index of the dimer where the exciton is at the moment
            {
              //ensure monomers at excCoup[excPos] are not part of excCoup[k]
              if ((excCoup[excPos.location].monA != excCoup[k].monA) && (excCoup[excPos.location].monA != excCoup[k].monB) && (excCoup[excPos.location].monB != excCoup[k].monA) && (excCoup[excPos.location].monB != excCoup[k].monB))
              {
                //check if excCoup[k] is close enough to excCoup[excPos]
                if (exciD::length(excCoup[excPos.location].position, excCoup[k].position) < 5.0)
                {
                  viablePartners.push_back(k);
                }
              }
            }// #if(excPos != k)

            if (excPos.state == 'c')
            {
              //same purpose as above but for location of possible second particle
              if (excPos.h_location != k)
              {
                if ((excCoup[excPos.h_location].monA != excCoup[k].monA) && (excCoup[excPos.h_location].monA != excCoup[k].monB) && (excCoup[excPos.h_location].monB != excCoup[k].monA) && (excCoup[excPos.h_location].monB != excCoup[k].monB))
                {
                  //check if excCoup[k] is close enough to excCoup[excPos]
                  if (exciD::length(excCoup[excPos.h_location].position, excCoup[k].position) < 5.0)
                  {
                    h_viablePartners.push_back(k);
                  }
                }
              }
            }

          }//dertermining of viable partners k

                //check if couplings of monomers in excPos exist to viable partner monomers, if not set value to zero
          for (std::size_t m = 0u; m < viablePartners.size(); m++)//loop over viable partners to find couplings between monomers in current posirtion and viable partners
          {
            std::vector<std::size_t> tmpG;
            std::cout << "  Partners: " << excCoup[viablePartners[m]].monA << " " << excCoup[viablePartners[m]].monB << " " << excCoup[viablePartners[m]].coupling << '\n';
            for (std::size_t l = 0u; l < excCoup.size(); l++)//loop over all dimerpairs which have couplings
            {
              if (excCoup[excPos.location].monA == excCoup[l].monA || excCoup[excPos.location].monA == excCoup[l].monB)//look if monomer A of the current location is part of the viewed Dimer
              {
                if (excCoup[excPos.location].monB != excCoup[l].monA && excCoup[excPos.location].monB != excCoup[l].monB)//ensure the other monomer is not also part of the viewed dimer
                {
                  if (excCoup[viablePartners[m]].monA == excCoup[l].monA || excCoup[viablePartners[m]].monA == excCoup[l].monB || excCoup[viablePartners[m]].monB == excCoup[l].monA || excCoup[viablePartners[m]].monB == excCoup[l].monB)//look if monomer of viablöe Parner is part of the viewed dimer
                  {
                    tmpG.push_back(l);
                  }
                }
              }
              else if (excCoup[excPos.location].monB == excCoup[l].monA || excCoup[excPos.location].monB == excCoup[l].monB)//look if monomer B of the current location is part of the viewed Dimer if monomer A is not
              {
                if (excCoup[excPos.location].monA != excCoup[l].monA && excCoup[excPos.location].monA != excCoup[l].monB)//ensure the other monomer is not also part of the viewed dimer
                {
                  //for (std::size_t m = 0u; m < viablePartners.size(); m++)//loop over viable partners to find couplings between monomers in current posirtion and viable partners
                  //{
                  if (excCoup[viablePartners[m]].monA == excCoup[l].monA || excCoup[viablePartners[m]].monA == excCoup[l].monB || excCoup[viablePartners[m]].monB == excCoup[l].monA || excCoup[viablePartners[m]].monB == excCoup[l].monB)//look if monomer of viablöe Parner is part of the viewed dimer
                  {
                    tmpG.push_back(l);
                  }
                  //}
                }
              }
            }// l
            exciD::Partners tmpH(viablePartners[m], tmpG);
            partnerConnections.push_back(tmpH);
          } //m

          if (excPos.state == 'c')//hole movement only of interesst if simulation of charges is done (steate=c)
          {
            for (std::size_t m = 0u; m < h_viablePartners.size(); m++)//loop over viable partners to find couplings between monomers in current posirtion and viable partners
            {
              std::vector<std::size_t> tmpG;
              std::cout << "  Partners: " << excCoup[h_viablePartners[m]].monA << " " << excCoup[h_viablePartners[m]].monB << " " << excCoup[h_viablePartners[m]].coupling << '\n';
              for (std::size_t l = 0u; l < excCoup.size(); l++)//loop over all dimerpairs which have couplings
              {
                if (excCoup[excPos.h_location].monA == excCoup[l].monA || excCoup[excPos.h_location].monA == excCoup[l].monB)//look if monomer A of the current location is part of the viewed Dimer
                {
                  if (excCoup[excPos.h_location].monB != excCoup[l].monA && excCoup[excPos.h_location].monB != excCoup[l].monB)//ensure the other monomer is not also part of the viewed dimer
                  {
                    if (excCoup[h_viablePartners[m]].monA == excCoup[l].monA || excCoup[h_viablePartners[m]].monA == excCoup[l].monB || excCoup[h_viablePartners[m]].monB == excCoup[l].monA || excCoup[h_viablePartners[m]].monB == excCoup[l].monB)//look if monomer of viablöe Parner is part of the viewed dimer
                    {
                      tmpG.push_back(l);
                    }
                  }
                }
                else if (excCoup[excPos.h_location].monB == excCoup[l].monA || excCoup[excPos.h_location].monB == excCoup[l].monB)//look if monomer B of the current location is part of the viewed Dimer if monomer A is not
                {
                  if (excCoup[excPos.h_location].monA != excCoup[l].monA && excCoup[excPos.h_location].monA != excCoup[l].monB)//ensure the other monomer is not also part of the viewed dimer
                  {
                    if (excCoup[h_viablePartners[m]].monA == excCoup[l].monA || excCoup[h_viablePartners[m]].monA == excCoup[l].monB || excCoup[h_viablePartners[m]].monB == excCoup[l].monA || excCoup[h_viablePartners[m]].monB == excCoup[l].monB)//look if monomer of viablöe Parner is part of the viewed dimer
                    {
                      tmpG.push_back(l);
                    }
                  }
                }
              }
              exciD::Partners tmpH(h_viablePartners[m], tmpG);
              h_partnerConnections.push_back(tmpH);
            }
          }

          //calculate avgCouplings & avgsecCouplings
          for (std::size_t n = 0u; n < partnerConnections.size(); n++)
          {
            partnerConnections[n].avgCoup = 0.0;//prevent visiting undefined behaviour land
            partnerConnections[n].avgsecCoup = 0.0;

            for (std::size_t o = 0u; o < partnerConnections[n].connect.size(); o++)
            {
              partnerConnections[n].avgCoup += excCoup[partnerConnections[n].connect[o]].coupling;
              
              if (excCoup[partnerConnections[n].connect[o]].seccoupling != 0.0)//if a opair has a second coupling its value is addet to avgsecCoup
              {
                partnerConnections[n].avgsecCoup += excCoup[partnerConnections[n].connect[o]].seccoupling;
              }
            }// o
            //partnerConnections[w].avgCoup /= partnerConnections[w].connect.size();
            partnerConnections[n].avgCoup /= 4; //unsure if the average coupling is gained by dividing the sum of relevant couplings by their number or the maximum (4) number of relevant couplings (assuming not added couplings are zero).
            
            if (partnerConnections[n].avgsecCoup != 0.0)//not every pair has a second coupling so bevoreS dividing its existence is checked
            {
              partnerConnections[n].avgsecCoup /= 4;
            }
          }// n

          if (excPos.state == 'c')//hole movement only of interesst if simulation of charges is done (steate=c)
          {
            for (std::size_t n = 0u; n < h_partnerConnections.size(); n++)
            {
              partnerConnections[n].avgCoup = 0.0;//prevent visiting undefined behaviour land
              partnerConnections[n].avgsecCoup = 0.0;

              for (std::size_t o = 0u; o < h_partnerConnections[n].connect.size(); o++)
              {
                h_partnerConnections[n].avgCoup += excCoup[h_partnerConnections[n].connect[o]].coupling;

                if (excCoup[h_partnerConnections[n].connect[o]].seccoupling != 0.0)//if a opair has a second coupling its value is addet to avgsecCoup
                {
                  h_partnerConnections[n].avgsecCoup += excCoup[h_partnerConnections[n].connect[o]].seccoupling;
                }
              }
              h_partnerConnections[n].avgCoup /= 4; //unsure if the average coupling is gained by dividing the sum of relevant couplings by their number or the maximum (4) number of relevant couplings (assuming not added couplings are zero).

              if (h_partnerConnections[n].avgsecCoup != 0.0)//not every pair has a second coupling so bevoreS dividing its existence is checked
              {
                h_partnerConnections[n].avgsecCoup /= 4;
              }
            }
          }

          //writing loop for calculated avg Couplings
          for (std::size_t n = 0u; n < partnerConnections.size(); n++)
          {
            std::cout << "PartnerIndex: " << partnerConnections[n].partnerIndex << " Monomers: " << excCoup[partnerConnections[n].partnerIndex].monA << " " << excCoup[partnerConnections[n].partnerIndex].monB << '\n';
            for (std::size_t o = 0u; o < partnerConnections[n].connect.size(); o++)
            {
              std::cout << " ConnectorIndex: " << partnerConnections[n].connect[o] << " Monomers: " << excCoup[partnerConnections[n].connect[o]].monA << " " << excCoup[partnerConnections[n].connect[o]].monB
                << " Coupling: " << excCoup[partnerConnections[n].connect[o]].coupling << " |" << " secCoupling: " << excCoup[partnerConnections[n].connect[o]].seccoupling << " |" << '\n';
            } //o
            std::cout << " Average Coupling: " << partnerConnections[n].avgCoup << " |" << "Average secCoupling: " << partnerConnections[n].avgsecCoup <<  '\n';
          }// n

          if (excPos.state == 'c')//hole movement only of uinteresst if simulation of charges is done (steate=c
          {
            for (std::size_t n = 0u; n < h_partnerConnections.size(); n++)
            {
              std::cout << "h_PartnerIndex: " << h_partnerConnections[n].partnerIndex << " Monomers: " << excCoup[h_partnerConnections[n].partnerIndex].monA << " " << excCoup[h_partnerConnections[n].partnerIndex].monB << '\n';
              for (std::size_t o = 0u; o < h_partnerConnections[n].connect.size(); o++)
              {
                std::cout << " h_ConnectorIndex: " << h_partnerConnections[n].connect[o] << " Monomers: " << excCoup[h_partnerConnections[n].connect[o]].monA << " " << excCoup[h_partnerConnections[n].connect[o]].monB
                  << " h_Coupling: " << excCoup[partnerConnections[n].connect[o]].coupling << " |" << " h_Coupling: " << excCoup[h_partnerConnections[n].connect[o]].seccoupling << " |" << '\n';
              } //o
              std::cout << " Average Coupling: " << h_partnerConnections[n].avgCoup << " |" << "Average secCoupling: " << h_partnerConnections[n].avgsecCoup << '\n';
            }// n
          }

          //algorithm for exciton movement starts here all before was preparation to know where movement to is possible

          std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << '\n';

          double rate_sum(0.0), rateFul_sum(0.0), coulombenergy, rate_KMC, tmp_ratesum(0.0);
          std::vector <double> raten;//used for exciton and electron rates
          std::vector <double> raten_hole;//used for hole rates
          double random_normal, random_normal1;
          double random_eq;
          bool heterodimer(false);
          random_normal1 = distributionN(engine);
//EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
          if (excPos.state == 'e')//Exciton state
          {
            for (std::size_t p = 0u; p < partnerConnections.size(); p++)
            {
              if (excCoup[partnerConnections[p].partnerIndex].monA < pscnumber && excCoup[partnerConnections[p].partnerIndex].monB < pscnumber)//only for homo p-type SC due to reorganisation energies
              {
                random_normal = distributionN(engine);//generating normal distributed random number

                rate_sum += marcus(partnerConnections[p].avgsecCoup, (random_normal - random_normal1), reorganisationsenergie_exciton);//rate for homoPSCpartner
              }

              else if (excCoup[partnerConnections[p].partnerIndex].monA > pscnumber && excCoup[partnerConnections[p].partnerIndex].monB > pscnumber)//CT only if BOTH patrner molecules are n-type SC
              {
                random_normal = distributionN(engine);//generating normal distributed random number

                coulombenergy = coulomb(excCoup[excPos.location].position, excCoup[partnerConnections[p].partnerIndex].position, 1);
                rate_sum += marcus(partnerConnections[p].avgCoup, (random_normal - random_normal1) + coulombenergy + triebkraft_ct, reorganisationsenergie_ct);//rate for heteroPSCpartner
              }
              else//to prevent heterodimersfrom participating as exciton location HOW TO HANDLE HETERO DIMERS? CT OR EXCITONDIFFUSION?
              {
                tmp_ratesum = rate_sum;
                rate_sum = 0.0;
                heterodimer = true;
              }
              raten.push_back(rate_sum);
              if (heterodimer) //set rate_sum back to previous value
              {
                heterodimer = false;
                rate_sum = tmp_ratesum;
              }
              std::cout << "Partner: " << partnerConnections[p].partnerIndex << " Rates: " << rate_sum << '\n';
            }// p
            rate_sum += k_rad;//accounting for fluorescence

            double random_real = distributionR(engine);

            rate_KMC = random_real * rate_sum;

            //calculate time needed for step
            time += (1 / rate_sum);

            //trapping
            random_real = distributionR(engine);
            if (random_real * (900e-1 + 1 / rate_sum) > (900e-1))
            {
              std::cout << "trapped" << '\n';
              trapped[i]++;
              excPos.state = 't';
              viablePartners.clear();//empties vector containing possible partners for step so it can be reused in next step
              partnerConnections.clear();
              break;
            }

            for (std::size_t q = 0u; q < viablePartners.size(); q++)
            {
              if (raten[q] > rate_KMC)
              {
                if (excCoup[viablePartners[q]].monA >= pscnumber && excCoup[viablePartners[q]].monB >= pscnumber)
                {
                  excPos.h_location = excPos.location; //set hole position to former exciton position for charge separation
                }

                excPos.location = viablePartners[q];//after chargeseparation excPos.location tracks electron position.

                if (excCoup[excPos.location].monA >= pscnumber && excCoup[excPos.location].monB >= pscnumber)
                {
                  excPos.state = 'c';
                  std::cout << "Chargeseparation." << '\n';
                  ex_diss[i]++;
                  excPos.state = 's';//till chargemovement is implemented
                }
                viablePartners.clear();//empties vector containing possible partners for step so it can be reused in next step
                partnerConnections.clear();
                break;
              }

              else if (raten.back() < rate_KMC)
              {
                std::cout << "Radiating decay." << '\n';
                radiating[i]++;
                excPos.state = 't';
                viablePartners.clear();//empties vector containing possible partners for step so it can be reused in next step
                partnerConnections.clear();
                break;
              }
            }
          }//state e end
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          else if (excPos.state == 'c')//charge separated state
          {
            if (nscnumber == 0)//if only pSC are availeable
            {
              std::cout << "Only Exciton movement in p-type semiconductor possible, due to lack of n-type semiconductors.";
              excPos.state = 's';
              break;
            }

            //hole rates in pSC
            for (std::size_t p=0u; p < partnerConnections.size(); p++)
            {
              if (excCoup[partnerConnections[p].partnerIndex].monA < pscnumber && excCoup[partnerConnections[p].partnerIndex].monB < pscnumber)//movement on pSC
              {
                random_normal = distributionN(engine);//generating normal distributed random number
                coulombenergy = coulomb(excCoup[excPos.h_location].position, excCoup[partnerConnections[p].partnerIndex].position, 3.4088) - coulomb(excCoup[excPos.h_location].position, excCoup[excPos.location].position, 3.4088);
                rate_sum += marcus(partnerConnections[p].avgCoup, (random_normal - random_normal1) + coulombenergy , reorganisationsenergie_charge);
              }
              else if (excCoup[partnerConnections[p].partnerIndex].monA > pscnumber && excCoup[partnerConnections[p].partnerIndex].monB > pscnumber && partnerConnections[p].partnerIndex == excPos.location)//movement to nSC --> recombination | only possible if electron present on nSC dimer
              {
                random_normal = distributionN(engine);//generating normal distributed random number
                coulombenergy = coulomb(excCoup[excPos.location].position, excCoup[partnerConnections[p].partnerIndex].position, 1);
                rate_sum += marcus(partnerConnections[p].avgCoup, (random_normal - random_normal1) + coulombenergy, reorganisationsenergie_rek);
              }   
              else if (excCoup[partnerConnections[p].partnerIndex].monA > pscnumber && excCoup[partnerConnections[p].partnerIndex].monB > pscnumber && partnerConnections[p].partnerIndex != excPos.location)
              { //recombination only possible if electron is pressent on nSC dimer => no hopping to nSC if no electron present
                raten[p] = 0.0;
              }
              else//to prevent heterodimersfrom participating as electron location HOW TO HANDLE HETERO DIMERS? holediffusion or recombination?
              {
                tmp_ratesum = rate_sum;
                rate_sum = 0.0;  
                heterodimer = true;
              }

              raten.push_back(rate_sum);

              if (heterodimer) //set rate_sum back to previous value
              {
                heterodimer = false;
                rate_sum = tmp_ratesum;
              }
              std::cout << "Partner: " << partnerConnections[p].partnerIndex << " h_Rates: " << rate_sum << '\n';
            }//p

            //electron rates in nSC
            for (std::size_t p = 0u; p < partnerConnections.size(); p++)
            {
              if (excCoup[partnerConnections[p].partnerIndex].monA > pscnumber && excCoup[partnerConnections[p].partnerIndex].monB > pscnumber)//movement on nSC
              {
                random_normal = distributionN(engine);//generating normal distributed random number
                coulombenergy = coulomb(excCoup[excPos.location].position, excCoup[partnerConnections[p].partnerIndex].position, 3.4088) - coulomb(excCoup[excPos.location].position, excCoup[excPos.h_location].position, 3.4088);
                rateFul_sum += marcus(partnerConnections[p].avgCoup, (random_normal - random_normal1) + coulombenergy , reorganisationsenergie_nSC);
              }
              else if (excCoup[partnerConnections[p].partnerIndex].monA < pscnumber && excCoup[partnerConnections[p].partnerIndex].monB < pscnumber && partnerConnections[p].partnerIndex == excPos.h_location)//movement to pSC --> recombination | only possible if electron present on nSC dimer
              {
                random_normal = distributionN(engine);//generating normal distributed random number
                coulombenergy = coulomb(excCoup[excPos.h_location].position, excCoup[partnerConnections[p].partnerIndex].position, 1);
                rateFul_sum += marcus(partnerConnections[p].avgCoup, (random_normal - random_normal1) + coulombenergy, reorganisationsenergie_rek);
              }
              else if (excCoup[partnerConnections[p].partnerIndex].monA < pscnumber && excCoup[partnerConnections[p].partnerIndex].monB < pscnumber && partnerConnections[p].partnerIndex != excPos.h_location)
              { //recombination only possible if hole is pressent on pSC dimer => no hopping to pSC if no hole present
                raten[p] = 0.0;
              }
              else//to prevent heterodimersfrom participating as electron location HOW TO HANDLE HETERO DIMERS? holediffusion or recombination?
              {
                tmp_ratesum = rateFul_sum;
                rateFul_sum = 0.0;
                heterodimer = true;
              }

              raten.push_back(rateFul_sum);

              if (heterodimer) //set rate_sum back to previous value
              {
                heterodimer = false;
                rateFul_sum = tmp_ratesum;
              }

              std::cout << "Partner: " << partnerConnections[p].partnerIndex << " e_Rates: " << rateFul_sum << '\n';
            }
            
            //decide hopping particle
            
            if((1 / rate_sum - time_p) < (1/rateFul_sum - time_n))
            {
              std::cout << "pSC hopps first. " << std::endl;

              //Update time
              if((1/rate_sum - time_p) > 0)
              {
                time   += (1/rate_sum - time_p);
                time_n += (1/rate_sum - time_p);
                time_p = 0.;
              }
              else if((1/rate_sum - time_p) < 0)
              {
                //time = time  //not necessary but intendet to help readeability
                //time_n = time_n
                time_p = 0.;
              }
              else
              {
                throw std::logic_error("Something went wrong with the time. Call the Doctor.");
              }

              //do hopping
               auto random_real = distributionR(engine);
               rate_KMC = random_real * rate_sum;

               for(std::size_t g = 0; g < viablePartners.size(); g++)
               {
                 if(raten[g] >rate_KMC)
                 {
                  if(excCoup[partnerConnections[g].partnerIndex].monA < pscnumber && excCoup[partnerConnections[g].partnerIndex].monB < pscnumber)
                  {
                    std::cout << "Chargetransport" << std::endl;
                    excPos.h_location = h_viablePartners[g];
                  //excPos.location = excPos.location; //electron stays in place
                  }

                  /*End criteria for simulation*/

                  switch (plane)
                  {
                    case 'x':
                      if((excCoup[excPos.h_location].position.x() - avg.x()) > (0.75* (excCoup[startPind[i]].position.x() - avg.x())))
                      {
                        ch_separation[i]++;
                      }
                      break;
                    
                    case 'y':
                      if((excCoup[excPos.h_location].position.y() - avg.y()) > (0.75* (excCoup[startPind[i]].position.y() - avg.y())))
                      {

                      }
                      break;
                    
                    case 'z':
                      if((excCoup[excPos.h_location].position.z() - avg.z()) > (0.75* (excCoup[startPind[i]].position.z() - avg.z())))
                      {

                      }
                      break;

                  }


                 }
               }

            }

            viablePartners.clear();//empties vector containing possible partners for step so it can be reused in next step
            partnerConnections.clear();
            break;
          }// state c end
//SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
          else if (excPos.state == 's')//separated state
          {
            std::cout << "Successful run." << '\n';
            viablePartners.clear();//empties vector containing possible partners for step so it can be reused in next step
            partnerConnections.clear();
            excPos.state = 'e';
            break;
          }//state s end
//TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
          else if (excPos.state == 't')//termination state
          {
            std::cout << "Broken." << '\n';
            viablePartners.clear();//empties vector containing possible partners for step so it can be reused in next step
            partnerConnections.clear();
            excPos.state = 'e';
            break;
          }// state t end
//??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
          else
          {
            throw std::logic_error("Something somewhere went terribly wrong and the simulation ended up in an unknown state.");
          }

          viablePartners.clear();//empties vector containing possible partners for step so it can be reused in next step
          partnerConnections.clear();

          std::cout << "Exciton Position after Movement:  " << excPos.location << " Monomer A: " << excCoup[excPos.location].monA << " Monomer B: " << excCoup[excPos.location].monB << '\n';
          std::cout << "Step: " << h << '\n';
          std::cout << "#################################################################################" << '\n';
        }//loop for steps h
        std::cout << "Try: " << j << '\n';
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << '\n';
      }//100 try loop j
    }//loop over startingpoints i
    for (std::size_t i = 0u; i < startPind.size(); i++)
    {
      std::cout << "Radiating decays for Starting point " << i << ": " << radiating[i] << "." << '\n';
      std::cout << "Trappings for Starting point " << i << ": " << trapped[i] << "." << '\n';
    }
  }//try

  catch (std::exception & e)
  {
    std::cout << "An exception occured. The execution of " << config::Programname << " failed. \n";
    std::cout << "Error: " << e.what() << '\n';
  }
}