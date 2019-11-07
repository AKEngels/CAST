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
  double len = sqrt((pointA.x() - pointB.x()) * (pointA.x() - pointB.x()) + (pointA.y() - pointB.y()) * (pointA.y() - pointB.y()) + (pointA.z() - pointB.z()) * (pointA.z() - pointB.z()));
  return len;
}

double exciD::marcus(double coupling, double drivingF, double reorganisation)
{
  double marc;
  marc = (coupling * coupling) / exciD::h_quer * sqrt(M_PI / (reorganisation * exciD::boltzmann_const * 298)) * exp(-(reorganisation + drivingF) * (reorganisation + drivingF) / (4 * exciD::boltzmann_const * 298 * reorganisation));//replace 298 by user defined temperatuer?
  return marc;
}

double exciD::coulomb(coords::Cartesian_Point aktPos, coords::Cartesian_Point targetPos, double e_relative)
{
  double e_0 = 8.854187e-12;
  double elementar = 1.60217662e-19;
  double coulomb = -elementar / (4 * M_PI * e_0 * e_relative * length(aktPos, targetPos) * 1e-10);
  return coulomb;
}

coords::Cartesian_Point exciD::structCenter(coords::Representation_3D com)//calculate average position of the given positions (more than two)
{
  coords::Cartesian_Point avg(0.0, 0.0, 0.0);
  for (std::size_t i = 0u; i < com.size(); i++)
  {
    avg += com[i];
  }
  avg /= com.size();
  return avg;
}

coords::Cartesian_Point exciD::min(coords::Representation_3D coords)
{
  coords::Cartesian_Point min;

  min = coords[0];//initialize with coordinates of first point

  for (std::size_t i = 1u; i < coords.size(); i++)
  {
    if (coords[i].x() < min.x())
    {
      min.x() = coords[i].x();
    }

    if (coords[i].y() < min.y())
    {
      min.y() = coords[i].y();
    }

    if (coords[i].z() < min.z())
    {
      min.z() = coords[i].z();
    }
  }
  return min;
}

coords::Cartesian_Point exciD::max(coords::Representation_3D coords)
{
  coords::Cartesian_Point max;
  max = coords[0];// initialize witzh coordinates of first point

  for (std::size_t i = 1u; i < coords.size(); i++)
  {
    if (coords[i].x() > max.x())
    {
      max.x() = coords[i].x();
    }

    if (coords[i].y() > max.y())
    {
      max.y() = coords[i].y();
    }

    if (coords[i].z() > max.z())
    {
      max.z() = coords[i].z();
    }
  }
  return max;
}

void exciD::dimexc(std::string masscenters, std::string couplings, std::size_t pscnumber, int nscnumber, char interfaceorientation, double startingPscaling, std::size_t nbrStatingpoins) {
  try {
	
	std::string couplingsname;
    double reorganisationsenergie_exciton = Config::get().exbreak.ReorgE_exc;//noch extra variablen in config.h und config.cc einfügen
    double reorganisationsenergie_charge = Config::get().exbreak.ReorgE_ch;
    double reorganisationsenergie_nSC = Config::get().exbreak.ReorgE_nSC;
    double reorganisationsenergie_ct = Config::get().exbreak.ReorgE_ct;
    double reorganisationsenergie_rek = Config::get().exbreak.ReorgE_rek;
    double triebkraft_ct = Config::get().exbreak.ct_triebkraft;
    double oszillatorstrength = Config::get().exbreak.oscillatorstrength;
    double wellenzahl = Config::get().exbreak.wellenzahl;
    double k_rad = wellenzahl * wellenzahl * oszillatorstrength; // fluoreszenz

    char plane = interfaceorientation;//don't forget to replace by userinput
    
    couplingsname = couplings;


    std::ifstream comf;
    std::ifstream coupf;
    std::size_t numbermon;
    std::vector<std::size_t> startPind, viablePartners, h_viablePartners;
    std::vector<exciD::Partners> partnerConnections, h_partnerConnections;
    coords::Representation_3D com, pSCcom;
    coords::Cartesian_Point avg, pSCavg;
    exciD::Exciton excPos;


    comf.open(masscenters);
    comf >> numbermon;
    //comf.ignore(256, '\n');//ignore all till newline 

    comf.ignore(std::numeric_limits < std::streamsize >::max(), '\n');

    coords::Cartesian_Point tmp;
    std::string unnecessaryindex;
    std::string line;

    //read centers of mass
    for (std::size_t i = 0; i < numbermon; i++)
    {
      comf >> unnecessaryindex >> tmp;
      com.push_back(tmp);
    }

    //check if the number of masscenters matches the expected value
    std::cout << numbermon << " " << com.size() << '\n';
    if (numbermon != com.size())
    {
      throw std::logic_error("Unclear number of centers of mass.");
    }

    coupf.open(couplingsname);

    std::string tmpcount;
    std::size_t tmpA, tmpB;
    double tmpC, tmpD;
    exciD::Couplings tmpE;
    std::vector<exciD::Couplings> excCoup;
    double avgCoup{ 0.0 };
    

std::cout << "Couplings are read from: " << couplings << '\n';

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
      excCoup[i].position = (exciD::avgDimCoM(com[excCoup[i].monA - 1], com[excCoup[i].monB - 1]));
    }

    for (size_t i = 0; i < pscnumber; i++)
    {
      pSCcom.push_back(com[i]);
    }

    avg = exciD::structCenter(com);

    pSCavg = exciD::structCenter(pSCcom);

    coords::Cartesian_Point minV = min(com);
    coords::Cartesian_Point maxV = max(com);

    std::random_device rd; //prepare rng
    std::default_random_engine engine(rd());
    std::normal_distribution<double> distributionN(0.0, 0.068584577);
    std::uniform_real_distribution<double> distributionR(0, 1); //beispiel für rng var: double rng = distributionN(engine);

    //choose startingpoints
    for (;startPind.size() < nbrStatingpoins;)
    {
      switch (plane)
      {
      case 'x':

        if (pSCavg.x() > avg.x())
        {
          for (std::size_t i = 0u; i < excCoup.size(); i++)
          {
            if (excCoup[i].position.x() - avg.x() > startingPscaling * (maxV.x() - avg.x()))
            {
              startPind.push_back(i);
            }
          }
        }
        else
        {
          for (std::size_t i = 0u; i < excCoup.size(); i++)
          {
            if (excCoup[i].position.x() - avg.x() < startingPscaling * (minV.x() - avg.x()))
            {
              startPind.push_back(i);
            }
          }
        }

        break;

      case 'y':

        if (pSCavg.y() > avg.y())
        {
          for (std::size_t i = 0u; i < excCoup.size(); i++)
          {
            if (excCoup[i].position.y() - avg.y() > startingPscaling * (maxV.y() - avg.y()))
            {
              startPind.push_back(i);
            }
          }
        }
        else
        {
          for (std::size_t i = 0u; i < excCoup.size(); i++)
          {
            if (excCoup[i].position.y() - avg.y() < startingPscaling * (minV.y() - avg.y()))
            {
              startPind.push_back(i);
            }
          }
        }

        break;

      case 'z':

        if (pSCavg.z() > avg.z())
        {
          for (std::size_t i = 0u; i < excCoup.size(); i++)
          {
            if (excCoup[i].position.z() - avg.z() > startingPscaling * (maxV.z() - avg.z()))
            {
              startPind.push_back(i);
            }
          }
        }
        else
        {
          for (std::size_t i = 0u; i < excCoup.size(); i++)
          {
            if (excCoup[i].position.z() - avg.z() < startingPscaling * (minV.z() - avg.z()))
            {
              startPind.push_back(i);
            }
          }
        }

        break;
      }

      if (startPind.size() < nbrStatingpoins)
      {
        startingPscaling -= 0.01;
        startPind.clear(); //To ensure no startingpoint is used more than once

      }
    }

    std::cout << "Used Startingpoint scaling factor: " << startingPscaling << '\n';
    std::cout << "Number of Startingpoints: : " << startPind.size() << '\n';
    std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << '\n';

    if (startPind.size() == 0)
    {
      throw std::logic_error("No Points to start the simulation were found.");
    }

    //loop for writing starting points
    std::ofstream startPout;
    startPout.open("Startingpoints.txt");

    for (std::size_t i = 0u; i < startPind.size(); i++)
    {
      startPout << "Startingpoint " << i << ": " << startPind[i] << " Monomer A " << excCoup[startPind[i]].monA << " Monomer B " << excCoup[startPind[i]].monB << '\n';
    }
    startPout.close();

    std::vector <int> trapped(startPind.size(), 0);//for counting the trapped excitons unable to reach the interface from each startingpoint
    std::vector <int> ex_diss(startPind.size(), 0);
    std::vector <int> radiating(startPind.size(), 0);
    std::vector <int> rekombined(startPind.size(), 0);
    std::vector <int> ch_separation(startPind.size(), 0);

    std::vector <std::vector<double>> time_ch(startPind.size(), std::vector<double>(100, 0.)), //vectors to keep time/velocities for different startingpoints and tries
      time_ex(startPind.size(), std::vector<double>(100, 0.)),
      vel_ch(startPind.size(), std::vector<double>(100, 0.)),
      vel_ex(startPind.size(), std::vector<double>(100, 0.));

    double time(0.0), time_p(0.0), time_n(0.0);
    int const nbrof_tries = 100;
    int const nbrof_steps = 2 * startPind.size() + 400;


    std::ofstream viabP;
    viabP.open("viablePartners.txt");
    //loop over all startingpoints 
    for (std::size_t i = 0u; i < startPind.size(); i++)
    {


      //loop ensures to start 100 times from every startingpoint
      for (std::size_t j = 0u; j < nbrof_tries; j++)//don't forget to set to 100 when all works fine
      {

        excPos.location = startPind[i];//startPind[i] is the startingpoint for the actual simulation, excPos is used to keep track of the position of the exciton during simulation

        std::cout << "Startingpoint " << i << ": " << startPind[i] << " Monomer A: " << excCoup[startPind[i]].monA << " Monomer B: " << excCoup[startPind[i]].monB << " Coupling: " << excCoup[excPos.location].coupling << '\n';

        for (int h = 0; h < nbrof_steps; h++)//steps in each try for testing hardcode
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
                  viabP << k << '\n';
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
            std::cout << " Average Coupling: " << partnerConnections[n].avgCoup << " |" << "Average secCoupling: " << partnerConnections[n].avgsecCoup << '\n';
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

          double rate_sum(0.0), rateFul_sum(0.0), coulombenergy, rate_KMC(0.0), tmp_ratesum(0.0);
          std::vector <double> raten;//used for exciton and electron rates
          std::vector <double> raten_hole;//used for hole rates
          double random_normal, random_normal1;
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

            std::cout << "Fluorescence Rate: " << k_rad << '\n';

            rate_KMC = random_real * rate_sum;
            std::cout << "KMC-Rate: " << rate_KMC << '\n';

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

            if (viablePartners.size() == 0)//if no partners are found the kmc point for radiating decay must be reached
            {
              viablePartners.push_back(0);
              raten.push_back(0.0);//so the vector does not leave its reange in the following loop
            }

            for (std::size_t q = 0u; q < viablePartners.size(); q++)
            {
              if (raten[q] > rate_KMC)
              {
                if (excCoup[viablePartners[q]].monA >= pscnumber && excCoup[viablePartners[q]].monB >= pscnumber)
                {
                  excPos.h_location_lastS = excPos.h_location;
                  excPos.h_location = excPos.location; //set hole position to former exciton position for charge separation
                }

                excPos.location_lastS = excPos.location;
                excPos.location = viablePartners[q];//after chargeseparation excPos.location tracks electron position.

                if (excCoup[excPos.location].monA >= pscnumber && excCoup[excPos.location].monB >= pscnumber)
                {
                  excPos.state = 'c';
                  std::cout << "Chargeseparation." << '\n';
                  ex_diss[i]++;
                  excPos.state = 's';//till chargemovement is implemented
                  vel_ex[i][j] = length(excCoup[startPind[i]].position, excCoup[excPos.h_location_lastS].position) / time;
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
            for (std::size_t p = 0u; p < partnerConnections.size(); p++)
            {
              if (excCoup[partnerConnections[p].partnerIndex].monA < pscnumber && excCoup[partnerConnections[p].partnerIndex].monB < pscnumber)//movement on pSC
              {
                random_normal = distributionN(engine);//generating normal distributed random number
                coulombenergy = coulomb(excCoup[excPos.h_location].position, excCoup[partnerConnections[p].partnerIndex].position, 3.4088) - coulomb(excCoup[excPos.h_location].position, excCoup[excPos.location].position, 3.4088);
                rate_sum += marcus(partnerConnections[p].avgCoup, (random_normal - random_normal1) + coulombenergy, reorganisationsenergie_charge);
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
                rateFul_sum += marcus(partnerConnections[p].avgCoup, (random_normal - random_normal1) + coulombenergy, reorganisationsenergie_nSC);
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

            if ((1 / rate_sum - time_p) < (1 / rateFul_sum - time_n))
            {
              std::cout << "pSC hopps first. " << std::endl;

              //Update time
              if ((1 / rate_sum - time_p) > 0)
              {
                time += (1 / rate_sum - time_p);
                time_n += (1 / rate_sum - time_p);
                time_p = 0.;
              }
              else if ((1 / rate_sum - time_p) < 0)
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

              for (std::size_t g = 0; g < h_viablePartners.size(); g++)
              {
                if (raten_hole[g] > rate_KMC)
                {
                  if (excCoup[partnerConnections[g].partnerIndex].monA < pscnumber && excCoup[partnerConnections[g].partnerIndex].monB < pscnumber)
                  {
                    std::cout << "Chargetransport" << std::endl;
                    excPos.location_lastS = excPos.location;
                    excPos.h_location = h_viablePartners[g];

                    /*End criteria for simulation*/

                    switch (plane)
                    {
                    case 'x':
                      if ((excCoup[excPos.h_location].position.x() - avg.x()) > (0.75 * (excCoup[startPind[i]].position.x() - avg.x())))
                      {
                        ch_separation[i]++;
                        time_ch[i][h] = time - time_ex[i][h];
                        vel_ch[i][h] = (excCoup[excPos.h_location].position.x() - avg.x()) / time_ch[i][h];

                        std::cout << "Chargetransport" << std::endl;
                        excPos.state = 's';
                      }
                      break;

                    case 'y':
                      if ((excCoup[excPos.h_location].position.y() - avg.y()) > (0.75 * (excCoup[startPind[i]].position.y() - avg.y())))
                      {
                        ch_separation[i]++;
                        time_ch[i][h] = time - time_ex[i][h];
                        vel_ch[i][h] = (excCoup[excPos.h_location].position.y() - avg.y()) / time_ch[i][h];

                        std::cout << "Chargetransport" << std::endl;
                        excPos.state = 's';
                      }
                      break;

                    case 'z':
                      if ((excCoup[excPos.h_location].position.z() - avg.z()) > (0.75 * (excCoup[startPind[i]].position.z() - avg.z())))
                      {
                        ch_separation[i]++;
                        time_ch[i][h] = time - time_ex[i][h];
                        vel_ch[i][h] = (excCoup[excPos.h_location].position.z() - avg.z()) / time_ch[i][h];

                        std::cout << "Chargetransport" << std::endl;
                        excPos.state = 's';
                      }
                      break;
                    }//switch end

                    //#########################################################################################################################
                    std::cout << "old pSC Monomer " << std::setw(5) << excPos.h_location_lastS << std::endl;
                    std::cout << "new pSC Monomer " << std::setw(5) << excPos.h_location << std::endl;
                    std::cout << "Coupling pSC " << std::setw(5) << std::setprecision(6) << std::fixed << avgCoup << std::endl;//Welches Coupling war nochmal electronen relevant? avg oder svgsec?
                    std::cout << "nSC " << std::setw(5) << excPos.location << std::setw(5) << excPos.location_lastS << std::endl;
                    break;

                  }
                  else if (excCoup[partnerConnections[g].partnerIndex].monA > pscnumber && excCoup[partnerConnections[g].partnerIndex].monB > pscnumber)//movement of electron back onto pSC
                  {
                    std::cout << "Recombination" << std::endl;
                    excPos.state = 't';
                    rekombined[i]++;
                    break;
                  }
                }
                else if (g == h_viablePartners.size())
                {
                  throw std::logic_error("WARNING: Errot during pSC chargetransfer, no suiteable hopping target was found.");
                  return;
                }
              }

            }//end of pSC hopping

            //psc hopping
            else if ((1 / rate_sum - time_p) > (1 / rateFul_sum - time_n))
            {
              std::cout << "Hopping on nSC first." << std::endl;

              if ((1 / rateFul_sum - time_n) > 0)
              {
                time = time + (1 / rateFul_sum - time_n);
                time_p = time_p + ((1 / rateFul_sum - time_n));
                time_n = 0;
              }
              else if ((1 / rateFul_sum - time_n) < 0)
              {
                //time = time;
                //time_p = time_p;
                time_n = 0;
              }
              else
              {
                throw std::logic_error("Something went wrong with the time. Call the Doctor.");
              }

              //execute nSC hopping
              auto random_real = distributionR(engine);
              rate_KMC = random_real * rateFul_sum;

              for (std::size_t g = 0; g < viablePartners.size(); g++)
              {
                if ((raten[g] > rate_KMC))
                {
                  if (excCoup[partnerConnections[g].partnerIndex].monA > pscnumber && excCoup[partnerConnections[g].partnerIndex].monB > pscnumber)
                  {
                    std::cout << "Chargetransport in nSC" << std::endl;

                    excPos.h_location_lastS = excPos.h_location;
                    excPos.location = viablePartners[g];

                    //#########################################################################################################################
                    std::cout << "old nSC Monomer " << std::setw(5) << excPos.location_lastS << std::endl;
                    std::cout << "new nSC Monomer " << std::setw(5) << excPos.location << std::endl;
                    std::cout << "Coupling nSC " << std::setw(5) << std::setprecision(6) << std::fixed << avgCoup << std::endl;//Welches Coupling war nochmal electronen relevant? avg oder svgsec?
                    std::cout << "pSC " << std::setw(5) << excPos.h_location << std::setw(5) << excPos.h_location_lastS << std::endl;

                    break;
                  }
                  else if (excCoup[partnerConnections[g].partnerIndex].monA < pscnumber && excCoup[partnerConnections[g].partnerIndex].monB < pscnumber)
                  {
                    std::cout << "Recombination" << std::endl;
                    excPos.state = 't';
                    rekombined[i]++;
                    break;
                  }
                }
                else if (g == viablePartners.size())
                {
                  throw std::logic_error("WARNING: ERROR during transport in nSC");
                  return;
                }
              }



            }//end nSC hopping

            viablePartners.clear();//empties vector containing possible partners for step so it can be reused in next step
            h_viablePartners.clear();
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
            std::cout << "Broken.-" << '\n';
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

    viabP.close();
    //___________________________________EVALUATION_______________________________________________________________________________________

    std::ofstream evaluation;
    evaluation.open("evaluation.txt");
    evaluation << std::setw(4) << "k" << std::setw(4) << "IX" << std::setw(9) << "Ex_Diss" << std::setw(9) << "Ch_Diss" << std::setw(9) << "Rek." << std::setw(9) << "Trapp." << std::setw(9) << "Fluor." << std::endl;

    double avg_ex(0.), avg_ch(0.), avg_rek(0.), avg_trap(0.), avg_rad(0.);

    for (std::size_t k = 0; k < startPind.size(); k++)
    {
      avg_ex += ex_diss[k];
      avg_ch += ch_separation[k];
      avg_rek += rekombined[k];
      avg_trap += trapped[k];
      avg_rad += radiating[k];

      evaluation << std::setw(4) << k << std::setw(4) << startPind[k] << std::setw(9) << std::setprecision(5) << ex_diss[k]
        << std::setw(9) << std::setprecision(5) << ch_separation[k] << std::setw(9) << std::setprecision(5) << rekombined[k]
        << std::setw(9) << std::setprecision(5) << trapped[k] << std::setw(9) << std::setprecision(5) << radiating[k] << std::endl;
    }

    evaluation << std::setw(9) << "Average " << std::setw(9) << std::setprecision(5) << std::fixed << avg_ex / startPind.size()
      << std::setw(9) << std::setprecision(5) << std::fixed << avg_ch / startPind.size()
      << std::setw(9) << std::setprecision(5) << std::fixed << avg_rek / startPind.size()
      << std::setw(9) << std::setprecision(5) << std::fixed << avg_trap / startPind.size()
      << std::setw(9) << std::setprecision(5) << std::fixed << avg_rad / startPind.size() << std::endl;

    evaluation << "Velocities " << std::endl;
    evaluation << std::setw(4) << "k" << std::setw(5) << "IX" << std::setw(11) << "Ex_vel" << std::setw(11) << "Ex_s_dev" << std::setw(11) << "Ch_vel" << std::setw(11) << "Ch_s_dev" << std::endl;

    //average veloceties
    std::vector <double> avg_ex_vel, avg_ch_vel;
    std::vector <double> standDevEX, standDevCH;

    for (std::size_t k = 0; k < startPind.size(); k++)
    {
      double tmp_avg_ch_vel(0.);
      double tmp_avg_ex_vel(0.);
      for (std::size_t j = 0; j < nbrof_tries; j++)// sum up all veloceties for exciton and holes
      {

        if (vel_ch[k][j] > 0.0001)
        {
          tmp_avg_ch_vel += vel_ch[k][j];
        }

        if (vel_ex[k][j] > 0.0001)
        {
          tmp_avg_ex_vel += vel_ex[k][j];
        }
      }//j
      avg_ch_vel.push_back(tmp_avg_ch_vel);
      avg_ex_vel.push_back(tmp_avg_ex_vel);
      // divide summed up velocities by number of relevant events happened
      if (ch_separation[k] > 0)
      {
        avg_ch_vel[k] /= ch_separation[k];
      }
      else if (ch_separation[k] == 0)
      {
        avg_ch_vel[k] = 0;
      }

      if (ex_diss[k] > 0)
      {
        avg_ex_vel[k] /= ex_diss[k];
      }
      else if (ex_diss[k] == 0)
      {
        avg_ex_vel[k] = 0;
      }

      //calculation of standart deviation
      double tmpstandDevEX(0.), tmpstandDevCH(0.);

      for (std::size_t j = 0; j < nbrof_tries; j++)
      {
        if (vel_ch[k][j] > 0.0001)
        {
          tmpstandDevCH += ((vel_ch[k][j] - avg_ch_vel[k]) * (vel_ch[k][j] - avg_ch_vel[k]));
        }

        if (vel_ex[k][j] > 0.0001)
        {
          tmpstandDevEX += ((vel_ex[k][j] - avg_ex_vel[k]) * (vel_ex[k][j] - avg_ex_vel[k]));
        }
      }//j
      standDevCH.push_back(tmpstandDevCH);
      standDevEX.push_back(tmpstandDevEX);

      if (ch_separation[k] > 1)
      {
        standDevCH[k] = sqrt(standDevCH[k] / ch_separation[k]);
      }
      else if (ch_separation[k] < 2)
      {
        standDevCH[k] = 0;
      }

      if (ex_diss[k] > 1)
      {
        standDevEX[k] = sqrt(standDevEX[k] / ex_diss[k]);
      }
      else if (ex_diss[k] < 2)
      {
        standDevEX[k] = 0;
      }

      evaluation << std::setw(4) << k << std::setw(5) << startPind[k] << std::setw(11) << std::setprecision(5) << std::fixed << avg_ex_vel[k] * 1e-9 <<
        std::setw(11) << std::setprecision(5) << std::fixed << standDevEX[k] * 1e-9 <<
        std::setw(11) << std::setprecision(5) << std::fixed << avg_ch_vel[k] * 1e-9 <<
        std::setw(11) << std::setprecision(5) << std::fixed << standDevCH[k] * 1e-9 << std::endl;
    }//k

    double mean_vel_ex(0.), mean_vel_ch(0.);
    for (std::size_t k = 0; k < startPind.size(); k++)
    {
      mean_vel_ex += avg_ex_vel[k];
      mean_vel_ch += avg_ch_vel[k];
    }

    mean_vel_ex /= startPind.size();
    mean_vel_ch /= startPind.size();

    evaluation << std::left << std::setw(7) << "Average    " << std::left << std::setw(22) << std::setprecision(5) << std::fixed << mean_vel_ex * 1e-9 <<
      std::left << std::setw(9) << std::setprecision(5) << std::fixed << mean_vel_ch * 1e-9 << std::endl;

    //distribution of chargecarrier and exciton-velocities
    std::ofstream distribution;
    distribution.open("Exciton_Distribution.txt");
    int nmbr;

    for (std::size_t i = 0; i < 20; i++)
    {
      nmbr = 0;

      for (std::size_t k = 0; k < startPind.size(); k++)
      {
        for (std::size_t j = 0; j < nbrof_tries; j++)
        {
          if ((vel_ex[k][j] > ((i + 1) * 50 * 1e9)))
          {
            nmbr++;
          }
        }
      }
      distribution << std::setw(9) << std::setprecision(5) << (i + 1) * 50 << std::setw(9) << nmbr / startPind.size() << std::endl;
    }
    distribution.close();

    distribution.open("Charge_Distribution.txt");

    for (std::size_t i = 0; i < 20; i++)
    {
      nmbr = 0;

      for (std::size_t k = 0; k < startPind.size(); k++)
      {
        for (std::size_t j = 0; j < nbrof_tries; j++)
        {
          if (vel_ch[k][j] > ((i + 1) * 50 * 1e9))
          {
            nmbr++;
          }
        }
      }
      distribution << std::setw(9) << std::setprecision(5) << (i + 1) * 50 << std::setw(9) << nmbr / (startPind.size()) << std::endl;
    }

    distribution.close();

  }//try
  catch (std::exception& e)
  {
    std::cout << "An exception occured. The execution of " << config::Programname << " failed. \n";
    std::cout << "Error: " << e.what() << '\n';
  }
}
