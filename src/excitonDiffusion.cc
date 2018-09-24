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

void exciD::dimexc(std::string masscenters, std::string couplings, double pscnumber, double nscnumber) {
  try {

    double reorganisationsenergie_exciton = Config::get().exbreak.ReorgE_exc;//noch extra variablen in config.h und config.cc einfügen
    double oszillatorstrength = Config::get().exbreak.oscillatorstrength;
    double wellenzahl = Config::get().exbreak.wellenzahl;
    double krad = wellenzahl * wellenzahl * oszillatorstrength; // fluoreszenz

    char plane = 'z';//don't forget to replace by userinput

    std::ifstream comf;
    std::ifstream coupf;
    std::size_t numbermon;
    std::vector<std::size_t> startPind, viablePartners;
    std::vector<exciD::Partners> partnerConnections;
    coords::Representation_3D com;
    coords::Cartesian_Point avg;
    exciD::Exciton excPos;

    comf.open(masscenters);
    comf >> numbermon;
    comf.ignore(256, '\n');//ignore all till newline

    coords::Cartesian_Point tmp;
    int unnecessaryindex;

    //read centers of mass
    while (!comf.eof())
    {
      comf >> unnecessaryindex >> tmp;
      com.push_back(tmp);
    }

    //check if the number of masscenters matches the expected value
    if (numbermon != com.size())
    {
      throw std::logic_error("Unclear number of centers of mass.");
    }

    coupf.open(couplings);

    std::size_t tmpA, tmpB;
    double tmpC;
    exciD::Couplings tmpE;
    std::vector<exciD::Couplings> excCoup;
    double avgCoup;

    //read monomer indices of dimers and corresponding exciton-coupling and save in vector of couplings
    while (!coupf.eof())
    {
      coupf >> tmpA >> tmpB >> tmpC;
      exciD::Couplings tmpE(tmpA, tmpB, tmpC);
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

        //loop over all startingpoints 
    for (std::size_t i = 0u; i < startPind.size(); i++)
    {


      //loop ensures to start 100 times from every startingpoint
      //      for (std::size_t j = 0u; j < 100; j++)
       //     {

      excPos.location = startPind[i];//startPind[i] is the startingpoint for the actual simulation, excPos is used to keep track of the position of the exciton during simulation

      std::cout << "Startingpoint " << i << ": " << startPind[i] << " Monomer A: " << excCoup[startPind[i]].monA << " Monomer B: " << excCoup[startPind[i]].monB << " Coupling: " << excCoup[excPos.location].coupling << '\n';

      //    for (std::size_t h = 0u; h < 100; h++)//steps in each try for testing hardcode at 100
       //   {

      for (std::size_t k = 0; k < excCoup.size(); k++)//loop over all dimerpairs for viable partners
      {

        //in this i-cause the viable partners to the actual location are determined and the dimer pairs necessary to calculate the average couplings between the dimers are determined.
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

      for (std::size_t n = 0u; n < partnerConnections.size(); n++)
      {
        partnerConnections[n].avgCoup = 0.0;

        for (std::size_t o = 0u; o < partnerConnections[n].connect.size(); o++)
        {
          partnerConnections[n].avgCoup += excCoup[partnerConnections[n].connect[o]].coupling;
        }// o
        //partnerConnections[w].avgCoup /= partnerConnections[w].connect.size();
        partnerConnections[n].avgCoup /= 4; //unsure if the average coupling is gained by dividing the sum of relevant couplings by their number or the maximum (4) number of relevant couplings (assuming not added couplings are zero).
      }// n

      for (std::size_t n = 0u; n < partnerConnections.size(); n++)
      {
        std::cout << "PartnerIndex: " << partnerConnections[n].partnerIndex << " Monomers: " << excCoup[partnerConnections[n].partnerIndex].monA << " " << excCoup[partnerConnections[n].partnerIndex].monB << '\n';
        for (std::size_t o = 0u; o < partnerConnections[n].connect.size(); o++)
        {
          std::cout << " ConnectorIndex: " << partnerConnections[n].connect[o] << " Monomers: " << excCoup[partnerConnections[n].connect[o]].monA << " " << excCoup[partnerConnections[n].connect[o]].monB 
            << " Coupling: " << excCoup[partnerConnections[n].connect[o]].coupling <<  " |" << '\n';
        } //o
        std::cout << " Average Coupling: " << partnerConnections[n].avgCoup << '\n';
      }// n

      //algorithm for exciton movement starts here all before was preparation to know where movement to is possible

      std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << '\n';

      if (excPos.state == 'e')//Exciton state
      {
        std::vector <double> raten;
       
        for (std::size_t p = 0u; p < viablePartners.size(); p++)
        {
          if (viablePartners[p])
        }// p

      }//state e
      else if (excPos.state == 'c')//charge separated state
      {

      }// state c
      else if(excPos.state == 's')//separated state
      {

      }//state s
      else if (excPos.state == 't')//termination state
      {

      }// state t
      else
      {
        throw std::logic_error("Something somewhere went terribly wrong and the simulation ended up in an unknown state.");
      }

      viablePartners.clear();//empties vector containing possible partners for step so it can be reused in next step
      partnerConnections.clear();


     //      }//loop for steps h
      //  }//100 try loop j
    }//loop over startingpoints i
  }//try

  catch (std::exception & e)
  {
    std::cout << "An exception occured. The execution of " << config::Programname << " failed. \n";
    std::cout << "Error: " << e.what() << '\n';
  }
}