#include "Couplings.h"

void couplings::coupling::kopplung()
{
  int gesanzahl_monomere = Config::get().couplings.nbr_nSC + Config::get().couplings.nbr_pSC;
  std::unique_ptr<coords::input::format> ci(coords::input::new_format());

  for (int i = 1; i < gesanzahl_monomere-1; i++)
  {
std::ofstream i_loop("i_loop.txt", std::ios::out);
i_loop << i ;
i_loop.close();

    for (int j = 2; j < gesanzahl_monomere; j++)
    {

std::ofstream j_loop("j_loop.txt", std::ios::out);
j_loop << j;
j_loop.close();

std::ofstream debug("debug.txt", std::ios::out);
std::ofstream test("test.txt", std::ios::out);

      std::stringstream idatname;
      idatname << "Dimerstrukt_" << i << "_" << j << ".xyz";

debug << idatname.str() << '\n';
    
      std::ifstream coord_test(idatname.str(), std::ios_base::in);

      if(coord_test) //there will be names for dimerpairs generated whom not exist so these errors shall be caught within the loop
      {
test << '4' << " " << i << " " << j << '\n';
debug << '4';



        coords::Coordinates dim_coords(ci->read(idatname.str()));

        coords::Coordinates *ptr_dim_coords = &dim_coords;

        if (i < Config::get().couplings.nbr_pSC && j < Config::get().couplings.nbr_pSC)//pSC homo-pair
        {
debug << '5';

test << '5' ;

test.close();
          pSC_homo_1.push_back(i);
          pSC_homo_2.push_back(j);

          INDO(dim_coords);

          double V_el_tmp = 0.5*(c_virtMO[1] - c_virtMO[0]);

test << '6' << '\n';
test << V_el_tmp << " " << c_virtMO[1] << " " << c_virtMO[0] << '\n';


          V_el.push_back(V_el_tmp);

test << V_el.back() << '\n';
test << '7';


          ZINDO(dim_coords);

test << '8';


          V_ex.push_back(0.5*(c_excitE[1] - c_excitE[0]));

        }//pSC homo-pair end

        if ( i >= Config::get().couplings.nbr_pSC && j >= Config::get().couplings.nbr_pSC) //nSC homo-pair
        {
debug << '6';

          nSC_homo_1.push_back(i);
          nSC_homo_2.push_back(j);

          INDO(dim_coords);

          V_hole.push_back(0.5*(c_occMO[0] - c_occMO[1]));

        }//nSC homo-pair end

        if (i < Config::get().couplings.nbr_pSC && j >= Config::get().couplings.nbr_pSC)//hetero-pair i pSC, j nSC  
        {
debug << '7';

          hetero_pSC.push_back(i);
          hetero_nSC.push_back(j);

          ZINDO(dim_coords);

          coords::Cartesian_Point monom1, monom2, dipol_ct;

          monom1 = dim_coords.center_of_mass_mol(0);//for molecules without static dipolemoment we use the masscenter for the dipolemoment
          monom2 = dim_coords.center_of_mass_mol(1);

          for (int i = 0; i < 3; i++) //calculation dipolemoment for dimer for unpolar monomers
          {
debug << '8';

            dipol_ct.x() = 0.5 * monom1.x() - 0.5 * monom2.x();
            dipol_ct.y() = 0.5 * monom1.y() - 0.5 * monom2.y();
            dipol_ct.z() = 0.5 * monom1.z() - 0.5 * monom2.z();
          }

          double dipolemoment = sqrt(dipol_ct.x()*dipol_ct.x() + dipol_ct.y()*dipol_ct.y() + dipol_ct.z()*dipol_ct.z());//length of total dipolmoment

          std::stringstream string_ct_relev_states(Config::get().couplings.ct_chara_all);
          std::vector<int> ct_relev_states;
          int ct_state (0);
          double projection(0), coupling (0), ct_square_coup_sum(0), rek_square_coup_sum(0);
          std::vector <double> ct_coupling, rek_coupling;
          double a_u (0.52917721067);//conversion factor

          while (string_ct_relev_states >> ct_state){ct_relev_states.push_back(ct_state); }//all ct_states relevant to the calculation are bundeled in a vector of ints


          //CALCULATION FOR CT-COUPLINGS########################################################################################################################
          for (int i = 0; i < c_ex_ex_trans.size(); i++)//loop over all ex_ex_dipoles
          {
debug << '9';
            if (c_state_j[i] == 1)//ensuring unly dipolemoments concering the first excited state are used
            {
              for (int j = 2; j < ct_relev_states.size(); j++)//loop over user defined relevant ct-states
              {
                if (c_state_i[i] == ct_relev_states[j])//only if the dipolemoment is concering a relevant state
                {
                  projection = dipol_ct.x() / dipolemoment * c_ex_ex_trans[i].x()
                             + dipol_ct.y() / dipolemoment * c_ex_ex_trans[i].y()
                             + dipol_ct.z() / dipolemoment * c_ex_ex_trans[i].z();

                  coupling = (projection * (c_excitE[ct_relev_states[j]-1] - c_excitE[0])) / sqrt((dipolemoment/a_u)*(dipolemoment/a_u) + 4* projection * projection);
                  ct_coupling.push_back(coupling);
                }//end if-clause for relevant states
              }//end loop over relevant ct-states
            }//end if-clause ensuring first excited state
          }//end loop over ex_ex_dipoles

          //CALCULATION FOR REK-COUPLINGS##########################################################################################################################
          for (int i = 0; i < c_gz_ex_trans.size(); i++)//loop over all gz_ex_dipoles
          {
debug << 'a';
            for (int j = 2; j < ct_relev_states.size(); j++)//loop over user defined relevant ct-states
            {
              if (c_gz_i_state[i] == ct_relev_states[j])//only if the dipolemoment is concering a relevant state
              {
                projection = dipol_ct.x() / dipolemoment * c_gz_ex_trans[i].x()
                           + dipol_ct.y() / dipolemoment * c_gz_ex_trans[i].y()
                           + dipol_ct.z() / dipolemoment * c_gz_ex_trans[i].z();

                coupling = (projection * (c_excitE[ct_relev_states[j] - 1] )) / sqrt((dipolemoment / a_u)*(dipolemoment / a_u) + 4 * projection * projection);
                rek_coupling.push_back(coupling);
              }//end if-clause for relevant states
            }//end loop over relevant ct-states
          }//end loop over ex_ex_dipoles

          for (int i = 0; i < ct_relev_states.size(); i++) //sum up squares of couplings between single states
          {
            ct_square_coup_sum  += ct_coupling[i]  * ct_coupling[i];
            rek_square_coup_sum += rek_coupling[i] * rek_coupling[i];
          }

          V_ct.push_back(sqrt(ct_square_coup_sum));//put the coupling for the dimer in the vector
          V_rek.push_back(sqrt(rek_square_coup_sum));
         
        }//hetero end
        delete ptr_dim_coords;
      }
      
debug << 'b';

      //WRITING CACULATED COUPLINGS#####################################################

      write();

debug << 'c';
debug.close();
    }

  }

}

void couplings::coupling::INDO(coords::Coordinates coords) //Funktion for INDO-Calculation for marcus-theorie couplings
{
  //Change gasussian parameters to the needed settings
  Config::set().energy.gaussian.method = "INDO";
  Config::set().energy.gaussian.basisset = " ";
  Config::set().energy.gaussian.spec = " ";

  coords.e();

  coords.get_catch_interface();

  c_occMO = coords.catch_interface->get_occMO();
  c_virtMO = coords.catch_interface->get_virtMO();
}

void couplings::coupling::ZINDO(coords::Coordinates coords)//Funktion for ZINDO-Calculation for marcus-theorie couplings
{
  //Change gasussian parameters to the needed settings
  Config::set().energy.gaussian.method = "ZINDO TD=(NStates=15, Singlets,AllTransitiondensities, ListWindow) gfinput IOP(6/7=3)";
  Config::set().energy.gaussian.basisset = " ";
  Config::set().energy.gaussian.spec = " ";

  coords.e();

  coords.get_catch_interface();

  c_excitE = coords.catch_interface->get_excitE();
  c_ex_ex_trans = coords.catch_interface->get_ex_ex_trans();
  c_gz_ex_trans = coords.catch_interface->get_gz_ex_trans();
  c_state_i = coords.catch_interface->get_state_i();
  c_state_j = coords.catch_interface->get_state_j();
  c_gz_i_state = coords.catch_interface->get_gz_i_state();

}

void couplings::coupling::write()
{
  std::ofstream all_couplings("couplings.txt", std::ios::out);

  for (int i = 0; i < pSC_homo_1.size(); i++)
  {
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << pSC_homo_1[i] << pSC_homo_2[i] << V_el[i] << V_ex[i] << '\n';
  }

  for (int i = 0; i < hetero_pSC.size(); i++)
  {
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << hetero_pSC[i] << hetero_nSC[i] << V_ct[i] << V_rek[i] << '\n';
  }

  for (int i = 0; i < nSC_homo_1.size(); i++)
  {
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << nSC_homo_1[i] << nSC_homo_2[i] << V_hole[i] << '\n';
  }

  all_couplings.close();
}
