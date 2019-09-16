#include "Couplings.h"

void couplings::coupling::calculateAndWriteToFile()
{
  using namespace constants;
  if (Config::get().general.energy_interface != config::interface_types::T::GAUSSIAN)
  {
    throw std::runtime_error("Wrong energy interface detected, the COUPLINGS task can only be used with the GAUSSIAN interface.");
  }

  std::size_t const gesanzahl_monomere = Config::get().couplings.nbr_nSC + Config::get().couplings.nbr_pSC;


  std::string inFilename_string;

  Config::set().energy.gaussian.basisset = " ";
  Config::set().energy.gaussian.spec = " ";



  for (std::size_t i = 1u; i < gesanzahl_monomere; i++)//Iterator for first monomer
  {

    for (std::size_t j = 2u; j <= gesanzahl_monomere; j++)//Iterator for second monomer
    {

      std::stringstream idatname;
      idatname << "Dimerstrukt_" << i << "_" << j << ".xyz";

      std::ofstream a;
      a.open("a.txt");
      a << i << "_" << j;
      a.close();

      std::ifstream coord_test(idatname.str(), std::ios_base::in);

      if (coord_test) //there will be names for dimerpairs generated that dont exist 
      {
        std::unique_ptr<coords::input::format> ci(coords::input::new_format());
        coords::Coordinates dim_coords(ci->read(idatname.str()));

        //CALCULATION FOR p-SC########################################################################################################################
        if (i <= Config::get().couplings.nbr_pSC && j <= Config::get().couplings.nbr_pSC)//pSC homo-pair
        {
          pSC_homo_1.push_back(i);
          pSC_homo_2.push_back(j);

          INDO(dim_coords, Config::get().couplings.pSCmethod_el, Config::get().couplings.pSCmultipl, Config::get().couplings.pSCcharge);

          V_hole.push_back(0.5 * (c_occMO[0] - c_occMO[1]) / au2kcal_mol);

          ZINDO(dim_coords, Config::get().couplings.pSCmethod_ex, Config::get().couplings.pSCmultipl, Config::get().couplings.pSCcharge);

          V_ex.push_back(0.5 * (c_excitE[1] - c_excitE[0]) / eV2kcal_mol);

        }//pSC homo-pair end


         //CALCULATION FOR n-SC########################################################################################################################
        if (i > Config::get().couplings.nbr_pSC && j > Config::get().couplings.nbr_pSC) //nSC homo-pair
        {

          nSC_homo_1.push_back(i);
          nSC_homo_2.push_back(j);

          INDO(dim_coords, Config::get().couplings.nSCmethod, Config::get().couplings.nSCmultipl, Config::get().couplings.nSCcharge);

          V_el.push_back(0.5 * (c_virtMO[1] - c_virtMO[0]) / au2kcal_mol);

        }//nSC homo-pair end

         //CALCULATION FOR HETERO-PAIR########################################################################################################################
        if (i <= Config::get().couplings.nbr_pSC && j > Config::get().couplings.nbr_pSC)//hetero-pair i pSC, j nSC  
        {

          hetero_pSC.push_back(i);
          hetero_nSC.push_back(j);

          ZINDO(dim_coords, Config::get().couplings.hetmethod, Config::get().couplings.hetmultipl, Config::get().couplings.hetcharge);

          coords::Cartesian_Point monom1, monom2, dipol_ct;

          monom1 = dim_coords.center_of_mass_mol(0);//for molecules without static dipolemoment we use the masscenter for the dipolemoment
          monom2 = dim_coords.center_of_mass_mol(1);



          for (std::size_t d = 0u; d < 3u; d++) //calculation dipolemoment for dimer for unpolar monomers
          {
            dipol_ct.x() = 0.5 * monom1.x() - 0.5 * monom2.x();
            dipol_ct.y() = 0.5 * monom1.y() - 0.5 * monom2.y();
            dipol_ct.z() = 0.5 * monom1.z() - 0.5 * monom2.z();
          }

          double const dipolemoment = sqrt(dipol_ct.x() * dipol_ct.x() + dipol_ct.y() * dipol_ct.y() + dipol_ct.z() * dipol_ct.z());//length of total dipolmoment

          std::stringstream string_ct_relev_states(Config::get().couplings.ct_chara_all);
          std::vector<std::size_t> ct_relev_states;
          std::size_t ct_state(0u);
          std::vector <double> ct_coupling, rek_coupling;
          const double a_u(0.52917721067);//conversion factor

          while (string_ct_relev_states >> ct_state)
          {
            ct_relev_states.push_back(ct_state);
          }//all ct_states relevant to the calculation are bundeled in a vector of ints

          //CALCULATION FOR CT-COUPLINGS########################################################################################################################
          for (auto c = 0u; c < c_ex_ex_trans.size(); c++)//loop over all ex_ex_dipoles
          {
            if (c_state_j[c] == 1)//ensuring unly dipolemoments concering the first excited state are used
            {
              for (auto d = 0u; d < ct_relev_states.size(); d++)//loop over user defined relevant ct-states
              {
                if (c_state_i[c] == (int)ct_relev_states[d])//only if the dipolemoment is concering a relevant state
                {
                  const double projection = dipol_ct.x() / dipolemoment * c_ex_ex_trans[c].x()
                    + dipol_ct.y() / dipolemoment * c_ex_ex_trans[c].y()
                    + dipol_ct.z() / dipolemoment * c_ex_ex_trans[c].z();

                  const double coupling = (projection * (c_excitE[0] - c_excitE[ct_relev_states[d] - 1]) / eV2kcal_mol) / sqrt((dipolemoment / a_u) * (dipolemoment / a_u) + 4 * projection * projection);//swaped the c_excitE's 26.02.19
                  ct_coupling.push_back(coupling);
                }//end if-clause for relevant states
              }//end loop over relevant ct-states
            }//end if-clause ensuring first excited state
          }//end loop over ex_ex_dipoles


           //CALCULATION FOR REK-COUPLINGS##########################################################################################################################
          for (auto c = 0u; c < c_gz_ex_trans.size(); c++)//loop over all gz_ex_dipoles
          {
            for (auto d = 0u; d < ct_relev_states.size(); d++)//loop over user defined relevant ct-states
            {
              if (c_gz_i_state[c] == (int)ct_relev_states[d])//only if the dipolemoment is concering a relevant state
              {
                const double projection = dipol_ct.x() / dipolemoment * c_gz_ex_trans[d].x()
                  + dipol_ct.y() / dipolemoment * c_gz_ex_trans[d].y()
                  + dipol_ct.z() / dipolemoment * c_gz_ex_trans[d].z();

                const double coupling = (projection * (c_excitE[ct_relev_states[d] - 1]) / eV2kcal_mol) / sqrt((dipolemoment / a_u) * (dipolemoment / a_u) + 4 * projection * projection);
                rek_coupling.push_back(coupling);
              }//end if-clause for relevant states
            }//end loop over relevant ct-states
          }//end loop over ex_ex_dipoles

          double ct_square_coup_sum(0), rek_square_coup_sum(0);
          for (std::size_t j = 0u; j < ct_relev_states.size(); j++) //sum up squares of couplings between single states
          {
            ct_square_coup_sum += ct_coupling[j] * ct_coupling[j];
            rek_square_coup_sum += rek_coupling[j] * rek_coupling[j];
          }

          V_ct.push_back(sqrt(ct_square_coup_sum));//put the coupling for the dimer in the vector
          V_rek.push_back(sqrt(rek_square_coup_sum));

        }//hetero end

      }//end if-coord_test

    }//end for j

  }//end for i
   //WRITING CACULATED COUPLINGS#####################################################
  write();
}

void couplings::coupling::INDO(coords::Coordinates& coords, std::string method, std::string multiplicity, std::string charge) //Funktion for INDO-Calculation for marcus-theorie couplings
{
  //Change gasussian parameters to the needed settings
  Config::set().energy.gaussian.method = method;
  Config::set().energy.gaussian.basisset = " ";
  Config::set().energy.gaussian.spec = " ";
  Config::set().energy.gaussian.charge = charge;
  Config::set().energy.gaussian.multipl = multiplicity;

  coords.e();

  /* coords.get_catch_interface();*/

  c_occMO = coords.catch_interface->get_occMO();
  c_virtMO = coords.catch_interface->get_virtMO();
}

void couplings::coupling::ZINDO(coords::Coordinates& coords, std::string method, std::string multiplicity, std::string charge)//Funktion for ZINDO-Calculation for marcus-theorie couplings
{
  //Change gasussian parameters to the needed settings
  Config::set().energy.gaussian.method = method;
  Config::set().energy.gaussian.basisset = " ";
  Config::set().energy.gaussian.spec = " ";
  Config::set().energy.gaussian.charge = charge;
  Config::set().energy.gaussian.multipl = multiplicity;

  coords.e();

  c_excitE = coords.catch_interface->get_excitE();
  c_ex_ex_trans = coords.catch_interface->get_ex_ex_trans();
  c_gz_ex_trans = coords.catch_interface->get_gz_ex_trans();
  c_state_i = coords.catch_interface->get_state_i();
  c_state_j = coords.catch_interface->get_state_j();
  c_gz_i_state = coords.catch_interface->get_gz_i_state();
}

void couplings::coupling::write()
{
  std::ofstream all_couplings("Couplings.txt", std::ios::out);
  //NEW OUTPUT FORMAT
  for (auto i = 0u; i < pSC_homo_1.size(); i++)
  {
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << pSC_homo_1[i] << " ";
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << pSC_homo_2[i] << " ";
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << V_hole[i] << " ";
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << V_ex[i] << '\n';
  }

  for (auto i = 0u; i < hetero_pSC.size(); i++)
  {
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << hetero_pSC[i] << " ";
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << hetero_nSC[i] << " ";
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << V_ct[i] << " ";
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << V_rek[i] << '\n';
  }

  for (auto i = 0u; i < nSC_homo_1.size(); i++)
  {
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << nSC_homo_1[i] << " ";
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << nSC_homo_2[i] << " ";
    all_couplings << std::setw(12) << std::setprecision(6) << std::fixed << std::left << V_el[i] << '\n';
  }

  all_couplings.close();
  // LEGACY OUTPUT FORMAT
  std::ofstream homo_ex("homodimer_exciton.txt", std::ios::out);
  for (auto i = 0u; i < pSC_homo_1.size(); i++)
  {
    homo_ex << std::setw(12) << std::setprecision(6) << std::fixed << std::left << pSC_homo_1[i] << " ";
    homo_ex << std::setw(12) << std::setprecision(6) << std::fixed << std::left << pSC_homo_2[i] << " ";
    homo_ex << std::setw(12) << std::setprecision(6) << std::fixed << std::left << V_ex[i] << '\n';
  }
  homo_ex.close();

  std::ofstream homo_ch("homodimer_ladung.txt", std::ios::out);
  for (auto i = 0u; i < pSC_homo_1.size(); i++)
  {
    homo_ch << std::setw(12) << std::setprecision(6) << std::fixed << std::left << pSC_homo_1[i] << " ";
    homo_ch << std::setw(12) << std::setprecision(6) << std::fixed << std::left << pSC_homo_2[i] << " ";
    homo_ch << std::setw(12) << std::setprecision(6) << std::fixed << std::left << V_hole[i] << '\n';
  }
  homo_ch.close();

  std::ofstream hetero("heterodimer.txt", std::ios::out);
  for (auto i = 0u; i < hetero_pSC.size(); i++)
  {
    hetero << std::setw(12) << std::setprecision(6) << std::fixed << std::left << hetero_pSC[i] << " ";
    hetero << std::setw(12) << std::setprecision(6) << std::fixed << std::left << hetero_nSC[i] << " ";
    hetero << std::setw(12) << std::setprecision(6) << std::fixed << std::left << V_ct[i] << " ";
    hetero << std::setw(12) << std::setprecision(6) << std::fixed << std::left << V_rek[i] << '\n';
  }
  hetero.close();

  std::ofstream nSC("nSC_homodimer.txt", std::ios::out);
  for (auto i = 0u; i < nSC_homo_1.size(); i++)
  {
    nSC << std::setw(12) << std::setprecision(6) << std::fixed << std::left << nSC_homo_1[i] << " ";
    nSC << std::setw(12) << std::setprecision(6) << std::fixed << std::left << nSC_homo_2[i] << " ";
    nSC << std::setw(12) << std::setprecision(6) << std::fixed << std::left << V_el[i] << '\n';
  }
}

