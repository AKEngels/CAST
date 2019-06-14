/**
CAST 3
Purpose: find and name aminoacids, determine OPLSAA atomtypes, especially for amino acids

@author Susanne Sauer
@version 1.0
*/
#include "coords_io.h"
#include "helperfunctions.h"

void coords::AtomtypeFinder::get_some_easy_atomtypes()
{
  for (auto i{ 0u }; i < atoms.size(); ++i)
  {
    auto &a = atoms.atom(i);
    if (a.symbol() == "Na" && a.bonds().size() == 0)         // sodium ion
    {
      a.set_energy_type(349);
      got_it[i] = true;
    }
    else if (a.symbol() == "O")                              // water molecule
    {
      if (a.bonds().size() == 2 && atoms.atom(a.bonds()[0]).symbol() == "H" && atoms.atom(a.bonds()[1]).symbol() == "H")
      {
        a.set_energy_type(63);
        atoms.atom(a.bonds()[0]).set_energy_type(64);
        atoms.atom(a.bonds()[1]).set_energy_type(64);

        got_it[i] = true;
        got_it[a.bonds()[0]] = true;
        got_it[a.bonds()[1]] = true;
      }
    }
  }
}

std::vector<coords::AminoAcid> coords::AtomtypeFinder::get_aminoacids()
{
  std::vector<coords::AminoAcid> amino_acids;

  for (auto i{ 0u }; i < atoms.size(); ++i)     // for all atoms
  {
    auto a = atoms.atom(i);
    if (a.symbol() == "O" && got_it[i] == false)      // O atom that is not in an amino acid yet
    {
      if (a.bonds().size() == 1 && atoms.atom(a.bonds()[0]).symbol() == "C")  // if only bound to C (=> carbonyle)
      {
        auto j = a.bonds()[0];
        auto b = atoms.atom(j);
        if (got_it[j] == false)                      // if carbonyle C is not in an amino acid yet (could be for a terminal amino acid)
        {
          std::vector<std::string> symbolvec_b = get_bonding_symbols(b, atoms);
          auto index = find_index("C", symbolvec_b);
          if (index < std::numeric_limits<int>::max())                        // if other C atom is bound to carbonyle C (=> C alpha)
          {
            auto k = b.bonds()[index];
            auto c = atoms.atom(k);
            if (c.bonds().size() == 4)
            {
              auto symbolvec_c = get_bonding_symbols(c, atoms);
              index = find_index("N", symbolvec_c);
              if (index < std::numeric_limits<int>::max())                    // if this C is bound to N 
              {
                auto l = c.bonds()[index];
                auto d = atoms.atom(l);

                // set all 4 atoms to got_it
                got_it[i] = true;      
                got_it[j] = true;
                got_it[k] = true;
                got_it[l] = true;

                // determine terminal state
                auto terminal = coords::terminalState::no;
                if (count_element("O", symbolvec_b) == 2) terminal = coords::terminalState::C;
                if (count_element("H", get_bonding_symbols(d, atoms)) > 1) terminal = coords::terminalState::N;

                // create amino acid and add it to vector
                coords::AminoAcid as({ i, j, k, l }, terminal);
                amino_acids.emplace_back(as);
              }
            }
          }
        }
      }
    }
  }
  complete_atoms_of_aminoacids(amino_acids);
  return amino_acids;
}

void coords::AtomtypeFinder::add_bonds_to_as(int index, AminoAcid &as)
{
  for (auto b : atoms.atom(index).bonds())
  {
    if (got_it[b] == false) 
    {
      if (atoms.atom(index).symbol() == "S" && atoms.atom(b).symbol() == "S") continue;  // cut disulfide bonds
      as.add_index(b);
      got_it[b] = true;
      add_bonds_to_as(b, as);
    }
  }
}

void coords::AtomtypeFinder::complete_atoms_of_aminoacids(std::vector<AminoAcid>& amino_acids)
{
  for (auto &as : amino_acids) {
    for (auto j = 0u; j < as.get_indices().size(); ++j) {
      auto i = as.get_indices()[j];
      add_bonds_to_as(i, as);
    }
  }
}

void coords::AminoAcid::get_chemical_formula(Atoms const &atoms)
{
  int C{ 0 }, H{ 0 }, N{ 0 }, O{ 0 }, S{ 0 }, other{ 0 };
  for (auto i : indices)
  {
    if (atoms.atom(i).symbol() == "C") C++;
    else if (atoms.atom(i).symbol() == "H") H++;
    else if (atoms.atom(i).symbol() == "N") N++;
    else if (atoms.atom(i).symbol() == "O") O++;
    else if (atoms.atom(i).symbol() == "S") S++;
    else other++;
  }
  chemical_formula = { C, H, N, O, S, other };
}

void coords::AminoAcid::get_name_from_chemical_formula()
{
  if (chemical_formula == std::vector<int>{ 2, 3, 1, 1, 0, 0 }) res_name = residueName::GLY;
  else if (chemical_formula == std::vector<int>{ 3, 5, 1, 1, 0, 0 }) res_name = residueName::ALA;
  else if (chemical_formula == std::vector<int>{ 6, 13, 4, 1, 0, 0 }) res_name = residueName::ARG; // protonated
  else if (chemical_formula == std::vector<int>{ 4, 6, 2, 2, 0, 0 }) res_name = residueName::ASN;
  else if (chemical_formula == std::vector<int>{ 4, 4, 1, 3, 0, 0 }) res_name = residueName::ASP;  // deprotonated
  else if (chemical_formula == std::vector<int>{ 3, 5, 1, 1, 1, 0 }) res_name = residueName::CYS;  // "normal" cysteine
  else if (chemical_formula == std::vector<int>{ 3, 4, 1, 1, 1, 0 }) res_name = residueName::CYM;  // deprotonated cysteine or disulfide bond, will later be distinguished
  else if (chemical_formula == std::vector<int>{ 5, 8, 2, 2, 0, 0 }) res_name = residueName::GLN;
  else if (chemical_formula == std::vector<int>{ 5, 6, 1, 3, 0, 0 }) res_name = residueName::GLU;  // deprotonated
  else if (chemical_formula == std::vector<int>{ 6, 7, 3, 1, 0, 0 }) res_name = residueName::HIS;  // one nitrogen protonated, will later be distinguished between HID and HIE
  else if (chemical_formula == std::vector<int>{ 6, 8, 3, 1, 0, 0 }) res_name = residueName::HIP;  // both nitrogen protonated
  else if (chemical_formula == std::vector<int>{ 6, 11, 1, 1, 0, 0 }) res_name = residueName::ILE; // can also be LEU at this point, will later be changed
  else if (chemical_formula == std::vector<int>{ 6, 13, 2, 1, 0, 0 }) res_name = residueName::LYS; // protonated
  else if (chemical_formula == std::vector<int>{ 5, 9, 1, 1, 1, 0 }) res_name = residueName::MET;
  else if (chemical_formula == std::vector<int>{ 9, 9, 1, 1, 0, 0 }) res_name = residueName::PHE;
  else if (chemical_formula == std::vector<int>{ 5, 7, 1, 1, 0, 0 }) res_name = residueName::PRO;  // no additional proton
  else if (chemical_formula == std::vector<int>{ 3, 5, 1, 2, 0, 0 }) res_name = residueName::SER;
  else if (chemical_formula == std::vector<int>{ 4, 7, 1, 2, 0, 0 }) res_name = residueName::THR;
  else if (chemical_formula == std::vector<int>{ 11, 10, 2, 1, 0, 0 }) res_name = residueName::TRP;
  else if (chemical_formula == std::vector<int>{ 9, 9, 1, 2, 0, 0 }) res_name = residueName::TYR;
  else if (chemical_formula == std::vector<int>{ 5, 9, 1, 1, 0, 0 }) res_name = residueName::VAL;
  else res_name = residueName::XXX;
}

void coords::AminoAcid::determine_aminoacid(Atoms const &atoms)
{
  get_chemical_formula(atoms);
  if (terminal == terminalState::no)
  {
    get_name_from_chemical_formula();
  }
  else if (terminal == terminalState::C)  
  {
    // deprotonated -> one oxygen more
    chemical_formula[3] = chemical_formula[3] - 1;
    get_name_from_chemical_formula();

    if (res_name == residueName::XXX)
    {
      // protonated -> one oxygen and one hydrogen more
      chemical_formula[1] = chemical_formula[1] - 1;
      get_name_from_chemical_formula();
      terminal = terminalState::C_prot;
    }
  }
  else if (terminal == terminalState::N)  
  {
    // not protonated -> one hydrogen more
    chemical_formula[1] = chemical_formula[1] - 1;
    get_name_from_chemical_formula();

    if (res_name == residueName::XXX)
    {
      // protonated -> two hydrogens more
      chemical_formula[1] = chemical_formula[1] - 1;
      get_name_from_chemical_formula();
    }
  }
  if (res_name == residueName::XXX) std::cout << "unknown amino acid: " << (*this) << "\n";
}

void coords::AminoAcid::assign_backbone_atom_types(Atoms &atoms)
{
  if (terminal == terminalState::no)
  {
    atoms.atom(indices[0]).set_energy_type(178);   // carbonyle O
    atoms.atom(indices[1]).set_energy_type(177);   // carbonyle C
    if (res_name == residueName::GLY) atoms.atom(indices[2]).set_energy_type(165);      // C alpha
    else if (res_name == residueName::PRO) atoms.atom(indices[2]).set_energy_type(188);
    else atoms.atom(indices[2]).set_energy_type(166);   
    if (res_name != residueName::PRO) atoms.atom(indices[3]).set_energy_type(180);   // backbone N
    else atoms.atom(indices[3]).set_energy_type(181);

    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.symbol() == "H" && is_in(indices[3], a.bonds())) a.set_energy_type(183);  // amide H
    }
  }
  if (terminal == terminalState::C || terminal == terminalState::C_prot)
  {
    if (terminal == terminalState::C) atoms.atom(indices[0]).set_energy_type(214);   // carbonyle O
    else   // protonated C-terminus
    {
      auto bonding_symbols = get_bonding_symbols(atoms.atom(indices[0]), atoms);
      if (bonding_symbols.size() == 1) atoms.atom(indices[0]).set_energy_type(210);  // O=C
      else if (bonding_symbols.size() == 2)         // OH
      {
        atoms.atom(indices[0]).set_energy_type(211);
        for (auto b : atoms.atom(indices[0]).bonds()) {
          if (atoms.atom(b).symbol() == "H") atoms.atom(b).set_energy_type(212);   // HO
        }
      }
      else std::cout << "Something strange happened in protonated C-terminus\n";
    }

    if (terminal == terminalState::C) atoms.atom(indices[1]).set_energy_type(213);   // carbonyle C
    else atoms.atom(indices[1]).set_energy_type(213); 

    if (res_name == residueName::GLY) atoms.atom(indices[2]).set_energy_type(226);       // C alpha
    else if (res_name == residueName::PRO) atoms.atom(indices[2]).set_energy_type(228);
    else atoms.atom(indices[2]).set_energy_type(225);
    if (res_name != residueName::PRO) atoms.atom(indices[3]).set_energy_type(180);   // backbone N
    else atoms.atom(indices[3]).set_energy_type(181);

    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);

      if (a.symbol() == "H" && is_in(indices[3], a.bonds())) a.set_energy_type(183);        // amide H

      else if (a.symbol() == "O" && is_in(indices[1], a.bonds()))   // second O
      {
        if (terminal == terminalState::C) a.set_energy_type(214);
        else
        {
          auto bonding_symbols = get_bonding_symbols(a, atoms);
          if (bonding_symbols.size() == 1) a.set_energy_type(210);  // O=C
          else if (bonding_symbols.size() == 2)         // OH
          {
            a.set_energy_type(211);
            for (auto b : a.bonds()) {
              if (atoms.atom(b).symbol() == "H") atoms.atom(b).set_energy_type(212);   // HO
            }
          }
          else std::cout << "Something strange happened in protonated C-terminus\n";
        }
      }
    }
  }
  if (terminal == terminalState::N)
  {
    atoms.atom(indices[0]).set_energy_type(178);   // carbonyle O
    atoms.atom(indices[1]).set_energy_type(177);   // carbonyle C
    if (res_name == residueName::GLY) atoms.atom(indices[2]).set_energy_type(235); // C alpha
    else if (res_name == residueName::PRO) atoms.atom(indices[2]).set_energy_type(238);
    else atoms.atom(indices[2]).set_energy_type(236);  

    if (res_name != residueName::PRO) atoms.atom(indices[3]).set_energy_type(230);   // backbone N
    else atoms.atom(indices[3]).set_energy_type(252);

    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.symbol() == "H" && is_in(indices[3], a.bonds()))
      {
        if (res_name != residueName::PRO) a.set_energy_type(233);  // amide H
        else a.set_energy_type(253);
      }
    }
  }
}

void coords::AminoAcid::assign_atom_types(Atoms &atoms)
{
  // for all known amino acids: assign atomtypes for backbone atoms
  if (res_name != residueName::XXX) assign_backbone_atom_types(atoms);    

  // assign atomtypes to aminoacid sidechains

  if (res_name == residueName::GLY || res_name == residueName::ALA || res_name == residueName::VAL || res_name == residueName::ILE)   // also contains LEU
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "H") a.set_energy_type(85);                         // sidechain H's
        else if (a.symbol() == "C")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (count_element("C", bonding_symbols) == 3 && count_element("H", bonding_symbols) == 1){    // CH
            a.set_energy_type(82);
          }
          else if (count_element("C", bonding_symbols) == 2 && count_element("H", bonding_symbols) == 2){  // CH2
            a.set_energy_type(81);
            if (is_in(indices[2], a.bonds())) res_name = residueName::LEU;      // set residue name to LEU if CH2 is bound to C alpha atom
          }
          else if (count_element("C", bonding_symbols) == 1 && count_element("H", bonding_symbols) == 3){  // CH3
            a.set_energy_type(80);
          }
          else std::cout << "Something went wrong in residue " << res_name << " with atom " << a.symbol() << ".\nNo atom type assigned.\n";
        }
        else std::cout << "strange atom in residue "<<res_name<<": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::CYS)    // normal Cys, with S-H
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "S") a.set_energy_type(142);
        else if (a.symbol() == "C") a.set_energy_type(148);
        else if (a.symbol() == "H")
        {
          auto bonding_partner = atoms.atom(a.bonds()[0]);
          if (bonding_partner.symbol() == "S") a.set_energy_type(146);
          else a.set_energy_type(85);
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::CYM)    // deprotonated, also includes CYX (disulfide)
  {
    for (auto i = 4u; i < indices.size(); ++i)       // first distinguish between CYM and CYX
    {
      auto &a = atoms.atom(indices[i]);
      if (a.symbol() == "S")
      {
        std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
        if (is_in("S", bonding_symbols)) res_name = residueName::CYX;           // set residue name to CYX if there is a disulfide bond
        break;
      }
    }
    if (res_name == residueName::CYM && Config::get().general.verbosity > 1)
    {
      std::cout << "Warning! Residue " << res_name << " can't be parametrized with OPLSAA. Taken parameters for CYS instead.\n";
    }

    for (auto i = 4u; i < indices.size(); ++i)      // assign atomtypes
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "H") a.set_energy_type(85);
        else if (a.symbol() == "S")
        {
          if (res_name == residueName::CYM) a.set_energy_type(142);        // this is for protonated state
          else if (res_name == residueName::CYX) a.set_energy_type(145);
        }
        else if (a.symbol() == "C")
        {
          if (res_name == residueName::CYM) a.set_energy_type(148);        
          else if (res_name == residueName::CYX) a.set_energy_type(156);
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::HIP)
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "N") a.set_energy_type(453);
        else if (a.symbol() == "H")
        {
          auto bonding_partner = atoms.atom(a.bonds()[0]);
          if (bonding_partner.symbol() == "N") a.set_energy_type(454);
          else if (bonding_partner.symbol() == "C")
          {
            if (bonding_partner.bonds().size() == 4) a.set_energy_type(85);
            else if (bonding_partner.bonds().size() == 3) a.set_energy_type(91);
            else std::cout << "Something went wrong with element " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
          }
          else std::cout << "Wrong bonding partner for " << a.symbol() << " in residue " << res_name << "\nNo atom type assigned.\n";
        }
        else if (a.symbol() == "C")
        {
          if (a.bonds().size() == 4) a.set_energy_type(446);
          else if (a.bonds().size() == 3)
          {
            std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
            if (count_element("N", bonding_symbols) == 2) a.set_energy_type(450);
            else if (count_element("N", bonding_symbols) == 1) a.set_energy_type(451);
            else std::cout << "Something went wrong with element " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
          }
          else std::cout << "Something went wrong with element " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::HIS)
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "N")
        {
          if (a.bonds().size() == 2) a.set_energy_type(452);
          else if (a.bonds().size() == 3) a.set_energy_type(444);
          else std::cout << "Wrong number of bonding partners for element " << a.symbol() << " in residue " << res_name << "\n";
        }
        else if (a.symbol() == "H")
        {
          auto bonding_partner = atoms.atom(a.bonds()[0]);
          if (bonding_partner.symbol() == "N") a.set_energy_type(445);
          else if (bonding_partner.symbol() == "C")
          {
            if (bonding_partner.bonds().size() == 4) a.set_energy_type(85);
            else if (bonding_partner.bonds().size() == 3) a.set_energy_type(91);
            else std::cout << "Something went wrong with element " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
          }
          else std::cout << "Wrong bonding partner for " << a.symbol() << " in residue " << res_name << "\nNo atom type assigned.\n";
        }
        else if (a.symbol() == "C")
        {
          if (a.bonds().size() == 4) a.set_energy_type(446);
          else if (a.bonds().size() == 3)
          {
            std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
            if (count_element("N", bonding_symbols) == 2) a.set_energy_type(447);
            else if (count_element("N", bonding_symbols) == 1 && count_element("C", bonding_symbols) == 2)   // C_gamma
            {
              int index_of_N = a.bonds()[find_index("N", bonding_symbols)];
              auto N_atom = atoms.atom(index_of_N);
              if (N_atom.bonds().size() == 2) {
                a.set_energy_type(448);
                res_name = residueName::HIE;
              }
              else if (N_atom.bonds().size() == 3) {
                a.set_energy_type(449);
                res_name = residueName::HID;
              }
              else std::cout << "Something went wrong with element " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
            }
            else if (count_element("N", bonding_symbols) == 1 && count_element("C", bonding_symbols) == 1 && count_element("H", bonding_symbols) == 1)  // C_delta
            {
              int index_of_N = a.bonds()[find_index("N", bonding_symbols)];
              auto N_atom = atoms.atom(index_of_N);
              if (N_atom.bonds().size() == 2) {
                a.set_energy_type(448);
                res_name = residueName::HID;
              }
              else if (N_atom.bonds().size() == 3) {
                a.set_energy_type(449);
                res_name = residueName::HIE;
              }
              else std::cout << "Something went wrong with element " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
            }
            else std::cout << "Something went wrong with element " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
          }
          else std::cout << "Wrong number of bonding partners for element " << a.symbol() << " in residue " << res_name << "\n";
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::PRO) 
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "H") a.set_energy_type(85);                        
        else if (a.symbol() == "C")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (is_in("N", bonding_symbols))                   // C_delta
          {
            if (terminal == terminalState::N) a.set_energy_type(239);
            else a.set_energy_type(187);
          }
          else a.set_energy_type(81);                        // C_beta and C_gamma
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::ASP || res_name == residueName::GLU)
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "H") a.set_energy_type(85);                        
        else if (a.symbol() == "C")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (is_in("O", bonding_symbols)) a.set_energy_type(213);
          else
          {
            for (auto b : a.bonds()) {
              if (atoms.atom(b).symbol() == "C")
              {
                std::vector<std::string> bbonds = get_bonding_symbols(atoms.atom(b), atoms);
                if (is_in("O", bbonds)) {
                  a.set_energy_type(216);     // CH2 next to COO
                  break;
                }
              }
              a.set_energy_type(81);          // CH2 in Glu, next to C_alpha
            } 
          }
        }
        else if (a.symbol() == "O") a.set_energy_type(214);
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::ASN || res_name == residueName::GLN)
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "O") a.set_energy_type(178);
        else if (a.symbol() == "N") a.set_energy_type(179);
        else if (a.symbol() == "C")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (is_in("N", bonding_symbols)) a.set_energy_type(177);
          else a.set_energy_type(81);
        }
        else if (a.symbol() == "H")
        {
          auto bonding_partner = atoms.atom(a.bonds()[0]);
          if (bonding_partner.symbol() == "N") a.set_energy_type(182);
          else if (bonding_partner.symbol() == "C") a.set_energy_type(85);
          else std::cout << "Wrong binding partner for " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::SER || res_name == residueName::THR)
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "O") a.set_energy_type(96);
        else if (a.symbol() == "H")
        {
          auto bonding_partner = atoms.atom(a.bonds()[0]);
          if (bonding_partner.symbol() == "O") a.set_energy_type(97);
          else if (bonding_partner.symbol() == "C") a.set_energy_type(85);
          else std::cout << "Something went wrong in residue " << res_name << " with atom " << a.symbol() << ".\nNo atom type assigned.\n";
        }
        else if (a.symbol() == "C")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (count_element("C", bonding_symbols) == 1 && count_element("H", bonding_symbols) == 3) a.set_energy_type(80);
          else if (count_element("C", bonding_symbols) == 1 && count_element("H", bonding_symbols) == 2 && count_element("O", bonding_symbols) == 1) a.set_energy_type(99);  // Ser
          else if (count_element("C", bonding_symbols) == 2 && count_element("H", bonding_symbols) == 1 && count_element("O", bonding_symbols) == 1) a.set_energy_type(100); // Thr
          else std::cout << "Wrong binding partners for " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::LYS)
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "N") a.set_energy_type(230);
        else if (a.symbol() == "C")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (is_in("N", bonding_symbols)) a.set_energy_type(236);
          else a.set_energy_type(81);
        }
        else if (a.symbol() == "H")
        {
          auto bonding_partner = atoms.atom(a.bonds()[0]);
          if (bonding_partner.symbol() == "N") a.set_energy_type(233);
          else if (bonding_partner.symbol() == "C") a.set_energy_type(85);
          else std::cout << "Wrong binding partner for " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::MET)
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "S") a.set_energy_type(144);
        else if (a.symbol() == "H") a.set_energy_type(85);
        else if (a.symbol() == "C")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (count_element("C", bonding_symbols) == 2 && count_element("H", bonding_symbols) == 2) a.set_energy_type(81);
          else if (count_element("S", bonding_symbols) == 1 && count_element("H", bonding_symbols) == 3) a.set_energy_type(151);
          else if (count_element("S", bonding_symbols) == 1 && count_element("H", bonding_symbols) == 2 && count_element("C", bonding_symbols) == 1) a.set_energy_type(152);
          else std::cout << "Wrong binding partners for " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::ARG)
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "N")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (count_element("C", bonding_symbols) == 2 && count_element("H", bonding_symbols) == 1) a.set_energy_type(246);
          else if (count_element("C", bonding_symbols) == 1 && count_element("H", bonding_symbols) == 2) a.set_energy_type(243);
          else std::cout << "Wrong binding partners for " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
        }
        else if (a.symbol() == "H")
        {
          auto bonding_partner = atoms.atom(a.bonds()[0]);
          if (bonding_partner.symbol() == "C") a.set_energy_type(85);
          else if (bonding_partner.symbol() == "N")
          {
            std::vector<std::string> bbonds = get_bonding_symbols(atoms.atom(a.bonds()[0]), atoms);
            if (count_element("C", bbonds) == 2 && count_element("H", bbonds) == 1) a.set_energy_type(247);
            else if (count_element("C", bbonds) == 1 && count_element("H", bbonds) == 2) a.set_energy_type(244);
            else std::cout << "Something went wrong with " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
          }
          else std::cout << "Wrong binding partners for " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
        }
        else if (a.symbol() == "C")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (bonding_symbols.size() == 3) a.set_energy_type(245);
          else if (count_element("C", bonding_symbols) == 1 && count_element("H", bonding_symbols) == 2 && count_element("N", bonding_symbols) == 1) a.set_energy_type(250); // C_delta
          else
          {
            for (auto b : a.bonds()) {
              if (atoms.atom(b).symbol() == "C")
              {
                std::vector<std::string> bbonds = get_bonding_symbols(atoms.atom(b), atoms);
                if (count_element("C", bbonds) == 2 && count_element("H", bbonds) == 2) {}   // both remaining atoms have neighbouring CH2 groups so they are ignored
                else if (count_element("C", bbonds) == 2 && count_element("H", bbonds) == 1 && count_element("N", bbonds) == 1) a.set_energy_type(81);   // C_beta
                else if (count_element("C", bbonds) == 1 && count_element("H", bbonds) == 2 && count_element("N", bbonds) == 1) a.set_energy_type(251);  // C_gamma 
                else std::cout << "Something went wrong with " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
              }
            }
          }
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }
  
  else if (res_name == residueName::PHE || res_name == residueName::TYR)
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "O") a.set_energy_type(109);
        else if (a.symbol() == "C")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (count_element("C", bonding_symbols) == 2 && count_element("H", bonding_symbols) == 2) a.set_energy_type(94);
          else if (count_element("C", bonding_symbols) == 2 && count_element("O", bonding_symbols) == 1) a.set_energy_type(108);
          else if ((count_element("C", bonding_symbols) == 2 && count_element("H", bonding_symbols) == 1) || count_element("C", bonding_symbols) == 3) a.set_energy_type(90);
          else std::cout << "Wrong binding partners for " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
        }
        else if (a.symbol() == "H")
        {
          auto bonding_partner = atoms.atom(a.bonds()[0]);
          if (bonding_partner.symbol() == "O") a.set_energy_type(110);
          else if (bonding_partner.symbol() == "C")
          {
            if (bonding_partner.bonds().size() == 3) a.set_energy_type(91);
            else a.set_energy_type(85);
          }
          else std::cout << "Wrong binding partner for " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }

  else if (res_name == residueName::TRP)
  {
    for (auto i = 4u; i < indices.size(); ++i)
    {
      auto &a = atoms.atom(indices[i]);
      if (a.energy_type() == 0)          // not yet assigned as backbone atom
      {
        if (a.symbol() == "N") a.set_energy_type(444);
        else if (a.symbol() == "H")
        {
          auto bonding_partner = atoms.atom(a.bonds()[0]);
          if (bonding_partner.symbol() == "N") a.set_energy_type(445);
          else if (bonding_partner.symbol() == "C")
          {
            if (bonding_partner.bonds().size() == 3) a.set_energy_type(91);
            else a.set_energy_type(85);
          }
          else std::cout << "Wrong bonding partner for " << a.symbol() << " in residue " << res_name << "\nNo atom type assigned.\n";
        }
        else if (a.symbol() == "C")
        {
          std::vector<std::string> bonding_symbols = get_bonding_symbols(a, atoms);
          if (bonding_symbols.size() == 4) a.set_energy_type(81);
          else   // 3 bonding partners
          {
            if (count_element("C", bonding_symbols) == 2 && count_element("N", bonding_symbols) == 1) a.set_energy_type(443);
            else if (count_element("C", bonding_symbols) == 2 && count_element("H", bonding_symbols) == 1) a.set_energy_type(90);
            else if (count_element("C", bonding_symbols) == 1 && count_element("H", bonding_symbols) == 1 && count_element("N", bonding_symbols) == 1) a.set_energy_type(455);
            else 
            {
              for (auto b : a.bonds()) {
                std::vector<std::string> bbonds = get_bonding_symbols(atoms.atom(b), atoms);
                if (count_element("C", bbonds) == 3) {}
                else if (count_element("C", bbonds) == 2 && count_element("H", bbonds) == 2) a.set_energy_type(441);
                else if (count_element("C", bbonds) == 2 && count_element("H", bbonds) == 1) a.set_energy_type(442);
                else if (count_element("C", bbonds) == 2 && count_element("N", bbonds) == 1) a.set_energy_type(442);
                else if (count_element("C", bbonds) == 1 && count_element("N", bbonds) == 1 && count_element("H", bbonds) == 1) a.set_energy_type(441);
                else std::cout << "Something went wrong with element " << a.symbol() << " in residue " << res_name << ".\nNo atom type assigned.\n";
              }
            }
          }
        }
        else std::cout << "strange atom in residue " << res_name << ": " << a.symbol() << "\n";
      }
    }
  }
}

void coords::AminoAcid::correct_residue_names(Atoms& atoms)       
{
  if (res_name == residueName::CYM) res_name = residueName::CYS;
  else if (res_name == residueName::HIP) res_name = residueName::HIS;

  else if (res_name == residueName::ILE)
  {
    auto c_alpha = atoms.atom(indices[2]);
    auto bonding_partners = c_alpha.bonds();  // bonding partners of C_alpha
    for (auto b : bonding_partners)
    {
      auto bond = atoms.atom(b);
      auto bonding_symbols = get_bonding_symbols(bond, atoms);
      if (count_element("H", bonding_symbols) == 2 && count_element("C", bonding_symbols) == 2) {
        res_name = residueName::LEU;    // if CH2 group next to C_alpha
      }
    }
  }
}

void coords::AtomtypeFinder::find_energy_types()
{
  get_some_easy_atomtypes();
  auto amino_acids = get_aminoacids();
  for (auto &as : amino_acids) as.determine_aminoacid(atoms);
  for (auto &as : amino_acids) as.assign_atom_types(atoms);

  if (Config::get().general.verbosity > 3)   // print aminoacid sequence
  {
    std::cout << "Amino acids:\n";
    for (auto as : amino_acids) std::cout << as << " ";
    std::cout << "\n";
  }

  bool first_not_assigned = true;          // print those atoms for which no atomtypes are assigned
  std::string indexstring{ "" };
  for (auto i = 0u; i < atoms.size(); ++i) 
  {
    if (atoms.atom(i).energy_type() == 0) 
    {
      if (first_not_assigned == true) {
        std::cout << "-------------------------------------------------\n";
        std::cout << "No atomtypes found for (to be copied into pymol):\nselect sele, id ";
        first_not_assigned = false;
      }
      indexstring += (std::to_string(i+1) + "+");
    }
  }
  indexstring.pop_back();   // remove '+'
  std::cout << indexstring << "\n";
}