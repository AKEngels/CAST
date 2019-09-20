/**
CAST 3
find_as.h
Purpose: header for finding amino acids and giving them atomtypes

@author Susanne Sauer
@version 1.0
*/
#pragma once

#include "coords.h"

/**main function for task FIND_AS
@param coords: Coordinates object
@param filename: name of the outputfile*/
void find_as(coords::Coordinates const & coords, std::string const& filename);

/**terminal states: not terminal, C-terminal (protonated or not), N-terminal*/
enum class terminalState { no, C, C_prot, N };
/**overloaded output operator for terminalState*/
inline std::ostream& operator<< (std::ostream& os, const terminalState& T)
{
  switch (T)
  {
  case terminalState::no: os << "not terminal"; break;
  case terminalState::C:  os << "C-terminal"; break;
  case terminalState::C_prot:  os << "C-terminal"; break;
  case terminalState::N:  os << "N-terminal"; break;
  }
  return os;
}

/**known aminoacids: 3-letter codes + some special names inspired by AMBER (http://ambermd.org/tutorials/advanced/tutorial1_orig/section1.htm)
CYX: cysteine in disulfide bridge
CYM: deprotonated cysteine
HID: histidine protonated at N_delta
HIE: histidine protonated at N_epsilon
HIP: histidine where both nitrogen atoms are protonated
XXX: just a wildcard for not known aminoacid*/
enum class residueName { ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, CYX, CYM, HID, HIE, HIP, XXX };
/**function to convert residueName to string*/
inline std::string res_to_string(const residueName& res)
{
  switch (res)
  {
  case residueName::ALA: return "ALA";
  case residueName::ARG: return "ARG";
  case residueName::ASN: return "ASN";
  case residueName::ASP: return "ASP";
  case residueName::CYS: return "CYS";
  case residueName::GLN: return "GLN";
  case residueName::GLU: return "GLU";
  case residueName::GLY: return "GLY";
  case residueName::HIS: return "HIS";
  case residueName::ILE: return "ILE";
  case residueName::LEU: return "LEU";
  case residueName::LYS: return "LYS";
  case residueName::MET: return "MET";
  case residueName::PHE: return "PHE";
  case residueName::PRO: return "PRO";
  case residueName::SER: return "SER";
  case residueName::THR: return "THR";
  case residueName::TRP: return "TRP";
  case residueName::TYR: return "TYR";
  case residueName::VAL: return "VAL";
  case residueName::CYX: return "CYX";
  case residueName::CYM: return "CYM";
  case residueName::HID: return "HID";
  case residueName::HIE: return "HIE";
  case residueName::HIP: return "HIP";
  case residueName::XXX: return "XXX";
  default: return "XXX";
  }
}
/**overloaded output operator for residueName*/
inline std::ostream& operator<< (std::ostream& os, const residueName& res)
{
  os << res_to_string(res);
  return os;
}

/**class for one amino acid*/
class AminoAcid
{
public:
  /**constructor
  @param i: indices of backbone atoms (order: carbonyle O, carbonyle C, C alpha, N)
  @param T: terminal state*/
  AminoAcid(std::vector<std::size_t> i, terminalState T) : indices(i), terminal(T) {};

  /**get all indices*/
  std::vector<std::size_t> get_indices() const { return indices; }
  /**add an index to indices
  @param i: index to be added*/
  void add_index(std::size_t const i) { indices.emplace_back(i); }

  /**determine residueName of aminoacid and saving it into res_name
  as this is only done by chemical formula there might be inaccuracies, i.e. for protonation states of HIS or the binding mode of CYS
  those will be corrected later when assigning atomtypes
  @param atoms: atom vector*/
  void determine_aminoacid(coords::Atoms const& atoms);
  /**assigns oplsaa atomtypes to atoms
  @param atoms: atom vector*/
  void assign_atom_types(coords::Atoms& atoms);
  /**function that assigns the correct three-letter codes to those aminoacids where determine_aminoacid() didn't
  i.e. HIP -> HIS, CYM -> CYS, ILE -> LEU (sometimes)
  @param atoms: atom vector*/
  void correct_residue_names(coords::Atoms& atoms);
  /**returns residue name as string*/
  std::string get_res_name() { return res_to_string(res_name); }

private:
  /**indices of all atoms belonging to amino acid
  first 4 indices are those of carbonyle O, carbonyle C, C alpha, amide N*/
  std::vector<std::size_t> indices;
  /**terminal state of amino acid*/
  terminalState terminal;
  /**residue name*/
  residueName res_name{ residueName::XXX };
  /**chemical formula: number of C, H, N, O, S, other in this order*/
  std::vector<int> chemical_formula;

  /**get chemical formula of aminoacid and save it into chemical_formula
  @param atoms: atom vector*/
  void get_chemical_formula(coords::Atoms const& atoms);
  /**assigns oplsaa atomtypes to backbone atoms (part of assign_atom_types())
  @param atoms: atom vector*/
  void assign_backbone_atom_types(coords::Atoms& atoms);
  /**determines residue name from chemical formula*/
  void get_name_from_chemical_formula();

  /**overloaded output operator for AminoAcid*/
  friend std::ostream& operator<< (std::ostream& os, const AminoAcid& as);
};

/**overloaded output operator for AminoAcid*/
inline std::ostream& operator<< (std::ostream& os, const AminoAcid& as)
{
  os << as.res_name;
  if (as.terminal != terminalState::no) os << "(" << as.terminal << ")";
  return os;
}

/**struct to get forcefield energy type from amino acids*/
class AtomtypeFinder
{

public:
  /**constructor
  sets size of got_it to number of atoms and sets all of them to false*/
  AtomtypeFinder(coords::Atoms& a) : atoms(a)
  {
    got_it.resize(atoms.size());
    for (auto&& g : got_it) g = false;
  };

  /**function that finds all possible atomtypes*/
  void find_energy_types();

  /**function that creates amino acids with backbone atoms and terminal state*/
  std::vector<AminoAcid> get_aminoacids();

  /**function to determine which atoms have already been recognized (normally if they are part of an amino acid)*/
  bool recognized_atom(std::size_t const a) const { return got_it[a]; };

private:
  /**reference to atoms
  will be changed inside this class (addition of atomtypes)*/
  coords::Atoms& atoms;

  /**vector that tells us if an atom either has already a forcefield type or is in an aminoacid*/
  std::vector<bool> got_it;

  /**function that finds atomtypes of some atoms that are quite easy to determine
  sets their value for got_it to true
  at the moment the atomtypes of Na ions and water molecules are found*/
  void get_some_easy_atomtypes();
  /**function that fills the rest of the atoms into the aminoacids*/
  void complete_atoms_of_aminoacids(std::vector<AminoAcid>& amino_acids);
  /**helperfunction for complete_atoms_of_aminoacids()
  is called recursively on every atom and adds all atoms that are bound to current atom to amino acid
  stops at disulfide bonds*/
  void add_bonds_to_as(int index, AminoAcid& as);
};