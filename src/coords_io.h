#pragma once 

#if defined _OPENMP
#include <omp.h>
#endif
#include <memory>
#include <string>
#include <iostream>
#include <sstream>
#include "scon_utility.h"
#include "coords.h"
#include "configuration.h"
#include "helperfunctions.h"
#pragma once

namespace coords
{
  namespace input
  {
    // format types
    struct types { enum T { TINKER, PDB, AMBER, XYZ }; };

    /**converts string into object of enum input::types::T*/
    inline types::T get_type(std::string const &str)
    {
      if (str.find("PDB") != str.npos) return types::PDB;
      else if (str.find("AMBER") != str.npos) return types::AMBER;
      else if (str.find("XYZ") != str.npos) return types::XYZ;
      else return types::TINKER;
    }
    
    /**general format class
    input formats like AMBER, TINKER, PDB or XYZ inherit from this*/
    class format
    {
    protected:
      format(void) { }
      coords::Ensemble_PES input_ensemble;
    public:
      virtual ~format(void) { }
      // Number of Atoms (N) and Structures (M)
      std::size_t atoms(void) const { return input_ensemble.size() > 0 ? input_ensemble.back().size() : 0U; }
      std::size_t size(void) const { return input_ensemble.size(); }
      /** Read Structure and return coordinates (this function must be overwritten for newly implemented inputformats)*/
      virtual Coordinates read(std::string) = 0;
      // Get structure i
      PES_Point structure(std::size_t const i = 0u) const { return input_ensemble[i]; }
      Ensemble_PES const & PES(void) const { return input_ensemble; }
      Ensemble_PES & PES(void) { return input_ensemble; }
      // Iterators
      Ensemble_PES::iterator begin(void) { return input_ensemble.begin(); }
      Ensemble_PES::iterator end(void) { return input_ensemble.end(); }
      Ensemble_PES::const_iterator begin(void) const { return input_ensemble.begin(); }
      Ensemble_PES::const_iterator end(void) const { return input_ensemble.end(); }
    };

    // new format creator
    input::format* new_format(void);
    //new format creator for interface creation
    input::format* additional_format(void);
    input::format* new_interf_format(void);

    namespace formats
    {
      // stuff for AMBER input

      /**number of sections in AMBER input file*/
      static const unsigned int sections_size = 91u;
      /**names of sections in AMBER input file*/
      static std::string const amber_sections[91] =
      { "TITLE", "POINTERS", "ATOM_NAME", "CHARGE", "ATOMIC_NUMBER",
        "MASS", "ATOM_TYPE_INDEX", "NUMBER_EXCLUDED_ATOMS", "NONBONDED_PARM_INDEX", "POLARIZABILITY",
        "RESIDUE_LABEL", "RESIDUE_POINTER", "BOND_FORCE_CONSTANT", "BOND_EQUIL_VALUE", "ANGLE_FORCE_CONSTANT",
        "ANGLE_EQUIL_VALUE", "DIHEDRAL_FORCE_CONSTANT", "DIHEDRAL_PERIODICITY", "DIHEDRAL_PHASE", "SCEE_SCALE_FACTOR",
        "SCNB_SCALE_FACTOR", "SOLTY", "LENNARD_JONES_ACOEF", "LENNARD_JONES_BCOEF", "BONDS_INC_HYDROGEN",
        "BONDS_WITHOUT_HYDROGEN", "ANGLES_INC_HYDROGEN", "ANGLES_WITHOUT_HYDROGEN", "DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN",
        "EXCLUDED_ATOMS_LIST", "HBOND_ACOEF", "HBOND_BCOEF", "HBCUT", "AMBER_ATOM_TYPE",
        "TREE_CHAIN_CLASSIFICATION", "JOIN_ARRAY", "IROTAT", /*"SOLVENT_POINTERS",*/ "ATOMS_PER_MOLECULE", // SOLVENT_POINTERS for now disable because it causes trouble with std::find and POINTERS
        "BOX_DIMENSIONS", "CAP_INFO", "CAP_INFO2", "RADIUS_SET", "RADII",
        "IPOL" };

      /**class for reading AMBER input*/
      class amber : public coords::input::format
      {
      public:
        /**function that reads in AMBER input*/
        Coordinates read(std::string);
      private:
        //struct sections;
        static const unsigned int sections_size = 91u;

        enum section_types
        {
          ILLEGAL = -1,
          TITLE, POINTERS, ATOM_NAME, CHARGE, ATOMIC_NUMBER,
          MASS, ATOM_TYPE_INDEX, NUMBER_EXCLUDED_ATOMS, NONBONDED_PARM_INDEX, POLARIZABILITY,
          RESIDUE_LABEL, RESIDUE_POINTER, BOND_FORCE_CONSTANT, BOND_EQUIL_VALUE, ANGLE_FORCE_CONSTANT,
          ANGLE_EQUIL_VALUE, DIHEDRAL_FORCE_CONSTANT, DIHEDRAL_PERIODICITY, DIHEDRAL_PHASE, SCEE_SCALE_FACTOR,
          SCNB_SCALE_FACTOR, SOLTY, LENNARD_JONES_ACOEF, LENNARD_JONES_BCOEF, BONDS_INC_HYDROGEN,
          BONDS_WITHOUT_HYDROGEN, ANGLES_INC_HYDROGEN, ANGLES_WITHOUT_HYDROGEN, DIHEDRALS_INC_HYDROGEN, DIHEDRALS_WITHOUT_HYDROGEN,
          EXCLUDED_ATOMS_LIST, HBOND_ACOEF, HBOND_BCOEF, HBCUT, AMBER_ATOM_TYPE,
          TREE_CHAIN_CLASSIFICATION, JOIN_ARRAY, IROTAT, SOLVENT_POINTERS, ATOMS_PER_MOLECULE,
          BOX_DIMENSIONS, CAP_INFO, CAP_INFO2, RADIUS_SET, RADII,
          IPOL
        };

        //Important Stuff
        std::vector<unsigned int> sectionsPresent;
        std::vector<unsigned int> bondsWithHydrogen, bondsWithoutHydrogen; // [0] binds to [1], [2] to [3], [4] to [5] and so on....
        Atoms atoms;
        Cartesian_Point position;
        std::string name;
        unsigned int numberOfAtoms;
        int pointers_raw[32];

      };

			/**input format XYZ*/
      class xyz : public coords::input::format
      {
      public:
				/**read from XYZ file*/
        Coordinates read(std::string);
      private:
				/**atoms*/
        Atoms atoms;
				/**positions*/
        Cartesian_Point position;

				// STUFF TO GET FORCEFIELD ENERGY TYPES

				/**terminal states: not terminal, C-terminal, N-terminal*/
				enum class terminalState { no, C, N };
        /**overloaded output operator for terminalState*/
        friend std::ostream & operator<< (std::ostream &os, const terminalState &T)
        {
          switch (T)
          {
            case terminalState::no: os << "not terminal"; break;
            case terminalState::C:  os << "C-terminal"; break;
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
        /**overloaded output operator for residueName*/
        friend std::ostream & operator<< (std::ostream &os, const residueName &res)
        {
          switch (res)
          {
            case residueName::ALA: os << "ALA"; break;
            case residueName::ARG: os << "ARG"; break;
            case residueName::ASN: os << "ASN"; break;
            case residueName::ASP: os << "ASP"; break;
            case residueName::CYS: os << "CYS"; break;
            case residueName::GLN: os << "GLN"; break;
            case residueName::GLU: os << "GLU"; break;
            case residueName::GLY: os << "GLY"; break;
            case residueName::HIS: os << "HIS"; break;
            case residueName::ILE: os << "ILE"; break;
            case residueName::LEU: os << "LEU"; break;
            case residueName::LYS: os << "LYS"; break;
            case residueName::MET: os << "MET"; break;
            case residueName::PHE: os << "PHE"; break;
            case residueName::PRO: os << "PRO"; break;
            case residueName::SER: os << "SER"; break;
            case residueName::THR: os << "THR"; break;
            case residueName::TRP: os << "TRP"; break;
            case residueName::TYR: os << "TYR"; break;
            case residueName::VAL: os << "VAL"; break;
            case residueName::CYX: os << "CYX"; break;
            case residueName::CYM: os << "CYM"; break;
            case residueName::HID: os << "HID"; break;
            case residueName::HIE: os << "HIE"; break;
            case residueName::HIP: os << "HIP"; break;
            case residueName::XXX: os << "XXX"; break;
            default: os << "unknown";
          }
          return os;
        }

				/**structure for one amino acid*/
				struct AminoAcid
				{
					/**constructor
					@param i: indices of backbone atoms (order: carbonyle O, carbonyle C, C alpha, N)
					@param T: terminal state*/
					AminoAcid(std::vector<std::size_t> i, terminalState T) : indices(i), terminal(T) {};

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
          void get_chemical_formula(Atoms const &atoms);
          /**determine residueName of aminoacid and saving it into res_name
          as this is only done by chemical formula there might be inaccuracies, i.e. for protonation states of HIS or the binding mode of CYS
          those will be corrected later when assigning atomtypes
          @param atoms: atom vector*/
          void determine_aminoacid(Atoms const &atoms);
          /**assigns oplsaa atomtypes to backbone atoms
          @param atoms: atom vector*/
          void assign_backbone_atom_types(Atoms &atoms);
          /**assigns oplsaa atomtypes to sidechain atoms
          @param atoms: atom vector*/
          void assign_atom_types(Atoms &atoms); 
          /**determines residue name from chemical formula*/
          void get_name_from_chemical_formula();
				};

        /**overloaded output operator for AminoAcid*/
        friend std::ostream & operator<< (std::ostream &os, const AminoAcid &as)
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
          AtomtypeFinder(Atoms &a) : atoms(a)
					{ 
						got_it.resize(atoms.size());
						for (auto &&g : got_it) g = false;
					};

          /**function that finds all possible atomtypes*/
          void find_energy_types();

        private:
					/**reference to atoms 
          will be changed inside this class (addition of atomtypes)*/
					Atoms &atoms;

					/**vector that tells us if an atom either has already a forcefield type or is in an aminoacid*/
					std::vector<bool> got_it;

          /**function that finds atomtypes of some atoms that are quite easy to determine
          sets their value for got_it to true
          at the moment the atomtypes of Na ions and water molecules are found*/
          void get_some_easy_atomtypes();
					/**function that creates amino acids with backbone atoms and terminal state*/
					std::vector<AminoAcid> get_aminoacids();
					/**function that fills the rest of the atoms into the aminoacids*/
					void complete_atoms_of_aminoacids(std::vector<AminoAcid> &amino_acids);
          /**helperfunction for complete_atoms_of_aminoacids()
          is called recursively on every atom and adds all atoms that are bound to current atom to amino acid
          stops at disulfide bonds*/
          void add_bonds_to_as(int index, AminoAcid &as);
				};

      };
      
      /**class for reading PDB as input*/
      class pdb : public coords::input::format
      {
      public:
        /**reads input*/
        Coordinates read(std::string);
      private:
        /**atoms (are filled during reading)*/
        Atoms atoms;
        /**positions (are filled during reading)*/
        Cartesian_Point position;
      };

      /*! Class to read from TINKER coordinate file (.arc)
       *
       * This class is used to read coordinates from
       * a tinker xyz file (.arc or sometimes .xyz).
       * The coordinates are then put into a
       * coords::Coordinates object.
       */
      class tinker : public coords::input::format
      {
      public:
        Coordinates read(std::string);
      private:
        struct line
        {
          Cartesian_Point position;
          std::vector<std::size_t> bonds;
          std::size_t index, tinker_type;
          std::string symbol;
          line(void)
            : bonds(7U, 0U), index(0U),
            tinker_type(0U) {}
          friend std::istream & operator>> (std::istream &is, line &l)
          {
            is >> l.index >> l.symbol >> l.position.x() >> l.position.y() >> l.position.z() >> l.tinker_type
              >> l.bonds[0u] >> l.bonds[1u] >> l.bonds[2u] >> l.bonds[3u] >> l.bonds[4u] >> l.bonds[5u] >> l.bonds[6u];
            return is;
          }
        };

      };
    }

  }


  /**namespace for structure output in different formats*/
  namespace output
  {
    /**motherclass of all output formats*/
    class format
    {
    protected:
      Coordinates const & ref;
      format& operator= (format const &);
      format(Coordinates const &coord_obj) : ref(coord_obj) {}
    public:
      /**this is the function where you have to look if you want to know how the format will look like
      has to be overwritten for all newly implemented formats*/
      virtual void to_stream(std::ostream&) const = 0;
    };

    inline std::ostream& operator<< (std::ostream &stream, format const & out_format)
    {
      out_format.to_stream(stream);
      return stream;
    }

    namespace formats
    {
      /**tinker format*/
      class tinker : public output::format
      {
        tinker& operator= (tinker const &);
      public:
        tinker(Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream(std::ostream&) const;
        static std::string filename(std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".arc")).get_unique_path();
        }
      };

      /**xyz format*/
      class xyz
        : public output::format
      {
        xyz& operator= (xyz const &);
      public:
        xyz(Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream(std::ostream&) const;
        static std::string filename(std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".xyz")).get_unique_path();
        }
      };

      class xyz_dftb
        : public output::format
      {
        xyz& operator= (xyz_dftb const &);
      public:
        xyz_dftb(Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream(std::ostream&) const;
        static std::string filename(std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".xyz")).get_unique_path();
        }
      };

      class xyz_gen
        : public output::format
      {
        xyz& operator= (xyz_gen const &);
      public:
        xyz_gen(Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream(std::ostream&) const;
        static std::string filename(std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".xyz")).get_unique_path();
        }
      };

      class xyz_mopac
        : public output::format
      {
        xyz_mopac& operator= (xyz_mopac const &);
      public:
        xyz_mopac(Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream(std::ostream&) const;
        static std::string filename(std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".xyz")).get_unique_path();
        }
      };

      class xyz_mopac7
        : public output::format
      {
        xyz_mopac7 operator= (xyz_mopac7 const &);
      public:
        xyz_mopac7(Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream(std::ostream&) const;

        static std::string filename(std::string postfix)
        {
          return std::string(Config::get().general.outputFilename).append(postfix).append(".xyz");
        }
      };
      class moldenxyz
        : public output::format
      {
        xyz& operator= (xyz const &);
      public:
        moldenxyz(Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream(std::ostream&) const;
        static std::string filename(std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".xyz")).get_unique_path();
        }
      };

      /**output of z matrix*/
      class zmatrix
        : public output::format
      {
        tinker& operator= (tinker const &);
      public:
        zmatrix(Coordinates const &coord_obj) : output::format(coord_obj) {}
        void to_stream(std::ostream&) const;
        static std::string filename(std::string postfix)
        {
          return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(".zm")).get_unique_path();
        }
      };

    }

    inline std::string filename(std::string postfix = "", std::string extension = "")
    {
      if (extension.length() > 0)
      {
        return scon::StringFilePath(std::string(Config::get().general.outputFilename).append(postfix).append(extension)).get_unique_path();
      }
      if (Config::get().general.output == config::output_types::TINKER) return formats::tinker::filename(postfix);
      else if (Config::get().general.output == config::output_types::XYZ) return formats::xyz::filename(postfix);
      else if (Config::get().general.output == config::output_types::MOLDEN) return formats::moldenxyz::filename(postfix);
      else if (Config::get().general.output == config::output_types::ZMATRIX) return formats::zmatrix::filename(postfix);
      return scon::StringFilePath(Config::get().general.outputFilename).get_unique_path();
    }

  }


}
