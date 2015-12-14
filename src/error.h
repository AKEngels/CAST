#ifndef MOSS_ERR_H

#define MOSS_ERR_H 83
#include <string>

enum MOS_ERR
{
  MOS_SUCCESS=-1,

  MOS_ERR_MEM=0,
    // thrown if allocation error occurs
    MOS_ERR_MEM_BAD_ALLOC,
    // if buffer allocation fails during resize
    MOS_ERR_MEM_RESIZE_ALLOC,
    // thrown if a fixed size array is to contain too many elements
    MOS_ERR_MEM_FIX_OVERFLOW,
    // thrown if an index is to be accessed that is out of range (if security is enabled)
    MOS_ERR_MEM_INDEX_OUT_OF_RANGE,

  // Class configuration errors
  MOS_ERR_CONFIG,
    // thrown if the program fails to open the configuration file
    MOS_ERR_CONFIG_FILE, 
    // thrown if the number of threads is manipulated through command line while compilation was done without openmp
    MOS_ERR_CONFIG_MISSING_OPENMP_SUPPORT,
    MOS_ERR_CONFIG_MISSING_STRUCTURE,
    MOS_ERR_CONFIG_MISSING_PARAMS,
    MOS_ERR_CONFIG_MISSING_intERFACE,
    MOS_ERR_CONFIG_MISSING_TASK,

  // Class coordinates errors
  MOS_ERR_COORD,
    // if structure contains no atoms
    MOS_ERR_NO_ATOMS,
    // thrown if the class fails to obtain parameters for a specific bond
    MOS_ERR_COORD_NOPARAM_BOND, 
    // thrown if the class fails to obtain parameters for a specific bond
    MOS_ERR_COORD_NOPARAM_ANGLE, 
    // Fixation of an Atom that is no in Range
    MOS_ERR_COORD_FIX_OUT_OF_RANGE,
    // thrown if the number of pairs exceed the maximum vector size
    MOS_ERR_COORD_PAIRS_VECTOR_EXCEED, 
    // thrown if the removal of any pair from the pairlist fails (too many pairs)
    MOS_ERR_COORD_PAIRS_OVERFLOW, 
    // thrown if scale pairs overflow occurs
    MOS_ERR_COORD_SCALES_OVERFLOW,
    // thrown if too many pairs are removed
    MOS_ERR_COORD_REMOVE_OVERFLOW,
    // thrown if the linkedcell class returns an element out of scope (> max elements)
    MOS_ERR_COORD_LINKEDCELL_SCOPE,
    // thrown if the coord file contains atoms that have an index which is higher than the number in the first line of the coord file
    MOS_ERR_COORD_ATOM_OVERFLOW,
    // thrown if the coord file contains too few atoms lines (number in the first line of the coord file too high or lines missing)
    MOS_ERR_COORD_ATOM_UNDERFLOW,
    // thrown if the atom index does not match the line number of the current atom (index1!=line1,...)
    MOS_ERR_COORD_ALIGNMENT,
    // thrown if an atom index of 0 or below is found
    MOS_ERR_COORD_INDEXATION,
    // thrown if the number of bonds of an atom does not match the paramteric number of bound atoms
    MOS_ERR_COORD_ATOM_NUMBER_BONDS,
    // the coord file has an insufficient number of lines (<=1) or was not openend successful
    MOS_ERR_COORD_MISSING_LINES,
    // if the index for extraction is out of range
    MOS_ERR_COORD_EXTRACT_INDEX_OOR,
    // thrown if the program was not able to evaluate the internal structure
    MOS_ERR_COORD_internal_INDEXATION,
    // if the internals struct has more or less elements than the main positions vector
    MOS_ERR_COORD_internal_WRONG_COUNT,
	  //! if switchdistance is larger than the cutoff radius
	  MOS_ERR_COORD_CUTOFF,

  //! Class fileManipulation errors
  MOS_ERR_FILE,
    // thrown if the file cannot be opened [from all constructors of fileManipulation in filemanipulation.cc]
    MOS_ERR_FILE_OPEN,
    // thrown if reading a file is not possible due to a stream that is not ready to read [from fileManipulation::readFile() in filemanipulation.cc]
    MOS_ERR_FILE_READ,
    // thrown if a line shall be accessed that is out of scope (> max num lines) [from fileManipulation::getLine() in filemanipulation.cc]
    MOS_ERR_FILE_INVALID_ACCESS,
    // thrown if writing to a file is not possible due to a stream that is not ready to write 
    MOS_ERR_FILE_WRITE,

  //! Energy class errors
  MOS_ERR_ENERGY,
    // The specified interface is not available
    MOS_ERR_ENERGY_UNKNOWN_INTERFACE,
    MOS_ERR_ENERGY_MOPAC_EXE_NOT_FOUND,
    MOS_ERR_ENERGY_MOPAC_DNF,
    MOS_ERR_ENERGY_MOPAC_MISSING_OUTPUT,
    MOS_ERR_ENERGY_MOPAC_OUTPUT_ERROR,
    // TERACHEM ERRORS
    MOS_ERR_ENERGY_TERACHEM_INIT_MPI,
    MOS_ERR_ENERGY_TERACHEM_SEND_MPI,
    MOS_ERR_ENERGY_TERACHEM_RECIEVE_MPI,
    MOS_ERR_ENERGY_TERACHEM_SCF_DNC,

  //! Class Linkedcell errors
  MOS_ERR_LINKEDCELL,
    // thrown if the number of elements is insufficient to use linkedcells
    MOS_ERR_LINKEDCELL_INSUFFICIENT_ELEMENTS, 
    // thrown if the position data is missing for initialization
    MOS_ERR_LINKEDCELL_INIT,
    // thrown if updating the linkedcell information is not possible due to a lack of data or if the initialization has not been done yet
    MOS_ERR_LINKEDCELL_UPDATE,
    MOS_ERR_LINKEDCELL_LINK,
    MOS_ERR_LINKEDCELL_PTR_NULL,
    
  //! Class parameters errors
  MOS_ERR_PARAMETER,
    MOS_ERR_PARAMETER_FILE,
    // thrown if the forcefield is unknown
    MOS_ERR_PARAMETER_UNKNOWN_FF,
    // thrown if a known FF could not be identified
    MOS_ERR_PARAMETER_NO_KNOWN_FF,
    // thrown if no vdw params can be found
    MOS_ERR_PARAMETER_NO_VDW,
    // thrown if no charge params can be found
    MOS_ERR_PARAMETER_NO_CHARGE,

  MOS_ERR_STARTOPT,
    MOS_ERR_STARTOPT_SOLVADD_FAIL,

  MOS_ERR_SOLVADD,
    MOS_ERR_SOLVADD_WATER_PARAM_BOND,
    MOS_ERR_SOLVADD_WATER_PARAM_ANGLE,
    MOS_ERR_SOLVADD_WATER_SET_OVERFLOW,
    MOS_ERR_SOLVADD_WATER_NUM_UNDERFLOW,

  MOS_ERR_RINGSEARCH,

  MOS_ERR_FOLD,
    
  MOS_ERR_VECT,
   MOS_ERR_VECT_OUT_OF_RANGE,

  MOS_ERR_VECT_int,

  MOS_ERR_MD,
    MOS_ERR_MD_OUTPUT_FILE,

  MOS_ERR_INPUT,
    MOS_ERR_INPUT_UNKNOWN_INPUT_TYPE,

  MOS_ERR_RDF,
    // thrown if a certain distance does not fit into the dists vector due to its size
    MOS_ERR_RDF_BOX_OVERFLOW,

  MOS_ERR_DIMER,
    MOS_ERR_DIMER_NUMBER,
    MOS_ERR_DIMER_DISTANCE,
    MOS_ERR_DIMER_MAX_ITERATIONS,

  MOS_ERR_PATH,
    MOS_ERR_PATH_RANGE
};

static const std::string ERR_STR_MAP[MOSS_ERR_H] =
{
  "ERR_MEM",
    // thrown if allocation error occurs
    "ERR_MEM_BAD_ALLOC",
    // if buffer allocation fails during resize
    "ERR_MEM_RESIZE_ALLOC",
    // thrown if a fixed size array is to contain too many elements
    "ERR_MEM_FIX_OVERFLOW",
    // thrown if an index is to be accessed that is out of range (if security is enabled)
    "ERR_MEM_INDEX_OUT_OF_RANGE",

  // Class configuration errors
  "ERR_CONFIG",
    // thrown if the program fails to open the configuration file
    "ERR_CONFIG_FILE", 
    // thrown if the number of threads is manipulated through command line while compilation was done without openmp
    "ERR_CONFIG_MISSING_OPENMP_SUPPORT",
    "ERR_CONFIG_MISSING_STRUCTURE",
    "ERR_CONFIG_MISSING_PARAMS",
    "ERR_CONFIG_MISSING_intERFACE",
    "ERR_CONFIG_MISSING_TASK",

  // Class coordinates errors
  "ERR_COORD",
    // if structure contains no atoms
    "ERR_NO_ATOMS",
    // thrown if the class fails to obtain parameters for a specific bond
    "ERR_COORD_NOPARAM_BOND", 
    // thrown if the class fails to obtain parameters for a specific bond
    "ERR_COORD_NOPARAM_ANGLE", 
    // Fixation of an Atom that is no in Range
    "ERR_COORD_FIX_OUT_OF_RANGE",
    // thrown if the number of pairs exceed the maximum vector size
    "ERR_COORD_PAIRS_VECTOR_EXCEED", 
    // thrown if the removal of any pair from the pairlist fails (too many pairs)
    "ERR_COORD_PAIRS_OVERFLOW", 
    // thrown if scale pairs overflow occurs
    "ERR_COORD_SCALES_OVERFLOW",
    // thrown if too many pairs are removed
    "ERR_COORD_REMOVE_OVERFLOW",
    // thrown if the linkedcell class returns an element out of scope (> max elements)
    "ERR_COORD_LINKEDCELL_SCOPE",
    // thrown if the coord file contains atoms that have an index which is higher than the number in the first line of the coord file
    "ERR_COORD_ATOM_OVERFLOW",
    // thrown if the coord file contains too few atoms lines (number in the first line of the coord file too high or lines missing)
    "ERR_COORD_ATOM_UNDERFLOW",
    // thrown if the atom index does not match the line number of the current atom (index1!=line1,...)
    "ERR_COORD_ALIGNMENT",
    // thrown if an atom index of 0 or below is found
    "ERR_COORD_INDEXATION",
    // thrown if the number of bonds of an atom does not match the paramteric number of bound atoms
    "ERR_COORD_ATOM_NUMBER_BONDS",
    // the coord file has an insufficient number of lines (<=1) or was not openend successful
    "ERR_COORD_MISSING_LINES",
    // if the index for extraction is out of range
    "ERR_COORD_EXTRACT_INDEX_OOR",
    // thrown if the program was not able to evaluate the internal structure
    "ERR_COORD_internal_INDEXATION",
    // if the internals struct has more or less elements than the main positions vector
    "ERR_COORD_internal_WRONG_COUNT",
	  //! if switchdistance is larger than the cutoff radius
	  "ERR_COORD_CUTOFF",

  //! Class fileManipulation errors
  "ERR_FILE",
    // thrown if the file cannot be opened [from all constructors of fileManipulation in filemanipulation.cc]
    "ERR_FILE_OPEN",
    // thrown if reading a file is not possible due to a stream that is not ready to read [from fileManipulation::readFile() in filemanipulation.cc]
    "ERR_FILE_READ",
    // thrown if a line shall be accessed that is out of scope (> max num lines) [from fileManipulation::getLine() in filemanipulation.cc]
    "ERR_FILE_INVALID_ACCESS",
    // thrown if writing to a file is not possible due to a stream that is not ready to write 
    "ERR_FILE_WRITE",

  //! Energy class errors
  "ERR_ENERGY",
    // The specified interface is not available
    "ERR_ENERGY_UNKNOWN_INTERFACE",
    "ERR_ENERGY_MOPAC_EXE_NOT_FOUND",
    "ERR_ENERGY_MOPAC_DNF",
    "ERR_ENERGY_MOPAC_MISSING_OUTPUT",
    "ERR_ENERGY_MOPAC_OUTPUT_ERROR",
    "ERR_ENERGY_TERACHEM_INIT_MPI",
    "ERR_ENERGY_TERACHEM_SEND_MPI",
    "ERR_ENERGY_TERACHEM_RECIEVE_MPI",
    "ERR_ENERGY_TERACHEM_SCF_DNC",

  //! Class Linkedcell errors
  "ERR_LINKEDCELL",
    // thrown if the number of elements is insufficient to use linkedcells
    "ERR_LINKEDCELL_INSUFFICIENT_ELEMENTS", 
    // thrown if the position data is missing for initialization
    "ERR_LINKEDCELL_INIT",
    // thrown if updating the linkedcell information is not possible due to a lack of data or if the initialization has not been done yet
    "ERR_LINKEDCELL_UPDATE",
    "ERR_LINKEDCELL_LINK",
    "ERR_LINKEDCELL_PTR_NULL",
    
  //! Class parameters errors
  "ERR_PARAMETER",
    "ERR_PARAMETER_FILE",
    // thrown if the forcefield is unknown
    "ERR_PARAMETER_UNKNOWN_FF",
    // thrown if a known FF could not be identified
    "ERR_PARAMETER_NO_KNOWN_FF",
    // thrown if no vdw params can be found
    "ERR_PARAMETER_NO_VDW",
    // thrown if no charge params can be found
    "ERR_PARAMETER_NO_CHARGE",

  "ERR_STARTOPT",
    "ERR_STARTOPT_SOLVADD_FAIL",

  "ERR_SOLVADD",
    "ERR_SOLVADD_WATER_PARAM_BOND",
    "ERR_SOLVADD_WATER_PARAM_ANGLE",
    "ERR_SOLVADD_WATER_SET_OVERFLOW",
    "ERR_SOLVADD_WATER_NUM_UNDERFLOW",

  "ERR_RINGSEARCH",

  "ERR_FOLD",
    
  "ERR_VECT",
   "ERR_VECT_OUT_OF_RANGE",

  "ERR_VECT_INT",

  "ERR_MD",
    "ERR_MD_OUTPUT_FILE",

  "ERR_INPUT",
    "ERR_INPUT_UNKNOWN_INPUT_TYPE",

  "ERR_RDF",
    // thrown if a certain distance does not fit into the dists vector due to its size
    "ERR_RDF_BOX_OVERFLOW",

  "ERR_DIMER",
    "ERR_DIMER_NUMBER",
    "ERR_DIMER_DISTANCE",
    "ERR_DIMER_MAX_ITERATIONS",

  "ERR_PATH",
    "ERR_PATH_RANGE"
};



#endif
