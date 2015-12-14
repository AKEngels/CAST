#ifndef ATOMIC_H

#define ATOMIC_H 87U

#include <string>
#include <cmath>
#include <cstddef>
  
namespace atomic 
{

  static const std::string symbolMap[ATOMIC_H] = 
  { "XX", 
   "H",                                                                                                                                                                                       "He", 
   "Li", "Be",                                                                                                                                                  "B",  "C",  "N",  "O",  "F",  "Ne", 
   "Na", "Mg",                                                                                                                                                  "Al", "Si", "P",  "S",  "Cl", "Ar",
   "K",  "Ca",                                                                                      "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
   "Rb", "Sr",                                                                                      "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
   "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",  "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn" 
  };

  static const std::size_t num_saturated_bonds[ATOMIC_H] =
  { 0,
    1,                                                                                           0,
    0, 0,                                                                         3, 4, 3, 2, 1, 0,
    0, 0,                                                                         3, 4, 3, 2, 1, 0,
    0, 0,                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0,
    0, 0,                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

  enum nameMap 
  {
    HYDROGEN=1,                                                                                                                                                                            HELIUM, 
    LITHIUM,   BERYLLIUM,                                                                                                           BORON,     CARBON,    NITROGEN,   OXYGEN,    FLUORINE, NEON, 
    SODIUM,    MAGNESIUM,                                                                                                           ALUMINIUM, SILICON,   PHOSPHORUS, SULFUR,    CHLORINE, ARGON,
    POTASSIUM, CALCIUM,   SCANDIUM, TITANIUM,  VANADIUM, CHROMIUM,   MANGANESE,  IRON,      COBALT,   NICKEL,    COPPER,  ZINC,     GALLIUM,   GERMANIUM, ARSENIC,    SELENIUM,  BROMINE,  KRYPTON,
    RUBIDIUM,  STRONTIUM, YTTRIUM,  ZIRCONIUM, NIOBIUM,  MOLYBDENUM, TECHNETIUM, RUTHENIUM, RHODIUM,  PALLADIUM, SILVER,  CADMIUM,  INDIUM,    TIN,       ANTIMONY,   TELLURIUM, IODINE,   XENON,
    CAESIUM,   BARIUM, LANTHANUM, 
                                    CERIUM,  PRASEODYMIUM, NEODYMIUM, PROMETHIUM, SAMARIUM, EUROPIUM, GADOLINIUM, TERBIN, DYSPROSIUM, HOLMIUM, ERBIUM, THULIUM, YTTERBIUM, LUTETIUM,
                                    HAFNIUM, TANTALUM,     WOLFRAM,   RHENIUM,    OSMIUM,   IRIDIUM,  PLATINUM,   GOLD,   MERCURY,    THALLIUM,  LEAD, BISMUTH, POLONIUM,   ASTATINE, RADON
  };

  static const double massMap[ATOMIC_H] = 
  { 0.0, 
    1.00794,                                                                                                                                                          4.0026, 
    6.941,     9.012182,                                                                                              10.811,  12.0107, 14.0067, 15.9994,  18.9984,  20.1797, 
   22.989770, 24.3050,                                                                                                26.9815, 28.0855, 30.9737, 32.065,   35.453,   39.948,
   39.0983,   40.078,  /* ... */ 44.9559, 47.867,  50.941,   51.996,  54.938,   55.845,  58.9332,  58.693,  63.546,   65.409,  69.723,  72.64,   74.9216, 78.96,    79.904,   83.798,
   85.4678,   87.62,  /* ... */  88.9059, 91.224,  92.906,   95.94,   98.9063, 101.07,  102.905,  106.42,  107.868,  112.411, 114.818, 118.710, 121.760, 127.60,   126.904,  131.293,
  132.90545, 137.327,  138.906, 140.116, 140.908,  144.24,  146.915,  150.36,  151.964,  157.25 , 158.925,  162.5,   164.93,  167.259, 168.934, 173.04,   174.967,
                                178.49,  180.9479, 183.84,  186.207,  190.23,  192.217,  195.078, 196.9665, 200.59,  204.383, 207.2,   208.908, 208.9824, 209.9871, 222.0176 
  };

  static const double radiusMap[ATOMIC_H] =
  {
    0.0,
    1.25 /*1*/, 2.0 /*2*/, 1.432 /*3*/, 2.0 /*4*/, 2.0 /*5*/,
    1.9 /*6*/, 1.7063 /*7*/, 1.535 /*8*/, 1.47 /*9*/, 1.39 /*10*/,
    1.992 /*11*/, 1.7 /*12*/, 2.0 /*13*/, 1.8 /*14*/, 1.87 /*15*/,
    1.775 /*16*/, 1.735 /*17*/, 1.7 /*18*/, 2.123 /*19*/, 1.817 /*20*/,
    2.0 /*21*/, 2.0 /*22*/, 2.0 /*23*/, 2.0 /*24*/, 2.0 /*25*/,
    2.0 /*26*/, 2.0 /*27*/, 2.0 /*28*/, 2.0 /*29*/, 2.0 /*30*/,
    2.0 /*31*/, 2.0 /*32*/, 2.0 /*33*/, 2.0 /*34*/, 1.9 /*35*/,
    1.812 /*36*/, 2.26 /*37*/, 2.0 /*38*/, 2.0 /*39*/, 1.9 /*40*/,
    1.812 /*41*/, 2.26 /*42*/, 2.0 /*43*/, 2.0 /*44*/, 1.9 /*45*/,
    1.812 /*46*/, 2.26 /*47*/, 2.0 /*48*/, 2.0 /*49*/, 1.9 /*50*/,
    1.812 /*51*/, 2.26 /*52*/, 2.0 /*53*/, 1.967 /*54*/, 2.507 /*55*/,
    2.188 /*56*/, 2.0 /*57*/, 2.0 /*58*/, 2.0 /*59*/, 2.0 /*60*/,
    2.0, 2.0, 2.0, 2.0, 2.0 /*61-65*/, 2.0, 2.0, 2.0, 2.0, 2.0 /*66-70*/,
    2.0, 2.0, 2.0, 2.0, 2.0 /*71-75*/, 2.0, 2.0, 2.0, 2.0, 2.0 /*76-80*/,
    2.0, 2.0, 2.0, 2.0, 2.0 /*81-85*/, 2.0 /*86*/
  };

  inline std::size_t atomic_number_by_symbol (std::string const & symbol)
  {
    std::size_t ret(0U);
    for (std::size_t i(0U); i<ATOMIC_H; ++i)
    {
      if (symbol.find(atomic::symbolMap[i]) == 0U) 
      {
        if (atomic::symbolMap[i].length() == 2U) return i;
        else ret = i;
      }
    }
    return ret;
  }

  inline std::size_t atomic_number_by_mass (double const value)
  {
    if (value < 0.0 || atomic::massMap[ATOMIC_H-1U] < value) return 0U;
    std::size_t low(0U), high(ATOMIC_H-1U), buffer(0U);
    while (low < high) 
    {
      buffer = (low+high)/2U;
      if (atomic::massMap[buffer] < value) low = buffer+1U;
      else high = buffer;
    }
    if (std::fabs(atomic::massMap[low] - value) < 1.0) return low;
    if ((low < ATOMIC_H-1U) && std::fabs(atomic::massMap[low+1U] - value) < 1.0) return low+1U;
    if ((low > 1U) && std::fabs(atomic::massMap[low-1U] - value) < 1.0) return low-1U;
    return 0;
  } 

  inline bool number_is_heteroatom (std::size_t const atomic_number) 
  {
    return atomic_number ==  7u || atomic_number ==  8u || atomic_number ==  9u || 
           atomic_number == 15u || atomic_number == 16u || atomic_number == 17u || 
           atomic_number == 34u || atomic_number == 35u || atomic_number == 53u;
  }

  inline bool saturated(std::size_t const atomic_number, std::size_t const bonds)
  {
    return num_saturated_bonds[atomic_number >= ATOMIC_H ? 0u : atomic_number] == bonds;
  }

}

#endif
