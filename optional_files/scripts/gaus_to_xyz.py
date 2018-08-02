### script to extract xyz coordinates of gaussian outputfile

inputfile = "ethen.log"
number_of_atoms = 6

"""function to get the name of the xyz file
inputlines: gaussian logfile read in with readlines()
linenumber: index of line that contains the title for the z-matrix (Z-MATRIX (ANGSTROMS AND DEGREES))

you should modify this function before using the script"""
def get_title(inputlines, linenumber):
    return inputlines[linenumber+1+number_of_atoms][78:81]


# list with atom symbols for every element in PSE
PSE = [ "XX", 
        "H",                                                                                                                                                                                       "He", 
        "Li", "Be",                                                                                                                                                  "B",  "C",  "N",  "O",  "F",  "Ne", 
        "Na", "Mg",                                                                                                                                                  "Al", "Si", "P",  "S",  "Cl", "Ar",
        "K",  "Ca",                                                                                      "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr",                                                                                      "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",  "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn" ]

# script is starting here

with open(inputfile) as inp: # read inputfile
    data = inp.readlines()

coordline = [0]*number_of_atoms
for i,line in enumerate(data):
    if "Z-MATRIX (ANGSTROMS AND DEGREES)" in line: # find structures (z-matrix)
        title = get_title(data,i)
        for j in range(number_of_atoms): # get line with atom number and xyz coordinates for every atom in structre
            coordline[j] = data[i+number_of_atoms+j+9].split()
            
        with open("{}.xyz".format(title), "w") as out:  # write xyz file
            out.write(str(number_of_atoms))
            out.write("\n\n")
            for a in range(number_of_atoms):
                atom_symbol = PSE[int(coordline[a][1])]
                out.write(str(atom_symbol))
                out.write("\t")
                out.write(coordline[a][3])
                out.write("\t")
                out.write(coordline[a][4])
                out.write("\t")
                out.write(coordline[a][5])
                out.write("\n")
