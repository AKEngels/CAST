### script to extract xyz coordinates of gaussian outputfile

inputfile = "water.log"
number_of_atoms = 3

"""function to get the name of the xyz file
inputlines: gaussian logfile read in with readlines()
linenumber: index of line that contains the title for the z-matrix (Z-MATRIX (ANGSTROMS AND DEGREES))

you should modify this function before using the script"""
def get_title(inputlines, linenumber):
    return inputlines[linenumber+2+number_of_atoms][54:59].replace(".","")

"""function to convert atomnumber (from PSE) to atom symbol
currently only hydrogen and oxygen implemented"""
def get_atom_symbol_from_atomic_number(an):
    if an == 1:
        return "H"
    elif an == 8:
        return "O"


# script is starting here

with open(inputfile) as inp: # read inputfile
    data = inp.readlines()

coordline = [0,0,0]
for i,line in enumerate(data):
    if "Z-MATRIX (ANGSTROMS AND DEGREES)" in line: # find structures (z-matrix)
        title = get_title(data,i)
        for j in range(number_of_atoms): # get line with atom number and xyz coordinates for every atom in structre
            coordline[j] = data[i+number_of_atoms+j+9].split()
            
        with open("{}.xyz".format(title), "w") as out:  # write xyz file
            out.write(str(number_of_atoms))
            out.write("\n\n")
            for a in range(number_of_atoms):
                atom_symbol = get_atom_symbol_from_atomic_number(int(coordline[a][1]))
                out.write(str(atom_symbol))
                out.write("\t")
                out.write(coordline[a][3])
                out.write("\t")
                out.write(coordline[a][4])
                out.write("\t")
                out.write(coordline[a][5])
                out.write("\n")
