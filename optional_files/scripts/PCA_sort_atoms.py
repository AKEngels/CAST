### With this script you can process the file 'pca_stocksdelta.dat' from task PCAgen.
###
### It reads the contribution of every atom for a given mode (starting with 0)
### and sorts the atoms by how much they contribute to this mode.
### The result is written into a file 'node<nodenumber>_atoms.csv'.
###
### Furthermore it prints the first X atoms in a format that you can copy into VMD
### where you can choose X as you desire.
###
### The last thing this script does is to write a tracefile of tinkerstructures
### where the last structure contains only the atom with the biggest contribution
### the structure before the two atoms with the biggest contribution and so on
### until the first structure contains all atoms that contribute to the mode.

MODE = 0
NUMBER_OF_ATOMS = 100
STRUCTURE = "ncov_full_19_21.arc"

import os

class Atom:

    def __init__(self, index, weight):
        self.index = index
        self.weight = weight

class AtomWithCoords:

    def __init__(self, element, x, y, z):
        self.element = element
        self.x = x
        self.y = y
        self.z = z

def find_modeline(lines):
    for i,line in enumerate(lines):
        if (line.startswith("PCA Mode {}".format(MODE))):
            return i

def create_atomlist(starting_index):
    line = lines[starting_index+1]

    atoms = []
    while line.startswith("---") == False:
        index = int(line.split()[1][:-1])
        weight = float(line.split()[2])
        atom = Atom(index, weight)
        atoms.append(atom)
        starting_index+=1
        line = lines[starting_index]

    return atoms

def write_to_file(atoms, filename):
    with open(filename, "w") as out:
        out.write("Atom,Weight\n")
        for a in atoms:
            out.write("{},{}\n".format(a.index, a.weight))

def get_vmd_indices(number, atoms):
    string = ""
    for i in range(number):
        string += "index " + str(atoms[i].index-1) + " or "
    print string[:-4]

def read_atoms_from_structure(filename):
    with open(filename) as inp:
        lines = inp.readlines()

    number_of_atoms = int(lines[0])

    atoms = []
    for i in range(number_of_atoms):
        line = lines[i+1]
        linelist = line.split()
        element = linelist[1]
        x = float(linelist[2])
        y = float(linelist[3])
        z = float(linelist[4])
        atom = AtomWithCoords(element, x, y, z)
        atoms.append(atom)
        i=i+1
    return atoms

def create_partial_structure(indices, structure, total_size):
    string = "{}\n".format(total_size)
    for i,index in enumerate(indices):
        atom = structure[index-1]
        string += "{}  {}  {}  {}  {}\n".format(i+1, atom.element, atom.x, atom.y, atom.z)
    number_of_dummys = total_size - len(indices)
    for i in range(number_of_dummys):
        string += "{}   H   0  0   0 \n".format(i+1+len(indices))
    return string

def create_strings_for_all_partial_structures(atoms, atoms_structure):
    tinkerstructures = []
    for i, a in enumerate(atoms):
        indices_for_structure = []
        for j in range(i+1):
            indices_for_structure.append(atoms[j].index)
        tinkerstructures.append(create_partial_structure(indices_for_structure, atoms_structure, len(atoms)))
    return tinkerstructures

def write_tinkerstructures_to_file(tinkerstructures):
    for i,t in enumerate(tinkerstructures):
        index = len(tinkerstructures) - i -1
        with open("mode{}.arc".format(MODE), "a") as out:
            out.write(tinkerstructures[index])

with open("pca_stocksdelta.dat") as inp:
    lines = inp.readlines()

i = find_modeline(lines)
atoms = create_atomlist(i)
atoms = sorted(atoms, key=lambda atom: atom.weight, reverse=True)
write_to_file(atoms, "mode{}_atoms.csv".format(MODE))
get_vmd_indices(NUMBER_OF_ATOMS, atoms)

os.remove("mode{}.arc".format(MODE))
atoms_structure = read_atoms_from_structure(STRUCTURE)
tinkerstructures = create_strings_for_all_partial_structures(atoms, atoms_structure)
write_tinkerstructures_to_file(tinkerstructures)




