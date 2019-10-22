### With this script you can process the file 'pca_stocksdelta.dat' from task PCAgen.
###
### It reads the contribution of every atom for a given node (starting with 0)
### and sorts the atoms by how much they contribute to this node.
### The result is written into a file 'node<nodenumber>_atoms.csv'.
###
### Furthermore it prints the first X atoms in a format that you can copy into VMD
### where you can choose X as you desire.

MODE = 0
NUMBER_OF_ATOMS = 100

class Atom:

    def __init__(self, index, weight):
        self.index = index
        self.weight = weight

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

with open("pca_stocksdelta.dat") as inp:
    lines = inp.readlines()

i = find_modeline(lines)
atoms = create_atomlist(i)
atoms = sorted(atoms, key=lambda atom: atom.weight, reverse=True)
write_to_file(atoms, "node{}_atoms.csv".format(MODE))
get_vmd_indices(NUMBER_OF_ATOMS, atoms)


