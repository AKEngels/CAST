### script that takes a tinkerstructure from a peptide and cuts in into single amino acids, writes them into xyz-files

STRUCTURENAME = "35_10.arc"

class Atom:

    def __init__(self, tinker_line):
        linelist = line.split()
        self.symbol = linelist[1]
        self.x = linelist[2]
        self.y = linelist[3]
        self.z = linelist[4]

        bondstr = linelist[6:]
        self.bonds = []
        for b in bondstr:
            self.bonds.append(int(b)-1)

    def bonds_as_str(self):
        res = ""
        for b in self.bonds:
            res += str(b+1) + " "
        return res

    def write_to_xyz(self):
        return "{}  {}  {}  {}\n".format(self.symbol, self.x, self.y, self.z) 

    def info(self):
        print "{} at {}, {}, {} with bonding partners {}".format(self.symbol, self.x, self.y, self.z, self.bonds_as_str())


### get atomic information from tinkerstructure 
with open(STRUCTURENAME) as tinkerfile:
    lines = tinkerfile.readlines()

atoms = []
for i, line in enumerate(lines):
    if i != 0:
        a = Atom(line)
        atoms.append(a)

### find out where the peptide bonds are (atoms in the order C, O, N, H, "cuts" between O and N)
indices_of_new_peptide = []
for i in range(len(atoms)-3):
    if atoms[i].symbol == "C" and atoms[i+1].symbol == "O" and atoms[i+2].symbol == "N" and atoms[i+3].symbol == "H":
        indices_of_new_peptide.append(i+2)

### get the indices of atoms that are in each segment and write segments into file "bla.txt"
segment_ranges = []
with open("bla.txt","w") as segfile:
    for i, index in enumerate(indices_of_new_peptide):
        if i == 0:
            segment_ranges.append(range(index))
        else:
            segment_ranges.append(range(indices_of_new_peptide[i-1], index))
        segfile.write("Segment {} : {} to {} ({} atoms):\n".format(i+1, segment_ranges[-1][0]+1, segment_ranges[-1][-1]+1, len(segment_ranges[-1])))
    segment_ranges.append(range(indices_of_new_peptide[-1], len(atoms)))
    segfile.write("Segment {} : {} to {} ({} atoms):\n".format(i+1, segment_ranges[-1][0]+1, segment_ranges[-1][-1]+1, len(segment_ranges[-1])))

### writes an xyz structure for every segment
for i, s in enumerate(segment_ranges):
    number_of_atoms = len(s)
    with open("segment_{}.xyz".format(i), "w") as output:
        output.write(str(number_of_atoms) + "\n\n")
        for a in s:
            atom = atoms[a]
            output.write(atom.write_to_xyz())
    
    

            
