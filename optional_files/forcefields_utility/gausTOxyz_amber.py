### This script converts a gaus.com file for an ONIOM calculation into an xyz-file.

### Which atoms (in tinker-numbering) are in which system (low, middle, high) is printed to the standard output. You can copy them into the CAST.txt inputfile.

### As an xyz-file doesn't contain all necessary information for forcefield calculations, a few additional files are produced:
### charges.txt: a chargefile with amber charges that you can use together with the amber forcefield in CAST
### atomtypes.csv: a file with the atomtypes in tinker form. you can copy this into the atomtypes column in a tinkerstructure CAST produces from the generated xyz-file by the task WRITE_TINKER
###                the file 'atomtypes.txt' which you find in this folder is used for that



class Atom:

    def __init__(self, x, y, z, system, symbol, atomtype, charge):
        self.x = x
        self.y = y
        self.z = z
        self.system = system
        self.symbol = symbol
        self.atomtype = atomtype
        self.charge = charge


def create_typefinder():
    with open("atomtypes.txt") as typefile:
        lines = typefile.readlines()

    ambertypes = []
    tinkertypes = []
    for line in lines:
        linelist = line.split()
        ambertype = linelist[4][1:]
        tinkertype = linelist[-1][:-2]
        ambertypes.append(ambertype)
        tinkertypes.append(tinkertype)

    return dict(zip(ambertypes, tinkertypes))


with open("gaus.com") as gaussianfile:
    lines = gaussianfile.readlines()

read = False
for i,line in enumerate(lines):
    if line.startswith("0 1"):
        start = i+1
        read = True
    if read and line == "\n":
        end = i
        break

atoms = []
for atomline in lines[start:end]:
    linelist = atomline.split()
    x = linelist[2]
    y = linelist[3]
    z = linelist[4]
    system = linelist[5]

    newlist = linelist[0].split("-", 2)
    symbol = newlist[0]
    atomtype = newlist[1]
    charge = newlist[2].split("(")[0]
    
    atom = Atom(x,y,z,system, symbol, atomtype, charge)
    atoms.append(atom)

number_of_atoms = len(atoms)
typedict = create_typefinder()

with open("output.xyz","w") as xyzfile:
    xyzfile.write(str(number_of_atoms) + "\n\n")
    for a in atoms:
        xyzfile.write("{} {}  {}  {}\n".format(a.symbol, a.x,a.y,a.z))
    
with open("charges.txt","w") as chargefile:
    for i,a in enumerate(atoms):
        charge = float(a.charge) * 18.2223
        chargefile.write("{}  {}  {}\n".format(i+1, a.atomtype, charge))

with open("atomtypes.csv","w") as typefile:
    typefile.write("\n")
    for a in atoms:
        try:
            tinkertype = typedict[a.atomtype]
        except KeyError:
            if a.atomtype == "IP":
                print "WARNING! We believe that IP is sodium. If it's not the atom type for this atom is wrong."
                tinkertype = str(2004)
            else:
                print "unknown atom type:", a.atomtype
                tinkertype = str(0)
                
        typefile.write(tinkertype + "\n")

low = []
medium = []
high = []
for i,a in enumerate(atoms):
    if a.system == "L":
        low.append(i+1)
    elif a.system == "M":
        medium.append(i+1)
    elif a.system == "H":
        high.append(i+1)
    else:
        print "something went wrong"

print "high system contains {} atoms:".format(len(high))
print high

print "middle system contains {} atoms:".format(len(medium))
print medium

print "low system contains {} atoms".format(len(low))
    
