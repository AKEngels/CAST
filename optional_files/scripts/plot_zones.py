### This script can help you to plot the temperature of a given region.
### For this it uses the CAST outputfile ..._MD_VELO.arc where the velocities of all atoms are saved.
### You can either define the region yourself by giving the atom indices for ATOMS or define an active center and a radius.
### In the second case the region consists of all atoms inside a sphere around the active site with the given radius.
### Then you must also give a structurefile where the atomic coordinates are taken from (in tinker format).
### In any case you have to give the number of steps you want to plot (plotting starts at the beginning of the MD) and the total number of atoms in the structure.
### The output is a file called temp.csv which contains the temperatures for every MD step.

### ATTENTION! For some reason I don't know the results are different from those received from CAST when directly calculating and plotting those temperatures!
###            So either in this script or in CAST there is an error.



### USER INPUT

# first part (if you want to take atoms from zones around a center)
STRUCTUREFILE = "21_13_equil.arc"
CENTER = [29.3636,31.0752,36.1577]
RADIUS = 3    # in angstrom

# second part
USE_FIRST = True
ATOMS = []    # does not need to be filled if first part is used
STEPS = 50    # number of MD steps plotted
VELOFILE = "output_MD_VELO.arc"
NUMBER_OF_ATOMS = 3176

# atomic masses (to be complemented if other elements are in the structure)
MASSES = {"H" : 1.0079, "C":12.0107, "O" :15.9994, "N" : 14.0067, "S" : 32.065, "Na": 22.9897}


### PROGRAM

import math    

# function that returns the indices of all atoms that are inside the given region (first part)
def get_atoms(filename, center, radius):

    # class that contains necessary information about an atom
    class Atom:
        def __init__(self, index, x, y, z):
            self.index = index
            self.x = x
            self.y = y
            self.z = z
        def calc_dist(self, point):  # calculate distance of the atom to a given point [x, y, z]
            return math.sqrt((point[0]-self.x)*(point[0]-self.x) + (point[1]-self.y)*(point[1]-self.y) + (point[2]-self.z)*(point[2]-self.z))
    
    with open(STRUCTUREFILE) as struc:
        lines = struc.readlines()
    atoms_in_zone = []
    for i, line in enumerate(lines):
        if i != 0:
            a = Atom(int(line.split()[0]), float(line.split()[2]), float(line.split()[3]), float(line.split()[4]))
            if a.calc_dist(center) < radius:
                atoms_in_zone.append(a.index)
    return atoms_in_zone


# calculates the temperature of a given region at a certain step from the atomic velocities in the velofile
# filename: name of the velofile (CAST output)
# atoms: atom indices that define the region
# step: MD step for which temperature should be calculated
# number_of_atoms: number of atoms in the system
# atom_masses: dictionary that maps the atomic masses to the element symbols
def calc_temp(filename, atoms, step, number_of_atoms, atom_masses):

    # translates atomic indices to the line number in the velofile taking into account MD step
    def translate_atoms_to_linenumbers(atoms, step, number_of_atoms):
        numbers = []
        for a in atoms:
            numbers.append(a + step*(number_of_atoms+1))
        return numbers

    # reads the velofile and returns only those lines that are needed to calculate the temperature,
    # i.e. those that correspond to the given atoms in the current MD step
    def get_current_lines(atoms, filename, step, number_of_atoms):
        count_lines = 0
        lines = []
        atom_lines = translate_atoms_to_linenumbers(atoms, step, number_of_atoms)
        with open(filename) as velo:
            while len(lines) < len(atoms):
                line = velo.readline()
                if count_lines in atom_lines:
                    lines.append(line)
                count_lines += 1
        return lines

    # calculates kinetic energy from the lines extracted in the former step
    def calc_Ekin(lines, atom_masses):
        e_kin = 0
        for i,line in enumerate(lines):
            print line
            v_x, v_y, v_z = float(lines[i].split()[2]), float(lines[i].split()[3]), float(lines[i].split()[4])
            v = math.sqrt(v_x*v_x + v_y*v_y + v_z*v_z)
            m = MASSES[lines[i].split()[1]]
            e_kin += 0.5*m*v*v / 418.4    # 418.4 is the same conversion factor as used in CAST
        return e_kin

    lines = get_current_lines(atoms, filename, step, number_of_atoms)
    ekin = calc_Ekin(lines, atom_masses)
    temp = 2*ekin / (3*len(lines)*1.9872066e-3)  # the last number is also a factor taken from CAST
    return temp
    
            
if USE_FIRST:
    ATOMS = get_atoms(STRUCTUREFILE, CENTER, RADIUS)

print ATOMS
with open("temp.csv", "w") as out:
    out.write("temp\n")
    for step in range(3):
        temp = calc_temp(VELOFILE, ATOMS, step, NUMBER_OF_ATOMS, MASSES)
        print temp
        out.write("{}\n".format(temp))
