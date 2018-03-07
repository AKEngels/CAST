### script to calculate the distance between 2 atoms in a PDB file
### just enter the atom indizes (starting with 1), the script will print the distance

ATOMS = [2276, 2277]

import math

def calc_dist(lines, atoms):
    atomline1 = lines[atoms[0]-1]
    x1 = float(atomline1[30:38])
    y1 = float(atomline1[38:46])
    z1 = float(atomline1[46:54])
    atomline2 = lines[atoms[1]-1]
    x2 = float(atomline2[30:38])
    y2 = float(atomline2[38:46])
    z2 = float(atomline2[46:54])

    res = math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
    return res

# open PDB file
with open("qmmm_opt.pdb") as pdbfile:
    lines = pdbfile.readlines()

print calc_dist(lines, ATOMS)
