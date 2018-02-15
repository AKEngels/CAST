"""
You can use this script to plot how the distance of chosen atom pairs changes during an MD simulation.
Copy this script into the folder where your CAST output file '*_MD_SNAP.arc' is, then change the FILENAME and the atom_number_pairs
(atom numbers starting with 1 as in tinkerstructure) and run it. If you want to can also change the number of FRAMES you want to plot.

Script could be sped up if it would search the lines from the output file only once and not once for every atom pair.
(relevant mainly if more than a few atom pairs are analyzed but then the picture would also be confusing)
"""

import math
import matplotlib.pyplot as plt

FILENAME = "output_MD_SNAP.arc"
ATOM_NUMEBR_PAIRS = [[347,3149], [2277,3151], [2277,2276]]  # pairs of atoms you want to analyze

"""class for every atom pair"""
class AtomPair(object):
    def __init__(self, atom_numbers):
        self.atom_numbers = atom_numbers  # atom numbers of pair
        self.distances = []               # distances for every step of MD
        self.symbols = ["dummy","dummy"]  # element symbols of the two atoms
        

with open(FILENAME) as mdfile: # read the output file
    lines = mdfile.readlines()

N = int(lines[0])         # get atom number
FRAMES = len(lines)/(N+1) # get number of frames (can be changed if you only want to plot the first ... frames)

atom_pairs = []

for a in ATOM_NUMEBR_PAIRS: # for every atom pair
    pair = AtomPair(a)

    for j in range(FRAMES): # for every frame in MD
        
        for i in range(N):  # search and save to corresponding lines
            line = lines[i+1+(j*(N+1))]
            if i == pair.atom_numbers[0]-1: 
                linelist1 = line.split()
            elif i == pair.atom_numbers[1]-1:
                linelist2 = line.split()

        x1 = float(linelist1[2])  # get atom coordinates
        y1 = float(linelist1[3])
        z1 = float(linelist1[4])
        x2 = float(linelist2[2])
        y2 = float(linelist2[3])
        z2 = float(linelist2[4])
        
        dist = math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) # calculate distance and save
        pair.distances.append(dist)

    pair.symbols[0] = linelist1[1]  # get element symbols from last frame 
    pair.symbols[1] = linelist2[1]

    atom_pairs.append(pair)

legends = []
for p in atom_pairs: # for every atom pair
    plt.plot(p.distances) # plot the distances
    legend = "{}{}-{}{}".format(p.symbols[0],p.atom_numbers[0],p.symbols[1],p.atom_numbers[1])
    legends.append(legend) # create and save the legend

plt.xlabel("frame")
plt.ylabel("distance [$\AA$]")
plt.legend(legends)
plt.savefig("distances.png")
        
