"""
You can use this script to plot how the angle between specified atoms changes during an MD simulation.
Copy this script into the folder where your CAST output file '*_MD_SNAP.arc' is, then change the FILENAME and the ANGLES
(atom numbers starting with 1 as in tinkerstructure) and run it. If you want to can also change the number of FRAMES you want to plot.

Script could be sped up if it would search the lines from the output file only once and not once for every angle.
(relevant mainly if more than a few angles are analyzed but then the picture would also be confusing)
"""

import math
import matplotlib.pyplot as plt

FILENAME = "output_MD_SNAP_1.arc"
ANGLES = [[2, 1, 3], [5, 4, 6]]  # pairs of atoms you want to analyze

"""class for every atom pair"""
class Angle(object):
    def __init__(self, atom_numbers):
        self.atom_numbers = atom_numbers  # atom numbers of pair
        self.symbols = None               # element symbols of the three atoms
        self.angles = []                  # angle for every MD step

    def read_atom_symbols(self, linelist_A, linelist_B, linelist_C):
        symbol_A = linelist_A[1]
        symbol_B = linelist_B[1]
        symbol_C = linelist_C[1]
        self.symbols = [symbol_A, symbol_B, symbol_C]

    def calc_current_angle(self, linelist_A, linelist_B, linelist_C):
        # get atomic coordinates
        x_a = float(linelist_A[2])  
        y_a = float(linelist_A[3])
        z_a = float(linelist_A[4])
        x_b = float(linelist_B[2])
        y_b = float(linelist_B[3])
        z_b = float(linelist_B[4])
        x_c = float(linelist_C[2])
        y_c = float(linelist_C[3])
        z_c = float(linelist_C[4])

        # calculate angle
        scalar_product = (x_a - x_b) * (x_c - x_b) + (y_a - y_b) * (y_c - y_b) + (z_a - z_b) * (z_c - z_b)
        length_of_AB = math.sqrt( (x_a-x_b)*(x_a-x_b) + (y_a-y_b)*(y_a-y_b) + (z_a-z_b)*(z_a-z_b) )
        length_of_BC = math.sqrt( (x_c-x_b)*(x_c-x_b) + (y_c-y_b)*(y_c-y_b) + (z_c-z_b)*(z_c-z_b) )
        cos = scalar_product / (length_of_AB * length_of_BC)
        alpha = math.acos(cos) * (180.0 / math.pi)  # convert to degree

        # add angle to list
        self.angles.append(alpha)
        

with open(FILENAME) as mdfile: # read the output file
    lines = mdfile.readlines()

N = int(lines[0])         # get atom number
FRAMES = len(lines)/(N+1) # get number of frames (can be changed if you only want to plot the first ... frames)

angles = []

for a in ANGLES: # for every atom pair
    angle = Angle(a)

    for j in range(FRAMES): # for every frame in MD
        
        for i in range(N):  # search and save to corresponding lines
            line = lines[i+1+(j*(N+1))]
            if i == angle.atom_numbers[0]-1: 
                linelist1 = line.split()
            elif i == angle.atom_numbers[1]-1:
                linelist2 = line.split()
            elif i == angle.atom_numbers[2]-1:
                linelist3 = line.split()

        if angle.symbols == None:  # read symbols
            angle.read_atom_symbols(linelist1, linelist2, linelist3)
            
        # calculate angle
        angle.calc_current_angle(linelist1, linelist2, linelist3) 

    angles.append(angle)

# write angles to .csv file
with open("angles.csv","w") as anglefile:
    anglefile.write("Schritte")
    for a in angles:
        anglefile.write(",{}{}-{}{}-{}{}".format(a.symbols[0],a.atom_numbers[0],a.symbols[1],a.atom_numbers[1],a.symbols[2],a.atom_numbers[2]))
    anglefile.write("\n")

    for f in range(FRAMES):
        anglefile.write("{}".format(f+1))
        for a in angles:
            anglefile.write(",{}".format(a.angles[f]))
        anglefile.write("\n")
