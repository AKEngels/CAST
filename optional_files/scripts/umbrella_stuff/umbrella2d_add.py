### script for graphical analysis of 2D Umbrella Sampling ###
### it plots one direction and adds up the other
### for the probability this is the correct way to determine the overall probability
### for the Free Energy the values it gives are not reasonable
### because the correct way to do this would be A = -k_B * T * (P1+P2+...)
### but the general shape of the curve can give a hint 

import math
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# USER INPUT
DIRECTION = "X"             # direction that should be plotted
STEPS = range(-180, 185, 5) # possible values for that direction

######################################################################

P = [0]*len(STEPS)             # list for saving probability
A = [0]*len(STEPS)             # list for saving sum of free energy

# read outputfile and fill matrices
with open("out.txt") as outfile:
    lines = outfile.readlines()

for line in lines:
    if line.startswith("#") == False and len(line) > 1:
        linelist = line.split()
        if (DIRECTION == "Y"):
            i = STEPS.index(float(linelist[1]))
            P[i] += float(linelist[3])
            A[i] += float(linelist[2])
        if (DIRECTION == "X"):
            i = STEPS.index(float(linelist[0]))
            P[i] += float(linelist[3])
            A[i] += float(linelist[2])

# plot probability
plt.plot(STEPS, P, linewidth=2, color='black')
plt.xlabel('reaction coordinate {}'.format(DIRECTION))
plt.ylabel('Probability P')
plt.savefig("Probability_{}.png".format(DIRECTION))
plt.close()

# plot Free Energy
plt.plot(STEPS, A, linewidth=2, color='black')
plt.xlabel('reaction coordinate {}'.format(DIRECTION))
plt.ylabel('Free Energy A [sum]')
plt.savefig("FreeEnergy_{}.png".format(DIRECTION))
plt.close()


