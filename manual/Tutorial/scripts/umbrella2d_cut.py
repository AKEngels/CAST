### script for graphical analysis of 2D Umbrella Sampling ###

import math
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# USER INPUT
DIRECTION = "X"   # direction that is fixed, can be X or Y
VALUE = 175       # value to which reaction coordinate is set


######################################################################

xi = []            # reaction coordinate that is not fixed
A = []             # list for saving free energy
P = []             # list for saving probability

# read outputfile and fill matrices
with open("out.txt") as outfile:
    lines = outfile.readlines()

for line in lines:
    if line.startswith("#") == False and len(line) > 1:
        linelist = line.split()
        x_value = float(linelist[0])
        y_value = float(linelist[1])
        if (DIRECTION == "Y" and y_value == VALUE):
            xi.append(x_value)
            A.append(float(linelist[2]))
            P.append(float(linelist[3]))
        if (DIRECTION == "X" and x_value == VALUE):
            xi.append(y_value)
            A.append(float(linelist[2]))
            P.append(float(linelist[3]))

if DIRECTION == 'X':
    plot_direction = 'Y'
elif DIRECTION == 'Y':
    plot_direction = 'X'

# plot Free Energy
plt.plot(xi, A, linewidth=2, color='black')
plt.xlim(min(xi),max(xi))
plt.xlabel('reaction coordinate {}'.format(plot_direction))
plt.ylabel('Free Energy A [kcal/mol]')
plt.savefig("FreeEnergy_{}{}.png".format(DIRECTION, VALUE))
plt.close()

# plot probability
plt.plot(xi, P, linewidth=2, color='black')
plt.xlim(min(xi),max(xi))
plt.xlabel('reaction coordinate {}'.format(plot_direction))
plt.ylabel('Probability P')
plt.savefig("Probability_{}{}.png".format(DIRECTION, VALUE))
plt.close()
