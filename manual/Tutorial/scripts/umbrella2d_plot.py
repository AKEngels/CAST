### script for graphical analysis of 2D Umbrella Sampling 
### plot heatmap and surface of Free Energy and Probability and saves csv files

import math
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# USER INPUT
X_MIN, X_MAX, X_STEP = -180, 180, 5
Y_MIN, Y_MAX, Y_STEP = -180, 180, 5

X_LABEL_STEP = 30  # has to be a multiple of X_STEP
Y_LABEL_STEP = 30  # has to be a multiple of Y_STEP

######################################################################

# function that writes a numpy array into a csv file
# RANGE_X and RANGE_Y are lists that contain the names of the columns and rows
def write_into_csv(numpy_array, RANGE_X, RANGE_Y, filename):
    with open(filename,'w') as csvfile:
        for x in RANGE_X:
            csvfile.write(',{}'.format(x))
        csvfile.write(",X\n")

        for i,y in enumerate(RANGE_Y):
            csvfile.write(str(y))
            for j in range(len(RANGE_X)):
                csvfile.write(',{}'.format(numpy_array[i][j]))
            csvfile.write(",\n")
            
        csvfile.write("Y")
        for j in range(len(RANGE_X)):
            csvfile.write(',')

# calculate some stuff
RANGE_X = range(X_MIN, X_MAX+X_STEP, X_STEP)
RANGE_Y = range(Y_MIN, Y_MAX+Y_STEP, Y_STEP)

FACTOR_X = X_LABEL_STEP / X_STEP
FACTOR_Y = Y_LABEL_STEP / Y_STEP 

LABELS_X = range(X_MIN, X_MAX+X_LABEL_STEP, X_LABEL_STEP)
LABELS_Y = range(Y_MIN, Y_MAX+Y_LABEL_STEP, Y_LABEL_STEP)

# create empty matrices for free energy and probability
free_energy = numpy.empty((len(RANGE_Y),len(RANGE_X)))
probability = numpy.empty((len(RANGE_Y),len(RANGE_X)))

# read outputfile and fill matrices
with open("out.txt") as outfile:
    lines = outfile.readlines()

for line in lines:
    if line.startswith("#") == False and len(line) > 1:
        linelist = line.split()
        x_value = float(linelist[0])
        y_value = float(linelist[1])
        A = float(linelist[2])
        P = float(linelist[3])

        x = RANGE_X.index(x_value)
        y = RANGE_Y.index(y_value)

        free_energy[y][x] = A
        probability[y][x] = P

# write free energy and probability into csv file
write_into_csv(free_energy,RANGE_X,RANGE_Y,"freeEnergy.csv")
write_into_csv(probability,RANGE_X,RANGE_Y,"Probability.csv")

# create a meshgrid from X and Y
X, Y = numpy.meshgrid(RANGE_X, RANGE_Y)

# plot free energy as surface
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xticklabels(LABELS_X)
ax.set_yticklabels(LABELS_Y)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Free Energy [kcal/mol]")
surf = ax.plot_surface(X, Y, numpy.asarray(free_energy), cmap='jet') # jet is colorscheme, can be changed
plt.savefig("FreeEnergy.png")
plt.close()

# plot probability as surface
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xticklabels(LABELS_X)
ax.set_yticklabels(LABELS_Y)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Probability")
surf = ax.plot_surface(X, Y, numpy.asarray(probability), cmap='jet') # jet is colorscheme, can be changed
plt.savefig("Probability.png")
plt.close()

# create heatmap with free energy
fig, ax = plt.subplots()
# customize ticks
ax.set_xticks(numpy.arange(0, len(LABELS_X)*FACTOR_X, FACTOR_X))
ax.set_yticks(numpy.arange(0, len(LABELS_Y)*FACTOR_Y, FACTOR_Y))
ax.set_xticklabels(LABELS_X)
ax.set_yticklabels(LABELS_Y)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
# save figure with colorbar
im = ax.imshow(free_energy)
cbar = fig.colorbar(im)
plt.savefig("FreeEnergy_heatmap.png")
plt.close()

# create heatmap with probability
fig, ax = plt.subplots()
# customize ticks
ax.set_xticks(numpy.arange(0, len(LABELS_X)*FACTOR_X, FACTOR_X))
ax.set_yticks(numpy.arange(0, len(LABELS_Y)*FACTOR_Y, FACTOR_Y))
ax.set_xticklabels(LABELS_X)
ax.set_yticklabels(LABELS_Y)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
# save figure with colorbar
im = ax.imshow(probability)
cbar = fig.colorbar(im)
plt.savefig("Probability_heatmap.png")
plt.close()
