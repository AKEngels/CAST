### This script is for plotting overlap of two PCA distributions in two dimensions. 
### It reads the files 'pca_modes_1.dat' and 'pca_modes_2.dat' from the CAST task PCAgen with 2 trajectories.
### You can choose which modes should be plotted (starting with 1).
### It marks the value of the first step which should be about the same in the two directories.

MODE = (1,2)

import math
import matplotlib.pyplot as plt

def find_line(mode, lines):
    for i, line in enumerate(lines):
        if line.startswith("Trajectory in PCA - Modes following (columns are frames, rows are modes)"):
            return i+mode+1

def get_values_from_line(line):
    linelist = line.split()
    valuelist = []
    for l in linelist:
        valuelist.append(float(l))
    return valuelist

with open("pca_modes_1.dat") as inp:
    lines_1 = inp.readlines()

with open("pca_modes_2.dat") as inp:
    lines_2 = inp.readlines()

linenumber_first = find_line(MODE[0], lines_1)
linenumber_second = find_line(MODE[1], lines_1)
values_first_1 = get_values_from_line(lines_1[linenumber_first])
values_second_1 = get_values_from_line(lines_1[linenumber_second])
plt.scatter(values_first_1,values_second_1,label="full",s=0.5,color="blue")

values_first_2 = get_values_from_line(lines_2[linenumber_first])
values_second_2 = get_values_from_line(lines_2[linenumber_second])
plt.scatter(values_first_2,values_second_2,label="new",s=0.5,color="red")

plt.xlim(min(values_first_1 + values_first_2), max(values_first_1 + values_first_2))
plt.ylim(min(values_second_1 + values_second_2), max(values_second_1 + values_second_2))
plt.scatter(values_first_1[0],values_second_1[0], label="starting point", color="green")

plt.legend()
plt.savefig("compare_{}_{}.png".format(MODE[0], MODE[1]))
plt.close()
