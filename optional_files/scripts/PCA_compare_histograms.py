### This script is for plotting overlap of two PCA distributions. 
### It reads the files 'pca_modes_1.dat' and 'pca_modes_2.dat' from the CAST task PCAgen with 2 trajectories.
### You can choose which mode should be plotted (starting with 1).
### It marks the value of the first step which should be about the same in the two directories.

MODE = 1

import math
import numpy
import matplotlib.pyplot as plt
    
def plot_histograms_and_calculate_overlap(data_1, data_2, filename):
    number = int(math.sqrt(len(data_1)))
    n, bins, patches = plt.hist([data_1, data_2], number, histtype='step', color = ["blue", "red"], label = ["full", "new"])
    plt.legend()
    plt.plot([data_1[0], data_1[0]], [0, max(numpy.concatenate(n))], linestyle="--")
    plt.savefig(filename)
    plt.close()
    minimal_area = []  # calculate overlap of distributions
    for j in range(number):
        minimal_area.append(min(n[0][j], n[1][j]))
    overlap = sum(minimal_area)/len(data_1)
    return overlap

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

linenumber = find_line(MODE, lines_1)
values_1 = get_values_from_line(lines_1[linenumber])
values_2 = get_values_from_line(lines_2[linenumber])
print plot_histograms_and_calculate_overlap(values_1, values_2, "mode_{}.png".format(MODE))
