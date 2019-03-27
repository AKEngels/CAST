### script for graphical analysis of Umbrella Sampling ###

import glob
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

# read outputfile
with open("out.txt") as outfile:
    lines = outfile.readlines()

# get values for reaction coordinate, free energy and probability
xi = []  # reaction coordinate
A = []   # free energy
P = []   # probability
for line in lines:
    if line.startswith("#") == False:
        linelist = line.split()
        xi.append(float(linelist[0]))
        A.append(float(linelist[1]))
        P.append(float(linelist[3]))
        
# determine minimum and maximum value for xi
MIN, MAX = min(xi), max(xi)

# plot Free Energy
plt.plot(xi, A, linewidth=2, color='black')
plt.xlim(min(xi),max(xi))
plt.xlabel('reaction coordinate')
plt.ylabel('Free Energy A [kcal/mol]')
plt.savefig("FreeEnergy.png")
plt.close()

# plot probability
plt.plot(xi, P, linewidth=2, color='black')
plt.xlim(min(xi),max(xi))
plt.xlabel('reaction coordinate')
plt.ylabel('Probability P')
plt.savefig("Probability.png")
plt.close()

# plot distributions from sampling
print "Start ploting distributions"
xi_lists = []
broken_files = []
for filename in glob.glob("umbrella_*.txt"):
    print "looking at file", filename
    with open(filename) as umbrella_file:
        lines = umbrella_file.readlines()
    xi = []
    for i,line in enumerate(lines):
        xi.append(float(line.split()[1]))
        if line.split()[1] == "-nan" or line.split()[1] == "nan":
            if filename not in broken_files:
                broken_files.append(filename)
    number = int(math.sqrt(len(xi)))
    n, bins, patches = plt.hist(xi, number, histtype='step')
    plt.savefig("distribution_{}.png".format(filename))
    plt.close()
    xi_lists.append(xi)

# all distributions in one file
plt.rcParams['figure.figsize'] = [11,8]  # 11x8 inches is next to DinA4
number = int(math.sqrt(len(xi_lists[0])*len(xi_lists)))
plt.xlim(MIN, MAX)
n, bins, patches = plt.hist(xi_lists, number, histtype='step')
plt.savefig("distribution.png", dpi=100)
print broken_files

    
