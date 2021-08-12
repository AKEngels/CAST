# Execute this file in the CAST calculation folder to plot 1D Histogram of PCA-Modes
# exec(open("/home/dustin/git/CAST/optional_files/scripts/pca_plot_1D.py").read())
#
#
import math
filestring = "pca_modes.dat" # 
targetModes = [1,2,3,4,5]
timeStepFS = 1
numTimeSteps = None
binsNumber = 30
#
#
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import itertools
import copy
#
inFile = open(filestring,"r+")
modes = []
found = False
for line in inFile:
    if line.find("Trajectory in PCA - Modes following") is not -1:
        found = True
        continue
    if found is True:
        if line.find("Additional Information following") is not -1:
            found = False
            continue
        modes.append(copy.deepcopy(line.split()))
print("IO done: " + str(modes[0]))
print("Commencing plotting...")
modes = modes[1:]
for i in range(len(modes)):
    for j in range(len(modes[i])):
        modes[i][j] = copy.deepcopy(float(modes[i][j])) #* math.sqrt(10)
for curMode in range(len(targetModes)):
    for curMode2 in range(len(targetModes)):
        if curMode2 > curMode:
            print("Plotting 2D PCA of Modes ",targetModes[curMode]," and ",targetModes[curMode2])
            fig, ax = plt.subplots()
            plt.hist2d(modes[targetModes[curMode]-1], modes[targetModes[curMode2]-1], bins=(binsNumber, binsNumber), cmap=plt.cm.YlGn)
            plt.colorbar()
            #plt.hist(modes[curMode], density=True, bins=binsNumber)
            outstring = "pca_2d_histo_values_mode_" + str(targetModes[curMode]) + "_" + str(targetModes[curMode2])
            plt.savefig(str(outstring +'.png'), transparent=True,bbox_inches='tight')
            plt.savefig(str(outstring +'.pgf'), transparent=True,bbox_inches='tight')
            plt.savefig(str(outstring +'.pdf'), transparent=True,bbox_inches='tight')
            plt.clf()
