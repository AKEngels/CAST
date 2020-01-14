# Execute this file in the CAST calculation folder to plot 1D Histogram of PCA-Modes
# exec(open("/home/dustin/git/CAST/optional_files/scripts/pca_plot_1D.py").read())
#
#
filestring = "pca_modes.dat" # 
targetModes = [0,1,2,3,4,5]
timeStepFS = 1
numTimeSteps = None
#
#
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import itertools
import copy
#
for curMode in targetModes:
    with open(filestring) as inFile:
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
        modes = modes[1:]
        for i in range(len(modes)):
            for j in range(len(modes[i])):
                modes[i][j] = copy.deepcopy(float(modes[i][j]))
        fig, ax = plt.subplots()
        plt.hist(modes[curMode], density=True, bins=30)
        outstring = "pca_1d_histo_values_mode_" + str(curMode)
        plt.savefig(str(outstring +'.png'), transparent=True,bbox_inches='tight')
        plt.savefig(str(outstring +'.pgf'), transparent=True,bbox_inches='tight')
        plt.savefig(str(outstring +'.pdf'), transparent=True,bbox_inches='tight')
        plt.clf()
        if numTimeSteps is None:
            ranger = list(range(len(modes[curMode])))
            #print(ranger)
            for i in range(len(ranger)):
                ranger[i] *= timeStepFS
            plt.plot(ranger,modes[curMode],color='blue')
            outstring = "pca_1d_timeline_values_mode_" + str(curMode)
            plt.savefig(str(outstring +'.png'), transparent=True,bbox_inches='tight')
            plt.savefig(str(outstring +'.pgf'), transparent=True,bbox_inches='tight')
            plt.savefig(str(outstring +'.pdf'), transparent=True,bbox_inches='tight')
            plt.clf()
        else:
            for j in numTimeSteps:
                ranger = range(len(modes[curMode][:j]))
                for i in range(len(ranger)):
                    ranger[i] *= timeStepFS
                plt.plot(ranger,modes[curMode][:j],color='blue')
                outstring = "pca_1d_timeline_values_mode_" + str(curMode) + "_" + str(j*timeStepFS) + "fs" 
                plt.savefig(str(outstring +'.png'), transparent=True,bbox_inches='tight')
                plt.savefig(str(outstring +'.pgf'), transparent=True,bbox_inches='tight')
                plt.savefig(str(outstring +'.pdf'), transparent=True,bbox_inches='tight')
                plt.clf()
