# Execute this file in the CAST calculation folder to plot the temperature of the simulation.
# exec(open("plotThermostat.py").read())
#
#
print("Begin plotting")
import sys
alpha = 0.01 # For exponentially weighted avergae, best not to change this
targetTemp = 300 #  Set the target temperature for the plot of temperature deviation
targTempStartFrame = 5000 # Set the STARTING FRAME of constant temperature (the frame when heating is done and temp is maintained constant)
filestring = "OUT" # 
if len(sys.argv) == 2:
    filestring = sys.argv[1]
castoptions = "CAST.txt"
#
#
import matplotlib.pyplot as plt
import itertools
import copy
instaTemps = []
#
with open(castoptions) as inFile:
    tempOptions = []
    for line in inFile:
        if line.find("MDheat") != -1:
            #print(line.split())
            tempOptions.append([line.split()[1],line.split()[2]])
    highestTemp = 0
    relatedFrame = 0
    for i in range(len(tempOptions)):
        curNum = int(tempOptions[i][1])
        if curNum > highestTemp:
            highestTemp = curNum
            relatedFrame = int(tempOptions[i][0])
    targetTemp = copy.deepcopy(highestTemp)
    targTempStartFrame = copy.deepcopy(relatedFrame)

#
with open(filestring) as inFile:
    for line in inFile:
        if filestring.find("TRACE") == -1: # File is CAST std::cout with hopefully verbosity high enough (>4) that thermostatting output is there
            if line.find("THERMOCONTROL - Full step:") is not -1:
                words = line.split()
                currentTemp = words[9]
                instaTemps.append(float(currentTemp[:-2]))
        else: #  Assuming we got a TRACE.csv file as input
            if line.find('It') == -1:
                words = line.split(',')
                currentTemp = words[1]
                instaTemps.append(float(currentTemp))
    fig, ax = plt.subplots()
    plt.plot(list(range(len(instaTemps))),instaTemps,color='blue')
    exponentialWeightedAvg = [0 for i in range(len(instaTemps))] # https://dsp.stackexchange.com/questions/20333/how-to-implement-a-moving-average-in-c-without-a-buffer
    for num, name in enumerate(instaTemps):
        if (num is 0):
            exponentialWeightedAvg[num] = copy.deepcopy(instaTemps[num])
        else:
            exponentialWeightedAvg[num] = alpha * instaTemps[num] + (1.0 - alpha) * exponentialWeightedAvg[num - 1]
    diff = [0 for i in range(len(instaTemps))]
    for num, name in enumerate(instaTemps):
        if (num > targTempStartFrame):
            diff[num] = copy.deepcopy(exponentialWeightedAvg[num]) - targetTemp
    plt.plot(list(range(len(instaTemps))),exponentialWeightedAvg,color='red')
    plt.savefig(str("thermostat_values" +'.png'), transparent=True,bbox_inches='tight')
    plt.savefig(str("thermostat_values" +'.pgf'), transparent=True,bbox_inches='tight')
    plt.savefig(str("thermostat_values" +'.pdf'), transparent=True,bbox_inches='tight')
    plt.clf()
    plt.plot(list(range(len(instaTemps)))[-targTempStartFrame:],diff[-targTempStartFrame:],color='red')
    plt.savefig(str("thermostat_values_diff" +'.png'), transparent=True,bbox_inches='tight')
    plt.savefig(str("thermostat_values_diff" +'.pgf'), transparent=True,bbox_inches='tight')
    plt.savefig(str("thermostat_values_diff" +'.pdf'), transparent=True,bbox_inches='tight')
    plt.clf()
print("Plotting Done")