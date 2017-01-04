"""This is a script to calculate the free energy out of the file "alchemical.txt" that is an output file of FEP calculations in CAST
should be used to control if the results of the calculation of dG out of the sample averages is reasonable"""

import math

boltz = 1.3806488*10**(-23) # boltzman constant
avogad = 6.022*10**23       # avogadro constant
conv = 4184.0               # conversion factor kcal to joule

de_ensemble = 0
go = False

runs = [[]]
dG_list = []

with open("alchemical.txt") as alch:   # read lines of this file that describe production run
    for line in alch:
        if line.startswith("Starting new data collection with values: "):
            go = True
        elif line.startswith("Free energy change for the current window: "):
            go = False
            runs.append([])
        elif go == True:
            runs[-1].append(line)
runs.pop()  # delete last (empty) element of this list

temp = float(runs[0][0][60:75])  # get temperature

for r in runs:  # for every production run
    de_ens = []
    for k in r:  # for every conformation in this run
        energy = float(k[75:90])  # find dE
        de_ens.append(math.exp(-1 / (boltz*temp)*conv*energy / avogad)) # calculate the e-function term
        
    de_ensemble = sum(de_ens)/len(de_ens) # calculate ensemble average

    dG = -1 * math.log(de_ensemble)*temp*boltz*avogad / conv  # calculate dG
    dG_list.append(dG)  # save dG
    print dG, sum(dG_list)  # print dG of current window and complete dG until the current window
