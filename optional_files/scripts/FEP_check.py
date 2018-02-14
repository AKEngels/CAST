"""This is a script to calculate the free energy out of the file "alchemical.txt" that is an output file of FEP calculations in CAST
can be used to control if the results of the calculation of dG out of the sample averages is reasonable"""

import math

boltz = 1.3806488*10**(-23) # boltzman constant
avogad = 6.022*10**23       # avogadro constant
conv = 4184.0               # conversion factor kcal to joule

de_ensemble = 0
de_ensemble_back = 0
go = False

runs = [[]]
dG_list = []
dG_list_back = []

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


for r in runs:  # for every production run
    de_ens = []
    de_ens_back = []
    temps = []
    for k in r:  # for every conformation in this run
        energy = float(k[105:120])  # find dE
        current_temp = float(k[90:105]) # find current temperature
        energy_back = float(k[135:150])
        #print energy, current_temp
        de_ens.append(math.exp(-1 / (boltz*current_temp)*conv*energy / avogad)) # calculate the e-function term
        de_ens_back.append(math.exp(1 / (boltz*current_temp)*conv*energy_back / avogad)) # calculate the e-function term
        temps.append(current_temp) # save current temperature
        
    de_ensemble = sum(de_ens)/len(de_ens) # calculate ensemble average
    de_ensemble_back = sum(de_ens_back)/len(de_ens_back) # calculate ensemble average
    temp = sum(temps)/len(temps)          # calculate average temperature

    dG = -1 * math.log(de_ensemble)*temp*boltz*avogad / conv  # calculate dG
    dG_back = 1 * math.log(de_ensemble_back)*temp*boltz*avogad / conv  # calculate dG
    dG_list.append(dG)  # save dG
    dG_list_back.append(dG_back)  # save dG
    print dG, sum(dG_list), dG_back, sum(dG_list_back)  # print dG of current window and complete dG until the current window
