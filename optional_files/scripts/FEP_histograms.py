"""run this script on the output file 'dE_pot.txt' of an FEP calculation
it plots the probability distributions P_i(dE_pot) und P_f(dE_pot) for every window and saves them as .png images
furhermore it calcultes the overlap of the forward and backwards distributions of every window which is a quality estimator"""

import math
import matplotlib.pyplot as plt

with open("dE_pot.txt") as dEpot:
    lines = dEpot.readlines()

pot_energies = []  # list of potential energies [E_pot(0), E_pot,back(dlambda), E_pot(dlambda), ..., E_pot(1-dlambda), E_pot,back(1)
for i,l in enumerate(lines):
    if (i%4 == 2 or i%4 == 3) and i != 2 and i != len(lines)-1:
        splitted_line = l.split(",")[1:-1]
        number_line = []
        for s in splitted_line:
            number = float(s)
            number_line.append(number)
        pot_energies.append(number_line)

value_number = len(pot_energies[0])   # number of single values in every histogram (= number of production steps in FEP)
number = int(math.sqrt(value_number)) # number of steps for histogram (change for finer or rougher diagram)

for i in range(len(pot_energies)/2):  # for every window plot distributions E_pot(lambda) und E_pot,back(lambda+dlambda)
    
    n, bins, patches = plt.hist([pot_energies[2*i], pot_energies[2*i+1]], number, histtype='step', color = ["blue", "red"], label = ["E_pot", "E_pot,back"])
    plt.legend()
    plt.savefig("{}.png".format(i+1))
    plt.close()

    minimal_area = []  # calculate overlap of distributions
    for j in range(number):
        minimal_area.append(min(n[0][j], n[1][j]))
    overlap = sum(minimal_area)/value_number
    print "Fenster {:2}: {:.1f} %".format(i+1, overlap*100) 
    
    
