import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt    
    
def plot_histograms_and_calculate_overlap(dE_pots, dE_pot_backs, window):
    try:
        number = int(math.sqrt(len(dE_pots)))
        n, bins, patches = plt.hist([dE_pots, dE_pot_backs], number, histtype='step', color = ["blue", "red"], label = ["E_pot", "E_pot,back"])
        plt.legend()
        plt.savefig("{}.png".format(window))
        plt.close()
        minimal_area = []  # calculate overlap of distributions
        for j in range(number):
            minimal_area.append(min(n[0][j], n[1][j]))
        overlap = sum(minimal_area)/len(dE_pots)
        return str(overlap)
    except:
        return "error"

