import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_dists(legends, distances):
    try:
        for d in distances: # for every atom pair
            plt.plot(d) # plot the distances

        plt.xlabel("frame")
        plt.ylabel("distance [$\AA$]")
        plt.legend(legends)
        plt.savefig("distances.png")
        
        return "Python here: All is wonderful!"
    except:
        print sys.exc_info()
        return "error"

