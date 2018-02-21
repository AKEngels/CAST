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
        plt.close()
        
        return "Python here: All is wonderful!"
    except:
        print sys.exc_info()
        return "error"

def plot_temp(temperatures):
    try:
        plt.plot(temperatures)
        plt.xlabel("frame")
        plt.ylabel("temperature [K]")
        plt.savefig("temperature.png")
        plt.close()
        
        return "Python here: All is wonderful!"
    except:
        print sys.exc_info()
        return "error"

def plot_zones(legends, temp_lists):
    try:
        for t in temp_lists:  # for every zone
            plt.plot(t) 

        plt.xlabel("frame")
        plt.ylabel("temperature [K]")
        plt.legend(legends)
        plt.savefig("zones.png")
        plt.close()
        
        return "Python here: All is wonderful!"
    except:
        print sys.exc_info()
        return "error"

