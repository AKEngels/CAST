### script to analyse a 2D umbrella calculation ###

import os
import shutil

# USER INPUT
PATH = "/apps/wham/wham-2d/wham-2d"   # path to WHAM program

WINDOWS_X = range(-180, 185, 5)       # values for first restraint
FORCE_CONSTANT_X = 0.05               # force constant of biasing potential for restraint 1
WINDOWS_Y = range(-180, 185, 5)       # values for second restraint
FORCE_CONSTANT_Y = 0.05               # force constant of biasing potential for restraint 2

# ... for WHAM
MIN_X = -182.5        # boundaries of historgram for restraint 1
MAX_X = 182.5
BINS_X = 73           # number of points in final PMF (= number of bins) for restraint 1

MIN_Y = -182.5        # boundaries of historgram for restraint 2
MAX_Y = 182.5
BINS_Y = 73           # number of points in final PMF (= number of bins) for restraint 2

TEMP = 300          # temperature for WHAM 
TOL = 0.00001       # tolerance for WHAM interations
AUTO_MASK = True    # should auto-mask feature be turned on?

COMBINATIONS = [[-180, -180], [-180, -175], [-175, -175], [-170, -175], [-170, -170], [-165, -170],
                [-165, -165], [-160, -165], [-160, -160], [-155, -160], [-155, -155], [-155, -150],
                [-150, -150], [-150, -145], [-150, -140], [-150, -135], [-150, -130], [-150, -125],
                [-150, -120], [-150, -115], [-150, -110], [-150, -105], [-150, -100], [-150, -95],
                [-150, -90], [-150, -85], [-150, -80], [-150, -75], [-150, -70], [-150, -65], [-150, -60],
                [-150, -55], [-150, -50], [-150, -45], [-145, -45], [-140, -45], [-135, -45], [-130, -45],
                [-125, -45], [-120, -45], [-115, -45], [-110, -45], [-105, -45], [-100, -45], [-95, -45],
                [-90, -45], [-85, -45], [-80, -45], [-75, -45], [-70, -45], [-65, -45], [-60, -45],
                [-55, -45], [-50, -45], [-45, -45], [-40, -45], [-40, -40], [-40, -35], [-40, -30],
                [-40, -25], [-40, -20], [-40, -15], [-40, -10], [-40, -5], [-40, 0], [-40, 5], [-40, 10],
                [-40, 15], [-40, 20], [-40, 25], [-40, 30], [-40, 35], [-40, 40], [-40, 45], [-40, 50],
                [-40, 55], [-40, 60], [-40, 65], [-40, 70], [-40, 75], [-40, 80], [-40, 85], [-40, 90],
                [-40, 95], [-40, 100], [-40, 105], [-40, 110], [-40, 115], [-40, 120], [-40, 125],
                [-40, 130], [-40, 135], [-40, 140], [-40, 145], [-40, 150], [-40, 155], [-40, 160],
                [-40, 165], [-40, 170], [-40, 175], [-40, 180], [-35, 180], [-30, 180], [-25, 180],
                [-20, 180], [-15, 180], [-10, 180], [-5, 180], [0, 180], [5, 180], [10, 180], [15, 180],
                [20, 180], [25, 180], [30, 180], [35, 180], [40, 180], [45, 180], [50, 180], [55, 180],
                [60, 180], [65, 180], [70, 180], [75, 180], [80, 180], [85, 180], [90, 180], [95, 180],
                [100, 180], [105, 180], [110, 180], [115, 180], [120, 180], [125, 180], [130, 180],
                [135, 180], [140, 180], [145, 180], [150, 180], [155, 180], [160, 180], [165, 180],
                [170, 180], [175, 180], [180, 180]]


# submits a calculation
# order is a string which should be executed
def submit_pbs(order):

    # parameters for PBS
    NUM_PROCESSORS = 4  # number of processors
    MEMORY = 4024       # working memory in MB
    WALLTIME = 999      # walltime in hours
    TITEL = "wham2d"    # jobtitle

    with open("pbs.sh","w") as pbs:
        pbs.write("#PBS -l nodes=1:ppn={}\n".format(NUM_PROCESSORS))
        pbs.write("#PBS -l mem={}MB\n".format(MEMORY))
        pbs.write("#PBS -l walltime={}:00:00\n".format(WALLTIME))
        pbs.write("#PBS -N {}\n\n".format(TITEL))

        pbs.write("cd $PBS_O_WORKDIR\n")
        pbs.write("cp -pr * $TMPDIR\n")
        pbs.write("cd $TMPDIR\n\n")

        pbs.write(order + "\n\n")
        pbs.write("cp -pr * $PBS_O_WORKDIR\n")

    os.popen("qsub -V pbs.sh")


# creates wham inputfile
def create_inputfile():
    for s1 in WINDOWS_X:
        for s2 in WINDOWS_Y:
            if [s1,s2] in COMBINATIONS:
                with open("in.txt","a") as inp:         
                    inp.write("umbrella_{}_{}.txt   {}  {}   {}  {}\n".format(
                        s1,s2,float(s1), float(s2),FORCE_CONSTANT_X, FORCE_CONSTANT_X))


create_inputfile()
order = ""
if AUTO_MASK == False:
    order = "{} Px=0 {} {} {} Py=0 {} {} {} {:.9f} {} 0 in.txt out.txt 0 > wham_out.txt".format(
        PATH, MIN_X, MAX_X, BINS_X, MIN_Y, MAX_Y, BINS_Y, TOL, TEMP)
else:
    order = "{} Px=0 {} {} {} Py=0 {} {} {} {:.9f} {} 0 in.txt out.txt 1 > wham_out.txt".format(
        PATH, MIN_X, MAX_X, BINS_X, MIN_Y, MAX_Y, BINS_Y, TOL, TEMP)
submit_pbs(order)
