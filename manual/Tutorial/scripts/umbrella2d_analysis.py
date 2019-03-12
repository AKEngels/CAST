### script to analyse a 2D umbrella calculation ###

import os
import shutil

# USER INPUT
PATH = "/apps/wham/wham-2d/wham-2d"   # path to WHAM program

WINDOWS_X = range(-180, 185, 5)       # values for first restraint (X)
FORCE_CONSTANT_X = 0.05               # force constant of biasing potential for restraint X
WINDOWS_Y = range(-180, 185, 5)       # values for second restraint (Y)
FORCE_CONSTANT_Y = 0.05               # force constant of biasing potential for restraint Y

# ... for WHAM
MIN_X = -182.5        # boundaries of historgram for restraint X
MAX_X = 182.5
BINS_X = 73           # number of points in final PMF (= number of bins) for restraint X

MIN_Y = -182.5        # boundaries of historgram for restraint Y
MAX_Y = 182.5
BINS_Y = 73           # number of points in final PMF (= number of bins) for restraint Y

TEMP = 300          # temperature for WHAM 
TOL = 0.00001       # tolerance for WHAM interations
AUTO_MASK = False   # should auto-mask feature be turned on?


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


# creates a folder 'analysis' and copies all 'umbrella.txt' files there
# to the names of the files the step is added
# furthermore a file 'in.txt' is created that is used as input for WHAM
# if files are missing this function returns false, otherwise true
def copy_and_and_create_inputfile():
    success = True
    os.mkdir("analysis") # create folder
    for s1 in WINDOWS_X:
        for s2 in WINDOWS_Y:
            print "looking at", s1, s2
            if os.path.isfile("f_{}/f_{}/umbrella.txt".format(s1,s2)):
                # copy and rename file
                shutil.copy("f_{}/f_{}/umbrella.txt".format(s1,s2),
                            "analysis/umbrella_{}_{}.txt".format(s1,s2))
                # write inputfile for WHAM
                with open("analysis/in.txt","a") as inp:         
                    inp.write("umbrella_{}_{}.txt   {}  {}   {}  {}\n".format(
                        s1,s2,float(s1), float(s2),FORCE_CONSTANT_X, FORCE_CONSTANT_Y))
            else:
                print "No file 'umbrella.txt' in folder", s
                sucess = False
    return success


if copy_and_and_create_inputfile():
    os.chdir("analysis")
    order = ""
    if AUTO_MASK == False:
        order = "{} Px=0 {} {} {} Py=0 {} {} {} {:.9f} {} 0 in.txt out.txt 0 > wham_output.txt".format(
            PATH, MIN_X, MAX_X, BINS_X, MIN_Y, MAX_Y, BINS_Y, TOL, TEMP)
    else:
        order = "{} Px=0 {} {} {} Py=0 {} {} {} {:.9f} {} 0 in.txt out.txt 1 > wham_output.txt".format(
            PATH, MIN_X, MAX_Y, BINS_X, MIN_Y, MAX_Y, BINS_Y, TOL, TEMP)
    submit_pbs(order)
