############################## DESCRIPTION ##############################################

# This is a script with which you can get a better control over the queue on UNICORN.
# You can give a number of jobs which you want to at maximum send to the queue and a time <t>.
# Every <t> seconds the script looks at qstat how many jobs you have currently running.
# It determines how many jobs can now be sent and submits them.

# You have to give a list of parameters with which the calculation is sent.
# Furthermore you have to give a function run_calc(param) which runs the calculation.
# For a better understanding you can imagine this program to just take the parameter list
# and run the function for every parameter in there (just in a more clever timing).


####################################### STUFF DEFINED BY USER ############################

MAX_JOBS = 25      # this is the maximum number of jobs I want to send to queue
DELAY = 300        # time in seconds for which the program waits before looking again

# list where every element is the parameter of one calculation
# elements of list can be anything, like numbers, lists or objects of a user-defined class
dim_1 = range(-180, 185, 5)
dim_2 = range(-180, 185, 5)
PARAMETER_LIST = []
for d1 in dim_1:
    for d2 in dim_2:
        PARAMETER_LIST.append([d1, d2])


# user defined function to start a job, depending on one element from PARAMETER_LIST
# this function is just an example (it performs an umbrella sampling with CAST)
# the user might delete the content of this function and write it new
def run_calc(param):
    
    # modules necessary for running the job
    import os
    import shutil
    
    # user variables for calculation
    MOLECULE = "pentan.arc"             # name of molecule file
    FORCEFIELD = "charmm22.prm"         # name of forcefield file
    
    PLACEHOLDER_1 = "RESTRAINT_1"       # placeholders in CAST.txt file that is replaced by step number
    PLACEHOLDER_2 = "RESTRAINT_2"
    
    # create folder for current window
    if os.path.isdir("f_{}".format(str(param[0]))) == False:
        os.mkdir("f_{}".format(str(param[0])))
    os.mkdir("f_{}/f_{}".format(param[0], param[1]))  

    # copy necessary files to that folder (USER INPUT)
    shutil.copy(MOLECULE,"f_{}/f_{}/{}".format(param[0], param[1],MOLECULE))
    shutil.copy(FORCEFIELD,"f_{}/f_{}/{}".format(param[0], param[1],FORCEFIELD))
    shutil.copy("cast.sh", "f_{}/f_{}/cast.sh".format(param[0], param[1]))
    shutil.copy("CAST.txt", "f_{}/f_{}/CAST.txt".format(param[0], param[1]))

    # optional: change jobname
    with open("f_{}/f_{}/cast.sh".format(param[0], param[1])) as inp:
        x = inp.read()
        x = x.replace("TITEL","calc{}_{}".format(param[0],param[1]))
    with open("f_{}/f_{}/cast.sh".format(param[0], param[1]),"w") as inp:
        inp.write(x)
    
    # important: set correct parameter in inputfile
    with open("f_{}/f_{}/CAST.txt".format(param[0], param[1])) as inp:
        x = inp.read()
        x = x.replace(PLACEHOLDER_1,str(float(param[0])))
        x = x.replace(PLACEHOLDER_2,str(float(param[1])))
    with open("f_{}/f_{}/CAST.txt".format(param[0], param[1]),"w") as inp:
        inp.write(x)

    # submit calculation
    os.chdir("f_{}".format(str(param[0])))
    os.chdir("f_{}".format(str(param[1])))
    os.popen("qsub -V cast.sh")
    os.chdir("..")
    os.chdir("..")


########################################### PROGRAM ########################################

import subprocess
import sys
import time

# function that executes qstat and counts how many jobs are already there
def look_for_running_jobs():
    qstat_str = subprocess.Popen(['qstat'], stdout=subprocess.PIPE).stdout.read()
    qstat_list = qstat_str.split('\n')
    if len(qstat_list) < 3:     # no calculations running
        return 0
    else:                       # subtract two lines for title and one empty line at the end
        return len(qstat_list) - 3
    
# function that submits a certain number of jobs
# starting_counter: index of the first calculation to be sent
# jobs_to_be_sent: total number of calculations that should be sent now
# returns the index of the next calculation that must be sent afterwards
def start_next_jobs(starting_counter, jobs_to_be_sent):
    counter = starting_counter
    for i in range(starting_counter, jobs_to_be_sent+starting_counter):
        print "Starting job", counter+1, "of", len(PARAMETER_LIST), ":", PARAMETER_LIST[counter]
        run_calc(PARAMETER_LIST[counter])
        counter += 1
        if counter == len(PARAMETER_LIST):
            sys.exit(0)
    return counter

# counter of the jobs that are already sent
counter = 0

# see how many jobs are running and how much can be sent now
current_jobs = look_for_running_jobs()
jobs_to_be_sent = MAX_JOBS - current_jobs
if jobs_to_be_sent < 0:
    raise ValueError("There are more jobs running than are allowed to be sent.")

# start the first jobs
print time.ctime(), ":", current_jobs, "jobs running,", jobs_to_be_sent, "will be sent"
counter = start_next_jobs(counter, jobs_to_be_sent)

# while there are still jobs left:
while True:
    time.sleep(DELAY)
    current_jobs = look_for_running_jobs()
    jobs_to_be_sent = MAX_JOBS - current_jobs
    if jobs_to_be_sent < 0:
        print time.ctime(), ":", current_jobs, "jobs running, no jobs will be sent"
    else:
        print time.ctime(), ":", current_jobs, "jobs running,", jobs_to_be_sent, "will be sent"
        counter = start_next_jobs(counter, jobs_to_be_sent)
