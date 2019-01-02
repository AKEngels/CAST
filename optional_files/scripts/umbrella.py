### script to start an umbrella calculation on ECPC, things users probably want to change are marked with USER INPUT ###

import os
import shutil

steps = range(-180, 185, 5)  # to which values should the restraint be set? (USER INPUT)

"""function that copies all necessary files for a calculation into folders and changes lines in inputfiles"""
def run_calc(step):

    # create folder for current window
    os.mkdir("{}".format(step))  

    # copy necessary files to that folder (USER INPUT)
    shutil.copy("butan.arc","{}/butan.arc".format(step))
    shutil.copy("oplsaa_mod2.prm","{}/oplsaa_mod2.prm".format(step))
    shutil.copy("cast.sh", "{}/cast.sh".format(step))
    shutil.copy("CAST.txt", "{}/CAST.txt".format(step))

    # optional: change jobname
    with open("{}/cast.sh".format(step)) as inp:
        x = inp.read()
        x = x.replace("TITEL",str(step))
    with open("{}/cast.sh".format(step),"w") as inp:
        inp.write(x)
    
    # important: set correct restraint in inputfile
    with open("{}/CAST.txt".format(step)) as inp:
        x = inp.read()
        x = x.replace("UMBRELLA_RESTRAINT",str(float(step)))
    with open("{}/CAST.txt".format(step),"w") as inp:
        inp.write(x)

    # submit calculation
    os.chdir(str(step))
    os.popen("qsub cast.sh")
    os.chdir("..")

for s in steps:
    run_calc(s)
