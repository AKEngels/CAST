### script to start an umbrella calculation on UNICORN ###

import os
import shutil

MOLECULE = "butan.arc"              # name of molecule file
FORCEFIELD = "charmm22.prm"         # name of forcefield file
STEPS = range(-180, 185, 5)         # to which values should the restraint be set?
PLACEHOLDER = "UMBRELLA_RESTRAINT"  # placeholder in CAST.txt file that is replaced by step number

"""function that copies all necessary files for a calculation into folders and changes lines in inputfiles"""
def run_calc(step):

    # create folder for current window
    os.mkdir("f_{}".format(step))  

    # copy necessary files to that folder (USER INPUT)
    shutil.copy(MOLECULE,"f_{}/{}".format(step,MOLECULE))
    shutil.copy(FORCEFIELD,"f_{}/{}".format(step,FORCEFIELD))
    shutil.copy("cast.sh", "f_{}/cast.sh".format(step))
    shutil.copy("CAST.txt", "f_{}/CAST.txt".format(step))

    # optional: change jobname
    with open("f_{}/cast.sh".format(step)) as inp:
        x = inp.read()
        x = x.replace("TITEL",str(step))
    with open("f_{}/cast.sh".format(step),"w") as inp:
        inp.write(x)
    
    # important: set correct restraint in inputfile
    with open("f_{}/CAST.txt".format(step)) as inp:
        x = inp.read()
        x = x.replace(PLACEHOLDER,str(float(step)))
    with open("f_{}/CAST.txt".format(step),"w") as inp:
        inp.write(x)

    # submit calculation
    os.chdir("f_{}".format(step))
    os.popen("qsub -V cast.sh")
    os.chdir("..")

for s in STEPS:
    run_calc(s)
