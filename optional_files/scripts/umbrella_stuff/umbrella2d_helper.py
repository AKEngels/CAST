# list of parameters (in this case torsional angle 2)
PARAMETER_LIST = range(-180, 185, 5)

# user defined function to start a job, depending on one element from PARAMETER_LIST
def run_calc(param):
    
    # modules necessary for running the job
    import os
    import shutil
    
    # user variables for calculation
    PLACEHOLDER = "RESTRAINT_2"       # placeholder in CAST.txt file that is replaced by step number
    
    # create folder for current window
    os.mkdir("f_{}".format(param))

    # copy necessary files to that folder (USER INPUT)
    shutil.copy("CAST.txt", "f_{}/CAST.txt".format(param))
    shutil.copy("pentan.arc", "f_{}/pentan.arc".format(param))
    shutil.copy("charmm22.prm", "f_{}/charmm22.prm".format(param))
    shutil.copy("/home/susanne/CAST/optional_files/build/CAST_linux_x64_release",
                "f_{}/CAST.exe".format(param))
    
    # important: set correct parameter in inputfile
    with open("f_{}/CAST.txt".format(param)) as inp:
        x = inp.read()
        x = x.replace(PLACEHOLDER,str(float(param)))
    with open("f_{}/CAST.txt".format(param),"w") as inp:
        inp.write(x)
            
    # submit calculation
    os.chdir("f_{}".format(param))
    os.popen("chmod +x CAST")
    os.popen("./CAST.exe | tee CAST_OUTPUT.txt")
    os.chdir("..")


########################################### PROGRAM ########################################

for p in PARAMETER_LIST:
    run_calc(p)
             
