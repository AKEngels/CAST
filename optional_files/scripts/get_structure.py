### script to cut out a structure from a MD trajectory ###

SNAPSHOT_NUMBER = 12 # aus VMD

with open("output_MD_SNAP.arc") as trajfile:
    inputfile = trajfile.readlines()

NUMBER_OF_ATOMS = int(inputfile[0])
print NUMBER_OF_ATOMS

number_of_leftout_lines = (NUMBER_OF_ATOMS+1)*SNAPSHOT_NUMBER

with open("structure_{}.arc".format(SNAPSHOT_NUMBER), "w") as structure:
    structure.writelines(inputfile[number_of_leftout_lines:number_of_leftout_lines+NUMBER_OF_ATOMS+1])
