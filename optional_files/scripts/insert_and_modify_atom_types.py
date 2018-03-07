### script that takes atom types from a file 'atom_types_bonded.txt' and asigns them to the tinkerstructure 'output.arc'
### 'atom_types_bonded.txt' needs a line break after the last atom type but this last line must be empty (i.e. no spaces or such)
### afterwards it changes the atom types of certain atoms

"""function to insert the atomtypes
lines: inputfile read in with readlines()
atomtypes: file with atomtypes read in with readlines()
starting_index: number of line where first atomtype is to be inserted (tinker atom index)"""
def insert_atomtypes(lines, atomtypes, starting_index):
    for i in range(len(atomtypes)):
        lines[i+starting_index] = "{}{:3}{}".format(lines[i+starting_index][:50], atomtypes[i][:-1], lines[i+starting_index][53:])

"""function to change an atom type
lines: inputfile read in with readlines()
changeline: line of the atom whose type is to be changed (tinker atom index)
new_atomtype: desired atom type"""
def change_atomtype(lines, changeline, new_atomtype):
    lines[changeline] = "{}{:3}{}".format(lines[changeline][:50], new_atomtype, lines[changeline][53:])


"""change some atom types (e.g. for a proton transfer from inhibitor to His)"""
def proton_transfer(lines):
    change_atomtype(lines, 3160, 454)
    change_atomtype(lines, 3118, 87)
    change_atomtype(lines, 3161, 89)

"""change some atom types (e.g. for breaking the bond between inhibitor and protein)"""
def break_bond(lines):
    change_atomtype(lines, 3117, 87)
    change_atomtype(lines, 3159, 89)
    

# open tinkerstructure
with open("output.arc") as arcfile:
    lines = arcfile.readlines()

# open atom types (every line the number of another atom type)
with open("atom_types_bonded.txt") as atfile:
    atomtypes_bonded = atfile.readlines()

# do stuff
insert_atomtypes(lines, atomtypes_bonded, 3088)
change_atomtype(lines, 2277, 85)
#proton_transfer(lines)
#break_bond(lines)

# write new tinkerstucture
with open("new.arc","w") as new_arc:
    new_arc.writelines(lines)
