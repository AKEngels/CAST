### script that takes atom types from a file 'atom_types_bonded.txt' and asigns them to the tinkerstructure 'output.arc'
### afterwards it changes the atom types of certain atoms
### improve: make a function where you give atom number and desired atom type that assigns that atom type to the atom

"""function to insert the atomtypes"""
def insert_atomtypes(lines, atomtypes):
    for i in range(81):
        lines[i+3087] = "{}{:3}{}".format(lines[i+3087][:50], atomtypes[i][:-1], lines[i+3087][53:])

"""change some atom types (e.g. for a proton transfer from inhibitor to His)"""
def proton_transfer(lines):
    lines[3160] = "{}{:3}{}".format(lines[3160][:50], 454, lines[3160][53:])
    lines[3118] = "{}{:3}{}".format(lines[3118][:50], 87, lines[3118][53:])
    lines[3161] = "{}{:3}{}".format(lines[3161][:50], 89, lines[3161][53:])

"""change some atom types (e.g. for breaking the bond between inhibitor and protein)"""
def break_bond(lines):
    lines[3117] = "{}{:3}{}".format(lines[3117][:50], 87, lines[3117][53:])
    lines[3159] = "{}{:3}{}".format(lines[3159][:50], 89, lines[3159][53:])

# open tinkerstructure
with open("output.arc") as arcfile:
    lines = arcfile.readlines()

# open atom types (every line the number of another atom type)
with open("atom_types_bonded.txt") as atfile:
    atomtypes_bonded = atfile.readlines()

# do stuff
insert_atomtypes(lines, atomtypes_bonded)
proton_transfer(lines)
break_bond(lines)

# write new tinkerstucture
with open("new.arc","w") as new_arc:
    new_arc.writelines(lines)
