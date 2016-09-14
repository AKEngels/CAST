CURRENT_NUMBER_OF_TINKERTYPES = 50
CURRENT_NUMBER_OF_ATOMTYPES = 3000

"""
This script can be used to convert an amber type parameter file (like gaff.dat) to a tinker type parameter file (like amber99.prm)

You have to copy the atoms from gaff.dat to a new file atoms.txt. Then comment out everything after "print ausgabe" and run the script.
Copy the output into the atom section of amber99.prm.

Then comment out "print ausgabe" and activate "print ausgabe_charge". Copy the output into the charge section of amber99.prm.

Then copy the vdw parameters from gaff.dat into a new file vdw.txt. Comment out "print ausgabe_charge", activate everything until "print ausgabe"
in the loop over vdwlist and copy the output into the vdw section of amber99.prm.

Then you can comment out everything you have activated in the step before (i.e. the conversion of the vdw parameters) and activate everything
until the next "print ausgabe". Copy the bonds from gaff.dat into a new file bindungen.txt, run the script and copy the output into the bond
section of amber99.prm.

Repeat the same for the angles (winkel.txt), torsions (torsions.txt) and improper torsions (impropers.txt).
"""

with open("atoms.txt") as atom_file:
    atoms = atom_file.readlines()

atomlist = []
for i, line in enumerate(atoms):
    atomlist.append(line[:2])
    ausgabe = 'atom{:11d}{:5d}    {}   "???????????????????????????"  ??{:10.3f}    ?  !! GAFF'.format(i+CURRENT_NUMBER_OF_ATOMTYPES, i+CURRENT_NUMBER_OF_TINKERTYPES, line[:2], float(line[3:8]))
    ausgabe_charge = "charge{:9d}              ??.????".format(i+CURRENT_NUMBER_OF_ATOMTYPES)
    #print ausgabe
    #print ausgabe_charge

endindex = len(atomlist) + CURRENT_NUMBER_OF_TINKERTYPES
    
with open("vdw.txt") as vdw_file:
    vdw = vdw_file.readlines()

vdwlist = []
for i, line in enumerate(vdw):
    x = atomlist.index(line[2:4])
    vdwlist.append([x, line])

vdwlist.sort()
for v in vdwlist:
    line = v[1]
    i = v[0]
    ausgabe = "vdw{:12d}{:21.4f}{:11.4f}   !! GAFF".format(i+CURRENT_NUMBER_OF_TINKERTYPES, float(line[14:20]), float(line[22:30]))
    #print ausgabe

with open("bindungen.txt") as bond_file:
    bonds = bond_file.readlines()

for i, line in enumerate(bonds):
    b1 = line[:2]
    b2 = line[3:5]
    z1 = atomlist.index(b1) + CURRENT_NUMBER_OF_TINKERTYPES
    z2 = atomlist.index(b2) + CURRENT_NUMBER_OF_TINKERTYPES
    ausgabe = "bond{:11d}{:5d}{:16.2f}{:11.4f}  !! GAFF".format(z1, z2, float(line[6:13]), float(line[15:22]))
    #print ausgabe

with open("winkel.txt") as winkel_file:
    winkel = winkel_file.readlines()

for i, line in enumerate(winkel):
    b1 = line[:2]
    b2 = line[3:5]
    b3 = line[6:8]
    z1 = atomlist.index(b1) + CURRENT_NUMBER_OF_TINKERTYPES
    z2 = atomlist.index(b2) + CURRENT_NUMBER_OF_TINKERTYPES
    z3 = atomlist.index(b3) + CURRENT_NUMBER_OF_TINKERTYPES
    ausgabe = "angle{:10d}{:5d}{:5d}{:11.2f}{:11.2f}  !! GAFF".format(z1, z2, z3, float(line[10:17]), float(line[22:29]))
    #print ausgabe

with open("torsions.txt") as torsion_file:
    torsions = torsion_file.readlines()

for i, line in enumerate(torsions):
    b1 = line[:2]
    b2 = line[3:5]
    b3 = line[6:8]
    b4 = line[9:11]
    z2 = atomlist.index(b2) + CURRENT_NUMBER_OF_TINKERTYPES
    if b3 != "cb":
        z3 = atomlist.index(b3) + CURRENT_NUMBER_OF_TINKERTYPES
        if b1 == "X " and b4 == "X ":
            z1_list = range(CURRENT_NUMBER_OF_TINKERTYPES, endindex)
            z4_list = range(CURRENT_NUMBER_OF_TINKERTYPES, endindex)
            for z1 in z1_list:
                for z4 in z4_list:
                    ausgabe = "torsion{:8d}{:5d}{:5d}{:5d}{:17.3f}{:7.1f}{:3d}".format(z1, z2, z3, z4, float(line[17:24]), float(line[31:38]), int(float(line[48:54])))
                    print ausgabe
        else:
            z1 = atomlist.index(b1) + CURRENT_NUMBER_OF_TINKERTYPES
            z4 = atomlist.index(b4) + CURRENT_NUMBER_OF_TINKERTYPES
            ausgabe = "torsion{:8d}{:5d}{:5d}{:5d}{:17.3f}{:7.1f}{:3d}".format(z1, z2, z3, z4, float(line[17:24]), float(line[31:38]), int(float(line[48:54])))
            print ausgabe

with open("impropers.txt") as improper_file:
    impropers = improper_file.readlines()

for i, line in enumerate(impropers):
    b1 = line[:2]
    b2 = line[3:5]
    b3 = line[6:8]
    b4 = line[9:11]
    z4 = atomlist.index(b4) + CURRENT_NUMBER_OF_TINKERTYPES
    z3 = atomlist.index(b3) + CURRENT_NUMBER_OF_TINKERTYPES
    if b1 == "X " and b2 == "X ":
        z1_list = range(CURRENT_NUMBER_OF_TINKERTYPES, endindex)
        z2_list = range(CURRENT_NUMBER_OF_TINKERTYPES, endindex)
        for z1 in z1_list:
            for z2 in z2_list:
                ausgabe = "imptors{:8d}{:5d}{:5d}{:5d}{:17.3f}{:7.1f}{:3d}  !! GAFF".format(z1, z2, z3, z4, float(line[20:24]), float(line[33:38]), int(float(line[47:CURRENT_NUMBER_OF_TINKERTYPES])))
                print ausgabe
    elif b1 == "X ":
        z2 = atomlist.index(b2)
        z1_list = range(CURRENT_NUMBER_OF_TINKERTYPES, endindex)
        for z1 in z1_list:
            ausgabe = "imptors{:8d}{:5d}{:5d}{:5d}{:17.3f}{:7.1f}{:3d}  !! GAFF".format(z1, z2, z3, z4, float(line[20:24]), float(line[33:38]), int(float(line[47:CURRENT_NUMBER_OF_TINKERTYPES])))
            print ausgabe
    else:
        z1 = atomlist.index(b1)
        z2 = atomlist.index(b2)
        ausgabe = "imptors{:8d}{:5d}{:5d}{:5d}{:17.3f}{:7.1f}{:3d}  !! GAFF".format(z1, z2, z3, z4, float(line[20:24]), float(line[33:38]), int(float(line[47:CURRENT_NUMBER_OF_TINKERTYPES])))
        print ausgabe
    
