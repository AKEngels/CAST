### script to modify a PDB file (give new atom types and change residue names)

"""changes a residue name
filelines: inputfile read in with readlines()
reslines: lines which contain the residues
oldname: old residue name
newname: new residue name"""
def change_resname(filelines, reslines, oldname, newname):
    for i in reslines:
        filelines[i] = filelines[i].replace(oldname, newname)

"""changes the atom type of an atom
filelines: inputfile read in with readlines()
changeline: number of the line where the atom should be changed
new_atomtype: name of the desired atom type"""
def change_atomtype(filelines, changeline, new_atomtype):
    filelines[changeline] = filelines[changeline][:13] + "{:2}".format(new_atomtype) + filelines[changeline][15:]


"""changes atom types of CYP"""
def rename_cyp(lines):
    change_atomtype(lines, 339, "N")
    change_atomtype(lines, 340, "H")
    change_atomtype(lines, 341, "CA")
    change_atomtype(lines, 347, "C")
    change_atomtype(lines, 348, "O")


# open PDB file
with open("qmmm_opt.pdb") as pdbfile:
    lines = pdbfile.readlines()

# do stuff
rename_cyp(lines)
change_resname(lines, range(2267,2285), "HIE","HIP")
change_resname(lines, range(339,349), "CYP", "CYM")

# write new PDB file
with open("new.pdb","w") as new_pdb:
    new_pdb.writelines(lines)

