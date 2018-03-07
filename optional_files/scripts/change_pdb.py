### script to modify a PDB file (give new atom types and change residue names)

"""changes a residue name
filelines: inputfile read in with readlines()
reslines: lines which contain the residues (in tinker/pdb numbering)
oldname: old residue name
newname: new residue name"""
def change_resname(filelines, reslines, oldname, newname):
    for i in reslines:
        filelines[i-1] = filelines[i-1].replace(oldname, newname)

"""changes the atom type of an atom
filelines: inputfile read in with readlines()
changeline: number of the line where the atom should be changed (in tinker/pdb numbering)
new_atomtype: name of the desired atom type"""
def change_atomtype(filelines, changeline, new_atomtype):
    filelines[changeline-1] = filelines[changeline-1][:13] + "{:2}".format(new_atomtype) + filelines[changeline-1][15:]


"""changes atom types of CYP"""
def rename_cyp(lines):
    change_atomtype(lines, 340, "N")
    change_atomtype(lines, 341, "H")
    change_atomtype(lines, 342, "CA")
    change_atomtype(lines, 348, "C")
    change_atomtype(lines, 349, "O")


# open PDB file
with open("qmmm_opt.pdb") as pdbfile:
    lines = pdbfile.readlines()

# do stuff
#rename_cyp(lines)
change_resname(lines, range(2268,2286), "HIP","HIE")
change_resname(lines, range(340,350), "CYM", "CYP")

# write new PDB file
with open("new.pdb","w") as new_pdb:
    new_pdb.writelines(lines)

