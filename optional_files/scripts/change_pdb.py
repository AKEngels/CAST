### script to modify a PDB file (give new atom types) and change residue names
### improve: make one function where you give atom number and desired atom type that changes the type of this atom
###          make one function where you give the atom numbers of a residue and the desired residue name that changes the name of this residue

"""changes atom types"""
def rename_cyp(lines):
    lines[339] = lines[339][:13] + "N " + lines[339][15:]
    lines[340] = lines[340][:13] + "H " + lines[340][15:]
    lines[341] = lines[341][:13] + "CA" + lines[341][15:]
    lines[347] = lines[347][:13] + "C " + lines[347][15:]
    lines[348] = lines[348][:13] + "O " + lines[348][15:]

"""changes residue name from HIE to HIP"""
def change_his(lines):
    for i in range(2267,2285):
        lines[i] = lines[i].replace("HIE","HIP")

"""changes residue name from CYP to CYM"""
def change_cyp(lines):
    for i in range(339,349):
        lines[i] = lines[i].replace("CYP","CYM")

# open PDB file
with open("qmmm_opt.pdb") as pdbfile:
    lines = pdbfile.readlines()

# do stuff
rename_cyp(lines)
change_his(lines)
change_cyp(lines)

# write new PDB file
with open("new.pdb","w") as new_pdb:
    new_pdb.writelines(lines)

