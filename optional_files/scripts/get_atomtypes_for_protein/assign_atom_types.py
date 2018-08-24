### script that takes a tinkerstructure and a file where the single residues are defined (like in bla.txt, the script cut_peptide.py can help you to create this file)
### and assigns atom types to those residues that it recognises. The are written into a new tinkerstructure new_struc.arc.

STRUCTURENAME = "19_15.arc"
SEGMENTFILE = "bla.txt"

# number are the numbers of atoms in the amino acid. If an amino acid is terminal this number might be different, but than the atom types might also be different. In the moment the script can't handle this.
RESIDUE_NAMES_AND_NUMBERS = {"ALA":10, "ARG":24, "ASN":14, "ASP":12, "CYS":11, "GLN":17, "GLU":15, "GLY":7, "HIS":16, "ILE":19, "LEU":19, "LYS":22, "MET":17, "PHE":20, "SER":11, "THR":14, "TRP":24, "TYR":21, "VAL":16, "CYX":10, "CYM":10, "CYP":10, "HID":17, "HIE":17, "HIP":18}

class Atom:

    def __init__(self, tinker_line):
        linelist = tinker_line.split()
        self.symbol = linelist[1]
        self.x = linelist[2]
        self.y = linelist[3]
        self.z = linelist[4]

        bondstr = linelist[6:]
        self.bonds = []
        for b in bondstr:
            self.bonds.append(int(b)-1)

    def bonds_as_str(self):
        res = ""
        for b in self.bonds:
            res += str(b+1) + " "
        return res

    def write_to_xyz(self):
        return "{}  {}  {}  {}\n".format(self.symbol, self.x, self.y, self.z) 

    def info(self):
        print "{} at {}, {}, {} with bonding partners {}".format(self.symbol, self.x, self.y, self.z, self.bonds_as_str())


class Segment:

    def __init__(self, line):
        linelist = line.split()
        start = int(linelist[3])
        stop = int(linelist[5])
        self.index_range = range(start-1, stop)
        self.number_of_atoms = len(self.index_range)
        self.res_name = linelist[8]
        self.test_atom_number()

    def test_atom_number(self):
        if self.res_name.upper() in RESIDUE_NAMES_AND_NUMBERS:
            if RESIDUE_NAMES_AND_NUMBERS[self.res_name.upper()] == self.number_of_atoms:
                return True
            else:
                print "Something went wrong. {} should have {} atoms but has {}.".format(self.res_name, RESIDUE_NAMES_AND_NUMBERS[self.res_name.upper()], self.number_of_atoms)
                return False
        else:
            print "Unknown residue {}\n You will have to assign the atom types manually.".format(self.res_name)
            return True

    # Here the import part of the program takes place. There is a function for every amino acid (or sometimes a group of amino acids) that assigns the atom types to the side chain atoms.

    def get_atom_types_gly(self, types):
        types[3] = 85
        types[4] = 85
        return types

    def get_atom_types_ala(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                types[a] = 80
            elif atom.symbol == "H":
                types[a] = 85
            else:
                print "strange element in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_cyx(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                types[a] = 156
            elif atom.symbol == "H":
                types[a] = 85
            elif atom.symbol == "S":
                types[a] = 145
            else:
                print "strange element {} in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_met(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if bonding_symbols.count("C") == 2 and bonding_symbols.count("H") == 2:
                    types[a] = 81
                elif bonding_symbols.count("S") == 1 and bonding_symbols.count("H") == 3:
                    types[a] = 151
                elif bonding_symbols.count("S") == 1 and bonding_symbols.count("H") == 2 and bonding_symbols.count("C") == 1:
                    types[a] = 152
                else:
                    print "wrong bonding partners for C in {}".format(self.res_name)
            elif atom.symbol == "H":
                types[a] = 85
            elif atom.symbol == "S":
                types[a] = 144
            else:
                print "strange element {} in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_ser_thr(self, atoms, types):

        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if bonding_symbols.count("C") == 1 and bonding_symbols.count("H") == 3:
                    types[a] = 80
                elif bonding_symbols.count("C") == 1 and bonding_symbols.count("H") == 2 and bonding_symbols.count("O") == 1:
                    types[a] = 115
                elif bonding_symbols.count("C") == 2 and bonding_symbols.count("H") == 1 and bonding_symbols.count("O") == 1:
                    types[a] = 99
                else:
                    print "wrong bonding partners for C in {}".format(self.res_name)
            elif atom.symbol == "H":
                bonding_partner = atoms[atom.bonds[0]]
                if bonding_partner.symbol == "O":
                    types[a] = 97
                elif bonding_partner.symbol == "C":
                    types[a] = 85
                else:
                    print "wrong bonding partner for H in {}".format(self.res_name)

            elif atom.symbol == "O":
                types[a] = 96
            else:
                print "strange element {} in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_val_leu_ile(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if bonding_symbols.count("C") == 3 and bonding_symbols.count("H") == 1:
                    types[a] = 82
                elif bonding_symbols.count("H") == 3 and bonding_symbols.count("C") == 1:
                    types[a] = 80
                elif bonding_symbols.count("H") == 2 and bonding_symbols.count("C") == 2:
                    types[a] = 81
                else:
                    print "wrong bonding partners for C in {}".format(self.res_name)
            elif atom.symbol == "H":
                types[a] = 85
            else:
                print "strange element in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_lys(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if "N" in bonding_symbols:
                    types[a] = 236
                else:
                    types[a] = 81
            elif atom.symbol == "N":
                types[a] = 230
            elif atom.symbol == "H":
                bonding_partner = atoms[atom.bonds[0]]
                if bonding_partner.symbol == "N":
                    types[a] = 233
                elif bonding_partner.symbol == "C":
                    types[a] = 85
                else:
                    print "wrong bonding partner for H in Lys"
            else:
                print "strange element in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_asn_gln(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if "N" in bonding_symbols:
                    types[a] = 177
                else:
                    types[a] = 81
            elif atom.symbol == "H":
                bonding_partner = atoms[atom.bonds[0]]
                if bonding_partner.symbol == "N":
                    types[a] = 182
                elif bonding_partner.symbol == "C":
                    types[a] = 85
                else:
                    print "wrong bonding partner for H in {}".format(self.res_name)
            elif atom.symbol == "O":
                types[a] = 178
            elif atom.symbol == "N":
                types[a] = 179
            else:
                print "strange element {} in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_asp_glu(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if "O" in bonding_symbols:
                    types[a] = 213
                else:
                    types[a] = 81
            elif atom.symbol == "H":
                types[a] = 85
            elif atom.symbol == "O":
                types[a] = 214
            else:
                print "strange element {} in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_phe_tyr(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if bonding_symbols.count("C") == 2 and bonding_symbols.count("H") == 2:
                    types[a] = 81
                elif bonding_symbols.count("C") == 2 and bonding_symbols.count("H") == 1 or bonding_symbols.count("C") == 3:
                    types[a] = 90
                elif bonding_symbols.count("C") == 2 and bonding_symbols.count("O") == 1:
                    types[a] = 108
                else:
                    print "wrong bonding partners for C in {}".format(self.res_name)
            elif atom.symbol == "H":
                bonding_partner = atoms[atom.bonds[0]]
                if bonding_partner.symbol == "O":
                    types[a] = 110
                elif bonding_partner.symbol == "C":
                    if len(bonding_partner.bonds) == 3:
                        types[a] = 91
                    else:
                        types[a] = 85
                else:
                    print "wrong bonding partner for H in {}".format(self.res_name)
            elif atom.symbol == "O":
                types[a] = 109
            else:
                print "strange element {} in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_hie(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if bonding_symbols.count("C") == 2 and bonding_symbols.count("H") == 2:
                    types[a] = 81
                elif bonding_symbols.count("C") == 2 and bonding_symbols.count("N") == 1:
                    types[a] = 448
                elif bonding_symbols.count("N") == 2 and bonding_symbols.count("H") == 1:
                    types[a] = 447
                elif bonding_symbols.count("C") == 1 and bonding_symbols.count("H") == 1 and bonding_symbols.count("N") == 1:
                    types[a] = 449
                else:
                    print "wrong bonding partners for C in {}".format(self.res_name)
            elif atom.symbol == "H":
                bonding_partner = atoms[atom.bonds[0]]
                if bonding_partner.symbol == "N":
                    types[a] = 445
                elif bonding_partner.symbol == "C":
                    if len(bonding_partner.bonds) == 3:
                        types[a] = 91
                    else:
                        types[a] = 85
                else:
                    print "wrong bonding partner for H in {}".format(self.res_name)
            elif atom.symbol == "N":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if len(bonding_symbols) == 3:
                    types[a] = 444
                elif len(bonding_symbols) == 2:
                    types[a] = 452
                else:
                    print "wrong number of bonding partners for N in {}".format(self.res_name)
            else:
                print "strange element {} in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_trp(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "C":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if len(bonding_symbols) == 4:
                    types[a] = 81
                elif len(bonding_symbols) == 3:
                    if bonding_symbols.count("C") == 2 and bonding_symbols.count("N") == 1:
                        types[a] = 443
                    elif bonding_symbols.count("C") == 2 and bonding_symbols.count("H") == 1:
                        types[a] = 90
                    elif bonding_symbols.count("C") == 1 and bonding_symbols.count("N") == 1 and bonding_symbols.count("H") == 1:
                        types[a] = 455
                    else:
                        for b in atom.bonds:
                            bbonds = []
                            for bb in atoms[b].bonds:
                                bbonds.append(atoms[bb].symbol)
                            if bbonds.count("C") == 3:
                                pass
                            elif bbonds.count("C") == 2 and bbonds.count("H") == 2:
                                types[a] = 441
                            elif bbonds.count("C") == 1 and bbonds.count("H") == 1 and bbonds.count("N") == 1:
                                types[a] = 441
                            elif bbonds.count("C") == 2 and bbonds.count("H") == 1:
                                types[a] = 442
                            elif bbonds.count("C") == 2 and bbonds.count("N") == 1:
                                types[a] = 441
                            else:
                                print "something went wrong with C in Trp"
                else:
                    print "wrong bonding partners for C in {}".format(self.res_name)
            elif atom.symbol == "H":
                bonding_partner = atoms[atom.bonds[0]]
                if bonding_partner.symbol == "N":
                    types[a] = 445
                elif bonding_partner.symbol == "C":
                    if len(bonding_partner.bonds) == 3:
                        types[a] = 91
                    else:
                        types[a] = 85
                else:
                    print "wrong bonding partner for H in {}".format(self.res_name)
            elif atom.symbol == "N":
                types[a] = 444
            else:
                print "strange element {} in residue {}".format(atom.symbol, self.res_name)
        return types

    def get_atom_types_arg(self, atoms, types):
        for a in range(3, self.number_of_atoms-2):
            atom = atoms[a + self.index_range[0]]
            if atom.symbol == "N":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if bonding_symbols.count("C") == 2 and bonding_symbols.count("H") == 1:
                    types[a] = 246
                elif bonding_symbols.count("C") == 1 and bonding_symbols.count("H") == 2:
                    types[a] = 243
                else:
                    print "wrong bonding partners for N in {}".format(self.res_name)
            elif atom.symbol == "H":
                bonding_partner = atoms[atom.bonds[0]]
                if bonding_partner.symbol == "C":
                    types[a] = 85
                elif bonding_partner.symbol == "N":
                    bbonds = []
                    for bb in bonding_partner.bonds:
                        bbonds.append(atoms[bb].symbol)
                    if bbonds.count("C") == 2 and bbonds.count("H") == 1:
                        types[a] = 247
                    elif bbonds.count("C") == 1 and bbonds.count("H") == 2:
                        types[a] = 244
                    else:
                        print "Something went wrong with H bound to N in Arg"
                else:
                    print "wrong bonding partner for H in {}".format(self.res_name)
            elif atom.symbol == "C":
                bonding_symbols = []
                for b in atom.bonds:
                    bonding_symbols.append(atoms[b].symbol)
                if len(bonding_symbols) == 3:
                    types[a] = 245
                elif bonding_symbols.count("C") == 1 and bonding_symbols.count("H") == 2 and bonding_symbols.count("N") == 1:
                    types[a] = 250
                else:
                    for b in atom.bonds:
                        if atoms[b].symbol == "C":
                            bbonds = []
                            for bb in atoms[b].bonds:
                                bbonds.append(atoms[bb].symbol)
                            if bbonds.count("C") == 2 and bbonds.count("H") == 2:
                                pass
                            elif bbonds.count("C") == 2 and bbonds.count("H") == 1 and bbonds.count("N") == 1:
                                types[a] = 81
                            elif bbonds.count("C") == 1 and bbonds.count("H") == 2 and bbonds.count("N") == 1:
                                types[a] = 251
                            else:
                                print "something went wrong with C_beta or C_gamma in Arg"
            else:
                print "strange element {} in residue {}".format(atom.symbol, self.res_name)
        return types

    def find_atom_types(self, atoms):
        """function to assign atom types to a residue"""
        
        types = [0] * self.number_of_atoms  # as many types as there are atoms (first filled with 0)

        if self.res_name.upper() not in RESIDUE_NAMES_AND_NUMBERS:
            return types

        else:   # residue is recognised

            # backbone atomtypes
            types[0] = 180
            types[1] = 183
            if self.res_name == "Gly":
                types[2] = 165
            else:
                types[2] = 166
            types[-1] = 178
            types[-2] = 177

            # side chain atomtypes
            if self.res_name == "Gly":
                types = self.get_atom_types_gly(types)
            elif self.res_name == "Ala":
                types = self.get_atom_types_ala(atoms, types)
            elif self.res_name == "Cyx":
                types = self.get_atom_types_cyx(atoms, types)
            elif self.res_name == "Ser" or self.res_name == "Thr":
                types = self.get_atom_types_ser_thr(atoms, types)
            elif self.res_name == "Val" or self.res_name == "Leu" or self.res_name == "Ile":
                types = self.get_atom_types_val_leu_ile(atoms, types)
            elif self.res_name == "Lys":
                types = self.get_atom_types_lys(atoms, types)
            elif self.res_name == "Asn" or self.res_name == "Gln":
                types = self.get_atom_types_asn_gln(atoms, types)
            elif self.res_name == "Asp" or self.res_name == "Glu":
                types = self.get_atom_types_asp_glu(atoms, types)
            elif self.res_name == "Met":
                types = self.get_atom_types_met(atoms, types)
            elif self.res_name == "Phe" or self.res_name == "Tyr":
                types = self.get_atom_types_phe_tyr(atoms, types)
            elif self.res_name == "Arg":
                types = self.get_atom_types_arg(atoms, types)
            elif self.res_name == "Trp":
                types = self.get_atom_types_trp(atoms, types)
            elif self.res_name == "Hie":
                types = self.get_atom_types_hie(atoms, types)
            else:
                print "residue {} not implemented yet".format(self.res_name)

        return types


    def info(self):
        print "Segment from {} to {} ({} atoms) is a {}".format(self.index_range[0]+1, self.index_range[-1]+1, self.number_of_atoms, self.res_name)
        


### get atomic information from tinkerstructure 
with open(STRUCTURENAME) as tinkerfile:
    lines = tinkerfile.readlines()

atoms = []
for i, line in enumerate(lines):
    if i != 0:
        a = Atom(line)
        atoms.append(a)

### get information about the residues from the segmentfile
with open(SEGMENTFILE) as segfile:
    lines = segfile.readlines()

segments = []
for line in lines:
    s = Segment(line)
    segments.append(s)

### assign atomtypes
atomtypes = []
for s in segments:
    atomtypes += s.find_atom_types(atoms)

### write new tinkerstructure
with open("new_struc.arc", "w") as output:
    output.write(str(len(atoms)) + "\n")
    for i, a in enumerate(atoms):
        output.write("{}  {}  {} {} {}  {}  {}\n".format(i+1, a.symbol, a.x, a.y, a.z, atomtypes[i], a.bonds_as_str()))



    
    

            
