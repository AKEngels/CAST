PRMTOP = "min2.prmtop"    #name of the prmtop-file

def get_flag(flag_name, signnumber):
    """searches the prmtop file for the flag (without %FLAG)
    and returns a list of the corresponding property of every atom
    
    signnumber is the length of the property-string for every atom
    (e.g. 4 for AMBER_ATOM_TYPE)
    """
    go = False
    flag_str = ""
    with open (PRMTOP) as prmtop_file:
        start = 10000000000000 # random big number
        for linenumber,line in enumerate(prmtop_file):
            if line.startswith("%FLAG "+flag_name):
                start = linenumber + 2
            if linenumber == start:
                go = True
            if linenumber > start and line.startswith("%"):
                end = linenumber
                break
            if go:
                flag_str = flag_str + line
    flag_str = flag_str.replace("\n","")
    flag = []
    i = 0
    while i < len(flag_str):
        new_atom = flag_str[i:i+signnumber]
        flag.append(new_atom)
        i = i+signnumber
    return flag

atom_types = get_flag("AMBER_ATOM_TYPE",4) # find atom types
charges = get_flag("CHARGE",16)            # find charges

if len(atom_types) != len(charges):
    print "ERROR"
else:    # writes output-file with atomnumber (compatible with tinker), atomtype and charge
    number_of_atoms = len(atom_types)
    with open("charges.txt", "w") as output:  
        for i in range(number_of_atoms):
            charge_koeff = float(charges[i][:12])
            charge_exp = int(charges[i][13:])
            charge = (charge_koeff * 10**charge_exp)/ 18.2223
            output.write("{}    {}    {}\n".format(i+1,atom_types[i], charge))
    

    
                                                  
                                                  
            

    
