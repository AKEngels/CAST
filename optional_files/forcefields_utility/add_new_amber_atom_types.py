"""This is a script to create the if loop in function amberUtil::toTinkerType out of the atom section of an tinker type amber forcefield"""

with open("atoms.txt") as parameterfile:
    atoms_str = parameterfile.readlines()

saved_types = []
for a in atoms_str:
    number = a[11:15]
    symbol = a[24:27]
    if symbol not in saved_types:
        print 'else if (in == "{}") return {}u;'.format(symbol, number)
        saved_types.append(symbol)
