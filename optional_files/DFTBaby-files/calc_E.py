import sys
from DFTB2 import DFTB2
from XYZ import read_xyz

SCF_OPTIONLIST = ['scf_conv', 'start_from_previous', 'level_shift', 'density_mixer', 'fock_interpolator', 'mixing_threshold', 'HOMO_LUMO_tol', 'maxiter', 'linear_mixing_coefficient']

def extract_scf_options(options):
    """extracts SCF options from total options"""
    scf_keylist = []
    scf_valuelist = []
    for i,o in enumerate(options.keys()):
        if o in SCF_OPTIONLIST:
            scf_keylist.append(o)
            scf_valuelist.append(options.values()[i])
    return dict(zip(scf_keylist, scf_valuelist))

def read_options(filename):
    """reads options from optionfile"""

    optionlist = []
    valuelist_str = []
    valuelist = []

    with open(filename) as optfile:
        for line in optfile.readlines():
            if line.find("=") != -1:
                option = line.split("=")
                optionlist.append(remove_spaces(option[0]))
                valuelist_str.append(remove_spaces(option[1][:-1]))

    for v in valuelist_str:
        valuelist.append(convert_string(v))

    return dict(zip(optionlist, valuelist))


def remove_spaces(string):
    """removes unnecessary spaces from a string

    input: string
    return: string without the unnecessary spaces"""

    newstring = ""
    for i, a in enumerate(string):
        if i == 0 and a == " ":  # spaces at the beginning
            pass
        elif a == " " and i == len(string) - 1:  # spaces at the end
            return remove_spaces(newstring)
        elif a == " " and string[i - 1] == " ":  # two spaces back to back
            pass
        else:
            newstring = newstring + a
    return newstring

def convert_string(string):
    """converts string to float or int"""
    try:
        string = float(string)
        if round(string, 0) == string:
            string = int(string)
    except:
        pass
    return string


def calc_energies(xyzfile, optionfile):
    """calculates DFTB energies for a molecule in the xyz_file
    with the options given in the optionfile"""

    outputfile = open("scf_output_dftb.txt", "a")  # redirect output to file
    sys.stdout = outputfile

    options = read_options(optionfile)  # read options
    atomlist = read_xyz(xyzfile)[0]  # read structure

    dftb2 = DFTB2(atomlist, **options)  # create dftb object
    dftb2.setGeometry(atomlist)

    scf_options = extract_scf_options(options) # calculate energy
    dftb2.getEnergy(**scf_options)
    energies = list(dftb2.getEnergies())  # get partial energies

    if dftb2.long_range_correction == 1: # add long range correction to partial energies
        energies.append(dftb2.E_HF_x)

    return str(energies)




