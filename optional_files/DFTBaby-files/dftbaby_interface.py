import sys
from DFTB.DFTB2 import DFTB2
from DFTB.XYZ import read_xyz
from DFTB.LR_TDDFTB import LR_TDDFTB
from DFTB.ExcGradients import Gradients

SCF_OPTIONLIST = ['scf_conv', 'start_from_previous', 'level_shift', 'density_mixer', 'fock_interpolator', 'mixing_threshold', 'HOMO_LUMO_tol', 'maxiter', 'linear_mixing_coefficient']
TD_INIT_OPTIONLIST = ["parameter_set","point_charges_xyz", "initial_charge_guess", "save_converged_charges", "verbose", "distance_cutoff", "long_range_correction", "long_range_radius", "long_range_T", "long_range_switching", "lc_implementation", "tune_range_radius", "save_tuning_curve", "nr_unpaired_electrons", "use_symmetry", "fluctuation_functions", "mulliken_dipoles", "dispersion_correction", "qmmm_partitioning", "qmmm_embedding", "periodic_force_field", "cavity_radius", "cavity_force_constant", "scratch_dir", "cpks_solver"]
GRAD_OPTIONS = ['gradient_file', 'gradient_check', 'gradient_state']


def extract_options(options, optionlist):
    """extracts options for a given function (defined by optionlist) from total options"""
    scf_keylist = []
    scf_valuelist = []
    for i,o in enumerate(options.keys()):
        if o in optionlist:
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

    outputfile = open("output_dftb.txt", "a")  # redirect output to file
    sys.stdout = outputfile

    options = read_options(optionfile)  # read options
    atomlist = read_xyz(xyzfile)[0]  # read structure

    dftb2 = DFTB2(atomlist, **options)  # create dftb object
    dftb2.setGeometry(atomlist)

    scf_options = extract_options(options, SCF_OPTIONLIST) # calculate energy
    dftb2.getEnergy(**scf_options)
    energies = list(dftb2.getEnergies())  # get partial energies

    if dftb2.long_range_correction == 1: # add long range correction to partial energies
        energies.append(dftb2.E_HF_x)

    return str(energies)


def calc_gradients(xyzfile, optionfile):
    """calculates DFTB energies and gradients for a molecule in the xyz_file
    with the options given in the optionfile"""

    outputfile = open("output_dftb.txt", "a")  # redirect output to file
    sys.stdout = outputfile

    options = read_options(optionfile)  # read options
    atomlist = read_xyz(xyzfile)[0]  # read structure

    init_options = extract_options(options, TD_INIT_OPTIONLIST)
    scf_options = extract_options(options, SCF_OPTIONLIST)
    grad_options = extract_options(options, GRAD_OPTIONS)

    tddftb = LR_TDDFTB(atomlist, **init_options)  # create object
    tddftb.setGeometry(atomlist)
    tddftb.getEnergies(**scf_options)  # calculate energies

    grad = Gradients(tddftb)            # calculate gradients
    grad.getGradients(**grad_options)

    energies = list(tddftb.dftb2.getEnergies())  # get partial energies

    if tddftb.dftb2.long_range_correction == 1: # add long range correction to partial energies
        energies.append(tddftb.dftb2.E_HF_x)

    return str(energies)




