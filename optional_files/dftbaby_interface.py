import sys
import numpy
from scipy import optimize
from copy import copy

from DFTB import utils, XYZ
from DFTB.DFTB2 import DFTB2, annotated_hessian
from DFTB.LR_TDDFTB import LR_TDDFTB
from DFTB.ExcGradients import Gradients
from DFTB.PES import PotentialEnergySurfaces
from DFTB.Dynamics import HarmonicApproximation


SCF_OPTIONLIST = ['scf_conv', 'start_from_previous', 'level_shift', 'density_mixer', 'fock_interpolator', 'mixing_threshold', 'HOMO_LUMO_tol', 'maxiter', 'linear_mixing_coefficient']
TD_INIT_OPTIONLIST = ["parameter_set", "point_charges_xyz", "initial_charge_guess", "save_converged_charges", "verbose", "distance_cutoff", "long_range_correction", "long_range_radius", "long_range_T", "long_range_switching", "lc_implementation", "tune_range_radius", "save_tuning_curve", "nr_unpaired_electrons", "use_symmetry", "fluctuation_functions", "mulliken_dipoles", "dispersion_correction", "qmmm_partitioning", "qmmm_embedding", "periodic_force_field", "cavity_radius", "cavity_force_constant", "scratch_dir", "cpks_solver"]
TD_OPTIONLIST = ["nr_active_occ", "nr_active_virt", "select_lm", "oszis", "response_method", "multiplicity", "nstates", "diag_ifact", "diag_conv", "diag_maxiter", "diag_check", "diag_L2threshold", "diag_selected_irreps", "ct_correction"] + SCF_OPTIONLIST
GRAD_OPTIONS = ['gradient_file', 'gradient_check', 'gradient_state']


class MyPES(PotentialEnergySurfaces):
    """overwritten __init__ function for PotentialEnergySurface
    because in original parse_args() is used
    that doesn't work together with C++"""
  
    def __init__(self, atomlist, options, Nst=2, **kwds):
        usage = "Type --help to show all options for DFTB"
        parser = utils.OptionParserFuncWrapper(
        [DFTB2.__init__,
         DFTB2.runSCC,
         LR_TDDFTB.getEnergies
         ], usage)
        self.options = options
        td_init_options = extract_options(self.options, TD_INIT_OPTIONLIST)
        self.atomlist = atomlist
        self.tddftb = LR_TDDFTB(atomlist, **td_init_options)
        self.grads = Gradients(self.tddftb)
        self.tddftb.setGeometry(atomlist, charge=kwds.get("charge", 0.0))
        self.scf_options = extract_options(self.options, SCF_OPTIONLIST)
        self.options = copy(self.scf_options)
        self.Nst = Nst
        #        # always use iterative diagonalizer for lowest Nst-1 excited states
        self.options["nstates"] = Nst - 1
        # save geometry, orbitals and TD-DFT coefficients from
        # last calculation
        self.last_calculation = None
        # save transition dipoles from last calculation
        self.tdip_old = None


def extract_options(options, optionlist):
    """extracts options for a given function (defined by optionlist) from total options"""
    scf_keylist = []
    scf_valuelist = []
    for i, o in enumerate(options.keys()):
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
    atomlist = XYZ.read_xyz(xyzfile)[0]  # read structure
    kwds = XYZ.extract_keywords_xyz(xyzfile)
    
    try:
        dftb2 = DFTB2(atomlist, **options)  # create dftb object
        dftb2.setGeometry(atomlist, charge=kwds.get("charge", 0.0))

        scf_options = extract_options(options, SCF_OPTIONLIST)  # calculate energy
        dftb2.getEnergy(**scf_options)
        energies = list(dftb2.getEnergies())  # get partial energies
    except:
        return "error"

    if dftb2.long_range_correction == 1:  # add long range correction to partial energies
        energies.append(dftb2.E_HF_x)

    return str(energies)


def calc_gradients(xyzfile, optionfile):
    """calculates DFTB energies and gradients for a molecule in the xyz_file
    with the options given in the optionfile"""

    outputfile = open("output_dftb.txt", "a")  # redirect output to file
    sys.stdout = outputfile

    options = read_options(optionfile)  # read options
    atomlist = XYZ.read_xyz(xyzfile)[0]  # read structure

    init_options = extract_options(options, TD_INIT_OPTIONLIST)
    td_options = extract_options(options, TD_OPTIONLIST)
    grad_options = extract_options(options, GRAD_OPTIONS)
    kwds = XYZ.extract_keywords_xyz(xyzfile)
    
    try:
        tddftb = LR_TDDFTB(atomlist, **init_options)  # create object
        tddftb.setGeometry(atomlist, charge=kwds.get("charge", 0.0))
        tddftb.getEnergies(**td_options)  # calculate energies

        grad = Gradients(tddftb)            # calculate gradients
        grad.getGradients(**grad_options)

        energies = list(tddftb.dftb2.getEnergies())  # get partial energies
    except:
        return "error"

    if tddftb.dftb2.long_range_correction == 1: # add long range correction to partial energies
        energies.append(tddftb.dftb2.E_HF_x)

    return str(energies)


def opt(xyzfile, optionfile):
    """performs an optimization"""
  
    outputfile = open("output_dftb.txt", "a")  # redirect output to file
    sys.stdout = outputfile
  
    I = 0  # index of electronic state (ground state)
  
    atomlist = XYZ.read_xyz(xyzfile)[0]   # read atomlist
    kwds = XYZ.extract_keywords_xyz(xyzfile)  # read keywords from xyz-file (charge)
    options = read_options(optionfile)  # read options
    scf_options = extract_options(options, SCF_OPTIONLIST)  # get scf-options
    
    try:
        # optimization (taken from optimize.py)
        pes = MyPES(atomlist, options, Nst=max(I + 1, 2), **kwds)

        x0 = XYZ.atomlist2vector(atomlist)  #convert geometry to a vector

        def f(x):
            save_xyz(x)  # also save geometries from line searches

            if I == 0 and type(pes.tddftb.XmY) != type(None):
            # only ground state is needed. However, at the start
            # a single TD-DFT calculation is performed to initialize
            # all variables (e.g. X-Y), so that the program does not
            # complain about non-existing variables.
                enI, gradI = pes.getEnergyAndGradient_S0(x)
            else:
                energies, gradI = pes.getEnergiesAndGradient(x, I)
                enI = energies[I]
            print "E = %2.7f" % (enI)
            return enI, gradI

        xyz_trace = xyzfile.replace(".xyz", "_trace.xyz")

        # This is a callback function that is executed by numpy for each optimization step.
        # It appends the current geometry to an xyz-file.
        def save_xyz(x, mode="a"):
            atomlist_opt = XYZ.vector2atomlist(x, atomlist)
            XYZ.write_xyz(xyz_trace, [atomlist_opt],
                  title="charge=%s" % kwds.get("charge", 0),
                  mode=mode)

        save_xyz(x0, mode="w")  # write original geometry

        Nat = len(atomlist)
        min_options = {'gtol': 1.0e-7, 'norm': 2}
        # The "BFGS" method is probably better than "CG", but the line search in BFGS is expensive.
        res = optimize.minimize(f, x0, method="CG", jac=True, callback=save_xyz, options=min_options)
        # res = optimize.minimize(f, x0, method="BFGS", jac=True, callback=save_xyz, options=options)
        xopt = res.x
        save_xyz(xopt)

        print "Intermediate geometries written into file {}".format(xyz_trace)
    
    
        # write optimized geometry into file
        atomlist_opt = XYZ.vector2atomlist(xopt, atomlist)
        xyz_opt = xyzfile.replace(".xyz", "_opt.xyz")
        XYZ.write_xyz(xyz_opt, [atomlist_opt],
                  title="charge=%s" % kwds.get("charge", 0),
                  mode="w")

        # calculate energy for optimized geometry
        dftb2 = DFTB2(atomlist_opt, **options)  # create dftb object
        dftb2.setGeometry(atomlist_opt, charge=kwds.get("charge", 0.0))
    
        dftb2.getEnergy(**scf_options)
        energies = list(dftb2.getEnergies())  # get partial energies
        
    except:
        return "error"

    if dftb2.long_range_correction == 1:  # add long range correction to partial energies
        energies.append(dftb2.E_HF_x)

    return str(energies)


def hessian(xyzfile, optionfile):
    """calculates hessian matrix"""

    outputfile = open("output_dftb.txt", "a")  # redirect output to file
    sys.stdout = outputfile

    I = 0  # index of electronic state (ground state)

    atomlist = XYZ.read_xyz(xyzfile)[0]  # read xyz file
    kwds = XYZ.extract_keywords_xyz(xyzfile)  # read keywords (charge)
    options = read_options(optionfile)  # read options
    scf_options = extract_options(options, SCF_OPTIONLIST)  # get scf-options
    
    try:
        pes = MyPES(atomlist, options, Nst=max(I + 1, 2), **kwds)  # create PES
  
        atomvec = XYZ.atomlist2vector(atomlist)  # convert atomlist to vector
  
        # FIND ENERGY MINIMUM
        # f is the objective function that should be minimized
        # it returns (f(x), f'(x))
        def f(x):
            if I == 0 and type(pes.tddftb.XmY) != type(None):
                # only ground state is needed. However, at the start
                # a single TD-DFT calculation is performed to initialize
                # all variables (e.g. X-Y), so that the program does not
                # complain about non-existing variables.
                enI, gradI = pes.getEnergyAndGradient_S0(x)
            else:
                energies, gradI = pes.getEnergiesAndGradient(x, I)
                enI = energies[I]
            return enI, gradI

        minoptions = {'gtol': 1.0e-7, 'norm': 2}
        # somehow numerical_hessian does not work without doing this mimimization before
        res = optimize.minimize(f, atomvec, method="CG", jac=True, options=minoptions)
  
        # COMPUTE HESSIAN AND VIBRATIONAL MODES
        # The hessian is calculated by numerical differentiation of the
        # analytical gradients
        def grad(x):
            if I == 0:
                enI, gradI = pes.getEnergyAndGradient_S0(x)
            else:
                energies, gradI = pes.getEnergiesAndGradient(x, I)
            return gradI
    
        print "Computing Hessian"   # calculate hessians from gradients
        hess = HarmonicApproximation.numerical_hessian_G(grad, atomvec)
    
        string = ""   # create string that is to be written into file
        for line in hess:
            for column in line:
                string += str(column) + " "
            string = string[:-1] + "\n"
    
        with open("hessian.txt", "w") as hessianfile:  # write hessian matrix to file
            hessianfile.write(string)
            # this would look nicer but is not as exact
            #hessianfile.write(annotated_hessian(atomlist, hess))

        # calculate energy for optimized geometry
        dftb2 = DFTB2(atomlist, **options)  # create dftb object
        dftb2.setGeometry(atomlist, charge=kwds.get("charge", 0.0))
    
        dftb2.getEnergy(**scf_options)
        energies = list(dftb2.getEnergies())  # get partial energies
        
    except:
        return "error"

    if dftb2.long_range_correction == 1:  # add long range correction to partial energies
        energies.append(dftb2.E_HF_x)

    return str(energies)