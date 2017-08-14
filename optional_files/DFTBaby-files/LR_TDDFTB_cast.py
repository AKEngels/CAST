#!/usr/bin/env python
"""
linear response time dependent tight binding density functional theory (LR-TDDFTB)
gives excited state energies and gradients
"""
from numpy import array, where, argsort, zeros, dot, tensordot, transpose, diag, reshape, sqrt, arange, outer, sum, meshgrid
import numpy as np
from numpy.linalg import eigh
import numpy.linalg as la
from scipy.special import erf
from os import path
import string
import sys
sys.path.append("..")

from DFTB.DFTB2_cast import DFTB2
from DFTB.AtomicData import hartree_to_eV, hartree_to_nm, atom_names
from DFTB import utils_cast as utils
from DFTB.Timer import GlobalTimer as T
from DFTB import Solver
from DFTB.Analyse.SymmetryAssignment import SymmetryAssignmentEx_new as SymmetryAssignmentEx
from DFTB import XYZ
from DFTB.ExcGradients import Gradients
from DFTB.LongRangeCT import LongRangeCTCorrection_Neugebauer as LongRangeCTCorrection

class LR_TDDFTB:
    XpY = None
    XmY = None
    def __init__(self, atomlist, **opts):
        self.dftb2 = DFTB2(atomlist, **opts)
    def setGeometry(self,atomlist, **opts):
        self.dftb2.setGeometry(atomlist, **opts)
    def getEnergies(self, \
                        nr_active_occ=None, nr_active_virt=None, select_lm=None, \
                        oszis="Dipoles", response_method="Casida", multiplicity="S", \
                        nstates=None, diag_ifact=1, diag_conv=1.0e-5, diag_maxiter=50, diag_check=0, \
                        diag_L2threshold=0.0, \
                        diag_selected_irreps=None, \
                        ct_correction=0, \
                        **scf_opts):
        """
        Parameters:
        ==========
        Linear-Response TD-DFTB.nr_active_occ: 
              integer, number of active occupied orbitals \
              if set to None, all occupied orbitals are selected 
        Linear-Response TD-DFTB.nr_active_virt: 
              integer, number of active virtual orbitals \
              if set to None, all virtual orbitals are included
        Linear-Response TD-DFTB.select_lm:
              Molecular orbitals are included in the active space if they contain atomic orbitals with angular quantum numbers (l,m). This options allows to select only the pz-orbitals in a planar conjugated molecule to mimic a PPP calculation, e.g. --select_lm='[(1,0)]'. If symmetry is enabled, the molecule will be reoriented and the pz-orbitals might actually be the px or py orbitals depending on the point group.
        Linear-Response TD-DFTB.oszis: 
              method by which oscillator strengths are calculated, \
              "Mulliken" -> from Mulliken transition charges, \
              "Dipoles" -> from transition dipole moments (not available with 'mio' or 'hotbit' parameters)
        Linear-Response TD-DFTB.response_method: 
              "Casida" (for RPA) or "Tamm-Dancoff"
        Linear-Response TD-DFTB.multiplicity:
              "S" for Singlets or "T" for Triplets
        Davidson-like Diagonalization.nstates: Solve iteratively for the lowest 'nstates' eigenvalues of the non-Hermitian TD-DFTB equations. If this options if not set, the full TD-DFTB matrix is constructed in memory and diagonalized.
        Davidson-like Diagonalization.diag_ifact: ifact*nstates singly excited states are created as initial guesses. This factor should be increased if the resulting states do not have the desired symmetry.
        Davidson-like Diagonalization.diag_conv: convergence threshold for eigenvalues
        Davidson-like Diagonalization.diag_maxiter: maximum number of iterations
        Davidson-like Diagonalization.diag_check: compare iterative solution with the full diagonalization (only for debugging purposes)
        Davidson-like Diagonalization.diag_L2threshold: solution vectors whose Lambda2 value is below this threshold are removed.
        Davidson-like Diagonalization.diag_selected_irreps: If symmetry is enabled, expansion vectors belonging to irreps listed in diag_selected_irreps are included preferentially. Note the use of quotes and double quotes, e.g. --diag_selected_irreps="['B1U','B2U']". 
        Long-Range Charge-Transfer.ct_correction: 0 -> no correction, 1 -> shifts low-lying spurious CT states to higher energies as described in J.Chem.Phys. 124, 213102 (2006)
        """

        self.en0 = self.dftb2.getEnergy(**scf_opts)
        if select_lm != None:
            self.setActiveLMOrbitals(nr_active_occ,nr_active_virt,select_lm=select_lm)
        else:
            self.setActiveOrbitals(nr_active_occ,nr_active_virt)
        self.setMultiplicity(multiplicity)
        self.nstates=nstates
        if nstates != None:
            assert nstates > 0, "nstates==0, no excited states can be calculated!"
        self._solve_lr_equations(\
                response_method=response_method, multiplicity=multiplicity, \
                nstates=nstates, diag_ifact=diag_ifact, \
                diag_conv=diag_conv, diag_maxiter=diag_maxiter, \
                diag_check=diag_check, diag_L2threshold=diag_L2threshold, \
                diag_selected_irreps=diag_selected_irreps, \
                ct_correction=ct_correction)
        if oszis == "Dipoles" and self.dftb2.parameter_set == "homegrown":
            if self.dftb2.verbose > 0:
                print "compute transition dipoles from Slater-Koster tables"
            self._TransitionDipoles()
        else:
            if self.dftb2.verbose > 0:
                print "compute transition dipoles using Mulliken approximation (ZDO)"
            self._TransitionDipolesMulliken()
        self.getOscillatorStrengths()
        #
        self._Lambda2_diagnostics()
        self._SymmetryOfExcitations()
        if nstates != None:
            self.determineActiveSpace()
        if self.dftb2.verbose > 0:
            self._writeExcitations()
        return self.Omega, self.oscillator_strength, self.en0


    def energy_func(self, x, atomlist, I):
        """helper function for numerical differentiation of excitation energies"""
        # fetch options from command line
        (options,args) = parser.parse_args(self.getEnergies)
        (scf_options,args) = parser.parse_args(self.dftb2.runSCC)
        options.update(scf_options)
        atpos = XYZ.vector2atomlist(x, atomlist)
        self.setGeometry(atpos)
        self.getEnergies(**options)
        return self.Omega[I]
    def setActiveLMOrbitals(self, nr_active_occ=1000, nr_active_virt=1000, select_lm=[(1,0)]):
        """
        Only those molecular orbitals are placed into the active space that contain atomic
        orbitals with angular and magnetic quantum numbers (l,m) from the set in select_lm.
        """
        atomlist = self.dftb2.getGeometry()
        valorbs = self.dftb2.getValorbs()
        orbs = self.dftb2.getKSCoefficients()
        lm_indeces = [] # indeces of atomic orbitals with quantum numbers (l,m) from the selected set
        iao = 0 # this counter runs over all atomic orbitals, iao=0,...,nao-1
        print "The active space contains only atomic orbitals with (l,m) in %s" % str(select_lm)
        # iterate over all atomic centers
        for i,(Zi,posi) in enumerate(atomlist):
            # iterate over all atomic orbitals at center i
            for (ni,li,mi) in valorbs[Zi]:
                # index iao is added only if (li,mi) is of the right type
                if (li,mi) in select_lm:
                    lm_indeces.append( iao )
                iao += 1
        lm_indeces = np.array(lm_indeces, dtype=int)
        nao,nmo = orbs.shape
        # indeces of MOs that contain orbitals with (l,m) quantum numbers
        is_lm = np.array([False for imo in range(0, nmo)])
        for imo in range(0, nmo):
            # The MO with index imo is added to the list if its norm in the lm-subspace
            # exceeds a certain threshold.
            lm_norm = np.sum(abs(orbs[lm_indeces,imo])**2)
            if lm_norm > 1.0e-15:
                is_lm[imo] = True
        print "SELECTED ORBITALS"
        print "NMO(SELECTED) = %d   NMO(TOTAL) = %d" % (len(np.where(is_lm == True)[0]),nmo)
        f = self.dftb2.getOccupation()
        orbe = self.dftb2.getKSEnergies()
        # occupied and unoccupied lm-orbitals
        self.occupied_orbs = where((f > 0.1) & (is_lm == True))[0]
        self.virtual_orbs = where((f <= 0.1) & (is_lm == True))[0]
        # highest nr_active_occ occupied orbitals
        self.active_occupied_orbs = self.occupied_orbs[argsort(orbe[self.occupied_orbs])[::-1][:nr_active_occ][::-1]]
        # lowest nr_active_virt virtual orbitals
        self.active_virtual_orbs = self.virtual_orbs[argsort(orbe[self.virtual_orbs])[:nr_active_virt]]
        self.nocc = len(self.active_occupied_orbs)
        self.nvirt = len(self.active_virtual_orbs)
    def setActiveOrbitals(self, nr_active_occ=30, nr_active_virt=30):
        """ 
        Select nr_active_occ orbitals below and nr_active_virt orbitals 
        above the HOMO-LUMO gap which can contribute as single excitations.

        Parameters:
        ==========
        nr_active_occ: integer, number of active occupied orbitals
          if set to None, all occupied orbitals are selected 
        nr_active_virt: integer, number of active virtual orbitals
          if set to None, all virtual orbitals are included
        """
        f = self.dftb2.getOccupation()
        orbe = self.dftb2.getKSEnergies()
        self.occupied_orbs = where(f > 0.1)[0]
        self.virtual_orbs = where(f <= 0.1)[0]
        # highest nr_active_occ occupied orbitals
        self.active_occupied_orbs = self.occupied_orbs[argsort(orbe[self.occupied_orbs])[::-1][:nr_active_occ][::-1]]
        # lowest nr_active_virt virtual orbitals
        self.active_virtual_orbs = self.virtual_orbs[argsort(orbe[self.virtual_orbs])[:nr_active_virt]]
        self.nocc = len(self.active_occupied_orbs)
        self.nvirt = len(self.active_virtual_orbs)
#        if self.dftb2.verbose > 0:
#            print "occupied orbitals : %s" % self.occupied_orbs
#            print "virtual orbitals: %s" % self.virtual_orbs
#            print "active occupied orbitals : %s" % self.active_occupied_orbs
#            print "active virtual  orbitals : %s" % self.active_virtual_orbs
    def getActiveOrbitals(self):
        return self.active_occupied_orbs, self.active_virtual_orbs
    def determineActiveSpace(self, threshold=1.0e-5):
        """
        determine how large the active space has to be so that the lowest N states have no amplitudes
        in the inactive space
        """
        Nst = self.Cij.shape[0]
        occ_mins, virt_maxs = [], []
        for n in range(0, Nst):
            Cn = self.Cij[n,:,:] # amplitudes of TD-DFT 'wavefunction' for excitation i->a
            Cn2 = np.ravel(Cn*Cn)
            sort_indx = np.argsort(-Cn2)
            Cn2 = Cn2[sort_indx]
            Sn = np.cumsum(Cn2)
            # find excitations that contribute to 1-threshold of the total probability 
            active = sort_indx[np.where(Sn <= 1.0-threshold)]
            if len(active) == 0:
                if Sn[0] > 1.0-threshold:
                    # state is totally dominated by a single excitation
                    active = sort_indx[0]
            occ,virt = np.unravel_index(active, Cn.shape)
            # lowest occupied orbital and highest virtual in active space
            occ_min = self.active_occupied_orbs[occ.min()]
            virt_max = self.active_virtual_orbs[virt.max()]
            occ_mins.append(occ_min)
            virt_maxs.append(virt_max)
        # determine the number of occupied orbitals below the HOMO and the number of virtuals above
        # the LUMO that are necessary to describe the lowest Nst states
        occ_min = min(occ_mins)
        virt_max = max(virt_maxs)
        HOMO = self.occupied_orbs.max()
        nr_active_occ = HOMO - occ_min + 1
        nr_active_virt = virt_max - HOMO
        if self.dftb2.verbose > 0:
            print ""
            print "Active Space"
            print "============"
            for n in range(0, Nst):
                print "State %3d: contains excitations from window  occ %3d - virt %3d" % (n+1,occ_mins[n], virt_maxs[n])
            print ""
            print "The suggested active space for describing the lowest %d states consists of" % Nst
            print "the highest %d occupied and the lowest %d virtual orbitals." % (nr_active_occ, nr_active_virt)
            print "This active space accounts for %2.7f of the total probability for each state." % (1.0-threshold)
            print "You should check that the active space contains all neighbouring (almost) degenerate orbitals."
            print "Dimension of active space: %d" % (nr_active_occ * nr_active_virt)
            print "Dimension of  full  space: %d" % (len(self.occupied_orbs)*len(self.virtual_orbs))
            print ""
        # return the suggested active space
        return nr_active_occ, nr_active_virt
    def hasActiveSpace(self):
        """check whether an active space is set or whether the full space of excitations is used"""
        if (len(self.occupied_orbs) == len(self.active_occupied_orbs)) \
                and \
           (len(self.virtual_orbs) == len(self.active_virtual_orbs)):
            # full space
            return False
        else:
            # limited active space
            return True
    def setMultiplicity(self, multiplicity):
        assert multiplicity in ["S", "T"], "multiplicity has to be either 'S' for Singlet or 'T' for Triplet"
        self.multiplicity = multiplicity
    @T.timer
    def _MullikenTransitionCharges(self):
        """
        point charge approximation of transition densities according to formula (14)
        in Heringer, Niehaus  J Comput Chem 28: 2589-2601 (2007)
        """
        atomlist = self.dftb2.getGeometry()
        valorbs = self.dftb2.getValorbs()
        orbs = self.dftb2.getKSCoefficients()
        S = self.dftb2.getOverlapMatrix()
        """
        # PYTHON CODE
        # transition charges between occupied and virtual orbitals
        Nat = len(atomlist)
        dim_o = len(self.active_occupied_orbs)
        dim_v = len(self.active_virtual_orbs)
        self.qtrans_ov = zeros((Nat, dim_o, dim_v))
        # transition charges between occupied orbitals
        self.qtrans_oo = zeros((Nat, dim_o, dim_o))
        # transition charges between virtual orbitals
        self.qtrans_vv = zeros((Nat, dim_v, dim_v))
        
        Sc = dot(S,orbs)
        mu = 0
        # TODO: implement this in Fortran (using OpenMP parallelization)
        for A,(ZA,posA) in enumerate(atomlist):
            for (nA,lA,mA) in valorbs[ZA]:
                # occupied - virtuals
                for i,occi in enumerate(self.active_occupied_orbs):
                    for a,virta in enumerate(self.active_virtual_orbs):
                        self.qtrans_ov[A,i,a] += 0.5*(  orbs[mu,occi]*Sc[mu,virta] \
                                                   + orbs[mu,virta]*Sc[mu,occi])
                # occupied - occupied
                for i,occi in enumerate(self.active_occupied_orbs):
                    for j,occj in enumerate(self.active_occupied_orbs):
                        self.qtrans_oo[A,i,j] += 0.5*(  orbs[mu,occi]*Sc[mu,occj] \
                                                   + orbs[mu,occj]*Sc[mu,occi])
                # virtual - virtual
                for a,virta in enumerate(self.active_virtual_orbs):
                    for b,virtb in enumerate(self.active_virtual_orbs):
                        self.qtrans_vv[A,a,b] += 0.5*(  orbs[mu,virta]*Sc[mu,virtb] \
                                                   + orbs[mu,virtb]*Sc[mu,virta])         
                mu += 1
        #
        """
        # FASTER FORTRAN CODE
        from DFTB.extensions import tddftb
        orbs_occ = orbs[:,self.active_occupied_orbs]
        orbs_virt = orbs[:,self.active_virtual_orbs]
        q_oo, q_ov, q_vv = tddftb.tddftb.trans_charges(orbs_occ, orbs_virt, S, self.dftb2.orbsPerAtom)
        """
        # COMPARE PYTHON AND FORTRAN RESULTS
        err = np.sum(abs(q_oo - self.qtrans_oo))
        assert err < 1.0e-10
        err = np.sum(abs(q_ov - self.qtrans_ov))
        assert err < 1.0e-10
        err = np.sum(abs(q_vv - self.qtrans_vv))
        assert err < 1.0e-10
        #
        """
        self.qtrans_oo, self.qtrans_ov, self.qtrans_vv = q_oo, q_ov, q_vv
        if self.dftb2.mulliken_dipoles == 0:
            return self.qtrans_ov, self.qtrans_oo, self.qtrans_vv
        #
        ############ MULTIPOLES #################
        D = self.dftb2.D
        # transition multipoles between occupied and virtual orbitals
        Nat = len(atomlist)
        dim_o = len(self.active_occupied_orbs)
        dim_v = len(self.active_virtual_orbs)
        qm_ov = zeros((4*Nat, dim_o, dim_v))
        # transition multipoles between occupied orbitals
        qm_oo = zeros((4*Nat, dim_o, dim_o))
        # transition multipoles between virtual orbitals
        qm_vv = zeros((4*Nat, dim_v, dim_v))
        
        Sc = np.dot(S,orbs)
        Dc = np.tensordot(D, orbs, axes=(1,0))
        mu = 0
        # TODO: implement this in Fortran (using OpenMP parallelization)
        for A,(ZA,posA) in enumerate(atomlist):
            for (nA,lA,mA) in valorbs[ZA]:
                # occupied - virtuals
                for i,occi in enumerate(self.active_occupied_orbs):
                    for a,virta in enumerate(self.active_virtual_orbs):
                        # transition monopoles
                        qm_ov[4*A,i,a] += 0.5*(  orbs[mu,occi]*Sc[mu,virta] \
                                                   + orbs[mu,virta]*Sc[mu,occi])
                        # transition dipoles
                        qm_ov[4*A+1:4*(A+1),i,a] += 0.5*(  orbs[mu,occi]*Dc[mu,:,virta] \
                                                   + orbs[mu,virta]*Dc[mu,:,occi])
                # occupied - occupied
                for i,occi in enumerate(self.active_occupied_orbs):
                    for j,occj in enumerate(self.active_occupied_orbs):
                        # transition monopoles
                        qm_oo[4*A,i,j] += 0.5*(  orbs[mu,occi]*Sc[mu,occj] \
                                                   + orbs[mu,occj]*Sc[mu,occi])
                        # transition dipoles
                        qm_oo[4*A+1:4*(A+1),i,j] += 0.5*(  orbs[mu,occi]*Dc[mu,:,occj] \
                                                   + orbs[mu,occj]*Dc[mu,:,occi])
                # virtual - virtual
                for a,virta in enumerate(self.active_virtual_orbs):
                    for b,virtb in enumerate(self.active_virtual_orbs):
                        # transition monopoles
                        qm_vv[4*A,a,b] += 0.5*(  orbs[mu,virta]*Sc[mu,virtb] \
                                                   + orbs[mu,virtb]*Sc[mu,virta])         
                        # transition dipoles
                        qm_vv[4*A+1:4*(A+1),a,b] += 0.5*(  orbs[mu,virta]*Dc[mu,:,virtb] \
                                                   + orbs[mu,virtb]*Dc[mu,:,virta])      
        
                mu += 1
        #
        self.qtrans_ov_multipoles = qm_ov
        self.qtrans_oo_multipoles = qm_oo
        self.qtrans_vv_multipoles = qm_vv

        #########################################

        return self.qtrans_ov, self.qtrans_oo, self.qtrans_vv
    def getTransitionCharges(self):
        return self.qtrans_oo, self.qtrans_vv, self.qtrans_ov
    def _getOrbitalEnDifferences(self):
        """
        omega_ia = en_a - en_i 
        """
        orbe = self.dftb2.getKSEnergies()
        # energy differences between occupied and virtual Kohn-Sham orbitals
        omega = zeros((self.nocc, self.nvirt))
        for i,occi in enumerate(self.active_occupied_orbs):
            for a,virta in enumerate(self.active_virtual_orbs):
                omega[i,a] = orbe[virta] - orbe[occi]
        return omega
    def _getOrbitalOccDifferences(self):
        """
        f_ia = f_i - f_a
        """
        f = self.dftb2.getOccupation()
        # occupation differences between occupied and virtual Kohn-Sham orbitals
        df = zeros((self.nocc,self.nvirt))
        for i,occi in enumerate(self.active_occupied_orbs):
            for a,virta in enumerate(self.active_virtual_orbs):
                df[i,a] = f[occi] - f[virta]
        return df
    #############################################
    # TDDFT linear response equations           #
    #############################################
    @T.timer
    def _solve_lr_equations(self, response_method="Casida", multiplicity="Singlet", nstates=None, diag_ifact=1, diag_conv=1.0e-10, diag_maxiter=10, diag_check=1, diag_L2threshold=0.0, diag_selected_irreps=[], ct_correction=1):
        ###
        # transition charges
        if self.dftb2.verbose > 0:
            print "compute transition charges"
        qtrans_ov, qtrans_oo, qtrans_vv = self._MullikenTransitionCharges()
        # compute Oia for Lambda diagnostics Oia = int int phi_i(r)^2 phi_a(r)^2 dr dr
        self._Lambda2_calcOia()
        # gamma matrices at short and long range
        if self.dftb2.verbose > 0:
            print "build gamma matrices"
        gamma = self.dftb2.getGamma()
        if self.dftb2.long_range_correction == 1:
            gamma_lr = self.dftb2.getGamma_lr()
        else:
            gamma_lr = 0.0*gamma
        ###
        # KS orbital energy differences
        self.omega = self._getOrbitalEnDifferences()
        ### long-range CT correction shift orbital energies
        # ea - ei up, if the overlap between i and a vanishes
        self.ct_correction = ct_correction
        if ct_correction == 1:
            assert self.dftb2.long_range_correction == 0, "CT correction --ct_correction=1 cannot be combined with long-range correction --long_range_correction=1!"
            CTcorrection = LongRangeCTCorrection(self)
            omega_shift = CTcorrection.diagonal_shifts(self.omega)
        else:
            omega_shift = 0*self.omega
        #
        # occupation number differences
        df = self._getOrbitalOccDifferences()
        if response_method == "Casida":
            if nstates == None:
                if self.dftb2.verbose > 0:
                    print "construct and diagonalize full TD-DFTB matrix with dimension %d" % (self.nocc*self.nvirt)
                # construct the full TD-DFTB matrix in memory
                # and diagonalize the Hermitian eigenvalue problem
                self.Omega, self.Cij, self.XmY, self.XpY = \
                    Solver.Casida(gamma, gamma_lr,\
                                  qtrans_oo, qtrans_vv, qtrans_ov,\
                                  self.omega, omega_shift, df, self.nocc, self.nvirt)
            else:
                # solve the eigenvalue problem
                # for the lowest nstates using an iterative
                # approach 
                if self.dftb2.verbose > 0:
                    print "solve eigenvalue problem iteratively"
                #
                if self.dftb2.use_symmetry > 0 and diag_selected_irreps != None:
                    selector = SymmetryAssignmentEx(self, self.dftb2.symmetry_group)
                    selector.select_irreps(diag_selected_irreps)
                else:
                    selector = None


                if self.dftb2.long_range_correction == 0:
                    ######### MULTIPOLES #########
                    if self.dftb2.mulliken_dipoles == 1:
                        qtrans_ov = self.qtrans_ov_multipoles
                        gamma = self.dftb2.gamma_multipoles
                    ##############################
                    self.Omega, self.Cij, self.XmY, self.XpY = \
                        Solver.HermitianDavidson(gamma, qtrans_ov, \
                                  self.omega, omega_shift, self.nocc, self.nvirt, \
                                  self.XmY, self.XpY, \
                                  self.Oia, \
                                  selector=selector, \
                                  nstates=nstates, ifact=diag_ifact, \
                                  conv=diag_conv, maxiter=diag_maxiter, \
                                  L2_thresh=diag_L2threshold, \
                                  verbose=self.dftb2.verbose)
                #
                else:
                    self.Omega, self.Cij, self.XmY, self.XpY = \
                        Solver.nonHermitianDavidson(gamma, gamma_lr,\
                                  qtrans_oo, qtrans_vv, qtrans_ov,\
                                  self.omega, self.nocc, self.nvirt, \
                                  self.XmY, self.XpY, \
                                  selector=selector,\
                                  nstates=nstates, ifact=diag_ifact, \
                                  conv=diag_conv, maxiter=diag_maxiter, \
                                  lc=self.dftb2.long_range_correction, \
                                  verbose=self.dftb2.verbose)
                if diag_check == 1:
                    ### DEBUG
                    print "CHECK ITERATIVE DIAGONALIZATION"
                    # only for debugging purposes, check that Davidson-like
                    # diagonalization and full diagonalization agree
                    # check that Davidson-like diagonalization yields the same results
                    Omega_full, C_full, XmY_full, XpY_full = \
                        Solver.Casida(gamma, gamma_lr,\
                                  qtrans_oo, qtrans_vv, qtrans_ov,\
                                  self.omega, omega_shift, df, self.nocc, self.nvirt)

                    err = np.sum(abs(self.Omega - Omega_full[:nstates]))
                    assert err < 1.0e-10, "eigenvalues differ, err = %s\nEn(iter) = %s\nEn(full) = %s" % (err, self.Omega, Omega_full[:nstates])
                    print "EIGENVALUES AGREE"
                    err = 0
                    for n in range(0, nstates):
                        # arbitrary global sign
                        XmYerr = min(np.sum(abs(self.XmY[n,:,:] - XmY_full[n,:,:])),
                                     np.sum(abs(self.XmY[n,:,:] + XmY_full[n,:,:])))
                        XpYerr = min(np.sum(abs(self.XpY[n,:,:] - XpY_full[n,:,:])),
                                     np.sum(abs(self.XpY[n,:,:] + XpY_full[n,:,:])))
                        Cerr = min(np.sum(abs(self.Cij[n,:,:] - C_full[n,:,:])),
                                   np.sum(abs(self.Cij[n,:,:] + C_full[n,:,:])))
                        err += XmYerr + XpYerr # + Cerr
                    assert err < 1.0e-10, "Davidson-like diagonalization and full diagonalization disagree, error = %s! You might want to tighten the convergence using --diag_conv=1.0e-13. Symmetry should be disabled." % err
                    print "EIGENVECTORS AGREE"
                    ####

        elif response_method == "Tamm-Dancoff":
            self.Omega, self.Cij = \
                Solver.TDA(gamma, gamma_lr,\
                       qtrans_oo, qtrans_vv, qtrans_ov,\
                       self.omega, df, self.nocc, self.nvirt)
        else:
            raise Exception("response_method should be \"Casida\" or \"Tamm-Dancoff\"")
    def getXY(self):
        if hasattr(self, "XmY"):
            return self.XmY, self.XpY
        else:
            raise NotImplemented("X-Y and X+Y are only available if the full TD-DFT equations are solved (Casida's equation)")
    def getExcEnergies(self):
        return self.Omega
    ########################################################
    # Lambda diagnostics
    ########################################################
    def _Lambda2_calcOia(self):
        """
        see J. Chem. Phys. 128, 044118 2008

        assumes that the member variables self.qtrans_oo and self.qtrans_vv 
        have been set
        """
        gaussOmega = self.dftb2.getGaussianOverlap()
#        print "GaussOmega"
#        print "=========="
#        print gaussOmega
        atomlist = self.dftb2.getGeometry()
        Nat = len(atomlist)
        dim_o = len(self.active_occupied_orbs)
        dim_v = len(self.active_virtual_orbs)
        Oia = zeros((dim_o, dim_v))
        for i,occi in enumerate(self.active_occupied_orbs):
            Oii = dot(self.qtrans_oo[:,i,i], dot(gaussOmega, self.qtrans_oo[:,i,i]))
            for a,virta in enumerate(self.active_virtual_orbs):
                Oaa = dot(self.qtrans_vv[:,a,a], dot(gaussOmega, self.qtrans_vv[:,a,a]))
                Oia[i,a] = dot(self.qtrans_oo[:,i,i], dot(gaussOmega, self.qtrans_vv[:,a,a]))
                Oia[i,a] /= sqrt(Oii*Oaa)
        self.Oia = Oia
    def _Lambda2_diagnostics(self):
        # compute Lambda2
        # assumes Oia has been calculated
        self.Lambda2 = sum(sum(pow(self.Cij,2) * self.Oia, axis=-1), axis=-1)
        return self.Lambda2
#
    def _TransitionDipolesMulliken(self):
        atomlist = self.dftb2.getGeometry()
        fvec = np.zeros((3,len(self.Omega)))
        for A,(ZA,posA) in enumerate(atomlist):
            tdipmom = np.tensordot(self.qtrans_ov[A,:,:],self.Cij,axes=([0,1],[1,2]))
            fvec += np.outer(posA, tdipmom)
        self.tdipX, self.tdipY, self.tdipZ = fvec[0,:], fvec[1,:], fvec[2,:]
    def _TransitionDipoles(self):
        DipMat = self.dftb2.getDipoleMatrix()
#        print DipMat
        orbs = self.dftb2.getKSCoefficients()
        # transition dipole moments between KS molecular orbitals
        self.tdipMO_X = dot(orbs.transpose(), dot(DipMat[:,:,0], orbs))
        self.tdipMO_Y = dot(orbs.transpose(), dot(DipMat[:,:,1], orbs))
        self.tdipMO_Z = dot(orbs.transpose(), dot(DipMat[:,:,2], orbs))
        # select transition dipole matrix elements between active orbitals
        actocc, actvirt = meshgrid(self.active_virtual_orbs, self.active_occupied_orbs)
        tdipActive_X = self.tdipMO_X[actocc, actvirt]
        tdipActive_Y = self.tdipMO_Y[actocc, actvirt]
        tdipActive_Z = self.tdipMO_Z[actocc, actvirt]
        # transition dipole matrix elements between ground state and each excited state
        # t[I] =  <Psi0|x|PsiI>
        self.tdipX = np.tensordot(self.Cij, tdipActive_X, axes=([1,2],[0,1]))
        self.tdipY = np.tensordot(self.Cij, tdipActive_Y, axes=([1,2],[0,1]))
        self.tdipZ = np.tensordot(self.Cij, tdipActive_Z, axes=([1,2],[0,1]))
        return (self.tdipX, self.tdipY, self.tdipZ)
    def getTransitionDipolesExc(self, Nst):
        """
        compute the transition dipoles between the lowest Nst states
        (including the ground state)
        """
        # combine x,y,z components of transition dipoles between active MOs
        occ, virt = self.active_occupied_orbs, self.active_virtual_orbs
        nocc, nvirt = len(occ), len(virt)
        # occupied-occupied  <i|r|j>
        tdip_oo = np.zeros((nocc,nocc,3))

        tdip_oo[:,:,0] = self.tdipMO_X[occ,:][:,occ]
        tdip_oo[:,:,1] = self.tdipMO_Y[occ,:][:,occ]
        tdip_oo[:,:,2] = self.tdipMO_Z[occ,:][:,occ]
        # occupied - virtual  <i|r|a>
        tdip_ov = np.zeros((nocc,nvirt,3))
        tdip_ov[:,:,0] = self.tdipMO_X[occ,:][:,virt]
        tdip_ov[:,:,1] = self.tdipMO_Y[occ,:][:,virt]
        tdip_ov[:,:,2] = self.tdipMO_Z[occ,:][:,virt]
        # virtual - virtual  <a|r|b>
        tdip_vv = np.zeros((nvirt,nvirt,3))
        tdip_vv[:,:,0] = self.tdipMO_X[virt,:][:,virt]
        tdip_vv[:,:,1] = self.tdipMO_Y[virt,:][:,virt]
        tdip_vv[:,:,2] = self.tdipMO_Z[virt,:][:,virt]
        
        # number of electrons
        N = self.dftb2.Nelec_val
        assert N % 2 == 0
        # ground state dipole
        # d0 = sum_k <k|r|k>
        d0 = np.trace(tdip_oo, axis1=0, axis2=1)
        tdip = np.zeros((Nst, Nst, 3))
        # excitation coefficients C[A,i,a]
        C = self.Cij
        # transition dipoles between ground and excited states
        for A in range(1, Nst):
            tdip[0,A,:] = np.tensordot(C[A-1,:,:], tdip_ov, axes=([0,1],[0,1]))
            tdip[A,0,:] = tdip[0,A,:]
            # check that we get the same result
            assert np.sum(abs(tdip[0,A,0] - self.tdipX[A-1])) < 1.0e-10
            assert np.sum(abs(tdip[0,A,1] - self.tdipY[A-1])) < 1.0e-10
            assert np.sum(abs(tdip[0,A,2] - self.tdipZ[A-1])) < 1.0e-10
        # 
        # transition dipoles between excited states
        for A in range(1,Nst):
            for B in range(A, Nst):
                # P0 = sum_i sum_a C_(ia)^A C_(ia)^B
                P0 = np.tensordot(C[A-1,:,:], C[B-1,:,:], axes=([0,1],[0,1]))
                # Pij[i,j] = sum_a C_(ia)^A C_(ja)^B
                Pij = np.tensordot(C[A-1,:,:], C[B-1,:,:], axes=(1,1))
                # Pab[a,b] = sum_i C_(ia)^A C_(ib)^B
                Pab = np.tensordot(C[A-1,:,:], C[B-1,:,:], axes=(0,0))
                tdip[A,B,:] = d0 * P0 \
                                  - np.tensordot(Pij, tdip_oo, axes=([0,1],[0,1])) \
                                  + np.tensordot(Pab, tdip_vv, axes=([0,1],[0,1]))
                tdip[B,A,:] = tdip[A,B,:]
        return tdip
    def getOscillatorStrengths(self):
        # UGLY, overwrites Mulliken approximation
        if self.multiplicity == "S":
            self.oscillator_strength = 4.0/3.0 * self.Omega * (abs(self.tdipX)**2 + abs(self.tdipY)**2 + abs(self.tdipZ)**2)
        else:
            # triplet states are dark
            self.oscillator_strength = 0.0 * self.Omega
        return self.oscillator_strength
    def _dominant_excitations(self, I, dominant_probability=0.99):
        """
        Find those single-particle excitations that together make up for
        at least 99% of the probability of state I

        Parameters:
        ===========
        I: index to excited state (starting from 0 ~ 1st excited state)
        dominant_probability:   sum_(dominant ij) |C^ij|^2 > dominant_probability

        Returns:
        ========
        ex_occ : list of occupied orbitals
        ex_virt: list of virtual orbitals
        C: excitation coefficients

        The i-th most important excitation is ex_occ[i] -> ex_virt[i] with strength |C[i]|^2
        """
        # sort excitations by contribution in decending order
        iocc_sort, avirt_sort = utils.argsort_2d(-abs(self.Cij[I,:,:])**2)
        prob_sum = 0.0
        ex_occ = []
        ex_virt = []
        C = []
        for (i,a) in zip(iocc_sort, avirt_sort):
            # add excitations until threshold is reached
            if prob_sum >= dominant_probability:
                break
            dp = abs(self.Cij[I,i,a])**2
            if dp > 1.0e-10:
                ex_occ.append(self.active_occupied_orbs[i])
                ex_virt.append(self.active_virtual_orbs[a])
                C.append(self.Cij[I,i,a])
            prob_sum += dp
        return ex_occ, ex_virt, C
    def _SymmetryOfExcitations(self):
        """
        try to determine the symmetry of the excited states
        """
        self.Irreps = ["" for I in range(0, len(self.Omega))]
        if self.dftb2.use_symmetry > 0:
            print "Assigning symmetries to excited states based on transformation properties of the transition densities"

            symas = SymmetryAssignmentEx(self, self.dftb2.symmetry_group)
            # labels of irreducible representations for each state
            if self.nstates == None:
                sym_max_states = self.dftb2.use_symmetry
            else:
                # If only the lowest few roots are calculated, determine
                # the symmetry of all of them 
                sym_max_states = self.nstates
            for I in range(0, min(sym_max_states, len(self.Omega))):
                self.Irreps[I] = symas.assign_symmetry(I)
    def _writeExcitations(self, print_max_states=1000000, nr_dominant_ex=2, print_MO_tdip=False):
        """
        write out excitation energies, oscillator strengths and the <nr_dominant_ex> dominant single particle
        excitations for the lowest <max_states> excited states.
        """
        import string
        if self.dftb2.use_symmetry > 0:
            symstr = self.dftb2.symmetry_group.name().rjust(3)
        else:
            symstr = "   "
        txt = " ground state energy: %s %s\n\n" % (string.rjust("%.7f hartree" % self.en0,15), \
                               string.rjust("%.7f eV" % (self.en0*hartree_to_eV), 15))
        txt += "           Excited States\n"
        txt += "           ==============\n"
        txt += "     Spn %s     exc. en. /hartree      exc. en. / eV   exc. en. / nm    osc. strength             dominant contrib.               en_KS / eV      transition dipole [TDx TDy TDz]   /\\ diagn.\n" % symstr
        for I in range(0, min(print_max_states, len(self.Omega))):
            # most dominant excitations
            domex = self._dominant_excitations(I,dominant_probability=0.99)
            domex_str = 7*" "
            domKS_str = 6*" "
            for count,(iocc,avirt,Cia) in enumerate(zip(*domex)):
                if count == 0:
                    orbe = self.dftb2.getKSEnergies()
                    # differences of KS orbital energies for dominant excitation
                    omegaKS = orbe[avirt] - orbe[iocc]
                    domKS_str = "%2.4f" % (omegaKS*hartree_to_eV)
                if print_MO_tdip == True:
                    # print transition dipoles <o|r|v> for dominant excitations
                    MOtdip = ",[%+2.7f %+2.7f %+2.7f]" % (self.tdipMO_X[iocc,avirt], self.tdipMO_Y[iocc,avirt], self.tdipMO_Z[iocc,avirt])
                else:
                    MOtdip = ""
                domex_str += "%d->%d(%+.3e%s)" % (iocc, avirt, Cia, MOtdip)
                if count < min(nr_dominant_ex, len(domex[0])) - 1:
                    domex_str += ","
                if count == nr_dominant_ex-1:
                    break
            txt += "%s %s %s: %s %s %s  %s %s %s %s %s\n" % (string.rjust(str(I+1),5), \
                               self.multiplicity, \
                               string.rjust(self.Irreps[I], 3), \
                               string.rjust("%.7f hartree" % self.Omega[I],20), \
                               string.rjust("%.7f eV" % (self.Omega[I]*hartree_to_eV), 17), \
                               string.rjust("%.7f nm" % (hartree_to_nm / self.Omega[I]), 17), \
                               string.rjust("%.7f" % self.oscillator_strength[I], 12), \
                               string.rjust("%s" % domex_str, 45), \
                               string.rjust("%s" % domKS_str, 8), \
                               string.rjust("[%+2.7f %+2.7f %+2.7f]" % (self.tdipX[I], self.tdipY[I], self.tdipZ[I]), 40), \
                               string.rjust("%.4f" % self.Lambda2[I], 7))
        print txt
    def saveAbsorptionSpectrum(self, spectrum_file="absorption_spectrum.dat"):
        """
        Write a table with absorption spectrum to a file.
        For each excitation energy the oscillator strength is given.
        
        Parameters:
        ==========
        Output.spectrum_file: write table with excitation energies and oscillator strengths to this file. If the states were classified by symmetry a seperate file is written for each irrep.
        """
        import string
        max_states = 10000000
        unique_irreps = np.unique(self.Irreps)
        # If the states were classified by symmetry a separate spectrum is written for each irrep
        for irrep in unique_irreps:
            spectrum_file_irrep = path.expanduser(path.expandvars(spectrum_file))
            if irrep != "":
                spectrum_file_irrep += ".%s" % irrep
            fh = open(spectrum_file_irrep, 'w')
            if irrep != "":
                print>>fh, "# states with %s symmetry" % irrep
            print>>fh, "# Absorption spectrum"
            print>>fh, "# exc. energy / hartree  osc.strength"
            for I in range(0, min(max_states, len(self.Omega))):
                if self.Irreps[I] == irrep:
                    print>>fh, "%s   %s" % (string.rjust("%.7f" % self.Omega[I],20), string.rjust("%.7f" % self.oscillator_strength[I],20))
            fh.close()
    def _getExcitationCoefficients(self):
        """
        Excited states are linear combination of single excitations (d^+ creation and d annihilation operator):
           Phi^I = sum_ia C^I_ia da^+ di Phi^0

        Returns:
        ========
        C: C[I,i,a] is the coefficient for the (occupied i)->(virtual a) excitation in excited state I 
        """
        return self.Cij
    def TransitionDensityMatrix(self, I):
        """
        compute transition density matrix for state I according to

        P^I_(a,b) = sum_(o) sum_(v) C[I,o,v] * orbs[a,o] * orbs[b,v]
        """
        orbs_occ  = self.dftb2.orbs[:,self.active_occupied_orbs]
        orbs_virt =self.dftb2.orbs[:,self.active_virtual_orbs]
        PtransI = np.dot(orbs_occ, np.dot(self.Cij[I,:,:], orbs_virt.transpose()))
        return PtransI
    def ParticleHoleCharges(self, I):
        """
        Compute the density difference between the excited state I and the ground
        state coarse grained to the atomic positions.
        """
        atomlist = self.dftb2.getGeometry()
        Nat = len(atomlist)
        particle_charges = np.zeros(Nat)
        hole_charges = np.zeros(Nat)
        for A in range(0, Nat):
            for o in range(0, len(self.active_occupied_orbs)):
                particle_charges[A] += np.dot(self.Cij[I,o,:], np.dot(self.qtrans_vv[A,:,:], self.Cij[I,o,:]))
            for v in range(0, len(self.active_virtual_orbs)):
                hole_charges[A] -= np.dot(self.Cij[I,:,v], np.dot(self.qtrans_oo[A,:,:], self.Cij[I,:,v]))
        # distance between particle and hole
        p_pos = np.zeros(3)
        h_pos = np.zeros(3)
        for A,(ZA,posA) in enumerate(atomlist):
            posA = np.array(posA)
            p_pos += posA * particle_charges[A]
            h_pos += posA * hole_charges[A]
        p_pos /= float(np.sum(particle_charges))
        h_pos /= float(np.sum(hole_charges))
        ph_dist = la.norm(p_pos - h_pos)
        if self.dftb2.verbose > 0:
            txt = ""
            txt += "Exciton: State %d\n" % I
            txt += "===================\n"
            txt += "Particle position: %s bohr\n" % p_pos
            txt += "Hole     position: %s bohr\n" % h_pos
            txt += "P-H distance     : %3.7f bohr\n" % ph_dist
            txt += "Particle-Hole Charges\n"
            txt += "======================\n"
            for A,(ZA,posA) in enumerate(atomlist):
                txt += "%s: %3.7f (p) \t %3.7f (h)\n" % (string.center("%s-%s" % (atom_names[ZA-1], A),12), particle_charges[A], hole_charges[A])
            print txt

        return particle_charges, hole_charges
    def analyseParticleHole(self, particle_hole_states="[]", particle_hole_dir=None):
        """
        Charge-Transfer-Analysis.particle_hole_states: list of excited states (starting from 1) for which the positions of the particle and the hole should be determined.
        Charge-Transfer-Analysis.particle_hole_dir: The charge density of the particle and the hole are written to this file. The density is coarse grained to atomic resolution, so for each atom the portion of the delocalized particle or hole sitting on that atom is save. The atoms are listed in the same order as in the xyz file.
        """
        for I in particle_hole_states:
            particle_charges, hole_charges = self.ParticleHoleCharges(I-1)
            data = np.vstack((particle_charges, hole_charges)).transpose()
            if particle_hole_dir != None:
                ph_file = path.join(path.expanduser(path.expandvars(particle_hole_dir)), "particle_hole_charges_%d.chrg" % I)
                np.savetxt(ph_file, data)
    # GRAPHICAL ANALYSIS
    def graphical_analysis(self, graphical=0):
        """
        Graphical-Analysis.graphical: If graphical > 0, a graphical user interface is opened for analysing the results of the calculations. Requires mayavi and a running X-server.
        """
        if graphical == 0:
            return
        print "STARTING GRAPHICAL USER INTERFACE"
        from DFTB.Analyse.mayavi import GUI
        GUI.start_graphical_analysis(self)

def main(xyz_file, conf_file):
    from DFTB.XYZ import read_xyz, extract_keywords_xyz
    from DFTB.DFTB2 import DFTB2
    from DFTB.Analyse.Cube import CubeExporterEx
    from DFTB.Molden import MoldenExporter

    outputfile = open("scf_output_dftb.txt", "a")
    sys.stdout = outputfile

    usage = "Usage: %s <xyz-file>\n" % xyz_file
    usage += "   --help option will give more information\n"

    parser = utils.OptionParserFuncWrapper(conf_file, [ \
        DFTB2.__init__, DFTB2.runSCC, \
        LR_TDDFTB.getEnergies, LR_TDDFTB.saveAbsorptionSpectrum, LR_TDDFTB.analyseParticleHole, \
        LR_TDDFTB.graphical_analysis, \
        CubeExporterEx.exportCubes, MoldenExporter.export, \
        Gradients.getGradients], \
        usage)

    atomlist = read_xyz(xyz_file)[0]
    kwds = extract_keywords_xyz(xyz_file)

    options = parser.options
    scf_options = parser.scf_options

    init_options = parser.init_options
    tddftb = LR_TDDFTB(atomlist, **init_options)

    options = parser.tde_options
    options.update(scf_options)
    tddftb.setGeometry(atomlist, charge=kwds.get("charge", 0.0))
    x = tddftb.getEnergies(**options)

    grad = Gradients(tddftb)
    options = parser.gradient_options
    grad.getGradients(**options)

    return str(x[-1])  # return energies in Hartree





if __name__ == "__main__":
    main(sys.argv[1], "dftbaby.cfg")