[![Build Status](https://img.shields.io/travis/com/AKEngels/CAST?style=plastic)](https://travis-ci.com/AKEngels/CAST) 
[![Coverage Status](https://img.shields.io/coveralls/github/AKEngels/CAST?style=plastic)](https://coveralls.io/github/AKEngels/CAST?branch=devel&service=github) (for branch devel)



# CAST (Conformational Analysis and Search Tool)
This is the official repository for **CAST** (**C**onformational **A**nalysis and **S**earch **T**ool).

The main feature of CAST is the calculation of energies and gradients of molecular systems from several energy interfaces. Those energies and gradients then are used for different tasks. Those tasks include local and global optimization (by **Monte Carlo** or **TabuSearch**), molecular dynamics simulations (**MD**) and based on them calculation of Entropies and Free Energie differences (**FEP** for alchemical transformations and **Umbrella Sampling** for geometrical transformations). Further features of CAST include searching for reaction paths (e. g. by **NEB**) or adding water shells to a given molecule. For a detailed description of all tasks see our [manual](https://github.com/AKEngels/CAST/blob/devel/manual/castmanual.tex). As source for energy and gradients for all tasks you can either choose one of the internal interfaces (i. e. the **forcefields** OPLSAA, AMBER, CHARMM and AMOEBA). Or you can use an external interface which means an external program is run and the energetic information is received from the output of this program. External programs with which CAST can communicate this way are for example the **quantumchemical** programs Gaussian, Psi4 and Orca or the **semiempirical** programs DFTB+ and Mopac (a complete list can also be found in our manual). It is also possible to choose a **QM/MM** interface which means that two or three of these interfaces are combined either by an additive or a subtractive QM/MM scheme.

## Installation

Please look at the installation instructions [here](https://github.com/AKEngels/CAST/wiki/How-to-build-CAST).

## First steps

For first steps please look at our [tutorial](https://github.com/AKEngels/CAST/blob/devel/manual/Tutorial/LaTeX/Tutorial.tex).

## Further resources

* [Tutorial](https://github.com/AKEngels/CAST/blob/devel/manual/Tutorial/LaTeX/Tutorial.tex) (LaTeX file)
* [Manual](https://github.com/AKEngels/CAST/blob/devel/manual/castmanual.tex) (LaTeX file)
* [Wiki](https://github.com/AKEngels/CAST/wiki)

## For developers

Please have a look at our [contributing guidelines](https://github.com/AKEngels/CAST/blob/devel/CONTRIBUTING.md)


