[![Build Status](https://travis-ci.com/AKEngels/CAST.svg?branch=devel)](https://travis-ci.com/AKEngels/CAST) 
[![Coverage Status](https://coveralls.io/repos/github/AKEngels/CAST/badge.svg?branch=devel&service=github)](https://coveralls.io/github/AKEngels/CAST?branch=devel&service=github) (for branch devel)


# CAST (Conformational Analysis and Search Tool)
This is the official repository for **CAST** (**C**onformational **A**nalysis and **S**earch **T**ool).

The main feature of CAST is the calculation of energies and gradients of molecular systems from several energy interfaces. Those energies and gradients then are used for different tasks. Those tasks include local and global optimization (by Monte Carlo or TabuSearch), molecular dynamics simulations and based on them calculation of Entropies and Free Energie differences (FEP for alchemical transformations and Umbrella Sampling for geometrical transformations). Further features of CAST include searching for reaction paths or adding water shells to a given molecule. For a detailed description of all tasks see our manual. As source for energy and gradients for all tasks you can either choose one of the internal interfaces (i. e. the forcefields OPLSAA, AMBER, CHARMM and AMOEBA). Or you can use an external interface which means an external program is run and the energetic information is received from the output of this program. External programs with which CAST can communicate this way are for example Gaussian, Psi4, Orca, DFTB+ or Mopac (a complete list can also be found in our manual). It is also possible to choose a QM/MM interface which means that two or three of these interfaces are combined either by an additive or a subtractive QM/MM scheme.

## Installation

### Windows

Prerequisites:
* Git, e.g. [Git for Windows](https://gitforwindows.org/)
* [Premake5](https://premake.github.io/) (yes, the alpha version)
* [Visual Studio](https://visualstudio.microsoft.com/de/downloads/)

Installation:
* open Git CMD, go into the folder where you want to download the CAST folder
* clone the repositoryby typing ``git clone https://github.com/AKEngels/CAST.git``
* type ``git submodule init`` and ``git submodule update`` to get the submodules Eigen and Googletest
* type ``cd optional_files`` to go into that folder
* copy the ``premak5.exe`` into this folder (or add the folder where this exe is to the path variable)
* type ``premake5 vs2017`` to create a visual studio project
* in Windows Explorer go into folder ``optional_files/project`` and double click on ``CAST.sln`` to open CAST in Visual Studio
* choose configuration "Release" and your system architecture (probably "x64"), then compile CAST (in German version of Visual Studio: Erstellen -> CAST erstellen)
* the program is called ``CAST_win_x64_release.exe`` and dropped into folder ``optional_files/build``

### Linux


Prerequisites:
* Git (normally already installed)
* [Premake5](https://premake.github.io/) (yes, the alpha version)
* g++ compiler (normally already installed)

Installation:
* open terminal, go into the folder where you want to download the CAST folder
* clone the repositoryby typing ``git clone https://github.com/AKEngels/CAST.git``
* type ``git submodule init`` and ``git submodule update`` to get the submodules Eigen and Googletest
* type ``cd optional_files`` to go into that folder
* type ``premake5 gmake`` to create a makefile (for this premake needs to be in your path environment variable)
* type ``cd project`` to enter project folder
* type ``make config=release_x64`` to compile CAST
* the program is called ``CAST_win_x64_release`` and dropped into folder ``optional_files/build``

### Different configurations

Besides the release configuration there are also a few others:

* The debug configuration is for debugging puposes. It is only recommended to use this if you want to develop in CAST.
* The configurations which start with "Armadillo_" include the [Armadillo C++ linear algebra library]( http://arma.sourceforge.net/ ) which will speed up matrix operations considerably. This enhances the performance of the tasks ENTROPY and PCA specifically. To enable armadillo, compile CAST with the flag USE_ARMADILLO enabled. You (yes, you yourself) will have to make sure that the necessary include-files are included and the lapack / blas linear algebra libraries are found by the compiler. You might get some clues from the makefile which is found in the optional_files/project/ folder. While CAST officially supports armadillo and you should report bugs occuring in these routines, we are afraid we cannot offer support on how to compile CAST with external dependencies enabled. But don't worry, it is not complicated! On the [armadillo webpage](http://arma.sourceforge.net/ ) comprehensive information on how to compile programs with armadillo may be found.  Good luck :)
* The configurations that start with "Python_" include python functions by the python API. When this is enabled the plotting features for tasks MD and FEP are activated. Furthermore the energy interface to DFTBaby is implemented with the python API so this interface is only available in these configurations. In order to be able to use it Python 2.7 with all necessary modules has to be installed on you machine. Information which modules that are can be found in the manual.
* Configurations that include "Testing" are testing configurations so they can't be used to run calculations but only run the tests. In order to run the tests just compile this configuration and then use a terminal to go into the build-folder where you run the program ``CAST_win_x64_testing``. If everything in the output is green all is fine. **Attention!** Tests do not work if there is a CAST inputfile in the folder where you run the tests!

## First steps

For first steps please look at our [tutorial](https://github.com/AKEngels/CAST/blob/devel/manual/Tutorial/LaTeX/Tutorial.tex).

## Further resources

* [Tutorial](https://github.com/AKEngels/CAST/blob/devel/manual/Tutorial/LaTeX/Tutorial.tex)
* [Manual](https://github.com/AKEngels/CAST/blob/devel/manual/castmanual.tex)
* [Wiki](https://github.com/AKEngels/CAST/wiki)




