[![Build Status](https://travis-ci.com/AKEngels/CAST.svg?branch=devel)](https://travis-ci.com/AKEngels/CAST) (for branch devel)

# CAST (Conformational Analysis and Search Tool)
This is the official repository for the
CAST
Conformational Analysis and Search Tool.

CAST reads its runtime options from the inputfile "CAST.txt". CAST's features are split into tasks that are outlined in the manual. It depends on the Eigen library that can be installed by using git submodules. If you want to run the tests googletest, another submodule, also has to be included.

Optionally, CAST is fit to include the "Armadillo C++ linear algebra library" ( http://arma.sourceforge.net/ ) which will speed up matrix operations considerably. This enhances the performance of the tasks ENTROPY and PCA specifically. To enable armadillo, compile CAST with the flag USE_ARMADILLO enabled. You (yes, you yourself) will have to make sure that the necessary include-files are included and the lapack / blas linear algebra libraries are found by the compiler. You might get some clues from the makefile which is found in the optional_files/project/ folder. While CAST officially supports armadillo and you should report bugs occuring in these routines, we are afraid we cannot offer support on how to compile CAST with external dependencies enabled. But don't worry, it is not complicated! On the armadillo webpage ( http://arma.sourceforge.net/ ) comprehensive information on how to compile programs with armadillo may be found.  Good luck :)

Another optional dependency is Python. You can find information how to work with it in our wiki: https://github.com/AKEngels/CAST/wiki/CAST-and-Python


Notes for (new) developers:

On coding CAST:
- Try to make your code understandible, use variable names that reflect the usage and purpose and so on.

On using the git repositry:
- Make commits atomic (as small as possible). This helps tracking bugs.
- Write useful commit messages (see commit history for examples)
- If your desired commit breaks backwards compatability or simply does not compile: COMMIT TO A BRANCH!
- The master branch should be the current "nightly" build and should always compile.
- Use issues and also resolve some issues if you have time. :)
- Don't hesitate to comment on commits or issues if you have questions or suggestions
