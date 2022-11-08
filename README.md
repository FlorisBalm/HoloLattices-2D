# HoloLattices-2D

C and python code for to compute holographic lattices with two transverse dimensions, for the Gubser-Rocha EMD Lattice, along with the thermoelectric DC conductivity.

This code is functional, and includes the Stokes methods to compute DC thermoelectric conductivity

[![DOI](https://zenodo.org/badge/561700485.svg)](https://zenodo.org/badge/latestdoi/561700485)

# Instructions

To compile and run, first have PETSc installed.

From the source directory, do the following: 

`mkdir build coefs/output_inlined equations/output HeadersSecondOrder`

`python3 scripts/ExpressionCoefficients2ndOrder.py`

`cd coefs; python3 CPP_From_Coef_Inline_Trivial.py`

`cd ../equations; python3 EOM_To_CPP.py`

`bash CompileScripts/BuildEquationObjectFiles.batch`

`bash CompileScripts/BuildBinaryDouble.batch`

These files contain code specifically to compile for the Snellius supercomputer of the Netherlands, but will work fine as batch files.

# Stokes

the stokes code has an additional dependency, namely the [HolographicLattice](https://github.com/FlorisBalm/HoloLattices), see also [![DOI](https://zenodo.org/badge/504278109.svg)](https://zenodo.org/badge/latestdoi/504278109)
