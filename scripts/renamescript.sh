#!/bin/sh
# Call this with ./renamescript.sh <files>  and it will do some in-place replacing to get stuff working. It makes a backup B file at the end so you can do

input=$@

for x in ${input}; do
    # Replace gf... notation with f
    sed -i"B" -E "s/g(f|x|y|z)[itb]/\1/gi; s/J/PETSC_i/g;s/Sinh/sinh/g;s/Sqrt/PetscSqrtScalar/g;s/Power/PetscPowScalar/g;s/(Sinh?)/Petsc\1Scalar/g;s/(Cosh?)/Petsc\1Scalar/g;s/(Tanh?)/Petsc\1Scalar/g;s/Pi/M_PI/g;s/Lx/lx/g;s/Ly/ly/g;s/source\[0\]/source1/g;s/source\[1\]/source2/g;s/Ax/ax/g;s/Ay/ay/g;s/NPeriodsx/nperiodsx/gi;s/NPeriodsY/nperiodsy/gi;s/PhaseX/phasex/gi; s/Phasey/phasey/gi" $x;
done

mkdir -p backup
mv *B backup

