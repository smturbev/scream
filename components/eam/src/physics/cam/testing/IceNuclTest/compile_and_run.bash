#!/bin/bash

module purge
module load intel

pushd obj
rm -f *.o *.mod
ifort -c ../src/shr_kind_mod.F90 
ifort -c ../src/phys_utils_minimal.f90
ifort -c ../src/IceNucleationOriginal.f90
popd
ifort -Iobj -o TestIceNucleation src/TestIceNucleation.f90 obj/*.o

nice ./TestIceNucleation
