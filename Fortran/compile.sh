#!/bin/sh
# Set source directory:
export SRC=..

# Set compiler (usually gfortran) and compiler flags:
export FC=gfortran
# Use -O0 to actually catch the bug in gfortran
#export FFLAGS_g="-g -Og -fcheck=bounds,array-temps" # Full debugging
export FFLAGS_g="-g -Og" # Add debugging info with some optimization
export FFLAGS_O="-Ofast" # Optimize for speed. Can try -O3

# If you have the Intel Fortran compiler
#export FC=ifort; export FFLAGS_g="-g -debug all -check all -traceback"; export FFLAGS_O="-O3"

export FFLAGS=${FFLAGS_g} # Choose _g (debug) or _O (optimize)

# If you want to limit the amount of error messages you see, use:
# ${FC} -c source.f90 2>&1 | head -10

# Compile Fortran "library":
${FC} ${FFLAGS} -c ${SRC}/FLib/Precision.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/r8lib.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/Tridiagonal.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/MatrixExponential.f90
${FC} ${FFLAGS} -c ${SRC}/FLib/MatrixInverse.f90

# Compile linear operators:
${FC} ${FFLAGS} -c ${SRC}/LinearOperator.f90
${FC} ${FFLAGS} -c ${SRC}/BandedMatrix.f90

# Make executable:
echo ${FC} ${FFLAGS} -o TestBandedMatrix.x ${SRC}/TestBandedMatrix.f90 BandedMatrix.o LinearOperator.o Tridiagonal.o MatrixExponential.o MatrixInverse.o r8lib.o Precision.o -llapack
${FC} ${FFLAGS} -o TestBandedMatrix.x ${SRC}/TestBandedMatrix.f90 BandedMatrix.o LinearOperator.o Tridiagonal.o MatrixExponential.o MatrixInverse.o r8lib.o Precision.o -llapack
echo "Compiled executable TestBandedMatrix.x"

# Compile ODE solvers:
${FC} ${FFLAGS} -c ${SRC}/ODEIntegrand.f90

echo ${FC} ${FFLAGS} -o TestODEIntegrator.x ${SRC}/TestODEIntegrator.f90 ODEIntegrand.o BandedMatrix.o LinearOperator.o Tridiagonal.o MatrixExponential.o MatrixInverse.o r8lib.o Precision.o -llapack
${FC} ${FFLAGS} -o TestODEIntegrator.x ${SRC}/TestODEIntegrator.f90 ODEIntegrand.o BandedMatrix.o LinearOperator.o Tridiagonal.o MatrixExponential.o MatrixInverse.o r8lib.o Precision.o -llapack
echo "Compiled executable TestODEIntegrator.x"

