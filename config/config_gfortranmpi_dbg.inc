# Configuration file for compiling software on Ubuntu linux system
# using optimised gfortran compiler

# Set EXEDIR 
EXEDIR = ../exec

# Fortran compiler
F90 = mpifort
FC = mpifort

# Fortran compiler flags
F90FLAGS = -fexceptions -fbacktrace --check=all -g -O -cpp -D WITH_MPI

FFLAGS = -fdefault-real-8 -fdefault-double-8 -O 

FSFLAG = -O

# Debug option
DBG =

# command to remove files
RM = /bin/rm -f

# extra cleanup command
EXTRA_CLEAN_CMD = :

# Extra files to clean
EXTRA_CLEAN_FILES = *.T

# Module extension
MOD=.mod
