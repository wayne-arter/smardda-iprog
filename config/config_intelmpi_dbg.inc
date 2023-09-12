# Configuration file for compiling software on scientific linux system
# using optimised intel compiler

# Set EXEDIR 
EXEDIR = ../exec

# Fortran compiler
F90 ?= mpif90
FC = ifort

# Fortran compiler flags
F90FLAGS = debug -debug-parameters all -check all -trace -g -O -fpp -D WITH_MPI

FFLAGS = -r8 -save -O

FSFLAG = -save -O

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
