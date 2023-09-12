# Configuration file for compiling software on scientific linux system
# using debugging intel compiler

# Set EXEDIR 
EXEDIR = ../exec

# Fortran compiler
F90 ?= ifort
FC = ifort

# Fortran compiler flags
#F90FLAGS = -O -std03 -unroll -fp-model source -diag-disable 8291,8577
F90FLAGS = -debug -debug-parameters all -trace -fpe0 -fp-stack-check -check -fpp
ifeq ($(findstring mpi,$(F90)),mpi)
   F90FLAGS += -DWITH_MPI
endif

FFLAGS = -r8 -save -debug -debug-parameters all -trace -fpe0 -fp-stack-check -check

FSFLAG = -save -debug -debug-parameters all -trace -fpe0 -fp-stack-check -check

# Debug option
DBG = g

# command to remove files
RM = /bin/rm -f

# extra cleanup command
EXTRA_CLEAN_CMD = :

# Extra files to clean
EXTRA_CLEAN_FILES = *.T

# Module extension
MOD=.mod
