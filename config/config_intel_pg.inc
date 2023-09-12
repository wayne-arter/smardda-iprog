# Configuration file for compiling software on scientific linux system
# using optimised intel compiler

# Set EXEDIR 
EXEDIR = ../exec

# Fortran compiler
F90 ?= ifort
FC = ifort

# Fortran compiler flags
F90FLAGS = -O3 -C -pg -fpp
ifeq ($(findstring mpi,$(F90)),mpi)
   F90FLAGS += -DWITH_MPI
endif

FFLAGS = -r8 -save -O -pg

FSFLAG = -save -O -pg

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
