# Configuration file for compiling software on Ubuntu linux system
# using debugging gfortran compiler

# Set EXEDIR 
EXEDIR = ../exec

# Fortran compiler
F90 ?= gfortran
FC = gfortran

# Fortran compiler flags
F90FLAGS = -g -fbacktrace -fbounds-check -fexceptions -cpp
ifeq ($(findstring mpi,$(F90)),mpi)
   F90FLAGS += -DWITH_MPI
endif

FFLAGS = -fdefault-real-8 -fdefault-double-8 -g -fbacktrace -fbounds-check -fexceptions

FSFLAG = -g -fbacktrace -fbounds-check -fexceptions

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
