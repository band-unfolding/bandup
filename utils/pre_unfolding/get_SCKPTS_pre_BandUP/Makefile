# Adapted from: 
#     Makefile-Fortran: A Makefile for a simple Fortran project
#     Copyright (C) 2009 Davide Cesari, dcesari <at> arpa.emr.it
#     Available at http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html
#     Distributed according to the GNU GPL license.
#     This makefile works with the GNU make command, the one find on
#     GNU/Linux systems and often called gmake on non-GNU systems
# ======================================================================
# Declarations
# ======================================================================
# The compiler
FC = ifort

# Needed flags
FCFLAGS = -assume byterecl 

# flags for debugging or for maximum performance, comment (or edit) as necessary
FCFLAGS += -openmp
FCFLAGS += -O2 -no-ipo # ipo was causing a few problems.

# Aliases
current_dir = $(shell pwd)
BandUP_src_folder = $(current_dir)/../../../src/
spglib_folder = $(BandUP_src_folder)/external/spglib-1.5.2

# libraries needed for linking
LIBS = $(spglib_folder)/src/.libs/libsymspg.a

# Folders containing needed sources
vpath %.f90 $(spglib_folder)/example
vpath %.f90 $(BandUP_src_folder)
# List of executables to be built within the package
SRC = main_get_SCKPTS_pre_BandUP.f90
OBJ = ${SRC:.f90=.o}
EXE = get_SCKPTS_pre_BandUP.x

# "make" builds all
all: $(EXE)

math_mod.o: spglib_f08.o
band_unfolding_routines_mod.o: math_mod.o
read_wavecar_mod.o: math_mod.o general_io_mod.o
io_routines_mod.o: math_mod.o strings_mod.o general_io_mod.o read_wavecar_mod.o
$(OBJ): spglib_f08.o math_mod.o io_routines_mod.o read_wavecar_mod.o band_unfolding_routines_mod.o
$(EXE): spglib_f08.o math_mod.o band_unfolding_routines_mod.o general_io_mod.o read_wavecar_mod.o strings_mod.o io_routines_mod.o $(OBJ)
# ======================================================================
# General rules
# ======================================================================
$(EXE): $(OBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(EXE)

