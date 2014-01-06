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
# Required flags: "-assume byterecl" is required for the compilation of the read_wavecar_mod module
FCFLAGS = -assume byterecl 
# Shared memory parallelization
FCFLAGS += -openmp -parallel
# Flags for maximum performance, portability, etc. Comment (or modify) as necessary.
#FCFLAGS += -static
FCFLAGS += -O2 -ipo
#FCFLAGS += -O2 -ipo -xHost
#FCFLAGS += -O2 -xAVX # Uncomment if you have an Intel processor supporting the AVX instruction set
# Enable/disable warnings and remarks
FCFLAGS += -warn unused
FCFLAGS += -diag-disable remark

# List of executables to be built within the package
SRC = main_BandUP.f90
OBJ = ${SRC:.f90=.o}
EXE = BandUP.x

# "make" builds all
all: $(EXE)

band_unfolding_routines_mod.o: math_mod.o
read_wavecar_mod.o: math_mod.o general_io_mod.o
io_routines_mod.o: math_mod.o strings_mod.o general_io_mod.o read_wavecar_mod.o
$(OBJ): math_mod.o io_routines_mod.o read_wavecar_mod.o band_unfolding_routines_mod.o
$(EXE): math_mod.o band_unfolding_routines_mod.o general_io_mod.o read_wavecar_mod.o strings_mod.o io_routines_mod.o $(OBJ)

# ======================================================================
# General rules
# ======================================================================
$(EXE): $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

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
