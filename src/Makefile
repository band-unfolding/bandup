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
GIT_HASH_LATEST_COMMIT := $(shell git log --format='%H' -n 1 2> /dev/null)
GIT_USED_BRANCH := $(shell git rev-parse --abbrev-ref HEAD 2> /dev/null)
GIT_AUTHOR_DATE := $(shell git log --format='%aI' -n 1 2> /dev/null)
GIT_TAG := $(shell git describe --tags --dirty 2> /dev/null) 
# The compiler
FC = ifort
# Required flags: 
# "-assume byterecl" is required for the read_wavecar_mod module (ifort)
FCFLAGS_for_ifort = -assume byterecl 
# gfortran with GNU extensions (default, but I prefer to set it explicitly)
FCFLAGS_for_gfortran = -std=gnu -pedantic
# Shared memory parallelization
FCFLAGS_for_ifort += -openmp -parallel
FCFLAGS_for_gfortran += -fopenmp
# Flags for maximum performance, portability, etc. Comment (or modify) as necessary.
FCFLAGS_for_ifort += -O2 -no-ipo # ipo was causing a few problems.
FCFLAGS_for_gfortran += -O2
# Enable/disable warnings and remarks
FCFLAGS_for_ifort += -diag-disable remark
FCFLAGS_for_gfortran += -w
# Setting correct compiler options and enabling preprocessor
ifeq ($(FC),ifort)
  FCFLAGS=$(FCFLAGS_for_ifort)
  FPP = fpp 
else ifeq ($(FC),gfortran)
  FCFLAGS=$(FCFLAGS_for_gfortran)
  FPP = cpp 
endif
# Preprocessor stuff
ifneq ($(GIT_TAG), )
    FPPFLAGS += -D__COMMIT_TAG__=\"$(strip $(GIT_TAG))\"
    FPPFLAGS += -D__COMMIT_HASH__=\"$(strip $(GIT_HASH_LATEST_COMMIT))\"
    FPPFLAGS += -D__COMMIT_DATE__=\"$(strip $(GIT_AUTHOR_DATE))\"
    FPPFLAGS += -D__USED_BRANCH__=\"$(strip $(GIT_USED_BRANCH))\"
endif

# Aliases
current_dir = $(shell pwd)
BandUP_root_dir = $(shell dirname $(current_dir))
doc_dir = $(BandUP_root_dir)/doc
source_doc_dir = $(doc_dir)/source_code
BandUP_src_folder = $(current_dir)
BandUP_python_dir = $(BandUP_src_folder)/python_interface
BandUP_external_mods_dir = $(BandUP_src_folder)/external
cla_folder = $(BandUP_src_folder)/external/claf90
spglib_folder = $(BandUP_src_folder)/external/spglib-1.5.2
vasp_related_modules_folder = $(BandUP_src_folder)/vasp_related_modules
qe_related_modules_folder = $(BandUP_src_folder)/qe_related_modules
abinit_related_modules_folder = $(BandUP_src_folder)/abinit_related_modules
castep_related_modules_folder = $(BandUP_src_folder)/castep_related_modules
check2xsf_modules_dir = $(BandUP_src_folder)/external/check2xsf2_modules_for_BandUP
espresso_folder = $(BandUP_src_folder)/external/espresso-5.1_modules_for_BandUP
QE_modules=$(espresso_folder)/Modules
QE_iotk=$(espresso_folder)/iotk/src
QE_flib=$(espresso_folder)/flib
QE_clib=$(espresso_folder)/clib
# Checking whether CASTEP interface dir exists
ifneq ($(wildcard $(castep_related_modules_folder)), )
    ifneq ($(wildcard $(check2xsf_modules_dir)), )
        CASTEP_SUPPORT = true
        # More preprocessor stuff
        FPPFLAGS += -D__CASTEP_SUPPORT__=\"$(CASTEP_SUPPORT)\"
        FPPFLAGS += -D__CASTEP_SUPPORT__=\"$(CASTEP_SUPPORT)\"
    endif
endif

# libraries needed for linking
LIBS = $(spglib_folder)/src/.libs/libsymspg.a
LIBS += $(QE_modules)/libqemod.a $(QE_flib)/ptools.a $(QE_flib)/flib.a \
        $(QE_clib)/clib.a $(QE_iotk)/libiotk.a
# Additional module search paths
FCFLAGS += -I$(espresso_folder)/Modules -I$(espresso_folder)/iotk/src
ifneq ($(CASTEP_SUPPORT), )
    CFLAGS += -I$(castep_related_modules_folder)
endif

# Folders containing needed sources
vpath %.f90 $(cla_folder)
vpath %.f90 $(spglib_folder)/example
vpath %.f90 $(BandUP_src_folder)
vpath %.f90 $(vasp_related_modules_folder)
vpath %.f90 $(qe_related_modules_folder)
vpath %.f90 $(abinit_related_modules_folder)
ifneq ($(CASTEP_SUPPORT), )
    vpath %.f90 $(castep_related_modules_folder)
    vpath %.c $(castep_related_modules_folder)
    vpath %.c $(check2xsf_modules_dir)
endif

# List of executables to be built within the package
BANDUP_SRC = main_BandUP.f90 
BANDUP_GET_SCKPTS_SRC = main_get_SCKPTS_pre_BandUP.f90
BANDUP_SRC_OBJ = $(patsubst %.f90, %.o, $(BANDUP_SRC))
BANDUP_GET_SCKPTS_SRC_OBJ = $(patsubst %.f90, %.o, $(BANDUP_GET_SCKPTS_SRC))
BANDUP_EXE = BandUP.x
BANDUP_GET_SCKPTS_EXE = get_SCKPTS_pre_BandUP.x
PYTHON_INTERFACE_SYMLINK = $(BandUP_root_dir)/bandup
# First rule: "make" will build this rule if no target is specified
bandup: all
all: $(BANDUP_GET_SCKPTS_EXE) $(BANDUP_EXE) $(PYTHON_INTERFACE_SYMLINK)

# Objects
ifneq ($(CASTEP_SUPPORT), )
    C2XSF_C_OBJS = check2xsf_mod_for_bandup.o check_read_for_bandup.o \
                   periodic_table_c2xsf.o basis_c2xsf.o super_c2xsf.o gpfa_c2xsf.o \
                   interpolate_c2xsf.o line_write_c2xsf.o ident_sym_c2xsf.o \
                   primitive_c2xsf.o ksym_c2xsf.o sort_atoms_c2xsf.o
    CASTEP_READ_WF_OBJ = read_castep_wavefunctions_mod.o
endif
OBJ := cla.o kinds.o constants_and_types_mod.o cla_wrappers_mod.o spglib_f08.o \
       math_mod.o band_unfolding_routines_mod.o general_io_mod.o read_wavecar_mod.o \
       read_vasp_files_mod.o read_qe_wavefunctions_mod.o \
       read_abinit_wavefunctions_mod.o $(CASTEP_READ_WF_OBJ) \
       write_vasp_files_mod.o strings_mod.o io_routines_mod.o $(C2XSF_C_OBJS)

# Rules to make the executables individually
main: $(BANDUP_EXE)
prebandup: $(BANDUP_GET_SCKPTS_EXE)
python: $(PYTHON_INTERFACE_SYMLINK)
# Rules for dependencies
cla.o: kinds.o
math_mod.o: constants_and_types_mod.o spglib_f08.o
general_io_mod.o: constants_and_types_mod.o math_mod.o
cla_wrappers_mod.o: general_io_mod.o kinds.o cla.o constants_and_types_mod.o
read_vasp_files_mod.o: constants_and_types_mod.o math_mod.o general_io_mod.o \
                       strings_mod.o
write_vasp_files_mod.o: constants_and_types_mod.o math_mod.o general_io_mod.o \
                        strings_mod.o
read_wavecar_mod.o: constants_and_types_mod.o math_mod.o general_io_mod.o
read_qe_wavefunctions_mod.o: constants_and_types_mod.o cla_wrappers_mod.o math_mod.o \
                             general_io_mod.o
read_abinit_wavefunctions_mod.o: constants_and_types_mod.o math_mod.o general_io_mod.o
read_castep_wavefunctions_mod.o: constants_and_types_mod.o
io_routines_mod.o: constants_and_types_mod.o math_mod.o cla_wrappers_mod.o \
                   strings_mod.o general_io_mod.o read_vasp_files_mod.o \
                   write_vasp_files_mod.o read_wavecar_mod.o \
                   read_qe_wavefunctions_mod.o read_abinit_wavefunctions_mod.o \
                   $(CASTEP_READ_WF_OBJ)
band_unfolding_routines_mod.o: constants_and_types_mod.o cla_wrappers_mod.o math_mod.o \
                               general_io_mod.o io_routines_mod.o
$(BANDUP_SRC_OBJ) $(BANDUP_GET_SCKPTS_SRC_OBJ): $(OBJ) 
# ======================================================================
# General rules
# ======================================================================
$(BANDUP_EXE): $(BANDUP_SRC_OBJ) $(OBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

$(BANDUP_GET_SCKPTS_EXE): $(BANDUP_GET_SCKPTS_SRC_OBJ) $(OBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

$(PYTHON_INTERFACE_SYMLINK): | prebandup main
	ln -s $(BandUP_python_dir)/bandup $@

%.o: %.f90
	$(FC) $(FCFLAGS) -$(FPP) $(FPPFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -$(FPP) $(FPPFLAGS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -c $<

# Utility targets
.PHONY: all sourcedoc clean veryclean print

sourcedoc:
	export BANDUPDIR=$(BandUP_root_dir);\
	export DOXIGENOUTDIR=$(source_doc_dir);\
	mkdir -p $(source_doc_dir);\
	doxygen
clean:
	rm -f *.o *.mod *.MOD
	find $(BandUP_python_dir) -name "*.pyc" -type f -delete
	find $(BandUP_python_dir) -name "__pycache__" -type d -delete
veryclean: clean
	rm -f *~ $(BANDUP_EXE) $(BANDUP_GET_SCKPTS_EXE) $(PYTHON_INTERFACE_SYMLINK)
	rm -rf $(doc_dir)
print-%:
	@echo '$*=$($*)'
