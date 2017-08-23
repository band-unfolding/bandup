# Aliases
rootdir = $(shell pwd)
cla = claf90
spglib = spglib-1.5.2
check2xsf = check2xsf2_modules_for_BandUP
espresso = espresso-5.1_modules_for_BandUP
spglib_makefile = $(spglib)/Makefile
espresso_makefile = $(espresso)/Makefile

# Rules
external_modules: all
all: $(cla) $(spglib) $(check2xsf) $(espresso)

$(cla):
	tar -xvzf $@.tgz
$(spglib_makefile):
	tar -xvzf $(spglib).tgz
	cd $(spglib);\
	./configure FC=$(FC) CC=$(CC) CFLAGS=$(OMP_CFLAG) FCFLAGS=$(OMP_FCFLAG);\
	cd $(rootdir)
$(spglib): $(spglib_makefile)
	cd $@;\
	$(MAKE) -j --max-load 2.5;\
	cd $(rootdir)
$(check2xsf):
ifneq ($(wildcard $(check2xsf).tgz), )
	tar -xvzf $@.tgz
endif
$(espresso_makefile):
	tar -xvzf $(espresso).tgz
	cd $(espresso);\
	./configure --disable-parallel F77=$(FC) F90=$(FC) FC=$(FC) CC=$(CC);\
	cd $(rootdir)
$(espresso): $(espresso_makefile)
	cd $@;\
	$(MAKE) -j --max-load 2.5 bandup;\
	cd $(rootdir)

# Utility targets
.PHONY: all external_modules clean veryclean print

clean:
	find $(rootdir) -name "*.o" -type f -delete
	find $(rootdir) -name "*.mod" -type f -delete
	find $(rootdir) -name "*.MOD" -type f -delete
veryclean: clean
	rm -rf $(cla) $(spglib) $(check2xsf) $(espresso)
print-%:
	@echo '$*=$($*)'
