
ifeq ($(USER_FC),$(null))
  FC := ifort
else
  FC := $(USER_FC)
endif

srcdir = src
blddir = bld

OBJS := $(srcdir)/system_mod.F90 $(srcdir)/lusolvec_mod.F90 $(srcdir)/adgaquad_types_mod.F90 $(srcdir)/adgaquad_mod.F90 $(srcdir)/miess_mod.F90 $(srcdir)/fractal_meanfield_mod.F90 $(srcdir)/main.F90


fractaloptics.exe: $(OBJS) 
	$(FC) -o $@  $(OBJS) 
	$(RM) -r -f $(blddir)
	mkdir $(blddir)
	mv *mod $(blddir)
clean: 
	$(RM) -f *.mod *.o fractaloptics.exe 
	$(RM) -r -f $(blddir)

