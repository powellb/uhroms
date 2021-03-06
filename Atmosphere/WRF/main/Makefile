#

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f

MODULES = 

OBJS    =

LIBPATHS = 

include ../configure.wrf

$(SOLVER)_wrf : wrf.o
	$(RANLIB) libwrflib.a
	$(AR) libWRF.a wrf.o libwrflib.a
	$(RANLIB) libWRF.a
#	$(LD) -o wrf.exe $(LDFLAGS) wrf.o libwrflib.a $(LIB)
#jcw

$(SOLVER)_wrf_ESMFApp : wrf_ESMFMod.o wrf_ESMFApp.o wrf_SST_ESMF.o
	$(RANLIB) libwrflib.a
	$(LD) -o wrf_ESMFApp.exe $(LDFLAGS) wrf_ESMFApp.o wrf_ESMFMod.o libwrflib.a $(LIB)
	$(LD) -o wrf_SST_ESMF.exe $(LDFLAGS) wrf_SST_ESMF.o wrf_ESMFMod.o libwrflib.a $(LIB)

$(SOLVER)_ideal : module_initialize ideal.o
	$(RANLIB) libwrflib.a
	$(LD) -o ideal.exe $(LDFLAGS) ideal.o ../dyn_$(SOLVER)/module_initialize_$(IDEAL_CASE).o libwrflib.a $(LIB)

$(SOLVER)_real : module_initialize real_$(SOLVER).o ndown_$(SOLVER).o nup_$(SOLVER).o
	$(RANLIB) libwrflib.a
	$(LD) -o real.exe $(LDFLAGS) real_$(SOLVER).o ../dyn_$(SOLVER)/module_initialize_$(IDEAL_CASE).o libwrflib.a $(LIB)
	$(LD) -o ndown.exe $(LDFLAGS) ndown_$(SOLVER).o  ../dyn_$(SOLVER)/module_initialize_$(IDEAL_CASE).o libwrflib.a $(LIB)
	$(LD) -o nup.exe $(LDFLAGS) nup_$(SOLVER).o  ../dyn_$(SOLVER)/module_initialize_$(IDEAL_CASE).o libwrflib.a $(LIB)

convert_em : convert_em.o
	$(RANLIB) libwrflib.a
	$(LD) -o convert_em.exe $(LDFLAGS) convert_em.o libwrflib.a $(LIB)

convert_nmm : convert_nmm.o
	$(RANLIB) libwrflib.a
	$(FC) -o convert_nmm.exe $(LDFLAGS) convert_nmm.o libwrflib.a $(LIB)

real_nmm : real_nmm.o
	( cd ../dyn_nmm ;  $(MAKE) module_initialize_real.o )
	$(RANLIB) libwrflib.a
	$(FC) -o real_nmm.exe $(LDFLAGS) real_nmm.o ../dyn_nmm/module_initialize_real.o ../share/module_optional_si_input.o ../share/input_wrf.o ../share/module_io_domain.o libwrflib.a $(LIB)

module_initialize :
	( cd ../dyn_$(SOLVER) ;  $(MAKE) module_initialize_$(IDEAL_CASE).o )

## prevent real being compiled for OMP -- only for regtesting
#$(SOLVER)_real : module_initialize real_$(SOLVER).o
#	$(RANLIB) libwrflib.a
#	if [ -z "$(OMP)" ] ; then $(FC) -o real.exe $(LDFLAGS) real_$(SOLVER).o ../dyn_$(SOLVER)/module_initialize_$(IDEAL_CASE).o libwrflib.a $(LIB) ; fi
#
## prevent module_initialize being compiled for OMP --remove after IBM debugging
#module_initialize :
#	if [ -z "$(OMP)" ] ; then ( cd ../dyn_$(SOLVER) ;  $(MAKE) module_initialize_$(IDEAL_CASE).o ) ; fi
# end of regtest changes

clean:
	@ echo 'use the clean script'

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

convert_nmm.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../share/module_bc.o \
	../share/module_io_domain.o \
	$(ESMF_MOD_DEPENDENCE)

convert_em.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../share/module_bc.o \
	../share/module_io_domain.o \
	$(ESMF_MOD_DEPENDENCE)

ideal.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../share/module_io_domain.o \
	../dyn_$(SOLVER)/$(CASE_MODULE) \
	$(ESMF_MOD_DEPENDENCE)

ndown_em.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../frame/module_wrf_error.o \
	../frame/module_integrate.o \
	../share/module_bc.o \
	../share/module_io_domain.o \
	../share/module_get_file_names.o \
	../share/module_soil_pre.o \
	../dyn_em/module_initialize_real.o \
	../dyn_em/module_big_step_utilities_em.o \
	$(ESMF_MOD_DEPENDENCE)

nup_em.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../frame/module_wrf_error.o \
	../frame/module_integrate.o \
	../share/module_bc.o \
	../share/module_io_domain.o \
	../share/module_get_file_names.o \
	../share/module_soil_pre.o \
	../dyn_em/module_initialize_real.o \
	../dyn_em/module_big_step_utilities_em.o \
	$(ESMF_MOD_DEPENDENCE)

# this already built above :../dyn_em/module_initialize.real.o \
real_em.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../dyn_em/module_si_io_em.o \
	../dyn_em/module_big_step_utilities_em.o \
	../share/module_io_domain.o \
	../share/module_date_time.o \
	../share/module_optional_si_input.o \
	../share/module_bc_time_utilities.o \
	$(ESMF_MOD_DEPENDENCE)
#	../chem/module_input_chem_data.o \
#	../chem/module_input_chem_bioemiss.o \

wrf.o:  ../share/module_wrf_top.o

wrf_ESMFMod.o:  ../share/module_wrf_top.o

wrf_ESMFApp.o:  wrf_ESMFMod.o

wrf_SST_ESMF.o:  wrf_ESMFMod.o

# DO NOT DELETE
